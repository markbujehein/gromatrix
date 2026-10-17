"""[G] Run the MD protocol of every built system: EM -> NVT -> NPT -> production.

    gromatrix-run --batch-dir B [--only ID ...] [--threads 8] [--gpu-id 0] [--deffnm md] [--maxh H]
    gromatrix-run --batch-dir B --background            # detach; progress in B/07-run-logs/
    gromatrix-run --batch-dir B --slurm --time 24:00:00 [--partition gpu --gpus-per-task 1]

Inputs are step F's 06-systems/<id>/ (topol.top, index.ndx, mdp/, solv_ions.gro).
The phases chain exactly like run.sh: -r gives POSRES reference coordinates,
-t carries velocities and box/thermostat state forward. -maxwarn is never used.

Resume rules, checked per phase before anything runs:
  * done     the .gro exists, the newest log segment finished normally, and the
             .cpt the next phase continues from exists      -> gated, then skipped
  * partial  a .cpt exists                                    -> mdrun -cpi continues it
  * new      otherwise                                        -> grompp, then mdrun
Once a phase actually runs, every later phase starts fresh, so no stale output
from an earlier attempt is ever chained. (Step F rebuilds a changed system from
scratch, which also discards its MD output; it refuses to without --force.)

Gates, applied to every phase including skipped ones (so an EM left by step F's
smoke test is checked too): EM must converge (Fmax < emtol, Epot not positive);
NVT/NPT must finish with no fatal error and no LINCS warning in any log segment;
production must finish without a fatal error (LINCS warnings are reported). A
phase that fails its gate keeps <phase>.failed.{gro,log} as evidence and is
cleared, so a rerun starts it fresh. A failed system stops; the others go on.

Time and signals: --maxh is one budget for the whole invocation; each mdrun gets
what is left, so a Slurm task stops cleanly with a checkpoint before its
walltime. SIGTERM / Ctrl-C make mdrun checkpoint and stop the batch, also when
they arrive between phases.

Output: progress lines {"step": "run", "system", "phase", "status", "progress_pct",
"overall_pct"} while mdrun runs, then one summary line {"ok": ...}. State per
system goes to 06-systems/<id>/run.json, the batch run to B/07-run-logs/run.json.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import signal
import socket
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import checks, layout
from . import mdp as mdpmod
from .forcefields import load_lock
from .gmx import ToolError, gmx_cmd, gmxlib_env, locate_ff_dir, run, tool_env
from .net import atomic_write


@dataclass(frozen=True)
class Phase:
    name: str  # also the mdp file stem
    coords: str | None  # -c input; None = the previous phase's output
    restrained: bool  # -r reference coordinates (POSRES)
    continues: bool  # -t previous checkpoint (velocities, box, thermostat state)


PHASES = (
    Phase("em", "solv_ions.gro", True, False),
    Phase("nvt", None, True, False),
    Phase("npt", None, True, True),
    Phase("md", None, False, True),
)
# Phases whose checkpoint the next phase reads with -t: "done" also needs their .cpt.
FEEDS_CPT = {p.name for p, nxt in zip(PHASES, PHASES[1:]) if nxt.continues}
DYNAMICS = ("nvt", "npt", "md")
STEP_RE = re.compile(r"^\s+Step\s+Time\s*\n\s+(\d+)\s+", re.M)
MIN_MAXH = 0.01  # h; below this, starting another mdrun is pointless


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def emit(rec: dict) -> None:
    print(json.dumps(rec), flush=True)


def load_json(path: Path) -> dict:
    """run.json may be missing or cut short by a kill; neither may break a rerun."""
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return {}


# ------------------------------------------------------------------ phase state
def deffnm(phase: Phase, prod: str) -> str:
    return prod if phase.name == "md" else phase.name


def last_step(text: str) -> int | None:
    m = STEP_RE.findall(text)
    return int(m[-1]) if m else None


def phase_state(wd: Path, phase: Phase, name: str) -> str:
    scan = checks.scan_md_log(wd / f"{name}.log")
    cpt = (wd / f"{name}.cpt").exists()
    if (wd / f"{name}.gro").exists() and scan["finished"] and not scan["fatal"]:
        if cpt or phase.name not in FEEDS_CPT:
            return "done"
    return "partial" if cpt else "new"


def clear_phase(wd: Path, name: str, keep_tpr: bool = False) -> None:
    files = [f"{name}.{e}" for e in ("cpt", "gro", "log", "edr", "xtc", "trr", "mdout.mdp")] + [f"{name}_prev.cpt"]
    if not keep_tpr:
        files.append(f"{name}.tpr")
    for f in files:
        (wd / f).unlink(missing_ok=True)


def nsteps_of(wd: Path, phase: Phase) -> int:
    try:
        return int(float(mdpmod.read_mdp(wd / "mdp" / f"{phase.name}.mdp").get("nsteps", 0)))
    except (OSError, ValueError):
        return 0


# ------------------------------------------------------------------ execution
class Interrupted(RuntimeError):
    pass


class Runner:
    """Runs the protocol for one system directory."""

    def __init__(self, wd: Path, a: argparse.Namespace, env: dict, deadline: float | None):
        self.wd, self.a, self.env, self.deadline = wd, a, env, deadline
        self.state_path = wd / "run.json"
        self.log = wd / "run.log"
        self.child: subprocess.Popen | None = None
        self.stop = False
        self.nsteps = {p.name: nsteps_of(wd, p) for p in PHASES}
        self.total = sum(self.nsteps[k] for k in DYNAMICS) or 1

    def g(self, *sub: str) -> list[str]:
        return gmx_cmd(self.a.gmx, *sub)

    def request_stop(self) -> None:
        """Signal handler side: stop at the next opportunity; mdrun checkpoints on SIGTERM."""
        self.stop = True
        if self.child and self.child.poll() is None:
            self.child.send_signal(signal.SIGTERM)

    def check_stop(self) -> None:
        if self.stop:
            raise Interrupted("stopped by a signal; rerun gromatrix-run to resume")

    def overall(self, name: str, step: int) -> float:
        done = 0
        for p in PHASES:
            if p.name == name:
                break
            done += self.nsteps[p.name] if p.name in DYNAMICS else 0
        return round(100 * (done + (step if name in DYNAMICS else 0)) / self.total, 1)

    def grompp(self, phase: Phase, prev: str | None) -> None:
        self.check_stop()
        name = deffnm(phase, self.a.deffnm)
        c = phase.coords or f"{prev}.gro"
        cmd = ["grompp", "-f", f"mdp/{phase.name}.mdp", "-c", c, "-p", "topol.top", "-n", "index.ndx"]
        if phase.restrained:
            cmd += ["-r", c]
        if phase.continues:
            cmd += ["-t", f"{prev}.cpt"]
        run(self.g(*cmd, "-o", f"{name}.tpr", "-po", f"{name}.mdout.mdp"), self.wd, self.log, env_extra=self.env)

    def maxh(self) -> float | None:
        if self.deadline is None:
            return None
        left = (self.deadline - time.monotonic()) / 3600
        if left < MIN_MAXH:
            raise Interrupted("--maxh time budget used up; rerun (or resubmit) to resume")
        return round(left, 4)

    def mdrun_cmd(self, name: str, resume: bool) -> list[str]:
        a = self.a
        cmd = self.g("mdrun", "-deffnm", name)
        if resume:
            cmd += ["-cpi", f"{name}.cpt"]
        if a.threads:
            cmd += ["-nt", str(a.threads)]
        if a.gpu_id is not None:
            cmd += ["-gpu_id", str(a.gpu_id)]
        maxh = self.maxh()
        if maxh is not None:
            cmd += ["-maxh", str(maxh)]
        return cmd + shlex.split(a.mdrun_args)

    def mdrun(self, phase: Phase, resume: bool) -> None:
        self.check_stop()
        name = deffnm(phase, self.a.deffnm)
        phase_log = self.wd / f"{name}.log"
        cmd = self.mdrun_cmd(name, resume)
        env = tool_env(cmd, self.env)
        start = len(checks.read_log(phase_log)) if resume else 0  # -cpi appends; judge only the new text
        n = self.nsteps[phase.name]
        with open(self.log, "a") as fh:
            fh.write(f"\n$ {shlex.join(cmd)}\n")
            fh.flush()
            self.child = subprocess.Popen(cmd, cwd=self.wd, stdout=fh, stderr=subprocess.STDOUT, env=env)
            while True:
                try:
                    rc = self.child.wait(timeout=self.a.poll)
                    break
                except KeyboardInterrupt:
                    self.request_stop()
                except subprocess.TimeoutExpired:
                    step = min(last_step(checks.last_segment(checks.read_log(phase_log))) or 0, n)
                    emit(
                        {
                            "step": "run",
                            "system": self.wd.name,
                            "phase": phase.name,
                            "status": "running",
                            "progress_pct": round(100 * step / n, 1) if n else None,
                            "overall_pct": self.overall(phase.name, step),
                        }
                    )
        self.child = None
        new_text = checks.read_log(phase_log)[start:]
        if self.stop or rc in (-signal.SIGTERM, -signal.SIGINT):
            raise Interrupted(f"{phase.name} stopped by a signal; rerun gromatrix-run to resume")
        if rc != 0:  # checked before the log: a crash must never look like a resumable stop
            tail = "\n".join(self.log.read_text(errors="replace").strip().splitlines()[-20:])
            raise ToolError(f"mdrun {phase.name} exited with {rc} (see {self.log.name}):\n{tail}")
        if checks.STOPPED_EARLY.search(new_text):
            raise Interrupted(f"{phase.name} stopped before nsteps (-maxh or signal); rerun to resume")

    def gate(self, phase: Phase) -> list[str]:
        """Raise CheckError on a hard failure; return warnings."""
        name = deffnm(phase, self.a.deffnm)
        log = self.wd / f"{name}.log"
        if phase.name == "em":
            tol = mdpmod.emtol(self.wd / "mdp")
            em = checks.parse_em_log(log, tol)
            if not em["converged"]:
                raise checks.CheckError(
                    f"EM did not converge (Fmax={em['fmax']}, Epot={em['epot']}, emtol={tol}); "
                    "look for clashes or a bad ligand pose"
                )
            return []
        scan = checks.scan_md_log(log)
        if scan["fatal"] or not scan["finished"]:
            raise checks.CheckError(f"{phase.name}: mdrun did not finish cleanly (see {log.name})")
        if scan["lincs_warnings"] and phase.name != "md":
            raise checks.CheckError(
                f"{phase.name}: {scan['lincs_warnings']} LINCS warnings during restrained equilibration"
            )
        return [f"{phase.name}: {scan['lincs_warnings']} LINCS warnings"] if scan["lincs_warnings"] else []

    def quarantine(self, name: str) -> None:
        """Keep the evidence of a gated-out phase, then clear it so a rerun starts it fresh."""
        for ext in ("gro", "log"):
            src = self.wd / f"{name}.{ext}"
            if src.exists():
                src.replace(self.wd / f"{name}.failed.{ext}")
        clear_phase(self.wd, name)

    def save(self, state: dict) -> None:
        state["updated"] = now()
        atomic_write(self.state_path, json.dumps(state, indent=2))

    def execute(self, inputs_sha: str | None) -> dict:
        state = {
            "id": self.wd.name,
            "status": "running",
            "inputs_sha": inputs_sha,
            "deffnm": self.a.deffnm,
            "host": socket.gethostname(),
            "pid": os.getpid(),
            "started": now(),
            "phases": {},
            "warnings": [],
        }
        self.save(state)
        prev, ran = None, False
        try:
            for phase in PHASES:
                name = deffnm(phase, self.a.deffnm)
                st = "new" if ran else phase_state(self.wd, phase, name)
                if st == "new":
                    # A fresh start never inherits a stale log or checkpoint. Step F's em.tpr
                    # passed a zero-warning grompp, so it is reused.
                    clear_phase(self.wd, name, keep_tpr=phase.name == "em" and not ran)
                    if not (self.wd / f"{name}.tpr").exists():
                        self.grompp(phase, prev)
                if st != "done":
                    state["phases"][phase.name] = {"status": "running", "resumed": st == "partial", "started": now()}
                    self.save(state)
                    emit(
                        {
                            "step": "run",
                            "system": self.wd.name,
                            "phase": phase.name,
                            "status": "resuming" if st == "partial" else "started",
                            "progress_pct": 0.0,
                            "overall_pct": self.overall(phase.name, 0),
                        }
                    )
                    self.mdrun(phase, resume=st == "partial")
                    ran = True
                try:
                    state["warnings"] += self.gate(phase)
                except checks.CheckError:
                    self.quarantine(name)
                    raise
                state["phases"][phase.name] = {
                    "status": "done",
                    "finished": now(),
                    **({"skipped": True} if st == "done" else {}),
                }
                prev = name
                emit(
                    {
                        "step": "run",
                        "system": self.wd.name,
                        "phase": phase.name,
                        "status": "skipped" if st == "done" else "done",
                        "progress_pct": 100.0,
                        "overall_pct": self.overall(phase.name, self.nsteps[phase.name]),
                    }
                )
            state["status"] = "done"
            state["trajectory"] = f"{self.a.deffnm}.xtc"
        except Interrupted as e:
            state["status"], state["reason"] = "interrupted", str(e)
        except KeyboardInterrupt:  # Ctrl-C outside mdrun (e.g. during grompp)
            self.request_stop()
            state["status"], state["reason"] = "interrupted", "interrupted by the user; rerun to resume"
        except Exception as e:  # ToolError, CheckError: this system stops, the batch goes on
            state["status"], state["reason"] = "failed", f"{type(e).__name__}: {str(e)[:600]}"
        for p in state["phases"].values():
            if p.get("status") == "running":
                p["status"] = state["status"]
        state.pop("pid", None)
        self.save(state)
        return state


# ------------------------------------------------------------------ selection
def runnable(b: Path, only: list[str] | None, force: bool) -> tuple[list[tuple[Path, dict]], list[dict]]:
    """Built systems that may run. A failed smoke test blocks a system unless --force."""
    picked, skipped = [], []
    for stamp_p in sorted((b / layout.SYSTEMS).glob("*/system.json")):
        wd = stamp_p.parent
        if only and wd.name not in only:
            continue
        stamp = load_json(stamp_p)
        if stamp.get("status") != "ok":
            skipped.append({"id": wd.name, "reason": "not built"})
        elif stamp.get("smoke", {}).get("ok") is False and not force:
            skipped.append({"id": wd.name, "reason": "smoke test failed; fix it or pass --force"})
        else:
            picked.append((wd, stamp))
    if only:
        known = {wd.name for wd, _ in picked} | {s["id"] for s in skipped}
        skipped += [{"id": i, "reason": "no such system"} for i in only if i not in known]
    return picked, skipped


def owner_alive(st: dict) -> bool | None:
    """Is the runner recorded in a 'running' run.json alive? None = on another host (unknown)."""
    if st.get("status") != "running" or not st.get("pid"):
        return False
    if st.get("host") != socket.gethostname():
        return None
    try:
        os.kill(int(st["pid"]), 0)
    except (OSError, ValueError):
        return False
    return True


def busy(wd: Path) -> str | None:
    """Best-effort guard against two runners on one system (not a lock: a check, then a write)."""
    st = load_json(wd / "run.json")
    alive = owner_alive(st)
    if alive is None:
        return f"running on {st.get('host')} (pid {st.get('pid')})"
    if alive and int(st["pid"]) != os.getpid():
        return f"running here (pid {st['pid']})"
    return None


def ff_env(b: Path, gmx: str) -> dict:
    try:
        return gmxlib_env(locate_ff_dir(load_lock(b).gmx_name, gmx, [b.resolve()]))
    except (ToolError, OSError, ValueError, KeyError):
        return {}  # systems carry an ff symlink (step F), so grompp still resolves the includes


# ------------------------------------------------------------------ slurm / background
def _hours(walltime: str) -> float:
    """Slurm --time -> hours. Formats: M, M:S, H:M:S, D-H, D-H:M, D-H:M:S."""
    days, sep, rest = walltime.partition("-") if "-" in walltime else ("0", "", walltime)
    f = [float(x) for x in rest.split(":")]
    if sep:  # after "D-" the fields are hours[:minutes[:seconds]]
        f += [0.0] * (3 - len(f))
        h, m, s = f
    elif len(f) == 3:
        h, m, s = f
    else:  # minutes[:seconds]
        h, m, s = 0.0, f[0], f[1] if len(f) > 1 else 0.0
    return 24 * float(days) + h + m / 60 + s / 3600


SBATCH = """#!/usr/bin/env bash
# Generated by gromatrix (gromatrix-run --slurm). One array task per system.
# The whole task shares one -maxh budget just under the walltime, so a task that
# runs out of time stops cleanly at a checkpoint: resubmit this script and every
# task resumes where it stopped.
#SBATCH --job-name={job}
#SBATCH --array=0-{last}{throttle}
#SBATCH --time={time}
#SBATCH --cpus-per-task={cpus}
#SBATCH --output={logdir}/%x_%A_%a.out
#SBATCH --signal=B:TERM@300
{extra}
set -euo pipefail
{prelude}
SYSTEMS=({systems})
ID=${{SYSTEMS[$SLURM_ARRAY_TASK_ID]}}
exec {exe} --batch-dir {batch} --only "$ID" --threads "${{SLURM_CPUS_PER_TASK}}" \\
    --deffnm {deffnm} --maxh {maxh} --gmx {gmx}{more}
"""


def write_slurm(b: Path, ids: list[str], a: argparse.Namespace) -> Path:
    extra = []
    if a.partition:
        extra.append(f"#SBATCH --partition={a.partition}")
    if a.account:
        extra.append(f"#SBATCH --account={a.account}")
    if a.gpus_per_task:
        extra.append(f"#SBATCH --gpus-per-task={a.gpus_per_task}")
    more = (f" --mdrun-args {shlex.quote(a.mdrun_args)}" if a.mdrun_args else "") + (" --force" if a.force else "")
    out = b / layout.SLURM
    out.mkdir(parents=True, exist_ok=True)
    script = out / "run_array.sbatch"
    script.write_text(
        SBATCH.format(
            job=a.job_name,
            last=len(ids) - 1,
            throttle=f"%{a.max_parallel}" if a.max_parallel else "",
            time=a.time,
            cpus=a.threads or 8,
            logdir=out.resolve(),
            extra="\n".join(extra),
            prelude="\n".join(a.prelude or []),
            systems=" ".join(shlex.quote(i) for i in ids),
            exe="gromatrix-run",
            batch=shlex.quote(str(b.resolve())),
            deffnm=shlex.quote(a.deffnm),
            # Leave 10 minutes (and >= 5 %) for mdrun to write its final checkpoint.
            maxh=round(max(_hours(a.time) * 0.95 - 1 / 6, 0.05), 3),
            gmx=shlex.quote(a.gmx),
            more=more,
        )
    )
    script.chmod(0o755)
    return script


def background(b: Path, argv: list[str]) -> dict:
    logs = b / layout.RUN_LOGS
    logs.mkdir(parents=True, exist_ok=True)
    log = logs / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.jsonl"
    args = [x for x in argv if x != "--background"]  # exact: the parser rejects abbreviations
    with open(log, "w") as fh:
        proc = subprocess.Popen(
            [sys.executable, "-m", "gromatrix.run", *args],
            stdout=fh,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            start_new_session=True,  # survives the terminal / agent session that started it
        )
    batch = {"mode": "background", "pid": proc.pid, "host": socket.gethostname(), "log": str(log), "started": now()}
    atomic_write(b / layout.BATCH_RUN, json.dumps(batch, indent=2))
    return {"ok": True, **batch, "next": f"gromatrix-status --batch-dir {b}  (or tail -f {log})"}


# ------------------------------------------------------------------ CLI
def cmd_run(a: argparse.Namespace, argv: list[str]) -> dict:
    b: Path = a.batch_dir
    if not (b / layout.SYSTEMS).is_dir():
        raise FileNotFoundError(f"{b}/systems not found: build the systems first (gromatrix-gmx build)")
    picked, skipped = runnable(b, a.only, a.force)
    if not picked:
        return {"ok": False, "error": "no runnable systems", "skipped": skipped}
    if a.slurm:
        script = write_slurm(b, [wd.name for wd, _ in picked], a)
        return {
            "ok": True,
            "mode": "slurm",
            "script": str(script),
            "systems": [wd.name for wd, _ in picked],
            "skipped": skipped,
            "next": f"sbatch {script}",
        }
    if a.background:
        return background(b, argv)

    env = ff_env(b, a.gmx)
    deadline = time.monotonic() + a.maxh * 3600 if a.maxh else None
    current: list[Runner] = []

    def on_term(signum, _frame):
        for r in current:
            r.request_stop()

    old = signal.signal(signal.SIGTERM, on_term)
    results = []
    try:
        for wd, stamp in picked:
            why = busy(wd)
            if why:
                results.append({"id": wd.name, "status": "busy", "reason": why})
                continue
            r = Runner(wd, a, env, deadline)
            current[:] = [r]
            st = r.execute(stamp.get("inputs_sha"))
            results.append(
                {"id": wd.name, "status": st["status"], **({"reason": st["reason"]} if "reason" in st else {})}
            )
            if st["status"] == "interrupted":
                break  # a signal or the time budget stops the whole batch, not just this system
    finally:
        signal.signal(signal.SIGTERM, old)
    summary = {
        "mode": "foreground",
        "only": a.only,
        "finished": now(),
        "systems": {r["id"]: r["status"] for r in results},
    }
    (b / layout.BATCH_RUN).parent.mkdir(parents=True, exist_ok=True)
    atomic_write(b / layout.BATCH_RUN, json.dumps(summary, indent=2))
    return {
        "ok": all(r["status"] == "done" for r in results),
        "done": [r["id"] for r in results if r["status"] == "done"],
        "failed": [r for r in results if r["status"] == "failed"],
        "interrupted": [r["id"] for r in results if r["status"] == "interrupted"],
        "not_started": [wd.name for wd, _ in picked if wd.name not in {r["id"] for r in results}],
        "busy": [r for r in results if r["status"] == "busy"],
        "skipped": skipped,
    }


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="gromatrix-run",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--only", nargs="+", help="system ids")
    ap.add_argument("--deffnm", default="md", help="file name of the production run (EM/NVT/NPT names are fixed)")
    ap.add_argument("--threads", type=int, help="mdrun -nt (Slurm: --cpus-per-task)")
    ap.add_argument("--gpu-id", help="mdrun -gpu_id, e.g. 0 or 01 (not used under Slurm: the scheduler binds GPUs)")
    ap.add_argument("--maxh", type=float, help="hours for this whole invocation; each mdrun gets what is left")
    ap.add_argument("--mdrun-args", default="", help='extra mdrun flags, e.g. "-pin on -nb gpu"')
    ap.add_argument("--gmx", default="gmx")
    ap.add_argument("--force", action="store_true", help="also run systems whose smoke test failed")
    ap.add_argument("--poll", type=float, default=30.0, help="seconds between progress lines")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--background", action="store_true", help="detach; follow B/07-run-logs/*.jsonl")
    mode.add_argument("--slurm", action="store_true", help="write B/07-slurm/run_array.sbatch instead of running")
    s = ap.add_argument_group("slurm")
    s.add_argument("--time", default="24:00:00", help="walltime per array task (Slurm --time format)")
    s.add_argument("--partition")
    s.add_argument("--account")
    s.add_argument("--gpus-per-task", type=int)
    s.add_argument("--max-parallel", type=int, help="array throttle (%%N)")
    s.add_argument("--job-name", default="gromatrix")
    s.add_argument("--prelude", action="append", help='shell line before the run, e.g. "module load gromacs"')
    a = ap.parse_args(argv)
    try:
        out = cmd_run(a, argv)
    except Exception as e:  # every failure is reported as the one-line JSON contract, never a traceback
        out = {"ok": False, "error": f"{type(e).__name__}: {e}"}
    print(json.dumps(out), flush=True)
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
