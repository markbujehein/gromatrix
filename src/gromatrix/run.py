"""[G] Run the MD protocol of every built system: EM -> NVT -> NPT -> production.

    gmxprep-run --batch-dir B [--only ID ...] [--threads 8] [--gpu-id 0] [--deffnm md]
    gmxprep-run --batch-dir B --background            # detach; progress in B/run_logs/
    gmxprep-run --batch-dir B --slurm --time 24:00:00 [--partition gpu --gpus-per-task 1]

Inputs are step F's systems/<id>/ (topol.top, index.ndx, mdp/, solv_ions.gro).
The phases chain exactly like run.sh: -r gives POSRES reference coordinates,
-t carries velocities and box/thermostat state forward. -maxwarn is never used.

Resume rules, checked per phase before anything runs:
  * done     the phase's .gro exists and the last run segment of its .log
             finished normally            -> skipped
  * partial  a .cpt exists                -> mdrun -cpi continues it
  * new      otherwise                    -> grompp, then mdrun
Once a phase actually runs, every later phase starts fresh, so no stale output
from an earlier attempt is ever chained. A system rebuilt by step F (new
inputs_sha) restarts from EM.

Gates: EM must converge (Fmax < emtol, Epot < 0); NVT/NPT must finish with no
LINCS warnings and no fatal error; production must finish without a fatal error
(LINCS warnings there are reported). A failed system stops; the others go on.

Output: progress lines {"step": "run", "system", "phase", "status", "progress_pct",
"overall_pct"} while mdrun runs, then one summary line {"ok": ...}. State per
system goes to systems/<id>/run.json, the batch run to B/run.json.
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
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import checks
from . import mdp as mdpmod
from .forcefields import load_lock
from .gmx import ToolError, gmxlib_env, locate_ff_dir, run


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
DYNAMICS = {"nvt", "npt", "md"}
# A log segment that ends like this did not run to nsteps (mdrun still prints "Finished mdrun").
STOPPED_EARLY = re.compile(r"Run time exceeded|Received the (?:TERM|INT|USR1) signal|will terminate the run", re.I)
STEP_RE = re.compile(r"^\s+Step\s+Time\s*\n\s+(\d+)\s+", re.M)


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def emit(rec: dict) -> None:
    print(json.dumps(rec), flush=True)


# ------------------------------------------------------------------ phase state
def deffnm(phase: Phase, prod: str) -> str:
    return prod if phase.name == "md" else phase.name


def last_segment(log: Path) -> str:
    """mdrun -cpi appends to the log; only the newest run segment says how it ended."""
    text = log.read_text(errors="replace") if log.exists() else ""
    i = text.rfind("Log file opened")
    return text[i:] if i >= 0 else text


def last_step(text: str) -> int | None:
    m = STEP_RE.findall(text)
    return int(m[-1]) if m else None


def phase_state(wd: Path, name: str) -> str:
    seg = last_segment(wd / f"{name}.log")
    finished = "Finished mdrun" in seg and not STOPPED_EARLY.search(seg) and not re.search(r"Fatal error", seg)
    if (wd / f"{name}.gro").exists() and finished:
        return "done"
    if (wd / f"{name}.cpt").exists():
        return "partial"
    return "new"


def clear_phase(wd: Path, name: str) -> None:
    for ext in ("tpr", "cpt", "gro", "log", "edr", "xtc", "trr", "_prev.cpt", "mdout.mdp"):
        p = wd / (f"{name}{ext}" if ext.startswith("_") else f"{name}.{ext}")
        p.unlink(missing_ok=True)


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

    def __init__(self, wd: Path, a: argparse.Namespace, env: dict):
        self.wd, self.a, self.env = wd, a, env
        self.state_path = wd / "run.json"
        self.log = wd / "run.log"
        self.child: subprocess.Popen | None = None
        steps = {p.name: nsteps_of(wd, p) for p in PHASES}
        self.weights = {k: v for k, v in steps.items() if k in DYNAMICS}
        self.total = sum(self.weights.values()) or 1
        self.nsteps = steps

    def g(self, *sub: str) -> list[str]:
        return [self.a.gmx, "-quiet", "-nobackup", *sub]

    def overall(self, name: str, step: int) -> float:
        done = 0
        for p in PHASES:
            if p.name == name:
                break
            done += self.weights.get(p.name, 0)
        return round(100 * (done + (step if name in DYNAMICS else 0)) / self.total, 1)

    def grompp(self, phase: Phase, prev: str | None) -> None:
        name = deffnm(phase, self.a.deffnm)
        c = phase.coords or f"{prev}.gro"
        cmd = ["grompp", "-f", f"mdp/{phase.name}.mdp", "-c", c, "-p", "topol.top", "-n", "index.ndx"]
        if phase.restrained:
            cmd += ["-r", c]
        if phase.continues:
            cmd += ["-t", f"{prev}.cpt"]
        run(self.g(*cmd, "-o", f"{name}.tpr", "-po", f"{name}.mdout.mdp"), self.wd, self.log, env_extra=self.env)

    def mdrun_cmd(self, name: str, resume: bool) -> list[str]:
        a = self.a
        cmd = self.g("mdrun", "-deffnm", name)
        if resume:
            cmd += ["-cpi", f"{name}.cpt"]
        if a.threads:
            cmd += ["-nt", str(a.threads)]
        if a.gpu_id is not None:
            cmd += ["-gpu_id", str(a.gpu_id)]
        if a.maxh:
            cmd += ["-maxh", str(a.maxh)]
        return cmd + shlex.split(a.mdrun_args)

    def mdrun(self, phase: Phase, resume: bool) -> None:
        name = deffnm(phase, self.a.deffnm)
        cmd = self.mdrun_cmd(name, resume)
        sysid = self.wd.name
        with open(self.log, "a") as fh:
            fh.write(f"\n$ {shlex.join(cmd)}\n")
            fh.flush()
            try:
                self.child = subprocess.Popen(
                    cmd, cwd=self.wd, stdout=fh, stderr=subprocess.STDOUT, env={**os.environ, **self.env}
                )
            except FileNotFoundError as e:
                raise ToolError(f"'{cmd[0]}' not found on PATH") from e
            n = self.nsteps.get(phase.name) or 0
            while True:
                try:
                    rc = self.child.wait(timeout=self.a.poll)
                    break
                except KeyboardInterrupt:
                    self.terminate()
                    self.child.wait()
                    self.child = None
                    raise Interrupted("interrupted by the user; rerun gmxprep-run to resume") from None
                except subprocess.TimeoutExpired:
                    step = last_step(last_segment(self.wd / f"{name}.log")) or 0
                    emit(
                        {
                            "step": "run",
                            "system": sysid,
                            "phase": phase.name,
                            "status": "running",
                            "progress_pct": round(100 * min(step, n) / n, 1) if n else None,
                            "overall_pct": self.overall(phase.name, min(step, n)),
                        }
                    )
        self.child = None
        if rc in (-signal.SIGTERM, -signal.SIGINT) or STOPPED_EARLY.search(last_segment(self.wd / f"{name}.log")):
            raise Interrupted(f"{phase.name} stopped before nsteps (signal or -maxh); rerun gmxprep-run to resume")
        if rc != 0:
            tail = "\n".join(self.log.read_text(errors="replace").strip().splitlines()[-20:])
            raise ToolError(f"mdrun {phase.name} exited with {rc} (see {self.log}):\n{tail}")

    def gate(self, phase: Phase) -> list[str]:
        """Raise on a hard failure; return warnings."""
        name = deffnm(phase, self.a.deffnm)
        log = self.wd / f"{name}.log"
        if phase.name == "em":
            emtol = float(mdpmod.read_mdp(self.wd / "mdp" / "em.mdp").get("emtol", 1000))
            em = checks.parse_em_log(log, emtol)
            if not em["converged"]:
                raise checks.CheckError(
                    f"EM did not converge (Fmax={em['fmax']}, Epot={em['epot']}, emtol={emtol}); "
                    "look for clashes or a bad ligand pose"
                )
            return []
        seg = last_segment(log)
        lincs = len(re.findall(r"LINCS WARNING", seg, re.I))
        if re.search(r"Fatal error", seg) or "Finished mdrun" not in seg:
            raise checks.CheckError(f"{phase.name}: mdrun did not finish cleanly (see {log.name})")
        if lincs and phase.name != "md":
            raise checks.CheckError(f"{phase.name}: {lincs} LINCS warnings during restrained equilibration")
        return [f"{phase.name}: {lincs} LINCS warnings"] if lincs else []

    def save(self, state: dict) -> None:
        state["updated"] = now()
        self.state_path.write_text(json.dumps(state, indent=2))

    def execute(self, inputs_sha: str | None) -> dict:
        old = json.loads(self.state_path.read_text()) if self.state_path.exists() else {}
        if old and old.get("inputs_sha") != inputs_sha:
            for p in PHASES:  # step F rebuilt the system: no earlier output is valid
                clear_phase(self.wd, deffnm(p, old.get("deffnm") or self.a.deffnm))
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
                st = "new" if ran else phase_state(self.wd, name)
                if ran:
                    clear_phase(self.wd, name)
                if st == "done":
                    state["phases"][phase.name] = {"status": "done", "skipped": True}
                    emit({"step": "run", "system": self.wd.name, "phase": phase.name, "status": "skipped"})
                    prev = name
                    continue
                if st == "new":
                    if not (phase.name == "em" and (self.wd / "em.tpr").exists()):
                        self.grompp(phase, prev)  # step F already wrote a gated em.tpr
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
                try:
                    state["warnings"] += self.gate(phase)
                except checks.CheckError:
                    # Keep the evidence, but never let a gated-out phase look complete on resume.
                    out_gro = self.wd / f"{name}.gro"
                    if out_gro.exists():
                        out_gro.replace(self.wd / f"{name}.failed.gro")
                    raise
                state["phases"][phase.name] = {"status": "done", "finished": now()}
                ran, prev = True, name
                emit(
                    {
                        "step": "run",
                        "system": self.wd.name,
                        "phase": phase.name,
                        "status": "done",
                        "progress_pct": 100.0,
                        "overall_pct": self.overall(phase.name, self.nsteps.get(phase.name, 0)),
                    }
                )
            state["status"] = "done"
            state["trajectory"] = f"{self.a.deffnm}.xtc"
        except Interrupted as e:
            state["status"], state["reason"] = "interrupted", str(e)
        except Exception as e:  # ToolError, CheckError: this system stops, the batch goes on
            state["status"], state["reason"] = "failed", f"{type(e).__name__}: {str(e)[:600]}"
        for p in state["phases"].values():
            if p.get("status") == "running":
                p["status"] = state["status"]
        state.pop("pid", None)
        self.save(state)
        return state

    def terminate(self) -> None:
        """Forward SIGTERM so mdrun writes a checkpoint and stops at the next step."""
        if self.child and self.child.poll() is None:
            self.child.send_signal(signal.SIGTERM)


# ------------------------------------------------------------------ selection
def runnable(b: Path, only: list[str] | None, force: bool) -> tuple[list[tuple[Path, dict]], list[dict]]:
    """Built systems that may run. A failed smoke test blocks a system unless --force."""
    picked, skipped = [], []
    for stamp_p in sorted((b / "systems").glob("*/system.json")):
        wd = stamp_p.parent
        if only and wd.name not in only:
            continue
        stamp = json.loads(stamp_p.read_text())
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


def busy(wd: Path) -> str | None:
    """Another live runner owns this system (same host, live pid)?"""
    p = wd / "run.json"
    if not p.exists():
        return None
    st = json.loads(p.read_text())
    if st.get("status") != "running" or not st.get("pid"):
        return None
    if st.get("host") != socket.gethostname():
        return f"running on {st.get('host')} (pid {st.get('pid')})"
    try:
        os.kill(int(st["pid"]), 0)
    except (OSError, ValueError):
        return None  # stale: the runner died without cleaning up
    return None if int(st["pid"]) == os.getpid() else f"running here (pid {st['pid']})"


def ff_env(b: Path, gmx: str) -> dict:
    try:
        return gmxlib_env(locate_ff_dir(load_lock(b).gmx_name, gmx, [b.resolve()]))
    except (ToolError, OSError, ValueError, KeyError):
        return {}  # systems carry an ff symlink (step F), so grompp still resolves the includes


# ------------------------------------------------------------------ slurm / background
def _hours(walltime: str) -> float:
    """Slurm time 'D-HH:MM:SS', 'HH:MM:SS', 'MM:SS' or minutes -> hours."""
    days, _, rest = walltime.rpartition("-")
    parts = [float(x) for x in rest.split(":")]
    if len(parts) == 1:
        h = parts[0] / 60
    else:
        parts = [0.0] * (3 - len(parts)) + parts
        h = parts[0] + parts[1] / 60 + parts[2] / 3600
    return h + 24 * float(days or 0)


SBATCH = """#!/usr/bin/env bash
# Generated by gromatrix (gmxprep-run --slurm). One array task per system.
# A task that hits the walltime stops cleanly (-maxh), leaving a checkpoint:
# resubmit this script and every task resumes where it stopped.
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
    --deffnm {deffnm} --maxh {maxh} --gmx {gmx}{mdrun_args}
"""


def write_slurm(b: Path, ids: list[str], a: argparse.Namespace) -> Path:
    extra = []
    if a.partition:
        extra.append(f"#SBATCH --partition={a.partition}")
    if a.account:
        extra.append(f"#SBATCH --account={a.account}")
    if a.gpus_per_task:
        extra.append(f"#SBATCH --gpus-per-task={a.gpus_per_task}")
    out = b / "slurm"
    out.mkdir(exist_ok=True)
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
            exe="gmxprep-run",
            batch=shlex.quote(str(b.resolve())),
            deffnm=shlex.quote(a.deffnm),
            # Leave 10 minutes (and >= 5 %) for mdrun to write its final checkpoint.
            maxh=round(max(_hours(a.time) * 0.95 - 1 / 6, 0.05), 3),
            gmx=shlex.quote(a.gmx),
            mdrun_args=f" --mdrun-args {shlex.quote(a.mdrun_args)}" if a.mdrun_args else "",
        )
    )
    script.chmod(0o755)
    return script


def background(b: Path, argv: list[str]) -> dict:
    logs = b / "run_logs"
    logs.mkdir(exist_ok=True)
    log = logs / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.jsonl"
    args = [a for a in argv if a != "--background"]
    with open(log, "w") as fh:
        proc = subprocess.Popen(
            [sys.executable, "-m", "gromatrix.run", *args],
            stdout=fh,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            start_new_session=True,  # survives the terminal / agent session that started it
        )
    batch = {"mode": "background", "pid": proc.pid, "host": socket.gethostname(), "log": str(log), "started": now()}
    (b / "run.json").write_text(json.dumps(batch, indent=2))
    return {"ok": True, **batch, "next": f"gmxprep-status --batch-dir {b}  (or tail -f {log})"}


# ------------------------------------------------------------------ CLI
def cmd_run(a: argparse.Namespace, argv: list[str]) -> dict:
    b: Path = a.batch_dir
    if not (b / "systems").is_dir():
        raise FileNotFoundError(f"{b}/systems not found: build the systems first (gmxprep-gmx build)")
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
    current: list[Runner] = []

    def on_term(signum, _frame):
        for r in current:
            r.terminate()

    old = signal.signal(signal.SIGTERM, on_term)
    results = []
    try:
        for wd, stamp in picked:
            why = busy(wd)
            if why:
                results.append({"id": wd.name, "status": "busy", "reason": why})
                continue
            r = Runner(wd, a, env)
            current[:] = [r]
            st = r.execute(stamp.get("inputs_sha"))
            results.append(
                {"id": wd.name, "status": st["status"], **({"reason": st["reason"]} if "reason" in st else {})}
            )
            if st["status"] == "interrupted":
                break  # a signal stops the whole batch, not just this system
    finally:
        signal.signal(signal.SIGTERM, old)
    summary = {
        "mode": "foreground",
        "finished": now(),
        "systems": {r["id"]: r["status"] for r in results},
    }
    if not a.only:
        (b / "run.json").write_text(json.dumps(summary, indent=2))
    return {
        "ok": all(r["status"] == "done" for r in results),
        "done": [r["id"] for r in results if r["status"] == "done"],
        "failed": [r for r in results if r["status"] == "failed"],
        "interrupted": [r["id"] for r in results if r["status"] == "interrupted"],
        "busy": [r for r in results if r["status"] == "busy"],
        "skipped": skipped,
    }


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="gmxprep-run", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--only", nargs="+", help="system ids")
    ap.add_argument("--deffnm", default="md", help="file name of the production run (EM/NVT/NPT names are fixed)")
    ap.add_argument("--threads", type=int, help="mdrun -nt (Slurm: --cpus-per-task)")
    ap.add_argument("--gpu-id", help="mdrun -gpu_id, e.g. 0 or 01 (not used under Slurm: the scheduler binds GPUs)")
    ap.add_argument("--maxh", type=float, help="mdrun -maxh: stop cleanly with a checkpoint after this many hours")
    ap.add_argument("--mdrun-args", default="", help='extra mdrun flags, e.g. "-pin on -nb gpu"')
    ap.add_argument("--gmx", default="gmx")
    ap.add_argument("--force", action="store_true", help="also run systems whose smoke test failed")
    ap.add_argument("--poll", type=float, default=30.0, help="seconds between progress lines")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--background", action="store_true", help="detach; follow B/run_logs/*.jsonl")
    mode.add_argument("--slurm", action="store_true", help="write B/slurm/run_array.sbatch instead of running")
    s = ap.add_argument_group("slurm")
    s.add_argument("--time", default="24:00:00", help="walltime per array task")
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
