"""End-to-end runner: steps A-H from one manifest.

    gmxprep-pipeline --manifest M [--dry-run] [--run-md] [--analyze]
                     [--skip-relax] [--skip-afdb] [--skip-smoke] [--from C] [--to F] [--keep-going]

The manifest is the step-F manifest (enzymes / targets / complexes / ligands)
plus pipeline keys. Paths are relative to the manifest file.

    batch_dir: work/batch1          # required (or --batch-dir)
    forcefield: amber               # locked before anything family-dependent
    allow_incomplete: false         # passed to step F; see below
    steps:                          # each key becomes that command's --flag
      fetch:     {ids: [1ABC, 2DEF]}        # present = start from PDB IDs (A-E run)
      extract:   {targets: LIG}
      relax:     {}                         # false = skip (as --skip-relax)
      afcompare: {}                         # false = skip (as --skip-afdb)
      build:     {md_ns: 100}
      smoke:     {nvt_steps: 5000}          # false = skip (as --skip-smoke)
      run:       {threads: 8, gpu_id: 0}    # only with --run-md
      analyze:   {}                         # only with --analyze

Without `enzymes:`, step F's manifest is generated from the step-C units after
A-E have run (relaxed files where they passed, never *.failed.pdb) and written
to B/manifest.generated.yaml. `ligands:` (SMILES + charge per HET) always comes
from the user. Crystal structures only contain native pairs, so a generated
matrix is usually incomplete: step F then stops unless `allow_incomplete: true`.

Each step runs in-process through its own CLI entry point, so the forcefield
lock, caches and gates are exactly those of the stage commands. After each step
one line {"step": "pipeline", "stage", "command", "status", "result"} is
printed; progress lines of step G pass through. The first step that reports
"ok": false stops the pipeline (fail-fast) unless --keep-going. The plan and
results are written to B/pipeline.json.
"""

from __future__ import annotations

import argparse
import importlib
import io
import json
import re
import shlex
import sys
from contextlib import redirect_stdout
from dataclasses import dataclass, field
from pathlib import Path

import yaml

LETTERS = "ABCDEFGH"
GENERATED = "<B/manifest.generated.yaml>"
PATH_KEYS = {"pdb", "sdf", "str"}


@dataclass
class Step:
    stage: str  # A-H, or "lock"
    command: str  # console script name
    module: str
    argv: list[str]
    skip: str | None = None  # reason, when the step is planned but not run
    deferred: bool = False  # argv completed at run time (generated manifest)
    result: dict = field(default_factory=dict)

    @property
    def label(self) -> str:
        return self.argv[0] if self.module == "gromatrix.gmx_prep" else self.command.removeprefix("gmxprep-")

    def shell(self) -> str:
        return shlex.join([self.command, *self.argv])


# ------------------------------------------------------------------ manifest
def to_flags(opts: dict | None) -> list[str]:
    """{max_loop: 8, ids: [A, B], keep_additives: true} -> --max-loop 8 --ids A B --keep-additives"""
    out: list[str] = []
    for k, v in (opts or {}).items():
        flag = "--" + str(k).replace("_", "-")
        if v is True:
            out.append(flag)
        elif v is False or v is None:
            continue
        elif isinstance(v, (list, tuple)):
            out += [flag, *map(str, v)]
        else:
            out += [flag, str(v)]
    return out


def load(path: Path, batch_override: Path | None) -> tuple[dict, Path]:
    raw = yaml.safe_load(Path(path).read_text()) or {}
    root = Path(path).resolve().parent
    b = batch_override or (root / raw["batch_dir"] if raw.get("batch_dir") else None)
    if b is None:
        raise ValueError("the manifest needs batch_dir (or pass --batch-dir)")
    return raw, Path(b).resolve()


def _abs(entry: dict, root: Path) -> dict:
    return {k: (str((root / v).resolve()) if k in PATH_KEYS and isinstance(v, str) else v) for k, v in entry.items()}


def derive_manifest(raw: dict, root: Path, b: Path) -> tuple[Path, list[str]]:
    """Write B/manifest.generated.yaml from units/*/unit.json (+ relax.json)."""
    notes: list[str] = []
    user_targets = [_abs(t, root) for t in raw.get("targets") or []]
    by_key = {}
    for t in user_targets:  # how step C names each target in its file names
        key = t.get("resname") if t.get("kind", "ligand") == "ligand" else "".join(sorted(t.get("chains") or []))
        by_key[key] = t
    enzymes, targets, complexes = [], {}, []
    for unit_p in sorted((b / "units").glob("*/unit.json")):
        unit = json.loads(unit_p.read_text())
        if unit.get("error"):
            notes.append(f"{unit_p.parent.name}: no unit selected ({unit['error'][:120]}); left out")
            continue
        pid = unit_p.parent.name
        relax = {}
        if (b / "relaxed" / pid / "relax.json").exists():
            relax = json.loads((b / "relaxed" / pid / "relax.json").read_text()).get("files", {})

        def pick(fname: str, pid=pid, relax=relax) -> str | None:
            st = relax.get(fname, {}).get("status")
            if st == "ok":
                return str(b / "relaxed" / pid / fname)
            if st == "failed":
                notes.append(f"{pid}/{fname}: failed relaxation ({relax[fname].get('failure')}); left out")
                return None
            return str(b / "units" / pid / fname)

        files = unit.get("files", {})
        enz = pick(files["enzyme"]) if "enzyme" in files else None
        if enz is None:
            continue
        enzymes.append({"id": pid, "pdb": enz})
        chain_tags = {"".join(sorted(unit.get("target_chains") or []))} - {""}
        for key, fname in sorted(files.items()):
            m = re.fullmatch(r"complex_([A-Za-z0-9]+?)(?:_alt([A-Za-z0-9]))?", key)
            if not m:
                continue
            tag, alt = m.groups()
            cpdb = pick(fname)
            if cpdb is None:
                continue
            t = by_key.get(tag)
            if t is None:  # not named by the user: derive it from step C's targets
                t = (
                    {"id": f"chains_{tag}", "kind": "protein", "chains": list(tag)}
                    if tag in chain_tags
                    else {"id": tag, "kind": "ligand", "resname": tag}
                )
                by_key[tag] = t
            if t["id"] not in targets:
                targets[t["id"]] = dict(t)
            tt = targets[t["id"]]
            if not tt.get("pdb") and not tt.get("sdf"):
                pose = files.get(f"target_{tag}" + (f"_alt{alt}" if alt else ""))
                tpdb = pick(pose) if pose else None
                if tpdb:
                    tt["pdb"] = tpdb
            complexes.append(
                {"enzyme": pid, "target": t["id"], "pdb": cpdb, **({"variant": f"alt{alt}"} if alt else {})}
            )
    if not enzymes:
        raise ValueError(f"no usable step-C units under {b}/units; cannot generate the step-F manifest")
    if not targets:
        raise ValueError("the units contain no target complexes; set steps.extract.targets (or target_chains)")
    man = {
        "forcefield": raw.get("forcefield"),
        "enzymes": enzymes,
        "targets": list(targets.values()),
        "complexes": complexes,
        "ligands": {res: _abs(d or {}, root) for res, d in (raw.get("ligands") or {}).items()},
    }
    out = b / "manifest.generated.yaml"
    out.write_text(
        "# Generated by gmxprep-pipeline from units/*/unit.json; edit the source manifest, not this file.\n"
        + yaml.safe_dump(man, sort_keys=False)
    )
    return out, notes


# ------------------------------------------------------------------ plan
def plan(raw: dict, root: Path, b: Path, a: argparse.Namespace) -> list[Step]:
    steps_cfg = raw.get("steps") or {}
    bd = ["--batch-dir", str(b)]

    def cfg(name):
        return to_flags(steps_cfg.get(name) if isinstance(steps_cfg.get(name), dict) else {})

    def off(name):
        return steps_cfg.get(name) is False

    ff = a.forcefield or raw.get("forcefield")
    lo, hi = LETTERS.index(a.start), LETTERS.index(a.stop)
    fetch_mode = "fetch" in steps_cfg
    user_f = bool(raw.get("enzymes"))
    f_manifest = str(a.manifest.resolve()) if user_f else GENERATED
    inc = ["--allow-incomplete"] if raw.get("allow_incomplete") else []

    out = []
    if (b / "forcefield.lock.json").exists() and not ff:
        pass  # an existing lock is used as is
    elif not ff:
        raise ValueError("no forcefield: set `forcefield:` in the manifest or pass --forcefield")
    else:
        out.append(Step("lock", "gmxprep-gmx", "gromatrix.gmx_prep", ["lock", *bd, "--forcefield", ff]))
    no_struct = None if fetch_mode else "no steps.fetch: the manifest brings its own structures"
    out += [
        Step("A", "gmxprep-fetch", "gromatrix.fetch", [*bd, *cfg("fetch")], no_struct),
        Step("B", "gmxprep-fix", "gromatrix.fix", [*bd, *cfg("fix")], no_struct),
        Step("C", "gmxprep-extract", "gromatrix.extract", [*bd, *cfg("extract")], no_struct),
        Step(
            "D",
            "gmxprep-relax",
            "gromatrix.relax",
            [*bd, *cfg("relax")],
            no_struct or ("skipped (--skip-relax)" if a.skip_relax or off("relax") else None),
        ),
        Step(
            "E",
            "gmxprep-afcompare",
            "gromatrix.afdb",
            [*bd, *cfg("afcompare")],
            no_struct or ("skipped (--skip-afdb)" if a.skip_afdb or off("afcompare") else None),
        ),
        Step(
            "F",
            "gmxprep-gmx",
            "gromatrix.gmx_prep",
            ["build", *bd, "--manifest", f_manifest, *inc, *cfg("build")],
            deferred=not user_f,
        ),
        Step("F", "gmxprep-gmx", "gromatrix.gmx_prep", ["verify", *bd]),
        Step(
            "F",
            "gmxprep-gmx",
            "gromatrix.gmx_prep",
            ["smoke", *bd, *cfg("smoke")],
            "skipped (--skip-smoke)" if a.skip_smoke or off("smoke") else None,
        ),
        Step(
            "G",
            "gmxprep-run",
            "gromatrix.run",
            [*bd, *cfg("run")],
            None if a.run_md else "not requested (--run-md)",
        ),
        Step(
            "H",
            "gmxprep-analyze",
            "gromatrix.analyze",
            [*bd, *cfg("analyze")],
            None if a.analyze else "not requested (--analyze)",
        ),
    ]
    for s in out:
        if s.stage != "lock" and not (lo <= LETTERS.index(s.stage) <= hi) and not s.skip:
            s.skip = f"outside --from {a.start} --to {a.stop}"
    if not fetch_mode and not user_f and any(s.stage == "F" and not s.skip for s in out):
        if not (b / "units").is_dir():
            raise ValueError("the manifest has neither steps.fetch nor enzymes/complexes, and B/units is empty")
    return out


# ------------------------------------------------------------------ execution
class _Tee(io.TextIOBase):
    """Captures a step's stdout; progress lines (they carry "step") pass straight through."""

    def __init__(self, sink):
        self.sink, self.buf, self.lines = sink, "", []

    def write(self, s: str) -> int:
        self.buf += s
        while "\n" in self.buf:
            line, self.buf = self.buf.split("\n", 1)
            self.lines.append(line)
            if line.startswith("{") and '"step"' in line:
                self.sink.write(line + "\n")
                self.sink.flush()
        return len(s)


def call(step: Step) -> dict:
    mod = importlib.import_module(step.module)
    tee = _Tee(sys.stdout)
    try:
        with redirect_stdout(tee):
            rc = mod.main(list(step.argv))
    except SystemExit as e:  # argparse rejected a flag from steps.<name>
        return {"ok": False, "error": f"{step.command}: bad arguments {step.argv} (exit {e.code})"}
    tee.write("\n")
    for line in reversed(tee.lines):
        if line.startswith("{") and '"step"' not in line:
            try:
                out = json.loads(line)
                out.setdefault("ok", rc == 0)
                return out
            except ValueError:
                continue
    return {"ok": rc == 0, "error": None if rc == 0 else f"{step.command} printed no result line"}


def execute(steps: list[Step], raw: dict, root: Path, b: Path, a: argparse.Namespace) -> dict:
    b.mkdir(parents=True, exist_ok=True)
    failed, notes = [], []
    for s in steps:
        rec = {"step": "pipeline", "stage": s.stage, "command": s.shell()}
        if s.skip:
            s.result = {"status": "skipped", "reason": s.skip}
            print(json.dumps({**rec, "status": "skipped", "reason": s.skip}), flush=True)
            continue
        if failed and not a.keep_going:
            s.result = {"status": "not run", "reason": f"stopped after failure in {failed[0]}"}
            continue
        if s.deferred:
            try:
                path, n = derive_manifest(raw, root, b)
            except ValueError as e:
                s.result = {"status": "failed", "ok": False, "error": str(e)}
                failed.append(s.stage)
                print(json.dumps({**rec, "status": "failed", "result": s.result}), flush=True)
                continue
            notes += n
            s.argv = [str(path) if x == GENERATED else x for x in s.argv]
            rec["command"] = s.shell()
        out = call(s)
        s.result = {"status": "ok" if out.get("ok") else "failed", **out}
        if not out.get("ok"):
            failed.append(s.stage)
        print(json.dumps({**rec, "status": s.result["status"], "result": out}), flush=True)
    report = {
        "ok": not failed,
        "batch": str(b),
        "failed_stage": failed[0] if failed else None,
        "steps": [{"stage": s.stage, "command": s.shell(), **s.result} for s in steps],
        "notes": notes,
        "next": f"gmxprep-status --batch-dir {b}",
    }
    (b / "pipeline.json").write_text(json.dumps(report, indent=1, default=str))
    summary = {k: v for k, v in report.items() if k != "steps"}
    summary["stages"] = [f"{s.stage} {s.label}: {s.result['status']}" for s in steps]
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gmxprep-pipeline", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--batch-dir", type=Path, help="overrides batch_dir in the manifest")
    ap.add_argument("--forcefield", help="overrides forcefield in the manifest (must match an existing lock)")
    ap.add_argument("--dry-run", action="store_true", help="print the plan (exact commands) and exit")
    ap.add_argument("--run-md", action="store_true", help="include step G (gmxprep-run)")
    ap.add_argument("--analyze", action="store_true", help="include step H (gmxprep-analyze)")
    ap.add_argument("--skip-relax", action="store_true", help="skip step D")
    ap.add_argument("--skip-afdb", action="store_true", help="skip step E")
    ap.add_argument("--skip-smoke", action="store_true", help="skip the step-F smoke test")
    ap.add_argument("--from", dest="start", default="A", choices=list(LETTERS))
    ap.add_argument("--to", dest="stop", default="H", choices=list(LETTERS))
    ap.add_argument("--keep-going", action="store_true", help="continue after a failed step")
    a = ap.parse_args(argv)
    try:
        raw, b = load(a.manifest, a.batch_dir)
        root = a.manifest.resolve().parent
        steps = plan(raw, root, b, a)
        if a.dry_run:
            out = {
                "ok": True,
                "batch": str(b),
                "plan": [
                    {"stage": s.stage, "command": s.shell(), **({"skip": s.skip} if s.skip else {})} for s in steps
                ],
            }
        else:
            out = execute(steps, raw, root, b, a)
    except Exception as e:  # every failure is reported as the one-line JSON contract, never a traceback
        out = {"ok": False, "error": f"{type(e).__name__}: {e}"}
    print(json.dumps(out), flush=True)
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
