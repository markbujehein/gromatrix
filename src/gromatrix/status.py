"""Batch status: the hand-off contract between the gromatrix agent skills.

    gmxprep-status --batch-dir B

Reads only the artefacts each step writes (no heavy imports, no network) and
prints one JSON line: the forcefield lock, a done/partial/pending verdict for
each step A-H, per-structure progress through A-E, per-system counts for F-H,
and `next`, the skill and command that should run next. Every skill starts by
running this instead of re-reading files.

State lives entirely in the batch directory:
    forcefield.lock.json         lock            (set once, read by D and F)
    raw/<ID>.cif, cache.json     A  fetch
    fixed/<ID>.pdb, .meta.json   B  fix
    units/<ID>/unit.json         C  extract      (selected unit, targets, files, altloc plan)
    relaxed/<ID>/relax.json      D  relax        (per-file ok / failed + label)
    af/<ID>.json, af_report.csv  E  afcompare
    systems/<ID>/system.json     F  gmx build    (status, smoke result)
    systems/<ID>/run.json        G  run          (per-phase status; done/running/failed/interrupted)
    systems/<ID>/analysis/analysis.json  H  analyze  (metrics, flags), analysis_summary.csv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ORCHESTRATOR = "gromatrix"
STRUCTURE_SKILL, VALIDATE_SKILL, GMX_SKILL = "gromatrix-structure-prep", "gromatrix-validate", "gromatrix-gmx-builder"
RUN_SKILL, ANALYSIS_SKILL = "gromatrix-execution", "gromatrix-analysis"


def _json(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def batch_status(b: Path) -> dict:
    lock = _json(b / "forcefield.lock.json")
    ids = sorted(
        {p.stem.upper() for p in (b / "raw").glob("*") if p.suffix in (".cif", ".pdb")}
        | {p.name for p in (b / "units").glob("*") if p.is_dir()}
    )
    per_id: dict[str, dict] = {}
    for pid in ids:
        entry: dict = {
            "fetched": any((b / "raw" / f"{pid}.{e}").exists() for e in ("cif", "pdb")),
            "fixed": (b / "fixed" / f"{pid}.meta.json").exists(),
        }
        unit = _json(b / "units" / pid / "unit.json")
        if unit:
            entry["unit"] = unit.get("selected") or ("error" if unit.get("error") else None)
            entry["targets"] = unit.get("targets", [])
            entry["matrix_files"] = (
                sorted(unit.get("files", {}).values()) if isinstance(unit.get("files"), dict) else []
            )
        relax = _json(b / "relaxed" / pid / "relax.json")
        if relax:
            st = [f.get("status") for f in relax.get("files", {}).values()]
            entry["relax"] = {"ok": st.count("ok"), "failed": st.count("failed"), "skipped": st.count("skipped")}
        af = _json(b / "af" / f"{pid}.json")
        if af:
            entry["af"] = {
                "scored": sum(1 for r in af.values() if r.get("tm_score") is not None),
                "flagged": sum(1 for r in af.values() if r.get("status") == "flagged"),
            }
        per_id[pid] = entry

    built = [p.parent for p in (b / "systems").glob("*/system.json") if (_json(p) or {}).get("status") == "ok"]
    stamps = {wd.name: _json(wd / "system.json") or {} for wd in built}
    runs = {wd.name: _run_status(_json(wd / "run.json") or {}) for wd in built}
    analyses = {wd.name: (_json(wd / "analysis" / "analysis.json") or {}) for wd in built}
    sys_summary = {
        "built": len(built),
        "failed": sum((p / "error.txt").exists() for p in (b / "systems").glob("*/")),
        "smoke_ok": sum(bool(s.get("smoke", {}).get("ok")) for s in stamps.values()),
        "smoke_failed": sum(s.get("smoke", {}).get("ok") is False for s in stamps.values()),
        "run": {st: sorted(k for k, v in runs.items() if v == st) for st in RUN_STATES if st in runs.values()},
        "run_done": sum(v == "done" for v in runs.values()),
        "analysed": sum(r.get("status") == "ok" for r in analyses.values()),
        "analysis_failed": sorted(k for k, r in analyses.items() if r.get("status") == "failed"),
        "analysis_flags": {k: r["flags"] for k, r in analyses.items() if r.get("flags")},
    }
    return {
        "ok": True,
        "batch": str(b),
        "forcefield": lock.get("key") if lock else None,
        "steps": _steps(per_id, sys_summary),
        "structures": per_id,
        "systems": sys_summary,
        "next": _next(b, lock, per_id, sys_summary),
    }


RUN_STATES = ("done", "running", "interrupted", "failed")


def _run_status(st: dict) -> str | None:
    """A 'running' run.json whose runner died on this host (kill -9, reboot) is 'interrupted'."""
    from .run import owner_alive

    if st.get("status") == "running" and owner_alive(st) is False:
        return "interrupted"
    return st.get("status")


def _verdict(n_done: int, n_total: int) -> str:
    if n_total == 0:
        return "pending"
    return "done" if n_done == n_total else ("partial" if n_done else "pending")


def _steps(per_id: dict, systems: dict) -> dict:
    """A-H at a glance. A-E count structures, F-H count built systems.
    A batch built from the user's own PDBs has no structures: A-E are n/a."""
    n = len(per_id)
    out = {}
    for letter, key in zip("ABCDE", ("fetched", "fixed", "unit", "relax", "af")):
        out[letter] = (
            "n/a" if not n and systems["built"] else _verdict(sum(bool(v.get(key)) for v in per_id.values()), n)
        )
    nb = systems["built"]
    out["F"] = "partial" if systems["failed"] and nb else ("done" if nb else "pending")
    out["G"] = _verdict(systems["run_done"], nb - systems["smoke_failed"])  # G never runs smoke failures
    out["H"] = _verdict(systems["analysed"], systems["run_done"]) if systems["run_done"] else "pending"
    return out


def _next(b: Path, lock: dict | None, per_id: dict, systems: dict) -> dict:
    if lock is None:
        return {
            "skill": ORCHESTRATOR,
            "command": f"gmxprep-gmx lock --batch-dir {b} --forcefield <amber|charmm|...>",
            "why": "no forcefield lock yet; decide it before anything family-dependent (D, F)",
        }
    if not per_id and not systems["built"]:
        return {"skill": STRUCTURE_SKILL, "command": f"gmxprep-fetch --batch-dir {b} --ids <IDs>", "why": "empty batch"}
    todo = lambda key: sorted(k for k, v in per_id.items() if not v.get(key))  # noqa: E731
    if systems["built"]:
        return _next_md(b, systems)
    if todo("fixed"):
        return {
            "skill": STRUCTURE_SKILL,
            "command": f"gmxprep-fix --batch-dir {b}",
            "why": f"unfixed: {todo('fixed')[:5]}",
        }
    if todo("unit"):
        return {
            "skill": STRUCTURE_SKILL,
            "command": f"gmxprep-extract --batch-dir {b} --targets <RES>",
            "why": f"no unit selected: {todo('unit')[:5]}",
        }
    if todo("relax"):
        return {
            "skill": VALIDATE_SKILL,
            "command": f"gmxprep-relax --batch-dir {b}",
            "why": f"not relaxed: {todo('relax')[:5]}",
        }
    if todo("af"):
        return {
            "skill": VALIDATE_SKILL,
            "command": f"gmxprep-afcompare --batch-dir {b}",
            "why": f"no AlphaFold comparison: {todo('af')[:5]}",
        }
    if not systems["built"]:
        return {
            "skill": GMX_SKILL,
            "command": f"gmxprep-gmx build --batch-dir {b} --manifest <manifest.yaml>",
            "why": "structures validated; no GROMACS systems built yet",
        }
    return _next_md(b, systems)


def _next_md(b: Path, systems: dict) -> dict:
    """F (smoke) -> G -> H, once systems exist."""
    runs = systems["run"]
    if systems["smoke_ok"] + systems["smoke_failed"] < systems["built"] and not systems["run_done"]:
        return {"skill": GMX_SKILL, "command": f"gmxprep-gmx smoke --batch-dir {b}", "why": "systems not smoke-tested"}
    if systems["smoke_failed"] and systems["smoke_failed"] == systems["built"]:
        return {"skill": GMX_SKILL, "command": f"gmxprep-gmx verify --batch-dir {b}", "why": "every smoke test failed"}
    if runs.get("failed"):
        return {
            "skill": RUN_SKILL,
            "command": f"gmxprep-run --batch-dir {b} --only {' '.join(runs['failed'][:5])}",
            "why": f"MD failed for {runs['failed'][:5]}: read systems/<id>/run.json, fix, rerun",
        }
    if systems["run_done"] + systems["smoke_failed"] < systems["built"]:
        if runs.get("running"):
            return {
                "skill": RUN_SKILL,
                "command": f"gmxprep-status --batch-dir {b}",
                "why": f"MD running for {runs['running'][:5]}; check again later",
            }
        return {
            "skill": RUN_SKILL,
            "command": f"gmxprep-run --batch-dir {b}",
            "why": "systems ready; MD not run or not finished (gmxprep-run resumes from checkpoints)",
        }
    if systems["analysis_failed"]:
        return {
            "skill": ANALYSIS_SKILL,
            "command": f"gmxprep-analyze --batch-dir {b} --only {' '.join(systems['analysis_failed'][:5])}",
            "why": f"analysis failed for {systems['analysis_failed'][:5]}",
        }
    if systems["analysed"] < systems["run_done"]:
        return {
            "skill": ANALYSIS_SKILL,
            "command": f"gmxprep-analyze --batch-dir {b}",
            "why": "trajectories not analysed yet",
        }
    return {
        "skill": None,
        "command": f"gmxprep-gmx verify --batch-dir {b}",
        "why": "pipeline complete (A-H); verify and write the final report from analysis_summary.csv",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gmxprep-status", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    a = ap.parse_args(argv)
    if not a.batch_dir.exists():
        print(
            json.dumps(
                {
                    "ok": False,
                    "error": f"{a.batch_dir} does not exist",
                    "next": {"skill": STRUCTURE_SKILL, "command": f"gmxprep-fetch --batch-dir {a.batch_dir} ..."},
                }
            )
        )
        return 1
    print(json.dumps(batch_status(a.batch_dir)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
