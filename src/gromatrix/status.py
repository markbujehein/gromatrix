"""Batch status: the hand-off contract between the gromatrix agent skills.

    gmxprep-status --batch-dir B

Reads only the artefacts each step writes (no heavy imports, no network) and
prints one JSON line: the forcefield lock, per-structure progress through steps
A-E, step-F system counts, and `next`, the skill and command that should run
next. Every skill starts by running this instead of re-reading files.

State lives entirely in the batch directory:
    forcefield.lock.json         lock            (set once, read by D and F)
    raw/<ID>.cif, cache.json     A  fetch
    fixed/<ID>.pdb, .meta.json   B  fix
    units/<ID>/unit.json         C  extract      (selected unit, targets, files, altloc plan)
    relaxed/<ID>/relax.json      D  relax        (per-file ok / failed + label)
    af/<ID>.json, af_report.csv  E  afcompare
    systems/<ID>/system.json     F  gmx build    (status, smoke result)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ORCHESTRATOR = "gromatrix"
STRUCTURE_SKILL, VALIDATE_SKILL, GMX_SKILL = "gromatrix-structure-prep", "gromatrix-validate", "gromatrix-gmx-builder"


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

    systems = [_json(p) or {} for p in (b / "systems").glob("*/system.json")]
    sys_summary = {
        "built": sum(s.get("status") == "ok" for s in systems),
        "failed": sum((b / "systems" / p.name / "error.txt").exists() for p in (b / "systems").glob("*/")),
        "smoke_ok": sum(bool(s.get("smoke", {}).get("ok")) for s in systems),
        "smoke_failed": sum(s.get("smoke", {}).get("ok") is False for s in systems),
    }
    return {
        "ok": True,
        "batch": str(b),
        "forcefield": lock.get("key") if lock else None,
        "structures": per_id,
        "systems": sys_summary,
        "next": _next(b, lock, per_id, sys_summary),
    }


def _next(b: Path, lock: dict | None, per_id: dict, systems: dict) -> dict:
    if lock is None:
        return {
            "skill": ORCHESTRATOR,
            "command": f"gmxprep-gmx lock --batch-dir {b} --forcefield <amber|charmm|...>",
            "why": "no forcefield lock yet; decide it before anything family-dependent (D, F)",
        }
    if not per_id:
        return {"skill": STRUCTURE_SKILL, "command": f"gmxprep-fetch --batch-dir {b} --ids <IDs>", "why": "empty batch"}
    todo = lambda key: sorted(k for k, v in per_id.items() if not v.get(key))  # noqa: E731
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
    if systems["smoke_ok"] + systems["smoke_failed"] < systems["built"]:
        return {"skill": GMX_SKILL, "command": f"gmxprep-gmx smoke --batch-dir {b}", "why": "systems not smoke-tested"}
    return {
        "skill": None,
        "command": f"gmxprep-gmx verify --batch-dir {b}",
        "why": "pipeline complete; verify and report",
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
