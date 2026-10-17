"""Follow the MD of a batch: one line per phase change per system, exit when every system is finished.

    gromatrix-watch --batch-dir B [--interval 30] [--once]

Reads 06-systems/<id>/run.json (and, for phases it does not list yet, the phase logs through
checks.scan_md_log). Prints one JSON line per change:
    {"system": "E1", "phase": "npt", "status": "running", "run": "running"}
Exit code: 0 when every system is done, 1 when any failed or was interrupted, 2 when there is nothing
to watch. `--once` prints the current state and exits. Plain stdout, so any background wait works.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from . import checks, layout, result

PHASE_NAMES = ("em", "nvt", "npt", "md")  # the phase order of run.PHASES (kept here: run imports are heavy)
FINAL = ("done", "failed", "interrupted")


def _json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def _runnable(system_json: Path) -> bool:
    """Built systems step G will run: a failed smoke test is skipped (unless forced, then run.json exists)."""
    st = _json(system_json)
    smoke_failed = st.get("smoke", {}).get("ok") is False
    return st.get("status") == "ok" and (not smoke_failed or (system_json.parent / "run.json").exists())


def system_state(wd: Path) -> tuple[str, dict[str, str]]:
    """(run status, {phase: status}) of one system."""
    from .run import owner_alive

    st = _json(wd / "run.json")
    run = st.get("status") or "pending"
    if run == "running" and owner_alive(st) is False:
        run = "interrupted"  # runner died (kill -9, reboot) before it could write a final status
    deffnm = st.get("deffnm") or "md"
    phases = {}
    for name in PHASE_NAMES:
        entry = (st.get("phases") or {}).get(name)
        log = wd / f"{deffnm if name == 'md' else name}.log"
        if entry and entry.get("status"):
            phases[name] = entry["status"]
        elif log.exists() and run != "pending":
            scan = checks.scan_md_log(log)
            phases[name] = "failed" if scan["fatal"] else ("done" if scan["finished"] else "running")
    return run, phases


def watch(batch: Path, interval: float, once: bool) -> int:
    seen: dict[tuple[str, str], str] = {}
    while True:
        wds = sorted(p.parent for p in (batch / layout.SYSTEMS).glob("*/system.json") if _runnable(p))
        if not wds:
            result.finish(result.fail("input_missing", f"no built systems in {batch / layout.SYSTEMS}"), flush=True)
            return 2
        runs = {}
        for wd in wds:
            runs[wd.name], phases = system_state(wd)
            for phase, status in phases.items():
                if seen.get((wd.name, phase)) != status:
                    seen[(wd.name, phase)] = status
                    rec = {"system": wd.name, "phase": phase, "status": status, "run": runs[wd.name]}
                    print(json.dumps(rec), flush=True)
        if once or all(r in FINAL for r in runs.values()):
            return 1 if any(r in ("failed", "interrupted") for r in runs.values()) else 0
        time.sleep(interval)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gromatrix-watch", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--interval", type=float, default=30.0, help="seconds between checks")
    ap.add_argument("--once", action="store_true", help="print the current state and exit")
    a = ap.parse_args(argv)
    try:
        return watch(a.batch_dir, a.interval, a.once)
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    sys.exit(main())
