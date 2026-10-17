"""SessionStart (startup, resume, clear, compact): inject `gromatrix-status` as additional context.

Runs only when the working directory looks like a batch (00-config/) or has a pipeline.yaml, and
only if `gromatrix-status` is on PATH. Silent no-op otherwise.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

from _common import read_input

LIMIT = 9500
KEEP = ("ok", "batch", "forcefield", "steps", "systems", "next")


def main() -> None:
    cwd = Path(read_input().get("cwd") or ".")
    if not ((cwd / "00-config").is_dir() or (cwd / "pipeline.yaml").is_file()):
        return
    exe = shutil.which("gromatrix-status")
    if not exe:
        return
    try:
        r = subprocess.run([exe, "--batch-dir", str(cwd)], capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return
    lines = r.stdout.strip().splitlines()
    line = lines[-1] if lines else ""
    if not line.startswith("{"):
        return
    if len(line) > LIMIT:  # keep the parts that decide what to do next
        try:
            full = json.loads(line)
            line = json.dumps({k: full.get(k) for k in KEEP})[:LIMIT]
        except ValueError:
            line = line[:LIMIT]
    ctx = f"gromatrix-status for {cwd} (from the SessionStart hook; may be stale later): {line}"
    print(json.dumps({"hookSpecificOutput": {"hookEventName": "SessionStart", "additionalContext": ctx}}))


if __name__ == "__main__":
    main()
