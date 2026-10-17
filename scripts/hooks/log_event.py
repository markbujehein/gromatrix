"""PostToolUse / PostToolUseFailure (Bash, gromatrix-*): append one JSONL line to <batch>/00-config/events.jsonl.

Fields: ts, command, and exit / duration_ms when the hook input carries them. Never fails the tool call.
"""

from __future__ import annotations

import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

from _common import batch_dir, read_input, segments


def main() -> None:
    data = read_input()
    cmd = str((data.get("tool_input") or {}).get("command", ""))
    cwd = str(data.get("cwd") or ".")
    for toks in segments(cmd):
        if not toks[0].startswith("gromatrix-"):
            continue
        b = batch_dir(toks, cwd) or Path(cwd)
        if not (b / "00-config").is_dir():  # never create a batch folder as a side effect
            return
        rec: dict = {"ts": datetime.now(timezone.utc).isoformat(timespec="seconds"), "command": cmd[:500]}
        resp = data.get("tool_response")
        resp = resp if isinstance(resp, dict) else {}
        code = resp.get("exit_code", resp.get("exitCode"))
        if code is None and data.get("hook_event_name") == "PostToolUseFailure":
            m = re.search(r"[Ee]xit code (\d+)", str(data.get("error", "")))
            code = int(m.group(1)) if m else None
        if code is not None:
            rec["exit"] = code
        if isinstance(data.get("duration_ms"), (int, float)):
            rec["duration_ms"] = data["duration_ms"]
        with open(b / "00-config" / "events.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        return


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"log_event: {e}", file=sys.stderr)
    sys.exit(0)
