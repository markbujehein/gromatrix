"""PreToolUse (Read, Bash): keep binary MD files and whole logs out of the context.

Denies reading *.xtc *.trr *.edr *.cpt, and dumping md.log / run.log (and the em/nvt/npt logs) whole.
"""

from __future__ import annotations

import re
from pathlib import Path

from _common import deny, read_input, segments

BINARY = (".xtc", ".trr", ".edr", ".cpt")
LOG = re.compile(r"^(run|em|nvt|npt|md)\.log$")
READERS = {"cat", "less", "more", "head", "tail", "xxd", "od", "hexdump", "strings", "bat", "grep", "rg", "nl", "most"}
DUMPERS = {"cat", "less", "more", "bat", "nl", "most"}  # read the whole file; head/tail/grep stay allowed
HINT = "Use `gromatrix-status --batch-dir B` or read the run.json; for a log, `tail -n 40 <log>`."


def main() -> None:
    data = read_input()
    ti = data.get("tool_input") or {}
    if data.get("tool_name") == "Read" or "file_path" in ti:
        name = Path(str(ti.get("file_path", ""))).name
        if name.endswith(BINARY):
            deny(f"{name} is binary and large; it is not for the context window. {HINT}")
        if LOG.match(name) and not ti.get("limit"):
            deny(f"Reading all of {name} floods the context. {HINT} (or pass a small `limit`).")
        return
    for toks in segments(str(ti.get("command", ""))):
        prog = Path(toks[0]).name
        for arg in (t for t in toks[1:] if not t.startswith("-")):
            name = Path(arg).name
            if prog in READERS and name.endswith(BINARY):
                deny(f"{prog} {name}: binary MD file. {HINT}")
            if prog in DUMPERS and LOG.match(name):
                deny(f"{prog} {name} dumps the whole log. {HINT}")


if __name__ == "__main__":
    main()
