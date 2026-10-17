"""PreToolUse (Bash): refuse commands that bypass the pipeline's gates.

- `-maxwarn N` with N > 0 (grompp warnings are findings, not noise);
- `gmx mdrun` (use gromatrix-run via the md-runner subagent: gates, resume, run.json);
- `rm -r` on a path inside a numbered batch folder (00-config ... 09-reports);
- `gromatrix-gmx build --force` when a targeted system already has md.xtc or run.json.

Second layer only: gromatrix-gmx build itself refuses to rebuild a system with run.json without --force.
"""

from __future__ import annotations

import re
from pathlib import Path

from _common import LAYOUT_DIR, batch_dir, deny, option_values, read_input, segments

GMX = {"gmx", "gmx_mpi", "gmx_d", "gmx_mpi_d"}


def _maxwarn(toks: list[str]) -> str | None:
    for i, t in enumerate(toks):
        m = re.match(r"^--?maxwarn(?:=(.*))?$", t)
        if m:
            val = m.group(1) if m.group(1) is not None else (toks[i + 1] if i + 1 < len(toks) else "")
            if val != "0":
                return "-maxwarn above 0 hides grompp warnings. Fix the cause (see skills/gmx-builder) instead."
    return None


def _rm(toks: list[str]) -> str | None:
    opts = [t for t in toks[1:] if t.startswith("-")]
    recursive = any(t == "--recursive" or (not t.startswith("--") and re.search(r"[rR]", t)) for t in opts)
    if recursive:
        for p in (t for t in toks[1:] if not t.startswith("-")):
            if LAYOUT_DIR.search(p.replace("\\", "/")):
                return f"rm -r on {p}: batch folders hold the pipeline state. Ask the user, or rerun the step with --force."
    return None


def _forced_build(toks: list[str], cwd: str) -> str | None:
    if "--force" not in toks:
        return None
    b = batch_dir(toks, cwd)
    if b is None:
        return None
    systems = b / "06-systems"
    only = option_values(toks, "--only")
    dirs = [systems / i for i in only] if only else [d for d in systems.glob("*") if d.is_dir()]
    hit = [d.name for d in dirs if (d / "md.xtc").exists() or (d / "run.json").exists()]
    if hit:
        return f"build --force would delete MD results of {hit[:5]}. Ask the user first; narrow it with --only."
    return None


def main() -> None:
    data = read_input()
    cmd = str((data.get("tool_input") or {}).get("command", ""))
    cwd = str(data.get("cwd") or ".")
    for toks in segments(cmd):
        prog = Path(toks[0]).name
        reason = None
        if prog in GMX or prog == "mdrun" or prog.startswith("gromatrix-"):
            reason = _maxwarn(toks)
        if reason is None and prog in GMX and toks[1:2] == ["mdrun"]:
            reason = "Direct gmx mdrun skips the gates and run.json. Use gromatrix-run (md-runner subagent)."
        if reason is None and prog == "rm":
            reason = _rm(toks)
        if reason is None and prog == "gromatrix-gmx" and toks[1:2] == ["build"]:
            reason = _forced_build(toks, cwd)
        if reason:
            deny(reason)


if __name__ == "__main__":
    main()
