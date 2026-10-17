"""Shared helpers for the gromatrix plugin hooks (stdlib only).

Hooks are a second layer: the CLI gates (gromatrix-run, gromatrix-gmx build, ...) stay authoritative
and must keep working without hooks. A hook that crashes must never block the user, so scripts treat
unreadable input as "allow".
"""

from __future__ import annotations

import json
import re
import shlex
import sys
from pathlib import Path

WRAPPERS = {"sudo", "nohup", "time", "env", "command", "exec", "nice"}
# Mirrors the folder names in src/gromatrix/layout.py (tests/test_hooks.py keeps them in sync).
LAYOUT_DIR = re.compile(
    r"(^|/)0[0-9]-(config|raw|fixed|units|relaxed|alphafold(-cache)?|params|systems|run-logs|slurm|analysis|reports)(/|$)"
)
PIXI_VALUE_OPTS = {"-e", "--environment", "--manifest-path", "--color"}


def read_input() -> dict:
    try:
        data = json.load(sys.stdin)
    except ValueError:
        return {}
    return data if isinstance(data, dict) else {}


def deny(reason: str) -> None:
    """Block the tool call: exit 2, the reason on stderr is shown to Claude."""
    print(reason, file=sys.stderr)
    sys.exit(2)


def segments(command: str) -> list[list[str]]:
    """Split a shell command into simple commands (on ; & | newline) and tokenise each.

    Leading VAR=x assignments and wrappers (sudo, env, nohup, `pixi run -e md`) are dropped so
    tokens[0] is the program. Best effort, not a shell parser: quoted strings stay one token.
    """
    out = []
    for part in re.split(r"&&|\|\||[;&|\n]", command):
        try:
            toks = shlex.split(part)
        except ValueError:
            toks = part.split()
        toks = _strip_wrappers(toks)
        if toks:
            out.append(toks)
    return out


def _strip_wrappers(toks: list[str]) -> list[str]:
    while toks:
        if re.match(r"^[A-Za-z_]\w*=", toks[0]) or toks[0] in WRAPPERS:
            toks = toks[1:]
        elif toks[0] == "pixi" and toks[1:2] == ["run"]:
            toks = toks[2:]
            while toks and toks[0].startswith("-"):
                toks = toks[2:] if toks[0] in PIXI_VALUE_OPTS else toks[1:]
        else:
            break
    return toks


def option_values(toks: list[str], name: str) -> list[str]:
    """Values after `--name V...` (until the next option) or `--name=V`."""
    vals, grab = [], False
    for t in toks:
        if t == name:
            grab = True
        elif t.startswith(name + "="):
            vals.append(t.split("=", 1)[1])
            grab = False
        elif t.startswith("-"):
            grab = False
        elif grab:
            vals.append(t)
    return vals


def batch_dir(toks: list[str], cwd: str) -> Path | None:
    vals = option_values(toks, "--batch-dir")
    return Path(cwd or ".").joinpath(vals[0]) if vals else None
