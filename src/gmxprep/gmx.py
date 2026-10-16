"""Thin, logged wrappers around external programs (gmx, acpype, obabel, cgenff)."""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
from functools import cache
from pathlib import Path


class ToolError(RuntimeError):
    pass


def run(cmd: list[str], cwd: Path, log: Path, stdin: str | None = None, env_extra: dict | None = None) -> str:
    """Run `cmd`, appending full output to `log`. Raises ToolError with the log tail."""
    if shutil.which(cmd[0]) is None and not Path(cmd[0]).exists():
        raise ToolError(f"'{cmd[0]}' not found on PATH")
    env = os.environ.copy()
    env.update(env_extra or {})
    proc = subprocess.run(cmd, cwd=cwd, input=stdin, capture_output=True, text=True, env=env)
    with open(log, "a") as fh:
        fh.write(f"\n$ {shlex.join(cmd)}\n{proc.stdout}{proc.stderr}")
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout + proc.stderr).strip().splitlines()[-25:])
        raise ToolError(f"{cmd[0]} {cmd[1] if len(cmd) > 1 else ''} failed (see {log}):\n{tail}")
    return proc.stdout + proc.stderr


@cache
def data_prefix(gmx: str) -> Path | None:
    try:
        out = subprocess.run([gmx, "-quiet", "--version"], capture_output=True, text=True).stdout
    except FileNotFoundError:
        return None
    m = re.search(r"Data prefix:\s*(\S+)", out)
    return Path(m.group(1)) if m else None


def locate_ff_dir(gmx_name: str, gmx: str, extra: list[Path]) -> Path:
    """Find <gmx_name>.ff in the batch dir, cwd, $GMXLIB, then the GROMACS share dir."""
    candidates = list(extra) + [Path.cwd()]
    candidates += [Path(p) for p in os.environ.get("GMXLIB", "").split(os.pathsep) if p]
    prefix = data_prefix(gmx)
    if prefix:
        candidates.append(prefix / "share" / "gromacs" / "top")
    for base in candidates:
        d = base / f"{gmx_name}.ff"
        if (d / "forcefield.itp").exists():
            return d.resolve()
    raise ToolError(
        f"{gmx_name}.ff not found in {[str(c) for c in candidates]}; "
        "download the port and place it in the batch dir or on GMXLIB"
    )


def gmxlib_env(ff_dir: Path) -> dict:
    parts = [str(ff_dir.parent)] + [p for p in os.environ.get("GMXLIB", "").split(os.pathsep) if p]
    return {"GMXLIB": os.pathsep.join(dict.fromkeys(parts))}


def ion_atom_names(ff_dir: Path) -> dict[str, str]:
    """moleculetype -> atom name, parsed from <ff>/ions.itp."""
    names, section, mol = {}, None, None
    for raw in (ff_dir / "ions.itp").read_text().splitlines():
        line = raw.split(";")[0].strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("["):
            section = line.strip("[] ").lower()
            continue
        if section == "moleculetype":
            mol = line.split()[0]
        elif section == "atoms" and mol and mol not in names:
            names[mol] = line.split()[4]
    return names
