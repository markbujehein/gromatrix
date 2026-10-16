"""Pre-flight and post-run checks for step F, adapted from patterns in SQUIP
(ESS-MimerAI/SQUIP: scripts/verify_systems.py, trajectory_validation/*.sh) and
made strict and automatic.

* topology <-> coordinate consistency: walks [ molecules ] through every
  included moleculetype and checks atom counts and residue names against the
  .gro, so failures come with a clear message instead of a terse grompp error
* net charge: integer before genion (catches broken ligand or hand-edited
  charges), zero after genion
* index groups by name (Solute/Solvent/per-ligand), written in Python rather
  than by piping group numbers into make_ndx/genion
* mdrun log parsing: EM convergence (Fmax < emtol), LINCS warnings,
  fatal errors, "Finished mdrun"
"""

from __future__ import annotations

import re
from pathlib import Path

from .pdbio import read_gro

SECTION_RE = re.compile(r"^\s*\[\s*([\w-]+)\s*\]")
INCLUDE_RE = re.compile(r'^\s*#include\s+"([^"]+)"')


class CheckError(RuntimeError):
    pass


def _resolve(inc: str, here: Path, search: list[Path]) -> Path | None:
    for base in [here] + list(search):
        p = base / inc
        if p.exists():
            return p
    return None


def moleculetypes(top: Path, search: list[Path], _seen: set | None = None) -> dict[str, dict]:
    """name -> {natoms, charge, resnames} for every moleculetype reachable from `top`.
    Every #include is followed regardless of #ifdef. Only [ atoms ] is read."""
    _seen = _seen if _seen is not None else set()
    top = Path(top).resolve()
    if top in _seen:
        return {}
    _seen.add(top)
    out: dict[str, dict] = {}
    section, cur = None, None
    for raw in top.read_text(errors="replace").splitlines():
        m = INCLUDE_RE.match(raw)
        if m:
            p = _resolve(m.group(1), top.parent, search)
            if p is not None:
                out.update(moleculetypes(p, search, _seen))
            continue
        line = raw.split(";")[0].strip()
        if not line or line.startswith("#"):
            continue
        m = SECTION_RE.match(line)
        if m:
            section = m.group(1).lower()
            continue
        if section == "moleculetype":
            cur = line.split()[0]
            out[cur] = {"natoms": 0, "charge": 0.0, "resnames": set()}
            section = "moleculetype-done"
        elif section == "atoms" and cur:
            f = line.split()
            out[cur]["natoms"] += 1
            out[cur]["resnames"].add(f[3])
            if len(f) > 6:
                out[cur]["charge"] += float(f[6])
    return out


def molecules(top: Path) -> list[tuple[str, int]]:
    out, section = [], None
    for raw in Path(top).read_text().splitlines():
        line = raw.split(";")[0].strip()
        m = SECTION_RE.match(line)
        if m:
            section = m.group(1).lower()
            continue
        if section == "molecules" and line:
            name, n = line.split()[:2]
            out.append((name, int(n)))
    return out


def top_vs_gro(top: Path, gro: Path, search: list[Path]) -> dict:
    """Walk [ molecules ] against the coordinate file. Returns
    {natoms_top, natoms_gro, charge, problems}."""
    types = moleculetypes(top, search)
    mols = molecules(top)
    atoms, _ = read_gro(gro)
    problems, pos, charge = [], 0, 0.0
    for name, n in mols:
        mt = types.get(name)
        if mt is None:
            problems.append(f"[ molecules ] lists '{name}' but no moleculetype of that name is included")
            continue
        charge += mt["charge"] * n
        block = atoms[pos : pos + mt["natoms"] * n]
        stray = {a.resname for a in block} - mt["resnames"]
        # gro residue names are cut to 5 characters
        stray = {r for r in stray if not any(x.startswith(r) for x in mt["resnames"])}
        if stray and not problems:
            problems.append(
                f"atoms {pos + 1}-{pos + len(block)} should be '{name}' "
                f"({sorted(mt['resnames'])}) but the .gro has residues {sorted(stray)[:5]}"
            )
        pos += mt["natoms"] * n
    if pos != len(atoms):
        problems.append(f"topology describes {pos} atoms, {gro.name} has {len(atoms)}")
    return {"natoms_top": pos, "natoms_gro": len(atoms), "charge": round(charge, 4), "problems": problems}


def assert_integer_charge(q: float, where: str, tol: float = 1e-3) -> int:
    if abs(q - round(q)) > tol:
        raise CheckError(
            f"{where}: net charge {q:.4f} is not an integer - a ligand or residue "
            "charge set is broken (check params/<RES>/meta.json and any edited .rtp)"
        )
    return int(round(q))


def write_index(
    gro: Path, top: Path, search: list[Path], out: Path, protein_mol: str | None, ligand_mols: list[str]
) -> dict[str, int]:
    """Named groups computed from [ molecules ]: System, Protein, one per
    ligand, Solute (everything before the first SOL) and Solvent (SOL plus
    the ions genion appended)."""
    types = moleculetypes(top, search)
    groups: dict[str, list[int]] = {"System": [], "Protein": [], "Solute": [], "Solvent": []}
    for lig in ligand_mols:
        groups[f"LIG_{lig}"] = []
    pos, in_solvent = 0, False
    for name, n in molecules(top):
        span = list(range(pos + 1, pos + 1 + types[name]["natoms"] * n))
        pos += len(span)
        in_solvent = in_solvent or name == "SOL"
        groups["System"] += span
        groups["Solvent" if in_solvent else "Solute"] += span
        if name == protein_mol:
            groups["Protein"] += span
        if name in ligand_mols:
            groups[f"LIG_{name}"] += span
    write_ndx(groups, out)
    return {g: len(v) for g, v in groups.items() if v}


def write_ndx(groups: dict[str, list[int]], out: Path) -> None:
    """GROMACS index file; empty groups are left out."""
    lines = []
    for g, idx in groups.items():
        if idx:
            lines.append(f"[ {g} ]")
            lines += [" ".join(f"{i:>6d}" for i in idx[k : k + 15]) for k in range(0, len(idx), 15)]
    Path(out).write_text("\n".join(lines) + "\n")


def read_ndx(path: Path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    cur = None
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith("["):
            cur = line.strip("[] ")
            groups[cur] = []
        elif cur and line:
            groups[cur] += [int(x) for x in line.split()]
    return groups


# ------------------------------------------------------------------ logs
# mdrun prints "Finished mdrun" even when -maxh or a signal stopped it before nsteps.
STOPPED_EARLY = re.compile(r"Run time exceeded|Received the (?:TERM|INT|USR1) signal|will terminate the run", re.I)


def read_log(log: Path) -> str:
    return Path(log).read_text(errors="replace") if Path(log).exists() else ""


def last_segment(text: str) -> str:
    """mdrun -cpi appends to the log; only the newest run segment says how it ended."""
    i = text.rfind("Log file opened")
    return text[i:] if i >= 0 else text


def parse_em_log(log: Path, emtol: float) -> dict:
    text = last_segment(read_log(log))
    fmax = re.search(r"Maximum force\s*=\s*([-\d.eE+]+)", text)
    epot = re.search(r"Potential Energy\s*=\s*([-\d.eE+]+)", text)
    res = {"fmax": float(fmax.group(1)) if fmax else None, "epot": float(epot.group(1)) if epot else None}
    res["converged"] = res["fmax"] is not None and res["fmax"] < emtol
    res["machine_precision"] = "machine precision" in text
    if res["epot"] is not None and res["epot"] > 0:
        res["converged"] = False  # positive potential energy after EM = overlapping atoms
    return res


def scan_md_log(log: Path) -> dict:
    """How the newest run segment ended; LINCS warnings are counted over every segment,
    so resuming a run (mdrun -cpi) never hides a warning from an earlier segment."""
    text = read_log(log)
    seg = last_segment(text)
    return {
        "finished": "Finished mdrun" in seg and not STOPPED_EARLY.search(seg),
        "stopped_early": bool(STOPPED_EARLY.search(seg)),
        "lincs_warnings": len(re.findall(r"LINCS WARNING", text, re.I)),
        "fatal": bool(re.search(r"Fatal error", seg)),
    }
