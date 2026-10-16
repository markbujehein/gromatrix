"""Topology text surgery: split/merge parameter sections, normalise ligand
moleculetypes, extract pdb2gmx moleculetypes and assemble topol.top."""
from __future__ import annotations

import re
from pathlib import Path

SECTION_RE = re.compile(r"^\s*\[\s*([\w-]+)\s*\]")
# Sections that define types and must precede every [ moleculetype ].
PARAM_SECTIONS = ("atomtypes", "bondtypes", "constrainttypes", "pairtypes", "angletypes",
                  "dihedraltypes", "nonbond_params", "cmaptypes")


class TopologyError(RuntimeError):
    pass


def _norm(line: str) -> str:
    return " ".join(line.split(";")[0].split())


def split_params(text: str) -> tuple[dict[str, list[str]], str]:
    """Separate type sections (returned per section, comments stripped) from the
    moleculetype text. A [ defaults ] block is dropped: forcefield.itp owns it."""
    params: dict[str, list[str]] = {}
    mol_lines: list[str] = []
    section = None
    for line in text.splitlines():
        m = SECTION_RE.match(line)
        if m:
            section = m.group(1).lower()
            if section not in PARAM_SECTIONS and section != "defaults":
                mol_lines.append(line)
            continue
        if section == "defaults":
            continue
        if section in PARAM_SECTIONS:
            if _norm(line):
                params.setdefault(section, []).append(_norm(line))
        else:
            mol_lines.append(line)
    return params, "\n".join(mol_lines).strip() + "\n"


def merge_params(blocks: list[dict[str, list[str]]]) -> str:
    """Merge type sections from several ligands. Exact duplicates are dropped;
    an atomtype redefined with different parameters is an error, because it
    would silently change another molecule's physics."""
    merged: dict[str, list[str]] = {s: [] for s in PARAM_SECTIONS}
    seen: dict[str, set[str]] = {s: set() for s in PARAM_SECTIONS}
    atomtype_def: dict[str, str] = {}
    for block in blocks:
        for section, lines in block.items():
            for line in lines:
                if section == "atomtypes":
                    name = line.split()[0]
                    if name in atomtype_def and atomtype_def[name] != line:
                        raise TopologyError(f"conflicting atomtype '{name}':\n  {atomtype_def[name]}\n  {line}")
                    atomtype_def[name] = line
                if line not in seen[section]:
                    seen[section].add(line)
                    merged[section].append(line)
    out = []
    for section in PARAM_SECTIONS:
        if merged[section]:
            out.append(f"[ {section} ]")
            out.extend(merged[section])
            out.append("")
    return "\n".join(out)


def iter_section(mol_text: str, name: str):
    section = None
    for line in mol_text.splitlines():
        m = SECTION_RE.match(line)
        if m:
            section = m.group(1).lower()
            continue
        if section == name and _norm(line) and not line.lstrip().startswith("#"):
            yield line


def itp_atoms(mol_text: str) -> list[dict]:
    atoms = []
    for line in iter_section(mol_text, "atoms"):
        f = _norm(line).split()
        atoms.append({"nr": int(f[0]), "type": f[1], "resname": f[3], "name": f[4],
                      "charge": float(f[6]), "mass": float(f[7]) if len(f) > 7 else 0.0})
    return atoms


def moleculetype_name(mol_text: str) -> str:
    return _norm(next(iter_section(mol_text, "moleculetype"))).split()[0]


def normalise_moleculetype(mol_text: str, resname: str) -> str:
    """Rename the moleculetype and every [ atoms ] residue column to `resname`."""
    out, section, renamed = [], None, False
    for line in mol_text.splitlines():
        m = SECTION_RE.match(line)
        if m:
            section = m.group(1).lower()
            out.append(line)
            continue
        body = _norm(line)
        if section == "moleculetype" and body and not renamed and not line.lstrip().startswith(";"):
            out.append(f"{resname:<10s} {body.split()[1]}")
            renamed = True
        elif section == "atoms" and body and not line.lstrip().startswith((";", "#")):
            f = body.split()
            f[3] = resname
            out.append("  ".join(f))
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def posre_block(mol_text: str, define: str = "POSRES", k: int = 1000) -> str:
    rows = [f"{a['nr']:6d}     1  {k}  {k}  {k}" for a in itp_atoms(mol_text) if a["mass"] > 1.5]
    return (f"\n#ifdef {define}\n[ position_restraints ]\n;  ai  funct  fcx  fcy  fcz\n"
            + "\n".join(rows) + "\n#endif\n")


def extract_pdb2gmx_moleculetype(top_text: str) -> str:
    """Pull the (single, -merge all) moleculetype plus its POSRES include out of a
    pdb2gmx .top so it can be included like any other .itp."""
    lines = top_text.splitlines()
    try:
        start = next(i for i, l in enumerate(lines) if SECTION_RE.match(l) and "moleculetype" in l)
    except StopIteration:
        raise TopologyError("no [ moleculetype ] in pdb2gmx topology")
    end = len(lines)
    for i in range(start + 1, len(lines)):
        l = lines[i]
        if l.startswith(("; Include water topology", "; Include topology for ions")) \
                or re.match(r'#include ".*\.ff/(tip|spc|opc|ions)', l) or (SECTION_RE.match(l) and "system" in l):
            end = i
            break
    return "\n".join(lines[start:end]).rstrip() + "\n"


def compress_molecules(molecules: list[tuple[str, int]]) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    for name, n in molecules:
        if out and out[-1][0] == name:
            out[-1] = (name, out[-1][1] + n)
        else:
            out.append((name, n))
    return out


def write_topol(path: Path, ff_gmx_name: str, water_itp: str, param_itp: str | None,
                moltype_itps: list[str], molecules: list[tuple[str, int]], title: str) -> None:
    inc = [f'#include "{ff_gmx_name}.ff/forcefield.itp"']
    if param_itp:
        inc.append(f'#include "{param_itp}"')
    inc += [f'#include "{i}"' for i in moltype_itps]
    inc += [f'#include "{ff_gmx_name}.ff/{water_itp}"', f'#include "{ff_gmx_name}.ff/ions.itp"']
    mols = "\n".join(f"{n:<20s} {c}" for n, c in compress_molecules(molecules))
    Path(path).write_text("; generated by gmxprep - do not edit by hand\n" + "\n".join(inc)
                          + f"\n\n[ system ]\n{title}\n\n[ molecules ]\n{mols}\n")
