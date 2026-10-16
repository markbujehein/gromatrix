"""Alternate-conformer (altloc) handling shared by steps B and C.

Policy, deterministic and never averaging:

1. Never average coordinates. Alternate conformers are separate partial-occupancy
   states. Their mean is not a physical structure: wherever conformers differ by a
   ring flip or rotamer, averaging collapses bonds (e.g. aromatic C-C to ~0.7 Å).
2. Primary conformer per residue = highest mean occupancy; ties go to the
   alphabetically first ID. This is what goes into fixed/<ID>.pdb.
3. Every conformer is kept in meta.json. For each pair we report the raw RMSD
   and a symmetry-corrected RMSD: atoms are re-matched within each element
   (Hungarian assignment), so a flip of a symmetric ring, which is a
   relabelling and not a new pose, is not counted as a difference.
4. A conformer is *distinct* when its symmetry-corrected RMSD to the primary
   is >= the threshold (default 0.5 Å, about the coordinate error at 2-3 Å
   resolution). A conformer whose bond lengths differ from the primary's by
   > 0.2 Å is flagged inconsistent (a broken or averaged model) and never used.
5. Step C branches the matrix (complex_<T>_alt<L>.pdb) only for target
   ligands with distinct, consistent conformers. Residues sharing altloc
   letter L (e.g. coupled side chains) switch to L together, following the
   wwPDB convention that one letter is one state. A branch whose pose does not
   touch the unit is rejected: some deposits store a second copy of a ligand,
   bound at another site, as an altloc of the same residue.
"""

from __future__ import annotations

from collections import OrderedDict, defaultdict
from pathlib import Path

import numpy as np

from .pdbio import Atom

BOND_TOL = 0.2  # Å, max bond-length disagreement between conformers
BOND_CUT = {"C": 1.9, "N": 1.9, "O": 1.9, "S": 2.0, "P": 2.0, "SE": 2.1, "CL": 2.0, "BR": 2.1, "I": 2.3, "F": 1.6}


def _clean(v: str) -> str:
    return "" if v in (".", "?") else v


def parse_all_atoms(path: Path) -> list[tuple[Atom, str]]:
    """Every atom of the first model as (Atom, altloc_id); altloc '' = shared."""
    path = Path(path)
    out: list[tuple[Atom, str]] = []
    if path.suffix.lower() in (".cif", ".mmcif"):
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict

        d = MMCIF2Dict(str(path))
        cols = [
            "group_PDB",
            "auth_atom_id",
            "auth_comp_id",
            "auth_asym_id",
            "auth_seq_id",
            "pdbx_PDB_ins_code",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "type_symbol",
            "B_iso_or_equiv",
            "occupancy",
            "label_alt_id",
            "pdbx_PDB_model_num",
        ]
        n = len(d.get("_atom_site.Cartn_x", []))
        required = {"auth_atom_id", "auth_comp_id", "auth_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z"}
        missing = [c for c in required if f"_atom_site.{c}" not in d]
        if missing:
            raise ValueError(f"{path.name}: mmCIF atom_site lacks required columns {missing}")
        # Optional columns may be absent in minimal mmCIFs: pad with '?' instead of zipping to nothing.
        defaults = {
            "group_PDB": "ATOM",
            "pdbx_PDB_ins_code": "?",
            "type_symbol": "?",
            "B_iso_or_equiv": "0",
            "occupancy": "1",
            "label_alt_id": ".",
            "pdbx_PDB_model_num": "1",
        }
        data = [d.get(f"_atom_site.{c}") or [defaults.get(c, "?")] * n for c in cols]
        data = [v if isinstance(v, list) else [v] for v in data]
        first = None
        for rec, name, comp, ch, seq, ins, x, y, z, el, b, occ, alt, model in zip(*data):
            first = first or model
            if model != first:
                break
            out.append(
                (
                    Atom(
                        rec,
                        name.strip('"'),
                        comp,
                        ch,
                        int(seq),
                        _clean(ins),
                        float(x),
                        float(y),
                        float(z),
                        (el if el not in (".", "?") else "".join(c for c in name if c.isalpha())[:1]).upper(),
                        float(b) if b not in (".", "?") else 0.0,
                        float(occ) if occ not in (".", "?") else 1.0,
                    ),
                    _clean(alt),
                )
            )
        return out
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        el = (line[76:78].strip() if len(line) >= 78 else "") or "".join(c for c in name if c.isalpha())[:1]
        out.append(
            (
                Atom(
                    line[:6].strip(),
                    name,
                    line[17:21].strip(),
                    line[21].strip(),
                    int(line[22:26]),
                    line[26].strip(),
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                    el.upper(),
                    float(line[60:66].strip() or 0),
                    float(line[54:60].strip() or 1),
                ),
                line[16].strip(),
            )
        )
    return out


def _res_key(a: Atom):
    """Residue *position*. The residue name is excluded so that microheterogeneity
    (e.g. altloc A = SER, B = THR at one position) is one group with two conformers."""
    return (a.chain, a.resseq, a.icode)


def resolve(records: list[tuple[Atom, str]]) -> tuple[list[Atom], list[dict]]:
    """Choose the primary conformer per residue position. Returns (primary atoms in
    file order, one group dict per position that had altlocs; hydrogens excluded)."""
    by_res: dict[tuple, list[tuple[Atom, str]]] = OrderedDict()
    for a, alt in records:
        by_res.setdefault(_res_key(a), []).append((a, alt))
    primary_atoms, groups = [], []
    for key, recs in by_res.items():
        alts = sorted({alt for _, alt in recs if alt})
        if not alts:
            primary_atoms += [a for a, _ in recs]
            continue
        occ = {L: float(np.mean([a.occ for a, alt in recs if alt == L])) for L in alts}
        primary = sorted(alts, key=lambda L: (-round(occ[L], 3), L))[0]
        primary_name = next(a.resname for a, alt in recs if alt == primary)
        # Shared atoms belong to every conformer, but only if they carry that conformer's
        # residue name; with microheterogeneity there are no shared atoms.
        primary_atoms += [a for a, alt in recs if alt == primary or (alt == "" and a.resname == primary_name)]
        heavy = [(a, alt) for a, alt in recs if a.element not in ("H", "D")]
        resname = {L: next(a.resname for a, alt in recs if alt == L) for L in alts}
        groups.append(
            {
                "chain": key[0],
                "resname": primary_name,
                "resseq": key[1],
                "icode": key[2],
                "primary": primary,
                "microheterogeneity": len(set(resname.values())) > 1,
                "alts": {
                    L: {
                        "occ": round(occ[L], 3),
                        "resname": resname[L],
                        "atoms": {a.name: [a.x, a.y, a.z] for a, alt in heavy if alt in ("", L)},
                        "elements": {a.name: a.element for a, alt in heavy if alt in ("", L)},
                    }
                    for L in alts
                },
            }
        )
    return primary_atoms, groups


def reduce_to_primary(path: Path, groups: list[dict], out_dir: Path) -> Path:
    """Write a copy of the input holding only the primary conformer of every residue,
    with altloc IDs cleared. All other records (SEQRES, assemblies, ...) are kept.
    PDBFixer then builds missing atoms and loops on the conformer we keep. Without
    this it builds on whichever altloc OpenMM met first."""
    primary = {(g["chain"], g["resseq"], g["icode"]): g["primary"] for g in groups}
    path, out_dir = Path(path), Path(out_dir)
    out = out_dir / path.name
    if path.suffix.lower() in (".cif", ".mmcif"):
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        from Bio.PDB.mmcifio import MMCIFIO

        d = MMCIF2Dict(str(path))
        alt = d.get("_atom_site.label_alt_id", [])
        keys = zip(
            d["_atom_site.auth_asym_id"],
            d["_atom_site.auth_seq_id"],
            d.get("_atom_site.pdbx_PDB_ins_code") or ["?"] * len(alt),
        )
        keep = [
            a in (".", "?") or a == primary.get((ch, int(seq), _clean(ins))) for a, (ch, seq, ins) in zip(alt, keys)
        ]
        for k in [k for k in d if k.startswith("_atom_site.")]:
            if isinstance(d[k], list) and len(d[k]) == len(keep):
                d[k] = [v for v, ok in zip(d[k], keep) if ok]
        if "_atom_site.label_alt_id" in d:
            d["_atom_site.label_alt_id"] = ["."] * sum(keep)
        io = MMCIFIO()
        io.set_dict(d)
        io.save(str(out))
        return out
    lines = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM", "ANISOU")) and len(line) > 26 and line[16] != " ":
            key = (line[21].strip(), int(line[22:26]), line[26].strip())
            if line[16] != primary.get(key):
                continue
            line = line[:16] + " " + line[17:]
        lines.append(line)
    out.write_text("\n".join(lines) + "\n")
    return out


# ------------------------------------------------------------------ geometry
def sym_rmsd(xa: np.ndarray, xb: np.ndarray, elements: list[str]) -> float:
    """RMSD after optimally re-matching atoms within each element. A lower
    bound that removes ring-flip / atom-label permutations."""
    from scipy.optimize import linear_sum_assignment

    perm = np.arange(len(xa))
    for e in set(elements):
        idx = np.array([i for i, el in enumerate(elements) if el == e])
        d = np.linalg.norm(xa[idx][:, None] - xb[idx][None], axis=2)
        r, c = linear_sum_assignment(d)
        perm[idx[r]] = idx[c]
    return float(np.sqrt(((xa - xb[perm]) ** 2).sum(1).mean()))


def bonds_from(x: np.ndarray, elements: list[str]) -> list[tuple[int, int]]:
    out = []
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            cut = max(BOND_CUT.get(elements[i], 1.9), BOND_CUT.get(elements[j], 1.9))
            if np.linalg.norm(x[i] - x[j]) < cut:
                out.append((i, j))
    return out


def analyse(group: dict, distinct_rmsd: float) -> dict:
    """Annotate a group with per-conformer RMSDs (raw and symmetry-corrected)
    against the primary, bond-geometry consistency and distinctness."""
    p = group["primary"]
    ref = group["alts"][p]
    report = {}
    for L, conf in group["alts"].items():
        names = [n for n in ref["atoms"] if n in conf["atoms"]]
        if not names:
            continue
        el = [ref["elements"][n] for n in names]
        xa = np.array([ref["atoms"][n] for n in names])
        xb = np.array([conf["atoms"][n] for n in names])
        bonds = bonds_from(xa, el)
        dev = max((abs(np.linalg.norm(xa[i] - xa[j]) - np.linalg.norm(xb[i] - xb[j])) for i, j in bonds), default=0.0)
        s = sym_rmsd(xa, xb, el)
        report[L] = {
            "occ": conf["occ"],
            "rmsd": round(float(np.sqrt(((xa - xb) ** 2).sum(1).mean())), 3),
            "sym_rmsd": round(s, 3),
            "max_bond_dev": round(float(dev), 3),
            "geometry_ok": bool(dev <= BOND_TOL),
            "distinct": bool(L != p and s >= distinct_rmsd),
        }
    group["analysis"] = report
    group["distinct"] = any(r["distinct"] and r["geometry_ok"] for r in report.values())
    return group


def summarise(groups: list[dict], polymer_res: set[str]) -> list[str]:
    """Human-readable one-liners for HET groups only (protein altlocs are routine)."""
    lines = []
    for g in groups:
        if g["resname"] in polymer_res:
            continue
        others = {L: r for L, r in g["analysis"].items() if L != g["primary"]}
        desc = ", ".join(
            f"{L}: occ {r['occ']}, RMSD {r['rmsd']} Å (sym {r['sym_rmsd']} Å)"
            + ("" if r["geometry_ok"] else ", INCONSISTENT GEOMETRY - excluded")
            for L, r in others.items()
        )
        verdict = (
            "distinct poses -> step C can branch"
            if g["distinct"]
            else "equivalent within coordinate error -> primary used"
        )
        lines.append(
            f"{g['resname']}{g['resseq']}{g['icode']}:{g['chain']} primary {g['primary']} "
            f"(occ {g['analysis'][g['primary']]['occ']}); {desc}; {verdict}"
        )
    return lines


def index(groups: list[dict]) -> dict[tuple, dict]:
    """(chain, primary resname, resseq, icode) -> group, matching extract.Res.key."""
    return {(g["chain"], g["resname"], g["resseq"], g["icode"]): g for g in groups}


def letters_for(groups_by_key: dict, keys: list[tuple], usable_only: bool = True) -> list[str]:
    """Altloc letters available across the given residues (consistent geometry only)."""
    out = set()
    for k in keys:
        g = groups_by_key.get(k)
        if not g:
            continue
        for L, r in g.get("analysis", {}).items():
            if not usable_only or r["geometry_ok"]:
                out.add(L)
    return sorted(out)


def alt_coords(group: dict, letter: str) -> dict[str, list[float]] | None:
    conf = group["alts"].get(letter)
    if conf is None or not group.get("analysis", {}).get(letter, {}).get("geometry_ok", True):
        return None
    return conf["atoms"]


def per_residue_counts(records) -> dict:
    counts = defaultdict(int)
    for a, alt in records:
        if alt:
            counts[_res_key(a)] += 1
    return counts
