"""Deposition metadata shared by steps B and C: entities, biological assemblies,
symmetry operators and label->auth chain mapping, from mmCIF or PDB.

Everything is expressed in *author* chain IDs and residue numbers, the
convention used by the PDB files the rest of the pipeline reads and writes.
"""

from __future__ import annotations

import itertools
import re
from collections import OrderedDict
from pathlib import Path

from .pdbio import PROTEIN_RES, Atom, read_pdb

NUCLEIC_RES = {"DA", "DC", "DG", "DT", "DU", "A", "C", "G", "U", "I", "DI"}
POLYMER_RES = PROTEIN_RES | NUCLEIC_RES
IDENTITY = {"R": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], "t": [0.0, 0.0, 0.0]}


def _list(d: dict, key: str) -> list[str]:
    v = d.get(key, [])
    return v if isinstance(v, list) else [v]


def _clean(v: str) -> str:
    return "" if v in (".", "?") else v


# ------------------------------------------------------------------ operators
def parse_oper_expression(expr: str) -> list[tuple[str, ...]]:
    """'1' -> [('1',)]; '(1-3)' -> [('1',),('2',),('3',)];
    '(1,2)(5)' -> [('1','5'),('2','5')]. Operators in a tuple are applied
    right to left, following the PDBx convention for products."""
    groups = re.findall(r"\(([^)]*)\)", expr) or [expr]
    expanded = []
    for g in groups:
        ids = []
        for item in g.split(","):
            item = item.strip()
            m = re.fullmatch(r"(\d+)-(\d+)", item)
            ids += [str(i) for i in range(int(m.group(1)), int(m.group(2)) + 1)] if m else [item]
        expanded.append(ids)
    return [tuple(p) for p in itertools.product(*expanded)]


def compose(ops: dict, seq: tuple[str, ...]) -> dict:
    """Matrix product for an operator sequence (rightmost applied first)."""
    R = [row[:] for row in IDENTITY["R"]]
    t = IDENTITY["t"][:]
    for oid in reversed(seq):
        o = ops[oid]
        # x' = Ro (R x + t) + to
        R = [[sum(o["R"][i][k] * R[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        t = [sum(o["R"][i][k] * t[k] for k in range(3)) + o["t"][i] for i in range(3)]
    return {"R": R, "t": t}


def is_identity(op: dict, tol: float = 1e-4) -> bool:
    return all(abs(op["R"][i][j] - IDENTITY["R"][i][j]) < tol for i in range(3) for j in range(3)) and all(
        abs(x) < tol for x in op["t"]
    )


# ------------------------------------------------------------------ mmCIF
def _mmcif_meta(path: Path) -> dict:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    d = MMCIF2Dict(str(path))
    label_asym = _list(d, "_atom_site.label_asym_id")
    auth_asym = _list(d, "_atom_site.auth_asym_id")
    label2auth: dict[str, str] = OrderedDict()
    label_residues: dict[str, list] = {}
    keys = zip(
        label_asym,
        auth_asym,
        _list(d, "_atom_site.auth_comp_id") or _list(d, "_atom_site.label_comp_id"),
        _list(d, "_atom_site.auth_seq_id"),
        _list(d, "_atom_site.pdbx_PDB_ins_code"),
        _list(d, "_atom_site.label_entity_id"),
    )
    entity_of_label = {}
    for la, aa, comp, seq, ins, ent in keys:
        label2auth.setdefault(la, aa)
        entity_of_label.setdefault(la, ent)
        res = [aa, comp, int(seq), _clean(ins)]
        lst = label_residues.setdefault(la, [])
        if not lst or lst[-1] != res:
            lst.append(res)

    entities = {}
    for eid, etype, desc in zip(_list(d, "_entity.id"), _list(d, "_entity.type"), _list(d, "_entity.pdbx_description")):
        entities[eid] = {"type": etype, "description": _clean(desc), "chains": []}
    for eid, strands in zip(_list(d, "_entity_poly.entity_id"), _list(d, "_entity_poly.pdbx_strand_id")):
        entities[eid]["chains"] = [s.strip() for s in strands.split(",")]
    for eid, comp in zip(_list(d, "_pdbx_entity_nonpoly.entity_id"), _list(d, "_pdbx_entity_nonpoly.comp_id")):
        entities[eid]["resname"] = comp

    ops = {}
    for i, oid in enumerate(_list(d, "_pdbx_struct_oper_list.id")):
        R = [[float(_list(d, f"_pdbx_struct_oper_list.matrix[{r}][{c}]")[i]) for c in (1, 2, 3)] for r in (1, 2, 3)]
        t = [float(_list(d, f"_pdbx_struct_oper_list.vector[{r}]")[i]) for r in (1, 2, 3)]
        ops[oid] = {"R": R, "t": t}

    details = dict(zip(_list(d, "_pdbx_struct_assembly.id"), _list(d, "_pdbx_struct_assembly.details")))
    assemblies: dict[str, dict] = OrderedDict()
    for aid, expr, asyms in zip(
        _list(d, "_pdbx_struct_assembly_gen.assembly_id"),
        _list(d, "_pdbx_struct_assembly_gen.oper_expression"),
        _list(d, "_pdbx_struct_assembly_gen.asym_id_list"),
    ):
        poly, het = [], []
        for la in asyms.split(","):
            la = la.strip()
            if la not in label2auth:
                continue
            if entities.get(entity_of_label[la], {}).get("type") == "polymer":
                if label2auth[la] not in poly:
                    poly.append(label2auth[la])
            else:
                het += label_residues[la]
        det = details.get(aid, "")
        a = assemblies.setdefault(aid, {"id": aid, "details": det, "author_defined": "author" in det, "gens": []})
        a["gens"].append({"polymer_chains": poly, "het_residues": het, "opers": parse_oper_expression(expr)})

    return {
        "format": "mmcif",
        "entities": entities,
        "operators": ops,
        "assemblies": list(assemblies.values()),
        "label2auth": dict(label2auth),
        "hetnam": {e["resname"]: e["description"] for e in entities.values() if "resname" in e},
    }


# ------------------------------------------------------------------ PDB
def _pdb_meta(path: Path) -> dict:
    lines = Path(path).read_text().splitlines()
    entities, compnd = {}, " ".join(l[10:].strip() for l in lines if l.startswith("COMPND"))
    for block in re.split(r"MOL_ID:\s*", compnd)[1:]:
        mid = block.split(";")[0].strip()
        mol = re.search(r"MOLECULE:\s*([^;]*)", block)
        ch = re.search(r"CHAIN:\s*([^;]*)", block)
        entities[mid] = {
            "type": "polymer",
            "description": mol.group(1).strip() if mol else "",
            "chains": [c.strip() for c in ch.group(1).split(",")] if ch else [],
        }
    hetnam: dict[str, str] = {}
    for l in lines:
        if l.startswith("HETNAM"):
            code = l[11:14].strip()
            hetnam[code] = (hetnam.get(code, "") + " " + l[15:].strip()).strip()

    ops, assemblies, cur = {}, [], None
    pending_chains: list[str] = []
    for l in lines:
        if not l.startswith("REMARK 350"):
            continue
        body = l[10:].strip()
        if body.startswith("BIOMOLECULE:"):
            cur = {"id": body.split(":")[1].strip(), "details": "", "author_defined": False, "gens": []}
            assemblies.append(cur)
        elif cur is not None and body.startswith("AUTHOR DETERMINED"):
            cur["author_defined"] = True
            cur["details"] = "author_defined_assembly"
        elif cur is not None and ("APPLY THE FOLLOWING TO CHAINS:" in body or body.startswith("AND CHAINS:")):
            chains = [c.strip() for c in body.split(":", 1)[1].split(",") if c.strip()]
            if "APPLY" in body:
                pending_chains = chains
                cur["gens"].append({"polymer_chains": pending_chains, "het_residues": None, "opers": []})
            else:
                pending_chains.extend(chains)
        elif cur is not None and body.startswith("BIOMT"):
            f = body.split()
            row, oid = int(f[0][-1]) - 1, f"{cur['id']}.{f[1]}"
            op = ops.setdefault(oid, {"R": [[0.0] * 3 for _ in range(3)], "t": [0.0] * 3})
            op["R"][row] = [float(x) for x in f[2:5]]
            op["t"][row] = float(f[5])
            if row == 0:
                cur["gens"][-1]["opers"].append((oid,))
    # PDB convention: HETATM groups belong to an assembly through their chain ID.
    hets = [a for a in read_pdb(path) if a.resname not in POLYMER_RES]
    for asm in assemblies:
        if not asm["details"]:
            asm["details"] = "software_defined_assembly"
        for g in asm["gens"]:
            keys = []
            for a in hets:
                k = [a.chain, a.resname, a.resseq, a.icode]
                if a.chain in g["polymer_chains"] and k not in keys:
                    keys.append(k)
            g["het_residues"] = keys
    return {
        "format": "pdb",
        "entities": entities,
        "operators": ops,
        "assemblies": assemblies,
        "label2auth": {},
        "hetnam": hetnam,
    }


def read_metadata(path: Path) -> dict:
    path = Path(path)
    meta = _mmcif_meta(path) if path.suffix.lower() in (".cif", ".mmcif") else _pdb_meta(path)
    meta["source"] = str(path)
    meta["chain_entity"] = {
        c: eid for eid, e in meta["entities"].items() if e["type"] == "polymer" for c in e["chains"]
    }
    return meta


def read_raw_atoms(path: Path, distinct_rmsd: float = 0.5) -> tuple[list[Atom], dict, list[dict]]:
    """Primary-conformer atoms (first model, author numbering, B-factors),
    statistics, and one analysed group per residue with alternate conformers
    (see altlocs.py for the selection policy)."""
    from . import altlocs

    records = altlocs.parse_all_atoms(Path(path))
    atoms, groups = altlocs.resolve(records)
    groups = [altlocs.analyse(g, distinct_rmsd) for g in groups]
    stats = {
        "altloc_atoms": sum(1 for _, alt in records if alt),
        "hydrogens": sum(a.element in ("H", "D") for a in atoms),
        "altloc_residues": len(groups),
        "altloc_hets": [f"{g['resname']}{g['resseq']}{g['chain']}" for g in groups if g["resname"] not in POLYMER_RES],
    }
    return atoms, stats, groups


def entities_by_sequence(atoms: list[Atom]) -> dict[str, str]:
    """Fallback entity assignment (e.g. docked models without COMPND):
    chains with identical polymer sequences share an entity."""
    seqs: dict[str, list[str]] = OrderedDict()
    last = {}
    for a in atoms:
        if a.resname in POLYMER_RES and last.get(a.chain) != (a.resseq, a.icode):
            seqs.setdefault(a.chain, []).append(a.resname)
            last[a.chain] = (a.resseq, a.icode)
    ent, out = {}, {}
    for ch, s in seqs.items():
        out[ch] = ent.setdefault(tuple(s), str(len(ent) + 1))
    return out
