"""Steps B/C. The 1XVP benchmark needs OpenMM + PDBFixer; the rest is pure NumPy/SciPy."""
import argparse
import gzip
import json
import shutil
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from gmxprep.extract import assign_hets, au_units, build_assembly, chain_contacts, extract_one, load_au
from gmxprep.pdbio import Atom, write_pdb
from gmxprep.structure import compose, is_identity, parse_oper_expression

DATA = Path(__file__).parent / "data"


def test_oper_expression():
    assert parse_oper_expression("1") == [("1",)]
    assert parse_oper_expression("(1-3)") == [("1",), ("2",), ("3",)]
    assert parse_oper_expression("(1,2)(5)") == [("1", "5"), ("2", "5")]


def test_compose_applies_right_to_left():
    ops = {"1": {"R": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "t": [1, 0, 0]},
           "2": {"R": [[0, -1, 0], [1, 0, 0], [0, 0, 1]], "t": [0, 0, 0]}}  # 90° about z
    op = compose(ops, ("1", "2"))  # rotate, then translate
    x = np.array(op["R"]) @ np.array([1.0, 0, 0]) + np.array(op["t"])
    assert np.allclose(x, [1, 1, 0]) and not is_identity(op)


def _chain(chain, x0, n=4, resname="ALA"):
    """A straight poly-Ala strand; consecutive C-N distance 1.33 Å."""
    atoms = []
    for i in range(n):
        x = x0 + 3.8 * i
        for name, dx, dy in (("N", 0.0, 0.0), ("CA", 1.2, 0.5), ("C", 2.47, 0.0), ("O", 2.5, 1.2), ("CB", 1.2, 1.9)):
            atoms.append(Atom("ATOM", name, resname, chain, i + 1, "", x + dx, dy, 0.0, name[0], 20.0, 1.0))
    return atoms


def _lig(chain, x, y, resname="LIG", seq=500):
    return [Atom("HETATM", f"C{k}", resname, chain, seq, "", x + 1.4 * k, y, 0.0, "C", 30.0, 1.0) for k in range(3)]


def replace_y(atom, dy):
    return replace(atom, y=atom.y + dy)


def _residues(atoms):
    with tempfile.TemporaryDirectory() as d:
        write_pdb(Path(d) / "x.pdb", atoms)
        return load_au(Path(d) / "x.pdb")


def test_ligand_only_in_second_copy_selects_second_copy(tmp_path):
    # Two copies of a receptor (entity 1: A,C) + partner (entity 2: B,D); ligand bound to D only.
    atoms = (_chain("A", 0) + _chain("C", 100) + [replace_y(a, 4.0) for a in _chain("B", 0)]
             + [replace_y(a, 4.0) for a in _chain("D", 100)] + _lig("D", 103, 8.0))
    fixed = tmp_path / "fixed"
    fixed.mkdir()
    write_pdb(fixed / "TOY.pdb", atoms)
    (fixed / "TOY.meta.json").write_text(json.dumps(
        {"assemblies": [], "operators": {}, "chain_entity": {"A": "1", "C": "1", "B": "2", "D": "2"}}))
    a = argparse.Namespace(targets="LIG", target_chains=None, expect=None, assembly="auto", ligand_cutoff=4.5,
                           contact_cutoff=5.0, min_contacts=2, max_chains=12, unit=None, allow_incomplete=False,
                           altloc="auto")
    out = extract_one(fixed / "TOY.pdb", tmp_path / "units", a)
    assert out["selected"] == "au-CD" and out["inventory"] == {"LIG": 1}
    files = set(out["files"])
    assert {"enzyme.pdb", "target_LIG.pdb", "complex_LIG.pdb"} <= files
    enzyme = (tmp_path / "units" / "TOY" / "enzyme.pdb").read_text()
    assert "LIG" not in enzyme and "ALA C" in enzyme and "ALA A" not in enzyme


def test_assign_hets_prefers_most_contacts():
    res = _residues(_chain("A", 0) + [replace_y(a, 10.0) for a in _chain("B", 0)] + _lig("B", 1.0, 7.5))
    assignment, orphans = assign_hets(res, 4.5)
    (i, info), = assignment.items()
    assert info["chain"] == "B" and not orphans


def test_assembly_copies_get_fresh_chain_ids():
    au = _residues(_chain("A", 0) + _lig("A", 1.0, 3.0))
    asm = {"id": "1", "author_defined": True,
           "gens": [{"polymer_chains": ["A"], "het_residues": [["A", "LIG", 500, ""]], "opers": [("1",), ("2",)]}]}
    ops = {"1": {"R": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "t": [0, 0, 0]},
           "2": {"R": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "t": [0, 50, 0]}}
    cand = build_assembly(asm, au, ops)
    assert cand.chains == ["A", "B"]
    lig = [r for r in cand.residues if r.resname == "LIG"]
    assert sorted(r.chain for r in lig) == ["A", "B"]  # ligand copy follows its chain copy


def test_au_units_require_connected_interfaces():
    au = _residues(_chain("A", 0) + [replace_y(a, 4.0) for a in _chain("B", 0)]
                   + _chain("C", 100) + [replace_y(a, 4.0) for a in _chain("D", 100)])
    contacts = chain_contacts(au, 5.0)
    units = au_units(au, {"A": "1", "C": "1", "B": "2", "D": "2"}, contacts, 2)
    assert sorted(units) == [("A", "B"), ("C", "D")]  # A+D and C+B are not in contact


@pytest.mark.slow
def test_1xvp_benchmark(tmp_path):
    pytest.importorskip("pdbfixer")
    from gmxprep import extract, fix

    raw = tmp_path / "raw"
    raw.mkdir()
    with gzip.open(DATA / "1XVP.cif.gz", "rb") as src, open(raw / "1XVP.cif", "wb") as dst:
        shutil.copyfileobj(src, dst)
    assert fix.main(["--batch-dir", str(tmp_path)]) == 0
    assert extract.main(["--batch-dir", str(tmp_path), "--targets", "CID"]) == 0
    rep = json.loads((tmp_path / "units" / "1XVP" / "unit.json").read_text())
    assert rep["selected"] == "assembly-2" and rep["chains"] == ["C", "D", "G", "H"]
    assert rep["metrics"]["inventory"] == {"F15": 1, "CID": 1}
    a1 = next(c for c in rep["candidates"] if c["id"] == "assembly-1")
    assert not a1["complete"] and a1["missing"] == {"CID": 1}
    enzyme = (tmp_path / "units" / "1XVP" / "enzyme.pdb").read_text()
    assert " CID " not in enzyme and " F15 C" in enzyme
    assert extract.main(["--batch-dir", str(tmp_path), "--out", "u2", "--targets", "CID", "--assembly", "1"]) == 1

    # Alternate conformers: CITCO A/B are a chlorophenyl ring flip -> equivalent, no branch by default ...
    meta = json.loads((tmp_path / "fixed" / "1XVP.meta.json").read_text())
    (cid,) = [g for g in meta["altlocs"] if g["resname"] == "CID"]
    assert cid["analysis"]["B"]["rmsd"] > 1.0 and cid["analysis"]["B"]["sym_rmsd"] < 0.5 and not cid["distinct"]
    assert rep["altlocs"]["CID"]["branch"] is False and "complex_CID" in rep["files"]
    # ... but an explicit request branches the matrix into two clash-free complexes.
    assert extract.main(["--batch-dir", str(tmp_path), "--out", "u3", "--targets", "CID", "--altloc", "branch"]) == 0
    rep3 = json.loads((tmp_path / "u3" / "1XVP" / "unit.json").read_text())
    assert {"complex_CID_altA", "complex_CID_altB"} <= set(rep3["files"])
    assert all(b["clash"]["clashes"] == 0 for b in rep3["altlocs"]["CID"]["branches"].values())


def _alt_group(xa, xb, elements, names):
    return {"chain": "A", "resname": "LIG", "resseq": 500, "icode": "", "primary": "A",
            "alts": {"A": {"occ": 0.5, "atoms": dict(zip(names, xa.tolist())), "elements": dict(zip(names, elements))},
                     "B": {"occ": 0.5, "atoms": dict(zip(names, xb.tolist())), "elements": dict(zip(names, elements))}}}


def test_altloc_policy_ring_flip_shift_and_average():
    from gmxprep import altlocs
    ang = np.radians(np.arange(6) * 60)
    ring = np.c_[1.39 * np.cos(ang), 1.39 * np.sin(ang), np.zeros(6)]
    names, el = [f"C{i}" for i in range(6)], ["C"] * 6
    flipped = ring[[0, 5, 4, 3, 2, 1]]                     # same positions, permuted labels
    g = altlocs.analyse(_alt_group(ring, flipped, el, names), 0.5)
    assert g["analysis"]["B"]["rmsd"] > 1.0 and g["analysis"]["B"]["sym_rmsd"] < 1e-6 and not g["distinct"]
    shifted = ring + [2.0, 0, 0]                           # a genuinely different pose
    g = altlocs.analyse(_alt_group(ring, shifted, el, names), 0.5)
    assert g["distinct"] and g["analysis"]["B"]["geometry_ok"]
    averaged = (ring + flipped) / 2                        # what "conformer C = mean(A,B)" would be
    g = altlocs.analyse(_alt_group(ring, averaged, el, names), 0.5)
    assert not g["analysis"]["B"]["geometry_ok"] and not g["distinct"]


def test_resolve_picks_highest_occupancy():
    from gmxprep import altlocs
    a = lambda occ, x: Atom("HETATM", "C1", "LIG", "A", 1, "", x, 0.0, 0.0, "C", 10.0, occ)
    atoms, groups = altlocs.resolve([(a(0.3, 0.0), "A"), (a(0.7, 3.0), "B")])
    assert groups[0]["primary"] == "B" and atoms[0].x == 3.0
