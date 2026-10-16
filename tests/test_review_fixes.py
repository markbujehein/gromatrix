"""Regression tests for defects found by the code-review agents (one test per fix)."""

import numpy as np
import pytest

from gromatrix import altlocs
from gromatrix.afdb import lddt_ca
from gromatrix.fix import openmm_uses_label_ids
from gromatrix.gmx_prep import load_manifest
from gromatrix.net import _retry_after
from gromatrix.pdbio import Atom, classify, read_pdb
from gromatrix.topology import TopologyError, merge_params


def test_chain_remap_mirrors_openmm_rule():
    # ligands/waters with their own label chains -> OpenMM reads label IDs -> remap needed
    assert openmm_uses_label_ids({"label2auth": {"A": "A", "B": "B", "C": "A", "D": "B"}})
    # polymer-only file with swapped label/auth IDs -> OpenMM already uses auth IDs -> no remap
    assert not openmm_uses_label_ids({"label2auth": {"A": "B", "B": "A"}})
    assert not openmm_uses_label_ids({})


def test_blank_occupancy_and_bfactor_columns(tmp_path):
    line = "ATOM      1  CA  ALA A   1      11.104   6.134  -6.504" + " " * 12 + "           C"
    p = tmp_path / "blank.pdb"
    p.write_text(line + "\n")
    ((atom, alt),) = altlocs.parse_all_atoms(p)
    assert atom.occ == 1.0 and atom.b == 0.0 and alt == ""


def test_microheterogeneity_is_one_position_with_two_conformers():
    def a(name, resname, alt, occ, x):
        return (Atom("ATOM", name, resname, "A", 5, "", x, 0.0, 0.0, name[0], 10.0, occ), alt)

    recs = [
        a("CA", "SER", "A", 0.4, 0.0),
        a("OG", "SER", "A", 0.4, 1.0),
        a("CA", "THR", "B", 0.6, 0.1),
        a("OG1", "THR", "B", 0.6, 1.1),
        a("CG2", "THR", "B", 0.6, 2.0),
    ]
    atoms, groups = altlocs.resolve(recs)
    assert len(groups) == 1 and groups[0]["microheterogeneity"] and groups[0]["primary"] == "B"
    assert {x.resname for x in atoms} == {"THR"} and len(atoms) == 3


def test_four_char_resname_round_trip_and_overflow(tmp_path):
    p = tmp_path / "w.pdb"
    p.write_text(Atom("HETATM", "OH2", "TIP3", "W", 12, "", 1, 2, 3, "O").pdb_line(1) + "\n")
    (w,) = read_pdb(p)
    assert (w.resname, w.chain, w.resseq) == ("TIP3", "W", 12)
    with pytest.raises(ValueError):
        Atom("ATOM", "CA", "ALA", "A", 10000, "", 0, 0, 0, "C").pdb_line(1)


def test_bonded_parameter_conflicts_are_detected_direction_independently():
    a = {"dihedraltypes": ["c3 ca ca c3 9 180.0 10.0 2", "c3 ca ca c3 9 0.0 1.0 3"], "bondtypes": ["c3 ca 1 0.15 3e5"]}
    same = {
        "dihedraltypes": ["c3 ca ca c3 9 0.0 1.0 3", "c3 ca ca c3 9 180.0 10.0 2"],
        "bondtypes": ["ca c3 1 0.15 3e5"],
    }
    other_terms = {"dihedraltypes": ["c3 ca ca c3 9 180.0 10.0 2", "c3 ca ca c3 9 0.0 0.5 1"]}
    assert merge_params([a, same]).count("c3 ca ca c3") == 2
    with pytest.raises(TopologyError):
        merge_params([a, other_terms])


def test_lddt_penalises_unmodelled_reference_residues():
    t = np.arange(60) * np.radians(100)
    ref = np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(60)]
    half = ref.copy()
    half[30:] = np.nan
    assert lddt_ca(ref, ref) == pytest.approx(1.0) and lddt_ca(half, ref) < 0.6


def test_unselected_protein_chains_are_dropped_not_ligands():
    atoms = [Atom("ATOM", "CA", "ALA", c, 1, "", 0, 0, 0, "C") for c in "AB"]
    split = classify(atoms, {"A"})
    assert len(split.protein) == 1 and split.ligands == []


def test_duplicate_ligand_target_resnames_rejected(tmp_path):
    m = tmp_path / "m.yaml"
    m.write_text("""
enzymes: [{id: E1, pdb: e.pdb}]
targets: [{id: T1, kind: ligand, resname: LIG, pdb: a.pdb}, {id: T2, kind: ligand, resname: LIG, pdb: b.pdb}]
complexes: [{enzyme: E1, target: T1, pdb: c1.pdb}, {enzyme: E1, target: T2, pdb: c2.pdb}]
""")
    with pytest.raises(ValueError, match="LIG"):
        load_manifest(m, allow_incomplete=False)


def test_retry_after_parsing_is_bounded():
    assert _retry_after("5", 2) == 5 and _retry_after("3600", 2) == 120 and _retry_after("junk", 2) == 2
    assert _retry_after("Wed, 21 Oct 2015 07:28:00 GMT", 2) == 0  # a date in the past
