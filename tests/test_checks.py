from pathlib import Path

import pytest

from gmxprep import checks
from gmxprep.gmx_prep import load_manifest
from gmxprep.pdbio import GroAtom, write_gro


def _setup(tmp_path):
    ff = tmp_path / "toy.ff"
    ff.mkdir()
    (ff / "forcefield.itp").write_text("[ defaults ]\n1 2 yes 0.5 0.8333\n")
    (ff / "tip3p.itp").write_text("[ moleculetype ]\nSOL 2\n[ atoms ]\n1 OW 1 SOL OW 1 -0.834 16\n"
                                  "2 HW 1 SOL HW1 1 0.417 1\n3 HW 1 SOL HW2 1 0.417 1\n")
    (ff / "ions.itp").write_text("[ moleculetype ]\nNA 1\n[ atoms ]\n1 Na 1 NA NA 1 1 23\n"
                                 "[ moleculetype ]\nCL 1\n[ atoms ]\n1 Cl 1 CL CL 1 -1 35\n")
    sysd = tmp_path / "sys"
    sysd.mkdir()
    (sysd / "LIG.itp").write_text("[ moleculetype ]\nLIG 3\n[ atoms ]\n1 c3 1 LIG C1 1 -0.5 12\n2 c3 1 LIG C2 2 -0.5 12\n"
                                  "#ifdef POSRES\n[ position_restraints ]\n1 1 1000 1000 1000\n#endif\n")
    return ff, sysd


def _top(sysd, mols):
    body = "\n".join(f"{n} {c}" for n, c in mols)
    (sysd / "topol.top").write_text('#include "toy.ff/forcefield.itp"\n#include "LIG.itp"\n'
                                    '#include "toy.ff/tip3p.itp"\n#include "toy.ff/ions.itp"\n'
                                    f"[ system ]\nx\n[ molecules ]\n{body}\n")


def _gro(sysd, spec):
    atoms, r = [], 0
    for resname, names in spec:
        r += 1
        atoms += [GroAtom(r, resname, n, 0.1 * len(atoms), 0, 0) for n in names]
    write_gro(sysd / "c.gro", atoms, "t", "   3.0   3.0   3.0")


def test_top_vs_gro_counts_charge_and_index(tmp_path):
    ff, sysd = _setup(tmp_path)
    _top(sysd, [("LIG", 1), ("SOL", 2), ("NA", 1)])
    _gro(sysd, [("LIG", ["C1", "C2"]), ("SOL", ["OW", "HW1", "HW2"]), ("SOL", ["OW", "HW1", "HW2"]), ("NA", ["NA"])])
    r = checks.top_vs_gro(sysd / "topol.top", sysd / "c.gro", [tmp_path])
    assert r["problems"] == [] and r["natoms_top"] == 9 and r["charge"] == pytest.approx(0.0)
    groups = checks.write_index(sysd / "c.gro", sysd / "topol.top", [tmp_path], sysd / "index.ndx", None, ["LIG"])
    assert groups == {"System": 9, "Solute": 2, "Solvent": 7, "LIG_LIG": 2}
    assert "[ Solvent ]" in (sysd / "index.ndx").read_text()


def test_top_vs_gro_detects_mismatch_and_bad_charge(tmp_path):
    ff, sysd = _setup(tmp_path)
    _top(sysd, [("LIG", 1), ("SOL", 3)])  # topology claims 3 waters
    _gro(sysd, [("LIG", ["C1", "C2"]), ("SOL", ["OW", "HW1", "HW2"])])
    r = checks.top_vs_gro(sysd / "topol.top", sysd / "c.gro", [tmp_path])
    assert any("11 atoms" in p for p in r["problems"])
    with pytest.raises(checks.CheckError):
        checks.assert_integer_charge(-0.9, "x")
    assert checks.assert_integer_charge(-1.0004, "x") == -1


def test_log_parsers(tmp_path):
    (tmp_path / "em.log").write_text("Steepest Descents converged to Fmax < 1000 in 812 steps\n"
                                     "Potential Energy  = -3.4567e+05\nMaximum force     =  8.21e+02 on atom 77\n")
    em = checks.parse_em_log(tmp_path / "em.log", 1000)
    assert em["converged"] and em["fmax"] == pytest.approx(821)
    (tmp_path / "bad.log").write_text("Potential Energy  = 1.2e+07\nMaximum force     =  5.0e+02\n")
    assert not checks.parse_em_log(tmp_path / "bad.log", 1000)["converged"]  # positive Epot = overlap
    (tmp_path / "md.log").write_text("Step 10 LINCS WARNING\nFinished mdrun on rank 0\n")
    assert checks.scan_md_log(tmp_path / "md.log") == {"finished": True, "lincs_warnings": 1, "fatal": False}


def test_manifest_variants_from_altloc_branches(tmp_path):
    m = tmp_path / "m.yaml"
    m.write_text("""
enzymes: [{id: E1, pdb: e.pdb}]
targets: [{id: T1, kind: ligand, resname: CID, pdb: t.pdb}]
complexes:
  - {enzyme: E1, target: T1, variant: altA, pdb: cA.pdb}
  - {enzyme: E1, target: T1, variant: altB, pdb: cB.pdb}
""")
    ids = [s.id for s in load_manifest(m, allow_incomplete=False).systems]
    assert ids == ["E1", "T1", "E1__T1__altA", "E1__T1__altB"]
