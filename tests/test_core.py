import pytest

from gmxprep import topology as topo
from gmxprep.forcefields import ForcefieldMismatch, acquire_lock, load_lock, resolve
from gmxprep.gmx_prep import load_manifest, main
from gmxprep.mdp import protocols
from gmxprep.pdbio import Atom, classify

ACPYPE_ITP = """
[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb
 c3       c3          0.00000  0.00000   A     3.39771e-01   4.51035e-01 ; 1.91  0.1078
 hc       hc          0.00000  0.00000   A     2.60018e-01   8.70272e-02 ; 1.46  0.0208

[ moleculetype ]
;name            nrexcl
 template         3

[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass
     1   c3     1   UNL    C1    1    -0.094100     12.01000
     2   hc     1   UNL    H1    2     0.031367      1.00800
     3   hc     1   UNL    H2    3     0.031367      1.00800
     4   hc     1   UNL    H3    4     0.031366      1.00800
"""


def test_split_normalise_posre():
    params, mol = topo.split_params(ACPYPE_ITP)
    assert len(params["atomtypes"]) == 2 and "atomtypes" not in mol
    mol = topo.normalise_moleculetype(mol, "LIG")
    assert topo.moleculetype_name(mol) == "LIG"
    atoms = topo.itp_atoms(mol)
    assert {a["resname"] for a in atoms} == {"LIG"}
    assert abs(sum(a["charge"] for a in atoms)) < 0.01
    assert topo.posre_block(mol).count("1000  1000  1000") == 1  # heavy atoms only


def test_merge_dedupes_and_rejects_conflicts():
    p1, _ = topo.split_params(ACPYPE_ITP)
    merged = topo.merge_params([p1, p1])
    assert merged.count("c3 c3") == 1
    p2 = {"atomtypes": [p1["atomtypes"][0].replace("4.51035e-01", "9.9e-01")]}
    with pytest.raises(topo.TopologyError):
        topo.merge_params([p1, p2])


def test_extract_pdb2gmx_moleculetype():
    top = """; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"

[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
     1         N3      1    MET      N      1     0.1592      14.01

; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include water topology
#include "amber99sb-ildn.ff/tip3p.itp"
"""
    mol = topo.extract_pdb2gmx_moleculetype(top)
    assert topo.moleculetype_name(mol) == "Protein" and "posre.itp" in mol and "tip3p" not in mol


def test_compress_molecules():
    assert topo.compress_molecules([("NAD", 1), ("NAD", 1), ("ZN", 1), ("NAD", 1)]) == [
        ("NAD", 2),
        ("ZN", 1),
        ("NAD", 1),
    ]


def test_lock_prevents_mixing(tmp_path):
    acquire_lock(tmp_path, "amber")
    assert load_lock(tmp_path).key == "amber99sb-ildn"
    assert acquire_lock(tmp_path, "amber99sb-ildn").key == "amber99sb-ildn"
    with pytest.raises(ForcefieldMismatch):
        acquire_lock(tmp_path, "charmm")


def test_cli_lock_mismatch_exit_code(tmp_path, capsys):
    assert main(["lock", "--batch-dir", str(tmp_path), "--forcefield", "charmm"]) == 0
    assert main(["lock", "--batch-dir", str(tmp_path), "--forcefield", "amber"]) == 1
    assert "ForcefieldMismatch" in capsys.readouterr().out.splitlines()[-1]


def test_charmm_mdp_uses_force_switch():
    em = protocols(resolve("charmm"))["em"]
    assert em["vdw-modifier"] == "Force-switch" and em["rvdw"] == "1.2" and em["DispCorr"] == "no"
    assert protocols(resolve("amber"))["md"]["DispCorr"] == "EnerPres"


def _atom(rec, name, res, chain, seq, el):
    return Atom(rec, name, res, chain, seq, "", 0.0, 0.0, 0.0, el)


def test_classify_drops_water_and_h():
    atoms = [
        _atom("ATOM", "N", "ALA", "A", 1, "N"),
        _atom("ATOM", "H", "ALA", "A", 1, "H"),
        _atom("HETATM", "O", "HOH", "A", 2, "O"),
        _atom("HETATM", "ZN", "ZN", "A", 3, "ZN"),
        _atom("HETATM", "C1", "LIG", "A", 4, "C"),
        _atom("HETATM", "C2", "LIG", "A", 4, "C"),
    ]
    s = classify(atoms)
    assert len(s.protein) == 1 and len(s.ions) == 1 and len(s.ligands) == 1 and len(s.ligands[0]) == 2


def test_manifest_requires_full_matrix(tmp_path):
    m = tmp_path / "m.yaml"
    m.write_text("""
enzymes: [{id: E1, pdb: e1.pdb}, {id: E2, pdb: e2.pdb}]
targets: [{id: T1, kind: ligand, resname: LIG, pdb: t1.pdb}]
complexes: [{enzyme: E1, target: T1, pdb: c.pdb}]
""")
    with pytest.raises(ValueError, match="E2__T1"):
        load_manifest(m, allow_incomplete=False)
    man = load_manifest(m, allow_incomplete=True)
    assert [s.id for s in man.systems] == ["E1", "E2", "T1", "E1__T1"]


def test_pose_maps_onto_template_order(tmp_path):
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from gmxprep.ligands import LigandParams, LigandSpec, build_template, pose_coords

    # acetate pose as heavy-atom PDB (no H, no bond orders), shuffled atom order for the 2nd pose
    ref = Chem.AddHs(Chem.MolFromSmiles("CC(=O)[O-]"))
    AllChem.EmbedMolecule(ref, randomSeed=1)
    heavy = [a for a in ref.GetAtoms() if a.GetAtomicNum() > 1]

    def block(order):
        lines = []
        for i, a in enumerate(order, 1):
            p = ref.GetConformer().GetAtomPosition(a.GetIdx())
            lines.append(
                f"HETATM{i:5d}  {a.GetSymbol()}{i:<2d} ACT A   1    {p.x:8.3f}{p.y:8.3f}{p.z:8.3f}"
                f"  1.00  0.00           {a.GetSymbol()}"
            )
        return "\n".join(lines) + "\nEND\n"

    spec = LigandSpec("ACT", smiles="CC(=O)[O-]", charge=-1)
    mol, warns = build_template(spec, block(heavy), tmp_path, False, tmp_path / "log")
    assert mol.GetNumAtoms() == 7 and Chem.GetFormalCharge(mol) == -1 and not warns
    names = [a.GetProp("_Name") for a in mol.GetAtoms()]
    assert len(set(names)) == 7
    Chem.MolToMolFile(mol, str(tmp_path / "template.sdf"))
    params = LigandParams("ACT", tmp_path, tmp_path / "ACT.itp", {}, names, -1, [])
    xyz = pose_coords(params, block(list(reversed(heavy))))
    assert len(xyz) == 7
    # heavy atoms must land on the same coordinates regardless of input order
    tmpl = mol.GetConformer()
    for i, a in enumerate(mol.GetAtoms()):
        if a.GetAtomicNum() > 1:
            assert all(abs(c - t / 10) < 2e-3 for c, t in zip(xyz[i], tmpl.GetAtomPosition(i)))
