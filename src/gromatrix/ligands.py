"""Ligand / cofactor parameterisation, routed by the batch forcefield.

    AMBER family  -> acpype (antechamber, GAFF2, AM1-BCC)
    CHARMM family -> CGenFF (.str from `cgenff` or user-supplied) -> cgenff_charmm2gmx.py

One *template* per residue name is built and parameterised once per batch and
cached in <batch>/06-params/<RES>/. Every pose of that residue in every system is
mapped onto the template's atom order, so the enzyme-alone, target-alone and
complex systems share byte-identical ligand topologies.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

from . import topology as topo
from .forcefields import ForceField, assert_same
from .gmx import ToolError, run
from .result import CodedError

ORGANIC = {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"}
CGENFF_PENALTY_LIMIT = 50.0


class LigandError(CodedError):
    code = "ligand_invalid"


@dataclass
class LigandSpec:
    resname: str
    smiles: str | None = None
    charge: int | None = None
    multiplicity: int = 1
    sdf: str | None = None  # curated template with explicit H (authoritative)
    str_file: str | None = None  # pre-computed CGenFF stream file

    @classmethod
    def from_manifest(cls, resname: str, d: dict | None, root: Path) -> LigandSpec:
        d = d or {}

        def path(k):
            return str((root / d[k]).resolve()) if d.get(k) else None

        return cls(resname, d.get("smiles"), d.get("charge"), d.get("multiplicity", 1), path("sdf"), path("str"))


@dataclass
class LigandParams:
    resname: str
    directory: Path
    itp: Path  # moleculetype only
    params: dict  # type sections, see topology.split_params
    atom_names: list[str]
    net_charge: float
    warnings: list[str]


# ------------------------------------------------------------------ RDKit
def _rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as e:  # pragma: no cover
        raise LigandError("RDKit is required for ligand handling (conda install -c conda-forge rdkit)") from e
    return Chem, AllChem


def _pose_from_pdb_block(block: str):
    Chem, _ = _rdkit()
    mol = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False, proximityBonding=True)
    if mol is None:
        raise LigandError("RDKit could not read ligand pose")
    return mol


def _assign_from_reference(ref_heavy, pose_heavy):
    """Copy bond orders and formal charges from `ref_heavy` onto a bond-order-less pose."""
    Chem, AllChem = _rdkit()
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(ref_heavy, pose_heavy)
    except ValueError as e:
        raise LigandError(f"pose does not match reference chemistry: {e}") from e
    match = mol.GetSubstructMatch(ref_heavy)
    for ti, pi in enumerate(match):
        a = mol.GetAtomWithIdx(pi)
        a.SetFormalCharge(ref_heavy.GetAtomWithIdx(ti).GetFormalCharge())
        a.SetNoImplicit(False)
        a.SetNumExplicitHs(0)
    Chem.SanitizeMol(mol)
    return mol


def _heavy(mol):
    Chem, _ = _rdkit()
    return Chem.RemoveHs(mol, sanitize=False)


def build_template(
    spec: LigandSpec, pose_block: str | None, workdir: Path, allow_auto: bool, log: Path
) -> tuple[object, list[str]]:
    """Return (RDKit mol with explicit H and 3D coords, warnings)."""
    Chem, AllChem = _rdkit()
    warnings: list[str] = []
    if spec.sdf:
        mol = Chem.MolFromMolFile(spec.sdf, removeHs=False)
        if mol is None:
            raise LigandError(f"cannot read {spec.sdf}")
    elif spec.smiles:
        ref = Chem.MolFromSmiles(spec.smiles)
        if ref is None:
            raise LigandError(f"bad SMILES for {spec.resname}")
        if pose_block:
            mol = Chem.AddHs(_assign_from_reference(ref, _heavy(_pose_from_pdb_block(pose_block))), addCoords=True)
        else:
            mol = Chem.AddHs(ref)
            AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
            AllChem.MMFFOptimizeMolecule(mol)
    elif allow_auto and pose_block:
        warnings.append(f"{spec.resname}: protonation guessed by Open Babel at pH 7.4 - verify!")
        (workdir / "pose.pdb").write_text(pose_block)
        run(["obabel", "pose.pdb", "-O", "auto.sdf", "-p", "7.4"], workdir, log)
        mol = Chem.MolFromMolFile(str(workdir / "auto.sdf"), removeHs=False)
    else:
        raise LigandError(
            f"no SMILES/SDF for '{spec.resname}' in manifest `ligands:`; add one (or pass --allow-auto-protonation)"
        )

    bad = {a.GetSymbol() for a in mol.GetAtoms()} - ORGANIC
    if bad:
        raise LigandError(
            f"{spec.resname} contains {sorted(bad)}; GAFF2/CGenFF cannot describe "
            "metal centres - use MCPB.py (AMBER) or a bonded/cationic-dummy model"
        )
    charge = Chem.GetFormalCharge(mol)
    if spec.charge is not None and charge != spec.charge:
        raise LigandError(f"{spec.resname}: template formal charge {charge} != manifest charge {spec.charge}")

    # Unique <=4-char atom names; keep crystallographic names for heavy atoms when available.
    pdb_names = []
    if pose_block and not spec.sdf:
        pdb_names = [l[12:16].strip() for l in pose_block.splitlines() if l.startswith(("ATOM", "HETATM"))]
    used, counters = set(), {}
    for a in mol.GetAtoms():
        name = None
        info = a.GetPDBResidueInfo()
        if info is not None:
            name = info.GetName().strip()
        elif a.GetAtomicNum() > 1 and a.GetIdx() < len(pdb_names):
            name = pdb_names[a.GetIdx()]
        while not name or name in used or len(name) > 4:
            counters[a.GetSymbol()] = counters.get(a.GetSymbol(), 0) + 1
            name = f"{a.GetSymbol().upper()}{counters[a.GetSymbol()]}"
        used.add(name)
        a.SetProp("_Name", name)
    return mol, warnings


def _write_mol2(mol, resname: str, workdir: Path, log: Path) -> Path:
    Chem, _ = _rdkit()
    Chem.MolToMolFile(mol, str(workdir / "template.sdf"))
    run(["obabel", "template.sdf", "-O", "raw.mol2"], workdir, log)
    names = [a.GetProp("_Name") for a in mol.GetAtoms()]
    out, section, i = [], None, 0
    lines = (workdir / "raw.mol2").read_text().splitlines()
    for n, line in enumerate(lines):
        if line.startswith("@<TRIPOS>"):
            section = line
            out.append(line)
            continue
        if section == "@<TRIPOS>MOLECULE" and lines[n - 1] == section:
            line = resname
        elif section == "@<TRIPOS>ATOM" and line.strip():
            f = line.split()
            f[1], f[7] = names[i], resname
            i += 1
            line = (
                f"{int(f[0]):7d} {f[1]:<4s} {float(f[2]):10.4f} {float(f[3]):10.4f} {float(f[4]):10.4f} "
                f"{f[5]:<6s} {f[6]:>3s} {f[7]:<4s} {float(f[8]):9.4f}"
            )
        elif section == "@<TRIPOS>SUBSTRUCTURE" and line.strip():
            f = line.split()
            f[1] = resname
            line = " ".join(f)
        out.append(line)
    if i != mol.GetNumAtoms():
        raise LigandError("mol2 atom count mismatch after Open Babel conversion")
    path = workdir / "template.mol2"
    path.write_text("\n".join(out) + "\n")
    return path


# ------------------------------------------------------------- engines
def _run_acpype(spec: LigandSpec, ff: ForceField, mol2: Path, charge: int, log: Path) -> str:
    wd = mol2.parent
    run(
        [
            "acpype",
            "-i",
            mol2.name,
            "-b",
            spec.resname,
            "-c",
            ff.charge_method,
            "-a",
            ff.ligand_ff,
            "-n",
            str(charge),
            "-m",
            str(spec.multiplicity),
            "-o",
            "gmx",
        ],
        wd,
        log,
    )
    itp = wd / f"{spec.resname}.acpype" / f"{spec.resname}_GMX.itp"
    if not itp.exists():
        raise ToolError(f"acpype produced no {itp}")
    return itp.read_text()


def _cgenff_penalties(str_text: str) -> tuple[float, float] | None:
    """(param, charge) penalties from a CGenFF stream file, or None if not found."""
    m = re.search(r"param\s+penalty\s*=\s*([\d.]+)\s*;\s*charge\s+penalty\s*=\s*([\d.]+)", str_text, re.I)
    return (float(m.group(1)), float(m.group(2))) if m else None


def _run_cgenff(
    spec: LigandSpec, ff_dir: Path, mol2: Path, log: Path, accept_high_penalty: bool, warnings: list[str]
) -> str:
    wd = mol2.parent
    str_path = wd / f"{spec.resname}.str"
    if spec.str_file:
        str_path.write_text(Path(spec.str_file).read_text())
    else:
        str_path.write_text(_cgenff_stdout(os.environ.get("CGENFF_BIN", "cgenff"), mol2, wd, log))
    penalties = _cgenff_penalties(str_path.read_text())
    if penalties is None:
        # An unparseable header must not silently pass the quality gate.
        if not accept_high_penalty:
            raise LigandError(
                f"{spec.resname}: no 'param penalty= ; charge penalty=' line in {str_path.name}; "
                "cannot assess CGenFF quality (pass --accept-high-penalty to proceed anyway)"
            )
        warnings.append(f"{spec.resname}: CGenFF penalties unknown (unparseable .str header)")
        penalties = (0.0, 0.0)
    pp, cp = penalties
    if max(pp, cp) > CGENFF_PENALTY_LIMIT and not accept_high_penalty:
        raise LigandError(
            f"{spec.resname}: CGenFF penalties param={pp} charge={cp} > "
            f"{CGENFF_PENALTY_LIMIT}; reparameterise (ffTK) or pass --accept-high-penalty"
        )
    if max(pp, cp) > 10:
        warnings.append(f"{spec.resname}: CGenFF penalties param={pp} charge={cp} - validate")
    conv = os.environ.get("CGENFF2GMX", "cgenff_charmm2gmx_py3_nx2.py")
    # Run from 06-params/<RES>, so a bare script name must be resolved on PATH first.
    conv = shutil.which(conv) or conv
    if not Path(conv).is_file():
        raise ToolError(f"CGenFF converter '{conv}' not found; set CGENFF2GMX to its absolute path")
    run([sys.executable, conv, spec.resname, mol2.name, str_path.name, str(ff_dir)], wd, log)
    stem = spec.resname.lower()
    itp, prm = wd / f"{stem}.itp", wd / f"{stem}.prm"
    if not itp.exists():
        raise ToolError(f"cgenff_charmm2gmx produced no {itp}")
    return (prm.read_text() if prm.exists() else "") + "\n" + itp.read_text()


def _cgenff_stdout(exe: str, mol2: Path, wd: Path, log: Path) -> str:
    import subprocess

    proc = subprocess.run([exe, mol2.name], cwd=wd, capture_output=True, text=True)
    with open(log, "a") as fh:
        fh.write(f"\n$ {exe} {mol2.name}\n{proc.stderr}")
    if proc.returncode != 0 or "RESI" not in proc.stdout:
        raise ToolError(f"cgenff failed for {mol2}; supply a server-generated .str via `str:` in the manifest")
    return proc.stdout


# ------------------------------------------------------------- public API
def parameterise(
    spec: LigandSpec,
    pose_block: str | None,
    ff: ForceField,
    ff_dir: Path,
    params_root: Path,
    allow_auto: bool = False,
    accept_high_penalty: bool = False,
    refresh: bool = False,
) -> LigandParams:
    wd = params_root / spec.resname
    meta_path = wd / "meta.json"
    fingerprint = _spec_fingerprint(spec)
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        assert_same(ff, meta["forcefield"], f"cached parameters {wd}")
        if meta.get("spec_fingerprint") != fingerprint:
            if not refresh:
                raise LigandError(
                    f"{spec.resname}: the manifest's SMILES/charge/SDF/str changed since {wd} was built; "
                    "rebuild with --force to re-parameterise (every system using it is rebuilt)"
                )
            shutil.rmtree(wd)
        else:
            return _load(wd, meta)

    wd.mkdir(parents=True, exist_ok=True)
    log = wd / "param.log"
    mol, warnings = build_template(spec, pose_block, wd, allow_auto, log)
    Chem, _ = _rdkit()
    charge = Chem.GetFormalCharge(mol)
    mol2 = _write_mol2(mol, spec.resname, wd, log)

    if ff.ligand_engine == "acpype":
        raw = _run_acpype(spec, ff, mol2, charge, log)
    elif ff.ligand_engine == "cgenff":
        raw = _run_cgenff(spec, ff_dir, mol2, log, accept_high_penalty, warnings)
    else:
        raise LigandError(f"no ligand engine '{ff.ligand_engine}'")

    params, mol_text = topo.split_params(raw)
    mol_text = topo.normalise_moleculetype(mol_text, spec.resname)
    atoms = topo.itp_atoms(mol_text)
    if len(atoms) != mol.GetNumAtoms():
        raise LigandError(f"{spec.resname}: itp has {len(atoms)} atoms, template {mol.GetNumAtoms()}")
    # Poses are written in template order with itp names, so the orders must agree exactly
    # (a server-generated .str can list atoms in a different order than our mol2).
    template_names = [a.GetProp("_Name") for a in mol.GetAtoms()]
    if [a["name"] for a in atoms] != template_names:
        raise LigandError(
            f"{spec.resname}: itp atom order differs from the template; regenerate the .str "
            "from 06-params/<RES>/template.mol2 so coordinates map onto the right atoms"
        )
    q = sum(a["charge"] for a in atoms)
    if abs(q - charge) > 1e-3:  # same tolerance as step F's integer-charge gate
        raise LigandError(f"{spec.resname}: itp net charge {q:.3f} != formal charge {charge}")
    (wd / f"{spec.resname}.itp").write_text(mol_text + topo.posre_block(mol_text))
    (wd / "params.json").write_text(json.dumps(params))

    meta = {
        "resname": spec.resname,
        "forcefield": ff.identity(),
        "engine": ff.ligand_engine,
        "net_charge": charge,
        "atom_names": [a["name"] for a in atoms],
        "warnings": warnings,
        "template_sha256": hashlib.sha256((wd / "template.sdf").read_bytes()).hexdigest(),
        "spec": spec.__dict__,
        "spec_fingerprint": fingerprint,
    }
    meta_path.write_text(json.dumps(meta, indent=2))
    return _load(wd, meta)


def _spec_fingerprint(spec: LigandSpec) -> str:
    """Hash of everything that defines a ligand's parameters, including file contents."""
    h = hashlib.sha256(
        json.dumps(
            {k: v for k, v in spec.__dict__.items() if k not in ("sdf", "str_file")}, sort_keys=True, default=str
        ).encode()
    )
    for f in (spec.sdf, spec.str_file):
        h.update(Path(f).read_bytes() if f and Path(f).exists() else b"-")
    return h.hexdigest()


def _load(wd: Path, meta: dict) -> LigandParams:
    return LigandParams(
        meta["resname"],
        wd,
        wd / f"{meta['resname']}.itp",
        json.loads((wd / "params.json").read_text()),
        meta["atom_names"],
        meta["net_charge"],
        meta["warnings"],
    )


def pose_coords(params: LigandParams, pose_block: str | None = None, pose_sdf: str | None = None):
    """Coordinates (nm) of one pose, in the template's (== itp) atom order."""
    Chem, _ = _rdkit()
    template = Chem.MolFromMolFile(str(params.directory / "template.sdf"), removeHs=False)
    if pose_sdf:
        pose = Chem.MolFromMolFile(pose_sdf, removeHs=False)
        match = pose.GetSubstructMatch(template) if pose.GetNumAtoms() == template.GetNumAtoms() else ()
        if match:
            conf = pose.GetConformer()
            return [tuple(c / 10 for c in conf.GetAtomPosition(i)) for i in match]
        pose_heavy = _heavy(pose)
    elif pose_block:
        pose_heavy = _heavy(_pose_from_pdb_block(pose_block))
    else:
        conf = template.GetConformer()
        return [tuple(c / 10 for c in conf.GetAtomPosition(i)) for i in range(template.GetNumAtoms())]
    full = Chem.AddHs(_assign_from_reference(_heavy(template), pose_heavy), addCoords=True)
    match = full.GetSubstructMatch(template)
    if len(match) != template.GetNumAtoms():
        raise LigandError(
            f"{params.resname}: pose could not be mapped onto parameter template (different tautomer/protonation?)"
        )
    conf = full.GetConformer()
    return [tuple(c / 10 for c in conf.GetAtomPosition(i)) for i in match]
