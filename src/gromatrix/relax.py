"""[D] OpenMM structural-integrity relaxation of the extracted matrix files.

    gromatrix-relax --batch-dir B [--in units] [--out relaxed] [--ids 1ABC ...]
                  [--files enzyme.pdb complex_*.pdb target_*.pdb] [--forcefield amber|charmm]
                  [--solvent auto|gbn2|rf] [--restraint-k 1000] [--max-iter 2000] [--free-iter 1000]
                  [--max-rmsd 2.5] [--max-drop 100] [--pocket-cutoff 5.0] [--platform auto|CPU|CUDA|OpenCL]

This is an integrity check, not equilibration. Per file:
  1. Polymer chains get missing terminal atoms (OXT at TER breaks) and
     hydrogens at --ph, using the locked forcefield family (amber14-all or
     charmm36). Solvent model:
       gbn2  GBn2 implicit solvent (AMBER only). About 0.7 s per evaluation for
             8k atoms on 4 CPU cores, so the default on a GPU only.
       rf    no explicit solvent; 1.6 nm non-periodic cutoff, which OpenMM treats
             with a reaction field (eps 78.3). About 25 ms per evaluation; CPU default.
     HETs and ions are *not*
     parameterised here (that happens once, in step F). They are held out of
     the system and written back unchanged.
  2. Stage 1: minimise with harmonic restraints on all heavy atoms (--restraint-k,
     kJ/mol/nm^2), which removes hydrogen and side-chain strain.
     Stage 2: minimise with restraints only on residues within --pocket-cutoff
     of a HET, so the empty ligand site cannot collapse while the rest relaxes freely.
  3. Gates. Each failure is reported with a label:
       explosion            non-finite energy or coordinates
       implosion/explosion  heavy-atom RMSD > --max-rmsd AND a dramatic energy drop
                            (more negative than -(--max-drop) kJ/mol per atom)
       distortion           heavy-atom RMSD > --max-rmsd (moderate energy change)
       chain break          a Calpha-Calpha distance stretched beyond 4.2 Å
       unresolved clash     final potential energy still positive
       pocket collapse      relaxed protein within 2.0 Å of a held-out HET
Relaxed files keep the input name under <out>/<ID>/. Hydrogens are stripped
again, because step F rebuilds them with pdb2gmx.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import math
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np

from . import layout
from .forcefields import ForcefieldMismatch, acquire_lock, lock_path, resolve
from .pdbio import ION_RES, PROTEIN_RES, Atom, read_pdb, write_pdb

FF_FILES = {
    "amber": {"gbn2": ["amber14-all.xml", "implicit/gbn2.xml"], "rf": ["amber14-all.xml"]},
    "charmm": {"rf": ["charmm36.xml"]},
}
CUTOFF_NM = 1.6
CA_BREAK = 4.2


def _kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    a0, b0 = a - a.mean(0), b - b.mean(0)
    u, s, vt = np.linalg.svd(a0.T @ b0)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    return float(np.sqrt((((a0 @ rot) - b0) ** 2).sum(1).mean()))


def _platform(name: str):
    import openmm

    if name != "auto":
        return openmm.Platform.getPlatformByName(name)
    for n in ("CUDA", "OpenCL", "CPU", "Reference"):
        try:
            return openmm.Platform.getPlatformByName(n)
        except Exception:
            continue


def _ca_breaks(atoms: list[Atom]) -> set[tuple]:
    ca = [a for a in atoms if a.name == "CA" and a.resname in PROTEIN_RES]
    out = set()
    for p, q in zip(ca, ca[1:]):
        if p.chain == q.chain and math.dist((p.x, p.y, p.z), (q.x, q.y, q.z)) > CA_BREAK:
            out.add((p.chain, p.resseq, q.resseq))
    return out


def relax_file(path: Path, out_path: Path, family: str, a) -> dict:
    import openmm
    import openmm.app as app
    from openmm import unit
    from pdbfixer import PDBFixer
    from scipy.spatial import cKDTree

    atoms = read_pdb(path)
    poly = [x for x in atoms if x.resname in PROTEIN_RES and x.element not in ("H", "D")]
    hets = [x for x in atoms if x.resname not in PROTEIN_RES]
    if not poly:
        return {"status": "skipped", "reason": "no polymer (ligand geometry is checked in step F)"}
    platform = _platform(a.platform)
    solvent = a.solvent
    if solvent == "auto":
        solvent = "gbn2" if platform.getName() in ("CUDA", "OpenCL") else "rf"
    if family == "charmm":
        solvent = "rf"  # OpenMM's GB parameter files are AMBER-specific
    ff = app.ForceField(*FF_FILES[family][solvent])

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "in.pdb"
        write_pdb(tmp, poly)  # TER at chain ends and real breaks
        fixer = PDBFixer(filename=str(tmp))
        fixer.missingResidues = {}  # never build residues here; step B owns that
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
    model = app.Modeller(fixer.topology, fixer.positions)
    model.addHydrogens(ff, pH=a.ph)
    system = ff.createSystem(
        model.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=CUTOFF_NM * unit.nanometer,
        constraints=None,
    )

    for force in system.getForces():
        force.setForceGroup(1)  # physical terms; restraints added below go to group 2
    pos0 = np.array(model.positions.value_in_unit(unit.nanometer))
    topo_atoms = list(model.topology.atoms())
    heavy = [at.index for at in topo_atoms if at.element is not None and at.element.symbol != "H"]
    pocket = set()
    het_xyz = np.array([[h.x, h.y, h.z] for h in hets if h.element not in ("H", "D")]) / 10.0
    # The pocket-collapse gate ignores ions: metal sites are legitimately coordinated at about 2 Å.
    organic = [h for h in hets if h.element not in ("H", "D") and h.resname not in ION_RES]
    org_xyz = np.array([[h.x, h.y, h.z] for h in organic]) if organic else np.zeros((0, 3))
    if len(het_xyz):
        near = cKDTree(het_xyz).query(pos0[heavy])[0] < a.pocket_cutoff / 10.0
        pocket_res = {topo_atoms[i].residue for i, close in zip(heavy, near) if close}
        pocket = {i for i in heavy if topo_atoms[i].residue in pocket_res}

    for gname, members in (("k_all", heavy), ("k_pocket", sorted(pocket))):
        f = openmm.CustomExternalForce(f"0.5*{gname}*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        f.setForceGroup(2)
        f.addGlobalParameter(gname, 0.0)
        for p in ("x0", "y0", "z0"):
            f.addPerParticleParameter(p)
        for i in members:
            f.addParticle(i, pos0[i].tolist())
        system.addForce(f)

    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
    ctx.setPositions(model.positions)
    n = model.topology.getNumAtoms()

    def restraint_free():
        state = ctx.getState(getEnergy=True, groups={1})
        return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    t0 = time.time()
    e0 = restraint_free()
    ctx.setParameter("k_all", a.restraint_k)
    openmm.LocalEnergyMinimizer.minimize(ctx, 10.0, a.max_iter)
    ctx.setParameter("k_all", 0.0)
    ctx.setParameter("k_pocket", a.restraint_k)
    openmm.LocalEnergyMinimizer.minimize(ctx, 10.0, a.free_iter)
    e2 = restraint_free()
    pos2 = np.array(ctx.getState(getPositions=True).getPositions().value_in_unit(unit.nanometer))

    res = {
        "status": "ok",
        "n_atoms": n,
        "n_heavy": len(heavy),
        "pocket_atoms_restrained": len(pocket),
        "solvent": solvent,
        "platform": ctx.getPlatform().getName(),
        "e_initial_kj_per_atom": round(e0 / n, 2),
        "e_final_kj_per_atom": round(e2 / n, 2),
        "seconds": round(time.time() - t0, 1),
    }
    if not (np.isfinite(pos2).all() and math.isfinite(e0) and math.isfinite(e2)):
        res.update(status="failed", failure="explosion", reason="non-finite energy or coordinates")
        out_path.unlink(missing_ok=True)  # never leave an older, passing output behind
        return res

    drop = (e2 - e0) / n
    res["rmsd_heavy"] = round(_kabsch_rmsd(pos0[heavy], pos2[heavy]) * 10, 3)
    res["max_displacement"] = round(float(np.linalg.norm(pos2[heavy] - pos0[heavy], axis=1).max()) * 10, 2)
    res["drop_kj_per_atom"] = round(drop, 2)

    # Back to heavy-atom Atoms (author IDs), B-factors carried over; HETs unchanged.
    bfac = {(x.chain, x.resseq, x.icode, x.name): x.b for x in poly}
    relaxed = []
    for at, xyz in zip(model.topology.atoms(), pos2):
        if at.element is None or at.element.symbol == "H":
            continue
        r = at.residue
        key = (r.chain.id, int(r.id), (r.insertionCode or "").strip(), at.name)
        relaxed.append(
            Atom(
                "ATOM",
                at.name,
                r.name,
                r.chain.id,
                int(r.id),
                key[2],
                *(xyz * 10).tolist(),
                at.element.symbol.upper(),
                bfac.get(key, 0.0),
                1.0 if key in bfac else 0.0,
            )
        )
    new_breaks = _ca_breaks(relaxed) - _ca_breaks(poly)
    if len(org_xyz):
        start_xyz = np.array([[x.x, x.y, x.z] for x in poly])
        rel_xyz = np.array([[x.x, x.y, x.z] for x in relaxed])
        res["pocket_min_dist_initial"] = round(float(cKDTree(start_xyz).query(org_xyz)[0].min()), 2)
        res["pocket_min_dist"] = round(float(cKDTree(rel_xyz).query(org_xyz)[0].min()), 2)

    failures = []
    if res["rmsd_heavy"] > a.max_rmsd:
        failures.append("implosion/explosion" if drop < -a.max_drop else "distortion")
    if new_breaks:
        failures.append("chain break")
        res["ca_breaks"] = sorted(new_breaks)[:10]
    if e2 > 0:
        failures.append("unresolved clash")
    # Collapse = the protein moved *into* the held-out ligand's space: closer than 2.0 Å
    # and at least 0.5 Å closer than in the input. Contacts that are already short in the
    # crystal (covalent cofactors, Schiff bases) are not failures.
    if (
        "pocket_min_dist" in res
        and res["pocket_min_dist"] < 2.0
        and res["pocket_min_dist"] < res["pocket_min_dist_initial"] - 0.5
    ):
        failures.append("pocket collapse")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    failed_path = out_path.with_suffix(".failed.pdb")
    if failures:
        res.update(status="failed", failure=", ".join(failures))
        # Failed structures are kept for inspection under a name no downstream glob matches.
        out_path.unlink(missing_ok=True)
        write_pdb(failed_path, relaxed + [h for h in hets if h.element not in ("H", "D")])
        res["output"] = failed_path.name
    else:
        failed_path.unlink(missing_ok=True)
        write_pdb(out_path, relaxed + [h for h in hets if h.element not in ("H", "D")])
    return res


def _family(a) -> str:
    """Forcefield family from the batch lock. An explicit --forcefield must agree
    with an existing lock (ForcefieldMismatch otherwise); amber is used only when
    no lock exists and nothing was requested."""
    if lock_path(a.batch_dir).exists():
        return acquire_lock(a.batch_dir, a.forcefield).family
    return resolve(a.forcefield).family if a.forcefield else "amber"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gromatrix-relax", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--in", dest="indir", default=layout.UNITS)
    ap.add_argument("--out", default=layout.RELAXED)
    ap.add_argument("--ids", nargs="+")
    ap.add_argument("--files", nargs="+", default=["enzyme.pdb", "complex_*.pdb", "target_*.pdb"])
    ap.add_argument("--forcefield", help="amber|charmm; default: the batch lock's family (amber if unlocked)")
    ap.add_argument(
        "--solvent",
        default="auto",
        choices=["auto", "gbn2", "rf"],
        help="auto = gbn2 on CUDA/OpenCL, rf on CPU; CHARMM always uses rf",
    )
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--restraint-k", type=float, default=1000.0, help="kJ/mol/nm^2 on heavy atoms")
    ap.add_argument("--max-iter", type=int, default=2000, help="stage 1 (restrained) iterations")
    ap.add_argument("--free-iter", type=int, default=1000, help="stage 2 (pocket-restrained) iterations")
    ap.add_argument("--max-rmsd", type=float, default=2.5, help="Å, heavy-atom RMSD gate")
    ap.add_argument(
        "--max-drop",
        type=float,
        default=100.0,
        help="kJ/mol/atom; a drop below -this together with RMSD > --max-rmsd = implosion/explosion",
    )
    ap.add_argument("--pocket-cutoff", type=float, default=5.0, help="Å around HETs kept restrained in stage 2")
    ap.add_argument("--platform", default="auto")
    a = ap.parse_args(argv)

    try:
        family = _family(a)
    except (ForcefieldMismatch, KeyError) as e:
        print(json.dumps({"ok": False, "error": f"{type(e).__name__}: {e}"}))
        return 1
    units = (
        sorted(p for p in (a.batch_dir / a.indir).iterdir() if p.is_dir()) if (a.batch_dir / a.indir).exists() else []
    )
    if a.ids:
        units = [u for u in units if u.name.upper() in {i.upper() for i in a.ids}]
    if not units:
        print(json.dumps({"ok": False, "error": f"no unit directories in {a.batch_dir / a.indir}"}))
        return 1
    results, failed = [], []
    for u in units:
        files = sorted(
            {
                f
                for pat in a.files
                for f in u.glob("*.pdb")
                if fnmatch.fnmatch(f.name, pat) and not f.name.endswith(".failed.pdb")
            }
        )
        report = {}
        for f in files:
            try:
                r = relax_file(f, a.batch_dir / a.out / u.name / f.name, family, a)
            except Exception as e:
                r = {"status": "failed", "failure": "error", "reason": f"{type(e).__name__}: {str(e)[:300]}"}
                (a.batch_dir / a.out / u.name).mkdir(parents=True, exist_ok=True)
                (a.batch_dir / a.out / u.name / f"{f.stem}.error.txt").write_text(traceback.format_exc())
            report[f.name] = r
            brief = {
                k: r.get(k)
                for k in ("status", "failure", "rmsd_heavy", "drop_kj_per_atom", "pocket_min_dist")
                if r.get(k) is not None
            }
            (failed if r["status"] == "failed" else results).append({"id": u.name, "file": f.name, **brief})
        (a.batch_dir / a.out / u.name).mkdir(parents=True, exist_ok=True)
        (a.batch_dir / a.out / u.name / "relax.json").write_text(
            json.dumps(
                {
                    "forcefield_family": family,
                    "thresholds": {"max_rmsd": a.max_rmsd, "max_drop": a.max_drop},
                    "files": report,
                },
                indent=1,
            )
        )
    print(json.dumps({"ok": not failed, "family": family, "relaxed": results, "failed": failed}))
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
