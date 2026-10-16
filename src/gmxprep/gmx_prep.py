"""[F] GROMACS preparation for an enzyme x target simulation matrix.

    gmx_prep lock   --batch-dir B --forcefield amber
    gmx_prep plan   --batch-dir B --manifest M
    gmx_prep build  --batch-dir B --manifest M [--only ID ...] [--box dodecahedron --pad 1.2 ...]
    gmx_prep verify --batch-dir B
    gmx_prep smoke  --batch-dir B [--only ID ...] [--nvt-steps 5000]

Each command prints one compact JSON object on stdout; tool chatter goes to
<batch>/systems/<id>/build.log and <batch>/params/<RES>/param.log.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import traceback
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from . import checks
from . import mdp as mdpmod
from . import topology as topo
from .forcefields import ForceField, ForcefieldMismatch, acquire_lock, assert_same, load_lock, resolve
from .gmx import ToolError, data_prefix, gmxlib_env, ion_atom_names, locate_ff_dir, run
from .ligands import LigandSpec, parameterise, pose_coords
from .pdbio import ION_RES, GroAtom, classify, read_gro, read_pdb, residue_block, write_gro, write_pdb


# ----------------------------------------------------------------- manifest
@dataclass
class SystemDef:
    id: str
    role: str  # enzyme | target | complex
    pdb: Path | None
    sdf: Path | None = None  # explicit-H pose for the target ligand
    target_resname: str | None = None
    protein_chains: set[str] | None = None
    enzyme: str | None = None
    target: str | None = None


@dataclass
class Manifest:
    root: Path
    forcefield: str | None
    systems: list[SystemDef]
    ligands: dict[str, LigandSpec]
    missing_pairs: list[str] = field(default_factory=list)


def load_manifest(path: Path, allow_incomplete: bool) -> Manifest:
    path = Path(path).resolve()
    root = path.parent
    raw = yaml.safe_load(path.read_text())

    def p(v):
        return (root / v).resolve() if v else None

    enzymes = {e["id"]: e for e in raw.get("enzymes", [])}
    targets = {t["id"]: t for t in raw.get("targets", [])}
    if not enzymes or not targets:
        raise ValueError("manifest needs at least one enzyme and one target")
    if set(enzymes) & set(targets):
        raise ValueError(f"ids used for both enzyme and target: {set(enzymes) & set(targets)}")

    # A pair may have several variants (e.g. altloc branches from step C: variant: altA / altB).
    complexes = {(c["enzyme"], c["target"], c.get("variant")): c for c in raw.get("complexes", [])}
    for e, t, _v in complexes:
        if e not in enzymes or t not in targets:
            raise ValueError(f"complex {e}+{t} references an unknown id")
    pairs = {(e, t) for e, t, _v in complexes}
    missing = [f"{e}__{t}" for e in enzymes for t in targets if (e, t) not in pairs]
    if missing and not allow_incomplete:
        raise ValueError(
            f"matrix incomplete - no docked complex for {missing}; "
            "every enzyme must be paired with every target (or pass --allow-incomplete)"
        )

    systems = [SystemDef(eid, "enzyme", p(e["pdb"])) for eid, e in enzymes.items()]
    for tid, t in targets.items():
        kind = t.get("kind", "ligand")
        if kind == "protein":
            systems.append(
                SystemDef(tid, "target", p(t["pdb"]), protein_chains=set(t["chains"]) if t.get("chains") else None)
            )
            continue
        src_pdb, src_sdf = p(t.get("pdb")), p(t.get("sdf"))
        if not src_pdb and not src_sdf:  # fall back to the first docked pose
            src_pdb = next((p(c["pdb"]) for (e, tt, _v), c in complexes.items() if tt == tid), None)
        systems.append(SystemDef(tid, "target", src_pdb, src_sdf, target_resname=t["resname"]))
    for (e, t, v), c in complexes.items():
        tdef = targets[t]
        tres = tdef.get("resname") if tdef.get("kind", "ligand") == "ligand" else None
        systems.append(
            SystemDef(
                f"{e}__{t}" + (f"__{v}" if v else ""),
                "complex",
                p(c["pdb"]),
                p(c.get("sdf")),
                target_resname=tres,
                enzyme=e,
                target=t,
            )
        )

    ligands = {res: LigandSpec.from_manifest(res, d, root) for res, d in (raw.get("ligands") or {}).items()}
    for t in targets.values():  # an explicit-H target SDF can serve as its own template
        if t.get("kind", "ligand") == "ligand" and t.get("sdf"):
            spec = ligands.setdefault(t["resname"], LigandSpec(t["resname"]))
            if not spec.sdf and not spec.smiles:
                spec.sdf = str(p(t["sdf"]))
    return Manifest(root, raw.get("forcefield"), systems, ligands, missing)


# ----------------------------------------------------------------- build
@dataclass
class Ctx:
    batch: Path
    ff: ForceField
    ff_dir: Path
    gmx: str
    args: argparse.Namespace
    ligand_specs: dict[str, LigandSpec]

    def g(self, *sub: str) -> list[str]:
        return [self.gmx, "-quiet", "-nobackup", *sub]

    @property
    def env(self) -> dict:
        return gmxlib_env(self.ff_dir)


def _components(sd: SystemDef):
    """-> (protein atoms, [(resname, pdb_block|None, sdf|None)], ion atoms)"""
    if sd.pdb is None:
        return [], [(sd.target_resname, None, str(sd.sdf))], []
    atoms = read_pdb(sd.pdb)
    if sd.role == "target" and sd.target_resname:  # ligand-alone system: keep only that residue
        atoms = [a for a in atoms if a.resname == sd.target_resname]
        if not atoms:
            raise ValueError(f"{sd.pdb} has no residue {sd.target_resname}")
    split = classify(atoms, sd.protein_chains)
    ligs = []
    for res in split.ligands:
        rn = res[0].resname
        sdf = str(sd.sdf) if sd.sdf and rn == sd.target_resname else None
        ligs.append((rn, residue_block(res), sdf))
    return split.protein, ligs, split.ions


def inputs_hash(sd: SystemDef, ctx: Ctx) -> str:
    """Fingerprint of everything a system is built from. A rebuild is skipped only
    when it is unchanged; this is SQUIP/TaskBlaster-style caching, but keyed on content."""
    a, h = ctx.args, hashlib.sha256()
    for f in (sd.pdb, sd.sdf):
        if f:
            h.update(Path(f).read_bytes())
    for res, spec in sorted(ctx.ligand_specs.items()):
        h.update(json.dumps([res, spec.__dict__], sort_keys=True, default=str).encode())
    h.update(
        json.dumps(
            [ctx.ff.identity(), a.box, a.pad, a.conc, a.temp, a.md_ns, sorted(sd.protein_chains or [])], sort_keys=True
        ).encode()
    )
    return h.hexdigest()[:16]


RUN_SH = """#!/usr/bin/env bash
# Generated by gmxprep. Coordinates and checkpoints are chained correctly:
# -r supplies POSRES reference coordinates, -t carries velocities and
# box/thermostat state forward. -maxwarn is never used.
set -euo pipefail
cd "$(dirname "$0")"
GMX=${GMX:-%(gmx)s}
$GMX grompp -f mdp/em.mdp  -c solv_ions.gro -r solv_ions.gro -p topol.top -n index.ndx -o em.tpr
$GMX mdrun -deffnm em
$GMX grompp -f mdp/nvt.mdp -c em.gro  -r em.gro  -p topol.top -n index.ndx -o nvt.tpr
$GMX mdrun -deffnm nvt
$GMX grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
$GMX mdrun -deffnm npt
$GMX grompp -f mdp/md.mdp  -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
$GMX mdrun -deffnm md
"""


def _link_ff(wd: Path, ctx: Ctx) -> None:
    """Make the system directory self-contained: a forcefield port that is not in
    the GROMACS share dir is symlinked next to topol.top, so a plain `gmx grompp`
    run later (without GMXLIB) still resolves the includes. SQUIP solves this by
    rewriting includes to absolute paths; a relative link survives moving the batch."""
    prefix = data_prefix(ctx.gmx)
    if prefix and ctx.ff_dir.parent == (prefix / "share" / "gromacs" / "top").resolve():
        return
    link = wd / ctx.ff_dir.name
    if not link.exists():
        link.symlink_to(Path(os.path.relpath(ctx.ff_dir, wd)))


def build_system(sd: SystemDef, ctx: Ctx) -> dict:
    ff, a = ctx.ff, ctx.args
    wd = ctx.batch / "systems" / sd.id
    stamp_path = wd / "system.json"
    digest = inputs_hash(sd, ctx)
    if stamp_path.exists() and not a.force:
        stamp = json.loads(stamp_path.read_text())
        if stamp.get("status") == "ok" and stamp.get("inputs_sha") == digest:
            assert_same(ff, stamp["forcefield"], f"system {sd.id}")
            return stamp
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)
    _link_ff(wd, ctx)
    log = wd / "build.log"
    warnings: list[str] = []

    protein, ligs, ions = _components(sd)
    gro: list[GroAtom] = []
    itps: list[str] = []
    molecules: list[tuple[str, int]] = []
    param_blocks: list[dict] = []

    # 1. protein (all chains merged into one moleculetype so inter-chain S-S bonds survive)
    if protein:
        write_pdb(wd / "protein.pdb", protein)
        run(
            ctx.g(
                "pdb2gmx",
                "-f",
                "protein.pdb",
                "-o",
                "protein.gro",
                "-p",
                "protein.top",
                # -water none: water comes from our own topol.top. This also avoids the
                # missing `-water tip4pew` enum of older GROMACS that SQUIP works around.
                "-i",
                "posre.itp",
                "-ff",
                ff.gmx_name,
                "-water",
                "none",
                "-ignh",
                "-merge",
                "all",
            ),
            wd,
            log,
            env_extra=ctx.env,
        )
        mol = topo.extract_pdb2gmx_moleculetype((wd / "protein.top").read_text())
        (wd / "protein.itp").write_text(mol)
        itps.append("protein.itp")
        molecules.append((topo.moleculetype_name(mol), 1))
        gro, _ = read_gro(wd / "protein.gro")

    # 2. cofactors / ligands (shared, cached parameters; per-pose coordinates)
    seen: dict[str, object] = {}
    for rn, block, sdf in ligs:
        spec = ctx.ligand_specs.get(rn) or LigandSpec(rn)
        params = parameterise(
            spec, block, ff, ctx.ff_dir, ctx.batch / "params", a.allow_auto_protonation, a.accept_high_penalty
        )
        if rn not in seen:
            seen[rn] = params
            shutil.copy(params.itp, wd / params.itp.name)
            itps.append(params.itp.name)
            param_blocks.append(params.params)
            warnings += params.warnings
        xyz = pose_coords(params, block, sdf)
        resnr = (gro[-1].resnr + 1) if gro else 1
        gro += [GroAtom(resnr, rn, n, *c) for n, c in zip(params.atom_names, xyz)]
        molecules.append((rn, 1))

    # 3. crystallographic ions, mapped to this forcefield's ion moleculetypes
    if ions:
        names = ion_atom_names(ctx.ff_dir)
        for ion in ions:
            elem = ION_RES[ion.resname]
            moltype = ff.ion_moltypes.get(elem)
            if moltype is None:
                raise ValueError(
                    f"{ff.key} has no non-bonded model for ion {ion.resname}; "
                    "treat it explicitly (e.g. MCPB.py) before building"
                )
            resnr = (gro[-1].resnr + 1) if gro else 1
            gro.append(GroAtom(resnr, moltype, names.get(moltype, moltype), ion.x / 10, ion.y / 10, ion.z / 10))
            molecules.append((moltype, 1))

    # 4. topology
    param_itp = None
    merged = topo.merge_params(param_blocks)
    if merged:
        (wd / "ligand_params.itp").write_text(merged)
        param_itp = "ligand_params.itp"
    topo.write_topol(
        wd / "topol.top", ff.gmx_name, ff.water.itp, param_itp, itps, molecules, f"{sd.id} ({sd.role}) {ff.key}"
    )
    write_gro(wd / "conf.gro", gro, f"{sd.id} {ff.key}")

    # 5. box, solvent, ions; then a zero-warning grompp as the integrity gate
    mdpmod.write_all(ff, wd / "mdp", a.temp, a.md_ns)
    run(
        ctx.g("editconf", "-f", "conf.gro", "-o", "boxed.gro", "-bt", a.box, "-d", str(a.pad), "-c"),
        wd,
        log,
        env_extra=ctx.env,
    )
    run(
        ctx.g("solvate", "-cp", "boxed.gro", "-cs", ff.water.solvent_gro, "-o", "solv.gro", "-p", "topol.top"),
        wd,
        log,
        env_extra=ctx.env,
    )

    # Pre-flight (before any grompp): topology vs coordinates, and an integer net charge.
    search = [ctx.ff_dir.parent, wd]
    pre = checks.top_vs_gro(wd / "topol.top", wd / "solv.gro", search)
    if pre["problems"]:
        raise checks.CheckError("topology/coordinate mismatch after solvate: " + "; ".join(pre["problems"]))
    q0 = checks.assert_integer_charge(pre["charge"], f"{sd.id} solute")
    # Only a charged system may pass the ions-stage grompp with one (net-charge) warning.
    run(
        ctx.g(
            "grompp",
            "-f",
            "mdp/ions.mdp",
            "-c",
            "solv.gro",
            "-p",
            "topol.top",
            "-o",
            "ions.tpr",
            "-maxwarn",
            "1" if q0 else "0",
        ),
        wd,
        log,
        env_extra=ctx.env,
    )
    # Select the solvent group by name, never by number (SQUIP pipes `echo 15`, which depends on group order).
    run(
        ctx.g(
            "genion",
            "-s",
            "ions.tpr",
            "-o",
            "solv_ions.gro",
            "-p",
            "topol.top",
            "-pname",
            ff.cation,
            "-nname",
            ff.anion,
            "-neutral",
            "-conc",
            str(a.conc),
        ),
        wd,
        log,
        stdin="SOL\n",
        env_extra=ctx.env,
    )
    post = checks.top_vs_gro(wd / "topol.top", wd / "solv_ions.gro", search)
    if post["problems"] or abs(post["charge"]) > 1e-3:
        raise checks.CheckError(f"after genion: charge {post['charge']} / {post['problems']}")
    protein_mol = molecules[0][0] if protein else None
    groups = checks.write_index(
        wd / "solv_ions.gro", wd / "topol.top", search, wd / "index.ndx", protein_mol, sorted(seen)
    )
    run(
        ctx.g(
            "grompp",
            "-f",
            "mdp/em.mdp",
            "-c",
            "solv_ions.gro",
            "-r",
            "solv_ions.gro",
            "-p",
            "topol.top",
            "-n",
            "index.ndx",
            "-o",
            "em.tpr",
            "-maxwarn",
            str(a.maxwarn),
        ),
        wd,
        log,
        env_extra=ctx.env,
    )
    (wd / "run.sh").write_text(RUN_SH % {"gmx": a.gmx})
    (wd / "run.sh").chmod(0o755)

    lig_counts: dict[str, int] = {}
    for rn, *_ in ligs:
        lig_counts[rn] = lig_counts.get(rn, 0) + 1
    stamp = {
        "id": sd.id,
        "role": sd.role,
        "enzyme": sd.enzyme,
        "target": sd.target,
        "status": "ok",
        "forcefield": ff.identity(),
        "inputs_sha": digest,
        "atoms": post["natoms_gro"],
        "solute_charge": q0,
        "index_groups": groups,
        "components": {"protein_heavy_atoms": len(protein), "ligands": lig_counts, "ions": [i.resname for i in ions]},
        "warnings": warnings,
        "next": "gmxprep-gmx smoke (EM + short NVT check), then ./run.sh",
    }
    stamp_path.write_text(json.dumps(stamp, indent=2))
    return stamp


# ----------------------------------------------------------------- commands
def cmd_lock(a) -> dict:
    ff = acquire_lock(a.batch_dir, a.forcefield)
    return {"ok": True, "locked": ff.identity()}


def cmd_plan(a) -> dict:
    m = load_manifest(a.manifest, a.allow_incomplete)
    return {
        "ok": True,
        "systems": [{"id": s.id, "role": s.role} for s in m.systems],
        "missing_pairs": m.missing_pairs,
        "ligand_specs": sorted(m.ligands),
    }


def cmd_build(a) -> dict:
    m = load_manifest(a.manifest, a.allow_incomplete)
    if a.forcefield and m.forcefield and resolve(a.forcefield).key != resolve(m.forcefield).key:
        raise ForcefieldMismatch(f"--forcefield {a.forcefield} disagrees with manifest forcefield {m.forcefield}")
    ff = acquire_lock(a.batch_dir, a.forcefield or m.forcefield)
    ff_dir = locate_ff_dir(ff.gmx_name, a.gmx, [a.batch_dir.resolve()])
    ctx = Ctx(a.batch_dir.resolve(), ff, ff_dir, a.gmx, a, m.ligands)
    results = []
    for sd in m.systems:
        if a.only and sd.id not in a.only:
            continue
        try:
            st = build_system(sd, ctx)
            results.append({"id": sd.id, "role": sd.role, "status": "ok", "warnings": st["warnings"]})
        except ForcefieldMismatch:
            raise
        except Exception as e:  # keep going; report per system
            (ctx.batch / "systems" / sd.id).mkdir(parents=True, exist_ok=True)
            (ctx.batch / "systems" / sd.id / "error.txt").write_text(traceback.format_exc())
            results.append({"id": sd.id, "role": sd.role, "status": "failed", "error": str(e)[:400]})
    (ctx.batch / "matrix.json").write_text(
        json.dumps({"forcefield": ff.identity(), "missing_pairs": m.missing_pairs, "systems": results}, indent=2)
    )
    return {
        "ok": all(r["status"] == "ok" for r in results),
        "forcefield": ff.key,
        "built": sum(r["status"] == "ok" for r in results),
        "failed": [r for r in results if r["status"] != "ok"],
        "warnings": sorted({w for r in results for w in r.get("warnings", [])}),
    }


def cmd_verify(a) -> dict:
    ff = load_lock(a.batch_dir)
    problems = []
    try:
        search = [locate_ff_dir(ff.gmx_name, a.gmx, [a.batch_dir.resolve()]).parent]
    except ToolError:
        search = [a.batch_dir.resolve()]
    for meta in sorted((a.batch_dir / "params").glob("*/meta.json")):
        try:
            assert_same(ff, json.loads(meta.read_text())["forcefield"], str(meta.parent))
        except ForcefieldMismatch as e:
            problems.append(str(e))
    systems = sorted((a.batch_dir / "systems").glob("*/"))
    for sd in systems:
        stamp_p = sd / "system.json"
        if not stamp_p.exists():
            problems.append(f"{sd.name}: not built")
            continue
        stamp = json.loads(stamp_p.read_text())
        try:
            assert_same(ff, stamp["forcefield"], sd.name)
        except ForcefieldMismatch as e:
            problems.append(str(e))
        top = (sd / "topol.top").read_text()
        if f'#include "{ff.gmx_name}.ff/forcefield.itp"' not in top or top.count(".ff/forcefield.itp") != 1:
            problems.append(f"{sd.name}: topol.top does not include exactly {ff.gmx_name}.ff/forcefield.itp")
        for f in (sd / "mdp").glob("*.mdp"):
            got = mdpmod.read_mdp(f)
            bad = {k: got.get(k) for k, v in ff.nonbonded.items() if got.get(k) != str(v)}
            if bad:
                problems.append(f"{sd.name}/{f.name}: non-bonded settings differ from {ff.key}: {bad}")
        if not (sd / "em.tpr").exists():
            problems.append(f"{sd.name}: em.tpr missing (grompp gate not passed)")
        if (sd / "solv_ions.gro").exists():
            chk = checks.top_vs_gro(sd / "topol.top", sd / "solv_ions.gro", search)
            problems += [f"{sd.name}: {p}" for p in chk["problems"]]
            if abs(chk["charge"]) > 1e-3:
                problems.append(f"{sd.name}: net charge {chk['charge']} (expected 0 after genion)")
        smoke = stamp.get("smoke")
        if smoke and not smoke.get("ok"):
            problems.append(f"{sd.name}: smoke test failed: {smoke.get('reason')}")
    matrix = a.batch_dir / "matrix.json"
    if matrix.exists():
        mj = json.loads(matrix.read_text())
        if mj.get("missing_pairs"):
            problems.append(f"matrix incomplete: {mj['missing_pairs']}")
    return {"ok": not problems, "forcefield": ff.key, "systems": len(systems), "problems": problems}


def _smoke_one(wd: Path, a, env: dict) -> dict:
    """EM to convergence, then a short restrained NVT. This is SQUIP's test_mode
    idea, with pass/fail decided from the logs instead of by eye."""

    def g(*x):
        return [a.gmx, "-quiet", "-nobackup", *x]

    log = wd / "smoke.log"
    emtol = float(mdpmod.read_mdp(wd / "mdp" / "em.mdp").get("emtol", 1000))
    run(g("mdrun", "-deffnm", "em", *a.mdrun_args.split()), wd, log, env_extra=env)
    em = checks.parse_em_log(wd / "em.log", emtol)
    res = {"em": em, "ok": em["converged"]}
    if not em["converged"]:
        res["reason"] = (
            f"EM did not reach Fmax < {emtol} (Fmax={em['fmax']}, Epot={em['epot']})"
            " - look for clashes or a bad ligand pose"
        )
        return res
    if a.nvt_steps:
        params = mdpmod.read_mdp(wd / "mdp" / "nvt.mdp")
        params.update({"nsteps": str(a.nvt_steps), "nstxout-compressed": "0"})
        (wd / "mdp" / "smoke_nvt.mdp").write_text("".join(f"{k:<24s} = {v}\n" for k, v in params.items()))
        run(
            g(
                "grompp",
                "-f",
                "mdp/smoke_nvt.mdp",
                "-c",
                "em.gro",
                "-r",
                "em.gro",
                "-p",
                "topol.top",
                "-n",
                "index.ndx",
                "-o",
                "smoke_nvt.tpr",
            ),
            wd,
            log,
            env_extra=env,
        )
        run(g("mdrun", "-deffnm", "smoke_nvt", *a.mdrun_args.split()), wd, log, env_extra=env)
        md = checks.scan_md_log(wd / "smoke_nvt.log")
        res["nvt"] = md
        # Stricter than SQUIP's "< 10 LINCS warnings": a healthy 10 ps restrained NVT has none.
        res["ok"] = md["finished"] and not md["fatal"] and md["lincs_warnings"] == 0
        if not res["ok"]:
            res["reason"] = f"NVT smoke: {md}"
    return res


def cmd_smoke(a) -> dict:
    ff = load_lock(a.batch_dir)
    env = gmxlib_env(locate_ff_dir(ff.gmx_name, a.gmx, [a.batch_dir.resolve()]))
    out = []
    for stamp_p in sorted((a.batch_dir / "systems").glob("*/system.json")):
        wd = stamp_p.parent
        if a.only and wd.name not in a.only:
            continue
        stamp = json.loads(stamp_p.read_text())
        if stamp.get("status") != "ok":
            continue
        try:
            res = _smoke_one(wd, a, env)
        except ToolError as e:
            res = {"ok": False, "reason": str(e)[:300]}
        stamp["smoke"] = res
        stamp_p.write_text(json.dumps(stamp, indent=2))
        out.append(
            {
                "id": wd.name,
                "ok": res["ok"],
                **({"reason": res["reason"]} if not res["ok"] else {"fmax": res["em"]["fmax"]}),
            }
        )
    return {"ok": bool(out) and all(r["ok"] for r in out), "systems": out}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gmx_prep", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, manifest=False):
        p.add_argument("--batch-dir", type=Path, required=True)
        if manifest:
            p.add_argument("--manifest", type=Path, required=True)
            p.add_argument(
                "--allow-incomplete", action="store_true", help="permit enzyme/target pairs without a docked complex"
            )
        return p

    p = common(sub.add_parser("lock", help="create or check the batch forcefield lock"))
    p.add_argument("--forcefield", required=True)
    common(sub.add_parser("plan", help="list the systems the manifest expands to"), manifest=True)
    p = common(sub.add_parser("build", help="build topologies, boxes, solvent, ions, em.tpr"), manifest=True)
    p.add_argument("--forcefield", help="must match the lock (and manifest) if given")
    p.add_argument("--only", nargs="+", help="system ids to (re)build")
    p.add_argument("--force", action="store_true", help="rebuild systems already marked ok")
    p.add_argument("--box", default="dodecahedron", choices=["dodecahedron", "cubic", "octahedron", "triclinic"])
    p.add_argument("--pad", type=float, default=1.2, help="solute-box distance, nm")
    p.add_argument("--conc", type=float, default=0.15, help="salt, mol/L")
    p.add_argument("--temp", type=float, default=300.0)
    p.add_argument("--md-ns", type=float, default=100.0)
    p.add_argument("--maxwarn", type=int, default=0, help="for the final em grompp gate")
    p.add_argument("--gmx", default="gmx")
    p.add_argument("--allow-auto-protonation", action="store_true")
    p.add_argument("--accept-high-penalty", action="store_true")
    p = common(sub.add_parser("verify", help="prove the whole batch uses one forcefield and is self-consistent"))
    p.add_argument("--gmx", default="gmx")
    p = common(sub.add_parser("smoke", help="EM + short NVT per system; record pass/fail in system.json"))
    p.add_argument("--only", nargs="+")
    p.add_argument("--nvt-steps", type=int, default=5000, help="0 = EM only")
    p.add_argument("--mdrun-args", default="", help='extra mdrun flags, e.g. "-nt 8 -gpu_id 0"')
    p.add_argument("--gmx", default="gmx")

    a = ap.parse_args(argv)
    try:
        out = {"lock": cmd_lock, "plan": cmd_plan, "build": cmd_build, "verify": cmd_verify, "smoke": cmd_smoke}[a.cmd](
            a
        )
    except (ForcefieldMismatch, ValueError, KeyError, FileNotFoundError) as e:
        out = {"ok": False, "error": f"{type(e).__name__}: {e}"}
    print(json.dumps(out))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
