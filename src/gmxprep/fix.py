"""[B] Crystallographic clean-up with PDBFixer.

    gmxprep-fix --batch-dir B [--in raw] [--out fixed] [--ids 1ABC ...]
                [--max-loop 8] [--waters none|ligand-shell|shell|all] [--water-shell 3.5]
                [--keep-het all|none|RES,RES] [--keep-additives] [--keep-nonstandard]

For every raw/<ID>.cif|.pdb this writes fixed/<ID>.pdb plus fixed/<ID>.meta.json,
which carries assemblies, entities and the gap report on to step C.

What it does, in order:
  1. maps OpenMM's chain IDs back to *author* chain IDs. OpenMM labels mmCIF
     chains by label_asym_id, which would move every ligand, ion and water onto
     its own label chain (a ligand bound to author chain D can appear as chain K).
  2. replaces non-standard residues (MSE->MET, ...) unless --keep-nonstandard.
  3. finds missing residues. Only internal gaps of <= --max-loop residues are
     modelled. Termini and longer loops are left open, reported, and split by
     TER records so no tool bonds across them.
  4. adds missing heavy atoms (incl. OXT). Never adds hydrogens.
  5. strips every H/D atom, hydronium, unknown HET (UNX/UNL), common
     crystallisation additives, and waters according to --waters.
  6. alternate conformers: the primary conformer of every residue is the one with
     the highest occupancy (ties: first letter), and its coordinates overwrite
     whatever altloc OpenMM picked. All conformers are analysed (raw and
     symmetry-corrected RMSD, bond-geometry consistency) and stored in
     meta.json for step C. Coordinates are never averaged (see altlocs.py).
  7. writes author-numbered PDB with the original B-factors. Atoms built by
     PDBFixer get occupancy 0.00 so they can be told apart later.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np
from openmm.unit import angstrom

from . import altlocs
from .pdbio import ION_RES, WATER_RES, read_pdb, residues, write_pdb
from .structure import POLYMER_RES, read_metadata, read_raw_atoms

# Crystallisation / cryo additives that are almost never biologically relevant.
ADDITIVES = {
    "GOL",
    "EDO",
    "PEG",
    "PGE",
    "PG4",
    "1PE",
    "P6G",
    "MPD",
    "DMS",
    "SO4",
    "PO4",
    "ACT",
    "FMT",
    "TRS",
    "EPE",
    "MES",
    "BME",
    "IPA",
    "EOH",
    "MOH",
    "NO3",
    "SCN",
    "IOD",
    "BR",
    "CIT",
    "TAR",
    "MLI",
    "IMD",
    "BTB",
    "ACY",
    "PEO",
    "DTT",
    "AZI",
    "NH4",
    "CAC",
}
HYDRONIUM = {"H3O", "OH3"}
UNKNOWN_HET = {"UNX", "UNL"}


def _is_h(atom) -> bool:
    el = atom.element.symbol if atom.element is not None else atom.name[:1]
    return el in ("H", "D")


def openmm_uses_label_ids(meta: dict) -> bool:
    """Mirror OpenMM's PDBx reader (openmm/app/pdbxfile.py). It uses auth_asym_id
    for chain IDs unless there are strictly more distinct label_asym_ids than
    auth_asym_ids, as with ligands and waters that have their own label chains."""
    l2a = meta.get("label2auth") or {}
    return bool(l2a) and len(l2a) > len(set(l2a.values()))


def _author_chain_ids(fixer, meta: dict) -> bool:
    """Rename chains to author IDs when OpenMM read label IDs. Returns True if remapped.
    The decision mirrors OpenMM's own rule. Guessing from chain names can swap
    chains when label and author IDs overlap (e.g. label A -> auth B, label B -> auth A)."""
    if not openmm_uses_label_ids(meta):
        return False
    l2a = meta["label2auth"]
    chains = list(fixer.topology.chains())
    unknown = [c.id for c in chains if c.id not in l2a]
    if unknown:
        raise ValueError(f"OpenMM chain IDs {unknown} are not label_asym_ids of this file")
    for c in chains:
        c.id = l2a[c.id]
    return True


def _gaps(fixer, max_loop: int) -> list[dict]:
    """Decide which missing segments PDBFixer may build; drop the rest."""
    chains = list(fixer.topology.chains())
    out = []
    for (ci, ri), names in sorted(fixer.missingResidues.items()):
        res = list(chains[ci].residues())
        terminal = ri == 0 or ri == len(res)
        modelled = not terminal and len(names) <= max_loop
        if not modelled:
            del fixer.missingResidues[(ci, ri)]
        out.append(
            {
                "chain": chains[ci].id,
                "after": res[ri - 1].id if ri > 0 else None,
                "before": res[ri].id if ri < len(res) else None,
                "n": len(names),
                "residues": names,
                "terminal": terminal,
                "modelled": modelled,
            }
        )
    return out


def _select_deletions(modeller, args) -> tuple[list, dict]:
    keep_list = None if args.keep_het in ("all", "none") else {r.strip().upper() for r in args.keep_het.split(",")}
    removed = {"hydrogens": 0, "waters": 0, "additives": [], "hydronium": 0, "unknown": [], "het_filtered": []}
    pos = np.asarray(modeller.positions.value_in_unit(angstrom))
    delete, kept_ligand_idx, polymer_idx, waters = [], [], [], []

    for res in modeller.topology.residues():
        name = res.name.upper()
        heavy = [a for a in res.atoms() if not _is_h(a)]
        if name in WATER_RES:
            waters.append(res)
            continue
        if name in HYDRONIUM:
            delete.append(res)
            removed["hydronium"] += 1
            continue
        if name in POLYMER_RES:
            polymer_idx += [a.index for a in heavy]
        elif name in UNKNOWN_HET:
            delete.append(res)
            removed["unknown"].append(f"{name}{res.id}{res.chain.id}")
            continue
        elif args.keep_het == "none" or (keep_list is not None and name not in keep_list):
            delete.append(res)
            removed["het_filtered"].append(f"{name}{res.id}{res.chain.id}")
            continue
        elif name in ADDITIVES and not args.keep_additives and not (keep_list and name in keep_list):
            delete.append(res)
            removed["additives"].append(f"{name}{res.id}{res.chain.id}")
            continue
        elif name not in ION_RES:
            kept_ligand_idx += [a.index for a in heavy]
        for a in res.atoms():
            if _is_h(a):
                delete.append(a)
                removed["hydrogens"] += 1

    # Waters: keep none / the ligand-site shell / the first shell of the solute / all.
    if waters and args.waters != "all":
        anchor = {"ligand-shell": kept_ligand_idx, "shell": kept_ligand_idx + polymer_idx}.get(args.waters, [])
        tree = None
        if anchor:
            from scipy.spatial import cKDTree

            tree = cKDTree(pos[anchor])
        for w in waters:
            o = [a for a in w.atoms() if not _is_h(a)]
            if tree is None or not o or tree.query(pos[o[0].index])[0] > args.water_shell:
                delete.append(w)
                removed["waters"] += 1
            else:
                delete += [a for a in w.atoms() if _is_h(a)]
    else:
        for w in waters:
            delete += [a for a in w.atoms() if _is_h(a)]
    return delete, removed


def fix_one(path: Path, out_dir: Path, args) -> dict:
    from openmm.app import Modeller, PDBFile
    from pdbfixer import PDBFixer

    pid = path.stem.upper()
    warnings: list[str] = []
    meta = read_metadata(path)
    raw_atoms, raw_stats, alt_groups = read_raw_atoms(path, args.altloc_rmsd)

    # PDBFixer sees only the primary conformer, so missing atoms and loops are built
    # on the conformer we keep (OpenMM would otherwise keep the first-listed altloc).
    tmpdir = tempfile.TemporaryDirectory()
    src = altlocs.reduce_to_primary(path, alt_groups, Path(tmpdir.name)) if alt_groups else path
    fixer = PDBFixer(filename=str(src))
    fixer.findMissingResidues()  # before the remap: PDBFixer matches SEQRES by its own chain IDs
    remapped = _author_chain_ids(fixer, meta)
    gaps = _gaps(fixer, args.max_loop)
    fixer.findNonstandardResidues()
    nonstd = [f"{r.name}{r.id}{r.chain.id}->{new}" for r, new in fixer.nonstandardResidues]
    if nonstd and not args.keep_nonstandard:
        fixer.replaceNonstandardResidues()
    elif nonstd:
        warnings.append(f"kept non-standard residues {nonstd}; pdb2gmx will reject them")
    fixer.findMissingAtoms()
    n_missing_atoms = sum(len(v) for v in fixer.missingAtoms.values())
    n_terminals = sum(len(v) for v in fixer.missingTerminals.values())
    try:
        fixer.addMissingAtoms(seed=0)
    except TypeError:  # older PDBFixer without the seed argument
        fixer.addMissingAtoms()

    tmpdir.cleanup()
    modeller = Modeller(fixer.topology, fixer.positions)
    delete, removed = _select_deletions(modeller, args)
    modeller.delete(delete)

    if any(len(c.id) > 1 for c in modeller.topology.chains()):
        raise ValueError(
            "multi-character chain IDs cannot be written as PDB; extract a smaller unit from the mmCIF first"
        )
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp:
        PDBFile.writeFile(modeller.topology, modeller.positions, tmp, keepIds=True)
    atoms = read_pdb(Path(tmp.name))
    Path(tmp.name).unlink()

    # Re-attach crystallographic B-factors; PDBFixer-built atoms get occupancy 0.
    # The key includes the residue name, and positions whose numbering became
    # ambiguous (PDBFixer numbers built loop residues backwards from the next
    # residue, which can collide with inconsistent author numbering) are not restored.
    ref = {(a.chain, a.resseq, a.icode, a.resname, a.name): a for a in raw_atoms}
    seen_res, dup = set(), set()
    for res in residues(atoms):
        k = (res[0].chain, res[0].resseq, res[0].icode)
        (dup if k in seen_res else seen_res).add(k)
    if dup:
        warnings.append(
            f"residue numbering collides after loop building at {sorted(dup)[:5]}; "
            "coordinates there are PDBFixer's, check the loop"
        )
    restored, built = [], 0
    for a in atoms:
        r = None if (a.chain, a.resseq, a.icode) in dup else ref.get((a.chain, a.resseq, a.icode, a.resname, a.name))
        if r is None:
            built += 1
            restored.append(a.__class__(**{**a.__dict__, "b": 0.0, "occ": 0.0}))
        else:
            restored.append(a.__class__(**{**a.__dict__, "x": r.x, "y": r.y, "z": r.z, "b": r.b, "occ": r.occ}))
    out_pdb = out_dir / f"{pid}.pdb"
    write_pdb(out_pdb, restored)

    open_gaps = [g for g in gaps if not g["modelled"] and not g["terminal"]]
    if open_gaps:
        warnings.append(
            "internal gaps left open (chain split with TER, charged termini in MD): "
            + ", ".join(f"{g['chain']}:{g['after']}-{g['before']} ({g['n']} res)" for g in open_gaps)
        )
    for line in altlocs.summarise(alt_groups, POLYMER_RES):
        warnings.append("alternate conformers: " + line)
    if removed["unknown"]:
        warnings.append(f"removed unknown HET groups {removed['unknown']}")
    if args.waters != "none":
        warnings.append("crystal waters kept; step F currently drops them when building topologies")
    if not remapped and meta["format"] == "mmcif":
        warnings.append("could not map label->author chain IDs; check chain IDs in the output")

    report = {
        "input": str(path),
        "author_chain_remap": remapped,
        "raw": raw_stats,
        "gaps": gaps,
        "nonstandard": nonstd,
        "missing_heavy_atoms": n_missing_atoms,
        "missing_terminal_atoms": n_terminals,
        "built_atoms": built,
        "removed": removed,
        "waters_policy": args.waters,
        "warnings": warnings,
    }
    meta["fix"] = report
    meta["altlocs"] = alt_groups
    (out_dir / f"{pid}.meta.json").write_text(json.dumps(meta, indent=1))
    return {
        "id": pid,
        "out": str(out_pdb),
        "built_atoms": built,
        "modelled_gaps": sum(g["modelled"] for g in gaps),
        "open_gaps": len(open_gaps),
        "removed_h": removed["hydrogens"],
        "removed_waters": removed["waters"],
        "removed_additives": len(removed["additives"]),
        "nonstandard": len(nonstd),
        "altloc_residues": len(alt_groups),
        "distinct_het_conformers": [
            f"{g['resname']}{g['resseq']}{g['chain']}"
            for g in alt_groups
            if g["distinct"] and g["resname"] not in POLYMER_RES
        ],
        "warnings": warnings,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gmxprep-fix", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--in", dest="indir", default="raw")
    ap.add_argument("--out", default="fixed")
    ap.add_argument("--ids", nargs="+", help="restrict to these IDs (file stems)")
    ap.add_argument("--max-loop", type=int, default=8, help="longest internal gap PDBFixer may build (residues)")
    ap.add_argument("--waters", default="none", choices=["none", "ligand-shell", "shell", "all"])
    ap.add_argument("--water-shell", type=float, default=3.5, help="Å, for ligand-shell/shell")
    ap.add_argument("--keep-het", default="all", help="all | none | comma list of residue names")
    ap.add_argument("--keep-additives", action="store_true", help="do not strip GOL/SO4/EDO/...")
    ap.add_argument("--keep-nonstandard", action="store_true")
    ap.add_argument(
        "--altloc-rmsd",
        type=float,
        default=0.5,
        help="Å; symmetry-corrected RMSD at which alternate conformers count as distinct poses",
    )
    a = ap.parse_args(argv)

    indir, outdir = a.batch_dir / a.indir, a.batch_dir / a.out
    outdir.mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in indir.glob("*") if p.suffix.lower() in (".cif", ".mmcif", ".pdb", ".ent"))
    if a.ids:
        want = {i.upper() for i in a.ids}
        files = [p for p in files if p.stem.upper() in want]
    if not files:
        print(json.dumps({"ok": False, "error": f"no structures in {indir}"}))
        return 1
    results, failed = [], []
    for p in files:
        try:
            results.append(fix_one(p, outdir, a))
        except Exception as e:
            (outdir / f"{p.stem.upper()}.error.txt").write_text(traceback.format_exc())
            failed.append({"id": p.stem.upper(), "error": f"{type(e).__name__}: {str(e)[:300]}"})
    print(json.dumps({"ok": not failed, "fixed": results, "failed": failed}))
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
