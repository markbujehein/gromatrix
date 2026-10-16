---
name: gromatrix-structure-prep
description: Steps A–C of the gromatrix pipeline. Fetch structures from the RCSB (by ID, text or UniProt search), clean them with PDBFixer (author chain IDs, missing atoms, short loops, hydrogens and waters stripped, alternate conformers analysed), and extract the right sub-unit from any asymmetric unit or biological assembly, sliced into enzyme / target / complex PDB files. Use this when the user wants PDB entries downloaded or cleaned, needs a biological unit or dimer chosen, has ligands bound in only some copies of a multi-chain structure, asks about alternate conformers (altlocs), or when `gmxprep-status` names this skill as `next`.
---

# gromatrix: structure preparation (steps A–C)

Commands: `gmxprep-fetch`, `gmxprep-fix`, `gmxprep-extract`. Full flags are in
`references/cli.md`; read it only when you need a flag that isn't shown here. Error and warning
texts are explained in `references/troubleshooting.md`.

## Before you start

```bash
gmxprep-status --batch-dir B     # what already exists; skip steps whose outputs are present
```

You need the **target** (HET residue name, or protein chain) and any **cofactors** that must
stay with the enzyme. If the orchestrator skill hasn't settled these with the user, ask now:
step C cannot guarantee a sub-unit contains a ligand nobody declared. Deposits often bind a
ligand in only one of several copies, or label it with the wrong chain.

## Run

```bash
gmxprep-fetch   --batch-dir B --ids 1ABC 2XYZ            # or --uniprot P12345 --method xray / --query TEXT
gmxprep-fix     --batch-dir B                            # all of raw/ -> fixed/
gmxprep-extract --batch-dir B --targets LIG --expect COF:1   # LIG@1ABC scopes a name to one entry
```

Check each JSON line before moving on:

| Command | Look at | Act when |
|---|---|---|
| fetch | `failed`, `counts` | `not found / obsolete`: tell the user. Reruns should report `cached`. |
| fix | per entry: `open_gaps`, `distinct_het_conformers`, `warnings` | Report open internal gaps near the binding site, and conformer verdicts |
| extract | per entry: `selected`, `inventory`, `warnings`, `files` | The target is missing from `inventory`, or there is any warning about a relabelled HET or a skipped altloc branch |

Open `units/<ID>/unit.json → candidates` only if the selection surprises you or there are
warnings. Report which unit was chosen and why the runner-up lost.

## Alternate conformers

Step B keeps the highest-occupancy conformer and stores the rest. Step C branches the matrix
(`complex_<T>_alt<L>.pdb`) only when a target's conformers are genuinely distinct
(symmetry-corrected RMSD ≥ 0.5 Å), clash-free, and inside the unit. Coordinates are never
averaged. Use `--altloc branch` to force branches, `--altloc primary` to suppress them.

## Hand-off

Done when `gmxprep-status` shows every entry with `unit` set and the expected `matrix_files`.
The next skill reads `units/<ID>/` directly. Tell the orchestrator:
- which entries are apo (no target);
- which entries branched (one manifest complex per branch, `variant: alt<L>`);
- any target that is a protein chain (`--target-chains`).
