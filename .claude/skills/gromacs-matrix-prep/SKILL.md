---
name: gromacs-matrix-prep
description: Prepare PDB structures of enzymes and docked targets (ligands or proteins) for GROMACS molecular dynamics as a full simulation matrix (enzyme alone, target alone, enzyme+target complex for every pair) with one forcefield locked per batch. Fetches from RCSB, fixes structures with PDBFixer, extracts subunits with their cofactors, relaxes with OpenMM, compares with AlphaFold (RMSD, TM-score, lDDT), then builds topologies, boxes, solvent and ions using AMBER (GAFF2/acpype) or CHARMM36m (CGenFF). Use this whenever the user wants to set up MD/GROMACS runs from PDB IDs or docked poses, parameterise a ligand or cofactor for GROMACS, compare enzyme-substrate systems, check forcefield consistency across simulations, or asks about pdb2gmx/acpype/CGenFF topology problems in this repository, even if they don't say "matrix".
---

# GROMACS matrix preparation

You drive tested CLI scripts. Do not re-implement their logic inline: the scripts carry the
forcefield lock, the shared ligand-parameter cache and the validation gates that keep the matrix
scientifically comparable, and ad-hoc code bypasses all three.

Every command prints a single JSON line (`{"ok": true|false, ...}`) and exits non-zero on failure.
Read that line; open log files only when it reports a failure.

## 0. Setup (once per machine)

```bash
conda env create -f environment.yml && conda activate gmxprep   # from the repo root
python -m pytest -q tests                                       # sanity check
```

If a `gmxprep-*` command is missing, that module is not implemented yet. Tell the user which step is
unavailable rather than improvising a replacement script.

## 1. Settle the decisions that are hard to undo

Before running anything, confirm with the user (or infer from their request):

1. **Batch directory**: one batch = one comparison matrix.
2. **Forcefield**: default `amber` (amber99sb-ildn + GAFF2/AM1-BCC + TIP3P). Use `charmm`
   (CHARMM36m + CGenFF) only if they ask or their earlier work uses it. For the trade-offs and the
   newer AMBER options, read `references/forcefields.md`.
3. **Which ligands/cofactors** need parameters, and their **protonation state** (SMILES with explicit
   charges, or an SDF with hydrogens). This is the scientific choice the user most often forgets.

Then lock the batch. The lock is permanent for this batch: a different forcefield needs a new batch dir.

```bash
gmxprep-gmx lock --batch-dir B --forcefield amber
```

## 2. Pipeline

Run the steps in order. Each step reads the previous step's directory. Full flags are in `references/cli.md`.

| Step | Command | Check before moving on |
|---|---|---|
| A fetch | `gmxprep-fetch --batch-dir B --ids 1TQN 3QM4` | every ID has a file in `B/raw/` |
| B fix | `gmxprep-fix --batch-dir B` | `open_gaps` and `warnings` (alternate ligand conformers, unmodelled loops) |
| C extract | `gmxprep-extract --batch-dir B --targets RES` | `selected`, `inventory` hold the target and cofactors |
| D relax | `gmxprep-relax --batch-dir B --in units --out relaxed` | no unit over `--max-rmsd` |
| E AF compare | `gmxprep-afcompare --batch-dir B --in relaxed --out af_report.csv` | flag TM < 0.8 or lDDT < 0.7 to the user |
| F plan | `gmxprep-gmx plan --batch-dir B --manifest M` | system list = enzymes + targets + all pairs |
| F build | `gmxprep-gmx build --batch-dir B --manifest M` | `"ok": true`, then read `warnings` |
| F verify | `gmxprep-gmx verify --batch-dir B` | `"ok": true` and `problems` is empty |
| F smoke | `gmxprep-gmx smoke --batch-dir B` (needs GROMACS run time, ~minutes per system) | every system `ok`. Otherwise read `reason`. |

**Step C needs to know the target.** Ask which HET residue (or chain) is the target, and which
cofactors must stay with the enzyme (`--expect`). Without that, C cannot guarantee that the chosen
sub-unit contains them. Some deposits bind a ligand in only one copy (1XVP: CITCO in chain D only).
After C, read `units/<ID>/unit.json` → `candidates` only when the selection surprises you or the
JSON reports warnings. Report the selected unit and why the runner-up lost.

**Alternate conformers.** Step B reports them, and Step C decides. Conformers are never averaged:
the mean of two partial-occupancy states is not a physical structure. If C branched the matrix
(`complex_<T>_altA/B.pdb`), add one manifest complex per branch with `variant: altA` etc. Report
the symmetry-corrected RMSD to the user when they ask why something was or wasn't branched.

Steps A–E are useful but optional: if the user already has relaxed, docked PDBs, go straight to F.

## 3. Writing the manifest for step F

Copy `examples/manifest.yaml` and edit it. Paths are relative to the manifest. Rules the build enforces:

- Every enzyme × target pair needs a `complexes` entry (the docked pose). Missing pairs abort the build.
  `--allow-incomplete` exists, but ask the user first because it breaks the matrix.
- `targets[].kind` is `ligand` (needs `resname`) or `protein` (needs `chains` present in the complex).
- Every non-ion HET residue name in any input needs a `ligands:` entry with `smiles` + `charge` or an
  explicit-H `sdf`. The first pose seen becomes the batch template; every other pose is mapped onto it,
  so all systems share one topology per ligand.
- An explicit-H docked pose for a complex can be given as `complexes[].sdf`; otherwise hydrogens are rebuilt.

## 4. Reading build results

- `failed[]`: each entry has the system id and the error. The full traceback is in
  `B/systems/<id>/error.txt` and tool output is in `build.log`. Fix the cause, then rerun with
  `--only <id>`. Systems already built are skipped unless you pass `--force`.
- `warnings[]` must be shown to the user. In particular: auto-guessed protonation, and CGenFF
  penalties from 10 to 50.
- Hard stops are deliberate. Do not paper over them with flags. Explain the stop to the user:
  - `ForcefieldMismatch`: someone tried to mix forcefields. Use a new batch dir.
  - metal-containing ligand (Fe, Cu, Zn bound to an organic residue, e.g. HEM): needs MCPB.py or a
    dedicated metal model, which is outside this pipeline.
  - CGenFF penalty > 50: needs ffTK/reparameterisation. `--accept-high-penalty` is for exploratory runs only.
  - grompp warnings at the `em.tpr` gate: read them. Raising `--maxwarn` hides real errors.

For specific error messages, read `references/troubleshooting.md`.

## 5. Report back

When F succeeds, give the user:

- forcefield key, number of systems built, and the matrix as a small table (enzyme rows × target
  columns, plus the alone systems)
- all warnings, verbatim
- the AlphaFold comparison outliers from step E, if it was run
- the smoke result per system, and the next command: `B/systems/<id>/run.sh` (em→nvt→npt→md with
  correct `-r`/`-t` chaining)

Do not start long `mdrun` jobs unless the user asks.
