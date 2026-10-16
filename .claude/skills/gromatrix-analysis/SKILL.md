---
name: gromatrix-analysis
description: Step H of the gromatrix pipeline. Analyses every finished production trajectory of a gromatrix batch, covering PBC removal and fitting with gmx trjconv, protein Cα RMSD, ligand/target heavy-atom RMSD in the protein frame, per-residue Cα RMSF, radius of gyration and target-enzyme hydrogen-bond persistence (MDAnalysis), with flags for ligand drift and protein instability, publication-ready plots and one comparison table for the whole receptor × target matrix. Use this when the user wants MD results analysed or compared, asks whether a ligand stayed bound, which H-bonds persist or how stable a complex was, or when `gmxprep-status` names this skill as `next`.
---

# gromatrix: trajectory analysis (step H)

Command: `gmxprep-analyze`. Full flags and the `analysis.json` schema are in `references/cli.md`.
How to read the metrics and flags, and their pitfalls, is in `references/interpretation.md`.

## Before you start

```bash
gmxprep-status --batch-dir B
```

Step H analyses systems whose `run.json` says `done` and that have a production trajectory
(`md.xtc`, or the run's `--deffnm`). Systems still running or failed are skipped and listed.
Needs GROMACS (`trjconv`, `convert-tpr`, `editconf`) plus MDAnalysis and matplotlib.

## Run

```bash
gmxprep-analyze --batch-dir B                    # every finished system
gmxprep-analyze --batch-dir B --stride 5 --begin 10000   # every 5th frame, skip the first 10 ns
```

Pass the same `--deffnm` as `gmxprep-run` if the user changed it. A system whose `analysis.json` is
newer than its trajectory and was made with the same settings is not recomputed (`--force`
overrides). `analysis_summary.csv` is always rebuilt from every system's results, so `--only` never
drops rows.

## What it does

1. **PBC.** The molecules are made whole and centred (`-pbc mol -center` on Protein, or on the
   ligand in a ligand-only system), then fitted (`-fit rot+trans` on the backbone). Index groups
   are selected by **name**, from step F's `index.ndx` plus Backbone and C-alpha. Only the solute
   is kept.
2. **Metrics**, relative to the first production frame:
   - Cα RMSD;
   - target heavy-atom RMSD *without refitting the target*, which measures pose retention;
   - Cα RMSF per residue;
   - radius of gyration;
   - target ↔ enzyme H-bonds (D–A < 3.0 Å, D–H–A > 150°), as an occupancy per donor/acceptor pair
     and a count per frame.
3. **Which atoms are the target** comes from step F's `system.json`: `target_resname` for a
   ligand, `target_residues` for a protein target that pdb2gmx merged into the enzyme molecule.

## Read the results

- `B/analysis_summary.csv`: one row per system. Use it for the matrix table in the final report:
  last-half mean RMSDs, mean H-bond count, top H-bond and its occupancy, flags.
- `systems/<id>/analysis/analysis.json` holds the full series. The plots are `rmsd.png`,
  `rmsf.png`, `rg.png` and `hbonds.png`.
- **Flags** (report them verbatim; they are findings, not errors):
  - `ligand_drift`: the pose was not retained;
  - `protein_unstable`: the protein drifted from its starting structure;
  - `no_target_hbonds`: a complex with no persistent H-bond to the enzyme.

  See `references/interpretation.md` before drawing conclusions. One replica of a given length
  is anecdotal; say so.

## Hand-off

Done when `gmxprep-status` shows `steps.H = "done"`. Return to the orchestrator with:
- `analysis_summary.csv`;
- the flags per system;
- the plot paths;
- any systems skipped or failed, and why.
