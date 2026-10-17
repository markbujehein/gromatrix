---
name: analysis
description: "Analyses finished production trajectories and explains the results (step H): periodic-boundary fixing, protein Calpha RMSD, target RMSD in the protein frame, per-residue RMSF, radius of gyration, target-enzyme H-bond persistence, drift and instability flags, plots with a README, and a matrix summary table. Use this skill ONLY inside the gromatrix:md-analyst subagent, or when the user asks to interpret or compare MD results, whether a ligand stayed bound, which H-bonds persist or how stable a complex was; the orchestrator must never run gmxprep-analyze itself. Covers gmxprep-analyze, gromatrix-report, analysis.json and honest interpretation. Not for running MD."
---

# gromatrix: trajectory analysis (step H)

**This step runs inside the `gromatrix:md-analyst` subagent.** It returns one JSON payload
(`${CLAUDE_PLUGIN_ROOT}/references/subagents.md`). Every `analysis/` folder gets a `README.md`
explaining each plot, and `gromatrix-report --batch-dir B` writes `09-reports/README.md` with
tables and figures for every step. Write the final report in the format of
`${CLAUDE_PLUGIN_ROOT}/examples/analysis_summary_example.md`; biophysical limits are in
`${CLAUDE_PLUGIN_ROOT}/references/domain-rules.md`.

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
overrides). `08-analysis/analysis_summary.csv` is always rebuilt from every system's results, so `--only` never
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

- `B/08-analysis/analysis_summary.csv`: one row per system. Use it for the matrix table in the final report:
  last-half mean RMSDs, mean H-bond count, top H-bond and its occupancy, flags.
- `06-systems/<id>/analysis/analysis.json` holds the full series. The plots are `rmsd.png`,
  `rmsf.png`, `rg.png` and `hbonds.png`.
- **Flags** (report them verbatim; they are findings, not errors):
  - `ligand_drift`: the pose was not retained;
  - `protein_unstable`: the protein drifted from its starting structure;
  - `no_target_hbonds`: a complex with no persistent H-bond to the enzyme.

  See `references/interpretation.md` before drawing conclusions. One replica of a given length
  is anecdotal; say so.

## Hand-off

Done when `gmxprep-status` shows `steps.H = "done"`. Return to the orchestrator with:
- `08-analysis/analysis_summary.csv`;
- the flags per system;
- the plot paths;
- any systems skipped or failed, and why.
