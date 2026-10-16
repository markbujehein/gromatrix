---
name: gromatrix
description: Orchestrator for the gromatrix toolkit, which takes PDB/mmCIF structures or docked poses end to end through a GROMACS simulation matrix (receptor alone, target alone, every receptor+target complex) with one forcefield locked per batch, covering preparation, MD runs and trajectory analysis. Use this first whenever the user wants MD/GROMACS systems prepared, simulated or analysed from PDB IDs or docked complexes, wants to compare several enzymes against several ligands/targets, wants the whole pipeline run in one go, asks "where is my batch at" or "what's next", or mentions the gromatrix pipeline as a whole. It decides the forcefield, then routes to the stage skills gromatrix-structure-prep (fetch/fix/extract), gromatrix-validate (OpenMM relax + AlphaFold comparison), gromatrix-gmx-builder (GROMACS topologies), gromatrix-execution (MD runs) and gromatrix-analysis (trajectory metrics), in that order.
---

# gromatrix orchestrator

gromatrix is a set of CLI commands covering steps A–H. This skill decides the order and the
hard-to-undo choices. The five stage skills each document only their own commands. Load a stage
skill only when the batch is at that stage.

## How state is handed between skills

All state lives in one **batch directory**. No skill keeps anything in conversation memory that
the next skill needs. Begin every session, and every stage, with:

```bash
gmxprep-status --batch-dir B
```

It prints one JSON line containing:
- `forcefield`: the lock;
- `steps`: `{A…H: done | partial | pending | n/a}`, the batch at a glance;
- `structures.<ID>`: `fetched`, `fixed`, `unit`, `targets`, `matrix_files`, `relax`, `af`;
- `systems`: built/failed/smoke counts, `run` (ids per MD status), `run_done`, `analysed`,
  `analysis_flags`;
- **`next`**: `{skill, command, why}`.

Follow `next`: load the named skill and run its command. Run `gmxprep-status` again after each
stage. The status command reads only these artefacts, so the files *are* the contract:

| Stage (skill) | Writes | The next stage reads |
|---|---|---|
| lock (this skill) | `forcefield.lock.json` | D (family), F (everything) |
| A–C (`gromatrix-structure-prep`) | `raw/`, `fixed/<ID>.pdb` + `.meta.json`, `units/<ID>/unit.json` + `enzyme.pdb` / `target_<T>.pdb` / `complex_<T>[_alt<L>].pdb` | D, E, and the F manifest |
| D–E (`gromatrix-validate`) | `relaxed/<ID>/*.pdb` + `relax.json`, `af/<ID>.json`, `af_report.csv` | F (use relaxed files that passed) |
| F (`gromatrix-gmx-builder`) | `params/<RES>/`, `systems/<ID>/` (`system.json`, `topol.top`, `index.ndx`, `mdp/`, `em.tpr`), `matrix.json` | G |
| G (`gromatrix-execution`) | `systems/<ID>/run.json`, `{em,nvt,npt,md}.{tpr,cpt,gro,log,edr}`, `md.xtc`; `slurm/run_array.sbatch` | H |
| H (`gromatrix-analysis`) | `systems/<ID>/analysis/` (`analysis.json`, `*.png`), `analysis_summary.csv` | the final report |

## Decisions to settle with the user before any command

1. **Batch directory.** One batch is one comparison matrix.
2. **Forcefield.** Default `amber` (amber99sb-ildn + GAFF2/AM1-BCC + TIP3P). Use `charmm`
   (CHARMM36m + CGenFF) if the user asks or their earlier work uses it. Lock it before step D:
   `gmxprep-gmx lock --batch-dir B --forcefield amber`. The lock is permanent: a different
   forcefield needs a new batch directory.
3. **Targets and cofactors.** Which HET residue name (or protein chain) is the target, which
   cofactors must stay bound, and each ligand's protonation state (SMILES with charges, or an SDF
   with explicit hydrogens).
4. **Matrix completeness.** Crystal structures give only native pairs. A full receptor × target
   matrix needs docked complexes for the cross pairs. Running the native pairs only
   (`--allow-incomplete`) needs the user's explicit consent.
5. **Compute** (only if MD will run): this machine or Slurm, threads/GPU, and production length
   (`--md-ns` at build time, default 100 ns per system).

## Order

```
lock → gromatrix-structure-prep (A fetch, B fix, C extract)
     → gromatrix-validate      (D relax, E AlphaFold)          [optional, recommended]
     → gromatrix-gmx-builder   (F manifest, build, verify, smoke)
     → gromatrix-execution     (G EM → NVT → NPT → production)  [hours to days of compute]
     → gromatrix-analysis      (H RMSD, RMSF, Rg, H-bonds)
```

Skip A–C when the user already has clean, docked PDBs; go straight to the gmx builder. Skip D–E
only if the user says so. G is the expensive step: confirm compute and length before starting it.

## One-command mode

When the user wants everything run unattended, write a pipeline manifest from the template in
`references/pipeline-manifest.md` (`batch_dir`, `forcefield`, `steps.<name>` flags, `ligands`,
optional `targets`) and preview it first:

```bash
gmxprep-pipeline --manifest M --dry-run                  # exact commands, nothing runs
gmxprep-pipeline --manifest M --run-md --analyze         # A–H; stops at the first failed step
```

It runs the same stage commands in order, so every gate still applies. Without `enzymes:` it
writes `B/manifest.generated.yaml` from the step-C units, using relaxed files that passed.
Useful flags:
- `--skip-relax`, `--skip-afdb`, `--skip-smoke`;
- `--from`/`--to` with step letters;
- `--keep-going`, to carry on past failures.

When it stops, `failed_stage` names the step. Load that stage's skill to diagnose it, fix the
cause, and rerun with `--from <letter>`; finished work is reused. For long MD on a cluster, use
the pipeline `--to F` and then `gmxprep-run --slurm` (gromatrix-execution).

## Rules for every stage

- Run the commands. Don't re-implement their logic: they carry the forcefield lock, the shared
  ligand-parameter cache and the validation gates.
- Every command prints one JSON line (`{"ok": ...}`). Read it. Open logs only on failure.
- A hard stop (forcefield mismatch, incomplete matrix, metal-containing ligand, high CGenFF
  penalty, failed relax gate) is a finding to report, not something to bypass with a flag.
- If `gmxprep-status` is not found, the package is not installed where commands run. Give the
  user the install. No clone is needed:
  `conda create -n gromatrix -c conda-forge python=3.12 openmm pdbfixer rdkit biopython scipy pyyaml`,
  then `pip install "gromatrix @ git+https://github.com/markbujehein/gromatrix"`. Step F also needs
  `gromacs acpype openbabel` from conda-forge. Inside a clone, `pixi shell -e md` does all of this.
  Never improvise a replacement script.

## Final report to the user

Report the forcefield key, then the matrix as a table (receptor rows × target columns, plus the
alone systems), with the status of each system. Use the furthest step reached: built, MD done,
or analysed. Once analysed, fill the cells from `B/analysis_summary.csv`: last-half target RMSD,
mean H-bond count, top H-bond and its occupancy. Then, verbatim:
- every warning from C (unit choice, altloc branches, relabelled ligands),
- every D failure label,
- every E flag,
- every F warning,
- every G failure and production LINCS warning,
- every H flag (`ligand_drift`, `protein_unstable`, `no_target_hbonds`).

State the sampling honestly: one replica of N ns per system. Finish with what is left, which is
`gmxprep-status`'s `next`.
