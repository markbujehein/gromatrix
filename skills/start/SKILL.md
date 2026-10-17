---
name: start
description: "Entry point and router for the gromatrix pipeline, which takes PDB/mmCIF structures or docked poses through GROMACS preparation, MD and trajectory analysis (steps A-H) for a receptor by target matrix under one locked forcefield. Use this skill FIRST whenever the user wants MD systems prepared, simulated or analysed from PDB IDs or docked complexes, wants several enzymes compared against several ligands, asks where a batch is or what comes next, or wants the whole workflow in one go. Covers the forcefield lock, target and cofactor choices, the local/Slurm/remote compute choice, routing to the five stage skills in order, delegating the long and noisy steps G and H to subagents, and the final report format."
---

# gromatrix: orchestrator

gromatrix is a set of CLI commands covering steps A–H. This skill decides the order and the
hard-to-undo choices. The five stage skills each document only their own commands; load one only
when the batch is at that stage. Shared reference material (read on demand, not up front):
`${CLAUDE_PLUGIN_ROOT}/references/` holds `batch-layout.md` (what is in each folder),
`domain-rules.md` (biophysical rules), `compute-modes.md` and `subagents.md`.

## State lives in the batch directory

Begin every session, and every stage, with:

```bash
gmxprep-status --batch-dir B
```

It prints one JSON line: `forcefield`, `compute`, `steps` (`A…H`: done / partial / pending / n/a),
`structures`, `systems`, and **`next`** = `{skill, agent?, command, why}`. Follow `next`: load the
named skill and run its command; if `next.agent` is set, delegate (below). Run `gmxprep-status`
again after each stage. The files are the contract; no skill keeps state in conversation memory.

The batch folders are numbered in pipeline order: `00-config`, `01-raw` … `04-relaxed`,
`05-alphafold`, `06-params`, `06-systems`, `07-run-logs`, `08-analysis`, `09-reports`.

## Decisions to settle with the user before any command

1. **Batch directory.** One batch is one comparison matrix. Ask where; never pick a generic
   folder name like `work`.
2. **Forcefield.** Default `amber` (amber99sb-ildn + GAFF2/AM1-BCC + TIP3P); `charmm` (CHARMM36m +
   CGenFF) if the user asks. Lock it before step D: `gmxprep-gmx lock --batch-dir B --forcefield amber`.
   The lock is permanent: another forcefield needs a new batch directory.
3. **Targets and cofactors.** Which HET residue (or protein chain) is the target, which cofactors
   stay bound, and each ligand's protonation state (SMILES with charges, or an SDF with explicit H).
4. **Matrix completeness.** Crystal structures give only native pairs; a full matrix needs docked
   cross pairs. `--allow-incomplete` needs the user's explicit consent.
5. **Compute** (before G, not before A–F, which need only CPU and GROMACS): run
   `gromatrix-compute --batch-dir B detect`, then ask where the MD runs and record it with
   `gromatrix-compute --batch-dir B set --mode local|slurm|remote …`. No GPU is fine for A–F and
   for short local runs; see `${CLAUDE_PLUGIN_ROOT}/references/compute-modes.md`. Production length
   is fixed at build time (`--md-ns`, default 100 ns): confirm it against the system size and the
   chosen compute before building.

## Order

```
lock → structure-prep (A fetch, B fix, C extract)
     → validate       (D relax, E AlphaFold)             [optional, recommended]
     → gmx-builder    (F manifest, build, verify, smoke)
     → execution      (G EM → NVT → NPT → production)    [subagent: gromatrix:md-runner]
     → analysis       (H RMSD, RMSF, Rg, H-bonds)        [subagent: gromatrix:md-analyst]
```

Skip A–C when the user has clean docked PDBs; skip D–E only if the user says so. Stage skills are
invoked as `gromatrix:structure-prep`, `gromatrix:validate`, `gromatrix:gmx-builder`,
`gromatrix:execution`, `gromatrix:analysis`.

## Steps G and H run in subagents, never in this context

**The rule:** never run `gmxprep-run` (G) or `gmxprep-analyze` (H) here, including through
`gmxprep-pipeline --run-md` or `--analyze`. They produce thousands of log lines and run for hours;
the subagent absorbs that and returns one small JSON payload.

Spawn the agent named in `next.agent`, in the background, with the Agent tool:

```
Agent(
  description = "Run MD for <batch name>",
  subagent_type = "gromatrix:md-runner",      # H: "gromatrix:md-analyst"
  run_in_background = true,
  prompt = "batch_dir: /abs/path/B\naction: start   # start | resume | check\nReturn only the JSON payload."
)
```

Then carry on with the user or stop. You are notified when it finishes; don't poll it, don't read
its transcript or output file. Its final message is the payload (schemas in
`${CLAUDE_PLUGIN_ROOT}/references/subagents.md`): `ok`, per-system status, the failed phase and a
one-line reason, LINCS warning counts, and the path of the evidence. If `ok` is false, report the
reason verbatim and load `gromatrix:execution` (or `gromatrix:analysis`) to decide the fix. Only
open a log yourself for the one failing system, and read just its tail.

Where the agent tool is unavailable (some chat clients): start G with `gmxprep-run … --background`
(flags from `gromatrix-compute run-args`) and check `gmxprep-status` at long intervals; never tail logs.

## One-command mode

For unattended preparation, write a pipeline manifest from `references/pipeline-manifest.md` and
preview it first:

```bash
gmxprep-pipeline --manifest M --dry-run     # exact commands, nothing runs
gmxprep-pipeline --manifest M --to F        # A–F; stops at the first failed step
```

It runs the same stage commands, so every gate still applies. It ends by writing
`B/09-reports/README.md`. When it stops, `failed_stage` names the step: load that stage skill,
fix the cause and rerun with `--from <letter>`. Then continue with G and H through the subagents.
Flags: `--skip-relax`, `--skip-afdb`, `--skip-smoke`, `--keep-going`.

## Rules for every stage

- Run the commands; don't re-implement their logic: they carry the forcefield lock, the shared
  ligand-parameter cache and the validation gates.
- Every command prints one JSON line (`{"ok": ...}`). Read it. Open logs only on failure.
- A hard stop (forcefield mismatch, incomplete matrix, metal-containing ligand, high CGenFF penalty,
  failed relax gate) is a finding to report, not something to bypass with a flag.
- Don't invent data. Results come from the files; if a step wasn't run, say so.
- If `gmxprep-status` is not found, the package isn't installed where commands run. Inside a clone,
  `pixi install -e md` then `pixi run -e md …`. Elsewhere:
  `conda create -n gromatrix -c conda-forge python=3.12 openmm pdbfixer rdkit biopython scipy pyyaml matplotlib mdanalysis gromacs acpype openbabel`,
  then `pip install "gromatrix @ git+https://github.com/markbujehein/gromatrix"`. Never improvise a replacement script.

## Final report

Run `gromatrix-report --batch-dir B` (the subagent does this after H), then write the report in the
format of `${CLAUDE_PLUGIN_ROOT}/examples/analysis_summary_example.md`: forcefield key, the matrix as
a table, the sampling statement, and every warning and flag verbatim (C unit choices and altloc
branches, D failures, E flags, F warnings, G failures and LINCS warnings, H flags). The example's
numbers are synthetic: never copy them. State the sampling honestly (one replica of N ns per
system) and finish with what is left, which is `gmxprep-status`'s `next`.
