---
name: execution
description: "Runs the MD protocol (energy minimisation, NVT, NPT, production) for built GROMACS systems on this machine, detached in the background, or as a Slurm array job, with checkpoint resume and per-phase gates (step G). Use this skill ONLY inside the gromatrix:md-runner subagent, or when the user asks to diagnose, resume or script an MD run; the orchestrator must never run gmxprep-run itself. Covers gmxprep-run, the compute choice from gromatrix-compute (CPU-only, local, Slurm, remote), resume rules, EM-convergence and LINCS failures, Slurm scripts. Not for building systems or analysing trajectories."
---

# gromatrix: MD execution (step G)

**This step runs inside the `gromatrix:md-runner` subagent.** The orchestrator delegates here so the
log noise stays out of its context; the subagent returns one JSON payload
(`${CLAUDE_PLUGIN_ROOT}/references/subagents.md`). Where MD runs, and with how many threads or which
GPU, comes from `00-config/compute.json` (`gromatrix-compute --batch-dir B show`), not from asking
again. No GPU is fine: omit `--gpu-id`. Compute modes: `${CLAUDE_PLUGIN_ROOT}/references/compute-modes.md`.

Command: `gmxprep-run`. Full flags and output fields are in `references/cli.md`. Failure reasons
and their fixes are in `references/troubleshooting.md`.

## Before you start

```bash
gmxprep-status --batch-dir B
```

Step G needs built systems (`systems.built > 0`). It runs only systems whose `system.json` is
`ok`. A system whose smoke test **failed** is skipped unless `--force`. Don't force it: the smoke
failure is a finding for the gmx builder.

Settle with the user before running:
1. **Where it runs.** This machine (foreground or `--background`) or a cluster (`--slurm`).
2. **Hardware.** `--threads N` and, locally, `--gpu-id`. Under Slurm the scheduler assigns GPUs;
   use `--gpus-per-task`.
3. **Length.** Set at build time by step F (`--md-ns`, default 100 ns). Step G does not change it.

## Run

```bash
gmxprep-run --batch-dir B --threads 8 --gpu-id 0             # foreground, streams progress
gmxprep-run --batch-dir B --threads 8 --background           # detached; returns the pid and a log path
gmxprep-run --batch-dir B --slurm --time 24:00:00 --threads 16 --gpus-per-task 1 \
            --partition gpu --prelude "module load gromacs"  # writes B/07-slurm/run_array.sbatch
```

- **Slurm.** Hand the user `sbatch B/07-slurm/run_array.sbatch`; the script never submits itself.
  Each array task runs one system with a `--maxh` budget just under the walltime, shared by all its
  phases. A task that runs out of time stops cleanly at a checkpoint, and resubmitting the same
  script resumes every system.
- **Progress.** Lines look like `{"step":"run","system":…,"phase":"npt","status":"running",
  "progress_pct":45,"overall_pct":61}`. The last line is the summary `{"ok":…}`. In background mode,
  read the log, or run `gmxprep-status`, which reads each `06-systems/<id>/run.json`.

## Resume

Run the same command again. Each phase is either
- **done**: its gate is checked again, then it is skipped. That includes an EM left by step F's
  smoke test;
- **partial** (it has a `.cpt`): continued with `mdrun -cpi`;
- **new**: run from scratch with `grompp`, then `mdrun`.

Once one phase reruns, every phase after it starts fresh. A phase that failed its gate is cleared
(its evidence is kept as `<phase>.failed.gro/.log`), so the rerun starts it fresh. Step F refuses
to rebuild a system that has MD results unless `--force`, because rebuilding deletes them. Never
delete checkpoints by hand to "restart cleanly". Never hand-edit `run.json`.

## Read the results

`done` / `failed` / `interrupted` / `busy` per system; the reason is in `06-systems/<id>/run.json`,
the tool output is in `run.log`, and each phase has its own `<phase>.log`.
- **`interrupted`** (signal, Ctrl-C, `-maxh`): not an error. Rerun to resume.
- **`failed`**: a gate stopped this system and the others carried on. The fix is upstream: a bad
  pose, a clash or wrong parameters. See `references/troubleshooting.md`. Do not add `-maxwarn`,
  and don't shorten the timestep to make warnings go away.
- **`busy`**: another `gmxprep-run` owns the system. Don't start a second one. This is a check,
  not a lock, so never start two runners on the same systems on purpose.
- **`not_started`**: the batch stopped (signal or `--maxh`) before reaching these systems.

## Hand-off

Done when `gmxprep-status` shows `steps.G = "done"`, i.e. every runnable system's `run.json` says
`done`. Report any production LINCS warnings (in `run.json` `warnings`). `next` then names
`gromatrix:analysis`.
