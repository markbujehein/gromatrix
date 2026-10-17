# CLI reference: step G

## gmxprep-run

```
gmxprep-run --batch-dir B [--only ID ...] [--deffnm md] [--threads N] [--gpu-id ID] [--maxh H]
            [--mdrun-args "..."] [--gmx gmx] [--force] [--poll 30]
            [--background | --slurm [--time 24:00:00] [--partition P] [--account A]
                                    [--gpus-per-task N] [--max-parallel N] [--job-name gromatrix]
                                    [--prelude "module load gromacs" ...]]
```

| Flag | Meaning |
|---|---|
| `--only` | system ids (directory names under `06-systems/`) |
| `--deffnm` | file name of the **production** run (`md` → `md.tpr/xtc/edr/log/gro/cpt`). EM/NVT/NPT are always `em`, `nvt`, `npt`, because the later phases read them |
| `--threads` | `mdrun -nt`; with `--slurm`, also `--cpus-per-task` |
| `--gpu-id` | `mdrun -gpu_id` (local runs only) |
| `--maxh` | time budget in hours for the whole invocation. Each `mdrun` gets `-maxh` = what is left, so the run stops cleanly with a checkpoint before the budget ends. `--slurm` derives it: walltime × 0.95 − 10 min |
| `--mdrun-args` | extra mdrun flags, e.g. `"-pin on -nb gpu -bonded gpu -update gpu"` |
| `--force` | also run systems whose step-F smoke test failed |
| `--poll` | seconds between progress lines |

### Protocol (mdp files from step F, `06-systems/<id>/mdp/`)

| Phase | grompp inputs | Notes |
|---|---|---|
| em | `-c solv_ions.gro -r solv_ions.gro` | reuses step F's gated `em.tpr` when present; an EM already run by step F's smoke test counts as done (and is gated) |
| nvt | `-c em.gro -r em.gro` | POSRES; new velocities |
| npt | `-c nvt.gro -r nvt.gro -t nvt.cpt` | POSRES; C-rescale barostat |
| md | `-c npt.gro -t npt.cpt` | unrestrained production |

All runs use `-n index.ndx -p topol.top`; `-maxwarn` is never passed.

### Gates

| Phase | Fails when |
|---|---|
| em | Fmax ≥ emtol, or a positive Epot |
| nvt, npt | fatal error, no `Finished mdrun`, or **any** LINCS warning in any log segment |
| md | fatal error or no `Finished mdrun` (LINCS warnings are recorded in `warnings`) |

Gates also run on phases that are skipped as done. When a gate fails, the phase's `.gro` and `.log`
are kept as `<phase>.failed.gro/.log` and the phase is cleared, checkpoint included, so a rerun
starts it fresh.

A phase counts as done when its `.gro` exists, its newest log segment finished without a `-maxh`
or signal stop, and the checkpoint the next phase reads (`nvt.cpt`, `npt.cpt`) exists. An `mdrun`
that exits non-zero is `failed`, never `interrupted`, whatever an older log segment says.

### Output

- Progress lines: `{"step":"run","system","phase","status":"started|resuming|running|done|skipped","progress_pct","overall_pct"}`.
  `overall_pct` weighs the NVT/NPT/production phases by `nsteps`.
- Summary line: `{"ok","done":[…],"failed":[{id,status,reason}],"interrupted":[…],"not_started":[…],"busy":[…],"skipped":[{id,reason}]}`.
  `ok` is true when every system that ran reached `done`. A signal or the end of the `--maxh` budget
  stops the batch; the remaining systems are listed in `not_started`.
- `--background`: `{"ok","mode":"background","pid","host","log"}`. Progress goes to `B/run_logs/run_<time>.jsonl`.
- `--slurm`: `{"ok","mode":"slurm","script","systems","next":"sbatch …"}`.
- `06-systems/<id>/run.json`: `status`, `inputs_sha`, `deffnm`, `host`, `pid` (while running),
  `phases.<name>.status` (`done|running|failed|interrupted`, plus `resumed`/`skipped`), `warnings`,
  `reason`, `trajectory`.
- `B/07-run-logs/run.json`: the last invocation (foreground summary including `only`, or the background pid/log).
- `gmxprep-status` reports a `running` system whose runner died on this host as `interrupted`.
