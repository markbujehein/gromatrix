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
| `--only` | system ids (directory names under `systems/`) |
| `--deffnm` | file name of the **production** run (`md` → `md.tpr/xtc/edr/log/gro/cpt`). EM/NVT/NPT are always `em`, `nvt`, `npt`, because the later phases read them |
| `--threads` | `mdrun -nt`; with `--slurm`, also `--cpus-per-task` |
| `--gpu-id` | `mdrun -gpu_id` (local runs only) |
| `--maxh` | `mdrun -maxh`: stop cleanly with a checkpoint after H hours. `--slurm` derives it: walltime × 0.95 − 10 min |
| `--mdrun-args` | extra mdrun flags, e.g. `"-pin on -nb gpu -bonded gpu -update gpu"` |
| `--force` | also run systems whose step-F smoke test failed |
| `--poll` | seconds between progress lines |

### Protocol (mdp files from step F, `systems/<id>/mdp/`)

| Phase | grompp inputs | Notes |
|---|---|---|
| em | `-c solv_ions.gro -r solv_ions.gro` | reuses step F's gated `em.tpr` when present |
| nvt | `-c em.gro -r em.gro` | POSRES; new velocities |
| npt | `-c nvt.gro -r nvt.gro -t nvt.cpt` | POSRES; C-rescale barostat |
| md | `-c npt.gro -t npt.cpt` | unrestrained production |

All runs use `-n index.ndx -p topol.top`; `-maxwarn` is never passed.

### Gates

| Phase | Fails when |
|---|---|
| em | Fmax ≥ emtol, or Epot ≥ 0 |
| nvt, npt | fatal error, no `Finished mdrun`, or **any** LINCS warning |
| md | fatal error or no `Finished mdrun` (LINCS warnings are recorded in `warnings`) |

When a gate fails, the phase's `.gro` is renamed to `<phase>.failed.gro`, so a rerun never takes
the phase as complete.

### Output

- Progress lines: `{"step":"run","system","phase","status":"started|resuming|running|done|skipped","progress_pct","overall_pct"}`.
  `overall_pct` weighs the NVT/NPT/production phases by `nsteps`.
- Summary line: `{"ok","done":[…],"failed":[{id,status,reason}],"interrupted":[…],"busy":[…],"skipped":[{id,reason}]}`.
  `ok` is true when every system that ran reached `done`.
- `--background`: `{"ok","mode":"background","pid","host","log"}`. Progress goes to `B/run_logs/run_<time>.jsonl`.
- `--slurm`: `{"ok","mode":"slurm","script","systems","next":"sbatch …"}`.
- `systems/<id>/run.json`: `status`, `inputs_sha`, `deffnm`, `host`, `pid` (while running),
  `phases.<name>.status` (`done|running|failed|interrupted`, plus `resumed`/`skipped`), `warnings`,
  `reason`, `trajectory`.
- `B/run.json`: the last batch-wide run (foreground summary, or the background pid/log).
