# Where the MD runs

Steps A-F (structure preparation, validation, GROMACS builds and the smoke test) need only CPU and
GROMACS and take minutes. Only production MD (G) is expensive. Choose once per batch:

```bash
gromatrix-compute --batch-dir B detect      # cpus, gpu, slurm, gmx on this machine
gromatrix-compute --batch-dir B set --mode local  --threads 8 [--gpu-id 0] [--maxh 20]
gromatrix-compute --batch-dir B set --mode slurm  --time 24:00:00 --partition gpu --gpus-per-task 1
gromatrix-compute --batch-dir B set --mode remote --time 48:00:00
gromatrix-compute --batch-dir B run-args    # exact flags for gmxprep-run
```

| Mode | What happens | Use when |
|---|---|---|
| `local` | `gmxprep-run --background` on this machine; omit `--gpu-id` for CPU-only | short runs (equilibration, smoke-grade), or a local GPU |
| `slurm` | the Slurm array script is written here and submitted with `sbatch` from this machine | this machine is a cluster login node |
| `remote` | the array script is written, nothing is submitted; the batch is copied elsewhere | the cluster is a separate machine |

## Cost guide (rule of thumb, measure on your machine)
A 65k-atom system on 12 CPU cores does very roughly 10 ns/day, so 100 ns takes about 10 days. The
100 ps NVT + 100 ps NPT + 0.1 ns production protocol of a smoke-grade run takes about 35 minutes.
A single modern GPU is typically 10-30 times faster. The production length is fixed at build time
(`gmxprep-gmx build --md-ns N`); rebuild with `--force --only <id>` to change it.

## Remote hand-off
1. Finish A-F locally and run `gmxprep-gmx smoke`.
2. `gromatrix-compute set --mode remote ...`, then the md-runner writes `07-slurm/run_array.sbatch`.
3. Copy the batch to the cluster (it needs the same pixi or conda environment with GROMACS):
   `rsync -a --exclude 01-raw --exclude 05-alphafold-cache B/ cluster:path/B/`
4. On the cluster: `sbatch 07-slurm/run_array.sbatch`. Resubmitting the same script resumes from checkpoints.
5. Copy `06-systems/` back and let the md-analyst run step H locally.
Absolute paths recorded inside old logs do not matter; `gmxprep-run` resolves everything from the batch directory.
