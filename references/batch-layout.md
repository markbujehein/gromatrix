# Batch layout: the contract between steps

All state lives in one batch directory `B`. Folders are numbered in pipeline order; folders that
belong to one stage share a number. Names come from `gromatrix.layout`, the only place that defines
them. `gmxprep-status --batch-dir B` reads exactly these files, so the files are the hand-off.

| Folder or file | Written by | Read by | Holds |
|---|---|---|---|
| `00-config/forcefield.lock.json` | lock | D, F | the forcefield identity, set once, permanent |
| `00-config/compute.json` | `gromatrix-compute set` | G agent, status | `mode` local / slurm / remote, threads, gpu_id, Slurm options |
| `00-config/manifest.generated.yaml`, `pipeline.json` | `gmxprep-pipeline` | F, humans | derived step-F manifest; per-stage result of a pipeline run |
| `01-raw/` | A fetch | B, status | `<ID>.cif`, `cache.json` (sha256, revision) |
| `02-fixed/` | B fix | C | `<ID>.pdb`, `<ID>.meta.json` (gaps, removals, altlocs) |
| `03-units/<ID>/` | C extract | D, E, F | `unit.json`, `enzyme.pdb`, `target_<T>.pdb`, `complex_<T>[_alt<L>].pdb` |
| `04-relaxed/<ID>/` | D relax | E, F | relaxed `*.pdb`, `*.failed.pdb`, `relax.json` |
| `05-alphafold/` | E afcompare | status, report | `<ID>.json`, `af_report.csv`, `<ID>/AF_*.pdb` |
| `05-alphafold-cache/` | E | E | downloaded AlphaFold DB files |
| `06-params/<RES>/` | F build | F | shared ligand parameters (`meta.json`, `param.log`) |
| `06-systems/<ID>/` | F, G, H | G, H, status | `system.json`, `topol.top`, `index.ndx`, `mdp/`, `em.tpr`, `run.json`, `md.xtc`, `analysis/` |
| `06-systems/matrix.json` | F build | verify | which receptor x target pairs exist |
| `07-run-logs/`, `07-slurm/` | G | humans | background progress `*.jsonl`; `run_array.sbatch` |
| `08-analysis/analysis_summary.csv` | H | report | one row per analysed system |
| `09-reports/README.md`, `figures/` | `gromatrix-report` | humans | tables and figures for steps B-H |

## Command contract

Every `gmxprep-*`, `gromatrix-report` and `gromatrix-compute` command prints exactly one JSON
object on stdout, `{"ok": true|false, ...}`, and exits 0 or 1. Read that line; open logs only on
failure. System ids are `<enzyme>`, `<target>` and `<enzyme>__<target>`. The status command's
`next` is `{skill, agent?, command, why}`; `agent` is set for G and H.

Per-system `analysis/` also holds a `README.md` that explains each plot.
