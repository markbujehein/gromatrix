"""Names of every directory and file a batch directory holds: the one place to change them.

Top-level folders are numbered so a plain `ls` reads in pipeline order. Folders that belong to the
same stage share a number (05-alphafold and 05-alphafold-cache; 06-params and 06-systems;
07-run-logs and 07-slurm).

    00-config/          forcefield lock, pipeline manifests and result      (lock, pipeline)
    01-raw/             downloaded RCSB files                               (A fetch)
    02-fixed/           PDBFixer-cleaned structures                         (B fix)
    03-units/           chosen biological unit, sliced per enzyme / target  (C extract)
    04-relaxed/         OpenMM-relaxed copies + relax.json                  (D relax)
    05-alphafold/       AlphaFold comparison JSON, CSV and superposed models (E afcompare)
    05-alphafold-cache/ downloaded AlphaFold DB files
    06-params/          ligand parameters shared by all systems             (F build)
    06-systems/         one folder per GROMACS system; MD (G) and analysis (H) write here too
    07-run-logs/        progress logs of background runs                    (G run)
    07-slurm/           generated Slurm array script                        (G run --slurm)
    08-analysis/        matrix-wide summary table                           (H analyze)
    09-reports/         human-readable README and figures for every step    (gromatrix-report)

Per-system files (system.json, run.json, analysis/ ...) keep their names inside 06-systems/<ID>/.
"""

from __future__ import annotations

from pathlib import Path

LAYOUT_VERSION = 2  # 1 = the unnumbered names (raw/, fixed/, ...) of v0.1/v0.2-dev batches

CONFIG = "00-config"
RAW = "01-raw"
FIXED = "02-fixed"
UNITS = "03-units"
RELAXED = "04-relaxed"
ALPHAFOLD = "05-alphafold"
ALPHAFOLD_CACHE = "05-alphafold-cache"
PARAMS = "06-params"
SYSTEMS = "06-systems"
RUN_LOGS = "07-run-logs"
SLURM = "07-slurm"
ANALYSIS = "08-analysis"
REPORTS = "09-reports"

# Single files, relative to the batch directory.
LOCK = f"{CONFIG}/forcefield.lock.json"
COMPUTE = f"{CONFIG}/compute.json"
GENERATED_MANIFEST = f"{CONFIG}/manifest.generated.yaml"
PIPELINE_RESULT = f"{CONFIG}/pipeline.json"
AF_REPORT = f"{ALPHAFOLD}/af_report.csv"
MATRIX = f"{SYSTEMS}/matrix.json"
BATCH_RUN = f"{RUN_LOGS}/run.json"
ANALYSIS_SUMMARY = f"{ANALYSIS}/analysis_summary.csv"


def at(batch_dir: Path | str, *parts: str) -> Path:
    """`at(B, UNITS, "1ABC")` -> B/03-units/1ABC."""
    return Path(batch_dir).joinpath(*parts)


def system_dir(batch_dir: Path | str, system_id: str) -> Path:
    return at(batch_dir, SYSTEMS, system_id)
