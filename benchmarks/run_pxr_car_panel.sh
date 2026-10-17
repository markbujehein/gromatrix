#!/usr/bin/env bash
# Reproduce the PXR/CAR validation panel (all entries, steps A-E).
# Usage: pixi run bash benchmarks/run_pxr_car_panel.sh [batch_dir]
# Expect roughly: A 2-5 min (network), B ~45 min, C <1 min, D ~40 min CPU (minutes on a GPU), E ~3 min.
set -euo pipefail
B=${1:-runs/pxr_car_panel}
mkdir -p "$B"
step() { local name=$1; shift; local s=$SECONDS; "$@"; echo "[timing] $name: $((SECONDS - s)) s" | tee -a "$B/timings.txt" >&2; }
step A gromatrix-fetch --batch-dir "$B" --uniprot O75469 Q14994 --method xray --max 500 > "$B/fetch.json" || true
step B gromatrix-fix --batch-dir "$B" > "$B/fix.json" || true
step C python benchmarks/pick_targets.py "$B"      # picks each entry's largest drug-like HET, runs step C
step D gromatrix-relax --batch-dir "$B" --files enzyme.pdb "complex_*.pdb" > "$B/relax.json" || true
step E gromatrix-afcompare --batch-dir "$B" > "$B/af.json" || true
echo "results in $B (fix.json, extract.json, relax.json, 05-alphafold/af_report.csv, timings.txt)"
