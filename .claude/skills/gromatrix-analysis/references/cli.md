# CLI reference: step H

## gmxprep-analyze

```
gmxprep-analyze --batch-dir B [--only ID ...] [--deffnm md] [--gmx gmx] [--stride 1] [--begin PS] [--force]
                [--hbond-distance 3.0] [--hbond-angle 150] [--hbond-min-occupancy 0.1]
                [--ligand-rmsd-max 3.0] [--ca-rmsd-max 4.0]
```

| Flag | Meaning |
|---|---|
| `--stride` | keep every Nth frame (`trjconv -skip`) |
| `--begin` | ps; drop frames before this time (`trjconv -b`), e.g. to exclude relaxation |
| `--hbond-distance` / `--hbond-angle` | geometric H-bond criteria: donor–acceptor Å, donor–H–acceptor ° |
| `--hbond-min-occupancy` | pairs below this fraction of frames are counted, not listed |
| `--ligand-rmsd-max` | Å; `ligand_drift` flag threshold (last-half mean target RMSD) |
| `--ca-rmsd-max` | Å; `protein_unstable` flag threshold (last-half mean Cα RMSD) |

### GROMACS calls (in `systems/<id>/`, logged to `analysis/trjconv.log`)

```
gmx convert-tpr -s md.tpr -n index.ndx -o analysis/solute.tpr                           <<< Solute
gmx trjconv -s md.tpr -f md.xtc -n index.ndx -o analysis/whole.xtc -pbc mol -center     <<< Protein|LIG_<RES>, Solute
gmx trjconv -s analysis/solute.tpr -f analysis/whole.xtc -n analysis/solute.ndx \
            -o analysis/fit.xtc -fit rot+trans                                          <<< Backbone|LIG_<RES>, Solute
gmx editconf -f analysis/solute.tpr -o analysis/solute.gro
```

`Solute` is always a prefix of the system (step F writes it that way), so every solute group
keeps its atom numbers in `solute.ndx`.

### analysis.json

| Key | Content |
|---|---|
| `role`, `enzyme`, `target`, `selections` | system identity and the MDAnalysis selections used |
| `frames`, `time_ns` | frame count and times |
| `ca_rmsd`, `target_rmsd`, `rg` | per-frame series, Å |
| `rmsf` | `{resindex, resid, resname, rmsf, label}` per Cα; `label` is the author residue (`A:100A`) from step F |
| `hbonds` | `{count: per frame, pairs: [{donor, acceptor, occupancy, target_is_donor}], pairs_below_threshold}` |
| `summary.<series>` | `{mean, std, max, final, last_half_mean}` |
| `flags`, `plots`, `settings` | flags with numbers, png names, criteria used |

Which series exist depends on the system's role:

| Role | Cα RMSD / RMSF / Rg | target RMSD | H-bonds |
|---|---|---|---|
| enzyme | ✓ | – | – |
| target (ligand) | – | ✓ (fitted on itself) | – |
| target (protein) | ✓ | – | – |
| complex | ✓ (enzyme part) | ✓ (protein frame) | ✓ |

### Summary output

`{"ok","analysed","skipped":[{id,reason}],"failed":[{id,reason}],"flags":{id:[…]},"summary_csv","summary_rows"}`.
`analysis_summary.csv` holds one row per system with a successful `analysis.json`, across the whole
batch, not just this invocation. A stored analysis is reused only if its `settings` (deffnm, stride,
begin, H-bond criteria, flag thresholds) match the current flags.
`B/analysis_summary.csv` columns: `id, role, enzyme, target, status, frames, time_ns,
ca_rmsd_last_half, target_rmsd_last_half, rg_mean, hbond_mean, top_hbond, top_hbond_occupancy, flags`.
