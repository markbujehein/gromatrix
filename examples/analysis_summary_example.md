# Example final report (format target)

> **SYNTHETIC EXAMPLE.** Every name and number below is invented to show the structure and tone of a
> gromatrix final report. Never copy them into a real report. Fill each field from the batch's own
> files: `00-config/forcefield.lock.json`, `matrix.json`, `analysis_summary.csv`, each system's
> `run.json` and `analysis.json`.

# gromatrix report: `demo-batch`

## Summary

Four systems were built and simulated (receptor ENZ1, target LIGX, their complex, and the apo
receptor) for 20 ns each. The complex kept LIGX in the pocket (target RMSD 1.4 Å, last half) and the
protein stayed folded (Cα RMSD 1.9 Å). One flag was raised. This is a single replica per system.

## Setup

| | |
|---|---|
| Forcefield | `amber99sb-ildn` (GAFF2 / AM1-BCC, TIP3P), locked in `00-config/` |
| Structures | ENZ1: PDB `1ABC`, assembly 2 (chains A, B), chosen because it holds both ligands |
| Compute | local, 12 CPU threads, no GPU |
| Sampling | **one replica of 20 ns per system** (production), after 100 ps NVT and 100 ps NPT |
| Progress | built: 4/4, MD done: 4/4, analysed: 4/4 |

## Matrix

Rows are receptors, columns are targets. Each cell: last-half target RMSD (Å), mean H-bonds to the
enzyme, and the top persistent H-bond with its occupancy.

| receptor \ target | LIGX |
|---|---|
| ENZ1 | 1.4 Å · 0.9 H-bonds · `LIGX:N2` → `ASN165:OD1` (41%) |

Alone systems: ENZ1 apo, Cα RMSD 2.2 Å · LIGX in solution, RMSD 1.8 Å.

## Per-system results

| system | role | ns | Cα RMSD (Å) | target RMSD (Å) | Rg (Å) | H-bonds | flags |
|---|---|---|---|---|---|---|---|
| ENZ1 | enzyme | 20.0 | 2.2 | - | 23.6 | - | none |
| LIGX | target | 20.0 | - | 1.8 | - | - | none |
| ENZ1__LIGX | complex | 20.0 | 1.9 | 1.4 | 23.5 | 0.9 | none |
| ENZ1__LIGX_altB | complex | 20.0 | 3.1 | 3.4 | 23.9 | 0.2 | `ligand_drift` |

## Findings and warnings (verbatim, by step)

- **C unit choice:** `1ABC: selected assembly-2 (chains A, B); runner-up assembly-1 held one ligand only.`
- **C altlocs:** `LIGX: conformers distinct (RMSD 1.9 Å); matrix branched into altA and altB.`
- **D relax:** no failures.
- **E AlphaFold:** no flags (chain A: TM 0.91, lDDT 0.86).
- **F build:** `LIGX: protonation taken from the supplied SMILES (charge -1).`
- **G run:** no failures; production LINCS warnings: 0 in every system.
- **H flags:** `ENZ1__LIGX_altB: ligand_drift: target RMSD 3.4 Å over the last half (> 3.0 Å) - pose not retained`

## Reading of the results

- The primary pose (ENZ1__LIGX) stayed in the pocket and kept one persistent hydrogen bond. The
  alternative conformer (altB) left its starting pose; with a single 20 ns replica this suggests
  altB is less stable, not that it is wrong.
- LIGX makes few hydrogen bonds; contacts are mostly hydrophobic, which is expected for this ligand class.
- **Limits:** 20 ns and one replica do not support a binding-affinity or ranking claim. Cα RMSD
  plateaus were reached after about 5 ns in every system.

## What is left

`gmxprep-status` reports the batch complete. Suggested next steps: repeat with three replicas, extend
the complex to 100 ns, and compare the two conformers with free-energy methods if a ranking is needed.

Figures and per-step tables: `09-reports/README.md`. Per-system plots and their guides:
`06-systems/<id>/analysis/README.md`.
