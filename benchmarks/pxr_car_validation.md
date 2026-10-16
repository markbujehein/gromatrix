# Validation: human PXR / CAR crystal-structure panel

gromatrix is target-agnostic. This panel is a **stress test**, not a use case: nuclear-receptor
ligand-binding-domain structures are a demanding real-world mix. They contain multi-copy asymmetric
units, heterodimers with coactivator peptides, fusion constructs, engineered binders, long
disordered inserts, ligands with alternate conformers, and inconsistent deposition metadata.

**Status: partial run.** All 82 entries were identified and downloaded. Steps B, C and E ran on
the first 16 entries (by ID). Step D was not run on the panel. The full run (all entries, all
steps) is deferred to local hardware with [`run_pxr_car_panel.sh`](run_pxr_car_panel.sh). This
report will be updated with those results.

## Method

| | |
|---|---|
| Panel definition | `gmxprep-fetch --uniprot O75469 Q14994 --method xray`: every X-ray entry containing human PXR (NR1I2, O75469) and/or human CAR (NR1I3, Q14994). **82 entries.** |
| Target per entry | Largest drug-like HET (by heavy atoms) after excluding water, ions and crystallisation additives ([`pick_targets.py`](pick_targets.py)). Entries with none are processed as apo (`enzyme.pdb` only). This heuristic exists only for the benchmark; real studies declare `--targets`. |
| Settings | All defaults: `--max-loop 8`, `--ligand-cutoff 4.5`, `--altloc auto`, `--altloc-rmsd 0.5`, `--plddt-min 70`, `--min-length 30`. |
| Hardware | Cloud container, 4 CPU cores, no GPU. |
| Outputs | Per-entry results in [`pxr_car_results/units.csv`](pxr_car_results/units.csv) and [`pxr_car_results/af_report.csv`](pxr_car_results/af_report.csv). |

## Success rates

| Step | Ran on | Succeeded | Failed | Notes |
|---|---|---|---|---|
| A fetch | 82 | 82 | 0 | mmCIF for all entries. The cache was verified separately: a rerun is 100% `cached`, and a corrupted file is detected by hash and `repaired`. |
| B fix | 16 | 16 | 0 | Run stopped after 16/82 on request (slow step); see timings. |
| C extract | 16 | 16 | 0 | Every entry resolved to an author-defined assembly. Of the 16, 12 have a target ligand and 4 are apo. |
| D relax | 0 | — | — | Not run on the panel. The offline 1XVP regression in CI covers it: RMSD 1.06–1.12 Å, ΔE −4.7 kJ/mol/atom, pocket clearance 2.5 Å, no failure labels. |
| E AlphaFold | 16 entries / 32 chains | 25 chains scored | 0 | 6 coactivator peptides (< 30 residues) skipped by design; 1 engineered binder has no UniProt entry. |

The 16 processed entries are 14 PXR and 2 CAR/RXRα structures. Step C picked units of 1 chain
(8 entries), 2 chains (4) and 4 chains (4). 11/16 entries had more than one candidate unit to
choose from.

## Execution time (4-core CPU)

| Step | Per entry | Basis |
|---|---|---|
| B fix | median **13 s**, mean 34 s, max 133 s (4S0S) | wall time between consecutive outputs, 15 intervals; dominated by PDBFixer loop building |
| C extract | **0.27 s** | 4.3 s for 16 entries |
| D relax | ~10–13 s per file (enzyme / complex) | 1XVP regression, reaction-field solvent; GBn2 is ~25× slower on CPU and is the default only on GPUs |
| E AlphaFold | **1.8 s** (1.2 s with a warm cache) | 28.6 s for 16 entries including SIFTS/AFDB downloads; 19 s cached |

## AlphaFold comparison (25 chains)

| Metric | Mean | Median | Range |
|---|---|---|---|
| TM-score | **0.949** | 0.964 | 0.896 – 0.991 |
| Cα-lDDT | **0.938** | 0.950 | 0.874 – 0.977 |
| RMSD, confident core (pLDDT ≥ 70) | **0.99 Å** | 0.92 Å | 0.58 – 1.43 Å |
| RMSD, all aligned residues | 1.80 Å | 1.42 Å | 0.59 – 3.85 Å |
| Sequence identity | 0.996 | 1.000 | ≥ 0.982 |
| Coverage of experimental chain | 0.983 | 1.000 | ≥ 0.943 |

| UniProt | Chains | TM-score | lDDT | Core RMSD |
|---|---|---|---|---|
| PXR (O75469) | 17 | 0.956 | 0.954 | 1.07 Å |
| CAR (Q14994) | 2 | 0.972 | 0.915 | 0.84 Å |
| RXRα (P19793, heterodimer partner) | 6 | 0.921 | 0.899 | 0.81 Å |

No chain was flagged (thresholds: TM < 0.8, lDDT < 0.7, identity < 0.9). The gap between core
and all-residue RMSD is expected: it comes from flexible termini and the PXR insert, where AF
pLDDT is low. The TM-score implementation was cross-checked against TM-align (`tmtools`) and
agrees within 0.003.

## Edge cases handled

The first two were **defects found by this benchmark**. Both were fixed generically, with
regression tests, before this report was written.

1. **Ligand labelled with the wrong chain (1XV9).** CI2 is deposited as chain D (`D1001`) but
   sits in chain B's pocket: 121 contacts with B, none with D.
   - *Before:* HET assembly membership trusted the deposited label, so no assembly was complete
     and the whole 8-chain asymmetric unit was selected.
   - *Fix:* membership now follows the chain the HET physically binds, and relabelled HETs
     produce a warning.
   - *After:* author assembly 1 (A, B, E, F), complete.
2. **A second ligand copy stored as an "alternate conformer" (1XV9).** CI2 altloc B is 48 Å from
   altloc A, in the other heterodimer's pocket. It is a copy, not a conformer.
   - *Before:* the clash check alone accepted a branch that floated the ligand outside the unit.
   - *Fix:* a branch must contact the unit. altB is now rejected with an explanation and only
     altA is written.
3. **Three genuine binding orientations (1ILH, SR12813).** Occupancies 0.4/0.4/0.2,
   symmetry-corrected RMSD 4.3 and 3.8 Å from the primary. The matrix branches into
   `complex_SRL_altA/B/C`, all clash-free (closest contacts 2.4–2.9 Å).
4. **Pseudo-conformers from a ring flip (1XVP, CITCO).** Raw RMSD 1.03 Å, symmetry-corrected
   0.38 Å, within coordinate error. The primary conformer is used with no branch. Averaging
   would have collapsed ring bonds to 0.70 Å; the pipeline never averages.
5. **Ligand present in only one of two copies (1XVP).** Completeness-first ranking selects
   assembly 2 (C, D, G, H), the only one holding both ligands. It does so with and without
   assembly records.
6. **Fusion constructs (4J5X).** Chains C/D are RXRα fused to an SRC-1 motif. SIFTS lists both
   accessions; the first-listed one (the 23-residue motif) gave TM 0.09. *Fix:* choose the
   accession covering most of the chain. Now TM 0.90 against RXRα.
7. **No SIFTS record at PDBe (4S0S).** Falls back to the RCSB Data API mapping (TM 0.99). The
   co-crystallised engineered Adnectin binder correctly reports "no UniProt mapping" instead of
   failing.
8. **mmCIF label vs author chain IDs.** 15/16 entries needed the label→author remap. Without it,
   OpenMM/PDBFixer would place every ligand on its own label chain.
9. **Long disordered inserts.** The PXR β-insert (residues ~177–210, 13–32 missing) was left
   open with TER records in all 14 PXR entries (35 open gaps), so no tool bonds across them. Short
   loops of 2–8 residues were modelled (15 gaps), and 58 terminal gaps were left untouched.
   1,241 missing heavy atoms were built (occupancy 0 marks them).
10. **Coactivator peptides and partners.** Six 10–15-residue SRC-1/coactivator chains travel
    with their receptor in the selected unit but are excluded from fold metrics by the length gate.

Not yet exercised on this panel: Step D failure labels (explosion, distortion, pocket collapse).
These need the local Step D run.

## Reproduce / complete

```bash
pixi run bash benchmarks/run_pxr_car_panel.sh runs/pxr_car_panel   # all entries, steps A-E, writes timings.txt
```
