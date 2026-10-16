# CLI reference: steps D–E

All commands take `--batch-dir`, print one JSON line on stdout and exit non-zero on failure.

## gmxprep-relax (step D)

```
gmxprep-relax --batch-dir B [--in units] [--out relaxed] [--ids ...]
              [--files enzyme.pdb complex_*.pdb target_*.pdb] [--forcefield amber|charmm]
              [--solvent auto|gbn2|rf] [--restraint-k 1000] [--max-iter 2000] [--free-iter 1000]
              [--max-rmsd 2.5] [--max-drop 100] [--pocket-cutoff 5.0] [--ph 7.0] [--platform auto|CPU|CUDA|OpenCL]
```

Integrity check, not equilibration. Polymer chains are protonated at `--ph` with the lock's family
(amber14-all / charmm36). HETs and ions are held out, because they are parameterised once in F,
and copied back unchanged. Stage 1: all heavy atoms restrained. Stage 2: only residues within
`--pocket-cutoff` of a HET restrained, so empty ligand sites cannot collapse.
`--solvent auto` = GBn2 on a GPU, otherwise `rf` (1.6 nm cutoff with reaction field, ε 78.3).
GBn2 on CPU costs about 0.7 s per evaluation for 8k atoms, versus 25 ms for rf.

Failure labels in the JSON: `explosion` (non-finite), `implosion/explosion` (RMSD > `--max-rmsd` **and**
energy drop below −`--max-drop` kJ/mol/atom), `distortion` (RMSD alone), `chain break` (Cα–Cα > 4.2 Å
created), `unresolved clash` (final E > 0), `pocket collapse` (< 2.0 Å to a held-out HET). Ligand-only
files are `skipped`. Per-file details: `relaxed/<ID>/relax.json`.

## gmxprep-afcompare (step E, module `afdb.py`)

```
gmxprep-afcompare --batch-dir B [--in relaxed] [--ids ...] [--files enzyme.pdb target_*.pdb]
                  [--uniprot A:P12345,B:Q67890[@ID]] [--plddt-min 70] [--min-length 30]
                  [--tm-min 0.8] [--lddt-min 0.7] [--out af_report.csv] [--afdb-dir afdb]
```

Chain → UniProt via SIFTS (author chain; assembly copies traced through `units/<ID>/unit.json`) or
`--uniprot`. The model comes from the AFDB API, with its version recorded and files hash-cached in
`afdb/`. Residues are matched by sequence alignment, not numbering. Reported per chain: `rmsd_all`,
`rmsd_core` (pLDDT ≥ `--plddt-min`), `tm_score` (normalised by experimental length; within 0.003 of
TM-align in validation, see `benchmarks/`), `lddt_ca` (local, superposition-free), identity, coverage, mean pLDDT. Flags:
TM < `--tm-min`, lDDT < `--lddt-min`, identity < 0.9. Chains shorter than `--min-length` are listed
but not scored. Outputs: `af_report.csv`, `af/<ID>.json`, superposed models `af/<ID>/AF_<chain>_<acc>.pdb`.

