# CLI reference: steps A–C

All commands take `--batch-dir`, print one JSON line on stdout and exit non-zero on failure.

## gmxprep-fetch (step A)

```
gmxprep-fetch --batch-dir B (--ids 1ABC ... | --query TEXT | --uniprot P12345 ...)
              [--max 25] [--method xray|em|nmr] [--max-resolution 3.0]
              [--format cif|pdb|both] [--out raw] [--check-updates]
```

| Flag | Meaning |
|---|---|
| `--uniprot` | every entry with a chain mapped to these accessions, best resolution first |
| `--query` | RCSB full-text search, by relevance |
| `--format` | `cif` (default; needed for assemblies in C), `pdb`, `both`. Large entries have no PDB file and report `unavailable`. |
| `--check-updates` | ask the RCSB Data API for each entry's revision and refetch revised entries |

Cache: `01-raw/cache.json` stores sha256/size/URL/revision. Files are reused only when present and
hash-identical (`cached`); truncated or edited files are refetched (`repaired`). Unknown or obsolete IDs
land in `failed` as `not found / obsolete`.

## gmxprep-fix (step B)

```
gmxprep-fix --batch-dir B [--in raw] [--out fixed] [--ids 1ABC ...] [--max-loop 8]
            [--waters none|ligand-shell|shell|all] [--water-shell 3.5]
            [--keep-het all|none|RES,RES] [--keep-additives] [--keep-nonstandard]
```

| Flag | Meaning |
|---|---|
| `--max-loop` | Longest *internal* gap PDBFixer may build. Longer gaps and all termini are left open and reported, and the PDB gets a TER there. |
| `--waters` | `none` (default); `ligand-shell` = waters within `--water-shell` Å of a ligand; `shell` = first shell of protein and ligands; `all`. Step F currently drops crystal waters regardless. |
| `--keep-het` | `all` (default, minus additives), `none`, or an allow-list. Listed names are kept even if they are additives. |
| `--keep-additives` | Keep GOL/EDO/SO4/PEG/... |
| `--keep-nonstandard` | Do not convert MSE→MET etc. (pdb2gmx will then fail) |
| `--altloc-rmsd` | Å (default 0.5). Symmetry-corrected RMSD at which alternate conformers count as distinct poses |

Always done: label→author chain-ID remap, missing heavy atoms (incl. OXT), removal of all H/D,
hydronium and UNX/UNL. Alternate conformers: the highest-occupancy conformer goes into the PDB (ties: first
letter). All conformers, with raw and symmetry-corrected RMSD and bond-geometry checks, go to
`meta.json → altlocs`. Coordinates are never averaged.
Outputs: `02-fixed/<ID>.pdb` (author numbering, original B-factors, built atoms have occupancy 0.00)
and `02-fixed/<ID>.meta.json` (entities, assemblies, operators, gap and removal report). JSON line
per structure: `built_atoms`, `modelled_gaps`, `open_gaps`, `removed_*`, `warnings`.

## gmxprep-extract (step C)

```
gmxprep-extract --batch-dir B [--in fixed] [--out units] [--ids ...]
                [--targets RES[,RES][@ID]] [--target-chains X[,Y][@ID]] [--expect RES[:n],...]
                [--assembly auto|none|N] [--ligand-cutoff 4.5] [--contact-cutoff 5.0]
                [--min-contacts 10] [--max-chains 12] [--unit CANDIDATE] [--allow-incomplete]
```

| Flag | Meaning |
|---|---|
| `--targets` | Ligand residue names that form the "target" of the matrix. `LIG@1ABC` scopes a name to one ID. Hard requirement for unit selection. |
| `--target-chains` | Author chains of a protein/peptide target (sliced like a ligand target) |
| `--expect` | Further hard requirements, e.g. `COF:1,HEM:2` (cofactors that must travel with the enzyme) |
| `--assembly` | `auto` = deposited assemblies + AU stoichiometric units; `none` = AU only (docked models); `N` = only assembly N |
| `--ligand-cutoff` (alias `--het-cutoff`) | Å; a HET is assigned to the chain with the most heavy-atom pairs within this distance |
| `--contact-cutoff`, `--min-contacts` | Chain–chain interface edge = ≥ `min-contacts` residue pairs within `contact-cutoff` Å |
| `--max-chains` | Discard candidates larger than this |
| `--unit` | Override the choice with a candidate id from `unit.json` (e.g. `assembly-1`, `au-ABEF`) |
| `--allow-incomplete` | Proceed although no unit holds every required HET (gives a warning) |
| `--altloc` | `auto` (default): branch a target only if it has a *distinct* conformer; `branch`: always branch multi-conformer targets; `primary`: never. Branches give `complex_<T>_alt<L>.pdb` / `target_<T>_alt<L>.pdb`, with coupled residues switched to the same letter, and are clash-checked. A pose closer than 1.5 Å to its environment is dropped. |

Ranking: complete → HET-type coverage → author-defined → fewer chains → larger interfaces →
fewer unmodelled residues → lower mean Cα B.
Outputs in `03-units/<ID>/`: `unit.pdb`, `enzyme.pdb` (polymer + cofactors, no targets),
`target_<RES>.pdb` (one copy, the best-contacting instance), `complex_<RES>.pdb` (enzyme + all
instances in the unit), `target_<chains>.pdb`/`complex_<chains>.pdb` for protein targets, and
`unit.json` (every candidate's metrics, every ligand assignment, warnings).

