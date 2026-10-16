# CLI reference

All commands take `--batch-dir`. stdout is one JSON object. Exit code 0 means `"ok": true`.

## gmxprep-gmx (step F, implemented)

```
gmxprep-gmx lock   --batch-dir B --forcefield NAME
gmxprep-gmx plan   --batch-dir B --manifest M [--allow-incomplete]
gmxprep-gmx build  --batch-dir B --manifest M [--forcefield NAME] [--only ID ...] [--force]
                   [--box dodecahedron|cubic|octahedron|triclinic] [--pad 1.2] [--conc 0.15]
                   [--temp 300] [--md-ns 100] [--maxwarn 0] [--gmx gmx]
                   [--allow-incomplete] [--allow-auto-protonation] [--accept-high-penalty]
gmxprep-gmx verify --batch-dir B [--gmx gmx]
gmxprep-gmx smoke  --batch-dir B [--only ID ...] [--nvt-steps 5000] [--mdrun-args "-nt 8"]
```

| Flag | Meaning |
|---|---|
| `--forcefield` | `amber`→amber99sb-ildn, `amber99sb-ildn-tip4pew`, `amber14sb`, `amber19sb` (OPC), `charmm`→charmm36m. Must equal the lock. |
| `--only` | system ids such as `E1`, `T1`, `E1__T1` |
| `--force` | rebuild systems already stamped `ok` |
| `--pad` | minimum solute–box distance (nm). Keep it the same across a matrix. |
| `--conc` | NaCl (KCl not supported) in mol/L on top of neutralisation |
| `--maxwarn` | applies only to the final `em.tpr` gate. Leave it at 0. |
| `--gmx` | e.g. `gmx_mpi` |
| `--allow-auto-protonation` | Open Babel pH 7.4 guess for ligands without SMILES/SDF. Results include a warning. |
| `--accept-high-penalty` | allow CGenFF penalties > 50 |

Environment variables for the CHARMM route: `CGENFF_BIN` (path to `cgenff`), `CGENFF2GMX`
(path to `cgenff_charmm2gmx_py3_nx2.py`). Forcefield ports that are not bundled
(`amber14sb_OL15.ff`, `amber19sb.ff`, `charmm36-jul2022.ff`) go in the batch dir or on `GMXLIB`.

Build gates (any failure stops that system): topology↔coordinate atom walk after `solvate`, integer
solute charge, zero net charge after `genion`, and a zero-warning `grompp` of `em.tpr`. Systems are
skipped on rebuild only when their `inputs_sha` (source files, ligand specs, build flags) is unchanged.
Manifest complexes may carry `variant:` (e.g. `altA`), giving system ids `E__T__altA`.

`smoke` runs EM per system (pass = Fmax < emtol and negative Epot), then `--nvt-steps` of restrained
NVT (pass = finished, no fatal error, zero LINCS warnings), and records the result in `system.json`.
`verify` also re-runs the atom walk and charge check, and fails on a recorded smoke failure.

### Outputs per system `B/systems/<id>/`

`topol.top`, `protein.itp`, `<RES>.itp`, `ligand_params.itp`, `posre.itp`, `conf.gro`, `solv_ions.gro`,
`index.ndx` (System, Protein, Solute, Solvent, LIG_<RES>), `em.tpr`, `mdp/{ions,em,nvt,npt,md}.mdp`,
`run.sh` (em→nvt→npt→md with correct `-r`/`-t` chaining), `<ff>.ff` symlink for non-bundled ports,
`system.json` (stamp), `build.log`, `error.txt` (on failure).

Ligand cache `B/params/<RES>/`: `template.sdf` (atom order is the order used in the topology),
`template.mol2`, `<RES>.itp`, `params.json`, `meta.json`, `param.log`.

## gmxprep-fix (step B, implemented)

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
Outputs: `fixed/<ID>.pdb` (author numbering, original B-factors, built atoms have occupancy 0.00)
and `fixed/<ID>.meta.json` (entities, assemblies, operators, gap and removal report). JSON line
per structure: `built_atoms`, `modelled_gaps`, `open_gaps`, `removed_*`, `warnings`.

## gmxprep-extract (step C, implemented)

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
Outputs in `units/<ID>/`: `unit.pdb`, `enzyme.pdb` (polymer + cofactors, no targets),
`target_<RES>.pdb` (one copy, the best-contacting instance), `complex_<RES>.pdb` (enzyme + all
instances in the unit), `target_<chains>.pdb`/`complex_<chains>.pdb` for protein targets, and
`unit.json` (every candidate's metrics, every ligand assignment, warnings).

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

Cache: `raw/cache.json` stores sha256/size/URL/revision. Files are reused only when present and
hash-identical (`cached`); truncated or edited files are refetched (`repaired`). Unknown or obsolete IDs
land in `failed` as `not found / obsolete`.

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
