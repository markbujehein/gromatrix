# gromatrix

Takes PDB files of enzymes docked with targets and turns them into a GROMACS simulation **matrix**
(each enzyme alone, each target alone, and every enzyme + target complex). One forcefield is locked
per batch and checked at every step.

## Repository layout

```
gromatrix/
├── pyproject.toml / environment.yml
├── examples/manifest.yaml            # matrix definition consumed by [F]
├── src/gmxprep/
│   ├── forcefields.py   # FF registry (amber99sb-ildn, amber14sb, amber19sb, charmm36m) + batch lock
│   ├── pdbio.py         # dependency-free PDB/GRO I/O, residue classification
│   ├── gmx.py           # logged subprocess wrappers, .ff dir discovery, ions.itp parsing
│   ├── topology.py      # [ atomtypes ] split/merge, moleculetype normalisation, topol.top writer
│   ├── ligands.py       # ligand router: acpype/GAFF2 (AMBER) | CGenFF (CHARMM), cache, pose mapping
│   ├── checks.py        # [F] safety checks: top/gro atom walk, charge, named index, EM/MD log parsing
│   ├── mdp.py           # em/nvt/npt/md mdp with FF-specific non-bonded settings
│   ├── gmx_prep.py      # [F] CLI: lock / plan / build / verify
│   ├── fetch.py         # [A] planned
│   ├── structure.py     # mmCIF/PDB metadata: entities, assemblies, operators, label->auth chains
│   ├── altlocs.py       # alternate-conformer policy (primary by occupancy, symmetry-aware RMSD, no averaging)
│   ├── fix.py           # [B] PDBFixer clean-up
│   ├── extract.py       # [C] assembly/AU candidates, ligand mapping, completeness scoring, matrix slicing
│   ├── relax.py         # [D] planned
│   └── afcompare.py     # [E] planned
├── tests/test_core.py, test_extract.py, data/1XVP.cif.gz (benchmark)
└── .claude/skills/gromatrix/
    ├── SKILL.md                      # operational router (what to run, in what order)
    └── references/                   # loaded only when needed
        ├── cli.md  forcefields.md  troubleshooting.md
```

## Batch directory (created by the pipeline)

```
<batch>/
├── forcefield.lock.json          # single source of truth for the FF
├── raw/ fixed/ units/ relaxed/   # outputs of A-D
├── af_report.csv                 # E
├── params/<RES>/                 # one parameter set per ligand/cofactor, shared by all systems
├── systems/<ID>/                 # E1, T1, E1__T1 ... each with topol.top, solv_ions.gro, em.tpr, mdp/
└── matrix.json
```

## CLI contracts

Every command takes `--batch-dir`, prints **one JSON line** to stdout (`{"ok": ..., ...}`) and exits
non-zero on failure; verbose tool output goes to log files. Steps C-F read the forcefield from the lock.

| Step | Command | Key arguments |
|---|---|---|
| lock | `gmxprep-gmx lock` | `--forcefield {amber,amber99sb-ildn,amber14sb,amber19sb,charmm,charmm36m}` |
| A | `gmxprep-fetch` | `--ids 1TQN 3QM4` or `--query TEXT [--max N]`, `--format {cif,pdb}`, `--assembly 1` |
| B | `gmxprep-fix` | `--in raw --out fixed`, `--max-loop 8`, `--waters {none,ligand-shell,shell,all}`, `--keep-het {all,none,LIST}` |
| C | `gmxprep-extract` | `--targets CID`, `--expect F15:1`, `--assembly {auto,none,N}`, `--ligand-cutoff 4.5`, `--min-contacts 10` |
| D | `gmxprep-relax` | `--in units --out relaxed`, `--restraint-k 10`, `--max-rmsd 0.5`, `--platform CUDA` |
| E | `gmxprep-afcompare` | `--in relaxed --out af_report.csv`, `--uniprot-map map.tsv`, `--plddt-min 70` |
| F | `gmxprep-gmx plan / build / verify / smoke` | `--manifest`, `--only`, `--box`, `--pad`, `--conc`, `--temp`, `--md-ns`, `--gmx` |

Design notes for the planned modules:

- **A** uses `search.rcsb.org/rcsbsearch/v2/query` for text search and `files.rcsb.org/download/{id}.cif`;
  mmCIF is the default because the biological assembly records are needed by C.
- **B** and **C** are implemented; see the docstrings of `fix.py` and `extract.py` and the
  1XVP section below.
- **D** uses OpenMM with the family-matched FF (`amber14-all.xml` or `charmm36.xml`, OBC2/GBn2 implicit
  solvent), harmonic heavy-atom restraints, ligands from `openmmforcefields` (GAFF2 or restrained), and
  fails a structure whose heavy-atom RMSD exceeds `--max-rmsd`.
- **E** maps chains to UniProt via SIFTS, fetches `alphafold.ebi.ac.uk/api/prediction/{acc}`, and reports
  Cα RMSD after superposition on aligned residues, TM-score (`tmtools`) and lDDT, with low-pLDDT regions excluded.

## [C] and the 1XVP benchmark

1XVP (CAR/RXR heterodimer with SRC1 peptides) has two copies in the asymmetric unit, but CITCO
(`CID`) is bound to CAR chain D only. Author assembly 1 (A,B,E,F) therefore has no CITCO.
`gmxprep-extract --targets CID` builds every assembly and every stoichiometric AU unit, maps each
HET to the chain it contacts most, and ranks completeness first. It selects assembly 2
(C,D,G,H; F15 on C, CID on D), and gets the same answer (`au-CDGH`) with assembly records
ignored. The deposit also models CITCO as two alternate conformers. Step B keeps A and reports it.

### Alternate conformers (1XVP CITCO)

1XVP has no conformer C, only A and B (27 atoms each, occupancy 0.50, on CID D1003). Their raw RMSD
is 1.03 Å, but almost all of it is a 180° flip of the para-chlorophenyl ring, which only relabels
C7↔C8 and C12↔C13. The symmetry-corrected RMSD is 0.38 Å, within the coordinate error at 2.6 Å
resolution. Averaging A and B would shrink ring C–C bonds from 1.44 to 0.70 Å. Policy: highest
occupancy → primary; every conformer analysed and kept; the matrix branches (`complex_<T>_alt<L>`)
only for distinct poses (`--altloc auto`) or on request (`--altloc branch`); never average.

## Patterns adopted from SQUIP (ESS-MimerAI/SQUIP)

| SQUIP pattern | Here |
|---|---|
| `test_mode` short runs before production | `gmxprep-gmx smoke`: EM + restrained NVT with logs parsed to pass/fail |
| `check_warnings.sh` (LINCS, fatal, "Finished mdrun") | `checks.scan_md_log`, zero LINCS tolerance in the smoke test |
| `verify_systems.py` (`[ molecules ]` vs `.gro`, water model) | `checks.top_vs_gro`: full atom walk + integer/zero charge gates in build and verify |
| Absolute `#include` rewriting for portability | relative `<ff>.ff` symlink in each system dir |
| `-water select` + stdin for TIP4P-Ew | `pdb2gmx -water none`, water from our own `topol.top`; `amber99sb-ildn-tip4pew` registry entry |
| TaskBlaster per-task caching | content-hash `inputs_sha` resume |
| `grompp -t prev.cpt` chaining | generated `run.sh` with `-r`/`-t` |
| *Avoided:* blanket `-maxwarn 1`, `echo 15 \| genion`, per-stage hand-written mdp drift, Parrinello-Rahman straight after NVT, hand-adjusted residue charges | `-maxwarn 0` (ions stage tolerates 1 only if charged), groups selected by name, single mdp generator, C-rescale, integer-charge gate |

## [F] in one paragraph

`build` expands the manifest into systems, refuses an incomplete matrix, and for each system: runs
`pdb2gmx -ignh -merge all` on the protein chains; parameterises every non-ion HET residue **once per batch**
(acpype/GAFF2/AM1-BCC for AMBER, CGenFF + `cgenff_charmm2gmx` for CHARMM, with penalty gating); maps
each pose onto the cached template so every system shares identical ligand topologies; maps crystal ions
to the FF's `ions.itp` moleculetypes; merges ligand `[ atomtypes ]` ahead of all moleculetypes; writes
`topol.top`; then `editconf`, `solvate`, `genion` and a **zero-warning `grompp`** of `em.tpr` as the
integrity gate. `verify` re-checks the lock, every parameter cache, every `topol.top` include and every
mdp's non-bonded block.

## Quick start

```bash
conda env create -f environment.yml && conda activate gmxprep
gmxprep-gmx lock  --batch-dir runs/cyp --forcefield amber
gmxprep-gmx plan  --batch-dir runs/cyp --manifest examples/manifest.yaml
gmxprep-gmx build --batch-dir runs/cyp --manifest examples/manifest.yaml
gmxprep-gmx verify --batch-dir runs/cyp
```
