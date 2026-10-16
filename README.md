# gromatrix

**gromatrix** is a general-purpose, target-agnostic pipeline that turns experimental or docked
structures into a GROMACS simulation **matrix**: every receptor alone, every target alone, and every
receptor + target complex. Targets can be small molecules, cofactors or other proteins. It works
from any PDB/mmCIF input: single chains, hetero-oligomers, multi-copy asymmetric units, and
symmetry-generated biological assemblies. One forcefield (AMBER or CHARMM family) is locked per batch
and enforced at every step, so all systems in a comparison are built with identical physics.

Nothing in the algorithms refers to a particular protein, chain, ligand or residue number. Every
decision is driven by the deposition metadata, geometry and the user's declared targets.

## Problems it solves, for any structure

| Problem in real deposits | How gromatrix handles it (step) |
|---|---|
| The asymmetric unit is not the biological unit; copies differ in which ligands they hold | Enumerates every deposited assembly *and* every stoichiometric, interface-connected unit of the AU; ranks them by completeness against the declared targets (C) |
| Ligands/cofactors must stay with the subunit they bind | Distance-based HET→chain assignment (most heavy-atom contacts within a cutoff) (C) |
| mmCIF `label_asym_id` vs author chain IDs (OpenMM/PDBFixer use label IDs) | Every chain mapped back to author IDs before anything else (B) |
| Alternate conformers | Primary by occupancy; symmetry-aware RMSD decides whether conformers are distinct; distinct poses branch the matrix; coordinates are never averaged (B, C) |
| Missing atoms, short loops, long disordered gaps | Heavy atoms and loops ≤ `--max-loop` built; longer gaps left open and split with TER so no tool bonds across them (B) |
| Crystallisation additives, hydrogens, waters | Stripped by configurable policy; H rebuilt later by the forcefield's own tools (B, F) |
| Is the structure physically sane before MD? | Two-stage restrained OpenMM minimisation with explosion/distortion/chain-break/clash/pocket-collapse gates (D) |
| Is the experimental unit consistent with the predicted fold? | AlphaFold DB comparison by sequence alignment: RMSD (all/confident core), TM-score, Cα-lDDT (E) |
| Mixing forcefields across a comparison | Batch-wide forcefield lock checked by every command and stamped on every artefact (F) |
| Ligand topologies differing between apo/holo/complex systems | One parameter set per ligand per batch (acpype/GAFF2 or CGenFF); every pose mapped onto it (F) |

## Pipeline

| Step | Command | What it does |
|---|---|---|
| lock | `gmxprep-gmx lock` | fix the forcefield for the batch (`amber`, `amber99sb-ildn-tip4pew`, `amber14sb`, `amber19sb`, `charmm`) |
| A | `gmxprep-fetch` | RCSB download by ID, full-text or UniProt search; sha256-verified cache |
| B | `gmxprep-fix` | PDBFixer clean-up, author chain IDs, gap policy, altloc analysis |
| C | `gmxprep-extract` | choose the sub-unit and slice `enzyme` / `target_<T>` / `complex_<T>` |
| D | `gmxprep-relax` | OpenMM integrity relaxation with failure labels |
| E | `gmxprep-afcompare` | AlphaFold DB comparison per protein chain |
| F | `gmxprep-gmx plan / build / verify / smoke` | topologies, boxes, solvent, ions, `em.tpr`, `run.sh`, smoke test |

Every command takes `--batch-dir`, prints **one JSON line** (`{"ok": ..., ...}`) and exits non-zero on
failure. Verbose tool output goes to log files. Full flags: `.claude/skills/gromatrix/references/cli.md`.
An agent skill (`.claude/skills/gromatrix/SKILL.md`) routes an AI agent through the same commands.

## Validation

The algorithms are general; these datasets are used to **stress-test** them:

* **Regression fixture: 1XVP** (`tests/data/1XVP.cif.gz`). A multi-copy asymmetric unit in which
  one ligand is bound to only one of two copies, the ligand is modelled with two alternate
  conformers, and five deposited assemblies compete. It exercises completeness-first unit
  selection, label→author chain remapping and the altloc policy in a single small file. It runs
  offline in the test suite.
* **Benchmark panel: 82 human PXR/CAR crystal structures.** Every X-ray entry of human PXR
  (UniProt O75469) and/or human CAR (Q14994), run through steps A–E. Results, timings and the
  anomalies encountered are in [`benchmarks/pxr_car_validation.md`](benchmarks/pxr_car_validation.md).

## Repository layout

```
gromatrix/
├── pyproject.toml / environment.yml
├── examples/manifest.yaml            # matrix definition consumed by step F
├── benchmarks/                       # validation reports on real structure panels
├── src/gmxprep/
│   ├── net.py           # stdlib HTTP with retries + sha256 download cache
│   ├── fetch.py         # [A] RCSB fetch: IDs, full-text or UniProt search
│   ├── structure.py     # mmCIF/PDB metadata: entities, assemblies, operators, label->auth chains
│   ├── altlocs.py       # alternate-conformer policy (primary by occupancy, symmetry-aware RMSD, no averaging)
│   ├── fix.py           # [B] PDBFixer clean-up
│   ├── extract.py       # [C] assembly/AU candidates, ligand mapping, completeness scoring, matrix slicing
│   ├── relax.py         # [D] OpenMM two-stage restrained minimisation + integrity gates
│   ├── afdb.py          # [E] AFDB comparison: RMSD (all/core), TM-score, Cα-lDDT
│   ├── forcefields.py   # FF registry + batch lock
│   ├── ligands.py       # ligand router: acpype/GAFF2 (AMBER) | CGenFF (CHARMM), cache, pose mapping
│   ├── topology.py      # [ atomtypes ] split/merge, moleculetype normalisation, topol.top writer
│   ├── checks.py        # [F] safety checks: top/gro atom walk, charge, named index, EM/MD log parsing
│   ├── mdp.py           # em/nvt/npt/md mdp with FF-specific non-bonded settings
│   ├── gmx.py           # logged subprocess wrappers, .ff dir discovery, ions.itp parsing
│   ├── pdbio.py         # dependency-free PDB/GRO I/O, residue classification
│   └── gmx_prep.py      # [F] CLI: lock / plan / build / verify / smoke
├── tests/                # unit tests + offline end-to-end regression (tests/data/)
└── .claude/skills/gromatrix/         # agent skill: SKILL.md router + references/
```

## Batch directory (created by the pipeline)

```
<batch>/
├── forcefield.lock.json          # single source of truth for the forcefield
├── raw/ fixed/ units/ relaxed/   # outputs of A-D (one entry per structure ID)
├── af/ af_report.csv afdb/       # E: per-chain metrics, superposed models, cached AF models
├── params/<RES>/                 # one parameter set per ligand/cofactor, shared by all systems
├── systems/<ID>/                 # E1, T1, E1__T1 ...: topol.top, solv_ions.gro, index.ndx, em.tpr, run.sh
└── matrix.json
```

## Quick start

```bash
conda env create -f environment.yml && conda activate gmxprep
B=runs/my_study
gmxprep-gmx lock     --batch-dir $B --forcefield amber
gmxprep-fetch        --batch-dir $B --ids <PDB IDs>            # or --uniprot <ACC> --method xray
gmxprep-fix          --batch-dir $B
gmxprep-extract      --batch-dir $B --targets <LIGAND RESNAME> # e.g. the drug; cofactors via --expect
gmxprep-relax        --batch-dir $B
gmxprep-afcompare    --batch-dir $B
# write a manifest pairing units/<ID>/enzyme.pdb, target_<T>.pdb, complex_<T>.pdb (see examples/)
gmxprep-gmx build    --batch-dir $B --manifest manifest.yaml
gmxprep-gmx verify   --batch-dir $B && gmxprep-gmx smoke --batch-dir $B
```

## Development

The environment is managed with [pixi](https://pixi.sh). It reads `pixi.toml` and the committed
`pixi.lock`, so everyone (and CI) gets byte-identical conda-forge packages:

```bash
curl -fsSL https://pixi.sh/install.sh | bash   # or: brew install pixi / winget install prefix-dev.pixi
pixi run test        # fast unit tests (Python 3.12, OpenMM, PDBFixer, RDKit from conda-forge)
pixi run test-all    # full suite incl. the offline end-to-end regression, with coverage
pixi run lint        # ruff check + ruff format --check (tiny env, no science stack)
pixi run fmt         # auto-format and auto-fix before committing
pixi run -e py310 test-all   # the other Python edges CI checks on pushes (py310, py313)
pixi shell -e md     # adds GROMACS / acpype / Open Babel for real step-F builds
```

Run `pixi run lint && pixi run test` before pushing. CI (`.github/workflows/ci.yml`) runs the same
tasks:
* Ruff, plus the full test suite with coverage, on every pull request, using Python 3.12.
* Additionally Python 3.10 and 3.13 on pushes to `main`/`dev`.
* Environments are cached from the lockfile, so jobs install in seconds. After changing
  dependencies, run `pixi lock` and commit `pixi.lock`; CI uses `--locked` and fails if the lock
  is stale.

## Design notes

* **Step F in one paragraph.** `build` expands the manifest into systems and refuses an incomplete
  matrix. For each system it:
  1. runs `pdb2gmx -ignh -merge all` on the protein chains;
  2. parameterises every non-ion HET once per batch, with acpype/GAFF2/AM1-BCC for AMBER or
     CGenFF + `cgenff_charmm2gmx` for CHARMM (penalty-gated);
  3. maps each pose onto the cached template;
  4. maps crystal ions to the forcefield's `ions.itp`;
  5. merges ligand `[ atomtypes ]`;
  6. writes `topol.top`;
  7. runs `editconf`, `solvate` and `genion`, with atom-walk and charge gates before and after;
  8. writes a named `index.ndx`, a **zero-warning `grompp`** of `em.tpr`, and `run.sh`.
* **Patterns adopted from SQUIP** (ESS-MimerAI/SQUIP), each in a stricter form:
  * test-mode short runs became `smoke`;
  * log scanning (LINCS, fatal errors, "Finished mdrun") became `checks.scan_md_log`;
  * `[ molecules ]` vs `.gro` checks became a full atom walk with charge gates;
  * absolute includes became a relative forcefield symlink;
  * TaskBlaster caching became an input-hash resume;
  * `-t` checkpoint chaining is built into `run.sh`.

  Avoided: blanket `-maxwarn 1`, genion group numbers, hand-adjusted residue charges.
