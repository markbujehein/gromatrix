# gromatrix

[![CI](https://github.com/markbujehein/gromatrix/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/markbujehein/gromatrix/actions/workflows/ci.yml)
[![coverage](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/markbujehein/gromatrix/badges/coverage.json)](https://github.com/markbujehein/gromatrix/actions/workflows/ci.yml)
[![release](https://img.shields.io/github/v/release/markbujehein/gromatrix?sort=semver)](https://github.com/markbujehein/gromatrix/releases/latest)
[![python](https://img.shields.io/badge/python-3.10%20%7C%203.12%20%7C%203.13-blue)](pixi.toml)
[![pixi](https://img.shields.io/badge/env-pixi-yellow)](https://pixi.sh)
[![ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

> **Status: alpha (0.x).** The CLI, file formats and skills may still change between minor
> versions. Steps A–E are tested on an 82-structure benchmark. Steps F–H are covered by unit tests
> against simulated GROMACS output but have not yet run against real GROMACS; confirm on your own
> system before relying on them.

**gromatrix** is a general-purpose, target-agnostic pipeline that turns experimental or docked
structures into a GROMACS simulation **matrix**: every receptor alone, every target alone, and every
receptor + target complex. Targets can be small molecules, cofactors or other proteins. It works
from any PDB/mmCIF input: single chains, hetero-oligomers, multi-copy asymmetric units, and
symmetry-generated biological assemblies. One forcefield (AMBER or CHARMM family) is locked per batch
and enforced at every step, so all systems in a comparison are built with identical physics.

It covers the whole path from a PDB ID to MD results: preparation (steps A–F), running the MD
protocol with checkpoint resume locally or on Slurm (G), and trajectory analysis with a
comparison table across the matrix (H). Each step is one command; `gmxprep-pipeline` chains them
all from one manifest.

Nothing in the algorithms refers to a particular protein, chain, ligand or residue number. Every
decision is driven by the deposition metadata, geometry and the user's declared targets.

It ships with a **Claude agent toolkit**: six skills that let Claude drive the whole pipeline.
They install into any Claude client **without cloning this repository** (see below).

## Install the Claude skills (no clone needed)

Every [GitHub Release](https://github.com/markbujehein/gromatrix/releases/latest) carries one
ready-to-upload `.zip` per skill, plus `SHA256SUMS`:

| Download | Role |
|---|---|
| [`gromatrix.zip`](https://github.com/markbujehein/gromatrix/releases/latest/download/gromatrix.zip) | Orchestrator: decides the forcefield and targets, chains the stages, writes the final report. **Install this one first.** |
| [`gromatrix-structure-prep.zip`](https://github.com/markbujehein/gromatrix/releases/latest/download/gromatrix-structure-prep.zip) | Steps A–C: fetch, clean, sub-unit extraction |
| [`gromatrix-validate.zip`](https://github.com/markbujehein/gromatrix/releases/latest/download/gromatrix-validate.zip) | Steps D–E: OpenMM integrity relaxation, AlphaFold comparison |
| [`gromatrix-gmx-builder.zip`](https://github.com/markbujehein/gromatrix/releases/latest/download/gromatrix-gmx-builder.zip) | Step F: GROMACS topologies, boxes, ions, smoke test |
| [`gromatrix-execution.zip`](https://github.com/markbujehein/gromatrix/releases/latest/download/gromatrix-execution.zip) | Step G: EM → NVT → NPT → production; resume, gates, background or Slurm |
| [`gromatrix-analysis.zip`](https://github.com/markbujehein/gromatrix/releases/latest/download/gromatrix-analysis.zip) | Step H: PBC removal, RMSD, RMSF, Rg, H-bond persistence, plots, matrix table |

Install all six. The orchestrator hands off to the others by name.

**claude.ai (web, desktop, Pro / Team / Enterprise).** Go to *Customize → Skills → Upload a
skill* and upload each `.zip` as downloaded. Don't unzip them; each archive already contains
the `<skill-name>/SKILL.md` folder layout that claude.ai expects.

**Claude Code (CLI).** Unzip into your personal skills directory, which makes them available in
every project:

```bash
cd ~/.claude/skills   # create it first if needed: mkdir -p ~/.claude/skills
for s in gromatrix gromatrix-structure-prep gromatrix-validate gromatrix-gmx-builder \
         gromatrix-execution gromatrix-analysis; do
  curl -fsSLO "https://github.com/markbujehein/gromatrix/releases/latest/download/$s.zip"
  unzip -oq "$s.zip" && rm "$s.zip"
done
```

Optionally verify a download against the release's `SHA256SUMS` with `sha256sum -c SHA256SUMS`.
The archives are built reproducibly.

**The skills drive the `gromatrix` command-line tools, so install the package where the commands
will run.** That needs no clone either:

```bash
conda create -n gromatrix -c conda-forge python=3.12 openmm pdbfixer rdkit biopython scipy pyyaml
conda activate gromatrix
pip install "gromatrix @ git+https://github.com/markbujehein/gromatrix@v0.2.0"
gmxprep-status --batch-dir test_batch   # prints {"ok": false, ...}: the commands are installed
```

What runs where:
* **Steps A–E** need network access to RCSB / PDBe / AlphaFold DB.
* **Step F** additionally needs GROMACS, acpype or CGenFF (conda-forge: `gromacs acpype openbabel`).
* **Step G** needs GROMACS on the machine that runs MD (ideally with a GPU), or a Slurm cluster.
* **Step H** needs GROMACS plus `mdanalysis matplotlib-base` (conda-forge), or
  `pip install "gromatrix[analysis] @ git+…"`.

In practice, Claude Code on your workstation or cluster runs the whole pipeline. On claude.ai the
skills still guide planning, manifests and interpretation of results, and run commands wherever
code execution can install the package.

### Why six skills instead of one

A single large prompt for eight pipeline stages wastes context and blurs responsibilities. The
toolkit follows **progressive disclosure**:
* **Small entry point.** The orchestrator is a ~80-line entry point that makes the decisions that
  are hard to undo and routes to a stage.
* **Narrow stage skills.** Each stage skill documents only its own commands; a test enforces it.
  Detailed flags and troubleshooting live in its `references/` and are read only when needed.
* **Stateless hand-off.** The skills hand off through the batch directory, not through
  conversation memory. Every stage starts with `gmxprep-status`, which reports progress and names
  the `next` skill and command. A fresh session, or a different agent, resumes any batch at the
  right stage.

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
| Long runs die: walltime, node failure, Ctrl-C | Per-phase resume from `.cpt`; `-maxh` derived from the Slurm walltime; SIGTERM forwarded so mdrun checkpoints; stale downstream output never chained (G) |
| A blown-up equilibration silently feeds production | EM convergence and zero-LINCS gates on NVT/NPT before production starts (G) |
| PBC artefacts and group-number scripts in analysis | `trjconv -pbc mol -center` then `-fit`, groups selected by name; target atoms taken from the build record, not guessed (H) |
| Comparing many systems by eye | One `analysis_summary.csv` row per system, with drift/instability flags (H) |

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
| G | `gmxprep-run` | EM → NVT → NPT → production with resume and gates; foreground, `--background` or `--slurm` |
| H | `gmxprep-analyze` | PBC removal + fit, Cα/target RMSD, RMSF, Rg, H-bond persistence, plots, `analysis_summary.csv` |
| A–H | `gmxprep-pipeline` | runs the steps above in order from one manifest; `--dry-run`, `--skip-*`, `--from/--to`, fail-fast |

Every command takes `--batch-dir` (the pipeline takes `--manifest`), prints **one JSON line**
(`{"ok": ..., ...}`; `gmxprep-run` and `gmxprep-pipeline` stream progress lines before it) and exits
non-zero on failure. Verbose tool output goes to log files. `gmxprep-status --batch-dir B` summarises a batch
and names the next step. Full flags for each step are in the `references/cli.md` of the matching
agent skill (below).

## Agent skill reference

Inside a clone, Claude Code discovers the skills in `.claude/skills/` automatically. Outside one,
install the release zips (above).

| Skill | Scope | Commands it documents |
|---|---|---|
| `gromatrix` | Orchestrator: settles forcefield, targets, protonation and compute with the user; chains the stages or runs the pipeline; final report | `gmxprep-status`, `gmxprep-gmx lock`, `gmxprep-pipeline` |
| `gromatrix-structure-prep` | Steps A–C: fetch, fix, sub-unit extraction, altlocs | `gmxprep-fetch`, `gmxprep-fix`, `gmxprep-extract` |
| `gromatrix-validate` | Steps D–E: OpenMM integrity relaxation, AlphaFold comparison | `gmxprep-relax`, `gmxprep-afcompare` |
| `gromatrix-gmx-builder` | Step F: manifest, GROMACS build, verify, smoke test | `gmxprep-gmx` |
| `gromatrix-execution` | Step G: running, resuming and submitting the MD protocol | `gmxprep-run` |
| `gromatrix-analysis` | Step H: trajectory metrics, flags, plots, interpretation | `gmxprep-analyze` |

`gmxprep-status` reads only the artefacts the steps write:
* the forcefield lock;
* `fixed/*.meta.json` and `units/*/unit.json`;
* `relaxed/*/relax.json` and `af/*.json`;
* `systems/*/system.json`, `systems/*/run.json` and `systems/*/analysis/analysis.json`.

From these it reports a `steps` verdict for A–H (`done` / `partial` / `pending` / `n/a`) and the
`next` skill and command. That keeps the hand-off between skills deterministic.

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
├── pyproject.toml / pixi.toml / pixi.lock
├── examples/manifest.yaml            # matrix definition consumed by step F
├── examples/pipeline.yaml            # one-command A-H run (gmxprep-pipeline)
├── benchmarks/                       # validation reports on real structure panels
├── scripts/package_skills.py         # validates + zips each skill for upload (used by the release workflow)
├── src/gromatrix/
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
│   ├── gmx_prep.py      # [F] CLI: lock / plan / build / verify / smoke
│   ├── run.py           # [G] MD protocol: phase table, resume, gates, progress, background, Slurm array
│   ├── analyze.py       # [H] trjconv PBC/fit, MDAnalysis metrics, plots, matrix summary
│   ├── pipeline.py      # A-H runner: manifest -> plan -> in-process steps, generated step-F manifest
│   └── status.py        # gmxprep-status: A-H progress + next step (skill hand-off contract)
├── tests/                # unit tests + offline end-to-end regression (tests/data/)
└── .claude/skills/                   # agent skills: gromatrix (orchestrator), gromatrix-structure-prep,
                                      # gromatrix-validate, gromatrix-gmx-builder, gromatrix-execution,
                                      # gromatrix-analysis (each SKILL.md + references/)
```

## Batch directory (created by the pipeline)

```
<batch>/
├── forcefield.lock.json          # single source of truth for the forcefield
├── raw/ fixed/ units/ relaxed/   # outputs of A-D (one entry per structure ID)
├── af/ af_report.csv afdb/       # E: per-chain metrics, superposed models, cached AF models
├── params/<RES>/                 # one parameter set per ligand/cofactor, shared by all systems
├── systems/<ID>/                 # E1, T1, E1__T1 ...: topol.top, solv_ions.gro, index.ndx, em.tpr, run.sh
│   ├── run.json, run.log         # G: per-phase state; {em,nvt,npt,md}.{tpr,cpt,gro,log,edr}, md.xtc
│   └── analysis/                 # H: analysis.json, rmsd.png, rmsf.png, rg.png, hbonds.png
├── slurm/run_array.sbatch        # G --slurm
├── matrix.json                   # F
├── analysis_summary.csv          # H: one row per system
└── manifest.generated.yaml, pipeline.json   # gmxprep-pipeline
```

## Quick start

```bash
pixi shell -e md          # from a clone; or the conda + pip install shown above
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
gmxprep-run          --batch-dir $B --threads 8 --gpu-id 0      # or --background, or --slurm --time 24:00:00
gmxprep-analyze      --batch-dir $B                             # -> $B/analysis_summary.csv
gmxprep-status       --batch-dir $B                             # where am I, what's next
```

Or all of it from one manifest ([`examples/pipeline.yaml`](examples/pipeline.yaml)):

```bash
gmxprep-pipeline --manifest pipeline.yaml --dry-run            # the exact commands; nothing runs
gmxprep-pipeline --manifest pipeline.yaml --run-md --analyze   # A-H, stops at the first failed step
gmxprep-pipeline --manifest pipeline.yaml --from G --run-md    # resume after fixing a failure
```

Rerunning is always safe:
* steps A–F skip work whose inputs haven't changed;
* G resumes from checkpoints;
* H skips analyses that are up to date.

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
pixi shell -e md     # adds GROMACS / acpype / Open Babel for real step F-H runs
```

Run `pixi run lint && pixi run test` before pushing. CI (`.github/workflows/ci.yml`) runs the same
tasks:
* Ruff, plus the full test suite with coverage, on every pull request, using Python 3.12.
* Additionally Python 3.10 and 3.13 on pushes to `main`/`dev`.
* Environments are cached from the lockfile, so jobs install in seconds. After changing
  dependencies, run `pixi lock` and commit `pixi.lock`; CI uses `--locked` and fails if the lock
  is stale.
* The `skills` job validates and zips every agent skill on each push and pull request, so a broken
  skill (e.g. invalid YAML frontmatter) fails CI instead of failing on upload.

**Releasing.** Set the version in `pyproject.toml` / `pixi.toml` / `src/gromatrix/__init__.py`, merge
to `main`, then tag:

```bash
pixi run -e skills package-skills          # optional local check: dist/skills/*.zip + SHA256SUMS
git tag v0.2.0 && git push origin v0.2.0   # .github/workflows/release.yml publishes the zips
```

The release workflow refuses a tag that doesn't match the package version.

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
  * `-t` checkpoint chaining is built into `run.sh` and `gmxprep-run`.

  Avoided: blanket `-maxwarn 1`, genion group numbers, hand-adjusted residue charges.
* **Step G.** The protocol is a four-row phase table: each phase's input coordinates, whether it
  is restrained (`-r`) and whether it continues a checkpoint (`-t`).
  * **Resume.** Each phase is classified *done*, *partial* or *new* from the last segment of its
    log, so an appended log from a `-maxh` stop is never mistaken for completion. Once a phase
    reruns, every later phase starts fresh.
  * **Interrupts.** `SIGTERM` is forwarded to mdrun, so it writes a checkpoint. The Slurm array
    script sets `--signal=B:TERM@300` and `-maxh` just under the walltime, so resubmitting
    continues every system.
  * **Gates.** A failed gate renames the phase's `.gro`, so the phase is never resumed as
    complete.
  * **Ownership.** `run.json` records the host and pid, so two runners never work on one system.
* **Step H.**
  * **PBC.** A two-pass `trjconv`: `-pbc mol -center`, then `-fit rot+trans`, writing the solute
    only. `Solute` is a prefix of every system, so the step-F index groups stay valid for the
    solute files, and selections are built from those groups rather than residue-name heuristics.
  * **Topology.** MDAnalysis reads a `.gro` (`editconf` of the solute `.tpr`), so the analysis
    doesn't depend on the GROMACS `.tpr` version.
  * **Pose retention.** Target RMSD is computed in the protein frame, without refitting the
    target.
  * **H-bonds.** Geometric (3.0 Å / 150°), with donors paired to hydrogens by distance; each
    donor–acceptor pair gets an occupancy.
* **Pipeline.** Each step runs in-process through its own CLI `main()`, so the pipeline adds no
  logic of its own: the same locks, caches and gates apply. Every `steps.<name>` key in the
  manifest becomes that command's flag.
