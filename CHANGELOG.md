# Changelog

## v0.2.0 (alpha, unreleased)

### Added
- **Claude Code plugin.** The repository is now an installable plugin (`/plugin marketplace add
  markbujehein/gromatrix`): `.claude-plugin/`, `skills/`, `agents/`, `references/`, `examples/`.
  - Skills renamed without the prefix (`gromatrix:start`, `structure-prep`, `validate`, `gmx-builder`,
    `execution`, `analysis`) with descriptions in the form what / when / covers.
  - Subagents `gromatrix:md-runner` (G) and `gromatrix:md-analyst` (H): the orchestrator never runs
    `gromatrix-run` or `gromatrix-analyze` itself; the agents return one JSON payload.
  - `gromatrix-status` `next` now carries `agent` for G and H and plugin-namespaced skill names.
  - Replaces the standalone skill zips and `scripts/package_skills.py` with `scripts/check_plugin.py`.
- **Plugin hooks** (`hooks/hooks.json`, `scripts/hooks/`): subagent-only `gromatrix-run` / `-analyze`,
  a Bash guard (`-maxwarn`, `gmx mdrun`, `rm -r` of batch folders, `build --force` over MD results),
  a read guard (trajectories, energy files, checkpoints, whole logs), a SubagentStop payload check,
  `gromatrix-status` at SessionStart, and `00-config/events.jsonl`. A second layer; CLI gates stay
  authoritative. `scripts/check_plugin.py` validates `hooks.json`.
- **Codecov**: CI uploads `coverage.xml` (`codecov.yml`, non-blocking); the README badge points to it and
  replaces the `badges`-branch job.
- **`gromatrix-watch`**: one line per MD phase change per system; exits when all systems are finished.
- **`gromatrix-compute`**: stores where MD runs (local, Slurm, remote) in `00-config/compute.json`;
  CPU-only is fine; the orchestrator asks once before step G.
- **`gromatrix-report`**: `09-reports/README.md` with tables and figures for steps B-H, and a README in
  every system's `analysis/` folder explaining each plot.
- MIT `LICENSE` and `THIRD_PARTY_NOTICES.md`.
- **G `gromatrix-run`**: runs EM → NVT → NPT → production for every built system.
  - Per-phase resume from checkpoints (`mdrun -cpi`).
  - Gates: EM convergence, plus zero LINCS warnings in NVT/NPT.
  - Live JSON progress.
  - Runs in the foreground, with `--background`, or as a Slurm array (`--slurm`, with `-maxh`
    derived from the walltime and SIGTERM checkpointing).
- **H `gromatrix-analyze`**: two-pass `trjconv` PBC removal and fit.
  - Metrics: Cα RMSD, target RMSD in the protein frame, RMSF, Rg, and target–enzyme H-bond
    occupancy.
  - Flags: `ligand_drift`, `protein_unstable`, `no_target_hbonds`.
  - Plots and `analysis_summary.csv` across the matrix.
- **`gromatrix-pipeline`**: steps A–H from one manifest (`examples/pipeline.yaml`).
  - Options: `--dry-run`, `--skip-relax` / `--skip-afdb` / `--skip-smoke`, `--from` / `--to`,
    and fail-fast with `--keep-going`.
  - Without `enzymes:`, the step-F manifest is generated from the step-C units.
- `gromatrix-status` reports `steps` A–H and routes to the new skills. A batch built from the
  user's own PDBs no longer points back to `gromatrix-fetch`.
- Skills `gromatrix-execution` (G) and `gromatrix-analysis` (H). The orchestrator now covers
  end-to-end steering and one-command mode. A test enforces that stage skills only run their own
  commands.
- Step F records `target_resname`, `target_residues` and `protein_residues` in `system.json`.
- pixi `analysis` feature (MDAnalysis, matplotlib) in the test and md environments; pyproject
  extra `analysis`.

### Changed
- **Breaking: every `gmxprep-*` command is now `gromatrix-*`** (`gromatrix-fetch`, `-fix`, `-extract`, `-relax`,
  `-afcompare`, `-gmx`, `-run`, `-analyze`, `-pipeline`, `-status`), one namespace with `gromatrix-report` and
  `gromatrix-compute`. No aliases: the package is unreleased.
- **Breaking: numbered batch folders** (`00-config` ... `09-reports`), all defined in `gromatrix.layout`.
  Batches from earlier development builds must be recreated or their folders renamed.
- Step F runs `pdb2gmx -renum`, so insertion codes stay distinct residues. It records
  `protein_residue_labels` (author labels) for step H. A build-recipe version in `inputs_sha`
  rebuilds existing systems once.
- Step F refuses to rebuild a system that has MD results (`run.json`) without `--force`, because
  a rebuild deletes them.
- `checks.scan_md_log` reads the newest log segment for how a run ended, and counts LINCS
  warnings over all segments. It is shared by smoke and step G.
- The import package is now `gromatrix`, as the distribution always was; it was `gmxprep`. The
  `gmxprep-*` command names were unchanged at that point (all commands are `gromatrix-*` as of v0.2.0).

## v0.1.0 (alpha, never tagged)

First alpha: steps A–F, status and the agent toolkit.

### Pipeline (`gromatrix-*` commands)
- **A `gromatrix-fetch`**: RCSB download by ID, full-text or UniProt search. Downloads are cached and
  verified by sha256, corrupted files are repaired, and entry revisions are tracked
  (`--check-updates`).
- **B `gromatrix-fix`**: PDBFixer clean-up.
  - Mirrors OpenMM's label→author chain-ID rule.
  - Builds on a primary-conformer copy of the input.
  - Builds short loops; long gaps are split with TER.
  - Strips hydrogens, waters and additives by policy.
- **C `gromatrix-extract`**: completeness-first selection of the sub-unit across deposited
  assemblies and stoichiometric asymmetric-unit units.
  - Ligands follow the chain they physically bind.
  - Alternate conformers branch only when genuinely distinct.
  - Output: enzyme / target / complex files.
- **D `gromatrix-relax`**: two-stage restrained OpenMM minimisation with explosion, distortion,
  chain-break, clash and pocket-collapse gates.
- **E `gromatrix-afcompare`**: AlphaFold DB comparison per protein chain: RMSD (all residues and
  confident core), TM-score (agrees with TM-align within 0.003) and Cα-lDDT.
- **F `gromatrix-gmx`**: GROMACS matrix builder with a batch-wide forcefield lock (AMBER / CHARMM
  families).
  - One shared parameter set per ligand (acpype/GAFF2 or CGenFF), with conflict detection.
  - Atom-count and charge gates; a zero-warning `em.tpr`; named index groups; `run.sh`.
  - `verify` and `smoke` subcommands.
- **`gromatrix-status`**: batch progress plus the next step; the hand-off contract between agent
  skills.

### Claude agent toolkit
- Four skills: the `gromatrix` orchestrator plus `gromatrix-structure-prep`,
  `gromatrix-validate` and `gromatrix-gmx-builder`.
- Shipped as reproducible `.zip` release assets with `SHA256SUMS`, installable on claude.ai or in
  Claude Code without cloning.

### Engineering
- pixi environments with a committed lockfile.
- GitHub Actions CI: Ruff, pytest with coverage on Python 3.10, 3.12 and 3.13, and skill
  validation.
- SHA-pinned actions, least-privilege permissions, and a coverage badge with no third-party
  service.
- Tag-triggered release workflow.

### Validation
- Offline end-to-end regression (1XVP) runs in CI.
- 82-entry human PXR/CAR benchmark: [`benchmarks/pxr_car_validation.md`](benchmarks/pxr_car_validation.md).
  It is a partial run: steps B, C and E on 16 entries, with step D pending on local hardware.

### Known limitations
- Step F has not yet been run against real GROMACS / acpype / CGenFF installations in CI. Confirm
  locally with `pixi shell -e md`.
- Covalent ligands and metal-containing cofactors are refused, not parameterised.
- Crystal waters are not carried into step F.
