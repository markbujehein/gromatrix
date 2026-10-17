# gromatrix architecture

This is the engineering reference behind the [README](../README.md). It covers the problems each
step solves, the code layout, the batch directory that carries all state, and the design decisions
per step. The CLI flags for each step are in the matching skill's `references/cli.md`
(`skills/<skill>/references/`); shared references (batch layout, biophysical rules, compute modes,
subagent contract) are in the top-level `references/`.

## Contracts every step follows

- **One command per step.** Each command takes `--batch-dir`, prints exactly one JSON result line
  (`{"ok": …}`) and exits non-zero on failure. `gmxprep-run` and `gmxprep-pipeline` stream
  progress lines before it. Exceptions are caught in `main()` and reported as `{"ok": false,
  "error": …}`, never as a traceback.
- **State lives in files.** Every step reads and writes only the batch directory.
  `gmxprep-status` reconstructs progress from those files alone, so any session or agent can
  resume a batch.
- **Content-keyed caching.** Downloads are sha256-verified. Ligand parameters are fingerprinted
  from their SMILES, charge and files. A built system's `inputs_sha` covers every input file and
  build flag. Reruns skip unchanged work.
- **Gates, not warnings.** A failed physical or bookkeeping check stops that structure or system
  with a reason: forcefield mismatch, non-integer charge, parameter conflict, a relaxation that
  blew up, an EM that did not converge, LINCS warnings during equilibration. No flag silently
  lowers a gate. `-maxwarn` is never raised to get past a warning. The one exception is the
  ions-stage `grompp` of a charged system: its single warning is the net charge that `genion` is
  about to neutralise.

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
| Comparing many systems by eye | One `08-analysis/analysis_summary.csv` row per system, with drift/instability flags (H) |

## Repository layout

```
gromatrix/
├── pyproject.toml / pixi.toml / pixi.lock
├── examples/manifest.yaml            # matrix definition consumed by step F
├── examples/pipeline.yaml            # one-command A-H run (gmxprep-pipeline)
├── benchmarks/                       # validation reports on real structure panels
├── scripts/check_plugin.py           # validates plugin.json, marketplace.json, skills/ and agents/ (CI)
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
│   ├── layout.py        # every batch folder and file name (numbered 00-config ... 09-reports)
│   ├── report.py        # gromatrix-report: Markdown tables + figures for B-H; per-system plot README
│   ├── plotstyle.py     # shared matplotlib style
│   ├── compute.py       # gromatrix-compute: local / slurm / remote choice for step G
│   ├── pipeline.py      # A-H runner: manifest -> plan -> in-process steps, generated step-F manifest
│   └── status.py        # gmxprep-status: A-H progress + next step (skill hand-off contract)
├── tests/                # unit tests + offline end-to-end regression (tests/data/)
├── .claude-plugin/       # plugin.json + marketplace.json: the repository is a Claude Code plugin
├── skills/               # start (orchestrator), structure-prep, validate, gmx-builder, execution, analysis
├── agents/               # md-runner (step G), md-analyst (step H): subagents returning one JSON payload
├── references/           # batch layout, biophysical rules, compute modes, subagent contract
└── examples/             # manifests and the final-report format example
```

## Batch directory (created by the pipeline)

```
<batch>/                            # numbered in pipeline order; related folders share a number
├── 00-config/                    # forcefield.lock.json (the forcefield), compute.json, pipeline files
├── 01-raw/ 02-fixed/             # A fetch, B fix (one entry per structure ID)
├── 03-units/ 04-relaxed/         # C extract, D relax
├── 05-alphafold/                 # E: <ID>.json, af_report.csv, superposed models (05-alphafold-cache/: AF files)
├── 06-params/<RES>/              # one parameter set per ligand/cofactor, shared by all systems
├── 06-systems/<ID>/              # E1, T1, E1__T1 ...: topol.top, solv_ions.gro, index.ndx, em.tpr, run.sh
│   ├── run.json, run.log         # G: per-phase state; {em,nvt,npt,md}.{tpr,cpt,gro,log,edr}, md.xtc
│   └── analysis/                 # H: analysis.json, rmsd/rmsf/rg/hbonds.png and a README explaining them
├── 06-systems/matrix.json        # F
├── 07-run-logs/ 07-slurm/        # G: background progress; run_array.sbatch for --slurm
├── 08-analysis/                  # H: analysis_summary.csv, one row per system
└── 09-reports/                   # gromatrix-report: README.md with tables and figures for B-H
All names are defined once, in src/gromatrix/layout.py.
```

## Design notes per step

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
  * log scanning (LINCS, fatal errors, "Finished mdrun") became `checks.scan_md_log`, which is
    aware of log segments and shared by smoke and step G;
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
  * **Ownership.** `run.json` records the host and pid. A second runner skips a system that a live
    runner owns. This is a check, not a lock. A runner that died on this host is reported as
    `interrupted`.
  * **Time budget.** `--maxh` covers the whole invocation, and each `mdrun` gets what is left.
  * **Gates.** Gates also run on phases skipped as done, so an EM left by the smoke test is checked
    too. NVT/NPT LINCS warnings are counted over every log segment. A phase that fails its gate
    keeps `<phase>.failed.gro/.log` as evidence and is cleared, checkpoint included, so it is never
    resumed as complete.
  * **Safety.** Step F refuses to rebuild a system that has MD results unless `--force`.
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
* **Pipeline.** Each step runs in-process through its own CLI `main()`, so the same locks, caches
  and gates apply. The pipeline's own logic is planning and generating the step-F manifest from
  the step-C units. Every `steps.<name>` key in the manifest becomes that command's flag.
