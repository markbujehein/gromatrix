---
name: gmx-builder
description: "Builds GROMACS systems for a receptor by target matrix under one locked forcefield (step F): pdb2gmx topologies, shared ligand and cofactor parameters (acpype GAFF2 for AMBER, CGenFF for CHARMM), ions, box, solvent, a zero-warning em.tpr, verify and smoke test. Use this skill ONLY when structures are ready (or the user brings docked PDBs) and no systems are built yet, or a pdb2gmx, grompp, acpype or CGenFF error must be fixed, or gmxprep-status names gromatrix:gmx-builder as next. Covers the manifest, gmxprep-gmx lock, plan, build, verify and smoke, and forcefield consistency. Not for running or analysing MD."
---

# gromatrix: GROMACS builder (step F)

Shared biophysical rules: `${CLAUDE_PLUGIN_ROOT}/references/domain-rules.md`. Compute needs are in `${CLAUDE_PLUGIN_ROOT}/references/compute-modes.md`: steps A–F need no GPU.

Command: `gmxprep-gmx lock | plan | build | verify | smoke`. Full flags and outputs are in
`references/cli.md`. Forcefield trade-offs and ligand routing are in `references/forcefields.md`.
Error texts are explained in `references/troubleshooting.md`.

## Before you start

```bash
gmxprep-status --batch-dir B
```

The lock must exist. If it doesn't, the orchestrator decides the forcefield first. Inputs are the
step-C files in `03-units/<ID>/`, or the step-D files in `04-relaxed/<ID>/` where those passed. Never
use `*.failed.pdb`. The user may also bring their own docked PDBs and skip A–E.

## Write the manifest

Copy `${CLAUDE_PLUGIN_ROOT}/examples/manifest.yaml`. Paths are relative to the manifest file.

- `enzymes:` → `03-units/<ID>/enzyme.pdb` (or the relaxed copy).
- `targets:` → `kind: ligand` with `resname` (pose from `target_<T>.pdb` or an explicit-H `sdf`),
  or `kind: protein` with `chains`. Every ligand target needs a **unique** resname.
- `complexes:` → one per enzyme × target pair (`complex_<T>.pdb`). For altloc branches add one
  entry per branch with `variant: altA`, `variant: altB`, ...
- `ligands:` → every non-ion HET residue name (targets *and* cofactors) needs `smiles` + `charge`,
  or an explicit-H `sdf`. This fixes the protonation state used by every system.

```bash
gmxprep-gmx plan  --batch-dir B --manifest M    # system list = receptors + targets + every pair
gmxprep-gmx build --batch-dir B --manifest M
gmxprep-gmx verify --batch-dir B
gmxprep-gmx smoke  --batch-dir B                # EM + short NVT; needs GROMACS run time
```

A pair without a complex aborts `plan`/`build`. Use `--allow-incomplete` only with the user's
consent, because it breaks the matrix.

## Read the results

- `build → failed[]`: the traceback is in `06-systems/<id>/error.txt` and tool output in `build.log`.
  Fix the cause and rebuild with `--only <id>`. Systems with unchanged inputs are skipped.
- Changing a ligand's SMILES, charge or file makes the cached parameters stale. The build then
  stops and asks for `--force`. That's intended: every system sharing the ligand must be rebuilt.
- Hard stops are deliberate: `ForcefieldMismatch`, a metal-containing ligand (needs MCPB.py), a
  CGenFF penalty > 50, conflicting ligand parameters, a non-integer charge. Report them. Don't
  raise `--maxwarn`.
- `smoke`: a system passes when EM converges (Fmax < emtol, negative Epot) and a 10 ps NVT run
  has zero LINCS warnings.

## Hand-off

Done when `verify` is `"ok": true`, and `smoke` passed if it was run. Return to the orchestrator;
`gmxprep-status` then names `gromatrix:execution` (MD runs) as `next`. Each system also has a
dependency-free `run.sh` (em → nvt → npt → md), for users who run without gromatrix installed.
`system.json` records `target_resname` / `target_residues` for step H.
