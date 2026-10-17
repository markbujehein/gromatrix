---
name: validate
description: "Checks that extracted structures are physically sane before MD: two-stage restrained OpenMM relaxation with failure gates, then comparison of each protein chain with its AlphaFold DB prediction (steps D-E). Use this skill ONLY when gromatrix-status shows steps A-C done and D or E next, or the user asks whether prepared structures are sane or how they compare with AlphaFold. Covers gromatrix-relax (explosion, distortion, chain-break, clash and pocket-collapse gates) and gromatrix-afcompare (RMSD, TM-score, Calpha-lDDT, flags). Not for fetching structures, GROMACS builds or MD."
---

# gromatrix: validation (steps D–E)

Shared biophysical rules: `${CLAUDE_PLUGIN_ROOT}/references/domain-rules.md`.

Commands: `gromatrix-relax`, `gromatrix-afcompare`. Full flags are in `references/cli.md`. Failure
labels and flags are explained in `references/troubleshooting.md`.

## Before you start

```bash
gromatrix-status --batch-dir B     # needs 03-units/<ID>/ from gromatrix:structure-prep, and the forcefield lock
```

Step D reads the forcefield **family** from the lock, so it relaxes with amber14 or charmm36 to
match the later GROMACS build. If `forcefield` is `null`, go back to the orchestrator to lock it
first. An explicit `--forcefield` that disagrees with the lock is refused.

## Run

```bash
gromatrix-relax     --batch-dir B          # enzyme.pdb, complex_*.pdb, target_*.pdb; ligand-only files are skipped
gromatrix-afcompare --batch-dir B          # every protein chain ≥ 30 residues vs its AFDB model
```

`gromatrix-relax` picks its solvent model itself: GBn2 implicit solvent on a GPU, a
reaction-field cutoff on CPU (about 25× faster). Expect about 10–15 s per file on a 4-core CPU.

## Read the results

- **relax** `failed[]`: every entry has a label (`explosion`, `implosion/explosion`, `distortion`,
  `chain break`, `unresolved clash`, `pocket collapse`). Report them verbatim. The fix is almost
  always upstream (a bad altloc branch, a mis-built loop), so send the user back to
  structure-prep. Do not loosen `--max-rmsd`.
- **afcompare** `flagged[]` (TM < 0.8, lDDT < 0.7, identity < 0.9):
  - A low `rmsd_core` with a high `rmsd_all` is normal: flexible regions have low pLDDT.
  - A low TM-score with a high lDDT means a domain or hinge rearrangement, not a local error.
  - `no UniProt mapping` is expected for engineered binders.

## Hand-off

Done when `gromatrix-status` shows `relax` and `af` for every entry. Tell the orchestrator:
- which files passed. The gmx builder uses `04-relaxed/<ID>/<file>` for those.
- which files failed. They stay as `*.failed.pdb` and must not go into the manifest.
- the AlphaFold flags, for the final report.
