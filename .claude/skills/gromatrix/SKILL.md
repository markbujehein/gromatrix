---
name: gromatrix
description: Orchestrator for the gromatrix toolkit, which turns PDB/mmCIF structures or docked poses into a GROMACS simulation matrix (receptor alone, target alone, every receptor+target complex) with one forcefield locked per batch. Use this first whenever the user wants MD/GROMACS systems prepared from PDB IDs or docked complexes, wants to compare several enzymes against several ligands/targets, asks "where is my batch at" or "what's next", or mentions the gromatrix pipeline as a whole. It decides the forcefield, then routes to the stage skills gromatrix-structure-prep (fetch/fix/extract), gromatrix-validate (OpenMM relax + AlphaFold comparison) and gromatrix-gmx-builder (GROMACS topologies), in that order.
---

# gromatrix orchestrator

gromatrix is a set of CLI commands. This skill decides the order and the hard-to-undo choices.
The three stage skills each document only their own commands. Load a stage skill only when the
batch is at that stage.

## How state is handed between skills

All state lives in one **batch directory**. No skill keeps anything in conversation memory that
the next skill needs. Begin every session, and every stage, with:

```bash
gmxprep-status --batch-dir B
```

It prints one JSON line containing:
- `forcefield`: the lock;
- `structures.<ID>`: `fetched`, `fixed`, `unit`, `targets`, `matrix_files`, `relax`, `af`;
- `systems`: built, failed and smoke counts;
- **`next`**: `{skill, command, why}`.

Follow `next`: load the named skill and run its command. Run `gmxprep-status` again after each
stage. The status command reads only these artefacts, so the files *are* the contract:

| Stage (skill) | Writes | The next stage reads |
|---|---|---|
| lock (this skill) | `forcefield.lock.json` | D (family), F (everything) |
| A–C (`gromatrix-structure-prep`) | `raw/`, `fixed/<ID>.pdb` + `.meta.json`, `units/<ID>/unit.json` + `enzyme.pdb` / `target_<T>.pdb` / `complex_<T>[_alt<L>].pdb` | D, E, and the F manifest |
| D–E (`gromatrix-validate`) | `relaxed/<ID>/*.pdb` + `relax.json`, `af/<ID>.json`, `af_report.csv` | F (use relaxed files that passed) |
| F (`gromatrix-gmx-builder`) | `params/<RES>/`, `systems/<ID>/` (`system.json`, `topol.top`, `em.tpr`, `run.sh`), `matrix.json` | the user |

## Decisions to settle with the user before any command

1. **Batch directory.** One batch is one comparison matrix.
2. **Forcefield.** Default `amber` (amber99sb-ildn + GAFF2/AM1-BCC + TIP3P). Use `charmm`
   (CHARMM36m + CGenFF) if the user asks or their earlier work uses it. Lock it before step D:
   `gmxprep-gmx lock --batch-dir B --forcefield amber`. The lock is permanent: a different
   forcefield needs a new batch directory.
3. **Targets and cofactors.** Which HET residue name (or protein chain) is the target, which
   cofactors must stay bound, and each ligand's protonation state (SMILES with charges, or an SDF
   with explicit hydrogens).

## Order

```
lock → gromatrix-structure-prep (A fetch, B fix, C extract)
     → gromatrix-validate (D relax, E AlphaFold)   [optional, recommended]
     → gromatrix-gmx-builder (F manifest, build, verify, smoke)
```

Skip A–C when the user already has clean, docked PDBs; go straight to the gmx builder. Skip D–E
only if the user says so.

## Rules for every stage

- Run the commands. Don't re-implement their logic: they carry the forcefield lock, the shared
  ligand-parameter cache and the validation gates.
- Every command prints one JSON line (`{"ok": ...}`). Read it. Open logs only on failure.
- A hard stop (forcefield mismatch, incomplete matrix, metal-containing ligand, high CGenFF
  penalty, failed relax gate) is a finding to report, not something to bypass with a flag.
- If a command is missing, the package is not installed (`pixi install` or `pip install -e .`).
  Don't improvise a replacement script.

## Final report to the user

Report the forcefield key, then the matrix as a table (receptor rows × target columns, plus the
alone systems), with the status of each system. Then, verbatim:
- every warning from C (unit choice, altloc branches, relabelled ligands),
- every D failure label,
- every E flag,
- every F warning.

Finish with the next command per system: `B/systems/<id>/run.sh`.
