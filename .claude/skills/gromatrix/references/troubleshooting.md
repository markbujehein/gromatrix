# Troubleshooting (by error text)

## Steps B and C

| Error / warning contains | Cause | Fix |
|---|---|---|
| `no candidate unit contains all required ligands` | The required HET is bound in no single assembly or stoichiometric unit, e.g. it sits only at a crystal contact. | Open `units/<ID>/unit.json` and compare the candidates' `inventory` and `missing`. Check the HET really binds the protein (`--ligand-cutoff`). As a last resort, use `--unit <id> --allow-incomplete` and dock the ligand into the chosen unit. |
| `dropped HET groups touching no chain` | Orphan ligand more than `--ligand-cutoff` Å from every chain in the unit (crystal packing, or a ligand bound to a symmetry mate) | Usually correct to drop. If it is the target, the deposit's assembly is wrong for it: try `--assembly none`. |
| `selected unit lacks HET types present elsewhere` | Coverage below 1. Another copy carries a ligand this unit lacks (e.g. a ligand bound in only one copy, without `--targets`). | Say which ligands are targets via `--targets`/`--expect`. Selection then guarantees them. |
| `selected the whole asymmetric unit` | No smaller connected stoichiometric unit was complete | Check `unit.json`. Often one ligand is bound across two copies. Simulating the AU is valid but costs more. |
| `no candidate units (assembly 'N' not found?)` | `--assembly N` does not exist | Pick an id listed in `fixed/<ID>.meta.json` |
| `every candidate exceeds --max-chains` | Large assemblies (capsids, chaperonins) | Raise `--max-chains` or use `--assembly none` |
| `needs more than 62 chain IDs` | Assembly too large for PDB format | Use a smaller assembly |
| `alternate conformers: ... equivalent within coordinate error` | Conformers differ by a label permutation (ring flip) or less than `--altloc-rmsd` (e.g. raw 1.0 Å but symmetry-corrected 0.4 Å for a flipped para-substituted ring) | Nothing to do: the primary conformer is used. To compare both anyway: `gmxprep-extract --altloc branch`. |
| `alternate conformers: ... distinct poses` | The deposit has two genuinely different binding modes | Step C branches the matrix automatically. Add each `complex_<T>_alt<L>.pdb` as a manifest complex with `variant: alt<L>`. Never average the poses. |
| `INCONSISTENT GEOMETRY - excluded` | A conformer's bond lengths differ from the primary's by > 0.2 Å (broken or averaged model) | That conformer is never used. Report it to the user. |
| `severe clash with environment` (unit.json → altlocs) | A branch pose overlaps the protein in that state | The branch is skipped. Usually an occupancy-coupled side chain carries a different letter. |
| `internal gaps left open` | A missing loop longer than `--max-loop` | Charged termini at the gap are an approximation. For loops near the site of interest, model them properly (MODELLER/AF) instead of raising `--max-loop`. |
| `could not map label->author chain IDs` | Unusual mmCIF | Inspect chain IDs in `fixed/<ID>.pdb`. Ligands may sit on the wrong chain letter. |
| `multi-character chain IDs` | Very large deposit | Not supported in PDB output |

## Step F

| Error contains | Cause | Fix |
|---|---|---|
| `topology/coordinate mismatch after solvate` | `[ molecules ]` order or counts do not match `solv.gro`, e.g. ligand coordinates written in the wrong place, or a stale file | Read the message: it names the first atom range that disagrees. Rebuild with `--force`. |
| `net charge ... is not an integer` | A ligand or residue charge set does not sum to an integer (bad acpype/CGenFF output, hand-edited `.rtp`) | Inspect `params/<RES>/meta.json`. Never "fix" this by editing one atom's charge (as SQUIP did for its ZGLY). Reparameterise instead. |
| `after genion: charge` | genion could not neutralise (no SOL to replace, tiny box) | Increase `--pad` |
| smoke `EM did not reach Fmax` / positive Epot | Atom overlap, usually a ligand pose clashing after protonation or a bad altloc branch | Check `em.log`'s "Maximum force ... on atom N" and map N with `index.ndx` |
| smoke `lincs_warnings` > 0 | Unstable geometry in the first picoseconds | Same as above. Do not raise the warning tolerance. |

| Error contains | Cause | Fix |
|---|---|---|
| `Residue 'XXX' not found in residue topology database` (pdb2gmx) | A non-standard residue was classified as protein, or step B was skipped | Run `gmxprep-fix`, or list the residue under `ligands:` if it is a real HET group |
| `atom X in residue Y was not found in rtp entry` | Missing heavy atoms | Rerun step B (PDBFixer `addMissingAtoms`) |
| `Long bonds and/or missing atoms` | Chain break with no TER | Step B builds gaps up to `--max-loop` and writes TER at longer ones. If this appears, the input skipped steps B/C. |
| `.ff not found` | Forcefield port is not installed | Put `<name>.ff/` in the batch dir or on `GMXLIB` |
| `pose could not be mapped onto parameter template` | The pose's tautomer or heavy atoms differ from the SMILES | Fix the SMILES, or give an explicit-H `sdf` for that complex |
| `template formal charge N != manifest charge M` | SMILES charge and declared charge disagree | Correct one of them. Do not delete the `charge` field. |
| `itp net charge ... != formal charge` | acpype/sqm did not converge or picked the wrong charge | Check `params/<RES>/param.log`. Often a bad starting geometry: supply an `sdf`. |
| `contains ['FE']` / metal | Organometallic cofactor | Outside scope (MCPB.py). Tell the user. |
| `CGenFF penalties ... > 50` | Poor analogy | ffTK reparameterisation, or a different ligand protonation/tautomer |
| `conflicting atomtype` | Two ligands define one type with different parameters | Usually a mixed batch or a stale cache. Delete `params/<RES>` for the offending ligand and rebuild with `--force`. |
| `No such moleculetype SOD`/`NA` | Ion naming does not match the forcefield | Means a forcefield mismatch. Run `verify`. |
| grompp `Too many warnings` at the em gate | Read the warning in `build.log` | Common real causes: atom-name mismatch (stale `conf.gro`, rebuild with `--force`), or a net charge left after genion (ligand charge wrong) |
| `ForcefieldMismatch` | Lock, manifest or cache disagree | Never edit the lock. Start a new batch dir for a different forcefield. |
