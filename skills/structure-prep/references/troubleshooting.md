# Troubleshooting: steps A–C (by message text)


| Error / warning contains | Cause | Fix |
|---|---|---|
| `no candidate unit contains all required ligands` | The required HET is bound in no single assembly or stoichiometric unit, e.g. it sits only at a crystal contact. | Open `03-units/<ID>/unit.json` and compare the candidates' `inventory` and `missing`. Check the HET really binds the protein (`--ligand-cutoff`). As a last resort, use `--unit <id> --allow-incomplete` and dock the ligand into the chosen unit. |
| `dropped HET groups touching no chain` | Orphan ligand more than `--ligand-cutoff` Å from every chain in the unit (crystal packing, or a ligand bound to a symmetry mate) | Usually correct to drop. If it is the target, the deposit's assembly is wrong for it: try `--assembly none`. |
| `selected unit lacks HET types present elsewhere` | Coverage below 1. Another copy carries a ligand this unit lacks (e.g. a ligand bound in only one copy, without `--targets`). | Say which ligands are targets via `--targets`/`--expect`. Selection then guarantees them. |
| `selected the whole asymmetric unit` | No smaller connected stoichiometric unit was complete | Check `unit.json`. Often one ligand is bound across two copies. Simulating the AU is valid but costs more. |
| `no candidate units (assembly 'N' not found?)` | `--assembly N` does not exist | Pick an id listed in `02-fixed/<ID>.meta.json` |
| `every candidate exceeds --max-chains` | Large assemblies (capsids, chaperonins) | Raise `--max-chains` or use `--assembly none` |
| `needs more than 62 chain IDs` | Assembly too large for PDB format | Use a smaller assembly |
| `alternate conformers: ... equivalent within coordinate error` | Conformers differ by a label permutation (ring flip) or less than `--altloc-rmsd` (e.g. raw 1.0 Å but symmetry-corrected 0.4 Å for a flipped para-substituted ring) | Nothing to do: the primary conformer is used. To compare both anyway: `gmxprep-extract --altloc branch`. |
| `alternate conformers: ... distinct poses` | The deposit has two genuinely different binding modes | Step C branches the matrix automatically. Add each `complex_<T>_alt<L>.pdb` as a manifest complex with `variant: alt<L>`. Never average the poses. |
| `INCONSISTENT GEOMETRY - excluded` | A conformer's bond lengths differ from the primary's by > 0.2 Å (broken or averaged model) | That conformer is never used. Report it to the user. |
| `severe clash with environment` (unit.json → altlocs) | A branch pose overlaps the protein in that state | The branch is skipped. Usually an occupancy-coupled side chain carries a different letter. |
| `internal gaps left open` | A missing loop longer than `--max-loop` | Charged termini at the gap are an approximation. For loops near the site of interest, model them properly (MODELLER/AF) instead of raising `--max-loop`. |
| `could not map label->author chain IDs` | Unusual mmCIF | Inspect chain IDs in `02-fixed/<ID>.pdb`. Ligands may sit on the wrong chain letter. |
| `HETs whose deposited chain label differs from the chain they bind` | The deposit assigns a ligand to one chain, but it sits in another chain's pocket | Informational: assembly membership already follows the physical binding chain. Check the selected unit holds the intended pocket. |
| `altloc L not used - pose does not touch this unit` | The deposit stores a second copy of the ligand, bound elsewhere, as an altloc of the same residue | Correct to drop for this unit. If that other site is the one you want, select the unit containing it (`--unit`). |
| `multi-character chain IDs` | Very large deposit | Not supported in PDB output |

