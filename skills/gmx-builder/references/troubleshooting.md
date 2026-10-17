# Troubleshooting: step F (by error text)


| Error contains | Cause | Fix |
|---|---|---|
| `topology/coordinate mismatch after solvate` | `[ molecules ]` order or counts do not match `solv.gro`, e.g. ligand coordinates written in the wrong place, or a stale file | Read the message: it names the first atom range that disagrees. Rebuild with `--force`. |
| `net charge ... is not an integer` | A ligand or residue charge set does not sum to an integer (bad acpype/CGenFF output, hand-edited `.rtp`) | Inspect `06-params/<RES>/meta.json`. Never "fix" this by editing one atom's charge (as SQUIP did for its ZGLY). Reparameterise instead. |
| `after genion: charge` | genion could not neutralise (no SOL to replace, tiny box) | Increase `--pad` |
| smoke `EM did not reach Fmax` / positive Epot | Atom overlap, usually a ligand pose clashing after protonation or a bad altloc branch | Check `em.log`'s "Maximum force ... on atom N" and map N with `index.ndx` |
| smoke `lincs_warnings` > 0 | Unstable geometry in the first picoseconds | Same as above. Do not raise the warning tolerance. |

| Error contains | Cause | Fix |
|---|---|---|
| `Residue 'XXX' not found in residue topology database` (pdb2gmx) | A non-standard residue was classified as protein, or step B was skipped | Hand back to the `gromatrix:structure-prep` skill (step B fixes non-standard residues), or list the residue under `ligands:` if it is a real HET group |
| `atom X in residue Y was not found in rtp entry` | Missing heavy atoms | Hand back to `gromatrix:structure-prep`: step B builds missing heavy atoms |
| `Long bonds and/or missing atoms` | Chain break with no TER | The input skipped steps B/C, which split long gaps with TER. Hand back to `gromatrix:structure-prep`. |
| `.ff not found` | Forcefield port is not installed | Put `<name>.ff/` in the batch dir or on `GMXLIB` |
| `pose could not be mapped onto parameter template` | The pose's tautomer or heavy atoms differ from the SMILES | Fix the SMILES, or give an explicit-H `sdf` for that complex |
| `template formal charge N != manifest charge M` | SMILES charge and declared charge disagree | Correct one of them. Do not delete the `charge` field. |
| `itp net charge ... != formal charge` | acpype/sqm did not converge or picked the wrong charge | Check `06-params/<RES>/param.log`. Often a bad starting geometry: supply an `sdf`. |
| `contains ['FE']` / metal | Organometallic cofactor | Outside scope (MCPB.py). Tell the user. |
| `CGenFF penalties ... > 50` | Poor analogy | ffTK reparameterisation, or a different ligand protonation/tautomer |
| `conflicting atomtype` | Two ligands define one type with different parameters | Usually a mixed batch or a stale cache. Delete `06-params/<RES>` for the offending ligand and rebuild with `--force`. |
| `No such moleculetype SOD`/`NA` | Ion naming does not match the forcefield | Means a forcefield mismatch. Run `verify`. |
| grompp `Too many warnings` at the em gate | Read the warning in `build.log` | Common real causes: atom-name mismatch (stale `conf.gro`, rebuild with `--force`), or a net charge left after genion (ligand charge wrong) |
| `ForcefieldMismatch` | Lock, manifest or cache disagree | Never edit the lock. Start a new batch dir for a different forcefield. |

