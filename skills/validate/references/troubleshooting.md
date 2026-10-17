# Troubleshooting: steps D–E

## Step E (AlphaFold comparison)

| Message | Cause | Fix |
|---|---|---|
| `no UniProt mapping` | Engineered or synthetic chain (e.g. designed binders), or neither PDBe nor RCSB maps it | Expected for non-natural chains. For natural ones pass `--uniprot CHAIN:ACC`. |
| `sequence identity < 0.9` flag | Wrong accession (e.g. a fusion construct), or heavy engineering | Fusion chains use the accession covering most of the chain. Override with `--uniprot` if that is still wrong. |


## Step D (relax failure labels)

| Label | Meaning | What to do |
|---|---|---|
| `explosion` | Non-finite energy or coordinates | Almost always overlapping atoms in the input. Look at alternate-conformer branches and modelled loops in step B/C output. |
| `implosion/explosion` | Heavy-atom RMSD > `--max-rmsd` with a large energy drop (< −`--max-drop` kJ/mol/atom) | Starting geometry far from a minimum (bad loop, wrong altloc). Fix the input in B/C. Do not loosen the gate. |
| `distortion` | RMSD > `--max-rmsd` without a dramatic energy drop | Often a long flexible tail or a domain held only by crystal contacts. Inspect `04-relaxed/<ID>/<file>.failed.pdb`. Consider a different unit in step C. |
| `chain break` | A Cα–Cα distance stretched beyond 4.2 Å | A gap not split by TER, or a mis-modelled loop. Rerun step B with a smaller `--max-loop`. |
| `unresolved clash` | Final potential energy still positive | Overlap the minimiser could not resolve. Usually a ligand-adjacent side chain from a different altloc. |
| `pocket collapse` | Protein moved ≥ 0.5 Å closer to a held-out ligand and is now < 2.0 Å from it | The empty site closed in. Increase `--pocket-cutoff`, or check that the ligand belongs to this unit. |

Failed structures are written as `*.failed.pdb`, so no later step picks them up. Step E skips
them automatically.
