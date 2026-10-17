# Biophysical rules the workflow relies on

Treat these as hard limits when interpreting or reporting. They are conventions of careful MD
practice, not options to be tuned until the answer looks good.

## Preparation (A-F)
- **One forcefield per matrix.** AMBER: amber99sb-ildn, GAFF2 with AM1-BCC charges, TIP3P. CHARMM:
  CHARMM36m, CGenFF, TIP3P-CHARMM. Mixing within a comparison makes differences meaningless; the lock
  enforces it.
- **Protonation is a decision.** Every non-ion HET needs a SMILES with explicit charge (or an H-complete
  SDF). Carboxylates in a pocket are normally ionised; state any assumption in the report.
- **Alternate conformers.** The highest-occupancy conformer is the primary one. Branch the matrix only when
  conformers are genuinely distinct (symmetry-corrected RMSD of at least 0.5 Å), clash-free and inside the unit.
  Never average coordinates. A ring flip is not a distinct conformer.
- **Biological unit.** Prefer the author-defined assembly holding every required ligand. A ligand present
  in only one copy of a multi-chain structure is a finding to report.
- **Relax gates are not tunable.** Do not raise `--max-rmsd` to make a structure pass; fix the cause upstream.
- **Metals and CGenFF.** Metal-containing ligands need MCPB.py and stop the build. A CGenFF penalty above 50
  stops it. Never raise `-maxwarn`.

## Simulation (G)
- Protocol: EM (must converge, negative Epot), NVT and NPT with position restraints, then unrestrained
  production. Any LINCS warning in NVT/NPT fails the phase; production LINCS warnings are reported.
- A failed gate is evidence of a bad pose, clash or parameters. Do not shorten the timestep or add
  `-maxwarn` to hide it.

## Interpretation (H)
- RMSD of the protein Calpha near 1-3 Å and a target RMSD under about 3 Å in the protein frame is the
  usual picture of a stable bound pose. A rising curve means not yet equilibrated or drifting.
- **Sampling.** One replica of N ns is a single sample. Under about 10 ns it shows equilibration and
  pose retention at most; it cannot support a binding-affinity or stability ranking. Say so every time.
- Hydrophobic ligands may form no persistent H-bond; `no_target_hbonds` is a finding, not an error.
- Flags (`ligand_drift`, `protein_unstable`, `no_target_hbonds`) are reported verbatim and never suppressed.
- Never present synthetic example numbers, or numbers from a system that was not simulated, as results.
