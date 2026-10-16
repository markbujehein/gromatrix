# Interpreting step-H results

## RMSD

- **Cα RMSD** is measured from the first *production* frame, after NVT/NPT. A plateau of
  1–3 Å is typical for a folded domain. A steady climb means the run hasn't converged; long
  flexible termini or loops inflate it. Check `rmsf.png` before calling a protein unstable.
- **Target RMSD** is computed in the protein's frame, without refitting the ligand, so it
  includes translation and rotation in the pocket. Under 2 Å: the pose is retained. 2–3 Å: it
  rearranged inside the pocket. Above 3 Å (`ligand_drift`): it moved out or flipped. Check the
  H-bond panel and look at the trajectory before calling it unbound.
- Symmetric ligands (a flipped phenyl, a carboxylate) can show a ~1–2 Å RMSD step without any
  real change.

## RMSF

The x axis is the sequential residue index of the simulated protein (all chains, in pdb2gmx order),
not the author numbering. `analysis.json` `rmsf.resid` holds the residue numbers. Termini and
rebuilt loops (step B, `--max-loop`) are usually the peaks.

## H-bonds

- These are geometric (D–A < 3.0 Å, D–H–A > 150°), with no energy criterion. Occupancy is the
  fraction of frames in which the specific donor→acceptor pair is present.
- More than 50 % is a persistent interaction worth naming; 10–50 % is transient; under 10 % is
  counted but not listed.
- Donors are paired with hydrogens by distance (≤ 1.2 Å), so no bond information is needed.
  Atom names come from the topology: `LIG12:O1` is atom O1 of residue LIG with GROMACS residue
  number 12.

## Radius of gyration

A stable Rg (within about 0.5 Å) means the fold is compact. A drift upwards can mean partial
unfolding, or a domain moving apart. Look at it alongside Cα RMSD.

## Comparing systems in the matrix

- Compare `*_last_half` values. The first half includes relaxation away from the crystal
  structure.
- With one replica per system, differences of under ~0.5 Å RMSD or ~10 % occupancy are within
  noise. Recommend replicas (more seeds: rebuild or rerun NVT with a new `gen_seed`) before
  ranking ligands.
- These metrics describe stability, not affinity. Binding free energies need FEP/MM-PBSA, which
  is outside gromatrix.
