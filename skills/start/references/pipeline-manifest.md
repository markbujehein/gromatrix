# Pipeline manifest (gromatrix-pipeline)

A template for one-command runs. Copy it, change `batch_dir`, `steps.fetch.ids`,
`steps.extract.targets` and `ligands`, then preview it with `--dry-run`. Every key under `steps:`
becomes that command's flag (`max_loop: 8` → `--max-loop 8`; `true` → a bare flag; `false` skips
the optional steps relax / afcompare / smoke).

```yaml
# gromatrix-pipeline manifest: steps A-H from PDB IDs. Paths are relative to this file.
#   gromatrix-pipeline --manifest examples/pipeline.yaml --dry-run          # print the exact commands
#   gromatrix-pipeline --manifest examples/pipeline.yaml --run-md --analyze
batch_dir: ../gromatrix-work/rxr_car
forcefield: amber

# Step C keeps only native (crystal) pairs, so the generated step-F matrix is usually
# incomplete. Set this only if comparing native complexes is what you want.
allow_incomplete: true

steps:                         # each key becomes that command's --flag
  fetch:     {ids: [1XVP]}     # present = start from PDB IDs; omit to bring your own PDBs (below)
  fix:       {}
  extract:   {targets: CID}
  relax:     {}                # false to skip (same as --skip-relax)
  afcompare: {}                # false to skip (same as --skip-afdb)
  build:     {md_ns: 100}
  smoke:     {nvt_steps: 5000} # false to skip
  run:       {threads: 8, gpu_id: 0}
  analyze:   {stride: 2}

# Optional: name targets (else derived from step C: ligand resnames, chains_<X> for chains).
targets:
  - {id: CITCO, kind: ligand, resname: CID}

# Required for every non-ion HET (targets and cofactors): fixes the protonation state.
ligands:
  CID: {smiles: "Clc1ccc(cc1)-c1nc2sccn2c1/C=N/OCc1ccc(Cl)c(Cl)c1", charge: 0}   # CITCO
  # Every other HET kept in the selected unit (for 1XVP, the F15 in the RXR partner's pocket)
  # needs an entry too, or step F stops and names it.

# To skip A-E with your own docked PDBs, drop steps.fetch and add enzymes/targets/complexes
# exactly as in examples/manifest.yaml.
```

The manifest is the step-F manifest plus pipeline keys. Skip A–E with your own docked PDBs by
dropping `steps.fetch` and adding `enzymes` / `targets` / `complexes` (the gmx builder documents
that format). Paths are relative to the manifest file.
