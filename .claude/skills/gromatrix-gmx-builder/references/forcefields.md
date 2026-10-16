# Forcefield choices

## Registry

| key | protein FF (.ff dir) | water | ligands | bundled with GROMACS |
|---|---|---|---|---|
| `amber99sb-ildn` (alias `amber`) | amber99sb-ildn | TIP3P | GAFF2 + AM1-BCC via acpype | yes |
| `amber99sb-ildn-tip4pew` | amber99sb-ildn | TIP4P-Ew (4-site) | GAFF2 + AM1-BCC | yes |
| `amber14sb` | amber14sb_OL15 port | TIP3P | GAFF2 + AM1-BCC | no |
| `amber19sb` | amber19sb port | OPC (4-site, solvated from tip4p.gro) | GAFF2 + AM1-BCC | no |
| `charmm36m` (alias `charmm`) | charmm36-jul2022 port | CHARMM TIP3P | CGenFF via cgenff_charmm2gmx | no |

How to choose:
- Default to `amber`: it is bundled, robust and widely benchmarked with GAFF2 ligands.
- `amber19sb` + OPC is the most accurate modern AMBER protein model. It costs more (4-site water)
  and needs the community port.
- `charmm36m` suits disordered regions, membranes and lipids, or groups that already use CHARMM/CGenFF.

To add a forcefield, add one entry in `src/gromatrix/forcefields.py`. It must declare the ligand
engine, the ion moleculetypes and the mdp non-bonded block.

## Why non-bonded settings are tied to the forcefield

- AMBER: 1.0 nm plain cut-offs, PME, `DispCorr = EnerPres`.
- CHARMM36(m): 1.2 nm cut-off with a 1.0–1.2 nm `Force-switch`, `DispCorr = no`.

Using the AMBER settings with CHARMM (or the reverse) quietly changes densities and lipid/protein
behaviour. `verify` checks every mdp against the lock.

## Ligand routing

| family | engine | charges | gate |
|---|---|---|---|
| amber | `acpype -a gaff2 -c bcc` | AM1-BCC (sqm), conformer from the first pose | itp net charge must equal the template's formal charge |
| charmm | `cgenff` (or the server's `.str`) → `cgenff_charmm2gmx` | by analogy | penalty > 50 stops; 10–50 gives a warning |

Both engines are given the same explicit-H template, built from the manifest SMILES/SDF plus the
first crystallographic pose. Protonation states are therefore identical in every system and under
either forcefield.

Not supported by GAFF2 or CGenFF, so the build stops: organometallic cofactors (heme Fe, Cu sites,
Fe–S clusters) and covalently bound ligands. Free ions (Na, Cl, K, Mg, Ca, Zn) are mapped to the
forcefield's own `ions.itp`. A Zn ion in a structural site is treated as non-bonded; mention this
to the user.
