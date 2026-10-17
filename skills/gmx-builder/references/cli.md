# CLI reference: step F (`gromatrix-gmx`)

All commands take `--batch-dir`, print one JSON line on stdout and exit non-zero on failure.

## gromatrix-gmx 

```
gromatrix-gmx lock   --batch-dir B --forcefield NAME
gromatrix-gmx plan   --batch-dir B --manifest M [--allow-incomplete]
gromatrix-gmx build  --batch-dir B --manifest M [--forcefield NAME] [--only ID ...] [--force]
                   [--box dodecahedron|cubic|octahedron|triclinic] [--pad 1.2] [--conc 0.15]
                   [--temp 300] [--md-ns 100] [--maxwarn 0] [--gmx gmx]
                   [--allow-incomplete] [--allow-auto-protonation] [--accept-high-penalty]
gromatrix-gmx verify --batch-dir B [--gmx gmx]
gromatrix-gmx smoke  --batch-dir B [--only ID ...] [--nvt-steps 5000] [--mdrun-args "-nt 8"]
```

| Flag | Meaning |
|---|---|
| `--forcefield` | `amber`→amber99sb-ildn, `amber99sb-ildn-tip4pew`, `amber14sb`, `amber19sb` (OPC), `charmm`→charmm36m. Must equal the lock. |
| `--only` | system ids such as `E1`, `T1`, `E1__T1` |
| `--force` | rebuild systems already stamped `ok`; also required to rebuild a system that has MD results (`run.json`), because rebuilding deletes them |
| `--pad` | minimum solute–box distance (nm). Keep it the same across a matrix. |
| `--conc` | NaCl (KCl not supported) in mol/L on top of neutralisation |
| `--maxwarn` | applies only to the final `em.tpr` gate. Leave it at 0. |
| `--gmx` | e.g. `gmx_mpi` |
| `--allow-auto-protonation` | Open Babel pH 7.4 guess for ligands without SMILES/SDF. Results include a warning. |
| `--accept-high-penalty` | allow CGenFF penalties > 50 |

Environment variables for the CHARMM route: `CGENFF_BIN` (path to `cgenff`), `CGENFF2GMX`
(path to `cgenff_charmm2gmx_py3_nx2.py`). Forcefield ports that are not bundled
(`amber14sb_OL15.ff`, `amber19sb.ff`, `charmm36-jul2022.ff`) go in the batch dir or on `GMXLIB`.

Build gates (any failure stops that system): topology↔coordinate atom walk after `solvate`, integer
solute charge, zero net charge after `genion`, and a zero-warning `grompp` of `em.tpr`. Systems are
skipped on rebuild only when their `inputs_sha` (source files, ligand specs, build flags, build recipe
version) is unchanged. `pdb2gmx` runs with `-renum`, so residues are numbered consecutively and
insertion codes stay distinct residues.
Manifest complexes may carry `variant:` (e.g. `altA`), giving system ids `E__T__altA`.

`smoke` runs EM per system (pass = Fmax < emtol and negative Epot), then `--nvt-steps` of restrained
NVT (pass = finished, no fatal error, zero LINCS warnings), and records the result in `system.json`.
`verify` also re-runs the atom walk and charge check, and fails on a recorded smoke failure.

### Outputs per system `B/systems/<id>/`

`topol.top`, `protein.itp`, `<RES>.itp`, `ligand_params.itp`, `posre.itp`, `conf.gro`, `solv_ions.gro`,
`index.ndx` (System, Protein, Solute, Solvent, LIG_<RES>), `em.tpr`, `mdp/{ions,em,nvt,npt,md}.mdp`,
`run.sh` (em→nvt→npt→md with correct `-r`/`-t` chaining), `<ff>.ff` symlink for non-bundled ports,
`system.json` (stamp, with `target_resname`, `target_residues` and `protein_residue_labels` for
step H), `build.log`, `error.txt` (on failure).

Ligand cache `B/params/<RES>/`: `template.sdf` (atom order is the order used in the topology),
`template.mol2`, `<RES>.itp`, `params.json`, `meta.json`, `param.log`.

