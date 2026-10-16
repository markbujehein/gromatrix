# Troubleshooting: step G

| `reason` / symptom | Cause | What to do |
|---|---|---|
| `EM did not converge (Fmax=…)` | Atoms overlap: a bad docked pose, a ligand placed into the protein, or a missing residue rebuilt through the ligand | Check the complex in the relaxed/units file. Fix the pose upstream, rebuild with `gmxprep-gmx build --only <id> --force`, rerun. |
| `nvt: N LINCS warnings during restrained equilibration` | The structure blows up under dynamics: bad ligand parameters (check the CGenFF penalty or the acpype log), a bad protonation state, or an EM that "converged" in a strained geometry | Read `nvt.log` around the first warning; the atom numbers point to the residue. Fix the parameters or the pose; never raise the LINCS limits. |
| `npt: …` failure after NVT passed | Density or box problem: the box is too small for the solute, or vacuum bubbles | Rebuild with a larger `--pad` in step F. |
| `mdrun … exited with 1` | GROMACS error (GPU not available, a bad `-gpu_id`, an out-of-range thread count, a corrupt checkpoint) | Read the tail in `reason` and `run.log`. GPU problems: drop `--gpu-id` or fix `--mdrun-args`. |
| `stopped before nsteps (signal or -maxh)` | Walltime reached, Ctrl-C, SIGTERM from the scheduler | Normal. Rerun or resubmit; the phase continues from its `.cpt`. |
| `busy: running on <host>` | `run.json` says a runner on another host owns the system | Wait. A dead runner on *this* host is detected automatically. For another host, run `--only <id>` there, or ask the user before deleting `run.json`. |
| `smoke test failed; fix it or pass --force` | Step F's smoke test failed | Go back to `gromatrix-gmx-builder`. |
| `nvt`/`npt` rerun unexpectedly | Its checkpoint was deleted by hand: a phase only counts as done while the checkpoint the next phase reads exists | Nothing to do: the phase reruns, and every later phase starts fresh. |
| `--maxh time budget used up` | The budget ended between phases | Normal. Rerun or resubmit. |

The production phase writes `md.xtc` with the `nstxout-compressed` set by step F (every 10 ps).
Step H reads that file.
