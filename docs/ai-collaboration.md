# How gromatrix was built: the human pilot and the AI engine

gromatrix is an open demonstration of AI-assisted scientific software engineering. One person and
Claude (Anthropic's model, working through Claude Code) built it with a deliberate split of
responsibility. This page records how that worked, with evidence from the git history and the
code. It also records where it didn't work first time.

## Roles

| | Human: domain architect, biophysics lead, pipeline strategist, quality assurer | AI (Claude): primary developer, implementer, test author, skill integrator |
|---|---|---|
| **What to build** | The simulation-matrix concept; the forcefield rules (AMBER preferred, CHARMM supported, no mixing within a batch); ligand routing (acpype for AMBER, CGenFF for CHARMM); the order of steps A–H; which benchmark proves it | Turned each brief into module boundaries, CLI contracts and file formats |
| **Domain judgement** | Named the hard cases: 1XVP's ligand present in only one copy, a single conformer picked arbitrarily, over-fitting to one benchmark; set physical thresholds (e.g. the 2.5 Å relax gate) | Investigated each case in the data, quantified it, and proposed and implemented a general rule |
| **Code** | Reviewed the results and decided scope at each review round | Wrote all of `src/`, `tests/`, the CI workflows, `scripts/` and the six agent skills |
| **Quality** | Required the 82-structure benchmark, independent review agents, pixi-based CI, honest reporting of what is untested | Ran the benchmark and two rounds of review agents, verified every finding before fixing it, and added a regression test per fix |
| **Release** | Owns versioning (held at 0.x alpha), merges and tags | Prepares changelogs, release workflow and PR text; never merges or tags on its own |

Every commit is in the maintainer's name and carries a `Co-Authored-By: Claude` trailer, plus a
link to the session that produced it.

## The agentic control plane

An agent driving an eight-step scientific pipeline fails in two ways: it loses track of state
across a long session, or its context fills with documentation for steps it isn't running. Both
were designed out.

- **State lives on disk, not in the conversation.**
  - Every step writes its result to the batch directory (`unit.json`, `relax.json`,
    `system.json`, `run.json`, `analysis.json`).
  - `gmxprep-status` reads only those files. It returns a verdict per step A–H plus
    `next: {skill, command, why}`.
  - A fresh session, or a different agent, resumes any batch at the right step.
- **Progressive disclosure.**
  - The `gromatrix` orchestrator skill holds only the hard-to-undo decisions: forcefield, targets
    and protonation, matrix completeness, compute. It also sets the step order.
  - Five stage skills each document only their own commands, with the details in
    `references/`.
  - A test (`tests/test_plugin.py`) fails if a stage skill's code blocks run another
    stage's commands.
- **Commands carry the logic; skills carry the judgement.**
  - Skills tell the agent to run commands and read their one-line JSON. They tell it never to
    re-implement a step, and never to work around a gate with a flag.
  - The forcefield lock, parameter caches and validation gates live in code, where an agent can't
    talk its way past them.
- **Distribution is checked like code.** `scripts/check_plugin.py` applies the same frontmatter
  rules as the Claude skill upload, validates the plugin manifest and agents, and CI runs it on every push. It
  has twice caught an invalid-YAML skill description before any user saw it.

## Human-in-the-loop domain interventions

These are the points where the human's biophysical judgement changed the design. Each is quoted
from the request, then followed into the code.

### 1. "Randomly picking one conformer is biophysically unsound"

- **What happened.** The first implementation of step C met the CITCO ligand in 1XVP. It was
  modelled as two alternate conformers at 0.5 occupancy each, and the code simply kept conformer
  A. The human rejected that, and also asked whether a third, averaged conformer existed, and
  whether averaged coordinates could be simulated at all.
- **What the investigation found.**
  - There is no conformer C.
  - A and B differ by 1.03 Å RMSD, but only 0.38 Å once symmetry-equivalent atoms are matched. It
    is a ring flip within coordinate error, not a second pose.
  - Averaging A and B would collapse aromatic C–C bonds to 0.70 Å. No forcefield can parameterise
    that geometry.
- **What it became** (`src/gromatrix/altlocs.py`, `extract.py`):
  - **Primary conformer:** the highest occupancy.
  - **Symmetry-corrected comparison:** conformers are compared by RMSD after per-element Hungarian
    matching (`scipy.optimize.linear_sum_assignment`).
  - **Distinct poses:** a conformer is a separate pose at a symmetry-corrected RMSD of 0.5 Å or
    more (`--altloc-rmsd`).
  - **Geometry check:** a conformer is unusable if its bond lengths disagree by more than 0.2 Å.
  - **Branching:** distinct poses branch the matrix into `complex_<T>_alt<L>.pdb`, but only if they
    touch the unit and don't clash below 1.5 Å.
  - **No averaging:** coordinates are never averaged.
- **Benchmark payoff.**
  - 1ILH: three genuine SR12813 orientations, branched into three complexes.
  - 1XV9: a ligand copy 48 Å away, mislabelled as "altloc B", rejected with an explanation.

### 2. The 1XVP problem: a target that exists in only one copy

- **The brief.** The human specified this failure before any extraction code existed. 1XVP's
  asymmetric unit holds two RXRα/CAR heterodimers, but CITCO is bound only in chain D. *"If a
  naive extraction script simply grabs the first heterodimer (chains A+B), the resulting
  simulation matrix will completely lack the CITCO target ligand."*
- **What it became** (`src/gromatrix/extract.py`):
  - **Candidates:** every deposited assembly, plus every stoichiometric, interface-connected unit
    of the asymmetric unit, is a candidate.
  - **Ligand assignment:** each HET is assigned to the chain it contacts most, using a
    `scipy.spatial.cKDTree` within 4.5 Å.
  - **Ranking is completeness-first:** a candidate that holds every declared target always wins.
    After that: target coverage; author-defined assemblies; fewer chains; larger interfaces; fewer
    unmodelled residues; lower Cα B-factor.
  - **Result on 1XVP:** assembly 2 (C, D, G, H), the only unit with both ligands, selected with
    and without assembly records.
- **The benchmark sharpened it further.** In 1XV9 the ligand is deposited with chain D's label
  but sits in chain B's pocket (121 contacts with B, none with D). Trusting the label made every
  assembly "incomplete". Assembly membership now follows physical contact.

### 3. "Mixing forcefields within a comparison matrix is strictly forbidden"

- **The requirement.** It was the human's first and firmest one. A comparison between an
  enzyme alone, a target alone and their complex is meaningless if the three were built with
  different physics.
- **What it became:**
  - **The lock:** `forcefields.py` writes one `forcefield.lock.json` per batch. Every
    family-dependent step reads and checks it: relax uses the lock's family, build checks the
    lock, and verify re-checks every system.
  - **Stamping:** every ligand parameter set and every built system carries the forcefield
    identity, and `gmxprep-gmx verify` fails on any mismatch.
  - **One template per ligand:** each ligand is parameterised once per batch (acpype/GAFF2/AM1-BCC
    or CGenFF, penalty-gated at 50), and every pose is mapped onto that template's atom order. So
    "ligand alone" and "ligand in complex" are the same molecule.
  - **Parameter conflicts:** an AI review round later found that two ligands defining the same
    bonded type differently would be silently mixed by `grompp`. That is now a hard stop, compared
    independently of atom-order direction.

### Corrections of scope and framing

Not every intervention was biophysics. The human also corrected the AI's engineering choices:

- **Over-fitting.** After the 1XVP work, the human required a target-agnostic rewrite: 1XVP and
  the 82 PXR/CAR structures are benchmark data, and core code must contain no chain IDs, ligand
  names or residue numbers from them. The audit is commit `56fbaa3`.
- **Honest benchmarking.** "Skip slow tasks; we'll later confirm on a local system." The
  benchmark report says exactly which steps ran on how many entries, rather than extrapolating.
- **Naming.** The import package had kept its working name, `gmxprep`, after the project was
  renamed. The human caught it, and it is now `gromatrix` (`bee70e6`).
- **Versioning.** A 1.0.0 version prepared by the AI was rejected as premature. The project stays
  at 0.x alpha until steps F–H have run on real GROMACS.
- **Tooling.** pixi (not conda/micromamba) for environments and CI; SHA-pinned actions and
  least-privilege permissions.

## Where the AI caught its own mistakes

The human asked for independent review agents (`/code-review`, `/feature-dev`). Two rounds ran,
each as several reviewers with different focuses. Every finding was re-verified before it was
fixed, and each fix has a regression test.

| Round | Example findings (all fixed) |
|---|---|
| Steps A–F (`6d6f394`) | The OpenMM chain-ID remap could swap chains in polymer-only mmCIFs. Fixed by mirroring OpenMM's own label/auth rule. The pocket-collapse gate failed correct metalloenzymes: the deposited Cu–His (1PLC) and Zn–His (1CA2) distances are 1.91 Å, and the PLP–Lys Schiff base (1ASN) is 1.37 Å. The gate is now relative to the input and ignores ions. Bonded-parameter conflicts between ligands were silently mixed. |
| Steps G–H (`8889db6`) | A phase that failed its LINCS gate could pass on resume, because the checkpoint was kept. An unconverged EM left by the smoke test skipped the EM gate. Slurm `D-HH` walltimes were mis-parsed. `--only` truncated the analysis summary. |

## Quality assurance, in numbers

These numbers come from the repository at the time of writing:

- **Code.** About 6,700 lines in `src/gromatrix/`, 1,600 lines of tests and 950 lines of agent
  skills across six skills.
- **Tests.** 82 pass in CI on Python 3.10, 3.12 and 3.13. They include:
  - an offline end-to-end regression on 1XVP: fetch fixture → fix → extract → relax;
  - a fake `gmx` that reproduces GROMACS log, checkpoint and `-maxh` behaviour, for step G;
  - a synthetic trajectory with known answers, for step H.
- **Benchmark.** 82 human PXR/CAR crystal structures (partial run; see the
  [report](../benchmarks/pxr_car_validation.md)).
- **What is not verified yet.** Steps F–H have not run against a real GROMACS / acpype / CGenFF
  installation; their logic is tested against simulated tool output. The benchmark ran steps B,
  C and E on 16 of 82 entries. Both are stated in the README and the changelog. Neither is to be
  claimed as done until the local runs exist.
