# Third-party software

gromatrix itself is MIT-licensed. It builds on the following software, each under its own licence.

**Imported as Python libraries** (permissive or LGPL; importing LGPL libraries from MIT code is allowed):
NumPy, SciPy (BSD-3-Clause); PyYAML, PDBFixer (MIT); Biopython (Biopython licence, BSD-style);
RDKit (BSD-3-Clause); matplotlib (PSF-style); OpenMM and MDAnalysis (LGPL-3.0-or-later).

**Run as separate programs, never bundled or imported** (so their licences do not apply to gromatrix):
GROMACS (LGPL-2.1-or-later); acpype and AmberTools (GPL-3.0-or-later and others); Open Babel
(GPL-2.0-only); CGenFF is accessed externally under its own terms.

Rules that keep this separation intact: do not ship these programs inside a gromatrix release, and do not
`import` acpype or openbabel from gromatrix code. Install them from conda-forge or pixi (`pixi install -e md`).
