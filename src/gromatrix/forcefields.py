"""Forcefield registry and the batch-level forcefield lock.

Every batch directory carries a ``00-config/forcefield.lock.json``. The first command that
touches a batch writes it; every later command compares its requested forcefield
against it and aborts on any difference. Ligand parameter caches and per-system
stamps record the same identity, so ``gmx_prep verify`` can prove that an entire
comparison matrix shares one forcefield.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

from . import layout

LOCK_NAME = layout.LOCK


class ForcefieldMismatch(RuntimeError):
    """Raised when a batch would mix forcefields."""


@dataclass(frozen=True)
class WaterModel:
    name: str  # value for `gmx pdb2gmx -water`
    itp: str  # topology file inside the .ff directory
    solvent_gro: str  # box for `gmx solvate -cs`


@dataclass(frozen=True)
class ForceField:
    key: str
    family: str  # "amber" | "charmm"
    gmx_name: str  # directory name without the .ff suffix
    bundled: bool  # ships with GROMACS; otherwise a port must be on GMXLIB / in the batch dir
    water: WaterModel
    ligand_engine: str  # "acpype" | "cgenff"
    ligand_ff: str  # "gaff2" | "cgenff"
    charge_method: str  # "bcc" (AM1-BCC) | "cgenff" (analogy)
    cation: str  # genion -pname
    anion: str  # genion -nname
    ion_moltypes: dict = field(default_factory=dict)  # element -> moleculetype in <ff>/ions.itp
    nonbonded: dict = field(default_factory=dict)  # mdp keys fixed by the forcefield

    def identity(self) -> dict:
        """The fields that must be identical across a batch."""
        return {
            "key": self.key,
            "family": self.family,
            "gmx_name": self.gmx_name,
            "water": self.water.name,
            "ligand_engine": self.ligand_engine,
            "ligand_ff": self.ligand_ff,
            "charge_method": self.charge_method,
        }


TIP3P = WaterModel("tip3p", "tip3p.itp", "spc216.gro")
OPC = WaterModel("opc", "opc.itp", "tip4p.gro")
TIP4PEW = WaterModel("tip4pew", "tip4pew.itp", "tip4p.gro")  # 4-site; SQUIP's AMBER choice

AMBER_IONS = {"NA": "NA", "CL": "CL", "K": "K", "MG": "MG", "CA": "CA", "ZN": "ZN"}
CHARMM_IONS = {"NA": "SOD", "CL": "CLA", "K": "POT", "MG": "MG", "CA": "CAL", "ZN": "ZN2"}

AMBER_NB = {
    "cutoff-scheme": "Verlet",
    "coulombtype": "PME",
    "rcoulomb": "1.0",
    "vdwtype": "Cut-off",
    "rvdw": "1.0",
    "DispCorr": "EnerPres",
}
# CHARMM36(m) is parameterised with a 1.0-1.2 nm force switch and no dispersion correction.
CHARMM_NB = {
    "cutoff-scheme": "Verlet",
    "coulombtype": "PME",
    "rcoulomb": "1.2",
    "vdwtype": "Cut-off",
    "vdw-modifier": "Force-switch",
    "rvdw-switch": "1.0",
    "rvdw": "1.2",
    "DispCorr": "no",
}


def _amber(key, gmx_name, bundled, water=TIP3P):
    return ForceField(
        key, "amber", gmx_name, bundled, water, "acpype", "gaff2", "bcc", "NA", "CL", AMBER_IONS, AMBER_NB
    )


def _charmm(key, gmx_name):
    return ForceField(
        key, "charmm", gmx_name, False, TIP3P, "cgenff", "cgenff", "cgenff", "SOD", "CLA", CHARMM_IONS, CHARMM_NB
    )


REGISTRY: dict[str, ForceField] = {
    ff.key: ff
    for ff in (
        _amber("amber99sb-ildn", "amber99sb-ildn", bundled=True),
        _amber("amber99sb-ildn-tip4pew", "amber99sb-ildn", bundled=True, water=TIP4PEW),
        # Community ports: place the .ff directory in the batch dir or on GMXLIB.
        _amber("amber14sb", "amber14sb_OL15", bundled=False),
        _amber("amber19sb", "amber19sb", bundled=False, water=OPC),
        _charmm("charmm36m", "charmm36-jul2022"),
    )
}

ALIASES = {"amber": "amber99sb-ildn", "charmm": "charmm36m", "charmm36": "charmm36m"}


def resolve(name: str) -> ForceField:
    key = ALIASES.get(name.lower(), name.lower())
    if key not in REGISTRY:
        raise KeyError(f"unknown forcefield '{name}'; choose from {sorted(REGISTRY) + sorted(ALIASES)}")
    return REGISTRY[key]


def lock_path(batch_dir: Path) -> Path:
    return Path(batch_dir) / LOCK_NAME


def load_lock(batch_dir: Path) -> ForceField:
    path = lock_path(batch_dir)
    if not path.exists():
        raise ForcefieldMismatch(f"{path} missing; run `gmx_prep lock` first")
    data = json.loads(path.read_text())
    ff = resolve(data["key"])
    if ff.identity() != data:
        raise ForcefieldMismatch(f"{path} does not match registry entry for {ff.key}: {data}")
    return ff


def acquire_lock(batch_dir: Path, requested: str | None) -> ForceField:
    """Create the lock, or check `requested` against an existing one.

    `requested=None` means "use whatever the batch is locked to".
    """
    path = lock_path(batch_dir)
    if path.exists():
        locked = load_lock(batch_dir)
        if requested is not None and resolve(requested).identity() != locked.identity():
            raise ForcefieldMismatch(
                f"batch {batch_dir} is locked to '{locked.key}' but '{requested}' was requested; "
                "forcefields must not be mixed within a comparison matrix"
            )
        return locked
    if requested is None:
        raise ForcefieldMismatch(f"{path} missing and no --forcefield given")
    ff = resolve(requested)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(ff.identity(), indent=2) + "\n")
    return ff


def assert_same(ff: ForceField, stamp: dict, where: str) -> None:
    if stamp != ff.identity():
        raise ForcefieldMismatch(f"{where} was built with {stamp.get('key')} but batch is locked to {ff.key}")


def describe(ff: ForceField) -> dict:
    d = asdict(ff)
    d["water"] = asdict(ff.water)
    return d
