"""Minimal PDB / GRO I/O and residue classification (no heavy dependencies)."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

PROTEIN_RES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "HID",
    "HIE",
    "HIP",
    "HSD",
    "HSE",
    "HSP",
    "CYX",
    "CYM",
    "ASH",
    "GLH",
    "LYN",
    "ACE",
    "NME",
    "NH2",
}
WATER_RES = {"HOH", "WAT", "SOL", "TIP3", "DOD", "H2O"}
# Residue name -> element for monatomic ions.
ION_RES = {
    "NA": "NA",
    "CL": "CL",
    "K": "K",
    "MG": "MG",
    "CA": "CA",
    "ZN": "ZN",
    "SOD": "NA",
    "CLA": "CL",
    "POT": "K",
    "CAL": "CA",
    "ZN2": "ZN",
    "MN": "MN",
    "FE": "FE",
    "FE2": "FE",
    "CU": "CU",
    "CU1": "CU",
    "CO": "CO",
    "NI": "NI",
    "CD": "CD",
}


@dataclass(frozen=True)
class Atom:
    record: str
    name: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    element: str
    b: float = 0.0
    occ: float = 1.0

    @property
    def res_key(self):
        return (self.chain, self.resseq, self.icode, self.resname)

    def pdb_line(self, serial: int) -> str:
        """Fixed-column PDB ATOM/HETATM record. A 4-character residue name uses
        column 21 (the GROMACS/CHARMM convention, which read_pdb understands);
        anything that cannot fit its columns raises instead of shifting them."""
        if len(self.name) > 4 or len(self.resname) > 4 or len(self.chain) > 1 or len(self.icode) > 1:
            raise ValueError(f"{self.resname}{self.resseq}{self.chain}:{self.name} does not fit PDB columns")
        if not -999 <= self.resseq <= 9999:
            raise ValueError(
                f"residue number {self.resseq} does not fit PDB columns 23-26; use mmCIF input or renumber the chain"
            )
        name = self.name if len(self.name) == 4 else f" {self.name:<3}"
        resname = f"{self.resname:<4}" if len(self.resname) == 4 else f"{self.resname:>3} "
        return (
            f"{self.record:<6}{serial % 100000:5d} {name:<4} {resname}{self.chain:1}"
            f"{self.resseq:4d}{self.icode:1}   {self.x:8.3f}{self.y:8.3f}{self.z:8.3f}"
            f"{self.occ:6.2f}{self.b:6.2f}          {self.element:>2}"
        )


def read_pdb(path: Path) -> list[Atom]:
    atoms, model_seen = [], False
    for line in Path(path).read_text().splitlines():
        if line.startswith("MODEL"):
            if model_seen:
                break  # first model only
            model_seen = True
        if line.startswith("ENDMDL") and model_seen:
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[16] not in (" ", "A"):
            continue  # keep altloc A only
        name = line[12:16].strip()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = "".join(c for c in name if c.isalpha())[:1].upper()
        occ = float(line[54:60]) if line[54:60].strip() else 1.0
        b = float(line[60:66]) if line[60:66].strip() else 0.0
        atoms.append(
            Atom(
                line[:6].strip(),
                name,
                line[17:21].strip(),
                line[21].strip(),
                int(line[22:26]),
                line[26].strip(),
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
                element,
                b,
                occ,
            )
        )
    return atoms


PEPTIDE_BREAK = 2.0  # Å; C(i)-N(i+1) longer than this is a chain break


def _linked(prev: list[Atom], cur: list[Atom]) -> bool:
    """True if `cur` continues the polymer chain of `prev` (same chain, bonded)."""
    p, c = prev[0], cur[0]
    if p.chain != c.chain or p.resname not in PROTEIN_RES or c.resname not in PROTEIN_RES:
        return False
    cp = next((a for a in prev if a.name == "C"), None)
    nc = next((a for a in cur if a.name == "N"), None)
    if cp is None or nc is None:
        return False
    return ((cp.x - nc.x) ** 2 + (cp.y - nc.y) ** 2 + (cp.z - nc.z) ** 2) ** 0.5 <= PEPTIDE_BREAK


def write_pdb(path: Path, atoms: list[Atom], ter: bool = True) -> None:
    """Write atoms; with `ter`, emit TER at chain ends and backbone breaks so
    pdb2gmx (-chainsep id_or_ter) and OpenMM never bond across a gap."""
    lines, serial, prev = [], 0, None
    for res in residues(atoms):
        if ter and prev is not None and prev[0].resname in PROTEIN_RES and not _linked(prev, res):
            lines.append("TER")
        for a in res:
            serial += 1
            lines.append(a.pdb_line(serial))
        prev = res
    if ter and prev is not None and prev[0].resname in PROTEIN_RES:
        lines.append("TER")
    Path(path).write_text("\n".join(lines + ["END"]) + "\n")


def residues(atoms: list[Atom]) -> list[list[Atom]]:
    out, cur, key = [], [], None
    for a in atoms:
        if a.res_key != key:
            if cur:
                out.append(cur)
            cur, key = [], a.res_key
        cur.append(a)
    if cur:
        out.append(cur)
    return out


def is_hydrogen(a: Atom) -> bool:
    return a.element in ("H", "D")


@dataclass
class Split:
    protein: list[Atom]
    ions: list[Atom]  # one atom per ion
    ligands: list[list[Atom]]  # one heavy-atom list per residue instance


def classify(atoms: list[Atom], protein_chains: set[str] | None = None) -> Split:
    """Partition atoms. Water and all hydrogens are dropped (pdb2gmx and the
    ligand router re-add hydrogens consistently with the forcefield)."""
    split = Split([], [], [])
    for res in residues(atoms):
        r = res[0]
        if r.resname in WATER_RES:
            continue
        heavy = [a for a in res if not is_hydrogen(a)]
        if r.resname in PROTEIN_RES:
            # Residues of unselected protein chains are dropped. They are never ligands.
            if protein_chains is None or r.chain in protein_chains:
                split.protein.extend(heavy)
        elif r.resname in ION_RES and len(heavy) == 1:
            split.ions.append(heavy[0])
        elif heavy:
            split.ligands.append(heavy)
    return split


def residue_block(res: list[Atom]) -> str:
    return "\n".join(a.pdb_line(i + 1) for i, a in enumerate(res)) + "\nEND\n"


# --------------------------------------------------------------------- GRO
@dataclass(frozen=True)
class GroAtom:
    resnr: int
    resname: str
    name: str
    x: float  # nm
    y: float
    z: float


def read_gro(path: Path) -> tuple[list[GroAtom], str]:
    lines = Path(path).read_text().splitlines()
    n = int(lines[1])
    atoms = [
        GroAtom(int(l[0:5]), l[5:10].strip(), l[10:15].strip(), float(l[20:28]), float(l[28:36]), float(l[36:44]))
        for l in lines[2 : 2 + n]
    ]
    return atoms, lines[2 + n]


def write_gro(path: Path, atoms: list[GroAtom], title: str, box: str = "   0.00000   0.00000   0.00000") -> None:
    out = [title, f"{len(atoms):5d}"]
    for i, a in enumerate(atoms, 1):
        out.append(
            f"{a.resnr % 100000:5d}{a.resname[:5]:<5}{a.name[:5]:>5}{i % 100000:5d}{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}"
        )
    out.append(box)
    Path(path).write_text("\n".join(out) + "\n")


def renumber(atoms: list[GroAtom], start: int) -> list[GroAtom]:
    """Shift residue numbers so the block starts at `start`; returns new atoms."""
    if not atoms:
        return []
    off = start - atoms[0].resnr
    return [replace(a, resnr=a.resnr + off) for a in atoms]
