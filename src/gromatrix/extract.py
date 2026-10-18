"""[C] Sub-unit extraction with spatial ligand mapping and completeness scoring.

    gromatrix-extract --batch-dir B [--in fixed] [--out units] [--ids 1ABC ...]
                    [--targets LIG[,RES][@ID]] [--target-chains X[,Y][@ID]] [--expect RES[:n],...]
                    [--assembly auto|none|N] [--ligand-cutoff 4.5] [--contact-cutoff 5.0]
                    [--min-contacts 10] [--max-chains 12] [--unit CANDIDATE_ID] [--allow-incomplete]

Algorithm (per structure):
  1. Candidates. Every biological assembly in the deposit is rebuilt by applying
     its operators. Copies made by non-identity operators get fresh chain IDs,
     so there are no collisions. Independently, the asymmetric unit is cut into
     stoichiometric units: entity copy numbers divided by their GCD, keeping
     only combinations whose chains form one connected interface graph.
  2. Ligand mapping. Inside each candidate, every non-polymer residue (ligand,
     cofactor, ion, kept water) is assigned to the polymer chain with the most
     heavy-atom pairs within --ligand-cutoff. Residues touching no chain are orphans.
  3. Completeness. Hard requirement: every --targets / --expect residue is present.
     Soft score: fraction of all HET types in the deposit that the unit carries.
  4. Ranking: complete > coverage > author-defined > fewer chains >
     larger interfaces > fewer missing residues > lower mean Cα B-factor.
     So when copies in the asymmetric unit differ in which ligands they
     hold, the copy that holds the declared targets wins, whatever its order
     in the file. HET membership in assemblies follows the chain each HET
     physically binds, not its deposited chain label.
  5. Matrix slicing. Writes enzyme.pdb (+ cofactors), target_<T>.pdb and
     complex_<T>.pdb with author numbering, TER records at breaks, and
     unit.json documenting every candidate and every ligand assignment.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import string
import sys
import traceback
from collections import Counter, OrderedDict
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from . import altlocs, layout, result
from .pdbio import ION_RES, WATER_RES, Atom, read_pdb, residues, write_pdb
from .structure import POLYMER_RES, compose, entities_by_sequence, is_identity

CHAIN_POOL = string.ascii_uppercase + string.digits + string.ascii_lowercase
MAX_COMBOS = 20000


class ExtractionError(result.CodedError):
    code = "extract_failed"


@dataclass
class Res:
    chain: str  # chain ID inside the candidate (after collision renaming)
    origin: str  # author chain in the deposit
    oper: str  # operator label; "1" = identity / asymmetric unit
    resname: str
    resseq: int
    icode: str
    kind: str  # polymer | ligand | ion | water
    atoms: list[Atom]
    xyz: np.ndarray
    R: np.ndarray | None = None  # operator that placed this copy (None = identity)
    t: np.ndarray | None = None

    @property
    def key(self):
        return (self.origin, self.resname, self.resseq, self.icode)

    def label(self) -> str:
        return f"{self.resname}{self.resseq}{self.icode}:{self.chain}"


@dataclass
class Candidate:
    id: str
    source: str  # assembly | au-unit | au
    author_defined: bool
    residues: list[Res]
    assignment: dict[int, dict] = field(default_factory=dict)  # residue index -> {chain, contacts, min_dist}
    orphans: list[int] = field(default_factory=list)
    metrics: dict = field(default_factory=dict)

    @property
    def chains(self) -> list[str]:
        return list(OrderedDict.fromkeys(r.chain for r in self.residues if r.kind == "polymer"))


def _kind(resname: str, n_heavy: int) -> str:
    if resname in POLYMER_RES:
        return "polymer"
    if resname in WATER_RES:
        return "water"
    if resname in ION_RES and n_heavy == 1:
        return "ion"
    return "ligand"


def load_au(path: Path) -> list[Res]:
    out = []
    for res in residues([a for a in read_pdb(path) if a.element not in ("H", "D")]):
        a0 = res[0]
        out.append(
            Res(
                a0.chain,
                a0.chain,
                "1",
                a0.resname,
                a0.resseq,
                a0.icode,
                _kind(a0.resname, len(res)),
                res,
                np.array([[a.x, a.y, a.z] for a in res]),
            )
        )
    return out


# ------------------------------------------------------------------ geometry
def chain_contacts(res: list[Res], cutoff: float) -> Counter:
    """(chainA, chainB) -> number of residue pairs with heavy atoms within `cutoff`."""
    poly = [(i, r) for i, r in enumerate(res) if r.kind == "polymer"]
    if not poly:
        return Counter()
    xyz = np.concatenate([r.xyz for _, r in poly])
    owner = np.concatenate([np.full(len(r.atoms), i) for i, r in poly])
    pairs = cKDTree(xyz).query_pairs(cutoff, output_type="ndarray")
    if len(pairs) == 0:
        return Counter()
    ri, rj = owner[pairs[:, 0]], owner[pairs[:, 1]]
    seen, out = set(), Counter()
    for a, b in zip(ri, rj):
        ca, cb = res[a].chain, res[b].chain
        if ca == cb or (a, b) in seen:
            continue
        seen.add((a, b))
        seen.add((b, a))
        out[tuple(sorted((ca, cb)))] += 1
    return out


def assign_hets(res: list[Res], cutoff: float) -> tuple[dict[int, dict], list[int]]:
    """Map each non-polymer residue to the polymer chain it contacts most."""
    poly = [r for r in res if r.kind == "polymer"]
    if not poly:
        return {}, [i for i, r in enumerate(res) if r.kind != "polymer"]
    xyz = np.concatenate([r.xyz for r in poly])
    chain_of = np.concatenate([np.full(len(r.atoms), r.chain, dtype=object) for r in poly])
    tree = cKDTree(xyz)
    assignment, orphans = {}, []
    for i, r in enumerate(res):
        if r.kind == "polymer":
            continue
        hits, mind = Counter(), {}
        for p, nbs in zip(r.xyz, tree.query_ball_point(r.xyz, cutoff)):
            if not nbs:
                continue
            for c, d in zip(chain_of[nbs], np.linalg.norm(xyz[nbs] - p, axis=1)):
                hits[c] += 1
                mind[c] = min(mind.get(c, cutoff), float(d))
        if not hits:
            orphans.append(i)
            continue
        best = max(hits, key=lambda c: (hits[c], -mind[c]))
        assignment[i] = {
            "chain": best,
            "contacts": int(hits[best]),
            "min_dist": round(mind[best], 2),
            "other_chains": {c: int(n) for c, n in hits.items() if c != best},
        }
    return assignment, orphans


# ------------------------------------------------------------------ candidates
def build_assembly(asm: dict, au: list[Res], ops: dict, het_chain: dict | None = None) -> Candidate:
    """Apply assembly operators. Identity copies keep their author chain IDs;
    every other copy gets an unused ID, so chain IDs never collide.

    A HET joins an operator copy when the chain it physically binds (`het_chain`,
    from the AU-wide spatial assignment) is in that copy. Deposited chain labels
    for HETs are not trusted, because a ligand labelled with one chain can sit in another
    chain's pocket. The deposit's own membership list is used only for HETs that touch no chain."""
    het_chain = het_chain or {}
    plan = []  # (residue, op_label, matrix)
    for gen in asm["gens"]:
        poly = set(gen["polymer_chains"])
        hets = {tuple(h) for h in gen["het_residues"]}
        for seq in gen["opers"]:
            op = compose(ops, seq)
            label = "x".join(seq)
            for r in au:
                if r.kind == "polymer":
                    member = r.origin in poly
                elif r.key in het_chain:
                    member = het_chain[r.key] in poly
                else:
                    member = r.key in hets
                if member:
                    plan.append((r, label, op))
    ids: dict[tuple[str, str], str] = {}
    for r, label, op in plan:  # pass 1: identity copies keep author IDs
        if is_identity(op) and (r.origin, label) not in ids and r.origin not in ids.values():
            ids[(r.origin, label)] = r.origin
    reserved = {r.chain for r in au}
    for r, label, _op in plan:  # pass 2: everything else gets a fresh ID
        if (r.origin, label) not in ids:
            used = set(ids.values()) | reserved
            free = next((c for c in CHAIN_POOL if c not in used), None)
            if free is None:
                raise ExtractionError(f"assembly {asm['id']} needs more than {len(CHAIN_POOL)} chain IDs")
            ids[(r.origin, label)] = free
    out = []
    for r, label, op in plan:
        cid = ids[(r.origin, label)]
        xyz = r.xyz @ np.array(op["R"]).T + np.array(op["t"])
        atoms = [replace(a, chain=cid, x=float(p[0]), y=float(p[1]), z=float(p[2])) for a, p in zip(r.atoms, xyz)]
        out.append(replace(r, chain=cid, oper=label, atoms=atoms, xyz=xyz, R=np.array(op["R"]), t=np.array(op["t"])))
    return Candidate(f"assembly-{asm['id']}", "assembly", bool(asm.get("author_defined")), out)


def au_units(au: list[Res], chain_entity: dict, contacts: Counter, min_contacts: int) -> list[tuple[str, ...]]:
    """Stoichiometric, interface-connected chain sets cut from the asymmetric unit."""
    chains = list(OrderedDict.fromkeys(r.chain for r in au if r.kind == "polymer"))
    fallback = entities_by_sequence([a for r in au for a in r.atoms])
    by_ent: dict[str, list[str]] = OrderedDict()
    for c in chains:
        by_ent.setdefault(chain_entity.get(c) or f"seq{fallback.get(c)}", []).append(c)
    g = math.gcd(*[len(v) for v in by_ent.values()])
    if g <= 1:
        return []
    need = {e: len(v) // g for e, v in by_ent.items()}
    per_entity = [list(itertools.combinations(v, need[e])) for e, v in by_ent.items()]
    if math.prod(len(p) for p in per_entity) > MAX_COMBOS:
        return _greedy_units(by_ent, need, contacts, min_contacts)
    units = []
    for combo in itertools.product(*per_entity):
        members = tuple(c for part in combo for c in part)
        if _connected(members, contacts, min_contacts):
            units.append(members)
    return units


def _connected(members, contacts, min_contacts) -> bool:
    members = list(members)
    seen, stack = {members[0]}, [members[0]]
    while stack:
        c = stack.pop()
        for d in members:
            if d not in seen and contacts.get(tuple(sorted((c, d))), 0) >= min_contacts:
                seen.add(d)
                stack.append(d)
    return len(seen) == len(members)


def _greedy_units(by_ent, need, contacts, min_contacts):
    seed_ent = min(by_ent, key=lambda e: len(by_ent[e]))
    units = []
    for seed in by_ent[seed_ent]:
        unit = [seed]
        for e, chains in by_ent.items():
            pool = [c for c in chains if c not in unit]
            for _ in range(need[e] - (1 if e == seed_ent else 0)):
                best = max(pool, key=lambda c: sum(contacts.get(tuple(sorted((c, u))), 0) for u in unit))
                unit.append(best)
                pool.remove(best)
        if _connected(unit, contacts, min_contacts):
            units.append(tuple(unit))
    return units


def subset(au: list[Res], members: tuple[str, ...], assignment: dict, name: str, source: str) -> Candidate:
    keep = []
    for i, r in enumerate(au):
        if r.kind == "polymer" and r.chain in members:
            keep.append(r)
        elif r.kind != "polymer" and i in assignment and assignment[i]["chain"] in members:
            keep.append(r)
    return Candidate(name, source, False, keep)


# ------------------------------------------------------------------ scoring
def score(c: Candidate, req: dict, target_chains: set[str], all_types: set[str], gaps: list[dict], a) -> None:
    c.assignment, c.orphans = assign_hets(c.residues, a.ligand_cutoff)
    inv = Counter(c.residues[i].resname for i in c.assignment if c.residues[i].kind in ("ligand", "ion"))
    origins = {r.origin for r in c.residues if r.kind == "polymer"}
    missing = {k: n - inv.get(k, 0) for k, n in req.items() if inv.get(k, 0) < n}
    missing_chains = sorted(target_chains - origins)
    contacts = chain_contacts(c.residues, a.contact_cutoff)
    ca_b = [at.b for r in c.residues if r.kind == "polymer" for at in r.atoms if at.name == "CA" and at.occ > 0]
    open_res = sum(g["n"] for g in gaps if not g["modelled"] and g["chain"] in origins)
    coverage = len(set(inv) & all_types) / len(all_types) if all_types else 1.0
    c.metrics = {
        "complete": not missing and not missing_chains,
        "missing": missing,
        "missing_target_chains": missing_chains,
        "coverage": round(coverage, 3),
        "inventory": dict(inv),
        "n_chains": len(c.chains),
        "interface_contacts": int(sum(contacts.values())),
        "unmodelled_residues": open_res,
        "mean_ca_b": round(float(np.mean(ca_b)), 1) if ca_b else None,
        "orphans": [c.residues[i].label() for i in c.orphans if c.residues[i].kind != "water"],
    }


def rank_key(c: Candidate):
    m = c.metrics
    return (
        not m["complete"],
        -m["coverage"],
        not c.author_defined,
        m["n_chains"],
        -m["interface_contacts"],
        m["unmodelled_residues"],
        m["mean_ca_b"] or 0.0,
    )


# ------------------------------------------------------------------ slicing
def _relabel_het(c: Candidate, i: int, used: set) -> Res:
    """Give a cofactor the chain ID of the chain it binds, unless that would
    duplicate an existing (chain, resseq, icode)."""
    r = c.residues[i]
    target = c.assignment[i]["chain"]
    if target == r.chain or (target, r.resseq, r.icode) in used:
        return r
    used.add((target, r.resseq, r.icode))
    return replace(r, chain=target, atoms=[replace(at, chain=target) for at in r.atoms])


def with_altloc(r: Res, groups: dict, letter: str) -> Res:
    """Copy of `r` with conformer `letter` where the deposit has it; shared atoms
    and residues without that letter are unchanged (wwPDB: one letter, one state)."""
    g = groups.get(r.key)
    coords = altlocs.alt_coords(g, letter) if g else None
    if not coords:
        return r
    atoms = []
    for a in r.atoms:
        if a.name in coords:
            x = np.array(coords[a.name])
            if r.R is not None:
                x = r.R @ x + r.t
            a = replace(a, x=float(x[0]), y=float(x[1]), z=float(x[2]))
        atoms.append(a)
    return replace(r, atoms=atoms, xyz=np.array([[a.x, a.y, a.z] for a in atoms]))


def clash_report(lig: Res, others: list[Res], cutoff: float = 2.2, contact: float = 4.5) -> dict:
    """Clashes (< `cutoff` Å) and contacts (< `contact` Å) between a ligand pose and its surroundings."""
    env = [r for r in others if r is not lig and r.kind != "water"]
    if not env:
        return {"min_dist": None, "clashes": 0, "contacts": 0}
    tree = cKDTree(np.concatenate([r.xyz for r in env]))
    d, _ = tree.query(lig.xyz)
    return {
        "min_dist": round(float(d.min()), 2),
        "clashes": int((d < cutoff).sum()),
        "contacts": int((d < contact).sum()),
    }


def plan_altlocs(c: Candidate, targets: list[str], groups: dict, mode: str) -> dict[str, dict]:
    """Per target residue name: branch over altloc letters or use the primary.
    auto   - branch only if some instance has a distinct, geometrically consistent conformer
    branch - branch whenever an instance has >= 2 usable conformers
    primary- never branch"""
    plan = {}
    for t in targets:
        keys = [r.key for r in c.residues if r.resname == t]
        gs = [groups[k] for k in keys if k in groups]
        letters = altlocs.letters_for(groups, keys)
        if not gs:
            plan[t] = {"branch": False, "reason": "single conformer in deposit"}
        elif mode == "primary":
            plan[t] = {"branch": False, "reason": "--altloc primary", "letters": letters}
        elif mode == "auto" and not any(g["distinct"] for g in gs):
            plan[t] = {
                "branch": False,
                "letters": letters,
                "reason": "conformers equivalent within --altloc-rmsd (symmetry-corrected); primary used",
                "sym_rmsd": {
                    f"{g['resname']}{g['resseq']}{g['chain']}": {L: r["sym_rmsd"] for L, r in g["analysis"].items()}
                    for g in gs
                },
            }
        elif len(letters) < 2:
            plan[t] = {
                "branch": False,
                "reason": "fewer than two geometrically consistent conformers",
                "letters": letters,
            }
        else:
            plan[t] = {"branch": True, "letters": letters, "reason": f"--altloc {mode}"}
    return plan


def slice_matrix(
    c: Candidate,
    targets: list[str],
    target_chains: set[str],
    outdir: Path,
    groups: dict | None = None,
    alt_plan: dict | None = None,
) -> dict:
    groups, alt_plan = groups or {}, alt_plan or {}
    used = {(r.chain, r.resseq, r.icode) for r in c.residues}
    tchains = {r.chain for r in c.residues if r.kind == "polymer" and r.origin in target_chains}
    poly = [r for r in c.residues if r.kind == "polymer"]
    hets = [(i, _relabel_het(c, i, used)) for i in c.assignment]

    enzyme_poly = [r for r in poly if r.chain not in tchains]
    enzyme_het = [r for i, r in hets if r.resname not in targets and c.assignment[i]["chain"] not in tchains]
    order = {ch: n for n, ch in enumerate(c.chains)}

    def write(name, res_list):
        res_list = sorted(res_list, key=lambda r: (r.kind != "polymer", order.get(r.chain, 99)))
        write_pdb(outdir / name, [a for r in res_list for a in r.atoms])
        return name

    files = {
        "unit": write("unit.pdb", poly + [r for _, r in hets]),
        "enzyme": write("enzyme.pdb", enzyme_poly + enzyme_het),
    }
    for t in targets:
        inst = [(i, r) for i, r in hets if r.resname == t]
        if not inst:
            continue
        lead_i, lead = max(inst, key=lambda ir: c.assignment[ir[0]]["contacts"])
        plan = alt_plan.get(t, {})
        if not plan.get("branch"):
            files[f"target_{t}"] = write(f"target_{t}.pdb", [lead])
            files[f"complex_{t}"] = write(f"complex_{t}.pdb", enzyme_poly + enzyme_het + [r for _, r in inst])
            continue
        plan["branches"] = {}
        for L in plan["letters"]:
            env = [with_altloc(r, groups, L) for r in enzyme_poly + enzyme_het]
            ligs = [with_altloc(r, groups, L) for _, r in inst]
            # Every instance that changes under letter L is checked, not just the lead copy.
            checks = [(orig, new) for (_, orig), new in zip(inst, ligs) if new is not orig] or [(lead, ligs[0])]
            reports = [clash_report(new, env + [x for x in ligs if x is not new]) for _, new in checks]
            plan["branches"][L] = {"clash": min(reports, key=lambda r: (r["min_dist"] is None, r["min_dist"] or 0))}
            if any(r["min_dist"] is not None and r["min_dist"] < 1.5 for r in reports):
                plan["branches"][L]["skipped"] = "severe clash with environment (< 1.5 Å)"
                continue
            if any(r["contacts"] == 0 for r in reports):
                # Some deposits store a second copy of the ligand, bound at another site,
                # as an altloc of the same residue. That is not a conformer of this pose.
                plan["branches"][L]["skipped"] = (
                    "pose does not touch this unit: an alternate copy bound elsewhere in the deposit, not a conformer"
                )
                continue
            lead_L = ligs[[r for _, r in inst].index(lead)]
            files[f"target_{t}_alt{L}"] = write(f"target_{t}_alt{L}.pdb", [lead_L])
            files[f"complex_{t}_alt{L}"] = write(f"complex_{t}_alt{L}.pdb", env + ligs)
        if not any(k.startswith(f"complex_{t}_alt") for k in files):
            # Every branch was rejected: never drop the target from the matrix; use the primary pose.
            plan["fallback"] = "all branches rejected; primary conformer used"
            files[f"target_{t}"] = write(f"target_{t}.pdb", [lead])
            files[f"complex_{t}"] = write(f"complex_{t}.pdb", enzyme_poly + enzyme_het + [r for _, r in inst])
    if tchains:
        tag = "".join(sorted(tchains))
        t_res = [r for r in poly if r.chain in tchains] + [
            r for i, r in hets if c.assignment[i]["chain"] in tchains and r.resname not in targets
        ]
        files[f"target_{tag}"] = write(f"target_{tag}.pdb", t_res)
        files[f"complex_{tag}"] = write(f"complex_{tag}.pdb", enzyme_poly + enzyme_het + t_res)
    return files


# ------------------------------------------------------------------ driver
def _scoped(spec: str | None, pid: str) -> list[str]:
    """'LIG,COF@1ABC,XYZ@2ABC' -> names that apply to `pid`."""
    out = []
    for item in (spec or "").split(","):
        item = item.strip()
        if not item:
            continue
        name, _, scope = item.partition("@")
        if not scope or scope.upper() == pid:
            out.append(name.upper() if len(name) > 1 else name)
    return out


def extract_one(pdb: Path, out_root: Path, a) -> dict:
    pid = pdb.stem.upper()
    meta_p = pdb.with_suffix(".meta.json")
    meta = (
        json.loads(meta_p.read_text()) if meta_p.exists() else {"assemblies": [], "operators": {}, "chain_entity": {}}
    )
    gaps = meta.get("fix", {}).get("gaps", [])
    warnings: list[str] = []
    if not meta_p.exists():
        warnings.append("no .meta.json (run gromatrix-fix first): assemblies unknown, using AU units only")

    au = load_au(pdb)
    targets = _scoped(a.targets, pid)
    target_chains = set(_scoped(a.target_chains, pid))
    req: dict[str, int] = {t: 1 for t in targets}
    for item in _scoped(a.expect, pid):
        name, _, n = item.partition(":")
        req[name.upper()] = max(req.get(name.upper(), 0), int(n or 1))
    # Soft coverage counts ligand types only: a stray crystallisation ion on one copy
    # must not make a larger unit outrank a complete smaller one.
    all_types = {r.resname for r in au if r.kind == "ligand"}

    # Candidates
    cands: list[Candidate] = []
    au_assign, _ = assign_hets(au, a.ligand_cutoff)
    het_chain = {au[i].key: v["chain"] for i, v in au_assign.items()}
    relabelled = [
        f"{au[i].label()} binds {v['chain']}"
        for i, v in au_assign.items()
        if au[i].kind in ("ligand", "ion") and v["chain"] != au[i].origin
    ]
    if relabelled:
        warnings.append(f"HETs whose deposited chain label differs from the chain they bind: {relabelled}")
    if a.assembly != "none":
        for asm in meta.get("assemblies", []):
            if a.assembly in ("auto", asm["id"]):
                cands.append(build_assembly(asm, au, meta["operators"], het_chain))
    if a.assembly in ("auto", "none"):
        contacts = chain_contacts(au, a.contact_cutoff)
        for members in au_units(au, meta.get("chain_entity", {}), contacts, a.min_contacts):
            cands.append(subset(au, members, au_assign, "au-" + "".join(members), "au-unit"))
        cands.append(Candidate("au", "au", False, list(au)))
    if not cands:
        raise ExtractionError(f"no candidate units (assembly '{a.assembly}' not found?)")

    # De-duplicate identical chain sets (an assembly that equals an AU unit keeps the assembly label).
    uniq, seen = [], set()
    for c in cands:
        sig = frozenset((r.origin, r.oper if r.oper != "1" else "1") for r in c.residues if r.kind == "polymer")
        if sig not in seen and len(c.chains) <= a.max_chains:
            seen.add(sig)
            uniq.append(c)
    if not uniq:
        raise ExtractionError(f"every candidate exceeds --max-chains {a.max_chains}")
    for c in uniq:
        score(c, req, target_chains, all_types, gaps, a)
    ranked = sorted(uniq, key=rank_key)

    table = [
        {
            "id": c.id,
            "source": c.source,
            "author_defined": c.author_defined,
            "chains": c.chains,
            **{
                k: c.metrics[k]
                for k in (
                    "complete",
                    "missing",
                    "coverage",
                    "inventory",
                    "n_chains",
                    "interface_contacts",
                    "unmodelled_residues",
                    "mean_ca_b",
                )
            },
        }
        for c in ranked
    ]
    outdir = out_root / pid
    outdir.mkdir(parents=True, exist_ok=True)
    for stale in [*outdir.glob("*.pdb"), outdir / "unit.json"]:  # never mix files from an earlier run
        stale.unlink(missing_ok=True)
    chosen = next((c for c in ranked if c.id == a.unit), None) if a.unit else ranked[0]
    if chosen is None:
        raise ExtractionError(f"--unit {a.unit} not among candidates {[c.id for c in ranked]}")
    m = chosen.metrics
    if not m["complete"]:
        msg = (
            f"no candidate unit contains all required ligands/chains; best '{chosen.id}' lacks "
            f"{m['missing'] or ''}{m['missing_target_chains'] or ''}"
        )
        if not a.allow_incomplete:
            (outdir / "unit.json").write_text(json.dumps({"id": pid, "error": msg, "candidates": table}, indent=1))
            raise ExtractionError(msg + " (see unit.json; --allow-incomplete to proceed)")
        warnings.append(msg)
    if m["coverage"] < 1.0:
        lacking = sorted(all_types - set(m["inventory"]))
        warnings.append(f"selected unit lacks HET types present elsewhere in the deposit: {lacking}")
    if chosen.source == "au" and len(ranked) > 1:
        warnings.append("selected the whole asymmetric unit: no smaller unit was complete")
    if m["orphans"]:
        warnings.append(f"dropped HET groups touching no chain within {a.ligand_cutoff} Å: {m['orphans']}")
    if not targets and not target_chains:
        warnings.append("no --targets/--target-chains: wrote unit.pdb and enzyme.pdb only (all HETs as cofactors)")

    groups = altlocs.index(meta.get("altlocs", []))
    alt_plan = plan_altlocs(chosen, targets, groups, a.altloc)
    files = slice_matrix(chosen, targets, target_chains, outdir, groups, alt_plan)
    for t, p in alt_plan.items():
        if p.get("branch"):
            made = [L for L, b in p["branches"].items() if "skipped" not in b]
            for L, b in p["branches"].items():
                if "skipped" in b:
                    warnings.append(f"{t}: altloc {L} not used - {b['skipped']}")
            if len(made) > 1:
                warnings.append(
                    f"{t}: matrix branched over conformers {made} ({p['reason']}) into "
                    f"complex_{t}_alt<L>.pdb - simulate each as its own complex"
                )
            elif made:
                warnings.append(
                    f"{t}: only altloc {made[0]} is a valid pose here; written as complex_{t}_alt{made[0]}.pdb"
                )
        elif p.get("letters"):
            warnings.append(f"{t}: {p['reason']}")

    report = {
        "id": pid,
        "selected": chosen.id,
        "source": chosen.source,
        "author_defined": chosen.author_defined,
        "chains": chosen.chains,
        "chain_origin": {r.chain: f"{r.origin}/op{r.oper}" for r in chosen.residues if r.kind == "polymer"},
        "targets": targets,
        "target_chains": sorted(target_chains),
        "required": req,
        "metrics": m,
        "ligand_assignment": [
            {"residue": chosen.residues[i].label(), **v}
            for i, v in chosen.assignment.items()
            if chosen.residues[i].kind != "water"
        ],
        "candidates": table,
        "altlocs": alt_plan,
        "files": files,
        "warnings": warnings,
    }
    (outdir / "unit.json").write_text(json.dumps(report, indent=1))
    return {
        "id": pid,
        "selected": chosen.id,
        "chains": chosen.chains,
        "inventory": m["inventory"],
        "runner_up": ranked[1].id if len(ranked) > 1 and ranked[1] is not chosen else None,
        "files": sorted(files.values()),
        "warnings": warnings,
    }


@result.guarded
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gromatrix-extract", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--in", dest="indir", default=layout.FIXED)
    ap.add_argument("--out", default=layout.UNITS)
    ap.add_argument("--ids", nargs="+")
    ap.add_argument("--targets", help="target ligand residue names, e.g. LIG or LIG@1ABC,XYZ@2ABC")
    ap.add_argument("--target-chains", help="author chains of a protein/peptide target, e.g. E@1ABC")
    ap.add_argument("--expect", help="additional required HETs, RES[:count], e.g. COF:1,HEM:2")
    ap.add_argument("--assembly", default="auto", help="auto | none (AU-derived units only) | assembly id")
    ap.add_argument("--ligand-cutoff", "--het-cutoff", type=float, default=4.5, help="Å, HET-to-chain mapping")
    ap.add_argument("--contact-cutoff", type=float, default=5.0, help="Å, chain-chain interface residues")
    ap.add_argument("--min-contacts", type=int, default=10, help="residue pairs for a chain-chain edge")
    ap.add_argument("--max-chains", type=int, default=12)
    ap.add_argument("--unit", help="force a candidate id from unit.json, e.g. assembly-1 or au-ABEF")
    ap.add_argument("--allow-incomplete", action="store_true")
    ap.add_argument(
        "--altloc",
        default="auto",
        choices=["auto", "branch", "primary"],
        help="target ligands with alternate conformers: branch only distinct ones (auto), "
        "always branch, or always use the primary conformer",
    )
    result.add_full_flag(ap)
    a = ap.parse_args(argv)

    indir, outdir = a.batch_dir / a.indir, a.batch_dir / a.out
    files = sorted(indir.glob("*.pdb"))
    if a.ids:
        want = {i.upper() for i in a.ids}
        files = [p for p in files if p.stem.upper() in want]
    if not files:
        return result.finish(result.fail("input_missing", f"no fixed structures in {indir}"))
    results, failed = [], []
    for p in files:
        try:
            results.append(extract_one(p, outdir, a))
        except Exception as e:
            (outdir / p.stem.upper()).mkdir(parents=True, exist_ok=True)
            (outdir / p.stem.upper() / "error.txt").write_text(traceback.format_exc())
            err = f"{type(e).__name__}: {str(e)[:400]}"
            failed.append({"id": p.stem.upper(), "code": result.code_of(e), "error": err})
    out = {"ok": not failed, "units": results, "failed": failed}
    return result.finish(out, command="extract", out_dir=a.batch_dir / layout.CONFIG, full=a.full)


if __name__ == "__main__":
    sys.exit(main())
