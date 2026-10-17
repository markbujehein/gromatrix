"""[E] Compare relaxed experimental sub-units with AlphaFold DB predictions.

    gromatrix-afcompare --batch-dir B [--in relaxed] [--ids 1ABC ...]
                      [--files enzyme.pdb target_*.pdb] [--uniprot A:P12345,B:Q67890[@ID]]
                      [--plddt-min 70] [--min-length 30] [--tm-min 0.8] [--lddt-min 0.7]
                      [--out 05-alphafold/af_report.csv] [--afdb-dir 05-alphafold-cache]

Per protein chain:
  1. UniProt accession: --uniprot if given, otherwise from SIFTS (PDBe API; RCSB
     Data API when PDBe has no record) using the chain's *author* ID. For fusion
     constructs the accession covering most of the chain is used. Chains renamed
     during assembly building are traced back through 03-units/<ID>/unit.json.
  2. AFDB model: the prediction API gives the current model URL (model version
     recorded). Files are cached by sha256 in <batch>/<afdb-dir>/.
  3. Residue correspondence by global sequence alignment of the experimental
     chain to the model (BLOSUM62, free end gaps), not by residue numbering.
  4. Metrics over aligned Calpha atoms:
       rmsd_all   Kabsch RMSD over every aligned residue
       rmsd_core  Kabsch RMSD over residues with AF pLDDT >= --plddt-min
       tm_score   TM-score normalised by the experimental chain length, using the
                  TM-score program's procedure (superpositions seeded on fragments
                  of length L, L/2, L/4 ..., iterated on residues within d0)
       lddt_ca    superposition-free Calpha lDDT (15 Å radius; 0.5/1/2/4 Å thresholds)
                  with the experimental structure as reference
     Flags: tm_score < --tm-min or lddt_ca < --lddt-min. Chains shorter than
     --min-length (e.g. coactivator peptides) are reported but not scored.
     These metrics are not meaningful for chains that short.
  5. The AF model, superposed onto the experimental chain, is written to
     <batch>/05-alphafold/<ID>/AF_<chain>_<UniProt>.pdb for visual inspection.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

from . import layout
from .net import Cache, HTTPError, atomic_write, get_json
from .pdbio import PROTEIN_RES, read_pdb, write_pdb

AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pid}"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pid}"
RCSB_ENTITY = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{eid}"
THREE2ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "HSD": "H",
    "HSE": "H",
    "HSP": "H",
    "CYX": "C",
    "CYM": "C",
    "ASH": "D",
    "GLH": "E",
    "LYN": "K",
    "MSE": "M",
}


# ------------------------------------------------------------------ geometry
def kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t minimising |(p @ R.T + t) - q|."""
    pc, qc = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, qc - pc @ r.T


def rmsd_after_fit(p: np.ndarray, q: np.ndarray) -> float:
    r, t = kabsch(p, q)
    return float(np.sqrt((((p @ r.T + t) - q) ** 2).sum(1).mean()))


def tm_score(model: np.ndarray, ref: np.ndarray, l_norm: int) -> tuple[float, np.ndarray, np.ndarray]:
    """TM-score of aligned Calpha pairs, normalised by l_norm. Returns (tm, R, t)
    for the best superposition (model -> ref)."""
    n = len(model)
    d0 = max(0.5, 1.24 * (max(l_norm, 19) - 15) ** (1 / 3) - 1.8)
    d_search = min(max(d0, 4.5), 8.0)
    best = (-1.0, np.eye(3), np.zeros(3))
    lengths, L = [], n
    while L >= 4:
        lengths.append(L)
        L //= 2
    for seg in lengths:
        step = max(1, seg // 2 if seg < n else n)
        for start in range(0, n - seg + 1, step):
            idx = np.arange(start, start + seg)
            for _ in range(20):
                r, t = kabsch(model[idx], ref[idx])
                d = np.linalg.norm(model @ r.T + t - ref, axis=1)
                tm = float((1.0 / (1.0 + (d / d0) ** 2)).sum() / l_norm)
                if tm > best[0]:
                    best = (tm, r, t)
                cut = d_search
                new = np.where(d < cut)[0]
                while len(new) < 3 and cut < 50:
                    cut += 0.5
                    new = np.where(d < cut)[0]
                if np.array_equal(new, idx):
                    break
                idx = new
    return best


def lddt_ca(model: np.ndarray, ref: np.ndarray, radius: float = 15.0) -> float:
    """Superposition-free Calpha lDDT over the *reference* residues. Rows of `model`
    that are NaN (reference residues the model lacks) count as not preserved, as in
    the standard definition, so a fragment covering half the chain cannot score 1.0."""
    dr = np.linalg.norm(ref[:, None] - ref[None], axis=2)
    dm = np.linalg.norm(model[:, None] - model[None], axis=2)
    mask = (dr < radius) & ~np.eye(len(ref), dtype=bool)
    if not mask.any():
        return float("nan")
    diff = np.abs(dr - dm)[mask]  # NaN wherever a residue is missing from the model
    with np.errstate(invalid="ignore"):
        return float(np.mean([(diff < t).mean() for t in (0.5, 1.0, 2.0, 4.0)]))


# ------------------------------------------------------------------ data
def chain_ca(atoms) -> dict[str, list]:
    out: dict[str, list] = {}
    for a in atoms:
        if a.name == "CA" and a.resname in PROTEIN_RES:
            out.setdefault(a.chain, []).append(a)
    return out


def align(seq_exp: str, seq_af: str) -> list[tuple[int, int]]:
    from Bio import Align
    from Bio.Align import substitution_matrices

    al = Align.PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10.0, -0.5
    al.end_gap_score = 0.0
    best = al.align(seq_exp, seq_af)[0]
    pairs = []
    for (a0, a1), (b0, b1) in zip(*best.aligned):
        pairs += list(zip(range(a0, a1), range(b0, b1)))
    return pairs


def sifts_map(pid: str, cache_dir: Path) -> dict[str, str]:
    """author chain -> UniProt accession. Fusion constructs map one chain to several
    accessions, so the accession covering the most residues of that chain wins.
    If PDBe has no SIFTS record, the RCSB Data API (also SIFTS-derived) is used instead."""
    f = cache_dir / f"sifts_{pid.lower()}.json"
    if f.exists():
        return json.loads(f.read_text())
    try:
        data = get_json(SIFTS.format(pid=pid.lower()))
        cover: dict[str, dict[str, int]] = {}
        for acc, v in data.get(pid.lower(), {}).get("UniProt", {}).items():
            for m in v["mappings"]:
                c = cover.setdefault(m["chain_id"], {})
                c[acc] = c.get(acc, 0) + (m["unp_end"] - m["unp_start"] + 1)
        out = {ch: max(c, key=c.get) for ch, c in cover.items()}
    except HTTPError as e:
        if e.status != 404:
            raise
        out = rcsb_map(pid)
    atomic_write(f, json.dumps(out))
    return out


def rcsb_map(pid: str) -> dict[str, str]:
    out = {}
    entry = get_json(RCSB_ENTRY.format(pid=pid.upper()))
    for eid in entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []):
        ids = get_json(RCSB_ENTITY.format(pid=pid.upper(), eid=eid)).get(
            "rcsb_polymer_entity_container_identifiers", {}
        )
        refs = [r for r in ids.get("reference_sequence_identifiers") or [] if r.get("database_name") == "UniProt"]
        if refs:
            best = max(refs, key=lambda r: r.get("entity_sequence_coverage") or 0)["database_accession"]
            for ch in ids.get("auth_asym_ids", []):
                out[ch] = best
    return out


def af_model(acc: str, cache: Cache, seq_exp: str):
    """Fetch the AFDB entry (fragment) that aligns best to the experimental sequence."""
    entries = get_json(AF_API.format(acc=acc)) or []
    entries = [e for e in entries if e.get("uniprotAccession") == acc] or entries
    best = None
    for e in entries:
        name = e["pdbUrl"].rsplit("/", 1)[-1]
        path, _ = cache.fetch(name, e["pdbUrl"], extra={"uniprot": acc, "version": e.get("latestVersion")})
        ca = chain_ca(read_pdb(path))
        model_ca = next(iter(ca.values()))
        seq = "".join(THREE2ONE.get(a.resname, "X") for a in model_ca)
        # Rank by identical aligned residues: with free end gaps, the raw pair count
        # rewards long unrelated sequences.
        n = sum(seq_exp[i] == seq[j] for i, j in align(seq_exp, seq))
        if best is None or n > best[0]:
            best = (n, e, path, model_ca, seq)
    if best is None:
        raise LookupError(f"no AlphaFold DB entry for {acc}")
    return best[1:]


def compare_chain(chain: str, exp_ca: list, acc: str, cache: Cache, a, out_dir: Path) -> dict:
    seq_exp = "".join(THREE2ONE.get(x.resname, "X") for x in exp_ca)
    entry, model_path, model_ca, seq_af = af_model(acc, cache, seq_exp)
    pairs = align(seq_exp, seq_af)
    ident = sum(seq_exp[i] == seq_af[j] for i, j in pairs) / max(1, len(pairs))
    full_ref = np.array([[x.x, x.y, x.z] for x in exp_ca])
    full_mod = np.full_like(full_ref, np.nan)
    for i, j in pairs:
        full_mod[i] = (model_ca[j].x, model_ca[j].y, model_ca[j].z)
    ref = np.array([[exp_ca[i].x, exp_ca[i].y, exp_ca[i].z] for i, _ in pairs])
    mod = np.array([[model_ca[j].x, model_ca[j].y, model_ca[j].z] for _, j in pairs])
    plddt = np.array([model_ca[j].b for _, j in pairs])
    core = plddt >= a.plddt_min
    tm, r, t = tm_score(mod, ref, len(exp_ca))
    res = {
        "af_entry": entry.get("entryId"),
        "af_version": entry.get("latestVersion"),
        "L_exp": len(exp_ca),
        "aligned": len(pairs),
        "coverage": round(len(pairs) / len(exp_ca), 3),
        "identity": round(ident, 3),
        "mean_plddt_aligned": round(float(plddt.mean()), 1),
        "n_core": int(core.sum()),
        "rmsd_all": round(rmsd_after_fit(mod, ref), 2),
        "rmsd_core": round(rmsd_after_fit(mod[core], ref[core]), 2) if core.sum() >= 3 else None,
        "tm_score": round(tm, 3),
        "lddt_ca": round(lddt_ca(full_mod, full_ref), 3),
    }
    flags = []
    if ident < 0.9:
        flags.append(f"sequence identity {ident:.2f}: check the UniProt mapping or engineered mutations")
    if tm < a.tm_min:
        flags.append(f"TM-score {tm:.2f} < {a.tm_min}")
    if res["lddt_ca"] < a.lddt_min:
        flags.append(f"lDDT {res['lddt_ca']:.2f} < {a.lddt_min}")
    res["flags"] = flags
    # superposed AF model for inspection
    full = read_pdb(model_path)
    xyz = np.array([[x.x, x.y, x.z] for x in full]) @ r.T + t
    out_dir.mkdir(parents=True, exist_ok=True)
    write_pdb(
        out_dir / f"AF_{chain}_{acc}.pdb",
        [replace(x, chain=chain, x=float(p[0]), y=float(p[1]), z=float(p[2])) for x, p in zip(full, xyz)],
    )
    return res


def _scoped_map(spec: str | None, pid: str) -> dict[str, str]:
    out = {}
    for item in (spec or "").split(","):
        item = item.strip()
        if not item:
            continue
        body, _, scope = item.partition("@")
        if scope and scope.upper() != pid:
            continue
        ch, _, acc = body.partition(":")
        out[ch] = acc.upper()
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gromatrix-afcompare", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--in", dest="indir", default=layout.RELAXED)
    ap.add_argument("--ids", nargs="+")
    ap.add_argument("--files", nargs="+", default=["enzyme.pdb", "target_*.pdb"])
    ap.add_argument("--uniprot", help="explicit chain:accession map, e.g. A:P12345,B:Q67890@1ABC")
    ap.add_argument("--plddt-min", type=float, default=70.0)
    ap.add_argument("--min-length", type=int, default=30)
    ap.add_argument("--tm-min", type=float, default=0.8)
    ap.add_argument("--lddt-min", type=float, default=0.7)
    ap.add_argument("--out", default=layout.AF_REPORT)
    ap.add_argument("--afdb-dir", default=layout.ALPHAFOLD_CACHE)
    a = ap.parse_args(argv)

    root = a.batch_dir / a.indir
    units = sorted(p for p in root.iterdir() if p.is_dir()) if root.exists() else []
    if a.ids:
        units = [u for u in units if u.name.upper() in {i.upper() for i in a.ids}]
    if not units:
        print(json.dumps({"ok": False, "error": f"no structures in {root}"}))
        return 1
    cache = Cache(a.batch_dir / a.afdb_dir)
    rows, flagged, failed = [], [], []
    for u in units:
        pid = u.name.upper()
        origin = {}
        uj = a.batch_dir / layout.UNITS / u.name / "unit.json"
        if uj.exists():
            origin = {c: o.split("/")[0] for c, o in json.loads(uj.read_text()).get("chain_origin", {}).items()}
        manual = _scoped_map(a.uniprot, pid)
        chains: dict[str, list] = {}
        # Only score structures that passed step D: skip *.failed.pdb and anything relax.json marks not ok.
        relax_status = {}
        if (u / "relax.json").exists():
            relax_status = {k: v.get("status") for k, v in json.loads((u / "relax.json").read_text())["files"].items()}
        for f in sorted(u.glob("*.pdb")):
            if f.name.endswith(".failed.pdb") or relax_status.get(f.name, "ok") not in ("ok", "skipped"):
                continue
            if any(fnmatch.fnmatch(f.name, pat) for pat in a.files):
                for ch, ca in chain_ca(read_pdb(f)).items():
                    chains.setdefault(ch, ca)
        try:
            auto = {} if manual and set(manual) >= set(chains) else sifts_map(pid, cache.dir)
        except HTTPError as e:
            auto = {}
            failed.append({"id": pid, "error": f"SIFTS lookup failed ({e}); pass --uniprot"})
        details = {}
        for ch, ca in chains.items():
            acc = manual.get(ch) or auto.get(origin.get(ch, ch))
            row = {"id": pid, "chain": ch, "origin_chain": origin.get(ch, ch), "uniprot": acc}
            if not acc:
                row["status"] = "no UniProt mapping"
            elif len(ca) < a.min_length:
                row["status"] = f"skipped: {len(ca)} residues < --min-length"
            else:
                try:
                    row.update(compare_chain(ch, ca, acc, cache, a, a.batch_dir / layout.ALPHAFOLD / u.name))
                    row["status"] = "flagged" if row["flags"] else "ok"
                except (HTTPError, LookupError) as e:
                    row["status"] = f"error: {e}"
                    failed.append({"id": pid, "chain": ch, "error": str(e)[:200]})
            details[ch] = row
            rows.append(row)
            if row["status"] == "flagged":
                flagged.append({k: row[k] for k in ("id", "chain", "uniprot", "tm_score", "lddt_ca", "flags")})
        (a.batch_dir / layout.ALPHAFOLD).mkdir(parents=True, exist_ok=True)
        (a.batch_dir / layout.ALPHAFOLD / f"{pid}.json").write_text(json.dumps(details, indent=1))

    cols = [
        "id",
        "chain",
        "origin_chain",
        "uniprot",
        "status",
        "af_entry",
        "af_version",
        "L_exp",
        "aligned",
        "coverage",
        "identity",
        "mean_plddt_aligned",
        "n_core",
        "rmsd_all",
        "rmsd_core",
        "tm_score",
        "lddt_ca",
    ]
    (a.batch_dir / a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.batch_dir / a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    scored = [r for r in rows if r.get("tm_score") is not None]
    print(
        json.dumps(
            {
                "ok": not failed,
                "report": str(a.batch_dir / a.out),
                "chains_scored": len(scored),
                "summary": [
                    {k: r.get(k) for k in ("id", "chain", "uniprot", "rmsd_core", "tm_score", "lddt_ca", "status")}
                    for r in rows
                ],
                "flagged": flagged,
                "failed": failed,
            }
        )
    )
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
