"""[H] Trajectory analysis of every finished production run.

    gmxprep-analyze --batch-dir B [--only ID ...] [--deffnm md] [--gmx gmx] [--stride 1] [--begin PS] [--force]
                    [--hbond-distance 3.0] [--hbond-angle 150] [--hbond-min-occupancy 0.1]
                    [--ligand-rmsd-max 3.0] [--ca-rmsd-max 4.0]

Per system (06-systems/<id>/analysis/):
  1. PBC removal with GROMACS, by index-group *name* (never by number):
       convert-tpr   solute-only run input (the Solute group is a prefix of the
                     system, so every index group keeps its atom numbers)
       trjconv       -pbc mol -center   (centre on Protein, or the ligand)
       trjconv       -fit rot+trans     (on Backbone, or the ligand alone)
       editconf      solute.gro, the topology MDAnalysis reads (no .tpr version
                     dependence)
  2. Metrics with MDAnalysis, relative to the first production frame:
       Cα RMSD (after fitting Cα) and target heavy-atom RMSD in that protein
       frame (a pose that drifts out of the pocket shows up here); per-residue
       Cα RMSF; protein radius of gyration; target <-> enzyme hydrogen bonds
       (geometric: D-A < 3.0 Å, D-H-A > 150°) with the occupancy of every
       donor/acceptor pair and the count per frame.
  3. analysis.json, rmsd.png, rmsf.png, rg.png and hbonds.png, and one row per
     system in B/08-analysis/analysis_summary.csv for comparing the matrix.

What the target is comes from step F's system.json: `target_resname` for a
ligand, `target_residues` for a protein target folded into the enzyme molecule.

Flags (reported, never fatal): ligand_drift (last-half mean target RMSD above
--ligand-rmsd-max), protein_unstable (last-half mean Cα RMSD above
--ca-rmsd-max), no_target_hbonds (complexes only).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

from . import checks, layout
from .gmx import gmx_cmd, run
from .pdbio import read_gro
from .plotstyle import C1, C2, INK, MUTED
from .plotstyle import figure as _figure
from .plotstyle import save as _save
from .plotstyle import style as _style
from .report import analysis_readme

BACKBONE = {"N", "CA", "C"}


# ------------------------------------------------------------------ index groups
def solute_groups(wd: Path) -> dict[str, list[int]]:
    """Step F's groups restricted to the solute, plus Backbone and C-alpha of the
    protein (computed from atom names in solv_ions.gro, same order as the .tpr)."""
    groups = checks.read_ndx(wd / "index.ndx")
    solute = groups.get("Solute") or []
    if solute != list(range(1, len(solute) + 1)):
        raise ValueError("Solute group is not a prefix of the system; the index was not written by step F")
    out = {k: v for k, v in groups.items() if k not in ("System", "Solvent")}
    prot = set(groups.get("Protein", []))
    if prot:
        atoms, _ = read_gro(wd / "solv_ions.gro")
        out["Backbone"] = [i for i in sorted(prot) if atoms[i - 1].name in BACKBONE]
        out["C-alpha"] = [i for i in sorted(prot) if atoms[i - 1].name == "CA"]
    return out


def pick_groups(groups: dict[str, list[int]], target_resname: str | None) -> tuple[str, str]:
    """(centre group, fit group)."""
    lig = f"LIG_{target_resname}" if target_resname and f"LIG_{target_resname}" in groups else None
    lig = lig or next((g for g in groups if g.startswith("LIG_")), None)
    if groups.get("Protein"):
        return "Protein", "Backbone"
    if lig:
        return lig, lig
    return "Solute", "Solute"


# ------------------------------------------------------------------ PBC treatment
def prepare_trajectory(
    wd: Path, a: argparse.Namespace, groups: dict[str, list[int]], target_resname: str | None
) -> tuple[Path, Path]:
    """-> (solute.gro, fit.xtc) inside wd/analysis."""
    out = wd / "analysis"
    out.mkdir(exist_ok=True)
    log = out / "trjconv.log"
    checks.write_ndx(groups, out / "solute.ndx")
    centre, fit = pick_groups(groups, target_resname)
    tpr, xtc = f"{a.deffnm}.tpr", f"{a.deffnm}.xtc"

    def g(*x):
        return gmx_cmd(a.gmx, *x)

    # Full-system index for the first pass (it reads the full .tpr and .xtc).
    run(g("convert-tpr", "-s", tpr, "-n", "index.ndx", "-o", "analysis/solute.tpr"), wd, log, stdin="Solute\n")
    stride = ["-skip", str(a.stride)] if a.stride > 1 else []
    begin = ["-b", str(a.begin)] if a.begin else []
    run(
        g("trjconv", "-s", tpr, "-f", xtc, "-n", "index.ndx", "-o", "analysis/whole.xtc", "-pbc", "mol", "-center",
          *stride, *begin),
        wd,
        log,
        stdin=f"{centre}\nSolute\n",
    )  # fmt: skip
    run(
        g("trjconv", "-s", "analysis/solute.tpr", "-f", "analysis/whole.xtc", "-n", "analysis/solute.ndx",
          "-o", "analysis/fit.xtc", "-fit", "rot+trans"),
        wd,
        log,
        stdin=f"{fit}\nSolute\n",
    )  # fmt: skip
    run(g("editconf", "-f", "analysis/solute.tpr", "-o", "analysis/solute.gro"), wd, log)
    (out / "whole.xtc").unlink(missing_ok=True)
    return out / "solute.gro", out / "fit.xtc"


# ------------------------------------------------------------------ metrics
def _stats(x) -> dict:
    x = np.asarray(x, float)
    if not len(x):
        return {}
    half = x[len(x) // 2 :]
    return {
        "mean": round(float(x.mean()), 3),
        "std": round(float(x.std()), 3),
        "max": round(float(x.max()), 3),
        "final": round(float(x[-1]), 3),
        "last_half_mean": round(float(half.mean()), 3),
    }


def index_sel(idx: list[int]) -> str | None:
    """1-based GROMACS atom numbers -> an MDAnalysis `index` selection (0-based, inclusive)."""
    if not idx:
        return None
    idx = sorted(idx)
    parts, start = [], idx[0]
    for prev, cur in zip(idx, idx[1:] + [None]):
        if cur != prev + 1:
            parts.append(f"index {start - 1}:{prev - 1}")
            start = cur
    return " or ".join(parts)


def selections(stamp: dict, groups: dict[str, list[int]]) -> tuple[str | None, str | None]:
    """(target, enzyme) MDAnalysis selections. Atoms come from step F's named index
    groups (independent of residue naming); what the target is comes from system.json."""
    role = stamp.get("role")
    tres = stamp.get("target_resname")
    ranges = stamp.get("target_residues") or []
    protein = index_sel(groups.get("Protein", []))
    if role == "enzyme":
        return None, protein
    if tres:  # ligand target (alone or in a complex)
        lig = index_sel(groups.get(f"LIG_{tres}", [])) or f"resname {tres}"
        return lig, (protein if role == "complex" else None)
    if ranges:  # protein target folded into the enzyme moleculetype
        tsel = " or ".join(f"resindex {s}:{e - 1}" for s, e in ranges)
        return f"({protein}) and ({tsel})", f"({protein}) and not ({tsel})"
    return None, protein  # protein target alone: protein metrics only


def hbond_table(u, target: str, partner: str, d_a: float, angle: float) -> tuple[np.ndarray, list[dict]]:
    """Geometric target<->partner H-bonds. Donors are paired with hydrogens by distance,
    so no bond information (and no .tpr) is needed."""
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

    polar = "(name N* or name O* or name S*)"
    both = f"(({target}) or ({partner}))"
    hba = HydrogenBondAnalysis(
        u,
        donors_sel=f"{both} and {polar}",
        hydrogens_sel=f"{both} and name H*",
        acceptors_sel=f"{both} and {polar}",
        between=[target, partner],
        d_a_cutoff=d_a,
        d_h_a_angle_cutoff=angle,
    )
    hba.run()
    hb = hba.results.hbonds  # frame, donor, hydrogen, acceptor, distance, angle
    n = u.trajectory.n_frames
    per_frame = np.bincount(hb[:, 0].astype(int), minlength=n) if len(hb) else np.zeros(n, int)
    pairs: dict[tuple[int, int], set[int]] = {}
    for f, d, _h, acc, *_ in hb:
        pairs.setdefault((int(d), int(acc)), set()).add(int(f))

    def label(i):
        at = u.atoms[i]
        return f"{at.resname}{at.resid}:{at.name}"

    table = sorted(
        (
            {
                "donor": label(d),
                "acceptor": label(acc),
                "occupancy": round(len(fr) / n, 3),
                "target_is_donor": bool(u.atoms[[d]].select_atoms(target).n_atoms),
            }
            for (d, acc), fr in pairs.items()
        ),
        key=lambda r: -r["occupancy"],
    )
    return per_frame, table


def compute_metrics(
    top: Path,
    traj: Path,
    target: str | None,
    enzyme: str | None,
    a: argparse.Namespace,
    labels: list[str] | None = None,
) -> dict:
    """`labels`: author residue labels in protein residue order (step F's system.json)."""
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms

    u = mda.Universe(str(top), str(traj))
    n = u.trajectory.n_frames
    if n < 2:
        raise ValueError(f"{traj.name} has {n} frame(s); nothing to analyse")
    times = np.array([ts.time for ts in u.trajectory]) / 1000.0  # ns
    res: dict = {"frames": n, "time_ns": [round(float(t), 4) for t in times]}

    ca_sel = f"({enzyme}) and name CA" if enzyme else None
    tgt_heavy = f"({target}) and not name H*" if target else None
    if ca_sel and not u.select_atoms(ca_sel).n_atoms:
        ca_sel = None
    if tgt_heavy and not u.select_atoms(tgt_heavy).n_atoms:
        raise ValueError(f"target selection '{target}' matches no atoms in {top.name}")

    if ca_sel:
        groups = [tgt_heavy] if tgt_heavy else []
        r = rms.RMSD(u, select=ca_sel, groupselections=groups, ref_frame=0).run()
        res["ca_rmsd"] = [round(float(x), 3) for x in r.results.rmsd[:, 2]]
        if tgt_heavy:  # in the protein's frame: no refit on the ligand
            res["target_rmsd"] = [round(float(x), 3) for x in r.results.rmsd[:, 3]]
        ca = u.select_atoms(ca_sel)
        f = rms.RMSF(ca).run()
        res["rmsf"] = {
            "resindex": [int(x) for x in ca.resindices],
            "resid": [int(x) for x in ca.resids],
            "resname": [str(x) for x in ca.resnames],
            "rmsf": [round(float(x), 3) for x in f.results.rmsf],
        }
        if labels:  # pdb2gmx -renum: residue index i is the i-th author residue
            res["rmsf"]["label"] = [labels[i] if i < len(labels) else None for i in res["rmsf"]["resindex"]]
        prot = u.select_atoms(enzyme)
        res["rg"] = [round(float(prot.radius_of_gyration()), 3) for _ in u.trajectory]
    elif tgt_heavy:  # target alone: internal RMSD after fitting the target itself
        r = rms.RMSD(u, select=tgt_heavy, ref_frame=0).run()
        res["target_rmsd"] = [round(float(x), 3) for x in r.results.rmsd[:, 2]]

    if target and enzyme:
        per_frame, table = hbond_table(u, target, enzyme, a.hbond_distance, a.hbond_angle)
        res["hbonds"] = {
            "count": [int(x) for x in per_frame],
            "pairs": [p for p in table if p["occupancy"] >= a.hbond_min_occupancy],
            "pairs_below_threshold": sum(p["occupancy"] < a.hbond_min_occupancy for p in table),
        }

    res["summary"] = {k: _stats(res[k]) for k in ("ca_rmsd", "target_rmsd", "rg") if k in res}
    if "hbonds" in res:
        res["summary"]["hbond_count"] = _stats(res["hbonds"]["count"])
    return res


def flags(res: dict, role: str, a: argparse.Namespace) -> list[str]:
    s, out = res["summary"], []
    if s.get("target_rmsd", {}).get("last_half_mean", 0) > a.ligand_rmsd_max:
        out.append(
            f"ligand_drift: target RMSD {s['target_rmsd']['last_half_mean']} Å over the last half "
            f"(> {a.ligand_rmsd_max} Å) - pose not retained"
        )
    if s.get("ca_rmsd", {}).get("last_half_mean", 0) > a.ca_rmsd_max:
        out.append(f"protein_unstable: Cα RMSD {s['ca_rmsd']['last_half_mean']} Å over the last half")
    if role == "complex" and "hbonds" in res and not res["hbonds"]["pairs"]:
        out.append(f"no_target_hbonds: no target-enzyme H-bond above {a.hbond_min_occupancy:.0%} occupancy")
    return out


# ------------------------------------------------------------------ plots
def _hbond_plot(hb: dict, t: list[float], sysid: str, path: Path) -> str:
    fig, (ax1, ax2) = _figure(2, (9.6, 3.8), gridspec_kw={"width_ratios": [1.2, 1]})
    _style(ax1, f"{sysid}: target-enzyme H-bonds per frame", "time (ns)", "H-bonds", grid="y")
    _style(ax2, "persistence (donor → acceptor)", "occupancy (% of frames)", None, grid="x")
    ax1.plot(t, hb["count"], color=C1, lw=1.0)
    ax1.set_ylim(bottom=0)
    top = hb["pairs"][:10]
    if top:
        labels = [f"{p['donor']} → {p['acceptor']}" for p in top][::-1]
        occ = [100 * p["occupancy"] for p in top][::-1]
        ax2.barh(labels, occ, color=C1, height=0.6)
        for y, v in enumerate(occ):
            ax2.text(v + 1, y, f"{v:.0f}%", va="center", fontsize=7, color=MUTED)
        ax2.set_xlim(0, 110)
        ax2.set_ylim(-0.75, max(len(top), 5) - 0.25)  # constant bar thickness for 1-10 pairs
    else:
        ax2.text(0.5, 0.5, "no persistent H-bonds", ha="center", va="center", color=MUTED, fontsize=9)
        ax2.set_yticks([])
    return _save(fig, path)


def plots(res: dict, sysid: str, out: Path) -> list[str]:
    t = res["time_ns"]
    made = []
    if "ca_rmsd" in res or "target_rmsd" in res:
        fig, ax = _figure()
        _style(ax, f"{sysid}: RMSD from the first production frame", "time (ns)", "RMSD (Å)")
        if "ca_rmsd" in res:
            ax.plot(t, res["ca_rmsd"], color=C1, lw=1.2, label="protein Cα")
        if "target_rmsd" in res:
            ax.plot(t, res["target_rmsd"], color=C2, lw=1.2, label="target heavy atoms")
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK)
        made.append(_save(fig, out / "rmsd.png"))
    if "rmsf" in res:
        r = res["rmsf"]
        fig, ax = _figure(size=(7.2, 3.4))
        _style(ax, f"{sysid}: per-residue Cα RMSF", "residue (sequential index)", "RMSF (Å)")
        ax.plot(np.array(r["resindex"]) + 1, r["rmsf"], color=C1, lw=1.0)
        ax.set_ylim(bottom=0)
        made.append(_save(fig, out / "rmsf.png"))
    if "rg" in res:
        fig, ax = _figure()
        _style(ax, f"{sysid}: protein radius of gyration", "time (ns)", "Rg (Å)")
        ax.plot(t, res["rg"], color=C1, lw=1.2)
        made.append(_save(fig, out / "rg.png"))
    if "hbonds" in res:
        made.append(_hbond_plot(res["hbonds"], t, sysid, out / "hbonds.png"))
    return made


# ------------------------------------------------------------------ driver
SUMMARY_COLS = [
    "id", "role", "enzyme", "target", "status", "frames", "time_ns",
    "ca_rmsd_last_half", "target_rmsd_last_half", "rg_mean", "hbond_mean", "top_hbond", "top_hbond_occupancy",
    "flags",
]  # fmt: skip


def settings(a: argparse.Namespace) -> dict:
    """Everything a stored analysis depends on; a cached result is reused only if these match."""
    return {
        "deffnm": a.deffnm,
        "stride": a.stride,
        "begin_ps": a.begin,
        "hbond": {"d_a_A": a.hbond_distance, "d_h_a_deg": a.hbond_angle, "min_occupancy": a.hbond_min_occupancy},
        "flags": {"ligand_rmsd_max": a.ligand_rmsd_max, "ca_rmsd_max": a.ca_rmsd_max},
    }


def analyse_system(wd: Path, a: argparse.Namespace) -> dict:
    stamp = json.loads((wd / "system.json").read_text())
    out = wd / "analysis"
    report_p = out / "analysis.json"
    xtc = wd / f"{a.deffnm}.xtc"
    run_state = load_json(wd / "run.json")
    if not xtc.exists() or run_state.get("status") not in (None, "done"):
        return {"id": wd.name, "status": "skipped", "reason": "production run not finished (gmxprep-run)"}
    old = load_json(report_p)
    if (
        not a.force
        and old.get("status") == "ok"
        and old.get("settings") == settings(a)
        and report_p.stat().st_mtime > xtc.stat().st_mtime
    ):
        return old
    groups = solute_groups(wd)
    target, enzyme = selections(stamp, groups)
    top, traj = prepare_trajectory(wd, a, groups, stamp.get("target_resname"))
    res = compute_metrics(top, traj, target, enzyme, a, stamp.get("protein_residue_labels"))
    report = {
        "id": wd.name,
        "status": "ok",
        "role": stamp.get("role"),
        "enzyme": stamp.get("enzyme"),
        "target": stamp.get("target"),
        "selections": {"target": target, "enzyme": enzyme},
        "settings": settings(a),
        **res,
    }
    report["flags"] = flags(res, stamp.get("role"), a)
    report["plots"] = plots(res, wd.name, out)
    report_p.write_text(json.dumps(report, indent=1))
    (out / "README.md").write_text(analysis_readme(report))
    return report


def load_json(path: Path) -> dict:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return {}


def write_summary(b: Path) -> int:
    """analysis_summary.csv from every system's analysis.json, so `--only X` never drops
    the other systems' rows. Returns the number of rows."""
    rows = []
    for p in sorted((b / layout.SYSTEMS).glob("*/analysis/analysis.json")):
        r = load_json(p)
        if r.get("status") == "ok":
            rows.append(summary_row(r))
    if rows:
        (b / layout.ANALYSIS_SUMMARY).parent.mkdir(parents=True, exist_ok=True)
        with open(b / layout.ANALYSIS_SUMMARY, "w", newline="") as fh:
            w = csv.DictWriter(fh, SUMMARY_COLS)
            w.writeheader()
            w.writerows(rows)
    return len(rows)


def summary_row(r: dict) -> dict:
    s = r.get("summary", {})
    top = (r.get("hbonds", {}).get("pairs") or [None])[0]
    return {
        "id": r["id"],
        "role": r.get("role"),
        "enzyme": r.get("enzyme"),
        "target": r.get("target"),
        "status": r["status"],
        "frames": r.get("frames"),
        "time_ns": r["time_ns"][-1] if r.get("time_ns") else None,
        "ca_rmsd_last_half": s.get("ca_rmsd", {}).get("last_half_mean"),
        "target_rmsd_last_half": s.get("target_rmsd", {}).get("last_half_mean"),
        "rg_mean": s.get("rg", {}).get("mean"),
        "hbond_mean": s.get("hbond_count", {}).get("mean"),
        "top_hbond": f"{top['donor']}->{top['acceptor']}" if top else None,
        "top_hbond_occupancy": top["occupancy"] if top else None,
        "flags": "; ".join(f.split(":")[0] for f in r.get("flags", [])),
    }


def cmd_analyze(a: argparse.Namespace) -> dict:
    b: Path = a.batch_dir
    results = []
    for stamp_p in sorted((b / layout.SYSTEMS).glob("*/system.json")):
        wd = stamp_p.parent
        if a.only and wd.name not in a.only:
            continue
        try:
            r = analyse_system(wd, a)
        except Exception as e:  # one bad system (tool, MDAnalysis, numpy) never stops the batch
            r = {"id": wd.name, "status": "failed", "reason": f"{type(e).__name__}: {str(e)[:400]}"}
            (wd / "analysis").mkdir(exist_ok=True)
            (wd / "analysis" / "analysis.json").write_text(json.dumps(r, indent=1))
        results.append(r)
    if not results:
        return {"ok": False, "error": f"no built systems under {b}/systems"}
    done = [r for r in results if r["status"] == "ok"]
    n_rows = write_summary(b)
    failed = [{"id": r["id"], "reason": r["reason"]} for r in results if r["status"] == "failed"]
    return {
        "ok": not failed and bool(done),
        "analysed": len(done),
        "skipped": [{"id": r["id"], "reason": r["reason"]} for r in results if r["status"] == "skipped"],
        "failed": failed,
        "flags": {r["id"]: r["flags"] for r in done if r.get("flags")},
        "summary_csv": str(b / layout.ANALYSIS_SUMMARY) if n_rows else None,
        "summary_rows": n_rows,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gmxprep-analyze", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--only", nargs="+")
    ap.add_argument("--deffnm", default="md", help="production run name used by gmxprep-run")
    ap.add_argument("--gmx", default="gmx")
    ap.add_argument("--stride", type=int, default=1, help="use every Nth frame")
    ap.add_argument("--begin", type=float, help="ps; discard frames before this time")
    ap.add_argument("--force", action="store_true", help="re-analyse even if analysis.json is current")
    ap.add_argument("--hbond-distance", type=float, default=3.0, help="Å, donor-acceptor")
    ap.add_argument("--hbond-angle", type=float, default=150.0, help="degrees, donor-H-acceptor")
    ap.add_argument("--hbond-min-occupancy", type=float, default=0.1, help="fraction of frames to report a pair")
    ap.add_argument("--ligand-rmsd-max", type=float, default=3.0, help="Å, ligand_drift flag")
    ap.add_argument("--ca-rmsd-max", type=float, default=4.0, help="Å, protein_unstable flag")
    a = ap.parse_args(argv)
    try:
        out = cmd_analyze(a)
    except Exception as e:  # every failure is reported as the one-line JSON contract, never a traceback
        out = {"ok": False, "error": f"{type(e).__name__}: {e}"}
    print(json.dumps(out), flush=True)
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
