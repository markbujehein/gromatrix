"""Human-readable reports for a batch: Markdown tables and figures built from the step artefacts.

    gromatrix-report --batch-dir B [--no-figures]

Reads only the files each step already wrote (the same contract as `gromatrix-status`), so it can
run at any point and again after every step. It writes

    B/09-reports/README.md      one section per step B-H: what happened, as tables, with figures
    B/09-reports/figures/*.png  relax gates, AlphaFold agreement, system sizes, MD timeline,
                                RMSD overlay and the receptor x target matrix

Per-system plot folders (B/06-systems/<ID>/analysis/) get their own README.md from
`analysis_readme`, which `gromatrix-analyze` calls after drawing the plots. Figures need matplotlib
(the `analysis` extra); without it the Markdown is still written.

Prints one JSON line {"ok": ..., "report": path, "figures": [...]}.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

from . import layout
from .plotstyle import BAD, C1, C2, GRID, INK, MUTED, OK, have_matplotlib
from .plotstyle import figure as _figure
from .plotstyle import save as _save
from .plotstyle import style as _style

# Thresholds the step defaults use; drawn as reference lines only (the stored JSON is authoritative).
TM_MIN, LDDT_MIN, MAX_RMSD = 0.8, 0.7, 2.5
SHORT_RUN_NS = 10.0


# ------------------------------------------------------------------ small helpers
def _json(path: Path) -> dict | None:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return None


def _num(x, nd: int = 2) -> str:
    return "-" if x is None else f"{x:.{nd}f}" if isinstance(x, float) else str(x)


def _table(headers: list[str], rows: list[list]) -> str:
    if not rows:
        return "_nothing to show yet._\n"
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join("" if c is None else str(c).replace("|", "/") for c in r) + " |" for r in rows]
    return "\n".join(out) + "\n"


def _warn_block(warnings: list[tuple[str, str]]) -> str:
    return "".join(f"- **{w[0]}**: {w[1]}\n" for w in warnings) + ("\n" if warnings else "")


def _figure_md(name: str, caption: str, how: str) -> str:
    return f"![{caption}](figures/{name})\n\n*{caption}* **How to read it:** {how}\n\n"


# ------------------------------------------------------------------ per-system plot README
PLOT_GUIDE = {
    "rmsd.png": (
        "RMSD from the first production frame",
        "Blue is the protein Cα RMSD after fitting the Cα atoms. Orange is the target's heavy-atom "
        "RMSD measured in the protein's frame, *without* refitting on the target, so it counts both "
        "internal change and the target moving or leaving the pocket. A curve that rises and then "
        "flattens means the system has settled; a steady climb means it is still relaxing or the "
        "target is drifting. Protein Cα near 1-3 Å and target under about 3 Å is the usual picture "
        "for a stable bound pose; the flags in `analysis.json` use 4 Å and 3 Å.",
    ),
    "rmsf.png": (
        "Per-residue Cα RMSF",
        "Root-mean-square fluctuation of each Cα around its mean position over the trajectory. Peaks "
        "are flexible loops and chain termini; low values are the rigid helix and sheet core and "
        "usually the binding pocket. The x axis is the sequential residue index; the author labels "
        "(chain:number) are in `analysis.json` under `rmsf.label`. Compare peaks between systems to "
        "see what the target stiffens or loosens.",
    ),
    "rg.png": (
        "Protein radius of gyration",
        "How compact the protein is over time. A flat line is a stable fold; a slow rise means "
        "swelling or opening, a drop means collapse. Changes of a few tenths of an ångström are noise.",
    ),
    "hbonds.png": (
        "Target-enzyme hydrogen bonds",
        "Left: number of target-enzyme H-bonds in each frame. Right: how persistent the individual "
        "donor to acceptor pairs are, as the percentage of frames in which each exists (only pairs "
        "above the occupancy threshold are drawn). A bar near 100% is a stable anchor point; "
        "several bars at 20-50% mean a bond that forms and breaks. Hydrophobic ligands may show none.",
    ),
}


def analysis_readme(rep: dict) -> str:
    """README.md for a system's analysis/ folder, from its analysis.json content."""
    s = rep.get("summary", {})
    st = rep.get("settings", {})
    ns = (rep.get("time_ns") or [0])[-1]
    lines = [
        f"# Analysis of `{rep.get('id')}`",
        "",
        f"Role: **{rep.get('role')}**. Enzyme: {rep.get('enzyme') or '-'}. Target: {rep.get('target') or '-'}. "
        f"{rep.get('frames')} frames covering {ns:g} ns of production (stride {st.get('stride')}).",
        "",
    ]
    if ns < SHORT_RUN_NS:
        lines += [
            f"> **Short trajectory ({ns:g} ns).** This shows equilibration and pose retention at most. "
            "It cannot support a stability or binding claim; one replica is a single sample.",
            "",
        ]
    if rep.get("flags"):
        lines += ["## Flags", ""] + [f"- {f}" for f in rep["flags"]] + [""]
    rows = []
    for key, label, unit in (
        ("ca_rmsd", "Protein Cα RMSD", "Å"),
        ("target_rmsd", "Target RMSD", "Å"),
        ("rg", "Radius of gyration", "Å"),
        ("hbond_count", "Target-enzyme H-bonds", "per frame"),
    ):
        if key in s:
            v = s[key]
            rows.append([label, unit, v["last_half_mean"], v["mean"], v["std"], v["max"]])
    if rows:
        lines += ["## Numbers", "", _table(["metric", "unit", "last-half mean", "mean", "std", "max"], rows)]
    lines += ["## Plots in this folder", ""]
    for name in rep.get("plots", []):
        title, how = PLOT_GUIDE.get(name, (name, ""))
        lines += [f"### `{name}`: {title}", "", how, ""]
    top = (rep.get("hbonds") or {}).get("pairs") or []
    if top:
        lines += [
            "## Persistent H-bonds",
            "",
            _table(
                ["donor", "acceptor", "occupancy"],
                [[p["donor"], p["acceptor"], f"{100 * p['occupancy']:.0f}%"] for p in top[:10]],
            ),
        ]
    hb = st.get("hbond", {})
    lines += [
        "## Other files",
        "",
        "- `analysis.json`: every number behind the plots (time series, RMSF, H-bond table, summary, flags).",
        "- `fit.xtc`, `solute.gro`, `solute.tpr`: the periodic-boundary-corrected, fitted solute trajectory, "
        "for opening in a viewer.",
        "",
        f"H-bond criteria: donor-acceptor distance ≤ {hb.get('d_a_A')} Å, angle ≥ {hb.get('d_h_a_deg')}°, "
        f"reported above {hb.get('min_occupancy')} occupancy. RMSD reference is the first production frame.",
        "",
    ]
    return "\n".join(lines)


# ------------------------------------------------------------------ sections B-H
def sec_fix(b: Path) -> tuple[str, list[str]]:
    rows, warns = [], []
    for p in sorted((b / layout.FIXED).glob("*.meta.json")):
        m = (_json(p) or {}).get("fix", {})
        pid = p.name.split(".")[0]
        gaps = m.get("gaps", [])
        rem = m.get("removed", {})
        rows.append(
            [
                pid,
                m.get("built_atoms"),
                f"{sum(g.get('modelled', False) for g in gaps)}/{sum(not g.get('modelled', False) for g in gaps)}",
                len(rem.get("waters", [])) if isinstance(rem.get("waters"), list) else rem.get("waters"),
                rem.get("hydrogens"),
                (m.get("raw") or {}).get("altloc_residues"),
            ]
        )
        warns += [(pid, w) for w in m.get("warnings", [])]
    t = _table(["ID", "atoms built", "gaps modelled/open", "waters removed", "H removed", "altloc residues"], rows)
    return t + "\n" + _warn_block(warns), []


def sec_units(b: Path) -> tuple[str, list[str]]:
    rows, warns = [], []
    for p in sorted((b / layout.UNITS).glob("*/unit.json")):
        u = _json(p) or {}
        pid = p.parent.name
        rows.append(
            [
                pid,
                u.get("selected") or "ERROR",
                ",".join(u.get("chains", [])),
                ",".join(u.get("targets", [])),
                ", ".join(sorted((u.get("files") or {}).values())) if isinstance(u.get("files"), dict) else "",
            ]
        )
        warns += [(pid, w) for w in u.get("warnings", [])]
    t = _table(["ID", "selected unit", "chains", "targets", "files for the matrix"], rows)
    return t + "\n" + _warn_block(warns), []


def _relax_rows(b: Path) -> tuple[list[dict], float]:
    out, gate = [], MAX_RMSD
    for p in sorted((b / layout.RELAXED).glob("*/relax.json")):
        rj = _json(p) or {}
        gate = (rj.get("thresholds") or {}).get("max_rmsd") or gate
        for fname, r in (rj.get("files") or {}).items():
            out.append({"label": f"{p.parent.name}/{fname.removesuffix('.pdb')}", **r})
    return out, gate


def sec_relax(b: Path, figdir: Path | None) -> tuple[str, list[str]]:
    rows, gate = _relax_rows(b)
    t = _table(
        ["structure/file", "status", "heavy RMSD (Å)", "ΔE (kJ/mol/atom)", "closest pocket contact (Å)", "failure"],
        [
            [
                r["label"],
                r.get("status"),
                _num(r.get("rmsd_heavy")),
                _num(r.get("drop_kj_per_atom"), 1),
                _num(r.get("pocket_min_dist")),
                r.get("failure") or "",
            ]
            for r in rows
        ],
    )
    figs = []
    scored = [r for r in rows if r.get("rmsd_heavy") is not None]
    if figdir and scored:
        fig, (a1, a2) = _figure(2, (9.6, 3.6))
        names = [r["label"] for r in scored]
        cols = [OK if r.get("status") == "ok" else BAD for r in scored]
        _style(a1, "Heavy-atom RMSD after relaxation", "Å", None, grid="x")
        a1.barh(names, [r["rmsd_heavy"] for r in scored], color=cols, height=0.6)
        a1.axvline(gate, color=MUTED, lw=1, ls="--")
        _style(a2, "Energy change", "kJ/mol per atom (negative = relaxed)", None, grid="x")
        a2.barh(names, [r.get("drop_kj_per_atom") or 0 for r in scored], color=C1, height=0.6)
        a2.set_yticks([])
        figs.append(_save(fig, figdir / "04-relax.png"))
        t += "\n" + _figure_md(
            "04-relax.png",
            "Step D relaxation gates",
            f"each bar is one relaxed file. Green passed. The dashed line is the {gate} Å RMSD gate: the "
            "structure moved too far if the bar crosses it. The right panel should be negative and modest: a "
            "very large drop means the input was strained.",
        )
    return t, figs


def sec_af(b: Path, figdir: Path | None) -> tuple[str, list[str]]:
    rows = []
    for p in sorted((b / layout.ALPHAFOLD).glob("*.json")):
        for ch, r in (_json(p) or {}).items():
            rows.append({"chain": ch, "id": p.stem, **r})
    t = _table(
        ["ID", "chain", "UniProt", "status", "TM-score", "Cα-lDDT", "core RMSD (Å)", "mean pLDDT"],
        [
            [
                r["id"],
                r["chain"],
                r.get("uniprot"),
                r.get("status"),
                _num(r.get("tm_score")),
                _num(r.get("lddt_ca")),
                _num(r.get("rmsd_core")),
                _num(r.get("mean_plddt_aligned"), 1),
            ]
            for r in rows
        ],
    )
    figs = []
    pts = [r for r in rows if r.get("tm_score") is not None and r.get("lddt_ca") is not None]
    if figdir and pts:
        fig, ax = _figure(size=(5.2, 4.4))
        _style(ax, "Experiment vs AlphaFold prediction", "TM-score", "Cα-lDDT")
        ax.scatter([r["tm_score"] for r in pts], [r["lddt_ca"] for r in pts], color=C1, s=28)
        for r in pts:
            ax.annotate(
                f"{r['id']}:{r['chain']}",
                (r["tm_score"], r["lddt_ca"]),
                fontsize=7,
                color=MUTED,
                xytext=(4, 3),
                textcoords="offset points",
            )
        ax.axvline(TM_MIN, color=GRID, lw=1.2)
        ax.axhline(LDDT_MIN, color=GRID, lw=1.2)
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.02)
        figs.append(_save(fig, figdir / "05-alphafold.png"))
        t += "\n" + _figure_md(
            "05-alphafold.png",
            "Step E agreement with AlphaFold",
            f"each point is one protein chain. Top right is close agreement. The grey lines are the TM {TM_MIN} and "
            f"lDDT {LDDT_MIN} flag thresholds; a point below or left of them is flagged. Disagreement is not an "
            "error by itself: flexible regions and bound-state conformations legitimately differ.",
        )
    return t, figs


def _systems(b: Path) -> list[dict]:
    out = []
    for p in sorted((b / layout.SYSTEMS).glob("*/system.json")):
        s = _json(p) or {}
        s["_dir"] = p.parent
        s.setdefault("id", p.parent.name)
        out.append(s)
    return out


def sec_build(b: Path, figdir: Path | None) -> tuple[str, list[str]]:
    sysl = _systems(b)
    rows, warns = [], []
    for s in sysl:
        sm = s.get("smoke") or {}
        comp = s.get("components") or {}
        rows.append(
            [
                s["id"],
                s.get("role"),
                s.get("status"),
                s.get("atoms"),
                _num(s.get("solute_charge"), 1),
                ", ".join(f"{k}x{v}" for k, v in (comp.get("ligands") or {}).items()) or "-",
                _num((sm.get("em") or {}).get("fmax"), 0),
                (sm.get("nvt") or {}).get("lincs_warnings", "-"),
                "ok" if sm.get("ok") else ("-" if not sm else "FAILED"),
            ]
        )
        warns += [(s["id"], w) for w in s.get("warnings", [])]
    t = _table(
        ["system", "role", "build", "atoms", "solute charge", "ligands", "EM Fmax", "NVT LINCS warnings", "smoke"], rows
    )
    figs = []
    sized = [s for s in sysl if s.get("atoms")]
    if figdir and sized:
        fig, ax = _figure(size=(6.4, 0.6 + 0.5 * len(sized)))
        _style(ax, "System size", "atoms (including solvent and ions)", None, grid="x")
        ax.barh([s["id"] for s in sized], [s["atoms"] for s in sized], color=C1, height=0.6)
        figs.append(_save(fig, figdir / "06-systems.png"))
        t += "\n" + _figure_md(
            "06-systems.png",
            "Step F system sizes",
            "atoms per built system. Cost of MD scales with this number, so it predicts the "
            "ns/day you will get and is the first thing to check when planning the run length.",
        )
    return t + "\n" + _warn_block(warns), figs


def _ts(x):
    try:
        return datetime.fromisoformat(x)
    except (TypeError, ValueError):
        return None


def sec_run(b: Path, figdir: Path | None) -> tuple[str, list[str]]:
    rows, bars, warns = [], [], []
    for s in _systems(b):
        r = _json(s["_dir"] / "run.json")
        if not r:
            continue
        ph = r.get("phases") or {}
        phase_status = [(ph.get(k) or {}).get("status", "-") for k in ("em", "nvt", "npt", "md")]
        rows.append([s["id"], r.get("status"), *phase_status])
        warns += [(s["id"], w) for w in r.get("warnings", [])]
        if r.get("reason"):
            warns.append((s["id"], r["reason"]))
        prev = _ts(r.get("started"))  # run.json records each phase's finish; a phase starts when the last ended
        for k in ("em", "nvt", "npt", "md"):
            p = ph.get(k) or {}
            a, z = _ts(p.get("started")) or prev, _ts(p.get("finished"))
            if a and z:
                bars.append((s["id"], k, (z - a).total_seconds() / 60))
            prev = z or prev
    t = _table(["system", "run", "EM", "NVT", "NPT", "production"], rows)
    figs = []
    if figdir and bars:
        ids = list(dict.fromkeys(x[0] for x in bars))
        fig, ax = _figure(size=(6.4, 1.8 + 0.5 * len(ids)))
        _style(ax, "Wall-clock time per MD phase", "minutes", None, grid="x")
        shades = {"em": GRID, "nvt": C1, "npt": C2, "md": INK}
        for y, sid in enumerate(ids):
            left = 0.0
            for k in ("em", "nvt", "npt", "md"):
                dur = sum(d for i, p, d in bars if i == sid and p == k)
                if dur:
                    ax.barh(y, dur, left=left, color=shades[k], height=0.6, label=k if y == 0 else None)
                    left += dur
        ax.set_yticks(range(len(ids)), ids)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK, ncols=4, loc="upper center", bbox_to_anchor=(0.5, -0.6))
        figs.append(_save(fig, figdir / "07-run.png"))
        t += "\n" + _figure_md(
            "07-run.png",
            "Step G time per phase",
            "stacked bars, one per system, in the order EM, NVT, NPT, production. Use the production "
            "bar to extrapolate how long a longer run would take on this machine.",
        )
    return t + "\n" + _warn_block(warns), figs


def _analyses(b: Path) -> list[dict]:
    out = []
    for s in _systems(b):
        a = _json(s["_dir"] / "analysis" / "analysis.json")
        if a and a.get("status") == "ok":
            out.append(a)
    return out


def sec_analysis(b: Path, figdir: Path | None) -> tuple[str, list[str]]:
    reps = _analyses(b)
    rows = []
    for a in reps:
        s = a.get("summary", {})
        top = ((a.get("hbonds") or {}).get("pairs") or [None])[0]
        rows.append(
            [
                a["id"],
                a.get("role"),
                (a.get("time_ns") or [0])[-1],
                _num((s.get("ca_rmsd") or {}).get("last_half_mean")),
                _num((s.get("target_rmsd") or {}).get("last_half_mean")),
                _num((s.get("rg") or {}).get("mean")),
                _num((s.get("hbond_count") or {}).get("mean")),
                f"{top['donor']}->{top['acceptor']} ({100 * top['occupancy']:.0f}%)" if top else "-",
                "; ".join(f.split(":")[0] for f in a.get("flags", [])) or "-",
            ]
        )
    t = _table(
        ["system", "role", "ns", "Cα RMSD (Å)", "target RMSD (Å)", "Rg (Å)", "H-bonds", "top H-bond", "flags"], rows
    )
    short = [a["id"] for a in reps if (a.get("time_ns") or [0])[-1] < SHORT_RUN_NS]
    if short:
        t += (
            f"\n> **Short trajectories** ({', '.join(short)}): under {SHORT_RUN_NS:g} ns. "
            "Equilibration and pose retention only.\n"
        )
    figs = []
    if figdir and reps:
        have_ca = [a for a in reps if a.get("ca_rmsd")]
        have_tg = [a for a in reps if a.get("target_rmsd")]
        if have_ca or have_tg:
            fig, (a1, a2) = _figure(2, (9.6, 3.6))
            _style(a1, "Protein Cα RMSD", "time (ns)", "Å")
            _style(a2, "Target RMSD (protein frame)", "time (ns)", None)
            for i, a in enumerate(have_ca):
                a1.plot(a["time_ns"], a["ca_rmsd"], lw=1.2, label=a["id"], color=[C1, C2, INK, MUTED][i % 4])
            for i, a in enumerate(have_tg):
                a2.plot(a["time_ns"], a["target_rmsd"], lw=1.2, label=a["id"], color=[C1, C2, INK, MUTED][i % 4])
            for ax in (a1, a2):
                ax.set_ylim(bottom=0)
                if ax.get_legend_handles_labels()[0]:
                    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
            figs.append(_save(fig, figdir / "08-rmsd-overlay.png"))
            t += "\n" + _figure_md(
                "08-rmsd-overlay.png",
                "Step H RMSD, all systems",
                "every system on the same axes. Compare plateau heights between systems with and without "
                "the target: a lower Cα plateau with the target suggests stabilisation, a higher target "
                "curve suggests a pose that is moving.",
            )
        cells = {
            (a.get("enzyme"), a.get("target")): (a.get("summary", {}).get("target_rmsd") or {}).get("last_half_mean")
            for a in reps
            if a.get("role") == "complex"
        }
        encs = sorted({k[0] for k in cells})
        tgs = sorted({k[1] for k in cells})
        if cells and (len(encs) > 1 or len(tgs) > 1):
            import numpy as np

            grid = np.array([[cells.get((e, g), np.nan) or np.nan for g in tgs] for e in encs], float)
            fig, ax = _figure(size=(1.6 + 1.1 * len(tgs), 1.2 + 0.7 * len(encs)))
            _style(ax, "Target RMSD, last half (Å)", "target", "receptor", grid="none")
            im = ax.imshow(grid, cmap="Blues", vmin=0)
            ax.set_xticks(range(len(tgs)), tgs)
            ax.set_yticks(range(len(encs)), encs)
            for i in range(len(encs)):
                for j in range(len(tgs)):
                    if not np.isnan(grid[i, j]):
                        ax.text(j, i, f"{grid[i, j]:.1f}", ha="center", va="center", fontsize=8, color=INK)
            fig.colorbar(im, ax=ax, shrink=0.8)
            figs.append(_save(fig, figdir / "08-matrix.png"))
            t += "\n" + _figure_md(
                "08-matrix.png",
                "Receptor x target matrix",
                "darker cells are poses that moved further from where they started. Blank cells were not simulated.",
            )
    return t, figs


# ------------------------------------------------------------------ driver
def build_report(b: Path, figures: bool = True) -> dict:
    from .status import batch_status

    st = batch_status(b)
    out = b / layout.REPORTS
    figdir = out / "figures" if figures and have_matplotlib() else None
    out.mkdir(parents=True, exist_ok=True)
    if figdir:
        figdir.mkdir(exist_ok=True)
    sections = [
        ("B", "Fix: cleaned structures", lambda: sec_fix(b)),
        ("C", "Extract: biological unit", lambda: sec_units(b)),
        ("D", "Relax: structural integrity", lambda: sec_relax(b, figdir)),
        ("E", "AlphaFold comparison", lambda: sec_af(b, figdir)),
        ("F", "Build: GROMACS systems", lambda: sec_build(b, figdir)),
        ("G", "Run: molecular dynamics", lambda: sec_run(b, figdir)),
        ("H", "Analyze: trajectories", lambda: sec_analysis(b, figdir)),
    ]
    ff = st.get("forcefield") or "not locked"
    nxt = st.get("next") or {}
    name = b.resolve().parent.name if b.resolve().name == "batch" else b.resolve().name
    md = [
        f"# Batch report: `{name}`",
        "",
        f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} from the files in this batch. "
        "Re-run `gromatrix-report` after any step; nothing here is edited by hand.",
        "",
        f"- Forcefield: **{ff}**",
        "- Steps: " + ", ".join(f"{k}={v}" for k, v in st["steps"].items()),
        f"- Next: {nxt.get('why', '-')}",
        "",
        "Folders: `01-raw` to `04-relaxed` are the structure steps, `05-alphafold*` the AlphaFold comparison, "
        "`06-params` and `06-systems` the GROMACS builds and MD, `07-*` run logs, `08-analysis` the summary "
        "table. Each system's plots and their own README are in `06-systems/<ID>/analysis/`.",
        "",
    ]
    if figures and not figdir:
        md += ["> Figures skipped: matplotlib is not installed (`pip install matplotlib`, or `pixi run -e md`).", ""]
    figs: list[str] = []
    done: list[str] = []
    for letter, title, fn in sections:
        text, f = fn()
        md += [f"## Step {letter}: {title}", "", text, ""]
        figs += f
        done.append(letter)
    path = out / "README.md"
    path.write_text("\n".join(md))
    return {"ok": True, "report": str(path), "figures": figs, "sections": done}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gromatrix-report", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--no-figures", action="store_true", help="Markdown only")
    a = ap.parse_args(argv)
    if not a.batch_dir.is_dir():
        print(json.dumps({"ok": False, "error": f"{a.batch_dir} is not a directory"}))
        return 1
    try:
        res = build_report(a.batch_dir, figures=not a.no_figures)
    except Exception as e:  # report generation must never be mistaken for a pipeline failure
        print(json.dumps({"ok": False, "error": f"{type(e).__name__}: {e}"}))
        return 1
    print(json.dumps(res))
    return 0


if __name__ == "__main__":
    sys.exit(main())
