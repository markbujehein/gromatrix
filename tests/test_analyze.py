"""Step H on a synthetic trajectory with known answers: a rigid 6-residue protein, a
ligand that H-bonds residue 3 for the first half and leaves the pocket for the second."""

import argparse
import csv
import json

import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")

from gromatrix import analyze, checks  # noqa: E402

NRES, NFRAMES = 6, 10
NAMES = ["N", "H", "CA", "C", "O"]


def _universe():
    n_prot = NRES * len(NAMES)
    n = n_prot + 3
    u = mda.Universe.empty(n, n_residues=NRES + 1, atom_resindex=[i // 5 for i in range(n_prot)] + [NRES] * 3,
                           trajectory=True)  # fmt: skip
    u.add_TopologyAttr("name", NAMES * NRES + ["C1", "O1", "H1"])
    u.add_TopologyAttr("resname", ["ALA"] * NRES + ["LIG"])
    u.add_TopologyAttr("resid", list(range(1, NRES + 2)))
    xyz = np.zeros((n, 3))
    for r in range(NRES):  # residue r along x, a zig-zag in y so Rg is not degenerate
        base = np.array([3.8 * r, 1.5 * (r % 2), 0.0])
        xyz[5 * r : 5 * r + 5] = base + np.array([[0, 0, 0], [0, -1.0, 0], [1.2, 0, 0], [2.4, 0, 0], [2.4, 1.2, 0]])
    o3 = xyz[5 * 2 + 4]  # residue 3 carbonyl O
    xyz[n_prot : n_prot + 3] = o3 + np.array([[0, 4.0, 0], [0, 2.8, 0], [0, 1.85, 0]])  # C1, O1-H1...O
    return u, xyz, n_prot


@pytest.fixture
def traj(tmp_path):
    u, xyz, n_prot = _universe()
    gro, xtc = tmp_path / "solute.gro", tmp_path / "fit.xtc"
    u.atoms.positions = xyz
    u.dimensions = [60, 60, 60, 90, 90, 90]
    u.atoms.write(str(gro))
    rng = np.random.default_rng(0)
    with mda.Writer(str(xtc), u.atoms.n_atoms) as w:
        for f in range(NFRAMES):
            pos = xyz + rng.normal(0, 0.05, xyz.shape)
            if f >= NFRAMES // 2:
                pos[n_prot:] += [0, 6.0, 0]  # ligand leaves: no H-bond, ~6 Å RMSD
            u.atoms.positions = pos
            u.trajectory.ts.time = 100.0 * f  # ps
            w.write(u.atoms)
    return gro, xtc, n_prot


def _args(**kw):
    d = dict(
        hbond_distance=3.0, hbond_angle=150.0, hbond_min_occupancy=0.1, ligand_rmsd_max=3.0, ca_rmsd_max=4.0,
        stride=1, begin=None, deffnm="md", gmx="gmx", force=False, only=None,
    )  # fmt: skip
    d.update(kw)
    return argparse.Namespace(**d)


def test_metrics_on_known_trajectory(traj):
    gro, xtc, n_prot = traj
    res = analyze.compute_metrics(gro, xtc, f"index {n_prot}:{n_prot + 2}", f"index 0:{n_prot - 1}", _args())
    assert res["frames"] == NFRAMES and res["time_ns"][-1] == pytest.approx(0.9)
    assert max(res["ca_rmsd"]) < 0.3  # rigid protein, 0.05 Å noise
    lig = np.array(res["target_rmsd"])
    assert lig[: NFRAMES // 2].max() < 0.5 and lig[NFRAMES // 2 :].min() > 5.0
    assert len(res["rmsf"]["rmsf"]) == NRES and max(res["rmsf"]["rmsf"]) < 0.3
    hb = res["hbonds"]
    assert hb["count"][: NFRAMES // 2] == [1] * 5 and hb["count"][NFRAMES // 2 :] == [0] * 5
    (pair,) = hb["pairs"]
    assert pair == {"donor": "LIG7:O1", "acceptor": "ALA3:O", "occupancy": 0.5, "target_is_donor": True}
    assert res["summary"]["target_rmsd"]["last_half_mean"] > 5.0
    flags = analyze.flags(res, "complex", _args())
    assert [f.split(":")[0] for f in flags] == ["ligand_drift"]


def test_selections_from_index_groups():
    groups = {"Protein": list(range(1, 31)), "LIG_LIG": [31, 32, 33], "Solute": list(range(1, 34))}
    assert analyze.index_sel([1, 2, 3, 7, 8]) == "index 0:2 or index 6:7"
    lig = {"role": "complex", "target_resname": "LIG"}
    assert analyze.selections(lig, groups) == ("index 30:32", "index 0:29")
    assert analyze.selections({"role": "target", "target_resname": "LIG"}, groups) == ("index 30:32", None)
    assert analyze.selections({"role": "enzyme"}, groups) == (None, "index 0:29")
    pep = {"role": "complex", "target_residues": [[4, 6]]}
    assert analyze.selections(pep, groups) == (
        "(index 0:29) and (resindex 4:5)",
        "(index 0:29) and not (resindex 4:5)",
    )
    assert analyze.pick_groups({"Protein": [1], "LIG_X": [2]}, "X") == ("Protein", "Backbone")
    assert analyze.pick_groups({"LIG_X": [2]}, "X") == ("LIG_X", "LIG_X")


def test_cli_end_to_end_with_prepared_trajectory(tmp_path, traj, monkeypatch, capsys):
    gro, xtc, n_prot = traj
    b = tmp_path / "B"
    wd = b / "systems" / "E1__LIG"
    wd.mkdir(parents=True)
    (wd / "system.json").write_text(
        json.dumps({"id": "E1__LIG", "role": "complex", "enzyme": "E1", "target": "LIG", "status": "ok",
                    "target_resname": "LIG", "target_residues": [],
                    "protein_residue_labels": [f"A:{100 + i}" for i in range(NRES - 1)] + ["A:104A"]})
    )  # fmt: skip
    (wd / "run.json").write_text(json.dumps({"status": "done"}))
    (wd / "md.xtc").write_text("x")
    (wd / "solv_ions.gro").write_text(gro.read_text())
    sol = list(range(1, n_prot + 4))
    checks.write_ndx(
        {"System": sol + [n_prot + 4], "Protein": sol[:n_prot], "LIG_LIG": sol[n_prot:], "Solute": sol,
         "Solvent": [n_prot + 4]},
        wd / "index.ndx",
    )  # fmt: skip
    groups = analyze.solute_groups(wd)
    assert "System" not in groups and groups["C-alpha"] == [3 + 5 * r for r in range(NRES)]
    assert len(groups["Backbone"]) == 3 * NRES

    def fake_prepare(wd_, a, groups_, tres):
        (wd_ / "analysis").mkdir(exist_ok=True)
        return gro, xtc

    monkeypatch.setattr(analyze, "prepare_trajectory", fake_prepare)
    assert analyze.main(["--batch-dir", str(b)]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["analysed"] == 1 and "ligand_drift" in out["flags"]["E1__LIG"][0]
    rep = json.loads((wd / "analysis" / "analysis.json").read_text())
    assert rep["selections"] == {"target": f"index {n_prot}:{n_prot + 2}", "enzyme": f"index 0:{n_prot - 1}"}
    assert set(rep["plots"]) == {"rmsd.png", "rmsf.png", "rg.png", "hbonds.png"}
    assert rep["rmsf"]["label"][-1] == "A:104A"  # author labels, insertion code kept
    for p in rep["plots"]:
        assert (wd / "analysis" / p).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    (row,) = csv.DictReader(open(b / "analysis_summary.csv"))
    assert row["id"] == "E1__LIG" and row["top_hbond"] == "LIG7:O1->ALA3:O" and row["flags"] == "ligand_drift"

    # An up-to-date analysis is reused; changed settings recompute it.
    calls = []
    monkeypatch.setattr(analyze, "prepare_trajectory", lambda *x: calls.append(1) or fake_prepare(*x))
    assert analyze.main(["--batch-dir", str(b)]) == 0 and calls == []
    assert analyze.main(["--batch-dir", str(b), "--ligand-rmsd-max", "8"]) == 0 and calls == [1]
    capsys.readouterr()
    assert not json.loads((wd / "analysis" / "analysis.json").read_text())["flags"]

    # --only on another system keeps this system's row in the summary.
    other = b / "systems" / "E1"
    other.mkdir()
    (other / "system.json").write_text(json.dumps({"id": "E1", "role": "enzyme", "status": "ok"}))
    (other / "run.json").write_text(json.dumps({"status": "done"}))
    (other / "md.xtc").write_text("x")
    monkeypatch.setattr(analyze, "solute_groups", lambda wd_: {"Protein": [1]})
    monkeypatch.setattr(analyze, "compute_metrics", lambda *x: (_ for _ in ()).throw(RuntimeError("numpy says no")))
    assert analyze.main(["--batch-dir", str(b), "--only", "E1"]) == 1  # unexpected error: reported, not raised
    out = json.loads(capsys.readouterr().out)
    assert out["failed"][0]["reason"].startswith("RuntimeError") and out["summary_rows"] == 1
    assert [r["id"] for r in csv.DictReader(open(b / "analysis_summary.csv"))] == ["E1__LIG"]

    # Unfinished runs are skipped.
    (wd / "run.json").write_text(json.dumps({"status": "running"}))
    assert analyze.main(["--batch-dir", str(b), "--only", "E1__LIG"]) == 1
    assert json.loads(capsys.readouterr().out)["skipped"][0]["id"] == "E1__LIG"
