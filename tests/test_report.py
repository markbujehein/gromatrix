"""gromatrix-report and the layout module: offline, built from hand-written step artefacts."""

import json

from gromatrix import layout, report


def _w(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def _batch(tmp_path):
    b = tmp_path / "B"
    _w(b / layout.LOCK, {"key": "amber99sb-ildn"})
    _w(b / layout.RELAXED / "1ABC" / "relax.json", {
        "thresholds": {"max_rmsd": 2.5},
        "files": {"enzyme.pdb": {"status": "ok", "rmsd_heavy": 1.1, "drop_kj_per_atom": -30.0, "pocket_min_dist": 2.6}},
    })
    _w(b / layout.ALPHAFOLD / "1ABC.json", {"A": {"status": "ok", "uniprot": "P1", "tm_score": 0.9, "lddt_ca": 0.8}})
    _w(b / layout.SYSTEMS / "E1" / "system.json", {
        "id": "E1", "role": "enzyme", "status": "ok", "atoms": 5000, "solute_charge": -2.0,
        "components": {"ligands": {"LIG": 1}}, "smoke": {"ok": True, "em": {"fmax": 900.0}, "nvt": {"lincs_warnings": 0}},
    })
    _w(b / layout.SYSTEMS / "E1" / "run.json", {
        "status": "done", "started": "2026-01-01T00:00:00+00:00",
        "phases": {"em": {"status": "done", "finished": "2026-01-01T00:01:00+00:00"},
                   "md": {"status": "done", "finished": "2026-01-01T00:11:00+00:00"}},
    })
    _w(b / layout.SYSTEMS / "E1" / "analysis" / "analysis.json", {
        "id": "E1", "status": "ok", "role": "enzyme", "enzyme": "E1", "frames": 3, "time_ns": [0.0, 0.5, 1.0],
        "ca_rmsd": [0.0, 1.0, 1.2], "rg": [20.0, 20.1, 20.0], "plots": ["rmsd.png", "rg.png"], "flags": [],
        "settings": {"stride": 1, "hbond": {"d_a_A": 3.0, "d_h_a_deg": 150, "min_occupancy": 0.1}},
        "summary": {"ca_rmsd": {"mean": 0.7, "std": 0.5, "max": 1.2, "final": 1.2, "last_half_mean": 1.1},
                    "rg": {"mean": 20.0, "std": 0.05, "max": 20.1, "final": 20.0, "last_half_mean": 20.0}},
    })
    return b


def test_layout_names_sort_in_pipeline_order():
    top = [layout.CONFIG, layout.RAW, layout.FIXED, layout.UNITS, layout.RELAXED, layout.ALPHAFOLD,
           layout.PARAMS, layout.RUN_LOGS, layout.ANALYSIS, layout.REPORTS]
    assert top == sorted(top)
    assert all(n[:2].isdigit() for n in top + [layout.SYSTEMS, layout.SLURM, layout.ALPHAFOLD_CACHE])
    assert layout.LOCK.startswith(layout.CONFIG + "/") and layout.system_dir("B", "E1").as_posix() == "B/06-systems/E1"


def test_report_text_covers_every_step(tmp_path):
    b = _batch(tmp_path)
    res = report.build_report(b, figures=False)
    assert res["ok"] and res["figures"] == [] and res["sections"] == list("BCDEFGH")
    md = (b / layout.REPORTS / "README.md").read_text()
    for needle in ("amber99sb-ildn", "1ABC/enzyme", "P1", "E1", "5000", "Short trajectories"):
        assert needle in md


def test_analysis_readme_explains_each_plot_and_warns_when_short():
    a = json.loads(json.dumps({
        "id": "S", "role": "complex", "enzyme": "E", "target": "T", "frames": 2, "time_ns": [0.0, 0.1],
        "plots": ["rmsd.png", "hbonds.png"], "flags": ["ligand_drift: x"],
        "settings": {"stride": 2, "hbond": {"d_a_A": 3.5, "d_h_a_deg": 150, "min_occupancy": 0.1}}, "summary": {},
    }))
    text = report.analysis_readme(a)
    assert "### `rmsd.png`" in text and "### `hbonds.png`" in text and "rg.png" not in text
    assert "Short trajectory" in text and "ligand_drift" in text


def test_report_cli_reports_missing_batch(tmp_path, capsys):
    assert report.main(["--batch-dir", str(tmp_path / "nope")]) == 1
    assert json.loads(capsys.readouterr().out)["ok"] is False


def test_compute_choice_round_trips_and_cpu_only_has_no_gpu_flag(tmp_path, capsys):
    from gromatrix import compute

    b = str(tmp_path)
    assert compute.main(["--batch-dir", b, "run-args"]) == 1  # nothing chosen yet
    capsys.readouterr()
    assert compute.main(["--batch-dir", b, "set", "--mode", "local", "--threads", "6"]) == 0
    capsys.readouterr()
    compute.main(["--batch-dir", b, "run-args"])
    out = json.loads(capsys.readouterr().out)
    assert out["args"] == ["--threads", "6", "--background"] and out["submit"] is False
    compute.main(["--batch-dir", b, "set", "--mode", "remote", "--time", "48:00:00", "--partition", "gpu"])
    capsys.readouterr()
    compute.main(["--batch-dir", b, "run-args"])
    out = json.loads(capsys.readouterr().out)
    assert out["args"][:3] == ["--slurm", "--time", "48:00:00"] and "--partition" in out["args"] and out["submit"] is False
