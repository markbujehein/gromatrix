"""gmxprep-pipeline: planning, the generated step-F manifest, fail-fast."""

import json
from pathlib import Path

import pytest
import yaml

from gromatrix import pipeline
from gromatrix.gmx_prep import load_manifest, target_residues
from gromatrix.pdbio import Atom


def _manifest(tmp_path, **raw):
    raw.setdefault("batch_dir", "B")
    raw.setdefault("forcefield", "amber")
    p = tmp_path / "pipe.yaml"
    p.write_text(yaml.safe_dump(raw))
    return p


def _plan(capsys, p, *extra):
    assert pipeline.main(["--manifest", str(p), "--dry-run", *extra]) == 0
    return json.loads(capsys.readouterr().out)["plan"]


def test_to_flags():
    assert pipeline.to_flags({"ids": ["1ABC", "2DEF"], "max_loop": 8, "keep_additives": True, "x": False}) == [
        "--ids", "1ABC", "2DEF", "--max-loop", "8", "--keep-additives",
    ]  # fmt: skip


def test_dry_run_plan_from_pdb_ids(tmp_path, capsys):
    p = _manifest(tmp_path, steps={"fetch": {"ids": ["1XVP"]}, "extract": {"targets": "CID"}, "run": {"threads": 8}})
    plan = _plan(capsys, p, "--run-md", "--skip-afdb")
    b = str((tmp_path / "B").resolve())
    cmds = {(s["stage"], s["command"].split()[0 if s["stage"] != "F" else 1]): s for s in plan}
    assert [s["stage"] for s in plan] == ["lock", "A", "B", "C", "D", "E", "F", "F", "F", "G", "H"]
    assert plan[0]["command"] == f"gmxprep-gmx lock --batch-dir {b} --forcefield amber"
    assert cmds[("A", "gmxprep-fetch")]["command"].endswith("--ids 1XVP")
    assert cmds[("C", "gmxprep-extract")]["command"].endswith("--targets CID")
    assert "skip" not in cmds[("D", "gmxprep-relax")] and "skip-afdb" in cmds[("E", "gmxprep-afcompare")]["skip"]
    assert "<B/00-config/manifest.generated.yaml>" in plan[6]["command"]  # derived after C/D
    assert plan[9]["command"].endswith("--threads 8") and "skip" not in plan[9]
    assert "--analyze" in plan[10]["skip"]


def test_dry_run_plan_from_user_manifest_and_range(tmp_path, capsys):
    p = _manifest(
        tmp_path,
        enzymes=[{"id": "E", "pdb": "e.pdb"}],
        targets=[{"id": "T", "kind": "ligand", "resname": "LIG"}],
        complexes=[{"enzyme": "E", "target": "T", "pdb": "c.pdb"}],
        steps={"smoke": False},
    )
    plan = _plan(capsys, p, "--from", "F", "--to", "G", "--run-md", "--analyze")
    by = {s["stage"]: s for s in plan}
    assert all("no steps.fetch" in by[x]["skip"] for x in "ABCDE")
    build = next(s for s in plan if " build " in s["command"])
    assert str(p.resolve()) in build["command"] and "skip" not in build
    assert any("smoke" in s["command"] and "--skip-smoke" in s.get("skip", "") for s in plan)
    assert "outside --from F --to G" in by["H"]["skip"]


def test_plan_errors_are_one_json_line(tmp_path, capsys):
    p = tmp_path / "m.yaml"
    p.write_text("forcefield: amber\n")
    assert pipeline.main(["--manifest", str(p), "--dry-run"]) == 1
    assert "batch_dir" in json.loads(capsys.readouterr().out)["error"]


def _unit(b, pid, files, targets=("LIG",), chains=()):
    d = b / "03-units" / pid
    d.mkdir(parents=True)
    for f in files.values():
        (d / f).write_text("ATOM\n")
    (d / "unit.json").write_text(
        json.dumps({"id": pid, "selected": "assembly-1", "targets": list(targets), "target_chains": list(chains),
                    "files": files})
    )  # fmt: skip


def test_generated_manifest_uses_relaxed_files_and_altloc_variants(tmp_path):
    b = tmp_path / "B"
    _unit(b, "1AAA", {"enzyme": "enzyme.pdb", "target_LIG": "target_LIG.pdb", "complex_LIG": "complex_LIG.pdb"})
    _unit(
        b,
        "2BBB",
        {"enzyme": "enzyme.pdb", "target_LIG_altA": "target_LIG_altA.pdb", "complex_LIG_altA": "complex_LIG_altA.pdb",
         "target_LIG_altB": "target_LIG_altB.pdb", "complex_LIG_altB": "complex_LIG_altB.pdb"},
    )  # fmt: skip
    _unit(b, "3CCC", {"enzyme": "enzyme.pdb", "target_P": "target_P.pdb", "complex_P": "complex_P.pdb"}, (), ("P",))
    rel = b / "04-relaxed" / "1AAA"
    rel.mkdir(parents=True)
    (rel / "relax.json").write_text(
        json.dumps({"files": {"enzyme.pdb": {"status": "ok"}, "complex_LIG.pdb": {"status": "ok"}}})
    )
    rel2 = b / "04-relaxed" / "2BBB"
    rel2.mkdir(parents=True)
    (rel2 / "relax.json").write_text(
        json.dumps({"files": {"complex_LIG_altB.pdb": {"status": "failed", "failure": "explosion"}}})
    )
    (tmp_path / "lig.sdf").write_text("x")
    raw = {"forcefield": "amber", "ligands": {"LIG": {"sdf": "lig.sdf", "charge": 0}}}
    path, notes = pipeline.derive_manifest(raw, tmp_path, b)
    man = yaml.safe_load(path.read_text())

    assert [e["id"] for e in man["enzymes"]] == ["1AAA", "2BBB", "3CCC"]
    assert man["enzymes"][0]["pdb"] == str(b / "04-relaxed" / "1AAA" / "enzyme.pdb")
    assert man["enzymes"][1]["pdb"] == str(b / "03-units" / "2BBB" / "enzyme.pdb")  # not relaxed: step-C file
    cx = {(c["enzyme"], c["target"], c.get("variant")) for c in man["complexes"]}
    assert cx == {("1AAA", "LIG", None), ("2BBB", "LIG", "altA"), ("3CCC", "chains_P", None)}
    assert any("altB" in n and "explosion" in n for n in notes)  # failed relax: never used
    t = {x["id"]: x for x in man["targets"]}
    assert t["LIG"]["resname"] == "LIG" and t["LIG"]["pdb"].endswith("units/1AAA/target_LIG.pdb")
    assert t["chains_P"] == {"id": "chains_P", "kind": "protein", "chains": ["P"], "pdb": t["chains_P"]["pdb"]}
    assert man["ligands"]["LIG"]["sdf"] == str((tmp_path / "lig.sdf").resolve())

    # Step F reads it; native pairs only, so the matrix is incomplete until the user allows it.
    with pytest.raises(ValueError, match="matrix incomplete"):
        load_manifest(path, allow_incomplete=False)
    m = load_manifest(path, allow_incomplete=True)
    assert {s.id for s in m.systems if s.role == "complex"} == {"1AAA__LIG", "2BBB__LIG__altA", "3CCC__chains_P"}
    assert next(s for s in m.systems if s.id == "3CCC__chains_P").target_chains == {"P"}


def test_user_named_targets_are_kept(tmp_path):
    b = tmp_path / "B"
    _unit(b, "1AAA", {"enzyme": "enzyme.pdb", "target_CID": "target_CID.pdb", "complex_CID": "complex_CID.pdb"})
    raw = {"targets": [{"id": "citco", "kind": "ligand", "resname": "CID", "sdf": "citco.sdf"}]}
    man = yaml.safe_load(pipeline.derive_manifest(raw, tmp_path, b)[0].read_text())
    assert man["complexes"] == [{"enzyme": "1AAA", "target": "citco", "pdb": str(b / "03-units/1AAA/complex_CID.pdb")}]
    assert man["targets"][0]["sdf"] == str((tmp_path / "citco.sdf").resolve()) and "pdb" not in man["targets"][0]


def test_execute_fails_fast_and_keep_going(tmp_path, capsys, monkeypatch):
    p = _manifest(
        tmp_path,
        enzymes=[{"id": "E", "pdb": "e.pdb"}],
        targets=[{"id": "T", "kind": "ligand", "resname": "LIG"}],
        complexes=[{"enzyme": "E", "target": "T", "pdb": "c.pdb"}],
    )
    calls = []

    def fake_call(step):
        calls.append(step.label)
        return {"ok": step.label != "build", "error": "boom" if step.label == "build" else None}

    monkeypatch.setattr(pipeline, "call", fake_call)
    assert pipeline.main(["--manifest", str(p), "--run-md"]) == 1
    out = capsys.readouterr().out.splitlines()
    final = json.loads(out[-1])
    assert calls == ["lock", "build"] and final["failed_stage"] == "F"
    assert "F verify: not run" in final["stages"] and "G run: not run" in final["stages"]
    rep = json.loads((tmp_path / "B" / "00-config" / "pipeline.json").read_text())
    assert rep["steps"][1]["status"] == "skipped"  # A: no steps.fetch

    calls.clear()
    assert pipeline.main(["--manifest", str(p), "--run-md", "--keep-going"]) == 1
    assert calls == ["lock", "build", "verify", "smoke", "run"]


def test_call_parses_result_and_passes_progress(tmp_path, capsys, monkeypatch):
    import types

    mod = types.ModuleType("fake_step")

    def main(argv):
        print(json.dumps({"step": "run", "phase": "nvt", "status": "running", "progress_pct": 45}))
        print(json.dumps({"ok": True, "done": ["S1"]}))
        return 0

    mod.main = main
    monkeypatch.setitem(__import__("sys").modules, "fake_step", mod)
    out = pipeline.call(pipeline.Step("G", "gmxprep-run", "fake_step", []))
    assert out == {"ok": True, "done": ["S1"]}
    assert json.loads(capsys.readouterr().out)["progress_pct"] == 45  # forwarded live


def test_target_residues_ranges():
    def atoms(chain, n):
        return [Atom("ATOM", "CA", "ALA", chain, i, "", 0.0, 0.0, 0.0, "C") for i in range(1, n + 1)]

    prot = atoms("A", 3) + atoms("P", 2) + atoms("B", 1) + atoms("Q", 2)
    assert target_residues(prot, {"P", "Q"}) == [[3, 5], [6, 8]]
    assert target_residues(prot, None) == []


def test_examples_manifest_is_valid_and_shipped_with_the_skill():
    root = Path(__file__).resolve().parents[1]
    ex = (root / "examples" / "pipeline.yaml").read_text()
    raw = yaml.safe_load(ex)
    assert raw["batch_dir"] and "fetch" in raw["steps"]
    ref = (root / ".claude" / "skills" / "gromatrix" / "references" / "pipeline-manifest.md").read_text()
    assert ex in ref  # the skill zip carries the same template (no clone needed)


def test_result_line_is_the_json_with_ok():
    assert pipeline._result('{"ok": true, "x": 1}') == {"ok": True, "x": 1}
    assert pipeline._result('{"step": "run", "phase": "md"}') is None
    assert pipeline._result('{"ok": false, "error": "a \\"step\\" failed"}')["ok"] is False
    assert pipeline._result("not json") is None
