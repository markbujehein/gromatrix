"""The status command is the hand-off contract between the agent skills."""

import json

from gromatrix.status import batch_status


def _w(path, data=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data) if not isinstance(data, str) else data)


def test_next_skill_follows_pipeline_state(tmp_path):
    b = tmp_path
    assert batch_status(b)["next"]["skill"] == "gromatrix"  # decide + lock forcefield
    _w(b / "forcefield.lock.json", {"key": "amber99sb-ildn"})
    _w(b / "raw" / "1ABC.cif", "data_1ABC\n")
    assert batch_status(b)["next"]["skill"] == "gromatrix-structure-prep"  # fix
    _w(b / "fixed" / "1ABC.meta.json", {})
    _w(
        b / "units" / "1ABC" / "unit.json",
        {"selected": "assembly-1", "targets": ["LIG"], "files": {"enzyme": "enzyme.pdb"}},
    )
    st = batch_status(b)
    assert st["structures"]["1ABC"]["unit"] == "assembly-1" and st["next"]["skill"] == "gromatrix-validate"
    _w(b / "relaxed" / "1ABC" / "relax.json", {"files": {"enzyme.pdb": {"status": "ok"}}})
    _w(b / "af" / "1ABC.json", {"A": {"tm_score": 0.95, "status": "ok"}})
    assert batch_status(b)["next"]["skill"] == "gromatrix-gmx-builder"  # build
    _w(b / "systems" / "E1" / "system.json", {"status": "ok"})
    assert batch_status(b)["next"]["command"].startswith("gmxprep-gmx smoke")
    _w(b / "systems" / "E1" / "system.json", {"status": "ok", "smoke": {"ok": True}})
    st = batch_status(b)
    assert st["systems"]["smoke_ok"] == 1 and st["next"]["skill"] == "gromatrix-execution"  # G
    assert st["steps"] == dict(zip("ABCDEFGH", ["done"] * 6 + ["pending"] * 2))
    import os
    import socket

    _w(b / "systems" / "E1" / "run.json", {"status": "running", "host": socket.gethostname(), "pid": os.getpid()})
    assert "running" in batch_status(b)["next"]["why"]
    _w(b / "systems" / "E1" / "run.json", {"status": "failed"})
    assert batch_status(b)["next"]["command"].endswith("--only E1")
    _w(b / "systems" / "E1" / "run.json", {"status": "done"})
    st = batch_status(b)
    assert st["next"]["skill"] == "gromatrix-analysis" and st["steps"]["G"] == "done"  # H
    _w(b / "systems" / "E1" / "analysis" / "analysis.json", {"status": "ok", "flags": ["ligand_drift: ..."]})
    final = batch_status(b)
    assert final["next"]["skill"] is None and final["steps"]["H"] == "done"
    assert final["systems"]["analysis_flags"] == {"E1": ["ligand_drift: ..."]}


def test_own_pdbs_batch_skips_structure_steps(tmp_path):
    """Systems built from the user's own docked PDBs: A-E are n/a, not 'pending'."""
    _w(tmp_path / "forcefield.lock.json", {"key": "amber99sb-ildn"})
    _w(tmp_path / "systems" / "X" / "system.json", {"status": "ok", "smoke": {"ok": True}})
    st = batch_status(tmp_path)
    assert st["steps"]["A"] == "n/a" and st["next"]["skill"] == "gromatrix-execution"


def test_dead_runner_and_smoke_failures(tmp_path):
    import os
    import socket

    _w(tmp_path / "forcefield.lock.json", {"key": "amber99sb-ildn"})
    _w(tmp_path / "systems" / "A" / "system.json", {"status": "ok", "smoke": {"ok": True}})
    _w(tmp_path / "systems" / "B" / "system.json", {"status": "ok", "smoke": {"ok": False}})
    _w(tmp_path / "systems" / "A" / "run.json", {"status": "done"})
    st = batch_status(tmp_path)
    assert st["steps"]["G"] == "done"  # B is never run by G: its smoke test failed
    # A runner killed with -9 leaves 'running' behind; status must not wait for it forever.
    _w(tmp_path / "systems" / "A" / "run.json", {"status": "running", "host": socket.gethostname(), "pid": 2**22 + 7})
    st = batch_status(tmp_path)
    assert st["systems"]["run"] == {"interrupted": ["A"]} and st["next"]["command"].startswith("gmxprep-run")
    _w(tmp_path / "systems" / "A" / "run.json", {"status": "running", "host": socket.gethostname(), "pid": os.getpid()})
    assert "running" in batch_status(tmp_path)["next"]["why"]
