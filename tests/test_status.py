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
    _w(b / "systems" / "E1" / "system.json", {"status": "ok", "smoke": {"ok": True}})
    final = batch_status(b)
    assert final["systems"]["smoke_ok"] == 1 and final["next"]["skill"] is None
