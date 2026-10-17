import json
from pathlib import Path

import pytest

from gromatrix import result


def test_order_puts_envelope_keys_first():
    rec = {"zzz": 1, "applied": {"a": 1}, "error": "boom", "ok": False, "aaa": 2, "code": "unexpected"}
    assert list(result.order(rec)) == ["ok", "code", "error", "applied", "zzz", "aaa"]


def test_order_does_not_mutate_input():
    rec = {"b": 1, "ok": True}
    result.order(rec)
    assert list(rec) == ["b", "ok"]


def test_fail_fills_next_action_from_registry():
    out = result.fail("unexpected", "boom")
    assert out["ok"] is False and out["code"] == "unexpected" and out["error"] == "boom"
    assert out["next_action"] == result.CODES["unexpected"].next_action


def test_fail_keeps_a_caller_supplied_next_action():
    assert result.fail("unexpected", "boom", next_action="do X")["next_action"] == "do X"


def test_fail_unknown_code_raises():
    with pytest.raises(KeyError):
        result.fail("not_a_registered_code", "boom")


def test_from_exception_without_code_is_unexpected():
    out = result.from_exception(ValueError("bad"))
    assert out["code"] == "unexpected" and out["error"] == "ValueError: bad"


def test_from_exception_uses_a_registered_class_code():
    class Custom(RuntimeError):
        code = "unexpected"

    assert result.from_exception(Custom("x"))["code"] == "unexpected"


def test_from_exception_ignores_an_unregistered_code():
    class Bad(RuntimeError):
        code = "nope"

    assert result.from_exception(Bad("x"))["code"] == "unexpected"


def test_finish_prints_one_line_and_returns_exit_status(capsys):
    assert result.finish({"zzz": 1, "ok": True}) == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 1 and list(json.loads(lines[0])) == ["ok", "zzz"]
    assert result.finish({"ok": False}) == 1


def test_ok_false_without_error_is_still_valid_json(capsys):
    result.finish({"ok": False, "code": "unexpected"})
    assert json.loads(capsys.readouterr().out) == {"ok": False, "code": "unexpected"}


def test_coded_error_default_and_per_raise_override():
    assert result.CodedError("x").code == "unexpected"
    assert result.CodedError("x", code="batch_missing").code == "batch_missing"
    assert str(result.CodedError("msg")) == "msg"


@pytest.mark.parametrize(
    "path",
    [
        "gromatrix.checks.CheckError",
        "gromatrix.ligands.LigandError",
        "gromatrix.gmx.ToolError",
        "gromatrix.forcefields.ForcefieldMismatch",
        "gromatrix.extract.ExtractionError",
        "gromatrix.topology.TopologyError",
        "gromatrix.net.HTTPError",
        "gromatrix.run.Interrupted",
    ],
)
def test_every_exception_class_has_its_own_registered_code(path):
    import importlib

    mod, name = path.rsplit(".", 1)
    cls = getattr(importlib.import_module(mod), name)
    assert issubclass(cls, result.CodedError)
    assert cls.code in result.CODES and cls.code != "unexpected"


def test_status_on_a_missing_batch_keeps_ok_error_next_and_adds_a_code(tmp_path, capsys):
    from gromatrix import status

    assert status.main(["--batch-dir", str(tmp_path / "nope")]) == 1
    out = json.loads(capsys.readouterr().out)
    assert out["ok"] is False and "does not exist" in out["error"] and "next" in out
    assert out["code"] == "batch_missing" and list(out)[:3] == ["ok", "code", "error"]


def test_fetch_with_a_malformed_id_reports_fetch_bad_id(tmp_path, capsys):
    from gromatrix import fetch

    assert fetch.main(["--batch-dir", str(tmp_path), "--ids", "not-an-id"]) == 1
    assert json.loads(capsys.readouterr().out)["code"] == "fetch_bad_id"


def test_non_integer_charge_is_check_failed():
    from gromatrix import checks

    with pytest.raises(checks.CheckError) as e:
        checks.assert_integer_charge(0.4, "system A")
    assert result.from_exception(e.value)["code"] == "check_failed"


BIG = [{"id": f"S{i}", "note": "x" * 80} for i in range(100)]  # well over BRIEF_LIMIT as JSON


def test_stub_replaces_only_listed_oversize_keys(tmp_path):
    rec = {"ok": True, "units": BIG, "other": BIG, "small": [1]}
    out = result.stub_large(rec, "extract", tmp_path)
    assert out["units"]["stub"] is True and out["units"]["chars"] == len(json.dumps(BIG))
    assert out["other"] == BIG and out["small"] == [1]


def test_stub_file_holds_the_full_value(tmp_path):
    out = result.stub_large({"ok": True, "units": BIG}, "extract", tmp_path)
    saved = json.loads(Path(out["units"]["file"]).read_text())
    assert saved == {"units": BIG} and Path(out["units"]["file"]).parent == tmp_path


def test_small_value_is_not_stubbed(tmp_path):
    rec = {"ok": True, "units": [{"id": "A"}]}
    assert result.stub_large(rec, "extract", tmp_path) == rec
    assert not list(tmp_path.iterdir())


def test_protected_keys_are_never_stubbed(tmp_path, monkeypatch):
    monkeypatch.setitem(result.STUBBABLE, "demo", ("error", "warnings", "units"))
    rec = {"ok": False, "error": "e" * 9000, "warnings": BIG, "units": BIG}
    out = result.stub_large(rec, "demo", tmp_path)
    assert out["error"] == rec["error"] and out["warnings"] == BIG and out["units"]["stub"] is True


def test_command_without_a_stub_list_is_untouched(tmp_path):
    rec = {"ok": True, "anything": BIG}
    assert result.stub_large(rec, "status", tmp_path) == rec


def test_finish_stubs_unless_full(tmp_path, capsys):
    result.finish({"ok": True, "units": BIG}, command="extract", out_dir=tmp_path)
    assert json.loads(capsys.readouterr().out)["units"]["stub"] is True
    result.finish({"ok": True, "units": BIG}, command="extract", out_dir=tmp_path, full=True)
    assert json.loads(capsys.readouterr().out)["units"] == BIG


def test_add_full_flag():
    import argparse

    ap = argparse.ArgumentParser()
    result.add_full_flag(ap)
    assert ap.parse_args(["--full"]).full is True and ap.parse_args([]).full is False


def test_receipt_marks_changed_options_and_drops_missing_values():
    out = result.receipt(
        {"a": 1, "b": 2, "c": None, "e": 7},
        {"a": 1, "b": 3, "c": 5, "d": 9, "e": None},
        {"x": 1, "y": None},
    )
    assert out == {
        "options": {
            "a": {"passed": 1, "used": 1, "changed": False},
            "b": {"passed": 2, "used": 3, "changed": True},
            "c": {"passed": None, "used": 5, "changed": False},
        },
        "facts": {"x": 1},
    }


def test_build_facts_come_only_from_the_stamps():
    from gromatrix import gmx_prep

    stamps = {
        "A": {"atoms": 100, "solute_charge": -2.0, "components": {"ions": ["NA", "CL"]}},
        "B": {"atoms": 50},
    }
    assert gmx_prep.build_facts(stamps) == {
        "atoms": {"A": 100, "B": 50},
        "solute_charge": {"A": -2.0},
        "input_ions": {"A": 2},
    }


def test_guarded_turns_a_crash_into_an_unexpected_result_line(capsys):
    @result.guarded
    def main(argv=None):
        raise PermissionError("read-only batch")

    assert main([]) == 1
    out = json.loads(capsys.readouterr().out)
    assert out["ok"] is False and out["code"] == "unexpected" and "PermissionError" in out["error"]


def test_guarded_passes_the_return_value_through(capsys):
    @result.guarded
    def main(argv=None):
        return 2

    assert main([]) == 2 and capsys.readouterr().out == ""


def test_fetch_survives_a_malformed_reply(tmp_path, monkeypatch, capsys):
    from gromatrix import fetch

    def boom(*a, **k):
        raise ValueError("Expecting value: line 1 column 1")

    monkeypatch.setattr(fetch, "search", boom)
    assert fetch.main(["--batch-dir", str(tmp_path), "--query", "x"]) == 1
    assert json.loads(capsys.readouterr().out)["code"] == "unexpected"


def test_a_failed_stub_write_keeps_the_full_value(tmp_path):
    blocker = tmp_path / "not-a-dir"
    blocker.write_text("")
    rec = {"ok": True, "units": BIG}
    assert result.stub_large(rec, "extract", blocker) == rec


def test_each_stub_call_writes_its_own_file(tmp_path):
    a = result.stub_large({"ok": True, "units": BIG}, "extract", tmp_path)["units"]["file"]
    b = result.stub_large({"ok": True, "units": BIG}, "extract", tmp_path)["units"]["file"]
    assert a != b and Path(a).exists() and Path(b).exists()
