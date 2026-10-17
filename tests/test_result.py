import json

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
