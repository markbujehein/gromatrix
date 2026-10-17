"""Plugin hook scripts: JSON on stdin, exit code and stderr/stdout out. Hooks are a second layer; the CLI
gates stay authoritative, so these tests only pin the hook contract."""

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from gromatrix import layout
from gromatrix.watch import main as watch_main

ROOT = Path(__file__).resolve().parents[1]
HOOKS = ROOT / "scripts" / "hooks"


def run_hook(script, payload, path_env=None):
    env = None
    if path_env is not None:
        env = {"PATH": str(path_env), "PYTHONIOENCODING": "utf-8"}
    return subprocess.run(
        [sys.executable, str(HOOKS / script)],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
    )


def bash(cmd, **extra):
    return {"tool_name": "Bash", "tool_input": {"command": cmd}, **extra}


@pytest.mark.parametrize(
    "cmd,agent_id,denied",
    [
        ("gromatrix-run --batch-dir B", None, True),
        ("gromatrix-run --batch-dir B", "", True),
        ("gromatrix-run --batch-dir B", "abc123", False),
        ("pixi run -e md gromatrix-analyze --batch-dir B", None, True),
        ("cd B && gromatrix-analyze --batch-dir B", None, True),
        ("echo gromatrix-run is for subagents", None, False),
        ("gromatrix-status --batch-dir B", None, False),
    ],
)
def test_subagent_only(cmd, agent_id, denied):
    extra = {} if agent_id is None else {"agent_id": agent_id}
    r = run_hook("subagent_only.py", bash(cmd, **extra))
    assert (r.returncode == 2) == denied
    if denied:
        assert "gromatrix:md-" in r.stderr


@pytest.mark.parametrize(
    "cmd",
    [
        "gmx grompp -f em.mdp -maxwarn 1",
        "gmx grompp -f em.mdp -maxwarn=2",
        "gromatrix-gmx build --maxwarn 3",
        "gmx mdrun -deffnm md",
        "pixi run -e md gmx mdrun -deffnm md",
        "rm -rf B/06-systems/E1",
        "rm -fr 03-units",
        "rm -r --force ../b/00-config/",
    ],
)
def test_bash_guard_denies(cmd):
    r = run_hook("bash_guard.py", bash(cmd, cwd=str(ROOT)))
    assert r.returncode == 2 and r.stderr.strip()


@pytest.mark.parametrize(
    "cmd",
    [
        "gmx grompp -f em.mdp -maxwarn 0",
        "gmx grompp -f em.mdp",
        "gmx check -f md.xtc",
        "rm -rf build/tmp",
        "rm B/00-config/events.jsonl",
        "echo gmx mdrun",
        "gromatrix-gmx build --batch-dir nonexistent --force",
    ],
)
def test_bash_guard_allows(cmd):
    assert run_hook("bash_guard.py", bash(cmd, cwd=str(ROOT))).returncode == 0


def test_layout_regex_covers_every_layout_folder():
    folders = [v for v in vars(layout).values() if isinstance(v, str) and re.match(r"^0\d-[\w-]+$", v)]
    sys.path.insert(0, str(HOOKS))
    try:
        from _common import LAYOUT_DIR
    finally:
        sys.path.remove(str(HOOKS))
    assert len(folders) == 13
    assert all(LAYOUT_DIR.search(f"B/{f}/x") for f in folders)


def test_forced_build_guard(tmp_path):
    (tmp_path / "06-systems" / "E1").mkdir(parents=True)
    (tmp_path / "06-systems" / "E2").mkdir()
    (tmp_path / "06-systems" / "E1" / "run.json").write_text("{}")
    base = "gromatrix-gmx build --batch-dir B --manifest m.yaml --force"
    payload = lambda c: bash(c, cwd=str(tmp_path))  # noqa: E731
    assert run_hook("bash_guard.py", payload(base.replace("B", "."))).returncode == 2
    assert run_hook("bash_guard.py", payload(base.replace("B", ".") + " --only E2")).returncode == 0
    assert run_hook("bash_guard.py", payload(base.replace("B", ".") + " --only E1")).returncode == 2
    assert run_hook("bash_guard.py", payload(base.replace("B", ".").replace(" --force", ""))).returncode == 0


@pytest.mark.parametrize(
    "tool_input,denied",
    [
        ({"file_path": "/b/06-systems/E1/md.xtc"}, True),
        ({"file_path": "/b/06-systems/E1/md.edr"}, True),
        ({"file_path": "/b/06-systems/E1/md.log"}, True),
        ({"file_path": "/b/06-systems/E1/md.log", "limit": 40}, False),
        ({"file_path": "/b/06-systems/E1/run.json"}, False),
        ({"command": "cat 06-systems/E1/md.log"}, True),
        ({"command": "less run.log"}, True),
        ({"command": "head -c 100 md.trr"}, True),
        ({"command": "tail -n 40 06-systems/E1/md.log"}, False),
        ({"command": "gmx check -f md.xtc"}, False),
        ({"command": "ls -l md.xtc"}, False),
    ],
)
def test_read_guard(tool_input, denied):
    tool = "Bash" if "command" in tool_input else "Read"
    r = run_hook("read_guard.py", {"tool_name": tool, "tool_input": tool_input})
    assert (r.returncode == 2) == denied
    if denied:
        assert "gromatrix-status" in r.stderr


RUNNER_OK = {"ok": True, "step": "G", "mode": "local", "systems": {"E1": "done"}, "failed": []}
ANALYST_OK = {"ok": True, "step": "H", "analysed": [], "skipped": [], "failed": [], "summary_csv": "a", "report": "b"}


@pytest.mark.parametrize(
    "agent,msg,blocked",
    [
        ("gromatrix:md-runner", json.dumps(RUNNER_OK), False),
        ("gromatrix:md-runner", "```json\n" + json.dumps(RUNNER_OK) + "\n```", False),
        ("gromatrix:md-runner", json.dumps({"ok": False, "reason": "no compute choice"}), False),
        ("gromatrix:md-runner", json.dumps({"ok": False}), True),
        ("gromatrix:md-runner", "Done! " + json.dumps(RUNNER_OK), True),
        ("gromatrix:md-runner", json.dumps({"ok": True, "step": "G"}), True),
        ("gromatrix:md-runner", json.dumps({**RUNNER_OK, "step": "H"}), True),
        ("gromatrix:md-runner", "[1, 2]", True),
        ("gromatrix:md-analyst", json.dumps(ANALYST_OK), False),
        ("gromatrix:md-analyst", json.dumps(RUNNER_OK), True),
        ("other:agent", "free text", False),
    ],
)
def test_subagent_payload(agent, msg, blocked):
    r = run_hook("subagent_payload.py", {"agent_type": agent, "last_assistant_message": msg})
    assert r.returncode == 0
    if blocked:
        out = json.loads(r.stdout)
        assert out["decision"] == "block" and agent in out["reason"]
    else:
        assert r.stdout == ""


def test_subagent_payload_honours_stop_hook_active():
    r = run_hook(
        "subagent_payload.py",
        {"agent_type": "gromatrix:md-runner", "last_assistant_message": "prose", "stop_hook_active": True},
    )
    assert r.returncode == 0 and r.stdout == ""


def test_session_status_injects_one_line(tmp_path):
    (tmp_path / "00-config").mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake = bin_dir / "gromatrix-status"
    fake.write_text('#!/bin/sh\necho \'{"ok": true, "next": {"skill": "gromatrix:start"}}\'\n')
    fake.chmod(0o755)
    r = run_hook("session_status.py", {"cwd": str(tmp_path), "source": "compact"}, path_env=f"{bin_dir}:/usr/bin:/bin")
    out = json.loads(r.stdout)["hookSpecificOutput"]
    assert out["hookEventName"] == "SessionStart" and "gromatrix:start" in out["additionalContext"]
    assert len(out["additionalContext"]) < 10000


def test_session_status_noops(tmp_path):
    # not a batch directory
    assert run_hook("session_status.py", {"cwd": str(tmp_path)}).stdout == ""
    # a batch directory but the command is not installed
    (tmp_path / "00-config").mkdir()
    r = run_hook("session_status.py", {"cwd": str(tmp_path)}, path_env="/nonexistent")
    assert r.returncode == 0 and r.stdout == ""


def test_log_event_appends_and_never_fails(tmp_path):
    (tmp_path / "00-config").mkdir()
    payload = bash(
        f"gromatrix-fix --batch-dir {tmp_path}",
        cwd="/",
        tool_response={"exit_code": 0},
        duration_ms=1200,
    )
    assert run_hook("log_event.py", payload).returncode == 0
    assert run_hook("log_event.py", payload).returncode == 0
    lines = (tmp_path / "00-config" / "events.jsonl").read_text().splitlines()
    rec = json.loads(lines[0])
    assert len(lines) == 2 and rec["exit"] == 0 and rec["duration_ms"] == 1200 and "ts" in rec
    # no batch, garbage input: still exit 0 and nothing created
    other = tmp_path / "other"
    other.mkdir()
    assert run_hook("log_event.py", bash("gromatrix-fix --batch-dir x", cwd=str(other))).returncode == 0
    assert list(other.iterdir()) == []
    r = subprocess.run([sys.executable, str(HOOKS / "log_event.py")], input="not json", capture_output=True, text=True)
    assert r.returncode == 0


def test_hooks_json_scripts_exist():
    hooks = json.loads((ROOT / "hooks" / "hooks.json").read_text())["hooks"]
    commands = [h["command"] for g in hooks.values() for grp in g for h in grp["hooks"]]
    assert commands
    for c in commands:
        rel = re.search(r"\$\{CLAUDE_PLUGIN_ROOT\}/([\w./-]+)", c).group(1)
        assert (ROOT / rel).is_file()


def _w(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def test_watch_prints_phase_changes_and_exits(tmp_path, capsys):
    sysd = tmp_path / "06-systems"
    _w(sysd / "E1" / "system.json", {"status": "ok"})
    _w(
        sysd / "E1" / "run.json",
        {"status": "done", "phases": {p: {"status": "done"} for p in ("em", "nvt", "npt", "md")}},
    )
    _w(sysd / "E2" / "system.json", {"status": "ok"})
    _w(
        sysd / "E2" / "run.json",
        {"status": "failed", "phases": {"em": {"status": "done"}, "nvt": {"status": "failed"}}},
    )
    code = watch_main(["--batch-dir", str(tmp_path), "--interval", "0.01"])
    lines = [json.loads(x) for x in capsys.readouterr().out.splitlines()]
    assert code == 1  # E2 failed; the loop still ended because every system is final
    assert len(lines) == 6
    assert {"system": "E2", "phase": "nvt", "status": "failed", "run": "failed"} in lines


def test_watch_all_done_and_empty(tmp_path, capsys):
    assert watch_main(["--batch-dir", str(tmp_path), "--once"]) == 2
    _w(tmp_path / "06-systems" / "E1" / "system.json", {"status": "ok"})
    _w(tmp_path / "06-systems" / "E1" / "run.json", {"status": "done", "phases": {"em": {"status": "done"}}})
    capsys.readouterr()
    assert watch_main(["--batch-dir", str(tmp_path)]) == 0


def test_watch_falls_back_to_phase_logs(tmp_path, capsys):
    wd = tmp_path / "06-systems" / "E1"
    _w(wd / "system.json", {"status": "ok"})
    _w(wd / "run.json", {"status": "failed"})
    (wd / "em.log").write_text("Fatal error:\nbad\n")
    watch_main(["--batch-dir", str(tmp_path), "--once"])
    assert json.loads(capsys.readouterr().out.splitlines()[0])["status"] == "failed"
