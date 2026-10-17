"""scripts/check_plugin.py and the shape of the plugin: skills, agents, manifests, delegation rules."""

import importlib.util
import json
import re
from pathlib import Path

import pytest

pytest.importorskip("yaml")

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_plugin", ROOT / "scripts" / "check_plugin.py")
cp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cp)


def _skill(tmp_path, name, frontmatter):
    d = tmp_path / name
    (d / "references").mkdir(parents=True)
    (d / "SKILL.md").write_text(f"---\n{frontmatter}\n---\n# body\n")
    (d / "references" / "cli.md").write_text("flags\n")
    return d


def test_the_repository_is_a_valid_plugin(capsys):
    assert cp.main(["--root", str(ROOT)]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["agents"] == ["md-analyst", "md-runner"]
    assert out["skills"] == ["analysis", "execution", "gmx-builder", "start", "structure-prep", "validate"]


# The commands each stage skill may tell the agent to run (plus gromatrix-status, the hand-off).
STAGE_COMMANDS = {
    "structure-prep": {"gromatrix-fetch", "gromatrix-fix", "gromatrix-extract"},
    "validate": {"gromatrix-relax", "gromatrix-afcompare"},
    "gmx-builder": {"gromatrix-gmx"},
    "execution": {"gromatrix-run"},
    "analysis": {"gromatrix-analyze"},
}


def test_stage_skills_only_run_their_own_commands():
    """Progressive disclosure: a stage skill's code blocks never run another stage's commands."""
    for skill, own in STAGE_COMMANDS.items():
        d = ROOT / "skills" / skill
        for md in [d / "SKILL.md", *sorted((d / "references").glob("*.md"))]:
            blocks = re.findall(r"```.*?```", md.read_text(), re.S)
            used = set(re.findall(r"\bgromatrix-[a-z]+", "\n".join(blocks)))
            assert used - {"gromatrix-report", "gromatrix-compute"} <= own | {"gromatrix-status"}, f"{md.name}: {used}"


def test_orchestrator_never_runs_steps_g_and_h_itself():
    """The long and noisy steps are delegated: no command block of the orchestrator runs them."""
    text = (ROOT / "skills" / "start" / "SKILL.md").read_text()
    blocks = "\n".join(re.findall(r"```.*?```", text, re.S))
    assert "gromatrix-run" not in blocks and "gromatrix-analyze" not in blocks and "--run-md" not in blocks
    assert 'subagent_type = "gromatrix:md-runner"' in text and "gromatrix:md-analyst" in text
    assert "never run `gromatrix-run`" in text


def test_agents_hand_back_json_only_and_match_status_names():
    from gromatrix import status

    for name, const in (("md-runner", status.MD_RUNNER), ("md-analyst", status.MD_ANALYST)):
        body = (ROOT / "agents" / f"{name}.md").read_text()
        assert const == f"gromatrix:{name}" and "only that JSON" in body
    assert "-maxwarn" in (ROOT / "agents" / "md-runner.md").read_text()


def test_manifest_version_tracks_pyproject():
    import tomllib

    plugin = json.loads((ROOT / ".claude-plugin" / "plugin.json").read_text())
    py = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert plugin["version"] == py["version"] and plugin["license"] == py["license"] == "MIT"


def test_validation_catches_common_mistakes(tmp_path):
    ok = _skill(tmp_path, "good-skill", "name: good-skill\ndescription: Does a thing. Use when asked.")
    assert cp.check_skill(ok) == []
    colon = _skill(tmp_path, "colon", "name: colon\ndescription: Builds systems: topologies and boxes")
    assert "invalid YAML" in cp.check_skill(colon)[0]
    mismatch = _skill(tmp_path, "folder", "name: other\ndescription: x")
    assert any("folder name" in e for e in cp.check_skill(mismatch))
    brackets = _skill(tmp_path, "brackets", "name: brackets\ndescription: use <this>")
    assert any("angle brackets" in e for e in cp.check_skill(brackets))
    agent = tmp_path / "runner.md"
    agent.write_text("---\nname: other\ndescription: x\n---\nbody\n")
    assert any("file name" in e for e in cp.check_agent(agent)) and any("tools" in e for e in cp.check_agent(agent))
