"""scripts/package_skills.py: validation rules and reproducible archives."""

import importlib.util
import re
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("package_skills", ROOT / "scripts" / "package_skills.py")
ps = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ps)


def _skill(tmp_path, name, frontmatter):
    d = tmp_path / name
    (d / "references").mkdir(parents=True)
    (d / "SKILL.md").write_text(f"---\n{frontmatter}\n---\n# body\n")
    (d / "references" / "cli.md").write_text("flags\n")
    return d


def test_repository_skills_are_valid():
    skills = [p for p in (ROOT / ".claude" / "skills").iterdir() if (p / "SKILL.md").is_file()]
    assert len(skills) == 6 and all(ps.validate(s) == [] for s in skills)


# The commands each stage skill may tell the agent to run (plus gmxprep-status, the hand-off).
STAGE_COMMANDS = {
    "gromatrix-structure-prep": {"gmxprep-fetch", "gmxprep-fix", "gmxprep-extract"},
    "gromatrix-validate": {"gmxprep-relax", "gmxprep-afcompare"},
    "gromatrix-gmx-builder": {"gmxprep-gmx"},
    "gromatrix-execution": {"gmxprep-run"},
    "gromatrix-analysis": {"gmxprep-analyze"},
}


def test_stage_skills_only_run_their_own_commands():
    """Progressive disclosure: a stage skill's code blocks never run another stage's commands."""
    for skill, own in STAGE_COMMANDS.items():
        d = ROOT / ".claude" / "skills" / skill
        for md in [d / "SKILL.md", *sorted((d / "references").glob("*.md"))]:
            blocks = re.findall(r"```.*?```", md.read_text(), re.S)
            used = set(re.findall(r"\bgmxprep-[a-z]+", "\n".join(blocks)))
            assert used <= own | {"gmxprep-status"}, f"{md.relative_to(ROOT)} runs {used - own}"


def test_validation_catches_common_mistakes(tmp_path):
    ok = _skill(tmp_path, "good-skill", "name: good-skill\ndescription: Does a thing. Use when asked.")
    assert ps.validate(ok) == []
    colon = _skill(tmp_path, "colon", "name: colon\ndescription: Builds systems: topologies and boxes")
    assert "invalid YAML" in ps.validate(colon)[0]
    mismatch = _skill(tmp_path, "folder", "name: other\ndescription: x")
    assert any("folder name" in e for e in ps.validate(mismatch))
    brackets = _skill(tmp_path, "brackets", "name: brackets\ndescription: use <this>")
    assert any("angle brackets" in e for e in ps.validate(brackets))


def test_archives_are_reproducible_and_upload_shaped(tmp_path):
    s = _skill(tmp_path / "src", "my-skill", "name: my-skill\ndescription: Does a thing.")
    a = ps.package(s, tmp_path).read_bytes()
    b = ps.package(s, tmp_path).read_bytes()
    assert a == b
    names = zipfile.ZipFile(tmp_path / "my-skill.zip").namelist()
    assert names == ["my-skill/SKILL.md", "my-skill/references/cli.md"]
