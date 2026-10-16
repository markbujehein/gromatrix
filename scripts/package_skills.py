"""Package every agent skill in .claude/skills/ as a standalone, uploadable .zip.

    python scripts/package_skills.py [--skills-dir .claude/skills] [--out dist/skills] [--check]

Each archive holds one top-level folder, `<skill-name>/SKILL.md` plus its references/.
That is the layout claude.ai (Customize > Skills > Upload) and the Skills API expect, and
it unzips straight into ~/.claude/skills/ for Claude Code.

Every skill is validated first, with the same rules as Anthropic's skill-creator:
- exactly one SKILL.md;
- YAML frontmatter with only the allowed keys;
- a kebab-case `name` equal to the folder name, at most 64 characters;
- a `description` of at most 1024 characters with no angle brackets.

Archives are byte-reproducible: sorted entries, fixed timestamps and permissions. So the
SHA256SUMS written next to them only changes when skill content changes.

--check validates without writing anything (for CI on pull requests).
Exit codes: 0 = ok, 1 = a skill failed validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import zipfile
from pathlib import Path

import yaml

ALLOWED_KEYS = {"name", "description", "license", "allowed-tools", "metadata", "compatibility"}
NAME_RE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
EXCLUDE_PARTS = {"__pycache__", "node_modules", ".DS_Store", "evals"}
FIXED_TIME = (1980, 1, 1, 0, 0, 0)  # earliest ZIP timestamp; makes archives reproducible


def validate(skill: Path) -> list[str]:
    """Return a list of problems (empty = valid)."""
    errors = []
    md = skill / "SKILL.md"
    if not md.is_file():
        return [f"{skill.name}: SKILL.md not found"]
    extra = [p for p in skill.rglob("SKILL.md") if p != md and not EXCLUDE_PARTS & set(p.parts)]
    if extra:
        errors.append(f"{skill.name}: nested SKILL.md files are rejected on upload: {[str(p) for p in extra]}")
    m = re.match(r"^---\n(.*?)\n---", md.read_text(), re.DOTALL)
    if not m:
        return errors + [f"{skill.name}: missing or malformed YAML frontmatter"]
    try:
        fm = yaml.safe_load(m.group(1))
    except yaml.YAMLError as e:  # e.g. an unquoted ': ' inside the description
        return errors + [f"{skill.name}: invalid YAML frontmatter: {e}"]
    if not isinstance(fm, dict):
        return errors + [f"{skill.name}: frontmatter must be a mapping"]
    if set(fm) - ALLOWED_KEYS:
        errors.append(f"{skill.name}: unexpected frontmatter keys {sorted(set(fm) - ALLOWED_KEYS)}")
    name, desc = str(fm.get("name", "")).strip(), fm.get("description")
    if not NAME_RE.match(name) or len(name) > 64:
        errors.append(f"{skill.name}: name {name!r} must be kebab-case and at most 64 characters")
    if name != skill.name:
        errors.append(f"{skill.name}: frontmatter name {name!r} must equal the folder name")
    if not isinstance(desc, str) or not desc.strip():
        errors.append(f"{skill.name}: description missing")
    elif len(desc.strip()) > 1024 or "<" in desc or ">" in desc:
        errors.append(f"{skill.name}: description must be at most 1024 characters with no angle brackets")
    return errors


def package(skill: Path, out_dir: Path) -> Path:
    archive = out_dir / f"{skill.name}.zip"
    files = sorted(p for p in skill.rglob("*") if p.is_file() and not EXCLUDE_PARTS & set(p.parts))
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in files:
            info = zipfile.ZipInfo(str(f.relative_to(skill.parent)), date_time=FIXED_TIME)
            info.external_attr = 0o644 << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            zf.writestr(info, f.read_bytes())
    return archive


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skills-dir", type=Path, default=Path(".claude/skills"))
    ap.add_argument("--out", type=Path, default=Path("dist/skills"))
    ap.add_argument("--check", action="store_true", help="validate only; write nothing")
    a = ap.parse_args(argv)

    skills = sorted(p for p in a.skills_dir.iterdir() if (p / "SKILL.md").is_file()) if a.skills_dir.is_dir() else []
    if not skills:
        print(json.dumps({"ok": False, "error": f"no skills found in {a.skills_dir}"}))
        return 1
    errors = [e for s in skills for e in validate(s)]
    if errors or a.check:
        print(json.dumps({"ok": not errors, "skills": [s.name for s in skills], "errors": errors}))
        return 1 if errors else 0

    a.out.mkdir(parents=True, exist_ok=True)
    built = []
    for s in skills:
        z = package(s, a.out)
        built.append(
            {
                "skill": s.name,
                "zip": z.name,
                "bytes": z.stat().st_size,
                "sha256": hashlib.sha256(z.read_bytes()).hexdigest(),
            }
        )
    (a.out / "SHA256SUMS").write_text("".join(f"{b['sha256']}  {b['zip']}\n" for b in built))
    print(json.dumps({"ok": True, "out": str(a.out), "packages": built}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
