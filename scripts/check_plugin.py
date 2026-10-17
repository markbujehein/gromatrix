"""Validate the gromatrix Claude Code plugin (the repository root is the plugin).

    python scripts/check_plugin.py [--root .]

Checks, with the same rules `claude plugin validate` and the skill uploader apply:
- .claude-plugin/plugin.json: valid JSON, kebab-case `name`, a `version` equal to pyproject.toml's;
- .claude-plugin/marketplace.json: name, owner and a plugin entry whose source resolves to this root;
- every skills/<name>/SKILL.md: YAML frontmatter with allowed keys only, `name` equal to the folder
  name, a `description` of at most 1024 characters with no angle brackets, no nested SKILL.md;
- every agents/<name>.md: frontmatter `name` equal to the file name, a description, a tools list;
- hooks/hooks.json: valid JSON, command handlers only, every ${CLAUDE_PLUGIN_ROOT} script exists;
- every ${CLAUDE_PLUGIN_ROOT}/<path> mentioned in a skill or agent exists.

Exit codes: 0 = ok, 1 = a problem was found. Prints one JSON line.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import tomllib
import yaml

SKILL_KEYS = {"name", "description", "license", "allowed-tools", "metadata", "compatibility"}
AGENT_KEYS = {"name", "description", "model", "effort", "maxTurns", "tools", "disallowedTools", "skills",
              "memory", "background", "isolation", "color"}  # fmt: skip
NAME_RE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
ROOT_REF = re.compile(r"\$\{CLAUDE_PLUGIN_ROOT\}/([\w./-]+)")


def frontmatter(path: Path) -> tuple[dict | None, str | None]:
    m = re.match(r"^---\n(.*?)\n---", path.read_text(), re.DOTALL)
    if not m:
        return None, f"{path}: missing or malformed YAML frontmatter"
    try:
        fm = yaml.safe_load(m.group(1))
    except yaml.YAMLError as e:  # e.g. an unquoted ': ' inside the description
        return None, f"{path}: invalid YAML frontmatter: {e}"
    return (fm, None) if isinstance(fm, dict) else (None, f"{path}: frontmatter must be a mapping")


def check_skill(skill: Path) -> list[str]:
    md = skill / "SKILL.md"
    if not md.is_file():
        return [f"{skill.name}: SKILL.md not found"]
    errors = []
    extra = [p for p in skill.rglob("SKILL.md") if p != md]
    if extra:
        errors.append(f"{skill.name}: nested SKILL.md files are rejected: {[str(p) for p in extra]}")
    fm, err = frontmatter(md)
    if err:
        return errors + [err]
    if set(fm) - SKILL_KEYS:
        errors.append(f"{skill.name}: unexpected frontmatter keys {sorted(set(fm) - SKILL_KEYS)}")
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


def check_agent(agent: Path) -> list[str]:
    fm, err = frontmatter(agent)
    if err:
        return [err]
    errors = []
    if set(fm) - AGENT_KEYS:
        errors.append(f"{agent.name}: unsupported frontmatter keys {sorted(set(fm) - AGENT_KEYS)}")
    if fm.get("name") != agent.stem:
        errors.append(f"{agent.name}: frontmatter name {fm.get('name')!r} must equal the file name")
    if not str(fm.get("description", "")).strip():
        errors.append(f"{agent.name}: description missing")
    if not fm.get("tools"):
        errors.append(f"{agent.name}: list the tools the agent may use")
    return errors


def check_manifests(root: Path) -> list[str]:
    errors = []
    try:
        plugin = json.loads((root / ".claude-plugin" / "plugin.json").read_text())
        market = json.loads((root / ".claude-plugin" / "marketplace.json").read_text())
        version = tomllib.loads((root / "pyproject.toml").read_text())["project"]["version"]
    except (OSError, ValueError, KeyError) as e:
        return [f"plugin manifests unreadable: {e}"]
    if not NAME_RE.match(str(plugin.get("name", ""))):
        errors.append("plugin.json: name must be kebab-case")
    if plugin.get("version") != version:
        errors.append(f"plugin.json version {plugin.get('version')!r} != pyproject.toml {version!r}")
    if not plugin.get("license"):
        errors.append("plugin.json: license missing")
    if not (market.get("name") and market.get("owner", {}).get("name")):
        errors.append("marketplace.json: name and owner.name are required")
    entries = {e.get("name"): e for e in market.get("plugins", [])}
    if plugin.get("name") not in entries:
        errors.append("marketplace.json: no entry for the plugin")
    elif entries[plugin["name"]].get("source") not in ("./", "."):
        errors.append("marketplace.json: the plugin entry must use source './' (this repository)")
    return errors


def check_hooks(root: Path) -> list[str]:
    """hooks/hooks.json: valid JSON, known structure, and every ${CLAUDE_PLUGIN_ROOT} script exists."""
    path = root / "hooks" / "hooks.json"
    if not path.is_file():
        return []
    try:
        hooks = json.loads(path.read_text())["hooks"]
    except (ValueError, KeyError) as e:
        return [f"hooks/hooks.json: invalid JSON or no top-level 'hooks' key: {e}"]
    errors = []
    for event, groups in hooks.items():
        for group in groups:
            for h in group.get("hooks", []):
                if h.get("type") != "command" or not h.get("command"):
                    errors.append(f"hooks/hooks.json: {event}: handlers need type 'command' and a command")
                    continue
                for ref in ROOT_REF.findall(h["command"]):
                    if not (root / ref).is_file():
                        errors.append(f"hooks/hooks.json: {event}: ${{CLAUDE_PLUGIN_ROOT}}/{ref} does not exist")
    return errors


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path("."))
    a = ap.parse_args(argv)
    root = a.root.resolve()
    skills = sorted(p for p in (root / "skills").glob("*") if p.is_dir())
    agents = sorted((root / "agents").glob("*.md"))
    errors = check_manifests(root) + check_hooks(root)
    errors += [e for s in skills for e in check_skill(s)]
    errors += [e for g in agents for e in check_agent(g)]
    for md in [*root.glob("skills/**/*.md"), *agents]:
        for ref in ROOT_REF.findall(md.read_text()):
            if not (root / ref.rstrip(".")).exists():
                errors.append(f"{md.relative_to(root)}: ${{CLAUDE_PLUGIN_ROOT}}/{ref} does not exist")
    if not skills:
        errors.append(f"no skills found in {root / 'skills'}")
    result = {
        "ok": not errors,
        "skills": [s.name for s in skills],
        "agents": [g.stem for g in agents],
        "errors": errors,
    }
    print(json.dumps(result))
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
