"""SubagentStop (gromatrix:md-runner, gromatrix:md-analyst): the last message must be one JSON object.

Required keys come from agents/*.md and references/subagents.md. A failed run (`ok: false`) needs
`ok` and a `reason` or `failed`. Blocks with {"decision": "block", "reason": ...}; honours
`stop_hook_active` so a subagent that cannot comply is not held forever.
"""

from __future__ import annotations

import json
import re

from _common import read_input

REQUIRED = {
    "gromatrix:md-runner": ("G", ("ok", "step", "mode", "systems", "failed")),
    "gromatrix:md-analyst": ("H", ("ok", "step", "analysed", "skipped", "failed", "summary_csv", "report")),
}


def problem(agent: str, msg: str) -> str | None:
    text = re.sub(r"^```(?:json)?\s*|\s*```$", "", (msg or "").strip())  # tolerate one code fence
    try:
        obj = json.loads(text)
    except ValueError:
        return "the final message is not a single JSON object (no prose, no log excerpts)"
    if not isinstance(obj, dict):
        return "the final message must be a JSON object"
    step, keys = REQUIRED[agent]
    if obj.get("ok") is False:
        return None if ("reason" in obj or "failed" in obj) else "ok is false: add `reason` (or `failed`)"
    missing = [k for k in keys if k not in obj]
    if missing:
        return f"missing keys {missing}"
    if obj["step"] != step:
        return f"step must be {step!r}"
    return None


def main() -> None:
    data = read_input()
    agent = data.get("agent_type")
    if agent not in REQUIRED or data.get("stop_hook_active"):
        return
    why = problem(agent, str(data.get("last_assistant_message") or ""))
    if why:
        reason = f"{agent} payload invalid: {why}. Reply again with exactly the JSON object from your instructions."
        print(json.dumps({"decision": "block", "reason": reason}))


if __name__ == "__main__":
    main()
