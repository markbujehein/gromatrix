"""PreToolUse (Bash): gromatrix-run and gromatrix-analyze belong to the subagents.

The hooks.json `if` filters already select these commands; this script re-checks the command word so a
quoted mention (`echo gromatrix-run`) is not blocked.

UNVERIFIED: the hook docs do not state that `agent_id` is absent on the main thread. We treat a
missing or empty `agent_id` as the main thread. If Claude Code ever sets it there, this hook would
block the subagents too; disable it (delete its handlers in hooks/hooks.json). The CLI is unaffected.
"""

from __future__ import annotations

from _common import deny, read_input, segments

AGENTS = {"gromatrix-run": "gromatrix:md-runner", "gromatrix-analyze": "gromatrix:md-analyst"}


def main() -> None:
    data = read_input()
    if data.get("agent_id"):
        return
    cmd = str((data.get("tool_input") or {}).get("command", ""))
    for toks in segments(cmd):
        if toks[0] in AGENTS:
            deny(
                f"{toks[0]} runs only inside a subagent, so logs stay out of the main context. "
                f"Spawn Agent(subagent_type='{AGENTS[toks[0]]}', run_in_background=true, "
                "prompt='batch_dir: <abs path>\\naction: start'); see references/subagents.md."
            )


if __name__ == "__main__":
    main()
