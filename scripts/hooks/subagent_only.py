"""PreToolUse (Bash): gromatrix-run and gromatrix-analyze belong to the subagents.

The hooks.json `if` filters already select these commands; this script re-checks the command word so a
quoted mention (`echo gromatrix-run`) is not blocked.

Verified in a live headless session (2026-10-17): `agent_id` (and `agent_type`) is absent from the hook
input on the main thread and present inside a subagent. A missing or empty `agent_id` is therefore the main
thread. The CLI is unaffected if the hook is removed.

Not blocked from the main thread, because they cannot flood the context: `--background` (detached, the
documented fallback where the Agent tool is unavailable), `--slurm` (only writes the sbatch script) and
`--help`.
"""

from __future__ import annotations

from _common import deny, read_input, segments

AGENTS = {"gromatrix-run": "gromatrix:md-runner", "gromatrix-analyze": "gromatrix:md-analyst"}
QUIET_FLAGS = {"--background", "--slurm", "--help", "-h"}


def main() -> None:
    data = read_input()
    if data.get("agent_id"):
        return
    cmd = str((data.get("tool_input") or {}).get("command", ""))
    for toks in segments(cmd):
        if toks[0] in AGENTS and not QUIET_FLAGS & set(toks):
            deny(
                f"{toks[0]} runs only inside a subagent, so logs stay out of the main context. "
                f"Spawn Agent(subagent_type='{AGENTS[toks[0]]}', run_in_background=true, "
                "prompt='batch_dir: <abs path>\\naction: start'); see references/subagents.md."
            )


if __name__ == "__main__":
    main()
