# Subagents for the long and noisy steps (G and H)

The main context never runs `gmxprep-run` or `gmxprep-analyze`. Two plugin subagents do, and return
one small JSON payload. This keeps hours of log output and polling out of the main context window,
and lets the user keep working while the MD runs.

## Spawning

```
Agent(
  description = "Run MD for <batch name>",
  subagent_type = "gromatrix:md-runner",       # or "gromatrix:md-analyst"
  run_in_background = true,
  prompt = "batch_dir: /abs/path/B\naction: start\nReturn only the JSON payload."
)
```

`action` is `start` (new run), `resume` (continue from checkpoints) or `check` (status only, no new run).
The orchestrator gets a notification when the agent finishes; do not poll, do not read the agent's
transcript. `gmxprep-status` names the right agent in `next.agent`.

## Payloads

`gromatrix:md-runner`:
```json
{"ok": true, "step": "G", "mode": "local", "systems": {"E1": "done", "E1__T1": "failed"},
 "failed": [{"id": "E1__T1", "phase": "nvt", "reason": "LINCS warnings in NVT", "evidence": "06-systems/E1__T1/nvt.log"}],
 "lincs_warnings": {"E1": 0}, "wall_minutes": 214, "slurm": null, "next": "gmxprep-status ..."}
```
For `slurm` and `remote` modes, `systems` stays empty, `slurm` holds `{script, submitted, job_id}`.

`gromatrix:md-analyst`:
```json
{"ok": true, "step": "H", "analysed": ["E1", "E1__T1"], "skipped": [], "failed": [],
 "flags": {"E1__T1": ["no_target_hbonds"]}, "summary_csv": "08-analysis/analysis_summary.csv",
 "report": "09-reports/README.md", "time_ns": {"E1": 20.0}}
```

## Rules for the agents
- Return only the JSON payload as the final message; never paste log excerpts.
- On failure set `ok: false`, name the system, phase and a one-line reason, and give the evidence path.
- Never edit `run.json`, delete checkpoints, add `-maxwarn` or loosen a gate.
- Poll with long waits (a background wait loop on the run's pid), not many short checks.
