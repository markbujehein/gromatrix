# Result envelope: error codes, fixed output shape, applied receipt

Status: draft for review. Date: 2026-10-17.

## Intent

Give every `gromatrix-*` command one predictable result line, so an agent (or a person) can branch on
a stable `code`, read the next step from `next_action`, and see what the command actually used
without opening artefacts. Ideas adapted from MDClaw (`matsunagalab/mdclaw`: guardrail codes,
`_envelope.py`, `_receipt.py`). Success criteria:

- Every failure carries a `code` from one registry and a `next_action`.
- The first keys of every result are in the same order.
- Build (F) and run (G) results say what was applied.
- Nothing that reads `ok` or `error` today breaks (additive; decided 2026-10-17).

Out of scope: renaming `ok`/`error`, a persistent DAG, new chemistry checks.

## Current state (read from the code)

- Each of 13 modules prints its own `json.dumps({"ok": …})` (about 44 sites). No shared emitter.
- Exceptions are caught in `main()` and printed as `{"ok": false, "error": …}`.
  Exception classes: `CheckError`, `LigandError`, `ToolError`, `ForcefieldMismatch`,
  `ExtractionError`, `TopologyError`, `HTTPError`, `Interrupted`.
- `docs/architecture.md` states the contract: one JSON line, non-zero exit on failure.

## Design

### `src/gromatrix/result.py` (new)

- `CODES: dict[str, Code]`, where `Code` holds `meaning` and `next_action` (one sentence each).
  Names are `snake_case`, prefixed by area: `ff_mismatch`, `ligand_charge_noninteger`,
  `em_not_converged`, `relax_gate_failed`, `matrix_incomplete`.
- `emit(ok, code=None, **fields)` prints one JSON line with keys in this order:
  `ok, code, error, next_action, applied, …rest`. Missing keys are omitted. On failure `next_action`
  defaults to the registry entry for `code`; a caller may override it.
- `fail(code, error, **fields)` is `emit(False, code, error=error, **fields)` and returns exit status 1.
- Unknown `code` raises at emit time (a bug, caught by the registry test below).

### (2) Codes on exceptions

Each existing exception class gets a `code` attribute (default for the class, overridable per
raise). Each `main()` catches the class and calls `fail(exc.code, str(exc))`. A bare `Exception` in
`main()` becomes `code="unexpected"`; this is the only generic code.

**Registry test** (`tests/test_result_codes.py`): scans `src/gromatrix` for `code="…"` literals and
`CODE = "…"` class attributes, then fails if a literal is not in `CODES` or an entry in `CODES` is
never raised. Written by the user with Claude reviewing, since it is a good learning exercise.

### (3) Envelope

- Key order is fixed by `emit`. Success results keep their current keys after the leading ones.
- Brief mode stubs only a named list of known-large keys per command (`STUBBABLE`: `units` in
  extract, `fixed` in fix, `systems` in gmx build; others added only when measured). A stubbed value
  is replaced by `{"stub": true, "file": "<path>", "chars": N}` and the full value is written to
  `<batch>/…/result-<command>.json`, and only when its JSON exceeds `BRIEF_LIMIT`. Everything else
  prints in full. Protected keys are never stubbed: `ok, code, error, next_action, applied,
  warnings, failed`. `--full` prints everything (escape hatch).
- `BRIEF_LIMIT` is one constant, provisionally 4000 (MDClaw's number). Commit 2 includes a throwaway
  script (not committed) that prints each command's result size on mock/1xvp data; the final value
  and its reason go in `docs/memo.md`.
- Commands that stream progress lines (`gromatrix-run`, `gromatrix-pipeline`, `gromatrix-watch`)
  keep streaming; only their final line goes through `emit`.

### (4) Applied receipt (F and G only)

`applied` is a dict built from data the command already holds (`system.json`, `run.json`):

- `options`: each option the caller passed next to the value used, and a `changed` flag when they differ.
- `facts`: forcefield key, net charge, atom count, water model, ion counts, box, for F;
  stage reached, steps, LINCS warning count, for G.

No new computation. If a fact is not available the key is absent (never guessed).

### (6) Docs

- `docs/adr/0001-forcefield-lock-is-permanent.md`, `0002-batch-directory-is-the-contract.md`,
  `0003-additive-result-envelope.md`: context, decision, consequence, a few lines each.
- `docs/memo.md`: dated entries, newest first; a correction names the entry it overturns.
- `docs/architecture.md` "Contracts every step follows" gains the envelope paragraph and links the ADRs.

## Order of work (one commit each)

1. `result.py` + exception codes + migrate error paths (2).
2. Migrate success paths + brief stubs (3).
3. `applied` for F and G (4).
4. ADRs, memo, architecture paragraph (6).

After each commit: `pixi run check` (the user confirmed on 2026-10-17 that running it is allowed).
No iterative fix loops: one considered pass, then report what the check says.

## Compatibility and risks

- `ok` and `error` are unchanged, so skills, hooks (`scripts/hooks/*`), and tests keep working.
  `tests/test_hooks.py` and `scripts/hooks/subagent_payload.py` read result lines; verify with the suite.
- Skills' `references/cli.md` need the new keys documented (mentioned in step 4's docs commit).
- Checked by grep 2026-10-17: the only exact result-dict assertion is `tests/test_pipeline.py:176`,
  which uses a fake step printing its own JSON, so it should survive. Unverified until `pixi run check`.
- Risk: a stub hides data a skill reads. Mitigation: protected keys, plus `--full`.

## Open questions

- `BRIEF_LIMIT` is provisional until measured (see Envelope). `gromatrix-work/` is not in this checkout.
- Whether `gromatrix-status` should also be stubbed (its `next` field is already the contract).
