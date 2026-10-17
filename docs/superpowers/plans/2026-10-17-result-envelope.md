# Result envelope Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every `gromatrix-*` command prints one result line with a stable `code`, a `next_action`, a fixed key order, optional stubs for large values, and (for F and G) an `applied` receipt.

**Architecture:** One new module `src/gromatrix/result.py` owns the code registry, key ordering, stubbing and printing. Each `main()` keeps its own `try/except` and hands its dict to `result.finish`. Exceptions carry a default `code`. Purely additive: `ok` and `error` are unchanged.

**Tech Stack:** Python 3.10+, argparse, pytest (+ hypothesis available), pixi (`pixi run check`).

**Spec:** `docs/superpowers/specs/2026-10-17-result-envelope-design.md`

## Global Constraints

- Additive only: never rename or remove `ok` or `error`; skills, hooks (`scripts/hooks/*`) and tests keep working.
- Key order of every result: `ok, code, error, next_action, applied`, then the command's own keys.
- Code names are `snake_case` with an area prefix (`ff_`, `ligand_`, `em_`, `relax_`, `matrix_`, `fetch_`, `run_`); the one generic code is `unexpected`.
- `BRIEF_LIMIT = 4000` is provisional; it is one constant. Stubbing applies only to keys named in `STUBBABLE`.
- Protected keys are never stubbed: `ok, code, error, next_action, applied, warnings, failed`.
- `--full` prints everything unstubbed.
- Streaming commands (`gromatrix-run`, `-pipeline`, `-watch`) keep streaming; only the final line goes through `finish`.
- Verification is `pixi run check` after each task (allowed per `/home/mark/Projects/CLAUDE.md`). One considered pass per task, no fix loops: if it fails, report the output and stop.
- Line length 120 (ruff); one commit per task; commit messages end with the attribution line from the session reminder.

## Review Focus

- A command that raises an exception whose class has no `code` prints `code: "unexpected"`, never a traceback.
- A caller passing an unregistered `code` fails loudly in tests, not silently in production output.
- A stubbed value must be recoverable: the `file` named in the stub exists and holds the full JSON.
- A result with `ok: false` and no `error` text still prints valid JSON with `code`.
- `gromatrix-status`'s `next` field is untouched by the envelope (it is already the contract).

---

### Task 1: `result.py` core (registry, ordering, finish)

**Files:**
- Create: `src/gromatrix/result.py`
- Create: `tests/test_result.py`

**Interfaces:**
- Produces:
  - `@dataclass(frozen=True) class Code: meaning: str; next_action: str`
  - `CODES: dict[str, Code]` (starts with `unexpected`; later tasks add entries)
  - `ENVELOPE_ORDER: tuple[str, ...] = ("ok", "code", "error", "next_action", "applied")`
  - `order(rec: dict) -> dict` (leading keys first, missing ones omitted, rest in original order)
  - `fail(code: str, error: str, **fields) -> dict` (builds `{"ok": False, "code", "error", "next_action", **fields}`; `next_action` from the registry unless passed; raises `KeyError` for an unknown code)
  - `from_exception(exc: BaseException) -> dict` (uses `exc.code` if present and registered, else `unexpected`; `error` is `f"{type(exc).__name__}: {exc}"`)
  - `finish(rec: dict, *, flush: bool = False) -> int` (prints `json.dumps(order(rec))`, returns 0 if `rec.get("ok")` else 1)

- [ ] **Step 1: Write the failing tests** in `tests/test_result.py`: `test_order_puts_envelope_keys_first`, `test_fail_fills_next_action_from_registry`, `test_fail_unknown_code_raises`, `test_from_exception_without_code_is_unexpected`, `test_finish_prints_one_line_and_returns_exit_status` (use `capsys`; assert exactly one line, `json.loads` it, return 1 for `ok: False`), `test_ok_false_without_error_still_valid_json`.
- [ ] **Step 2: Run** `pixi run -e md pytest tests/test_result.py -q`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement** the signatures above in `result.py`. `order` must not mutate its input.
- [ ] **Step 4: Run** the same command. Expected: all pass.
- [ ] **Step 5: Run** `pixi run check`; then commit `result: shared emit/finish helper and code registry`.

### Task 2: Codes on exceptions and every `main()` (spec item 2)

**Files:**
- Modify: `src/gromatrix/{checks,ligands,gmx,forcefields,extract,topology,net,run}.py` (exception classes)
- Modify: `src/gromatrix/{analyze,afdb,compute,extract,fetch,fix,gmx_prep,pipeline,relax,report,run,status,watch}.py` (print sites)
- Modify: `src/gromatrix/result.py` (register codes)
- Create: `tests/test_result_codes.py` (**the user writes this one**; Claude reviews)

**Interfaces:**
- Consumes: `result.fail`, `result.from_exception`, `result.finish`, `CODES`.
- Produces: class attribute `code: str` on `CheckError`, `LigandError`, `ToolError`, `ForcefieldMismatch`, `ExtractionError`, `TopologyError`, `HTTPError`, `Interrupted`, with these initial registry entries (code, next_action in one sentence each):
  `check_failed`, `ligand_invalid`, `tool_failed`, `ff_mismatch`, `extract_failed`, `topology_invalid`, `fetch_http`, `run_interrupted`, plus `em_not_converged`, `relax_gate_failed`, `matrix_incomplete`, `batch_missing`, `unexpected`.
  Read the matching skill reference (`skills/*/references/`) for the exact recovery sentence; do not invent recovery steps the skills do not state.

- [ ] **Step 1 (user): write `tests/test_result_codes.py`** with `test_every_raised_code_is_registered` (scan `src/gromatrix/*.py` for `code="..."` literals and `code = "..."` class attributes) and `test_every_registered_code_is_raised` (every `CODES` key except `unexpected` appears in the scan). Run it; expected FAIL until step 3.
- [ ] **Step 2: Add `code` class attributes** to the eight exception classes.
- [ ] **Step 3: Replace the print sites.** The common pattern, `except Exception as e: out = {"ok": False, "error": f"{type(e).__name__}: {e}"}` then `print(json.dumps(out)); return 0 if out.get("ok") else 1`, becomes `out = result.from_exception(e)` and `return result.finish(out)`. Sites that print `{"ok": False, "error": …}` directly (`fetch.py:163,167`, `extract.py:685`, `fix.py:349`, `afdb.py:316`, `compute.py:116`, `watch.py:66`, `status.py` batch-missing) use `result.fail(<code>, <message>)`. `run.py`'s `emit()` for progress lines is unchanged.
- [ ] **Step 4: Run** `pixi run check`. Expected: lint, plugin validation and tests pass, including the registry test.
- [ ] **Step 5: Commit** `result: stable error codes on every exception and failure path`.

### Task 3: Envelope on success paths and brief stubs (spec item 3)

**Files:**
- Modify: `src/gromatrix/result.py`
- Modify: the same command modules (success prints → `finish`)
- Modify: `tests/test_result.py`

**Interfaces:**
- Consumes: `finish`, `order`.
- Produces:
  - `BRIEF_LIMIT: int = 4000`
  - `STUBBABLE: dict[str, tuple[str, ...]] = {"extract": ("units",), "fix": ("fixed",), "gmx": ("systems",)}`
  - `PROTECTED: frozenset[str]` (the five envelope keys + `warnings`, `failed`)
  - `stub_large(rec: dict, command: str, out_dir: Path, limit: int = BRIEF_LIMIT) -> dict` (replaces each listed key whose `json.dumps` length exceeds `limit` with `{"stub": True, "file": str(path), "chars": n}` and writes the full value to `out_dir / f"result-{command}.json"` via the project's `atomic_write`)
  - `add_full_flag(ap: argparse.ArgumentParser) -> None` (adds `--full`, help "print large values in full")
  - `finish(rec, *, command: str | None = None, out_dir: Path | None = None, full: bool = False, flush: bool = False) -> int` (stubs only when `command` and `out_dir` are given and `full` is false)

- [ ] **Step 1: Write failing tests:** `test_stub_replaces_only_listed_oversize_keys`, `test_stub_file_holds_full_value` (Review Focus: recoverable), `test_protected_keys_never_stubbed`, `test_full_flag_disables_stubbing`, `test_small_value_not_stubbed`.
- [ ] **Step 2: Run** `pixi run -e md pytest tests/test_result.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the signatures above. Reuse `atomic_write` from `gromatrix.net` (same helper `afdb.py` uses).
- [ ] **Step 4: Migrate success prints** in the 13 modules to `finish(out, command=..., out_dir=batch_dir / layout.CONFIG, full=a.full)`; add `result.add_full_flag(ap)` to the three stubbable commands only. Streaming commands: only the last line.
- [ ] **Step 5: Measure.** Write a throwaway script in the scratchpad (not committed) that runs each command on mock or 1xvp data and prints result sizes. If no such data exists in this checkout (`gromatrix-work/` is absent), state "unmeasured" and keep 4000.
- [ ] **Step 6: Run** `pixi run check`; commit `result: fixed key order and opt-out stubs for large values`.

### Task 4: `applied` receipt for F and G (spec item 4)

**Files:**
- Modify: `src/gromatrix/result.py`, `src/gromatrix/gmx_prep.py` (`cmd_build`), `src/gromatrix/run.py` (`cmd_run`)
- Modify: `tests/test_result.py`, `tests/test_run.py`

**Interfaces:**
- Produces:
  - `receipt(options: dict, used: dict, facts: dict) -> dict` returning `{"options": {name: {"passed": v, "used": u, "changed": bool}}, "facts": {k: v}}`; `changed` is true when `passed != used`; names with `used is None` are omitted; fact keys with value `None` are dropped.
  - `cmd_build` result gains `applied` with options `box, pad, conc, temp, md_ns, maxwarn` and facts read from each built `system.json` (forcefield key, net charge, atom count, water model, ion counts, box). Use only keys `system.json` already stores; read the file to see which.
  - `cmd_run` result gains `applied` with stage reached, steps, LINCS warning count from `run.json`.

- [ ] **Step 1: Write failing tests:** `test_receipt_marks_changed_options`, `test_receipt_drops_missing_facts`, and in `tests/test_run.py` `test_run_result_has_applied_stage` (follow the existing fake-gmx tests in that file).
- [ ] **Step 2: Run** the tests; expected FAIL.
- [ ] **Step 3: Implement** `receipt`, then wire it into `cmd_build` and `cmd_run`. Missing facts stay absent; never guess.
- [ ] **Step 4: Run** `pixi run check`; commit `result: applied receipt for build and run`.

### Task 5: Docs (spec item 6)

**Files:**
- Create: `docs/adr/0001-forcefield-lock-is-permanent.md`, `docs/adr/0002-batch-directory-is-the-contract.md`, `docs/adr/0003-additive-result-envelope.md`, `docs/memo.md`
- Modify: `docs/architecture.md` ("Contracts every step follows"), `skills/*/references/cli.md` (new keys), `CHANGELOG.md`

**Interfaces:** none (documentation).

- [ ] **Step 1: Write the three ADRs**, each: title, context (2–3 lines), decision, consequence. Source the facts from `docs/architecture.md`, the forcefield lock code and spec; no claims that are not in them.
- [ ] **Step 2: Write `docs/memo.md`**: newest first; entries for (a) 2026-10-17 envelope design with the measured or "unmeasured" `BRIEF_LIMIT` result from Task 3, (b) MDClaw review (what was adopted and what was not).
- [ ] **Step 3: Update** the architecture paragraph (envelope order, codes, `--full`, `applied`), link the ADRs, document the new keys in each skill's `references/cli.md`, and add a CHANGELOG entry ("Added: result `code`, `next_action`, `applied`; `ok`/`error` unchanged").
- [ ] **Step 4: Run** `pixi run check` (the plugin validator checks skills and references); commit `docs: ADRs, memo and envelope contract`.

---

## Self-review

- **Spec coverage:** codes + registry test (Task 2), envelope order, brief stubs, `--full`, `BRIEF_LIMIT` measurement (Task 3), `applied` for F and G (Task 4), ADRs, memo, architecture, CHANGELOG (Task 5), compatibility (Global Constraints). Open question on `gromatrix-status` stubbing: left out on purpose (status untouched).
- **Type consistency:** `finish`'s signature grows in Task 3 with keyword-only additions; Task 2 callers pass only `rec` and stay valid.
- **Unverified until run:** `tests/test_pipeline.py:176` surviving the change; the exact keys `system.json` stores (Task 4 reads it first); that `gromatrix-work/` data exists to measure sizes.
