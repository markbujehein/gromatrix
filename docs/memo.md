# Memo

Dated record of decisions and findings, newest first. A correction names the entry it overturns.

## 2026-10-17: result envelope (codes, fixed order, receipt)

- Adopted three ideas from MDClaw (`matsunagalab/mdclaw`, read, not run): stable error codes,
  a fixed result-key order with stubs for large values, an `applied` receipt. Not adopted: its
  persistent DAG engine, the renaming of `ok`/`error`, the size of its test suite. See ADR 0003.
- `BRIEF_LIMIT = 4000` is MDClaw's figure. It was **not measured** on gromatrix output: the
  checkout had no real batch directory. Stubs apply only to the named lists (`STUBBABLE`).
  Measure on a real 1XVP batch before changing it.
- `cmd_build`'s `applied` receipt is covered by `build_facts` unit tests only; the command needs
  GROMACS and has no end-to-end test. Unverified until a real build runs.
- Python 3.10 support dropped (`requires-python = ">=3.11"`): `tomllib` is standard from 3.11.
  Ruff target raised to py311; the plugin hook scripts keep `timezone.utc` because they run on the
  user's system `python3`.
