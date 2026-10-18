# The result line

Every `gromatrix-*` command prints one JSON line, and exits non-zero when `ok` is false. The
leading keys are always in this order; absent ones are omitted:

| Key | Meaning |
|---|---|
| `ok` | true when nothing failed |
| `code` | on failure, a stable name from `gromatrix.result.CODES` (for example `ff_mismatch`, `em_not_converged`, `batch_missing`); branch on this, not on `error` |
| `error` | the human-readable text |
| `next_action` | on failure, the one-sentence recovery from the registry |
| `applied` | build (F) and run (G) only: `facts`, plus `options` (what you passed, what was used, `changed`) where the command can tell them apart: build compares the forcefield, run reports none (it uses every option as given) |

When some items of a batch fail (`fetch`, `fix`, `extract`, `relax`, `afdb`, build), the top level has
`ok: false` and no `code`: each entry of `failed` carries its own `code` instead.

Then the command's own keys. `extract`, `fix`, `fetch` and `relax` replace a very long list
(`units`, `fixed`, `fetched`, `relaxed`) by `{"stub": true, "file": …, "chars": N}` and write the
full list to that file in `00-config/`; add `--full` to print it inline. `gromatrix-status` keeps
its own shape, its `next` field is the contract.
