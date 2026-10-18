# The result envelope is additive

Every command prints one JSON result line. Failures now also carry a stable `code` (a key of
`gromatrix.result.CODES`) and a `next_action`; the leading keys are always in the same order
(`ok, code, error, next_action, applied`); build and run add an `applied` receipt; a few large
lists (`units`, `fixed`, `fetched`, `relaxed`) are replaced by a stub that names a file unless
`--full` is given.

`ok` and `error` keep their meaning, so skills, hooks and tests that read them still work. The
alternative, renaming to MDClaw's `success`/`message`, would repeat the breaking change of the
`gmxprep-` to `gromatrix-` rename. The cost is that codes are a second vocabulary to maintain:
a registry in `result.py` is the only place that defines them.
