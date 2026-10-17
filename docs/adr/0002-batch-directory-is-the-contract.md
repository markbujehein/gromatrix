# The batch directory is the contract between steps

All state lives in one batch directory with numbered folders (`references/batch-layout.md`,
`gromatrix.layout`). Each step reads and writes only that directory, and `gromatrix-status`
reconstructs progress from the files alone, so any session or agent can resume a batch.

No step keeps anything in conversation memory that the next one needs. The cost is that a file
format change breaks readers, so names live in one module (`layout.py`) and the layout carries a
version (`LAYOUT_VERSION`).
