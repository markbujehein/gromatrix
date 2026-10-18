# The forcefield lock is permanent

A batch directory carries `00-config/forcefield.lock.json`. The first command that touches a batch
writes it, and every later command compares its requested forcefield with it and aborts on any
difference. Ligand-parameter caches and per-system stamps record the same identity, so
`gromatrix-gmx verify` can prove that a whole comparison matrix shares one forcefield.

Mixing forcefields inside a matrix would make the systems incomparable without any visible error,
so a mismatch is a hard stop (`ff_mismatch`), not a warning. The cost is that changing forcefield
means a new batch directory; the lock is never edited in place.
