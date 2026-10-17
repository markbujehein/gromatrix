"""The one result line every gromatrix command prints, and the registry of failure codes.

Additive on top of the historic `{"ok": ..., "error": ...}` contract: `ok` and `error` keep their
meaning. A failure also carries a stable `code` (a key of `CODES`) and a `next_action`, and the
leading keys always come in the same order, so an agent can branch on `code` instead of parsing
prose. Design: docs/superpowers/specs/2026-10-17-result-envelope-design.md.
"""

from __future__ import annotations

import json
from dataclasses import dataclass


@dataclass(frozen=True)
class Code:
    meaning: str
    next_action: str


CODES: dict[str, Code] = {
    "unexpected": Code(
        "an error no command anticipated",
        "read the error text; if it is not a usage mistake, report it with the command and its arguments",
    ),
}

ENVELOPE_ORDER = ("ok", "code", "error", "next_action", "applied")


def order(rec: dict) -> dict:
    """Envelope keys first, in ENVELOPE_ORDER; every other key keeps its order. Input is not changed."""
    head = {k: rec[k] for k in ENVELOPE_ORDER if k in rec}
    return {**head, **{k: v for k, v in rec.items() if k not in head}}


def fail(code: str, error: str, **fields) -> dict:
    """A failed result. `next_action` comes from the registry unless the caller passes its own."""
    entry = CODES[code]  # an unregistered code is a bug: KeyError
    return {"ok": False, "code": code, "error": error, "next_action": entry.next_action, **fields}


def from_exception(exc: BaseException) -> dict:
    """The failed result for an exception: its class `code` if registered, else `unexpected`."""
    code = getattr(exc, "code", None)
    return fail(code if code in CODES else "unexpected", f"{type(exc).__name__}: {exc}")


def finish(rec: dict, *, flush: bool = False) -> int:
    """Print the result as one JSON line and return the process exit status."""
    print(json.dumps(order(rec)), flush=flush)
    return 0 if rec.get("ok") else 1
