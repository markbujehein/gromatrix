"""The failure codes raised in the source and the registry in gromatrix.result must be the same set."""

import re
from pathlib import Path

from gromatrix.result import CODES

SRC = Path(__file__).resolve().parent.parent / "src" / "gromatrix"
PATTERNS = (
    re.compile(r"^\s+code = \"([a-z_]+)\"", re.M),  # class attribute: code = "x"
    re.compile(r"\bcode=\"([a-z_]+)\""),  # per-raise override: Error(msg, code="x")
    re.compile(r"\bfail\(\s*\"([a-z_]+)\""),  # result.fail("x", ...)
)


def raised_codes() -> dict[str, set[str]]:
    """code -> files that raise it (result.py itself defines the registry, so it is not scanned)."""
    found: dict[str, set[str]] = {}
    for path in SRC.glob("*.py"):
        if path.name == "result.py":
            continue
        text = path.read_text()
        for pattern in PATTERNS:
            for code in pattern.findall(text):
                found.setdefault(code, set()).add(path.name)
    return found


def test_every_raised_code_is_registered():
    unknown = {c: sorted(f) for c, f in raised_codes().items() if c not in CODES}
    assert not unknown, f"codes used in the source but missing from result.CODES: {unknown}"


def test_every_registered_code_is_raised():
    unused = sorted(set(CODES) - set(raised_codes()) - {"unexpected"})
    assert not unused, f"codes in result.CODES that no source file raises: {unused}"
