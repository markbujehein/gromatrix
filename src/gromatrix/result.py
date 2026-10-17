"""The one result line every gromatrix command prints, and the registry of failure codes.

Additive on top of the historic `{"ok": ..., "error": ...}` contract: `ok` and `error` keep their
meaning. A failure also carries a stable `code` (a key of `CODES`) and a `next_action`, and the
leading keys always come in the same order, so an agent can branch on `code` instead of parsing
prose. Design: docs/superpowers/specs/2026-10-17-result-envelope-design.md.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Code:
    meaning: str
    next_action: str


class CodedError(RuntimeError):
    """Base of the exceptions gromatrix raises on purpose: `code` names the failure in `CODES`.

    A subclass sets its default `code`; one raise can override it: `CodedError(msg, code="batch_missing")`.
    """

    code = "unexpected"

    def __init__(self, *args, code: str | None = None):
        super().__init__(*args)
        if code is not None:
            self.code = code


CODES: dict[str, Code] = {
    "unexpected": Code(
        "an error no command anticipated",
        "read the error text; if it is not a usage mistake, report it with the command and its arguments",
    ),
    "check_failed": Code(
        "a physical or bookkeeping gate failed (charge, LINCS, unfinished mdrun)",
        "read the error: it names the input to fix; fix it and rerun the command, never raise a tolerance",
    ),
    "em_not_converged": Code(
        "energy minimisation did not converge",
        "look for clashes or a bad ligand pose, repair the structure or ligand, rebuild the system and rerun",
    ),
    "ligand_invalid": Code(
        "a ligand failed parameterisation or a chemistry check",
        "fix the ligand SMILES, SDF or charge in the manifest and rerun the build",
    ),
    "tool_failed": Code(
        "an external tool (gmx, acpype, cgenff, mdrun) failed or is not installed",
        "open the log named in the error; install the tool if it is missing",
    ),
    "ff_mismatch": Code(
        "a command would mix forcefields in one batch",
        "use the batch's locked forcefield, or start a new batch directory for a different one",
    ),
    "extract_failed": Code(
        "no biological unit could be chosen or built",
        "see unit.json; adjust --unit or --max-chains, and use --allow-incomplete only with the user's consent",
    ),
    "topology_invalid": Code(
        "a GROMACS topology could not be read or is inconsistent",
        "regenerate the topology; report the error if it persists",
    ),
    "fetch_http": Code(
        "a download failed (network error or HTTP status)",
        "check the ID and the network; 404 usually means an obsolete or wrong ID",
    ),
    "fetch_bad_id": Code("a PDB ID is malformed", "pass 4-character PDB IDs (digit plus three alphanumerics)"),
    "run_interrupted": Code(
        "mdrun stopped early (signal or --maxh); the state is resumable",
        "rerun gromatrix-run with the same arguments: it resumes from the checkpoint",
    ),
    "matrix_incomplete": Code(
        "an enzyme x target pair has no docked complex",
        "supply the docked complexes, or run native pairs only with --allow-incomplete (needs the user's consent)",
    ),
    "batch_missing": Code(
        "the batch directory or the systems it should hold do not exist",
        "check --batch-dir; run gromatrix-status to see the next step",
    ),
    "input_missing": Code(
        "the input folder for this step is empty or absent",
        "run the previous step first; gromatrix-status shows which",
    ),
    "bad_arguments": Code(
        "a step of the pipeline rejected its command-line arguments",
        "run the named command with --help and correct the step's arguments in the manifest",
    ),
    "compute_unset": Code(
        "no compute configuration has been saved",
        "run gromatrix-compute set --mode ... (ask the user for the mode)",
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


# Brief mode. Only these keys of these commands are ever replaced by a stub, and only when their JSON is
# longer than BRIEF_LIMIT characters (4000 is MDClaw's figure; unmeasured on gromatrix output).
BRIEF_LIMIT = 4000
STUBBABLE: dict[str, tuple[str, ...]] = {
    "extract": ("units",),
    "fix": ("fixed",),
    "fetch": ("fetched",),
    "relax": ("relaxed",),
}
PROTECTED = frozenset({*ENVELOPE_ORDER, "warnings", "failed"})


def receipt(options: dict, used: dict, facts: dict) -> dict:
    """What a command did with its options (`changed` = it used a different value than it was given) and
    the stage's key figures. Options whose used value is unknown, and facts that are None, are left out."""
    return {
        "options": {
            k: {"passed": v, "used": used[k], "changed": v is not None and v != used[k]}
            for k, v in options.items()
            if used.get(k) is not None
        },
        "facts": {k: v for k, v in facts.items() if v is not None},
    }


def stub_large(rec: dict, command: str, out_dir: Path, limit: int = BRIEF_LIMIT) -> dict:
    """Replace each oversize STUBBABLE value by {"stub", "file", "chars"}; the full values go to a file."""
    from .net import atomic_write  # net imports this module, so import it late

    big = {
        k: len(json.dumps(rec[k]))
        for k in STUBBABLE.get(command, ())
        if k in rec and k not in PROTECTED and len(json.dumps(rec[k])) > limit
    }
    if not big:
        return rec
    path = Path(out_dir) / f"result-{command}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_write(path, json.dumps({k: rec[k] for k in big}, indent=1))
    return {k: {"stub": True, "file": str(path), "chars": big[k]} if k in big else v for k, v in rec.items()}


def add_full_flag(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--full", action="store_true", help="print large values in full instead of a stub + file")


def finish(
    rec: dict,
    *,
    command: str | None = None,
    out_dir: Path | None = None,
    full: bool = False,
    flush: bool = False,
) -> int:
    """Print the result as one JSON line and return the process exit status.

    With `command` and `out_dir`, oversize values of that command's STUBBABLE keys become stubs unless `full`.
    """
    if command and out_dir and not full:
        rec = stub_large(rec, command, out_dir)
    print(json.dumps(order(rec)), flush=flush)
    return 0 if rec.get("ok") else 1
