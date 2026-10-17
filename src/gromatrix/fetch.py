"""[A] Fetch structures from the RCSB PDB.

    gromatrix-fetch --batch-dir B (--ids 1ABC 2XYZ ... | --query TEXT | --uniprot P12345 ...)
                  [--max 25] [--method xray] [--max-resolution 3.0]
                  [--format cif|pdb|both] [--out raw] [--check-updates]

* Downloads gzip-compressed mmCIF (default) and/or legacy PDB files from
  files.rcsb.org into <batch>/<out>/<ID>.cif|.pdb.
* Cache: <out>/cache.json stores sha256, size, source URL and the entry's
  revision. A file is reused only when present *and* hash-identical, so
  re-running never re-downloads, while truncated or edited files are
  replaced. --check-updates also asks the RCSB Data API for the current
  revision and refetches entries revised since they were cached.
* Search: --query runs an RCSB full-text search; --uniprot finds every entry
  containing a chain mapped to those UniProt accessions, optionally
  filtered by experimental method and
  resolution, best resolution first. Both fetch at most --max entries.
* Large structures have no legacy PDB file. With --format pdb they fail
  with a clear message; with `both` the mmCIF is still fetched.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from . import layout, result
from .net import Cache, HTTPError, get_json

FILES = "https://files.rcsb.org/download/{id}.{ext}.gz"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{id}"
UNIPROT_ATTR = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
# Classic 4-character IDs and the extended wwPDB form pdb_0000xxxx.
ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$|^pdb_0000[0-9][a-z0-9]{3}$", re.I)


def search(text: str, limit: int) -> list[str]:
    q = {
        "query": {"type": "terminal", "service": "full_text", "parameters": {"value": text}},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }
    res = get_json(SEARCH, q)
    return [hit["identifier"] for hit in (res or {}).get("result_set", [])]


def search_uniprot(accs: list[str], limit: int, method: str | None, max_res: float | None) -> list[str]:
    nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": UNIPROT_ATTR,
                "operator": "in",
                "value": accs,
            },
        }
    ]
    if method:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "exptl.method",
                    "operator": "exact_match",
                    "value": {"xray": "X-RAY DIFFRACTION", "em": "ELECTRON MICROSCOPY", "nmr": "SOLUTION NMR"}[method],
                },
            }
        )
    if max_res:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": max_res,
                },
            }
        )
    q = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
    }
    res = get_json(SEARCH, q)
    return [hit["identifier"] for hit in (res or {}).get("result_set", [])]


def revision(pid: str) -> dict:
    info = get_json(ENTRY.format(id=pid)).get("rcsb_accession_info", {})
    return {k: info.get(k) for k in ("revision_date", "major_revision", "minor_revision", "initial_release_date")}


def _revision_or_none(pid: str) -> dict | None:
    """Revision metadata is bookkeeping: an unavailable Data API must not fail a download."""
    try:
        return revision(pid)
    except HTTPError:
        return None


def fetch_one(pid: str, fmt: str, cache: Cache, check_updates: bool) -> dict:
    pid = pid.upper()
    out = {"id": pid, "files": {}}
    for ext in ("cif", "pdb") if fmt == "both" else (fmt,):
        name = f"{pid}.{ext}"
        stale, current = False, None
        if check_updates and cache.valid(name):
            current = _revision_or_none(pid)
            stale = current is not None and current != cache.entries[name].get("revision")
        need_download = stale or not cache.valid(name)
        # The revision is recorded with every download, so a later --check-updates has a baseline.
        extra = {"revision": current or _revision_or_none(pid)} if need_download else {}
        try:
            path, status = cache.fetch(name, FILES.format(id=pid, ext=ext), gunzip=True, extra=extra, force=stale)
        except HTTPError as e:
            if e.status == 404 and ext == "pdb":
                out["files"][ext] = {"status": "unavailable", "reason": "no legacy PDB file (large entry); use mmCIF"}
                continue
            raise
        out["files"][ext] = {"status": status, "sha256": cache.entries[name]["sha256"][:16], "path": str(path)}
    return out


@result.guarded
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="gromatrix-fetch", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--batch-dir", type=Path, required=True)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--ids", nargs="+")
    src.add_argument("--query", help="RCSB full-text search, e.g. 'kinase inhibitor complex'")
    src.add_argument("--uniprot", nargs="+", help="entries containing these UniProt accessions")
    ap.add_argument("--max", type=int, default=25, help="max entries for --query/--uniprot")
    ap.add_argument("--method", choices=["xray", "em", "nmr"], help="filter for --uniprot")
    ap.add_argument("--max-resolution", type=float, help="Å, filter for --uniprot")
    ap.add_argument("--format", default="cif", choices=["cif", "pdb", "both"])
    ap.add_argument("--out", default=layout.RAW)
    ap.add_argument("--check-updates", action="store_true", help="refetch entries revised since caching")
    result.add_full_flag(ap)
    a = ap.parse_args(argv)

    try:
        if a.ids:
            ids = [i.upper() for i in a.ids]
        elif a.uniprot:
            ids = search_uniprot([u.upper() for u in a.uniprot], a.max, a.method, a.max_resolution)
        else:
            ids = search(a.query, a.max)
    except HTTPError as e:
        return result.finish(result.fail(e.code, str(e)))
    bad = [i for i in ids if not ID_RE.match(i)]
    if bad:
        return result.finish(result.fail("fetch_bad_id", f"not PDB IDs: {bad}"))
    cache = Cache(a.batch_dir / a.out)
    done, failed = [], []
    for pid in dict.fromkeys(ids):
        try:
            r = fetch_one(pid, a.format, cache, a.check_updates)
            if all(f["status"] == "unavailable" for f in r["files"].values()):
                failed.append({"id": pid, "error": r["files"]})
            else:
                done.append(r)
        except HTTPError as e:
            failed.append({"id": pid, "error": "not found / obsolete" if e.status == 404 else str(e)})
    summary = {
        "ok": not failed,
        "query": a.query or a.uniprot,
        "fetched": done,
        "failed": failed,
        "counts": {
            s: sum(f["status"] == s for d in done for f in d["files"].values())
            for s in ("downloaded", "cached", "repaired", "updated")
        },
    }
    return result.finish(summary, command="fetch", out_dir=a.batch_dir / layout.CONFIG, full=a.full)


if __name__ == "__main__":
    sys.exit(main())
