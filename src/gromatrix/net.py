"""Minimal HTTP client (stdlib only) with retries, plus a content-addressed
download cache used by steps A and E."""

from __future__ import annotations

import gzip
import hashlib
import http.client
import json
import os
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

USER_AGENT = "gromatrix/1.0 (+https://github.com/markbujehein/gromatrix)"
RETRY_STATUS = {429, 500, 502, 503, 504}
MAX_WAIT = 120.0  # s; never let a server's Retry-After stall a batch for long


def _retry_after(value: str | None, default: float) -> float:
    """Retry-After is either delta-seconds or an HTTP-date (RFC 9110)."""
    if not value:
        return default
    try:
        wait = float(value)
    except ValueError:
        from email.utils import parsedate_to_datetime

        try:
            wait = (parsedate_to_datetime(value) - datetime.now(timezone.utc)).total_seconds()
        except (TypeError, ValueError):
            wait = default
    return min(max(wait, 0.0), MAX_WAIT)


def atomic_write(path: Path, text: str) -> None:
    """Write via a temp file + rename, so an interrupted write never corrupts `path`."""
    tmp = Path(path).with_suffix(Path(path).suffix + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


class HTTPError(RuntimeError):
    def __init__(self, url: str, status: int | None, msg: str):
        super().__init__(f"{status or 'network'} for {url}: {msg}")
        self.url, self.status = url, status


def get(
    url: str, data: bytes | None = None, retries: int = 4, timeout: float = 60.0, headers: dict | None = None
) -> bytes:
    """GET (or POST when `data` is given) with exponential backoff on 429/5xx
    and network errors. 404 and other 4xx fail immediately."""
    hdrs = {"User-Agent": USER_AGENT, **(headers or {})}
    if data is not None:
        hdrs.setdefault("Content-Type", "application/json")
    delay = 2.0
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            if e.code not in RETRY_STATUS or attempt == retries:
                raise HTTPError(url, e.code, e.reason) from e
            wait = _retry_after(e.headers.get("Retry-After"), delay)
        # http.client.HTTPException covers IncompleteRead (body truncated mid-transfer).
        except (urllib.error.URLError, http.client.HTTPException, TimeoutError, ConnectionError) as e:
            if attempt == retries:
                raise HTTPError(url, None, str(e)) from e
            wait = delay
        time.sleep(wait)
        delay *= 2
    raise AssertionError("unreachable")


def get_json(url: str, payload: dict | None = None, **kw):
    body = get(url, json.dumps(payload).encode() if payload is not None else None, **kw)
    return json.loads(body) if body else None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Cache:
    """`<dir>/cache.json` records url, sha256, size and fetch time for every file.
    A file is reused only if it exists *and* its hash still matches, so a
    truncated or edited file is fetched again instead of being trusted."""

    def __init__(self, directory: Path):
        self.dir = Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.path = self.dir / "cache.json"
        self.entries: dict = json.loads(self.path.read_text()) if self.path.exists() else {}

    def valid(self, name: str) -> bool:
        e, f = self.entries.get(name), self.dir / name
        return bool(e) and f.exists() and sha256(f) == e["sha256"]

    def fetch(
        self, name: str, url: str, gunzip: bool = False, extra: dict | None = None, force: bool = False
    ) -> tuple[Path, str]:
        """-> (path, 'cached'|'downloaded'|'repaired'|'updated'). With `force`, a valid
        file is refetched, and replaced only once the new copy has arrived intact."""
        f = self.dir / name
        if self.valid(name) and not force:
            return f, "cached"
        status = (
            "updated"
            if force and self.valid(name)
            else ("repaired" if f.exists() or name in self.entries else "downloaded")
        )
        blob = get(url)
        if gunzip:
            try:
                blob = gzip.decompress(blob)
            except (OSError, EOFError) as e:  # BadGzipFile is an OSError; truncated stream -> EOFError
                raise HTTPError(url, None, f"corrupt gzip payload: {e}") from e
        tmp = f.with_suffix(f.suffix + ".part")
        tmp.write_bytes(blob)
        tmp.replace(f)  # atomic: never leave a half-written structure behind
        self.entries[name] = {
            "url": url,
            "sha256": sha256(f),
            "bytes": f.stat().st_size,
            "fetched": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            **(extra or {}),
        }
        self.save()
        return f, status

    def save(self) -> None:
        atomic_write(self.path, json.dumps(self.entries, indent=1, sort_keys=True))
