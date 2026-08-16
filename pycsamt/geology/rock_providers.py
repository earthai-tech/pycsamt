# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Pluggable rock-property sources for :class:`~pycsamt.geology.lithology.RockDatabase`.

No public, machine-readable service currently maps a rock or lithology
name directly to a resistivity range (see :mod:`pycsamt.geology.rock_library`
for why the built-in table is a literature compilation rather than a live
fetch). This module instead gives
:meth:`~pycsamt.geology.lithology.RockDatabase.from_url` and
:meth:`~pycsamt.geology.lithology.RockDatabase.from_provider` a small,
source-agnostic contract, so a project- or organisation-controlled
endpoint -- an internal REST API, a JSON file on a shared drive or object
store, or a public service, should one appear -- can be plugged in without
touching :mod:`pycsamt.geology.lithology`.

Every provider implements :class:`RockPropertyProvider`: a zero-argument
``fetch()`` returning ``(entries, metadata)``, a list of
:class:`~pycsamt.geology.lithology.RockEntry` plus a small provenance
dictionary recording where they came from. :class:`RemoteRockPropertyProvider`
caches successful fetches under ``~/.pycsamt/rock_db`` (override with the
``PYCSAMT_ROCKDB_CACHE`` environment variable or the ``cache_dir``
argument, matching the convention already used by
:mod:`pycsamt.ai._zoo`'s model cache) and, on any failure -- network error,
timeout, or a response that does not match the expected schema -- falls
back first to a stale cache entry and then to
:class:`LocalRockPropertyProvider`, rather than raising, unless
``fallback=False``.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any, Protocol, runtime_checkable
from urllib.error import URLError
from urllib.request import urlopen

from ..log.logger import get_logger
from .lithology import RockDatabase, RockEntry

__all__ = [
    "RockPropertyProvider",
    "LocalRockPropertyProvider",
    "RemoteRockPropertyProvider",
    "RockProviderFetchError",
]

_logger = get_logger(__name__)

_REQUIRED_FIELDS = ("name", "rho_min", "rho_max")
_OPTIONAL_FIELDS = ("color", "description", "code", "source")


class RockProviderFetchError(RuntimeError):
    """Raised by a provider when a fetch fails and ``fallback=False``.

    Examples
    --------
    >>> isinstance(RockProviderFetchError("network down"), RuntimeError)
    True
    """


@runtime_checkable
class RockPropertyProvider(Protocol):
    """Minimal contract every rock-property source must satisfy."""

    def fetch(self) -> tuple[list[RockEntry], dict[str, Any]]:
        """Return ``(entries, metadata)`` for a :class:`RockDatabase`."""
        ...


class LocalRockPropertyProvider:
    """Provider wrapping the bundled table or a local CSV file.

    This is what :class:`RemoteRockPropertyProvider` falls back to, and is
    also usable directly wherever a :class:`RockPropertyProvider` is
    expected but no remote source is involved.

    Parameters
    ----------
    csv_path : path-like, optional
        When given, entries come from
        :meth:`~pycsamt.geology.lithology.RockDatabase.from_csv`. When
        omitted, entries come from
        :meth:`~pycsamt.geology.lithology.RockDatabase.default`.
    """

    def __init__(self, csv_path: str | Path | None = None) -> None:
        self.csv_path = Path(csv_path) if csv_path is not None else None

    def fetch(self) -> tuple[list[RockEntry], dict[str, Any]]:
        db = (
            RockDatabase.from_csv(self.csv_path)
            if self.csv_path is not None
            else RockDatabase.default()
        )
        return list(db.entries), dict(db.metadata)


class RemoteRockPropertyProvider:
    """Fetch rock entries as JSON from *url*, with caching and fallback.

    Parameters
    ----------
    url : str
        Any URL :func:`urllib.request.urlopen` can open (``http(s)://``,
        ``file://``, ...), serving a JSON array of objects with at least
        ``name``, ``rho_min``, ``rho_max``.
    cache_dir : path-like, optional
        Override the local cache location. Defaults to
        ``$PYCSAMT_ROCKDB_CACHE`` or ``~/.pycsamt/rock_db``.
    ttl_seconds : float
        Reuse a cached response younger than this many seconds instead of
        re-fetching.
    timeout : float
        Network timeout in seconds for the fetch itself.
    force : bool
        Re-fetch even if a fresh cache entry exists.
    fallback : bool
        On failure, fall back to a stale cache entry (if any) and then to
        :class:`LocalRockPropertyProvider`, instead of raising
        :class:`RockProviderFetchError`.
    """

    def __init__(
        self,
        url: str,
        *,
        cache_dir: str | Path | None = None,
        ttl_seconds: float = 86400.0,
        timeout: float = 10.0,
        force: bool = False,
        fallback: bool = True,
    ) -> None:
        self.url = url
        self.cache_dir = cache_dir
        self.ttl_seconds = float(ttl_seconds)
        self.timeout = float(timeout)
        self.force = bool(force)
        self.fallback = bool(fallback)

    def fetch(self) -> tuple[list[RockEntry], dict[str, Any]]:
        cache_path = _cache_path_for_url(self.url, self._resolve_cache_dir())

        if not self.force:
            cached = _read_cache(cache_path, max_age_seconds=self.ttl_seconds)
            if cached is not None:
                entries = _entries_from_payload(cached["entries"])
                return entries, {
                    "origin": "url",
                    "url": self.url,
                    "cache_path": str(cache_path),
                    "cache_hit": True,
                    "fetched_at": cached["fetched_at"],
                }

        try:
            payload = _fetch_json(self.url, timeout=self.timeout)
            entries = _entries_from_payload(payload)
        except Exception as exc:  # network, timeout, or bad schema
            _logger.warning(
                "RockDatabase.from_url: fetch from %r failed (%s); "
                "falling back.",
                self.url,
                exc,
            )
            return self._on_failure(cache_path, exc)

        fetched_at = time.time()
        _write_cache(cache_path, payload, fetched_at)
        return entries, {
            "origin": "url",
            "url": self.url,
            "cache_path": str(cache_path),
            "cache_hit": False,
            "fetched_at": fetched_at,
        }

    # ------------------------------------------------------------------

    def _resolve_cache_dir(self) -> Path:
        if self.cache_dir is not None:
            return Path(self.cache_dir)
        env = os.environ.get("PYCSAMT_ROCKDB_CACHE")
        if env:
            return Path(env)
        return Path.home() / ".pycsamt" / "rock_db"

    def _on_failure(
        self,
        cache_path: Path,
        exc: Exception,
    ) -> tuple[list[RockEntry], dict[str, Any]]:
        if not self.fallback:
            raise RockProviderFetchError(
                f"Failed to fetch rock database from {self.url!r}: {exc}"
            ) from exc

        stale = _read_cache(cache_path, max_age_seconds=float("inf"))
        if stale is not None:
            entries = _entries_from_payload(stale["entries"])
            return entries, {
                "origin": "url-stale-cache",
                "url": self.url,
                "cache_path": str(cache_path),
                "cache_hit": True,
                "fetched_at": stale["fetched_at"],
                "error": str(exc),
            }

        entries, metadata = LocalRockPropertyProvider().fetch()
        metadata = {
            **metadata,
            "origin": "default-fallback",
            "attempted_url": self.url,
            "error": str(exc),
        }
        return entries, metadata


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _cache_path_for_url(url: str, cache_dir: Path) -> Path:
    digest = hashlib.sha256(url.encode("utf-8")).hexdigest()[:32]
    return cache_dir / f"{digest}.json"


def _read_cache(path: Path, *, max_age_seconds: float) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    fetched_at = payload.get("fetched_at")
    entries = payload.get("entries")
    if not isinstance(fetched_at, (int, float)) or not isinstance(
        entries, list
    ):
        return None
    if time.time() - float(fetched_at) > max_age_seconds:
        return None
    return {"fetched_at": float(fetched_at), "entries": entries}


def _write_cache(path: Path, entries: list[dict], fetched_at: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(
        json.dumps({"fetched_at": fetched_at, "entries": entries}),
        encoding="utf-8",
    )
    os.replace(tmp, path)


def _fetch_json(url: str, *, timeout: float) -> list[dict]:
    try:
        with urlopen(url, timeout=timeout) as response:  # noqa: S310
            raw = response.read()
    except URLError as exc:
        raise RockProviderFetchError(
            f"Could not reach {url!r}: {exc}"
        ) from exc
    try:
        payload = json.loads(raw.decode("utf-8"))
    except ValueError as exc:
        raise RockProviderFetchError(
            f"Response from {url!r} is not valid JSON: {exc}"
        ) from exc
    if not isinstance(payload, list):
        raise RockProviderFetchError(
            f"Response from {url!r} must be a JSON array of rock entries."
        )
    return payload


def _entries_from_payload(payload: list[dict]) -> list[RockEntry]:
    entries: list[RockEntry] = []
    for i, row in enumerate(payload):
        if not isinstance(row, dict):
            raise RockProviderFetchError(
                f"Entry {i} is not a JSON object: {row!r}"
            )
        missing = [f for f in _REQUIRED_FIELDS if f not in row]
        if missing:
            raise RockProviderFetchError(
                f"Entry {i} is missing required field(s) {missing}: {row!r}"
            )
        entries.append(
            RockEntry(
                name=str(row["name"]),
                rho_min=float(row["rho_min"]),
                rho_max=float(row["rho_max"]),
                color=str(row.get("color", "#AAAAAA")),
                description=str(row.get("description", "")),
                code=int(row.get("code", i + 1)),
                source=str(row.get("source", "")),
            )
        )
    return entries
