"""Content-addressed disk cache for pipeline step outputs.

A cache hit for step *i* requires the exact same fingerprint of the sites
flowing into it, the exact same step code, and the exact same params as a
previous run — a Make/Docker-layer-style chain hash. This is what makes a
crashed-and-rerun ``Pipeline.run(..., cache=True)`` call "resume": every
already-completed step replays from cache instead of recomputing, and only
the step that was interrupted (and anything after it) actually runs.

Deliberately simpler than :class:`pycsamt.forward.maxwell.cache.MaxwellResultCache`
— no cross-process locking, no size-based eviction. This is a disposable,
regeneratable local cache, not scientifically load-bearing output, so that
extra complexity isn't worth it here. See :doc:`/user_guide/pipeline/caching`
for the caveats this implies (non-deterministic steps are not safe to cache).
"""

from __future__ import annotations

import hashlib
import json
import os
import pickle
import uuid
import warnings
from pathlib import Path
from typing import Any

import numpy as np

__all__ = ["StepCache", "fingerprint_sites", "chain_key"]

_MISS = object()  # never pickled — always the in-memory singleton, safe to check via `is`


class _DiagnosticOkType:
    """Marker stored in place of Sites for a cached diagnostic-step hit.

    Plain ``object()`` sentinels don't survive a pickle round-trip by
    identity (unpickling constructs a *new* instance) — a real bug this
    module's own test suite caught. ``__reduce__`` fixes that by telling
    pickle to always resolve back to the single module-level instance,
    the same trick :data:`None`/:data:`NotImplemented` rely on.
    """

    def __reduce__(self):
        return (_diagnostic_ok_singleton, ())

    def __repr__(self) -> str:
        return "DIAGNOSTIC_OK"


def _diagnostic_ok_singleton() -> _DiagnosticOkType:
    return DIAGNOSTIC_OK


DIAGNOSTIC_OK = _DiagnosticOkType()


def fingerprint_sites(sites: Any) -> str:
    """Deterministic content hash of a Sites-like collection.

    Tries the real ``Sites`` shape first — iterates stations, hashing each
    one's name plus its ``freq``/``z`` arrays, in iteration order. Order is
    part of the fingerprint on purpose: under-caching on a harmless reorder
    is a safer failure mode than treating two differently-ordered (and
    potentially differently-behaving) inputs as identical.

    Falls back to hashing the pickled object for anything that doesn't
    match that shape (plain test doubles, custom plugin-defined sites-like
    objects) — this keeps the cache usable without special-casing every
    possible "sites" representation.
    """
    h = hashlib.sha256()
    try:
        touched = False
        for site in sites:
            touched = True
            h.update(str(getattr(site, "name", "")).encode("utf-8"))
            for attr in ("freq", "z"):
                arr = getattr(site, attr, None)
                if arr is not None:
                    h.update(np.asarray(arr).tobytes())
        if not touched and len(sites) > 0:  # non-empty but not iterable-as-sites
            raise TypeError
    except Exception:
        h = hashlib.sha256(pickle.dumps(sites))
    return h.hexdigest()


def chain_key(upstream_fp: str, code: str, params: dict) -> str:
    """sha256 of (upstream fingerprint, step code, normalized params).

    ``default=str`` in the JSON dump means an unhashable/odd param type
    (a callable, a numpy array) never raises — it just gets stringified,
    which is deterministic even if imperfect (two different objects with
    identical ``repr`` would collide). Documented, not silently hidden.
    """
    payload = json.dumps(
        {"upstream": upstream_fp, "code": code, "params": params},
        sort_keys=True,
        default=str,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


class StepCache:
    """Sharded, content-addressed disk cache for pipeline step outputs.

    Layout: ``root/<key[:2]>/<key>.joblib``. Writes are atomic (temp file in
    the same directory, then ``os.replace``), so a killed process never
    leaves a half-written entry visible. A corrupt or unreadable entry is
    treated as a miss — warns, never crashes the run — the same "one bad
    entry must not break everything else" principle
    :func:`pycsamt.pipeline.discover_plugins` already established.
    """

    def __init__(self, root: str | Path | None = None) -> None:
        self.root = Path(root) if root is not None else default_cache_root()

    def _shard_dir(self, key: str) -> Path:
        return self.root / key[:2]

    def _entry_path(self, key: str) -> Path:
        return self._shard_dir(key) / f"{key}.joblib"

    def get(self, key: str) -> Any:
        """Return the cached value for *key*, or the ``_MISS`` sentinel."""
        path = self._entry_path(key)
        if not path.exists():
            return _MISS
        import joblib

        try:
            return joblib.load(path)
        except Exception as exc:
            warnings.warn(
                f"pyCSAMT pipeline cache entry {key!r} is unreadable "
                f"({type(exc).__name__}: {exc}) — treating as a cache miss.",
                stacklevel=2,
            )
            return _MISS

    def put(self, key: str, value: Any) -> None:
        """Atomically store *value* under *key*."""
        import joblib

        shard = self._shard_dir(key)
        shard.mkdir(parents=True, exist_ok=True)
        final_path = self._entry_path(key)
        token = uuid.uuid4().hex
        tmp_path = shard / f".{key}.{token}.tmp"
        try:
            joblib.dump(value, tmp_path)
            os.replace(tmp_path, final_path)
        finally:
            try:
                tmp_path.unlink()
            except FileNotFoundError:
                pass

    def clear(self) -> None:
        """Remove every entry (and the root directory itself) from disk."""
        import shutil

        shutil.rmtree(self.root, ignore_errors=True)


def default_cache_root() -> Path:
    return Path.home() / ".pycsamt" / "pipeline_cache"
