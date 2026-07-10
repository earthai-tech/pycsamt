# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.survey
==================

Persistent active-survey context and Sites cache for the CLI.

Every CLI command that needs EDI data calls :func:`resolve_survey`.
That function applies a strict resolution priority so that users only
have to specify their survey once per session:

Resolution priority
-------------------
1. **Explicit positional path** (e.g. ``EDI_DIR`` argument) — always
   re-parsed fresh, cache is updated.
2. **--survey / -S flag** — resolved path is looked up in the cache;
   rebuilt only when the source is newer than the cached copy.
3. **Active context** (``~/.pycsamt/context.json``) — same cache
   logic; rebuilt when the source directory has changed.
4. **UsageError** — with a clear hint to run ``pycsamt survey set``.

Cache layout
------------
::

    ~/.pycsamt/
      context.json           active survey path + summary metadata
      cache/
        <12-char-hash>/
          sites.pkl          pickled Sites object
          meta.json          path, mtime at cache time, station list

Cache invalidation
------------------
The source mtime is stored at cache-build time.  On every load the
current mtime of the survey source (file or directory) is compared; if
it is newer the cache is silently rebuilt.  A corrupt pickle also
triggers a silent rebuild.  ``pycsamt survey rebuild`` forces it.
"""

from __future__ import annotations

import hashlib
import json
import pickle
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import click

if TYPE_CHECKING:
    from pycsamt.site.base import Sites

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PYCSAMT_DIR  = Path.home() / ".pycsamt"
_CONTEXT_FILE = _PYCSAMT_DIR / "context.json"
_CACHE_ROOT   = _PYCSAMT_DIR / "cache"


# ---------------------------------------------------------------------------
# Cache key
# ---------------------------------------------------------------------------

def _cache_key(path: Path) -> str:
    """Return a 12-char hex key derived from the resolved absolute path."""
    return hashlib.md5(str(path.resolve()).encode()).hexdigest()[:12]


def _source_mtime(path: Path) -> float:
    """Return the most recent mtime across a file or a directory of EDIs."""
    if path.is_file():
        return path.stat().st_mtime
    # For a directory, take the newest mtime among .edi files (or the dir itself)
    edis = list(path.glob("*.edi"))
    if not edis:
        return path.stat().st_mtime
    return max(p.stat().st_mtime for p in edis)


# ---------------------------------------------------------------------------
# SurveyContext — reads/writes ~/.pycsamt/context.json
# ---------------------------------------------------------------------------

class SurveyContext:
    """Thin wrapper around ``~/.pycsamt/context.json``.

    Only stores lightweight metadata (path, station names, timestamps).
    The actual ``Sites`` object lives in the cache directory.
    """

    def __init__(self, data: dict[str, Any]) -> None:
        self._data = data

    # -- properties ----------------------------------------------------------

    @property
    def survey_path(self) -> Path:
        return Path(self._data["survey_path"])

    @property
    def cache_key(self) -> str:
        return self._data["cache_key"]

    @property
    def n_stations(self) -> int:
        return self._data.get("n_stations", 0)

    @property
    def station_names(self) -> list[str]:
        return self._data.get("station_names", [])

    @property
    def set_at(self) -> str:
        return self._data.get("set_at", "—")

    def to_dict(self) -> dict[str, Any]:
        return dict(self._data)

    # -- class methods -------------------------------------------------------

    @classmethod
    def load(cls) -> SurveyContext | None:
        """Return the active context or ``None`` if none is set."""
        if not _CONTEXT_FILE.exists():
            return None
        try:
            data = json.loads(_CONTEXT_FILE.read_text(encoding="utf-8"))
            return cls(data)
        except Exception:  # noqa: BLE001
            return None

    @classmethod
    def save(cls, path: Path, sites: Sites) -> SurveyContext:
        """Write a new context for *path* after *sites* has been built."""
        _PYCSAMT_DIR.mkdir(parents=True, exist_ok=True)
        key = _cache_key(path)
        try:
            names = [s.name for s in sites]
        except Exception:  # noqa: BLE001
            names = []
        data: dict[str, Any] = {
            "survey_path":   str(path.resolve()),
            "cache_key":     key,
            "set_at":        datetime.now().isoformat(timespec="seconds"),
            "n_stations":    len(sites),
            "station_names": names,
        }
        _CONTEXT_FILE.write_text(
            json.dumps(data, indent=2), encoding="utf-8"
        )
        return cls(data)

    @classmethod
    def clear(cls) -> None:
        """Remove the active context file (does not delete the cache)."""
        if _CONTEXT_FILE.exists():
            _CONTEXT_FILE.unlink()


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def _cache_dir(key: str) -> Path:
    return _CACHE_ROOT / key


def _cache_pkl(key: str) -> Path:
    return _cache_dir(key) / "sites.pkl"


def _cache_meta(key: str) -> Path:
    return _cache_dir(key) / "meta.json"


def _write_cache(key: str, path: Path, sites: Sites) -> None:
    """Pickle *sites* and write metadata under ``~/.pycsamt/cache/<key>/``."""
    d = _cache_dir(key)
    d.mkdir(parents=True, exist_ok=True)
    with open(_cache_pkl(key), "wb") as fh:
        pickle.dump(sites, fh, protocol=pickle.HIGHEST_PROTOCOL)
    meta = {
        "survey_path":       str(path.resolve()),
        "cached_at":         datetime.now().isoformat(timespec="seconds"),
        "cached_at_mtime":   _source_mtime(path),
        "n_stations":        len(sites),
    }
    _cache_meta(key).write_text(json.dumps(meta, indent=2), encoding="utf-8")


def _read_cache(key: str) -> Sites | None:
    """Return a cached ``Sites`` or ``None`` on any failure."""
    pkl = _cache_pkl(key)
    if not pkl.exists():
        return None
    try:
        with open(pkl, "rb") as fh:
            return pickle.load(fh)  # noqa: S301
    except Exception:  # noqa: BLE001
        return None


def _cache_is_valid(key: str, path: Path) -> bool:
    """Return ``True`` when the cache exists and the source has not changed."""
    meta_path = _cache_meta(key)
    if not meta_path.exists() or not _cache_pkl(key).exists():
        return False
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        return _source_mtime(path) <= meta["cached_at_mtime"]
    except Exception:  # noqa: BLE001
        return False


def _purge_cache(key: str) -> None:
    """Delete all cache files for *key*."""
    import shutil  # noqa: PLC0415
    d = _cache_dir(key)
    if d.exists():
        shutil.rmtree(d)


# ---------------------------------------------------------------------------
# Core: build Sites from source
# ---------------------------------------------------------------------------

def _build_sites(path: Path, verbose: int = 0) -> Sites:
    """Call ``ensure_sites`` and return a ``Sites`` object for *path*."""
    from pycsamt.emtools._core import (
        ensure_sites,  # noqa: PLC0415
    )

    if verbose >= 1:
        click.echo(f"  Parsing EDI files in {path} …", err=True)
    return ensure_sites(path, verbose=max(0, verbose - 1))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def set_survey(path: Path, *, force: bool = False, verbose: int = 0) -> Sites:
    """Build (or load cached) ``Sites`` for *path* and write the context.

    Parameters
    ----------
    path:
        EDI file or directory.  Resolved to an absolute path before hashing.
    force:
        When ``True`` rebuild the cache even if it is still valid.
    verbose:
        Verbosity forwarded to ``ensure_sites``.

    Returns
    -------
    Sites
    """
    path = path.resolve()
    key  = _cache_key(path)

    if not force and _cache_is_valid(key, path):
        sites = _read_cache(key)
        if sites is not None:
            if verbose >= 1:
                click.echo(
                    f"  Survey already cached "
                    f"({len(sites)} stations).  Use --force to rebuild.",
                    err=True,
                )
            SurveyContext.save(path, sites)
            return sites

    sites = _build_sites(path, verbose=verbose)
    _write_cache(key, path, sites)
    SurveyContext.save(path, sites)
    return sites


def resolve_survey(
    explicit: Path | None,
    *,
    fresh: bool = False,
    verbose: int = 0,
) -> Sites:
    """Resolve the active survey to a ``Sites`` object.

    Resolution priority
    -------------------
    1. *explicit* path  → always fresh, cache updated.
    2. Active context (``~/.pycsamt/context.json``) → cache used when valid.
    3. :class:`click.UsageError` with a hint.

    Parameters
    ----------
    explicit:
        Path passed directly by the user (positional arg or ``--survey``).
        ``None`` triggers fallback to the active context.
    fresh:
        Force re-parse even when a valid cache exists.
    verbose:
        Verbosity level.
    """
    if explicit is not None:
        path = explicit.resolve()
        key  = _cache_key(path)

        if not fresh and _cache_is_valid(key, path):
            sites = _read_cache(key)
            if sites is not None:
                if verbose >= 1:
                    click.echo(
                        f"  Using cached survey for {path.name} "
                        f"({len(sites)} stations).",
                        err=True,
                    )
                return sites

        # Build fresh and update cache
        sites = _build_sites(path, verbose=verbose)
        _write_cache(key, path, sites)
        return sites

    # --- fall back to active context ---
    ctx = SurveyContext.load()
    if ctx is None:
        raise click.UsageError(
            "No active survey is set.\n\n"
            "  Set one with:\n"
            "    pycsamt survey set <edi_dir>\n\n"
            "  Or pass the survey explicitly:\n"
            "    pycsamt <command> --survey <edi_dir>"
        )

    path = ctx.survey_path
    key  = ctx.cache_key

    if not fresh and _cache_is_valid(key, path):
        sites = _read_cache(key)
        if sites is not None:
            if verbose >= 1:
                click.echo(
                    f"  Active survey: {path}  "
                    f"({ctx.n_stations} stations, cached).",
                    err=True,
                )
            return sites

    # Cache stale or missing — rebuild
    if verbose >= 1:
        click.echo(
            f"  Active survey cache stale, rebuilding from {path} …",
            err=True,
        )
    sites = _build_sites(path, verbose=verbose)
    _write_cache(key, path, sites)
    SurveyContext.save(path, sites)
    return sites


def survey_summary(verbose: int = 0) -> dict[str, Any] | None:
    """Return a summary dict for the active survey, or ``None`` if unset."""
    ctx = SurveyContext.load()
    if ctx is None:
        return None
    key  = ctx.cache_key
    meta_path = _cache_meta(key)
    summary: dict[str, Any] = {
        "survey_path": str(ctx.survey_path),
        "set_at":      ctx.set_at,
        "n_stations":  ctx.n_stations,
        "station_names": ctx.station_names,
        "cache_key":   key,
        "cache_valid": _cache_is_valid(key, ctx.survey_path),
        "cache_path":  str(_cache_dir(key)),
    }
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            summary["cached_at"] = meta.get("cached_at", "—")
        except Exception:  # noqa: BLE001
            pass
    return summary


__all__ = [
    "SurveyContext",
    "resolve_survey",
    "set_survey",
    "survey_summary",
]
