# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
web/cache.py — shared diskcache instance for the Dash web application.

The cache stores heavy Python objects (Sites, etc.) keyed by browser
session ID so that each tab gets its own isolated data context.

Usage::

    from pycsamt.app.web.cache import cache, SESSION_TTL
    cache.set(session_id, sites, expire=SESSION_TTL)
    sites = cache.get(session_id)   # None if expired or not yet loaded
"""

from __future__ import annotations

import os
from typing import Any

try:
    import diskcache
    _BASE_DIR = os.environ.get(
        "PYCSAMT_DATA",
        os.path.join(os.path.expanduser("~"), ".pycsamt"),
    )
    _DIR = os.path.join(_BASE_DIR, "web_cache")
    cache: diskcache.Cache = diskcache.Cache(_DIR)
except Exception:
    cache = None  # type: ignore[assignment]

SESSION_TTL: int = 86_400  # seconds — cached Sites expire after 24 hours

# In-memory fallback — used when diskcache is not installed.
# Works fine for single-process development servers; install diskcache for
# multi-worker / background-callback production use.
_MEM: dict[str, Any] = {}

# True only when the full Dash background-callback chain works:
# diskcache installed + DiskcacheManager instantiable (requires multiprocess).
# Callbacks check this before adding `background=True` to avoid
# MissingLongCallbackManagerError when multiprocess is absent.
HAS_BG_MANAGER: bool = False
bg_manager = None  # type: ignore[assignment]
if cache is not None:
    try:
        from dash import DiskcacheManager as _DM  # noqa: F401
        bg_manager = _DM(cache)
        HAS_BG_MANAGER = True
    except Exception:
        pass


# ──────────────────────────────────────────────────────────────────────────────
# Convenience wrappers used by every callbacks submodule
# ──────────────────────────────────────────────────────────────────────────────

def cache_set(session_id: str, sites: Any) -> None:
    """Store *sites* keyed by *session_id* (disk cache when available, else RAM)."""
    if not session_id:
        return
    _MEM[session_id] = sites          # always keep in-memory copy
    if cache is not None:
        try:
            cache.set(session_id, sites, expire=SESSION_TTL)
        except Exception:
            pass


def cache_get(session_id: str) -> Any:
    """Return the Sites object for *session_id*, or ``None`` on miss/error."""
    if not session_id:
        return None
    if cache is not None:
        try:
            result = cache.get(session_id)
            if result is not None:
                return result
        except Exception:
            pass
    return _MEM.get(session_id)


def cache_merge_sites(session_id: str, new_sites: Any) -> Any:
    """
    Merge *new_sites* into the Sites object already cached for *session_id*.

    Deduplicates by station name: if a station appears in both, the
    version from *new_sites* wins.  Returns the merged Sites object
    (also stored in cache as a side-effect).
    """
    existing = cache_get(session_id)
    if existing is None:
        cache_set(session_id, new_sites)
        return new_sites
    try:
        from pycsamt.site.base import Sites
        # Build a name→edi dict for deduplication (new overwrites old on clash)
        merged: dict = {}
        for edi in existing.edic:
            name = getattr(edi, "station", None) or str(id(edi))
            merged[name] = edi
        for edi in new_sites.edic:
            name = getattr(edi, "station", None) or str(id(edi))
            merged[name] = edi
        combined = Sites(list(merged.values()))
        cache_set(session_id, combined)
        return combined
    except Exception:
        # If merging fails for any reason, keep the new data so the user
        # is not left with a completely broken state.
        cache_set(session_id, new_sites)
        return new_sites


def has_diskcache() -> bool:
    """True when the full background-callback chain is usable (diskcache + multiprocess)."""
    return HAS_BG_MANAGER


# ── Forward-response cache (keyed by session + dimension) ────────────────────

def cache_set_fwd(session_id: str, dim: str, response: Any) -> None:
    """Store a ForwardResponse object for the given session and dimension ('2d' or '3d')."""
    if not session_id:
        return
    key = f"{session_id}_fwd_{dim}"
    _MEM[key] = response
    if cache is not None:
        try:
            cache.set(key, response, expire=SESSION_TTL)
        except Exception:
            pass


def cache_get_fwd(session_id: str, dim: str) -> Any:
    """Return the ForwardResponse for *session_id*/*dim*, or ``None`` on miss."""
    if not session_id:
        return None
    key = f"{session_id}_fwd_{dim}"
    if cache is not None:
        try:
            result = cache.get(key)
            if result is not None:
                return result
        except Exception:
            pass
    return _MEM.get(key)


# ── Inversion-result cache (consumed by 3D map / interpretation pages) ───────

def cache_set_inversion_result(session_id: str, result: Any) -> None:
    """Store the latest inversion result for *session_id*."""
    if not session_id:
        return
    key = f"{session_id}_inversion_result"
    _MEM[key] = result
    if cache is not None:
        try:
            cache.set(key, result, expire=SESSION_TTL)
        except Exception:
            pass


def cache_get_inversion_result(session_id: str) -> Any:
    """Return the latest cached inversion result for *session_id*."""
    if not session_id:
        return None
    key = f"{session_id}_inversion_result"
    if cache is not None:
        try:
            result = cache.get(key)
            if result is not None:
                return result
        except Exception:
            pass
    return _MEM.get(key)


# ── Frequency-editor cache (edited Sites after confidence/band/regrid edits) ──

def cache_set_freq_edit(session_id: str, sites: Any) -> None:
    """Store edited Sites from the frequency editor for *session_id*."""
    if not session_id:
        return
    key = f"{session_id}_freq_edit"
    _MEM[key] = sites
    if cache is not None:
        try:
            cache.set(key, sites, expire=SESSION_TTL)
        except Exception:
            pass


def cache_get_freq_edit(session_id: str) -> Any:
    """Return the frequency-edited Sites for *session_id*, or ``None``."""
    if not session_id:
        return None
    key = f"{session_id}_freq_edit"
    if cache is not None:
        try:
            result = cache.get(key)
            if result is not None:
                return result
        except Exception:
            pass
    return _MEM.get(key)
