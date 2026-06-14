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
    _DIR = os.path.join(os.path.expanduser("~"), ".pycsamt", "web_cache")
    cache: diskcache.Cache = diskcache.Cache(_DIR)
except ImportError:
    cache = None  # type: ignore[assignment]

SESSION_TTL: int = 86_400  # seconds — cached Sites expire after 24 hours

# In-memory fallback — used when diskcache is not installed.
# Works fine for single-process development servers; install diskcache for
# multi-worker / background-callback production use.
_MEM: dict[str, Any] = {}


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


def has_diskcache() -> bool:
    """True when diskcache is available (background callbacks can run safely)."""
    return cache is not None
