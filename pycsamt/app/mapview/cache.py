# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Server-side session storage for the map-view platform.

The dcc stores keep only JSON-serialisable station records.  The heavy
:class:`~pycsamt.map.MapView` (which holds parsed EDI objects) lives
here, keyed by browser session ID — mirroring ``app/web/cache.py``.

A process-global *seed* slot lets :func:`pycsamt.map.MapView.launch`
hand an already-built view to the freshly started server so the first
session opens pre-populated.
"""

from __future__ import annotations

from typing import Any

# Reuse the web app's diskcache instance / RAM fallback so both apps
# share one storage backend and TTL policy.
try:
    from pycsamt.app.web.cache import (
        cache_get as _web_get,
        cache_set as _web_set,
    )
except Exception:  # pragma: no cover - web extra missing
    _MEM: dict[str, Any] = {}

    def _web_set(session_id: str, obj: Any) -> None:
        if session_id:
            _MEM[session_id] = obj

    def _web_get(session_id: str) -> Any:
        return _MEM.get(session_id) if session_id else None


_VIEW_PREFIX = "mapview::"
_SEED: dict[str, Any] = {"view": None}


def set_view(session_id: str, view: Any) -> None:
    """Store the :class:`MapView` for *session_id*."""
    if session_id:
        _web_set(_VIEW_PREFIX + session_id, view)


def get_view(session_id: str) -> Any:
    """Return the cached :class:`MapView`, or ``None`` on miss."""
    if not session_id:
        return None
    return _web_get(_VIEW_PREFIX + session_id)


def set_seed(view: Any) -> None:
    """Stash a view for the first session to adopt (launch handoff)."""
    _SEED["view"] = view


def take_seed() -> Any:
    """Return and clear the pending seed view (one-shot)."""
    view = _SEED.get("view")
    _SEED["view"] = None
    return view
