# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/session.py — save / restore the MapView workbench session.

A session snapshot captures *view state*, not raw EDI data: active view
(map / 3-D), inspector controls, active lines, masked stations, theme,
and the loaded survey's station metadata (for display only). The
station metadata alone can't rebuild the server-side ``MapView`` object
cached in ``cache.py`` (that requires the actual EDI bytes, which never
leave the temp directory of the browser tab that uploaded them), so a
restored session always needs a fresh "Load lines" afterwards to bring
the map/3-D canvas back — the same trade-off the sibling web app makes
(see ``pycsamt.app.web.callbacks._session``, which this module mirrors).

Responsibilities
----------------
- Open / close the Session offcanvas.
- Auto-snapshot the view-state stores to localStorage on every change.
- Download session JSON (BTN_SESSION_SAVE).
- Restore session from an uploaded JSON file (SESSION_UL).
- Restore from the browser's localStorage snapshot (BTN_SESSION_RESTORE).
- Clear the localStorage snapshot (BTN_SESSION_CLEAR).
"""

from __future__ import annotations

import base64
import json
from datetime import datetime, timezone

from dash import Input, Output, State, dcc, html, no_update

from .._ids import IDs

_VERSION = "1.0"


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _ok_icon(msg: str) -> html.Span:
    return html.Span(
        [
            html.I(
                className="bi bi-check-circle-fill me-1",
                style={"color": "var(--mv-accent)"},
            ),
            msg,
        ],
    )


def _err_icon(msg: str) -> html.Span:
    return html.Span(
        [
            html.I(
                className="bi bi-x-circle-fill me-1",
                style={"color": "#e05252"},
            ),
            msg,
        ],
    )


def _build_snapshot(
    view, controls, lines, line_filter, masked, theme, store_data, note
) -> dict:
    return {
        "version": _VERSION,
        "app": "mapview",
        "saved_at": _now_iso(),
        "note": note or "",
        "view": view,
        "controls": controls,
        "lines": lines,
        "line_filter": line_filter,
        "masked": masked,
        "theme": theme,
        "store_data": store_data,
    }


def _validate_snapshot(snap: dict) -> None:
    """Raise ValueError if the snapshot is missing mandatory keys."""
    if not isinstance(snap, dict):
        raise ValueError("Not a valid session file.")
    ver = str(snap.get("version", ""))
    if not ver.startswith("1"):
        raise ValueError(f"Unsupported session version '{ver}'.")
    if "view" not in snap:
        raise ValueError("Missing 'view' key — file may be corrupt.")


def register_session(app) -> None:

    # ── 1. Open / close session offcanvas ─────────────────────────────────
    @app.callback(
        Output(IDs.CANVAS_SESSION, "is_open"),
        Input(IDs.BTN_SESSION, "n_clicks"),
        State(IDs.CANVAS_SESSION, "is_open"),
        prevent_initial_call=True,
    )
    def _toggle(_n, is_open):
        return not is_open

    # ── 2. Auto-snapshot to localStorage on every view-state change ───────
    @app.callback(
        Output(IDs.SESSION_SNAPSHOT, "data"),
        Output(IDs.SESSION_AUTOSAVE, "children"),
        Input(IDs.STORE_VIEW, "data"),
        Input(IDs.STORE_CONTROLS, "data"),
        Input(IDs.STORE_LINES, "data"),
        Input(IDs.STORE_LINE_FILTER, "data"),
        Input(IDs.STORE_MASKED, "data"),
        Input(IDs.STORE_THEME, "data"),
        Input(IDs.STORE_DATA, "data"),
        State(IDs.SESSION_NOTE, "value"),
        prevent_initial_call=True,
    )
    def _auto_snapshot(
        view, controls, lines, line_filter, masked, theme, store_data, note
    ):
        if not store_data or not store_data.get("n_stations"):
            return no_update, no_update
        snap = _build_snapshot(
            view,
            controls,
            lines,
            line_filter,
            masked,
            theme,
            store_data,
            note,
        )
        n_sta = store_data.get("n_stations", 0)
        chip = html.Span(
            [
                html.I(
                    className="bi bi-check-circle-fill me-1",
                    style={"color": "var(--mv-accent)", "fontSize": "10px"},
                ),
                f"Auto-saved · {n_sta} stations",
            ],
        )
        return snap, chip

    # ── 3. Download session JSON ────────────────────────────────────────────
    @app.callback(
        Output(IDs.SESSION_DL, "data"),
        Output(IDs.SESSION_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_SESSION_SAVE, "n_clicks"),
        State(IDs.STORE_VIEW, "data"),
        State(IDs.STORE_CONTROLS, "data"),
        State(IDs.STORE_LINES, "data"),
        State(IDs.STORE_LINE_FILTER, "data"),
        State(IDs.STORE_MASKED, "data"),
        State(IDs.STORE_THEME, "data"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.SESSION_NOTE, "value"),
        prevent_initial_call=True,
    )
    def _download(
        n, view, controls, lines, line_filter, masked, theme, store_data, note
    ):
        if not n:
            return no_update, no_update
        if not store_data or not store_data.get("n_stations"):
            return no_update, _err_icon("No survey loaded — nothing to save.")
        snap = _build_snapshot(
            view,
            controls,
            lines,
            line_filter,
            masked,
            theme,
            store_data,
            note,
        )
        n_sta = store_data.get("n_stations", 0)
        fname = (
            f"mapview_session_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        )
        return (
            dcc.send_string(json.dumps(snap, indent=2), fname),
            _ok_icon(f"Session downloaded  ({n_sta} stations)."),
        )

    # ── 4. Restore from an uploaded JSON file ──────────────────────────────
    @app.callback(
        Output(IDs.STORE_VIEW, "data", allow_duplicate=True),
        Output(IDs.STORE_CONTROLS, "data", allow_duplicate=True),
        Output(IDs.STORE_LINES, "data", allow_duplicate=True),
        Output(IDs.STORE_LINE_FILTER, "data", allow_duplicate=True),
        Output(IDs.STORE_MASKED, "data", allow_duplicate=True),
        Output(IDs.STORE_THEME, "data", allow_duplicate=True),
        Output(IDs.SESSION_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.SESSION_UL, "contents"),
        State(IDs.SESSION_UL, "filename"),
        prevent_initial_call=True,
    )
    def _upload_restore(contents, filename):
        if not contents:
            return (no_update,) * 7
        try:
            _, enc = contents.split(",", 1)
            snap = json.loads(base64.b64decode(enc).decode())
            _validate_snapshot(snap)
            n_sta = (snap.get("store_data") or {}).get("n_stations", "?")
            msg = _ok_icon(
                f"Session restored from {filename}  ({n_sta} stations). "
                "Reload EDI files to re-enable the map/3-D canvas."
            )
            return (
                snap.get("view") or "map",
                snap.get("controls") or {},
                snap.get("lines") or {},
                snap.get("line_filter"),
                snap.get("masked") or [],
                snap.get("theme") or "light",
                msg,
            )
        except Exception as exc:
            return (no_update,) * 6 + (
                _err_icon(f"Could not parse {filename}: {exc}"),
            )

    # ── 5. Restore from the browser's localStorage snapshot ───────────────
    @app.callback(
        Output(IDs.STORE_VIEW, "data", allow_duplicate=True),
        Output(IDs.STORE_CONTROLS, "data", allow_duplicate=True),
        Output(IDs.STORE_LINES, "data", allow_duplicate=True),
        Output(IDs.STORE_LINE_FILTER, "data", allow_duplicate=True),
        Output(IDs.STORE_MASKED, "data", allow_duplicate=True),
        Output(IDs.STORE_THEME, "data", allow_duplicate=True),
        Output(IDs.SESSION_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_SESSION_RESTORE, "n_clicks"),
        State(IDs.SESSION_SNAPSHOT, "data"),
        prevent_initial_call=True,
    )
    def _browser_restore(n, snap):
        if not n:
            return (no_update,) * 7
        if not snap:
            return (no_update,) * 6 + (
                _err_icon("No browser snapshot found — save a session first."),
            )
        try:
            _validate_snapshot(snap)
            saved_at = snap.get("saved_at", "unknown time")
            n_sta = (snap.get("store_data") or {}).get("n_stations", "?")
            note = snap.get("note", "")
            detail = f" · {note}" if note else ""
            msg = _ok_icon(
                f"Restored snapshot from {saved_at}{detail} "
                f"({n_sta} stations). "
                "Reload EDI files to re-enable the map/3-D canvas."
            )
            return (
                snap.get("view") or "map",
                snap.get("controls") or {},
                snap.get("lines") or {},
                snap.get("line_filter"),
                snap.get("masked") or [],
                snap.get("theme") or "light",
                msg,
            )
        except Exception as exc:
            return (no_update,) * 6 + (_err_icon(str(exc)),)

    # ── 6. Clear the localStorage snapshot ─────────────────────────────────
    @app.callback(
        Output(IDs.SESSION_SNAPSHOT, "data", allow_duplicate=True),
        Output(IDs.SESSION_AUTOSAVE, "children", allow_duplicate=True),
        Output(IDs.SESSION_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_SESSION_CLEAR, "n_clicks"),
        prevent_initial_call=True,
    )
    def _clear(n):
        if not n:
            return (no_update,) * 3
        return None, "", _ok_icon("Browser snapshot cleared.")
