# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/_session.py

Responsibilities
----------------
- Initialise the per-tab session UUID on first load.
- Open / close the Session offcanvas.
- Auto-snapshot workflow stores to localStorage on every STORE_DATA change.
- Download session JSON (BTN_SESSION_SAVE).
- Restore session from uploaded JSON (SESSION_UL).
- Restore from browser localStorage snapshot (BTN_SESSION_RESTORE).
- Clear snapshot (BTN_SESSION_CLEAR).
"""

from __future__ import annotations

import base64
import json
import uuid
from datetime import datetime, timezone

from dash import (
    Input,
    Output,
    State,
    dcc,
    html,
    no_update,
)

from pycsamt.app.web.layout import IDs

# ── Keys included in every session snapshot ────────────────────────────────────
_SESSION_KEYS = {
    "store_data": IDs.STORE_DATA,
    "active_lines": IDs.STORE_ACTIVE_LINES,
    "corr_store": IDs.CORR_STORE,
    "elevation_raw": "tool-elev-raw-store",
    "elevation_corrected": "tool-elev-corrected-store",
}

_VERSION = "2.0"
_PYCSAMT_VERSION = "v2"


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _ok_icon(msg: str) -> html.Span:
    return html.Span(
        [
            html.I(
                className="bi bi-check-circle-fill me-1",
                style={"color": "var(--green)"},
            ),
            msg,
        ],
        className="small",
    )


def _err_icon(msg: str) -> html.Span:
    return html.Span(
        [
            html.I(
                className="bi bi-x-circle-fill me-1",
                style={"color": "var(--red)"},
            ),
            msg,
        ],
        className="small",
    )


def register_session(app) -> None:

    # ── 1. Initialise session UUID ────────────────────────────────────────────
    @app.callback(
        Output(IDs.SESSION_ID, "data"),
        Input(IDs.SESSION_ID, "data"),
        prevent_initial_call=False,
    )
    def _init_session_id(existing_id):
        if existing_id:
            return no_update
        return str(uuid.uuid4())

    # ── 2. Open / close session offcanvas ─────────────────────────────────────
    @app.callback(
        Output(IDs.SESSION_CANVAS, "is_open"),
        Input(IDs.BTN_SESSION_OPEN, "n_clicks"),
        State(IDs.SESSION_CANVAS, "is_open"),
        prevent_initial_call=True,
    )
    def _toggle_session_canvas(n, is_open):
        return not is_open if n else no_update

    # ── 3. Auto-snapshot to localStorage whenever STORE_DATA changes ──────────
    @app.callback(
        Output(IDs.SESSION_SNAPSHOT, "data"),
        Output(IDs.SESSION_AUTOSAVE, "children"),
        Input(IDs.STORE_DATA, "data"),
        Input("tool-elev-raw-store", "data"),
        Input("tool-elev-corrected-store", "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.CORR_STORE, "data"),
        State(IDs.SESSION_NOTE, "value"),
        prevent_initial_call=True,
    )
    def _auto_snapshot(
        store_data, elev_raw, elev_corr, active_lines, corr_store, note
    ):
        if not store_data:
            return no_update, no_update
        snapshot = {
            "version": _VERSION,
            "pycsamt": _PYCSAMT_VERSION,
            "saved_at": _now_iso(),
            "note": note or "",
            "store_data": store_data,
            "active_lines": active_lines,
            "corr_store": corr_store,
            "elevation_raw": elev_raw,
            "elevation_corrected": elev_corr,
        }
        n_sta = (store_data or {}).get("n_stations", 0)
        chip = html.Span(
            [
                html.I(
                    className="bi bi-check-circle-fill me-1",
                    style={"color": "var(--green)", "fontSize": "10px"},
                ),
                f"Auto-saved · {n_sta} stations",
            ],
            className="session-autosave-chip",
        )
        return snapshot, chip

    # ── 4. Download session JSON ──────────────────────────────────────────────
    @app.callback(
        Output(IDs.SESSION_DL, "data"),
        Output(IDs.SESSION_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_SESSION_SAVE, "n_clicks"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.CORR_STORE, "data"),
        State("tool-elev-raw-store", "data"),
        State("tool-elev-corrected-store", "data"),
        State(IDs.SESSION_NOTE, "value"),
        prevent_initial_call=True,
    )
    def _session_download(
        n, store_data, active_lines, corr, elev_raw, elev_corr, note
    ):
        if not n:
            return no_update, no_update
        if not store_data:
            return no_update, _err_icon("No survey loaded — nothing to save.")
        snapshot = {
            "version": _VERSION,
            "pycsamt": _PYCSAMT_VERSION,
            "saved_at": _now_iso(),
            "note": note or "",
            "store_data": store_data,
            "active_lines": active_lines,
            "corr_store": corr,
            "elevation_raw": elev_raw,
            "elevation_corrected": elev_corr,
        }
        n_sta = (store_data or {}).get("n_stations", 0)
        fname = (
            f"pycsamt_session_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        )
        return (
            dcc.send_string(json.dumps(snapshot, indent=2), fname),
            _ok_icon(f"Session downloaded  ({n_sta} stations)."),
        )

    # ── 5. Upload + restore session from JSON file ────────────────────────────
    @app.callback(
        Output(IDs.STORE_DATA, "data", allow_duplicate=True),
        Output(IDs.STORE_ACTIVE_LINES, "data", allow_duplicate=True),
        Output(IDs.CORR_STORE, "data", allow_duplicate=True),
        Output("tool-elev-raw-store", "data", allow_duplicate=True),
        Output("tool-elev-corrected-store", "data", allow_duplicate=True),
        Output(IDs.SESSION_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.SESSION_UL, "contents"),
        State(IDs.SESSION_UL, "filename"),
        prevent_initial_call=True,
    )
    def _session_upload_restore(contents, filename):
        if not contents:
            return (no_update,) * 6
        try:
            _, enc = contents.split(",", 1)
            snap = json.loads(base64.b64decode(enc).decode())
            _validate_snapshot(snap)
            n_sta = (snap.get("store_data") or {}).get("n_stations", "?")
            msg = _ok_icon(
                f"Session restored from {filename}  ({n_sta} stations). "
                "Reload EDI files to re-enable plots."
            )
            return (
                snap.get("store_data"),
                snap.get("active_lines"),
                snap.get("corr_store") or {"steps": []},
                snap.get("elevation_raw"),
                snap.get("elevation_corrected"),
                msg,
            )
        except Exception as exc:
            return (no_update,) * 5 + (
                _err_icon(f"Could not parse {filename}: {exc}"),
            )

    # ── 6. Restore from localStorage snapshot (BTN_SESSION_RESTORE) ──────────
    @app.callback(
        Output(IDs.STORE_DATA, "data", allow_duplicate=True),
        Output(IDs.STORE_ACTIVE_LINES, "data", allow_duplicate=True),
        Output(IDs.CORR_STORE, "data", allow_duplicate=True),
        Output("tool-elev-raw-store", "data", allow_duplicate=True),
        Output("tool-elev-corrected-store", "data", allow_duplicate=True),
        Output(IDs.SESSION_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_SESSION_RESTORE, "n_clicks"),
        State(IDs.SESSION_SNAPSHOT, "data"),
        prevent_initial_call=True,
    )
    def _session_browser_restore(n, snap):
        if not n:
            return (no_update,) * 6
        if not snap:
            return (no_update,) * 5 + (
                _err_icon(
                    "No browser snapshot found — save a session first."
                ),
            )
        try:
            _validate_snapshot(snap)
            saved_at = snap.get("saved_at", "unknown time")
            n_sta = (snap.get("store_data") or {}).get("n_stations", "?")
            note = snap.get("note", "")
            detail = f" · {note}" if note else ""
            msg = _ok_icon(
                f"Restored snapshot from {saved_at}{detail} ({n_sta} stations). "
                "Reload EDI files to re-enable plots."
            )
            return (
                snap.get("store_data"),
                snap.get("active_lines"),
                snap.get("corr_store") or {"steps": []},
                snap.get("elevation_raw"),
                snap.get("elevation_corrected"),
                msg,
            )
        except Exception as exc:
            return (no_update,) * 5 + (_err_icon(str(exc)),)

    # ── 7. Clear snapshot ─────────────────────────────────────────────────────
    @app.callback(
        Output(IDs.SESSION_SNAPSHOT, "data", allow_duplicate=True),
        Output(IDs.SESSION_AUTOSAVE, "children", allow_duplicate=True),
        Output(IDs.SESSION_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_SESSION_CLEAR, "n_clicks"),
        prevent_initial_call=True,
    )
    def _session_clear(n):
        if not n:
            return (no_update,) * 3
        return None, "", _ok_icon("Browser snapshot cleared.")


def _validate_snapshot(snap: dict) -> None:
    """Raise ValueError if the snapshot is missing mandatory keys."""
    if not isinstance(snap, dict):
        raise ValueError("Not a valid session file.")
    ver = snap.get("version", "")
    if not ver.startswith("2"):
        raise ValueError(f"Unsupported session version '{ver}'.")
    if "store_data" not in snap:
        raise ValueError("Missing 'store_data' key — file may be corrupt.")
