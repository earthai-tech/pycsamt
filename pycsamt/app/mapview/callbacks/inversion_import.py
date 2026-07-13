# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Import a ModEM 3-D inversion result folder into the map-view session.

The folder is picked with the same native-OS-dialog pattern as the EDI
"Load lines" tab (see ``assets/modem_loader.js``, modeled on
``assets/edi_loader.js``): a ``webkitdirectory`` input is injected over
the Browse button, the (filtered) file contents are staged client-side
into ``IDs.INV_FOLDER_STORE``, and this module decodes them into a temp
directory before loading through
:meth:`pycsamt.map.MapView.from_inversion_results` — into the same
session store as EDI loading, so the 3-D fence/depth views work
identically regardless of data source.
"""

from __future__ import annotations

import base64
import os
import shutil
import tempfile
from pathlib import Path

from dash import Input, Output, State, no_update

from .._ids import IDs
from .._render import merge_views, store_from_view
from ..cache import get_view, set_view

_MODE_APPEND = "append"


def register_inversion_import(app) -> None:
    _register_confirm(app)


def _decode_to_tempdir(filenames: list[str], contents: list[str]) -> str:
    """Write staged base64 file contents to a flat temp directory."""
    tmpdir = tempfile.mkdtemp(prefix="pycsamt_modem_upload_")
    for name, content in zip(filenames, contents):
        try:
            _, b64 = content.split(",", 1)
            raw = base64.b64decode(b64)
        except Exception:  # noqa: BLE001 - skip unreadable entries
            continue
        with open(os.path.join(tmpdir, Path(name).name), "wb") as fh:
            fh.write(raw)
    return tmpdir


def _register_confirm(app) -> None:
    @app.callback(
        Output(IDs.STORE_DATA, "data", allow_duplicate=True),
        Output(IDs.INV_STATUS, "children"),
        Output(IDs.MODAL_LOAD, "is_open", allow_duplicate=True),
        Output(IDs.DATA_BADGE_TEXT, "children", allow_duplicate=True),
        Output(IDs.DATA_BADGE, "className", allow_duplicate=True),
        Input(IDs.BTN_INV_CONFIRM, "n_clicks"),
        State(IDs.INV_FOLDER_STORE, "data"),
        State(IDs.CK_INV_KNOWN_STA, "value"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.LOAD_MODE_STORE, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def confirm(_n, staged, use_known, session_id, mode, theme):
        filenames = (staged or {}).get("filenames") or []
        contents = (staged or {}).get("contents") or []
        if not filenames:
            return (
                no_update,
                "⚠ Browse to a ModEM results folder first.",
                no_update,
                no_update,
                no_update,
            )
        if not session_id:
            return (
                no_update,
                "⚠ Session not initialised — refresh.",
                no_update,
                no_update,
                no_update,
            )

        theme = theme or "light"
        old = get_view(session_id)
        known_stations = None
        if use_known and old is not None:
            known_stations = old.data.stations

        tmpdir = None
        try:
            from pycsamt.map import MapView

            tmpdir = _decode_to_tempdir(filenames, contents)
            view = MapView.from_inversion_results(
                tmpdir,
                known_stations=known_stations,
                theme=theme,
            )
            if view.n_stations == 0:
                return (
                    no_update,
                    "⚠ No ModEM stations could be parsed.",
                    no_update,
                    no_update,
                    no_update,
                )

            if mode == _MODE_APPEND and old is not None:
                view = merge_views(old, view)
            set_view(session_id, view)
            store = store_from_view(view)
            n_s, n_l = store["n_stations"], store["n_lines"]
            feedback = f"✓ Imported {n_s} station(s) from {n_l} line(s)."
            badge = f"{n_s} stations · {n_l} line(s)"
            return store, feedback, False, badge, "mv-data-badge visible"
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            return (
                no_update,
                f"✗ Error: {exc}",
                no_update,
                no_update,
                no_update,
            )
        finally:
            if tmpdir:
                shutil.rmtree(tmpdir, ignore_errors=True)
