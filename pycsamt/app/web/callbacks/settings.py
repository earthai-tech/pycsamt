# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/settings.py — Settings off-canvas panel.

Responsibilities
----------------
- Open/close the Settings offcanvas (BTN_SETTINGS click).
- Persist settings in SETTINGS_STORE (localStorage).
- Save profile as downloadable JSON (BTN_SETTINGS_SAVE).
- Load a profile JSON uploaded via SETTINGS_UPLOAD.
- Reset all settings to defaults (BTN_SETTINGS_RESET).
"""
from __future__ import annotations

import base64
import json

from dash import (
    Input, Output, State,
    callback_context as ctx,
    no_update,
    dcc,
)

from pycsamt.app.web.layout import IDs

_DEFAULTS = {
    "interp":     "bilinear",
    "skindepth":  "cagniard",
    "dem":        "open-elevation",
    "cmap":       "viridis",
}


def register_settings(app) -> None:

    # 1. Toggle offcanvas open / close
    @app.callback(
        Output(IDs.SETTINGS_CANVAS, "is_open"),
        Input(IDs.BTN_SETTINGS, "n_clicks"),
        State(IDs.SETTINGS_CANVAS, "is_open"),
        prevent_initial_call=True,
    )
    def _toggle(n, is_open):
        if n:
            return not is_open
        return no_update

    # 2. Populate form fields from store on open
    @app.callback(
        Output("settings-interp",    "value"),
        Output("settings-skindepth", "value"),
        Output("settings-dem",       "value"),
        Output("settings-cmap",      "value"),
        Input(IDs.SETTINGS_CANVAS,   "is_open"),
        State(IDs.SETTINGS_STORE,    "data"),
        prevent_initial_call=True,
    )
    def _populate(is_open, stored):
        if not is_open:
            return no_update, no_update, no_update, no_update
        cfg = {**_DEFAULTS, **(stored or {})}
        return cfg["interp"], cfg["skindepth"], cfg["dem"], cfg["cmap"]

    # 3. Save form → store on any value change
    @app.callback(
        Output(IDs.SETTINGS_STORE, "data"),
        Input("settings-interp",    "value"),
        Input("settings-skindepth", "value"),
        Input("settings-dem",       "value"),
        Input("settings-cmap",      "value"),
        prevent_initial_call=True,
    )
    def _persist(interp, skindepth, dem, cmap):
        return {
            "interp":    interp    or _DEFAULTS["interp"],
            "skindepth": skindepth or _DEFAULTS["skindepth"],
            "dem":       dem       or _DEFAULTS["dem"],
            "cmap":      cmap      or _DEFAULTS["cmap"],
        }

    # 4. Reset to defaults
    @app.callback(
        Output("settings-interp",    "value", allow_duplicate=True),
        Output("settings-skindepth", "value", allow_duplicate=True),
        Output("settings-dem",       "value", allow_duplicate=True),
        Output("settings-cmap",      "value", allow_duplicate=True),
        Output(IDs.SETTINGS_STORE,   "data",  allow_duplicate=True),
        Output(IDs.SETTINGS_FEEDBACK,"children"),
        Input(IDs.BTN_SETTINGS_RESET,"n_clicks"),
        prevent_initial_call=True,
    )
    def _reset(n):
        if not n:
            return no_update, no_update, no_update, no_update, no_update, no_update
        d = _DEFAULTS
        return (
            d["interp"], d["skindepth"], d["dem"], d["cmap"],
            d,
            "✓ All settings reset to defaults.",
        )

    # 5. Save profile — trigger a browser download via dcc.Download
    app.layout.children  # ensure layout is ready (no-op, just safe guard)

    @app.callback(
        Output("settings-download", "data"),
        Output(IDs.SETTINGS_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_SETTINGS_SAVE, "n_clicks"),
        State(IDs.SETTINGS_STORE,    "data"),
        prevent_initial_call=True,
    )
    def _save(n, stored):
        if not n:
            return no_update, no_update
        cfg = {**_DEFAULTS, **(stored or {})}
        content = json.dumps(cfg, indent=2)
        return (
            dict(content=content, filename="pycsamt_settings.json"),
            "✓ Profile saved.",
        )

    # 6. Load profile from uploaded JSON
    @app.callback(
        Output("settings-interp",    "value", allow_duplicate=True),
        Output("settings-skindepth", "value", allow_duplicate=True),
        Output("settings-dem",       "value", allow_duplicate=True),
        Output("settings-cmap",      "value", allow_duplicate=True),
        Output(IDs.SETTINGS_STORE,   "data",  allow_duplicate=True),
        Output(IDs.SETTINGS_FEEDBACK,"children", allow_duplicate=True),
        Input(IDs.SETTINGS_UPLOAD,   "contents"),
        State(IDs.SETTINGS_UPLOAD,   "filename"),
        prevent_initial_call=True,
    )
    def _load(contents, filename):
        if not contents:
            return no_update, no_update, no_update, no_update, no_update, no_update
        try:
            _header, encoded = contents.split(",", 1)
            raw = base64.b64decode(encoded).decode("utf-8")
            cfg = {**_DEFAULTS, **json.loads(raw)}
            return (
                cfg["interp"], cfg["skindepth"], cfg["dem"], cfg["cmap"],
                cfg,
                f"✓ Profile loaded from {filename}.",
            )
        except Exception as exc:
            return (
                no_update, no_update, no_update, no_update, no_update,
                f"✗ Could not parse {filename}: {exc}",
            )
