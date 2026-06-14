# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/data.py — modal toggle, EDI loading, station table update.

Two load triggers are handled in a single callback:
  • dcc.Upload (drag-and-drop / file picker) → decodes base64, writes tmpdir
  • BTN_LOAD_CONFIRM (server-path text input) → scans filesystem

Note: load_data is intentionally synchronous so it works even when diskcache
is not installed.  Install diskcache for background agent execution:
    pip install diskcache
"""

from __future__ import annotations

from dash import ctx, Input, Output, State, no_update

from pycsamt.app.desktop.controllers.data_controller import DataController
from pycsamt.app.web.cache import cache_set
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import decode_upload_to_tempdir, find_edi_files_by_line


def register_data(app) -> None:
    _register_modal_toggle(app)
    _register_load_data(app)
    _register_table_update(app)


# ──────────────────────────────────────────────────────────────────────────────

def _register_modal_toggle(app) -> None:
    @app.callback(
        Output(IDs.MODAL_LOAD,          "is_open"),
        Input(IDs.BTN_LOAD,             "n_clicks"),
        Input("modal-close-btn",        "n_clicks"),
        Input(IDs.WELCOME_CTA,          "n_clicks"),
        Input(IDs.WELCOME_CARD_LOAD,    "n_clicks"),
        Input(IDs.WELCOME_CARD_PIPE,    "n_clicks"),
        Input(IDs.WELCOME_CARD_AGENTS,  "n_clicks"),
        Input(IDs.WELCOME_CARD_VIZ,     "n_clicks"),
        State(IDs.MODAL_LOAD,           "is_open"),
        prevent_initial_call=True,
    )
    def toggle_modal(_nb, _close, _cta, _cl, _cp, _ca, _cv, is_open):
        if ctx.triggered_id == "modal-close-btn":
            return False
        return True


def _register_load_data(app) -> None:
    @app.callback(
        Output(IDs.STORE_DATA,      "data"),
        Output(IDs.LOAD_FEEDBACK,   "children"),
        Output(IDs.MODAL_LOAD,      "is_open",  allow_duplicate=True),
        Output(IDs.UPLOAD_FILELIST, "children"),
        Input(IDs.UPLOAD,           "contents"),
        Input(IDs.BTN_LOAD_CONFIRM, "n_clicks"),
        State(IDs.UPLOAD,           "filename"),
        State(IDs.INPUT_PATH,       "value"),
        State(IDs.SESSION_ID,       "data"),
        prevent_initial_call=True,
    )
    def load_data(upload_contents, _n_clicks, upload_filenames, path, session_id):
        _SKIP = (no_update, no_update, no_update, no_update)

        if not session_id:
            return no_update, "⚠ Session not initialised — refresh the page.", no_update, no_update

        trigger = ctx.triggered_id

        # ── Upload branch ─────────────────────────────────────────────────
        if trigger == IDs.UPLOAD:
            if not upload_contents:
                return _SKIP

            import shutil

            paths, tmpdir = decode_upload_to_tempdir(upload_contents, upload_filenames)
            if not paths:
                return (
                    no_update,
                    "⚠ No recognised files in upload (.edi / .avg / .j).",
                    no_update,
                    no_update,
                )

            names   = upload_filenames if isinstance(upload_filenames, list) else [upload_filenames]
            preview = ", ".join(names[:5]) + (f" … and {len(names) - 5} more" if len(names) > 5 else "")

            try:
                ctrl  = DataController()
                sites = ctrl.load(paths)
                cache_set(session_id, sites)
                df    = ctrl.dataframe
                store = {
                    "station_records": df.to_dict("records"),
                    "data_dir":        "[uploaded]",
                    "n_stations":      len(df),
                    "n_lines":         1,
                    "line_counts":     {"uploaded": len(df)},
                }
                return (
                    store,
                    f"✓ Loaded {len(df)} station(s) from {len(paths)} file(s).",
                    False,
                    f"✓ {len(paths)} file(s): {preview}",
                )
            except Exception as exc:
                return no_update, f"✗ Error loading files: {exc}", no_update, ""
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)

        # ── Server-path branch ────────────────────────────────────────────
        if trigger == IDs.BTN_LOAD_CONFIRM:
            if not path or not path.strip():
                return no_update, "⚠ Enter a path first.", no_update, no_update

            paths, line_map = find_edi_files_by_line(path.strip())
            if not paths:
                return no_update, f"⚠ No EDI files found in: {path}", no_update, no_update

            path_to_line = {
                p: line
                for line, plist in line_map.items()
                for p in plist
            }
            n_lines = len(line_map)

            try:
                ctrl  = DataController()
                sites = ctrl.load(paths, path_to_line=path_to_line)
                cache_set(session_id, sites)
                df    = ctrl.dataframe

                # Build a compact per-line summary (max 5 lines shown)
                line_counts = {line: len(plist) for line, plist in line_map.items()}
                shown = list(line_counts.items())[:5]
                line_summary = ", ".join(f"{ln}: {cnt}" for ln, cnt in shown)
                if n_lines > 5:
                    line_summary += f" … +{n_lines - 5} more"

                store = {
                    "station_records": df.to_dict("records"),
                    "data_dir":        path.strip(),
                    "n_stations":      len(df),
                    "n_lines":         n_lines,
                    "line_counts":     line_counts,
                }
                n_s = len(df)
                feedback = (
                    f"✓ Loaded {n_s} station{'s' if n_s != 1 else ''} "
                    f"from {n_lines} line{'s' if n_lines != 1 else ''}"
                    + (f" ({line_summary})" if n_lines > 1 else "")
                    + "."
                )
                return store, feedback, False, no_update
            except Exception as exc:
                return no_update, f"✗ Error: {exc}", no_update, no_update

        return _SKIP


def _register_table_update(app) -> None:
    @app.callback(
        Output(IDs.STATION_TABLE,   "data"),
        Output(IDs.STATION_SUMMARY, "children"),
        Output(IDs.STATUS_LABEL,    "children"),
        Output(IDs.KPI_STATIONS,    "children"),
        Output(IDs.KPI_FREQ,        "children"),
        Output(IDs.KPI_PROFILES,    "children"),
        Output(IDs.KPI_SURVEY,      "children"),
        Input(IDs.STORE_DATA,       "data"),
    )
    def update_table(store_data):
        if not store_data:
            return [], "No stations loaded", "No data loaded", "—", "—", "—", "—"

        import os
        import statistics

        records = store_data.get("station_records", [])
        n       = store_data.get("n_stations", 0)
        label   = f"{n} station{'s' if n != 1 else ''} loaded"

        kpi_stations = str(n) if n else "—"

        freq_vals = [r.get("N_freq", 0) for r in records if r.get("N_freq")]
        kpi_freq  = str(int(statistics.median(freq_vals))) if freq_vals else "—"

        # Use stored n_lines when available; fall back to distinct lat clusters
        n_lines = store_data.get("n_lines", 0)
        if n_lines:
            kpi_profiles = str(n_lines)
        else:
            lats = sorted({round(r.get("Latitude", 0), 3) for r in records if r.get("Latitude")})
            kpi_profiles = str(len(lats)) if lats else "—"

        data_dir = store_data.get("data_dir", "")
        if data_dir == "[uploaded]":
            kpi_survey = "upload"
        elif data_dir:
            kpi_survey = os.path.basename(data_dir.rstrip("/\\")) or data_dir[:10]
        else:
            kpi_survey = "—"

        return records, label, f"✓ {label}", kpi_stations, kpi_freq, kpi_profiles, kpi_survey
