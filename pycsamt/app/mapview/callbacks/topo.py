# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Elevation-source callbacks: upload, fetch, and apply to the view."""

from __future__ import annotations

from dash import Input, Output, State, dcc, html, no_update

from .._ids import IDs
from .._render import store_from_view
from ..cache import get_view, set_view


def register_topo(app) -> None:
    _register_source_visibility(app)
    _register_parse_upload(app)
    _register_apply(app)
    _register_export(app)


def _register_source_visibility(app) -> None:
    app.clientside_callback(
        """
        function(src) {
            var s = src || 'stations';
            return [
                s === 'upload' ? {display:'block'} : {display:'none'},
                s === 'fetch'  ? {display:'block'} : {display:'none'}
            ];
        }
        """,
        Output(IDs.TOPO_UPLOAD_WRAP, "style"),
        Output(IDs.TOPO_FETCH_WRAP, "style"),
        Input(IDs.TOPO_SRC, "value"),
        prevent_initial_call=False,
    )


def _register_parse_upload(app) -> None:
    @app.callback(
        Output(IDs.TOPO_UPLOAD_STORE, "data"),
        Output(IDs.TOPO_STATUS, "children"),
        Input(IDs.TOPO_UPLOAD, "contents"),
        State(IDs.TOPO_UPLOAD, "filename"),
        prevent_initial_call=True,
    )
    def parse(contents, filename):
        if not contents:
            return no_update, no_update
        from pycsamt.map import parse_elevation_file

        elev = parse_elevation_file(contents, filename or "")
        if not elev:
            return {}, _msg(
                "No station/elevation columns found in file.", ok=False
            )
        return elev, _msg(f"Parsed {len(elev)} elevations from file.")


def _register_apply(app) -> None:
    @app.callback(
        Output(IDs.STORE_DATA, "data", allow_duplicate=True),
        Output(IDs.TOPO_STATUS, "children", allow_duplicate=True),
        Output(IDs.CTL_TOPO, "value"),
        Input(IDs.BTN_TOPO_APPLY, "n_clicks"),
        State(IDs.TOPO_SRC, "value"),
        State(IDs.TOPO_UPLOAD_STORE, "data"),
        State(IDs.TOPO_API, "value"),
        State(IDs.SESSION_ID, "data"),
        prevent_initial_call=True,
    )
    def apply(_n, src, upload_elev, api, session_id):
        view = get_view(session_id)
        if view is None:
            return no_update, _msg("Load lines first.", ok=False), no_update

        src = src or "stations"
        try:
            if src == "upload":
                elev_map = upload_elev or {}
                if not elev_map:
                    return (no_update,
                            _msg("Upload a file first.", ok=False),
                            no_update)
            elif src == "fetch":
                elev_map = view.fetch_elevations(api_name=api or "open_meteo")
                if not elev_map:
                    return (no_update,
                            _msg("Fetch returned no elevations "
                                 "(check coordinates / connection).",
                                 ok=False),
                            no_update)
            else:  # stations — already on the view, just confirm draping
                store = store_from_view(view)
                return store, _msg("Using station EDI elevations."), True

            new_view = view.with_elevations(elev_map)
            set_view(session_id, new_view)
            store = store_from_view(new_view)
            n = sum(
                1 for s in new_view.data.stations
                if str(s.id) in {str(k) for k in elev_map}
            )
            return store, _msg(f"Applied {n} elevations. Draping on."), True
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            return no_update, _msg(f"Error: {exc}", ok=False), no_update


def _register_export(app) -> None:
    @app.callback(
        Output(IDs.TOPO_EXPORT_DL, "data"),
        Output(IDs.TOPO_STATUS, "children", allow_duplicate=True),
        Input(IDs.BTN_TOPO_EXPORT, "n_clicks"),
        State(IDs.TOPO_EXPORT_FMT, "value"),
        State(IDs.SESSION_ID, "data"),
        prevent_initial_call=True,
    )
    def export(_n, fmt, session_id):
        import os
        import shutil
        import tempfile

        view = get_view(session_id)
        if view is None:
            return no_update, _msg("Load lines first.", ok=False)

        fmt = fmt or "csv"
        tmpdir = tempfile.mkdtemp(prefix="pycsamt_topo_export_")
        try:
            path = os.path.join(tmpdir, f"topography.{fmt}")
            view.export_topography(path, fmt=fmt)
            filename = f"topography.{fmt}"
            if fmt == "csv":
                with open(path, encoding="utf-8") as fh:
                    payload = dcc.send_string(fh.read(), filename)
            else:
                with open(path, "rb") as fh:
                    payload = dcc.send_bytes(fh.read(), filename)
            n = sum(1 for s in view.data.stations if s.elevation is not None)
            return payload, _msg(f"Exported {n} station elevation(s).")
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            return no_update, _msg(f"Error: {exc}", ok=False)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


def _msg(text: str, *, ok: bool = True):
    icon = "bi-check-circle" if ok else "bi-exclamation-triangle"
    color = "var(--mv-ok)" if ok else "var(--mv-err)"
    return html.Span(
        [html.I(className=f"bi {icon} me-1"), text],
        style={"color": color, "fontSize": "11px"},
    )
