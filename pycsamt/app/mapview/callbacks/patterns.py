# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Pattern packs (Geology Studio's "Patterns" tab) — browse/import a
swatch pack, crop a named swatch by hand out of a rasterized sheet, and
assign it to the selected Legend row's ``pattern_id``.

pycsamt never vendors a third-party pattern pack (e.g. the USGS FGDC
geologic-symbol set) -- only this import/persistence/crop machinery,
built on :mod:`pycsamt.geology.patterns`.
"""

from __future__ import annotations

from dash import ALL, Input, Output, State, ctx, html, no_update
from dash.exceptions import PreventUpdate

from pycsamt.app._patterns import (
    add_tile_from_shape,
    import_pack_path,
    import_pack_upload,
    list_packs,
    load_pack,
    sheet_figure,
    tile_grid_children,
)

from .._ids import IDs


def register_patterns(app) -> None:
    _register_pack_options(app)
    _register_import(app)
    _register_sheet_options(app)
    _register_sheet_figure(app)
    _register_save_crop(app)
    _register_swatch_grid(app)
    _register_swatch_select(app)
    _register_assign(app)


# ---------------------------------------------------------------------------
# pack selection
# ---------------------------------------------------------------------------


def _register_pack_options(app) -> None:
    @app.callback(
        Output(IDs.GEO_PACK_SELECT, "options"),
        Output(IDs.GEO_PACK_SELECT, "value"),
        Input(IDs.BTN_GEO_STUDIO, "n_clicks"),
        Input(IDs.GEO_PACKS_REFRESH_STORE, "data"),
        State(IDs.GEO_PACK_SELECT, "value"),
        prevent_initial_call=False,
    )
    def pack_options(_n, _refresh, current):
        packs = list_packs()
        options = [
            {"label": f"{p['name']} ({p['n_tiles']})", "value": p["id"]}
            for p in packs
        ]
        ids = {p["id"] for p in packs}
        value = current if current in ids else (packs[0]["id"] if packs else None)
        return options, value


# ---------------------------------------------------------------------------
# import (upload or local path)
# ---------------------------------------------------------------------------


def _register_import(app) -> None:
    @app.callback(
        Output(IDs.GEO_PACKS_REFRESH_STORE, "data", allow_duplicate=True),
        Output(IDs.GEO_PACK_SELECT, "value", allow_duplicate=True),
        Output(IDs.GEO_PACK_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_PACK_IMPORT_UPLOAD, "contents"),
        Input(IDs.GEO_BTN_IMPORT_PACK_PATH, "n_clicks"),
        State(IDs.GEO_PACK_IMPORT_UPLOAD, "filename"),
        State(IDs.GEO_PACK_PATH_INPUT, "value"),
        State(IDs.GEO_PACK_NAME_INPUT, "value"),
        State(IDs.GEO_PACKS_REFRESH_STORE, "data"),
        prevent_initial_call=True,
    )
    def do_import(contents, n_clicks, filename, path, name, refresh):
        trigger = ctx.triggered_id
        try:
            if trigger == IDs.GEO_PACK_IMPORT_UPLOAD and contents:
                summary = import_pack_upload(
                    contents, filename, name=name or None
                )
            elif trigger == IDs.GEO_BTN_IMPORT_PACK_PATH and n_clicks:
                summary = import_pack_path(path, name=name or None)
            else:
                raise PreventUpdate
        except PreventUpdate:
            raise
        except Exception as error:  # noqa: BLE001
            return no_update, no_update, _err(str(error))
        return (
            (refresh or 0) + 1,
            summary["id"],
            _ok(
                f"Imported {summary['name']} — "
                f"{summary['n_sheets']} sheet(s)"
            ),
        )


# ---------------------------------------------------------------------------
# sheet select + crop-ready figure
# ---------------------------------------------------------------------------


def _register_sheet_options(app) -> None:
    @app.callback(
        Output(IDs.GEO_PACK_SHEET_SELECT, "options"),
        Output(IDs.GEO_PACK_SHEET_SELECT, "value"),
        Input(IDs.GEO_PACK_SELECT, "value"),
        Input(IDs.GEO_PACKS_REFRESH_STORE, "data"),
        prevent_initial_call=False,
    )
    def sheet_options(pack_id, _refresh):
        if not pack_id:
            return [], None
        try:
            pack = load_pack(pack_id)
        except Exception:  # noqa: BLE001
            return [], None
        options = [
            {"label": f"{s.source_file} p.{s.page}", "value": s.id}
            for s in pack.sheets
        ]
        return options, (options[0]["value"] if options else None)


def _register_sheet_figure(app) -> None:
    @app.callback(
        Output(IDs.GEO_PACK_SHEET_GRAPH, "figure"),
        Input(IDs.GEO_PACK_SHEET_SELECT, "value"),
        State(IDs.GEO_PACK_SELECT, "value"),
        prevent_initial_call=False,
    )
    def figure(sheet_id, pack_id):
        import plotly.graph_objects as go

        if not pack_id or not sheet_id:
            return go.Figure()
        try:
            return sheet_figure(pack_id, sheet_id)
        except Exception:  # noqa: BLE001
            return go.Figure()


# ---------------------------------------------------------------------------
# save a crop as a named swatch
# ---------------------------------------------------------------------------


def _register_save_crop(app) -> None:
    @app.callback(
        Output(IDs.GEO_PACKS_REFRESH_STORE, "data", allow_duplicate=True),
        Output(IDs.GEO_PACK_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_BTN_SAVE_CROP, "n_clicks"),
        State(IDs.GEO_PACK_SHEET_GRAPH, "relayoutData"),
        State(IDs.GEO_PACK_SELECT, "value"),
        State(IDs.GEO_PACK_SHEET_SELECT, "value"),
        State(IDs.GEO_PACK_CROP_NAME, "value"),
        State(IDs.GEO_PACKS_REFRESH_STORE, "data"),
        prevent_initial_call=True,
    )
    def save_crop(n_clicks, relayout, pack_id, sheet_id, name, refresh):
        if not n_clicks:
            raise PreventUpdate
        if not pack_id or not sheet_id:
            return no_update, _err("pick a pack and a sheet first")
        shapes = (relayout or {}).get("shapes")
        if not shapes:
            return no_update, _err(
                "draw a rectangle on the sheet first (top-right toolbar)"
            )
        shape = shapes[-1]
        try:
            add_tile_from_shape(pack_id, sheet_id, shape, name or "")
        except Exception as error:  # noqa: BLE001
            return no_update, _err(str(error))
        return (refresh or 0) + 1, _ok(f"Saved swatch {name!r}")


# ---------------------------------------------------------------------------
# swatch grid + selection
# ---------------------------------------------------------------------------


def _register_swatch_grid(app) -> None:
    @app.callback(
        Output(IDs.GEO_PACK_SWATCH_GRID, "children"),
        Input(IDs.GEO_PACK_SELECT, "value"),
        Input(IDs.GEO_PACKS_REFRESH_STORE, "data"),
        Input(IDs.GEO_PACK_SELECTED_TILE_STORE, "data"),
        prevent_initial_call=False,
    )
    def grid(pack_id, _refresh, selected):
        if not pack_id:
            return []
        selected_tile = (
            (selected or {}).get("tile")
            if (selected or {}).get("pack") == pack_id
            else None
        )
        try:
            return tile_grid_children(pack_id, selected_tile_id=selected_tile)
        except Exception:  # noqa: BLE001
            return []


def _register_swatch_select(app) -> None:
    @app.callback(
        Output(IDs.GEO_PACK_SELECTED_TILE_STORE, "data"),
        Input({"type": "geo-swatch-btn", "pack": ALL, "tile": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def select(clicks):
        triggered = ctx.triggered_id
        if not triggered or not any(c for c in clicks if c):
            raise PreventUpdate
        return {"pack": triggered["pack"], "tile": triggered["tile"]}


# ---------------------------------------------------------------------------
# assign the selected swatch to the selected Legend row
# ---------------------------------------------------------------------------


def _register_assign(app) -> None:
    @app.callback(
        Output(IDs.GEO_TABLE_LEGEND, "data", allow_duplicate=True),
        Output(IDs.GEO_PACK_STATUS, "children", allow_duplicate=True),
        Input(IDs.GEO_BTN_ASSIGN_PATTERN, "n_clicks"),
        State(IDs.GEO_TABLE_LEGEND, "data"),
        State(IDs.GEO_TABLE_LEGEND, "selected_rows"),
        State(IDs.GEO_PACK_SELECTED_TILE_STORE, "data"),
        prevent_initial_call=True,
    )
    def assign(n_clicks, rows, selected_rows, tile_sel):
        if not n_clicks:
            raise PreventUpdate
        if not selected_rows:
            return no_update, _err(
                "select a row in the Legend tab's table first"
            )
        if not tile_sel or not tile_sel.get("tile"):
            return no_update, _err("select a swatch first")
        rows = list(rows or [])
        idx = selected_rows[0]
        if idx >= len(rows):
            return no_update, _err("selected row no longer exists")
        rows[idx] = {
            **rows[idx],
            "pattern_id": tile_sel["tile"],
            "pattern_source": tile_sel.get("pack") or "",
        }
        return rows, _ok(f"Assigned pattern to legend row {idx + 1}")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _ok(text: str):
    return html.Span(text, className="text-success")


def _err(text: str):
    return html.Span(text, className="text-danger")
