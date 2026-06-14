# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""callbacks/map.py — station map figure and map↔table selection sync."""

from __future__ import annotations

import pandas as pd
from dash import ctx, Input, Output, State, no_update

from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import build_station_map


def register_map(app) -> None:
    _register_map_update(app)
    _register_selection_sync(app)
    _register_status_station(app)


# ──────────────────────────────────────────────────────────────────────────────

def _register_map_update(app) -> None:
    @app.callback(
        Output(IDs.MAP_GRAPH,            "figure"),
        Input(IDs.STORE_DATA,            "data"),
        Input(IDs.MAP_OVERLAY,           "value"),
        Input(IDs.STORE_SELECTION,       "data"),
        Input(IDs.STORE_THEME,           "data"),
        Input(IDs.MAP_PAGE_LINE_FILTER,  "value"),
        Input(IDs.MAP_PAGE_MARKER_SIZE,  "value"),
        Input(IDs.MAP_PAGE_BASEMAP,      "value"),
    )
    def update_map(store_data, overlay, selection, theme,
                   line_filter, marker_size, basemap_style):
        dark = (theme or "dark") == "dark"
        if not store_data:
            return build_station_map(pd.DataFrame(), dark=dark)

        df = pd.DataFrame(store_data.get("station_records", []))
        if df.empty:
            return build_station_map(df, dark=dark)

        selected_id = (selection or {}).get("station_id")
        return build_station_map(
            df,
            selected_id=selected_id,
            overlay=overlay or "Index",
            dark=dark,
            line_filter=line_filter or None,
            marker_size=int(marker_size) if marker_size else 10,
            basemap_style=basemap_style or None,
        )


def _register_selection_sync(app) -> None:
    @app.callback(
        Output(IDs.STORE_SELECTION,  "data"),
        Input(IDs.MAP_GRAPH,         "clickData"),
        Input(IDs.STATION_TABLE,     "selected_rows"),
        State(IDs.STORE_DATA,        "data"),
        prevent_initial_call=True,
    )
    def sync_selection(click_data, selected_rows, store_data):
        trigger = ctx.triggered_id

        if trigger == IDs.MAP_GRAPH and click_data:
            points = click_data.get("points", [])
            if points:
                station_id = points[0].get("customdata") or points[0].get("text")
                return {"station_id": str(station_id)}

        if trigger == IDs.STATION_TABLE and selected_rows and store_data:
            records = store_data.get("station_records", [])
            if selected_rows[0] < len(records):
                return {"station_id": records[selected_rows[0]]["ID"]}

        return no_update


def _register_status_station(app) -> None:
    @app.callback(
        Output(IDs.STATUS_STATION, "children"),
        Input(IDs.STORE_SELECTION, "data"),
        State(IDs.STORE_DATA,      "data"),
    )
    def update_status_station(selection, store_data):
        if not selection:
            return "No station selected"
        station_id = selection.get("station_id", "")
        if not station_id:
            return "No station selected"

        label = f"Station: {station_id}"
        if store_data:
            records = store_data.get("station_records", [])
            match = next((r for r in records if str(r.get("ID")) == str(station_id)), None)
            if match:
                lat = match.get("Latitude", "")
                lon = match.get("Longitude", "")
                if lat and lon:
                    label = f"Station: {station_id}  ({lat:.4f}°N, {lon:.4f}°E)"
        return label
