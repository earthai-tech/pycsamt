# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/dashboard.py — Home-page survey dashboard.

Populates:
  - Per-line summary table
  - Data-quality indicators
  - Launch buttons (navigate to Map / Profile pages)
  - Profile page station dropdown
  - Map/Profile page "back to dashboard" navigation
"""
from __future__ import annotations

from dash import (
    Input, Output, State, clientside_callback,
    callback_context as ctx, no_update, html,
)
import dash_bootstrap_components as dbc

from pycsamt.app.web.layout import IDs


def register_dashboard(app) -> None:

    # ── 1. Populate dashboard when store data changes ─────────────────────────
    @app.callback(
        Output("dash-empty-state",   "style"),
        Output("dash-data-state",    "style"),
        Output(IDs.DASH_QUALITY,     "style"),
        Output("dash-launch-row",    "style"),
        Output(IDs.DASH_LINE_TABLE,  "children"),
        Output(IDs.DASH_QUALITY,     "children"),
        Output("dash-survey-name",   "children"),
        Output("dash-health-badge",  "children"),
        Output("dash-health-badge",  "className"),
        Input(IDs.STORE_DATA, "data"),
    )
    def update_dashboard(store_data):
        _hide = {"display": "none"}
        _show = {"display": "block"}

        if not store_data or not store_data.get("n_stations"):
            return _show, _hide, _hide, _hide, [], [], "", "", "dash-health-badge"

        n_sta   = store_data.get("n_stations", 0)
        n_lines = store_data.get("n_lines", 1)
        data_dir = store_data.get("data_dir", "")
        line_counts: dict = store_data.get("line_counts", {})
        records = store_data.get("station_records", [])

        # ── Per-line table ────────────────────────────────────────────────
        if line_counts and n_lines > 1:
            header = html.Div(
                [
                    html.Span("Line",      className="dlt-col dlt-head"),
                    html.Span("Stations",  className="dlt-col dlt-head dlt-num"),
                    html.Span("Coverage",  className="dlt-col dlt-head dlt-bar"),
                ],
                className="dlt-row dlt-header",
            )
            rows = [header]
            for line_name, cnt in line_counts.items():
                pct = int(cnt / n_sta * 100)
                rows.append(html.Div(
                    [
                        html.Span(str(line_name), className="dlt-col dlt-name",
                                  title=str(line_name)),
                        html.Span(str(cnt), className="dlt-col dlt-num"),
                        html.Div(
                            html.Div(className="dlt-fill",
                                     style={"width": f"{pct}%"}),
                            className="dlt-bar-track",
                        ),
                    ],
                    className="dlt-row",
                ))
            line_table_content = html.Div(rows, className="dlt-table")
        else:
            line_table_content = html.Span(
                f"{n_sta} station(s) loaded from a single folder.",
                className="text-muted small",
            )

        # ── Data quality indicators ───────────────────────────────────────
        n_tipper = sum(1 for r in records if r.get("Tipper"))
        n_valid  = sum(
            1 for r in records
            if r.get("Latitude") and r.get("Longitude")
        )
        n_nfreq  = [r.get("N_freq", 0) for r in records if r.get("N_freq")]
        avg_freq = int(sum(n_nfreq) / len(n_nfreq)) if n_nfreq else 0

        def _qc_chip(icon: str, text: str, ok: bool) -> html.Div:
            col = "text-success" if ok else "text-warning"
            return html.Div(
                [html.I(className=f"bi bi-{'check-circle' if ok else 'exclamation-triangle'} me-1 {col}"),
                 html.Span(text, style={"fontSize": "11px"})],
                className="dash-qc-chip",
            )

        quality_content = [
            _qc_chip("geo", f"{n_valid}/{n_sta} with coordinates", n_valid == n_sta),
            _qc_chip("freq", f"avg {avg_freq} frequencies/station", avg_freq >= 30),
            _qc_chip("tip", f"{n_tipper} stations with tipper", n_tipper > 0),
        ]

        # ── Survey name badge ─────────────────────────────────────────────
        import os
        dir_name = os.path.basename(data_dir.rstrip("/\\")) if data_dir else "survey"
        survey_badge = dir_name or "survey"

        # ── Health badge ──────────────────────────────────────────────────
        all_ok = (n_valid == n_sta) and (avg_freq >= 30)
        badge_text = "✓ All checks passed" if all_ok else "⚠ Review data quality"
        badge_cls  = "dash-health-badge good" if all_ok else "dash-health-badge warn"

        return (
            _hide, _show, _show, _show,
            line_table_content,
            quality_content,
            survey_badge,
            badge_text, badge_cls,
        )

    # ── 2. Profile page station dropdown — populate options ───────────────────
    @app.callback(
        Output(IDs.PROFILE_PAGE_STATION, "options"),
        Output(IDs.PROFILE_PAGE_STATION, "value"),
        Input(IDs.STORE_DATA, "data"),
    )
    def populate_station_dropdown(store_data):
        if not store_data:
            return [], None
        records = store_data.get("station_records", [])
        opts = [{"label": r["ID"], "value": r["ID"]} for r in records]
        first = opts[0]["value"] if opts else None
        return opts, first

    # ── 3. Profile page station dropdown → write STORE_SELECTION ─────────────
    @app.callback(
        Output(IDs.STORE_SELECTION, "data", allow_duplicate=True),
        Input(IDs.PROFILE_PAGE_STATION, "value"),
        prevent_initial_call=True,
    )
    def profile_station_selected(station_id):
        if not station_id:
            return no_update
        return {"station_id": str(station_id)}

    # ── 4. STORE_SELECTION → sync profile page dropdown ──────────────────────
    @app.callback(
        Output(IDs.PROFILE_PAGE_STATION, "value", allow_duplicate=True),
        Input(IDs.STORE_SELECTION, "data"),
        State(IDs.PROFILE_PAGE_STATION, "options"),
        prevent_initial_call=True,
    )
    def sync_profile_dropdown(selection, opts):
        if not selection or not opts:
            return no_update
        sid = selection.get("station_id", "")
        valid = {o["value"] for o in opts}
        return sid if sid in valid else no_update

    # ── 5. Map page: update station info card on map click ────────────────────
    @app.callback(
        Output(IDs.MAP_PAGE_INFO, "children"),
        Input(IDs.STORE_SELECTION, "data"),
        State(IDs.STORE_DATA,     "data"),
    )
    def update_map_station_info(selection, store_data):
        if not selection:
            return html.Span("Click a station on the map", className="text-muted small")
        sid = selection.get("station_id", "")
        if not sid or not store_data:
            return html.Span("Click a station on the map", className="text-muted small")
        records = store_data.get("station_records", [])
        match = next((r for r in records if str(r.get("ID")) == str(sid)), None)
        if not match:
            return html.Span(sid, className="fw-semibold small")
        lat  = match.get("Latitude", "—")
        lon  = match.get("Longitude", "—")
        elev = match.get("Elevation", "—")
        nf   = match.get("N_freq", "—")
        tip  = "Yes" if match.get("Tipper") else "No"
        line = match.get("Line", "—")
        return html.Div([
            html.Span(sid, className="map-info-name"),
            html.Div([
                html.Span("Line: ",  className="map-info-key"), html.Span(str(line)),
            ], className="map-info-row"),
            html.Div([
                html.Span("Lat: ",   className="map-info-key"),
                html.Span(f"{lat:.4f}°" if isinstance(lat, float) else str(lat)),
                html.Span("  Lon: ", className="map-info-key"),
                html.Span(f"{lon:.4f}°" if isinstance(lon, float) else str(lon)),
            ], className="map-info-row"),
            html.Div([
                html.Span("Elev: ",  className="map-info-key"),
                html.Span(f"{elev:.1f} m" if isinstance(elev, float) else str(elev)),
                html.Span("  N freq: ", className="map-info-key"),
                html.Span(str(nf)),
            ], className="map-info-row"),
            html.Div([
                html.Span("Tipper: ", className="map-info-key"),
                html.Span(tip),
            ], className="map-info-row"),
        ])

    # ── 6. Map page: populate survey-line filter checklist ───────────────────
    @app.callback(
        Output(IDs.MAP_PAGE_LINE_FILTER, "options"),
        Output(IDs.MAP_PAGE_LINE_FILTER, "value"),
        Output("map-line-filter-count",  "children"),
        Input(IDs.STORE_DATA, "data"),
    )
    def populate_map_line_filter(store_data):
        if not store_data:
            return [], [], "(all)"
        line_counts: dict = store_data.get("line_counts", {})
        if not line_counts:
            return [], [], "(all)"
        opts = [{"label": f"{ln} ({cnt})", "value": ln}
                for ln, cnt in line_counts.items()]
        values = [o["value"] for o in opts]
        count_label = f"({len(values)})" if values else "(all)"
        return opts, values, count_label

    # ── 7. Map page: stats strip ──────────────────────────────────────────────
    @app.callback(
        Output(IDs.MAP_PAGE_STATS, "children"),
        Input(IDs.STORE_DATA,            "data"),
        Input(IDs.MAP_PAGE_LINE_FILTER,  "value"),
    )
    def update_map_stats(store_data, line_filter):
        if not store_data:
            return ""
        records = store_data.get("station_records", [])
        if line_filter:
            records = [r for r in records if r.get("Line") in line_filter]
        n_sta   = len(records)
        n_lines = len({r.get("Line") for r in records if r.get("Line")})
        n_tip   = sum(1 for r in records if r.get("Tipper"))
        n_geo   = sum(1 for r in records if r.get("Latitude") and r.get("Longitude"))

        def _stat(label, value):
            return html.Div(
                [html.Span(str(value), className="map-stat-value"),
                 html.Span(label,      className="map-stat-label")],
                className="map-stat-item",
            )

        return [
            _stat("stations",          n_sta),
            _stat("lines",             n_lines),
            _stat("with tipper",       n_tip),
            _stat("with coordinates",  n_geo),
        ]

    # ── 9. Navigation: dashboard launch buttons → switch page ─────────────────
    clientside_callback(
        "function(n) { if (n) return 'map'; return window.dash_clientside.no_update; }",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input(IDs.DASH_LAUNCH_MAP, "n_clicks"),
        prevent_initial_call=True,
    )
    clientside_callback(
        "function(n) { if (n) return 'profile'; return window.dash_clientside.no_update; }",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input(IDs.DASH_LAUNCH_PROFILE, "n_clicks"),
        prevent_initial_call=True,
    )
    clientside_callback(
        "function(n) { if (n) return 'home'; return window.dash_clientside.no_update; }",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input("map-page-home-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    clientside_callback(
        "function(n) { if (n) return 'home'; return window.dash_clientside.no_update; }",
        Output(IDs.NAV_SECTION, "data", allow_duplicate=True),
        Input("profile-page-home-btn", "n_clicks"),
        prevent_initial_call=True,
    )
