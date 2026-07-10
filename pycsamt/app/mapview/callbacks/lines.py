# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Line pills (visibility), station list, selection, inspector, table."""

from __future__ import annotations

from dash import (
    ALL,
    Input,
    Output,
    State,
    ctx,
    html,
    no_update,
)

from .._ids import IDs

_PALETTE = [
    "#1e66f5", "#40a02b", "#fe640b", "#8839ef", "#179299",
    "#ea76cb", "#df8e1d", "#d20f39", "#04a5e5", "#7287fd",
]


def register_lines(app) -> None:
    _register_init_lines(app)
    _register_render_panel(app)
    _register_pill_toggle(app)
    _register_selection(app)
    _register_inspector(app)
    _register_table(app)


def _colors_for(lines):
    return {ln: _PALETTE[i % len(_PALETTE)] for i, ln in enumerate(lines)}


def _register_init_lines(app) -> None:
    @app.callback(
        Output(IDs.STORE_LINES, "data"),
        Input(IDs.STORE_DATA, "data"),
        State(IDs.STORE_LINES, "data"),
        prevent_initial_call=True,
    )
    def init_lines(store, prev):
        if not store:
            return {}
        all_lines = list(store.get("lines", []))
        prev_active = set((prev or {}).get("active", []))
        active = [ln for ln in all_lines if ln in prev_active] or all_lines
        return {
            "all": all_lines,
            "active": list(active),
            "colors": _colors_for(all_lines),
        }


def _register_render_panel(app) -> None:
    @app.callback(
        Output(IDs.LINE_PILLS, "children"),
        Output(IDs.STATION_LIST, "children"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_LINES, "data"),
        Input(IDs.STORE_SELECTION, "data"),
        Input(IDs.STORE_MASKED, "data"),
        prevent_initial_call=False,
    )
    def render(store, lines, selection, masked):
        if not store or not store.get("station_records"):
            return [], [html.Div("No lines loaded.", className="mv-empty")]

        all_lines = (lines or {}).get("all", [])
        active = set((lines or {}).get("active", all_lines))
        colors = (lines or {}).get("colors") or _colors_for(all_lines)
        selected = (selection or {}).get("station_id")
        masked_set = {str(m) for m in (masked or [])}

        pills = [
            html.Span(
                "All",
                id={"type": "mv-pill", "index": "__all__"},
                className="mv-pill mv-pill-all",
                n_clicks=0,
            )
        ]
        for ln in all_lines:
            col = colors.get(ln, _PALETTE[0])
            is_active = ln in active
            pills.append(html.Span(
                ["● ", ln],
                id={"type": "mv-pill", "index": ln},
                className="mv-pill" + ("" if is_active else " mv-pill-off"),
                n_clicks=0,
                style={"color": col, "borderColor": col},
            ))

        rows = []
        for rec in store.get("station_records", []):
            line = rec.get("Line", "") or "line"
            if line not in active:
                continue
            sid = str(rec.get("ID", ""))
            col = colors.get(line, "#6c7086")
            is_sel = sid == selected
            is_masked = sid in masked_set
            meta = _row_meta(rec)
            cls = "mv-sta-row"
            if is_sel:
                cls += " mv-sta-sel"
            if is_masked:
                cls += " mv-sta-masked"
            rows.append(html.Div(
                [
                    html.Div(
                        [
                            html.Span(line, className="mv-sta-badge",
                                      style={"color": col}),
                            html.Span(sid, className="mv-sta-id"),
                            html.Button(
                                html.I(className=(
                                    "bi bi-eye-slash" if is_masked
                                    else "bi bi-eye"
                                )),
                                id={"type": "mv-mask", "index": sid},
                                className="mv-mask-btn",
                                n_clicks=0,
                                title=("Unmask station" if is_masked
                                       else "Mask (exclude) station"),
                            ),
                        ],
                        className="mv-sta-top",
                    ),
                    html.Div(meta, className="mv-sta-meta"),
                ],
                id={"type": "mv-sta-row", "index": sid},
                className=cls,
                n_clicks=0,
                style={"borderLeftColor": col},
            ))
        if not rows:
            rows = [html.Div("No active lines.", className="mv-empty")]
        return pills, rows


def _row_meta(rec) -> str:
    parts = []
    try:
        if rec.get("Latitude") is not None:
            parts.append(f"{float(rec['Latitude']):+.4f}°")
        if rec.get("Longitude") is not None:
            parts.append(f"{float(rec['Longitude']):+.4f}°")
        if rec.get("Elevation") is not None:
            parts.append(f"{int(rec['Elevation']):,} m")
    except (TypeError, ValueError):
        pass
    return " · ".join(parts) or "—"


def _register_pill_toggle(app) -> None:
    @app.callback(
        Output(IDs.STORE_LINES, "data", allow_duplicate=True),
        Input({"type": "mv-pill", "index": ALL}, "n_clicks"),
        State(IDs.STORE_LINES, "data"),
        prevent_initial_call=True,
    )
    def toggle(_clicks, lines):
        trig = ctx.triggered_id
        if not isinstance(trig, dict):
            return no_update
        if not (ctx.triggered and ctx.triggered[0].get("value")):
            return no_update
        lines = dict(lines or {})
        all_lines = list(lines.get("all", []))
        active = set(lines.get("active", all_lines))
        idx = trig.get("index")
        if idx == "__all__":
            active = set(all_lines)
        elif idx in active:
            active.discard(idx)
        else:
            active.add(idx)
        if not active:  # never allow empty selection
            active = set(all_lines)
        lines["active"] = [ln for ln in all_lines if ln in active]
        return lines


def _register_selection(app) -> None:
    @app.callback(
        Output(IDs.STORE_SELECTION, "data"),
        Input({"type": "mv-sta-row", "index": ALL}, "n_clicks"),
        Input(IDs.CANVAS_GRAPH, "clickData"),
        prevent_initial_call=True,
    )
    def select(_row_clicks, click_data):
        trig = ctx.triggered_id
        if trig == IDs.CANVAS_GRAPH:
            if not click_data:
                return no_update
            pt = (click_data.get("points") or [{}])[0]
            sid = pt.get("text") or pt.get("customdata")
            if isinstance(sid, list):
                sid = sid[0] if sid else None
            return {"station_id": str(sid)} if sid else no_update
        if isinstance(trig, dict) and trig.get("type") == "mv-sta-row":
            if not (ctx.triggered and ctx.triggered[0].get("value")):
                return no_update
            return {"station_id": str(trig.get("index"))}
        return no_update


def _register_inspector(app) -> None:
    @app.callback(
        Output(IDs.STATION_INSPECT, "children"),
        Input(IDs.STORE_SELECTION, "data"),
        Input(IDs.STORE_CONTROLS, "data"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def inspect(selection, controls, store):
        sid = (selection or {}).get("station_id")
        if not sid or not store:
            return html.Div("Click a station on the map.",
                            className="mv-empty")
        rec = next(
            (r for r in store.get("station_records", [])
             if str(r.get("ID")) == str(sid)),
            None,
        )
        if rec is None:
            return html.Div("Station not found.", className="mv-empty")

        def field(label, value):
            return html.Div(
                [html.Span(label, className="mv-insp-k"),
                 html.Span(str(value), className="mv-insp-v")],
                className="mv-insp-row",
            )

        rows = [
            html.Div(str(rec.get("ID")), className="mv-insp-title"),
            field("Line", rec.get("Line", "—")),
        ]
        rows += _coordinate_fields(rec, controls or {}, field)
        rows.append(field("Elevation", rec.get("Elevation", "—")))
        return html.Div(rows)


def _coordinate_fields(rec, controls, field) -> list:
    """Lat/Lon always, plus easting/northing in the selected display CRS."""
    lat, lon = rec.get("Latitude"), rec.get("Longitude")
    base = [
        field("Latitude", _fmt_coord(lat, "°")),
        field("Longitude", _fmt_coord(lon, "°")),
    ]
    mode = (controls or {}).get("crs_mode", "geo")
    if mode in (None, "geo") or lat is None or lon is None:
        return base
    try:
        from pycsamt.app.mapview._render import project_to_crs

        east, north, code = project_to_crs(
            [float(lon)], [float(lat)], mode,
            (controls or {}).get("utm_zone"),
            (controls or {}).get("utm_hem"),
            (controls or {}).get("epsg"),
        )
        if east is None:
            return base
        return base + [
            field(f"Easting (EPSG:{code})", _fmt_coord(float(east[0]), " m")),
            field("Northing", _fmt_coord(float(north[0]), " m")),
        ]
    except (TypeError, ValueError):
        return base


def _fmt_coord(value, suffix) -> str:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "—"
    if suffix == "°":
        return f"{v:.5f}°"
    return f"{v:,.1f}{suffix}"


_BASE_COLS = ["ID", "Line", "Latitude", "Longitude", "Elevation", "Index"]


def _register_table(app) -> None:
    @app.callback(
        Output("mv-station-table", "data"),
        Output("mv-station-table", "columns"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_LINES, "data"),
        Input(IDs.STORE_CONTROLS, "data"),
        Input(IDs.STORE_MASKED, "data"),
        prevent_initial_call=True,
    )
    def fill_table(store, lines, controls, masked):
        cols = [{"name": c, "id": c} for c in _BASE_COLS]
        if not store:
            return [], cols
        active = set((lines or {}).get("active", []))
        masked_set = {str(m) for m in (masked or [])}
        records = [
            dict(r) for r in store.get("station_records", [])
            if (not active or (r.get("Line") or "line") in active)
            and str(r.get("ID")) not in masked_set
        ]
        mode = (controls or {}).get("crs_mode", "geo")
        if mode not in (None, "geo") and records:
            records, cols = _add_projected_columns(records, controls, cols)
        return records, cols


def _add_projected_columns(records, controls, cols):
    """Append Easting/Northing (display CRS) columns to the table."""
    from pycsamt.app.mapview._render import project_to_crs

    lons = [r.get("Longitude") for r in records]
    lats = [r.get("Latitude") for r in records]
    try:
        safe_lon = [float(v) if v is not None else float("nan") for v in lons]
        safe_lat = [float(v) if v is not None else float("nan") for v in lats]
        east, north, code = project_to_crs(
            safe_lon, safe_lat,
            controls.get("crs_mode"), controls.get("utm_zone"),
            controls.get("utm_hem"), controls.get("epsg"),
        )
    except (TypeError, ValueError):
        return records, cols
    if east is None:
        return records, cols
    for i, rec in enumerate(records):
        rec[f"E ({code})"] = round(float(east[i]), 1)
        rec["N"] = round(float(north[i]), 1)
    cols = cols + [
        {"name": f"E ({code})", "id": f"E ({code})"},
        {"name": "N", "id": "N"},
    ]
    return records, cols
