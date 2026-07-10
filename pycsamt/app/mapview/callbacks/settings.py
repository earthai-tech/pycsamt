# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Station / sites settings: line visibility + per-site masking.

Masked stations and inactive lines are removed from every view (map,
3-D and the data table) via ``STORE_MASKED`` / ``STORE_LINES``.
"""

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


def register_settings(app) -> None:
    _register_toggle_canvas(app)
    _register_mask_toggle(app)
    _register_mask_actions(app)
    _register_unmask_chips(app)
    _register_active_lines(app)
    _register_summary(app)


def _register_toggle_canvas(app) -> None:
    @app.callback(
        Output(IDs.CANVAS_SETTINGS, "is_open"),
        Input(IDs.BTN_SETTINGS, "n_clicks"),
        State(IDs.CANVAS_SETTINGS, "is_open"),
        prevent_initial_call=True,
    )
    def toggle(_n, is_open):
        return not is_open


def _register_mask_toggle(app) -> None:
    @app.callback(
        Output(IDs.STORE_MASKED, "data", allow_duplicate=True),
        Input({"type": "mv-mask", "index": ALL}, "n_clicks"),
        State(IDs.STORE_MASKED, "data"),
        prevent_initial_call=True,
    )
    def toggle_mask(_clicks, masked):
        trig = ctx.triggered_id
        if not isinstance(trig, dict):
            return no_update
        if not (ctx.triggered and ctx.triggered[0].get("value")):
            return no_update
        sid = str(trig.get("index"))
        masked = [str(m) for m in (masked or [])]
        if sid in masked:
            masked.remove(sid)
        else:
            masked.append(sid)
        return masked


def _register_mask_actions(app) -> None:
    @app.callback(
        Output(IDs.STORE_MASKED, "data", allow_duplicate=True),
        Input(IDs.BTN_CLEAR_MASKS, "n_clicks"),
        prevent_initial_call=True,
    )
    def clear(_n):
        return []

    @app.callback(
        Output(IDs.STORE_MASKED, "data", allow_duplicate=True),
        Input(IDs.BTN_MASK_HIDDEN, "n_clicks"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_LINES, "data"),
        State(IDs.STORE_MASKED, "data"),
        prevent_initial_call=True,
    )
    def mask_hidden(_n, store, lines, masked):
        if not store:
            return no_update
        active = set((lines or {}).get("active", []))
        masked = {str(m) for m in (masked or [])}
        for rec in store.get("station_records", []):
            line = rec.get("Line") or "line"
            if active and line not in active:
                masked.add(str(rec.get("ID")))
        return sorted(masked)


def _register_unmask_chips(app) -> None:
    @app.callback(
        Output(IDs.STORE_MASKED, "data", allow_duplicate=True),
        Input({"type": "mv-unmask", "index": ALL}, "n_clicks"),
        State(IDs.STORE_MASKED, "data"),
        prevent_initial_call=True,
    )
    def unmask(_clicks, masked):
        trig = ctx.triggered_id
        if not isinstance(trig, dict):
            return no_update
        if not (ctx.triggered and ctx.triggered[0].get("value")):
            return no_update
        sid = str(trig.get("index"))
        return [m for m in (masked or []) if str(m) != sid]


def _register_active_lines(app) -> None:
    @app.callback(
        Output(IDs.SET_LINES, "options"),
        Output(IDs.SET_LINES, "value"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_LINES, "data"),
        prevent_initial_call=False,
    )
    def populate(store, lines):
        all_lines = (lines or {}).get("all") or (store or {}).get("lines", [])
        active = (lines or {}).get("active", all_lines)
        opts = [{"label": ln, "value": ln} for ln in all_lines]
        return opts, list(active)

    @app.callback(
        Output(IDs.STORE_LINES, "data", allow_duplicate=True),
        Input(IDs.SET_LINES, "value"),
        State(IDs.STORE_LINES, "data"),
        prevent_initial_call=True,
    )
    def set_active(value, lines):
        lines = dict(lines or {})
        all_lines = list(lines.get("all", []))
        chosen = set(value or [])
        new_active = [ln for ln in all_lines if ln in chosen] or all_lines
        if new_active == lines.get("active"):
            return no_update
        lines["active"] = new_active
        return lines


def _register_summary(app) -> None:
    @app.callback(
        Output(IDs.SET_SUMMARY, "children"),
        Output(IDs.SET_MASKED_LIST, "children"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_LINES, "data"),
        Input(IDs.STORE_MASKED, "data"),
        prevent_initial_call=False,
    )
    def summary(store, lines, masked):
        if not store or not store.get("station_records"):
            return html.Span("No data loaded.", className="mv-crs-info"), ""
        masked = [str(m) for m in (masked or [])]
        active = set((lines or {}).get("active", []))
        recs = store.get("station_records", [])
        mset = set(masked)
        visible = sum(
            1 for r in recs
            if (not active or (r.get("Line") or "line") in active)
            and str(r.get("ID")) not in mset
        )
        summ = html.Div(
            [
                _stat_pill(str(visible), "visible"),
                _stat_pill(str(len(masked)), "masked"),
                _stat_pill(str(store.get("n_stations", 0)), "total"),
            ],
            className="mv-set-stats",
        )
        if not masked:
            chips = html.Div("No masked stations.", className="mv-empty")
        else:
            chips = html.Div(
                [
                    html.Span(
                        [
                            sid,
                            html.Button(
                                html.I(className="bi bi-x"),
                                id={"type": "mv-unmask", "index": sid},
                                className="mv-chip-x",
                                n_clicks=0,
                                title="Unmask",
                            ),
                        ],
                        className="mv-masked-chip",
                    )
                    for sid in masked
                ],
                className="mv-masked-chips",
            )
        return summ, chips


def _stat_pill(value, label):
    return html.Div(
        [html.Div(value, className="mv-set-stat-v"),
         html.Div(label, className="mv-set-stat-l")],
        className="mv-set-stat",
    )
