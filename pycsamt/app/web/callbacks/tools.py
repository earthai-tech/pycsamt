# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/tools.py — Tools off-canvas panel.

Responsibilities
----------------
- Open the Tools offcanvas when any item in the navbar Tools dropdown is clicked
  OR when a tool-picker rail button inside the offcanvas is clicked.
- Track the active tool in ACTIVE_TOOL store.
- Render the correct tool body panel inside the offcanvas.
- Execute each tool and stream results back into the panel.
"""
from __future__ import annotations

import io
import os

import dash_bootstrap_components as dbc
from dash import (
    ALL,
    Input,
    Output,
    State,
    dash_table,
    dcc,
    html,
    no_update,
    set_props,
)
from dash import (
    callback_context as ctx,
)
from dash.exceptions import PreventUpdate

from pycsamt.app.web.layout import (
    _TOOL_BY_ID,
    _TOOL_REGISTRY,
    IDs,
    _icon,
    _tool_panel,
)
from pycsamt.app.web.utils import filter_sites_by_lines

# ── Shared helpers ────────────────────────────────────────────────────────────

def _styled_table(df) -> dash_table.DataTable:
    """Return a compact, theme-neutral DataTable from a DataFrame."""
    from dash import dash_table
    return dash_table.DataTable(
        columns=[{"name": c, "id": c} for c in df.columns],
        data=df.to_dict("records"),
        page_size=12,
        sort_action="native",
        filter_action="native",
        style_table={"overflowX": "auto"},
        style_cell={
            "fontSize": "11px",
            "padding": "3px 6px",
            "border": "1px solid #313244",
            "backgroundColor": "#1e1e2e",
            "color": "#cdd6f4",
            "textAlign": "left",
        },
        style_header={
            "backgroundColor": "#181825",
            "fontWeight": "bold",
            "color": "#89b4fa",
            "border": "1px solid #45475a",
        },
        style_data_conditional=[{
            "if": {"row_index": "odd"},
            "backgroundColor": "#181825",
        }],
    )


def _freq_before_after_heatmaps(
    before, after, mpl_wrap_fn, style_ax_fn,
    bg, fg, grid, ax_bg, heatmap_fn,
) -> list:
    """Return a list of Div children showing before/after coverage heatmaps."""
    import matplotlib
    matplotlib.use("Agg")
    children = []

    def _counts(sites_obj):
        try:
            from pycsamt.emtools._core import (
                _get_z_block,
                _iter_items,
            )
            total = sum(
                1 for _, ed in enumerate(_iter_items(sites_obj))
                for Z, z, fr in [_get_z_block(ed)]
                if Z is not None
                for _ in [None]
            )
            freqs = set()
            for _i, ed in enumerate(_iter_items(sites_obj)):
                Z, z, fr = _get_z_block(ed)
                if Z is not None:
                    freqs.update(fr.tolist())
            return total, len(freqs)
        except Exception:
            return 0, 0

    n_sta_b, n_fr_b = _counts(before)
    n_sta_a, n_fr_a = _counts(after)
    children.append(dbc.Alert([
        html.I(className="bi bi-check-circle me-2"),
        f"Before: {n_sta_b} stations · {n_fr_b} unique freqs  →  "
        f"After: {n_sta_a} stations · {n_fr_a} unique freqs",
    ], color="success", className="py-2 px-3 small mb-2"))

    try:
        children.append(html.Div([
            html.P("Coverage quality — before", className="text-muted small mb-1"),
            mpl_wrap_fn(heatmap_fn, before, axis="period", verbose=0),
        ], className="mb-2"))
        children.append(html.Div([
            html.P("Coverage quality — after", className="text-muted small mb-1"),
            mpl_wrap_fn(heatmap_fn, after, axis="period", verbose=0),
        ], className="mb-3"))
    except Exception as exc:
        children.append(html.Span(str(exc), className="text-danger small"))
    return children


def _no_data_msg() -> html.Div:
    return html.Div(
        [html.I(className="bi bi-exclamation-triangle me-2 text-warning"),
         "Load survey data first."],
        className="tool-no-data",
    )


def _run_btn(btn_id: str, label: str, icon: str = "play-fill") -> html.Button:
    return html.Button(
        [html.I(className=f"bi bi-{icon} me-1"), label],
        id=btn_id,
        className="btn btn-primary btn-sm",
        n_clicks=0,
    )


def _action_bar(*children) -> html.Div:
    return html.Div(list(children), className="tool-action-bar")


def _plotly_config(filename: str = "pycsamt_plot") -> dict:
    return {
        "displayModeBar": True,
        "displaylogo": False,
        "responsive": True,
        "toImageButtonOptions": {
            "format": "png",
            "filename": filename,
            "height": None,
            "width": None,
            "scale": 2,
        },
        "modeBarButtonsToRemove": ["lasso2d", "select2d"],
    }


def _out_area(
    div_id: str,
    min_h: int = 160,
    cls: str = "mt-3 fig-area",
    children=None,
    last_output=None,
) -> html.Div:
    if children is None and last_output is not None:
        children = _restore_from_store(last_output)
    return html.Div(
        children=children,
        id=div_id,
        className=cls,
        style={"minHeight": f"{min_h}px"},
    )


def _restore_from_store(stored: dict | None):
    """Convert a serialised tool-outputs-store entry back to Dash components."""
    if not stored:
        return None
    t = stored.get("type")
    try:
        if t == "fig_json":
            import plotly.io as pio
            fig = pio.from_json(stored["fig"])
            return dcc.Graph(
                figure=fig,
                config=_plotly_config("pycsamt_restored"),
                className="tool-graph",
            )
        if t == "strike_payload":
            return _render_strike_result(stored["payload"], stored.get("theme"))
        if t == "validator_rows":
            rows = stored.get("rows", [])
            cols = stored.get("cols", [])
            summary = stored.get("summary", {})
            if not rows or not cols:
                return None
            chips = html.Div(
                [
                    _metric_chip("Scope",       summary.get("scope", "—")),
                    _metric_chip("Validated",   str(summary.get("total", len(rows)))),
                    _metric_chip("Pass",        str(summary.get("pass", 0))),
                    _metric_chip("Warn / Fail",
                                 f"{summary.get('warn', 0)} / {summary.get('fail', 0)}"),
                    html.Span(" (restored)", className="small text-muted ms-2"),
                ],
                className="tool-metric-row",
            )
            tbl = dash_table.DataTable(
                columns=[{"name": c, "id": c} for c in cols],
                data=rows,
                page_size=12,
                sort_action="native",
                filter_action="native",
                style_table={"overflowX": "auto", "maxHeight": "420px",
                             "overflowY": "auto"},
                style_cell={
                    "fontSize": "11px", "fontFamily": "Inter, sans-serif",
                    "padding": "6px", "backgroundColor": "var(--mantle)",
                    "color": "var(--text)", "border": "1px solid var(--surface0)",
                },
                style_header={"backgroundColor": "var(--surface0)", "fontWeight": "700"},
            )
            return html.Div([chips, tbl])
        if t == "text":
            return html.Span(stored.get("content", ""),
                             className=stored.get("cls", "small"))
        if t == "html":
            imgs = stored.get("imgs") or []
            if not imgs:
                return dbc.Alert(
                    "Last run complete — click Run to regenerate results.",
                    color="secondary", className="small py-2",
                )
            return html.Div([
                dbc.Alert(
                    "Results from last run (restored).",
                    color="secondary", className="small py-2 mb-2",
                ),
                *[html.Img(
                    src=f"data:image/png;base64,{b64}",
                    style={"width": "100%", "borderRadius": "6px",
                           "display": "block", "marginBottom": "12px"},
                ) for b64 in imgs],
            ])
        if t == "lm":
            model = stored.get("model") or {}
            imgs  = stored.get("imgs") or []
            solver = stored.get("solver", "")
            rhos   = model.get("resistivity") or []
            layer_desc = "  ·  ".join(
                f"L{i+1}: {r:.1f} Ω·m" for i, r in enumerate(rhos)
            )
            banner = dbc.Alert([
                html.B(f"{len(rhos)}-layer model"),
                html.Span(f"  [{model.get('name', '')}]  ",
                          className="text-muted") if model.get("name") else " ",
                html.Br(),
                html.Small(layer_desc or "—", className="text-muted"),
                html.Span(f"  ·  solver: {solver}", className="text-muted small"),
                html.Span(" (restored)", className="small text-muted ms-2"),
            ], color="info", className="py-2 px-3 small mb-2")
            img_divs = [
                html.Img(
                    src=f"data:image/png;base64,{b64}",
                    style={"width": "100%", "borderRadius": "6px",
                           "display": "block", "marginBottom": "12px"},
                ) for b64 in imgs
            ]
            return html.Div([banner, *img_divs])
    except Exception:
        return None
    return None


def _path_row(
    input_id: str,
    placeholder: str,
    *,
    file_btn_id: str | None = None,
    folder_btn_id: str | None = None,
    className: str = "mb-2",
) -> dbc.InputGroup:
    buttons = []
    if file_btn_id:
        buttons.append(
            dbc.Button(
                [html.I(className="bi bi-file-earmark me-1"), "File"],
                id=file_btn_id,
                color="secondary",
                outline=True,
                size="sm",
                title="Browse for a server-side file",
            )
        )
    if folder_btn_id:
        buttons.append(
            dbc.Button(
                [html.I(className="bi bi-folder2-open me-1"), "Folder"],
                id=folder_btn_id,
                color="secondary",
                outline=True,
                size="sm",
                title="Browse for a server-side folder",
            )
        )
    return dbc.InputGroup(
        [
            dbc.Input(id=input_id, placeholder=placeholder, size="sm"),
            *buttons,
        ],
        size="sm",
        className=className,
    )


def _active_line_names(store_data: dict | None, active_lines_store: dict | None) -> list[str]:
    records = (store_data or {}).get("station_records", [])
    active = list((active_lines_store or {}).get("active", []) or [])
    if active:
        return active
    return sorted({r.get("Line", "") for r in records if r.get("Line")})


def _station_options_for_lines(
    store_data: dict | None,
    active_lines_store: dict | None,
    selected_lines: list[str] | None = None,
) -> list[dict]:
    records = (store_data or {}).get("station_records", [])
    active_lines = _active_line_names(store_data, active_lines_store)
    active_set = set(active_lines)
    selected_set = set(selected_lines or [])
    opts = []
    for r in records:
        line = r.get("Line", "")
        if active_set and line and line not in active_set:
            continue
        if selected_set and line not in selected_set:
            continue
        sid = r.get("ID") or r.get("station") or r.get("name")
        if sid:
            label = f"{sid} · {line}" if line else sid
            opts.append({"label": label, "value": sid})
    return opts


def _station_line_map(store_data: dict | None) -> dict[str, str]:
    out = {}
    for r in (store_data or {}).get("station_records", []):
        sid = r.get("ID") or r.get("station") or r.get("name")
        if sid:
            out[str(sid)] = str(r.get("Line", "") or "Unassigned")
    return out


def _filter_sites_by_station_ids(sites, station_ids):
    if not station_ids:
        return sites
    try:
        from pycsamt.emtools._core import (
            _iter_items,
            _name,
            ensure_sites,
        )
        wanted = {str(s) for s in station_ids}
        items = [
            getattr(ed, "edi", ed)
            for i, ed in enumerate(_iter_items(sites))
            if _name(ed, i) in wanted
        ]
        if not items:
            return None
        return ensure_sites(
            items, recursive=False, on_dup="replace", strict=False, verbose=0
        )
    except Exception:
        return sites


def _scoped_sites(
    session_id,
    store_data,
    active_lines_store,
    selected_lines=None,
    selected_stations=None,
):
    from pycsamt.app.web.cache import cache_get

    sites = cache_get(session_id)
    if sites is None:
        return None, "Session expired — reload data."
    records = (store_data or {}).get("station_records", [])
    active = (active_lines_store or {}).get("active", [])
    if active:
        sites = filter_sites_by_lines(sites, records, active)
        if sites is None:
            return None, "No active survey lines. Enable at least one line."
    if selected_lines:
        sites = filter_sites_by_lines(sites, records, selected_lines)
        if sites is None:
            return None, "No stations remain for the selected line scope."
    if selected_stations:
        sites = _filter_sites_by_station_ids(sites, selected_stations)
        if sites is None:
            return None, "No stations remain for the selected station scope."
    return sites, None


# ── Tool body builders ────────────────────────────────────────────────────────

def _build_tool_body(
    tool_id: str,
    store_data: dict | None,
    active_lines_store: dict | None = None,
    run_store: dict | None = None,
    theme: str | None = "dark",
    last_output: dict | None = None,
) -> html.Div:
    n  = (store_data or {}).get("n_stations", 0)
    nd = _no_data_msg()

    builders = {
        "strike":        lambda: _strike_body(
            n, nd, store_data, active_lines_store, run_store, theme, last_output
        ),
        "validator":     lambda: _validator_body(
            n, nd, store_data, active_lines_store, last_output
        ),
        "converter":     lambda: _converter_body(),
        "coords":        lambda: _coords_body(),
        "dimensionality":lambda: _dimensionality_body(
            n, nd, store_data, active_lines_store, last_output
        ),
        "freq-editor":   lambda: _freq_editor_body(
            n, nd, store_data, active_lines_store, last_output
        ),
        "batch":         lambda: _batch_body(n, nd),
        "elevation":     lambda: _elevation_body(n, nd, store_data),
        "response":      lambda: _response_body(n, nd, store_data, last_output),
        "strike-profile":lambda: _strike_profile_body(n, nd, store_data, last_output),
        "phase-tensor":  lambda: _phase_tensor_body(n, nd, store_data, last_output),
        "layered-model": lambda: _layered_model_body(theme, last_output),
    }
    fn = builders.get(tool_id)
    return fn() if fn else _placeholder_body(tool_id)


# ── Implemented tools ─────────────────────────────────────────────────────────

def _strike_body(
    n: int,
    no_data: html.Div,
    store_data: dict | None,
    active_lines_store: dict | None,
    run_store: dict | None = None,
    theme: str | None = "dark",
    last_output: dict | None = None,
) -> html.Div:
    if not n:
        return no_data
    line_opts = [
        {"label": ln, "value": ln}
        for ln in _active_line_names(store_data, active_lines_store)
    ]
    station_opts = _station_options_for_lines(store_data, active_lines_store)
    return html.Div([
        html.P(
            "Estimate geoelectric strike from the currently loaded survey. "
            "By default the tool uses all active lines from the main Lines panel.",
            className="text-muted small mb-2",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label("Line scope", className="ctrl-label"),
                        dcc.Dropdown(
                            id="tool-strike-lines",
                            options=line_opts,
                            value=None,
                            multi=True,
                            placeholder="All active lines",
                            className="mb-2",
                        ),
                    ],
                    width=6,
                ),
                dbc.Col(
                    [
                        html.Label("Station scope", className="ctrl-label"),
                        dcc.Dropdown(
                            id="tool-strike-stations",
                            options=station_opts,
                            value=None,
                            multi=True,
                            placeholder="All stations",
                            className="mb-2",
                            optionHeight=26,
                        ),
                    ],
                    width=6,
                ),
            ],
            className="tool-control-row",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label("Strike method", className="ctrl-label"),
                        dcc.Dropdown(
                            id="tool-strike-method",
                            options=[
                                {"label": "Consensus", "value": "consensus"},
                                {"label": "Z sweep", "value": "sweep"},
                                {"label": "Phase tensor", "value": "pt"},
                            ],
                            value="consensus", clearable=False, className="mb-2",
                        ),
                    ],
                    width=6,
                ),
                dbc.Col(
                    [
                        html.Label("Period / frequency", className="ctrl-label"),
                        dcc.Dropdown(
                            id="tool-strike-freq",
                            options=[
                                {"label": "All frequencies", "value": "all"},
                                {"label": "≥ 100 Hz", "value": "high"},
                                {"label": "≤ 1 Hz", "value": "low"},
                            ],
                            value="all", clearable=False, className="mb-2",
                        ),
                    ],
                    width=6,
                ),
            ],
            className="tool-control-row",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label("Angle step (°)", className="ctrl-label"),
                        dbc.Input(
                            id="tool-strike-step",
                            type="number",
                            min=0.25,
                            max=10.0,
                            step=0.25,
                            value=1.0,
                            size="sm",
                        ),
                    ],
                    width=6,
                ),
                dbc.Col(
                    [
                        html.Label("Max rows", className="ctrl-label"),
                        dbc.Input(
                            id="tool-strike-maxrows",
                            type="number",
                            min=5,
                            max=200,
                            step=5,
                            value=40,
                            size="sm",
                        ),
                    ],
                    width=6,
                ),
            ],
            className="mb-2",
        ),
        _action_bar(_run_btn("tool-strike-run", "Run Strike Analysis")),
        _out_area(
            "tool-strike-out",
            min_h=360,
            children=_render_stored_tool_result("strike", run_store, theme),
            last_output=last_output,
        ),
    ])


def _validator_body(
    n: int,
    no_data: html.Div,
    store_data: dict | None,
    active_lines_store: dict | None,
    last_output: dict | None = None,
) -> html.Div:
    if not n:
        return no_data
    line_opts = [
        {"label": ln, "value": ln}
        for ln in _active_line_names(store_data, active_lines_store)
    ]
    station_opts = _station_options_for_lines(store_data, active_lines_store)
    return html.Div([
        html.P(
            "Validate station metadata, coordinates, impedance tensors, "
            "frequency grids, error arrays, and optional tipper availability.",
            className="text-muted small mb-2",
        ),
        html.Label("Line scope", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-valid-lines",
            options=line_opts,
            value=None,
            multi=True,
            placeholder="All active lines",
            className="mb-2",
        ),
        html.Label("Station scope", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-valid-stations",
            options=station_opts,
            value=None,
            multi=True,
            placeholder="All stations in selected lines",
            className="mb-2",
            optionHeight=26,
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label("Show", className="ctrl-label"),
                        dcc.Dropdown(
                            id="tool-valid-severity",
                            options=[
                                {"label": "All statuses", "value": "all"},
                                {"label": "Warnings + failures", "value": "issues"},
                                {"label": "Failures only", "value": "fail"},
                                {"label": "Pass only", "value": "pass"},
                            ],
                            value="all",
                            clearable=False,
                        ),
                    ],
                    width=7,
                ),
                dbc.Col(
                    [
                        html.Label("Min N freq", className="ctrl-label"),
                        dbc.Input(
                            id="tool-valid-minfreq",
                            type="number",
                            min=1,
                            max=10000,
                            step=1,
                            value=3,
                            size="sm",
                        ),
                    ],
                    width=5,
                ),
            ],
            className="mb-2",
        ),
        html.Label("Checks", className="ctrl-label"),
        dbc.Checklist(
            id="tool-valid-checks",
            options=[
                {"label": "Coordinates", "value": "coords"},
                {"label": "Z tensor", "value": "z"},
                {"label": "Z errors", "value": "errors"},
                {"label": "Frequency grid", "value": "freq"},
                {"label": "NaN rows", "value": "nan"},
                {"label": "Duplicates", "value": "duplicates"},
                {"label": "Tipper", "value": "tipper"},
            ],
            value=["coords", "z", "errors", "freq", "nan", "duplicates"],
            inline=False,
            className="tool-checklist mb-2",
        ),
        _action_bar(_run_btn("tool-valid-run", "Run Validation", "clipboard-check")),
        _out_area("tool-valid-out", min_h=360, cls="mt-3", last_output=last_output),
    ])


def _converter_body() -> html.Div:
    return html.Div([
        html.P(
            "Convert external processing outputs into pyCSAMT EDI files. "
            "AVG conversion can attach a Zonge .stn profile/topography file "
            "and convert UTM coordinates to latitude/longitude.",
            className="text-muted small mb-2",
        ),
        html.Label("Conversion", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-conv-type",
            options=[
                {"label": "AVG → EDI", "value": "AVG -> EDI"},
                {"label": "J → EDI", "value": "J -> EDI"},
                {"label": "Spectra EDI → impedance EDI", "value": "Spectra -> EDI"},
            ],
            value="AVG -> EDI",
            clearable=False,
            className="mb-2",
        ),
        html.Label("Source path", className="ctrl-label"),
        _path_row(
            "tool-conv-source",
            "/data/avg/K1.AVG or folder/",
            file_btn_id={
                "type": "tools-path-browse",
                "target": "tool-conv-source",
                "mode": "file",
                "title": "Select conversion source file",
            },
            folder_btn_id={
                "type": "tools-path-browse",
                "target": "tool-conv-source",
                "mode": "folder",
                "title": "Select conversion source folder",
            },
        ),
        html.Label("Output folder", className="ctrl-label"),
        _path_row(
            "tool-conv-output",
            "Optional folder for converted EDIs",
            folder_btn_id={
                "type": "tools-path-browse",
                "target": "tool-conv-output",
                "mode": "folder",
                "title": "Select converter output folder",
            },
        ),

        html.Div(
            [
                html.Label("AVG station profile / topography", className="ctrl-label"),
                _path_row(
                    "tool-conv-stn",
                    "/data/avg/K1.stn  (optional but recommended)",
                    file_btn_id={
                        "type": "tools-path-browse",
                        "target": "tool-conv-stn",
                        "mode": "file",
                        "title": "Select AVG station profile / topography file",
                    },
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.Label("UTM zone", className="ctrl-label"),
                                dbc.Input(
                                    id="tool-conv-utm",
                                    placeholder="49N",
                                    size="sm",
                                ),
                            ],
                            width=6,
                        ),
                        dbc.Col(
                            [
                                html.Label("EPSG", className="ctrl-label"),
                                dbc.Input(
                                    id="tool-conv-epsg",
                                    placeholder="32649",
                                    type="number",
                                    size="sm",
                                ),
                            ],
                            width=6,
                        ),
                    ],
                    className="mb-2",
                ),
                dbc.Switch(
                    id="tool-conv-convert-coords",
                    label="Convert station UTM/easting/northing to lat/lon",
                    value=True,
                    className="mb-2",
                ),
            ],
            id="tool-conv-avg-options",
            className="tool-option-box",
        ),

        html.Div(
            [
                html.Label("Spectra options", className="ctrl-label"),
                dbc.Row(
                    [
                        dbc.Col(
                            dbc.Input(
                                id="tool-conv-elabels",
                                placeholder="E labels: EX,EY",
                                value="EX,EY",
                                size="sm",
                            ),
                            width=6,
                        ),
                        dbc.Col(
                            dbc.Input(
                                id="tool-conv-hlabels",
                                placeholder="H labels: HX,HY",
                                value="HX,HY",
                                size="sm",
                            ),
                            width=6,
                        ),
                    ],
                    className="mb-2",
                ),
                dbc.Input(
                    id="tool-conv-suffix",
                    placeholder="Optional station suffix, e.g. _IMP",
                    size="sm",
                    className="mb-2",
                ),
                dbc.Checklist(
                    id="tool-conv-spectra-flags",
                    options=[
                        {"label": "Estimate errors", "value": "estimate_error"},
                        {"label": "Use remote reference", "value": "use_remote"},
                        {"label": "Skip files that fail", "value": "skip_errors"},
                    ],
                    value=["skip_errors"],
                    className="tool-checklist mb-2",
                ),
            ],
            id="tool-conv-spectra-options",
            className="tool-option-box",
            style={"display": "none"},
        ),

        html.Div(
            [
                html.Label("Core transform policy", className="ctrl-label"),
                dbc.Row(
                    [
                        dbc.Col(
                            dcc.Dropdown(
                                id="tool-conv-freq-order",
                                options=[
                                    {"label": "Descending frequency", "value": "descending"},
                                    {"label": "Ascending frequency", "value": "ascending"},
                                ],
                                value="descending",
                                clearable=False,
                            ),
                            width=7,
                        ),
                        dbc.Col(
                            dbc.Input(
                                id="tool-conv-freq-tol",
                                placeholder="freq tol",
                                type="number",
                                size="sm",
                            ),
                            width=5,
                        ),
                    ],
                    className="mb-2",
                ),
                dbc.Checklist(
                    id="tool-conv-core-flags",
                    options=[
                        {"label": "Compute Z from rho/phase if needed", "value": "compute_z"},
                        {"label": "Recompute rho/phase from Z", "value": "compute_rho_phi"},
                        {"label": "Load converted EDIs into current session", "value": "load_session"},
                    ],
                    value=["compute_z", "compute_rho_phi"],
                    className="tool-checklist mb-2",
                ),
                dbc.Input(
                    id="tool-conv-station-name",
                    placeholder="Optional station name override for single-file input",
                    size="sm",
                    className="mb-2",
                ),
            ],
            className="tool-option-box",
        ),
        _action_bar(_run_btn("tool-conv-run", "Convert to EDI", "download")),
        _out_area("tool-conv-out", min_h=340, cls="mt-2"),
    ])


def _coords_body() -> html.Div:

    def _lbl(t):
        return html.Div(t, className="ctrl-label mt-2 mb-1")

    # ── Mode pills ────────────────────────────────────────────────────────
    mode_bar = html.Div(
        [
            html.Button(
                [html.I(className="bi bi-cursor me-1"), "Single"],
                id="tool-coord-mode-single", n_clicks=0,
                className="coord-mode-btn active",
            ),
            html.Button(
                [html.I(className="bi bi-file-earmark-spreadsheet me-1"), "Batch CSV"],
                id="tool-coord-mode-batch", n_clicks=0,
                className="coord-mode-btn",
            ),
            html.Button(
                [html.I(className="bi bi-broadcast me-1"), "Survey"],
                id="tool-coord-mode-survey", n_clicks=0,
                className="coord-mode-btn",
            ),
        ],
        className="coord-mode-bar mb-3",
    )

    # ── Direction & common params ─────────────────────────────────────────
    direction_row = html.Div(
        [
            html.Div(
                [
                    _lbl("Conversion"),
                    dbc.Select(
                        id="tool-coord-direction",
                        options=[
                            {"label": "UTM → Lat / Lon",   "value": "utm2ll"},
                            {"label": "Lat / Lon → UTM",   "value": "ll2utm"},
                            {"label": "EPSG → EPSG",       "value": "epsg2epsg"},
                        ],
                        value="utm2ll",
                        size="sm",
                    ),
                ],
                style={"flex": "1"},
            ),
            html.Div(
                [
                    _lbl("Datum"),
                    dbc.Select(
                        id="tool-coord-datum",
                        options=[
                            {"label": "WGS84",  "value": "WGS84"},
                            {"label": "NAD83",  "value": "NAD83"},
                            {"label": "NAD27",  "value": "NAD27"},
                            {"label": "GDA94",  "value": "GDA94"},
                        ],
                        value="WGS84",
                        size="sm",
                    ),
                ],
                style={"flex": "1"},
            ),
        ],
        className="d-flex gap-2 mb-2",
    )

    # ── Single-point input panel (UTM → LL) ───────────────────────────────
    single_utm_panel = html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            _lbl("Easting (m)"),
                            dbc.Input(id="tool-coord-easting", type="number",
                                      placeholder="e.g. 500000", size="sm"),
                        ],
                        style={"flex": "1"},
                    ),
                    html.Div(
                        [
                            _lbl("Northing (m)"),
                            dbc.Input(id="tool-coord-northing", type="number",
                                      placeholder="e.g. 4500000", size="sm"),
                        ],
                        style={"flex": "1"},
                    ),
                ],
                className="d-flex gap-2 mb-2",
            ),
            html.Div(
                [
                    _lbl("UTM Zone"),
                    dbc.Input(id="tool-coord-zone", type="text",
                              placeholder="e.g. 49N", size="sm"),
                ],
                id="tool-coord-zone-wrap",
                className="mb-2",
            ),
        ],
        id="tool-coord-single-utm",
    )

    # ── Single-point input panel (LL → UTM) ───────────────────────────────
    single_ll_panel = html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            _lbl("Latitude (°)"),
                            dbc.Input(id="tool-coord-lat", type="number",
                                      placeholder="e.g. 40.7128",
                                      min=-90, max=90, size="sm"),
                        ],
                        style={"flex": "1"},
                    ),
                    html.Div(
                        [
                            _lbl("Longitude (°)"),
                            dbc.Input(id="tool-coord-lon", type="number",
                                      placeholder="e.g. -74.006",
                                      min=-180, max=180, size="sm"),
                        ],
                        style={"flex": "1"},
                    ),
                ],
                className="d-flex gap-2 mb-2",
            ),
            html.Div(
                [
                    _lbl("Target UTM Zone (optional)"),
                    dbc.Input(id="tool-coord-target-zone", type="text",
                              placeholder="Auto-detect or e.g. 18N", size="sm"),
                ],
                className="mb-2",
            ),
        ],
        id="tool-coord-single-ll",
        style={"display": "none"},
    )

    # ── EPSG → EPSG panel ─────────────────────────────────────────────────
    epsg_panel = html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            _lbl("Source EPSG"),
                            dbc.Input(id="tool-coord-epsg-src", type="number",
                                      placeholder="e.g. 32649", size="sm"),
                        ],
                        style={"flex": "1"},
                    ),
                    html.Div(
                        [
                            _lbl("Target EPSG"),
                            dbc.Input(id="tool-coord-epsg-dst", type="number",
                                      placeholder="e.g. 4326", size="sm"),
                        ],
                        style={"flex": "1"},
                    ),
                ],
                className="d-flex gap-2 mb-2",
            ),
        ],
        id="tool-coord-epsg-panel",
        style={"display": "none"},
    )

    # ── Batch CSV upload panel ─────────────────────────────────────────────
    batch_panel = html.Div(
        [
            dcc.Upload(
                id="tool-coord-upload",
                children=html.Div(
                    [
                        html.I(className="bi bi-file-earmark-arrow-up me-2",
                               style={"fontSize": "20px", "color": "var(--teal)"}),
                        html.Div("Drop CSV / TXT here or click to browse",
                                 style={"fontSize": "11px", "color": "var(--sub0)"}),
                        html.Div("Accepts .csv, .txt with header row",
                                 style={"fontSize": "10px", "color": "var(--ov0)"}),
                    ],
                    className="d-flex flex-column align-items-center py-3",
                ),
                accept=".csv,.txt",
                multiple=False,
                className="coord-upload-zone mb-2",
            ),
            html.Div(id="tool-coord-upload-info", className="coord-upload-info mb-2"),
            html.Div(
                [
                    html.Div(
                        [_lbl("Easting col"),
                         dbc.Select(id="tool-coord-col-e", options=[], size="sm")],
                        style={"flex": "1"},
                    ),
                    html.Div(
                        [_lbl("Northing col"),
                         dbc.Select(id="tool-coord-col-n", options=[], size="sm")],
                        style={"flex": "1"},
                    ),
                    html.Div(
                        [_lbl("ID col"),
                         dbc.Select(id="tool-coord-col-id",
                                    options=[{"label": "— none —", "value": ""}],
                                    value="", size="sm")],
                        style={"flex": "1"},
                    ),
                ],
                className="d-flex gap-1 mb-2",
            ),
            _lbl("UTM Zone (batch)"),
            dbc.Input(id="tool-coord-zone-batch", type="text",
                      placeholder="e.g. 49N (or leave blank to auto)", size="sm",
                      className="mb-2"),
        ],
        id="tool-coord-batch-panel",
        style={"display": "none"},
    )

    # ── Survey stations panel ─────────────────────────────────────────────
    survey_panel = html.Div(
        [
            html.Div(
                [html.I(className="bi bi-broadcast me-1 text-teal"),
                 " Convert coordinates for all currently loaded survey stations."],
                className="small text-muted mb-2",
            ),
            _lbl("UTM Zone"),
            dbc.Input(id="tool-coord-zone-survey", type="text",
                      placeholder="e.g. 49N (or leave blank to auto)", size="sm",
                      className="mb-2"),
        ],
        id="tool-coord-survey-panel",
        style={"display": "none"},
    )

    # ── Action bar ────────────────────────────────────────────────────────
    action_bar = html.Div(
        [
            html.Button(
                [html.I(className="bi bi-arrow-left-right me-1"), "Convert"],
                id="tool-coord-run", n_clicks=0,
                className="btn btn-primary btn-sm flex-grow-1",
            ),
            html.Button(
                [html.I(className="bi bi-map me-1"), "Map"],
                id="tool-coord-map-btn", n_clicks=0,
                className="btn btn-outline-secondary btn-sm",
                title="Preview converted points on map",
            ),
            html.Button(
                [html.I(className="bi bi-download me-1"), "Export"],
                id="tool-coord-export-btn", n_clicks=0,
                className="btn btn-outline-secondary btn-sm",
                title="Export results as CSV",
                disabled=True,
            ),
        ],
        className="d-flex gap-1 mt-2 mb-2",
    )

    # ── Progress bar ──────────────────────────────────────────────────────
    progress = html.Div(
        [
            dbc.Progress(
                id="tool-coord-progress",
                value=0, max=100, striped=True, animated=True,
                style={"height": "6px"},
                className="mb-1",
                color="info",
            ),
            html.Div(id="tool-coord-progress-label",
                     className="coord-progress-label"),
        ],
        id="tool-coord-progress-wrap",
        style={"display": "none"},
        className="mb-2",
    )

    return html.Div(
        [
            html.P(
                "Convert coordinates — single point, batch CSV, or all survey "
                "stations. Supports UTM ↔ Lat/Lon and arbitrary EPSG reprojections.",
                className="text-muted small mb-2",
            ),
            mode_bar,
            direction_row,
            single_utm_panel,
            single_ll_panel,
            epsg_panel,
            batch_panel,
            survey_panel,
            action_bar,
            progress,
            _out_area("tool-coord-out", min_h=80, cls="mt-1"),
        ]
    )


def _dimensionality_body(
    n: int,
    no_data: html.Div,
    store_data: dict | None = None,
    active_lines_store: dict | None = None,
    last_output: dict | None = None,
) -> html.Div:
    if not n:
        return no_data
    line_opts    = [{"label": ln, "value": ln}
                    for ln in _active_line_names(store_data, active_lines_store)]
    station_opts = _station_options_for_lines(store_data, active_lines_store)
    return html.Div([
        html.P(
            "Rule-based dimensionality classification (1D / 2D / 3D) from "
            "the Bahr phase-sensitive skew and phase-tensor ellipticity.",
            className="text-muted small mb-2",
        ),

        # ── Scope selectors ───────────────────────────────────────────────
        dbc.Row([
            dbc.Col([
                html.Label("Line scope", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-dim-lines",
                    options=line_opts, value=None, multi=True,
                    placeholder="All active lines", className="mb-2",
                ),
            ], width=6),
            dbc.Col([
                html.Label("Station scope", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-dim-stations",
                    options=station_opts, value=None, multi=True,
                    placeholder="All stations", className="mb-2", optionHeight=26,
                ),
            ], width=6),
        ], className="tool-control-row"),

        # ── Classifier thresholds ─────────────────────────────────────────
        dbc.Row([
            dbc.Col([
                html.Label("Skew threshold (°)", className="ctrl-label"),
                dbc.Input(
                    id="tool-dim-skew-th", type="number",
                    value=3.0, min=0.5, max=30.0, step=0.5,
                    debounce=True, className="mb-2",
                ),
            ], width=4),
            dbc.Col([
                html.Label("Ellipticity threshold", className="ctrl-label"),
                dbc.Input(
                    id="tool-dim-ellipt-th", type="number",
                    value=0.2, min=0.01, max=1.0, step=0.01,
                    debounce=True, className="mb-2",
                ),
            ], width=4),
            dbc.Col([
                html.Label("Map period (s)", className="ctrl-label"),
                dbc.Input(
                    id="tool-dim-map-period", type="number",
                    value=10.0, min=1e-4, step=1.0,
                    debounce=True, className="mb-2",
                ),
            ], width=4),
        ], className="tool-control-row"),

        # ── Period filter (diagnostic scatter / stats) ────────────────────
        dbc.Row([
            dbc.Col([
                html.Label("Period range — min (s)", className="ctrl-label"),
                dbc.Input(
                    id="tool-dim-period-lo", type="number",
                    placeholder="auto", min=1e-5, debounce=True, className="mb-2",
                ),
            ], width=4),
            dbc.Col([
                html.Label("Period range — max (s)", className="ctrl-label"),
                dbc.Input(
                    id="tool-dim-period-hi", type="number",
                    placeholder="auto", min=1e-5, debounce=True, className="mb-2",
                ),
            ], width=4),
            dbc.Col([
                html.Label("Occupancy bands", className="ctrl-label"),
                dbc.Input(
                    id="tool-dim-n-bands", type="number",
                    value=24, min=4, max=64, step=2,
                    debounce=True, className="mb-2",
                ),
            ], width=4),
        ], className="tool-control-row"),

        # ── View selector ─────────────────────────────────────────────────
        html.Label("View", className="ctrl-label"),
        dbc.RadioItems(
            id="tool-dim-view",
            options=[
                {"label": "Confidence grid",   "value": "psection"},
                {"label": "Occupancy area",    "value": "occupancy"},
                {"label": "Spatial map",       "value": "map"},
                {"label": "β / ellipticity",   "value": "scatter"},
                {"label": "All views",         "value": "all"},
            ],
            value="psection",
            inline=True,
            className="mb-2 tool-radio",
            inputClassName="me-1",
        ),

        _action_bar(_run_btn("tool-dim-run", "Classify", "grid")),
        _out_area("tool-dim-out", min_h=260, last_output=last_output),
    ])


def _freq_editor_body(
    n: int,
    no_data: html.Div,
    store_data: dict | None = None,
    active_lines_store: dict | None = None,
    last_output: dict | None = None,
) -> html.Div:
    if not n:
        return no_data
    line_opts    = [{"label": ln, "value": ln}
                    for ln in _active_line_names(store_data, active_lines_store)]
    station_opts = _station_options_for_lines(store_data, active_lines_store)

    _conf_methods = [
        {"label": "Composite",  "value": "composite"},
        {"label": "Coherence",  "value": "coherence"},
        {"label": "Rel. error", "value": "relerr"},
        {"label": "Phase",      "value": "phase"},
        {"label": "Spatial",    "value": "spatial"},
    ]
    _show = {"display": "block"}
    _hide = {"display": "none"}

    return html.Div([
        html.P(
            "Inspect, edit, and export frequency grids. "
            "Choose a data source, select a line/station scope, pick an action "
            "and its parameters, then run — or just use Diagnostics to preview coverage.",
            className="text-muted small mb-2",
        ),

        # ── Data source ───────────────────────────────────────────────────
        dbc.Row([
            dbc.Col([
                html.Label("Data source", className="ctrl-label"),
                dbc.RadioItems(
                    id="tool-freq-source",
                    options=[
                        {"label": "Loaded session", "value": "session"},
                        {"label": "Folder on server", "value": "folder"},
                    ],
                    value="session", inline=True,
                    className="mb-1 tool-radio", inputClassName="me-1",
                ),
            ], width=12),
        ]),
        html.Div(
            [
                html.Label("EDI folder (server-side path)", className="ctrl-label"),
                _path_row(
                    "tool-freq-folder",
                    placeholder="/data/edi/",
                    folder_btn_id={
                        "type": "tools-path-browse",
                        "target": "tool-freq-folder",
                        "mode": "folder",
                        "title": "Select EDI folder",
                    },
                ),
            ],
            id="tool-freq-folder-row",
            style=_hide,
        ),

        # ── Scope ─────────────────────────────────────────────────────────
        dbc.Row([
            dbc.Col([
                html.Label("Line scope", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-freq-lines",
                    options=line_opts, value=None, multi=True,
                    placeholder="All active lines", className="mb-2",
                ),
            ], width=6),
            dbc.Col([
                html.Label("Station scope", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-freq-stations",
                    options=station_opts, value=None, multi=True,
                    placeholder="All stations", className="mb-2", optionHeight=26,
                ),
            ], width=6),
        ], className="tool-control-row"),

        # ── Action ────────────────────────────────────────────────────────
        html.Label("Action", className="ctrl-label"),
        dbc.RadioItems(
            id="tool-freq-action",
            options=[
                {"label": "Diagnostics",      "value": "diagnostics"},
                {"label": "Confidence edit",  "value": "confidence"},
                {"label": "Band selection",   "value": "band"},
                {"label": "Regrid (log)",     "value": "regrid"},
                {"label": "Grid align",       "value": "align"},
                {"label": "Smooth / Decimate","value": "smooth"},
            ],
            value="diagnostics", inline=True,
            className="mb-2 tool-radio", inputClassName="me-1",
        ),

        # ── Confidence edit params ────────────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([
                    html.Label("Edit mode", className="ctrl-label"),
                    dbc.RadioItems(
                        id="tool-freq-mode",
                        options=[
                            {"label": "Recover (interpolate)", "value": "recover"},
                            {"label": "Drop",                  "value": "drop"},
                            {"label": "Mask (NaN)",            "value": "mask"},
                        ],
                        value="recover", inline=False,
                        className="mb-2 tool-radio", inputClassName="me-1",
                    ),
                ], width=5),
                dbc.Col([
                    html.Label("Confidence method", className="ctrl-label"),
                    dcc.Dropdown(
                        id="tool-freq-method",
                        options=_conf_methods, value="composite",
                        clearable=False, className="mb-2",
                    ),
                    html.Label("Edit also", className="ctrl-label"),
                    dbc.RadioItems(
                        id="tool-freq-also",
                        options=[
                            {"label": "Z + Tipper", "value": "both"},
                            {"label": "Z only",     "value": "z"},
                            {"label": "Tipper only","value": "tipper"},
                        ],
                        value="both", inline=True,
                        className="mb-2 tool-radio", inputClassName="me-1",
                    ),
                ], width=7),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Label("Threshold (drop/mask)", className="ctrl-label"),
                    dcc.Slider(
                        id="tool-freq-thresh",
                        min=0.0, max=1.0, step=0.05, value=0.50,
                        marks={0: "0", 0.25: "0.25", 0.5: "0.5",
                               0.75: "0.75", 1: "1"},
                        className="mb-2",
                    ),
                ], width=12),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Label("CI high (recover)", className="ctrl-label"),
                    dbc.Input(id="tool-freq-ci-hi", type="number",
                              value=0.9, min=0.5, max=1.0, step=0.05,
                              debounce=True, className="mb-2"),
                ], width=4),
                dbc.Col([
                    html.Label("CI low (recover)", className="ctrl-label"),
                    dbc.Input(id="tool-freq-ci-lo", type="number",
                              value=0.5, min=0.0, max=0.9, step=0.05,
                              debounce=True, className="mb-2"),
                ], width=4),
                dbc.Col([
                    html.Label("Reject below CI low", className="ctrl-label"),
                    dbc.RadioItems(
                        id="tool-freq-reject",
                        options=[
                            {"label": "Drop",  "value": "drop"},
                            {"label": "Mask",  "value": "mask"},
                            {"label": "Keep",  "value": "keep"},
                        ],
                        value="drop", inline=True,
                        className="mb-2 tool-radio", inputClassName="me-1",
                    ),
                ], width=4),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Label("Interpolation (recover)", className="ctrl-label"),
                    dbc.RadioItems(
                        id="tool-freq-interp",
                        options=[
                            {"label": "Linear (log-f)", "value": "linear"},
                            {"label": "Nearest",        "value": "nearest"},
                        ],
                        value="linear", inline=True,
                        className="mb-2 tool-radio", inputClassName="me-1",
                    ),
                ], width=12),
            ]),
        ], id="tool-freq-conf-params", style=_hide),

        # ── Band selection params ─────────────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([
                    html.Label("fmin (Hz)", className="ctrl-label"),
                    dbc.Input(id="tool-freq-fmin", type="number",
                              placeholder="auto", min=1e-6, debounce=True,
                              className="mb-2"),
                ], width=4),
                dbc.Col([
                    html.Label("fmax (Hz)", className="ctrl-label"),
                    dbc.Input(id="tool-freq-fmax", type="number",
                              placeholder="auto", min=1e-6, debounce=True,
                              className="mb-2"),
                ], width=4),
                dbc.Col([
                    html.Label("Keep — explicit Hz (comma-sep)", className="ctrl-label"),
                    dbc.Input(id="tool-freq-keep", type="text",
                              placeholder="e.g. 100,50,10,1",
                              className="mb-2"),
                ], width=4),
            ]),
        ], id="tool-freq-band-params", style=_hide),

        # ── Regrid params ─────────────────────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([
                    html.Label("fmin (Hz)", className="ctrl-label"),
                    dbc.Input(id="tool-freq-rg-fmin", type="number",
                              placeholder="auto", min=1e-6, debounce=True,
                              className="mb-2"),
                ], width=3),
                dbc.Col([
                    html.Label("fmax (Hz)", className="ctrl-label"),
                    dbc.Input(id="tool-freq-rg-fmax", type="number",
                              placeholder="auto", min=1e-6, debounce=True,
                              className="mb-2"),
                ], width=3),
                dbc.Col([
                    html.Label("Points per decade", className="ctrl-label"),
                    dbc.Input(id="tool-freq-rg-ppd", type="number",
                              value=6, min=2, max=30, step=1,
                              debounce=True, className="mb-2"),
                ], width=3),
                dbc.Col([
                    html.Label("Interpolation", className="ctrl-label"),
                    dcc.Dropdown(
                        id="tool-freq-rg-method",
                        options=[{"label": "Nearest", "value": "nearest"},
                                 {"label": "Linear",  "value": "linear"}],
                        value="nearest", clearable=False, className="mb-2",
                    ),
                ], width=3),
            ]),
        ], id="tool-freq-regrid-params", style=_hide),

        # ── Grid-align params ─────────────────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([
                    html.Label("Alignment mode", className="ctrl-label"),
                    dbc.RadioItems(
                        id="tool-freq-align-mode",
                        options=[
                            {"label": "Union (all freqs)",       "value": "union"},
                            {"label": "Intersection (common)",   "value": "intersection"},
                            {"label": "Reference station",       "value": "ref"},
                        ],
                        value="union", inline=False,
                        className="mb-2 tool-radio", inputClassName="me-1",
                    ),
                ], width=5),
                dbc.Col([
                    html.Label("Reference station", className="ctrl-label"),
                    dcc.Dropdown(
                        id="tool-freq-align-ref",
                        options=station_opts, value=None,
                        placeholder="Required for ref mode",
                        className="mb-2",
                    ),
                ], width=7),
            ]),
        ], id="tool-freq-align-params", style=_hide),

        # ── Smooth / Decimate params ──────────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([
                    html.Label("Operation", className="ctrl-label"),
                    dbc.RadioItems(
                        id="tool-freq-smooth-op",
                        options=[
                            {"label": "Moving average (smooth)", "value": "smooth"},
                            {"label": "Step decimate",           "value": "decimate"},
                            {"label": "Drop duplicates",         "value": "dedup"},
                        ],
                        value="smooth", inline=False,
                        className="mb-2 tool-radio", inputClassName="me-1",
                    ),
                ], width=6),
                dbc.Col([
                    html.Label("k — window (smooth)", className="ctrl-label"),
                    dbc.Input(id="tool-freq-smooth-k", type="number",
                              value=3, min=2, max=15, step=1,
                              debounce=True, className="mb-2"),
                    html.Label("Step (decimate)", className="ctrl-label mt-1"),
                    dbc.Input(id="tool-freq-smooth-step", type="number",
                              value=2, min=2, max=10, step=1,
                              debounce=True, className="mb-2"),
                ], width=6),
            ]),
        ], id="tool-freq-smooth-params", style=_hide),

        # ── Export ────────────────────────────────────────────────────────
        html.Div([
            html.Hr(style={"borderColor": "#45475a", "margin": "8px 0"}),
            html.Label("Export edited EDIs to disk", className="ctrl-label"),
            dbc.Row([
                dbc.Col([
                    _path_row(
                        "tool-freq-out-path",
                        placeholder="/data/edited_edi/",
                        folder_btn_id={
                            "type": "tools-path-browse",
                            "target": "tool-freq-out-path",
                            "mode": "folder",
                            "title": "Select EDI output folder",
                        },
                    ),
                ], width=8),
                dbc.Col([
                    html.Label("Filename template", className="ctrl-label"),
                    dbc.Input(
                        id="tool-freq-out-template",
                        value="{station}.edi",
                        debounce=True, className="mb-2",
                    ),
                ], width=4),
            ]),
            dbc.Row([
                dbc.Col(
                    dbc.Button(
                        [html.I(className="bi bi-download me-1"), "Export EDIs"],
                        id="tool-freq-export-btn",
                        color="success", size="sm",
                        outline=True, disabled=True,
                        className="mb-2",
                    ),
                    width="auto",
                ),
                dbc.Col(
                    html.Div(id="tool-freq-export-out", className="small"),
                    width=True,
                    className="align-self-center",
                ),
            ], align="center"),
        ]),

        _action_bar(_run_btn("tool-freq-run", "Run", "play-fill")),
        _out_area("tool-freq-out", min_h=260, last_output=last_output),
    ])


# ── New tool bodies ───────────────────────────────────────────────────────────

def _batch_body(n: int, no_data: html.Div) -> html.Div:
    if not n:
        return no_data
    return html.Div([
        html.P(
            "Generate publication-ready figures from the active loaded survey "
            "and save them to a server-side folder.",
            className="text-muted small",
        ),
        html.Label("Output folder (server-side)", className="ctrl-label"),
        _path_row(
            "tool-batch-path",
            "/data/figures/",
            folder_btn_id={
                "type": "tools-path-browse",
                "target": "tool-batch-path",
                "mode": "folder",
                "title": "Select batch export folder",
            },
        ),
        dbc.Alert(
            [
                html.I(className="bi bi-info-circle me-1"),
                "Uses active line filters. Muted lines are not exported.",
            ],
            color="info",
            className="tool-mini-alert mb-2",
        ),
        dbc.Row([
            dbc.Col([
                html.Label("Format", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-batch-fmt",
                    options=[
                        {"label": f, "value": f.lower()}
                        for f in ["PNG", "PDF", "SVG", "EPS", "TIFF"]
                    ],
                    value="png", clearable=False,
                ),
            ], width=6),
            dbc.Col([
                html.Label("DPI", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-batch-dpi",
                    options=[
                        {"label": str(d), "value": d}
                        for d in [72, 96, 150, 200, 300, 450, 600]
                    ],
                    value=150, clearable=False,
                ),
            ], width=6),
        ], className="mb-2"),
        html.Label("Include", className="ctrl-label"),
        dbc.Checklist(
            id="tool-batch-items",
            options=[
                {"label": "MT curves (ρ, φ)", "value": "curves"},
                {"label": "ρ pseudosection", "value": "rho_pseudo"},
                {"label": "φ pseudosection", "value": "phi_pseudo"},
                {"label": "Station map", "value": "map"},
                {"label": "Skew / dimensionality QC", "value": "strike"},
            ],
            value=["curves", "rho_pseudo", "phi_pseudo", "map"],
            inline=False,
            className="tool-checklist mb-2",
        ),
        dbc.Checklist(
            id="tool-batch-flags",
            options=[
                {"label": "Overwrite existing files", "value": "overwrite"},
                {"label": "Transparent background", "value": "transparent"},
                {"label": "Write manifest CSV", "value": "manifest"},
            ],
            value=["manifest"],
            className="tool-checklist mb-2",
        ),
        _action_bar(_run_btn("tool-batch-run", "Export All Figures", "images")),
        _out_area("tool-batch-out", min_h=260, cls="mt-2"),
    ])


def _elevation_body(
    n: int, no_data: html.Div, store_data: dict | None = None
) -> html.Div:
    if not n:
        return no_data

    def _lbl(t):
        return html.Div(t, className="ctrl-label")

    records  = (store_data or {}).get("station_records", [])
    line_set = sorted({r.get("Line", "") for r in records if r.get("Line")})
    line_opts = [{"label": ln, "value": ln} for ln in line_set]
    sta_opts  = []
    for r in records:
        sid = r.get("ID") or r.get("station") or r.get("name")
        if sid:
            ln = r.get("Line", "")
            sta_opts.append({"label": f"{sid} · {ln}" if ln else sid, "value": sid})

    # ── Mode pill bar ─────────────────────────────────────────────────────
    mode_bar = html.Div(
        [
            html.Button(
                [html.I(className="bi bi-cloud-download me-1"), "Fetch Online"],
                id="tool-elev-mode-fetch",
                className="elev-mode-btn active",
                n_clicks=0,
            ),
            html.Button(
                [html.I(className="bi bi-upload me-1"), "Upload File"],
                id="tool-elev-mode-upload",
                className="elev-mode-btn",
                n_clicks=0,
            ),
            html.Button(
                [html.I(className="bi bi-sliders me-1"), "Correct & Export"],
                id="tool-elev-mode-correct",
                className="elev-mode-btn",
                n_clicks=0,
            ),
        ],
        className="elev-mode-bar mb-3",
    )

    # ── FETCH PANEL ───────────────────────────────────────────────────────
    fetch_panel = html.Div(
        [
            html.Div(
                [
                    _lbl("DEM Source"),
                    dbc.Select(
                        id="tool-elev-source",
                        options=[
                            {"label": "Open-Elevation  (free, global)",
                             "value": "open-elevation"},
                            {"label": "Open-Meteo  (SRTM 90 m)",
                             "value": "open-meteo"},
                            {"label": "NASADEM  (SRTM 30 m, NASA)",
                             "value": "nasadem"},
                            {"label": "Open-Topo-Data  (multi-DEM)",
                             "value": "open-topo"},
                            {"label": "Session EDI headers  (offline)",
                             "value": "session"},
                        ],
                        value="open-elevation",
                        size="sm",
                    ),
                    html.Div(
                        "Requires internet except 'Session EDI headers'.",
                        className="elev-hint mt-1",
                    ),
                ],
                className="ctrl-card",
            ),
            html.Div(
                [
                    _lbl("Geodetic Datum"),
                    dbc.Select(
                        id="tool-elev-datum",
                        options=[
                            {"label": "WGS84  (EGM96 geoid)", "value": "WGS84"},
                            {"label": "NAD83",                "value": "NAD83"},
                            {"label": "GDA94",                "value": "GDA94"},
                        ],
                        value="WGS84",
                        size="sm",
                    ),
                ],
                className="ctrl-card",
            ),
            html.Div(
                [
                    _lbl("Station Scope"),
                    html.Div(
                        [
                            html.Button(
                                "All Stations",
                                id="tool-elev-scope-all",
                                className="elev-scope-btn active",
                                n_clicks=0,
                            ),
                            html.Button(
                                "Active Lines",
                                id="tool-elev-scope-lines",
                                className="elev-scope-btn",
                                n_clicks=0,
                            ),
                            html.Button(
                                "Selected",
                                id="tool-elev-scope-sel",
                                className="elev-scope-btn",
                                n_clicks=0,
                            ),
                        ],
                        className="elev-scope-bar mt-1 mb-2",
                    ),
                    html.Div(
                        [
                            _lbl("Filter Lines"),
                            dcc.Dropdown(
                                id="tool-elev-lines",
                                options=line_opts,
                                value=None,
                                multi=True,
                                placeholder="All active lines",
                                className="mb-2",
                                optionHeight=24,
                            ),
                            _lbl("Filter Stations"),
                            dcc.Dropdown(
                                id="tool-elev-stations",
                                options=sta_opts,
                                value=None,
                                multi=True,
                                placeholder="All stations in scope",
                                optionHeight=24,
                            ),
                        ],
                        id="tool-elev-scope-filters",
                        style={"display": "none"},
                    ),
                ],
                className="ctrl-card",
            ),
            html.Div(
                [
                    dbc.Button(
                        [
                            html.I(className="bi bi-cloud-download me-1"),
                            "Fetch Elevations",
                        ],
                        id="tool-elev-run",
                        color="primary",
                        size="sm",
                        className="w-100 mb-1",
                        n_clicks=0,
                    ),
                    dbc.Spinner(
                        html.Div(id="tool-elev-fetch-spin"),
                        size="sm",
                        color="primary",
                    ),
                ],
                className="ctrl-card",
            ),
        ],
        id="tool-elev-fetch-panel",
    )

    # ── UPLOAD PANEL ──────────────────────────────────────────────────────
    upload_panel = html.Div(
        [
            html.P(
                "Upload a CSV/TXT with elevation data. Columns are "
                "auto-detected; remap below if needed.",
                className="text-muted small mb-2",
            ),
            dcc.Upload(
                id="tool-elev-upload-input",
                children=html.Div(
                    [
                        html.I(className="bi bi-file-earmark-spreadsheet me-2"),
                        "Drag & drop CSV / TXT  or  ",
                        html.A("browse", style={"cursor": "pointer"}),
                    ]
                ),
                accept=".csv,.txt",
                className="elev-upload-zone mb-2",
            ),
            html.Div(id="tool-elev-upload-info", className="elev-upload-info mb-2"),
            html.Div(
                [
                    html.Div(
                        [
                            _lbl("Station / ID col"),
                            dcc.Dropdown(
                                id="tool-elev-col-id",
                                options=[],
                                value=None,
                                clearable=True,
                                placeholder="auto-detect",
                                className="mb-1",
                            ),
                        ]
                    ),
                    html.Div(
                        [
                            _lbl("Latitude col"),
                            dcc.Dropdown(
                                id="tool-elev-col-lat",
                                options=[],
                                value=None,
                                clearable=True,
                                placeholder="auto-detect",
                                className="mb-1",
                            ),
                        ]
                    ),
                    html.Div(
                        [
                            _lbl("Longitude col"),
                            dcc.Dropdown(
                                id="tool-elev-col-lon",
                                options=[],
                                value=None,
                                clearable=True,
                                placeholder="auto-detect",
                                className="mb-1",
                            ),
                        ]
                    ),
                    html.Div(
                        [
                            _lbl("Elevation col"),
                            dcc.Dropdown(
                                id="tool-elev-col-elev",
                                options=[],
                                value=None,
                                clearable=True,
                                placeholder="auto-detect",
                            ),
                        ]
                    ),
                ],
                id="tool-elev-col-map",
                className="ctrl-card",
            ),
            html.Div(
                [
                    dbc.Button(
                        [
                            html.I(className="bi bi-file-earmark-check me-1"),
                            "Load File",
                        ],
                        id="tool-elev-load-btn",
                        color="secondary",
                        size="sm",
                        className="w-100",
                        n_clicks=0,
                    ),
                ],
                className="ctrl-card",
            ),
        ],
        id="tool-elev-upload-panel",
        style={"display": "none"},
    )

    # ── CORRECT & EXPORT PANEL ────────────────────────────────────────────
    correct_panel = html.Div(
        [
            html.Div(
                [
                    _lbl("Corrections to apply"),
                    dbc.Checklist(
                        id="tool-elev-corr-ops",
                        options=[
                            {"label": "Smooth profile",     "value": "smooth"},
                            {"label": "Shift / Offset",     "value": "shift"},
                            {"label": "Fill missing / NaN", "value": "fill"},
                            {"label": "Remove outliers",    "value": "outlier"},
                        ],
                        value=["smooth"],
                        inline=False,
                        className="mb-2",
                    ),
                    # Smooth params ────────────────────────────────────────
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            _lbl("Algorithm"),
                                            dbc.Select(
                                                id="tool-elev-smooth-method",
                                                options=[
                                                    {
                                                        "label": "LOESS  (local polynomial)",
                                                        "value": "loess",
                                                    },
                                                    {
                                                        "label": "Moving average",
                                                        "value": "mean",
                                                    },
                                                    {
                                                        "label": "Savitzky-Golay  (scipy)",
                                                        "value": "savgol",
                                                    },
                                                    {
                                                        "label": "Gaussian  (scipy)",
                                                        "value": "gaussian",
                                                    },
                                                    {
                                                        "label": "Median filter  (scipy)",
                                                        "value": "median",
                                                    },
                                                ],
                                                value="loess",
                                                size="sm",
                                            ),
                                        ],
                                        style={"flex": "1"},
                                    ),
                                    html.Div(
                                        [
                                            _lbl("Window"),
                                            dbc.Input(
                                                id="tool-elev-smooth-window",
                                                type="number",
                                                value=5,
                                                min=3,
                                                max=31,
                                                step=2,
                                                size="sm",
                                            ),
                                        ],
                                        style={"flex": "0 0 70px"},
                                    ),
                                ],
                                className="d-flex gap-2 align-items-end",
                            ),
                        ],
                        id="tool-elev-smooth-params",
                        className="elev-corr-sub mb-2",
                    ),
                    # Shift params ─────────────────────────────────────────
                    html.Div(
                        [
                            _lbl("Δ Elevation offset (m)"),
                            dbc.Input(
                                id="tool-elev-shift-delta",
                                type="number",
                                value=0.0,
                                step=0.5,
                                size="sm",
                                placeholder="e.g. −12.5",
                            ),
                            html.Div(
                                "Constant added to every station elevation "
                                "(datum shift, GPS antenna height, etc.).",
                                className="elev-hint mt-1",
                            ),
                        ],
                        id="tool-elev-shift-params",
                        className="elev-corr-sub mb-2",
                        style={"display": "none"},
                    ),
                    # Fill params ──────────────────────────────────────────
                    html.Div(
                        [
                            dbc.Checklist(
                                id="tool-elev-fill-zeros",
                                options=[
                                    {
                                        "label": "Also treat zeros as missing",
                                        "value": "zeros",
                                    }
                                ],
                                value=[],
                                inline=False,
                            ),
                            html.Div(
                                "NaN / missing elevations filled by linear "
                                "interpolation from neighbours.",
                                className="elev-hint mt-1",
                            ),
                        ],
                        id="tool-elev-fill-params",
                        className="elev-corr-sub mb-2",
                        style={"display": "none"},
                    ),
                    # Outlier params ───────────────────────────────────────
                    html.Div(
                        [
                            _lbl("σ threshold  (std-dev multiplier)"),
                            dbc.Input(
                                id="tool-elev-outlier-sigma",
                                type="number",
                                value=3.0,
                                min=1.0,
                                max=10.0,
                                step=0.5,
                                size="sm",
                            ),
                            html.Div(
                                "Points beyond median ± σ × std are flagged "
                                "and replaced by interpolation.",
                                className="elev-hint mt-1",
                            ),
                        ],
                        id="tool-elev-outlier-params",
                        className="elev-corr-sub mb-2",
                        style={"display": "none"},
                    ),
                ],
                className="ctrl-card",
            ),
            html.Div(
                [
                    dbc.Button(
                        [html.I(className="bi bi-sliders me-1"), "Run Corrections"],
                        id="tool-elev-correct-btn",
                        color="info",
                        size="sm",
                        outline=True,
                        className="w-100 mb-1",
                        n_clicks=0,
                    ),
                    dbc.Spinner(
                        html.Div(id="tool-elev-correct-spin"),
                        size="sm",
                        color="info",
                    ),
                ],
                className="ctrl-card",
            ),
            html.Hr(className="tool-sep"),
            html.Div(
                [
                    _lbl("Export"),
                    html.Div(
                        [
                            dbc.Button(
                                [
                                    html.I(className="bi bi-filetype-csv me-1"),
                                    "CSV",
                                ],
                                id="tool-elev-export-csv",
                                color="success",
                                size="sm",
                                outline=True,
                                n_clicks=0,
                            ),
                            dbc.Button(
                                [
                                    html.I(
                                        className="bi bi-file-earmark-binary me-1"
                                    ),
                                    "HDF5 / NPZ",
                                ],
                                id="tool-elev-export-h5",
                                color="warning",
                                size="sm",
                                outline=True,
                                n_clicks=0,
                            ),
                            dbc.Button(
                                [
                                    html.I(className="bi bi-arrow-repeat me-1"),
                                    "Update EDI",
                                ],
                                id="tool-elev-update-session",
                                color="primary",
                                size="sm",
                                outline=True,
                                n_clicks=0,
                            ),
                        ],
                        className="d-flex gap-2 flex-wrap mt-1",
                    ),
                    html.Div(
                        id="tool-elev-export-status",
                        className="elev-hint mt-2",
                    ),
                ],
                className="ctrl-card",
            ),
            dcc.Download(id="tool-elev-dl-csv"),
            dcc.Download(id="tool-elev-dl-h5"),
        ],
        id="tool-elev-correct-panel",
        style={"display": "none"},
    )

    # ── Progress ──────────────────────────────────────────────────────────
    progress = html.Div(
        [
            dbc.Progress(
                id="tool-elev-progress",
                value=0,
                max=100,
                striped=True,
                animated=True,
                style={"height": "6px"},
                className="mb-1",
                color="info",
            ),
            html.Div(
                id="tool-elev-progress-label",
                className="elev-progress-label",
            ),
        ],
        id="tool-elev-progress-wrap",
        style={"display": "none"},
        className="mb-2",
    )

    return html.Div(
        [
            html.P(
                [
                    html.I(className="bi bi-layers me-1"),
                    "Elevation Enrichment — fetch, upload, correct, and export "
                    "elevation data for survey stations.",
                ],
                className="text-muted small mb-2",
            ),
            mode_bar,
            fetch_panel,
            upload_panel,
            correct_panel,
            progress,
            _out_area("tool-elev-out", min_h=100, cls="mt-2"),
        ]
    )


def _response_body(n: int, no_data: html.Div, store_data: dict | None, last_output: dict | None = None) -> html.Div:
    if not n:
        return no_data
    records  = (store_data or {}).get("station_records", [])
    sta_opts = [{"label": r["ID"], "value": r["ID"]} for r in records]
    first    = sta_opts[0]["value"] if sta_opts else None

    def _sep(label: str) -> html.Div:
        return html.Div(
            html.Span(label, style={"fontSize": "10px", "color": "var(--sub0)",
                                    "textTransform": "uppercase",
                                    "letterSpacing": "0.05em"}),
            style={"borderBottom": "1px solid var(--overlay0)",
                   "paddingBottom": "2px", "marginBottom": "5px",
                   "marginTop": "8px"},
        )

    return html.Div([
        # ── Station selection ───────────────────────────────────────────
        dbc.Row([
            dbc.Col([
                html.Label("Station", className="ctrl-label"),
                dcc.Dropdown(id="tool-resp-station", options=sta_opts, value=first,
                             clearable=False),
            ], width=6),
            dbc.Col([
                html.Label("Compare", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-resp-station2",
                    options=[{"label": "— none —", "value": ""}] + sta_opts,
                    value="", clearable=False,
                    placeholder="optional 2nd station",
                ),
            ], width=6),
        ], className="g-1 mb-2"),

        # ── Impedance components ─────────────────────────────────────────
        _sep("Impedance"),
        dbc.Checklist(
            id="tool-resp-comps",
            options=[
                {"label": "Zxy (TE)", "value": "xy"},
                {"label": "Zyx (TM)", "value": "yx"},
                {"label": "Zxx",      "value": "xx"},
                {"label": "Zyy",      "value": "yy"},
                {"label": "Det",      "value": "det"},
            ],
            value=["xy", "yx"],
            inline=True,
            inputStyle={"cursor": "pointer"},
            labelStyle={"fontSize": "12px", "cursor": "pointer",
                        "marginRight": "8px"},
            className="mb-1",
        ),

        # ── Tipper ───────────────────────────────────────────────────────
        _sep("Tipper"),
        dbc.Checklist(
            id="tool-resp-tipper",
            options=[
                {"label": "Tzx",        "value": "tx"},
                {"label": "Tzy",        "value": "ty"},
                {"label": "|T| (mag)",  "value": "mag"},
            ],
            value=[],
            inline=True,
            inputStyle={"cursor": "pointer"},
            labelStyle={"fontSize": "12px", "cursor": "pointer",
                        "marginRight": "8px"},
            className="mb-1",
        ),

        # ── Period range ─────────────────────────────────────────────────
        _sep("Period filter (s)"),
        dbc.Row([
            dbc.Col([
                html.Label("Min", className="ctrl-label"),
                dbc.Input(id="tool-resp-period-lo", type="number",
                          placeholder="auto", min=1e-6, step=None,
                          size="sm", debounce=True),
            ], width=6),
            dbc.Col([
                html.Label("Max", className="ctrl-label"),
                dbc.Input(id="tool-resp-period-hi", type="number",
                          placeholder="auto", min=1e-6, step=None,
                          size="sm", debounce=True),
            ], width=6),
        ], className="g-1 mb-2"),

        # ── Plot options ─────────────────────────────────────────────────
        _sep("Display"),
        dbc.Row([
            dbc.Col([
                html.Label("Plot", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-resp-ytype",
                    options=[
                        {"label": "ρ & φ",    "value": "rho_phi"},
                        {"label": "Re / Im Z", "value": "reim"},
                        {"label": "ρ only",   "value": "rho"},
                        {"label": "Phase only","value": "phi"},
                    ],
                    value="rho_phi", clearable=False,
                ),
            ], width=6),
            dbc.Col([
                html.Label("Period axis", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-resp-xscale",
                    options=[{"label": "Log",    "value": "log"},
                             {"label": "Linear", "value": "linear"}],
                    value="log", clearable=False,
                ),
            ], width=6),
        ], className="g-1 mb-2"),
        dbc.Switch(
            id="tool-resp-errbars",
            value=True,
            label="Error bars (from z_err)",
            style={"fontSize": "12px"},
            className="mb-2",
        ),

        _action_bar(_run_btn("tool-resp-run", "Plot Response")),
        _out_area("tool-resp-out", min_h=380, last_output=last_output),
    ])


def _strike_profile_body(n: int, no_data: html.Div, store_data: dict | None, last_output: dict | None = None) -> html.Div:
    if not n:
        return no_data
    line_counts = (store_data or {}).get("line_counts", {})
    line_opts   = [{"label": f"{ln}  ({cnt})", "value": ln}
                   for ln, cnt in line_counts.items()]
    if not line_opts:
        line_opts = [{"label": "All stations", "value": "__all__"}]
    first_line  = line_opts[0]["value"]

    def _sep(label: str) -> html.Div:
        return html.Div(
            html.Span(label, style={"fontSize": "10px", "color": "var(--sub0)",
                                    "textTransform": "uppercase",
                                    "letterSpacing": "0.05em"}),
            style={"borderBottom": "1px solid var(--overlay0)",
                   "paddingBottom": "2px", "marginBottom": "5px",
                   "marginTop": "8px"},
        )

    return html.Div([
        # ── Survey line ─────────────────────────────────────────────────────
        html.Label("Survey line(s)", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-sprof-line",
            options=line_opts,
            value=[first_line],
            multi=True,
            clearable=False,
            className="mb-2",
            placeholder="Select one or more lines…",
        ),

        # ── Method ──────────────────────────────────────────────────────────
        _sep("Estimation method"),
        dbc.RadioItems(
            id="tool-sprof-method",
            options=[
                {"label": "Sweep (tensor diagonal)",     "value": "sweep"},
                {"label": "Phase Tensor  (θ)",           "value": "phase_tensor"},
                {"label": "Consensus  (weighted blend)", "value": "consensus"},
            ],
            value="sweep",
            labelStyle={"fontSize": "12px", "cursor": "pointer",
                        "marginBottom": "2px"},
            className="mb-2",
        ),

        # ── Period band ──────────────────────────────────────────────────────
        _sep("Period band"),
        dcc.Dropdown(
            id="tool-sprof-band",
            options=[
                {"label": "All periods",            "value": "all"},
                {"label": "< 0.01 s  (shallow)",    "value": "shallow"},
                {"label": "0.01 – 1 s  (mid)",      "value": "mid"},
                {"label": "> 1 s  (deep)",           "value": "deep"},
                {"label": "Custom range…",           "value": "custom"},
            ],
            value="all",
            clearable=False,
            className="mb-1",
        ),
        dbc.Row([
            dbc.Col(dbc.Input(id="tool-sprof-period-lo", type="number",
                              placeholder="T min (s)", size="sm", debounce=True), width=6),
            dbc.Col(dbc.Input(id="tool-sprof-period-hi", type="number",
                              placeholder="T max (s)", size="sm", debounce=True), width=6),
        ], id="tool-sprof-custom-band", className="g-1 mb-2",
           style={"display": "none"}),

        # ── View ────────────────────────────────────────────────────────────
        _sep("View & overlays"),
        dbc.Row([
            dbc.Col([
                html.Label("View", className="ctrl-label"),
                dbc.RadioItems(
                    id="tool-sprof-display",
                    options=[
                        {"label": "Profile",  "value": "profile"},
                        {"label": "Heatmap",  "value": "heatmap"},
                    ],
                    value="profile",
                    inline=True,
                    labelStyle={"fontSize": "12px"},
                ),
            ]),
            dbc.Col([
                html.Label("Color by", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-sprof-colorby",
                    options=[
                        {"label": "Strike angle", "value": "angle"},
                        {"label": "IQR (±°)",     "value": "iqr"},
                        {"label": "Swift skew",   "value": "skew"},
                        {"label": "n periods",    "value": "n"},
                    ],
                    value="angle",
                    clearable=False,
                ),
            ]),
        ], className="g-2 mb-2"),
        dbc.Checklist(
            id="tool-sprof-overlays",
            options=[
                {"label": "IQR ribbon",        "value": "iqr"},
                {"label": "Swift skew subplot", "value": "skew"},
                {"label": "Normalise 0 – 90°", "value": "norm"},
                {"label": "Flag 3-D stations", "value": "flag3d"},
                {"label": "Stats table",       "value": "table"},
            ],
            value=["iqr", "table"],
            inline=False,
            inputStyle={"cursor": "pointer"},
            labelStyle={"fontSize": "12px", "cursor": "pointer",
                        "marginRight": "8px"},
            className="mb-2",
        ),
        dbc.Row([
            dbc.Col([
                html.Label("IQR flag (°)", className="ctrl-label"),
                dbc.Input(id="tool-sprof-iqr-thresh", type="number",
                          value=30, min=0, max=90, step=5, size="sm"),
            ], width=6),
            dbc.Col([
                html.Label("Skew flag", className="ctrl-label"),
                dbc.Input(id="tool-sprof-skew-thresh", type="number",
                          value=0.2, min=0.0, max=1.0, step=0.05, size="sm"),
            ], width=6),
        ], className="g-1 mb-2"),

        _action_bar(_run_btn("tool-sprof-run", "Compute Strike")),
        _out_area("tool-sprof-out", min_h=320, last_output=last_output),
    ])


def _phase_tensor_body(n: int, no_data: html.Div, store_data: dict | None = None,
                       last_output: dict | None = None) -> html.Div:
    if not n:
        return no_data

    line_counts = (store_data or {}).get("line_counts", {})
    line_opts   = [{"label": f"{ln}  ({cnt})", "value": ln}
                   for ln, cnt in line_counts.items()]
    if not line_opts:
        line_opts = [{"label": "All stations", "value": "__all__"}]

    def _sep(label: str) -> html.Div:
        return html.Div(
            html.Span(label, style={"fontSize": "10px", "color": "var(--sub0)",
                                    "textTransform": "uppercase",
                                    "letterSpacing": "0.05em"}),
            style={"borderBottom": "1px solid var(--overlay0)",
                   "paddingBottom": "2px", "marginBottom": "5px",
                   "marginTop": "8px"},
        )

    _CMAP_OPTS = [
        {"label": "RdBu_r  (diverging)",  "value": "RdBu_r"},
        {"label": "RdYlBu_r  (warm→cool)","value": "RdYlBu_r"},
        {"label": "coolwarm",             "value": "coolwarm"},
        {"label": "viridis  (sequential)","value": "viridis"},
        {"label": "plasma",               "value": "plasma"},
        {"label": "inferno",              "value": "inferno"},
        {"label": "jet",                  "value": "jet"},
    ]

    return html.Div([
        # ── Lines ──────────────────────────────────────────────────────────
        html.Label("Survey line(s)", className="ctrl-label"),
        dcc.Dropdown(
            id="tool-pt-lines",
            options=line_opts,
            value=[line_opts[0]["value"]],
            multi=True, clearable=False,
            className="mb-2",
        ),

        # ── View ───────────────────────────────────────────────────────────
        _sep("View"),
        dbc.RadioItems(
            id="tool-pt-view",
            options=[
                {"label": "Ellipses  (T × station, Caldwell style)", "value": "ellipses"},
                {"label": "Heatmap  (T × station, colour only)",     "value": "psection"},
                {"label": "Map  (geographic, single period)",        "value": "map"},
                {"label": "All three views",                         "value": "all"},
            ],
            value="ellipses",
            labelStyle={"fontSize": "12px", "cursor": "pointer",
                        "marginBottom": "2px"},
            className="mb-2",
        ),

        # ── Period ─────────────────────────────────────────────────────────
        _sep("Period"),
        dbc.Row([
            dbc.Col([
                html.Label("Target period (s)", className="ctrl-label"),
                dbc.Input(id="tool-pt-period", type="number", value=1.0,
                          min=1e-6, step=None, size="sm", debounce=True),
            ], width=12),
        ], className="mb-1"),
        dbc.Row([
            dbc.Col([
                html.Label("Section T min (s)", className="ctrl-label"),
                dbc.Input(id="tool-pt-period-lo", type="number",
                          placeholder="auto", min=1e-6, size="sm", debounce=True),
            ], width=6),
            dbc.Col([
                html.Label("Section T max (s)", className="ctrl-label"),
                dbc.Input(id="tool-pt-period-hi", type="number",
                          placeholder="auto", min=1e-6, size="sm", debounce=True),
            ], width=6),
        ], className="g-1 mb-2"),

        # ── Colour ─────────────────────────────────────────────────────────
        _sep("Colour"),
        dbc.Row([
            dbc.Col([
                html.Label("Colour by", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-pt-colorby",
                    options=[
                        {"label": "Skew β  (3-D indicator)",  "value": "skew"},
                        {"label": "Strike θ  (CCW from E)",   "value": "theta"},
                        {"label": "Ellipticity  (φ_min/φ_max)","value": "ellipt"},
                        {"label": "φ_max  (mean phase)",      "value": "phi_max"},
                        {"label": "φ_min  (min phase)",       "value": "phi_min"},
                        {"label": "Azimuth α",                "value": "alpha"},
                        {"label": "s1  (major axis)",         "value": "s1"},
                        {"label": "s2  (minor axis)",         "value": "s2"},
                    ],
                    value="skew", clearable=False,
                ),
            ], width=7),
            dbc.Col([
                html.Label("Colourmap", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-pt-cmap",
                    options=_CMAP_OPTS,
                    value="RdBu_r", clearable=False,
                ),
            ], width=5),
        ], className="g-1 mb-2"),
        dbc.Row([
            dbc.Col([
                html.Label("Clim lo (%)", className="ctrl-label"),
                dbc.Input(id="tool-pt-clim-lo", type="number",
                          value=5, min=0, max=50, step=1, size="sm"),
            ], width=6),
            dbc.Col([
                html.Label("Clim hi (%)", className="ctrl-label"),
                dbc.Input(id="tool-pt-clim-hi", type="number",
                          value=95, min=50, max=100, step=1, size="sm"),
            ], width=6),
        ], className="g-1 mb-2"),

        # ── Ellipse ────────────────────────────────────────────────────────
        _sep("Ellipse style"),
        dbc.Row([
            dbc.Col([
                html.Label("Normalise by", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-pt-norm",
                    options=[
                        {"label": "Cell  (relative size)", "value": "cell"},
                        {"label": "Unity  (45° = circle)", "value": "unity"},
                        {"label": "Absolute",              "value": "abs"},
                    ],
                    value="cell", clearable=False,
                ),
            ], width=7),
            dbc.Col([
                html.Label("Scale", className="ctrl-label"),
                dbc.Input(id="tool-pt-scale", type="number",
                          value=0.85, min=0.1, max=2.0, step=0.05, size="sm"),
            ], width=5),
        ], className="g-1 mb-1"),
        dbc.Row([
            dbc.Col([
                html.Label("Ellipse α", className="ctrl-label"),
                dbc.Input(id="tool-pt-alpha", type="number",
                          value=0.92, min=0.1, max=1.0, step=0.05, size="sm"),
            ], width=4),
            dbc.Col([
                html.Label("Edge lw", className="ctrl-label"),
                dbc.Input(id="tool-pt-lw", type="number",
                          value=0.3, min=0.0, max=2.0, step=0.1, size="sm"),
            ], width=4),
            dbc.Col([
                html.Label("Y-axis", className="ctrl-label"),
                dcc.Dropdown(
                    id="tool-pt-yaxis",
                    options=[
                        {"label": "log T ↑", "value": "logperiod_up"},
                        {"label": "log T ↓", "value": "logperiod_dn"},
                        {"label": "log f ↑", "value": "logfreq_up"},
                    ],
                    value="logperiod_up", clearable=False,
                ),
            ], width=4),
        ], className="g-1 mb-1"),
        dbc.Row([
            dbc.Col([
                dbc.Switch(id="tool-pt-ref-ellipse",
                           label="Reference circle",
                           value=True, style={"fontSize": "12px"}),
            ], width=6),
            dbc.Col([
                dbc.Switch(id="tool-pt-sym-clim",
                           label="Symmetric clim",
                           value=True, style={"fontSize": "12px"}),
            ], width=6),
        ], className="mb-2"),

        # ── Overlays ───────────────────────────────────────────────────────
        _sep("Overlays"),
        dbc.Row([
            dbc.Col([
                html.Label("Tipper", className="ctrl-label"),
                dbc.RadioItems(
                    id="tool-pt-tipper",
                    options=[
                        {"label": "Off",       "value": "none"},
                        {"label": "Parkinson", "value": "parkinson"},
                        {"label": "Wiese",     "value": "wiese"},
                    ],
                    value="none",
                    inline=True,
                    labelStyle={"fontSize": "11px"},
                ),
            ]),
        ], className="mb-1"),
        dbc.Row([
            dbc.Col([
                dbc.Switch(id="tool-pt-mark3d",
                           label="Flag 3-D stations",
                           value=True, style={"fontSize": "12px"}),
            ], width=7),
            dbc.Col([
                html.Label("Skew thresh.", className="ctrl-label"),
                dbc.Input(id="tool-pt-skew-thresh", type="number",
                          value=5.0, min=0.0, max=45.0, step=0.5, size="sm"),
            ], width=5),
        ], className="g-1 mb-2"),

        _action_bar(_run_btn("tool-pt-run", "Render Phase Tensor")),
        _out_area("tool-pt-out", min_h=380, last_output=last_output),
    ])


_LM_DEFAULT_TABLE = [
    {"#": 1, "ρ (Ω·m)": 100.0,   "Thickness (m)": 50.0},
    {"#": 2, "ρ (Ω·m)": 500.0,   "Thickness (m)": 200.0},
    {"#": 3, "ρ (Ω·m)": 50.0,    "Thickness (m)": 500.0},
    {"#": 4, "ρ (Ω·m)": 1000.0,  "Thickness (m)": None},
]

_GEOLOGY_PRESETS = [
    "basement", "coastal", "crystalline", "evaporite",
    "geothermal", "hydrothermal", "laterite", "marine",
    "mineralized", "permafrost", "porphyry", "sedimentary", "volcanic",
]


def _layered_model_body(theme: str | None = "dark",
                        last_output: dict | None = None) -> html.Div:
    _S = {"display": "block"}
    _H = {"display": "none"}
    dark = (theme or "dark") == "dark"
    _tbl_cell  = {"padding": "3px 6px", "textAlign": "right",
                  "backgroundColor": "#1e1e2e" if dark else "#e6e9ef",
                  "color": "#cdd6f4" if dark else "#4c4f69",
                  "border": "1px solid #313244" if dark else "1px solid #ccd0da",
                  "fontSize": "0.82rem"}
    _tbl_hdr   = {"backgroundColor": "#181825" if dark else "#dce0e8",
                  "color": "#89b4fa" if dark else "#1e66f5",
                  "fontWeight": "600", "fontSize": "0.82rem"}

    return html.Div([
        html.P(
            "Build a 1-D layered earth, run a MT / CSAMT / TEM forward solver, "
            "optionally add realistic noise, preview the model + synthetic response, "
            "and export for training or comparison.",
            className="text-muted small mb-2",
        ),

        # ═══ 1. MODEL CONSTRUCTOR ════════════════════════════════════════
        html.Label("Model constructor", className="ctrl-label"),
        dbc.RadioItems(
            id="tool-lm-constructor",
            options=[
                {"label": "Manual table",    "value": "manual"},
                {"label": "Random",          "value": "random"},
                {"label": "Blocky anomaly",  "value": "blocky"},
                {"label": "Smooth gradient", "value": "smooth"},
                {"label": "Geology preset",  "value": "geology"},
            ],
            value="manual", inline=True,
            className="mb-2 tool-radio", inputClassName="me-1",
        ),

        # ── Manual table ──────────────────────────────────────────────────
        html.Div([
            html.Label("Layers  (last row = halfspace — leave Thickness blank)",
                       className="ctrl-label"),
            dash_table.DataTable(
                id="tool-lm-table",
                columns=[
                    {"name": "#",             "id": "#",             "editable": False},
                    {"name": "ρ (Ω·m)",      "id": "ρ (Ω·m)",      "editable": True,
                     "type": "numeric"},
                    {"name": "Thickness (m)","id": "Thickness (m)", "editable": True,
                     "type": "numeric"},
                ],
                data=_LM_DEFAULT_TABLE,
                editable=True, row_deletable=True,
                style_table={"overflowX": "auto"},
                style_cell=_tbl_cell,
                style_header=_tbl_hdr,
                style_data_conditional=[{
                    "if": {"state": "active"},
                    "backgroundColor": "#313244" if dark else "#ccd0da",
                    "border": "1px solid #89b4fa",
                }],
            ),
            dbc.Row([
                dbc.Col(
                    html.Button(
                        [html.I(className="bi bi-plus-lg me-1"), "Add layer"],
                        id="tool-lm-add",
                        className="btn btn-sm btn-outline-secondary me-2",
                        n_clicks=0,
                    ),
                    width="auto",
                ),
                dbc.Col(
                    html.Label("Quick preset", className="ctrl-label mt-1"),
                    width="auto", className="align-self-center",
                ),
                dbc.Col(
                    dcc.Dropdown(
                        id="tool-lm-preset",
                        options=[
                            {"label": "Simple 3-layer", "value": "simple"},
                            {"label": "Continental",    "value": "continental"},
                            {"label": "Marine",         "value": "marine"},
                        ],
                        value=None, clearable=True,
                        placeholder="Load…", className="mb-0",
                        style={"minWidth": "160px"},
                    ),
                    width="auto",
                ),
            ], align="center", className="mb-2 g-1"),
        ], id="tool-lm-manual-params", style=_S),

        # ── Random ────────────────────────────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([html.Label("Layers",   className="ctrl-label"),
                         dbc.Input(id="tool-lm-r-nlayers", type="number",
                                   value=5, min=2, max=20, step=1,
                                   debounce=True, className="mb-2")], width=3),
                dbc.Col([html.Label("ρ min (Ω·m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-r-rho-min", type="number",
                                   value=1.0, min=0.01, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("ρ max (Ω·m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-r-rho-max", type="number",
                                   value=10000.0, min=1.0, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("Depth max (m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-r-depth-max", type="number",
                                   value=2000.0, min=10.0, debounce=True,
                                   className="mb-2")], width=3),
            ]),
            dbc.Row([
                dbc.Col([html.Label("Seed (blank = random)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-r-seed", type="number",
                                   placeholder="auto", debounce=True,
                                   className="mb-2")], width=4),
            ]),
        ], id="tool-lm-random-params", style=_H),

        # ── Blocky ────────────────────────────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([html.Label("Layers", className="ctrl-label"),
                         dbc.Input(id="tool-lm-b-nlayers", type="number",
                                   value=4, min=2, max=20, step=1,
                                   debounce=True, className="mb-2")], width=3),
                dbc.Col([html.Label("ρ background (Ω·m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-b-rho-bg", type="number",
                                   value=100.0, min=0.01, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("ρ anomaly (Ω·m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-b-rho-anom", type="number",
                                   value=5.0, min=0.01, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("Anomaly layer #", className="ctrl-label"),
                         dbc.Input(id="tool-lm-b-anom-layer", type="number",
                                   value=1, min=0, step=1, debounce=True,
                                   className="mb-2")], width=3),
            ]),
            dbc.Row([
                dbc.Col([html.Label("Depth max (m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-b-depth-max", type="number",
                                   value=1000.0, min=10.0, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("Equal thickness", className="ctrl-label mt-3"),
                         dbc.Switch(id="tool-lm-b-equal-th", value=True)], width=3),
                dbc.Col([html.Label("Seed", className="ctrl-label"),
                         dbc.Input(id="tool-lm-b-seed", type="number",
                                   placeholder="auto", debounce=True,
                                   className="mb-2")], width=3),
            ]),
        ], id="tool-lm-blocky-params", style=_H),

        # ── Smooth ────────────────────────────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([html.Label("Layers", className="ctrl-label"),
                         dbc.Input(id="tool-lm-s-nlayers", type="number",
                                   value=10, min=3, max=50, step=1,
                                   debounce=True, className="mb-2")], width=3),
                dbc.Col([html.Label("ρ surface (Ω·m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-s-rho-surf", type="number",
                                   value=100.0, min=0.01, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("ρ deep (Ω·m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-s-rho-deep", type="number",
                                   value=10.0, min=0.01, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("Depth max (m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-s-depth-max", type="number",
                                   value=5000.0, min=10.0, debounce=True,
                                   className="mb-2")], width=3),
            ]),
            dbc.Row([
                dbc.Col([html.Label("Perturbation (0–1)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-s-perturb", type="number",
                                   value=0.2, min=0.0, max=1.0, step=0.05,
                                   debounce=True, className="mb-2")], width=3),
                dbc.Col([html.Label("Seed", className="ctrl-label"),
                         dbc.Input(id="tool-lm-s-seed", type="number",
                                   placeholder="auto", debounce=True,
                                   className="mb-2")], width=3),
            ]),
        ], id="tool-lm-smooth-params", style=_H),

        # ── Geology preset ────────────────────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([
                    html.Label("Geology scenario", className="ctrl-label"),
                    dcc.Dropdown(
                        id="tool-lm-geology",
                        options=[{"label": g.capitalize(), "value": g}
                                 for g in _GEOLOGY_PRESETS],
                        value="sedimentary", clearable=False,
                        className="mb-2",
                    ),
                ], width=6),
                dbc.Col([
                    html.Label("Seed", className="ctrl-label"),
                    dbc.Input(id="tool-lm-geo-seed", type="number",
                              placeholder="auto", debounce=True,
                              className="mb-2"),
                ], width=3),
            ]),
        ], id="tool-lm-geology-params", style=_H),

        html.Hr(style={"borderColor": "#45475a", "margin": "6px 0"}),

        # ═══ 2. FORWARD SOLVER ═══════════════════════════════════════════
        html.Label("Forward solver", className="ctrl-label"),
        dbc.RadioItems(
            id="tool-lm-solver",
            options=[
                {"label": "MT 1D",     "value": "mt1d"},
                {"label": "CSAMT 1D",  "value": "csamt1d"},
                {"label": "TEM 1D",    "value": "tem1d"},
            ],
            value="mt1d", inline=True,
            className="mb-2 tool-radio", inputClassName="me-1",
        ),
        # MT / CSAMT frequency grid
        html.Div([
            dbc.Row([
                dbc.Col([html.Label("f min (Hz)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-freq-min", type="number",
                                   value=0.001, min=1e-6, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("f max (Hz)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-freq-max", type="number",
                                   value=10000.0, min=1e-4, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("N frequencies", className="ctrl-label"),
                         dbc.Input(id="tool-lm-n-freqs", type="number",
                                   value=30, min=5, max=200, step=1,
                                   debounce=True, className="mb-2")], width=3),
                # CSAMT extra (hidden by default)
                dbc.Col(html.Div([
                    html.Label("Source offset (m)", className="ctrl-label"),
                    dbc.Input(id="tool-lm-src-offset", type="number",
                              value=1000.0, min=1.0, debounce=True,
                              className="mb-2"),
                ], id="tool-lm-csamt-offset", style=_H), width=3),
            ]),
        ], id="tool-lm-freq-params", style=_S),
        # TEM time grid
        html.Div([
            dbc.Row([
                dbc.Col([html.Label("t min (s)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-time-min", type="number",
                                   value=1e-6, min=1e-9, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("t max (s)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-time-max", type="number",
                                   value=0.01, min=1e-6, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("N time gates", className="ctrl-label"),
                         dbc.Input(id="tool-lm-n-times", type="number",
                                   value=25, min=5, max=200, step=1,
                                   debounce=True, className="mb-2")], width=3),
                dbc.Col([html.Label("Loop radius (m)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-loop-radius", type="number",
                                   value=50.0, min=1.0, debounce=True,
                                   className="mb-2")], width=3),
            ]),
            dbc.Row([
                dbc.Col([html.Label("Dipole moment (A·m²)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-moment", type="number",
                                   value=1.0, min=0.01, debounce=True,
                                   className="mb-2")], width=3),
            ]),
        ], id="tool-lm-tem-params", style=_H),

        html.Hr(style={"borderColor": "#45475a", "margin": "6px 0"}),

        # ═══ 3. NOISE MODEL ══════════════════════════════════════════════
        html.Label("Noise model", className="ctrl-label"),
        dbc.RadioItems(
            id="tool-lm-noise",
            options=[
                {"label": "None",              "value": "none"},
                {"label": "Gaussian",          "value": "gaussian"},
                {"label": "Field-realistic",   "value": "field"},
                {"label": "Multiplicative",    "value": "mult"},
            ],
            value="none", inline=True,
            className="mb-2 tool-radio", inputClassName="me-1",
        ),
        html.Div([
            dbc.Row([
                dbc.Col([html.Label("Noise level (%)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-noise-level", type="number",
                                   value=5.0, min=0.1, max=50.0, step=0.5,
                                   debounce=True, className="mb-2")], width=4),
                dbc.Col([html.Label("Phase noise level (°, blank = same)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-noise-phase", type="number",
                                   placeholder="same as level", min=0.1,
                                   debounce=True, className="mb-2")], width=4),
                dbc.Col([html.Label("Seed", className="ctrl-label"),
                         dbc.Input(id="tool-lm-noise-seed", type="number",
                                   placeholder="auto", debounce=True,
                                   className="mb-2")], width=4),
            ]),
        ], id="tool-lm-gauss-params", style=_H),
        html.Div([
            dbc.Row([
                dbc.Col([html.Label("Base noise (%)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-field-base", type="number",
                                   value=2.0, min=0.1, max=50.0, step=0.5,
                                   debounce=True, className="mb-2")], width=3),
                dbc.Col([html.Label("Powerline freq (Hz)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-field-plfreq", type="number",
                                   value=50.0, min=1.0, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("Powerline level (%)", className="ctrl-label"),
                         dbc.Input(id="tool-lm-field-pllevel", type="number",
                                   value=30.0, min=0.0, debounce=True,
                                   className="mb-2")], width=3),
                dbc.Col([html.Label("MT dead-band", className="ctrl-label mt-3"),
                         dbc.Switch(id="tool-lm-field-dead-band", value=True)], width=3),
            ]),
        ], id="tool-lm-field-params", style=_H),
        html.Div([
            dbc.Row([
                dbc.Col([html.Label("σ log₁₀ρ", className="ctrl-label"),
                         dbc.Input(id="tool-lm-mult-sigma", type="number",
                                   value=0.05, min=0.001, max=1.0, step=0.005,
                                   debounce=True, className="mb-2")], width=4),
            ]),
        ], id="tool-lm-mult-params", style=_H),

        html.Hr(style={"borderColor": "#45475a", "margin": "6px 0"}),

        # ═══ 4. DISPLAY OPTIONS ══════════════════════════════════════════
        dbc.Row([
            dbc.Col([
                html.Label("View", className="ctrl-label"),
                dbc.RadioItems(
                    id="tool-lm-view",
                    options=[
                        {"label": "Model + response", "value": "both"},
                        {"label": "Model only",       "value": "model"},
                        {"label": "Response only",    "value": "response"},
                    ],
                    value="both", inline=True,
                    className="mb-2 tool-radio", inputClassName="me-1",
                ),
            ], width=7),
            dbc.Col([
                html.Label("Depth max plot (m)", className="ctrl-label"),
                dbc.Input(id="tool-lm-depth-plot", type="number",
                          placeholder="auto", min=1.0, debounce=True,
                          className="mb-2"),
            ], width=3),
            dbc.Col([
                html.Label("Log ρ axis", className="ctrl-label mt-3"),
                dbc.Switch(id="tool-lm-log-rho", value=True),
            ], width=2),
        ], className="tool-control-row"),

        html.Hr(style={"borderColor": "#45475a", "margin": "6px 0"}),

        # ═══ 5. EXPORT ═══════════════════════════════════════════════════
        dbc.Row([
            dbc.Col([
                html.Label("Export", className="ctrl-label"),
                dbc.RadioItems(
                    id="tool-lm-export-fmt",
                    options=[
                        {"label": "Model CSV",    "value": "model_csv"},
                        {"label": "Response CSV", "value": "resp_csv"},
                        {"label": "JSON (both)",  "value": "json"},
                    ],
                    value="model_csv", inline=True,
                    className="mb-2 tool-radio", inputClassName="me-1",
                ),
            ], width=9),
            dbc.Col(
                dbc.Button(
                    [html.I(className="bi bi-download me-1"), "Export"],
                    id="tool-lm-export-btn",
                    color="secondary", size="sm", outline=True,
                    disabled=True, className="mt-3",
                ),
                width=3,
            ),
        ], align="end", className="tool-control-row"),
        html.Div(id="tool-lm-export-out", className="small mb-2"),
        dcc.Download(id="tool-lm-download"),

        _action_bar(_run_btn("tool-lm-run", "Run Forward", "play-fill")),
        _out_area("tool-lm-out", min_h=300, last_output=last_output),
    ])


def _placeholder_body(tool_id: str) -> html.Div:
    _, _, label, _ = _TOOL_BY_ID.get(tool_id, (tool_id, "tools", tool_id, ""))
    return html.Div([
        html.Div([
            html.I(className="bi bi-wrench-adjustable",
                   style={"fontSize": "2rem", "color": "var(--blue)", "marginBottom": "8px"}),
            html.P(f"{label} is available in the desktop application.",
                   className="text-muted small mb-1"),
            html.P("Web integration coming in a future release.",
                   className="text-muted small"),
        ], className="text-center py-4"),
    ], className="tool-placeholder")


# ── Callbacks ─────────────────────────────────────────────────────────────────

def register_tools(app) -> None:

    # 1. Open offcanvas when a navbar dropdown item or picker rail button is clicked
    @app.callback(
        Output(IDs.TOOLS_CANVAS, "is_open"),
        Output(IDs.ACTIVE_TOOL,  "data"),
        Input({"type": "tool-menu-item",  "index": ALL}, "n_clicks"),
        Input({"type": "tool-picker-btn", "index": ALL}, "n_clicks"),
        State(IDs.TOOLS_CANVAS, "is_open"),
        State(IDs.ACTIVE_TOOL,  "data"),
        prevent_initial_call=True,
    )
    def _open_tools(menu_clicks, picker_clicks, is_open, current_tool):
        triggered = ctx.triggered
        if not triggered:
            return no_update, no_update
        if not any(t["value"] for t in triggered if t["value"]):
            return no_update, no_update
        trigger_id = ctx.triggered_id
        if trigger_id is None:
            return no_update, no_update
        new_tool = (trigger_id.get("index", current_tool)
                    if isinstance(trigger_id, dict) else current_tool)
        return True, new_tool

    # 2. Render tool panel (header + body) when active tool or loaded data changes.
    #    STORE_ACTIVE_LINES and STORE_THEME are State (not Input) so that navigating
    #    to the Lines panel or toggling the theme does NOT rebuild the body and
    #    wipe the user's last run result.  tools-run-store is also State so a
    #    strike run no longer triggers a spurious rebuild.
    @app.callback(
        Output("tools-panel-content", "children"),
        Output(IDs.TOOLS_CANVAS, "title"),
        Input(IDs.ACTIVE_TOOL,        "data"),
        Input(IDs.STORE_DATA,         "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State("tools-run-store",      "data"),
        State(IDs.STORE_THEME,        "data"),
        State("tool-outputs-store",   "data"),
        prevent_initial_call=False,
    )
    def _render_tool(tool_id, store_data, active_lines_store, run_store, theme, saved_outputs):
        tool_id = tool_id or "strike"
        last_output = (saved_outputs or {}).get(tool_id)
        body = _build_tool_body(
            tool_id, store_data, active_lines_store, run_store, theme,
            last_output=last_output,
        )
        _, icon, label, _ = _TOOL_BY_ID.get(
            tool_id, (tool_id, "tools", "Tools", "")
        )
        title = html.Span(
            [
                _icon("tools", size=16, cls="menu-item-icon"),
                html.Span("Tools", className="tool-title-root"),
                html.Span("·", className="tool-title-sep"),
                _icon(icon, size=15, cls="menu-item-icon"),
                html.Span(label, className="tool-title-active"),
            ],
            className="tools-title-smart",
        )
        return _tool_panel(tool_id, body=body), title

    # 3. Highlight active picker button
    @app.callback(
        Output({"type": "tool-picker-btn", "index": ALL}, "className"),
        Input(IDs.ACTIVE_TOOL, "data"),
        prevent_initial_call=False,
    )
    def _highlight_picker(active):
        active = active or "strike"
        ids = [t[0] for t in _TOOL_REGISTRY]
        return [
            "tool-picker-btn active" if tid == active else "tool-picker-btn"
            for tid in ids
        ]

    def _path_file_allowed(name: str) -> bool:
        suffixes = (
            ".avg", ".j", ".edi", ".stn", ".csv", ".txt", ".dat",
            ".json", ".yaml", ".yml", ".z", ".log",
        )
        return name.lower().endswith(suffixes)

    def _path_listing(path: str, mode: str) -> list:
        path = os.path.abspath(os.path.expanduser(path or "~"))
        try:
            entries = sorted(
                os.scandir(path),
                key=lambda e: (not e.is_dir(follow_symlinks=False), e.name.lower()),
            )
        except Exception as exc:
            return [
                html.Div(
                    f"Cannot read this folder: {exc}",
                    className="tools-path-empty",
                )
            ]

        rows = []
        for entry in entries:
            if entry.name.startswith("."):
                continue
            try:
                is_dir = entry.is_dir(follow_symlinks=False)
                is_file = entry.is_file(follow_symlinks=False)
            except OSError:
                continue
            if is_dir:
                rows.append(
                    dbc.ListGroupItem(
                        [
                            html.I(className="bi bi-folder-fill tools-path-folder me-2"),
                            html.Span(entry.name),
                        ],
                        id={"type": "tools-path-dir", "path": entry.path},
                        action=True,
                        n_clicks=0,
                        className="tools-path-item",
                    )
                )
            elif mode in ("file", "any") and is_file and _path_file_allowed(entry.name):
                rows.append(
                    dbc.ListGroupItem(
                        [
                            html.I(className="bi bi-file-earmark-text tools-path-file me-2"),
                            html.Span(entry.name),
                        ],
                        id={"type": "tools-path-file", "path": entry.path},
                        action=True,
                        n_clicks=0,
                        className="tools-path-item",
                    )
                )

        if not rows:
            text = (
                "No folders or supported files here."
                if mode in ("file", "any")
                else "No sub-folders here."
            )
            return [html.Div(text, className="tools-path-empty")]
        return [dbc.ListGroup(rows, flush=True)]

    @app.callback(
        Output("tools-path-modal", "is_open"),
        Output("tools-path-current", "data"),
        Output("tools-path-target", "data"),
        Output("tools-path-title", "children"),
        Input({"type": "tools-path-browse", "target": ALL, "mode": ALL, "title": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def _open_tools_path_browser(_clicks):
        if not any(_clicks or []):
            raise PreventUpdate
        trigger = ctx.triggered_id
        if not isinstance(trigger, dict):
            raise PreventUpdate
        target = {
            "target": trigger.get("target"),
            "mode": trigger.get("mode", "folder"),
            "title": trigger.get("title", "Browse Server Path"),
        }
        return True, os.path.expanduser("~"), target, target["title"]

    @app.callback(
        Output("tools-path-list", "children"),
        Output("tools-path-cwd", "children"),
        Input("tools-path-current", "data"),
        State("tools-path-target", "data"),
        prevent_initial_call=False,
    )
    def _refresh_tools_path_listing(path, target):
        mode = (target or {}).get("mode", "folder")
        path = os.path.abspath(os.path.expanduser(path or "~"))
        return _path_listing(path, mode), path

    @app.callback(
        Output("tools-path-current", "data", allow_duplicate=True),
        Input({"type": "tools-path-dir", "path": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def _tools_path_enter_dir(_clicks):
        if not any(_clicks or []):
            raise PreventUpdate
        trigger = ctx.triggered_id
        if isinstance(trigger, dict) and trigger.get("path"):
            return trigger["path"]
        raise PreventUpdate

    @app.callback(
        Output("tools-path-current", "data", allow_duplicate=True),
        Input("tools-path-up", "n_clicks"),
        State("tools-path-current", "data"),
        prevent_initial_call=True,
    )
    def _tools_path_parent(n_clicks, current):
        if not n_clicks:
            raise PreventUpdate
        current = os.path.abspath(os.path.expanduser(current or "~"))
        parent = os.path.dirname(current)
        if parent and parent != current:
            return parent
        raise PreventUpdate

    @app.callback(
        Output("tools-path-current", "data", allow_duplicate=True),
        Output("tools-path-new", "value"),
        Input("tools-path-mkdir", "n_clicks"),
        State("tools-path-new", "value"),
        State("tools-path-current", "data"),
        prevent_initial_call=True,
    )
    def _tools_path_mkdir(n_clicks, name, current):
        if not n_clicks or not name:
            raise PreventUpdate
        clean = str(name).strip().replace("/", "_").replace("\\", "_")
        if not clean:
            raise PreventUpdate
        current = os.path.abspath(os.path.expanduser(current or "~"))
        new_path = os.path.join(current, clean)
        os.makedirs(new_path, exist_ok=True)
        return new_path, ""

    @app.callback(
        Output("tools-path-modal", "is_open", allow_duplicate=True),
        Input("tools-path-cancel", "n_clicks"),
        prevent_initial_call=True,
    )
    def _tools_path_cancel(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        return False

    @app.callback(
        Output("tools-path-modal", "is_open", allow_duplicate=True),
        Input("tools-path-select-folder", "n_clicks"),
        State("tools-path-current", "data"),
        State("tools-path-target", "data"),
        prevent_initial_call=True,
    )
    def _tools_path_select_folder(n_clicks, current, target):
        if not n_clicks or not target:
            raise PreventUpdate
        target_id = target.get("target")
        if not target_id:
            raise PreventUpdate
        path = os.path.abspath(os.path.expanduser(current or "~"))
        set_props(target_id, {"value": path})
        return False

    @app.callback(
        Output("tools-path-modal", "is_open", allow_duplicate=True),
        Input({"type": "tools-path-file", "path": ALL}, "n_clicks"),
        State("tools-path-target", "data"),
        prevent_initial_call=True,
    )
    def _tools_path_select_file(_clicks, target):
        if not any(_clicks or []):
            raise PreventUpdate
        trigger = ctx.triggered_id
        if not isinstance(trigger, dict) or not target:
            raise PreventUpdate
        target_id = target.get("target")
        path = trigger.get("path")
        if not target_id or not path:
            raise PreventUpdate
        set_props(target_id, {"value": os.path.abspath(os.path.expanduser(path))})
        return False

    # ── Tool run callbacks ─────────────────────────────────────────────────────

    @app.callback(
        Output("tool-strike-stations", "options"),
        Output("tool-strike-stations", "value"),
        Input("tool-strike-lines", "value"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def _sync_strike_station_scope(selected_lines, store_data, active_lines_store):
        opts = _station_options_for_lines(
            store_data, active_lines_store, selected_lines
        )
        return opts, None

    @app.callback(
        Output("tool-valid-stations", "options"),
        Output("tool-valid-stations", "value"),
        Input("tool-valid-lines", "value"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def _sync_validator_station_scope(selected_lines, store_data, active_lines_store):
        opts = _station_options_for_lines(
            store_data, active_lines_store, selected_lines
        )
        return opts, None

    @app.callback(
        Output("tool-conv-avg-options", "style"),
        Output("tool-conv-spectra-options", "style"),
        Input("tool-conv-type", "value"),
        prevent_initial_call=False,
    )
    def _sync_converter_options(conv_type):
        text = str(conv_type or "")
        show_avg = {"display": "block"} if text.startswith("AVG") else {"display": "none"}
        show_spec = (
            {"display": "block"}
            if text.startswith("Spectra")
            else {"display": "none"}
        )
        return show_avg, show_spec

    # ── Coordinate transformer — clientside mode / direction toggles ─────────

    app.clientside_callback(
        """
        function(n0, n1, n2, cur) {
            const trig = dash_clientside.callback_context.triggered_id;
            let mode = cur || 'single';
            if (trig === 'tool-coord-mode-single') mode = 'single';
            if (trig === 'tool-coord-mode-batch')  mode = 'batch';
            if (trig === 'tool-coord-mode-survey') mode = 'survey';
            const show = m => mode === m ? {} : {display:'none'};
            const btn  = m => mode === m
                ? 'coord-mode-btn active' : 'coord-mode-btn';
            return [mode,
                    show('batch'), show('survey'),
                    btn('single'), btn('batch'), btn('survey')];
        }
        """,
        Output("tool-coord-mode-store",   "data"),
        Output("tool-coord-batch-panel",  "style"),
        Output("tool-coord-survey-panel", "style"),
        Output("tool-coord-mode-single",  "className"),
        Output("tool-coord-mode-batch",   "className"),
        Output("tool-coord-mode-survey",  "className"),
        Input("tool-coord-mode-single",   "n_clicks"),
        Input("tool-coord-mode-batch",    "n_clicks"),
        Input("tool-coord-mode-survey",   "n_clicks"),
        State("tool-coord-mode-store",    "data"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        """
        function(dir) {
            const utm2ll = dir === 'utm2ll';
            const ll2utm = dir === 'll2utm';
            const epsg   = dir === 'epsg2epsg';
            const show = v => v ? {} : {display:'none'};
            return [show(utm2ll), show(ll2utm), show(epsg)];
        }
        """,
        Output("tool-coord-single-utm",  "style"),
        Output("tool-coord-single-ll",   "style"),
        Output("tool-coord-epsg-panel",  "style"),
        Input("tool-coord-direction",    "value"),
        prevent_initial_call=False,
    )

    # ── Parse uploaded CSV → column options ───────────────────────────────────
    @app.callback(
        Output("tool-coord-upload-store", "data"),
        Output("tool-coord-upload-info",  "children"),
        Output("tool-coord-col-e",        "options"),
        Output("tool-coord-col-n",        "options"),
        Output("tool-coord-col-id",       "options"),
        Input("tool-coord-upload",        "contents"),
        State("tool-coord-upload",        "filename"),
        prevent_initial_call=True,
    )
    def _coord_parse_csv(contents, filename):
        if not contents:
            return no_update, no_update, no_update, no_update, no_update
        import base64

        import pandas as pd
        try:
            _, encoded = contents.split(",", 1)
            raw = base64.b64decode(encoded).decode("utf-8", errors="replace")
            sep = "," if "," in raw.split("\n")[0] else "\t" if "\t" in raw.split("\n")[0] else None
            df = pd.read_csv(io.StringIO(raw), sep=sep, engine="python", nrows=5000)
            cols = list(df.columns)
            opts = [{"label": c, "value": c} for c in cols]
            id_opts = [{"label": "— none —", "value": ""}] + opts
            # auto-detect columns
            def _guess(keywords):
                for c in cols:
                    if any(k in c.lower() for k in keywords):
                        return c
                return cols[0] if cols else ""
            e_guess = _guess(["east", "easting", "x", "lon", "e"])
            n_guess = _guess(["north", "northing", "y", "lat", "n"])
            info = html.Div(
                [html.I(className="bi bi-check-circle-fill text-success me-1"),
                 f"{filename} — {len(df)} rows, {len(cols)} columns"],
                className="small",
            )
            # store as list-of-dicts (JSON-serialisable)
            return df.to_dict("records"), info, \
                   [{"label": c, "value": c, "selected": c == e_guess} for c in cols], \
                   [{"label": c, "value": c, "selected": c == n_guess} for c in cols], \
                   id_opts
        except Exception as exc:
            return None, html.Span(f"✗ {exc}", className="text-danger small"), [], [], []

    # ── Main convert callback ─────────────────────────────────────────────────
    @app.callback(
        Output("tool-coord-out",           "children"),
        Output("tool-coord-result-store",  "data"),
        Output("tool-coord-progress",      "value"),
        Output("tool-coord-progress-label","children"),
        Output("tool-coord-progress-wrap", "style"),
        Output("tool-coord-export-btn",    "disabled"),
        Input("tool-coord-run",            "n_clicks"),
        State("tool-coord-mode-store",     "data"),
        State("tool-coord-direction",      "value"),
        State("tool-coord-datum",          "value"),
        # single UTM→LL
        State("tool-coord-easting",        "value"),
        State("tool-coord-northing",       "value"),
        State("tool-coord-zone",           "value"),
        # single LL→UTM
        State("tool-coord-lat",            "value"),
        State("tool-coord-lon",            "value"),
        State("tool-coord-target-zone",    "value"),
        # EPSG
        State("tool-coord-epsg-src",       "value"),
        State("tool-coord-epsg-dst",       "value"),
        # batch
        State("tool-coord-upload-store",   "data"),
        State("tool-coord-col-e",          "value"),
        State("tool-coord-col-n",          "value"),
        State("tool-coord-col-id",         "value"),
        State("tool-coord-zone-batch",     "value"),
        # survey
        State("tool-coord-zone-survey",    "value"),
        State(IDs.STORE_DATA,              "data"),
        prevent_initial_call=True,
    )
    def _run_coords(
        n, mode, direction, datum,
        easting, northing, zone,
        lat, lon, target_zone,
        epsg_src, epsg_dst,
        csv_data, col_e, col_n, col_id, zone_batch,
        zone_survey, survey_store,
    ):
        if not n:
            return (no_update,) * 6

        hide_prog  = {"display": "none"}
        show_prog  = {}

        import pandas as pd

        from pycsamt.gis.utils import (
            project_point_ll2utm,
            project_point_utm2ll,
        )

        # ── helpers ──────────────────────────────────────────────────────────
        def _result_table(rows: list[dict]) -> html.Div:
            if not rows:
                return html.Div("No results.", className="text-muted small")
            df_out = pd.DataFrame(rows)
            return html.Div([
                dash_table.DataTable(
                    data=df_out.to_dict("records"),
                    columns=[{"name": c, "id": c} for c in df_out.columns],
                    page_size=10,
                    style_table={"overflowX": "auto", "fontSize": "11px"},
                    style_cell={"padding": "4px 8px", "textAlign": "left",
                                "border": "1px solid var(--srf0)",
                                "background": "var(--crust)", "color": "var(--txt)"},
                    style_header={"fontWeight": "700", "background": "var(--srf0)",
                                  "color": "var(--ov0)", "border": "1px solid var(--srf1)"},
                    style_data_conditional=[
                        {"if": {"row_index": "odd"},
                         "background": "var(--mantle)"},
                    ],
                ),
                html.Div(f"✓ {len(rows)} points converted.",
                         className="small text-success mt-1"),
            ])

        def _single_result_card(label_vals: list[tuple]) -> html.Div:
            items = []
            for lbl, val in label_vals:
                items.append(html.Div([
                    html.Span(lbl, className="coord-result-label"),
                    html.Span(val, className="coord-result-value"),
                ], className="coord-result-row"))
            return html.Div(items, className="coord-result-card mt-2")

        try:
            # ── SINGLE MODE ──────────────────────────────────────────────────
            if mode == "single":
                if direction == "utm2ll":
                    if None in (easting, northing) or not zone:
                        return (_warn("Fill in Easting, Northing and UTM Zone."),
                                no_update, 0, "", hide_prog, True)
                    lat_r, lon_r = project_point_utm2ll(
                        float(easting), float(northing),
                        zone.strip(), datum=datum or "WGS84",
                    )
                    card = _single_result_card([
                        ("Latitude",  f"{lat_r:.8f} °"),
                        ("Longitude", f"{lon_r:.8f} °"),
                        ("DMS Lat",   _dd2dms(lat_r, "lat")),
                        ("DMS Lon",   _dd2dms(lon_r, "lon")),
                    ])
                    result = [{"Easting(m)": easting, "Northing(m)": northing,
                               "Zone": zone, "Latitude": round(lat_r, 8),
                               "Longitude": round(lon_r, 8)}]
                    return card, result, 100, "Done", show_prog, False

                elif direction == "ll2utm":
                    if None in (lat, lon):
                        return (_warn("Fill in Latitude and Longitude."),
                                no_update, 0, "", hide_prog, True)
                    tz = (target_zone or "").strip() or None
                    e_r, n_r, zone_r = project_point_ll2utm(
                        float(lat), float(lon),
                        datum=datum or "WGS84",
                        utm_zone=tz,
                    )
                    card = _single_result_card([
                        ("Easting (m)",  f"{e_r:.2f}"),
                        ("Northing (m)", f"{n_r:.2f}"),
                        ("UTM Zone",     str(zone_r or "auto")),
                    ])
                    result = [{"Latitude": lat, "Longitude": lon,
                               "Easting(m)": round(e_r, 2),
                               "Northing(m)": round(n_r, 2),
                               "Zone": zone_r}]
                    return card, result, 100, "Done", show_prog, False

                elif direction == "epsg2epsg":
                    if None in (easting, northing, epsg_src, epsg_dst):
                        return (_warn("Fill in X, Y and both EPSG codes."),
                                no_update, 0, "", hide_prog, True)
                    try:
                        from pyproj import Transformer
                        tr = Transformer.from_crs(
                            int(epsg_src), int(epsg_dst), always_xy=True)
                        xo, yo = tr.transform(float(easting), float(northing))
                        card = _single_result_card([
                            ("X out", f"{xo:.6f}"),
                            ("Y out", f"{yo:.6f}"),
                            ("CRS",   f"EPSG:{epsg_src} → EPSG:{epsg_dst}"),
                        ])
                        result = [{"X_in": easting, "Y_in": northing,
                                   "EPSG_src": epsg_src, "EPSG_dst": epsg_dst,
                                   "X_out": round(xo, 6), "Y_out": round(yo, 6)}]
                        return card, result, 100, "Done", show_prog, False
                    except ImportError:
                        return (_warn("pyproj is required for EPSG→EPSG conversion."),
                                no_update, 0, "", hide_prog, True)

            # ── BATCH CSV MODE ───────────────────────────────────────────────
            elif mode == "batch":
                if not csv_data:
                    return (_warn("Upload a CSV file first."),
                            no_update, 0, "", hide_prog, True)
                if not col_e or not col_n:
                    return (_warn("Select the Easting and Northing columns."),
                            no_update, 0, "", hide_prog, True)
                df = pd.DataFrame(csv_data)
                if col_e not in df.columns or col_n not in df.columns:
                    return (_warn(f"Columns '{col_e}' / '{col_n}' not found in data."),
                            no_update, 0, "", hide_prog, True)

                rows_out = []
                total = len(df)
                z_batch = (zone_batch or "").strip() or None

                for i, row in df.iterrows():
                    try:
                        e_v = float(row[col_e])
                        n_v = float(row[col_n])
                        if direction == "utm2ll":
                            zo = z_batch or str(row.get("zone", row.get("Zone", "49N")))
                            la, lo = project_point_utm2ll(
                                e_v, n_v, zo, datum=datum or "WGS84")
                            rec = {"Easting(m)": e_v, "Northing(m)": n_v,
                                   "Zone": zo, "Latitude": round(la, 8),
                                   "Longitude": round(lo, 8)}
                        elif direction == "ll2utm":
                            e_r, n_r, zone_r = project_point_ll2utm(
                                e_v, n_v, datum=datum or "WGS84",
                                utm_zone=z_batch)
                            rec = {"Lat_in": e_v, "Lon_in": n_v,
                                   "Easting(m)": round(e_r, 2),
                                   "Northing(m)": round(n_r, 2),
                                   "Zone": zone_r}
                        elif direction == "epsg2epsg":
                            if not epsg_src or not epsg_dst:
                                raise ValueError("Specify source and target EPSG.")
                            from pyproj import Transformer
                            tr = Transformer.from_crs(
                                int(epsg_src), int(epsg_dst), always_xy=True)
                            xo, yo = tr.transform(e_v, n_v)
                            rec = {"X_in": e_v, "Y_in": n_v,
                                   "X_out": round(xo, 6), "Y_out": round(yo, 6)}
                        else:
                            continue
                        if col_id and col_id in row.index:
                            rec = {"ID": row[col_id], **rec}
                        rows_out.append(rec)
                    except Exception:
                        rows_out.append({"Error": f"row {i}: conversion failed"})

                pct = 100 if total else 0
                label = f"{len(rows_out)}/{total} points"
                return (_result_table(rows_out), rows_out,
                        pct, label, show_prog, False)

            # ── SURVEY MODE ──────────────────────────────────────────────────
            elif mode == "survey":
                session = survey_store or {}
                stations = session.get("stations") or []
                if not stations:
                    return (_warn("No survey data loaded."),
                            no_update, 0, "", hide_prog, True)
                z_srv = (zone_survey or "").strip() or None
                rows_out = []
                for stn in stations:
                    name = stn.get("name", "?")
                    e_v  = stn.get("easting")  or stn.get("x")
                    n_v  = stn.get("northing") or stn.get("y")
                    if e_v is None or n_v is None:
                        rows_out.append({"Station": name, "Error": "no coords"})
                        continue
                    try:
                        zo = z_srv or stn.get("zone", "49N")
                        la, lo = project_point_utm2ll(
                            float(e_v), float(n_v), zo,
                            datum=datum or "WGS84")
                        rows_out.append({"Station": name,
                                         "Easting(m)": e_v, "Northing(m)": n_v,
                                         "Zone": zo,
                                         "Latitude": round(la, 8),
                                         "Longitude": round(lo, 8)})
                    except Exception as exc:
                        rows_out.append({"Station": name, "Error": str(exc)})
                pct = 100 if rows_out else 0
                return (_result_table(rows_out), rows_out,
                        pct, f"{len(rows_out)} stations", show_prog, False)

        except Exception as exc:
            return _err(exc), no_update, 0, "", hide_prog, True

        return no_update, no_update, 0, "", hide_prog, True

    # ── Coordinate transformer — Map preview ──────────────────────────────────
    @app.callback(
        Output("tool-coord-out", "children", allow_duplicate=True),
        Input("tool-coord-map-btn",    "n_clicks"),
        State("tool-coord-result-store","data"),
        prevent_initial_call=True,
    )
    def _coord_map(n, results):
        if not n or not results:
            return no_update
        import plotly.graph_objects as go
        rows = [r for r in results if "Latitude" in r and "Longitude" in r]
        if not rows:
            return html.Div("No Lat/Lon in results — run a UTM→LL conversion first.",
                            className="text-muted small")
        lats  = [r["Latitude"]  for r in rows]
        lons  = [r["Longitude"] for r in rows]
        names = [str(r.get("Station", r.get("ID", i)))
                 for i, r in enumerate(rows)]
        fig = go.Figure(go.Scattermapbox(
            lat=lats, lon=lons, mode="markers+text",
            marker=dict(size=8, color="#89b4fa"),
            text=names, textposition="top right",
            hovertemplate="<b>%{text}</b><br>Lat: %{lat:.6f}<br>Lon: %{lon:.6f}<extra></extra>",
        ))
        fig.update_layout(
            mapbox=dict(style="carto-darkmatter",
                        center=dict(lat=sum(lats)/len(lats),
                                    lon=sum(lons)/len(lons)),
                        zoom=8),
            margin=dict(l=0, r=0, t=0, b=0),
            height=300,
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        )
        return dcc.Graph(figure=fig, config={"displaylogo": False,
                                              "scrollZoom": True},
                         style={"height": "300px"})

    # ── Coordinate transformer — Export CSV ───────────────────────────────────
    @app.callback(
        Output("tool-coord-download",    "data"),
        Input("tool-coord-export-btn",   "n_clicks"),
        State("tool-coord-result-store", "data"),
        prevent_initial_call=True,
    )
    def _coord_export(n, results):
        if not n or not results:
            return no_update
        import pandas as pd
        df = pd.DataFrame(results)
        return dict(
            content=df.to_csv(index=False),
            filename="coordinates_converted.csv",
            type="text/csv",
        )

    # Layered model — show/hide constructor param sections
    @app.callback(
        Output("tool-lm-manual-params",  "style"),
        Output("tool-lm-random-params",  "style"),
        Output("tool-lm-blocky-params",  "style"),
        Output("tool-lm-smooth-params",  "style"),
        Output("tool-lm-geology-params", "style"),
        Input("tool-lm-constructor", "value"),
        prevent_initial_call=False,
    )
    def _lm_constructor_sync(ctor):
        S, H = {"display": "block"}, {"display": "none"}
        return (
            S if ctor == "manual"  else H,
            S if ctor == "random"  else H,
            S if ctor == "blocky"  else H,
            S if ctor == "smooth"  else H,
            S if ctor == "geology" else H,
        )

    # Layered model — show/hide solver param sections
    @app.callback(
        Output("tool-lm-freq-params",   "style"),
        Output("tool-lm-tem-params",    "style"),
        Output("tool-lm-csamt-offset",  "style"),
        Input("tool-lm-solver", "value"),
        prevent_initial_call=False,
    )
    def _lm_solver_sync(solver):
        S, H = {"display": "block"}, {"display": "none"}
        return (
            S if solver in ("mt1d", "csamt1d") else H,
            S if solver == "tem1d"             else H,
            S if solver == "csamt1d"           else H,
        )

    # Layered model — show/hide noise param sections
    @app.callback(
        Output("tool-lm-gauss-params", "style"),
        Output("tool-lm-field-params", "style"),
        Output("tool-lm-mult-params",  "style"),
        Input("tool-lm-noise", "value"),
        prevent_initial_call=False,
    )
    def _lm_noise_sync(noise):
        S, H = {"display": "block"}, {"display": "none"}
        return (
            S if noise == "gaussian" else H,
            S if noise == "field"    else H,
            S if noise == "mult"     else H,
        )

    # Layered model — load quick preset into manual table
    @app.callback(
        Output("tool-lm-table", "data"),
        Input("tool-lm-preset", "value"),
        prevent_initial_call=True,
    )
    def _lm_preset(preset):
        if not preset:
            return no_update
        _PRESETS = {
            "simple": [
                {"#": 1, "ρ (Ω·m)": 100.0,  "Thickness (m)": 50.0},
                {"#": 2, "ρ (Ω·m)": 500.0,  "Thickness (m)": 200.0},
                {"#": 3, "ρ (Ω·m)": 50.0,   "Thickness (m)": 500.0},
                {"#": 4, "ρ (Ω·m)": 1000.0, "Thickness (m)": None},
            ],
            "marine": [
                {"#": 1, "ρ (Ω·m)": 0.3,   "Thickness (m)": 1500.0},
                {"#": 2, "ρ (Ω·m)": 1.0,   "Thickness (m)": 500.0},
                {"#": 3, "ρ (Ω·m)": 10.0,  "Thickness (m)": 3000.0},
                {"#": 4, "ρ (Ω·m)": 100.0, "Thickness (m)": None},
            ],
            "continental": [
                {"#": 1, "ρ (Ω·m)": 200.0,  "Thickness (m)": 2000.0},
                {"#": 2, "ρ (Ω·m)": 50.0,   "Thickness (m)": 10000.0},
                {"#": 3, "ρ (Ω·m)": 500.0,  "Thickness (m)": 20000.0},
                {"#": 4, "ρ (Ω·m)": 1000.0, "Thickness (m)": None},
            ],
        }
        return _PRESETS.get(preset, no_update)

    # Layered model — add row to manual table
    @app.callback(
        Output("tool-lm-table", "data", allow_duplicate=True),
        Input("tool-lm-add", "n_clicks"),
        State("tool-lm-table", "data"),
        prevent_initial_call=True,
    )
    def _lm_add_row(n, rows):
        if not n:
            return no_update
        rows = list(rows or [])
        rows.append({"#": len(rows) + 1, "ρ (Ω·m)": 100.0, "Thickness (m)": 100.0})
        return rows

    # Layered model — main forward run
    @app.callback(
        Output("tool-lm-out",         "children"),
        Output("tool-outputs-store",  "data", allow_duplicate=True),
        Output("tool-lm-export-btn",  "disabled"),
        Input("tool-lm-run",          "n_clicks"),
        # constructor
        State("tool-lm-constructor",  "value"),
        State("tool-lm-table",        "data"),
        State("tool-lm-r-nlayers",    "value"),
        State("tool-lm-r-rho-min",    "value"),
        State("tool-lm-r-rho-max",    "value"),
        State("tool-lm-r-depth-max",  "value"),
        State("tool-lm-r-seed",       "value"),
        State("tool-lm-b-nlayers",    "value"),
        State("tool-lm-b-rho-bg",     "value"),
        State("tool-lm-b-rho-anom",   "value"),
        State("tool-lm-b-anom-layer", "value"),
        State("tool-lm-b-depth-max",  "value"),
        State("tool-lm-b-equal-th",   "value"),
        State("tool-lm-b-seed",       "value"),
        State("tool-lm-s-nlayers",    "value"),
        State("tool-lm-s-rho-surf",   "value"),
        State("tool-lm-s-rho-deep",   "value"),
        State("tool-lm-s-depth-max",  "value"),
        State("tool-lm-s-perturb",    "value"),
        State("tool-lm-s-seed",       "value"),
        State("tool-lm-geology",      "value"),
        State("tool-lm-geo-seed",     "value"),
        # solver
        State("tool-lm-solver",       "value"),
        State("tool-lm-freq-min",     "value"),
        State("tool-lm-freq-max",     "value"),
        State("tool-lm-n-freqs",      "value"),
        State("tool-lm-src-offset",   "value"),
        State("tool-lm-time-min",     "value"),
        State("tool-lm-time-max",     "value"),
        State("tool-lm-n-times",      "value"),
        State("tool-lm-loop-radius",  "value"),
        State("tool-lm-moment",       "value"),
        # noise
        State("tool-lm-noise",        "value"),
        State("tool-lm-noise-level",  "value"),
        State("tool-lm-noise-phase",  "value"),
        State("tool-lm-noise-seed",   "value"),
        State("tool-lm-field-base",   "value"),
        State("tool-lm-field-plfreq", "value"),
        State("tool-lm-field-pllevel","value"),
        State("tool-lm-field-dead-band","value"),
        State("tool-lm-mult-sigma",   "value"),
        # display
        State("tool-lm-view",         "value"),
        State("tool-lm-depth-plot",   "value"),
        State("tool-lm-log-rho",      "value"),
        State(IDs.STORE_THEME,        "data"),
        State("tool-outputs-store",   "data"),
        prevent_initial_call=True,
    )
    def _lm_run(
        n,
        ctor, tbl_rows,
        r_nl, r_rmin, r_rmax, r_dmax, r_seed,
        b_nl, b_rbg, b_ranom, b_alyr, b_dmax, b_equal_th, b_seed,
        s_nl, s_rsurf, s_rdeep, s_dmax, s_perturb, s_seed,
        geo_name, geo_seed,
        solver, freq_min, freq_max, n_freqs, src_offset,
        time_min, time_max, n_times, loop_radius, moment,
        noise_type, noise_level, noise_phase, noise_seed,
        field_base, field_plfreq, field_pllevel, field_dead,
        mult_sigma,
        view, depth_plot, log_rho,
        theme, saved_outputs,
    ):
        if not n:
            return no_update, no_update, no_update
        try:
            import base64
            import io

            import matplotlib
            import numpy as np
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            from pycsamt.forward import (
                CSAMT1DForward,
                FieldRealisticNoise,
                GaussianNoise,
                LayeredModel,
                MT1DForward,
                MultiplicativeNoise,
                TEM1DForward,
                add_noise,
                plot_model_1d,
                plot_response_1d,
                plot_response_and_model_1d,
            )

            dark   = (theme or "dark") == "dark"
            BG     = "#1e1e2e" if dark else "#eff1f5"
            FG     = "#cdd6f4" if dark else "#4c4f69"
            GRID   = "#313244" if dark else "#ccd0da"
            AX_BG  = "#181825" if dark else "#e6e9ef"

            _lm_imgs: list = []

            def _mpl_to_img(mpl_fig: plt.Figure) -> html.Img:
                buf = io.BytesIO()
                mpl_fig.savefig(buf, format="png", dpi=120,
                                bbox_inches="tight",
                                facecolor=mpl_fig.get_facecolor())
                buf.seek(0)
                b64 = base64.b64encode(buf.read()).decode()
                plt.close(mpl_fig)
                _lm_imgs.append(b64)
                return html.Img(
                    src=f"data:image/png;base64,{b64}",
                    style={"width": "100%", "borderRadius": "6px",
                           "display": "block"},
                )

            def _style_fig(mpl_fig: plt.Figure) -> None:
                mpl_fig.patch.set_facecolor(BG)
                for ax_ in mpl_fig.axes:
                    ax_.set_facecolor(AX_BG)
                    ax_.tick_params(colors=FG, labelsize=8)
                    for sp in ax_.spines.values():
                        sp.set_edgecolor(GRID)
                    ax_.title.set_color(FG)
                    ax_.xaxis.label.set_color(FG)
                    ax_.yaxis.label.set_color(FG)
                    ax_.grid(color=GRID, alpha=0.4)

            # ── 1. BUILD MODEL ────────────────────────────────────────────
            ctor = ctor or "manual"
            model: LayeredModel

            if ctor == "manual":
                rows = tbl_rows or _LM_DEFAULT_TABLE
                rhos, thicks = [], []
                for row in rows:
                    r = float(row.get("ρ (Ω·m)", 100) or 100)
                    t = row.get("Thickness (m)")
                    rhos.append(max(r, 1e-6))
                    if t not in (None, "", "None"):
                        thicks.append(float(t))
                if len(rhos) < 2:
                    return _warn("Add at least 2 layers."), no_update, True
                if len(thicks) != len(rhos) - 1:
                    thicks = thicks[:len(rhos) - 1]
                    while len(thicks) < len(rhos) - 1:
                        thicks.append(100.0)
                model = LayeredModel(
                    resistivity=np.array(rhos),
                    thickness=np.array(thicks),
                )

            elif ctor == "random":
                model = LayeredModel.random(
                    n_layers=int(r_nl or 5),
                    rho_range=(float(r_rmin or 1.0), float(r_rmax or 10000.0)),
                    depth_max=float(r_dmax or 2000.0),
                    seed=int(r_seed) if r_seed not in (None, "") else None,
                    name="random",
                )

            elif ctor == "blocky":
                model = LayeredModel.blocky(
                    n_layers=int(b_nl or 4),
                    rho_background=float(b_rbg or 100.0),
                    rho_anomaly=float(b_ranom or 5.0),
                    anomaly_layer=int(b_alyr or 1),
                    depth_max=float(b_dmax or 1000.0),
                    equal_thickness=bool(b_equal_th),
                    seed=int(b_seed) if b_seed not in (None, "") else None,
                    name="blocky",
                )

            elif ctor == "smooth":
                model = LayeredModel.smooth(
                    n_layers=int(s_nl or 10),
                    rho_surface=float(s_rsurf or 100.0),
                    rho_deep=float(s_rdeep or 10.0),
                    depth_max=float(s_dmax or 5000.0),
                    perturbation=float(s_perturb or 0.2),
                    seed=int(s_seed) if s_seed not in (None, "") else None,
                    name="smooth",
                )

            elif ctor == "geology":
                model = LayeredModel.from_geology(
                    geo_name or "sedimentary",
                    seed=int(geo_seed) if geo_seed not in (None, "") else None,
                )

            else:
                return _warn(f"Unknown constructor: {ctor}"), no_update, True

            # ── 2. RUN FORWARD SOLVER ─────────────────────────────────────
            solver = solver or "mt1d"
            response = None

            if solver in ("mt1d", "csamt1d"):
                freqs = np.logspace(
                    np.log10(float(freq_min or 0.001)),
                    np.log10(float(freq_max or 10000.0)),
                    int(n_freqs or 30),
                )
                if solver == "csamt1d":
                    fwd = CSAMT1DForward(
                        freqs,
                        source_offset=float(src_offset or 1000.0),
                    )
                else:
                    fwd = MT1DForward(freqs)
                response = fwd.run(model)

            elif solver == "tem1d":
                times = np.logspace(
                    np.log10(float(time_min or 1e-6)),
                    np.log10(float(time_max or 0.01)),
                    int(n_times or 25),
                )
                fwd = TEM1DForward(
                    times,
                    loop_radius=float(loop_radius or 50.0),
                    moment=float(moment or 1.0),
                )
                response = fwd.run(model)

            # ── 3. ADD NOISE ──────────────────────────────────────────────
            noise_type = noise_type or "none"
            noisy_resp = None
            if noise_type != "none" and response is not None:
                try:
                    if noise_type == "gaussian":
                        lvl = float(noise_level or 5.0) / 100.0
                        ph_lvl = (float(noise_phase) / 100.0
                                  if noise_phase not in (None, "") else None)
                        nm = GaussianNoise(
                            level=lvl,
                            apply_to="rho_phase",
                            phase_level=ph_lvl,
                        )
                        seed_v = int(noise_seed) if noise_seed not in (None,"") else None
                        noisy_resp = add_noise(
                            response, noise_model=nm, level=lvl, seed=seed_v,
                        )
                    elif noise_type == "field":
                        nm = FieldRealisticNoise(
                            base_level=float(field_base or 2.0) / 100.0,
                            powerline_freq=float(field_plfreq or 50.0),
                            powerline_level=float(field_pllevel or 30.0) / 100.0,
                            dead_band=bool(field_dead),
                        )
                        noisy_resp = add_noise(response, noise_model=nm,
                                               level=float(field_base or 2.0) / 100.0)
                    elif noise_type == "mult":
                        nm = MultiplicativeNoise(
                            sigma_log10=float(mult_sigma or 0.05),
                        )
                        noisy_resp = add_noise(response, noise_model=nm, level=0.05)
                except Exception:
                    noisy_resp = None

            # ── 4. PLOT ───────────────────────────────────────────────────
            view = view or "both"
            depth_max_v = float(depth_plot) if depth_plot not in (None, "") else None
            log_rho_v   = bool(log_rho)
            children: list = []

            # ── model layer summary banner ────────────────────────────────
            layer_desc = "  ·  ".join(
                f"L{i+1}: {r:.1f} Ω·m"
                for i, r in enumerate(model.resistivity)
            )
            children.append(dbc.Alert([
                html.B(f"{model.n_layers}-layer model"),
                html.Span(f"  [{model.name}]  " if model.name else "  ",
                          className="text-muted"),
                html.Br(),
                html.Small(layer_desc, className="text-muted"),
            ], color="info", className="py-2 px-3 small mb-2"))

            if view in ("model", "both"):
                try:
                    plot_models = [model]
                    labels      = ["model"]
                    if noisy_resp is not None:
                        pass  # model is the same, noise affects response only
                    mpl_fig, ax = plt.subplots(figsize=(4.5, 5.5))
                    mpl_fig.patch.set_facecolor(BG)
                    ax.set_facecolor(AX_BG)
                    plot_model_1d(
                        plot_models, labels,
                        log_rho=log_rho_v,
                        depth_max=depth_max_v,
                        ax=ax,
                    )
                    _style_fig(mpl_fig)
                    children.append(html.Div([
                        html.P("Resistivity–depth profile",
                               className="text-muted small mb-1"),
                        _mpl_to_img(mpl_fig),
                    ], className="mb-3"))
                except Exception as exc:
                    children.append(_err(exc))

            if view in ("response", "both") and response is not None:
                try:
                    if view == "both":
                        mpl_fig = plot_response_and_model_1d(
                            noisy_resp if noisy_resp is not None else response,
                            model=model,
                            title=f"{solver.upper()} synthetic response"
                                  + (" + noise" if noisy_resp else ""),
                            figsize=(11, 5),
                        )
                    else:
                        mpl_fig_r, _ = plt.subplots(figsize=(7, 5.5))
                        mpl_fig = plt.gcf()
                        plot_response_1d(
                            noisy_resp if noisy_resp is not None else response,
                            axes=None,
                        )
                        mpl_fig = plt.gcf()
                    _style_fig(mpl_fig)
                    resp_label = (
                        f"{solver.upper()} response"
                        + (" · clean + noisy overlay" if noisy_resp else "")
                    )
                    # overlay clean dashed if noisy was added
                    if noisy_resp is not None and view == "both":
                        pass  # plot_response_and_model_1d already takes the noisy one
                    children.append(html.Div([
                        html.P(resp_label, className="text-muted small mb-1"),
                        _mpl_to_img(mpl_fig),
                    ], className="mb-3"))
                except Exception as exc:
                    children.append(_err(exc))

            if response is None:
                children.append(_warn("Forward solver produced no response."))

            # serialise for export
            model_dict = {
                "resistivity": model.resistivity.tolist(),
                "thickness":   model.thickness.tolist(),
                "depth":       model.depth.tolist(),
                "name":        model.name,
            }
            resp_dict = {}
            if response is not None:
                if response.freqs is not None:
                    resp_dict["freqs"] = response.freqs.tolist()
                if response.times is not None:
                    resp_dict["times"] = response.times.tolist()
                if response.rho_a is not None:
                    resp_dict["rho_a"] = response.rho_a.tolist()
                if response.phase is not None:
                    resp_dict["phase"] = response.phase.tolist()
                if response.dBz_dt is not None:
                    resp_dict["dBz_dt"] = response.dBz_dt.tolist()
                resp_dict["method"] = response.method

            new_saved = dict(saved_outputs or {})
            new_saved["layered-model"] = {
                "type": "lm",
                "model": model_dict,
                "response": resp_dict,
                "solver": solver,
                "imgs": _lm_imgs,
            }
            return html.Div(children), new_saved, False

        except Exception as exc:
            return _err(exc), no_update, True

    # Layered model — export
    @app.callback(
        Output("tool-lm-download",    "data"),
        Output("tool-lm-export-out",  "children"),
        Input("tool-lm-export-btn",   "n_clicks"),
        State("tool-lm-export-fmt",   "value"),
        State("tool-outputs-store",   "data"),
        prevent_initial_call=True,
    )
    def _lm_export(n, fmt, saved_outputs):
        if not n:
            return no_update, no_update
        try:
            import json

            import numpy as np
            import pandas as pd
            payload = (saved_outputs or {}).get("layered-model", {})
            if payload.get("type") != "lm":
                return no_update, html.Span(
                    "⚠ Run the forward model first.", className="text-warning")
            model_d = payload.get("model", {})
            resp_d  = payload.get("response", {})

            if fmt == "model_csv":
                rho   = model_d.get("resistivity", [])
                thick = list(model_d.get("thickness", [])) + [None]
                depth = model_d.get("depth", [])
                df = pd.DataFrame({
                    "layer":       list(range(1, len(rho) + 1)),
                    "rho_ohm_m":   rho,
                    "thickness_m": thick,
                    "depth_m":     depth,
                })
                return dcc.send_string(df.to_csv(index=False),
                                       filename="layered_model.csv"), \
                       html.Span("✓ Model CSV downloaded.", className="text-success")

            elif fmt == "resp_csv":
                if not resp_d:
                    return no_update, html.Span(
                        "⚠ No response data.", className="text-warning")
                rows: dict = {}
                if "freqs" in resp_d:
                    rows["freq_hz"] = resp_d["freqs"]
                if "times" in resp_d:
                    rows["time_s"] = resp_d["times"]
                for key in ("rho_a", "phase", "dBz_dt"):
                    if key in resp_d:
                        arr = np.asarray(resp_d[key])
                        if arr.ndim == 2:
                            for j in range(arr.shape[1]):
                                rows[f"{key}_{j}"] = arr[:, j].tolist()
                        else:
                            rows[key] = arr.tolist()
                df = pd.DataFrame(rows)
                return dcc.send_string(df.to_csv(index=False),
                                       filename="forward_response.csv"), \
                       html.Span("✓ Response CSV downloaded.", className="text-success")

            else:
                blob = json.dumps({"model": model_d, "response": resp_d}, indent=2)
                return dcc.send_string(blob, filename="layered_model.json"), \
                       html.Span("✓ JSON downloaded.", className="text-success")

        except Exception as exc:
            return no_update, _err(exc)

    # Strike analysis
    @app.callback(
        Output("tool-strike-out",    "children"),
        Output("tools-run-store",    "data"),
        Output("tool-outputs-store", "data", allow_duplicate=True),
        Input("tool-strike-run", "n_clicks"),
        State("tool-strike-lines", "value"),
        State("tool-strike-stations", "value"),
        State("tool-strike-method", "value"),
        State("tool-strike-freq",   "value"),
        State("tool-strike-step",   "value"),
        State("tool-strike-maxrows", "value"),
        State(IDs.STORE_DATA,       "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.SESSION_ID,       "data"),
        State(IDs.STORE_THEME,      "data"),
        State("tools-run-store",    "data"),
        State("tool-outputs-store", "data"),
        prevent_initial_call=True,
    )
    def _run_strike(
        n,
        selected_lines,
        selected_stations,
        method,
        freq_band,
        angle_step,
        max_rows,
        store_data,
        active_lines_store,
        session_id,
        theme,
        run_store,
        saved_outputs,
    ):
        if not n:
            return no_update, no_update, no_update
        try:
            import numpy as np

            from pycsamt.emtools.strike import (
                estimate_strike_consensus,
                estimate_strike_phase_tensor,
                estimate_strike_sweep,
            )

            sites, msg = _scoped_sites(
                session_id,
                store_data,
                active_lines_store,
                selected_lines=selected_lines,
                selected_stations=selected_stations,
            )
            if msg:
                return _warn(msg), no_update

            band = _strike_band_from_choice(freq_band)
            step = float(angle_step or 1.0)
            step = max(0.25, min(10.0, step))
            angles = np.arange(-90.0, 90.0 + step, step)
            method = (method or "consensus").lower()
            if method == "sweep":
                df = estimate_strike_sweep(
                    sites, band=band, angles=angles,
                    recursive=False, strict=False, verbose=0,
                )
                method_label = "Z tensor rotation sweep"
            elif method == "pt":
                df = estimate_strike_phase_tensor(
                    sites, band=band, recursive=False, strict=False, verbose=0
                )
                method_label = "Phase tensor strike"
            else:
                df = estimate_strike_consensus(
                    sites, band=band, angles=angles,
                    recursive=False, strict=False, verbose=0,
                )
                method_label = "Consensus strike"

            if df is None or df.empty:
                return _warn(
                    "No strike estimates were produced for the selected scope. "
                    "Check that the selected stations have valid impedance tensors."
                ), no_update

            line_map = _station_line_map(store_data)
            df = df.copy()
            df["line"] = df["station"].map(line_map).fillna("Selected")
            df["ang_axial"] = ((df["ang"].astype(float) + 180.0) % 180.0)
            finite = df[np.isfinite(df["ang_axial"])]
            if finite.empty:
                return _warn("Strike estimates are all non-finite for this scope."), no_update

            table_df = finite.copy()
            for col in ("ang", "iqr", "lo", "hi"):
                if col in table_df:
                    table_df[col] = table_df[col].astype(float).round(3)
            table_df = table_df.rename(
                columns={
                    "station": "Station",
                    "line": "Line",
                    "ang": "Strike (°)",
                    "iqr": "IQR (°)",
                    "lo": "T min (s)",
                    "hi": "T max (s)",
                    "n": "N",
                }
            )
            keep_cols = [
                c for c in
                ["Station", "Line", "Strike (°)", "IQR (°)", "T min (s)", "T max (s)", "N"]
                if c in table_df.columns
            ]
            max_rows = int(max_rows or 40)
            max_rows = max(5, min(200, max_rows))
            med = float(np.nanmedian(finite["ang_axial"]))
            spread = float(np.nanmedian(finite["iqr"])) if "iqr" in finite else float("nan")
            scope_label = _strike_scope_label(
                store_data, active_lines_store, selected_lines, selected_stations
            )
            payload = {
                "tool": "strike",
                "method_label": method_label,
                "scope_label": scope_label,
                "n_stations": int(finite["station"].nunique()),
                "median_strike": med,
                "median_iqr": spread,
                "records": finite[["station", "line", "ang_axial"]].to_dict("records"),
                "table": table_df[keep_cols].head(max_rows).to_dict("records"),
                "columns": keep_cols,
                "page_size": min(10, max_rows),
            }
            run_store = dict(run_store or {})
            run_store["strike"] = payload
            new_saved = dict(saved_outputs or {})
            new_saved["strike"] = {"type": "strike_payload", "payload": payload, "theme": theme}
            return _render_strike_result(payload, theme), run_store, new_saved
        except Exception as exc:
            return _err(exc), no_update, no_update

    # EDI Validator
    @app.callback(
        Output("tool-valid-out",     "children"),
        Output("tool-outputs-store", "data", allow_duplicate=True),
        Input("tool-valid-run", "n_clicks"),
        State("tool-valid-lines", "value"),
        State("tool-valid-stations", "value"),
        State("tool-valid-severity", "value"),
        State("tool-valid-minfreq", "value"),
        State("tool-valid-checks", "value"),
        State(IDs.STORE_DATA,  "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.SESSION_ID,  "data"),
        State("tool-outputs-store", "data"),
        prevent_initial_call=True,
    )
    def _run_validator(
        n,
        selected_lines,
        selected_stations,
        severity,
        min_freq,
        checks,
        store_data,
        active_lines_store,
        session_id,
        saved_outputs,
    ):
        if not n:
            return no_update, no_update
        try:
            import pandas as pd

            sites, msg = _scoped_sites(
                session_id,
                store_data,
                active_lines_store,
                selected_lines=selected_lines,
                selected_stations=selected_stations,
            )
            if msg:
                return _warn(msg), no_update

            checks = checks or []
            min_freq = int(min_freq or 3)
            rows = _validate_sites_rows(
                sites,
                store_data=store_data,
                checks=set(checks),
                min_freq=min_freq,
            )
            if not rows:
                return _warn("No stations were available for validation."), no_update

            df = pd.DataFrame(rows)
            severity = (severity or "all").lower()
            if severity == "issues":
                shown = df[df["Status"].isin(["WARN", "FAIL"])]
            elif severity == "fail":
                shown = df[df["Status"].eq("FAIL")]
            elif severity == "pass":
                shown = df[df["Status"].eq("PASS")]
            else:
                shown = df

            counts = df["Status"].value_counts().to_dict()
            n_pass = int(counts.get("PASS", 0))
            n_warn = int(counts.get("WARN", 0))
            n_fail = int(counts.get("FAIL", 0))
            scope_label = _strike_scope_label(
                store_data, active_lines_store, selected_lines, selected_stations
            )
            columns = [
                "Station", "Line", "Lat", "Lon", "Z ok", "Errors", "Tipper",
                "N freq", "T min (s)", "T max (s)", "NaN rows", "Status",
                "Issues",
            ]
            shown = shown[[c for c in columns if c in shown.columns]]
            shown_records = shown.to_dict("records")
            style_rows = [
                {
                    "if": {"filter_query": '{Status} = "PASS"'},
                    "backgroundColor": "rgba(166, 227, 161, 0.12)",
                },
                {
                    "if": {"filter_query": '{Status} = "WARN"'},
                    "backgroundColor": "rgba(249, 226, 175, 0.14)",
                },
                {
                    "if": {"filter_query": '{Status} = "FAIL"'},
                    "backgroundColor": "rgba(243, 139, 168, 0.16)",
                },
            ]
            result = html.Div(
                [
                    html.Div(
                        [
                            _metric_chip("Scope", scope_label),
                            _metric_chip("Validated", str(len(df))),
                            _metric_chip("Pass", str(n_pass)),
                            _metric_chip("Warn / Fail", f"{n_warn} / {n_fail}"),
                        ],
                        className="tool-metric-row",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Span("Validation policy", className="tool-section-title"),
                                    html.Span(
                                        f"Min frequencies: {min_freq} · "
                                        f"Checks: {', '.join(checks) if checks else 'none'}",
                                        className="tool-section-subtitle",
                                    ),
                                ],
                                className="tool-section-head",
                            ),
                            _validator_status_bar(n_pass, n_warn, n_fail),
                        ],
                        className="tool-valid-summary",
                    ),
                    dash_table.DataTable(
                        columns=[{"name": c, "id": c} for c in shown.columns],
                        data=shown_records,
                        page_size=12,
                        sort_action="native",
                        filter_action="native",
                        tooltip_data=[
                            {
                                "Issues": {
                                    "value": str(row.get("Issues", "")),
                                    "type": "markdown",
                                }
                            }
                            for row in shown_records
                        ],
                        tooltip_duration=None,
                        style_table={"overflowX": "auto", "maxHeight": "420px", "overflowY": "auto"},
                        style_cell={
                            "fontSize": "11px",
                            "fontFamily": "Inter, sans-serif",
                            "padding": "6px",
                            "backgroundColor": "var(--mantle)",
                            "color": "var(--text)",
                            "border": "1px solid var(--surface0)",
                            "maxWidth": "220px",
                            "overflow": "hidden",
                            "textOverflow": "ellipsis",
                        },
                        style_header={
                            "backgroundColor": "var(--surface0)",
                            "fontWeight": "700",
                        },
                        style_data_conditional=style_rows,
                    ),
                ]
            )
            new_saved = dict(saved_outputs or {})
            new_saved["validator"] = {
                "type": "validator_rows",
                "rows": shown_records,
                "cols": shown.columns.tolist(),
                "summary": {
                    "scope": scope_label, "total": len(df),
                    "pass": n_pass, "warn": n_warn, "fail": n_fail,
                },
            }
            return result, new_saved
        except Exception as exc:
            return _err(exc), no_update

    # Format converter
    @app.callback(
        Output("tool-conv-out", "children"),
        Output(IDs.STORE_DATA, "data", allow_duplicate=True),
        Input("tool-conv-run", "n_clicks"),
        State("tool-conv-type", "value"),
        State("tool-conv-source", "value"),
        State("tool-conv-output", "value"),
        State("tool-conv-stn", "value"),
        State("tool-conv-utm", "value"),
        State("tool-conv-epsg", "value"),
        State("tool-conv-convert-coords", "value"),
        State("tool-conv-elabels", "value"),
        State("tool-conv-hlabels", "value"),
        State("tool-conv-suffix", "value"),
        State("tool-conv-spectra-flags", "value"),
        State("tool-conv-freq-order", "value"),
        State("tool-conv-freq-tol", "value"),
        State("tool-conv-core-flags", "value"),
        State("tool-conv-station-name", "value"),
        State(IDs.STORE_DATA,   "data"),
        State(IDs.SESSION_ID,   "data"),
        prevent_initial_call=True,
    )
    def _run_converter(
        n,
        conv_type,
        source_path,
        output_dir,
        stn_path,
        utm_zone,
        epsg,
        convert_coords,
        e_labels,
        h_labels,
        suffix,
        spectra_flags,
        freq_order,
        freq_tol,
        core_flags,
        station_name,
        store_data,
        session_id,
    ):
        if not n:
            return no_update, no_update
        if not source_path or not str(source_path).strip():
            return _warn("Enter a source AVG/J/Spectra path first."), no_update
        try:
            from pathlib import Path

            from pycsamt.app.desktop.controllers.advanced_controller import (
                ConversionController,
            )

            source = Path(str(source_path).strip()).expanduser()
            if not source.exists():
                return _warn(f"Source path does not exist: {source}"), no_update
            out = str(output_dir or "").strip()
            if out:
                Path(out).expanduser().mkdir(parents=True, exist_ok=True)

            options = _converter_options_from_values(
                output_dir=out,
                stn_path=stn_path,
                utm_zone=utm_zone,
                epsg=epsg,
                convert_coords=convert_coords,
                e_labels=e_labels,
                h_labels=h_labels,
                suffix=suffix,
                spectra_flags=spectra_flags,
                freq_order=freq_order,
                freq_tol=freq_tol,
                core_flags=core_flags,
                station_name=station_name,
            )
            ctrl = ConversionController()
            ctrl.set_source(str(conv_type or "AVG -> EDI"), str(source))
            collection, failures = ctrl.run(options)
            stats = ctrl.build_stats(collection, failures)

            load_session = "load_session" in set(core_flags or [])
            new_store = no_update
            if load_session and collection is not None:
                from pycsamt.app.web.cache import cache_set
                from pycsamt.emtools._core import ensure_sites

                sites = ensure_sites(
                    collection,
                    recursive=False,
                    on_dup="replace",
                    strict=False,
                    verbose=0,
                )
                cache_set(session_id, sites)
                new_store = _converted_sites_to_store(
                    sites,
                    old_store=store_data,
                    data_dir=str(out or source.parent),
                )

            return _converter_result_view(
                stats=stats,
                failures=failures,
                output_dir=out,
                loaded=load_session and new_store is not no_update,
            ), new_store
        except Exception as exc:
            return _err(exc), no_update

    # Dimensionality classifier — station scope sync
    @app.callback(
        Output("tool-dim-stations", "options"),
        Output("tool-dim-stations", "value"),
        Input("tool-dim-lines", "value"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def _sync_dim_station_scope(selected_lines, store_data, active_lines_store):
        opts = _station_options_for_lines(store_data, active_lines_store, selected_lines)
        return opts, None

    # Dimensionality classifier — main run
    @app.callback(
        Output("tool-dim-out",       "children"),
        Output("tool-outputs-store", "data", allow_duplicate=True),
        Input("tool-dim-run",         "n_clicks"),
        State("tool-dim-lines",       "value"),
        State("tool-dim-stations",    "value"),
        State("tool-dim-view",        "value"),
        State("tool-dim-skew-th",     "value"),
        State("tool-dim-ellipt-th",   "value"),
        State("tool-dim-period-lo",   "value"),
        State("tool-dim-period-hi",   "value"),
        State("tool-dim-map-period",  "value"),
        State("tool-dim-n-bands",     "value"),
        State(IDs.STORE_DATA,         "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.SESSION_ID,         "data"),
        State(IDs.STORE_THEME,        "data"),
        State("tool-outputs-store",   "data"),
        prevent_initial_call=True,
    )
    def _run_dimensionality(
        n, sel_lines, sel_stations, view,
        skew_th, ellipt_th, period_lo, period_hi,
        map_period, n_bands,
        store_data, active_lines_store, session_id, theme, saved_outputs,
    ):
        if not n:
            return no_update, no_update
        try:
            import base64
            import io

            import matplotlib
            import numpy as np
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import plotly.graph_objects as go

            from pycsamt.emtools.dimensionality import (
                classify_dimensionality,
                plot_dim_confidence_grid,
                plot_dim_map,
                plot_dim_occupancy_area,
            )

            sub, err = _scoped_sites(
                session_id, store_data, active_lines_store,
                selected_lines=sel_lines or [],
                selected_stations=sel_stations or [],
            )
            if sub is None:
                return _warn(err), no_update

            # ── resolved params ───────────────────────────────────────────
            view       = view or "psection"
            skew_th    = float(skew_th  or 3.0)
            ellipt_th  = float(ellipt_th or 0.2)
            map_period = float(map_period or 10.0)
            n_bands    = int(n_bands or 24)
            p_lo       = float(period_lo) if period_lo not in (None, "") else None
            p_hi       = float(period_hi) if period_hi not in (None, "") else None

            dark   = (theme or "dark") == "dark"
            BG     = "#1e1e2e" if dark else "#eff1f5"
            FG     = "#cdd6f4" if dark else "#4c4f69"
            GRID   = "#313244" if dark else "#ccd0da"
            AX_BG  = "#181825" if dark else "#e6e9ef"
            C_1D   = "#a6e3a1"
            C_2D   = "#89b4fa"
            C_3D   = "#f38ba8"

            _dim_imgs: list = []

            # ── common matplotlib helpers ──────────────────────────────────
            def _mpl_to_img(mpl_fig: plt.Figure) -> html.Img:
                buf = io.BytesIO()
                mpl_fig.savefig(buf, format="png", dpi=120,
                                bbox_inches="tight",
                                facecolor=mpl_fig.get_facecolor())
                buf.seek(0)
                b64 = base64.b64encode(buf.read()).decode()
                plt.close(mpl_fig)
                _dim_imgs.append(b64)
                return html.Img(
                    src=f"data:image/png;base64,{b64}",
                    style={"width": "100%", "borderRadius": "6px", "display": "block"},
                )

            def _style_ax(ax_: plt.Axes, fig_: plt.Figure) -> None:
                fig_.patch.set_facecolor(BG)
                ax_.set_facecolor(AX_BG)
                ax_.tick_params(colors=FG, labelsize=8)
                for sp in ax_.spines.values():
                    sp.set_edgecolor(GRID)
                ax_.title.set_color(FG)
                ax_.xaxis.label.set_color(FG)
                ax_.yaxis.label.set_color(FG)

            def _mpl_section(plot_fn, **kw) -> html.Img:
                mpl_fig, ax = plt.subplots(figsize=(10, 4.5))
                _style_ax(ax, mpl_fig)
                plot_fn(sub, skew_th=skew_th, ellipt_th=ellipt_th, ax=ax,
                        verbose=0, **kw)
                _style_ax(ax, mpl_fig)
                # style any colourbar axes appended by the plot function
                for ax_cb in mpl_fig.axes[1:]:
                    ax_cb.tick_params(colors=FG, labelsize=7)
                    ax_cb.yaxis.label.set_color(FG)
                    try:
                        ax_cb.spines["outline"].set_edgecolor(GRID)
                    except Exception:
                        pass
                return _mpl_to_img(mpl_fig)

            # ── classify (needed for stats + scatter) ─────────────────────
            df = classify_dimensionality(
                sub, skew_th=skew_th, ellipt_th=ellipt_th, verbose=0,
            )

            # apply optional period filter to the DataFrame for scatter/stats
            df_filt = df.copy()
            if p_lo is not None:
                df_filt = df_filt[df_filt["period"] >= p_lo]
            if p_hi is not None:
                df_filt = df_filt[df_filt["period"] <= p_hi]

            _DIM_LABEL = {0: "1D", 1: "2D", 2: "3D"}
            _DIM_COLOR = {0: C_1D, 1: C_2D, 2: C_3D}

            children: list = []

            # ── 1. Confidence grid pseudo-section ─────────────────────────
            if view in ("psection", "all"):
                try:
                    img = _mpl_section(plot_dim_confidence_grid)
                    children.append(html.Div([
                        html.P("Confidence grid — 1D(green) / 2D(blue) / 3D(red), "
                               "α = confidence",
                               className="text-muted small mb-1"),
                        img,
                    ], className="mb-3"))
                except Exception as exc:
                    children.append(_err(exc))

            # ── 2. Occupancy stacked-area ─────────────────────────────────
            if view in ("occupancy", "all"):
                try:
                    mpl_fig, ax = plt.subplots(figsize=(10, 3.2))
                    _style_ax(ax, mpl_fig)
                    plot_dim_occupancy_area(
                        sub, skew_th=skew_th, ellipt_th=ellipt_th,
                        n_bands=n_bands, ax=ax, verbose=0,
                    )
                    _style_ax(ax, mpl_fig)
                    legend = ax.get_legend()
                    if legend:
                        legend.get_frame().set_facecolor(AX_BG)
                        for t in legend.get_texts():
                            t.set_color(FG)
                    children.append(html.Div([
                        html.P("Dimensionality occupancy by period band",
                               className="text-muted small mb-1"),
                        _mpl_to_img(mpl_fig),
                    ], className="mb-3"))
                except Exception as exc:
                    children.append(_err(exc))

            # ── 3. Spatial map ────────────────────────────────────────────
            if view in ("map", "all"):
                try:
                    mpl_fig, ax = plt.subplots(figsize=(8, 5.5))
                    _style_ax(ax, mpl_fig)
                    plot_dim_map(
                        sub, period=map_period,
                        skew_th=skew_th, ellipt_th=ellipt_th,
                        ax=ax, verbose=0,
                    )
                    _style_ax(ax, mpl_fig)
                    legend = ax.get_legend()
                    if legend:
                        legend.get_frame().set_facecolor(AX_BG)
                        for t in legend.get_texts():
                            t.set_color(FG)
                    children.append(html.Div([
                        html.P(
                            f"Spatial map @ T = {map_period} s  "
                            "— ○ 1D  ■ 2D  △ 3D, size ∝ confidence",
                            className="text-muted small mb-1",
                        ),
                        _mpl_to_img(mpl_fig),
                    ], className="mb-3"))
                except Exception as exc:
                    children.append(_err(exc))

            # ── 4. β / ellipticity diagnostic scatter (Plotly) ────────────
            if view in ("scatter", "all"):
                try:
                    if df_filt.empty:
                        children.append(_warn("No data in selected period range."))
                    else:
                        scatter_traces = []
                        for d_val, grp in df_filt.groupby("dim"):
                            lbl = _DIM_LABEL.get(int(d_val), str(d_val))
                            col = _DIM_COLOR.get(int(d_val), "#aaa")
                            scatter_traces.append(go.Scatter(
                                x=grp["beta_abs"].tolist(),
                                y=grp["ellipt_abs"].tolist(),
                                mode="markers",
                                name=lbl,
                                marker={"color": col, "size": 5, "opacity": 0.75},
                                text=grp["station"].tolist(),
                                customdata=np.log10(
                                    np.maximum(grp["period"].to_numpy(), 1e-9)
                                ).tolist(),
                                hovertemplate=(
                                    "<b>%{text}</b><br>"
                                    "β = %{x:.2f}°<br>"
                                    "ε = %{y:.3f}<br>"
                                    "log₁₀T = %{customdata:.2f}<extra></extra>"
                                ),
                            ))
                        # threshold reference lines
                        x_max = float(df_filt["beta_abs"].max()) * 1.1 + 1.0
                        y_max = float(df_filt["ellipt_abs"].max()) * 1.1 + 0.05
                        scatter_traces += [
                            go.Scatter(
                                x=[skew_th, skew_th], y=[0, y_max],
                                mode="lines",
                                line={"color": "#f38ba8", "width": 1.2,
                                      "dash": "dash"},
                                name=f"skew_th={skew_th}°",
                                hoverinfo="skip",
                            ),
                            go.Scatter(
                                x=[0, x_max], y=[ellipt_th, ellipt_th],
                                mode="lines",
                                line={"color": "#fab387", "width": 1.2,
                                      "dash": "dot"},
                                name=f"ellipt_th={ellipt_th}",
                                hoverinfo="skip",
                            ),
                        ]
                        scat_fig = go.Figure(scatter_traces)
                        _style_fig(scat_fig, "Ellipticity ε", height=320)
                        scat_fig.update_xaxes(title_text="Phase-sensitive skew β (°)")
                        scat_fig.update_yaxes(title_text="Ellipticity ε")
                        p_note = ""
                        if p_lo or p_hi:
                            p_note = (f"  ·  Period range: "
                                      f"[{p_lo or 'auto'} – {p_hi or 'auto'}] s")
                        children.append(html.Div([
                            html.P(
                                f"β vs ε scatter — skew_th={skew_th}°, "
                                f"ellipt_th={ellipt_th}{p_note}",
                                className="text-muted small mb-1",
                            ),
                            dcc.Graph(
                                figure=scat_fig,
                                config=_plotly_config("pycsamt_dim_scatter"),
                                className="tool-graph",
                            ),
                        ], className="mb-3"))
                except Exception as exc:
                    children.append(_err(exc))

            # ── 5. Summary stats bar (always shown) ───────────────────────
            try:
                if not df_filt.empty:
                    counts     = df_filt["dim"].value_counts().sort_index()
                    bar_labels = [_DIM_LABEL.get(int(i), str(i))
                                  for i in counts.index]
                    bar_colors = [_DIM_COLOR.get(int(i), "#aaa")
                                  for i in counts.index]
                    bar_fig = go.Figure(go.Bar(
                        x=bar_labels, y=counts.values.tolist(),
                        marker_color=bar_colors,
                        text=counts.values.tolist(),
                        textposition="auto",
                        hovertemplate="%{x}: %{y} measurements<extra></extra>",
                    ))
                    _style_fig(bar_fig, "Count", height=200)
                    bar_fig.update_layout(
                        xaxis_title="Class",
                        yaxis_title="# measurements",
                        showlegend=False,
                    )
                    total = int(counts.sum())
                    pcts  = [f"{_DIM_LABEL.get(int(i), str(i))} "
                             f"{v/total*100:.0f}% ({v})"
                             for i, v in counts.items()]
                    children.append(html.Div([
                        html.P(
                            ["Summary: ", html.Span(" · ".join(pcts),
                                                    className="text-muted")],
                            className="small mb-1",
                        ),
                        dcc.Graph(
                            figure=bar_fig,
                            config=_plotly_config("pycsamt_dim_counts"),
                            className="tool-graph",
                        ),
                    ], className="mb-2"))
            except Exception:
                pass

            if not children:
                children.append(_warn("No output produced. Check your data and thresholds."))

            new_saved = dict(saved_outputs or {})
            new_saved["dimensionality"] = {"type": "html", "imgs": _dim_imgs}
            return html.Div(children), new_saved

        except Exception as exc:
            return _err(exc), no_update

    # Frequency editor — station scope sync
    @app.callback(
        Output("tool-freq-stations", "options"),
        Output("tool-freq-stations", "value"),
        Input("tool-freq-lines", "value"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def _sync_freq_station_scope(selected_lines, store_data, active_lines_store):
        opts = _station_options_for_lines(store_data, active_lines_store, selected_lines)
        return opts, None

    # Frequency editor — show/hide folder row when source changes
    @app.callback(
        Output("tool-freq-folder-row", "style"),
        Input("tool-freq-source", "value"),
        prevent_initial_call=False,
    )
    def _sync_freq_source(source):
        return {"display": "block"} if source == "folder" else {"display": "none"}

    # Frequency editor — show/hide param sections per action
    @app.callback(
        Output("tool-freq-conf-params",   "style"),
        Output("tool-freq-band-params",   "style"),
        Output("tool-freq-regrid-params", "style"),
        Output("tool-freq-align-params",  "style"),
        Output("tool-freq-smooth-params", "style"),
        Input("tool-freq-action", "value"),
        prevent_initial_call=False,
    )
    def _sync_freq_action_params(action):
        show = {"display": "block"}
        hide = {"display": "none"}
        return (
            show if action == "confidence" else hide,
            show if action == "band"       else hide,
            show if action == "regrid"     else hide,
            show if action == "align"      else hide,
            show if action == "smooth"     else hide,
        )

    # Frequency editor — main run
    @app.callback(
        Output("tool-freq-out",          "children"),
        Output("tool-outputs-store",     "data", allow_duplicate=True),
        Output("tool-freq-export-btn",   "disabled"),
        Input("tool-freq-run",           "n_clicks"),
        State("tool-freq-source",        "value"),
        State("tool-freq-folder",        "value"),
        State("tool-freq-lines",         "value"),
        State("tool-freq-stations",      "value"),
        State("tool-freq-action",        "value"),
        # confidence params
        State("tool-freq-mode",          "value"),
        State("tool-freq-method",        "value"),
        State("tool-freq-thresh",        "value"),
        State("tool-freq-ci-hi",         "value"),
        State("tool-freq-ci-lo",         "value"),
        State("tool-freq-also",          "value"),
        State("tool-freq-interp",        "value"),
        State("tool-freq-reject",        "value"),
        # band params
        State("tool-freq-fmin",          "value"),
        State("tool-freq-fmax",          "value"),
        State("tool-freq-keep",          "value"),
        # regrid params
        State("tool-freq-rg-fmin",       "value"),
        State("tool-freq-rg-fmax",       "value"),
        State("tool-freq-rg-ppd",        "value"),
        State("tool-freq-rg-method",     "value"),
        # align params
        State("tool-freq-align-mode",    "value"),
        State("tool-freq-align-ref",     "value"),
        # smooth params
        State("tool-freq-smooth-op",     "value"),
        State("tool-freq-smooth-k",      "value"),
        State("tool-freq-smooth-step",   "value"),
        # session
        State(IDs.STORE_DATA,            "data"),
        State(IDs.STORE_ACTIVE_LINES,    "data"),
        State(IDs.SESSION_ID,            "data"),
        State(IDs.STORE_THEME,           "data"),
        State("tool-outputs-store",      "data"),
        prevent_initial_call=True,
    )
    def _run_freq_editor(
        n,
        source, folder,
        sel_lines, sel_stations,
        action,
        # confidence
        mode, method, thresh, ci_hi, ci_lo, also, interp, reject,
        # band
        fmin, fmax, keep_str,
        # regrid
        rg_fmin, rg_fmax, rg_ppd, rg_method,
        # align
        align_mode, align_ref,
        # smooth
        smooth_op, smooth_k, smooth_step,
        # session
        store_data, active_lines_store, session_id, theme, saved_outputs,
    ):
        if not n:
            return no_update, no_update, no_update
        try:
            import base64
            import io

            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            from pycsamt.app.web.cache import (
                cache_set_freq_edit,
            )
            from pycsamt.emtools._core import ensure_sites
            from pycsamt.emtools.frequency import (
                align_grid,
                decimate_step,
                drop_duplicates,
                edit_frequencies_by_confidence,
                plot_apparent_depth_psection,
                plot_band_microstrips,
                plot_coverage_quality_heatmap,
                plot_frequency_edit_decisions,
                plot_frequency_edit_summary,
                regrid_logspace,
                select_band,
                smooth_mavg,
            )

            # ── resolve source ────────────────────────────────────────────
            if source == "folder" and folder:
                import pathlib
                p = pathlib.Path(folder)
                if not p.is_dir():
                    return _warn(f"Folder not found: {folder}"), no_update, True
                sub = ensure_sites(str(p), recursive=True, on_dup="replace",
                                   strict=False, verbose=0)
                if sub is None:
                    return _warn("No EDI files found in folder."), no_update, True
            else:
                sub, err = _scoped_sites(
                    session_id, store_data, active_lines_store,
                    selected_lines=sel_lines or [],
                    selected_stations=sel_stations or [],
                )
                if sub is None:
                    return _warn(err), no_update, True

            action    = action or "diagnostics"
            dark      = (theme or "dark") == "dark"
            BG        = "#1e1e2e" if dark else "#eff1f5"
            FG        = "#cdd6f4" if dark else "#4c4f69"
            GRID      = "#313244" if dark else "#ccd0da"
            AX_BG     = "#181825" if dark else "#e6e9ef"

            _freq_imgs: list = []

            # ── shared mpl helpers ────────────────────────────────────────
            def _mpl_to_img(fig_: plt.Figure) -> html.Img:
                buf = io.BytesIO()
                fig_.savefig(buf, format="png", dpi=120,
                             bbox_inches="tight",
                             facecolor=fig_.get_facecolor())
                buf.seek(0)
                b64 = base64.b64encode(buf.read()).decode()
                plt.close(fig_)
                _freq_imgs.append(b64)
                return html.Img(
                    src=f"data:image/png;base64,{b64}",
                    style={"width": "100%", "borderRadius": "6px",
                           "display": "block"},
                )

            def _style_ax(ax_: plt.Axes, fig_: plt.Figure) -> None:
                fig_.patch.set_facecolor(BG)
                ax_.set_facecolor(AX_BG)
                ax_.tick_params(colors=FG, labelsize=8)
                for sp in ax_.spines.values():
                    sp.set_edgecolor(GRID)
                ax_.title.set_color(FG)
                ax_.xaxis.label.set_color(FG)
                ax_.yaxis.label.set_color(FG)

            def _mpl_wrap(plot_fn, *args, figsize=(10, 4.2), **kw) -> html.Img:
                mpl_fig, ax = plt.subplots(figsize=figsize)
                _style_ax(ax, mpl_fig)
                try:
                    plot_fn(*args, ax=ax, **kw)
                except TypeError:
                    plot_fn(*args, **kw)
                _style_ax(ax, mpl_fig)
                leg = ax.get_legend()
                if leg:
                    leg.get_frame().set_facecolor(AX_BG)
                    for t in leg.get_texts():
                        t.set_color(FG)
                return _mpl_to_img(mpl_fig)

            children: list = []
            has_edit = False  # becomes True when edited Sites are cached

            # ── DIAGNOSTICS ───────────────────────────────────────────────
            if action == "diagnostics":
                try:
                    children.append(html.P(
                        "Coverage & quality diagnostics — no data modified.",
                        className="text-muted small mb-2",
                    ))
                    # 1. Coverage/quality heatmap
                    children.append(html.Div([
                        html.P("Frequency coverage quality heatmap",
                               className="text-muted small mb-1"),
                        _mpl_wrap(plot_coverage_quality_heatmap, sub,
                                  axis="period", verbose=0),
                    ], className="mb-3"))
                    # 2. Band micro-strips
                    children.append(html.Div([
                        html.P("Band micro-strips — data presence per station",
                               className="text-muted small mb-1"),
                        _mpl_wrap(plot_band_microstrips, sub,
                                  figsize=(10, 5), verbose=0),
                    ], className="mb-3"))
                    # 3. Apparent depth pseudo-section
                    children.append(html.Div([
                        html.P("Apparent investigation depth pseudo-section",
                               className="text-muted small mb-1"),
                        _mpl_wrap(plot_apparent_depth_psection, sub,
                                  axis_y="period", verbose=0),
                    ], className="mb-3"))
                except Exception as exc:
                    children.append(_err(exc))

            # ── CONFIDENCE EDIT ───────────────────────────────────────────
            elif action == "confidence":
                try:
                    thresh  = float(thresh  or 0.5)
                    ci_hi_v = float(ci_hi   or 0.9)
                    ci_lo_v = float(ci_lo   or 0.5)
                    result  = edit_frequencies_by_confidence(
                        sub,
                        mode          = mode     or "recover",
                        method        = method   or "composite",
                        threshold     = thresh,
                        ci_hi         = ci_hi_v,
                        ci_lo         = ci_lo_v,
                        also          = also     or "both",
                        interpolation = interp   or "linear",
                        reject        = reject   or "drop",
                        inplace       = False,
                        verbose       = 0,
                    )
                    edited = result.sites
                    cache_set_freq_edit(session_id, edited)
                    has_edit = True

                    # summary banner
                    children.append(dbc.Alert([
                        html.I(className="bi bi-check-circle me-2"),
                        f"Mode: {result.mode} · "
                        f"Dropped: {result.n_dropped} · "
                        f"Masked: {result.n_masked} · "
                        f"Recovered: {result.n_recovered}",
                    ], color="success", className="py-2 px-3 small mb-2"))

                    # before/after summary plot
                    try:
                        mpl_fig, ax = plt.subplots(figsize=(10, 4.0))
                        _style_ax(ax, mpl_fig)
                        plot_frequency_edit_summary(
                            sub, edited,
                            method=method or "composite",
                            ci_hi=ci_hi_v, ci_lo=ci_lo_v,
                            ax=ax,
                        )
                        _style_ax(ax, mpl_fig)
                        children.append(html.Div([
                            html.P("Before / after station-level summary",
                                   className="text-muted small mb-1"),
                            _mpl_to_img(mpl_fig),
                        ], className="mb-3"))
                    except Exception as exc:
                        children.append(_err(exc))

                    # decision grid
                    try:
                        mpl_fig2, ax2 = plt.subplots(figsize=(10, 5.0))
                        _style_ax(ax2, mpl_fig2)
                        plot_frequency_edit_decisions(
                            sub, edited,
                            method=method or "composite",
                            ci_hi=ci_hi_v, ci_lo=ci_lo_v,
                            ax=ax2,
                        )
                        _style_ax(ax2, mpl_fig2)
                        children.append(html.Div([
                            html.P("Frequency edit decision grid "
                                   "(keep / drop / mask / recover)",
                                   className="text-muted small mb-1"),
                            _mpl_to_img(mpl_fig2),
                        ], className="mb-3"))
                    except Exception as exc:
                        children.append(_err(exc))

                    # station report table
                    if not result.report.empty:
                        rpt = result.report.copy()
                        rpt.columns = [c.replace("_", " ").title()
                                       for c in rpt.columns]
                        rpt_num = rpt.select_dtypes("number")
                        rpt[rpt_num.columns] = rpt_num.round(3)
                        children.append(html.Div([
                            html.P("Station-level report",
                                   className="text-muted small mb-1"),
                            _styled_table(rpt),
                        ], className="mb-3"))

                except Exception as exc:
                    children.append(_err(exc))

            # ── BAND SELECTION ────────────────────────────────────────────
            elif action == "band":
                try:
                    keep_parsed = None
                    if keep_str:
                        try:
                            keep_parsed = [float(x.strip())
                                           for x in keep_str.split(",")
                                           if x.strip()]
                        except ValueError:
                            children.append(_warn(
                                "Could not parse 'Keep' Hz list — "
                                "use comma-separated numbers."))
                    edited = select_band(
                        sub,
                        fmin=float(fmin) if fmin not in (None, "") else None,
                        fmax=float(fmax) if fmax not in (None, "") else None,
                        keep=keep_parsed,
                        inplace=False, verbose=0,
                    )
                    cache_set_freq_edit(session_id, edited)
                    has_edit = True
                    children += _freq_before_after_heatmaps(
                        sub, edited, _mpl_wrap, _style_ax, BG, FG, GRID, AX_BG,
                        plot_coverage_quality_heatmap,
                    )
                except Exception as exc:
                    children.append(_err(exc))

            # ── REGRID ────────────────────────────────────────────────────
            elif action == "regrid":
                try:
                    edited = regrid_logspace(
                        sub,
                        fmin=float(rg_fmin) if rg_fmin not in (None,"") else None,
                        fmax=float(rg_fmax) if rg_fmax not in (None,"") else None,
                        per_decade=int(rg_ppd or 6),
                        method=rg_method or "nearest",
                        inplace=False, verbose=0,
                    )
                    cache_set_freq_edit(session_id, edited)
                    has_edit = True
                    children += _freq_before_after_heatmaps(
                        sub, edited, _mpl_wrap, _style_ax, BG, FG, GRID, AX_BG,
                        plot_coverage_quality_heatmap,
                    )
                except Exception as exc:
                    children.append(_err(exc))

            # ── GRID ALIGN ────────────────────────────────────────────────
            elif action == "align":
                try:
                    edited = align_grid(
                        sub,
                        mode=align_mode or "union",
                        ref_station=align_ref or None,
                        method="nearest",
                        inplace=False, verbose=0,
                    )
                    cache_set_freq_edit(session_id, edited)
                    has_edit = True
                    children += _freq_before_after_heatmaps(
                        sub, edited, _mpl_wrap, _style_ax, BG, FG, GRID, AX_BG,
                        plot_coverage_quality_heatmap,
                    )
                except Exception as exc:
                    children.append(_err(exc))

            # ── SMOOTH / DECIMATE ─────────────────────────────────────────
            elif action == "smooth":
                try:
                    op = smooth_op or "smooth"
                    if op == "smooth":
                        edited = smooth_mavg(
                            sub,
                            k=int(smooth_k or 3),
                            on="both", inplace=False, verbose=0,
                        )
                    elif op == "decimate":
                        edited = decimate_step(
                            sub,
                            step=int(smooth_step or 2),
                            inplace=False, verbose=0,
                        )
                    else:
                        edited = drop_duplicates(
                            sub, inplace=False, verbose=0,
                        )
                    cache_set_freq_edit(session_id, edited)
                    has_edit = True
                    children += _freq_before_after_heatmaps(
                        sub, edited, _mpl_wrap, _style_ax, BG, FG, GRID, AX_BG,
                        plot_coverage_quality_heatmap,
                    )
                except Exception as exc:
                    children.append(_err(exc))

            if not children:
                children.append(_warn("No output produced."))

            new_saved = dict(saved_outputs or {})
            new_saved["freq-editor"] = {"type": "html", "imgs": _freq_imgs}
            return html.Div(children), new_saved, not has_edit

        except Exception as exc:
            return _err(exc), no_update, True

    # Frequency editor — export edited EDIs
    @app.callback(
        Output("tool-freq-export-out", "children"),
        Input("tool-freq-export-btn",     "n_clicks"),
        State("tool-freq-out-path",       "value"),
        State("tool-freq-out-template",   "value"),
        State(IDs.SESSION_ID,             "data"),
        prevent_initial_call=True,
    )
    def _export_freq_edis(n, out_path, template, session_id):
        if not n:
            return no_update
        if not out_path:
            return html.Span("⚠ No output folder specified.", className="text-warning")
        try:
            import pathlib

            from pycsamt.app.web.cache import (
                cache_get_freq_edit,
            )
            edited = cache_get_freq_edit(session_id)
            if edited is None:
                return html.Span(
                    "⚠ No edited data in cache — run the editor first.",
                    className="text-warning",
                )
            p = pathlib.Path(out_path)
            p.mkdir(parents=True, exist_ok=True)
            tmpl = (template or "{station}.edi").strip() or "{station}.edi"
            written = edited.write(str(p), template=tmpl, exist_ok=True)
            n_ok = len(written)
            return html.Span(
                f"✓ Exported {n_ok} EDI file(s) → {p}",
                className="text-success",
            )
        except Exception as exc:
            return _err(exc)

    # Batch export
    @app.callback(
        Output("tool-batch-out", "children"),
        Input("tool-batch-run",   "n_clicks"),
        State("tool-batch-path",  "value"),
        State("tool-batch-fmt",   "value"),
        State("tool-batch-dpi",   "value"),
        State("tool-batch-items", "value"),
        State("tool-batch-flags", "value"),
        State(IDs.STORE_DATA,     "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.SESSION_ID,     "data"),
        prevent_initial_call=True,
    )
    def _run_batch(
        n, out_path, fmt, dpi, items, flags,
        store_data, active_lines_store, session_id,
    ):
        if not n:
            return no_update
        if not out_path or not out_path.strip():
            return _warn("Enter an output folder first.")
        if not items:
            return _warn("Select at least one figure type.")
        try:
            from pathlib import Path
            sites, err = _scoped_sites(
                session_id, store_data, active_lines_store,
            )
            if err:
                return _warn(err)
            dest  = Path(out_path.strip())
            dest.mkdir(parents=True, exist_ok=True)
            result = _batch_export_tool_figures(
                sites,
                dest=dest,
                fmt=str(fmt or "png").lower(),
                dpi=int(dpi or 150),
                items=list(items or []),
                flags=set(flags or []),
                store_data=store_data or {},
                active_lines_store=active_lines_store or {},
            )
            return _batch_result_view(result)
        except Exception as exc:
            return _err(exc)

    # ── Elevation enrichment — clientside: mode toggle ────────────────────
    app.clientside_callback(
        """
        function(n0, n1, n2, cur) {
            const trig = dash_clientside.callback_context.triggered_id;
            let mode = cur || 'fetch';
            if (trig === 'tool-elev-mode-fetch')   mode = 'fetch';
            if (trig === 'tool-elev-mode-upload')  mode = 'upload';
            if (trig === 'tool-elev-mode-correct') mode = 'correct';
            const show = m => mode === m ? {} : {display:'none'};
            const btn  = m => mode === m
                ? 'elev-mode-btn active' : 'elev-mode-btn';
            return [mode,
                    show('fetch'), show('upload'), show('correct'),
                    btn('fetch'),  btn('upload'),  btn('correct')];
        }
        """,
        Output("tool-elev-mode-store",      "data"),
        Output("tool-elev-fetch-panel",     "style"),
        Output("tool-elev-upload-panel",    "style"),
        Output("tool-elev-correct-panel",   "style"),
        Output("tool-elev-mode-fetch",      "className"),
        Output("tool-elev-mode-upload",     "className"),
        Output("tool-elev-mode-correct",    "className"),
        Input("tool-elev-mode-fetch",       "n_clicks"),
        Input("tool-elev-mode-upload",      "n_clicks"),
        Input("tool-elev-mode-correct",     "n_clicks"),
        State("tool-elev-mode-store",       "data"),
        prevent_initial_call=True,
    )

    # ── Clientside: scope toggle (All / Active Lines / Selected) ──────────
    app.clientside_callback(
        """
        function(n0, n1, n2, cur) {
            const trig = dash_clientside.callback_context.triggered_id;
            let scope = cur || 'all';
            if (trig === 'tool-elev-scope-all')   scope = 'all';
            if (trig === 'tool-elev-scope-lines')  scope = 'lines';
            if (trig === 'tool-elev-scope-sel')    scope = 'sel';
            const showF = scope !== 'all' ? {} : {display:'none'};
            const btn   = s => scope === s
                ? 'elev-scope-btn active' : 'elev-scope-btn';
            return [scope, showF, btn('all'), btn('lines'), btn('sel')];
        }
        """,
        Output("tool-elev-scope-store",    "data"),
        Output("tool-elev-scope-filters",  "style"),
        Output("tool-elev-scope-all",      "className"),
        Output("tool-elev-scope-lines",    "className"),
        Output("tool-elev-scope-sel",      "className"),
        Input("tool-elev-scope-all",       "n_clicks"),
        Input("tool-elev-scope-lines",     "n_clicks"),
        Input("tool-elev-scope-sel",       "n_clicks"),
        State("tool-elev-scope-store",     "data"),
        prevent_initial_call=True,
    )

    # ── Clientside: correction-ops checklist → show/hide sub-params ───────
    app.clientside_callback(
        """
        function(ops) {
            ops = ops || [];
            const show = k => ops.includes(k) ? {} : {display:'none'};
            return [show('smooth'), show('shift'), show('fill'), show('outlier')];
        }
        """,
        Output("tool-elev-smooth-params",  "style"),
        Output("tool-elev-shift-params",   "style"),
        Output("tool-elev-fill-params",    "style"),
        Output("tool-elev-outlier-params", "style"),
        Input("tool-elev-corr-ops",        "value"),
        prevent_initial_call=False,
    )

    # ── Update station dropdown when line filter changes ───────────────────
    @app.callback(
        Output("tool-elev-stations", "options"),
        Input("tool-elev-lines",     "value"),
        State(IDs.STORE_DATA,        "data"),
        prevent_initial_call=True,
    )
    def _elev_update_stations(sel_lines, store_data):
        records = (store_data or {}).get("station_records", [])
        sel_set = set(sel_lines or [])
        opts = []
        for r in records:
            sid = r.get("ID") or r.get("station") or r.get("name")
            ln  = r.get("Line", "")
            if not sid:
                continue
            if sel_set and ln not in sel_set:
                continue
            opts.append({"label": f"{sid} · {ln}" if ln else sid, "value": sid})
        return opts

    # ── Parse uploaded elevation CSV ───────────────────────────────────────
    @app.callback(
        Output("tool-elev-upload-store", "data"),
        Output("tool-elev-upload-info",  "children"),
        Output("tool-elev-col-id",       "options"),
        Output("tool-elev-col-lat",      "options"),
        Output("tool-elev-col-lon",      "options"),
        Output("tool-elev-col-elev",     "options"),
        Output("tool-elev-col-id",       "value"),
        Output("tool-elev-col-lat",      "value"),
        Output("tool-elev-col-lon",      "value"),
        Output("tool-elev-col-elev",     "value"),
        Input("tool-elev-upload-input",  "contents"),
        State("tool-elev-upload-input",  "filename"),
        prevent_initial_call=True,
    )
    def _elev_parse_upload(contents, filename):
        import base64
        import io

        import pandas as pd
        if not contents:
            return (no_update,) * 10
        try:
            _, enc = contents.split(",", 1)
            raw = base64.b64decode(enc).decode(errors="replace")
            sep = "," if "," in raw.split("\n")[0] else "\t"
            df  = pd.read_csv(io.StringIO(raw), sep=sep)
            cols = list(df.columns)
            opts = [{"label": c, "value": c} for c in cols]

            def _match(*keywords):
                for k in keywords:
                    for c in cols:
                        if k.lower() in c.lower():
                            return c
                return None

            v_id   = _match("station", "name", "id", "site")
            v_lat  = _match("lat", "latitude", "y")
            v_lon  = _match("lon", "long", "longitude", "x")
            v_elev = _match("elev", "alt", "height", "z", "dem")

            info = html.Span(
                [
                    html.I(className="bi bi-check-circle-fill me-1",
                           style={"color": "var(--green)"}),
                    f"{filename}  —  {len(df)} rows, "
                    f"{len(cols)} columns: {', '.join(cols[:6])}"
                    + ("…" if len(cols) > 6 else ""),
                ],
                className="text-success small",
            )
            data = {"records": df.to_dict("records"), "columns": cols}
            return (data, info,
                    opts, opts, opts, opts,
                    v_id, v_lat, v_lon, v_elev)
        except Exception as exc:
            err = html.Span(f"Parse error: {exc}", className="text-danger small")
            return None, err, [], [], [], [], None, None, None, None

    # ── Load file → push into raw-store ───────────────────────────────────
    @app.callback(
        Output("tool-elev-raw-store",       "data",     allow_duplicate=True),
        Output("tool-elev-out",             "children", allow_duplicate=True),
        Input("tool-elev-load-btn",         "n_clicks"),
        State("tool-elev-upload-store",     "data"),
        State("tool-elev-col-id",           "value"),
        State("tool-elev-col-lat",          "value"),
        State("tool-elev-col-lon",          "value"),
        State("tool-elev-col-elev",         "value"),
        prevent_initial_call=True,
    )
    def _elev_load_file(n, upload_data, col_id, col_lat, col_lon, col_elev):
        import pandas as pd
        if not n or not upload_data:
            return no_update, no_update
        try:
            df = pd.DataFrame(upload_data["records"])
            cols = upload_data["columns"]

            def _auto(keywords):
                for k in keywords:
                    for c in cols:
                        if k.lower() in c.lower():
                            return c
                return None

            cid   = col_id   or _auto(["station", "name", "id", "site"])
            clat  = col_lat  or _auto(["lat", "latitude", "y"])
            clon  = col_lon  or _auto(["lon", "long", "longitude", "x"])
            celev = col_elev or _auto(["elev", "alt", "height", "z", "dem"])

            if not celev:
                return no_update, _warn("Could not find an elevation column.")

            rows = []
            for _, r in df.iterrows():
                rows.append({
                    "station": str(r[cid])   if cid  and cid  in r else f"row{_}",
                    "lat":     float(r[clat])  if clat  and clat  in r else None,
                    "lon":     float(r[clon])  if clon  and clon  in r else None,
                    "elev":    float(r[celev]),
                })

            out = html.Div([
                html.Span(
                    [html.I(className="bi bi-check-circle-fill me-1",
                            style={"color": "var(--green)"}),
                     f" {len(rows)} station(s) loaded. "
                     "Switch to Correct & Export to proceed."],
                    className="small",
                ),
            ])
            return rows, out
        except Exception as exc:
            return no_update, _err(exc)

    # ── Fetch elevations from API or session ───────────────────────────────
    @app.callback(
        Output("tool-elev-out",            "children", allow_duplicate=True),
        Output("tool-elev-raw-store",      "data",     allow_duplicate=True),
        Output("tool-elev-progress",       "value"),
        Output("tool-elev-progress-label", "children"),
        Output("tool-elev-progress-wrap",  "style"),
        Input("tool-elev-run",             "n_clicks"),
        State("tool-elev-source",          "value"),
        State("tool-elev-datum",           "value"),
        State("tool-elev-scope-store",     "data"),
        State("tool-elev-lines",           "value"),
        State("tool-elev-stations",        "value"),
        State(IDs.STORE_DATA,              "data"),
        State(IDs.SESSION_ID,              "data"),
        prevent_initial_call=True,
    )
    def _elev_fetch(n, source, datum, scope, sel_lines, sel_stations,
                    store_data, session_id):
        import pandas as pd
        if not n:
            return (no_update,) * 5
        show = {}
        hide = {"display": "none"}
        try:
            from pycsamt.app.web.cache import cache_get
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data."), no_update, 0, "", hide

            records = list((store_data or {}).get("station_records", []) or [])
            sel_line_set = set(sel_lines or [])
            sel_sta_set  = set(sel_stations or [])

            # Filter records by scope
            if scope == "lines" and sel_line_set:
                records = [r for r in records if r.get("Line", "") in sel_line_set]
            elif scope == "sel":
                if sel_sta_set:
                    records = [r for r in records
                               if (r.get("ID") or r.get("station") or r.get("name", ""))
                               in sel_sta_set]
                elif sel_line_set:
                    records = [r for r in records if r.get("Line", "") in sel_line_set]

            if not records:
                return _warn("No stations matched the scope."), no_update, 0, "", hide

            rows_out = []
            total = len(records)

            def _get(r, *keys):
                """Try several key names (case variants) and return first non-None."""
                for k in keys:
                    v = r.get(k)
                    if v is not None:
                        try:
                            f = float(v)
                            if not (f != f):   # exclude NaN
                                return f
                        except (TypeError, ValueError):
                            pass
                return None

            if source == "session":
                # Pull elevation from already-loaded EDI headers
                for r in records:
                    sid = r.get("ID") or r.get("station") or r.get("name", "")
                    lat = _get(r, "Latitude",  "lat", "latitude")
                    lon = _get(r, "Longitude", "lon", "longitude")
                    elv = _get(r, "Elevation", "elev", "elevation")
                    rows_out.append({
                        "station": sid,
                        "lat":  lat,
                        "lon":  lon,
                        "elev": elv,
                    })
                n_ok = sum(1 for r in rows_out if r["elev"] is not None)
                label = f"{n_ok}/{total} stations have elevation in EDI headers"
            else:
                from pycsamt.gis.utils import (
                    get_elevation_from_api,
                )
                lats, lons, sids = [], [], []
                for r in records:
                    sid = r.get("ID") or r.get("station") or r.get("name", "")
                    lat = _get(r, "Latitude",  "lat", "latitude")
                    lon = _get(r, "Longitude", "lon", "longitude")
                    if lat is None or lon is None:
                        rows_out.append({
                            "station": sid, "lat": None,
                            "lon": None,   "elev": None,
                        })
                        continue
                    lats.append(lat)
                    lons.append(lon)
                    sids.append(sid)

                if lats:
                    _api_map = {
                        "open-elevation": "open_elevation",
                        "open-meteo":     "open_meteo",
                        "nasadem":        "nasadem",
                        "open-topo":      "open_topo_data",
                    }
                    api_name = _api_map.get(source)
                    import numpy as np
                    elevs = get_elevation_from_api(
                        np.array(lats), np.array(lons),
                        api_name=api_name,
                    )
                    for sid, lat, lon, elv in zip(sids, lats, lons,
                                                   np.atleast_1d(elevs)):
                        rows_out.append({
                            "station": sid,
                            "lat":     lat,
                            "lon":     lon,
                            "elev":    float(elv) if elv is not None else None,
                        })
                n_ok = sum(1 for r in rows_out if r.get("elev") is not None)
                label = f"{n_ok}/{total} stations fetched"

            df = pd.DataFrame(rows_out)
            out = _elev_summary_table(df, "Fetched elevation")
            return out, rows_out, 100, label, show

        except Exception as exc:
            return _err(exc), no_update, 0, "", hide

    # ── Apply corrections → profile preview ───────────────────────────────
    @app.callback(
        Output("tool-elev-out",            "children", allow_duplicate=True),
        Output("tool-elev-corrected-store","data"),
        Output("tool-elev-progress",       "value",    allow_duplicate=True),
        Output("tool-elev-progress-label", "children", allow_duplicate=True),
        Output("tool-elev-progress-wrap",  "style",    allow_duplicate=True),
        Input("tool-elev-correct-btn",     "n_clicks"),
        State("tool-elev-raw-store",       "data"),
        State("tool-elev-corr-ops",        "value"),
        State("tool-elev-smooth-method",   "value"),
        State("tool-elev-smooth-window",   "value"),
        State("tool-elev-shift-delta",     "value"),
        State("tool-elev-fill-zeros",      "value"),
        State("tool-elev-outlier-sigma",   "value"),
        prevent_initial_call=True,
    )
    def _elev_correct(n, raw_data, ops, sm_method, sm_window,
                      shift_delta, fill_zeros, outlier_sigma):
        import numpy as np
        import pandas as pd
        import plotly.graph_objects as go
        if not n or not raw_data:
            return (no_update,) * 5
        ops = ops or []
        show = {}
        hide = {"display": "none"}
        try:
            df = pd.DataFrame(raw_data)
            if "elev" not in df.columns:
                return _warn("No elevation column found in raw data."), \
                       no_update, 0, "", hide

            df["elev"] = pd.to_numeric(df["elev"], errors="coerce")
            # Always float64 — int arrays cannot hold NaN (outlier / fill ops)
            elev_raw  = df["elev"].values.astype(float)

            # Apply corrections in sensible order
            elev_work = elev_raw.copy()

            if "fill" in ops:
                # Linear interpolation of NaN (and zeros if toggled)
                fill_nan = pd.Series(elev_work, dtype=float)
                if "zeros" in (fill_zeros or []):
                    fill_nan = fill_nan.replace(0.0, float("nan"))
                fill_nan = fill_nan.interpolate(
                    method="linear", limit_direction="both"
                )
                elev_work = fill_nan.values.astype(float)

            if "outlier" in ops:
                sigma = float(outlier_sigma or 3.0)
                med = float(np.nanmedian(elev_work))
                std = float(np.nanstd(elev_work))
                mask = np.abs(elev_work - med) > sigma * std
                elev_work = elev_work.astype(float)   # ensure float before NaN assign
                elev_work[mask] = float("nan")
                s = pd.Series(elev_work, dtype=float)
                elev_work = s.interpolate(
                    method="linear", limit_direction="both"
                ).values.astype(float)
                n_out = int(mask.sum())
            else:
                n_out = 0

            if "smooth" in ops:
                method = sm_method or "loess"
                window = max(3, int(sm_window or 5))
                valid  = ~np.isnan(elev_work)
                x_idx  = np.arange(len(elev_work), dtype=float)
                if valid.sum() >= 3:
                    if method == "loess":
                        from pycsamt.gis.coord_correction import (
                            _loess_smooth,
                        )
                        elev_work[valid] = _loess_smooth(
                            x_idx[valid], elev_work[valid], span=window // 2
                        )
                    elif method == "mean":
                        from scipy.ndimage import (
                            uniform_filter1d,
                        )
                        elev_work = uniform_filter1d(
                            np.where(np.isnan(elev_work), 0, elev_work),
                            size=window,
                        )
                    elif method == "savgol":
                        from scipy.signal import savgol_filter
                        polyorder = min(3, window - 1)
                        elev_work = savgol_filter(
                            np.where(np.isnan(elev_work), 0, elev_work),
                            window_length=window,
                            polyorder=polyorder,
                        )
                    elif method == "gaussian":
                        from scipy.ndimage import (
                            gaussian_filter1d,
                        )
                        elev_work = gaussian_filter1d(
                            np.where(np.isnan(elev_work), 0, elev_work),
                            sigma=window / 3.0,
                        )
                    elif method == "median":
                        from scipy.ndimage import (
                            median_filter,
                        )
                        elev_work = median_filter(
                            np.where(np.isnan(elev_work), 0, elev_work),
                            size=window,
                        ).astype(float)

            if "shift" in ops:
                delta = float(shift_delta or 0.0)
                elev_work = elev_work + delta

            # Build corrected records
            df["elev_corrected"] = elev_work
            corrected_rows = df.to_dict("records")

            # Build profile figure
            x = np.arange(len(elev_raw))
            labels = df.get("station", pd.Series(range(len(df)))).tolist()

            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=x, y=elev_raw,
                mode="lines+markers",
                name="Raw elevation",
                line={"color": "#89b4fa", "width": 1.5, "dash": "dot"},
                marker={"size": 4},
                hovertemplate="%{text}<br>Raw: %{y:.1f} m<extra></extra>",
                text=labels,
            ))
            fig.add_trace(go.Scatter(
                x=x, y=elev_work,
                mode="lines+markers",
                name="Corrected",
                line={"color": "#a6e3a1", "width": 2},
                marker={"size": 5},
                hovertemplate="%{text}<br>Corrected: %{y:.1f} m<extra></extra>",
                text=labels,
            ))
            fig.update_layout(
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                margin={"l": 40, "r": 10, "t": 28, "b": 36},
                height=220,
                legend={"orientation": "h", "y": 1.12, "x": 0,
                        "font": {"size": 11}},
                xaxis={
                    "title": "Station index",
                    "gridcolor": "rgba(127,127,127,0.15)",
                    "tickfont": {"size": 10},
                },
                yaxis={
                    "title": "Elevation (m)",
                    "gridcolor": "rgba(127,127,127,0.15)",
                    "tickfont": {"size": 10},
                },
                font={"color": "#cdd6f4"},
            )

            n_changed = int(np.sum(np.abs(elev_work - elev_raw) > 0.01))
            parts = []
            if "outlier" in ops:
                parts.append(f"{n_out} outlier(s) removed")
            if "smooth" in ops:
                parts.append(f"smoothed ({sm_method}, w={sm_window})")
            if "shift" in ops and shift_delta:
                parts.append(f"shifted {shift_delta:+.1f} m")
            if "fill" in ops:
                parts.append("NaN filled")
            summary = "  ·  ".join(parts) if parts else "No ops applied"

            out = html.Div([
                html.Div(
                    [
                        html.I(className="bi bi-check-circle-fill me-1",
                               style={"color": "var(--green)"}),
                        f" {n_changed} point(s) modified  ·  {summary}",
                    ],
                    className="small mb-2",
                ),
                dcc.Graph(
                    figure=fig,
                    config=_plotly_config("elevation_profile"),
                    style={"height": "220px"},
                ),
            ])
            label = f"{len(df)} stations · {summary}"
            return out, corrected_rows, 100, label, show
        except Exception as exc:
            return _err(exc), no_update, 0, "", hide

    # ── Export CSV ────────────────────────────────────────────────────────
    @app.callback(
        Output("tool-elev-dl-csv",         "data"),
        Output("tool-elev-export-status",  "children", allow_duplicate=True),
        Input("tool-elev-export-csv",      "n_clicks"),
        State("tool-elev-corrected-store", "data"),
        State("tool-elev-raw-store",       "data"),
        prevent_initial_call=True,
    )
    def _elev_export_csv(n, corrected, raw):
        import pandas as pd
        if not n:
            return no_update, no_update
        src = corrected or raw
        if not src:
            return no_update, _warn("No elevation data — fetch or load first.")
        try:
            df = pd.DataFrame(src)
            col_order = [c for c in
                         ["station", "lat", "lon", "elev", "elev_corrected"]
                         if c in df.columns]
            df = df[col_order]
            csv_str = df.to_csv(index=False)
            status = html.Span(
                [html.I(className="bi bi-check-circle-fill me-1",
                        style={"color": "var(--green)"}),
                 f" CSV exported  ({len(df)} rows)."],
                className="small",
            )
            return dcc.send_string(csv_str, "pycsamt_elevation.csv"), status
        except Exception as exc:
            return no_update, _err(exc)

    # ── Export HDF5 / NPZ ─────────────────────────────────────────────────
    @app.callback(
        Output("tool-elev-dl-h5",          "data"),
        Output("tool-elev-export-status",  "children", allow_duplicate=True),
        Input("tool-elev-export-h5",       "n_clicks"),
        State("tool-elev-corrected-store", "data"),
        State("tool-elev-raw-store",       "data"),
        prevent_initial_call=True,
    )
    def _elev_export_h5(n, corrected, raw):
        import base64
        import io
        import os
        import tempfile

        import numpy as np
        import pandas as pd
        if not n:
            return no_update, no_update
        src = corrected or raw
        if not src:
            return no_update, _warn("No elevation data — fetch or load first.")
        try:
            df = pd.DataFrame(src)
            stations = df.get("station", pd.Series(range(len(df)))).astype(str).values
            lats   = pd.to_numeric(df.get("lat",  pd.Series([None]*len(df))),
                                   errors="coerce").values
            lons   = pd.to_numeric(df.get("lon",  pd.Series([None]*len(df))),
                                   errors="coerce").values
            elev_r = pd.to_numeric(df.get("elev", pd.Series([None]*len(df))),
                                   errors="coerce").values
            elev_c = pd.to_numeric(df.get("elev_corrected",
                                          pd.Series([None]*len(df))),
                                   errors="coerce").values

            try:
                import h5py
                with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as tf:
                    tmp_path = tf.name
                with h5py.File(tmp_path, "w") as hf:
                    grp = hf.create_group("elevation")
                    grp.create_dataset("station",        data=stations.astype("S"))
                    grp.create_dataset("lat",            data=lats)
                    grp.create_dataset("lon",            data=lons)
                    grp.create_dataset("elev_raw",       data=elev_r)
                    grp.create_dataset("elev_corrected", data=elev_c)
                    hf.attrs["source"] = "pycsamt elevation enrichment"
                with open(tmp_path, "rb") as fh:
                    b64 = base64.b64encode(fh.read()).decode()
                os.unlink(tmp_path)
                fname = "pycsamt_elevation.h5"
                mimetype = "application/x-hdf"
            except ImportError:
                buf = io.BytesIO()
                np.savez(buf,
                         station=stations, lat=lats, lon=lons,
                         elev_raw=elev_r, elev_corrected=elev_c)
                b64 = base64.b64encode(buf.getvalue()).decode()
                fname = "pycsamt_elevation.npz"
                mimetype = "application/octet-stream"

            status = html.Span(
                [html.I(className="bi bi-check-circle-fill me-1",
                        style={"color": "var(--green)"}),
                 f" {fname} exported ({len(df)} stations)."],
                className="small",
            )
            return (
                {"base64": True, "content": b64,
                 "filename": fname, "type": mimetype},
                status,
            )
        except Exception as exc:
            return no_update, _err(exc)

    # ── Update EDI session headers with corrected elevations ──────────────
    @app.callback(
        Output("tool-elev-export-status",  "children", allow_duplicate=True),
        Input("tool-elev-update-session",  "n_clicks"),
        State("tool-elev-corrected-store", "data"),
        State(IDs.SESSION_ID,              "data"),
        prevent_initial_call=True,
    )
    def _elev_update_session(n, corrected, session_id):
        if not n or not corrected:
            return no_update
        try:
            from pycsamt.app.web.cache import (
                cache_get,
                cache_set,
            )
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data.")

            updated = 0
            for row in corrected:
                sid  = str(row.get("station", ""))
                elev = row.get("elev_corrected") or row.get("elev")
                site = sites.get(sid)
                if elev is None or site is None:
                    continue
                edi  = getattr(site, "edi", site)
                if hasattr(edi, "Head"):
                    edi.Head.update({"elevation": float(elev)})
                elif hasattr(edi, "elevation"):
                    edi.elevation = float(elev)
                updated += 1

            cache_set(session_id, sites)
            return html.Span(
                [html.I(className="bi bi-check-circle-fill me-1",
                        style={"color": "var(--green)"}),
                 f" {updated} station(s) updated in session."],
                className="small",
            )
        except Exception as exc:
            return _err(exc)

    # ── Restore last elevation result when switching back to the tool ────
    @app.callback(
        Output("tool-elev-out",   "children", allow_duplicate=True),
        Output("tool-elev-progress-wrap",  "style",    allow_duplicate=True),
        Output("tool-elev-progress",       "value",    allow_duplicate=True),
        Output("tool-elev-progress-label", "children", allow_duplicate=True),
        Input(IDs.ACTIVE_TOOL,            "data"),
        State("tool-elev-raw-store",       "data"),
        State("tool-elev-corrected-store", "data"),
        prevent_initial_call=True,
    )
    def _restore_elev_output(tool_id, raw, corrected):
        import pandas as pd
        if tool_id != "elevation":
            return (no_update,) * 4
        src = corrected or raw
        if not src:
            return (no_update,) * 4
        df    = pd.DataFrame(src)
        label = ("Corrected elevation (restored)"
                 if corrected else "Fetched elevation (restored)")
        n_ok  = int(df["elev"].notna().sum()) if "elev" in df.columns else 0
        prog_label = f"{n_ok}/{len(df)} stations"
        return _elev_summary_table(df, label), {}, 100, prog_label

    # ── Restore last coord result when switching back to coords tool ──────
    @app.callback(
        Output("tool-coord-out", "children", allow_duplicate=True),
        Input(IDs.ACTIVE_TOOL,             "data"),
        State("tool-coord-result-store",   "data"),
        prevent_initial_call=True,
    )
    def _restore_coord_output(tool_id, result_data):
        import pandas as pd
        if tool_id != "coords" or not result_data:
            return no_update
        rows = result_data if isinstance(result_data, list) else []
        if not rows:
            return no_update
        df = pd.DataFrame(rows)
        tbl = dash_table.DataTable(
            data=df.to_dict("records"),
            columns=[{"name": c, "id": c} for c in df.columns],
            page_size=10,
            style_table={"overflowX": "auto", "fontSize": "11px"},
            style_cell={"padding": "3px 8px",
                        "border": "1px solid var(--srf1)",
                        "background": "var(--base)", "color": "var(--text)"},
            style_header={"background": "var(--srf0)", "fontWeight": "600",
                          "border": "1px solid var(--srf1)"},
        )
        return html.Div([
            html.Div(
                [html.I(className="bi bi-arrow-clockwise me-1",
                        style={"color": "var(--sky)"}),
                 f" {len(df)} result(s) restored from last run."],
                className="small mb-2",
            ),
            tbl,
        ])

    # Station response inspector
    @app.callback(
        Output("tool-resp-out",       "children"),
        Output("tool-outputs-store",  "data", allow_duplicate=True),
        Input("tool-resp-run",        "n_clicks"),
        State("tool-resp-station",    "value"),
        State("tool-resp-station2",   "value"),
        State("tool-resp-comps",      "value"),
        State("tool-resp-tipper",     "value"),
        State("tool-resp-ytype",      "value"),
        State("tool-resp-xscale",     "value"),
        State("tool-resp-period-lo",  "value"),
        State("tool-resp-period-hi",  "value"),
        State("tool-resp-errbars",    "value"),
        State(IDs.SESSION_ID,         "data"),
        State(IDs.STORE_THEME,        "data"),
        State("tool-outputs-store",   "data"),
        prevent_initial_call=True,
    )
    def _run_response(n, station, station2, comps, tipper_comps, ytype, xscale,
                      period_lo, period_hi, show_err, session_id, theme, saved_outputs):
        if not n or not station:
            return no_update, no_update
        try:
            import numpy as np
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots

            from pycsamt.app.web.cache import cache_get

            # Apparent resistivity from the STORED impedance, which pycsamt
            # keeps in field units (mV/km / nT). The field-unit formula is
            # ρ_a = 0.2 · |Z|² / f  (the 0.2 already folds in μ₀ and the unit
            # conversion). Do NOT use the SI form |Z|²/(ω·μ₀) here — applied to
            # field-unit Z it over-estimates ρ_a by 1/(0.2·2π·μ₀) ≈ 6.3·10⁵.
            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data."), no_update

            dark = (theme or "dark") == "dark"
            BG   = "#1e1e2e" if dark else "#eff1f5"
            FG   = "#cdd6f4" if dark else "#4c4f69"
            GRID = "#313244" if dark else "#ccd0da"

            COLORS = (
                {"xx": "#f38ba8", "xy": "#89b4fa", "yx": "#a6e3a1",
                 "yy": "#fab387", "det": "#cba6f7",
                 "tx": "#94e2d5", "ty": "#f5c2e7", "mag": "#f9e2af"}
                if dark else
                {"xx": "#d20f39", "xy": "#1e66f5", "yx": "#40a02b",
                 "yy": "#fe640b", "det": "#8839ef",
                 "tx": "#179299", "ty": "#ea76cb", "mag": "#df8e1d"}
            )
            _COMP_IDX = {"xx": (0, 0), "xy": (0, 1), "yx": (1, 0), "yy": (1, 1)}

            comps         = comps or ["xy", "yx"]
            tipper_comps  = tipper_comps or []
            has_tip_req   = bool(tipper_comps)

            # ── extract site data ──────────────────────────────────────────
            def _extract(s):
                freqs = np.asarray(s.freq, dtype=float).ravel()
                T     = 1.0 / np.where(freqs != 0, freqs, np.nan)

                try:
                    z  = np.asarray(s.z,     dtype=complex)
                except Exception:
                    z  = None
                try:
                    ze = np.asarray(s.z_err, dtype=float) if s.z_err is not None else None
                except Exception:
                    ze = None

                try:
                    rho = np.asarray(s.rho,   dtype=float)
                except Exception:
                    rho = None
                try:
                    phi = np.asarray(s.phase, dtype=float)
                except Exception:
                    phi = None

                if (rho is None or phi is None) and z is not None:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        rho = 0.2 * np.abs(z)**2 / freqs[:, None, None]
                        phi = np.degrees(np.angle(z))

                try:
                    tip = (np.asarray(s.tipper, dtype=complex)
                           if s.has_component("tipper") else None)
                except Exception:
                    tip = None

                return freqs, T, z, ze, rho, phi, tip

            # ── period mask ────────────────────────────────────────────────
            def _period_mask(T):
                mask = np.isfinite(T)
                if period_lo is not None:
                    try:
                        mask &= T >= float(period_lo)
                    except (TypeError, ValueError):
                        pass
                if period_hi is not None:
                    try:
                        mask &= T <= float(period_hi)
                    except (TypeError, ValueError):
                        pass
                return mask

            # ── error propagation ──────────────────────────────────────────
            def _rho_err_arr(z_c, ze_c, freqs_m):
                if ze_c is None:
                    return None
                with np.errstate(divide="ignore", invalid="ignore"):
                    # d(ρ_a)/d|Z| · σ_Z = 2·ρ_a·σ_Z/|Z| with ρ_a = 0.2|Z|²/f
                    #                    = 0.4·|Z|·σ_Z / f
                    return np.where(np.abs(z_c) > 0,
                                    0.4 * np.abs(z_c) * ze_c / freqs_m,
                                    np.nan)

            def _phi_err_arr(z_c, ze_c):
                if ze_c is None:
                    return None
                with np.errstate(divide="ignore", invalid="ignore"):
                    return np.where(np.abs(z_c) > 0,
                                    np.degrees(ze_c / np.abs(z_c)), np.nan)

            # ── gather stations ────────────────────────────────────────────
            s1 = sites.get(station)
            if s1 is None:
                return _warn(f"Station '{station}' not found in session."), no_update
            f1, T1, z1, ze1, rho1, phi1, tip1 = _extract(s1)
            m1 = _period_mask(T1)

            s2 = sites.get(station2) if station2 else None
            if s2 is not None:
                f2, T2, z2, ze2, rho2, phi2, tip2 = _extract(s2)
                m2 = _period_mask(T2)
            else:
                f2 = T2 = z2 = ze2 = rho2 = phi2 = tip2 = m2 = None

            # ── subplot layout ─────────────────────────────────────────────
            want_tip_plot  = has_tip_req and tip1 is not None
            n_imp_rows     = {"rho_phi": 2, "reim": 2, "rho": 1, "phi": 1}.get(ytype, 2)
            total_rows     = n_imp_rows + (1 if want_tip_plot else 0)

            row_heights = (
                [0.38, 0.38, 0.24] if total_rows == 3 else
                [0.5,  0.5]        if total_rows == 2 else
                None
            )
            subplot_titles: list[str] = []
            if ytype in ("rho_phi", "rho"):
                subplot_titles.append("Apparent Resistivity (Ω·m)")
            if ytype in ("rho_phi", "phi"):
                subplot_titles.append("Phase (°)")
            if ytype == "reim":
                subplot_titles += ["Re[Z]  (Ω)", "Im[Z]  (Ω)"]
            if want_tip_plot:
                subplot_titles.append("Tipper")

            fig = make_subplots(
                rows=total_rows, cols=1,
                shared_xaxes=True,
                vertical_spacing=0.07,
                row_heights=row_heights,
                subplot_titles=subplot_titles or None,
            )

            # ── trace builder ──────────────────────────────────────────────
            def _add_imp(T, z, ze, rho, phi, freqs, mask, dash="solid", sta_lbl=""):
                sfx = f" [{sta_lbl}]" if sta_lbl else ""
                for c in comps:
                    col = COLORS.get(c, FG)
                    if c == "det":
                        if z is None:
                            continue
                        z_c  = np.sqrt(z[mask, 0, 1] * z[mask, 1, 0])
                        r_c  = 0.2 * np.abs(z_c)**2 / freqs[mask]
                        p_c  = np.degrees(np.angle(z_c))
                        ze_c = None
                    else:
                        i, j = _COMP_IDX[c]
                        z_c  = z[mask, i, j]   if z   is not None else None
                        ze_c = ze[mask, i, j]  if ze  is not None else None
                        r_c  = rho[mask, i, j] if rho is not None else (
                               0.2 * np.abs(z_c)**2 / freqs[mask]
                               if z_c is not None else None)
                        p_c  = phi[mask, i, j] if phi is not None else (
                               np.degrees(np.angle(z_c)) if z_c is not None else None)

                    T_m = T[mask]

                    if ytype in ("rho_phi", "rho") and r_c is not None:
                        rho_v = np.where(r_c > 0, r_c, np.nan)
                        ekw: dict = {}
                        if show_err and z_c is not None and ze_c is not None:
                            ea = _rho_err_arr(z_c, ze_c, freqs[mask])
                            if ea is not None:
                                ekw = {"error_y": {"type": "data", "array": ea.tolist(),
                                                   "visible": True, "color": col,
                                                   "thickness": 1, "width": 3}}
                        fig.add_trace(go.Scatter(
                            x=T_m.tolist(), y=rho_v.tolist(),
                            name=f"Z{c}{sfx}", mode="lines+markers",
                            line={"color": col, "dash": dash, "width": 1.5},
                            marker={"size": 4}, **ekw,
                        ), row=1, col=1)

                    if ytype in ("rho_phi", "phi") and p_c is not None:
                        phi_row = 2 if ytype == "rho_phi" else 1
                        ekw = {}
                        if show_err and z_c is not None and ze_c is not None:
                            pa = _phi_err_arr(z_c, ze_c)
                            if pa is not None:
                                ekw = {"error_y": {"type": "data", "array": pa.tolist(),
                                                   "visible": True, "color": col,
                                                   "thickness": 1, "width": 3}}
                        fig.add_trace(go.Scatter(
                            x=T_m.tolist(), y=p_c.tolist(),
                            name=f"φ Z{c}{sfx}", mode="lines+markers",
                            line={"color": col, "dash": dash, "width": 1.5},
                            marker={"size": 4},
                            showlegend=(ytype != "rho_phi"),
                            **ekw,
                        ), row=phi_row, col=1)

                    if ytype == "reim" and z_c is not None:
                        fig.add_trace(go.Scatter(
                            x=T_m.tolist(), y=z_c.real.tolist(),
                            name=f"Re Z{c}{sfx}", mode="lines+markers",
                            line={"color": col, "dash": dash, "width": 1.5},
                            marker={"size": 4},
                        ), row=1, col=1)
                        fig.add_trace(go.Scatter(
                            x=T_m.tolist(), y=z_c.imag.tolist(),
                            name=f"Im Z{c}{sfx}", mode="lines+markers",
                            line={"color": col, "dash": "dot" if dash == "solid" else dash,
                                  "width": 1.5},
                            marker={"size": 4}, showlegend=False,
                        ), row=2, col=1)

            def _add_tip(T, tip, mask, dash="solid", sta_lbl=""):
                if tip is None:
                    return
                sfx = f" [{sta_lbl}]" if sta_lbl else ""
                T_m = T[mask]
                for tc in tipper_comps:
                    col = COLORS.get(tc, FG)
                    if   tc == "tx":
                        vals = tip[mask, 0, 0].real
                        lbl  = f"Re Tzx{sfx}"
                    elif tc == "ty":
                        vals = tip[mask, 0, 1].real
                        lbl  = f"Re Tzy{sfx}"
                    elif tc == "mag":
                        vals = np.sqrt(np.abs(tip[mask, 0, 0])**2
                                       + np.abs(tip[mask, 0, 1])**2)
                        lbl  = f"|T|{sfx}"
                    else:
                        continue
                    fig.add_trace(go.Scatter(
                        x=T_m.tolist(), y=vals.tolist(), name=lbl,
                        mode="lines+markers",
                        line={"color": col, "dash": dash, "width": 1.5},
                        marker={"size": 4},
                    ), row=n_imp_rows + 1, col=1)

            # ── populate figure ────────────────────────────────────────────
            _add_imp(T1, z1, ze1, rho1, phi1, f1, m1)
            if want_tip_plot:
                _add_tip(T1, tip1, m1)
            if s2 is not None:
                _add_imp(T2, z2, ze2, rho2, phi2, f2, m2, "dash", station2)
                if want_tip_plot:
                    _add_tip(T2, tip2, m2, "dash", station2)

            # ── axis labels & styling ──────────────────────────────────────
            for r in range(1, total_rows + 1):
                fig.update_xaxes(type=xscale, gridcolor=GRID, zerolinecolor=GRID,
                                 tickfont={"color": FG}, row=r, col=1)
                fig.update_yaxes(gridcolor=GRID, zerolinecolor=GRID,
                                 tickfont={"color": FG}, row=r, col=1)
            if ytype in ("rho_phi", "rho"):
                fig.update_yaxes(title_text="ρ (Ω·m)", type="log", row=1, col=1)
            if ytype in ("rho_phi", "phi"):
                fig.update_yaxes(title_text="Phase (°)",
                                 row=(2 if ytype == "rho_phi" else 1), col=1)
            if ytype == "reim":
                fig.update_yaxes(title_text="Re[Z]", row=1, col=1)
                fig.update_yaxes(title_text="Im[Z]", row=2, col=1)
            if want_tip_plot:
                fig.update_yaxes(title_text="T", row=n_imp_rows + 1, col=1)
            fig.update_xaxes(title_text="Period (s)", row=total_rows, col=1)

            fig.update_layout(
                plot_bgcolor=BG, paper_bgcolor=BG,
                font={"color": FG, "size": 11},
                margin={"l": 65, "r": 20, "t": 35, "b": 50},
                height=max(340, 170 * total_rows),
                legend={"orientation": "h", "y": 1.05, "font": {"size": 10}},
                hovermode="x unified",
            )

            # ── station metadata bar ───────────────────────────────────────
            parts = [f"Station: {station}"]
            for attr in ("lat", "lon", "elev"):
                v = getattr(s1, attr, None)
                if v is not None:
                    try:
                        parts.append(f"{attr}={float(v):.4f}")
                    except (TypeError, ValueError):
                        pass
            parts.append(f"n={int(m1.sum())} periods")
            if s2 is not None:
                parts.append(f"vs {station2}")

            meta_bar = html.Div(
                " · ".join(parts),
                style={"fontSize": "11px", "color": "var(--subtext0)",
                       "fontFamily": "monospace", "marginBottom": "4px"},
            )

            new_saved = dict(saved_outputs or {})
            new_saved["response"] = {"type": "fig_json", "fig": fig.to_json()}
            return (
                html.Div([
                    meta_bar,
                    dcc.Graph(
                        figure=fig,
                        config=_plotly_config("pycsamt_station_response"),
                        className="tool-graph",
                    ),
                ]),
                new_saved,
            )
        except Exception as exc:
            return _err(exc), no_update

    # Strike profile viewer
    @app.callback(
        Output("tool-sprof-out",          "children"),
        Output("tool-outputs-store",      "data", allow_duplicate=True),
        Output("tool-sprof-custom-band",  "style"),
        Input("tool-sprof-run",           "n_clicks"),
        Input("tool-sprof-band",          "value"),
        State("tool-sprof-line",          "value"),
        State("tool-sprof-method",        "value"),
        State("tool-sprof-period-lo",     "value"),
        State("tool-sprof-period-hi",     "value"),
        State("tool-sprof-display",       "value"),
        State("tool-sprof-overlays",      "value"),
        State("tool-sprof-colorby",       "value"),
        State("tool-sprof-iqr-thresh",    "value"),
        State("tool-sprof-skew-thresh",   "value"),
        State(IDs.STORE_DATA,             "data"),
        State(IDs.SESSION_ID,             "data"),
        State(IDs.STORE_THEME,            "data"),
        State("tool-outputs-store",       "data"),
        prevent_initial_call=True,
    )
    def _run_strike_profile(n_clicks, band_choice,
                            lines, method, t_lo, t_hi,
                            display, overlays, colorby,
                            iqr_thresh, skew_thresh,
                            store_data, session_id, theme, saved_outputs):

        # always update the custom-band visibility
        show_custom = {"display": "flex", "gap": "4px", "marginBottom": "8px"} \
                      if band_choice == "custom" \
                      else {"display": "none"}

        if not n_clicks or ctx.triggered_id == "tool-sprof-band":
            return no_update, no_update, show_custom

        try:
            import numpy as np
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots

            from pycsamt.app.web.cache import cache_get
            from pycsamt.emtools.strike import (
                estimate_strike_consensus,
                estimate_strike_phase_tensor,
                estimate_strike_sweep,
            )

            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data."), no_update, show_custom

            dark = (theme or "dark") == "dark"
            BG   = "#1e1e2e" if dark else "#eff1f5"
            FG   = "#cdd6f4" if dark else "#4c4f69"
            GRID = "#313244" if dark else "#ccd0da"

            overlays     = overlays or []
            do_iqr       = "iqr"    in overlays
            do_skew      = "skew"   in overlays
            do_norm      = "norm"   in overlays
            do_flag3d    = "flag3d" in overlays
            do_table     = "table"  in overlays
            iqr_thresh   = float(iqr_thresh  or 30)
            skew_thresh  = float(skew_thresh or 0.2)

            # ── period band ─────────────────────────────────────────────────
            _BAND = {
                "all":     None,
                "shallow": (None, 0.01),
                "mid":     (0.01, 1.0),
                "deep":    (1.0, None),
            }
            if band_choice == "custom":
                lo = float(t_lo) if t_lo is not None else None
                hi = float(t_hi) if t_hi is not None else None
                band = (lo, hi) if (lo is not None or hi is not None) else None
            else:
                band = _BAND.get(band_choice or "all")

            # ── filter sites to selected line(s) ────────────────────────────
            lines = lines or []
            if isinstance(lines, str):
                lines = [lines]
            records = (store_data or {}).get("station_records", [])
            keep    = {r["ID"] for r in records
                       if not lines or "__all__" in lines
                       or r.get("Line", "") in lines}
            filtered_edis = [s.edi for s in sites if s.name in keep] \
                            if keep else [s.edi for s in sites]
            if not filtered_edis:
                return _warn("No stations found for the selected line(s)."), \
                       no_update, show_custom

            from pycsamt.site.base import Sites as _Sites
            sub = _Sites(filtered_edis)

            # ── run estimator ────────────────────────────────────────────────
            _FN = {
                "sweep":        estimate_strike_sweep,
                "bahr":         estimate_strike_sweep,
                "phase_tensor": estimate_strike_phase_tensor,
                "consensus":    estimate_strike_consensus,
            }
            est_fn = _FN.get(method or "sweep", estimate_strike_sweep)
            df     = est_fn(sub, band=band)
            if df.empty:
                return _warn("No results — check period band or data coverage."), \
                       no_update, show_custom

            # normalise angles to 0–90° if requested
            if do_norm:
                df["ang"] = df["ang"].abs()

            # ── optionally get per-station Swift skew ────────────────────────
            skew_map: dict[str, float] = {}
            if do_skew or colorby == "skew" or do_flag3d:
                try:
                    from pycsamt.emtools.anisotropy import (
                        analyze_anisotropy,
                    )
                    adf = analyze_anisotropy(sub)
                    if not adf.empty and "swift_skew" in adf.columns:
                        # apply same period filter as band
                        if band:
                            lo_b = band[0] or 0.0
                            hi_b = band[1] or np.inf
                            m    = ((adf["period_s"] >= lo_b) &
                                    (adf["period_s"] <= hi_b))
                            adf  = adf[m]
                        skew_map = (
                            adf.groupby("station")["swift_skew"]
                            .median().to_dict()
                        )
                except Exception:
                    skew_map = {}

            # ── build colour array ───────────────────────────────────────────
            stations  = df["station"].tolist()
            angles    = df["ang"].tolist()
            iqr_vals  = df["iqr"].tolist()
            n_vals    = df["n"].tolist()
            skew_vals = [skew_map.get(s, np.nan) for s in stations]
            x_idx     = list(range(len(stations)))

            _CMAP_DARK  = "RdBu_r"
            _CMAP_LIGHT = "RdBu"

            if colorby == "iqr":
                c_arr  = iqr_vals
                c_cmap = "YlOrRd"
                c_lbl  = "IQR (°)"
            elif colorby == "skew":
                c_arr  = skew_vals
                c_cmap = "RdYlGn_r"
                c_lbl  = "Swift skew"
            elif colorby == "n":
                c_arr  = n_vals
                c_cmap = "Blues" if not dark else "Blues"
                c_lbl  = "n periods"
            else:
                c_arr  = angles
                c_cmap = _CMAP_DARK if dark else _CMAP_LIGHT
                c_lbl  = "Strike (°)"

            # ── HEATMAP view ─────────────────────────────────────────────────
            if display == "heatmap":
                try:
                    from pycsamt.emtools.anisotropy import (
                        analyze_anisotropy,
                    )
                    adf = analyze_anisotropy(sub)
                    if adf.empty or "strike_deg" not in adf.columns:
                        raise ValueError("no per-freq strike data")
                    if band:
                        lo_b = band[0] or 0.0
                        hi_b = band[1] or np.inf
                        adf  = adf[(adf["period_s"] >= lo_b) &
                                    (adf["period_s"] <= hi_b)]
                    piv = (adf.pivot_table(index="period_s",
                                           columns="station",
                                           values="strike_deg",
                                           aggfunc="median")
                              .sort_index())
                    col_order = [s for s in stations if s in piv.columns]
                    piv       = piv[col_order]
                    zv        = piv.values
                    if do_norm:
                        zv = np.abs(zv)
                    heat_fig = go.Figure(go.Heatmap(
                        z=zv,
                        x=col_order,
                        y=np.log10(piv.index.to_numpy()),
                        colorscale=_CMAP_DARK if dark else _CMAP_LIGHT,
                        zmid=0 if not do_norm else None,
                        colorbar={"title": {"text": "Strike (°)", "side": "right", "font": {"size": 10}}},
                        hovertemplate="Station: %{x}<br>log₁₀T: %{y:.2f}<br>"
                                      "Strike: %{z:.1f}°<extra></extra>",
                    ))
                    heat_fig.update_yaxes(
                        title="log₁₀ Period (s)",
                        tickvals=np.log10([0.001, 0.01, 0.1, 1, 10, 100, 1000]),
                        ticktext=["0.001","0.01","0.1","1","10","100","1000"],
                        gridcolor=GRID,
                    )
                    heat_fig.update_xaxes(tickangle=-45, gridcolor=GRID)
                    heat_fig.update_layout(
                        plot_bgcolor=BG, paper_bgcolor=BG,
                        font={"color": FG, "size": 10},
                        margin={"l": 65, "r": 20, "t": 25, "b": 80},
                        height=340,
                    )
                    new_saved = dict(saved_outputs or {})
                    new_saved["strike-profile"] = {
                        "type": "fig_json", "fig": heat_fig.to_json()}
                    return (
                        dcc.Graph(figure=heat_fig,
                                  config=_plotly_config("pycsamt_strike_heatmap"),
                                  className="tool-graph"),
                        new_saved, show_custom,
                    )
                except Exception as hm_exc:
                    return _err(hm_exc), no_update, show_custom

            # ── PROFILE view ──────────────────────────────────────────────────
            n_rows = 1 + (1 if do_skew and skew_map else 0)
            row_h  = [0.72, 0.28] if n_rows == 2 else None

            fig = make_subplots(
                rows=n_rows, cols=1,
                shared_xaxes=True,
                vertical_spacing=0.06,
                row_heights=row_h,
            )

            # ── IQR ribbon ────────────────────────────────────────────────────
            if do_iqr:
                q_hi = [a + q / 2 for a, q in zip(angles, iqr_vals)]
                q_lo = [a - q / 2 for a, q in zip(angles, iqr_vals)]
                rib_x = x_idx + x_idx[::-1]
                rib_y = q_hi + q_lo[::-1]
                fig.add_trace(go.Scatter(
                    x=rib_x, y=rib_y,
                    fill="toself",
                    fillcolor="rgba(137,180,250,0.12)" if dark
                              else "rgba(30,102,245,0.10)",
                    line={"color": "rgba(0,0,0,0)"},
                    showlegend=True, name="IQR ±½",
                    hoverinfo="skip",
                ), row=1, col=1)

            # ── 3-D background shading ────────────────────────────────────────
            if do_flag3d and skew_map:
                for xi, sta in enumerate(stations):
                    if skew_map.get(sta, 0) > skew_thresh:
                        fig.add_vrect(
                            x0=xi - 0.4, x1=xi + 0.4,
                            fillcolor="rgba(243,139,168,0.12)",
                            line_width=0,
                            annotation_text="3D?",
                            annotation_font_size=8,
                            annotation_position="top left",
                            row=1, col=1,
                        )

            # ── reference lines at 0°, ±45°, ±90° ────────────────────────────
            for yref in ([-90, -45, 0, 45, 90] if not do_norm else [0, 45, 90]):
                fig.add_hline(
                    y=yref, line_dash="dot",
                    line_color="rgba(108,112,134,0.35)",
                    line_width=1, row=1, col=1,
                )

            # ── main scatter — colour-mapped markers ──────────────────────────
            marker_kw: dict = {
                "size": [max(5, min(12, 5 + nv // 4)) for nv in n_vals],
                "colorscale": c_cmap,
                "color":      c_arr,
                "showscale":  True,
                "colorbar":   {"title": {"text": c_lbl, "side": "right", "font": {"size": 9}},
                               "thickness": 10, "len": 0.5, "y": 0.75},
                "line":       {"width": 0.5,
                               "color": "rgba(0,0,0,0.3)"},
            }
            # highlight high-IQR stations
            fig.add_trace(go.Scatter(
                x=x_idx, y=angles,
                mode="lines+markers",
                line={"width": 1.8,
                      "color": "#89b4fa" if dark else "#1e66f5"},
                marker=marker_kw,
                name=f"Strike ({method})",
                customdata=list(zip(stations, iqr_vals, n_vals, skew_vals)),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Strike: %{y:.1f}°<br>"
                    "IQR: ±%{customdata[1]:.1f}°<br>"
                    "n freq: %{customdata[2]}<br>"
                    "Skew: %{customdata[3]:.3f}"
                    "<extra></extra>"
                ),
            ), row=1, col=1)

            # high-IQR markers — red ring overlay
            hi_iqr_x = [xi for xi, q in zip(x_idx, iqr_vals) if q > iqr_thresh]
            hi_iqr_y = [a  for a,  q in zip(angles,iqr_vals) if q > iqr_thresh]
            if hi_iqr_x:
                fig.add_trace(go.Scatter(
                    x=hi_iqr_x, y=hi_iqr_y,
                    mode="markers",
                    marker={"symbol": "circle-open", "size": 14,
                            "color": "#f38ba8", "line": {"width": 2}},
                    name=f"IQR > {iqr_thresh:.0f}°",
                    hoverinfo="skip",
                ), row=1, col=1)

            fig.update_yaxes(
                title_text="Strike angle (°)",
                gridcolor=GRID, zerolinecolor=GRID,
                range=[-95, 95] if not do_norm else [-5, 95],
                row=1, col=1,
            )

            # ── Swift skew subplot ────────────────────────────────────────────
            if do_skew and skew_map:
                sk_y = [skew_map.get(s, np.nan) for s in stations]
                fig.add_hline(
                    y=skew_thresh, line_dash="dash",
                    line_color="#f38ba8", line_width=1.2,
                    annotation_text=f"Skew={skew_thresh}",
                    annotation_font_size=8,
                    row=2, col=1,
                )
                fig.add_trace(go.Bar(
                    x=x_idx, y=sk_y,
                    name="Swift skew",
                    marker_color=[
                        "#f38ba8" if (v or 0) > skew_thresh else
                        ("#89b4fa" if dark else "#1e66f5")
                        for v in sk_y
                    ],
                    opacity=0.75,
                ), row=2, col=1)
                fig.update_yaxes(
                    title_text="Swift skew",
                    gridcolor=GRID,
                    range=[0, max(max(v for v in sk_y if not np.isnan(v)), 0.6)],
                    row=2, col=1,
                )

            # ── x-axis (bottom row) ───────────────────────────────────────────
            fig.update_xaxes(
                tickvals=x_idx,
                ticktext=stations,
                tickangle=-45,
                gridcolor=GRID,
                zerolinecolor=GRID,
                row=n_rows, col=1,
            )
            fig.update_layout(
                plot_bgcolor=BG, paper_bgcolor=BG,
                font={"color": FG, "size": 10},
                margin={"l": 60, "r": 80, "t": 15, "b": 80},
                height=280 + (100 if do_skew and skew_map else 0),
                legend={"orientation": "h", "y": 1.04,
                        "font": {"size": 9}},
                hovermode="x unified",
                bargap=0.25,
            )

            # ── Stats table ───────────────────────────────────────────────────
            table_el = None
            if do_table:
                t_rows = []
                for sta, ang, iqr_v, n_v, sk_v in zip(
                        stations, angles, iqr_vals, n_vals, skew_vals):
                    flag3d = (not np.isnan(sk_v) and sk_v > skew_thresh)
                    iqr_hi = iqr_v > iqr_thresh
                    t_rows.append(html.Tr([
                        html.Td(sta, style={"fontFamily": "monospace",
                                            "fontSize": "11px"}),
                        html.Td(f"{ang:.1f}°",
                                style={"textAlign": "right", "fontSize": "11px",
                                       "fontWeight": "600"}),
                        html.Td(f"±{iqr_v:.1f}°",
                                style={"textAlign": "right", "fontSize": "11px",
                                       "color": "#f38ba8" if iqr_hi else "inherit"}),
                        html.Td(f"{sk_v:.3f}" if not np.isnan(sk_v) else "—",
                                style={"textAlign": "right", "fontSize": "11px",
                                       "color": "#f38ba8" if flag3d else "inherit"}),
                        html.Td(str(int(n_v)),
                                style={"textAlign": "right", "fontSize": "11px",
                                       "color": "var(--subtext0)"}),
                        html.Td(
                            html.Span("3D?",
                                      style={"background": "rgba(243,139,168,.18)",
                                             "color": "#f38ba8",
                                             "borderRadius": "4px",
                                             "padding": "1px 5px",
                                             "fontSize": "9px",
                                             "fontWeight": "700"})
                            if flag3d else "",
                        ),
                    ]))
                table_el = html.Div([
                    html.Div("Per-station summary",
                             style={"fontSize": "10px",
                                    "color": "var(--subtext0)",
                                    "textTransform": "uppercase",
                                    "letterSpacing": "0.05em",
                                    "marginBottom": "4px",
                                    "marginTop": "8px"}),
                    html.Div(
                        html.Table(
                            [html.Thead(html.Tr([
                                html.Th("Station"), html.Th("Strike"),
                                html.Th("IQR"), html.Th("Skew"),
                                html.Th("n"), html.Th(""),
                            ], style={"fontSize": "9px",
                                      "color": "var(--subtext0)",
                                      "textTransform": "uppercase",
                                      "letterSpacing": "0.05em"})),
                             html.Tbody(t_rows)],
                            style={"width": "100%",
                                   "borderCollapse": "collapse"},
                        ),
                        style={"overflowX": "auto",
                               "fontSize": "11px",
                               "color": FG},
                        className="small",
                    ),
                ])

            new_saved = dict(saved_outputs or {})
            new_saved["strike-profile"] = {
                "type": "fig_json", "fig": fig.to_json()}

            return (
                html.Div([
                    dcc.Graph(figure=fig,
                              config=_plotly_config("pycsamt_strike_profile"),
                              className="tool-graph"),
                    table_el or "",
                ]),
                new_saved,
                show_custom,
            )

        except Exception as exc:
            return _err(exc), no_update, show_custom

    # Phase tensor map / pseudo-section
    @app.callback(
        Output("tool-pt-out",        "children"),
        Output("tool-outputs-store", "data", allow_duplicate=True),
        Input("tool-pt-run",         "n_clicks"),
        State("tool-pt-lines",       "value"),
        State("tool-pt-view",        "value"),
        State("tool-pt-period",      "value"),
        State("tool-pt-period-lo",   "value"),
        State("tool-pt-period-hi",   "value"),
        State("tool-pt-colorby",     "value"),
        State("tool-pt-cmap",        "value"),
        State("tool-pt-clim-lo",     "value"),
        State("tool-pt-clim-hi",     "value"),
        State("tool-pt-norm",        "value"),
        State("tool-pt-scale",        "value"),
        State("tool-pt-alpha",        "value"),
        State("tool-pt-lw",           "value"),
        State("tool-pt-yaxis",        "value"),
        State("tool-pt-ref-ellipse",  "value"),
        State("tool-pt-sym-clim",     "value"),
        State("tool-pt-tipper",       "value"),
        State("tool-pt-mark3d",       "value"),
        State("tool-pt-skew-thresh",  "value"),
        State(IDs.STORE_DATA,         "data"),
        State(IDs.SESSION_ID,         "data"),
        State(IDs.STORE_THEME,        "data"),
        State("tool-outputs-store",   "data"),
        prevent_initial_call=True,
    )
    def _run_phase_tensor(n, lines, view, period, t_lo, t_hi,
                          colorby, cmap, clim_lo, clim_hi,
                          norm_by, scale, alpha_v, lw_v, yaxis_v,
                          ref_ellipse, sym_clim,
                          tipper_mode, mark3d, skew_thresh,
                          store_data, session_id, theme, saved_outputs):
        if not n:
            return no_update, no_update
        try:
            import base64
            import io

            import matplotlib
            import numpy as np
            import plotly.graph_objects as go
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            from pycsamt.app.web.cache import cache_get
            from pycsamt.emtools.tensor import (
                build_phase_tensor_table,
                plot_phase_tensor_map,
                plot_phase_tensor_psection,
            )

            sites = cache_get(session_id)
            if sites is None:
                return _warn("Session expired — reload data."), no_update

            dark    = (theme or "dark") == "dark"
            BG      = "#1e1e2e" if dark else "#eff1f5"
            FG      = "#cdd6f4" if dark else "#4c4f69"
            GRID    = "#313244" if dark else "#ccd0da"
            AX_BG   = "#181825" if dark else "#e6e9ef"
            EC      = "0.55"    if dark else "0.3"   # edge colour (greyscale frac)

            # ── filter sites by line ──────────────────────────────────────────
            lines   = lines or []
            if isinstance(lines, str):
                lines = [lines]
            records = (store_data or {}).get("station_records", [])
            keep    = {r["ID"] for r in records
                       if not lines or "__all__" in lines
                       or r.get("Line", "") in lines}
            filtered_edis = [s.edi for s in sites if s.name in keep] \
                            if keep else [s.edi for s in sites]
            if not filtered_edis:
                return _warn("No stations for the selected line(s)."), no_update

            from pycsamt.site.base import Sites as _Sites
            sub = _Sites(filtered_edis)

            period     = float(period or 1.0)
            colorby    = colorby or "skew"
            cmap       = cmap    or "RdBu_r"
            norm_by    = norm_by or "cell"
            scale_v    = float(scale       or 0.85)
            alpha_vv   = float(alpha_v     or 0.92)
            lw_vv      = float(lw_v        or 0.3)
            clim_pct   = (float(clim_lo or 5), float(clim_hi or 95))
            skew_thr   = float(skew_thresh or 5.0)
            show_tip   = tipper_mode not in (None, "none", False)
            tip_conv   = tipper_mode if tipper_mode not in (None, "none", False) \
                         else "parkinson"
            t_lo_v     = float(t_lo) if t_lo is not None else None
            t_hi_v     = float(t_hi) if t_hi is not None else None
            p_range    = (t_lo_v, t_hi_v) \
                         if (t_lo_v is not None or t_hi_v is not None) else None
            # y-axis decode: "logperiod_up" → axis_y="logperiod", period_up=True
            _ymap      = {
                "logperiod_up": ("logperiod", True),
                "logperiod_dn": ("logperiod", False),
                "logfreq_up":   ("logfreq",   True),
            }
            axis_y_str, period_up_v = _ymap.get(yaxis_v or "logperiod_up",
                                                  ("logperiod", True))
            ref_ell    = bool(ref_ellipse) if ref_ellipse is not None else True
            sym_clim_v = bool(sym_clim)    if sym_clim    is not None else True

            view = view or "ellipses"
            children: list = []

            # ── helper: matplotlib fig → base64 PNG ───────────────────────────
            def _mpl_to_img(mpl_fig: plt.Figure) -> html.Img:
                buf = io.BytesIO()
                mpl_fig.savefig(buf, format="png", dpi=120,
                                bbox_inches="tight",
                                facecolor=mpl_fig.get_facecolor())
                buf.seek(0)
                b64 = base64.b64encode(buf.read()).decode()
                plt.close(mpl_fig)
                return html.Img(
                    src=f"data:image/png;base64,{b64}",
                    style={"width": "100%", "borderRadius": "6px",
                           "display": "block"},
                )

            # ── helper: theme matplotlib axes ─────────────────────────────────
            def _style_ax(ax_: plt.Axes, fig_: plt.Figure) -> None:
                fig_.patch.set_facecolor(BG)
                ax_.set_facecolor(AX_BG)
                ax_.tick_params(colors=FG, labelsize=8)
                for sp in ax_.spines.values():
                    sp.set_edgecolor(GRID)
                ax_.title.set_color(FG)
                ax_.xaxis.label.set_color(FG)
                ax_.yaxis.label.set_color(FG)

            # ── ELLIPSE PSEUDO-SECTION (Caldwell 2004 style, matplotlib) ─────────
            if view in ("ellipses", "all"):
                try:
                    fig_w = 10.0
                    fig_h = max(4.5, min(7.0, 0.6 * len(filtered_edis)))
                    mpl_fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                    _style_ax(ax, mpl_fig)

                    # colourbar styling patch: make it theme-aware
                    _orig_colorbar = mpl_fig.colorbar
                    def _themed_cb(mappable, ax_=None, **kw):
                        cb = _orig_colorbar(mappable, ax=ax_, **kw)
                        cb.ax.tick_params(colors=FG, labelsize=7)
                        cb.ax.yaxis.label.set_color(FG)
                        cb.outline.set_edgecolor(GRID)
                        return cb
                    mpl_fig.colorbar = _themed_cb

                    plot_phase_tensor_psection(
                        sub,
                        period_range  = p_range,
                        axis_y        = axis_y_str,
                        period_up     = period_up_v,
                        c_by          = colorby,
                        cmap          = cmap,
                        clim_pct      = clim_pct,
                        symmetric_clim= sym_clim_v,
                        normalise_by  = norm_by,
                        scale         = scale_v,
                        alpha         = alpha_vv,
                        linewidth     = lw_vv,
                        edgecolor     = EC,
                        mark_3d       = bool(mark3d),
                        skew_threshold= skew_thr,
                        ref_ellipse   = ref_ell,
                        tick_label_rotation=45.0,
                        ax            = ax,
                    )
                    _style_ax(ax, mpl_fig)  # re-apply after ellipses drawn

                    # style the colorbar that was added by the function
                    for ax_cb in mpl_fig.axes[1:]:
                        ax_cb.tick_params(colors=FG, labelsize=7)
                        ax_cb.yaxis.label.set_color(FG)
                        for sp in ax_cb.spines.values():
                            sp.set_edgecolor(GRID)

                    children.append(html.Div([
                        html.Div(
                            "Phase-Tensor Ellipse Section  (Caldwell et al. 2004)  —  "
                            f"colour: {colorby}",
                            style={"fontSize": "10px",
                                   "color": "var(--subtext0)",
                                   "textTransform": "uppercase",
                                   "letterSpacing": "0.05em",
                                   "marginBottom": "4px"},
                        ),
                        _mpl_to_img(mpl_fig),
                    ]))
                except Exception as ell_exc:
                    children.append(_warn(
                        f"Ellipse section could not render: {ell_exc}"))

            # ── PSEUDO-SECTION (Plotly heatmap) ───────────────────────────────
            if view in ("psection", "all"):
                df = build_phase_tensor_table(sub)
                if df.empty:
                    children.append(_warn(
                        "No phase tensor data available for these stations."))
                else:
                    # apply period range filter
                    if p_range:
                        lo_p = p_range[0] or 0.0
                        hi_p = p_range[1] or np.inf
                        df   = df[(df["period"] >= lo_p) & (df["period"] <= hi_p)]

                    # choose colour column
                    _COL = {
                        "skew": "skew", "beta": "skew",
                        "theta": "theta", "alpha": "alpha",
                        "ellipt": "ellipt", "ellipticity": "ellipt",
                        "phi_max": "s1", "s1": "s1",
                        "phi_min": "s2", "s2": "s2",
                    }
                    col_key = _COL.get(colorby, "skew")
                    if col_key not in df.columns:
                        col_key = "skew"

                    _CLABEL = {
                        "skew": "Skew β (°)", "theta": "Strike θ (°)",
                        "alpha": "α (°)", "ellipt": "Ellipticity",
                        "s1": "φ_max", "s2": "φ_min",
                    }
                    clabel = _CLABEL.get(col_key, col_key)

                    # station order from records (profile order)
                    sta_order = [r["ID"] for r in records
                                 if not keep or r["ID"] in keep]
                    sta_order = [s for s in sta_order if s in df["station"].values]
                    if not sta_order:
                        sta_order = df["station"].unique().tolist()

                    piv = (df.pivot_table(index="period", columns="station",
                                          values=col_key, aggfunc="median")
                              .sort_index())[sta_order]

                    vals  = piv.values
                    flat  = vals[np.isfinite(vals)]
                    vlo   = float(np.nanpercentile(flat, clim_pct[0])) if flat.size else -10
                    vhi   = float(np.nanpercentile(flat, clim_pct[1])) if flat.size else  10
                    # symmetric for skew-like
                    if col_key in ("skew", "theta", "alpha") and vlo < 0 < vhi:
                        vm = max(abs(vlo), abs(vhi))
                        vlo, vhi = -vm, vm

                    # Plotly colourscale from matplotlib cmap
                    import matplotlib.cm as _mcm
                    import matplotlib.colors as _mco
                    _cm   = _mcm.get_cmap(cmap, 256)
                    _cs   = [[i / 255, _mco.to_hex(_cm(i))]
                             for i in range(256)]

                    log_periods = np.log10(piv.index.to_numpy())
                    heat = go.Heatmap(
                        z=vals,
                        x=sta_order,
                        y=log_periods,
                        colorscale=_cs,
                        zmin=vlo, zmax=vhi,
                        colorbar={"title": {"text": clabel, "side": "right",
                                             "font": {"size": 9, "color": FG}},
                                  "tickfont": {"size": 8, "color": FG},
                                  "thickness": 12, "len": 0.85},
                        hovertemplate=(
                            "Station: %{x}<br>"
                            "Period: 10^%{y:.2f} s<br>"
                            f"{clabel}: " + "%{z:.2f}<extra></extra>"
                        ),
                    )
                    ps_fig = go.Figure(heat)

                    # mark selected period as horizontal line
                    ps_fig.add_hline(
                        y=np.log10(period),
                        line_dash="dash", line_color="#f38ba8",
                        line_width=1.5,
                        annotation_text=f"T={period} s",
                        annotation_font={"size": 9, "color": "#f38ba8"},
                    )

                    # 3-D marker: flag stations where |skew| > thresh at target period
                    if mark3d and "skew" in df.columns:
                        t_near = df.copy()
                        t_near["dt"] = np.abs(np.log10(t_near["period"])
                                               - np.log10(period))
                        nearest = (t_near.sort_values("dt")
                                         .groupby("station")
                                         .first()
                                         .reset_index())
                        flagged = nearest[nearest["skew"].abs() > skew_thr]["station"].tolist()
                        for fsta in flagged:
                            if fsta in sta_order:
                                xi = sta_order.index(fsta)
                                ps_fig.add_vrect(
                                    x0=xi - 0.4, x1=xi + 0.4,
                                    fillcolor="rgba(243,139,168,0.10)",
                                    line_width=0,
                                )

                    # y tick labels in period (not log)
                    tick_exp = np.arange(
                        np.floor(log_periods.min()),
                        np.ceil(log_periods.max()) + 1,
                    )
                    ps_fig.update_layout(
                        plot_bgcolor=AX_BG, paper_bgcolor=BG,
                        font={"color": FG, "size": 10},
                        margin={"l": 65, "r": 80, "t": 25, "b": 80},
                        height=320,
                        hovermode="closest",
                    )
                    ps_fig.update_xaxes(
                        tickangle=-45, gridcolor=GRID, zerolinecolor=GRID,
                        tickfont={"size": 9},
                    )
                    ps_fig.update_yaxes(
                        title="Period (s)",
                        tickvals=tick_exp.tolist(),
                        ticktext=[f"10^{int(e)}" for e in tick_exp],
                        gridcolor=GRID, zerolinecolor=GRID,
                        tickfont={"size": 9},
                    )

                    children.append(
                        html.Div([
                            html.Div("Phase Tensor Pseudo-section",
                                     style={"fontSize": "10px",
                                            "color": "var(--subtext0)",
                                            "textTransform": "uppercase",
                                            "letterSpacing": "0.05em",
                                            "marginBottom": "4px"}),
                            dcc.Graph(
                                figure=ps_fig,
                                config=_plotly_config("pycsamt_pt_psection"),
                                className="tool-graph",
                            ),
                        ])
                    )

            # ── MAP (matplotlib → PNG) ────────────────────────────────────────
            if view in ("map", "all"):
                try:
                    mpl_fig, ax = plt.subplots(figsize=(7, 6))
                    _style_ax(ax, mpl_fig)
                    plot_phase_tensor_map(
                        sub,
                        period=period,
                        c_by=colorby,
                        cmap=cmap,
                        clim_pct=clim_pct,
                        normalise_by=norm_by,
                        scale=scale_v,
                        show_tipper=show_tip,
                        tipper_convention=tip_conv,
                        mark_3d=mark3d,
                        skew_threshold=skew_thr,
                        edgecolor=EC,
                        ax=ax,
                    )
                    _style_ax(ax, mpl_fig)
                    children.append(
                        html.Div([
                            html.Div(
                                f"Phase Tensor Map  —  T = {period} s  —  colour: {colorby}",
                                style={"fontSize": "10px",
                                       "color": "var(--subtext0)",
                                       "textTransform": "uppercase",
                                       "letterSpacing": "0.05em",
                                       "marginBottom": "4px",
                                       "marginTop": "10px" if view == "all" else "0"},
                            ),
                            _mpl_to_img(mpl_fig),
                        ])
                    )
                except Exception as map_exc:
                    children.append(_warn(
                        f"Map could not be rendered: {map_exc}. "
                        "Check that stations have lat/lon coordinates."))

            if not children:
                return _warn("Nothing to display — select a view and run."), no_update

            new_saved = dict(saved_outputs or {})
            new_saved["phase-tensor"] = {"type": "custom"}
            return html.Div(children), new_saved

        except Exception as exc:
            return _err(exc), no_update


# ── Private helpers ───────────────────────────────────────────────────────────

def _style_fig(fig, ytitle: str, height: int = 220) -> None:
    fig.update_layout(
        plot_bgcolor="#1e1e2e", paper_bgcolor="#1e1e2e",
        font={"color": "#cdd6f4"},
        yaxis={"title": ytitle, "gridcolor": "#313244"},
        xaxis={"gridcolor": "#313244"},
        margin={"l": 50, "r": 20, "t": 20, "b": 50},
        height=height,
    )


def _tool_theme(theme: str | None) -> dict:
    light = str(theme or "").lower() == "light"
    if light:
        return {
            "bg": "#ffffff",
            "paper": "#ffffff",
            "grid": "#d8dce6",
            "fg": "#4c4f69",
            "muted": "#6c6f85",
            "border": "#ccd0da",
            "header": "#e6e9ef",
            "blue": "#1e66f5",
            "green": "#40a02b",
            "yellow": "#df8e1d",
            "red": "#d20f39",
            "row": "#ffffff",
        }
    return {
        "bg": "#1e1e2e",
        "paper": "#1e1e2e",
        "grid": "#313244",
        "fg": "#cdd6f4",
        "muted": "#a6adc8",
        "border": "#45475a",
        "header": "#313244",
        "blue": "#89b4fa",
        "green": "#a6e3a1",
        "yellow": "#f9e2af",
        "red": "#f38ba8",
        "row": "#1e1e2e",
    }


def _tool_table_styles(theme: str | None) -> tuple[dict, dict, list]:
    p = _tool_theme(theme)
    cell = {
        "fontSize": "11px",
        "fontFamily": "Inter, sans-serif",
        "padding": "6px",
        "backgroundColor": p["row"],
        "color": p["fg"],
        "border": f"1px solid {p['border']}",
    }
    header = {
        "backgroundColor": p["header"],
        "color": p["fg"],
        "fontWeight": "700",
        "border": f"1px solid {p['border']}",
    }
    conditional = [
        {
            "if": {"row_index": "odd"},
            "backgroundColor": "#f6f7fb" if str(theme).lower() == "light" else "#242438",
        },
        {
            "if": {"state": "active"},
            "backgroundColor": "#e8f0fe" if str(theme).lower() == "light" else "#313244",
            "border": f"1px solid {p['blue']}",
        },
    ]
    return cell, header, conditional


def _safe_stem(label: str) -> str:
    import re

    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label or "figure")).strip("._")
    return stem[:96] or "figure"


def _station_count(sites) -> int:
    try:
        from pycsamt.emtools._core import _iter_items

        return len(list(_iter_items(sites)))
    except Exception:
        return 0


def _save_mpl_figure(fig, path, *, dpi: int, transparent: bool, overwrite: bool):
    if path.exists() and not overwrite:
        return {
            "Figure": path.stem,
            "Status": "SKIPPED",
            "File": str(path),
            "Message": "File exists; enable overwrite to replace it.",
        }
    fig.savefig(
        str(path),
        dpi=int(dpi),
        bbox_inches="tight",
        transparent=bool(transparent),
        facecolor="none" if transparent else fig.get_facecolor(),
    )
    return {
        "Figure": path.stem,
        "Status": "SAVED",
        "File": str(path),
        "Message": "",
    }


def _plot_station_map_mpl(store_data: dict | None, active_lines_store: dict | None):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    records = list((store_data or {}).get("station_records", []) or [])
    df = pd.DataFrame(records)
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.set_title("Station map")
    if df.empty:
        ax.text(0.5, 0.5, "No station records", ha="center", va="center")
        return fig

    active = list((active_lines_store or {}).get("active", []) or [])
    if active and "Line" in df.columns:
        df = df[df["Line"].isin(active)].copy()
    if df.empty:
        ax.text(0.5, 0.5, "No active station records", ha="center", va="center")
        return fig

    x = pd.to_numeric(df.get("Longitude"), errors="coerce")
    y = pd.to_numeric(df.get("Latitude"), errors="coerce")
    xlabel, ylabel = "Longitude", "Latitude"
    if not (np.isfinite(x).any() and np.isfinite(y).any()):
        x = pd.to_numeric(df.get("Easting"), errors="coerce")
        y = pd.to_numeric(df.get("Northing"), errors="coerce")
        xlabel, ylabel = "Easting (m)", "Northing (m)"
    if not (np.isfinite(x).any() and np.isfinite(y).any()):
        x = np.arange(len(df), dtype=float)
        y = np.zeros(len(df), dtype=float)
        xlabel, ylabel = "Station index", "Profile offset"

    df = df.copy()
    df["_x"] = x
    df["_y"] = y
    line_col = "Line" if "Line" in df.columns else None
    groups = df.groupby(line_col, dropna=False) if line_col else [("Survey", df)]
    for line, sub in groups:
        sx = pd.to_numeric(sub["_x"], errors="coerce").to_numpy(dtype=float)
        sy = pd.to_numeric(sub["_y"], errors="coerce").to_numpy(dtype=float)
        good = np.isfinite(sx) & np.isfinite(sy)
        if not good.any():
            continue
        ax.plot(sx[good], sy[good], "-", alpha=0.45, linewidth=1.0)
        ax.scatter(sx[good], sy[good], s=34, label=str(line), zorder=3)
        ids = sub.get("ID", sub.get("station", None))
        if ids is not None and good.sum() <= 80:
            for xi, yi, sid in zip(sx[good], sy[good], np.asarray(ids)[good]):
                ax.text(xi, yi, str(sid), fontsize=6, ha="left", va="bottom")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.35)
    if xlabel in ("Longitude", "Easting (m)"):
        ax.set_aspect("equal", adjustable="datalim")
    if line_col:
        ax.legend(loc="best", fontsize=7)
    return fig


def _batch_export_tool_figures(
    sites,
    *,
    dest,
    fmt: str,
    dpi: int,
    items: list[str],
    flags: set[str],
    store_data: dict,
    active_lines_store: dict,
) -> dict:
    import csv

    import matplotlib.pyplot as plt

    fmt = (fmt or "png").lower().lstrip(".")
    if fmt not in {"png", "pdf", "svg", "eps", "tiff", "tif"}:
        raise ValueError(f"Unsupported export format: {fmt!r}.")
    overwrite = "overwrite" in flags
    transparent = "transparent" in flags
    manifest = "manifest" in flags
    selected = set(items or [])
    if "pseudo" in selected:
        selected.update({"rho_pseudo", "phi_pseudo"})

    rows: list[dict] = []

    def _attempt(label: str, builder):
        path = dest / f"{_safe_stem(label)}.{fmt}"
        fig = None
        try:
            fig = builder()
            rows.append(
                _save_mpl_figure(
                    fig,
                    path,
                    dpi=dpi,
                    transparent=transparent,
                    overwrite=overwrite,
                )
            )
        except Exception as exc:
            rows.append({
                "Figure": _safe_stem(label),
                "Status": "FAILED",
                "File": str(path),
                "Message": str(exc),
            })
        finally:
            if fig is not None:
                plt.close(fig)

    if "curves" in selected:
        def _curves():
            from pycsamt.emtools.inspect import plot_rhoa_phi

            ax_r, ax_p = plot_rhoa_phi(
                sites,
                components=("xy", "yx"),
                axis="period",
                errorbar=False,
                figsize=(8.4, 6.4),
                verbose=0,
            )
            fig = ax_r.figure
            fig.suptitle("Apparent resistivity and phase")
            return fig
        _attempt("mt_curves_rho_phase", _curves)

    if "rho_pseudo" in selected:
        def _rho_ps():
            import matplotlib.pyplot as plt

            from pycsamt.emtools.inspect import pseudosection

            fig, ax = plt.subplots(figsize=(9.2, 5.0))
            pseudosection(
                sites,
                quantity="rho_xy",
                ax=ax,
                dark=False,
                verbose=0,
            )
            ax.set_title("Apparent resistivity pseudosection (Zxy)")
            return fig
        _attempt("rho_xy_pseudosection", _rho_ps)

    if "phi_pseudo" in selected:
        def _phi_ps():
            import matplotlib.pyplot as plt

            from pycsamt.emtools.inspect import pseudosection

            fig, ax = plt.subplots(figsize=(9.2, 5.0))
            pseudosection(
                sites,
                quantity="phi_xy",
                ax=ax,
                dark=False,
                verbose=0,
            )
            ax.set_title("Phase pseudosection (Zxy)")
            return fig
        _attempt("phase_xy_pseudosection", _phi_ps)

    if "map" in selected:
        _attempt(
            "station_map",
            lambda: _plot_station_map_mpl(store_data, active_lines_store),
        )

    if "strike" in selected:
        def _skew():
            import matplotlib.pyplot as plt

            from pycsamt.emtools.skew import (
                plot_skew_traffic_psection,
            )

            fig, ax = plt.subplots(figsize=(9.0, 5.0))
            plot_skew_traffic_psection(sites, ax=ax, verbose=0)
            ax.set_title("Skew traffic pseudosection")
            return fig

        def _dim():
            import matplotlib.pyplot as plt

            from pycsamt.emtools.tensor import (
                plot_dimensionality_psection,
            )

            fig, ax = plt.subplots(figsize=(9.0, 5.0))
            plot_dimensionality_psection(sites, ax=ax, verbose=0)
            ax.set_title("Dimensionality pseudosection")
            return fig

        _attempt("skew_traffic_pseudosection", _skew)
        _attempt("dimensionality_pseudosection", _dim)

    manifest_path = dest / "pycsamt_batch_export_manifest.csv"
    if manifest:
        with manifest_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["Figure", "Status", "File", "Message"],
            )
            writer.writeheader()
            writer.writerows(rows)

    return {
        "dest": str(dest),
        "fmt": fmt,
        "dpi": int(dpi),
        "n_stations": _station_count(sites),
        "saved": sum(1 for r in rows if r["Status"] == "SAVED"),
        "skipped": sum(1 for r in rows if r["Status"] == "SKIPPED"),
        "failed": sum(1 for r in rows if r["Status"] == "FAILED"),
        "rows": rows,
        "manifest": str(manifest_path) if manifest else "",
    }


def _batch_result_view(result: dict) -> html.Div:
    rows = result.get("rows", [])
    cell, header, conditional = _tool_table_styles("dark")
    conditional = list(conditional or []) + [
        {
            "if": {"filter_query": "{Status} = 'SAVED'"},
            "color": "#40a02b",
            "fontWeight": "700",
        },
        {
            "if": {"filter_query": "{Status} = 'FAILED'"},
            "color": "#d20f39",
            "fontWeight": "700",
        },
        {
            "if": {"filter_query": "{Status} = 'SKIPPED'"},
            "color": "#df8e1d",
            "fontWeight": "700",
        },
    ]
    return html.Div(
        [
            html.Div(
                [
                    _metric_chip("Saved", str(result.get("saved", 0))),
                    _metric_chip("Skipped", str(result.get("skipped", 0))),
                    _metric_chip("Failed", str(result.get("failed", 0))),
                    _metric_chip("Stations", str(result.get("n_stations", 0))),
                ],
                className="tool-metric-row",
            ),
            html.Div(
                [
                    html.Span("Batch export complete", className="tool-section-title"),
                    html.Span(
                        f"{result.get('fmt', '').upper()} · "
                        f"{result.get('dpi')} dpi · {result.get('dest')}",
                        className="tool-section-subtitle",
                    ),
                    html.Span(
                        f"Manifest: {result['manifest']}",
                        className="tool-section-subtitle",
                    ) if result.get("manifest") else None,
                ],
                className="tool-conv-summary",
            ),
            dash_table.DataTable(
                columns=[
                    {"name": c, "id": c}
                    for c in ["Figure", "Status", "File", "Message"]
                ],
                data=rows,
                page_size=8,
                sort_action="native",
                style_table={
                    "overflowX": "auto",
                    "maxHeight": "360px",
                    "overflowY": "auto",
                },
                style_cell=cell,
                style_header=header,
                style_data_conditional=conditional,
            ) if rows else _warn("No figures were generated."),
        ]
    )


def _render_stored_tool_result(tool_id: str, run_store: dict | None, theme: str | None):
    payload = (run_store or {}).get(tool_id)
    if not payload:
        return None
    if tool_id == "strike":
        return _render_strike_result(payload, theme)
    return None


def _render_strike_result(payload: dict, theme: str | None):
    import numpy as np
    import pandas as pd

    records = payload.get("records", [])
    df = pd.DataFrame(records)
    fig_rose = _strike_rose_figure(df, payload.get("method_label", "Strike"), theme)
    fig_box = _strike_box_figure(df, theme)
    spread = float(payload.get("median_iqr", float("nan")))
    cell, header, conditional = _tool_table_styles(theme)
    return html.Div(
        [
            html.Div(
                [
                    _metric_chip("Scope", str(payload.get("scope_label", "Selected"))),
                    _metric_chip("Stations", str(payload.get("n_stations", 0))),
                    _metric_chip(
                        "Median strike",
                        f"{float(payload.get('median_strike', 0.0)):.1f}°",
                    ),
                    _metric_chip(
                        "Median IQR",
                        "n/a" if not np.isfinite(spread) else f"{spread:.1f}°",
                    ),
                ],
                className="tool-metric-row",
            ),
            dcc.Graph(
                figure=fig_rose,
                config=_plotly_config("pycsamt_strike_rose"),
                className="tool-graph",
            ),
            dcc.Graph(
                figure=fig_box,
                config=_plotly_config("pycsamt_strike_distribution"),
                className="tool-graph",
            ),
            dash_table.DataTable(
                columns=[{"name": c, "id": c} for c in payload.get("columns", [])],
                data=payload.get("table", []),
                page_size=int(payload.get("page_size", 10) or 10),
                sort_action="native",
                style_table={"overflowX": "auto"},
                style_cell=cell,
                style_header=header,
                style_data_conditional=conditional,
            ),
        ]
    )


def _strike_band_from_choice(choice: str | None):
    choice = (choice or "all").lower()
    if choice == "high":
        return (0.0, 0.01)       # f >= 100 Hz
    if choice == "low":
        return (1.0, 1.0e12)     # f <= 1 Hz
    return None


def _metric_chip(label: str, value: str) -> html.Div:
    return html.Div(
        [
            html.Span(label, className="tool-metric-label"),
            html.Span(value, className="tool-metric-value"),
        ],
        className="tool-metric-chip",
    )


def _strike_scope_label(
    store_data,
    active_lines_store,
    selected_lines,
    selected_stations,
) -> str:
    if selected_stations:
        return f"{len(selected_stations)} selected station(s)"
    if selected_lines:
        return ", ".join(selected_lines[:3]) + (
            f" +{len(selected_lines) - 3}" if len(selected_lines) > 3 else ""
        )
    active = _active_line_names(store_data, active_lines_store)
    if active:
        return f"All active lines ({len(active)})"
    return "All loaded stations"


def _strike_rose_figure(df, title: str, theme: str | None = "dark"):
    import numpy as np
    import plotly.graph_objects as go

    p = _tool_theme(theme)
    fig = go.Figure()
    bins = np.arange(0.0, 190.0, 10.0)
    centers = (bins[:-1] + bins[1:]) / 2.0
    colors = [
        "#89b4fa", "#a6e3a1", "#f9e2af", "#f38ba8", "#cba6f7",
        "#94e2d5", "#fab387", "#eba0ac",
    ]
    grouped = list(df.groupby("line", dropna=False))
    for idx, (line, sub) in enumerate(grouped):
        vals = np.asarray(sub["ang_axial"], dtype=float)
        vals = vals[np.isfinite(vals)]
        hist, _ = np.histogram(vals, bins=bins)
        fig.add_trace(
            go.Barpolar(
                r=hist.tolist(),
                theta=centers.tolist(),
                width=[9.5] * len(centers),
                name=str(line),
                marker_color=colors[idx % len(colors)],
                marker_line_color="rgba(255,255,255,0.20)",
                marker_line_width=0.5,
                opacity=0.72,
            )
        )
    fig.update_layout(
        title={"text": title, "font": {"size": 13}},
        polar={
            "bgcolor": p["bg"],
            "angularaxis": {
                "direction": "clockwise",
                "rotation": 90,
                "gridcolor": p["grid"],
                "tickfont": {"size": 10},
            },
            "radialaxis": {
                "showticklabels": False,
                "gridcolor": p["grid"],
                "ticks": "",
            },
        },
        plot_bgcolor=p["bg"],
        paper_bgcolor=p["paper"],
        font={"color": p["fg"]},
        margin={"l": 20, "r": 20, "t": 42, "b": 20},
        legend={"orientation": "h", "y": -0.08},
        height=300,
    )
    return fig


def _strike_box_figure(df, theme: str | None = "dark"):
    import plotly.graph_objects as go

    p = _tool_theme(theme)
    fig = go.Figure()
    colors = [
        "#89b4fa", "#a6e3a1", "#f9e2af", "#f38ba8", "#cba6f7",
        "#94e2d5", "#fab387", "#eba0ac",
    ]
    for idx, (line, sub) in enumerate(df.groupby("line", dropna=False)):
        fig.add_trace(
            go.Box(
                y=sub["ang_axial"].astype(float).tolist(),
                x=[str(line)] * len(sub),
                name=str(line),
                boxpoints="all",
                jitter=0.35,
                pointpos=-1.6,
                marker={"color": colors[idx % len(colors)], "size": 5},
                line={"color": colors[idx % len(colors)]},
                showlegend=False,
            )
        )
    fig.update_layout(
        title={"text": "Station strike distribution by line", "font": {"size": 13}},
        plot_bgcolor=p["bg"],
        paper_bgcolor=p["paper"],
        font={"color": p["fg"]},
        margin={"l": 50, "r": 20, "t": 40, "b": 50},
        height=260,
        xaxis={"title": "Line", "gridcolor": p["grid"], "zerolinecolor": p["border"]},
        yaxis={
            "title": "Axial strike (0-180°)",
            "gridcolor": p["grid"],
            "zerolinecolor": p["border"],
        },
    )
    return fig


def _validator_status_bar(n_pass: int, n_warn: int, n_fail: int) -> html.Div:
    total = max(1, n_pass + n_warn + n_fail)
    parts = [
        ("PASS", n_pass, "tool-valid-pass"),
        ("WARN", n_warn, "tool-valid-warn"),
        ("FAIL", n_fail, "tool-valid-fail"),
    ]
    return html.Div(
        [
            html.Div(
                title=f"{label}: {count}",
                className=f"tool-valid-bar-seg {cls}",
                style={"width": f"{100.0 * count / total:.2f}%"},
            )
            for label, count, cls in parts
            if count > 0
        ],
        className="tool-valid-bar",
    )


def _validate_sites_rows(sites, *, store_data, checks: set, min_freq: int) -> list[dict]:
    import numpy as np

    from pycsamt.emtools._core import (
        _get_t_block,
        _iter_items,
        _name,
    )

    records = _station_line_map(store_data)
    rows = []
    names = []
    items = list(_iter_items(sites))
    for idx, ed in enumerate(items):
        raw = getattr(ed, "edi", ed)
        station = _name(ed, idx)
        names.append(station)
        issues: list[str] = []
        warnings: list[str] = []
        row = {
            "Station": station,
            "Line": records.get(station, "Selected"),
            "Lat": "MISSING",
            "Lon": "MISSING",
            "Z ok": "—",
            "Errors": "—",
            "Tipper": "—",
            "N freq": 0,
            "T min (s)": "—",
            "T max (s)": "—",
            "NaN rows": "—",
        }

        lat, lon = _extract_lat_lon(raw, ed)
        if "coords" in checks:
            if lat is None or lon is None:
                warnings.append("Missing latitude/longitude.")
            else:
                try:
                    lat_f = float(lat)
                    lon_f = float(lon)
                    row["Lat"] = f"{lat_f:.5f}"
                    row["Lon"] = f"{lon_f:.5f}"
                    if not (-90.0 <= lat_f <= 90.0):
                        issues.append(f"Latitude out of range: {lat_f}.")
                    if not (-180.0 <= lon_f <= 180.0):
                        issues.append(f"Longitude out of range: {lon_f}.")
                except Exception:
                    warnings.append("Latitude/longitude cannot be parsed as numbers.")
        elif lat is not None and lon is not None:
            try:
                row["Lat"] = f"{float(lat):.5f}"
                row["Lon"] = f"{float(lon):.5f}"
            except Exception:
                pass

        Z, z, fr, ze = _z_block_with_errors(ed)
        if "z" in checks:
            if Z is None or z is None:
                issues.append("Missing impedance tensor Z.")
                row["Z ok"] = "NO"
            else:
                z_arr = np.asarray(z)
                if z_arr.ndim != 3 or z_arr.shape[1:3] != (2, 2):
                    issues.append(f"Z tensor has unexpected shape {z_arr.shape}.")
                    row["Z ok"] = "NO"
                else:
                    xy = np.isfinite(z_arr[:, 0, 1])
                    yx = np.isfinite(z_arr[:, 1, 0])
                    n_xy = int(np.count_nonzero(xy))
                    n_yx = int(np.count_nonzero(yx))
                    if n_xy and n_yx:
                        row["Z ok"] = "YES"
                    elif n_xy or n_yx:
                        row["Z ok"] = "PARTIAL"
                        warnings.append(
                            f"Only one off-diagonal component has finite data "
                            f"(xy={n_xy}, yx={n_yx})."
                        )
                    else:
                        row["Z ok"] = "NO"
                        issues.append("No finite off-diagonal Zxy/Zyx samples.")

        if "errors" in checks:
            if ze is None:
                warnings.append("Missing impedance error estimates.")
                row["Errors"] = "NO"
            else:
                ze_arr = np.asarray(ze, dtype=float)
                finite_err = int(np.count_nonzero(np.isfinite(ze_arr)))
                if finite_err <= 0:
                    warnings.append("Impedance error array exists but has no finite values.")
                    row["Errors"] = "NO"
                else:
                    row["Errors"] = "YES"

        if fr is not None:
            try:
                fa = np.asarray(fr, dtype=float).ravel()
                finite_freq = fa[np.isfinite(fa) & (fa > 0.0)]
            except Exception:
                finite_freq = np.asarray([], dtype=float)
        else:
            finite_freq = np.asarray([], dtype=float)

        if "freq" in checks:
            row["N freq"] = int(finite_freq.size)
            if finite_freq.size == 0:
                issues.append("Missing positive frequency values.")
            else:
                periods = 1.0 / finite_freq
                row["T min (s)"] = f"{float(np.nanmin(periods)):.4g}"
                row["T max (s)"] = f"{float(np.nanmax(periods)):.4g}"
                if finite_freq.size < int(min_freq):
                    warnings.append(
                        f"Only {finite_freq.size} positive frequencies "
                        f"(minimum requested: {int(min_freq)})."
                    )
                if np.unique(finite_freq).size != finite_freq.size:
                    warnings.append("Frequency vector contains duplicate values.")
                if finite_freq.size > 1:
                    d = np.diff(finite_freq)
                    if not (np.all(d > 0) or np.all(d < 0)):
                        warnings.append("Frequency vector is not monotonic.")
        else:
            row["N freq"] = int(finite_freq.size)
            if finite_freq.size:
                periods = 1.0 / finite_freq
                row["T min (s)"] = f"{float(np.nanmin(periods)):.4g}"
                row["T max (s)"] = f"{float(np.nanmax(periods)):.4g}"

        if "nan" in checks and z is not None:
            z_arr = np.asarray(z)
            if z_arr.ndim == 3 and z_arr.shape[0]:
                flat = z_arr.reshape(z_arr.shape[0], -1)
                nan_rows = int(np.count_nonzero(np.all(~np.isfinite(flat), axis=1)))
                row["NaN rows"] = nan_rows
                if nan_rows:
                    warnings.append(f"{nan_rows} frequency row(s) have all-NaN Z.")

        if "tipper" in checks:
            T, tip, tf = _get_t_block(ed)
            if T is None or tip is None:
                row["Tipper"] = "NO"
                warnings.append("Missing tipper data.")
            else:
                tip_arr = np.asarray(tip)
                row["Tipper"] = "YES" if np.isfinite(tip_arr).any() else "NO"
                if row["Tipper"] == "NO":
                    warnings.append("Tipper data has no finite values.")

        status = "FAIL" if issues else ("WARN" if warnings else "PASS")
        all_issues = issues + warnings
        row["Status"] = status
        row["Issues"] = " ".join(all_issues) if all_issues else "OK"
        rows.append(row)

    if "duplicates" in checks and names:
        from collections import Counter

        dupes = {name for name, count in Counter(names).items() if count > 1}
        if dupes:
            for row in rows:
                if row["Station"] in dupes:
                    prefix = "Duplicate station identifier."
                    row["Issues"] = (
                        prefix if row["Issues"] == "OK"
                        else f"{prefix} {row['Issues']}"
                    )
                    if row["Status"] == "PASS":
                        row["Status"] = "WARN"
    return rows


def _converter_options_from_values(
    *,
    output_dir,
    stn_path,
    utm_zone,
    epsg,
    convert_coords,
    e_labels,
    h_labels,
    suffix,
    spectra_flags,
    freq_order,
    freq_tol,
    core_flags,
    station_name,
) -> dict:
    flags = set(core_flags or [])
    spectra_flags = set(spectra_flags or [])
    options = {
        "output_dir": str(output_dir or "").strip(),
        "stn_path": str(stn_path or "").strip(),
        "utm_zone": str(utm_zone or "").strip(),
        "epsg": str(epsg or "").strip(),
        "convert_stn_coords": bool(convert_coords),
        "freq_order": freq_order or "descending",
        "compute_z": "compute_z" in flags,
        "compute_rho_phi": "compute_rho_phi" in flags,
        "e_labels": str(e_labels or "EX,EY").strip(),
        "h_labels": str(h_labels or "HX,HY").strip(),
        "station_suffix": str(suffix or "").strip(),
        "estimate_error": "estimate_error" in spectra_flags,
        "use_remote": "use_remote" in spectra_flags,
        "skip_errors": "skip_errors" in spectra_flags,
        "station_name": str(station_name or "").strip(),
    }
    if freq_tol not in (None, ""):
        options["freq_tol"] = float(freq_tol)
    return options


def _converted_sites_to_store(sites, *, old_store: dict | None, data_dir: str) -> dict:
    from pycsamt.app.desktop.controllers.data_controller import (
        DataController,
    )

    ctrl = DataController()
    ctrl._sites = sites
    if old_store and old_store.get("station_records"):
        records = old_store["station_records"]
        ctrl._station_to_line = {
            r["ID"]: r.get("Line", "converted") for r in records if r.get("ID")
        }
    ctrl._df = ctrl._build_dataframe()
    df = ctrl.dataframe
    if "Line" in df.columns:
        line_counts = df.groupby("Line").size().to_dict()
        n_lines = int(df["Line"].nunique())
    else:
        line_counts = {"converted": len(df)}
        n_lines = 1
    return {
        "station_records": df.to_dict("records"),
        "data_dir": data_dir or "[converted]",
        "n_stations": len(df),
        "n_lines": n_lines,
        "line_counts": line_counts,
        "converted": True,
    }


def _converter_result_view(*, stats: dict, failures: list, output_dir: str, loaded: bool):
    rows = list((stats or {}).get("rows", []) or [])
    n_total = int((stats or {}).get("n_total", len(rows)) or 0)
    n_fail = int((stats or {}).get("n_failures", len(failures or [])) or 0)
    n_with_coords = sum(
        1 for r in rows
        if _is_finite_number(r.get("lat")) and _is_finite_number(r.get("lon"))
    )
    n_with_elev = sum(1 for r in rows if _is_finite_number(r.get("elev")))
    sum(1 for r in rows if bool(r.get("has_Z")))
    table_rows = []
    for r in rows:
        table_rows.append({
            "Station": r.get("station", ""),
            "N freq": r.get("n_freqs", 0),
            "f min": _fmt_float(r.get("f_min")),
            "f max": _fmt_float(r.get("f_max")),
            "Lat": _fmt_float(r.get("lat"), nd=5),
            "Lon": _fmt_float(r.get("lon"), nd=5),
            "Elev": _fmt_float(r.get("elev"), nd=2),
            "Z": "yes" if r.get("has_Z") else "no",
            "Tipper": "yes" if r.get("has_tipper") else "no",
        })
    fail_rows = [
        {
            "Source": getattr(f, "source", None) or f.get("source", ""),
            "Error": getattr(f, "error", None) or f.get("error", ""),
        }
        for f in (failures or [])
        if isinstance(f, dict) or getattr(f, "source", None) is not None
    ]
    return html.Div(
        [
            html.Div(
                [
                    _metric_chip("Converted", str(n_total)),
                    _metric_chip("Failures", str(n_fail)),
                    _metric_chip("Coordinates", f"{n_with_coords}/{n_total}"),
                    _metric_chip("Elevation", f"{n_with_elev}/{n_total}"),
                ],
                className="tool-metric-row",
            ),
            html.Div(
                [
                    html.Span("Conversion complete", className="tool-section-title"),
                    html.Span(
                        (
                            f"Written to {output_dir}" if output_dir
                            else "No output folder provided; result kept in memory."
                        )
                        + (" Loaded into session." if loaded else ""),
                        className="tool-section-subtitle",
                    ),
                ],
                className="tool-conv-summary",
            ),
            dash_table.DataTable(
                columns=[{"name": c, "id": c} for c in (
                    "Station", "N freq", "f min", "f max", "Lat", "Lon",
                    "Elev", "Z", "Tipper"
                )],
                data=table_rows,
                page_size=10,
                sort_action="native",
                style_table={"overflowX": "auto", "maxHeight": "360px", "overflowY": "auto"},
                style_cell={
                    "fontSize": "11px",
                    "fontFamily": "Inter, sans-serif",
                    "padding": "6px",
                    "backgroundColor": "var(--mantle)",
                    "color": "var(--text)",
                    "border": "1px solid var(--surface0)",
                },
                style_header={
                    "backgroundColor": "var(--surface0)",
                    "fontWeight": "700",
                },
            ) if table_rows else _warn("No EDI objects were produced."),
            html.Div(
                [
                    html.Span("Failures", className="tool-section-title"),
                    dash_table.DataTable(
                        columns=[{"name": c, "id": c} for c in ("Source", "Error")],
                        data=fail_rows,
                        page_size=5,
                        style_table={"overflowX": "auto"},
                        style_cell={
                            "fontSize": "11px",
                            "fontFamily": "Inter, sans-serif",
                            "padding": "6px",
                            "backgroundColor": "var(--mantle)",
                            "color": "var(--text)",
                            "border": "1px solid var(--surface0)",
                        },
                    ),
                ],
                className="mt-2",
            ) if fail_rows else None,
        ]
    )


def _fmt_float(value, nd: int = 4) -> str:
    try:
        import numpy as np
        v = float(value)
        if not np.isfinite(v):
            return "—"
        return f"{v:.{nd}g}" if nd <= 4 else f"{v:.{nd}f}"
    except Exception:
        return "—"


def _is_finite_number(value) -> bool:
    try:
        import numpy as np
        return bool(np.isfinite(float(value)))
    except Exception:
        return False


def _z_block_with_errors(ed):
    from pycsamt.emtools._core import _get_z_block

    out = _get_z_block(ed, with_errors=True)
    if len(out) == 4:
        return out
    Z, z, fr = out
    ze = getattr(Z, "z_err", None) if Z is not None else None
    return Z, z, fr, ze


def _extract_lat_lon(raw, wrapped=None):
    for obj in (raw, wrapped):
        if obj is None:
            continue
        lat = _first_existing_attr(obj, ("lat", "latitude"))
        lon = _first_existing_attr(obj, ("lon", "longitude"))
        if lat is not None or lon is not None:
            return lat, lon
        try:
            coords = getattr(obj, "coords", None)
            if coords is not None and len(coords) >= 2:
                return coords[0], coords[1]
        except Exception:
            pass
    for head_name in ("Header", "HEAD", "head"):
        head = getattr(raw, head_name, None)
        if head is None:
            continue
        lat = _first_existing_attr(head, ("lat", "latitude"))
        lon = _first_existing_attr(head, ("lon", "longitude"))
        if lat is not None or lon is not None:
            return lat, lon
    try:
        h = getattr(raw, "Head", None) or getattr(raw, "head", None)
        if isinstance(h, dict):
            return h.get("lat") or h.get("latitude"), h.get("lon") or h.get("longitude")
    except Exception:
        pass
    return None, None


def _first_existing_attr(obj, names):
    for name in names:
        try:
            value = getattr(obj, name)
        except Exception:
            value = None
        if value is not None:
            return value
    return None


def _err(exc: Exception) -> html.Span:
    return html.Span(f"✗ {exc}", className="text-danger small")


def _warn(msg: str) -> html.Span:
    return html.Span(f"⚠ {msg}", className="text-warning small")


def _dd2dms(dd: float, axis: str = "lat") -> str:
    """Decimal degrees → DMS string, e.g. '40°42′46.80″N'."""
    hemi = ("N" if dd >= 0 else "S") if axis == "lat" else ("E" if dd >= 0 else "W")
    dd = abs(dd)
    d = int(dd)
    m = int((dd - d) * 60)
    s = (dd - d - m / 60) * 3600
    return f"{d}°{m:02d}′{s:05.2f}″{hemi}"


def _coming_soon(name: str) -> html.Div:
    return html.Div([
        html.I(className="bi bi-clock-history me-2 text-muted"),
        html.Span(f"{name} backend not yet available in web mode.",
                  className="text-muted small"),
    ], className="mt-2")


def _elev_summary_table(df, title: str = "Elevation") -> html.Div:
    """Compact summary table for elevation fetch/load results."""
    import pandas as pd

    cols = [c for c in ["station", "lat", "lon", "elev", "elev_corrected"]
            if c in df.columns]
    display_df = df[cols].copy()
    for c in ["lat", "lon", "elev", "elev_corrected"]:
        if c in display_df.columns:
            display_df[c] = display_df[c].apply(
                lambda v: f"{v:.4f}" if pd.notna(v) else "—"
            )

    n_ok  = int(df["elev"].notna().sum()) if "elev" in df.columns else 0
    total = len(df)
    header = html.Div(
        [
            html.I(className="bi bi-check-circle-fill me-1",
                   style={"color": "var(--green)"}),
            f" {title}: {n_ok}/{total} station(s)  ·  "
            f"Switch to Correct & Export to process.",
        ],
        className="small mb-2",
    )
    tbl = dash_table.DataTable(
        data=display_df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in display_df.columns],
        page_size=8,
        style_table={"overflowX": "auto", "fontSize": "11px"},
        style_cell={"padding": "3px 8px", "border": "1px solid var(--srf1)",
                    "background": "var(--base)", "color": "var(--text)"},
        style_header={"background": "var(--srf0)", "fontWeight": "600",
                      "border": "1px solid var(--srf1)"},
    )
    return html.Div([header, tbl])
