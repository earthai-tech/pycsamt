# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Workbench layout for the pyCSAMT Map View platform.

Four zones: a top bar, a left view-rail + data/lines panel, a dominant
map canvas, a contextual right inspector, and a collapsible bottom data
dock.  A web-style Load modal handles multi-line EDI ingestion.
"""

from __future__ import annotations

import uuid

from dash import dcc, html
import dash_bootstrap_components as dbc

from ._ids import IDs

try:
    from pycsamt.map import BASEMAP_STYLES
except Exception:  # pragma: no cover - fallback if map extra missing
    BASEMAP_STYLES = (
        ("open-street-map", "Open Street Map"),
        ("esri-satellite", "ESRI Satellite"),
        ("esri-topo", "ESRI Topographic"),
        ("carto-positron", "Light (Carto)"),
        ("carto-darkmatter", "Dark (Carto)"),
    )

_COMPONENTS = ["xy", "yx", "xx", "yy", "det", "avg"]
_CMAPS = [
    "plasma", "viridis", "jet", "turbo", "magma",
    "cividis", "RdYlBu_r", "RdBu_r", "coolwarm",
]
_OVERLAYS = [
    ("index", "Station index"),
    ("elevation", "Elevation"),
    ("rho", "App. resistivity"),
    ("phase", "Phase"),
    ("skin_depth", "Skin depth"),
]
_MODES_3D = [
    ("fence", "Fence"),
    ("block", "Block"),
    ("depth", "Depth slices"),
    ("surface", "Iso-surface"),
]
_VIEWS = [
    (IDs.RAIL_MAP, "map", "map-view", "Map"),
    (IDs.RAIL_3D, "map3d", "3d", "3-D"),
]


# ── top bar ────────────────────────────────────────────


def _topbar() -> html.Div:
    return html.Div(
        id="mv-topbar",
        children=[
            html.Button(
                html.I(className="bi bi-layout-sidebar"),
                id=IDs.BTN_SIDEBAR,
                className="mv-icon-btn",
                title="Toggle data panel",
                n_clicks=0,
            ),
            html.Div(
                [
                    html.Img(
                        src="/mv-icons/pycsamt-v2-symbol.svg",
                        className="mv-brand-logo",
                    ),
                    html.Span(
                        [
                            html.Span("py", className="mv-brand-py"),
                            html.Span("CSAMT", className="mv-brand-csamt"),
                            html.Span(" Map View", className="mv-brand-sub"),
                        ],
                        className="mv-brand-text",
                    ),
                ],
                className="mv-brand",
            ),
            html.Button(
                [html.I(className="bi bi-folder2-open me-1"), "Load lines"],
                id=IDs.BTN_LOAD,
                className="mv-tbtn primary",
                n_clicks=0,
            ),
            html.Button(
                [html.I(className="bi bi-layer-forward me-1"), "Add line"],
                id=IDs.BTN_ADD_LINE,
                className="mv-tbtn",
                n_clicks=0,
            ),
            html.Div(
                [
                    html.I(className="bi bi-check-circle-fill me-1"),
                    html.Span("", id=IDs.DATA_BADGE_TEXT),
                ],
                id=IDs.DATA_BADGE,
                className="mv-data-badge",
            ),
            html.Div(className="mv-topbar-spacer"),
            html.Button(
                [html.I(className="bi bi-download me-1"), "Export"],
                id=IDs.BTN_EXPORT,
                className="mv-tbtn",
                n_clicks=0,
            ),
            html.Button(
                html.I(className="bi bi-moon-stars", id="mv-theme-icon"),
                id=IDs.BTN_THEME,
                className="mv-icon-btn",
                title="Toggle light / dark",
                n_clicks=0,
            ),
            html.Button(
                html.I(className="bi bi-question-circle"),
                id=IDs.BTN_HELP,
                className="mv-icon-btn",
                title="Help & About",
                n_clicks=0,
            ),
            html.Button(
                html.I(className="bi bi-layout-sidebar-reverse"),
                id=IDs.BTN_INSPECTOR,
                className="mv-icon-btn",
                title="Toggle inspector panel",
                n_clicks=0,
            ),
        ],
    )


# ── left rail + data panel ─────────────────────────────


def _view_rail() -> html.Div:
    return html.Div(
        id="mv-rail",
        className="mv-rail",
        children=[
            html.Button(
                [
                    html.Img(src=f"/mv-icons/{icon}.svg",
                             className="mv-rail-icon"),
                    html.Span(label),
                ],
                id=btn_id,
                className=(
                    "mv-rail-btn"
                    + (" mv-rail-active" if view == "map" else "")
                ),
                n_clicks=0,
                title=label,
            )
            for btn_id, view, icon, label in _VIEWS
        ],
    )


def _data_panel() -> html.Div:
    return html.Div(
        id="mv-datapanel",
        className="mv-datapanel",
        children=[
            html.Div("Lines", className="mv-panel-lbl"),
            html.Div(id=IDs.LINE_PILLS, className="mv-line-pills"),
            html.Div(
                id=IDs.STA_LOAD_BAR,
                className="mv-sta-load-bar mv-sta-load-bar--hidden",
            ),
            html.Div("Stations", className="mv-panel-lbl mt-2"),
            html.Div(
                html.Div(
                    "No lines loaded.",
                    className="mv-empty",
                ),
                id=IDs.STATION_LIST,
                className="mv-station-list",
            ),
        ],
    )


# ── right inspector ────────────────────────────────────


def _ctl_row(label, comp):
    return html.Div(
        [html.Label(label, className="mv-field-lbl"), comp],
        className="mv-field-row",
    )


def _num(id_, value, **kw):
    return dbc.Input(id=id_, type="number", value=value, size="sm",
                     debounce=True, **kw)


def _coord_system() -> html.Div:
    """Coordinate system: interpret station coords + reproject to lon/lat."""
    return html.Div(
        [
            _ctl_row(
                "Coordinate system",
                dbc.Select(
                    id=IDs.CTL_CRS_MODE,
                    options=[
                        {"label": "Geographic (lat/lon · EPSG:4326)",
                         "value": "geo"},
                        {"label": "UTM zone", "value": "utm"},
                        {"label": "Custom EPSG", "value": "custom"},
                    ],
                    value="geo",
                    size="sm",
                ),
            ),
            html.Div(
                dbc.InputGroup(
                    [
                        dbc.InputGroupText("Zone"),
                        dbc.Input(id=IDs.CTL_UTM_ZONE, type="number",
                                  value=50, min=1, max=60, step=1, size="sm"),
                        dbc.Select(
                            id=IDs.CTL_UTM_HEM,
                            options=[{"label": "N", "value": "N"},
                                     {"label": "S", "value": "S"}],
                            value="N", size="sm",
                        ),
                    ],
                    size="sm",
                ),
                id=IDs.GRP_UTM,
                style={"display": "none"},
                className="mb-1",
            ),
            html.Div(
                dbc.InputGroup(
                    [
                        dbc.InputGroupText("EPSG"),
                        dbc.Input(id=IDs.CTL_EPSG, type="text", value="32650",
                                  placeholder="e.g. 32650", size="sm"),
                    ],
                    size="sm",
                ),
                id=IDs.GRP_EPSG,
                style={"display": "none"},
                className="mb-1",
            ),
            html.Div(id=IDs.CRS_INFO, className="mv-crs-info"),
        ],
    )


def _map_controls() -> html.Div:
    """All Map-view controls (points over basemap + contour overlay)."""
    return html.Div(
        id=IDs.GRP_MAP,
        children=[
            _ctl_row(
                "Colour by",
                dbc.Select(
                    id=IDs.CTL_OVERLAY,
                    options=[{"label": lbl, "value": val}
                             for val, lbl in _OVERLAYS],
                    value="index",
                    size="sm",
                ),
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Frequency", className="mv-field-lbl"),
                            html.Span("—", id=IDs.CTL_FREQ_LABEL,
                                      className="mv-freq-val"),
                        ],
                        className="mv-field-head",
                    ),
                    dcc.Slider(
                        id=IDs.CTL_FREQUENCY,
                        min=0, max=1, step=1, value=0,
                        marks=None,
                        tooltip={"placement": "bottom"},
                    ),
                ],
                className="mv-field-row",
            ),
            _ctl_row(
                "Basemap",
                dbc.Select(
                    id=IDs.CTL_BASEMAP,
                    options=[{"label": lbl, "value": val}
                             for val, lbl in BASEMAP_STYLES],
                    value="esri-satellite",
                    size="sm",
                ),
            ),
            _coord_system(),
            html.Div(
                [
                    html.Div([html.Label("Marker size",
                                         className="mv-field-lbl"),
                              _num(IDs.CTL_MARKER_SIZE, 10, min=4, max=24,
                                   step=1)],
                             className="mv-col"),
                    html.Div([html.Label("Opacity %",
                                         className="mv-field-lbl"),
                              _num(IDs.CTL_OPACITY_MAP, 92, min=20, max=100,
                                   step=2)],
                             className="mv-col"),
                ],
                className="mv-two-col",
            ),
            html.Div(
                [
                    dbc.Switch(id=IDs.CTL_PROFILES, label="Profile lines",
                               value=True, className="mv-switch"),
                ],
                className="mv-switch-row",
            ),
            html.Hr(className="mv-hr"),
            dbc.Switch(id=IDs.CTL_CONTOUR_ENABLE,
                       label="Contour overlay (Surfer-style)",
                       value=False, className="mv-switch"),
            _ctl_row(
                "Contour style",
                dbc.Select(
                    id=IDs.CTL_CONTOUR_MODE,
                    options=[
                        {"label": "Filled + lines", "value": "filled+lines"},
                        {"label": "Filled", "value": "filled"},
                        {"label": "Lines", "value": "lines"},
                    ],
                    value="filled+lines",
                    size="sm",
                ),
            ),
            html.Div(
                [
                    html.Div([html.Label("Levels", className="mv-field-lbl"),
                              _num(IDs.CTL_CONTOUR_LEVELS, 12, min=2, max=40,
                                   step=1)],
                             className="mv-col"),
                    html.Div([html.Label("Smoothing σ",
                                         className="mv-field-lbl"),
                              _num(IDs.CTL_CONTOUR_SMOOTH, 1.0, min=0, max=5,
                                   step=0.5)],
                             className="mv-col"),
                ],
                className="mv-two-col",
            ),
            _ctl_row(
                "Interpolation",
                dbc.Select(
                    id=IDs.CTL_CONTOUR_INTERP,
                    options=[
                        {"label": "Cubic", "value": "cubic"},
                        {"label": "Linear", "value": "linear"},
                        {"label": "Nearest", "value": "nearest"},
                    ],
                    value="cubic",
                    size="sm",
                ),
            ),
            _ctl_row(
                "Grid resolution",
                dbc.Select(
                    id=IDs.CTL_CONTOUR_RES,
                    options=[
                        {"label": "Coarse (100²)", "value": "100"},
                        {"label": "Normal (150²)", "value": "150"},
                        {"label": "Fine (200²)", "value": "200"},
                        {"label": "Ultra (300²)", "value": "300"},
                    ],
                    value="150",
                    size="sm",
                ),
            ),
        ],
    )


def _three_d_group() -> html.Div:
    """3-D-only controls (mode + geometry + rendering)."""
    return html.Div(
        id=IDs.GRP_3D,
        style={"display": "none"},
        children=[
            html.Hr(className="mv-hr"),
            html.Div("3-D options", className="mv-panel-lbl"),
            _ctl_row(
                "Quantity",
                dbc.Select(
                    id=IDs.CTL_QUANTITY,
                    options=[
                        {"label": "Resistivity", "value": "rho"},
                        {"label": "Phase", "value": "phase"},
                    ],
                    value="rho",
                    size="sm",
                ),
            ),
            _ctl_row(
                "Mode",
                dbc.Select(
                    id=IDs.CTL_MODE3D,
                    options=[{"label": lbl, "value": val}
                             for val, lbl in _MODES_3D],
                    value="fence",
                    size="sm",
                ),
            ),
            html.Div(
                [
                    html.Label("Opacity", className="mv-field-lbl"),
                    dcc.Slider(
                        id=IDs.CTL_OPACITY,
                        min=0.1, max=1.0, step=0.05, value=0.85,
                        marks={0.1: "10%", 1.0: "100%"},
                        tooltip={"placement": "bottom"},
                    ),
                ],
                className="mv-field-row",
            ),
            html.Div(
                [
                    html.Div([_num(IDs.CTL_AZIMUTH, 0, min=-180, max=180,
                                   step=5)],
                             className="mv-col"),
                    html.Div([_num(IDs.CTL_SPACING, 1.0, min=0.1, step=0.1)],
                             className="mv-col"),
                ],
                className="mv-two-col",
            ),
            html.Div(
                [
                    html.Span("Azimuth °", className="mv-mini-lbl"),
                    html.Span("Line spacing", className="mv-mini-lbl"),
                ],
                className="mv-two-col mv-mini-row",
            ),
            _ctl_row(
                "Depth window (m)",
                html.Div(
                    [
                        _num(IDs.CTL_DEPTH_LO, None, min=0,
                             placeholder="min"),
                        _num(IDs.CTL_DEPTH_HI, None, min=0,
                             placeholder="max"),
                    ],
                    className="mv-two-col",
                ),
            ),
            html.Div(
                [
                    html.Div([html.Label("Depth slices",
                                         className="mv-field-lbl"),
                              _num(IDs.CTL_NSLICES, 8, min=2, max=40, step=1)],
                             className="mv-col"),
                    html.Div([html.Label("Iso-surfaces",
                                         className="mv-field-lbl"),
                              _num(IDs.CTL_SURFACES, 12, min=2, max=50,
                                   step=1)],
                             className="mv-col"),
                ],
                className="mv-two-col",
            ),
            _ctl_row(
                "Color scale",
                dbc.Select(
                    id=IDs.CTL_SCALE,
                    options=[{"label": "Log", "value": "log"},
                             {"label": "Linear", "value": "linear"}],
                    value="log",
                    size="sm",
                ),
            ),
            _ctl_row(
                "Colorbar range (Ω·m)",
                html.Div(
                    [
                        _num(IDs.CTL_VMIN, None, min=0, placeholder="min"),
                        _num(IDs.CTL_VMAX, None, min=0, placeholder="max"),
                    ],
                    className="mv-two-col",
                ),
            ),
            _ctl_row(
                "ρ visibility (Ω·m)",
                html.Div(
                    [
                        _num(IDs.CTL_RHO_LO, None, min=0, placeholder="min"),
                        _num(IDs.CTL_RHO_HI, None, min=0, placeholder="max"),
                    ],
                    className="mv-two-col",
                ),
            ),
            dbc.Switch(id=IDs.CTL_CONTOURS, label="Contours", value=False,
                       className="mv-switch"),
            dbc.Switch(id=IDs.CTL_TOPO, label="Drape topography", value=False,
                       className="mv-switch"),
            dbc.Switch(id=IDs.CTL_TERRAIN, label="Show terrain line",
                       value=False, className="mv-switch"),
            dbc.Switch(id=IDs.CTL_SHOW_STA, label="Show station markers",
                       value=False, className="mv-switch"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Symbol", className="mv-field-lbl"),
                            dbc.Select(
                                id=IDs.CTL_STA_SYMBOL,
                                options=[
                                    {"label": s.title(), "value": s}
                                    for s in ("diamond", "circle", "square",
                                              "cross", "x", "diamond-open",
                                              "circle-open")
                                ],
                                value="diamond",
                                size="sm",
                            ),
                        ],
                        className="mv-col",
                    ),
                    html.Div(
                        [
                            html.Label("Size", className="mv-field-lbl"),
                            _num(IDs.CTL_STA_SIZE, 4, min=2, max=20, step=1),
                        ],
                        className="mv-col",
                    ),
                    html.Div(
                        [
                            html.Label("Color", className="mv-field-lbl"),
                            dbc.Input(id=IDs.CTL_STA_COLOR, type="color",
                                      value="#1f2937",
                                      style={"height": "31px",
                                             "padding": "2px"}),
                        ],
                        className="mv-col",
                    ),
                ],
                className="mv-two-col",
            ),
            _elevation_source(),
        ],
    )


def _elevation_source() -> html.Div:
    """Elevation source: station EDIs, an uploaded file, or online fetch."""
    return html.Div(
        [
            html.Div("Elevation source", className="mv-panel-lbl mt-2"),
            dbc.Select(
                id=IDs.TOPO_SRC,
                options=[
                    {"label": "Station EDIs", "value": "stations"},
                    {"label": "Upload file", "value": "upload"},
                    {"label": "Fetch online", "value": "fetch"},
                ],
                value="stations",
                size="sm",
            ),
            html.Div(
                dcc.Upload(
                    id=IDs.TOPO_UPLOAD,
                    children=html.Div(
                        [html.I(className="bi bi-cloud-upload me-1"),
                         "Drop / pick .csv · .h5 · .npz"],
                    ),
                    className="mv-upload-drop mv-topo-drop",
                    multiple=False,
                    accept=".csv,.h5,.hdf5,.npz",
                ),
                id=IDs.TOPO_UPLOAD_WRAP,
                style={"display": "none"},
            ),
            html.Div(
                dbc.Select(
                    id=IDs.TOPO_API,
                    options=[
                        {"label": "Open-Meteo", "value": "open_meteo"},
                        {"label": "OpenTopoData",
                         "value": "open_topo_data"},
                    ],
                    value="open_meteo",
                    size="sm",
                ),
                id=IDs.TOPO_FETCH_WRAP,
                style={"display": "none"},
            ),
            dbc.Button(
                [html.I(className="bi bi-geo-alt me-1"), "Apply elevations"],
                id=IDs.BTN_TOPO_APPLY,
                color="secondary",
                size="sm",
                outline=True,
                className="w-100 mt-2",
                n_clicks=0,
            ),
            html.Div(id=IDs.TOPO_STATUS, className="mv-topo-status"),
        ],
    )


def _inspector() -> html.Div:
    return html.Div(
        id=IDs.INSPECTOR,
        className="mv-inspector",
        children=[
            html.Div("Inspector", className="mv-panel-lbl"),
            _map_controls(),
            _ctl_row(
                "Component",
                dbc.Select(
                    id=IDs.CTL_COMPONENT,
                    options=[
                        {"label": c.upper(), "value": c}
                        for c in _COMPONENTS
                    ],
                    value="xy",
                    size="sm",
                ),
            ),
            _ctl_row(
                "Colormap",
                dbc.Select(
                    id=IDs.CTL_CMAP,
                    options=[
                        {"label": c, "value": c} for c in _CMAPS
                    ],
                    value="plasma",
                    size="sm",
                ),
            ),
            html.Div(
                [
                    dbc.Switch(
                        id=IDs.CTL_LOG,
                        label="Log scale",
                        value=False,
                        className="mv-switch",
                    ),
                    dbc.Switch(
                        id=IDs.CTL_LABELS,
                        label="Show labels",
                        value=True,
                        className="mv-switch",
                    ),
                ],
                className="mv-switch-row",
            ),
            _three_d_group(),
            html.Hr(className="mv-hr"),
            html.Div("Station", className="mv-panel-lbl"),
            html.Div(
                html.Div("Click a station on the map.",
                         className="mv-empty"),
                id=IDs.STATION_INSPECT,
                className="mv-station-inspect",
            ),
        ],
    )


# ── canvas + bottom dock ───────────────────────────────


def _tb_btn(icon, label, btn_id, title=""):
    inner = [html.I(className=f"bi {icon}")]
    if label:
        inner.append(html.Span(label, className="mv-tb-label"))
    return html.Button(inner, id=btn_id, n_clicks=0,
                       className="mv-tb-btn", title=title or label)


def _canvas_toolbar() -> html.Div:
    return html.Div(
        [
            html.Span("Map view", id=IDs.CANVAS_TITLE,
                      className="mv-canvas-title"),
            html.Span(id=IDs.TB_INFO, className="mv-tb-info"),
            html.Div(className="mv-tb-sep"),
            _tb_btn("bi-arrows-fullscreen", "Fit", IDs.TB_FIT,
                    "Zoom to fit all stations"),
            html.Div(className="mv-tb-sep"),
            _tb_btn("bi-card-text", "Labels", IDs.TB_LABELS,
                    "Toggle station labels"),
            _tb_btn("bi-bezier", "Profiles", IDs.TB_PROFILES,
                    "Toggle survey-line polylines"),
            _tb_btn("bi-layers-half", "Contour", IDs.TB_CONTOUR,
                    "Toggle contour overlay"),
            html.Div(className="mv-tb-sep"),
            html.Span("Basemap", className="mv-tb-grp"),
            _tb_btn("bi-moon-stars-fill", "", IDs.TB_BM_DARK, "Dark (Carto)"),
            _tb_btn("bi-sun-fill", "", IDs.TB_BM_LIGHT, "Light (Carto)"),
            _tb_btn("bi-globe-americas", "", IDs.TB_BM_SAT,
                    "Satellite (ESRI)"),
            _tb_btn("bi-signpost-split", "", IDs.TB_BM_STREET,
                    "Street (ESRI)"),
            _tb_btn("bi-map-fill", "", IDs.TB_BM_TOPO, "Topographic (ESRI)"),
            html.Div(className="mv-tb-sep"),
            html.Span("Markers", className="mv-tb-grp"),
            _tb_btn("bi-dash-lg", "", IDs.TB_MARK_DEC, "Smaller markers"),
            html.Span("10", id=IDs.TB_MARK_VAL, className="mv-tb-markval"),
            _tb_btn("bi-plus-lg", "", IDs.TB_MARK_INC, "Larger markers"),
        ],
        className="mv-toolbar",
    )


def _canvas() -> html.Div:
    return html.Div(
        id="mv-canvas-wrap",
        className="mv-canvas-wrap",
        children=[
            _canvas_toolbar(),
            dcc.Loading(
                dcc.Graph(
                    id=IDs.CANVAS_GRAPH,
                    className="mv-graph",
                    config={
                        "displaylogo": False,
                        "responsive": True,
                        "toImageButtonOptions": {"format": "png", "scale": 2},
                    },
                    style={"height": "100%"},
                ),
                type="circle",
                color="#1e66f5",
            ),
        ],
    )


def _dock() -> html.Div:
    return html.Div(
        id=IDs.DOCK_TABLE,
        className="mv-dock",
        children=[
            html.Button(
                [
                    html.I(className="bi bi-table me-1"),
                    html.Span("Stations"),
                    html.I(className="bi bi-chevron-up ms-2",
                           id="mv-dock-chevron"),
                ],
                id=IDs.DOCK_TOGGLE,
                className="mv-dock-toggle",
                n_clicks=0,
            ),
            html.Div(
                dash_table_placeholder(),
                id=IDs.DOCK_BODY,
                className="mv-dock-body",
            ),
        ],
    )


def dash_table_placeholder():
    from dash import dash_table

    return dash_table.DataTable(
        id="mv-station-table",
        columns=[
            {"name": c, "id": c}
            for c in ["ID", "Line", "Latitude", "Longitude",
                      "Elevation", "Index"]
        ],
        data=[],
        page_size=8,
        style_as_list_view=True,
        style_table={"overflowX": "auto"},
        style_cell={
            "fontSize": "11px",
            "padding": "4px 8px",
            "fontFamily": "monospace",
        },
    )


# ── load modal (web-style) ─────────────────────────────


def _load_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle(
                    [html.I(className="bi bi-cloud-upload-fill me-2"),
                     "Load Survey Lines"],
                    id=IDs.MODAL_TITLE,
                ),
                close_button=True,
            ),
            dbc.ModalBody(
                [
                    # mode toggle
                    html.Div(
                        [
                            html.Button(
                                [html.I(className="bi bi-arrow-repeat me-1"),
                                 "Replace"],
                                id=IDs.MODE_BTN_REPLACE,
                                className="mv-mode-btn active",
                                n_clicks=0,
                            ),
                            html.Button(
                                [html.I(className="bi bi-layer-forward me-1"),
                                 "Add lines"],
                                id=IDs.MODE_BTN_APPEND,
                                className="mv-mode-btn",
                                n_clicks=0,
                            ),
                        ],
                        className="mv-mode-toggle",
                    ),
                    html.Div(
                        "Existing survey data will be replaced.",
                        id=IDs.MODE_HINT,
                        className="mv-mode-hint",
                    ),
                    # browse + drop
                    dbc.Button(
                        [html.I(className="bi bi-folder2-open me-2"),
                         html.Span("Browse folder",
                                   style={"fontWeight": "600"}),
                         html.Span(" — pick a survey directory",
                                   style={"fontSize": "11px",
                                          "opacity": ".75"})],
                        id=IDs.BTN_BROWSE,
                        color="primary",
                        className="w-100 mb-2",
                        n_clicks=0,
                    ),
                    html.Div(
                        [
                            html.Div("Or drag & drop files / folders",
                                     id=IDs.DROP_TITLE,
                                     className="mv-panel-lbl"),
                            dcc.Upload(
                                id=IDs.UPLOAD,
                                children=html.Div(
                                    [html.I(className="bi bi-file-earmark-plus me-2"),
                                     "Drop EDI files or click to pick"],
                                ),
                                className="mv-upload-drop",
                                multiple=True,
                                accept=".edi,.EDI,.avg,.j",
                            ),
                        ],
                        id=IDs.DROP_WRAP,
                    ),
                    # spinner overlay (driven by JS)
                    html.Div(
                        [
                            dbc.Spinner(size="sm", color="primary"),
                            html.Span("Reading files…", id=IDs.LOADER_MSG,
                                      className="ms-2"),
                        ],
                        id=IDs.LOADER_OVERLAY,
                        className="mv-loader-overlay",
                        style={"display": "none"},
                    ),
                    html.Div(
                        [
                            html.Span("", id=IDs.FILE_COUNT,
                                      className="mv-file-count"),
                            html.Span("", id=IDs.BROWSE_STATUS,
                                      className="mv-browse-status"),
                        ],
                        className="d-flex align-items-center gap-2 mt-1",
                    ),
                    # file manager + preflight
                    html.Div(id=IDs.FILE_MANAGER, className="mv-file-manager"),
                    html.Div(id=IDs.PREFLIGHT, className="mv-preflight"),
                    # progress bar
                    html.Div(
                        [
                            html.Div(
                                html.Div(id=IDs.PROG_FILL,
                                         className="mv-prog-fill"),
                                className="mv-prog-track",
                            ),
                            html.Div(
                                [
                                    html.Span("Loading…", id=IDs.PROG_LABEL,
                                              className="mv-prog-label"),
                                    html.Span("", id=IDs.PROG_SUBLABEL,
                                              className="mv-prog-sublabel"),
                                ],
                                className="mv-prog-meta",
                            ),
                        ],
                        id=IDs.PROGRESS_WRAP,
                        className="mv-progress-wrap",
                        style={"display": "none"},
                    ),
                ]
            ),
            dbc.ModalFooter(
                [
                    html.Div(
                        [
                            html.Div(id=IDs.DETECTED_SUMMARY,
                                     className="mv-detected-summary"),
                            html.Div(id=IDs.LOAD_FEEDBACK,
                                     className="mv-load-feedback"),
                        ],
                        className="mv-footer-info",
                    ),
                    html.Div(className="mv-topbar-spacer"),
                    dbc.Button(
                        [html.I(className="bi bi-check-lg me-1"),
                         "Load into view"],
                        id=IDs.BTN_LOAD_CONFIRM,
                        color="primary",
                        disabled=True,
                        n_clicks=0,
                    ),
                ]
            ),
        ],
        id=IDs.MODAL_LOAD,
        size="lg",
        is_open=False,
        centered=True,
        className="mv-modal",
    )


def _help_cap(icon, title, desc) -> html.Div:
    return html.Div(
        [
            html.Div(html.I(className=f"bi {icon}"),
                     className="mv-help-cap-icon"),
            html.Div(
                [
                    html.Div(title, className="mv-help-cap-title"),
                    html.Div(desc, className="mv-help-cap-desc"),
                ]
            ),
        ],
        className="mv-help-cap",
    )


def _help_modal() -> dbc.Modal:
    try:
        from importlib.metadata import version as _v
        _ver = _v("pycsamt")
    except Exception:
        _ver = "2.0"

    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle(
                    [
                        html.Img(src="/mv-icons/pycsamt-v2-symbol.svg",
                                 className="mv-help-logo"),
                        html.Span(
                            [
                                html.Span("py", className="mv-brand-py"),
                                html.Span("CSAMT", className="mv-brand-csamt"),
                                html.Span(" Map View", className="mv-brand-sub"),
                            ],
                            className="mv-help-title",
                        ),
                    ],
                    className="mv-help-titlewrap",
                ),
                close_button=True,
            ),
            dbc.ModalBody(
                [
                    html.P(
                        "An interactive map viewer for CSAMT / AMT / MT "
                        "surveys — station maps, depth-slice contour maps "
                        "over satellite imagery, and 3-D models.",
                        className="mv-help-lead",
                    ),
                    html.Div("What you can do", className="mv-help-section"),
                    html.Div(
                        [
                            _help_cap("bi-folder2-open", "Load survey lines",
                                      "EDI folders, drag-drop, or add lines "
                                      "incrementally to one survey."),
                            _help_cap("bi-geo-alt", "Station & profile maps",
                                      "Plot stations or depth-slice contours "
                                      "over satellite / street / topo basemaps."),
                            _help_cap("bi-box", "3-D models",
                                      "Fence, block, depth-slice and "
                                      "iso-surface views with topography."),
                            _help_cap("bi-download", "Export",
                                      "Save any view as an interactive HTML "
                                      "file; everything is also scriptable "
                                      "via pycsamt.map.MapView."),
                        ],
                        className="mv-help-caps",
                    ),
                    html.Div("Tips", className="mv-help-section"),
                    html.Ul(
                        [
                            html.Li("Use the left rail to switch view: "
                                    "Station, Profile, Pseudosection, 3-D."),
                            html.Li("Toggle the line pills to show/hide "
                                    "survey lines; click a station to inspect."),
                            html.Li("In 3-D, enable Drape topography and set "
                                    "an elevation source (EDIs, upload, fetch)."),
                            html.Li("Collapse the inspector (⟷ in the top bar) "
                                    "to give the map more room."),
                        ],
                        className="mv-help-tips",
                    ),
                ]
            ),
            dbc.ModalFooter(
                [
                    html.Span(f"pyCSAMT v{_ver} · LGPL-3.0 · L. Kouadio",
                              className="mv-help-meta"),
                    html.Div(className="mv-topbar-spacer"),
                    html.A([html.I(className="bi bi-book me-1"), "Docs"],
                           href="https://pycsamt.readthedocs.io",
                           target="_blank", className="mv-help-link"),
                    html.A([html.I(className="bi bi-github me-1"), "GitHub"],
                           href="https://github.com/earthai-tech/pycsamt",
                           target="_blank", className="mv-help-link"),
                    dbc.Button("Close", id=IDs.BTN_HELP_CLOSE, n_clicks=0,
                               size="sm"),
                ]
            ),
        ],
        id=IDs.MODAL_HELP,
        is_open=False,
        size="lg",
        centered=True,
        scrollable=True,
    )


# ── stores ─────────────────────────────────────────────


def _stores() -> list:
    return [
        dcc.Store(id=IDs.SESSION_ID, data=str(uuid.uuid4())),
        dcc.Store(id=IDs.STORE_DATA, data={}),
        dcc.Store(id=IDs.STORE_THEME, storage_type="local", data="light"),
        dcc.Store(id=IDs.STORE_VIEW, data="map"),
        dcc.Store(id=IDs.STORE_CONTROLS, data={}),
        dcc.Store(id=IDs.STORE_LINES, data={}),
        dcc.Store(id=IDs.STORE_LINE_FILTER, data=None),
        dcc.Store(id=IDs.STORE_SELECTION, data={}),
        dcc.Store(id=IDs.LOAD_MODE_STORE, data="replace"),
        dcc.Store(id=IDs.FOLDER_STORE, data={}),
        dcc.Store(id=IDs.UPLOAD_SELECTION, data=[]),
        dcc.Store(id=IDs.SOURCE_SELECTION, data="none"),
        dcc.Store(id=IDs.TOPO_UPLOAD_STORE, data={}),
        dcc.Store(id=IDs.STORE_FIT, data=0),
        dcc.Download(id=IDs.EXPORT_DL),
    ]


# ── full layout ────────────────────────────────────────


def create_layout() -> html.Div:
    return html.Div(
        id=IDs.ROOT,
        className="mv-root mv-root-light",
        children=[
            *_stores(),
            _load_modal(),
            _help_modal(),
            _topbar(),
            html.Div(
                id="mv-body",
                className="mv-body",
                children=[
                    _view_rail(),
                    _data_panel(),
                    html.Div(
                        id="mv-center",
                        className="mv-center",
                        children=[_canvas(), _dock()],
                    ),
                    _inspector(),
                ],
            ),
        ],
    )
