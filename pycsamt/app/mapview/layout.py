# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Workbench layout for the pyCSAMT Map View platform.

Four zones: a top bar, a left view-rail + data/lines panel, a dominant
map canvas, a contextual right inspector, and a collapsible bottom data
dock.  A web-style Load modal handles multi-line EDI ingestion.
"""

from __future__ import annotations

import uuid

import dash_bootstrap_components as dbc
from dash import dcc, html

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
    "plasma", "viridis", "jet", "jet_r", "turbo", "magma",
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
                        src="/mv-icons/pycsamt_logo.svg",
                        className="mv-brand-logo",
                    ),
                    html.Span(
                        "Map View",
                        className="mv-brand-sub",
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
                [html.I(className="bi bi-floppy me-1"), "Session"],
                id=IDs.BTN_SESSION,
                className="mv-tbtn",
                title="Save / restore this workbench session",
                n_clicks=0,
            ),
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
                [html.I(className="bi bi-sliders2 me-1"), "Sites"],
                id=IDs.BTN_SETTINGS,
                className="mv-tbtn",
                title="Station / sites settings",
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


def _depth_presets() -> html.Div:
    """Quick depth-window presets (like the web map view)."""
    def _b(label, bid):
        return dbc.Button(label, id=bid, size="sm", color="secondary",
                          outline=True, className="mv-preset-btn", n_clicks=0)
    return html.Div(
        [
            _b("Full", IDs.BTN_DEPTH_FULL),
            _b("500 m", IDs.BTN_DEPTH_500),
            _b("1 km", IDs.BTN_DEPTH_1K),
            _b("2 km", IDs.BTN_DEPTH_2K),
        ],
        className="mv-preset-bar",
    )


def _rho_presets() -> html.Div:
    """Quick resistivity-band presets (like the web map view).

    In fence/depth-slice mode this masks individual cells to
    transparent; in block/iso-surface mode it sets the visible
    isomin/isomax band (e.g. "Conduct." isolates the conductive
    zone, "Resist." isolates the resistive zone) — everything
    outside the picked band fades out rather than being removed,
    so the block still reads as one solid shape.
    """
    def _b(label, bid):
        return dbc.Button(label, id=bid, size="sm", color="secondary",
                          outline=True, className="mv-preset-btn", n_clicks=0)
    return html.Div(
        [
            _b("All", IDs.BTN_RHO_ALL),
            _b("Conduct.", IDs.BTN_RHO_COND),
            _b("Mid", IDs.BTN_RHO_MID),
            _b("Resist.", IDs.BTN_RHO_RES),
        ],
        className="mv-preset-bar",
    )


def _acc_item(title, icon, item_id, children):
    return dbc.AccordionItem(
        children,
        title=html.Span([html.I(className=f"bi {icon} me-2"), title]),
        item_id=item_id,
    )


def _three_d_group() -> html.Div:
    """3-D-only controls, grouped into collapsible accordion sections."""
    mode_sec = _acc_item("Mode & quantity", "bi-box", "3d-mode", [
        _ctl_row("Quantity", dbc.Select(
            id=IDs.CTL_QUANTITY,
            options=[{"label": "Resistivity", "value": "rho"},
                     {"label": "Phase", "value": "phase"}],
            value="rho", size="sm")),
        _ctl_row("Mode", dbc.Select(
            id=IDs.CTL_MODE3D,
            options=[{"label": lbl, "value": val} for val, lbl in _MODES_3D],
            value="fence", size="sm")),
        _ctl_row("Aspect / auto-fit", dbc.Select(
            id=IDs.CTL_ASPECT,
            options=[{"label": "True proportions (recommended)", "value": "data"},
                     {"label": "Equal cube (clear)", "value": "cube"},
                     {"label": "Stretch to box", "value": "manual"}],
            value="data", size="sm")),
        html.Div([
            html.Div([html.Label("Profile unit", className="mv-field-lbl"),
                      dbc.Select(id=IDs.CTL_X_UNIT,
                                 options=[{"label": "metres", "value": "m"},
                                          {"label": "km", "value": "km"}],
                                 value="m", size="sm")],
                     className="mv-col"),
            html.Div([html.Label("Depth unit", className="mv-field-lbl"),
                      dbc.Select(id=IDs.CTL_DEPTH_UNIT,
                                 options=[{"label": "metres", "value": "m"},
                                          {"label": "km", "value": "km"}],
                                 value="m", size="sm")],
                     className="mv-col"),
        ], className="mv-two-col"),
    ])

    depth_sec = _acc_item("Depth filter", "bi-arrows-collapse", "3d-depth", [
        html.Div("Quick window", className="mv-field-lbl"),
        _depth_presets(),
        _ctl_row("Depth window (m)", html.Div([
            _num(IDs.CTL_DEPTH_LO, None, min=0, placeholder="min"),
            _num(IDs.CTL_DEPTH_HI, None, min=0, placeholder="max"),
        ], className="mv-two-col")),
        html.Div([
            html.Div([html.Label("Depth slices", className="mv-field-lbl"),
                      _num(IDs.CTL_NSLICES, 8, min=2, max=40, step=1)],
                     className="mv-col"),
            html.Div([html.Label("Iso-surfaces", className="mv-field-lbl"),
                      _num(IDs.CTL_SURFACES, 12, min=2, max=50, step=1)],
                     className="mv-col"),
        ], className="mv-two-col"),
    ])

    rho_sec = _acc_item("Resistivity range", "bi-funnel", "3d-rho", [
        html.Div("ρ visibility band", className="mv-field-lbl"),
        _rho_presets(),
        _ctl_row("ρ visibility (Ω·m)", html.Div([
            _num(IDs.CTL_RHO_LO, None, min=0, placeholder="min"),
            _num(IDs.CTL_RHO_HI, None, min=0, placeholder="max"),
        ], className="mv-two-col")),
        html.Div(
            "Block / iso-surface: fades everything outside this band "
            "(isolate the conductive or resistive zone). "
            "Fence / depth-slice: hides individual cells outside it.",
            className="mv-help-hint",
            style={"fontSize": "10.5px", "opacity": ".7", "marginTop": "-4px"},
        ),
    ])

    geom_sec = _acc_item("Geometry", "bi-compass", "3d-geom", [
        html.Div([
            html.Div([html.Label("Azimuth °", className="mv-field-lbl"),
                      _num(IDs.CTL_AZIMUTH, 0, min=-180, max=180, step=5)],
                     className="mv-col"),
            html.Div([html.Label("Line spacing", className="mv-field-lbl"),
                      _num(IDs.CTL_SPACING, 1.0, min=0.1, step=0.1)],
                     className="mv-col"),
        ], className="mv-two-col"),
    ])

    sta_sec = _acc_item("Stations (3-D)", "bi-geo", "3d-sta", [
        dbc.Switch(id=IDs.CTL_SHOW_STA, label="Show station markers",
                   value=False, className="mv-switch"),
        dbc.Switch(id=IDs.CTL_STA_LABELS, label="Show station labels",
                   value=False, className="mv-switch"),
        html.Div([
            html.Div([html.Label("Symbol", className="mv-field-lbl"),
                      dbc.Select(id=IDs.CTL_STA_SYMBOL,
                                 options=[{"label": s.title(), "value": s}
                                          for s in ("diamond", "circle",
                                                    "square", "cross", "x",
                                                    "diamond-open",
                                                    "circle-open")],
                                 value="diamond", size="sm")],
                     className="mv-col"),
            html.Div([html.Label("Size", className="mv-field-lbl"),
                      _num(IDs.CTL_STA_SIZE, 4, min=2, max=20, step=1)],
                     className="mv-col"),
            html.Div([html.Label("Color", className="mv-field-lbl"),
                      dbc.Input(id=IDs.CTL_STA_COLOR, type="color",
                                value="#1f2937",
                                style={"height": "31px", "padding": "2px"})],
                     className="mv-col"),
        ], className="mv-two-col"),
    ])

    topo_sec = _acc_item("Topography", "bi-geo-alt", "3d-topo", [
        dbc.Switch(id=IDs.CTL_TOPO, label="Drape topography", value=True,
                   className="mv-switch"),
        dbc.Switch(id=IDs.CTL_TERRAIN, label="Show terrain line",
                   value=True, className="mv-switch"),
        _elevation_source(),
    ])

    appearance_sec = _acc_item("Appearance", "bi-palette", "3d-appearance", [
        html.Div([
            html.Label("Opacity", className="mv-field-lbl"),
            dcc.Slider(id=IDs.CTL_OPACITY, min=0.1, max=1.0, step=0.05,
                       value=0.85, marks={0.1: "10%", 1.0: "100%"},
                       tooltip={"placement": "bottom"}),
        ], className="mv-field-row"),
        _ctl_row("Color scale", dbc.Select(
            id=IDs.CTL_SCALE,
            options=[{"label": "Log", "value": "log"},
                     {"label": "Linear", "value": "linear"}],
            value="log", size="sm")),
        _ctl_row("Colorbar range (Ω·m)", html.Div([
            _num(IDs.CTL_VMIN, None, min=0, placeholder="min"),
            _num(IDs.CTL_VMAX, None, min=0, placeholder="max"),
        ], className="mv-two-col")),
        dbc.Switch(id=IDs.CTL_CONTOURS, label="Contours", value=False,
                   className="mv-switch"),
        html.Hr(className="mv-hr"),
        dbc.Switch(id=IDs.CTL_SMOOTH, label="Smooth sections (interpolate)",
                   value=True, className="mv-switch"),
        _ctl_row("Section resolution", dbc.Select(
            id=IDs.CTL_SECTION_RES,
            options=[{"label": "Coarse (60)", "value": "60"},
                     {"label": "Normal (100)", "value": "100"},
                     {"label": "Fine (160)", "value": "160"},
                     {"label": "Ultra (240)", "value": "240"}],
            value="100", size="sm")),
    ])

    return html.Div(
        id=IDs.GRP_3D,
        style={"display": "none"},
        children=[
            html.Hr(className="mv-hr"),
            html.Div("3-D options", className="mv-panel-lbl"),
            dbc.Accordion(
                # Most-frequently-tuned params first, purely cosmetic
                # appearance settings (color, opacity, …) last.
                [mode_sec, depth_sec, rho_sec, geom_sec, topo_sec,
                 sta_sec, appearance_sec],
                active_item=["3d-mode", "3d-depth"],
                always_open=True,
                flush=True,
                className="mv-3d-accordion",
            ),
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
            html.Hr(className="mv-hr"),
            html.Div(
                "Export topography", className="mv-panel-lbl mt-2"
            ),
            html.Div(
                "Save this survey's station elevations (e.g. real, "
                "field-surveyed EDI elevation) to re-apply later on "
                "an inversion result — which carries no elevation of "
                "its own — via “Upload file” above, even in a "
                "session that never reloads these EDIs.",
                className="mv-help-hint",
                style={"fontSize": "10.5px", "opacity": ".7"},
            ),
            html.Div(
                [
                    dbc.Select(
                        id=IDs.TOPO_EXPORT_FMT,
                        options=[
                            {"label": "CSV", "value": "csv"},
                            {"label": "HDF5", "value": "h5"},
                        ],
                        value="csv",
                        size="sm",
                    ),
                    dbc.Button(
                        [html.I(className="bi bi-download me-1"), "Export"],
                        id=IDs.BTN_TOPO_EXPORT,
                        color="secondary",
                        size="sm",
                        outline=True,
                        n_clicks=0,
                    ),
                ],
                className="mv-two-col mt-1",
            ),
            dcc.Download(id=IDs.TOPO_EXPORT_DL),
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
    """2-D map-view toolbar — Fit / labels / profiles / contour / basemap /
    markers. Hidden in 3-D mode (see ``_canvas_toolbar_3d`` and
    ``callbacks.toolbar._register_view_visibility``)."""
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
        id=IDs.TOOLBAR_2D,
        className="mv-toolbar",
    )


def _canvas_toolbar_3d() -> html.Div:
    """3-D map-view toolbar — reset view / render mode / depth window /
    topography. Hidden by default; shown only in 3-D mode (see
    ``callbacks.toolbar._register_view_visibility``). Mirrors, and writes
    into, the same stores as the Inspector's "3-D options" accordion, so
    the two never drift out of sync."""
    return html.Div(
        [
            html.Span("3-D view", className="mv-canvas-title"),
            html.Div(className="mv-tb-sep"),
            _tb_btn("bi-arrows-fullscreen", "Reset view", IDs.TB3D_RESET,
                    "Reset camera to fit the data"),
            html.Div(className="mv-tb-sep"),
            html.Span("Mode", className="mv-tb-grp"),
            _tb_btn("bi-bezier", "Fence", IDs.TB3D_MODE_FENCE,
                    "Vertical fence sections"),
            _tb_btn("bi-box", "Block", IDs.TB3D_MODE_BLOCK,
                    "Solid block volume"),
            _tb_btn("bi-layers-half", "Depth", IDs.TB3D_MODE_DEPTH,
                    "Horizontal depth slices"),
            _tb_btn("bi-circle-half", "Iso", IDs.TB3D_MODE_SURFACE,
                    "Iso-resistivity surfaces"),
            html.Div(className="mv-tb-sep"),
            html.Span("Depth", className="mv-tb-grp"),
            _tb_btn("bi-arrows-vertical", "Full", IDs.TB3D_DEPTH_FULL,
                    "Full depth range"),
            _tb_btn("bi-arrows-vertical", "500 m", IDs.TB3D_DEPTH_500,
                    "0 - 500 m"),
            _tb_btn("bi-arrows-vertical", "1 km", IDs.TB3D_DEPTH_1K,
                    "0 - 1000 m"),
            _tb_btn("bi-arrows-vertical", "2 km", IDs.TB3D_DEPTH_2K,
                    "0 - 2000 m"),
            html.Div(className="mv-tb-sep"),
            _tb_btn("bi-geo-alt", "Topo", IDs.TB3D_TOPO,
                    "Toggle topography drape"),
        ],
        id=IDs.TOOLBAR_3D,
        className="mv-toolbar",
        style={"display": "none"},
    )


def _welcome_card(icon, title, desc, tone) -> html.Div:
    return html.Div(
        [
            html.I(
                className=f"bi {icon} mv-welcome-card-icon",
                style={"color": f"var(--mv-{tone})"},
            ),
            html.Div(title, className="mv-welcome-card-title"),
            html.Div(desc, className="mv-welcome-card-desc"),
        ],
        className="mv-welcome-card",
    )


def _welcome() -> html.Div:
    particles = html.Div(
        [
            html.Span(className=f"mv-welcome-particle p{i}")
            for i in range(1, 10)
        ],
        className="mv-welcome-particles",
    )
    return html.Div(
        id=IDs.WELCOME,
        className="mv-welcome",
        children=[
            particles,
            html.Div(
                [
                    html.Img(
                        src="/mv-icons/pycsamt-v2-symbol.svg",
                        className="mv-welcome-logo",
                    ),
                    html.H1(
                        [
                            "Welcome to ",
                            html.Span("py", className="mv-brand-py"),
                            html.Span("CSAMT", className="mv-brand-csamt"),
                            html.Span(
                                " Map View",
                                className="mv-brand-sub",
                            ),
                        ],
                        className="mv-welcome-title",
                    ),
                    html.P(
                        "Where survey geometry becomes spatial insight",
                        className="mv-welcome-slogan",
                    ),
                    html.P(
                        "Basemap station maps · profile overlays · "
                        "3-D survey sections",
                        className="mv-welcome-subtitle",
                    ),
                    html.Div(
                        [
                            _welcome_card(
                                "bi-geo-alt-fill",
                                "Load Data",
                                "Plot lines on satellite, street, "
                                "topographic or light/dark maps.",
                                "py",
                            ),
                            _welcome_card(
                                "bi-layers-half",
                                "Map Layers",
                                "Switch colour fields, contours, "
                                "profiles and station labels.",
                                "ok",
                            ),
                            _welcome_card(
                                "bi-box",
                                "3-D Views",
                                "Inspect fence, block, depth-slice "
                                "and surface views.",
                                "csamt",
                            ),
                            _welcome_card(
                                "bi-download",
                                "Export",
                                "Export figures or reproduce them "
                                "with pycsamt.map.MapView.",
                                "accent",
                            ),
                        ],
                        className="mv-welcome-cards",
                    ),
                    html.Button(
                        [
                            html.I(
                                className=(
                                    "bi bi-cloud-upload-fill me-2"
                                )
                            ),
                            "Load Survey Lines - Start",
                        ],
                        id=IDs.WELCOME_CTA,
                        className="mv-welcome-cta",
                        n_clicks=0,
                    ),
                    html.P(
                        "Drop EDI, AVG or J files, or choose a folder "
                        "containing multiple survey lines.",
                        className="mv-welcome-hint",
                    ),
                ],
                className="mv-welcome-hero",
            )
        ],
    )


def _canvas() -> html.Div:
    return html.Div(
        id="mv-canvas-wrap",
        className="mv-canvas-wrap",
        children=[
            _canvas_toolbar(),
            _canvas_toolbar_3d(),
            dcc.Loading(
                html.Div(
                    [
                        dcc.Graph(
                            id=IDs.CANVAS_GRAPH,
                            className="mv-graph",
                            config={
                                "displaylogo": False,
                                "responsive": True,
                                "toImageButtonOptions": {
                                    "format": "png",
                                    "scale": 2,
                                },
                            },
                            style={"height": "100%"},
                        ),
                    ],
                    className="mv-graph-stage",
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
                [
                    html.Div(
                        [
                            html.Span(
                                "Station data",
                                className="mv-dock-body-title",
                            ),
                            html.Button(
                                [
                                    html.I(className="bi bi-chevron-down me-1"),
                                    "Hide",
                                ],
                                id=IDs.DOCK_CLOSE,
                                className="mv-dock-close",
                                n_clicks=0,
                            ),
                        ],
                        className="mv-dock-body-head",
                    ),
                    dash_table_placeholder(),
                ],
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
                dbc.Tabs(
                    [
                        dbc.Tab(
                            _load_modal_edi_body(),
                            label="EDI stations",
                            tab_id="load-tab-edi",
                        ),
                        dbc.Tab(
                            _load_modal_inversion_body(),
                            label="Inversion results",
                            tab_id="load-tab-inversion",
                        ),
                    ],
                    id=IDs.LOAD_TABS,
                    active_tab="load-tab-edi",
                ),
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


def _load_modal_edi_body() -> html.Div:
    return html.Div(
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
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Span(
                                        "Survey folders",
                                        className="mv-folder-filter-title",
                                    ),
                                    html.Span(
                                        "Select one or many; empty loads all.",
                                        className="mv-folder-filter-hint",
                                    ),
                                ],
                                className="mv-folder-filter-head",
                            ),
                            dcc.Checklist(
                                id=IDs.LOAD_LINE_FILTER,
                                options=[],
                                value=[],
                                className="mv-folder-filter-list",
                                labelClassName="mv-folder-filter-option",
                                inputClassName="mv-folder-filter-input",
                            ),
                        ],
                        id=IDs.LOAD_LINE_FILTER_WRAP,
                        className="mv-folder-filter",
                        style={"display": "none"},
                    ),
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
        ],
    )


def _load_modal_inversion_body() -> html.Div:
    """ModEM inversion-result import: native folder picker + confirm.

    Same browse-button pattern as the EDI tab (``edi_loader.js``): a
    transparent ``<input webkitdirectory>`` is injected over the
    button by ``modem_loader.js``, so the OS's native folder dialog
    opens instead of an in-app folder navigator.
    """
    return html.Div(
        [
            html.Div(
                "Import a ModEM 3-D inversion result folder. A single "
                "run inverts one volume for the whole survey — it is "
                "sliced into one panel per detected line.",
                className="mv-panel-lbl",
                style={"marginBottom": "8px"},
            ),
            dbc.Button(
                [html.I(className="bi bi-folder2-open me-2"),
                 html.Span("Browse folder",
                           style={"fontWeight": "600"}),
                 html.Span(" — pick a ModEM results directory",
                           style={"fontSize": "11px",
                                  "opacity": ".75"})],
                id=IDs.BTN_INV_BROWSE,
                color="primary",
                className="w-100 mb-2",
                n_clicks=0,
            ),
            html.Div(
                [
                    dbc.Spinner(size="sm", color="primary"),
                    html.Span("Reading files…", id=IDs.INV_LOADER_MSG,
                              className="ms-2"),
                ],
                id=IDs.INV_LOADER_OVERLAY,
                className="mv-loader-overlay",
                style={"display": "none"},
            ),
            html.Div(
                [
                    html.Span("", id=IDs.INV_FILE_COUNT,
                              className="mv-file-count"),
                    html.Span("", id=IDs.INV_BROWSE_STATUS,
                              className="mv-browse-status"),
                ],
                className="d-flex align-items-center gap-2 mt-1 mb-2",
            ),
            dbc.Switch(
                id=IDs.CK_INV_KNOWN_STA,
                label="Match coordinates from already-loaded EDI stations",
                value=True,
                className="mv-switch mb-2",
            ),
            dbc.Button(
                [html.I(className="bi bi-cpu me-1"), "Import inversion results"],
                id=IDs.BTN_INV_CONFIRM,
                color="primary",
                className="w-100",
                n_clicks=0,
            ),
            html.Div(id=IDs.INV_STATUS, className="mv-load-feedback mt-2"),
        ],
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


def _help_tip(icon, text) -> html.Li:
    return html.Li(
        [
            html.I(className=f"bi {icon} mv-help-tip-icon"),
            html.Span(text),
        ],
        className="mv-help-tip",
    )


def _settings_canvas() -> dbc.Offcanvas:
    """Station / sites settings — line visibility + per-site masking."""
    return dbc.Offcanvas(
        [
            html.Div(id=IDs.SET_SUMMARY, className="mv-set-summary"),
            html.Div("Active lines", className="mv-panel-lbl mt-2"),
            html.Div(
                "Toggle whole survey lines on/off across the map, "
                "3-D and table.",
                className="mv-crs-info",
            ),
            dbc.Checklist(
                id=IDs.SET_LINES,
                options=[],
                value=[],
                switch=True,
                className="mv-set-lines",
            ),
            html.Hr(className="mv-hr"),
            html.Div("Masked sites", className="mv-panel-lbl"),
            html.Div(
                "Masked stations are removed from every view (map, "
                "3-D and table). Toggle the eye on a station row, or "
                "unmask them below.",
                className="mv-crs-info",
            ),
            html.Div(
                [
                    dbc.Button(
                        [html.I(className="bi bi-eye me-1"), "Clear masks"],
                        id=IDs.BTN_CLEAR_MASKS,
                        color="secondary", outline=True, size="sm",
                        n_clicks=0,
                    ),
                    dbc.Button(
                        [html.I(className="bi bi-eye-slash me-1"),
                         "Mask hidden lines"],
                        id=IDs.BTN_MASK_HIDDEN,
                        color="secondary", outline=True, size="sm",
                        n_clicks=0,
                    ),
                ],
                className="d-flex gap-2 mb-2",
            ),
            html.Div(id=IDs.SET_MASKED_LIST, className="mv-masked-list"),
        ],
        id=IDs.CANVAS_SETTINGS,
        title="Station / Sites Settings",
        placement="end",
        is_open=False,
        style={"width": "340px"},
        className="mv-panel",
    )


def _session_canvas() -> dbc.Offcanvas:
    """Save / restore this workbench session.

    A session captures the current *view state* — active view (map /
    3-D), inspector controls, active lines, masked stations and theme —
    plus the loaded survey's station metadata for display. It does
    **not** capture the raw EDI bytes (they only ever live in the
    server-side view cache for this browser tab), so restoring a
    session always needs a "Load lines" reload afterwards to get the
    map/3-D canvas working again.
    """
    return dbc.Offcanvas(
        [
            html.Div(
                [
                    html.I(className="bi bi-info-circle me-1"),
                    "Sessions are auto-saved to your browser on every "
                    "change. Use Download / Upload to move a session to "
                    "another machine.",
                ],
                className="mv-crs-info mb-2",
            ),
            html.Div("Note", className="mv-panel-lbl"),
            dbc.Input(
                id=IDs.SESSION_NOTE,
                placeholder="e.g. WILLY AMT 2024 — pre-inversion",
                size="sm",
                className="mb-2",
            ),
            html.Hr(className="mv-hr"),
            html.Div("Save", className="mv-panel-lbl"),
            html.Div(
                [
                    dbc.Button(
                        [html.I(className="bi bi-download me-1"),
                         "Download JSON"],
                        id=IDs.BTN_SESSION_SAVE,
                        color="primary", size="sm", n_clicks=0,
                    ),
                    html.Div(id=IDs.SESSION_AUTOSAVE,
                             className="mv-topo-status"),
                ],
                className="d-flex align-items-center gap-2 mb-2",
            ),
            html.Details(
                [
                    html.Summary("What is included?"),
                    html.Ul(
                        [
                            html.Li("Active view (map / 3-D)"),
                            html.Li("Inspector controls (colours, mode, "
                                    "depth/ρ range, topo, basemap, …)"),
                            html.Li("Active lines & masked stations"),
                            html.Li("Theme"),
                            html.Li("Loaded station metadata (for display "
                                    "only — reload EDI files to re-enable "
                                    "the map/3-D canvas)"),
                        ],
                        className="mv-crs-info",
                        style={"paddingLeft": "16px"},
                    ),
                ],
                className="mb-2",
            ),
            html.Hr(className="mv-hr"),
            html.Div("Restore", className="mv-panel-lbl"),
            dcc.Upload(
                id=IDs.SESSION_UL,
                children=html.Div(
                    [
                        html.I(className="bi bi-file-earmark-arrow-up me-2"),
                        "Drop session JSON here or ",
                        html.A("browse"),
                    ],
                ),
                accept=".json",
                className="mv-upload-drop mv-topo-drop mb-2",
            ),
            dbc.Button(
                [html.I(className="bi bi-arrow-counterclockwise me-1"),
                 "Restore browser snapshot"],
                id=IDs.BTN_SESSION_RESTORE,
                color="secondary", outline=True, size="sm",
                className="w-100 mb-1", n_clicks=0,
            ),
            dbc.Button(
                [html.I(className="bi bi-trash me-1"), "Clear snapshot"],
                id=IDs.BTN_SESSION_CLEAR,
                color="danger", outline=True, size="sm",
                className="w-100", n_clicks=0,
            ),
            html.Div(id=IDs.SESSION_FEEDBACK, className="mv-topo-status mt-2"),
        ],
        id=IDs.CANVAS_SESSION,
        title="Workbench Session",
        placement="end",
        is_open=False,
        style={"width": "360px"},
        className="mv-panel",
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
                    html.Div("Core workflows", className="mv-help-section"),
                    html.Div(
                        [
                            _help_cap(
                                "bi-folder2-open",
                                "Load survey lines",
                                "Import EDI, AVG or J files by folder, "
                                "drag/drop, replace, or append.",
                            ),
                            _help_cap(
                                "bi-geo-alt",
                                "Station maps",
                                "Draw stations and line paths on selectable "
                                "satellite, street and topographic basemaps.",
                            ),
                            _help_cap(
                                "bi-layers-half",
                                "Contours & depth",
                                "Inspect resistivity, phase, elevation and "
                                "skin-depth layers by frequency.",
                            ),
                            _help_cap(
                                "bi-box",
                                "3-D model views",
                                "Move between fence, block, depth-slice and "
                                "iso-surface renderers.",
                            ),
                        ],
                        className="mv-help-caps",
                    ),
                    html.Div("Working cleanly", className="mv-help-section"),
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
                    html.Button(
                        "Close",
                        id=IDs.BTN_HELP_CLOSE,
                        className="mv-tbtn",
                        n_clicks=0,
                    ),
                ]
            ),
        ],
        id=IDs.MODAL_HELP,
        is_open=False,
        size="lg",
        centered=True,
        scrollable=True,
        className="mv-help-modal",
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
        dcc.Store(id=IDs.STORE_MASKED, data=[]),
        dcc.Store(id=IDs.LOAD_MODE_STORE, data="replace"),
        dcc.Store(id=IDs.FOLDER_STORE, data={}),
        dcc.Store(id=IDs.INV_FOLDER_STORE, data={}),
        dcc.Store(id=IDs.UPLOAD_SELECTION, data=[]),
        dcc.Store(id=IDs.SOURCE_SELECTION, data="none"),
        dcc.Store(id=IDs.TOPO_UPLOAD_STORE, data={}),
        dcc.Store(id=IDs.STORE_FIT, data=0),
        dcc.Download(id=IDs.EXPORT_DL),
        dcc.Store(id=IDs.SESSION_SNAPSHOT, storage_type="local", data=None),
        dcc.Download(id=IDs.SESSION_DL),
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
            _settings_canvas(),
            _session_canvas(),
            _topbar(),
            _welcome(),
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
