# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""3D Map page — fence diagram, volumetric block, and depth-slice visualizations."""

from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html

from pycsamt.app.web.layout import IDs, _command_bar

PAGE_ID = "map3d"

_CMAPS = [
    "RdYlBu_r",
    "RdBu_r",
    "seismic",
    "coolwarm",
    "jet",
    "viridis",
    "plasma",
    "inferno",
    "hot_r",
    "rainbow",
]

# ── Mode metadata ─────────────────────────────────────────────────────────────
_MAP3D_MODES = [
    ("fence", "Fence", "bi-layers"),
    ("block", "Block", "bi-box"),
    ("depth", "Slices", "bi-stack"),
]
_DEFAULT_MODE = "fence"

# ── Helpers ───────────────────────────────────────────────────────────────────


def _lbl(txt: str) -> html.Div:
    return html.Div(txt, className="fwd-feedback-mini mb-1")


def _ctrl_label(txt: str) -> html.Div:
    return html.Div(txt, className="ctrl-label")


def _filter_divider(txt: str) -> html.Div:
    return html.Div(txt, className="map3d-filter-section-label")


def _preset_bar(*btns) -> html.Div:
    return html.Div(list(btns), className="map3d-preset-bar")


def _preset(label: str, id_: str) -> html.Button:
    return html.Button(label, id=id_, className="map3d-preset-btn", n_clicks=0)


def _topo_panel() -> dbc.AccordionItem:
    """Topography accordion item — elevation source, terrain surface + z-shift."""
    return _settings_item(
        "Topography",
        "bi-geo-alt",
        [
            html.Div(
                "Drape terrain over the 3D model. "
                "Pick a source below, then re-generate 3D to apply.",
                className="fwd-feedback-mini mb-2",
                style={"color": "var(--sub0)", "fontStyle": "italic"},
            ),
            _lbl("Elevation source"),
            dbc.RadioItems(
                id=IDs.MAP3D_TOPO_SRC,
                options=[
                    {"label": "Off (no terrain)", "value": "none"},
                    {"label": "Station elevations", "value": "stations"},
                    {"label": "Elev. tool — corrected", "value": "elev-corr"},
                    {"label": "Elev. tool — raw fetch", "value": "elev-raw"},
                    {
                        "label": "Upload file (CSV / H5 / NPZ)",
                        "value": "upload",
                    },
                ],
                value="none",
                inputStyle={"cursor": "pointer"},
                labelStyle={
                    "fontSize": "12px",
                    "cursor": "pointer",
                    "padding": "2px 0",
                },
                className="mb-2",
            ),
            # ── Upload section (shown only when source = "upload") ──────────
            html.Div(
                id="map3d-topo-upload-wrap",
                children=[
                    dcc.Upload(
                        id=IDs.MAP3D_TOPO_UPLOAD,
                        children=html.Div(
                            [
                                html.I(className="bi bi-cloud-upload me-1"),
                                "Drop or ",
                                html.A(
                                    "browse",
                                    style={
                                        "textDecoration": "underline",
                                        "cursor": "pointer",
                                    },
                                ),
                                " a topography file",
                                html.Br(),
                                html.Span(
                                    ".csv  ·  .h5 / .hdf5  ·  .npz",
                                    style={
                                        "fontSize": "10px",
                                        "color": "var(--sub0)",
                                    },
                                ),
                            ]
                        ),
                        accept=".csv,.h5,.hdf5,.npz",
                        multiple=False,
                        style={
                            "border": "1px dashed var(--overlay0)",
                            "borderRadius": "6px",
                            "padding": "8px 6px",
                            "textAlign": "center",
                            "cursor": "pointer",
                            "fontSize": "12px",
                            "color": "var(--text)",
                            "marginBottom": "6px",
                        },
                    ),
                    html.Div(
                        id=IDs.MAP3D_TOPO_UPLOAD_INFO,
                        className="fwd-feedback-mini",
                        style={"minHeight": "16px"},
                    ),
                ],
                style={"display": "none"},  # hidden until source == "upload"
            ),
            html.Div(style={"marginTop": "8px"}),
            _lbl("Terrain opacity"),
            dcc.Slider(
                id=IDs.MAP3D_TOPO_OPACITY,
                min=0.1,
                max=1.0,
                step=0.05,
                value=0.70,
                marks={0.1: "10%", 0.5: "50%", 1.0: "100%"},
                tooltip={"always_visible": False},
                className="map3d-opacity-slider",
            ),
            html.Div(style={"marginTop": "10px"}),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Switch(
                            id=IDs.MAP3D_TOPO_STATIONS,
                            value=True,
                            label="Station markers on terrain",
                            className="map3d-switch",
                            style={"fontSize": "12px"},
                        ),
                        width=12,
                    ),
                ],
                className="mb-2",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Switch(
                            id=IDs.MAP3D_TOPO_APPLY,
                            value=False,
                            label="Shift model to follow terrain",
                            className="map3d-switch",
                            style={"fontSize": "12px"},
                        ),
                        width=12,
                    ),
                ]
            ),
            html.Div(
                "When enabled, each depth column is shifted by the surface "
                "elevation so the model follows real topography.",
                className="fwd-feedback-mini mt-1",
                style={"color": "var(--sub0)", "fontStyle": "italic"},
            ),
            # ── Station marker style ──────────────────────────────────────
            html.Hr(style={"borderColor": "var(--overlay0)", "margin": "10px 0"}),
            _lbl("Station marker"),
            html.Div(
                [
                    html.Span(
                        "Auto: ◇ open diamond (skin-depth) · "
                        "◆ filled diamond (inversion). "
                        "Plotly 3D does not support triangle markers.",
                        style={
                            "fontSize": "10px",
                            "color": "var(--sub0)",
                            "fontStyle": "italic",
                        },
                    ),
                ],
                className="mb-1",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            _lbl("Symbol"),
                            dbc.Select(
                                id=IDs.MAP3D_STA_SYMBOL,
                                options=[
                                    {
                                        "label": "Auto (by source)",
                                        "value": "auto",
                                    },
                                    {
                                        "label": "◇  Diamond open",
                                        "value": "diamond-open",
                                    },
                                    {
                                        "label": "◆  Diamond filled",
                                        "value": "diamond",
                                    },
                                    {
                                        "label": "○  Circle open",
                                        "value": "circle-open",
                                    },
                                    {
                                        "label": "●  Circle filled",
                                        "value": "circle",
                                    },
                                    {
                                        "label": "□  Square open",
                                        "value": "square-open",
                                    },
                                    {
                                        "label": "■  Square filled",
                                        "value": "square",
                                    },
                                    {"label": "+  Cross", "value": "cross"},
                                    {"label": "×  X mark", "value": "x"},
                                ],
                                value="auto",
                                size="sm",
                            ),
                        ],
                        width=7,
                    ),
                    dbc.Col(
                        [
                            _lbl("Size (px)"),
                            dbc.Input(
                                id=IDs.MAP3D_STA_SIZE,
                                type="number",
                                value=8,
                                min=3,
                                max=30,
                                step=1,
                                size="sm",
                                debounce=True,
                            ),
                        ],
                        width=5,
                    ),
                ],
                className="g-1 mb-2",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            _lbl("Colour"),
                            dbc.Input(
                                id=IDs.MAP3D_STA_COLOR,
                                type="color",
                                value="#000000",
                                style={
                                    "height": "28px",
                                    "padding": "1px 3px",
                                    "cursor": "pointer",
                                },
                            ),
                        ],
                        width=4,
                    ),
                    dbc.Col(
                        html.Div(
                            "Applies only when Show stations is on.",
                            className="fwd-feedback-mini mt-3",
                            style={
                                "color": "var(--sub0)",
                                "fontStyle": "italic",
                            },
                        ),
                        width=8,
                    ),
                ],
                className="g-1",
            ),
        ],
        "map3d-settings-topo",
    )


def _settings_item(title: str, icon: str, children, item_id: str) -> dbc.AccordionItem:
    return dbc.AccordionItem(
        children,
        title=html.Span(
            [
                html.I(className=f"bi {icon} me-2"),
                title,
            ],
            className="map3d-settings-title",
        ),
        item_id=item_id,
    )


# ── Sub-components ────────────────────────────────────────────────────────────


def _run_bar() -> html.Div:
    return html.Div(
        [
            dbc.Button(
                [html.I(className="bi bi-box me-2"), "Load & Generate 3D"],
                id=IDs.BTN_MAP3D_RUN,
                color="primary",
                size="sm",
                className="w-100 mb-2",
            ),
            html.Div(
                "Controls update live after first generation.",
                className="fwd-feedback-mini mb-1",
            ),
            dbc.Spinner(html.Div(id=IDs.MAP3D_SPINNER), size="sm", color="primary"),
        ],
        className="fwd-run-bar",
    )


def _mode_bar() -> html.Div:
    btns = []
    for slug, label, icon in _MAP3D_MODES:
        active = slug == _DEFAULT_MODE
        btns.append(
            html.Button(
                [html.I(className=f"bi {icon} me-1"), label],
                id=f"map3d-mode-btn-{slug}",
                className="fwd-dim-btn" + (" active" if active else ""),
                n_clicks=0,
            )
        )
    return html.Div(btns, className="fwd-dim-bar")


def _fence_panel() -> html.Div:
    return html.Div(
        [
            _ctrl_label("Fence Options"),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            _lbl("Line spacing (km)"),
                            dbc.Input(
                                id=IDs.MAP3D_LINE_SPACING,
                                type="number",
                                value=1.0,
                                min=0.1,
                                max=10.0,
                                step=0.1,
                                size="sm",
                                debounce=True,
                            ),
                        ]
                    ),
                    dbc.Col(
                        [
                            _lbl("Azimuth offset (°)"),
                            dbc.Input(
                                id=IDs.MAP3D_AZIMUTH,
                                type="number",
                                value=0,
                                min=-180,
                                max=180,
                                step=5,
                                size="sm",
                                debounce=True,
                            ),
                        ]
                    ),
                ],
                className="g-1 mb-2",
            ),
            _lbl("Panel visual thickness (m)"),
            dbc.Input(
                id=IDs.MAP3D_PANEL_THICK,
                type="number",
                value=50,
                min=1,
                max=500,
                step=10,
                size="sm",
                debounce=True,
            ),
        ],
        id="map3d-fence-panel",
        className="map3d-mode-panel",
    )


def _block_panel() -> html.Div:
    return html.Div(
        [
            _ctrl_label("Volume Block Options"),
            _lbl("Iso-surface count"),
            dbc.Input(
                id=IDs.MAP3D_N_SURFACES,
                type="number",
                value=10,
                min=2,
                max=50,
                step=1,
                size="sm",
                debounce=True,
                className="mb-2",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Span("Clip X", className="fwd-bound-label me-1"),
                            dcc.Slider(
                                id=IDs.MAP3D_CLIP_X,
                                min=0,
                                max=1,
                                step=0.05,
                                value=1,
                                marks=None,
                                tooltip={"always_visible": False},
                                className="map3d-clip-slider",
                            ),
                        ]
                    ),
                    html.Div(
                        [
                            html.Span("Clip Y", className="fwd-bound-label me-1"),
                            dcc.Slider(
                                id=IDs.MAP3D_CLIP_Y,
                                min=0,
                                max=1,
                                step=0.05,
                                value=1,
                                marks=None,
                                tooltip={"always_visible": False},
                                className="map3d-clip-slider",
                            ),
                        ]
                    ),
                    html.Div(
                        [
                            html.Span(
                                "Clip Z (depth)",
                                className="fwd-bound-label me-1",
                            ),
                            dcc.Slider(
                                id=IDs.MAP3D_CLIP_Z,
                                min=0,
                                max=1,
                                step=0.05,
                                value=1,
                                marks=None,
                                tooltip={"always_visible": False},
                                className="map3d-clip-slider",
                            ),
                        ]
                    ),
                ],
                className="mt-2",
            ),
        ],
        id="map3d-block-panel",
        className="map3d-mode-panel",
        style={"display": "none"},
    )


def _depth_panel() -> html.Div:
    return html.Div(
        [
            _ctrl_label("Depth Slice Options"),
            _lbl("Number of depth slices"),
            dbc.Input(
                id=IDs.MAP3D_N_SLICES,
                type="number",
                value=8,
                min=2,
                max=40,
                step=1,
                size="sm",
                debounce=True,
            ),
        ],
        id="map3d-depth-panel",
        className="map3d-mode-panel",
        style={"display": "none"},
    )


def _controls_scroll() -> html.Div:
    return html.Div(
        [
            dbc.Accordion(
                [
                    _settings_item(
                        "Data Source",
                        "bi-database",
                        [
                            dbc.RadioItems(
                                id=IDs.MAP3D_DATA_SRC,
                                options=[
                                    {
                                        "label": html.Span(
                                            [
                                                html.I(
                                                    className="bi bi-database me-2",
                                                    style={"color": "#89b4fa"},
                                                ),
                                                "Skin-depth pseudo",
                                            ],
                                            className="d-flex align-items-center",
                                        ),
                                        "value": "pseudo",
                                    },
                                    {
                                        "label": html.Span(
                                            [
                                                html.I(
                                                    className="bi bi-layers me-2",
                                                    style={"color": "#a6e3a1"},
                                                ),
                                                "Session inversion result",
                                            ],
                                            className="d-flex align-items-center",
                                        ),
                                        "value": "inversion",
                                    },
                                    {
                                        "label": html.Span(
                                            [
                                                html.I(
                                                    className="bi bi-diagram-3 me-2",
                                                    style={"color": "#fab387"},
                                                ),
                                                "Survey line profiles",
                                            ],
                                            className="d-flex align-items-center",
                                        ),
                                        "value": "profiles",
                                    },
                                ],
                                value="pseudo",
                                inputStyle={"cursor": "pointer"},
                                labelStyle={
                                    "fontSize": "12px",
                                    "cursor": "pointer",
                                    "padding": "3px 0",
                                },
                            ),
                            html.Div(
                                "Skin-depth uses loaded EDI stations. Survey line "
                                "profiles groups those stations by line metadata. "
                                "Session inversion uses the current inversion model "
                                "when one is available.",
                                className="fwd-feedback-mini mt-2",
                            ),
                        ],
                        "map3d-settings-source",
                    ),
                    _settings_item(
                        "View Filters",
                        "bi-funnel",
                        [
                            _filter_divider("Depth Window (m)"),
                            _preset_bar(
                                _preset("Full", "map3d-depth-preset-full"),
                                _preset("500 m", "map3d-depth-preset-500"),
                                _preset("1 km", "map3d-depth-preset-1k"),
                                _preset("2 km", "map3d-depth-preset-2k"),
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            _lbl("Min depth"),
                                            dbc.Input(
                                                id=IDs.MAP3D_DEPTH_LO,
                                                type="number",
                                                value=0,
                                                min=0,
                                                max=50000,
                                                step="any",
                                                size="sm",
                                                debounce=True,
                                            ),
                                        ]
                                    ),
                                    dbc.Col(
                                        [
                                            _lbl("Max depth"),
                                            dbc.Input(
                                                id=IDs.MAP3D_DEPTH_HI,
                                                type="number",
                                                value=2000,
                                                min=10,
                                                max=50000,
                                                step="any",
                                                size="sm",
                                                debounce=True,
                                            ),
                                        ]
                                    ),
                                ],
                                className="g-1 mb-3",
                            ),
                            _filter_divider("ρ Visibility Filter (Ω·m)"),
                            _preset_bar(
                                _preset("All", "map3d-rho-preset-all"),
                                _preset("Conduct.", "map3d-rho-preset-cond"),
                                _preset("Mid", "map3d-rho-preset-mid"),
                                _preset("Resist.", "map3d-rho-preset-res"),
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            _lbl("ρ min (Ω·m)"),
                                            dbc.Input(
                                                id=IDs.MAP3D_RHO_LO,
                                                type="number",
                                                value=1,
                                                min=0.001,
                                                step="any",
                                                size="sm",
                                                debounce=True,
                                            ),
                                        ]
                                    ),
                                    dbc.Col(
                                        [
                                            _lbl("ρ max (Ω·m)"),
                                            dbc.Input(
                                                id=IDs.MAP3D_RHO_HI,
                                                type="number",
                                                value=100_000,
                                                min=1,
                                                step="any",
                                                size="sm",
                                                debounce=True,
                                            ),
                                        ]
                                    ),
                                ],
                                className="g-1 mb-2",
                            ),
                            html.Div(
                                id=IDs.MAP3D_RANGE_HINT,
                                className="fwd-feedback-mini mt-1",
                                style={
                                    "color": "var(--sub0)",
                                    "fontStyle": "italic",
                                },
                            ),
                        ],
                        "map3d-settings-filters",
                    ),
                    _settings_item(
                        "Color Scale",
                        "bi-palette",
                        [
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            _lbl("Colormap"),
                                            dbc.Select(
                                                id=IDs.MAP3D_CMAP,
                                                options=[
                                                    {"label": c, "value": c}
                                                    for c in _CMAPS
                                                ],
                                                value="RdYlBu_r",
                                                size="sm",
                                            ),
                                        ],
                                        width=7,
                                    ),
                                    dbc.Col(
                                        [
                                            _lbl("Scale"),
                                            dbc.Select(
                                                id=IDs.MAP3D_SCALE,
                                                options=[
                                                    {
                                                        "label": "Log",
                                                        "value": "log",
                                                    },
                                                    {
                                                        "label": "Linear",
                                                        "value": "linear",
                                                    },
                                                ],
                                                value="log",
                                                size="sm",
                                            ),
                                        ],
                                        width=5,
                                    ),
                                ],
                                className="g-1 mb-2",
                            ),
                            _lbl("Colorbar range (Ω·m)"),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            _lbl("min"),
                                            dbc.Input(
                                                id=IDs.MAP3D_VMIN,
                                                type="number",
                                                value=1,
                                                min=0.001,
                                                step="any",
                                                size="sm",
                                                debounce=True,
                                            ),
                                        ]
                                    ),
                                    dbc.Col(
                                        [
                                            _lbl("max"),
                                            dbc.Input(
                                                id=IDs.MAP3D_VMAX,
                                                type="number",
                                                value=10_000,
                                                min=1,
                                                step="any",
                                                size="sm",
                                                debounce=True,
                                            ),
                                        ]
                                    ),
                                ],
                                className="g-1",
                            ),
                        ],
                        "map3d-settings-color",
                    ),
                    _settings_item(
                        "Rendering",
                        "bi-sliders",
                        [
                            _lbl("Surface / volume opacity"),
                            dcc.Slider(
                                id=IDs.MAP3D_OPACITY,
                                min=0.05,
                                max=1.0,
                                step=0.05,
                                value=0.85,
                                marks={0.05: "5%", 0.5: "50%", 1.0: "100%"},
                                tooltip={"always_visible": False},
                                className="map3d-opacity-slider",
                            ),
                            html.Div(style={"marginTop": "8px"}),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        dbc.Switch(
                                            id=IDs.MAP3D_CONTOURS,
                                            value=True,
                                            label="Contours",
                                            className="map3d-switch",
                                            style={"fontSize": "12px"},
                                        ),
                                        width=6,
                                    ),
                                    dbc.Col(
                                        [
                                            _lbl("N lines"),
                                            dbc.Input(
                                                id=IDs.MAP3D_N_CONTOURS,
                                                type="number",
                                                value=15,
                                                min=2,
                                                max=50,
                                                step=1,
                                                size="sm",
                                                debounce=True,
                                            ),
                                        ],
                                        width=6,
                                    ),
                                ],
                                className="g-1",
                            ),
                        ],
                        "map3d-settings-rendering",
                    ),
                    _settings_item(
                        "Mode Options",
                        "bi-boxes",
                        [
                            _fence_panel(),
                            _block_panel(),
                            _depth_panel(),
                            html.Div(
                                [
                                    _ctrl_label("Block Volume Band (log₁₀ Ω·m)"),
                                    html.Div(
                                        "Controls the visible resistivity interval for "
                                        "the semi-transparent volume.",
                                        className="fwd-feedback-mini mb-2",
                                        style={"color": "var(--sub0)"},
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    _lbl("Band lo"),
                                                    dbc.Input(
                                                        id=IDs.MAP3D_ISOVALUE_LO,
                                                        type="number",
                                                        value=0.5,
                                                        min=-1,
                                                        max=5,
                                                        step=0.1,
                                                        size="sm",
                                                        debounce=True,
                                                    ),
                                                ]
                                            ),
                                            dbc.Col(
                                                [
                                                    _lbl("Band hi"),
                                                    dbc.Input(
                                                        id=IDs.MAP3D_ISOVALUE_HI,
                                                        type="number",
                                                        value=2.5,
                                                        min=-1,
                                                        max=5,
                                                        step=0.1,
                                                        size="sm",
                                                        debounce=True,
                                                    ),
                                                ]
                                            ),
                                        ],
                                        className="g-1",
                                    ),
                                ],
                                id="map3d-iso-band-card",
                                className="map3d-mode-panel",
                                style={"display": "none"},
                            ),
                        ],
                        "map3d-settings-mode",
                    ),
                    _settings_item(
                        "Annotations",
                        "bi-type",
                        [
                            dbc.Switch(
                                id=IDs.MAP3D_SHOW_LABELS,
                                value=True,
                                label="Profile / line labels",
                                style={"fontSize": "12px"},
                            ),
                            html.Div(style={"marginTop": "6px"}),
                            _lbl("Figure title (optional)"),
                            dbc.Input(
                                id=IDs.MAP3D_TITLE,
                                type="text",
                                placeholder="e.g. Pandian Area - resistivity",
                                size="sm",
                                debounce=True,
                            ),
                        ],
                        "map3d-settings-annotations",
                    ),
                    _topo_panel(),
                    _settings_item(
                        "Export",
                        "bi-download",
                        [
                            dbc.Row(
                                [
                                    dbc.Col(
                                        dbc.Button(
                                            [
                                                html.I(className="bi bi-image me-1"),
                                                "PNG",
                                            ],
                                            id=IDs.BTN_MAP3D_EXPORT_PNG,
                                            color="secondary",
                                            size="sm",
                                            outline=True,
                                            className="w-100",
                                        )
                                    ),
                                    dbc.Col(
                                        dbc.Button(
                                            [
                                                html.I(
                                                    className="bi bi-filetype-html me-1"
                                                ),
                                                "HTML",
                                            ],
                                            id=IDs.BTN_MAP3D_EXPORT_HTML,
                                            color="secondary",
                                            size="sm",
                                            outline=True,
                                            className="w-100",
                                        )
                                    ),
                                ],
                                className="g-2",
                            ),
                            html.Div(
                                "PNG uses Plotly's toolbar. HTML exports an "
                                "interactive self-contained file.",
                                className="fwd-feedback-mini mt-2",
                            ),
                        ],
                        "map3d-settings-export",
                    ),
                ],
                active_item=["map3d-settings-source"],
                always_open=True,
                className="map3d-settings-accordion",
            ),
        ],
        className="fwd-ctrl-scroll",
    )


def layout() -> html.Div:
    store = dcc.Store(id=IDs.MAP3D_ACTIVE_MODE, data=_DEFAULT_MODE)
    grid_store = dcc.Store(id=IDs.MAP3D_GRID_STORE)

    sidebar = html.Div(
        [store, grid_store, _run_bar(), _mode_bar(), _controls_scroll()],
        className="analysis-controls fwd-sidebar",
    )

    info_bar = html.Div(
        id=IDs.MAP3D_INFO,
        className="map3d-info-strip",
        children=html.Span(
            "Select mode + data source, then click Load & Generate 3D.",
            style={"fontSize": "11px", "color": "#6c7086"},
        ),
    )

    graph_area = html.Div(
        [
            info_bar,
            dcc.Graph(
                id=IDs.MAP3D_GRAPH,
                figure=_empty_3d_fig(),
                config={
                    "displayModeBar": True,
                    "modeBarButtonsToAdd": [
                        "hoverClosest3d",
                        "orbitRotation",
                        "tableRotation",
                        "resetCameraLastSave3d",
                    ],
                    "toImageButtonOptions": {
                        "format": "png",
                        "filename": "map3d",
                        "height": 900,
                        "width": 1400,
                        "scale": 2,
                    },
                    "responsive": True,
                },
                style={"height": "100%", "width": "100%"},
                className="map3d-graph",
            ),
        ],
        className="analysis-output map3d-output",
        style={"display": "flex", "flexDirection": "column"},
    )

    return html.Div(
        [
            _command_bar(
                "map3d",
                "3D Resistivity Map",
                "Fence · Volume · Depth slices "
                "· ρ filter · Depth filter · Opacity · Contours · PNG / HTML",
            ),
            html.Div(
                [sidebar, graph_area],
                className="analysis-layout",
                style={"flex": "1"},
            ),
        ],
        style={
            "display": "flex",
            "flexDirection": "column",
            "height": "100%",
        },
    )


def _empty_3d_fig(dark: bool = True):
    import plotly.graph_objects as go

    text_color = "#6c7086" if dark else "#6c6f85"
    fig = go.Figure()
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                visible=False,
                showgrid=False,
                zeroline=False,
                showticklabels=False,
                title="",
                backgroundcolor="rgba(0,0,0,0)",
            ),
            yaxis=dict(
                visible=False,
                showgrid=False,
                zeroline=False,
                showticklabels=False,
                title="",
                backgroundcolor="rgba(0,0,0,0)",
            ),
            zaxis=dict(
                visible=False,
                showgrid=False,
                zeroline=False,
                showticklabels=False,
                title="",
                backgroundcolor="rgba(0,0,0,0)",
            ),
            bgcolor="rgba(0,0,0,0)",
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=0, r=0, t=0, b=0),
        annotations=[
            dict(
                text="Load & Generate 3D to begin.",
                x=0.5,
                y=0.5,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(color=text_color, size=13),
            )
        ],
    )
    return fig


def register_callbacks(app) -> None:
    pass  # handled by callbacks/map3d.py
