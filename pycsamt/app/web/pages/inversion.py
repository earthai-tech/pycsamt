# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Inversion page — Traditional (builtin/Occam2D) + AI Neural (1D/2D/3D).

Layout
------
Sidebar:
  • Sticky Run-Inversion bar (always visible)
  • Track selector (Traditional · AI Neural) as pill buttons
  • Scrollable per-track controls + shared freq range

View area (right):
  • Fixed output tab bar: Result · Convergence · Statistics · Log
  • Context info bar
  • Per-tab panels (persistent content; switching tabs preserves last result)
  • PyTorch install modal (overlay)
"""
from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html

try:
    import torch as _torch
    _TORCH_AVAIL = True
except ImportError:
    _TORCH_AVAIL = False

from pycsamt.app.web.layout import IDs, _command_bar
from pycsamt.app.web.utils import empty_src

PAGE_ID = "inversion"

# ── Selector option lists ─────────────────────────────────────────────────────

_METHOD_OPTS = [
    {"label": "MT  (Magnetotelluric)",      "value": "mt"},
    {"label": "CSAMT  (Controlled-source)", "value": "csamt"},
    {"label": "TDEM  (Time-domain EM)",     "value": "tdem"},
]
_BACKEND_OPTS = [
    {"label": "Built-in  (SciPy, no binaries) [Recommended]", "value": "builtin"},
    {"label": "Occam2D  (requires external binary)",           "value": "occam2d"},
    {"label": "ModEM    (requires external binary)",           "value": "modem"},
]
_REGULARIZE_OPTS = [
    {"label": "Smooth  (Occam-style L2)", "value": "smooth"},
    {"label": "Damped  (Tikhonov)",       "value": "damped"},
    {"label": "Blocky  (L1/sparse)",      "value": "blocky"},
    {"label": "None",                     "value": "none"},
]
_AI_DIM_OPTS = [
    {"label": "1-D Sounding  (ResNet/CNN/FCN)",   "value": "1d"},
    {"label": "2-D Profile   (U-Net)",            "value": "2d"},
    {"label": "3-D Volume    (Graph Conv. Net)",  "value": "3d"},
]
_AI_ARCH_1D_OPTS = [
    {"label": "ResNet1D  [recommended]", "value": "resnet"},
    {"label": "CNN1D",                   "value": "cnn1d"},
    {"label": "FCN1D",                   "value": "fcn"},
]
_AI_SOLVER_OPTS = [
    {"label": "MT 1-D",    "value": "mt1d"},
    {"label": "CSAMT 1-D", "value": "csamt1d"},
    {"label": "TEM 1-D",   "value": "tem1d"},
]

_DEFAULT_TRUE_LAYERS = [
    {"layer": 1, "resistivity": 50,   "thickness": 200},
    {"layer": 2, "resistivity": 500,  "thickness": 800},
    {"layer": 3, "resistivity": 20,   "thickness": 1000},
]

_COMP_OPTS = [
    {"label": "xy", "value": "xy"},
    {"label": "yx", "value": "yx"},
    {"label": "xx", "value": "xx"},
    {"label": "yy", "value": "yy"},
]
_DEVICE_OPTS = [
    {"label": "Auto-detect", "value": "auto"},
    {"label": "CPU",          "value": "cpu"},
    {"label": "CUDA (GPU)",   "value": "cuda"},
    {"label": "MPS (Apple)",  "value": "mps"},
]
_PINN_MODE_OPTS = [
    {"label": "TE",         "value": "te"},
    {"label": "TM",         "value": "tm"},
    {"label": "Both (avg)", "value": "both"},
]
_BACKEND_DL_OPTS = [
    {"label": "Auto-detect", "value": "auto"},
    {"label": "PyTorch",     "value": "torch"},
    {"label": "TensorFlow",  "value": "tensorflow"},
]
_UNET_D_OPTS = [
    {"label": "Auto",     "value": "auto"},
    {"label": "1 stage",  "value": "1"},
    {"label": "2 stages", "value": "2"},
    {"label": "3 stages", "value": "3"},
    {"label": "4 stages", "value": "4"},
]
_AI_COMP_OPTS = [
    {"label": "2  (TE mode)", "value": "2"},
    {"label": "4  (TE + TM)", "value": "4"},
]
_INV_TRACKS = [
    ("trad",   "Traditional", "bi-cpu"),
    ("ai",     "AI Neural",   "bi-robot"),
    ("pinn",   "PINN",        "bi-bezier2"),
    ("hybrid", "Hybrid",      "bi-node-plus"),
]
_INV_OUT_TABS = [
    ("result",   "Result",      "bi-image"),
    ("conv",     "Convergence", "bi-graph-down-arrow"),
    ("stats",    "Statistics",  "bi-table"),
    ("log",      "Log",         "bi-terminal"),
    ("data-fit", "Data Fit",    "bi-activity"),
]
_DEFAULT_TRACK = "trad"
_DEFAULT_OUT   = "result"

_IMG_STYLE = {
    "width": "100%", "height": "100%",
    "objectFit": "contain",
}


# ── Small helpers ─────────────────────────────────────────────────────────────

def _num(id_, val, **kw):
    return dbc.Input(
        id=id_, type="number",
        value=val, size="sm", **kw
    )


def _lbl(text):
    return html.Small(text, className="text-muted")


def _clabel(text):
    return html.Div(text, className="ctrl-label")


def _acc_item(
    title: str,
    icon: str,
    children,
    item_id: str,
) -> dbc.AccordionItem:
    return dbc.AccordionItem(
        children,
        title=html.Span([
            html.I(className=f"bi {icon} me-2"),
            title,
        ], className="inv-acc-title"),
        item_id=item_id,
    )


_GRP_LABEL_STYLE = {
    "fontSize": "9.5px",
    "color": "var(--sub0, #6c7086)",
    "textTransform": "uppercase",
    "letterSpacing": "0.07em",
    "marginBottom": "4px",
    "paddingLeft": "2px",
}


# ── True model DataTable (in sidebar) ────────────────────────────────────────

def _true_model_table() -> html.Div:
    tbl = dash_table.DataTable(
        id=IDs.INV_TRUE_TABLE,
        columns=[
            {"name": "#",              "id": "layer",       "editable": False},
            {"name": "ρ true (Ω·m)",  "id": "resistivity", "editable": True, "type": "numeric"},
            {"name": "Thickness (m)", "id": "thickness",   "editable": True, "type": "numeric"},
        ],
        data=_DEFAULT_TRUE_LAYERS,
        editable=True,
        row_deletable=True,
        style_table={"overflowY": "auto", "maxHeight": "160px"},
        style_header={
            "backgroundColor": "#11111b", "color": "#6c7086",
            "fontSize": "11px", "border": "none",
            "textTransform": "uppercase", "fontWeight": "600",
        },
        style_cell={
            "backgroundColor": "#181825", "color": "#cdd6f4",
            "border": "1px solid #313244", "fontSize": "11px", "padding": "4px 8px",
        },
        style_data_conditional=[
            {"if": {"state": "active"},
             "backgroundColor": "#313244", "border": "1px solid #89b4fa"},
        ],
    )
    return html.Div(tbl, className="inv-layer-table")


# ── Traditional inversion sidebar section ────────────────────────────────────

def _trad_controls():
    acc = dbc.Accordion(
        [
            _acc_item(
                "Method & Solver",
                "bi-cpu",
                [
                    _lbl("EM method"),
                    dbc.Select(
                        id=IDs.INV_METHOD,
                        options=_METHOD_OPTS,
                        value="mt",
                        size="sm",
                        className="mb-2",
                    ),
                    _lbl("Backend"),
                    dbc.Select(
                        id=IDs.INV_BACKEND,
                        options=_BACKEND_OPTS,
                        value="builtin",
                        size="sm",
                        className="mb-2",
                    ),
                    _lbl("Dimension"),
                    dbc.RadioItems(
                        id=IDs.INV_T_DIM,
                        options=[
                            {
                                "label": "1-D sounding",
                                "value": "1d",
                            },
                            {
                                "label": "2-D profile",
                                "value": "2d",
                            },
                        ],
                        value="1d",
                        inline=True,
                        className=(
                            "text-muted small mt-1"
                        ),
                    ),
                ],
                "trad-acc-method",
            ),

            _acc_item(
                "Data Source",
                "bi-database",
                [
                    dbc.RadioItems(
                        id=IDs.INV_DATA_SRC,
                        options=[
                            {
                                "label": "Synthetic"
                                " (define true model)",
                                "value": "synthetic",
                            },
                            {
                                "label": "Session"
                                " (EDIs / loaded sites)",
                                "value": "session",
                            },
                        ],
                        value="synthetic",
                        inline=False,
                        className="text-muted small",
                    ),
                    html.Hr(className="inv-acc-hr"),
                    # Synthetic panel (conditional)
                    html.Div([
                        _lbl("True Earth Model"),
                        html.Div(
                            "Synthetic forward data source",
                            className="inv-sub-hint mb-1",
                        ),
                        _true_model_table(),
                        dbc.Button(
                            [
                                html.I(
                                    className=(
                                        "bi bi-plus-lg"
                                        " me-1"
                                    )
                                ),
                                "Add Layer",
                            ],
                            id=IDs.BTN_INV_ADD_LAYER,
                            color="secondary",
                            size="sm",
                            className="mt-2",
                        ),
                        html.Div([
                            html.Div([
                                html.Span(
                                    [
                                        html.I(
                                            className=(
                                                "bi "
                                                "bi-infinity"
                                                " me-1"
                                            )
                                        ),
                                        "Halfspace",
                                    ],
                                    className="fwd-hs-label",
                                ),
                                html.Span(
                                    " half-space",
                                    style={
                                        "fontSize": "10px",
                                        "color": "var(--sub0)",
                                    },  # noqa: E501
                                ),
                            ], className="mb-1"),
                            dbc.InputGroup([
                                dbc.InputGroupText(
                                    "rho",
                                    style={
                                        "fontSize": "11px"
                                    },
                                ),
                                _num(
                                    IDs.INV_HALFSPACE_RHO,
                                    1000, min=0.01, step=10,
                                ),
                                dbc.InputGroupText(
                                    "Ohm*m",
                                    style={
                                        "fontSize": "11px"
                                    },
                                ),
                            ], size="sm"),
                        ], className="fwd-halfspace-row"),
                        dbc.Row([
                            dbc.Col([
                                _lbl("Noise level (sigma)"),
                                _num(
                                    IDs.INV_NOISE_LEVEL,
                                    0.05, min=0,
                                    max=0.5, step=0.01,
                                ),
                            ]),
                        ], className="g-1 mt-2"),
                        html.Div([
                            html.Hr(className="inv-acc-hr"),
                            _lbl("2-D Profile"),
                            dbc.Row([
                                dbc.Col([
                                    _lbl("N stations"),
                                    _num(
                                        IDs.INV_N_STATIONS,
                                        8, min=3,
                                        max=30, step=1,
                                    ),
                                ]),
                                dbc.Col([
                                    _lbl("Profile (m)"),
                                    _num(
                                        IDs.INV_PROFILE_LEN,
                                        5000, min=500,
                                        step=500,
                                    ),
                                ]),
                            ], className="g-1"),
                        ], id=IDs.INV_2D_PARAMS,
                           style={"display": "none"}),
                    ], id=IDs.INV_SYN_PANEL),
                    # Session panel (conditional)
                    html.Div([
                        html.I(
                            className=(
                                "bi bi-info-circle me-1"
                            ),
                            style={"color": "var(--blue)"},
                        ),
                        html.Span(
                            "Using sites from the"
                            " QC / Correction page.",
                            className="inv-sub-hint",
                        ),
                        html.Div(
                            "Load EDI files first if"
                            " no session data available.",
                            className="inv-sub-hint mt-1",
                            style={"color": "var(--sub0)"},
                        ),
                    ], id=IDs.INV_SESS_PANEL,
                       style={"display": "none"}),
                ],
                "trad-acc-data",
            ),

            _acc_item(
                "Configuration",
                "bi-sliders",
                [
                    dbc.Row([
                        dbc.Col([
                            _lbl("N layers (recover)"),
                            _num(
                                IDs.INV_N_LAYERS,
                                4, min=2, max=20, step=1,
                            ),
                        ]),
                        dbc.Col([
                            _lbl("Max iterations"),
                            _num(
                                IDs.INV_MAX_ITER,
                                60, min=10,
                                max=300, step=10,
                            ),
                        ]),
                    ], className="g-1 mb-2"),
                    dbc.Row([
                        dbc.Col([
                            _lbl("Error floor (%)"),
                            _num(
                                IDs.INV_ERROR_FLOOR,
                                0.05, min=0.001,
                                max=0.5, step=0.005,
                            ),
                        ]),
                    ], className="g-1 mb-2"),
                    _lbl("Regularization"),
                    dbc.Select(
                        id=IDs.INV_REGULARIZE,
                        options=_REGULARIZE_OPTS,
                        value="smooth",
                        size="sm",
                        className="mt-1 mb-2",
                    ),
                    dbc.Switch(
                        id=IDs.INV_INCLUDE_PHASE,
                        label="Include phase in objective",
                        value=True,
                        className="small",
                    ),
                ],
                "trad-acc-config",
            ),
        ],
        active_item=["trad-acc-method", "trad-acc-data"],
        always_open=True,
        className="inv-pinn-accordion",
        flush=True,
    )
    return html.Div(acc, id="inv-ctrl-trad")


# ── AI inversion sidebar section ──────────────────────────────────────────────

def _ai_controls():
    torch_badge = (
        dbc.Badge(
            "Installed",
            color="success",
            className="ms-1 py-0",
        )
        if _TORCH_AVAIL else
        dbc.Badge(
            "Not installed — click Run to install",
            color="warning",
            className="ms-1 py-0",
            style={"fontSize": "9px"},
        )
    )
    acc = dbc.Accordion(
        [
            _acc_item(
                "Setup",
                "bi-robot",
                [
                    html.Div(
                        [
                            html.I(
                                className=(
                                    "bi bi-cpu me-1"
                                ),
                                style={
                                    "color": "var(--green)"
                                },
                            ),
                            html.Span("PyTorch"),
                            torch_badge,
                        ],
                        className="inv-torch-badge mb-2",
                    ),
                    _lbl("Problem dimension"),
                    dbc.Select(
                        id=IDs.INV_AI_DIM,
                        options=_AI_DIM_OPTS,
                        value="1d",
                        size="sm",
                        className="mb-2",
                    ),
                    _lbl("Architecture"),
                    dbc.Select(
                        id=IDs.INV_AI_ARCH,
                        options=_AI_ARCH_1D_OPTS,
                        value="resnet",
                        size="sm",
                        className="mb-1",
                    ),
                    html.Div(
                        id="inv-ai-arch-note",
                        className="inv-sub-hint mb-2",
                    ),
                    _lbl("Forward solver (training data)"),
                    dbc.Select(
                        id=IDs.INV_AI_SOLVER,
                        options=_AI_SOLVER_OPTS,
                        value="mt1d",
                        size="sm",
                    ),
                ],
                "ai-acc-setup",
            ),

            _acc_item(
                "Network Config",
                "bi-diagram-3",
                [
                    # 2-D UNet panel (conditional)
                    html.Div([
                        html.Div(
                            [
                                html.I(
                                    className=(
                                        "bi bi-grid-3x3"
                                        "-gap me-1"
                                    ),
                                    style={
                                        "color": "var(--blue)"
                                    },
                                ),
                                html.Span(
                                    "2-D UNet",
                                    className=(
                                        "ctrl-label d-inline"
                                    ),
                                ),
                            ],
                            className="mb-2",
                        ),
                        dbc.Row([
                            dbc.Col([
                                _lbl("Input channels"),
                                dbc.Select(
                                    id=IDs.INV_AI_N_COMPONENTS,
                                    options=_AI_COMP_OPTS,
                                    value="2",
                                    size="sm",
                                ),
                            ]),
                            dbc.Col([
                                _lbl("Depth cells"),
                                _num(
                                    IDs.INV_AI_N_DEPTH,
                                    4, min=2, max=80, step=1,
                                ),
                            ]),
                        ], className="g-1 mb-1"),
                        dbc.Row([
                            dbc.Col([
                                _lbl("Stations/profile"),
                                _num(
                                    IDs.INV_AI_N_STATIONS,
                                    8, min=4,
                                    max=100, step=1,
                                ),
                            ]),
                            dbc.Col([
                                _lbl("UNet depth"),
                                dbc.Select(
                                    id=IDs.INV_AI_UNET_DEPTH,
                                    options=_UNET_D_OPTS,
                                    value="auto",
                                    size="sm",
                                ),
                            ]),
                        ], className="g-1 mb-1"),
                        dbc.Alert(
                            "Auto depth: max safe"
                            " pooling stages.",
                            color="dark",
                            className=(
                                "mb-0 py-1 px-2 inv-alert"
                            ),
                        ),
                    ], id=IDs.INV_AI_2D_PANEL,
                       style={"display": "none"}),

                    # 3-D GCN panel (conditional)
                    html.Div([
                        html.Div(
                            [
                                html.I(
                                    className=(
                                        "bi bi-diagram-3 me-1"
                                    ),
                                    style={
                                        "color":
                                        "var(--mauve)"
                                    },
                                ),
                                html.Span(
                                    "3-D GCN",
                                    className=(
                                        "ctrl-label"
                                        " d-inline"
                                    ),
                                ),
                            ],
                            className="mb-2",
                        ),
                        dbc.Row([
                            dbc.Col([
                                _lbl("Stations/survey"),
                                _num(
                                    IDs.INV_AI_N_STATIONS_3D,
                                    20, min=9,
                                    max=200, step=1,
                                ),
                            ]),
                            dbc.Col([
                                _lbl("Corr. length (m)"),
                                _num(
                                    "inv-ai-corr-length",
                                    2000, min=100,
                                    max=20000, step=100,
                                ),
                            ]),
                        ], className="g-1 mb-1"),
                        dbc.Alert(
                            "Adj. radius = 30th-pct"
                            " inter-station dist.",
                            color="dark",
                            className=(
                                "mb-0 py-1 px-2 inv-alert"
                            ),
                        ),
                    ], id=IDs.INV_AI_3D_PANEL,
                       style={"display": "none"}),

                    # Placeholder shown for 1-D
                    html.Div(
                        "Select 2-D or 3-D above"
                        " to see network settings.",
                        id="inv-ai-net-hint",
                        className="inv-sub-hint",
                    ),
                ],
                "ai-acc-net",
            ),

            _acc_item(
                "Training",
                "bi-graph-up-arrow",
                [
                    dbc.Row([
                        dbc.Col([
                            _lbl("N samples"),
                            _num(
                                IDs.INV_AI_N_SAMPLES,
                                1000, min=100,
                                max=20000, step=100,
                            ),
                        ]),
                        dbc.Col([
                            _lbl("N layers"),
                            _num(
                                IDs.INV_AI_N_LAYERS,
                                4, min=2, max=10, step=1,
                            ),
                        ]),
                    ], className="g-1 mb-2"),
                    dbc.Row([
                        dbc.Col([
                            _lbl("Epochs"),
                            _num(
                                IDs.INV_AI_EPOCHS,
                                40, min=5,
                                max=500, step=5,
                            ),
                        ]),
                        dbc.Col([
                            _lbl("Batch size"),
                            _num(
                                IDs.INV_AI_BATCH,
                                64, min=8,
                                max=1024, step=8,
                            ),
                        ]),
                    ], className="g-1 mb-2"),
                    dbc.Row([
                        dbc.Col([
                            _lbl("Learning rate"),
                            _num(
                                IDs.INV_AI_LR,
                                0.001, min=1e-5,
                                max=0.1, step=1e-4,
                            ),
                        ]),
                        dbc.Col([
                            _lbl("Noise level"),
                            _num(
                                IDs.INV_AI_NOISE,
                                0.05, min=0.0,
                                max=0.5, step=0.01,
                            ),
                        ]),
                    ], className="g-1"),
                ],
                "ai-acc-train",
            ),
        ],
        active_item=["ai-acc-setup"],
        always_open=True,
        className="inv-pinn-accordion",
        flush=True,
    )
    return html.Div(
        acc,
        id="inv-ctrl-ai",
        style={"display": "none"},
    )


# ── PINN inversion sidebar section ───────────────────────────────────────────

def _pinn_controls():
    acc = dbc.Accordion(
        [
            _acc_item(
                "Data & Setup",
                "bi-database",
                [
                    _lbl("DL backend"),
                    dbc.Select(
                        id=IDs.INV_BACKEND_SELECT,
                        options=_BACKEND_DL_OPTS,
                        value="auto",
                        size="sm",
                        className="mb-3",
                    ),
                    _lbl("Data source"),
                    dbc.RadioItems(
                        id=IDs.INV_PINN_DATA_SRC,
                        options=[
                            {
                                "label": "Session (EDIs)",
                                "value": "session",
                            },
                            {
                                "label": "Synthetic (fwd)",
                                "value": "synthetic",
                            },
                        ],
                        value="session",
                        inline=False,
                        className="text-muted small",
                    ),
                    html.Div([
                        _lbl("Line"),
                        dcc.Dropdown(
                            id=IDs.INV_PINN_LINE_SEL,
                            placeholder="All lines",
                            multi=True,
                            clearable=True,
                            className=(
                                "mt-1 inv-sta-dropdown"
                            ),
                        ),
                        _lbl("Stations"),
                        dcc.Dropdown(
                            id=IDs.INV_PINN_STATION_SEL,
                            placeholder=(
                                "All stations"
                            ),
                            multi=True,
                            className=(
                                "mt-1 inv-sta-dropdown"
                            ),
                        ),
                    ], className="mt-2"),
                    html.Hr(className="inv-acc-hr"),
                    _lbl("Problem dimension"),
                    dbc.RadioItems(
                        id=IDs.INV_PINN_DIM,
                        options=[
                            {
                                "label": "1-D  sounding",
                                "value": "1d",
                            },
                            {
                                "label": "2-D  joint profile",
                                "value": "2d",
                            },
                            {
                                "label": "3-D  quasi survey",
                                "value": "3d",
                            },
                        ],
                        value="1d",
                        inline=False,
                        className="text-muted small mt-1",
                    ),
                    # 1-D solver (hidden for 2D/3D)
                    html.Div([
                        html.Hr(className="inv-acc-hr"),
                        _lbl("Forward solver"),
                        dbc.Select(
                            id=IDs.INV_PINN_SOLVER,
                            options=[
                                {
                                    "label": "MT 1-D",
                                    "value": "mt1d",
                                },
                                {
                                    "label": "CSAMT 1-D",
                                    "value": "csamt1d",
                                },
                            ],
                            value="mt1d",
                            size="sm",
                        ),
                    ], id="inv-pinn-solver-card"),
                    # 2-D / 3-D mode (hidden for 1D)
                    html.Div([
                        html.Hr(className="inv-acc-hr"),
                        _lbl("EM polarisation"),
                        dbc.RadioItems(
                            id=IDs.INV_PINN_MODE,
                            options=_PINN_MODE_OPTS,
                            value="te",
                            inline=True,
                            className=(
                                "text-muted small mt-1"
                            ),
                        ),
                        dbc.Row([
                            dbc.Col([
                                _lbl("TE comp"),
                                dbc.Select(
                                    id=IDs.INV_PINN_COMP_TE,
                                    options=_COMP_OPTS,
                                    value="xy",
                                    size="sm",
                                ),
                            ]),
                            dbc.Col([
                                _lbl("TM comp"),
                                dbc.Select(
                                    id=IDs.INV_PINN_COMP_TM,
                                    options=_COMP_OPTS,
                                    value="yx",
                                    size="sm",
                                ),
                            ]),
                        ], className="g-1 mt-1"),
                    ], id=IDs.INV_PINN_2D_PANEL,
                       style={"display": "none"}),
                ],
                "pinn-acc-data",
            ),

            _acc_item(
                "Earth Model",
                "bi-layers",
                [
                    dbc.Row([
                        dbc.Col([
                            _lbl("N layers"),
                            _num(
                                IDs.INV_PINN_N_LAYERS,
                                5, min=2, max=30, step=1,
                            ),
                        ]),
                        dbc.Col([
                            _lbl("Depth max (m)"),
                            _num(
                                IDs.INV_PINN_DEPTH_MAX,
                                2000, min=100, step=100,
                            ),
                        ]),
                    ], className="g-1 mb-1"),
                    dbc.Row([
                        dbc.Col([
                            _lbl("N freqs (grid)"),
                            _num(
                                IDs.INV_PINN_N_FREQS,
                                32, min=8, max=128,
                                step=8,
                            ),
                        ]),
                    ], className="g-1"),
                ],
                "pinn-acc-model",
            ),

            _acc_item(
                "Optimiser",
                "bi-lightning-charge",
                [
                    dbc.Row([
                        dbc.Col([
                            _lbl("Epochs"),
                            _num(
                                IDs.INV_PINN_EPOCHS,
                                300, min=10,
                                max=2000, step=10,
                            ),
                        ]),
                        dbc.Col([
                            _lbl("Learning rate"),
                            _num(
                                IDs.INV_PINN_LR,
                                0.01, min=1e-5,
                                max=0.5, step=1e-4,
                            ),
                        ]),
                    ], className="g-1 mb-2"),
                    dbc.Row([
                        dbc.Col([
                            _lbl("Log every N"),
                            _num(
                                IDs.INV_PINN_LOG_EVERY,
                                50, min=1,
                                max=500, step=10,
                            ),
                        ]),
                        dbc.Col([
                            _lbl("Device"),
                            dbc.Select(
                                id=IDs.INV_PINN_DEVICE,
                                options=_DEVICE_OPTS,
                                value="auto",
                                size="sm",
                            ),
                        ]),
                    ], className="g-1"),
                ],
                "pinn-acc-optim",
            ),

            _acc_item(
                "Regularisation",
                "bi-sliders",
                [
                    dbc.Row([
                        dbc.Col([
                            _lbl("Smoothness lam_z"),
                            _num(
                                IDs.INV_PINN_LAM_Z,
                                0.01, min=0,
                                max=1, step=0.001,
                            ),
                        ]),
                    ], className="g-1 mb-1"),
                    html.Div([
                        dbc.Row([
                            dbc.Col([
                                _lbl("Lateral lam_x  (2-D)"),
                                _num(
                                    IDs.INV_PINN_LAM_X,
                                    0.005, min=0,
                                    max=1, step=0.001,
                                ),
                            ]),
                        ], className="g-1"),
                    ], id="inv-pinn-lamx-row",
                       style={"display": "none"}),
                    html.Div([
                        dbc.Row([
                            dbc.Col([
                                _lbl("Graph lam_g  (3-D)"),
                                _num(
                                    IDs.INV_PINN_LAM_G,
                                    0.003, min=0,
                                    max=1, step=0.001,
                                ),
                            ]),
                        ], className="g-1"),
                    ], id="inv-pinn-lamg-row",
                       style={"display": "none"}),
                ],
                "pinn-acc-reg",
            ),
        ],
        active_item=["pinn-acc-data"],
        always_open=True,
        className="inv-pinn-accordion",
        flush=True,
    )

    # 3-D network card: outside accordion so the
    # callback can hide the entire block cleanly.
    net_3d = html.Div([
        html.Div(
            [
                html.I(
                    className="bi bi-diagram-3 me-1",
                    style={"color": "var(--mauve)"},
                ),
                html.Span(
                    "3-D Station Network",
                    className="ctrl-label d-inline",
                ),
            ],
            className="mb-2",
        ),
        dbc.Row([
            dbc.Col([
                _lbl("Adj. radius (m)"),
                _num(
                    IDs.INV_PINN_RADIUS, 5000,
                    min=100, step=100,
                ),
            ]),
            dbc.Col([
                _lbl("Grid spacing (m)"),
                _num(
                    IDs.INV_PINN_SPC, 500,
                    min=10, step=10,
                ),
            ]),
        ], className="g-1 mb-1"),
        _lbl("Coordinate source"),
        dbc.RadioItems(
            id="inv-pinn-coord-src",
            options=[
                {
                    "label": "Auto (EDI lat/lon)",
                    "value": "auto",
                },
                {
                    "label": "Uniform grid",
                    "value": "uniform",
                },
            ],
            value="auto",
            inline=True,
            className="text-muted small mt-1",
        ),
    ], id=IDs.INV_PINN_3D_PANEL,
       className="ctrl-card mt-1",
       style={"display": "none"})

    return html.Div(
        [acc, net_3d],
        id="inv-ctrl-pinn",
        style={"display": "none"},
    )


# ── Hybrid inversion sidebar section ─────────────────────────────────────────

def _hybrid_controls():
    acc = dbc.Accordion(
        [
            _acc_item(
                "Data & Dimension",
                "bi-database",
                [
                    _lbl("Data source"),
                    dbc.RadioItems(
                        id=IDs.INV_HYB_DATA_SRC,
                        options=[
                            {
                                "label": "Session (EDIs)",
                                "value": "session",
                            },
                            {
                                "label": "Synthetic (fwd)",
                                "value": "synthetic",
                            },
                        ],
                        value="session",
                        inline=False,
                        className="text-muted small",
                    ),
                    html.Div([
                        _lbl("Line"),
                        dcc.Dropdown(
                            id=IDs.INV_HYB_LINE_SEL,
                            placeholder="All lines",
                            multi=True,
                            clearable=True,
                            className=(
                                "mt-1 inv-sta-dropdown"
                            ),
                        ),
                        _lbl("Stations"),
                        dcc.Dropdown(
                            id=IDs.INV_HYB_STATION_SEL,
                            placeholder=(
                                "All stations"
                            ),
                            multi=True,
                            className=(
                                "mt-1 inv-sta-dropdown"
                            ),
                        ),
                    ], className="mt-2"),
                    html.Hr(className="inv-acc-hr"),
                    _lbl("Problem dimension"),
                    dbc.RadioItems(
                        id=IDs.INV_HYB_DIM,
                        options=[
                            {
                                "label": "1-D  sounding",
                                "value": "1d",
                            },
                            {
                                "label": "2-D  profile",
                                "value": "2d",
                            },
                            {
                                "label": "3-D  survey",
                                "value": "3d",
                            },
                        ],
                        value="1d",
                        inline=False,
                        className=(
                            "text-muted small mt-1"
                        ),
                    ),
                ],
                "hyb-acc-data",
            ),

            _acc_item(
                "Stage 1 — AI Pre-model",
                "bi-cpu",
                [
                    html.Div(
                        "Upload a saved EMInverter"
                        " checkpoint (.npz):",
                        className="inv-sub-hint mb-2",
                    ),
                    dcc.Upload(
                        id=IDs.INV_HYB_CHECKPOINT,
                        children=html.Div([
                            html.I(
                                className=(
                                    "bi bi-file-earmark"
                                    "-zip me-1"
                                ),
                            ),
                            "Drop .npz / click",
                        ]),
                        accept=".npz",
                        className=(
                            "hyb-checkpoint-upload"
                        ),
                    ),
                    html.Div(
                        id=IDs.INV_HYB_MODEL_INFO,
                        className="inv-sub-hint mt-1",
                        children="No checkpoint loaded.",
                    ),
                    dbc.Button(
                        [
                            html.I(
                                className="bi bi-eye me-1"
                            ),
                            "Preview Stage-1",
                        ],
                        id=IDs.BTN_HYB_STAGE1_PREV,
                        color="secondary",
                        size="sm",
                        outline=True,
                        n_clicks=0,
                        className="mt-2",
                    ),
                ],
                "hyb-acc-stage1",
            ),

            _acc_item(
                "Stage 2 — PINN Refinement",
                "bi-bezier2",
                [
                    dbc.Row([
                        dbc.Col([
                            _lbl("Epochs"),
                            _num(
                                IDs.INV_HYB_EPOCHS,
                                150, min=10,
                                max=2000, step=10,
                            ),
                        ]),
                        dbc.Col([
                            _lbl("Learning rate"),
                            _num(
                                IDs.INV_HYB_LR,
                                0.005, min=1e-5,
                                max=0.5, step=1e-4,
                            ),
                        ]),
                    ], className="g-1 mb-2"),
                    dbc.Row([
                        dbc.Col([
                            _lbl("Device"),
                            dbc.Select(
                                id=IDs.INV_HYB_DEVICE,
                                options=_DEVICE_OPTS,
                                value="auto",
                                size="sm",
                            ),
                        ]),
                        dbc.Col([
                            _lbl("Smoothness lam_z"),
                            _num(
                                IDs.INV_HYB_LAM_Z,
                                0.005, min=0,
                                max=1, step=0.001,
                            ),
                        ]),
                    ], className="g-1 mb-1"),
                    html.Div([
                        dbc.Row([
                            dbc.Col([
                                _lbl(
                                    "Lateral lam_x  (2-D)"
                                ),
                                _num(
                                    IDs.INV_HYB_LAM_X,
                                    0.003, min=0,
                                    max=1, step=0.001,
                                ),
                            ]),
                        ], className="g-1"),
                    ], id="inv-hyb-lamx-row",
                       style={"display": "none"}),
                    html.Div([
                        dbc.Row([
                            dbc.Col([
                                _lbl(
                                    "Graph lam_g  (3-D)"
                                ),
                                _num(
                                    IDs.INV_HYB_LAM_G,
                                    0.003, min=0,
                                    max=1, step=0.001,
                                ),
                            ]),
                        ], className="g-1"),
                    ], id="inv-hyb-lamg-row",
                       style={"display": "none"}),
                ],
                "hyb-acc-stage2",
            ),

            _acc_item(
                "Settings",
                "bi-gear",
                [
                    dbc.Row([
                        dbc.Col([
                            _lbl("N layers"),
                            _num(
                                IDs.INV_HYB_N_LAYERS,
                                5, min=2,
                                max=30, step=1,
                            ),
                        ]),
                        dbc.Col([
                            _lbl("N freqs"),
                            _num(
                                IDs.INV_HYB_N_FREQS,
                                32, min=8,
                                max=128, step=8,
                            ),
                        ]),
                    ], className="g-1 mb-2"),
                    dbc.Row([
                        dbc.Col([
                            _lbl("EM mode"),
                            dbc.Select(
                                id=IDs.INV_HYB_MODE,
                                options=_PINN_MODE_OPTS,
                                value="te",
                                size="sm",
                            ),
                        ]),
                    ], className="g-1 mb-1"),
                    dbc.Row([
                        dbc.Col([
                            _lbl("TE comp"),
                            dbc.Select(
                                id=IDs.INV_HYB_COMP_TE,
                                options=_COMP_OPTS,
                                value="xy",
                                size="sm",
                            ),
                        ]),
                        dbc.Col([
                            _lbl("TM comp"),
                            dbc.Select(
                                id=IDs.INV_HYB_COMP_TM,
                                options=_COMP_OPTS,
                                value="yx",
                                size="sm",
                            ),
                        ]),
                    ], className="g-1"),
                ],
                "hyb-acc-settings",
            ),
        ],
        active_item=["hyb-acc-data"],
        always_open=True,
        className="inv-pinn-accordion",
        flush=True,
    )

    cond_panels = html.Div([
        # 2-D hint card (shown for 2D dim)
        html.Div(
            html.Div(
                "2-D: lateral smoothness is active.",
                className="inv-sub-hint",
            ),
            id=IDs.INV_HYB_2D_PANEL,
            className="ctrl-card",
            style={"display": "none"},
        ),
        # 3-D network card (shown for 3D dim)
        html.Div([
            html.Div(
                [
                    html.I(
                        className=(
                            "bi bi-diagram-3 me-1"
                        ),
                        style={
                            "color": "var(--mauve)"
                        },
                    ),
                    html.Span(
                        "3-D Network",
                        className=(
                            "ctrl-label d-inline"
                        ),
                    ),
                ],
                className="mb-2",
            ),
            dbc.Row([
                dbc.Col([
                    _lbl("Adj. radius (m)"),
                    _num(
                        IDs.INV_HYB_RADIUS, 5000,
                        min=100, step=100,
                    ),
                ]),
            ], className="g-1"),
        ], id=IDs.INV_HYB_3D_PANEL,
           className="ctrl-card mt-1",
           style={"display": "none"}),
    ])

    return html.Div(
        [acc, cond_panels],
        id="inv-ctrl-hybrid",
        style={"display": "none"},
    )


# ── PyTorch install modal (overlay) ───────────────────────────────────────────

def _install_modal() -> dbc.Modal:
    torch_badge = (
        dbc.Badge("Installed", color="success", className="ms-2")
        if _TORCH_AVAIL else
        dbc.Badge("Not installed", color="warning", className="ms-2")
    )
    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle([
                    html.I(className="bi bi-cpu me-2", style={"color": "#a6e3a1"}),
                    "PyTorch", torch_badge,
                ]),
                close_button=False,
            ),
            dbc.ModalBody([
                dbc.Alert(
                    [
                        html.I(className="bi bi-exclamation-triangle-fill me-2"),
                        "PyTorch is required for AI inversion but is not installed "
                        "in the current Python environment.",
                    ],
                    color="warning", className="mb-3",
                    style={"fontSize": "12px"},
                ),
                html.Div([
                    dbc.Row([
                        dbc.Col([
                            html.Div("Package", style={"fontSize": "10px", "color": "#6c7086"}),
                            html.Div("torch (CPU build)",
                                     style={"fontSize": "12px", "color": "#cdd6f4"}),
                        ]),
                        dbc.Col([
                            html.Div("Size", style={"fontSize": "10px", "color": "#6c7086"}),
                            html.Div("~750 MB",
                                     style={"fontSize": "12px", "color": "#cdd6f4"}),
                        ]),
                        dbc.Col([
                            html.Div("Est. time", style={"fontSize": "10px", "color": "#6c7086"}),
                            html.Div("2–5 min",
                                     style={"fontSize": "12px", "color": "#cdd6f4"}),
                        ]),
                    ], className="g-2 mb-3"),
                    html.Div(
                        "The CPU build is installed automatically into the current Python "
                        "environment. For GPU (CUDA), install manually with torch+cuXXX.",
                        style={"fontSize": "10px", "color": "#6c7086"},
                    ),
                ], style={"padding": "8px 12px", "background": "#1e1e2e",
                          "borderRadius": "6px", "marginBottom": "12px"}),

                # Progress area (hidden until install starts)
                html.Div([
                    html.Div(id=IDs.INV_INSTALL_LABEL,
                             style={"fontSize": "11px", "color": "#a6adc8",
                                    "marginBottom": "6px"}),
                    dbc.Progress(
                        id=IDs.INV_INSTALL_PROGRESS, value=0,
                        animated=True, striped=True, color="info",
                        style={"height": "20px", "marginBottom": "10px",
                               "borderRadius": "6px"},
                    ),
                    html.Pre(
                        id=IDs.INV_INSTALL_LOG, className="log-area",
                        style={"height": "180px", "margin": "0",
                               "fontSize": "10px", "lineHeight": "1.4"},
                    ),
                ], id=IDs.INV_INSTALL_PROG_AREA, style={"display": "none"}),

                dcc.Interval(
                    id=IDs.INV_INSTALL_INTERVAL,
                    interval=800, n_intervals=0, disabled=True,
                ),
            ]),
            dbc.ModalFooter([
                dbc.Button(
                    [html.I(className="bi bi-x-circle me-1"), "Cancel"],
                    id=IDs.BTN_INV_INSTALL_NO, color="secondary", size="sm",
                    className="me-auto",
                ),
                dbc.Button(
                    [html.I(className="bi bi-download me-1"), "Install PyTorch"],
                    id=IDs.BTN_INV_INSTALL_YES, color="primary", size="sm",
                ),
            ], id="inv-install-footer"),
        ],
        id=IDs.INV_INSTALL_MODAL,
        is_open=False, size="lg",
        backdrop="static", keyboard=False,
    )


# ── Page layout ───────────────────────────────────────────────────────────────

def layout() -> html.Div:

    # ── Left sidebar ──────────────────────────────────────────────────────────

    stores = html.Div([
        dcc.Store(id=IDs.INV_ACTIVE_TRACK,
                  data=_DEFAULT_TRACK),
        dcc.Store(id=IDs.INV_ACTIVE_OUT,
                  data=_DEFAULT_OUT),
        dcc.Store(id=IDs.INV_PINN_RESULT_STORE),
        dcc.Store(id=IDs.INV_HYB_RESULT_STORE),
    ])

    # Sticky run bar
    run_bar = html.Div([
        dbc.Button(
            [
                html.I(
                    className="bi bi-play-fill me-1"
                ),
                "Run Inversion",
            ],
            id=IDs.BTN_INV_RUN,
            color="primary",
            size="sm",
            className="w-100 mb-1",
            n_clicks=0,
        ),
        html.Div(
            [
                dbc.Spinner(
                    html.Div(id=IDs.INV_SPINNER),
                    size="sm",
                    color="primary",
                    spinnerClassName=(
                        "inv-run-spinner"
                    ),
                ),
                html.Div(
                    id=IDs.INV_RUN_MSG,
                    className="inv-run-msg",
                ),
            ],
            className="inv-spinner-row",
        ),
        html.Div(
            id=IDs.INV_FEEDBACK,
            className="fwd-feedback-mini",
        ),
    ], className="fwd-run-bar inv-run-bar")

    def _track_btn(tid, label, icon):
        active = " active" if tid == _DEFAULT_TRACK else ""
        return html.Button(
            [
                html.I(className=f"bi {icon} me-1"),
                label,
            ],
            id=f"inv-track-{tid}",
            className=f"fwd-dim-btn{active}",
            n_clicks=0,
        )

    track_bar = html.Div([
        html.Div(
            "Classical",
            className="inv-track-group-label",
        ),
        html.Div(
            [
                _track_btn(
                    "trad", "Traditional", "bi-cpu"
                ),
                _track_btn(
                    "ai", "AI Neural", "bi-robot"
                ),
            ],
            className="fwd-dim-bar",
        ),
        html.Div(
            "Physics-informed",
            className="inv-track-group-label mt-2",
        ),
        html.Div(
            [
                _track_btn(
                    "pinn", "PINN", "bi-bezier2"
                ),
                _track_btn(
                    "hybrid", "Hybrid", "bi-node-plus"
                ),
            ],
            className="fwd-dim-bar",
        ),
    ])

    # Shared: frequency range
    freq_card = html.Div([
        _clabel("Frequency range [log₁₀ Hz]"),
        dbc.Row([
            dbc.Col([_lbl("Min"), _num(IDs.INV_FREQ_MIN, -3, step=0.5)]),
            dbc.Col([_lbl("Max"), _num(IDs.INV_FREQ_MAX,  4, step=0.5)]),
            dbc.Col([_lbl("N pts"), _num(IDs.INV_N_FREQ, 20, min=5, max=80, step=5)]),
        ], className="g-1"),
    ], className="ctrl-card")

    ctrl_scroll = html.Div([
        track_bar,
        _trad_controls(),
        _ai_controls(),
        _pinn_controls(),
        _hybrid_controls(),
        freq_card,
    ], className="fwd-ctrl-scroll")

    sidebar = html.Div(
        [stores, run_bar, ctrl_scroll],
        className="analysis-controls fwd-sidebar",
    )

    # ── View area (right side) ────────────────────────────────────────────────

    out_tab_bar = html.Div([
        html.Button(
            [html.I(className=f"bi {icon} me-1"), label],
            id=f"inv-out-tab-{tid}",
            className=f"prof-tab-btn{' active' if tid == _DEFAULT_OUT else ''}",
            n_clicks=0,
        )
        for tid, label, icon in _INV_OUT_TABS
    ], className="prof-tab-bar")

    ctx_bar = html.Div(
        [
            html.Span(
                [html.I(className="bi bi-cpu me-1"), "Traditional"],
                className="adv-info-group",
            ),
            html.Span("·", className="adv-info-sep"),
            html.Span(
                "Configure the model then click Run Inversion",
                style={"color": "var(--sub0)", "fontSize": "10.5px"},
            ),
        ],
        id=IDs.INV_CTX_BAR,
        className="adv-info-bar",
    )

    # Per-tab output panels
    def _fig_panel(panel_id: str, img_id: str, visible: bool = False) -> html.Div:
        return html.Div(
            html.Div(
                html.Img(id=img_id, src=empty_src(), style=_IMG_STYLE),
                className="fig-img-wrap profile-page-fig",
            ),
            id=panel_id,
            className="prof-panel",
            style={"display": "flex" if visible else "none"},
        )

    panel_result = _fig_panel("inv-panel-result", IDs.IMG_INV, visible=True)

    panel_conv = _fig_panel("inv-panel-conv", IDs.IMG_INV_CONV)

    panel_stats = html.Div(
        html.Div(
            html.Div(
                "Run the inversion to see statistics here.",
                className="inv-stats-placeholder",
            ),
            id=IDs.INV_STATS,
            className="inv-stats-body",
        ),
        id="inv-panel-stats",
        className="prof-panel",
        style={"display": "none"},
    )

    panel_log = html.Div(
        html.Pre(
            id=IDs.INV_LOG,
            className="log-area inv-log-area",
        ),
        id="inv-panel-log",
        className="prof-panel",
        style={"display": "none"},
    )

    panel_data_fit = html.Div(
        html.Div(
            [
                html.Div(
                    "Run PINN or Hybrid inversion,"
                    " then select a station to view"
                    " observed vs predicted curves.",
                    className="inv-stats-placeholder",
                    id="inv-datafit-hint",
                ),
                html.Div(
                    dbc.Row([
                        dbc.Col([
                            html.Div([
                                html.I(
                                    className=(
                                        "bi bi-reception-4"
                                        " me-1"
                                    ),
                                    style={
                                        "color":
                                        "var(--mauve)"
                                    },
                                ),
                                _lbl("Line"),
                            ], className=(
                                "d-flex"
                                " align-items-center"
                                " gap-1 mb-1"
                            )),
                            dcc.Dropdown(
                                id=IDs.INV_DATAFIT_LINE,
                                placeholder="All lines",
                                clearable=True,
                                className=(
                                    "inv-sta-dropdown"
                                ),
                            ),
                        ], width=5),
                        dbc.Col([
                            html.Div([
                                html.I(
                                    className=(
                                        "bi bi-geo-alt"
                                        " me-1"
                                    ),
                                    style={
                                        "color":
                                        "var(--blue)"
                                    },
                                ),
                                _lbl("Station"),
                            ], className=(
                                "d-flex"
                                " align-items-center"
                                " gap-1 mb-1"
                            )),
                            dcc.Dropdown(
                                id=IDs.INV_DATAFIT_STA,
                                placeholder=(
                                    "Pick station..."
                                ),
                                className=(
                                    "inv-sta-dropdown"
                                ),
                            ),
                        ], width=7),
                    ], className="g-2 mt-2"),
                    id="inv-datafit-selectors",
                    style={"display": "none"},
                ),
                html.Img(
                    id=IDs.IMG_INV_DATA_FIT,
                    src=empty_src(),
                    style=_IMG_STYLE,
                    className="mt-2",
                ),
            ],
        ),
        id="inv-panel-data-fit",
        className="prof-panel",
        style={"display": "none"},
    )

    panels = html.Div(
        [
            panel_result, panel_conv,
            panel_stats, panel_log,
            panel_data_fit,
        ],
        className="prof-view-panel",
    )

    view_area = html.Div(
        [out_tab_bar, ctx_bar, panels],
        className="adv-view-area",
    )

    # ── Final assembly ────────────────────────────────────────────────────────

    return html.Div([
        _command_bar(
            "inversion", "MT Inversion",
            "Traditional · AI Neural · PINN · Hybrid"
            "  |  1D / 2D / 3D",
        ),
        html.Div(
            [sidebar, view_area],
            className="analysis-layout",
            style={"flex": "1"},
        ),
        _install_modal(),
    ], style={"display": "flex", "flexDirection": "column", "height": "100%"})


def register_callbacks(app) -> None:
    pass  # handled by callbacks/inversion.py
