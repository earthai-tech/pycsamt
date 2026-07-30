# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Inversion Results Viewer page.

Supports three inversion codes:
  • ModEM    — 3-D MT (Modular_NLCG_*.rho, *.log, *.cov)
  • Occam2D  — 2-D MT (STARTUP, ITER*, *.resp files)
  • MARE2DEM — 2.5-D MT/CSEM (*.emdata, *.resistivity)

Layout (mirrors 3D Map page structure)
---------------------------------------
Left panel (260 px, analysis-layout grid):
  ┌ sticky top ─────────────────────────────────────────────────┐
  │  [folder icon path input] [Browse…] [Load]                  │
  │  Loading spinner (pycsamt colours) + status line            │
  └─────────────────────────────────────────────────────────────┘
  ┌ fwd-ctrl-scroll (scrollable, overflow-y: auto) ─────────────┐
  │  Accordion                                                   │
  │    ▶ Solver        — radio (auto/ModEM/Occam2D/MARE2DEM)    │
  │    ▶ Display       — model sel · ρ range · depth · cmap     │
  │    ▶ ModEM extras  — iteration# · section/depthmap/profiles │
  │                      (hidden when solver ≠ modem)           │
  │    ▶ Occam2D extras— mode · (hidden when solver ≠ occam2d) │
  │    ▶ View Controls — context panel per active tab           │
  │    ▶ Stations      — markers · names · tolerance            │
  └─────────────────────────────────────────────────────────────┘
  ┌ sticky bottom ──────────────────────────────────────────────┐
  │  [Generate Plot]  [Export PNG]                              │
  └─────────────────────────────────────────────────────────────┘
"""

from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html

from pycsamt.app.web.layout import IDs, _command_bar
from pycsamt.app.web.utils import empty_src

PAGE_ID = "inv-results"

# ── Option lists ──────────────────────────────────────────────────────────────

_SOLVER_OPTS = [
    {"label": "Auto-detect", "value": "auto"},
    {"label": "ModEM (3-D)", "value": "modem"},
    {"label": "Occam2D", "value": "occam2d"},
    {"label": "MARE2DEM", "value": "mare2dem"},
]

_WHICH_OPTS = [
    {"label": "Final model", "value": "final"},
    {"label": "Initial model", "value": "initial"},
    {"label": "Iteration …", "value": "iteration"},
]

_DIR_OPTS = [
    {"label": "N-S  (fixed easting)", "value": "NS"},
    {"label": "E-W  (fixed northing)", "value": "EW"},
]

_COMP_OPTS = [
    {"label": "ZXX", "value": "ZXX"},
    {"label": "ZXY", "value": "ZXY"},
    {"label": "ZYX", "value": "ZYX"},
    {"label": "ZYY", "value": "ZYY"},
]

_CMAP_OPTS = [
    {"label": "jet_r", "value": "jet_r"},
    {"label": "viridis", "value": "viridis"},
    {"label": "plasma", "value": "plasma"},
    {"label": "RdBu_r", "value": "RdBu_r"},
    {"label": "coolwarm", "value": "coolwarm"},
    {"label": "bwr", "value": "bwr"},
    {"label": "seismic", "value": "seismic"},
]

_COV_CMAP_OPTS = [
    {"label": "Blues", "value": "Blues"},
    {"label": "Greens", "value": "Greens"},
    {"label": "Oranges", "value": "Oranges"},
    {"label": "Greys", "value": "Greys"},
]

# (slug, label, icon, solvers)
_ALL_TABS = [
    (
        "conv",
        "Convergence",
        "bi-graph-down-arrow",
        {"modem", "occam2d", "mare2dem"},
    ),
    ("section", "Section", "bi-layout-split", {"modem"}),
    ("depthmap", "Depth Map", "bi-layers", {"modem"}),
    ("profiles", "All Profiles", "bi-list-columns-reverse", {"modem"}),
    ("covariance", "Covariance", "bi-grid-3x3", {"modem"}),
    ("model", "Model", "bi-square-half", {"occam2d", "mare2dem"}),
    ("response", "Response", "bi-activity", {"modem", "occam2d", "mare2dem"}),
    ("pseudo", "Pseudo", "bi-bar-chart-steps", {"modem", "occam2d"}),
    ("sounding1d", "1D Sounding", "bi-reception-4", {"occam2d"}),
    ("sitemisfit", "Site Misfit", "bi-exclamation-diamond", {"occam2d"}),
    ("respgrid", "Resp. Grid", "bi-grid", {"occam2d"}),
    ("survey", "Survey Layout", "bi-geo-alt", {"mare2dem"}),
]

_DEFAULT_TAB = "conv"
_IMG_STYLE = {"width": "100%", "height": "100%", "objectFit": "contain"}


# ── Small helpers ─────────────────────────────────────────────────────────────


def _lbl(text: str) -> html.Small:
    return html.Small(text, className="text-muted")


def _clabel(text: str) -> html.Div:
    return html.Div(text, className="ctrl-label")


def _num(id_: str, val, **kw) -> dbc.Input:
    return dbc.Input(
        id=id_,
        type="number",
        value=val,
        size="sm",
        step="any",
        debounce=True,
        **kw,
    )


def _acc_title(icon: str, text: str) -> html.Span:
    return html.Span(
        [
            html.I(
                className=f"bi {icon} me-2",
                style={"color": "var(--mauve)", "fontSize": "13px"},
            ),
            text,
        ],
        className="map3d-settings-title",
    )


def _acc_item(
    item_id: str, icon: str, title: str, *children, start_open: bool = False
) -> dbc.AccordionItem:
    return dbc.AccordionItem(
        list(children),
        title=_acc_title(icon, title),
        item_id=item_id,
    )


# ── Per-tab context panels (shown/hidden by active tab) ───────────────────────


def _ctx_conv() -> html.Div:
    return html.Div(
        html.Div(
            "Convergence curve — no extra controls needed.",
            className="inv-sub-hint",
        ),
        id="invr-ctx-conv",
    )


def _ctx_section() -> html.Div:
    return html.Div(
        [
            _clabel("Profile Mode"),
            dbc.RadioItems(
                id=IDs.INVR_SECT_MODE,
                options=[
                    {"label": "Axis-aligned (NS / EW)", "value": "axial"},
                    {"label": "Arbitrary azimuth", "value": "arbitrary"},
                ],
                value="axial",
                inline=False,
                className="text-muted small",
            ),
            html.Div(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    _lbl("Direction"),
                                    dbc.Select(
                                        id=IDs.INVR_SECT_DIR,
                                        options=_DIR_OPTS,
                                        value="NS",
                                        size="sm",
                                    ),
                                ]
                            ),
                            dbc.Col(
                                [
                                    _lbl("Offset (m)"),
                                    _num(IDs.INVR_SECT_OFFSET, 0),
                                ]
                            ),
                        ],
                        className="g-1 mt-1",
                    ),
                ],
                id=IDs.INVR_SECT_AX_PANEL,
            ),
            html.Div(
                [
                    _clabel("Start point"),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    _lbl("Lat °N"),
                                    _num(
                                        IDs.INVR_SECT_START_LAT,
                                        None,
                                        placeholder="32.05",
                                    ),
                                ]
                            ),
                            dbc.Col(
                                [
                                    _lbl("Lon °E"),
                                    _num(
                                        IDs.INVR_SECT_START_LON,
                                        None,
                                        placeholder="119.08",
                                    ),
                                ]
                            ),
                        ],
                        className="g-1 mt-1",
                    ),
                    _clabel("End point"),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    _lbl("Lat °N"),
                                    _num(
                                        IDs.INVR_SECT_END_LAT,
                                        None,
                                        placeholder="32.22",
                                    ),
                                ]
                            ),
                            dbc.Col(
                                [
                                    _lbl("Lon °E"),
                                    _num(
                                        IDs.INVR_SECT_END_LON,
                                        None,
                                        placeholder="119.19",
                                    ),
                                ]
                            ),
                        ],
                        className="g-1 mt-1",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    _lbl("N samples"),
                                    _num(
                                        IDs.INVR_SECT_N_SAMP,
                                        200,
                                        min=20,
                                        max=1000,
                                    ),
                                ]
                            ),
                        ],
                        className="g-1 mt-1",
                    ),
                    _clabel("Geo origin (opt.)"),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    _lbl("Lat °N"),
                                    _num(
                                        IDs.INVR_SECT_ORIG_LAT,
                                        None,
                                        placeholder="32.129",
                                    ),
                                ]
                            ),
                            dbc.Col(
                                [
                                    _lbl("Lon °E"),
                                    _num(
                                        IDs.INVR_SECT_ORIG_LON,
                                        None,
                                        placeholder="119.125",
                                    ),
                                ]
                            ),
                        ],
                        className="g-1 mt-1",
                    ),
                ],
                id=IDs.INVR_SECT_ARB_PANEL,
                style={"display": "none"},
            ),
        ],
        id="invr-ctx-section",
        style={"display": "none"},
    )


def _ctx_depthmap() -> html.Div:
    return html.Div(
        [
            _clabel("Depths (m)"),
            dbc.Input(
                id=IDs.INVR_DM_DEPTHS,
                placeholder="50, 200, 500, 1000",
                value="50, 200, 500, 1000",
                type="text",
                size="sm",
                debounce=True,
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            _lbl("Grid cols"),
                            _num(IDs.INVR_DM_NCOLS, 2, min=1, max=4),
                        ]
                    ),
                ],
                className="g-1 mt-1",
            ),
            _clabel("Geo origin (opt.)"),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            _lbl("Lat °N"),
                            _num(IDs.INVR_DM_LAT, None, placeholder="32.129"),
                        ]
                    ),
                    dbc.Col(
                        [
                            _lbl("Lon °E"),
                            _num(IDs.INVR_DM_LON, None, placeholder="119.125"),
                        ]
                    ),
                ],
                className="g-1 mt-1",
            ),
        ],
        id="invr-ctx-depthmap",
        style={"display": "none"},
    )


def _ctx_profiles() -> html.Div:
    return html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            _lbl("Direction"),
                            dbc.Select(
                                id=IDs.INVR_AP_DIR,
                                options=_DIR_OPTS,
                                value="NS",
                                size="sm",
                            ),
                        ]
                    ),
                    dbc.Col(
                        [
                            _lbl("N profiles"),
                            _num(IDs.INVR_AP_N, 5, min=1, max=20),
                        ]
                    ),
                ],
                className="g-1",
            ),
            _clabel("Offsets (m, comma-sep · leave blank = auto)"),
            dbc.Input(
                id=IDs.INVR_AP_OFFSETS,
                placeholder="353, 155, -41, -242, -445",
                type="text",
                size="sm",
                debounce=True,
            ),
            _clabel("Profile names (opt.)"),
            dbc.Input(
                id=IDs.INVR_AP_NAMES,
                placeholder="L18, L22, L26, L30, L34",
                type="text",
                size="sm",
                debounce=True,
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            _lbl("Grid cols"),
                            _num(IDs.INVR_AP_NCOLS, 3, min=1, max=5),
                        ]
                    ),
                ],
                className="g-1 mt-1",
            ),
        ],
        id="invr-ctx-profiles",
        style={"display": "none"},
    )


def _ctx_covariance() -> html.Div:
    return html.Div(
        [
            dbc.Switch(
                id=IDs.INVR_COV_SMOOTH,
                label="Show smoothing panel",
                value=True,
                className="small",
            ),
            _lbl("Mask colourmap"),
            dbc.Select(
                id=IDs.INVR_COV_CMAP,
                options=_COV_CMAP_OPTS,
                value="Blues",
                size="sm",
                className="mt-1",
            ),
        ],
        id="invr-ctx-covariance",
        style={"display": "none"},
    )


def _ctx_model() -> html.Div:
    return html.Div(
        [
            html.Div(
                "Occam2D / MARE2DEM 2-D resistivity section.",
                className="inv-sub-hint",
            ),
        ],
        id="invr-ctx-model",
        style={"display": "none"},
    )


def _ctx_response() -> html.Div:
    return html.Div(
        [
            _clabel("Station"),
            dbc.Select(
                id=IDs.INVR_RESP_STATION,
                options=[],
                value=None,
                placeholder="Load a result first…",
                size="sm",
            ),
            _clabel("Component"),
            dbc.Select(
                id=IDs.INVR_RESP_COMP,
                options=_COMP_OPTS,
                value="ZXY",
                size="sm",
                className="mt-1",
            ),
        ],
        id="invr-ctx-response",
        style={"display": "none"},
    )


def _ctx_pseudo() -> html.Div:
    return html.Div(
        [
            _clabel("Component"),
            dbc.Select(
                id=IDs.INVR_PS_COMP,
                options=_COMP_OPTS,
                value="ZXY",
                size="sm",
            ),
        ],
        id="invr-ctx-pseudo",
        style={"display": "none"},
    )


def _ctx_sounding1d() -> html.Div:
    return html.Div(
        [
            _clabel("Station"),
            dbc.Select(
                id="invr-s1d-station",
                options=[],
                value=None,
                placeholder="Load an Occam2D result first…",
                size="sm",
            ),
        ],
        id="invr-ctx-sounding1d",
        style={"display": "none"},
    )


def _ctx_sitemisfit() -> html.Div:
    return html.Div(
        [
            html.Div(
                "Per-station RMS misfit map (Occam2D).",
                className="inv-sub-hint",
            ),
        ],
        id="invr-ctx-sitemisfit",
        style={"display": "none"},
    )


def _ctx_respgrid() -> html.Div:
    return html.Div(
        [
            _clabel("Component"),
            dbc.Select(
                id="invr-rg-comp", options=_COMP_OPTS, value="ZXY", size="sm"
            ),
        ],
        id="invr-ctx-respgrid",
        style={"display": "none"},
    )


def _ctx_survey() -> html.Div:
    return html.Div(
        [
            html.Div(
                "Transmitter + receiver geometry (MARE2DEM).",
                className="inv-sub-hint",
            ),
        ],
        id="invr-ctx-survey",
        style={"display": "none"},
    )


_CTX_BUILDERS = {
    "conv": _ctx_conv,
    "section": _ctx_section,
    "depthmap": _ctx_depthmap,
    "profiles": _ctx_profiles,
    "covariance": _ctx_covariance,
    "model": _ctx_model,
    "response": _ctx_response,
    "pseudo": _ctx_pseudo,
    "sounding1d": _ctx_sounding1d,
    "sitemisfit": _ctx_sitemisfit,
    "respgrid": _ctx_respgrid,
    "survey": _ctx_survey,
}


# ── Accordion content builders (one per logical section) ─────────────────────


def _acc_solver() -> list:
    """Solver selector + status strip."""
    return [
        dbc.RadioItems(
            id=IDs.INVR_SOLVER_RADIO,
            options=_SOLVER_OPTS,
            value="auto",
            inline=False,
            className="text-muted small",
        ),
        # Status line (updates after Load)
        html.Div(
            [
                html.I(className="bi bi-hourglass me-1"),
                "No results loaded yet",
            ],
            id=IDs.INVR_STATUS,
            className="invr-status-line mt-2",
        ),
    ]


def _acc_display() -> list:
    """Display settings shared across all solvers."""
    return [
        _clabel("ρ range (Ω·m)"),
        dbc.Row(
            [
                dbc.Col([_lbl("min"), _num(IDs.INVR_RHO_MIN, 1, min=0.001)]),
                dbc.Col([_lbl("max"), _num(IDs.INVR_RHO_MAX, 1000, min=0.01)]),
            ],
            className="g-1",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        _lbl("Depth max (m)"),
                        _num(IDs.INVR_DEPTH_MAX, 2000, min=10),
                    ]
                ),
            ],
            className="g-1 mt-1",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        _lbl("Colourmap"),
                        dbc.Select(
                            id=IDs.INVR_CMAP,
                            options=_CMAP_OPTS,
                            value="jet_r",
                            size="sm",
                        ),
                    ]
                ),
            ],
            className="g-1 mt-1",
        ),
    ]


def _acc_modem_extras() -> list:
    """ModEM-specific display settings (hidden until modem loaded)."""
    return [
        _clabel("Model snapshot"),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Select(
                            id=IDs.INVR_WHICH,
                            options=_WHICH_OPTS,
                            value="final",
                            size="sm",
                        ),
                    ],
                    width=7,
                ),
                dbc.Col(
                    [
                        _num(
                            IDs.INVR_ITERATION,
                            None,
                            min=0,
                            placeholder="iter #",
                        ),
                    ],
                    width=5,
                ),
            ],
            className="g-1",
        ),
    ]


def _acc_occam2d_extras() -> list:
    """Occam2D-specific display settings."""
    return [
        html.Div(
            [
                html.I(className="bi bi-info-circle me-1 text-muted"),
                "Model uses the final ITER file automatically.",
            ],
            className="inv-sub-hint",
        ),
        _clabel("Modes"),
        dbc.Checklist(
            id="invr-occ-modes",
            options=[
                {"label": "TE", "value": "TE"},
                {"label": "TM", "value": "TM"},
            ],
            value=["TE", "TM"],
            inline=True,
            className="text-muted small mt-1",
        ),
    ]


def _acc_mare2dem_extras() -> list:
    """MARE2DEM-specific display settings."""
    return [
        html.Div(
            [
                html.I(className="bi bi-info-circle me-1 text-muted"),
                "Loads resistivity from final *.resistivity file.",
            ],
            className="inv-sub-hint",
        ),
    ]


def _acc_view_controls() -> list:
    """Per-tab view controls (one hidden panel per tab)."""
    return [fn() for fn in _CTX_BUILDERS.values()]


def _acc_stations() -> list:
    """Station marker controls."""
    return [
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Switch(
                            id=IDs.INVR_SHOW_STA,
                            label="Show markers",
                            value=True,
                            className="small",
                        )
                    ]
                ),
                dbc.Col(
                    [
                        dbc.Switch(
                            id=IDs.INVR_SHOW_NAMES,
                            label="Show names",
                            value=True,
                            className="small",
                        )
                    ]
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        _lbl("Tolerance (m)"),
                        _num(IDs.INVR_STA_TOL, 300, min=0),
                    ]
                ),
            ],
            className="g-1 mt-1",
        ),
    ]


# ── Sidebar (full assembly) ───────────────────────────────────────────────────


def _load_bar() -> html.Div:
    """Sticky top bar: path input + browse button + load button + spinner."""
    return html.Div(
        [
            _clabel("Results folder"),
            dbc.InputGroup(
                [
                    dbc.Input(
                        id=IDs.INVR_PATH_INPUT,
                        placeholder="Paste path or click Browse…",
                        type="text",
                        size="sm",
                        debounce=False,
                        style={"fontFamily": "monospace", "fontSize": "11px"},
                    ),
                    dbc.Button(
                        [html.I(className="bi bi-folder2-open")],
                        id="btn-invr-browse",
                        color="secondary",
                        size="sm",
                        n_clicks=0,
                        title="Browse for folder",
                    ),
                ],
                size="sm",
            ),
            # pycsamt-coloured loading spinner (visible only during load)
            html.Div(
                [
                    html.Div(
                        className="invr-spinner",
                        id="invr-load-spinner",
                        style={"display": "none"},
                    ),
                    dbc.Button(
                        [
                            html.I(className="bi bi-database-fill-up me-1"),
                            "Load",
                        ],
                        id=IDs.BTN_INVR_LOAD,
                        color="primary",
                        size="sm",
                        className="w-100 mt-1",
                        n_clicks=0,
                    ),
                ]
            ),
        ],
        className="invr-load-bar",
    )


def _sidebar() -> html.Div:
    accordion = dbc.Accordion(
        [
            _acc_item(
                "invr-acc-solver",
                "bi-cpu",
                "Solver",
                *_acc_solver(),
                start_open=True,
            ),
            _acc_item(
                "invr-acc-display", "bi-sliders", "Display", *_acc_display()
            ),
            # Solver-specific extras — hidden/shown by clientside callback
            html.Div(
                _acc_item(
                    "invr-acc-modem",
                    "bi-boxes",
                    "ModEM Settings",
                    *_acc_modem_extras(),
                ),
                id="invr-acc-modem-wrap",
                style={"display": "none"},
            ),
            html.Div(
                _acc_item(
                    "invr-acc-occam2d",
                    "bi-bar-chart",
                    "Occam2D Settings",
                    *_acc_occam2d_extras(),
                ),
                id="invr-acc-occam2d-wrap",
                style={"display": "none"},
            ),
            html.Div(
                _acc_item(
                    "invr-acc-mare2dem",
                    "bi-water",
                    "MARE2DEM Settings",
                    *_acc_mare2dem_extras(),
                ),
                id="invr-acc-mare2dem-wrap",
                style={"display": "none"},
            ),
            _acc_item(
                "invr-acc-view",
                "bi-eye",
                "View Controls",
                html.Div(_acc_view_controls(), id=IDs.INVR_CTX_PANEL),
            ),
            _acc_item(
                "invr-acc-sta", "bi-pin-map", "Stations", *_acc_stations()
            ),
        ],
        id="invr-accordion",
        always_open=True,
        active_item=["invr-acc-solver", "invr-acc-display"],
        className="invr-accordion",
    )

    actions = html.Div(
        [
            dbc.Button(
                [html.I(className="bi bi-play-fill me-1"), "Generate Plot"],
                id=IDs.BTN_INVR_GEN,
                color="primary",
                size="sm",
                className="w-100 mb-1",
                n_clicks=0,
            ),
            dbc.Button(
                [html.I(className="bi bi-download me-1"), "Export PNG"],
                id=IDs.BTN_INVR_EXPORT,
                color="secondary",
                size="sm",
                outline=True,
                className="w-100",
                n_clicks=0,
            ),
            dcc.Download(id=IDs.INVR_DOWNLOAD),
        ],
        className="fwd-run-bar invr-action-bar",
    )

    return html.Div(
        [
            dcc.Store(id=IDs.INVR_STORE, storage_type="memory"),
            dcc.Store(id=IDs.INVR_ACTIVE_TAB, data=_DEFAULT_TAB),
            _load_bar(),
            html.Div(accordion, className="fwd-ctrl-scroll"),
            actions,
        ],
        className="analysis-controls fwd-sidebar",
    )


# ── Main view area ────────────────────────────────────────────────────────────


def _tab_button(
    slug: str, label: str, icon: str, active: bool = False
) -> html.Button:
    cls = f"prof-tab-btn{' active' if active else ''}"
    return html.Button(
        [html.I(className=f"bi {icon} me-1"), label],
        id=f"invr-tab-btn-{slug}",
        className=cls,
        n_clicks=0,
    )


def _view_area() -> html.Div:
    tab_bar = html.Div(
        [
            _tab_button(s, l, i, active=(s == _DEFAULT_TAB))
            for s, l, i, _ in _ALL_TABS
        ],
        id=IDs.INVR_TAB_BAR,
        className="prof-tab-bar",
    )

    info_strip = html.Div(
        [
            html.Span(
                [
                    html.I(className="bi bi-bar-chart-steps me-1"),
                    "No results loaded",
                ],
                className="adv-info-group",
            )
        ],
        id=IDs.INVR_INFO_STRIP,
        className="adv-info-bar",
    )

    # Two-layer loading: outer dcc.Loading for plot generation,
    # inner canvas shows matplotlib figure as base64 PNG.
    canvas = html.Div(
        dcc.Loading(
            html.Div(
                html.Img(
                    id=IDs.INVR_CANVAS,
                    src=empty_src(),
                    style=_IMG_STYLE,
                ),
                className="fig-img-wrap profile-page-fig",
            ),
            type="circle",
            color="#89b4fa",  # --blue
            overlay_style={
                "visibility": "visible",
                "background": "rgba(30,30,46,.55)",
            },
        ),
        className="prof-panel",
        style={"display": "flex", "flex": "1", "overflow": "hidden"},
    )

    return html.Div(
        [tab_bar, info_strip, canvas],
        style={
            "display": "flex",
            "flexDirection": "column",
            "flex": "1",
            "overflow": "hidden",
            "background": "var(--base)",
        },
    )


# ── Folder browser modal ──────────────────────────────────────────────────────


def _folder_modal() -> dbc.Modal:
    """In-app folder picker (same pattern as pipeline page)."""
    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle(
                    [
                        html.I(className="bi bi-folder2-open me-2"),
                        "Select Inversion Results Folder",
                    ]
                ),
                close_button=True,
            ),
            dbc.ModalBody(
                [
                    html.Div(
                        "Navigate into the folder that contains your inversion output files "
                        "(*.rho / STARTUP / *.emdata …), then click  Select this folder.",
                        style={
                            "fontSize": "11px",
                            "color": "var(--sub0)",
                            "marginBottom": "10px",
                        },
                    ),
                    html.Div(
                        [
                            html.Span(
                                [
                                    html.I(className="bi bi-hdd me-1"),
                                    "Current folder:",
                                ],
                                style={
                                    "fontSize": "11px",
                                    "color": "var(--sub0)",
                                    "marginRight": "6px",
                                },
                            ),
                            html.Code(
                                id="invr-browser-cwd",
                                style={
                                    "fontSize": "11px",
                                    "wordBreak": "break-all",
                                },
                            ),
                        ],
                        className="mb-2",
                    ),
                    dbc.Button(
                        [
                            html.I(className="bi bi-arrow-up me-1"),
                            "Parent folder",
                        ],
                        id="btn-invr-browser-up",
                        color="secondary",
                        outline=True,
                        size="sm",
                        className="mb-2",
                        n_clicks=0,
                    ),
                    html.Div(
                        id="invr-browser-list",
                        style={
                            "maxHeight": "340px",
                            "overflowY": "auto",
                            "border": "1px solid var(--srf0)",
                            "borderRadius": "6px",
                            "padding": "4px",
                        },
                    ),
                    dcc.Store(id="invr-browser-path", storage_type="memory"),
                ]
            ),
            dbc.ModalFooter(
                [
                    dbc.Button(
                        "Cancel",
                        id="btn-invr-browser-cancel",
                        color="secondary",
                        size="sm",
                        className="me-auto",
                        n_clicks=0,
                    ),
                    dbc.Button(
                        [
                            html.I(className="bi bi-check2-circle me-1"),
                            "Select this folder",
                        ],
                        id="btn-invr-browser-select",
                        color="success",
                        size="sm",
                        n_clicks=0,
                    ),
                ]
            ),
        ],
        id="invr-browser-modal",
        is_open=False,
        size="lg",
        backdrop="static",
    )


# ── Page entry point ──────────────────────────────────────────────────────────


def layout() -> html.Div:
    return html.Div(
        [
            _command_bar(
                "inv-results",
                "Inversion Results Viewer",
                "ModEM · Occam2D · MARE2DEM — browse, inspect, and export",
            ),
            html.Div(
                [_sidebar(), _view_area()],
                className="analysis-layout",
                style={"flex": "1"},
            ),
            _folder_modal(),
        ],
        style={
            "display": "flex",
            "flexDirection": "column",
            "height": "100%",
        },
    )
