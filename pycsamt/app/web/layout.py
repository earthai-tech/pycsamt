# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
web/layout.py — Dash component tree for the pycsamt web application.

Four-zone layout (v2 redesign):

  Navbar      — brand, breadcrumb, Load button, theme toggle, Docs / GitHub
  App Shell   — collapsible icon-rail sidebar + content body
  Content     — KPI strip (home only) · command bar · page content
  Status bar  — survey context, selected station, pyCSAMT version
"""

from __future__ import annotations

from dash import dcc, html
import dash_bootstrap_components as dbc
from dash import dash_table

from pycsamt.app.desktop.agent_registry import agent_names, agents_by_category
from pycsamt.app.web.utils import empty_src, _empty_map, _empty_plotly_heatmap

# ──────────────────────────────────────────────────────────────────────────────
# IDs — centralised to avoid typos across layout + callbacks
# ──────────────────────────────────────────────────────────────────────────────

class IDs:
    # Stores
    STORE_DATA      = "store-data"
    STORE_SELECTION = "store-selection"
    STORE_THEME     = "store-theme"
    STORE_SIDEBAR   = "store-sidebar"
    SESSION_ID      = "session-id"

    # Navbar
    BTN_LOAD        = "btn-load"
    BTN_THEME       = "btn-theme"

    # Sidebar
    BTN_SIDEBAR_TOGGLE = "btn-sidebar-toggle"
    SIDEBAR_DIV        = "sidebar-div"

    # KPI strip (home page)
    KPI_STATIONS    = "kpi-stations"
    KPI_FREQ        = "kpi-freq"
    KPI_PROFILES    = "kpi-profiles"
    KPI_SURVEY      = "kpi-survey"

    # Load modal
    MODAL_LOAD       = "modal-load"
    UPLOAD           = "upload-edi"
    UPLOAD_FILELIST  = "upload-filelist"
    INPUT_PATH       = "input-data-path"
    BTN_LOAD_CONFIRM = "btn-load-confirm"
    LOAD_FEEDBACK    = "load-feedback"
    LOAD_SPINNER     = "load-spinner"

    # Station panel
    STATION_TABLE   = "station-table"
    STATION_SUMMARY = "station-summary"

    # Map
    MAP_GRAPH       = "map-graph"
    MAP_OVERLAY     = "map-overlay-select"

    # Profile tabs
    PROFILE_TABS    = "profile-tabs"
    IMG_RHO_PHI     = "img-rho-phi"
    IMG_RHO_PS      = "img-rho-ps"
    IMG_PHASE_PS    = "img-phase-ps"
    IMG_TIPPER      = "img-tipper"
    IMG_PT          = "img-pt"

    # Period / frequency range slider
    FREQ_SLIDER     = "freq-slider"

    # QC tab (inside profile panel)
    QC_GROUP_SELECT = "qc-group-select"
    QC_PLOT_SELECT  = "qc-plot-select"
    BTN_QC_RUN      = "btn-qc-run"
    IMG_QC          = "img-qc"
    QC_SPINNER      = "qc-spinner"

    # Agent panel (home)
    AGENT_SELECT    = "agent-select"
    BTN_RUN_AGENT   = "btn-run-agent"
    AGENT_LOG       = "agent-log"
    AGENT_RESULT    = "agent-result-img"
    AGENT_SPINNER   = "agent-spinner"
    AGENT_SUMMARY   = "agent-summary"

    # Error toast
    TOAST_ERROR     = "toast-error"
    TOAST_BODY      = "toast-body"

    # Status footer
    STATUS_LABEL    = "status-label"
    STATUS_FREQ     = "status-freq"
    STATUS_STATION  = "status-station"

    # Interval (polling)
    INTERVAL        = "interval-poll"

    # Section (2D inversion result)
    IMG_SECTION     = "img-section"

    # ── Multi-page navigation ────────────────────────────────────────────
    NAV_SECTION     = "nav-section"
    CONTENT_AREA    = "content-area"

    # ── Advanced plots page ──────────────────────────────────────────────
    ADV_GROUP       = "adv-group"
    ADV_PLOT        = "adv-plot"
    BTN_ADV_RUN     = "btn-adv-run"
    IMG_ADV         = "img-adv"
    ADV_SPINNER     = "adv-spinner"

    # ── TDEM page ────────────────────────────────────────────────────────
    TDEM_CAT        = "tdem-cat"
    TDEM_PLOT       = "tdem-plot"
    TDEM_FOLDER     = "tdem-folder"
    BTN_TDEM_LOAD   = "btn-tdem-load"
    BTN_TDEM_RUN    = "btn-tdem-run"
    IMG_TDEM        = "img-tdem"
    TDEM_INFO       = "tdem-info"
    TDEM_SPINNER    = "tdem-spinner"

    # ── Correction page ──────────────────────────────────────────────────
    CORR_CAT        = "corr-cat"
    CORR_METHOD     = "corr-method"
    CORR_STACK      = "corr-stack"
    CORR_VIEW       = "corr-view"
    BTN_CORR_PREVIEW = "btn-corr-preview"
    BTN_CORR_APPLY  = "btn-corr-apply"
    BTN_CORR_UNDO   = "btn-corr-undo"
    BTN_CORR_RESET  = "btn-corr-reset"
    IMG_CORR_LEFT   = "img-corr-left"
    IMG_CORR_RIGHT  = "img-corr-right"
    IMG_CORR_OVERLAY = "img-corr-overlay"
    IMG_CORR_DIFF   = "img-corr-diff"
    CORR_SPINNER    = "corr-spinner"
    CORR_FEEDBACK   = "corr-feedback"

    # ── Pipeline page ────────────────────────────────────────────────────
    PIPE_STEP       = "pipe-step"
    PIPE_STEP_INFO  = "pipe-step-info"
    PIPE_STATUS     = "pipe-status"
    PIPE_METHOD     = "pipe-method"
    BTN_PIPE_RUN    = "btn-pipe-run"
    BTN_PIPE_ALL    = "btn-pipe-all"
    BTN_PIPE_SKIP   = "btn-pipe-skip"
    BTN_PIPE_RESET  = "btn-pipe-reset"
    PIPE_LOG        = "pipe-log"
    IMG_PIPE        = "img-pipe"
    PIPE_SPINNER    = "pipe-spinner"
    PIPE_STORE      = "pipe-store"

    # ── Forward modelling page ───────────────────────────────────────────
    FWD_DIM         = "fwd-dim"
    FWD_PRESET      = "fwd-preset"
    BTN_FWD_PRESET  = "btn-fwd-preset"
    BTN_FWD_RUN     = "btn-fwd-run"
    BTN_FWD_SAVE    = "btn-fwd-save"
    FWD_MODEL_NAME  = "fwd-model-name"
    IMG_FWD         = "img-fwd"
    FWD_SPINNER     = "fwd-spinner"
    FWD_FEEDBACK    = "fwd-feedback"
    FWD_LAYER_TABLE = "fwd-layer-table"
    BTN_FWD_ADD_LAYER = "btn-fwd-add-layer"

    # ── Inversion page ───────────────────────────────────────────────────
    INV_DIM         = "inv-dim"
    INV_SOLVER      = "inv-solver"
    INV_DATA_DIR    = "inv-data-dir"
    BTN_INV_RUN     = "btn-inv-run"
    INV_SPINNER     = "inv-spinner"
    INV_LOG         = "inv-log"
    IMG_INV         = "img-inv"

    # ── Interpretation page ──────────────────────────────────────────────
    INTERP_CAT      = "interp-cat"
    INTERP_PLOT     = "interp-plot"
    INTERP_DESC     = "interp-desc"
    BTN_INTERP_RUN  = "btn-interp-run"
    IMG_INTERP      = "img-interp"
    INTERP_SPINNER  = "interp-spinner"

    # ── Agents full page ─────────────────────────────────────────────────
    AGENTS_CAT      = "agents-cat"
    AGENTS_NAME     = "agents-name"
    AGENTS_PARAMS   = "agents-params"
    BTN_AGENTS_RUN  = "btn-agents-run"
    AGENTS_OUT      = "agents-out"
    IMG_AGENTS      = "img-agents"
    AGENTS_SPINNER  = "agents-spinner"
    AGENTS_STORE    = "agents-store"

    # ── Tools offcanvas ──────────────────────────────────────────────────
    TOOLS_CANVAS        = "tools-offcanvas"
    ACTIVE_TOOL         = "active-tool-store"
    BTN_TOOLS           = "btn-tools-menu"

    # ── Settings offcanvas ────────────────────────────────────────────────
    SETTINGS_CANVAS     = "settings-offcanvas"
    BTN_SETTINGS        = "btn-settings-menu"
    SETTINGS_STORE      = "settings-store"
    BTN_SETTINGS_RESET  = "btn-settings-reset"
    BTN_SETTINGS_SAVE   = "btn-settings-save"
    BTN_SETTINGS_LOAD   = "btn-settings-load"
    SETTINGS_UPLOAD     = "settings-profile-upload"
    SETTINGS_FEEDBACK   = "settings-feedback"

    # ── Help / About modal ────────────────────────────────────────────────
    ABOUT_MODAL         = "about-modal"
    BTN_ABOUT           = "btn-about"

    # ── Home dashboard ────────────────────────────────────────────────────
    DASH_LINE_TABLE     = "dash-line-table"
    DASH_QUALITY        = "dash-quality-area"
    DASH_LAUNCH_MAP     = "dash-launch-map"
    DASH_LAUNCH_PROFILE = "dash-launch-profile"

    # ── Map page ──────────────────────────────────────────────────────────
    MAP_PAGE_INFO        = "map-page-station-info"
    MAP_PAGE_BASEMAP     = "map-page-basemap"
    MAP_PAGE_LINE_FILTER = "map-page-line-filter"
    MAP_PAGE_MARKER_SIZE = "map-page-marker-size"
    MAP_PAGE_STATS       = "map-page-stats"

    # ── Profile page ──────────────────────────────────────────────────────
    PROFILE_PAGE_STATION = "profile-page-station"
    PROFILE_PAGE_ERRBAR  = "profile-page-errbar"
    PROFILE_PAGE_PREV    = "profile-page-prev"
    PROFILE_PAGE_NEXT    = "profile-page-next"

    # ── Welcome / landing overlay ─────────────────────────────────────────
    WELCOME_OVERLAY      = "welcome-overlay"
    WELCOME_CTA          = "welcome-cta-btn"
    WELCOME_CARD_LOAD    = "welcome-card-load"
    WELCOME_CARD_PIPE    = "welcome-card-pipeline"
    WELCOME_CARD_AGENTS  = "welcome-card-agents"
    WELCOME_CARD_VIZ     = "welcome-card-viz"


# ──────────────────────────────────────────────────────────────────────────────
# Icon helper
# ──────────────────────────────────────────────────────────────────────────────

def _icon(name: str, size: int = 15, cls: str = "panel-icon") -> html.Img:
    return html.Img(
        src=f"/icons/{name}.svg",
        className=cls,
        style={"width": f"{size}px", "height": f"{size}px"},
    )


# ──────────────────────────────────────────────────────────────────────────────
# Command Bar — shared across all inner pages
# ──────────────────────────────────────────────────────────────────────────────

def _command_bar(
    icon: str,
    title: str,
    subtitle: str = "",
    actions: list | None = None,
) -> html.Div:
    """
    Thin 36 px bar at the top of every page.

    Parameters
    ----------
    icon : str
        SVG icon name (without ``.svg``), served from ``/icons/<icon>.svg``.
    title : str
        Page title shown in bold.
    subtitle : str
        Short descriptor shown after a separator.
    actions : list
        Optional list of Dash components (buttons, dropdowns) on the right side.
    """
    left = html.Div(
        [
            html.Img(src=f"/icons/{icon}.svg", className="cmd-bar-icon"),
            html.Div(
                [
                    html.Span(title, className="cmd-bar-title"),
                    *(
                        [html.Span("·", className="cmd-bar-sep"),
                         html.Span(subtitle, className="cmd-bar-sub")]
                        if subtitle else []
                    ),
                ],
                className="cmd-bar-text",
            ),
        ],
        className="cmd-bar-left",
    )
    right = html.Div(actions or [], className="cmd-bar-actions")
    return html.Div([left, right], className="cmd-bar")


# ──────────────────────────────────────────────────────────────────────────────
# KPI Strip — home page survey summary bar
# ──────────────────────────────────────────────────────────────────────────────

def _kpi_strip() -> html.Div:
    def _card(bi_icon: str, label: str, value_id: str, default: str, color: str):
        return html.Div(
            [
                html.I(
                    className=f"bi bi-{bi_icon} kpi-icon",
                    style={"color": f"var(--{color})"},
                ),
                html.Div(
                    [
                        html.Span(default, id=value_id, className="kpi-value"),
                        html.Span(label,   className="kpi-label"),
                    ],
                    className="kpi-text",
                ),
            ],
            className="kpi-card",
        )

    return html.Div(
        [
            _card("broadcast",           "Stations",  IDs.KPI_STATIONS, "0 sites",    "blue"),
            _card("activity",            "Frequency", IDs.KPI_FREQ,     "—",          "sapphire"),
            _card("layout-three-columns","Profiles",  IDs.KPI_PROFILES, "0 profiles", "teal"),
            _card("geo-alt-fill",        "Survey",    IDs.KPI_SURVEY,   "No survey",  "peach"),
        ],
        className="kpi-strip",
    )


# ──────────────────────────────────────────────────────────────────────────────
# Overlay options (map)
# ──────────────────────────────────────────────────────────────────────────────

_OVERLAY_OPTIONS = [
    {"label": lbl, "value": lbl}
    for lbl in [
        "Index", "Apparent Resistivity", "Phase",
        "Quality Score", "Z-strike", "N_freq",
    ]
]


# ──────────────────────────────────────────────────────────────────────────────
# Component builders
# ──────────────────────────────────────────────────────────────────────────────

def _mi(icon_name: str, label: str, tool_id: str) -> dbc.DropdownMenuItem:
    """One item in the Tools dropdown — clicking stores the tool id."""
    return dbc.DropdownMenuItem(
        [_icon(icon_name, size=13, cls="menu-item-icon"), " ", label],
        id={"type": "tool-menu-item", "index": tool_id},
        n_clicks=0,
        className="navbar-menu-item",
    )


_TOOLS_MENU = dbc.DropdownMenu(
    label=html.Span(
        [_icon("tools", size=13, cls="menu-item-icon"), " Tools"],
        className="d-flex align-items-center gap-1",
    ),
    children=[
        # ── Analysis ──────────────────────────────────────────
        dbc.DropdownMenuItem("Analysis", header=True),
        _mi("strike-analyzer",        "Strike Analyzer",           "strike"),
        _mi("EDI-validator",          "EDI Validator",             "validator"),
        dbc.DropdownMenuItem(divider=True),
        # ── Conversion & Export ───────────────────────────────
        dbc.DropdownMenuItem("Conversion & Export", header=True),
        _mi("format-converter",       "Format Converter",          "converter"),
        _mi("batch-export",           "Batch Export Plots",        "batch"),
        dbc.DropdownMenuItem(divider=True),
        # ── Geospatial ────────────────────────────────────────
        dbc.DropdownMenuItem("Geospatial", header=True),
        _mi("coordinate-transformer", "Coordinate Transformer",    "coords"),
        _mi("elevation",              "Elevation Enrichment",      "elevation"),
        dbc.DropdownMenuItem(divider=True),
        # ── Visualisation ─────────────────────────────────────
        dbc.DropdownMenuItem("Visualisation", header=True),
        _mi("station-response",       "Station Response Inspector","response"),
        _mi("strike-profile",         "Strike Profile Viewer",     "strike-profile"),
        _mi("phase-tensor",           "Phase Tensor Map",          "phase-tensor"),
        dbc.DropdownMenuItem(divider=True),
        # ── Classification & Editing ──────────────────────────
        dbc.DropdownMenuItem("Classification & Editing", header=True),
        _mi("dimensionnality",        "Dimensionality Classifier", "dimensionality"),
        _mi("frequency-editor",       "Frequency Editor",          "freq-editor"),
        dbc.DropdownMenuItem(divider=True),
        # ── Modelling ─────────────────────────────────────────
        dbc.DropdownMenuItem("Modelling", header=True),
        _mi("layered-model",          "Layered Model Builder",     "layered-model"),
    ],
    id="tools-dropdown",
    nav=True,
    in_navbar=True,
    className="navbar-menu-dropdown",
    toggle_class_name="btn-navbar-menu",
    align_end=False,
)

_HELP_MENU = dbc.DropdownMenu(
    label=html.Span(
        [html.I(className="bi bi-question-circle"), " Help"],
        className="d-flex align-items-center gap-1",
    ),
    children=[
        dbc.DropdownMenuItem(
            [_icon("docs", size=13, cls="menu-item-icon"), " Documentation"],
            href="https://pycsamt.readthedocs.io",
            target="_blank",
            className="navbar-menu-item",
        ),
        dbc.DropdownMenuItem(
            [_icon("github", size=13, cls="menu-item-icon"), " pycsamt on GitHub"],
            href="https://github.com/earthai-tech/pycsamt",
            target="_blank",
            className="navbar-menu-item",
        ),
        dbc.DropdownMenuItem(divider=True),
        dbc.DropdownMenuItem(
            [_icon("help", size=13, cls="menu-item-icon"), " About pycsamt"],
            id=IDs.BTN_ABOUT,
            n_clicks=0,
            className="navbar-menu-item",
        ),
    ],
    id="help-dropdown",
    nav=True,
    in_navbar=True,
    className="navbar-menu-dropdown",
    toggle_class_name="btn-navbar-menu",
    align_end=True,
)


def _navbar() -> dbc.Navbar:
    brand = html.Div(
        html.Img(src="/icons/pycsamt_logo.svg", className="nav-logo",
                 style={"height": "26px"}),
        id="nav-brand",
        n_clicks=0,
        className="me-2 d-flex align-items-center nav-brand-btn",
        title="Go to Home",
    )

    page_indicator = html.Div(
        html.Span("Home", id="nav-page-chip", className="nav-page-chip"),
        className="ms-2 d-flex align-items-center",
    )

    # Centre menu group: Tools | Settings | Help
    menu_group = html.Div(
        [
            _TOOLS_MENU,
            html.Button(
                [_icon("tools", size=13, cls="menu-item-icon"),
                 html.Span(" Settings", style={"marginLeft": "4px"})],
                id=IDs.BTN_SETTINGS,
                className="btn-navbar-menu",
                title="Open settings panel",
            ),
            _HELP_MENU,
        ],
        className="d-flex align-items-center gap-1 mx-3",
    )

    actions = html.Div(
        [
            html.Div(className="nav-divider"),
            html.Button(
                [html.I(className="bi bi-cloud-upload-fill"),
                 html.Span(" Load Data", style={"marginLeft": "6px"})],
                id=IDs.BTN_LOAD,
                className="btn-navbar-load",
            ),
            html.Button(
                [html.I(className="bi bi-circle-half")],
                id=IDs.BTN_THEME,
                className="btn-navbar-icon ms-2",
                title="Toggle theme",
            ),
        ],
        className="d-flex align-items-center gap-1 ms-auto",
    )

    return dbc.Navbar(
        dbc.Container(
            [brand, page_indicator, menu_group, actions],
            fluid=True,
            className="d-flex align-items-center",
        ),
        color="dark",
        dark=True,
        className="mb-0",
    )


def _load_modal() -> dbc.Modal:
    # ── Section A: Folder / server-path (PRIMARY — works recursively) ──
    folder_section = html.Div([
        html.Div([
            html.Img(src="/icons/open.svg",
                     style={"width": "18px", "height": "18px",
                            "filter": "brightness(0) invert(.55)",
                            "marginRight": "7px", "verticalAlign": "middle"}),
            html.Span("Load from folder path", className="load-section-title"),
            html.Span(
                " · supports nested sub-folders as separate survey lines",
                className="load-section-hint",
            ),
        ], className="load-section-head"),
        dbc.InputGroup([
            dbc.Input(
                id=IDs.INPUT_PATH,
                placeholder="e.g.  /data/WILLY_DATA   or   /data/WILLY_DATA/L18PLT",
                type="text",
            ),
            dbc.Button(
                [html.I(className="bi bi-folder2-open me-1"), "Load"],
                id=IDs.BTN_LOAD_CONFIRM,
                color="success",
            ),
        ], className="mt-2"),
        html.Div([
            html.I(className="bi bi-info-circle me-1",
                   style={"color": "var(--blue)"}),
            html.Span(
                "Point to a parent folder (e.g. WILLY_DATA) to load all lines "
                "at once — every sub-directory becomes one survey line. "
                "Accepts .edi, .EDI, .avg and .j files recursively.",
                style={"fontSize": "11px"},
            ),
        ], className="load-info-row mt-2"),
    ], className="load-section mb-3")

    # ── Section B: Upload individual files (SECONDARY — browser only) ──
    upload_zone = dcc.Upload(
        id=IDs.UPLOAD,
        className="upload-drop",
        children=html.Div([
            html.Img(
                src="/icons/open.svg",
                style={"width": "30px", "height": "30px",
                       "filter": "brightness(0) invert(0.45)",
                       "display": "block", "margin": "0 auto 6px"},
            ),
            html.Span("Drag & drop individual EDI / AVG / J files",
                      style={"fontWeight": "500", "fontSize": "13px"}),
            html.Br(),
            html.Small("or click to browse — select multiple files with Ctrl/⌘",
                       className="text-muted"),
            html.Br(),
            html.Small(
                "Note: browser upload loads flat files only — "
                "use the folder path above to load multi-line surveys.",
                className="text-muted",
                style={"fontSize": "10px"},
            ),
        ]),
        style={"padding": "20px 18px"},
        multiple=True,
        accept=".edi,.EDI,.avg,.AVG,.j,.J",
    )

    upload_section = html.Div([
        html.Div([
            html.I(className="bi bi-cloud-upload me-1",
                   style={"color": "var(--sub0)", "fontSize": "15px",
                          "verticalAlign": "middle"}),
            html.Span("Or upload individual files", className="load-section-title"),
        ], className="load-section-head"),
        upload_zone,
        html.Div(id=IDs.UPLOAD_FILELIST, className="mt-1 text-muted small"),
    ], className="load-section")

    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle([_icon("load-edis", size=16), " Load Survey Data"])
            ),
            dbc.ModalBody([
                folder_section,
                html.Hr(className="my-2"),
                upload_section,
                html.Hr(className="my-2"),
                html.Div(id=IDs.LOAD_FEEDBACK, className="mt-2 small"),
                dbc.Spinner(html.Div(id=IDs.LOAD_SPINNER), size="sm", color="primary"),
            ]),
            dbc.ModalFooter(
                dbc.Button("Close", id="modal-close-btn", className="ms-auto")
            ),
        ],
        id=IDs.MODAL_LOAD,
        is_open=False,
        size="lg",
    )


def _station_panel() -> html.Div:
    cols = [{"name": c, "id": c}
            for c in ["ID", "Line", "Latitude", "Longitude", "Elevation", "N_freq"]]
    return html.Div(
        [
            html.Div(
                [
                    html.Span(
                        [html.Img(src="/icons/station-response.svg"),
                         "Stations"],
                        className="panel-card-title",
                    ),
                    html.Span(id=IDs.STATION_SUMMARY,
                              className="stat-chip",
                              children="0 sites"),
                ],
                className="panel-card-header",
            ),
            html.Div(
                dash_table.DataTable(
                    id=IDs.STATION_TABLE,
                    columns=cols,
                    data=[],
                    row_selectable="single",
                    page_size=25,
                    style_table={"overflowY": "auto", "flex": "1",
                                 "overflowX": "hidden"},
                    style_as_list_view=True,
                    style_header={
                        "backgroundColor": "#11111b", "color": "#6c7086",
                        "fontWeight": "800", "fontSize": "10px",
                        "border": "none", "textTransform": "uppercase",
                        "letterSpacing": "0.9px", "padding": "5px 8px",
                    },
                    style_cell={
                        "backgroundColor": "#181825", "color": "#cdd6f4",
                        "border": "none",
                        "borderBottom": "1px solid rgba(49,50,68,.55)",
                        "fontSize": "11.5px", "padding": "5px 8px",
                        "textOverflow": "ellipsis", "maxWidth": "90px",
                    },
                    style_data_conditional=[
                        {"if": {"state": "selected"},
                         "backgroundColor": "rgba(137,180,250,.12)",
                         "color": "#89b4fa"},
                        {"if": {"row_index": "odd"},
                         "backgroundColor": "rgba(17,17,27,.40)"},
                    ],
                    sort_action="native",
                    filter_action="native",
                    filter_options={"placeholder_text": "Filter…"},
                    tooltip_delay=0, tooltip_duration=None,
                ),
                className="panel-card-body",
                style={"padding": "0"},
            ),
        ],
        className="panel-card accent",
        style={"flex": "1", "overflow": "hidden"},
    )


def _map_panel() -> html.Div:
    return html.Div(
        [
            html.Div(
                [
                    html.Span(
                        [html.Img(src="/icons/map-view.svg"), "Station Map"],
                        className="panel-card-title",
                    ),
                    dbc.Select(
                        id=IDs.MAP_OVERLAY,
                        options=_OVERLAY_OPTIONS,
                        value="Index",
                        size="sm",
                        className="map-overlay-select",
                    ),
                ],
                className="panel-card-header",
            ),
            html.Div(
                dcc.Graph(
                    id=IDs.MAP_GRAPH,
                    figure=_empty_map(dark=True),
                    config={
                        "displayModeBar": True,
                        "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                        "scrollZoom": True,
                        "displaylogo": False,
                    },
                    style={"height": "100%", "minHeight": "300px"},
                    clear_on_unhover=True,
                ),
                className="panel-card-body",
                style={"padding": "0", "flex": "1", "minHeight": "0"},
            ),
        ],
        className="panel-card",
        style={"flex": "0 0 58%", "minHeight": "320px"},
    )


def _profile_panel() -> html.Div:
    _img_style = {"width": "100%", "maxHeight": "260px", "objectFit": "contain"}

    def _img_tab(img_id: str) -> html.Div:
        return html.Div(
            html.Img(id=img_id, src=empty_src(dark=True), style=_img_style),
            className="fig-img-wrap",
            style={"minHeight": "220px"},
        )

    def _graph_tab(graph_id: str) -> html.Div:
        return html.Div(
            dcc.Graph(
                id=graph_id,
                figure=_empty_plotly_heatmap(dark=True),
                config={
                    "displayModeBar": True,
                    "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                    "scrollZoom": True,
                    "displaylogo": False,
                    "toImageButtonOptions": {
                        "format": "png",
                        "filename": "pycsamt_pseudosection",
                        "scale": 2,
                    },
                },
                style={"height": "250px"},
            ),
        )

    _freq_slider = dcc.RangeSlider(
        id=IDs.FREQ_SLIDER,
        min=-4, max=4, step=0.25,
        value=[-4, 4],
        marks=None,
        tooltip={"placement": "bottom", "always_visible": False},
        disabled=True,
        className="px-1 profile-header-slider",
    )

    return html.Div(
        [
            html.Div(
                [
                    html.Span(
                        [html.Img(src="/icons/profile-view.svg"), "Profile View"],
                        className="panel-card-title",
                    ),
                    html.Span("log₁₀(T)",
                              className="freq-range-label",
                              style={"fontSize": "10px", "whiteSpace": "nowrap"}),
                    html.Div(_freq_slider, style={"flex": "1", "minWidth": "80px"}),
                ],
                className="panel-card-header",
            ),
            html.Div(
                dbc.Tabs(
                    [
                        dbc.Tab(_img_tab(IDs.IMG_RHO_PHI),    label="ρₐ / φ",      tab_id="tab-rho-phi"),
                        dbc.Tab(_graph_tab(IDs.IMG_RHO_PS),    label="ρₐ section",  tab_id="tab-rho-ps"),
                        dbc.Tab(_graph_tab(IDs.IMG_PHASE_PS),  label="φ section",   tab_id="tab-phase-ps"),
                        dbc.Tab(_img_tab(IDs.IMG_TIPPER),      label="Tipper",      tab_id="tab-tipper"),
                        dbc.Tab(_img_tab(IDs.IMG_PT),          label="Phase Tensor",tab_id="tab-pt"),
                    ],
                    id=IDs.PROFILE_TABS,
                    active_tab="tab-rho-ps",
                ),
                className="panel-card-body",
                style={"padding": "0"},
            ),
        ],
        className="panel-card",
        style={"flex": "1"},
    )


def _agent_panel() -> html.Div:
    """Home right-panel AI Agents — accordion-style sections."""
    _names = agent_names()

    def _section(header_text: str, bi_icon: str, body_children: list,
                 open_: bool = True) -> html.Div:
        return html.Div(
            [
                html.Div(
                    [
                        html.I(className=f"bi bi-{bi_icon} me-2"),
                        html.Span(header_text),
                        html.I(className="bi bi-chevron-down ms-auto"
                               + (" open" if open_ else "")),
                    ],
                    className="agent-section-header" + (" open" if open_ else ""),
                ),
                html.Div(body_children, className="agent-section-body"),
            ],
            className="agent-section",
        )

    agent_select = _section(
        "Agent", "robot",
        [
            html.Div("Agent", className="ctrl-label"),
            dbc.Select(
                id=IDs.AGENT_SELECT,
                options=[{"label": n, "value": n} for n in _names],
                value=_names[0] if _names else None,
                size="sm",
                className="mb-2",
            ),
            dbc.Button(
                [html.I(className="bi bi-play-fill me-1"), "Run"],
                id=IDs.BTN_RUN_AGENT,
                color="success", size="sm",
                className="w-100",
            ),
            html.Div(
                dbc.Spinner(html.Div(id=IDs.AGENT_SPINNER), size="sm", color="success"),
                className="mt-2 text-center",
            ),
        ],
        open_=True,
    )

    agent_log = _section(
        "Execution Log", "terminal",
        [
            html.Pre(
                "", id=IDs.AGENT_LOG,
                className="log-area",
                style={"height": "110px", "fontSize": "10.5px"},
            ),
        ],
        open_=True,
    )

    agent_result = _section(
        "Result", "image",
        [
            html.Img(
                id=IDs.AGENT_RESULT,
                src=empty_src(dark=True),
                className="agent-result-img",
                style={
                    "width": "100%", "maxHeight": "160px",
                    "objectFit": "contain",
                    "borderRadius": "4px",
                    "display": "block",
                },
            ),
            html.Div("", id=IDs.AGENT_SUMMARY,
                     className="agent-summary-txt",
                     style={"fontSize": "11px", "marginTop": "5px"}),
        ],
        open_=True,
    )

    return html.Div(
        [
            html.Div(
                html.Span(
                    [html.Img(src="/icons/agents.svg"), "AI Agents"],
                    className="panel-card-title",
                ),
                className="panel-card-header",
            ),
            html.Div(
                [agent_select, agent_log, agent_result],
                className="panel-card-body",
                style={"overflowY": "auto", "display": "flex",
                       "flexDirection": "column", "gap": "0"},
            ),
        ],
        className="panel-card",
        style={"flex": "1", "overflow": "hidden"},
    )


def _error_toast() -> dbc.Toast:
    return dbc.Toast(
        html.Div(id=IDs.TOAST_BODY),
        id=IDs.TOAST_ERROR,
        header="Error",
        is_open=False,
        dismissable=True,
        icon="danger",
        duration=8000,
        style={"position": "fixed", "top": 70, "right": 16,
               "width": 360, "zIndex": 1500},
    )



def _status_footer() -> html.Div:
    return html.Div(
        [
            html.Span(className="status-dot"),
            html.Span("Ready", id=IDs.STATUS_LABEL,
                      style={"fontSize": "11px"}),
            html.Span("|", className="status-sep"),
            html.Span("", id=IDs.STATUS_STATION,
                      style={"fontSize": "11px"}),
            html.Span("|", className="status-sep"),
            html.Span("", id=IDs.STATUS_FREQ,
                      className="status-freq-lbl",
                      style={"fontSize": "11px"}),
            html.Div(
                [
                    html.I(className="bi bi-circle-fill me-1 status-ver-dot",
                           style={"fontSize": "5px", "verticalAlign": "1px"}),
                    html.Span("pyCSAMT v2",
                              className="status-ver-txt",
                              style={"fontSize": "10px"}),
                ],
                style={"marginLeft": "auto", "display": "flex",
                       "alignItems": "center", "gap": "4px"},
            ),
        ],
        className="status-footer",
    )


# ──────────────────────────────────────────────────────────────────────────────
# Navigation entries & sidebar
# ──────────────────────────────────────────────────────────────────────────────

_NAV_ENTRIES = [
    ("home",           "summary",          "Home"),
    ("map",            "map-view",         "Map View"),
    ("profile",        "profile-view",     "Profile View"),
    ("qc",             "qc",               "Quality Control"),
    ("correction",     "sites-correction", "Correction"),
    ("advanced",       "advanced-tools",   "Advanced Plots"),
    ("tdem",           "tdem",             "TDEM"),
    ("pipeline",       "pipeline",         "Pipeline"),
    ("forward",        "forward",          "Forward Model"),
    ("inversion",      "inversion",        "Inversion"),
    ("interpretation", "interpret",        "Interpretation"),
    ("agents",         "agents",           "AI Agents"),
]

_NAV_GROUPS = [
    (None,         ["home"]),
    ("Survey",     ["map", "profile"]),
    ("Data",       ["qc", "correction"]),
    ("Analysis",   ["advanced", "tdem"]),
    ("Processing", ["pipeline"]),
    ("Modelling",  ["forward", "inversion"]),
    ("Results",    ["interpretation", "agents"]),
]

_NAV_LOOKUP = {pid: (icon, lbl) for pid, icon, lbl in _NAV_ENTRIES}


def _sidebar() -> html.Div:
    try:
        from pycsamt.metadata import __version__ as _ver
    except Exception:
        _ver = "2"

    children: list = []

    # Collapse / expand toggle button at the top
    children.append(
        html.Button(
            html.I(className="bi bi-chevron-left", id="sidebar-chevron"),
            id=IDs.BTN_SIDEBAR_TOGGLE,
            className="sidebar-toggle-btn",
            n_clicks=0,
            title="Collapse sidebar",
        )
    )

    for group_label, page_ids in _NAV_GROUPS:
        if group_label:
            children.append(
                html.Div(
                    html.Span(group_label, className="sidebar-group-label"),
                    className="sidebar-group",
                )
            )
        for page_id in page_ids:
            icon_name, label = _NAV_LOOKUP[page_id]
            children.append(
                html.Button(
                    [
                        html.Img(src=f"/icons/{icon_name}.svg",
                                 className="nav-btn-icon"),
                        html.Span(label, className="nav-btn-label"),
                    ],
                    id=f"nav-btn-{page_id}",
                    className="nav-btn" + (" active" if page_id == "home" else ""),
                    n_clicks=0,
                    title=label,
                )
            )

    children.append(
        html.Div(
            html.Span(f"v{_ver}", className="version-badge"),
            className="sidebar-footer",
        )
    )
    return html.Div(children, className="sidebar", id=IDs.SIDEBAR_DIV)


def _survey_dashboard() -> html.Div:
    """Center column of the home page: stats, quality, launch cards."""
    empty_state = html.Div(
        [
            html.I(className="bi bi-broadcast",
                   style={"fontSize": "2.5rem", "color": "var(--blue)",
                          "marginBottom": "12px"}),
            html.H6("No survey data loaded", className="dash-empty-title"),
            html.P("Load EDI files using the Load Data button above "
                   "to see survey statistics and launch the map and profile viewers.",
                   className="text-muted small text-center"),
        ],
        className="dash-empty-state",
        id="dash-empty-state",
    )

    line_table = html.Div(
        [
            html.Div(
                [
                    html.Span("Per-Line Summary", className="dash-section-label"),
                    html.Span(id="dash-survey-name", className="dash-survey-badge"),
                ],
                className="dash-section-head",
            ),
            html.Div(id=IDs.DASH_LINE_TABLE, className="dash-line-table-wrap"),
        ],
        id="dash-data-state",
        style={"display": "none"},
    )

    quality_row = html.Div(
        id=IDs.DASH_QUALITY,
        className="dash-quality-row",
        style={"display": "none"},
    )

    launch_row = html.Div(
        [
            # Map View card
            html.Div(
                [
                    html.Div(
                        html.I(className="bi bi-map", style={"fontSize": "1.6rem"}),
                        className="dash-card-icon dash-card-icon-map",
                    ),
                    html.Div(
                        [
                            html.Div("Map View", className="dash-card-title"),
                            html.Div(
                                "Interactive station map with survey line overlays, "
                                "tile basemaps, and per-station attributes.",
                                className="dash-card-desc",
                            ),
                        ],
                        className="dash-card-body",
                    ),
                    html.Div(
                        html.I(className="bi bi-arrow-right-circle"),
                        className="dash-card-arrow",
                    ),
                ],
                id=IDs.DASH_LAUNCH_MAP,
                n_clicks=0,
                className="dash-view-card dash-view-card-map",
            ),
            # Profile View card
            html.Div(
                [
                    html.Div(
                        html.I(className="bi bi-graph-up", style={"fontSize": "1.6rem"}),
                        className="dash-card-icon dash-card-icon-profile",
                    ),
                    html.Div(
                        [
                            html.Div("Profile View", className="dash-card-title"),
                            html.Div(
                                "Apparent resistivity, phase, tipper, phase tensor "
                                "and 2D section pseudosection tabs.",
                                className="dash-card-desc",
                            ),
                        ],
                        className="dash-card-body",
                    ),
                    html.Div(
                        html.I(className="bi bi-arrow-right-circle"),
                        className="dash-card-arrow",
                    ),
                ],
                id=IDs.DASH_LAUNCH_PROFILE,
                n_clicks=0,
                className="dash-view-card dash-view-card-profile",
            ),
        ],
        id="dash-launch-row",
        className="dash-launch-row",
        style={"display": "none"},
    )

    return html.Div(
        [
            html.Div(
                [
                    html.Span(
                        [html.I(className="bi bi-speedometer2 me-2"), "Survey Dashboard"],
                        className="panel-card-title",
                    ),
                    html.Span(id="dash-health-badge", className="dash-health-badge"),
                ],
                className="panel-card-header",
            ),
            html.Div(
                [empty_state, line_table, quality_row, launch_row],
                className="panel-card-body dash-body",
                style={"overflowY": "auto"},
            ),
        ],
        className="panel-card",
        style={"flex": "1"},
    )


def _home_content() -> html.Div:
    """Home page: KPI strip + 3-column workspace (station · dashboard · agents)."""
    return html.Div(
        [
            _kpi_strip(),
            html.Div(
                html.Div(
                    [
                        html.Div(_station_panel(),    className="home-col"),
                        html.Div(_survey_dashboard(), className="home-col"),
                        html.Div(_agent_panel(),      className="home-col"),
                    ],
                    className="home-grid",
                ),
                style={"flex": "1", "overflow": "hidden", "display": "flex",
                       "flexDirection": "column", "minHeight": "0"},
            ),
        ],
        id="page-home",
        className="page-content",
        style={"display": "flex"},
    )


def _map_page_content() -> html.Div:
    """Full-width Map View page: left sidebar controls + right Scattermapbox."""
    sidebar = html.Div(
        [
            # Overlay (colour-by)
            html.Div("Colour by", className="ctrl-label mt-2"),
            dbc.Select(
                id=IDs.MAP_OVERLAY,
                options=_OVERLAY_OPTIONS,
                value="Index",
                size="sm",
                className="mb-3",
            ),

            # Basemap style
            html.Div("Basemap", className="ctrl-label"),
            dbc.Select(
                id=IDs.MAP_PAGE_BASEMAP,
                options=[
                    {"label": "Dark (Carto)",      "value": "carto-darkmatter"},
                    {"label": "Light (Carto)",     "value": "carto-positron"},
                    {"label": "Open Street Map",   "value": "open-street-map"},
                    {"label": "Stamen Terrain",    "value": "stamen-terrain"},
                ],
                value="carto-darkmatter",
                size="sm",
                className="mb-3",
            ),

            # Marker size
            html.Div("Marker size", className="ctrl-label"),
            dcc.Slider(
                id=IDs.MAP_PAGE_MARKER_SIZE,
                min=6, max=22, step=2, value=10,
                marks={6: "6", 14: "14", 22: "22"},
                className="mb-3",
            ),

            html.Hr(className="tool-sep"),

            # Survey line filter
            html.Div(
                [html.Span("Survey lines", className="ctrl-label"),
                 html.Span("(all)", id="map-line-filter-count",
                           className="text-muted", style={"fontSize": "10px"})],
                className="d-flex justify-content-between align-items-center",
            ),
            dbc.Checklist(
                id=IDs.MAP_PAGE_LINE_FILTER,
                options=[],
                value=[],
                className="mb-3 map-line-filter",
            ),

            html.Hr(className="tool-sep"),

            # Selected station card
            html.Div("Selected station", className="ctrl-label"),
            html.Div(
                id=IDs.MAP_PAGE_INFO,
                className="map-page-station-card",
                children=[
                    html.Span("Click a station on the map",
                              className="text-muted small"),
                ],
            ),

            html.Hr(className="tool-sep"),

            # Quick-nav back to Home
            html.Button(
                [html.I(className="bi bi-house me-2"), "Back to Dashboard"],
                id="map-page-home-btn",
                className="btn btn-sm btn-outline-secondary w-100 mt-2",
                n_clicks=0,
            ),
        ],
        className="page-sidebar",
    )

    map_area = html.Div(
        [
            dcc.Graph(
                id=IDs.MAP_GRAPH,
                figure=_empty_map(dark=True),
                config={
                    "displayModeBar": True,
                    "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                    "scrollZoom": True,
                    "displaylogo": False,
                    "toImageButtonOptions": {
                        "format": "png", "filename": "pycsamt_map", "scale": 2,
                    },
                },
                style={"flex": "1", "minHeight": "0", "width": "100%"},
                clear_on_unhover=True,
            ),
            # Survey stats strip below the map
            html.Div(
                id=IDs.MAP_PAGE_STATS,
                className="map-stats-strip",
            ),
        ],
        className="page-main-area",
        style={"display": "flex", "flexDirection": "column"},
    )

    return html.Div(
        [sidebar, map_area],
        className="page-two-zone",
    )


def _profile_page_content() -> html.Div:
    """Full-width Profile View page: left sidebar controls + right tab panel."""
    _img_style = {"width": "100%", "height": "100%", "objectFit": "contain"}

    def _img_tab(img_id: str) -> html.Div:
        return html.Div(
            html.Img(id=img_id, src=empty_src(dark=True), style=_img_style),
            className="fig-img-wrap profile-page-fig",
        )

    def _graph_tab(graph_id: str) -> html.Div:
        return html.Div(
            dcc.Graph(
                id=graph_id,
                figure=_empty_plotly_heatmap(dark=True),
                config={
                    "displayModeBar": True,
                    "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                    "scrollZoom": True,
                    "displaylogo": False,
                    "toImageButtonOptions": {
                        "format": "png", "filename": "pycsamt_section", "scale": 2,
                    },
                },
                style={"height": "100%", "width": "100%"},
            ),
            className="profile-page-fig",
        )

    _freq_slider = dcc.RangeSlider(
        id=IDs.FREQ_SLIDER,
        min=-4, max=4, step=0.25,
        value=[-4, 4],
        marks=None,
        tooltip={"placement": "bottom", "always_visible": False},
        disabled=True,
        className="px-1",
    )

    sidebar = html.Div(
        [
            # Station selector with prev/next navigation
            html.Div("Station", className="ctrl-label mt-2"),
            html.Div(
                [
                    html.Button(
                        html.I(className="bi bi-chevron-left"),
                        id=IDs.PROFILE_PAGE_PREV,
                        className="btn btn-sm btn-outline-secondary profile-nav-btn",
                        n_clicks=0, title="Previous station",
                    ),
                    dcc.Dropdown(
                        id=IDs.PROFILE_PAGE_STATION,
                        options=[],
                        placeholder="Select station…",
                        clearable=False,
                        style={"flex": "1", "minWidth": "0"},
                    ),
                    html.Button(
                        html.I(className="bi bi-chevron-right"),
                        id=IDs.PROFILE_PAGE_NEXT,
                        className="btn btn-sm btn-outline-secondary profile-nav-btn",
                        n_clicks=0, title="Next station",
                    ),
                ],
                className="d-flex align-items-center gap-1 mb-3",
            ),

            html.Hr(className="tool-sep"),

            # Period range
            html.Div(
                [html.Span("Period range", className="ctrl-label"),
                 html.Span("log₁₀(T)", className="text-muted",
                           style={"fontSize": "10px"})],
                className="d-flex align-items-center justify-content-between mb-1",
            ),
            _freq_slider,

            html.Hr(className="tool-sep"),

            # Components
            html.Div("Components", className="ctrl-label"),
            dbc.Checklist(
                id="profile-page-comps",
                options=[
                    {"label": "Zxy (TE)", "value": "xy"},
                    {"label": "Zyx (TM)", "value": "yx"},
                    {"label": "Zxx",      "value": "xx"},
                    {"label": "Zyy",      "value": "yy"},
                ],
                value=["xy", "yx"],
                className="mb-2",
            ),

            # Error bars toggle
            dbc.Switch(
                id=IDs.PROFILE_PAGE_ERRBAR,
                label="Error bars",
                value=False,
                className="mb-3",
            ),

            html.Hr(className="tool-sep"),

            # Actions
            dbc.Button(
                [html.I(className="bi bi-arrow-repeat me-1"), "Refresh"],
                id="profile-page-refresh",
                color="primary", size="sm", className="w-100 mb-2",
                n_clicks=0,
            ),
            dbc.Button(
                [html.I(className="bi bi-download me-1"), "Export Tab"],
                id="profile-page-export",
                color="secondary", size="sm", className="w-100 mb-2",
                n_clicks=0,
            ),
            dcc.Download(id="profile-page-download"),
            html.Button(
                [html.I(className="bi bi-house me-2"), "Back to Dashboard"],
                id="profile-page-home-btn",
                className="btn btn-sm btn-outline-secondary w-100",
                n_clicks=0,
            ),
        ],
        className="page-sidebar",
    )

    tab_area = html.Div(
        dbc.Tabs(
            [
                dbc.Tab(_img_tab(IDs.IMG_RHO_PHI),    label="ρₐ / φ",       tab_id="tab-rho-phi"),
                dbc.Tab(_graph_tab(IDs.IMG_RHO_PS),    label="ρₐ section",   tab_id="tab-rho-ps"),
                dbc.Tab(_graph_tab(IDs.IMG_PHASE_PS),  label="φ section",    tab_id="tab-phase-ps"),
                dbc.Tab(_img_tab(IDs.IMG_TIPPER),      label="Tipper",       tab_id="tab-tipper"),
                dbc.Tab(_img_tab(IDs.IMG_PT),          label="Phase Tensor", tab_id="tab-pt"),
                dbc.Tab(_img_tab(IDs.IMG_SECTION),     label="2D Section",   tab_id="tab-section"),
            ],
            id=IDs.PROFILE_TABS,
            active_tab="tab-rho-ps",
            className="profile-page-tabs",
        ),
        className="page-main-area profile-page-main",
    )

    return html.Div(
        [sidebar, tab_area],
        className="page-two-zone",
    )


# ──────────────────────────────────────────────────────────────────────────────
# Global pop-out lightbox (managed entirely by assets/popout.js)
# ──────────────────────────────────────────────────────────────────────────────

def _lightbox() -> html.Div:
    """
    Full-screen lightbox shown when the user clicks a pop-out button on any
    figure.  Visibility and image ``src`` are controlled by ``popout.js``.
    """
    return html.Div(
        html.Div(
            [
                # ── Toolbar ────────────────────────────────────────────
                html.Div(
                    [
                        html.Span(
                            [html.I(className="bi bi-arrows-fullscreen me-1"),
                             "Esc  to redock"],
                            className="lb-hint",
                        ),
                        html.A(
                            [html.I(className="bi bi-download me-1"), "Save PNG"],
                            id="lightbox-download",
                            className="lb-btn",
                            href="#",
                            download="pycsamt_figure.png",
                        ),
                        html.Button(
                            [html.I(className="bi bi-arrows-angle-contract me-1"),
                             "Redock"],
                            id="lightbox-close",
                            className="lb-btn lb-close",
                            n_clicks=0,
                        ),
                    ],
                    className="lb-toolbar",
                ),
                # ── Figure ─────────────────────────────────────────────
                html.Img(id="lightbox-img", src="", className="lb-img"),
            ],
            className="lb-inner",
        ),
        id="lightbox-overlay",
        className="lightbox-overlay",
        style={"display": "none"},
    )


# ──────────────────────────────────────────────────────────────────────────────
# Welcome / landing overlay
# ──────────────────────────────────────────────────────────────────────────────

def _welcome_overlay() -> html.Div:
    """
    Full-viewport animated welcome screen displayed until the user loads data.
    Hidden via JS (adds class ``wlc-hiding`` then ``wlc-gone``) when
    STORE_DATA becomes populated.
    """

    # 14 deterministic floating particles (size_px, left%, top%, delay_s, dur_s)
    _DOTS = [
        (8,  8,  15, 0.0, 8), (5,  18, 72, 1.2, 6), (12, 28, 35, 0.4, 10),
        (4,  40, 80, 2.1, 7), (7,  52, 25, 0.8, 9), (10, 62, 60, 1.6,  8),
        (5,  73, 88, 0.2, 6), (8,  82, 42, 1.9, 9), (4,  90, 18, 0.6,  7),
        (6,  95, 65, 2.5, 8), (9,  35, 92, 1.1,10), (3,  48,  8, 0.9,  6),
        (6,  68, 78, 1.7, 8), (4,  15, 50, 2.3, 7),
    ]
    particles = html.Div([
        html.Div(className="wlc-particle", style={
            "width":             f"{sz}px",
            "height":            f"{sz}px",
            "left":              f"{lft}%",
            "top":               f"{top}%",
            "animationDelay":    f"{dly}s",
            "animationDuration": f"{dur}s",
        })
        for sz, lft, top, dly, dur in _DOTS
    ], className="wlc-particles")

    feature_cards = html.Div([
        html.Button([
            html.I(className="bi bi-cloud-upload-fill wlc-card-icon",
                   style={"color": "var(--blue)"}),
            html.Div("Load Data",             className="wlc-card-title"),
            html.Div("EDI · AVG · J · multi-line", className="wlc-card-desc"),
        ], className="wlc-card", id=IDs.WELCOME_CARD_LOAD, n_clicks=0),
        html.Button([
            html.I(className="bi bi-bezier2 wlc-card-icon",
                   style={"color": "var(--green)"}),
            html.Div("Pipeline",              className="wlc-card-title"),
            html.Div("8-step processing workflow",   className="wlc-card-desc"),
        ], className="wlc-card", id=IDs.WELCOME_CARD_PIPE, n_clicks=0),
        html.Button([
            html.I(className="bi bi-robot wlc-card-icon",
                   style={"color": "var(--mauve)"}),
            html.Div("AI Agents",             className="wlc-card-title"),
            html.Div("ML · DL · LLM solvers",       className="wlc-card-desc"),
        ], className="wlc-card", id=IDs.WELCOME_CARD_AGENTS, n_clicks=0),
        html.Button([
            html.I(className="bi bi-bar-chart-line-fill wlc-card-icon",
                   style={"color": "var(--teal)"}),
            html.Div("Visualize",             className="wlc-card-title"),
            html.Div("25 QC · 16 advanced plots",   className="wlc-card-desc"),
        ], className="wlc-card", id=IDs.WELCOME_CARD_VIZ, n_clicks=0),
    ], className="wlc-cards")

    return html.Div([
        particles,
        html.Div([
            html.Img(src="/icons/pycsamt_logo.svg", className="wlc-logo"),
            html.H1([
                "Welcome to ",
                html.Span("py",    className="wlc-py"),
                html.Span("CSAMT", className="wlc-csamt"),
            ], className="wlc-title"),
            html.P(
                "Where field data becomes geological insight",
                className="wlc-slogan",
            ),
            html.P(
                "MT/AMT processing · rich visualization · AI-powered inversion",
                className="wlc-subtitle",
            ),
            feature_cards,
            html.Button(
                [html.I(className="bi bi-cloud-upload-fill me-2"),
                 "Load EDI Data  —  Start"],
                id=IDs.WELCOME_CTA,
                className="wlc-cta",
                n_clicks=0,
            ),
            html.P(
                "Point to a parent folder to load all survey lines at once",
                className="wlc-hint",
            ),
        ], className="wlc-hero"),
    ], id=IDs.WELCOME_OVERLAY, className="wlc-overlay")


# ──────────────────────────────────────────────────────────────────────────────
# Tools off-canvas panel
# ──────────────────────────────────────────────────────────────────────────────

# Registry drives both the menu item list above and the panel content below
_TOOL_REGISTRY = [
    # (tool_id, icon_name, label, description)
    ("strike",        "strike-analyzer",        "Strike Analyzer",
     "Compute regional-strike angles and dimensionality from the impedance tensor."),
    ("validator",     "EDI-validator",           "EDI Validator",
     "Per-station data-quality checklist — flag bad impedances, missing frequencies."),
    ("converter",     "format-converter",        "Format Converter",
     "Export the loaded survey to EDI, CSV, or JSON."),
    ("batch",         "batch-export",            "Batch Export Plots",
     "Save every profile / pseudosection figure to a folder in one shot."),
    ("coords",        "coordinate-transformer",  "Coordinate Transformer",
     "Convert station coordinates between UTM and geographic Lat / Lon (pyproj)."),
    ("elevation",     "elevation",               "Elevation Enrichment",
     "Fetch SRTM elevation for every loaded station via open-elevation API."),
    ("response",      "station-response",        "Station Response Inspector",
     "Full impedance-tensor response curves for a selected station."),
    ("strike-profile","strike-profile",          "Strike Profile Viewer",
     "Strike angle vs. station-position line chart with IQR ribbon."),
    ("phase-tensor",  "phase-tensor",            "Phase Tensor Map",
     "Geographic map of phase-tensor ellipses at a chosen period."),
    ("dimensionality","dimensionnality",         "Dimensionality Classifier",
     "Label each station × frequency as 1-D / 2-D / 3-D automatically."),
    ("freq-editor",   "frequency-editor",        "Frequency Editor",
     "Confidence-based QC: drop, mask, or recover frequency bands per station."),
    ("layered-model", "layered-model",           "Layered Model Builder",
     "Build and preview a 1-D layered-earth resistivity model."),
]

_TOOL_BY_ID = {t[0]: t for t in _TOOL_REGISTRY}


def _tool_panel(tool_id: str, body=None) -> html.Div:
    """Render the body for one tool panel inside the tools offcanvas.

    Parameters
    ----------
    tool_id : str
        Key into ``_TOOL_BY_ID``.
    body : Dash component, optional
        Pre-built tool body.  When supplied the body is embedded directly;
        when omitted the inner div is left empty (populated later by a callback).
    """
    if tool_id not in _TOOL_BY_ID:
        return html.Div("Select a tool from the menu.", className="tool-empty")

    tid, icon, label, desc = _TOOL_BY_ID[tool_id]

    return html.Div([
        html.Div([
            _icon(icon, size=22, cls="tool-panel-icon"),
            html.Div([
                html.H6(label, className="tool-panel-title"),
                html.P(desc, className="tool-panel-desc"),
            ]),
        ], className="tool-panel-head"),
        html.Hr(className="tool-sep"),
        html.Div(body, className="tool-panel-body"),
    ], className="tool-panel")


def _tools_offcanvas() -> dbc.Offcanvas:
    return dbc.Offcanvas(
        [
            dcc.Store(id=IDs.ACTIVE_TOOL, data="strike"),
            # Tool picker: compact icon buttons for quick switching
            html.Div([
                html.Span("Select tool", className="offcanvas-section-lbl"),
                html.Div([
                    html.Button(
                        _icon(icon, size=14, cls="menu-item-icon"),
                        id={"type": "tool-picker-btn", "index": tid},
                        className="tool-picker-btn",
                        title=label,
                        n_clicks=0,
                    )
                    for tid, icon, label, _ in _TOOL_REGISTRY
                ], className="tool-picker-rail"),
            ], className="offcanvas-section"),
            # Active tool content — updated by callback
            html.Div(id="tools-panel-content", className="tool-panel-wrap"),
        ],
        id=IDs.TOOLS_CANVAS,
        title=html.Span([_icon("tools", size=16, cls="menu-item-icon"), " Tools"],
                        className="d-flex align-items-center gap-2"),
        placement="end",
        scrollable=True,
        is_open=False,
        className="pycsamt-offcanvas",
        style={"width": "480px"},
    )


# ──────────────────────────────────────────────────────────────────────────────
# Settings off-canvas panel
# ──────────────────────────────────────────────────────────────────────────────

def _settings_offcanvas() -> dbc.Offcanvas:
    def _section(title: str, children: list) -> html.Div:
        return html.Div([
            html.Span(title, className="offcanvas-section-lbl"),
            *children,
        ], className="offcanvas-section")

    api_section = _section("API Configuration", [
        html.Div([
            html.Label("Pseudosection interpolation", className="ctrl-label mt-2"),
            dbc.Select(
                id="settings-interp",
                options=[
                    {"label": "Bilinear",  "value": "bilinear"},
                    {"label": "Nearest",   "value": "nearest"},
                    {"label": "Cubic",     "value": "cubic"},
                ],
                value="bilinear",
                className="form-select-sm",
            ),
            html.Label("Skin-depth formula", className="ctrl-label mt-2"),
            dbc.Select(
                id="settings-skindepth",
                options=[
                    {"label": "Standard (Cagniard)",  "value": "cagniard"},
                    {"label": "True skin depth",      "value": "true"},
                ],
                value="cagniard",
                className="form-select-sm",
            ),
            html.Label("Topography DEM source", className="ctrl-label mt-2"),
            dbc.Select(
                id="settings-dem",
                options=[
                    {"label": "Open-Elevation API",    "value": "open-elevation"},
                    {"label": "SRTM local file",       "value": "srtm"},
                    {"label": "None (disable)",        "value": "none"},
                ],
                value="open-elevation",
                className="form-select-sm",
            ),
            html.Label("Default colormap", className="ctrl-label mt-2"),
            dbc.Select(
                id="settings-cmap",
                options=[{"label": c, "value": c}
                         for c in ["viridis", "plasma", "seismic",
                                   "jet", "RdBu_r", "rainbow"]],
                value="viridis",
                className="form-select-sm",
            ),
        ]),
    ])

    profile_section = _section("Configuration Profile", [
        html.Div([
            dbc.Button([_icon("save-session", size=12, cls="menu-item-icon"),
                        " Save Profile"],
                       id=IDs.BTN_SETTINGS_SAVE,
                       size="sm", color="secondary",
                       className="me-2 mt-2"),
            dcc.Upload(
                id=IDs.SETTINGS_UPLOAD,
                children=dbc.Button(
                    [_icon("open", size=12, cls="menu-item-icon"), " Load Profile"],
                    size="sm", color="secondary",
                ),
                accept=".json",
                multiple=False,
                className="d-inline-block mt-2",
            ),
        ]),
    ])

    reset_section = _section("Reset", [
        html.P("Restore all API singletons to package defaults.",
               className="text-muted small mt-1"),
        dbc.Button(
            [html.I(className="bi bi-arrow-counterclockwise me-1"),
             "Reset All to Defaults"],
            id=IDs.BTN_SETTINGS_RESET,
            size="sm",
            color="danger",
            outline=True,
            className="mt-1",
        ),
    ])

    return dbc.Offcanvas(
        [
            dcc.Store(id=IDs.SETTINGS_STORE, storage_type="local", data={}),
            api_section,
            html.Hr(className="my-2"),
            profile_section,
            html.Hr(className="my-2"),
            reset_section,
            html.Div(id=IDs.SETTINGS_FEEDBACK, className="mt-3 small"),
        ],
        id=IDs.SETTINGS_CANVAS,
        title=html.Span(
            [_icon("tools", size=16, cls="menu-item-icon"), " Settings"],
            className="d-flex align-items-center gap-2",
        ),
        placement="end",
        scrollable=True,
        is_open=False,
        className="pycsamt-offcanvas",
        style={"width": "360px"},
    )


# ──────────────────────────────────────────────────────────────────────────────
# About modal
# ──────────────────────────────────────────────────────────────────────────────

def _about_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle([_icon("help", size=16), " About pycsamt"]),
                close_button=True,
            ),
            dbc.ModalBody([
                html.Div([
                    html.Img(
                        src="/icons/pycsamt_logo.svg",
                        style={"height": "64px", "marginBottom": "12px"},
                    ),
                    html.H4([
                        html.Span("py",    style={"color": "var(--peach)",
                                                  "fontStyle": "italic"}),
                        html.Span("CSAMT", style={"color": "var(--blue)"}),
                        html.Span(" v2",   style={"color": "var(--ov0)",
                                                  "fontSize": "0.8em"}),
                    ], className="mb-1"),
                    html.P("Controlled-Source Audio Magnetotellurics · "
                           "Python Processing Suite",
                           className="text-muted small mb-3"),
                ], className="text-center"),

                html.Div([
                    html.Div([
                        html.Span("Author", className="about-key"),
                        html.Span("L. Kouadio  ·  etanoyau@gmail.com",
                                  className="about-val"),
                    ], className="about-row"),
                    html.Div([
                        html.Span("License", className="about-key"),
                        html.Span("LGPL-3.0", className="about-val"),
                    ], className="about-row"),
                    html.Div([
                        html.Span("Framework", className="about-key"),
                        html.Span("Dash 4  ·  DBC 2  ·  Plotly 6",
                                  className="about-val"),
                    ], className="about-row"),
                    html.Div([
                        html.Span("Stack", className="about-key"),
                        html.Span("NumPy · SciPy · Matplotlib · pyproj",
                                  className="about-val"),
                    ], className="about-row"),
                ], className="about-table"),
            ]),
            dbc.ModalFooter([
                html.A(
                    [_icon("docs", size=13, cls="menu-item-icon"), " Documentation"],
                    href="https://pycsamt.readthedocs.io",
                    target="_blank",
                    className="btn btn-outline-secondary btn-sm me-2",
                ),
                html.A(
                    [_icon("github", size=13, cls="menu-item-icon"), " GitHub"],
                    href="https://github.com/earthai-tech/pycsamt",
                    target="_blank",
                    className="btn btn-outline-secondary btn-sm me-2",
                ),
                dbc.Button("Close", id="about-close-btn",
                           color="secondary", size="sm"),
            ]),
        ],
        id=IDs.ABOUT_MODAL,
        is_open=False,
        size="md",
        centered=True,
    )


def layout() -> html.Div:
    """Return the complete Dash application layout."""
    from pycsamt.app.web.pages import (  # noqa: PLC0415
        qc_page, correction, advanced, tdem,
        pipeline, forward, inversion, interpretation, agents_page,
    )
    _extra_pages = [
        ("qc",             qc_page),
        ("correction",     correction),
        ("advanced",       advanced),
        ("tdem",           tdem),
        ("pipeline",       pipeline),
        ("forward",        forward),
        ("inversion",      inversion),
        ("interpretation", interpretation),
        ("agents",         agents_page),
    ]

    extra_divs = [
        html.Div(
            mod.layout(),
            id=f"page-{pid}",
            className="page-content",
            style={"display": "none"},
        )
        for pid, mod in _extra_pages
    ]

    # Inline pages (map + profile) — use layout IDs directly, no separate module
    extra_divs.insert(0, html.Div(
        _map_page_content(),
        id="page-map",
        className="page-content",
        style={"display": "none"},
    ))
    extra_divs.insert(1, html.Div(
        _profile_page_content(),
        id="page-profile",
        className="page-content",
        style={"display": "none"},
    ))

    return html.Div(
        [
            # Client-side state stores
            dcc.Store(id=IDs.SESSION_ID,      storage_type="session"),
            dcc.Store(id=IDs.STORE_DATA,      storage_type="memory"),
            dcc.Store(id=IDs.STORE_SELECTION, storage_type="memory"),
            dcc.Store(id=IDs.STORE_THEME,     storage_type="local",   data="dark"),
            dcc.Store(id=IDs.STORE_SIDEBAR,   storage_type="local",   data="expanded"),
            dcc.Store(id=IDs.NAV_SECTION,     storage_type="memory",  data="home"),

            # Polling interval (agent result streaming)
            dcc.Interval(id=IDs.INTERVAL, interval=2000,
                         n_intervals=0, disabled=True),

            # Download trigger (settings profile export)
            dcc.Download(id="settings-download"),

            # Navbar
            _navbar(),

            # Load data modal
            _load_modal(),

            # Error toast (floating, top-right)
            _error_toast(),

            # App shell: sidebar + content area
            html.Div(
                [
                    _sidebar(),
                    html.Div(
                        [_home_content()] + extra_divs,
                        id=IDs.CONTENT_AREA,
                        className="app-body",
                    ),
                ],
                className="app-shell",
            ),

            # Status footer
            _status_footer(),

            # Pop-out lightbox (managed by assets/popout.js)
            _lightbox(),

            # Tools, Settings, About panels
            _tools_offcanvas(),
            _settings_offcanvas(),
            _about_modal(),

            # Welcome / landing overlay (hidden once data is loaded)
            _welcome_overlay(),
        ],
        id="app-root",
        style={"minHeight": "100vh"},
    )
