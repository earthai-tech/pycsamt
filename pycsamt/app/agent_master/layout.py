# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Main layout for pyCSAMT Agent Master."""

from __future__ import annotations

from dash import dcc, html
import dash_bootstrap_components as dbc

from ._ids import IDs

_PROVIDER_MODELS = {
    "claude": [
        "claude-sonnet-4-6",
        "claude-opus-4-8",
        "claude-haiku-4-5-20251001",
    ],
    "openai": [
        "gpt-4o",
        "gpt-4o-mini",
        "gpt-4-turbo",
    ],
    "gemini": [
        "gemini-2.0-flash",
        "gemini-1.5-pro",
        "gemini-1.5-flash",
    ],
    "deepseek": [
        "deepseek-chat",
        "deepseek-reasoner",
    ],
}

_PROMPT_CHIPS = [
    "Load EDI data, run quality control and clean the data",
    "Run full AI-assisted 1-D neural inversion",
    "Run phase tensor and dimensionality analysis",
    "Prepare Occam2D inversion files and mesh",
    "Run full pipeline: QC, inversion, and report",
    "Decimate frequencies and select period range",
    "Run PINN physics-informed MT inversion",
    "Run hybrid two-stage AI + physics inversion",
]

_URLS = {
    "docs": (
        "https://pycsamt.readthedocs.io"
    ),
    "github": (
        "https://github.com/"
        "earthai-tech/pycsamt"
    ),
    "web": "http://localhost:8051",
}

# deterministic particles (px, left%, top%,
# delay_s, dur_s)
_DOTS = [
    ( 8,  8, 15, 0.0,  8),
    ( 5, 18, 72, 1.2,  6),
    (12, 28, 35, 0.4, 10),
    ( 4, 40, 80, 2.1,  7),
    ( 7, 52, 25, 0.8,  9),
    (10, 62, 60, 1.6,  8),
    ( 5, 73, 88, 0.2,  6),
    ( 8, 82, 42, 1.9,  9),
    ( 4, 90, 18, 0.6,  7),
    ( 6, 95, 65, 2.5,  8),
    ( 9, 35, 92, 1.1, 10),
    ( 3, 48,  8, 0.9,  6),
    ( 6, 68, 78, 1.7,  8),
    ( 4, 15, 50, 2.3,  7),
]


# ── Splash overlay ─────────────────────────────────

def _splash() -> html.Div:
    particles = html.Div(
        [
            html.Div(
                className="wlc-particle",
                style={
                    "width":  f"{sz}px",
                    "height": f"{sz}px",
                    "left":   f"{lft}%",
                    "top":    f"{top}%",
                    "animationDelay":    f"{dly}s",
                    "animationDuration": f"{dur}s",
                },
            )
            for sz, lft, top, dly, dur in _DOTS
        ],
        className="wlc-particles",
    )

    cards = html.Div(
        [
            html.Button(
                [
                    html.Img(
                        src="/am-icons/load-edis.svg",
                        className="wlc-card-icon",
                    ),
                    html.Div(
                        "Load EDI",
                        className="wlc-card-title",
                    ),
                    html.Div(
                        "EDI files · folders"
                        " · multi-line",
                        className="wlc-card-desc",
                    ),
                ],
                className="wlc-card",
                id=IDs.SPLASH_CARD_LOAD,
                n_clicks=0,
            ),
            html.Button(
                [
                    html.Img(
                        src="/am-icons/chat.svg",
                        className="wlc-card-icon",
                    ),
                    html.Div(
                        "Chat & Plan",
                        className="wlc-card-title",
                    ),
                    html.Div(
                        "Natural-language"
                        " workflow routing",
                        className="wlc-card-desc",
                    ),
                ],
                className="wlc-card",
                id=IDs.SPLASH_CARD_CHAT,
                n_clicks=0,
            ),
            html.Button(
                [
                    html.Img(
                        src="/am-icons/ai.svg",
                        className="wlc-card-icon",
                    ),
                    html.Div(
                        "AI Inversion",
                        className="wlc-card-title",
                    ),
                    html.Div(
                        "CNN · U-Net · GCN"
                        " inverters",
                        className="wlc-card-desc",
                    ),
                ],
                className="wlc-card",
                id=IDs.SPLASH_CARD_AI,
                n_clicks=0,
            ),
            html.Button(
                [
                    html.Img(
                        src="/am-icons/results.svg",
                        className="wlc-card-icon",
                    ),
                    html.Div(
                        "Reports",
                        className="wlc-card-title",
                    ),
                    html.Div(
                        "Provenance ·"
                        " reproducibility",
                        className="wlc-card-desc",
                    ),
                ],
                className="wlc-card",
                id=IDs.SPLASH_CARD_REPORT,
                n_clicks=0,
            ),
        ],
        className="wlc-cards",
    )

    return html.Div(
        [
            particles,
            html.Div(
                [
                    html.Img(
                        src=(
                            "/am-icons/"
                            "pycsamt-v2-symbol.svg"
                        ),
                        className="wlc-logo",
                    ),
                    html.H1(
                        [
                            html.Span(
                                "py",
                                className="wlc-py",
                            ),
                            html.Span(
                                "CSAMT",
                                className="wlc-csamt",
                            ),
                            " ",
                            html.Span(
                                "Agent Master",
                                className="wlc-agent",
                            ),
                        ],
                        className="wlc-title",
                    ),
                    html.P(
                        "Describe your workflow."
                        " I handle the rest.",
                        className="wlc-slogan",
                    ),
                    html.P(
                        "Natural-language AI · "
                        "MT/AMT/CSAMT processing · "
                        "AI-powered inversion · "
                        "full provenance",
                        className="wlc-subtitle",
                    ),
                    cards,
                    html.Button(
                        [
                            html.I(
                                className=(
                                    "bi bi-chat-dots-fill"
                                    " me-2"
                                )
                            ),
                            "Start Agent Master",
                        ],
                        id=IDs.SPLASH_CTA,
                        className="wlc-cta",
                        n_clicks=0,
                    ),
                    html.P(
                        "Load an EDI folder or type"
                        " a request to begin",
                        className="wlc-hint",
                    ),
                ],
                className="wlc-hero",
            ),
        ],
        id=IDs.SPLASH_OVERLAY,
        className="wlc-overlay",
    )


# ── Sidebar ────────────────────────────────────────

def _sidebar() -> html.Div:
    return html.Div(
        id=IDs.SIDEBAR,
        className="am-sidebar am-sidebar-open",
        children=[
            html.Div(
                [
                    html.Span(
                        "History",
                        className="am-sidebar-title",
                    ),
                    html.Button(
                        html.I(
                            className=(
                                "bi bi-layout-sidebar"
                            ),
                        ),
                        id="am-sidebar-close",
                        className="am-icon-btn",
                        title="Close sidebar",
                        n_clicks=0,
                    ),
                ],
                className="am-sidebar-header",
            ),
            html.Button(
                [
                    html.I(
                        className=(
                            "bi bi-plus-circle me-2"
                        ),
                    ),
                    "New Chat",
                ],
                id=IDs.BTN_NEW_CHAT,
                className="am-sidebar-new-chat",
                n_clicks=0,
            ),
            html.Div(
                "Sessions",
                className="am-sidebar-section-lbl",
            ),
            html.Div(
                html.Div(
                    "No saved sessions yet.",
                    className="am-sidebar-empty",
                ),
                id=IDs.SIDEBAR_HISTORY,
                className="am-sidebar-history",
            ),
            html.Div(
                "Pinned",
                className="am-sidebar-section-lbl",
            ),
            html.Div(
                html.Div(
                    "No pinned messages yet.",
                    className="am-sidebar-empty",
                ),
                id=IDs.SIDEBAR_PINS,
                className="am-sidebar-pins",
            ),
            html.Div(
                "Figures",
                className="am-sidebar-section-lbl",
            ),
            html.Div(
                html.Div(
                    "No figures yet.",
                    className="am-sidebar-empty",
                ),
                id=IDs.SIDEBAR_FIGS,
                className="am-sidebar-figs",
            ),
        ],
    )


# ── Top bar ────────────────────────────────────────

def _topbar() -> html.Div:
    return html.Div(
        id="am-topbar",
        children=[
            # Sidebar toggle
            html.Button(
                html.I(
                    className=(
                        "bi bi-layout-sidebar"
                    ),
                ),
                id=IDs.BTN_SIDEBAR,
                className="am-icon-btn",
                title="Toggle sidebar",
                n_clicks=0,
            ),
            # Brand
            html.Div(
                [
                    html.Img(
                        src=(
                            "/am-icons/"
                            "pycsamt-v2-symbol.svg"
                        ),
                        className="am-nav-logo",
                    ),
                    html.Span(
                        [
                            html.Span(
                                "py",
                                className="am-brand-py",
                            ),
                            html.Span(
                                "CSAMT",
                                className=(
                                    "am-brand-csamt"
                                ),
                            ),
                            " ",
                            html.Span(
                                "Agent",
                                className=(
                                    "am-brand-agent"
                                ),
                            ),
                        ],
                        className="am-brand-text",
                    ),
                ],
                className="am-brand",
            ),
            # Load EDI
            html.Button(
                [
                    html.I(
                        className=(
                            "bi bi-folder2-open me-1"
                        )
                    ),
                    "Load EDI",
                ],
                id=IDs.BTN_LOAD_EDI,
                className="am-tbtn primary",
                n_clicks=0,
            ),
            # EDI badge
            html.Div(
                [
                    html.I(
                        className=(
                            "bi bi-check-circle-fill"
                            " me-1"
                        )
                    ),
                    html.Span(
                        "",
                        id=IDs.EDI_BADGE_TEXT,
                    ),
                ],
                id=IDs.EDI_BADGE,
                className="am-edi-badge",
            ),
            # Save session
            html.Button(
                [
                    html.I(
                        className=(
                            "bi bi-cloud-check me-1"
                        )
                    ),
                    "Save",
                ],
                id=IDs.BTN_SAVE_SESSION,
                className="am-tbtn save",
                n_clicks=0,
            ),
            html.Div(
                id=IDs.SAVE_STATUS,
                className="am-save-status",
            ),
            html.Div(className="am-topbar-spacer"),
            # Theme toggle
            html.Button(
                html.I(
                    className="bi bi-moon-stars",
                    id="am-theme-icon",
                ),
                id=IDs.BTN_THEME,
                className="am-icon-btn",
                title="Toggle light / dark",
                n_clicks=0,
            ),
            # Settings
            html.Button(
                html.I(
                    className="bi bi-gear-fill"
                ),
                id=IDs.BTN_SETTINGS,
                className="am-icon-btn",
                title="Settings",
                n_clicks=0,
            ),
        ],
    )


# ── In-chat welcome (post-splash) ─────────────────

def _chat_welcome() -> html.Div:
    chips = [
        html.Button(
            c,
            className="am-chip",
            id={"type": "am-chip", "index": i},
            n_clicks=0,
        )
        for i, c in enumerate(_PROMPT_CHIPS)
    ]
    return html.Div(
        id=IDs.WELCOME,
        className="am-chat-welcome",
        children=[
            html.Img(
                src=(
                    "/am-icons/"
                    "pycsamt-v2-symbol.svg"
                ),
                className="am-welcome-logo",
            ),
            html.H3(
                "pyCSAMT Agent Master"
            ),
            html.P(
                "Load an EDI dataset and describe "
                "your workflow in natural language."
            ),
            html.Div(
                chips,
                className="am-prompt-chips",
            ),
        ],
    )


# ── "+" resource menu ──────────────────────────────

def _plus_menu() -> html.Div:

    def _ext(icon, label, href):
        return html.A(
            [
                html.I(
                    className=f"bi {icon} me-2"
                ),
                html.Span(label),
                html.I(
                    className=(
                        "bi bi-arrow-up-right"
                        " ms-auto"
                    ),
                ),
            ],
            href=href,
            target="_blank",
            rel="noopener noreferrer",
            className=(
                "am-plus-item am-plus-ext"
            ),
        )

    def _action(icon, label, btn_id):
        return html.Button(
            [
                html.I(
                    className=f"bi {icon} me-2"
                ),
                html.Span(label),
            ],
            id=btn_id,
            className="am-plus-item",
            n_clicks=0,
        )

    return html.Div(
        [
            html.Div(
                "Resources",
                className="am-plus-section",
            ),
            _ext(
                "bi-book-half",
                "Documentation",
                _URLS["docs"],
            ),
            _ext(
                "bi-github",
                "GitHub",
                _URLS["github"],
            ),
            _action(
                "bi-window-fullscreen",
                "Launch Web App",
                IDs.PLUS_WEB,
            ),
            html.Hr(
                className="am-plus-divider"
            ),
            html.Div(
                "Open panel",
                className="am-plus-section",
            ),
            _action(
                "bi-folder2-open",
                "Load EDI files",
                IDs.PLUS_LOAD,
            ),
            _action(
                "bi-gear-fill",
                "Settings",
                IDs.PLUS_SETTINGS,
            ),
            html.Hr(
                className="am-plus-divider"
            ),
            html.Div(
                "Insert",
                className="am-plus-section",
            ),
            _action(
                "bi-cursor-text",
                "Paste path template",
                IDs.PLUS_PATH,
            ),
        ],
        id=IDs.PLUS_POPOVER,
        className="am-plus-dropdown",
    )


# ── Input bar ──────────────────────────────────────

def _input_bar() -> html.Div:
    return html.Div(
        id="am-input-bar",
        children=[
            _plus_menu(),
            html.Div(
                [
                    html.Button(
                        html.I(
                            className=(
                                "bi bi-plus-lg"
                            )
                        ),
                        id=IDs.BTN_PLUS,
                        className="am-plus-btn",
                        title=(
                            "Resources & actions"
                        ),
                        n_clicks=0,
                    ),
                    html.Button(
                        html.I(
                            className=(
                                "bi bi-paperclip"
                            )
                        ),
                        id=IDs.BTN_ATTACH,
                        className="am-attach-btn",
                        title="Load EDI file",
                        n_clicks=0,
                    ),
                    dcc.Textarea(
                        id=IDs.INPUT,
                        placeholder=(
                            "Describe your "
                            "workflow..."
                        ),
                        style={},
                    ),
                    html.Button(
                        html.I(
                            className=(
                                "bi bi-arrow-up"
                            )
                        ),
                        id=IDs.BTN_SEND,
                        n_clicks=0,
                        disabled=False,
                    ),
                ],
                className="am-input-wrap",
            ),
        ],
    )


# ── EDI offcanvas ──────────────────────────────────

def _edi_canvas() -> dbc.Offcanvas:
    mode_opts = [
        {
            "label": html.Span(
                [
                    html.I(
                        className=(
                            "bi bi-folder2-open"
                            " me-2"
                        ),
                        style={"color": "#f9e2af"},
                    ),
                    "Folder names (default)",
                ],
                className="d-flex align-items-center",
            ),
            "value": "folder",
        },
        {
            "label": html.Span(
                [
                    html.I(
                        className="bi bi-stars me-2",
                        style={"color": "#89b4fa"},
                    ),
                    "Auto-detect from IDs",
                ],
                className="d-flex align-items-center",
            ),
            "value": "auto",
        },
        {
            "label": html.Span(
                [
                    html.I(
                        className=(
                            "bi bi-pencil-square"
                            " me-2"
                        ),
                        style={"color": "#a6e3a1"},
                    ),
                    "Edit / Rename",
                ],
                className="d-flex align-items-center",
            ),
            "value": "edit",
        },
    ]
    return dbc.Offcanvas(
        [
            # ── Source selection ──────────────────
            html.Div(
                "Select EDI source",
                className="am-section-lbl",
            ),
            # Big browse-folder button (primary)
            dbc.Button(
                [
                    html.I(
                        className=(
                            "bi bi-folder2-open"
                            " me-2"
                        )
                    ),
                    html.Span(
                        "Browse folder",
                        style={"fontWeight": "600"},
                    ),
                    html.Span(
                        " — navigate & select",
                        style={
                            "fontSize": "11px",
                            "opacity": ".75",
                        },
                    ),
                ],
                id=IDs.BTN_BROWSE,
                color="primary",
                outline=False,
                size="sm",
                className="w-100 mb-2",
                n_clicks=0,
            ),
            # Drop-zone for individual files
            dcc.Upload(
                id=IDs.UPLOAD_EDI,
                children=html.Div(
                    [
                        html.I(
                            className=(
                                "bi bi-file-earmark"
                                "-plus me-2"
                            )
                        ),
                        html.Span(
                            "Drop EDI files"
                            " or click to pick"
                        ),
                    ]
                ),
                className="am-upload-zone mb-2",
                multiple=True,
                accept=".edi,.EDI",
            ),
            # Passive path display / manual edit
            html.Div(
                "Selected path",
                className="am-section-lbl mt-1",
            ),
            dbc.Input(
                id=IDs.EDI_PATH_INPUT,
                placeholder=(
                    "Path appears here after"
                    " browsing — or type one"
                ),
                size="sm",
                className="mb-3",
                style={
                    "fontSize": "11px",
                    "fontFamily": "monospace",
                },
            ),
            html.Hr(className="my-2"),
            # ── Line assignment ───────────────────
            html.Div(
                "Line assignment mode",
                className="am-section-lbl",
            ),
            dbc.RadioItems(
                id=IDs.LINES_MODE,
                options=mode_opts,
                value="folder",
                className="mb-2",
                inputStyle={"cursor": "pointer"},
                labelStyle={
                    "cursor": "pointer",
                    "fontSize": "12px",
                    "padding": "2px 0",
                },
            ),
            dbc.Button(
                [
                    html.I(
                        className="bi bi-stars me-2"
                    ),
                    "Detect lines",
                ],
                id=IDs.BTN_DETECT_LINES,
                color="info",
                size="sm",
                outline=True,
                className="w-100 mb-1",
                style={"display": "none"},
                n_clicks=0,
            ),
            dbc.Button(
                [
                    html.I(
                        className=(
                            "bi bi-check2-circle"
                            " me-2"
                        )
                    ),
                    "Apply renames",
                ],
                id=IDs.BTN_APPLY_RENAME,
                color="success",
                size="sm",
                outline=True,
                className="w-100 mb-1",
                style={"display": "none"},
                n_clicks=0,
            ),
            html.Hr(className="my-2"),
            html.Div(
                id=IDs.LINES_PANEL,
                className="mb-3",
            ),
            html.Div(
                id=IDs.EDI_LOAD_STATUS,
                style={"fontSize": "12px"},
            ),
            html.Hr(className="my-2"),
            dbc.Button(
                [
                    html.I(
                        className=(
                            "bi bi-check-lg me-2"
                        )
                    ),
                    "Load into session",
                ],
                id=IDs.BTN_LOAD_CONFIRM,
                color="primary",
                size="sm",
                className="w-100",
                n_clicks=0,
            ),
        ],
        id=IDs.CANVAS_EDI,
        title="Load EDI Data",
        placement="start",
        is_open=False,
        style={"width": "360px"},
        className="am-panel",
    )


# ── Settings offcanvas ─────────────────────────────

def _settings_canvas() -> dbc.Offcanvas:
    def _key_tab(
        provider, key_id,
        models=None, model_id=None,
    ):
        children = [
            dbc.Label(
                "API Key",
                style={"fontSize": "12px"},
            ),
            dbc.InputGroup(
                [
                    dbc.Input(
                        id=key_id,
                        type="password",
                        placeholder=f"key ({provider})",
                        size="sm",
                        autocomplete="off",
                    ),
                    dbc.InputGroupText(
                        html.I(
                            className="bi bi-key"
                        ),
                        style={"fontSize": "13px"},
                    ),
                ],
                size="sm",
                className="mb-3",
            ),
        ]
        if models and model_id:
            children += [
                dbc.Label(
                    "Model",
                    style={"fontSize": "12px"},
                ),
                dbc.Select(
                    id=model_id,
                    options=[
                        {"label": m, "value": m}
                        for m in models
                    ],
                    value=models[0],
                    size="sm",
                    className="mb-3",
                ),
            ]
        return html.Div(
            children, className="pt-3"
        )

    tabs = dbc.Tabs(
        [
            dbc.Tab(
                _key_tab(
                    "Claude",
                    IDs.KEY_CLAUDE,
                    _PROVIDER_MODELS["claude"],
                    IDs.MODEL_CLAUDE,
                ),
                label="Claude",
                tab_id="claude",
            ),
            dbc.Tab(
                _key_tab(
                    "OpenAI",
                    IDs.KEY_OPENAI,
                    _PROVIDER_MODELS["openai"],
                    IDs.MODEL_OPENAI,
                ),
                label="OpenAI",
                tab_id="openai",
            ),
            dbc.Tab(
                _key_tab(
                    "Gemini",
                    IDs.KEY_GEMINI,
                    _PROVIDER_MODELS["gemini"],
                    IDs.MODEL_GEMINI,
                ),
                label="Gemini",
                tab_id="gemini",
            ),
            dbc.Tab(
                _key_tab(
                    "DeepSeek",
                    IDs.KEY_DEEPSEEK,
                    _PROVIDER_MODELS["deepseek"],
                    IDs.MODEL_DEEPSEEK,
                ),
                label="DeepSeek",
                tab_id="deepseek",
            ),
        ],
        id=IDs.PROVIDER_TABS,
        active_tab="claude",
        className="mb-2",
    )

    return dbc.Offcanvas(
        [
            html.Div(
                "LLM Providers",
                className="am-section-lbl",
            ),
            tabs,
            html.Div(
                id=IDs.KEYS_STATUS,
                style={
                    "fontSize": "11px",
                    "marginBottom": "6px",
                },
            ),
            dbc.Button(
                [
                    html.I(
                        className=(
                            "bi bi-floppy me-2"
                        )
                    ),
                    "Save API keys",
                ],
                id=IDs.BTN_SAVE_KEYS,
                color="primary",
                size="sm",
                className="w-100 mb-4",
                n_clicks=0,
            ),
            html.Hr(className="my-2"),
            html.Div(
                "Default figure export",
                className="am-section-lbl",
            ),
            dbc.Select(
                id=IDs.EXPORT_FORMAT,
                options=[
                    {"label": "PNG (300 dpi)",
                     "value": "png"},
                    {"label": "SVG (vector)",
                     "value": "svg"},
                    {"label": "EPS (vector)",
                     "value": "eps"},
                    {"label": "PDF",
                     "value": "pdf"},
                ],
                value="png",
                size="sm",
                className="mb-4",
            ),
            html.Hr(className="my-2"),
            html.Div(
                "Active LLM provider",
                className="am-section-lbl",
            ),
            dbc.RadioItems(
                id=IDs.ACTIVE_PROVIDER,
                options=[
                    {
                        "label": "Claude (Anthropic)",
                        "value": "claude",
                    },
                    {
                        "label": "OpenAI",
                        "value": "openai",
                    },
                    {
                        "label": "Gemini",
                        "value": "gemini",
                    },
                    {
                        "label": "DeepSeek",
                        "value": "deepseek",
                    },
                    {
                        "label": "Offline (rule-based)",
                        "value": "offline",
                    },
                ],
                value="offline",
                labelStyle={
                    "fontSize": "13px",
                    "padding": "2px 0",
                    "cursor": "pointer",
                },
                inputStyle={"cursor": "pointer"},
            ),
            html.Hr(className="my-2"),
            html.Div(
                "Results output directory",
                className="am-section-lbl",
            ),
            html.Div(
                [
                    dbc.Input(
                        id=IDs.OUTPUT_DIR,
                        placeholder=(
                            "Default: "
                            "pycsamt_workflow_output/"
                        ),
                        size="sm",
                        style={
                            "flex": "1",
                            "borderRadius":
                            "6px 0 0 6px",
                            "fontSize": "12px",
                        },
                        readonly=True,
                    ),
                    dbc.Button(
                        html.I(
                            className=(
                                "bi bi-folder2-open"
                            )
                        ),
                        id=IDs.BTN_OUTPUT_BROWSE,
                        size="sm",
                        color="secondary",
                        outline=True,
                        title="Browse / create folder",
                        style={
                            "borderRadius":
                            "0 6px 6px 0",
                            "borderLeft": "none",
                        },
                        n_clicks=0,
                    ),
                ],
                className="d-flex mb-1",
            ),
            html.Div(
                "Click the folder icon to browse"
                " or create an output folder.",
                style={
                    "fontSize": "10px",
                    "color": "var(--am-muted)",
                    "marginBottom": "4px",
                },
            ),
            html.Hr(className="my-2"),
            html.Div(
                "Line Registry (YAML)",
                className="am-section-lbl",
            ),
            dbc.Textarea(
                id=IDs.LINE_REGISTRY,
                placeholder=(
                    "L22PLT: /data/willy/L22PLT\n"
                    "L18PLT: /data/willy/L18PLT"
                ),
                rows=4,
                style={
                    "fontSize": "11px",
                    "fontFamily": "monospace",
                    "resize": "vertical",
                },
                className="mb-1",
            ),
            html.Div(
                "Map survey line names to EDI"
                " directories (YAML format)."
                " Agent Master resolves them"
                " automatically when you name"
                " a line in the chat.",
                style={
                    "fontSize": "10px",
                    "color": "var(--am-muted)",
                    "marginBottom": "4px",
                },
            ),
        ],
        id=IDs.CANVAS_SETTINGS,
        title="Settings",
        placement="end",
        is_open=False,
        style={"width": "340px"},
        className="am-panel",
    )


# ── Line selector modal ────────────────────────────

def _line_sel_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                [
                    html.I(
                        className=(
                            "bi bi-layers me-2"
                        )
                    ),
                    "Select line(s) to process",
                ],
                close_button=True,
            ),
            dbc.ModalBody(
                [
                    html.P(
                        "Multiple survey lines are"
                        " loaded. Choose which"
                        " line(s) to process:",
                        style={"fontSize": "13px"},
                    ),
                    dbc.Checklist(
                        id=IDs.LINE_SEL_LIST,
                        options=[],
                        value=[],
                        labelStyle={
                            "fontSize": "13px",
                            "padding": "2px 0",
                        },
                        className="mb-3",
                        switch=True,
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dbc.Button(
                                    "Run selected",
                                    id=(
                                        IDs.BTN_LINE_RUN_SEL
                                    ),
                                    color="primary",
                                    size="sm",
                                    className="w-100",
                                    n_clicks=0,
                                ),
                            ),
                            dbc.Col(
                                dbc.Button(
                                    "Run all lines",
                                    id=(
                                        IDs.BTN_LINE_RUN_ALL
                                    ),
                                    color="secondary",
                                    size="sm",
                                    className="w-100",
                                    n_clicks=0,
                                ),
                            ),
                        ],
                    ),
                    html.Div(
                        id=IDs.LINE_SEL_STATUS,
                        style={
                            "fontSize": "12px",
                            "marginTop": "6px",
                        },
                    ),
                ]
            ),
        ],
        id=IDs.MODAL_LINE_SEL,
        is_open=False,
        size="sm",
        className="am-modal",
    )


# ── Post-correction action modal ───────────────────

def _postproc_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                "Processing complete",
                close_button=True,
            ),
            dbc.ModalBody(
                [
                    html.P(
                        id=IDs.POSTPROC_MSG,
                        style={"fontSize": "13px"},
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dbc.Button(
                                    "Apply to session",
                                    id=(
                                        IDs.BTN_POSTPROC_APPLY
                                    ),
                                    color="primary",
                                    size="sm",
                                    className="w-100",
                                    n_clicks=0,
                                ),
                            ),
                            dbc.Col(
                                dbc.Button(
                                    "Export to folder",
                                    id=(
                                        IDs.BTN_POSTPROC_EXPORT
                                    ),
                                    color="secondary",
                                    size="sm",
                                    className="w-100",
                                    n_clicks=0,
                                ),
                            ),
                        ],
                        className="mb-3",
                    ),
                    dbc.Collapse(
                        [
                            dbc.Input(
                                id=IDs.POSTPROC_PATH,
                                placeholder=(
                                    "Enter export"
                                    " folder path..."
                                ),
                                size="sm",
                                className="mb-2",
                            ),
                            dbc.Button(
                                "Confirm export",
                                id=(
                                    IDs.BTN_POSTPROC_OK
                                ),
                                color="success",
                                size="sm",
                                n_clicks=0,
                            ),
                        ],
                        id=IDs.POSTPROC_COLLAPSE,
                        is_open=False,
                    ),
                    html.Div(
                        id=IDs.POSTPROC_STATUS,
                        style={
                            "fontSize": "12px",
                            "marginTop": "6px",
                        },
                    ),
                ]
            ),
        ],
        id=IDs.MODAL_POSTPROC,
        is_open=False,
        size="sm",
        className="am-modal",
    )


# ── Output folder browser modal ────────────────────

def _output_browse_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                [
                    html.I(
                        className=(
                            "bi bi-folder2-open me-2"
                        )
                    ),
                    "Select output folder",
                ],
                close_button=True,
            ),
            dbc.ModalBody(
                [
                    # breadcrumb current path
                    html.Div(
                        [
                            dbc.Button(
                                html.I(
                                    className=(
                                        "bi bi-arrow-up"
                                    )
                                ),
                                id=IDs.BTN_OUTPUT_UP,
                                size="sm",
                                color="light",
                                outline=True,
                                className="me-2",
                                n_clicks=0,
                                title="Go up",
                            ),
                            html.Code(
                                id=IDs.OUTPUT_BROWSE_PATH,
                                style={
                                    "fontSize": "11px",
                                    "wordBreak":
                                    "break-all",
                                },
                            ),
                        ],
                        className=(
                            "d-flex align-items-center"
                            " mb-2"
                        ),
                    ),
                    # directory listing
                    html.Div(
                        id=IDs.OUTPUT_BROWSE_LIST,
                        style={
                            "maxHeight": "260px",
                            "overflowY": "auto",
                            "border": (
                                "1px solid"
                                " var(--am-border)"
                            ),
                            "borderRadius": "6px",
                            "padding": "6px",
                        },
                    ),
                    html.Hr(className="my-2"),
                    # create new subfolder
                    html.Small(
                        "Create new subfolder here:",
                        className="am-section-lbl",
                    ),
                    html.Div(
                        [
                            dbc.Input(
                                id=IDs.OUTPUT_MKDIR_INPUT,
                                placeholder=(
                                    "new_folder_name"
                                ),
                                size="sm",
                                style={
                                    "flex": "1",
                                    "borderRadius":
                                    "6px 0 0 6px",
                                },
                                debounce=False,
                            ),
                            dbc.Button(
                                [
                                    html.I(
                                        className=(
                                            "bi bi-folder-plus"
                                            " me-1"
                                        )
                                    ),
                                    "Create",
                                ],
                                id=IDs.BTN_OUTPUT_MKDIR,
                                size="sm",
                                color="secondary",
                                style={
                                    "borderRadius":
                                    "0 6px 6px 0",
                                },
                                n_clicks=0,
                            ),
                        ],
                        className="d-flex mt-1 mb-1",
                    ),
                    html.Div(
                        id=IDs.OUTPUT_BROWSE_STATUS,
                        style={
                            "fontSize": "11px",
                            "minHeight": "16px",
                            "color": "var(--am-muted)",
                        },
                    ),
                ]
            ),
            dbc.ModalFooter(
                dbc.Button(
                    [
                        html.I(
                            className=(
                                "bi bi-check2-circle"
                                " me-1"
                            )
                        ),
                        "Use this folder",
                    ],
                    id=IDs.BTN_OUTPUT_CONFIRM,
                    color="primary",
                    size="sm",
                    n_clicks=0,
                ),
            ),
        ],
        id=IDs.MODAL_OUTPUT_BROWSE,
        is_open=False,
        size="md",
        className="am-modal",
    )


# ── Figure modal ───────────────────────────────────

def _fig_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                html.Span(
                    id=IDs.MODAL_FIG_TITLE,
                    style={
                        "fontSize": "14px",
                        "fontWeight": "600",
                    },
                ),
                className="am-modal-header",
                close_button=True,
            ),
            dbc.ModalBody(
                html.Img(
                    id=IDs.MODAL_FIG_IMG,
                    style={
                        "maxWidth": "100%",
                        "borderRadius": "6px",
                    },
                ),
                className="am-modal-body",
            ),
            dbc.ModalFooter(
                [
                    dcc.Store(
                        id=IDs.MODAL_FIG_KEY
                    ),
                    html.Small(
                        "Export:",
                        style={
                            "color": "var(--ov0)",
                            "marginRight": "6px",
                        },
                    ),
                    dbc.ButtonGroup(
                        [
                            dbc.Button(
                                fmt.upper(),
                                id={
                                    "type": "am-export",
                                    "fmt": fmt,
                                },
                                size="sm",
                                color="secondary",
                                outline=True,
                                n_clicks=0,
                            )
                            for fmt in [
                                "png", "svg",
                                "eps", "pdf",
                            ]
                        ],
                        size="sm",
                    ),
                    dcc.Download(
                        id=IDs.EXPORT_DL
                    ),
                ],
                className="am-modal-footer",
            ),
        ],
        id=IDs.MODAL_FIG,
        size="xl",
        is_open=False,
        centered=True,
        contentClassName="am-modal-content",
    )



# ── Smart param-collection modal ───────────────────

def _param_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                html.Div(
                    [
                        html.I(
                            id="am-param-icon",
                            className=(
                                "bi bi-sliders2 me-2"
                            ),
                        ),
                        html.Span(
                            "Configure Workflow",
                            id=IDs.PARAM_MODAL_TITLE,
                        ),
                    ],
                    className=(
                        "d-flex align-items-center"
                    ),
                ),
                className="am-modal-header",
                close_button=True,
            ),
            dbc.ModalBody(
                [
                    html.P(
                        id=IDs.PARAM_MODAL_DESC,
                        className="am-param-desc",
                    ),
                    html.Div(
                        id=IDs.PARAM_FORM_BODY,
                        className="am-param-form",
                    ),
                ],
                className="am-modal-body",
            ),
            dbc.ModalFooter(
                [
                    dbc.Button(
                        [
                            html.I(
                                className=(
                                    "bi bi-x-circle"
                                    " me-1"
                                )
                            ),
                            "Cancel",
                        ],
                        id=IDs.BTN_PARAM_CANCEL,
                        color="secondary",
                        outline=True,
                        size="sm",
                        n_clicks=0,
                    ),
                    dbc.Button(
                        [
                            html.I(
                                className=(
                                    "bi bi-play-fill"
                                    " me-1"
                                )
                            ),
                            "Run Workflow",
                        ],
                        id=IDs.BTN_PARAM_RUN,
                        color="success",
                        size="sm",
                        n_clicks=0,
                    ),
                ],
                className="am-modal-footer",
            ),
        ],
        id=IDs.MODAL_PARAMS,
        is_open=False,
        size="lg",
        backdrop="static",
        centered=True,
        scrollable=True,
        contentClassName="am-modal-content",
    )


# ── PINN / Hybrid parameter offcanvas ─────────────

def _inv_params_canvas() -> dbc.Offcanvas:
    def _label(text):
        return html.Label(
            text,
            className="am-field-label",
        )

    def _row(label, comp):
        return html.Div(
            [_label(label), comp],
            className="am-field-row",
            style={"marginBottom": "10px"},
        )

    return dbc.Offcanvas(
        id=IDs.CANVAS_INV_PARAMS,
        title="Inversion Parameters",
        placement="end",
        is_open=False,
        style={"width": "360px"},
        children=[
            # ── Mode ───────────────────────────────
            html.Div(
                [
                    _label("Inversion mode"),
                    dbc.RadioItems(
                        id=IDs.INV_MODE,
                        options=[
                            {"label": "PINN (physics-informed)",
                             "value": "pinn"},
                            {"label": "Hybrid (AI warm-start)",
                             "value": "hybrid"},
                        ],
                        value="pinn",
                        inline=True,
                        className="am-radio-group",
                    ),
                ],
                style={"marginBottom": "12px"},
            ),
            # ── Dimension ──────────────────────────
            html.Div(
                [
                    _label("Dimension"),
                    dbc.RadioItems(
                        id=IDs.INV_DIM,
                        options=[
                            {"label": "1-D", "value": 1},
                            {"label": "2-D", "value": 2},
                            {"label": "3-D", "value": 3},
                        ],
                        value=1,
                        inline=True,
                        className="am-radio-group",
                    ),
                ],
                style={"marginBottom": "12px"},
            ),
            # ── Solver (1-D only) ──────────────────
            _row(
                "Solver (1-D)",
                dbc.Select(
                    id=IDs.INV_SOLVER,
                    options=[
                        {"label": "mt1d", "value": "mt1d"},
                        {"label": "csamt1d",
                         "value": "csamt1d"},
                    ],
                    value="mt1d",
                ),
            ),
            html.Hr(style={"margin": "8px 0"}),
            # ── Earth model ───────────────────────
            html.Small(
                "Earth model",
                className="am-section-heading",
            ),
            _row(
                "Number of layers",
                dbc.Input(
                    id=IDs.INV_N_LAYERS,
                    type="number",
                    value=10,
                    min=2,
                    max=100,
                    step=1,
                ),
            ),
            _row(
                "Max depth (m)",
                dbc.Input(
                    id=IDs.INV_DEPTH_MAX,
                    type="number",
                    value=2000.0,
                    min=100,
                    step=50,
                ),
            ),
            html.Hr(style={"margin": "8px 0"}),
            # ── Optimiser ────────────────────────
            html.Small(
                "Optimiser",
                className="am-section-heading",
            ),
            _row(
                "Epochs / max iterations",
                dbc.Input(
                    id=IDs.INV_EPOCHS,
                    type="number",
                    value=500,
                    min=10,
                    step=10,
                ),
            ),
            _row(
                "Learning rate",
                dbc.Input(
                    id=IDs.INV_LR,
                    type="number",
                    value=0.01,
                    min=1e-5,
                    step=0.001,
                ),
            ),
            html.Hr(style={"margin": "8px 0"}),
            # ── Regularisation ────────────────────
            html.Small(
                "Regularisation",
                className="am-section-heading",
            ),
            _row(
                "Smoothness weight",
                dbc.Input(
                    id=IDs.INV_SMOOTH_W,
                    type="number",
                    value=0.01,
                    min=0,
                    step=0.001,
                ),
            ),
            # 2-D and 3-D panel
            html.Div(
                id=IDs.INV_PANEL_2D,
                style={"display": "none"},
                children=[
                    _row(
                        "Lateral weight (2D+)",
                        dbc.Input(
                            id=IDs.INV_LAT_W,
                            type="number",
                            value=0.005,
                            min=0,
                            step=0.001,
                        ),
                    ),
                ],
            ),
            # 3-D panel
            html.Div(
                id=IDs.INV_PANEL_3D,
                style={"display": "none"},
                children=[
                    _row(
                        "Graph weight (3D)",
                        dbc.Input(
                            id=IDs.INV_GRAPH_W,
                            type="number",
                            value=0.005,
                            min=0,
                            step=0.001,
                        ),
                    ),
                    _row(
                        "Station radius (m)",
                        dbc.Input(
                            id=IDs.INV_RADIUS,
                            type="number",
                            value=5000.0,
                            min=100,
                            step=100,
                        ),
                    ),
                ],
            ),
            # Hybrid checkpoint panel
            html.Div(
                id=IDs.INV_PANEL_HYBRID,
                style={"display": "none"},
                children=[
                    html.Hr(
                        style={"margin": "8px 0"}
                    ),
                    html.Small(
                        "Hybrid: AI inverter",
                        className="am-section-heading",
                    ),
                    _row(
                        "Checkpoint path"
                        " (optional)",
                        dbc.Input(
                            id=IDs.INV_CHECKPOINT,
                            type="text",
                            value="",
                            placeholder=(
                                "/path/to/checkpoint.pt"
                            ),
                        ),
                    ),
                ],
            ),
            html.Hr(style={"margin": "12px 0"}),
            # ── Confirm ───────────────────────────
            html.Div(
                id=IDs.INV_STATUS,
                className="am-keys-status",
                style={"minHeight": "20px"},
            ),
            dbc.Button(
                [
                    html.I(
                        className=(
                            "bi bi-check2-circle me-1"
                        )
                    ),
                    "Confirm & Apply",
                ],
                id=IDs.BTN_INV_CONFIRM,
                color="success",
                size="sm",
                className="w-100",
                n_clicks=0,
            ),
            html.Div(
                html.Small(
                    "Parameters are used on the next "
                    "PINN or hybrid inversion request.",
                    style={
                        "color": "var(--subtext0)",
                        "marginTop": "6px",
                    },
                )
            ),
        ],
    )


# ── Full layout ────────────────────────────────────

def create_layout() -> html.Div:
    return html.Div(
        id="am-root",
        className="am-root",
        children=[
            # ── data stores ────────────────────────
            dcc.Store(
                id=IDs.STORE_MESSAGES, data=[]
            ),
            dcc.Store(
                id=IDs.STORE_EDI, data={}
            ),
            dcc.Store(
                id=IDs.STORE_SESSION, data={}
            ),
            dcc.Store(
                id=IDs.STORE_SETTINGS,
                storage_type="local",
                data={
                    "provider": "offline",
                    "export_fmt": "png",
                },
            ),
            dcc.Store(
                id=IDs.STORE_JOB, data={}
            ),
            dcc.Store(
                id=IDs.STORE_THEME,
                storage_type="local",
                data="dark",
            ),
            dcc.Store(
                id=IDs.STORE_FIGS, data={}
            ),
            dcc.Store(
                id=IDs.LOAD_MODE_STORE,
                data="replace",
            ),
            # JS folder loader → Python callback
            dcc.Store(
                id=IDs.FOLDER_STORE, data={}
            ),
            # chat session history (localStorage)
            dcc.Store(
                id=IDs.STORE_HISTORY,
                storage_type="local",
                data=[],
            ),
            # pinned messages (localStorage)
            dcc.Store(
                id=IDs.STORE_PINS,
                storage_type="local",
                data=[],
            ),
            dcc.Store(id=IDs.PIN_SCROLL_DUMMY),
            # PINN / Hybrid inversion config
            dcc.Store(
                id=IDs.STORE_INV_CONFIG,
                storage_type="local",
                data={
                    "mode": "pinn",
                    "dim": 1,
                    "solver": "mt1d",
                    "n_layers": 10,
                    "depth_max": 2000.0,
                    "epochs": 500,
                    "lr": 0.01,
                    "smoothness_weight": 0.01,
                    "lateral_weight": 0.005,
                    "graph_weight": 0.005,
                    "radius": 5000.0,
                    "checkpoint": "",
                },
            ),
            # auto-scroll trigger (dummy target)
            dcc.Store(id="am-scroll-dummy"),
            # pending workflow config
            dcc.Store(
                id=IDs.STORE_PENDING, data={}
            ),
            # post-correction action store
            dcc.Store(
                id=IDs.STORE_POSTPROC, data={}
            ),
            # output folder browser nav state
            dcc.Store(
                id=IDs.OUTPUT_BROWSE_STORE,
                data={},
            ),
            # polling interval
            dcc.Interval(
                id=IDs.INTERVAL_POLL,
                interval=600,
                disabled=True,
                n_intervals=0,
            ),
            # ── overlays / drawers ─────────────────
            _splash(),
            _edi_canvas(),
            _settings_canvas(),
            _inv_params_canvas(),
            _fig_modal(),
            _param_modal(),
            _line_sel_modal(),
            _postproc_modal(),
            _output_browse_modal(),
            # ── main chrome ────────────────────────
            _topbar(),
            html.Div(
                id=IDs.AM_BODY,
                className="am-body",
                children=[
                    _sidebar(),
                    html.Div(
                        id=IDs.AM_MAIN,
                        className="am-main",
                        children=[
                            html.Div(
                                id=IDs.CHAT_WINDOW,
                                children=[
                                    _chat_welcome()
                                ],
                            ),
                            _input_bar(),
                        ],
                    ),
                ],
            ),
        ],
    )
