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

import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html

from pycsamt.app.desktop.agent_registry import (
    agent_names,
)
from pycsamt.app.web.utils import (
    _empty_map,
    _empty_plotly_heatmap,
    empty_src,
)

# ──────────────────────────────────────────────────────────────────────────────
# IDs — centralised to avoid typos across layout + callbacks
# ──────────────────────────────────────────────────────────────────────────────


class IDs:
    # Stores
    STORE_DATA = "store-data"
    STORE_SELECTION = "store-selection"
    STORE_THEME = "store-theme"
    STORE_SIDEBAR = "store-sidebar"
    SESSION_ID = "session-id"
    STORE_ACTIVE_LINES = "store-active-lines"  # {"active":[...], "all":[...]}

    # Navbar
    BTN_LOAD = "btn-load"
    BTN_ADD_LINE = "btn-add-line"  # opens load modal in append mode
    BTN_THEME = "btn-theme"
    BTN_LINES_OPEN = "btn-lines-open"  # opens Lines offcanvas

    # Lines offcanvas
    LINES_OFFCANVAS = "lines-offcanvas"
    LINES_PANEL_BODY = "lines-panel-body"  # dynamic content target
    BTN_LINES_ALL = "btn-lines-all"
    BTN_LINES_NONE = "btn-lines-none"
    LINES_DETECT_MODE = "lines-detect-mode"  # folder | auto | edit
    BTN_LINES_DETECT = "btn-lines-detect"  # run auto-detect
    BTN_LINES_RENAME_APPLY = "btn-lines-rename-apply"  # commit rename inputs

    # Sidebar
    BTN_SIDEBAR_TOGGLE = "btn-sidebar-toggle"
    SIDEBAR_DIV = "sidebar-div"

    # KPI strip (home page)
    KPI_STATIONS = "kpi-stations"
    KPI_FREQ = "kpi-freq"
    KPI_PROFILES = "kpi-profiles"
    KPI_SURVEY = "kpi-survey"

    # Load modal
    MODAL_LOAD = "modal-load"
    UPLOAD = "upload-edi"
    UPLOAD_FILELIST = "upload-filelist"
    INPUT_PATH = "input-data-path"
    BTN_LOAD_CONFIRM = "btn-load-confirm"
    LOAD_FEEDBACK = "load-feedback"
    LOAD_SPINNER = "load-spinner"
    # Folder-browse (JS-driven) — separate from dcc.Upload
    FOLDER_UPLOAD_STORE = "folder-upload-store"
    BTN_BROWSE_FOLDER = "btn-browse-folder"
    LOAD_FOLDER_FILTER_WRAP = "load-folder-filter-wrap"
    LOAD_FOLDER_FILTER = "load-folder-filter"

    # Station panel
    STATION_TABLE = "station-table"
    STATION_SUMMARY = "station-summary"
    STATION_LINE_PILLS = "station-line-pills"
    STATION_LINE_FILTER = "station-line-filter"
    STATION_LIST_BODY = "station-list-body"

    # Map
    MAP_GRAPH = "map-graph"
    MAP_OVERLAY = "map-overlay-select"

    # Profile tabs
    PROFILE_TABS = "profile-tabs"  # dbc.Tabs on home page mini-preview
    PROF_ACTIVE_TAB = (
        "prof-active-tab"  # dcc.Store: active tab on full Profile page
    )
    IMG_RHO_PHI = "img-rho-phi"
    IMG_RHO_PS = "img-rho-ps"
    IMG_PHASE_PS = "img-phase-ps"
    IMG_TIPPER = "img-tipper"
    IMG_PT = "img-pt"
    IMG_PUB = "img-pub"

    # Period / frequency range slider
    FREQ_SLIDER = "freq-slider"

    # QC page v2
    QC_GROUP_SELECT = "qc-group-select"  # kept for compat (hidden dbc.Select)
    QC_GROUP_SEL = "qc-group-sel"  # dcc.Store: active group name
    QC_PLOT_SELECT = "qc-plot-select"  # dbc.RadioItems: plot within group
    BTN_QC_RUN = "btn-qc-run"
    IMG_QC = "img-qc-plot"  # alias → same as IMG_QC_PLOT
    IMG_QC_PLOT = "img-qc-plot"  # main selected-plot image
    IMG_QC_OVERVIEW = "img-qc-overview"  # always-on quicklook image
    QC_ACTIVE_TAB = "qc-active-tab"  # dcc.Store: "plot"|"overview"
    QC_SPINNER = "qc-spinner"
    QC_FEEDBACK = "qc-feedback"
    QC_LINE_SEL = "qc-line-sel"  # Lines multi-dropdown
    QC_STN_SEL = "qc-stn-sel"  # Stations multi-dropdown
    QC_COMP_SEL = "qc-comp-sel"  # Component radio: xy|yx|both
    QC_PARAM_METRIC = "qc-param-metric"  # metric combo
    QC_PARAM_METHOD = "qc-param-method"  # method combo
    QC_PARAM_THRESH = "qc-param-thresh"  # threshold / ci_lo number
    QC_PARAM_NMAX = "qc-param-nmax"  # max stations integer
    QC_PARAM_MAINS = "qc-param-mains"  # mains Hz (50|60)
    QC_FIGSIZE = "qc-figsize"  # figsize radio: compact|wide|tall|pub
    QC_PARAMS_PANEL = "qc-params-panel"  # contextual params div
    QC_AUTO_TOGGLE = "qc-auto-toggle"  # checklist: ["auto"] = auto-render on

    # Agent panel (home)
    AGENT_SELECT = "agent-select"
    BTN_RUN_AGENT = "btn-run-agent"
    AGENT_LOG = "agent-log"
    AGENT_RESULT = "agent-result-img"
    AGENT_SPINNER = "agent-spinner"
    AGENT_SUMMARY = "agent-summary"

    # Error toast
    TOAST_ERROR = "toast-error"
    TOAST_BODY = "toast-body"
    # Success toast
    TOAST_SUCCESS = "toast-success"
    TOAST_SUCCESS_BODY = "toast-success-body"
    # PyTorch install modal
    INV_INSTALL_MODAL = "inv-install-modal"
    BTN_INV_INSTALL_YES = "btn-inv-install-yes"
    BTN_INV_INSTALL_NO = "btn-inv-install-no"
    INV_INSTALL_PROGRESS = "inv-install-progress"
    INV_INSTALL_PROG_AREA = "inv-install-prog-area"
    INV_INSTALL_LOG = "inv-install-log"
    INV_INSTALL_INTERVAL = "inv-install-interval"
    INV_INSTALL_LABEL = "inv-install-label"

    # Status footer
    STATUS_LABEL = "status-label"
    STATUS_FREQ = "status-freq"
    STATUS_STATION = "status-station"

    # Interval (polling)
    INTERVAL = "interval-poll"

    # Section (2D inversion result)
    IMG_SECTION = "img-section"

    # ── Multi-page navigation ────────────────────────────────────────────
    NAV_SECTION = "nav-section"
    CONTENT_AREA = "content-area"

    # ── Advanced plots page ──────────────────────────────────────────────
    ADV_GROUP = "adv-group"  # kept for compat (no longer rendered)
    ADV_PLOT = "adv-plot"
    BTN_ADV_RUN = "btn-adv-run"
    IMG_ADV = "img-adv"  # kept for compat
    ADV_SPINNER = "adv-spinner"
    ADV_ACTIVE_TAB = "adv-active-tab"  # dcc.Store: active tab group
    ADV_INFO_BAR = "adv-info-bar"  # context strip below tab bar
    # Per-group figure IDs (persistent figure per tab)
    IMG_ADV_STRIKE = "img-adv-strike"
    IMG_ADV_PT = "img-adv-pt"
    IMG_ADV_INDUCTION = "img-adv-induction"
    IMG_ADV_IMPEDANCE = "img-adv-impedance"
    IMG_ADV_DEPTH = "img-adv-depth"
    IMG_ADV_SURVEY = "img-adv-survey"
    # Advanced page — station/line filter + parameter panel
    ADV_LINE_FILTER = "adv-line-filter"
    ADV_STN_FILTER = "adv-stn-filter"
    ADV_STN_ALL = "adv-stn-all"
    ADV_STN_NONE = "adv-stn-none"
    ADV_DESCRIPTION = "adv-description"
    ADV_DATA_BAR = "adv-data-bar"
    ADV_BTN_TRAIN = "adv-btn-train"
    ADV_TRAIN_SPINNER = "adv-train-spinner"
    ADV_TRAIN_STATUS = "adv-train-status"
    # PT Strip / PT Strip Grid parameters (plot_phase_tensor_strip[_grid])
    ADV_PT_STATION = "adv-pt-station"  # single-station dropdown
    ADV_PT_PER_LINE = "adv-pt-per-line"  # stations-per-line count (grid)

    # ── TDEM page ────────────────────────────────────────────────────────
    TDEM_CAT = "tdem-cat"  # kept for backwards compat
    TDEM_PLOT = "tdem-plot"
    TDEM_FOLDER = "tdem-folder"
    BTN_TDEM_LOAD = "btn-tdem-load"
    BTN_TDEM_SAMPLE = "btn-tdem-sample"
    BTN_TDEM_RUN = "btn-tdem-run"
    IMG_TDEM = "img-tdem"  # kept for backwards compat
    TDEM_INFO = "tdem-info"
    TDEM_SPINNER = "tdem-spinner"
    STORE_TDEM_LOADED = "store-tdem-loaded"
    # ── TDEM tab system ──────────────────────────────────────────────
    TDEM_ACTIVE_TAB = "tdem-active-tab"
    TDEM_CTX_BAR = "tdem-ctx-bar"
    IMG_TDEM_DECAY = "img-tdem-decay"
    IMG_TDEM_SECTION = "img-tdem-section"
    IMG_TDEM_MAP = "img-tdem-map"
    IMG_TDEM_DASH = "img-tdem-dash"
    TDEM_FIGSIZE = "tdem-figsize"
    TDEM_CMAP = "tdem-cmap"
    TDEM_SOUNDING_IDX = "tdem-sounding-idx"
    TDEM_GATE_IDX = "tdem-gate-idx"
    # ── TDEM folder browse + export ──────────────────────────────────
    BTN_TDEM_BROWSE = "btn-tdem-browse"
    TDEM_EXPORT_FMT = "tdem-export-fmt"
    TDEM_EXPORT_DPI = "tdem-export-dpi"
    BTN_TDEM_EXPORT = "btn-tdem-export"
    TDEM_EXPORT_DL = "tdem-export-dl"

    # ── Correction page ──────────────────────────────────────────────────
    CORR_CAT = "corr-cat"
    CORR_METHOD = "corr-method"
    CORR_STACK = "corr-stack"
    CORR_VIEW = "corr-view"
    BTN_CORR_PREVIEW = "btn-corr-preview"
    BTN_CORR_APPLY = "btn-corr-apply"
    BTN_CORR_UNDO = "btn-corr-undo"
    BTN_CORR_RESET = "btn-corr-reset"
    IMG_CORR_LEFT = "img-corr-left"
    IMG_CORR_RIGHT = "img-corr-right"
    IMG_CORR_OVERLAY = "img-corr-overlay"
    IMG_CORR_DIFF = "img-corr-diff"
    CORR_SPINNER = "corr-spinner"
    CORR_FEEDBACK = "corr-feedback"
    CORR_STORE = "corr-state-store"
    CORR_PARAMS_PANEL = "corr-params-panel"
    CORR_DESC = "corr-method-desc"
    CORR_DIFF_STATS = "corr-diff-stats"  # Δρₐ stats strip on Diff tab
    BTN_CORR_EXPORT = "btn-corr-export"  # export corrected data
    CORR_EXPORT_FMT = "corr-export-fmt"  # edi | csv | xyz
    CORR_DOWNLOAD = "corr-download"  # dcc.Download target
    CORR_SCOPE_FREQ_LO = "corr-scope-freq-lo"  # frequency scope low
    CORR_SCOPE_FREQ_HI = "corr-scope-freq-hi"  # frequency scope high
    # v2 redesign additions
    CORR_ACTIVE_TAB = (
        "corr-active-tab"  # dcc.Store: "ba"|"rho-phi"|"overlay"|"diff"
    )
    IMG_CORR_RHO_PHI = "img-corr-rho-phi"  # per-station ρ/φ grid figure
    CORR_RHO_PHI_MODE = "corr-rho-phi-mode"  # RadioItems: raw|corrected|both
    CORR_RHO_PHI_STN = "corr-rho-phi-stn"  # max station count for ρ/φ grid
    CORR_RHO_PHI_STYLE = (
        "corr-rho-phi-style"  # view style: grid|pairs|freq-bands
    )
    BTN_CORR_DROP_FREQ = (
        "btn-corr-drop-freq"  # quick-apply drop bad frequencies
    )
    CORR_SNR_THRESH = "corr-snr-thresh"  # SNR threshold for drop-freq
    CORR_MIN_FRAC = "corr-min-frac"  # min station fraction for drop-freq
    CORR_LINE_SEL = "corr-line-sel"  # Dropdown: lines to display in ρ/φ view
    CORR_STN_SEL = "corr-stn-sel"  # Dropdown: stations to display in ρ/φ view
    CORR_COMP_SEL = "corr-comp-sel"  # RadioItems: xy | yx | both
    CORR_OVERLAY_STYLE = (
        "corr-overlay-style"  # Select: median-band|spaghetti|per-line
    )
    BTN_CORR_SMOOTH_PREVIEW = "btn-corr-smooth-preview"
    BTN_CORR_SMOOTH_APPLY = "btn-corr-smooth-apply"
    CORR_SMOOTH_DEGREE = "corr-smooth-degree"  # poly degree (0–8)
    CORR_SMOOTH_BLEND = "corr-smooth-blend"  # blend 0–1
    CORR_SMOOTH_COMP = "corr-smooth-comp"  # components combo
    CORR_SMOOTH_ROBUST = "corr-smooth-robust"  # robust checkbox
    # Manual frequency picker
    BTN_CORR_MFREQ_LOAD = (
        "btn-corr-mfreq-load"  # load freq table from session
    )
    BTN_CORR_MFREQ_ALL = "btn-corr-mfreq-all"  # select all
    BTN_CORR_MFREQ_NONE = "btn-corr-mfreq-none"  # deselect all
    BTN_CORR_MFREQ_PREVIEW = "btn-corr-mfreq-preview"
    BTN_CORR_MFREQ_APPLY = "btn-corr-mfreq-apply"
    CORR_MFREQ_LIST = "corr-mfreq-list"  # dcc.Checklist of frequencies
    CORR_MFREQ_STATUS = "corr-mfreq-status"  # status / count text

    # ── Pipeline page ────────────────────────────────────────────────────
    PIPE_STEP = "pipe-step"
    PIPE_STEP_INFO = "pipe-step-info"
    PIPE_STATUS = "pipe-status"
    PIPE_METHOD = "pipe-method"
    BTN_PIPE_RUN = "btn-pipe-run"
    BTN_PIPE_ALL = "btn-pipe-all"
    BTN_PIPE_SKIP = "btn-pipe-skip"
    BTN_PIPE_RESET = "btn-pipe-reset"
    PIPE_LOG = "pipe-log"
    IMG_PIPE = "img-pipe"
    PIPE_SPINNER = "pipe-spinner"
    PIPE_STORE = "pipe-store"  # {step_id: base64_src}
    PIPE_EXPORT_FOLDER = "pipe-export-folder"
    PIPE_EXPORT_CARD = "pipe-export-card"  # shown/hidden based on step
    PIPE_VIEW_STEP = "pipe-view-step"  # preview step selector
    # Export folder browser modal
    PIPE_BROWSE_MODAL = "pipe-browse-modal"
    PIPE_BROWSE_PATH = "pipe-browse-path"  # dcc.Store: current dir path
    PIPE_BROWSE_CWD = "pipe-browse-cwd"  # current path display label
    PIPE_BROWSE_LIST = "pipe-browse-list"  # directory listing body
    BTN_PIPE_BROWSE = "btn-pipe-browse"  # opens the modal
    BTN_PIPE_BROWSE_UP = "btn-pipe-browse-up"  # navigate to parent
    BTN_PIPE_BROWSE_SEL = "btn-pipe-browse-sel"  # confirm selection
    PIPE_BROWSE_NEW = "pipe-browse-new"  # new-folder name input
    BTN_PIPE_BROWSE_MK = "btn-pipe-browse-mk"  # mkdir button
    PIPE_ACTIVE_TAB = "pipe-active-tab"  # dcc.Store: "log"|"preview"|"status"
    PIPE_FEEDBACK = "pipe-feedback"  # text below spinner in run bar

    # ── Forward modelling page ───────────────────────────────────────────
    FWD_DIM = "fwd-dim"  # kept for save_model compat
    FWD_METHOD = "fwd-method"  # MT1D / CSAMT1D / TEM1D
    FWD_PRESET = "fwd-preset"
    BTN_FWD_PRESET = "btn-fwd-preset"
    BTN_FWD_RUN = "btn-fwd-run"
    BTN_FWD_SAVE = "btn-fwd-save"
    FWD_MODEL_NAME = "fwd-model-name"
    IMG_FWD = "img-fwd"
    FWD_SPINNER = "fwd-spinner"
    FWD_FEEDBACK = "fwd-feedback"
    FWD_LAYER_TABLE = "fwd-layer-table"
    BTN_FWD_ADD_LAYER = "btn-fwd-add-layer"
    FWD_HALFSPACE_RHO = "fwd-halfspace-rho"  # halfspace resistivity input
    FWD_FREQ_MIN = "fwd-freq-min"  # log10 Hz
    FWD_FREQ_MAX = "fwd-freq-max"  # log10 Hz
    FWD_N_FREQ = "fwd-n-freq"  # number of frequencies/times
    FWD_OFFSET = "fwd-offset"  # CSAMT source offset [m]
    FWD_LOOP_R = "fwd-loop-r"  # TEM loop radius [m]
    FWD_CSAMT_CARD = "fwd-csamt-card"  # show/hide CSAMT params
    FWD_TEM_CARD = "fwd-tem-card"  # show/hide TEM params
    FWD_FREQ_LABEL = "fwd-freq-label"  # "Frequency [Hz]" vs "Time [s]" label
    FWD_DIMENSION = (
        "fwd-dimension"  # kept for legacy compat (unused in new layout)
    )
    FWD_ACTIVE_TAB = "fwd-active-tab"  # dcc.Store: fwd-tab-1d / 2d / 3d
    FWD_CTX_BAR = "fwd-ctx-bar"  # context info bar in view area
    IMG_FWD_1D = "img-fwd-1d"  # per-tab persistent figures
    IMG_FWD_2D = "img-fwd-2d"
    IMG_FWD_3D = "img-fwd-3d"
    FWD_1D_PANEL = "fwd-1d-panel"  # view panel: 1-D figure
    FWD_2D_PANEL = "fwd-2d-panel"  # view panel: 2D figure
    FWD_3D_PANEL = "fwd-3d-panel"  # view panel: 3D figure
    # 2D MT config
    FWD2_MODEL_TYPE = (
        "fwd2-model-type"  # halfspace | anomaly | from_layers | random
    )
    FWD2_BG_RHO = "fwd2-bg-rho"
    FWD2_NX = "fwd2-nx"
    FWD2_NZ = "fwd2-nz"
    FWD2_X_MAX = "fwd2-x-max"
    FWD2_Z_MAX = "fwd2-z-max"
    FWD2_N_STATIONS = "fwd2-n-stations"
    FWD2_ANOMALY_RHO = "fwd2-anomaly-rho"
    FWD2_AX_LO = "fwd2-ax-lo"
    FWD2_AX_HI = "fwd2-ax-hi"
    FWD2_AZ_LO = "fwd2-az-lo"
    FWD2_AZ_HI = "fwd2-az-hi"
    FWD2_N_LAYERS = "fwd2-n-layers"
    FWD2_SEED = "fwd2-seed"
    FWD2_PLOT_TYPE = "fwd2-plot-type"
    FWD2_ANOMALY_CARD = "fwd2-anomaly-card"
    FWD2_RANDOM_CARD = "fwd2-random-card"
    # 3D MT config
    FWD3_MODEL_TYPE = "fwd3-model-type"  # halfspace | block | random
    FWD3_BG_RHO = "fwd3-bg-rho"
    FWD3_NX = "fwd3-nx"
    FWD3_NY = "fwd3-ny"
    FWD3_NZ = "fwd3-nz"
    FWD3_X_MAX = "fwd3-x-max"
    FWD3_Y_MAX = "fwd3-y-max"
    FWD3_Z_MAX = "fwd3-z-max"
    FWD3_NX_ST = "fwd3-nx-st"
    FWD3_NY_ST = "fwd3-ny-st"
    FWD3_ANOMALY_RHO = "fwd3-anomaly-rho"
    FWD3_AX_LO = "fwd3-ax-lo"
    FWD3_AX_HI = "fwd3-ax-hi"
    FWD3_AY_LO = "fwd3-ay-lo"
    FWD3_AY_HI = "fwd3-ay-hi"
    FWD3_AZ_LO = "fwd3-az-lo"
    FWD3_AZ_HI = "fwd3-az-hi"
    FWD3_N_LAYERS = "fwd3-n-layers"
    FWD3_SEED = "fwd3-seed"
    FWD3_PLOT_TYPE = "fwd3-plot-type"
    FWD3_COMPONENT = "fwd3-component"  # xy | yx
    FWD3_FREQ_IDX = "fwd3-freq-idx"
    FWD3_BLOCK_CARD = "fwd3-block-card"
    FWD3_RANDOM_CARD = "fwd3-random-card"

    # ── Inversion page ───────────────────────────────────────────────────
    INV_DIM = "inv-dim"  # kept for compat
    INV_SOLVER = "inv-solver"  # kept for compat
    INV_DATA_DIR = "inv-data-dir"  # kept for compat
    BTN_INV_RUN = "btn-inv-run"
    INV_SPINNER = "inv-spinner"
    INV_RUN_MSG = "inv-run-msg"
    INV_LOG = "inv-log"
    IMG_INV = "img-inv"
    # Extended inversion controls
    INV_TRACK = "inv-track"  # active tab: inv-tab-trad / inv-tab-ai
    INV_FEEDBACK = "inv-feedback"
    INV_METHOD = "inv-method"  # mt | csamt | tdem
    INV_T_DIM = "inv-t-dim"  # 1d | 2d (traditional)
    INV_BACKEND = "inv-backend"  # builtin | occam2d
    INV_DATA_SRC = "inv-data-src"  # synthetic | session
    INV_TRUE_TABLE = "inv-true-table"
    BTN_INV_ADD_LAYER = "btn-inv-add-layer"
    INV_HALFSPACE_RHO = "inv-halfspace-rho"
    INV_NOISE_LEVEL = "inv-noise-level"
    INV_N_STATIONS = "inv-n-stations"  # synthetic 2D station count
    INV_PROFILE_LEN = "inv-profile-len"  # synthetic 2D profile length [m]
    INV_N_LAYERS = "inv-n-layers"  # layers to recover
    INV_REGULARIZE = "inv-regularize"
    INV_MAX_ITER = "inv-max-iter"
    INV_ERROR_FLOOR = "inv-error-floor"
    INV_INCLUDE_PHASE = "inv-include-phase"
    INV_SYN_PANEL = "inv-syn-panel"  # show/hide synthetic true-model panel
    INV_SESS_PANEL = "inv-sess-panel"  # show/hide session-data info panel
    INV_2D_PARAMS = "inv-2d-params"  # show/hide 2D station/profile params
    INV_FREQ_MIN = "inv-freq-min"
    INV_FREQ_MAX = "inv-freq-max"
    INV_N_FREQ = "inv-n-freq"
    # AI inversion controls
    INV_AI_DIM = "inv-ai-dim"
    INV_AI_ARCH = "inv-ai-arch"
    INV_AI_SOLVER = "inv-ai-solver"
    INV_AI_N_LAYERS = "inv-ai-n-layers"
    INV_AI_N_SAMPLES = "inv-ai-n-samples"
    INV_AI_EPOCHS = "inv-ai-epochs"
    INV_AI_BATCH = "inv-ai-batch"
    INV_AI_LR = "inv-ai-lr"
    INV_AI_NOISE = "inv-ai-noise"
    # 2D / 3D extra controls
    INV_AI_N_COMPONENTS = "inv-ai-n-components"
    INV_AI_N_DEPTH = "inv-ai-n-depth"
    INV_AI_N_STATIONS = "inv-ai-n-stations"  # 2D panel
    INV_AI_N_STATIONS_3D = "inv-ai-n-stations-3d"  # 3D panel
    INV_AI_UNET_DEPTH = "inv-ai-unet-depth"
    INV_AI_2D_PANEL = "inv-ai-2d-panel"
    INV_AI_3D_PANEL = "inv-ai-3d-panel"
    # Inversion output
    INV_RESULT_TABS = "inv-result-tabs"  # kept for legacy compat
    IMG_INV_CONV = "img-inv-conv"
    INV_STATS = "inv-stats"
    # New fixed-tab system
    INV_ACTIVE_TRACK = "inv-active-track"  # dcc.Store: "trad" | "ai"
    INV_ACTIVE_OUT = (
        "inv-active-out"  # dcc.Store: "result"|"conv"|"stats"|"log"
    )
    INV_CTX_BAR = "inv-ctx-bar"  # context info bar (right panel)

    # ── PINN inversion controls ──────────────────────────────────────────
    INV_PINN_DIM = "inv-pinn-dim"
    INV_PINN_DATA_SRC = "inv-pinn-data-src"
    INV_PINN_LINE_SEL = "inv-pinn-line-sel"
    INV_PINN_STATION_SEL = "inv-pinn-station-sel"
    INV_PINN_SOLVER = "inv-pinn-solver"
    INV_PINN_MODE = "inv-pinn-mode"
    INV_PINN_N_LAYERS = "inv-pinn-n-layers"
    INV_PINN_DEPTH_MAX = "inv-pinn-depth-max"
    INV_PINN_N_FREQS = "inv-pinn-n-freqs"
    INV_PINN_EPOCHS = "inv-pinn-epochs"
    INV_PINN_LR = "inv-pinn-lr"
    INV_PINN_LOG_EVERY = "inv-pinn-log-every"
    INV_PINN_DEVICE = "inv-pinn-device"
    INV_PINN_LAM_Z = "inv-pinn-lam-z"
    INV_PINN_LAM_X = "inv-pinn-lam-x"
    INV_PINN_LAM_G = "inv-pinn-lam-g"
    INV_PINN_RADIUS = "inv-pinn-radius"
    INV_PINN_SPC = "inv-pinn-spc"
    INV_PINN_COMP_TE = "inv-pinn-comp-te"
    INV_PINN_COMP_TM = "inv-pinn-comp-tm"
    INV_PINN_2D_PANEL = "inv-pinn-2d-panel"
    INV_PINN_3D_PANEL = "inv-pinn-3d-panel"

    # ── Hybrid inversion controls ────────────────────────────────────────
    INV_HYB_DIM = "inv-hyb-dim"
    INV_HYB_DATA_SRC = "inv-hyb-data-src"
    INV_HYB_LINE_SEL = "inv-hyb-line-sel"
    INV_HYB_STATION_SEL = "inv-hyb-station-sel"
    INV_HYB_CHECKPOINT = "inv-hyb-checkpoint"
    INV_HYB_MODE = "inv-hyb-mode"
    INV_HYB_N_LAYERS = "inv-hyb-n-layers"
    INV_HYB_N_FREQS = "inv-hyb-n-freqs"
    INV_HYB_EPOCHS = "inv-hyb-epochs"
    INV_HYB_LR = "inv-hyb-lr"
    INV_HYB_DEVICE = "inv-hyb-device"
    INV_HYB_LAM_Z = "inv-hyb-lam-z"
    INV_HYB_LAM_X = "inv-hyb-lam-x"
    INV_HYB_LAM_G = "inv-hyb-lam-g"
    INV_HYB_RADIUS = "inv-hyb-radius"
    INV_HYB_COMP_TE = "inv-hyb-comp-te"
    INV_HYB_COMP_TM = "inv-hyb-comp-tm"
    INV_HYB_MODEL_INFO = "inv-hyb-model-info"
    INV_HYB_2D_PANEL = "inv-hyb-2d-panel"
    INV_HYB_3D_PANEL = "inv-hyb-3d-panel"
    BTN_HYB_STAGE1_PREV = "btn-hyb-stage1-prev"

    # ── PINN / Hybrid shared ─────────────────────────────────────────────
    INV_BACKEND_SELECT = "inv-backend-select"
    INV_DATAFIT_LINE = "inv-datafit-line"
    INV_DATAFIT_STA = "inv-datafit-sta"
    IMG_INV_DATA_FIT = "img-inv-data-fit"
    INV_PINN_RESULT_STORE = "inv-pinn-result-store"
    INV_HYB_RESULT_STORE = "inv-hyb-result-store"

    # ── Interpretation page ──────────────────────────────────────────────
    INTERP_CAT = "interp-cat"
    INTERP_PLOT = "interp-plot"
    INTERP_DESC = "interp-desc"
    BTN_INTERP_RUN = "btn-interp-run"
    IMG_INTERP = "img-interp"
    INTERP_SPINNER = "interp-spinner"
    INTERP_DATA_SRC = "interp-data-src"  # raw | corrected
    INTERP_FIGSIZE = "interp-figsize"  # sm/md/lg/wide
    INTERP_CMAP = "interp-cmap"  # colormap name
    INTERP_DEPTH_MAX = "interp-depth-max"  # max depth (section plots)
    INTERP_FREQ_LO = "interp-freq-lo"  # freq band low Hz
    INTERP_FREQ_HI = "interp-freq-hi"  # freq band high Hz
    BTN_INTERP_EXPORT = "btn-interp-export"  # download PNG
    INTERP_DOWNLOAD = "interp-download"  # dcc.Download target
    INTERP_INFO_STRIP = "interp-info-strip"  # status / stats bar
    BTN_INTERP_PREREQ = "btn-interp-prereq"  # run prerequisite estimation
    INTERP_PREREQ_OUT = "interp-prereq-out"  # prereq result message
    INTERP_ACTIVE_CAT = "interp-active-cat"  # dcc.Store: category slug
    INTERP_EXPORT_FMT = "interp-export-fmt"  # export format: png/svg/pdf/eps

    # ── 3D Map page ──────────────────────────────────────────────────────────
    MAP3D_GRAPH = "map3d-graph"  # dcc.Graph Plotly 3D
    MAP3D_MODE = "map3d-mode"  # fence | block | depth (legacy)
    MAP3D_ACTIVE_MODE = (
        "map3d-active-mode"  # dcc.Store: "fence"|"block"|"depth"
    )
    MAP3D_DATA_SRC = "map3d-data-src"  # inversion | pseudo | profiles
    MAP3D_CMAP = "map3d-cmap"
    MAP3D_SCALE = "map3d-scale"  # log | linear
    MAP3D_VMIN = "map3d-vmin"
    MAP3D_VMAX = "map3d-vmax"
    MAP3D_OPACITY = "map3d-opacity"  # 0.0 – 1.0
    MAP3D_CONTOURS = "map3d-contours"  # bool switch
    MAP3D_N_CONTOURS = "map3d-n-contours"
    MAP3D_MODE_PANEL = "map3d-mode-panel"  # dynamic sub-panel
    MAP3D_RHO_LO = "map3d-rho-lo"  # universal ρ filter min (Ω·m)
    MAP3D_RHO_HI = "map3d-rho-hi"  # universal ρ filter max (Ω·m)
    # Fence-specific
    MAP3D_LINE_SPACING = "map3d-line-spacing"
    MAP3D_AZIMUTH = "map3d-azimuth"
    MAP3D_PANEL_THICK = "map3d-panel-thick"
    # Block-specific
    MAP3D_ISOVALUE_LO = "map3d-isovalue-lo"
    MAP3D_ISOVALUE_HI = "map3d-isovalue-hi"
    MAP3D_N_SURFACES = "map3d-n-surfaces"
    MAP3D_CLIP_X = "map3d-clip-x"
    MAP3D_CLIP_Y = "map3d-clip-y"
    MAP3D_CLIP_Z = "map3d-clip-z"
    # Depth-slices-specific
    MAP3D_DEPTH_LO = "map3d-depth-lo"
    MAP3D_DEPTH_HI = "map3d-depth-hi"
    MAP3D_N_SLICES = "map3d-n-slices"
    # Annotations
    MAP3D_SHOW_LABELS = "map3d-show-labels"
    MAP3D_TITLE = "map3d-title"
    # Actions / download
    BTN_MAP3D_RUN = "btn-map3d-run"
    BTN_MAP3D_EXPORT_PNG = "btn-map3d-export-png"
    BTN_MAP3D_EXPORT_HTML = "btn-map3d-export-html"
    MAP3D_DOWNLOAD = "map3d-download"
    MAP3D_INFO = "map3d-info"  # status strip
    MAP3D_SPINNER = "map3d-spinner"
    MAP3D_GRID_STORE = (
        "map3d-grid-store"  # serialised profiles for dynamic re-render
    )
    MAP3D_RANGE_HINT = (
        "map3d-range-hint"  # data-range text shown under iso sliders
    )
    # Topography controls
    MAP3D_TOPO_SRC = (
        "map3d-topo-src"  # "none"|"stations"|"elev-corr"|"elev-raw"|"upload"
    )
    MAP3D_TOPO_OPACITY = (
        "map3d-topo-opacity"  # 0.0 – 1.0 terrain surface opacity
    )
    MAP3D_TOPO_STATIONS = (
        "map3d-topo-stations"  # bool — show station markers on terrain
    )
    MAP3D_TOPO_APPLY = (
        "map3d-topo-apply"  # bool — shift model cells to follow topo
    )
    MAP3D_TOPO_UPLOAD = "map3d-topo-upload"  # dcc.Upload for topo file
    MAP3D_TOPO_UPLOAD_STORE = (
        "map3d-topo-upload-store"  # parsed topo upload records
    )
    MAP3D_TOPO_UPLOAD_INFO = "map3d-topo-upload-info"  # upload status label
    # Station marker style (within topography/annotation controls)
    MAP3D_STA_SYMBOL = "map3d-sta-symbol"  # Plotly 3D marker symbol or "auto"
    MAP3D_STA_SIZE = "map3d-sta-size"  # marker size in px (int)
    MAP3D_STA_COLOR = "map3d-sta-color"  # marker colour (hex)

    # ── Agents full page ─────────────────────────────────────────────────
    AGENTS_CAT = "agents-cat"
    AGENTS_NAME = "agents-name"
    AGENTS_PARAMS = "agents-params"  # kept for compat (unused in layout)
    AGENTS_PARAM_FORM = "agents-param-form"  # live param form container
    AGENTS_DESC = "agents-desc"  # agent info block container
    BTN_AGENTS_RUN = "btn-agents-run"
    AGENTS_OUT = "agents-out"
    IMG_AGENTS = "img-agents"
    AGENTS_SPINNER = "agents-spinner"
    AGENTS_STORE = "agents-store"

    # ── Tools offcanvas ──────────────────────────────────────────────────
    TOOLS_CANVAS = "tools-offcanvas"
    ACTIVE_TOOL = "active-tool-store"
    BTN_TOOLS = "btn-tools-menu"

    # ── Settings offcanvas ────────────────────────────────────────────────
    SETTINGS_CANVAS = "settings-offcanvas"
    BTN_SETTINGS = "btn-settings-menu"
    SETTINGS_STORE = "settings-store"
    BTN_SETTINGS_RESET = "btn-settings-reset"
    BTN_SETTINGS_SAVE = "btn-settings-save"
    BTN_SETTINGS_LOAD = "btn-settings-load"
    SETTINGS_UPLOAD = "settings-profile-upload"
    SETTINGS_FEEDBACK = "settings-feedback"
    # AI / Agents sub-section
    SETTINGS_AI_PROVIDER = "settings-ai-provider-store"
    SETTINGS_AI_KEY = "settings-ai-key"
    SETTINGS_AI_KEY_SHOW = "settings-ai-key-show"
    SETTINGS_AI_MODEL = "settings-ai-model"
    SETTINGS_AI_TEST = "settings-ai-test"
    SETTINGS_AI_STATUS = "settings-ai-status"
    # ── Workflow session save / load ──────────────────────────────────────
    SESSION_CANVAS = "session-offcanvas"
    BTN_SESSION_OPEN = "btn-session-open"
    BTN_SESSION_SAVE = "btn-session-save"
    BTN_SESSION_RESTORE = "btn-session-restore"
    BTN_SESSION_CLEAR = "btn-session-clear"
    SESSION_SNAPSHOT = "session-snapshot"  # localStorage auto-save
    SESSION_DL = "session-download"
    SESSION_UL = "session-upload"
    SESSION_NOTE = "session-note"
    SESSION_FEEDBACK = "session-feedback"
    SESSION_AUTOSAVE = "session-autosave-indicator"

    # ── Help / About modal ────────────────────────────────────────────────
    ABOUT_MODAL = "about-modal"
    BTN_ABOUT = "btn-about"

    # ── Inversion Results Viewer page ─────────────────────────────────────
    INVR_STORE = "invr-store"  # {path, solver, meta}
    INVR_ACTIVE_TAB = "invr-active-tab"  # current view-tab slug
    INVR_PATH_INPUT = "invr-path-input"
    BTN_INVR_LOAD = "btn-invr-load"
    INVR_SOLVER_RADIO = "invr-solver-radio"
    INVR_STATUS = "invr-status"
    INVR_WHICH = "invr-which"  # final|initial|iteration
    INVR_ITERATION = "invr-iteration"
    INVR_RHO_MIN = "invr-rho-min"
    INVR_RHO_MAX = "invr-rho-max"
    INVR_DEPTH_MAX = "invr-depth-max"
    INVR_CMAP = "invr-cmap"
    # Section tab controls
    INVR_SECT_MODE = "invr-sect-mode"  # axial|arbitrary
    INVR_SECT_DIR = "invr-sect-dir"
    INVR_SECT_OFFSET = "invr-sect-offset"
    INVR_SECT_START_LAT = "invr-sect-start-lat"
    INVR_SECT_START_LON = "invr-sect-start-lon"
    INVR_SECT_END_LAT = "invr-sect-end-lat"
    INVR_SECT_END_LON = "invr-sect-end-lon"
    INVR_SECT_N_SAMP = "invr-sect-n-samples"
    INVR_SECT_ORIG_LAT = "invr-sect-orig-lat"
    INVR_SECT_ORIG_LON = "invr-sect-orig-lon"
    INVR_SECT_AX_PANEL = "invr-sect-ax-panel"
    INVR_SECT_ARB_PANEL = "invr-sect-arb-panel"
    # Depth Map tab controls
    INVR_DM_DEPTHS = "invr-dm-depths"
    INVR_DM_NCOLS = "invr-dm-ncols"
    INVR_DM_LAT = "invr-dm-lat"
    INVR_DM_LON = "invr-dm-lon"
    # All Profiles tab controls
    INVR_AP_DIR = "invr-ap-dir"
    INVR_AP_OFFSETS = "invr-ap-offsets"
    INVR_AP_NAMES = "invr-ap-names"
    INVR_AP_N = "invr-ap-n"
    INVR_AP_NCOLS = "invr-ap-ncols"
    # Covariance tab controls
    INVR_COV_SMOOTH = "invr-cov-smooth"
    INVR_COV_CMAP = "invr-cov-cmap"
    # Response tab controls
    INVR_RESP_STATION = "invr-resp-station"
    INVR_RESP_COMP = "invr-resp-comp"
    # Pseudo tab controls
    INVR_PS_COMP = "invr-ps-comp"
    # Stations section
    INVR_SHOW_STA = "invr-show-sta"
    INVR_SHOW_NAMES = "invr-show-names"
    INVR_STA_TOL = "invr-sta-tol"
    # Main view
    INVR_TAB_BAR = "invr-tab-bar"
    INVR_INFO_STRIP = "invr-info-strip"
    INVR_CANVAS = "invr-canvas"
    INVR_CTX_PANEL = "invr-ctx-panel"
    BTN_INVR_GEN = "btn-invr-gen"
    BTN_INVR_EXPORT = "btn-invr-export"
    INVR_DOWNLOAD = "invr-download"

    # ── Home dashboard ────────────────────────────────────────────────────
    DASH_LINE_TABLE = "dash-line-table"
    DASH_QUALITY = "dash-quality-area"
    DASH_LAUNCH_MAP = (
        "dash-launch-map"  # kept for back-compat; no longer in layout
    )
    DASH_LAUNCH_PROFILE = (
        "dash-launch-profile"  # kept for back-compat; no longer in layout
    )
    DASH_CHART_DONUT = "dash-chart-donut"
    DASH_CHART_FREQ = "dash-chart-freq"
    DASH_CHART_RADAR = "dash-chart-radar"
    DASH_CHART_GEO = "dash-chart-geo"
    DASH_CHART_RANK = "dash-chart-rank"
    DASH_CHART_HEATMAP = (
        "dash-chart-heatmap"  # slide 6: line × station coverage grid
    )
    DASH_CHART_VIOLIN = (
        "dash-chart-violin"  # slide 7: per-line freq distribution
    )

    # ── Home-page mini map (separate ID to avoid duplicate-ID bug) ───────────
    HOME_MAP_GRAPH = "home-map-graph"
    HOME_MAP_OVERLAY = "home-map-overlay"

    # ── Map page ──────────────────────────────────────────────────────────
    MAP_PAGE_INFO = "map-page-station-info"
    MAP_PAGE_REFRESH = "map-page-refresh"
    MAP_PAGE_BASEMAP = "map-page-basemap"
    MAP_PAGE_LINE_FILTER = "map-page-line-filter"
    MAP_PAGE_MARKER_SIZE = "map-page-marker-size"
    MAP_PAGE_STATS = "map-page-stats"
    # Map type & data
    MAP_PAGE_MAP_TYPE = "map-page-map-type"
    MAP_PAGE_FREQ_UNIT = "map-page-freq-unit"
    MAP_PAGE_FREQ_SEL = "map-page-freq-sel"
    MAP_PAGE_COMP = "map-page-comp"
    MAP_PAGE_DEPTH_TARGET = "map-page-depth-target"
    MAP_PAGE_DEPTH_INFO = "map-page-depth-info"
    # CRS
    MAP_PAGE_CRS_MODE = "map-page-crs-mode"
    MAP_PAGE_UTM_ZONE = "map-page-utm-zone"
    MAP_PAGE_UTM_HEM = "map-page-utm-hem"
    MAP_PAGE_EPSG = "map-page-epsg"
    MAP_PAGE_CRS_INFO = "map-page-crs-info"
    # Display
    MAP_PAGE_CMAP = "map-page-cmap"
    MAP_PAGE_OPACITY = "map-page-opacity"
    # Overlay toggles
    MAP_PAGE_SHOW_LABELS = "map-page-show-labels"
    MAP_PAGE_SHOW_PROFILES = "map-page-show-profiles"
    # Collapsible section wrappers (clientside show/hide)
    MAP_PAGE_FREQ_SEC = "map-page-freq-section"
    MAP_PAGE_DEPTH_SEC = "map-page-depth-section"
    MAP_PAGE_UTM_SEC = "map-page-utm-section"
    MAP_PAGE_EPSG_SEC = "map-page-epsg-section"
    # Inversion overlay
    MAP_PAGE_INV_SEC = "map-page-inv-section"
    MAP_PAGE_INV_DEPTH_M = "map-page-inv-depth-m"
    MAP_PAGE_INV_LINE_SEL = "map-page-inv-line-sel"
    MAP_PAGE_INV_DEPTHS = (
        "map-page-inv-depths"  # dropdown of available depths
    )
    # Rotation / bearing
    MAP_PAGE_BEARING = "map-page-bearing"

    # ── Recompute offcanvas panel ─────────────────────────────────────────
    BTN_RECOMPUTE_OPEN = "btn-recompute-open"
    RECOMPUTE_CANVAS = "recompute-offcanvas"
    RECOMPUTE_SOURCE = "recompute-source"
    RECOMPUTE_PATH = "recompute-path"
    RECOMPUTE_PATH_SEC = "recompute-path-section"
    RECOMPUTE_NAV_BADGE = "recompute-nav-badge"
    RECOMPUTE_SPIN_ICON = "recompute-spin-icon"
    RECOMPUTE_CONFIRM_SEC = "recompute-confirm-section"
    BTN_RECOMPUTE_YES = "btn-recompute-yes"
    BTN_RECOMPUTE_CANCEL = "btn-recompute-cancel-confirm"
    RECOMPUTE_RESPHASE = "recompute-resphase"
    RECOMPUTE_ROTATE_ANG = "recompute-rotate-angle"
    RECOMPUTE_COMPS = "recompute-comps"
    RECOMPUTE_FMIN = "recompute-fmin"
    RECOMPUTE_FMAX = "recompute-fmax"
    RECOMPUTE_FILL = "recompute-fill"
    BTN_RECOMPUTE = "btn-recompute"
    RECOMPUTE_PROG_SEC = "recompute-progress-section"
    RECOMPUTE_TICKER = "recompute-ticker"
    RECOMPUTE_PROGRESS = "recompute-progress-bar"
    RECOMPUTE_LOG = "recompute-log"
    RECOMPUTE_BADGE = "recompute-badge"
    RECOMPUTE_INTERVAL = "recompute-interval"
    RECOMPUTE_STORE = "recompute-state-store"
    RECOMPUTE_FEEDBACK = "recompute-feedback"

    # ── Map page contour overlay ─────────────────────────────────────────
    MAP_PAGE_CONTOUR_ENABLE = "map-page-contour-enable"
    MAP_PAGE_CONTOUR_SEC = "map-page-contour-section"
    MAP_PAGE_CONTOUR_MODE = "map-page-contour-mode"
    MAP_PAGE_CONTOUR_CMAP = "map-page-contour-cmap"
    MAP_PAGE_CONTOUR_LEVELS = "map-page-contour-levels"
    MAP_PAGE_CONTOUR_OPACITY = "map-page-contour-opacity"
    MAP_PAGE_CONTOUR_INTERP = "map-page-contour-interp"
    MAP_PAGE_CONTOUR_EXTRA = "map-page-contour-extra"
    MAP_PAGE_CONTOUR_SMOOTH = "map-page-contour-smooth"
    MAP_PAGE_CONTOUR_LWIDTH = "map-page-contour-lwidth"
    MAP_PAGE_CONTOUR_OPTS = "map-page-contour-opts"
    MAP_PAGE_CONTOUR_RES = "map-page-contour-res"

    # ── Profile page ──────────────────────────────────────────────────────
    PROFILE_PAGE_STATION = "profile-page-station"
    PROFILE_PAGE_ERRBAR = "profile-page-errbar"
    PROFILE_PAGE_PREV = "profile-page-prev"
    PROFILE_PAGE_NEXT = "profile-page-next"
    PROFILE_PAGE_PHASE_RANGE = "profile-page-phase-range"
    PROF_INFO_BAR = "prof-info-bar"

    # ── Welcome / landing overlay ─────────────────────────────────────────
    WELCOME_OVERLAY = "welcome-overlay"
    WELCOME_CTA = "welcome-cta-btn"
    WELCOME_CARD_LOAD = "welcome-card-load"
    WELCOME_CARD_PIPE = "welcome-card-pipeline"
    WELCOME_CARD_AGENTS = "welcome-card-agents"
    WELCOME_CARD_VIZ = "welcome-card-viz"


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
                        [
                            html.Span("·", className="cmd-bar-sep"),
                            html.Span(subtitle, className="cmd-bar-sub"),
                        ]
                        if subtitle
                        else []
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
    def _card(
        bi_icon: str, label: str, value_id: str, default: str, color: str
    ):
        return html.Div(
            [
                html.I(
                    className=f"bi bi-{bi_icon} kpi-icon",
                    style={"color": f"var(--{color})"},
                ),
                html.Div(
                    [
                        html.Span(
                            default, id=value_id, className="kpi-value"
                        ),
                        html.Span(label, className="kpi-label"),
                    ],
                    className="kpi-text",
                ),
            ],
            className="kpi-card",
        )

    return html.Div(
        [
            _card(
                "broadcast", "Stations", IDs.KPI_STATIONS, "0 sites", "blue"
            ),
            _card("activity", "Frequency", IDs.KPI_FREQ, "—", "sapphire"),
            _card(
                "layout-three-columns",
                "Profiles",
                IDs.KPI_PROFILES,
                "0 profiles",
                "teal",
            ),
            _card(
                "geo-alt-fill", "Survey", IDs.KPI_SURVEY, "No survey", "peach"
            ),
        ],
        className="kpi-strip",
    )


# ──────────────────────────────────────────────────────────────────────────────
# Overlay options (map)
# ──────────────────────────────────────────────────────────────────────────────

_OVERLAY_OPTIONS = [
    {"label": lbl, "value": lbl}
    for lbl in [
        "Index",
        "Apparent Resistivity",
        "Phase",
        "Quality Score",
        "Z-strike",
        "N_freq",
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
        _mi("strike-analyzer", "Strike Analyzer", "strike"),
        _mi("EDI-validator", "EDI Validator", "validator"),
        dbc.DropdownMenuItem(divider=True),
        # ── Conversion & Export ───────────────────────────────
        dbc.DropdownMenuItem("Conversion & Export", header=True),
        _mi("format-converter", "Format Converter", "converter"),
        _mi("batch-export", "Batch Export Plots", "batch"),
        dbc.DropdownMenuItem(divider=True),
        # ── Geospatial ────────────────────────────────────────
        dbc.DropdownMenuItem("Geospatial", header=True),
        _mi("coordinate-transformer", "Coordinate Transformer", "coords"),
        _mi("elevation", "Elevation Enrichment", "elevation"),
        dbc.DropdownMenuItem(divider=True),
        # ── Visualisation ─────────────────────────────────────
        dbc.DropdownMenuItem("Visualisation", header=True),
        _mi("station-response", "Station Response Inspector", "response"),
        _mi("strike-profile", "Strike Profile Viewer", "strike-profile"),
        _mi("phase-tensor", "Phase Tensor Map", "phase-tensor"),
        dbc.DropdownMenuItem(divider=True),
        # ── Classification & Editing ──────────────────────────
        dbc.DropdownMenuItem("Classification & Editing", header=True),
        _mi("dimensionnality", "Dimensionality Classifier", "dimensionality"),
        _mi("frequency-editor", "Frequency Editor", "freq-editor"),
        dbc.DropdownMenuItem(divider=True),
        # ── Modelling ─────────────────────────────────────────
        dbc.DropdownMenuItem("Modelling", header=True),
        _mi("layered-model", "Layered Model Builder", "layered-model"),
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
            [
                _icon("github", size=13, cls="menu-item-icon"),
                " pycsamt on GitHub",
            ],
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
        html.Img(
            src="/icons/pycsamt_logo.svg",
            className="nav-logo",
            style={"height": "26px"},
        ),
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
                [
                    _icon("tools", size=13, cls="menu-item-icon"),
                    html.Span(" Settings", style={"marginLeft": "4px"}),
                ],
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
                [
                    html.I(className="bi bi-cloud-upload-fill"),
                    html.Span(" Load Data", style={"marginLeft": "6px"}),
                ],
                id=IDs.BTN_LOAD,
                className="btn-navbar-load",
                title="Load (replace) survey data",
            ),
            html.Button(
                [
                    html.I(className="bi bi-layer-forward"),
                    html.Span(" +Lines", style={"marginLeft": "5px"}),
                ],
                id=IDs.BTN_ADD_LINE,
                className="btn-navbar-addline ms-1",
                title="Add new lines to current survey",
            ),
            html.Button(
                [
                    html.I(className="bi bi-arrow-repeat"),
                    html.Span(" Recompute", style={"marginLeft": "6px"}),
                ],
                id=IDs.BTN_RECOMPUTE_OPEN,
                className="btn-navbar-recompute ms-1",
                title="Recompute transfer functions",
            ),
            # Shown when active data is recomputed; hidden otherwise
            html.Span(
                [
                    html.I(className="bi bi-check-circle-fill me-1"),
                    "recomputed",
                ],
                id=IDs.RECOMPUTE_NAV_BADGE,
                className="rcp-nav-badge",
                style={"display": "none"},
            ),
            html.Button(
                [
                    html.I(className="bi bi-layers-half"),
                    html.Span(" Lines", style={"marginLeft": "5px"}),
                ],
                id=IDs.BTN_LINES_OPEN,
                className="btn-navbar-lines ms-1",
                title="Manage active survey lines",
            ),
            html.Button(
                [
                    html.I(className="bi bi-floppy me-1"),
                    html.Span("Session", className="btn-session-label"),
                ],
                id=IDs.BTN_SESSION_OPEN,
                className="btn-navbar-session ms-1",
                title="Save / restore workflow session",
                n_clicks=0,
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
    # ── Section A: Browse Folder (PRIMARY — JS-driven, recursive) ──────────
    browse_section = html.Div(
        [
            html.Div(
                [
                    html.I(
                        className="bi bi-folder2-open me-1",
                        style={
                            "color": "var(--blue)",
                            "fontSize": "15px",
                            "verticalAlign": "middle",
                        },
                    ),
                    html.Span(
                        "Browse folder", className="load-section-title"
                    ),
                    html.Span(
                        " · selects all .edi / .avg / .j files recursively",
                        className="load-section-hint",
                    ),
                ],
                className="load-section-head",
            ),
            # Browse button row
            html.Div(
                [
                    dbc.Button(
                        [
                            html.I(className="bi bi-folder2-open me-2"),
                            "Browse Folder",
                        ],
                        id=IDs.BTN_BROWSE_FOLDER,
                        color="primary",
                        outline=True,
                        size="sm",
                        className="w-100",
                        style={"fontWeight": "500"},
                    ),
                    # File count badge (populated by JS)
                    html.Span(
                        id="folder-file-count",
                        style={
                            "display": "none",
                            "marginLeft": "10px",
                            "fontSize": "11px",
                            "padding": "2px 8px",
                            "borderRadius": "10px",
                            "background": "var(--green)",
                            "color": "#1e1e2e",
                            "fontWeight": "600",
                        },
                    ),
                ],
                className="mt-2 d-flex align-items-center flex-wrap gap-2",
            ),
            # Status line (populated by JS)
            html.Div(
                id="folder-browse-status",
                className="small mt-1 folder-browse-status-info",
                style={"minHeight": "18px"},
            ),
            # Hidden text input kept for compatibility (cleared on browse)
            dbc.Input(id=IDs.INPUT_PATH, type="hidden", value=""),
        ],
        className="load-section load-choice-card",
    )

    # ── Section B: Drag-and-drop / individual files (dcc.Upload) ───────────
    # JS intercepts folder drops before dcc.Upload can reject them.
    upload_zone = dcc.Upload(
        id=IDs.UPLOAD,
        className="upload-drop",
        children=html.Div(
            [
                html.I(
                    className="bi bi-cloud-upload",
                    style={
                        "fontSize": "28px",
                        "display": "block",
                        "textAlign": "center",
                        "marginBottom": "6px",
                        "color": "var(--sub0)",
                    },
                ),
                html.Span(
                    "Drag survey folders or files here",
                    style={"fontWeight": "500", "fontSize": "13px"},
                ),
                html.Br(),
                html.Small(
                    "Drag folder or individual files.",
                    className="text-muted",
                ),
                html.Br(),
                html.Small(
                    "Accepts .edi · .avg · .j",
                    className="text-muted",
                    style={"fontSize": "10px"},
                ),
            ]
        ),
        style={"padding": "16px 14px"},
        multiple=True,
        accept=".edi,.EDI,.avg,.AVG,.j,.J",
    )

    upload_section = html.Div(
        [
            html.Div(
                [
                    html.I(
                        className="bi bi-upload me-1",
                        style={
                            "color": "var(--sub0)",
                            "fontSize": "14px",
                            "verticalAlign": "middle",
                        },
                    ),
                    html.Span(
                        id="upload-section-title",
                        children="Or drag & drop files / folders",
                        className="load-section-title",
                    ),
                ],
                className="load-section-head",
            ),
            # Drag zone — hidden once files are selected; file manager takes its place
            html.Div(id="upload-drop-wrap", children=[upload_zone]),
            # Spinner shown by JS while files are being read; hidden when list renders
            html.Div(
                [
                    html.Div(className="file-loader-ring"),
                    html.Span(
                        "Reading files…",
                        className="file-loader-label",
                        id="file-loader-msg",
                    ),
                ],
                id="file-loading-overlay",
                style={"display": "none"},
                className="file-loading-overlay",
            ),
            dcc.Store(id="load-upload-selection", data=[]),
            html.Div(
                id="load-upload-file-manager", className="load-file-manager"
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Span(
                                "Survey folders",
                                className="load-folder-filter-title",
                            ),
                            html.Span(
                                "Select one or many; empty loads all.",
                                className="load-folder-filter-hint",
                            ),
                        ],
                        className="load-folder-filter-head",
                    ),
                    dcc.Checklist(
                        id=IDs.LOAD_FOLDER_FILTER,
                        options=[],
                        value=[],
                        className="load-folder-filter-list",
                        labelClassName="load-folder-filter-option",
                        inputClassName="load-folder-filter-input",
                    ),
                ],
                id=IDs.LOAD_FOLDER_FILTER_WRAP,
                className="load-folder-filter",
                style={"display": "none"},
            ),
            html.Div(
                id=IDs.UPLOAD_FILELIST, className="mt-1 text-muted small"
            ),
        ],
        className="load-section load-choice-card",
    )

    # Mode toggle — shown at the very top of the modal body
    mode_toggle = html.Div(
        [
            html.Div(
                [
                    html.Button(
                        [
                            html.I(className="bi bi-cloud-upload-fill me-1"),
                            "Replace survey",
                        ],
                        id="load-mode-btn-replace",
                        className="load-mode-btn active",
                        n_clicks=0,
                    ),
                    html.Button(
                        [
                            html.I(className="bi bi-layer-forward me-1"),
                            "Add lines",
                        ],
                        id="load-mode-btn-append",
                        className="load-mode-btn",
                        n_clicks=0,
                    ),
                ],
                className="load-mode-toggle",
            ),
            # Hidden store that carries the selected mode into the load callback
            dcc.Store(id="load-mode-store", data="replace"),
            dcc.Store(id="load-source-selection", data="none"),
            # Contextual hint that changes with mode
            html.Div(id="load-mode-hint-text", className="load-mode-hint"),
        ],
        className="mb-3",
    )

    preflight_panel = html.Div(
        [
            html.Div(
                id="load-preflight-preview",
                className="load-preflight-empty",
                children=[
                    html.I(className="bi bi-info-circle me-1"),
                    "Choose a folder or drop files to preview stations and line grouping.",
                ],
            ),
        ],
        className="load-section load-preflight-panel",
    )

    details_accordion = dbc.Accordion(
        [
            dbc.AccordionItem(
                [
                    preflight_panel,
                ],
                title=html.Span(
                    [
                        html.I(className="bi bi-clipboard2-check me-2"),
                        "Preflight details",
                    ]
                ),
                item_id="load-preflight-details",
            ),
            dbc.AccordionItem(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Span(
                                        "Line grouping",
                                        className="load-adv-label",
                                    ),
                                    html.Span(
                                        "Folder name or loaded metadata",
                                        className="load-adv-value",
                                    ),
                                ],
                                className="load-adv-row",
                            ),
                            html.Div(
                                [
                                    html.Span(
                                        "Duplicate stations",
                                        className="load-adv-label",
                                    ),
                                    html.Span(
                                        "New data replaces duplicates on append",
                                        className="load-adv-value",
                                    ),
                                ],
                                className="load-adv-row",
                            ),
                            html.Div(
                                [
                                    html.Span(
                                        "Coordinate priority",
                                        className="load-adv-label",
                                    ),
                                    html.Span(
                                        "EDI header, companion files, computed fallback",
                                        className="load-adv-value",
                                    ),
                                ],
                                className="load-adv-row",
                            ),
                            html.Div(
                                [
                                    html.Span(
                                        "Post-load",
                                        className="load-adv-label",
                                    ),
                                    html.Span(
                                        "Line manager activates newly loaded lines",
                                        className="load-adv-value",
                                    ),
                                ],
                                className="load-adv-row",
                            ),
                        ],
                        className="load-advanced-body",
                    ),
                ],
                title=html.Span(
                    [
                        html.I(className="bi bi-sliders me-2"),
                        "Import policy",
                    ]
                ),
                item_id="load-advanced-policy",
            ),
        ],
        active_item=None,
        className="load-advanced-accordion",
    )

    choose_data = html.Div(
        [
            browse_section,
            upload_section,
        ],
        className="load-choose-grid",
    )

    detected_summary = html.Div(
        id="load-detected-summary",
        className="load-detected-summary",
        children=[
            html.I(className="bi bi-info-circle me-1"),
            "No data selected yet.",
        ],
    )

    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle(
                    html.Span(
                        id="load-modal-title",
                        children=[
                            _icon("load-edis", size=16),
                            " Load Survey Data",
                        ],
                    )
                )
            ),
            dbc.ModalBody(
                [
                    mode_toggle,
                    choose_data,
                    detected_summary,
                    details_accordion,
                    # ── Loading progress bar (shown while load_data runs) ──────
                    html.Div(
                        id="load-progress-wrap",
                        style={"display": "none"},
                        className="load-progress-wrap",
                        children=[
                            html.Div(
                                [
                                    html.Span(
                                        [
                                            html.I(
                                                className="bi bi-cpu-fill me-1"
                                            ),
                                            html.Span(
                                                id="load-prog-label",
                                                children="Processing stations…",
                                            ),
                                        ],
                                        className="load-prog-label",
                                    ),
                                    html.Span(
                                        id="load-prog-sublabel",
                                        className="load-prog-sublabel",
                                    ),
                                ],
                                className="load-prog-header",
                            ),
                            html.Div(
                                className="load-prog-track",
                                children=[
                                    html.Div(
                                        id="load-prog-fill",
                                        className="load-prog-fill",
                                    ),
                                ],
                            ),
                        ],
                    ),
                    html.Div(id=IDs.LOAD_FEEDBACK, className="mt-2 small"),
                    html.Div(id=IDs.LOAD_SPINNER, style={"display": "none"}),
                ]
            ),
            dbc.ModalFooter(
                [
                    dbc.Button(
                        "Close",
                        id="modal-close-btn",
                        color="secondary",
                        outline=True,
                        className="me-auto",
                    ),
                    dbc.Button(
                        [
                            html.I(className="bi bi-cloud-upload me-1"),
                            "Load Survey",
                        ],
                        id=IDs.BTN_LOAD_CONFIRM,
                        color="success",
                        className="load-primary-action",
                    ),
                ]
            ),
        ],
        id=IDs.MODAL_LOAD,
        is_open=False,
        size="lg",
    )


def _station_panel() -> html.Div:
    _hidden_cols = [
        {"name": c, "id": c}
        for c in [
            "ID",
            "Line",
            "Latitude",
            "Longitude",
            "Elevation",
            "N_freq",
        ]
    ]
    return html.Div(
        [
            # ── header ────────────────────────────────────────────
            html.Div(
                [
                    html.Span(
                        [
                            html.Img(src="/icons/station-response.svg"),
                            "Stations",
                        ],
                        className="panel-card-title",
                    ),
                    html.Span(
                        id=IDs.STATION_SUMMARY,
                        className="stat-chip",
                        children="0 sites",
                    ),
                ],
                className="panel-card-header",
            ),
            # ── station loading progress bar ───────────────────────
            html.Div(
                id="sta-load-bar",
                className="sta-load-bar sta-load-bar--hidden",
            ),
            # ── line filter pills ──────────────────────────────────
            html.Div(
                id=IDs.STATION_LINE_PILLS,
                className="sta-line-pills-row",
                children=[],
            ),
            # ── scrollable station list ────────────────────────────
            html.Div(
                [
                    html.Div(
                        id=IDs.STATION_LIST_BODY, className="sta-list-body"
                    ),
                    # Lemniscate spinner — shown while station list is being built
                    html.Div(
                        [
                            html.Div(className="lem-ring"),
                            html.Span(
                                "Building station list…",
                                className="sta-loader-label",
                                id="sta-loader-msg",
                            ),
                        ],
                        id="sta-loading-spinner",
                        className="sta-loading-spinner sta-loading-spinner--hidden",
                    ),
                ],
                className="panel-card-body",
                style={
                    "padding": "0",
                    "overflowY": "auto",
                    "flex": "1",
                    "minHeight": "0",
                    "position": "relative",
                },
            ),
            # ── hidden DataTable (keeps existing callbacks wired) ──
            dash_table.DataTable(
                id=IDs.STATION_TABLE,
                columns=_hidden_cols,
                data=[],
                row_selectable="single",
                style_table={"display": "none", "height": "0"},
            ),
        ],
        className="panel-card accent",
        style={
            "flex": "1",
            "overflow": "hidden",
            "display": "flex",
            "flexDirection": "column",
        },
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
                        id=IDs.HOME_MAP_OVERLAY,
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
                    id=IDs.HOME_MAP_GRAPH,
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
    _img_style = {
        "width": "100%",
        "maxHeight": "260px",
        "objectFit": "contain",
    }

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
        min=-4,
        max=4,
        step=0.25,
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
                        [
                            html.Img(src="/icons/profile-view.svg"),
                            "Profile View",
                        ],
                        className="panel-card-title",
                    ),
                    html.Span(
                        "log₁₀(T)",
                        className="freq-range-label",
                        style={"fontSize": "10px", "whiteSpace": "nowrap"},
                    ),
                    html.Div(
                        _freq_slider, style={"flex": "1", "minWidth": "80px"}
                    ),
                ],
                className="panel-card-header",
            ),
            html.Div(
                dbc.Tabs(
                    [
                        dbc.Tab(
                            _img_tab(IDs.IMG_RHO_PHI),
                            label="ρₐ / φ",
                            tab_id="tab-rho-phi",
                        ),
                        dbc.Tab(
                            _graph_tab(IDs.IMG_RHO_PS),
                            label="ρₐ section",
                            tab_id="tab-rho-ps",
                        ),
                        dbc.Tab(
                            _graph_tab(IDs.IMG_PHASE_PS),
                            label="φ section",
                            tab_id="tab-phase-ps",
                        ),
                        dbc.Tab(
                            _img_tab(IDs.IMG_TIPPER),
                            label="Tipper",
                            tab_id="tab-tipper",
                        ),
                        dbc.Tab(
                            _img_tab(IDs.IMG_PT),
                            label="Phase Tensor",
                            tab_id="tab-pt",
                        ),
                        dbc.Tab(
                            html.Div(
                                html.Img(
                                    id=IDs.IMG_PUB,
                                    src=empty_src(dark=True),
                                    style={
                                        "width": "100%",
                                        "maxHeight": "380px",
                                        "objectFit": "contain",
                                    },
                                ),
                                className="fig-img-wrap",
                                style={"minHeight": "340px"},
                            ),
                            label="➔ Publication",
                            tab_id="tab-pub",
                        ),
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
    """Home right-panel — Survey Insights (replaces the old plain agent runner).

    Three live sections:
      1. Data Health  — auto-refreshed from STORE_DATA (no run needed)
      2. Quick Workflows — 4 icon cards that navigate to Agents page
      3. Last Run — most recent result pulled from AGENTS_STORE

    The old AGENT_SELECT / BTN_RUN_AGENT / AGENT_LOG / AGENT_RESULT /
    AGENT_SUMMARY IDs are kept in hidden elements so the existing
    ``register_agents`` callback continues to function.
    """
    _names = agent_names()

    # ── Quick-workflow cards ───────────────────────────────────────────────
    _WORKFLOWS = [
        (
            "bi-shield-check",
            "QC Analysis",
            "Run quality control\non all stations",
            "qc",
        ),
        (
            "bi-activity",
            "Static Shift",
            "Detect & correct\nstatic shift",
            "static",
        ),
        (
            "bi-grid-3x3",
            "Skew / Dim.",
            "2-D phase-sensitive\nskew map",
            "skew",
        ),
        (
            "bi-layers",
            "Inversion Prep",
            "Build Occam2D input\nfiles for profiles",
            "inversion",
        ),
    ]

    workflow_cards = [
        html.Button(
            [
                html.I(className=f"bi {icon} si-wf-icon"),
                html.Span(label, className="si-wf-label"),
                html.Span(desc, className="si-wf-desc"),
            ],
            id={"type": "si-workflow", "index": key},
            className="si-wf-card",
            n_clicks=0,
            title=desc.replace("\n", " "),
        )
        for icon, label, desc, key in _WORKFLOWS
    ]

    # ── Hidden legacy IDs — keep existing agents callback wired ───────────
    hidden_legacy = html.Div(
        [
            dbc.Select(
                id=IDs.AGENT_SELECT,
                options=[{"label": n, "value": n} for n in _names],
                value=_names[0] if _names else None,
                size="sm",
            ),
            dbc.Button(
                "Run", id=IDs.BTN_RUN_AGENT, color="success", size="sm"
            ),
            dbc.Spinner(
                html.Div(id=IDs.AGENT_SPINNER), size="sm", color="success"
            ),
            html.Pre("", id=IDs.AGENT_LOG),
            html.Img(id=IDs.AGENT_RESULT, src=empty_src(dark=True)),
            html.Div("", id=IDs.AGENT_SUMMARY),
        ],
        style={"display": "none"},
        id="si-legacy-hidden",
    )

    return html.Div(
        [
            # ── Header ────────────────────────────────────────────────────
            html.Div(
                [
                    html.I(className="bi bi-stars me-2 si-header-icon"),
                    html.Span(
                        "Survey Insights", className="panel-card-title"
                    ),
                    html.Button(
                        [
                            html.I(className="bi bi-arrow-right-circle"),
                            " Agents",
                        ],
                        id="si-open-agents-btn",
                        className="si-open-btn",
                        n_clicks=0,
                        title="Open full Agents page",
                    ),
                ],
                className="panel-card-header si-header",
            ),
            html.Div(
                [
                    # ── 1. Data health ─────────────────────────────────────────
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Span(
                                        id="si-health-dot",
                                        className="si-dot si-dot-empty",
                                    ),
                                    html.Span(
                                        id="si-health-text",
                                        children="No survey loaded",
                                        className="si-health-text",
                                    ),
                                ],
                                className="si-health-row",
                            ),
                            html.Div(
                                id="si-health-chips",
                                className="si-health-chips",
                            ),
                        ],
                        className="si-section si-health-section",
                    ),
                    html.Hr(className="si-divider"),
                    # ── 2. Quick workflows ─────────────────────────────────────
                    html.Div(
                        [
                            html.Span(
                                "QUICK WORKFLOWS",
                                className="si-section-label",
                            ),
                            html.Div(workflow_cards, className="si-wf-grid"),
                        ],
                        className="si-section",
                    ),
                    html.Hr(className="si-divider"),
                    # ── 3. Last run ────────────────────────────────────────────
                    html.Div(
                        [
                            html.Span(
                                "LAST RUN", className="si-section-label"
                            ),
                            html.Div(
                                html.Span(
                                    "No runs yet", className="si-no-run"
                                ),
                                id="si-last-run",
                                className="si-last-run-body",
                            ),
                        ],
                        className="si-section",
                    ),
                ],
                className="panel-card-body si-body",
            ),
            hidden_legacy,
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
        style={
            "position": "fixed",
            "top": 70,
            "right": 16,
            "width": 360,
            "zIndex": 1500,
        },
    )


def _success_toast() -> dbc.Toast:
    return dbc.Toast(
        html.Div(id=IDs.TOAST_SUCCESS_BODY),
        id=IDs.TOAST_SUCCESS,
        header="Success",
        is_open=False,
        dismissable=True,
        icon="success",
        duration=6000,
        style={
            "position": "fixed",
            "top": 70,
            "right": 16,
            "width": 360,
            "zIndex": 1501,
        },
    )


def _status_footer() -> html.Div:
    return html.Div(
        [
            html.Span(className="status-dot"),
            html.Span(
                "Ready", id=IDs.STATUS_LABEL, style={"fontSize": "11px"}
            ),
            html.Span("|", className="status-sep"),
            html.Span("", id=IDs.STATUS_STATION, style={"fontSize": "11px"}),
            html.Span("|", className="status-sep"),
            html.Span(
                "",
                id=IDs.STATUS_FREQ,
                className="status-freq-lbl",
                style={"fontSize": "11px"},
            ),
            html.Div(
                [
                    html.I(
                        className="bi bi-circle-fill me-1 status-ver-dot",
                        style={"fontSize": "5px", "verticalAlign": "1px"},
                    ),
                    html.Span(
                        "pyCSAMT v2",
                        className="status-ver-txt",
                        style={"fontSize": "10px"},
                    ),
                ],
                style={
                    "marginLeft": "auto",
                    "display": "flex",
                    "alignItems": "center",
                    "gap": "4px",
                },
            ),
        ],
        className="status-footer",
    )


# ──────────────────────────────────────────────────────────────────────────────
# Navigation entries & sidebar
# ──────────────────────────────────────────────────────────────────────────────

_NAV_ENTRIES = [
    ("home", "summary", "Home"),
    ("map", "map-view", "Map View"),
    ("profile", "profile-view", "Profile View"),
    ("qc", "qc", "Quality Control"),
    ("correction", "sites-correction", "Correction"),
    ("advanced", "advanced-tools", "Advanced Plots"),
    ("tdem", "tdem", "TDEM"),
    ("pipeline", "pipeline", "Pipeline"),
    ("forward", "forward", "Forward Model"),
    ("inversion", "inversion", "Inversion"),
    ("interpretation", "interpret", "Interpretation"),
    ("map3d", "3d", "3D Map"),
    ("inv-results", "results", "Results View"),
    ("agents", "agents", "AI Agents"),
]

_NAV_GROUPS = [
    (None, ["home"]),
    ("Survey", ["map", "profile"]),
    ("Data", ["qc", "correction"]),
    ("Analysis", ["advanced", "tdem"]),
    ("Processing", ["pipeline"]),
    ("Modelling", ["forward", "inversion"]),
    ("Results", ["interpretation", "map3d", "inv-results", "agents"]),
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
                        html.Img(
                            src=f"/icons/{icon_name}.svg",
                            className="nav-btn-icon",
                        ),
                        html.Span(label, className="nav-btn-label"),
                    ],
                    id=f"nav-btn-{page_id}",
                    className="nav-btn"
                    + (" active" if page_id == "home" else ""),
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
            html.I(
                className="bi bi-broadcast",
                style={
                    "fontSize": "2.5rem",
                    "color": "var(--blue)",
                    "marginBottom": "12px",
                },
            ),
            html.H6("No survey data loaded", className="dash-empty-title"),
            html.P(
                "Load EDI files using the Load Data button above "
                "to see survey statistics and launch the map and profile viewers.",
                className="text-muted small text-center",
            ),
        ],
        className="dash-empty-state",
        id="dash-empty-state",
    )

    line_table = html.Div(
        [
            html.Div(
                [
                    html.Span(
                        "Per-Line Summary", className="dash-section-label"
                    ),
                    html.Span(
                        id="dash-survey-name", className="dash-survey-badge"
                    ),
                ],
                className="dash-section-head",
            ),
            html.Div(
                id=IDs.DASH_LINE_TABLE, className="dash-line-table-wrap"
            ),
        ],
        id="dash-data-state",
        style={"display": "none"},
    )

    quality_row = html.Div(
        id=IDs.DASH_QUALITY,
        className="dash-quality-row",
        style={"display": "none"},
    )

    _gcfg = {
        "displayModeBar": False,
        "displaylogo": False,
        "responsive": True,
    }
    _gempty = _empty_plotly_heatmap(dark=True)

    analytics_row = html.Div(
        [
            # Header: title + clickable dot indicators (1 per slide)
            html.Div(
                [
                    html.Span(
                        "Station Analytics", className="dash-section-label"
                    ),
                    html.Div(
                        [
                            html.Span(
                                className="crs-dot",
                                id="crs-dot-0",
                                title="Tipper Coverage",
                            ),
                            html.Span(
                                className="crs-dot",
                                id="crs-dot-1",
                                title="Frequencies / Station",
                            ),
                            html.Span(
                                className="crs-dot",
                                id="crs-dot-2",
                                title="Station Quality Radar",
                            ),
                            html.Span(
                                className="crs-dot",
                                id="crs-dot-3",
                                title="Survey Geometry",
                            ),
                            html.Span(
                                className="crs-dot",
                                id="crs-dot-4",
                                title="Station Ranking",
                            ),
                            html.Span(
                                className="crs-dot",
                                id="crs-dot-5",
                                title="Coverage Heatmap",
                            ),
                            html.Span(
                                className="crs-dot",
                                id="crs-dot-6",
                                title="Per-Line Distribution",
                            ),
                        ],
                        className="crs-dots",
                    ),
                ],
                className="dash-section-head",
                style={"marginBottom": "6px"},
            ),
            # Carousel viewport — 5 slides stacked absolutely
            html.Div(
                [
                    html.Div(
                        dcc.Graph(
                            id=IDs.DASH_CHART_DONUT,
                            figure=_gempty,
                            config=_gcfg,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="crs-slide",
                        id="crs-slide-0",
                    ),
                    html.Div(
                        dcc.Graph(
                            id=IDs.DASH_CHART_FREQ,
                            figure=_gempty,
                            config=_gcfg,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="crs-slide",
                        id="crs-slide-1",
                    ),
                    html.Div(
                        dcc.Graph(
                            id=IDs.DASH_CHART_RADAR,
                            figure=_gempty,
                            config=_gcfg,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="crs-slide",
                        id="crs-slide-2",
                    ),
                    html.Div(
                        dcc.Graph(
                            id=IDs.DASH_CHART_GEO,
                            figure=_gempty,
                            config=_gcfg,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="crs-slide",
                        id="crs-slide-3",
                    ),
                    html.Div(
                        dcc.Graph(
                            id=IDs.DASH_CHART_RANK,
                            figure=_gempty,
                            config=_gcfg,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="crs-slide",
                        id="crs-slide-4",
                    ),
                    html.Div(
                        dcc.Graph(
                            id=IDs.DASH_CHART_HEATMAP,
                            figure=_gempty,
                            config=_gcfg,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="crs-slide",
                        id="crs-slide-5",
                    ),
                    html.Div(
                        dcc.Graph(
                            id=IDs.DASH_CHART_VIOLIN,
                            figure=_gempty,
                            config=_gcfg,
                            style={"height": "100%", "width": "100%"},
                        ),
                        className="crs-slide",
                        id="crs-slide-6",
                    ),
                ],
                className="crs-viewport",
                id="crs-viewport",
            ),
            # Typewriter description text
            html.Div(
                [
                    html.Span(id="crs-typewriter"),
                    html.Span("▌", id="crs-cursor", className="crs-cursor"),
                ],
                className="crs-typewriter",
            ),
            # Slide counter caption
            html.Div(id="crs-label", className="crs-label"),
        ],
        id="dash-launch-row",
        style={"display": "none"},
    )

    return html.Div(
        [
            html.Div(
                [
                    html.Span(
                        [
                            html.I(className="bi bi-speedometer2 me-2"),
                            "Survey Dashboard",
                        ],
                        className="panel-card-title",
                    ),
                    html.Span(
                        id="dash-health-badge", className="dash-health-badge"
                    ),
                ],
                className="panel-card-header",
            ),
            html.Div(
                [empty_state, line_table, quality_row, analytics_row],
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
                        html.Div(_station_panel(), className="home-col"),
                        html.Div(_survey_dashboard(), className="home-col"),
                        html.Div(_agent_panel(), className="home-col"),
                    ],
                    className="home-grid",
                ),
                style={
                    "flex": "1",
                    "overflow": "hidden",
                    "display": "flex",
                    "flexDirection": "column",
                    "minHeight": "0",
                },
            ),
        ],
        id="page-home",
        className="page-content",
        style={"display": "flex"},
    )


_MAP_CMAPS = [
    # All names must be valid Plotly colorscale names (used in Scattermapbox.marker)
    "plasma",
    "viridis",
    "magma",
    "inferno",
    "jet",
    "rdbu_r",
    "balance",
    "spectral",
    "earth",
    "hot",
    "ylorrd",
    "thermal",
    "turbo",
]

_CONTOUR_CMAPS = [
    # Used by matplotlib contourf — any matplotlib name is valid here
    "jet",
    "rainbow",
    "plasma",
    "viridis",
    "magma",
    "RdBu_r",
    "seismic",
    "coolwarm",
    "terrain",
    "hot",
    "YlOrRd",
    "copper",
    "turbo",
    "gnuplot2",
]

_BASEMAP_OPTIONS = [
    # ── Vector tile styles ─────────────────────────────────────
    {"label": "Dark  (Carto)", "value": "carto-darkmatter"},
    {"label": "Light (Carto)", "value": "carto-positron"},
    {"label": "Open Street Map", "value": "open-street-map"},
    {"label": "Stamen Terrain", "value": "stamen-terrain"},
    {"label": "Stamen Toner", "value": "stamen-toner"},
    # ── ESRI raster tiles (free, no token) ────────────────────
    {"label": "ESRI Satellite", "value": "esri-satellite"},
    {"label": "ESRI Topographic", "value": "esri-topo"},
    {"label": "ESRI NatGeo", "value": "esri-natgeo"},
    {"label": "ESRI Ocean", "value": "esri-ocean"},
    {"label": "ESRI Street Map", "value": "esri-street"},
]


def _map_page_content() -> html.Div:
    """Full-width Map View page: left sidebar (accordion panels) + Scattermapbox."""

    def _sep():
        return html.Hr(className="tool-sep")

    def _lbl(text, cls=""):
        return html.Div(text, className=f"ctrl-label {cls}".strip())

    sidebar = html.Div(
        [
            # ── Run bar — FIXED at top ─────────────────────────────────────
            html.Div(
                [
                    html.Div(
                        [
                            html.I(
                                className="bi bi-map-fill me-2",
                                style={
                                    "color": "var(--teal)",
                                    "fontSize": "13px",
                                },
                            ),
                            html.Span("Map View", className="map-run-title"),
                        ],
                        className="d-flex align-items-center mb-2",
                    ),
                    html.Button(
                        [
                            html.I(className="bi bi-house me-2"),
                            "Back to Dashboard",
                        ],
                        id="map-page-home-btn",
                        className="btn btn-sm btn-outline-secondary w-100",
                        n_clicks=0,
                    ),
                ],
                className="map-run-bar",
            ),
            # ── Scrollable controls ────────────────────────────────────────
            html.Div(
                [
                    # Selected station — fixed info card, always visible
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.I(
                                        className="bi bi-geo-alt-fill me-1"
                                    ),
                                    "Selected Station",
                                ],
                                className="ctrl-label",
                            ),
                            html.Div(
                                id=IDs.MAP_PAGE_INFO,
                                className="map-page-station-card",
                                children=[
                                    html.Span(
                                        "Click a station on the map",
                                        className="text-muted small",
                                    )
                                ],
                            ),
                        ],
                        className="ctrl-card mb-1",
                    ),
                    # Accordion parameter sections
                    dbc.Accordion(
                        [
                            # ── Map Type ──────────────────────────────────────────
                            dbc.AccordionItem(
                                [
                                    dbc.Select(
                                        id=IDs.MAP_PAGE_MAP_TYPE,
                                        options=[
                                            {
                                                "label": "Station overview",
                                                "value": "station",
                                            },
                                            {
                                                "label": "Elevation",
                                                "value": "elevation",
                                            },
                                            {
                                                "label": "Apparent Resistivity",
                                                "value": "resistivity",
                                            },
                                            {
                                                "label": "Skin Depth",
                                                "value": "depth",
                                            },
                                            {
                                                "label": "Inversion — Profile",
                                                "value": "inv_profile",
                                            },
                                            {
                                                "label": "Inversion — Depth Slice",
                                                "value": "inv_depth_slice",
                                            },
                                        ],
                                        value="station",
                                        size="sm",
                                        className="mb-2",
                                    ),
                                    _lbl("Colour By"),
                                    dbc.Select(
                                        id=IDs.MAP_OVERLAY,
                                        options=_OVERLAY_OPTIONS,
                                        value="Index",
                                        size="sm",
                                        className="mb-2",
                                    ),
                                    html.Div(
                                        id=IDs.MAP_PAGE_FREQ_SEC,
                                        children=[
                                            html.Hr(className="tool-sep"),
                                            _lbl("Frequency / Period"),
                                            dbc.RadioItems(
                                                id=IDs.MAP_PAGE_FREQ_UNIT,
                                                options=[
                                                    {
                                                        "label": "Hz",
                                                        "value": "hz",
                                                    },
                                                    {
                                                        "label": "Period s",
                                                        "value": "period",
                                                    },
                                                ],
                                                value="hz",
                                                inline=True,
                                                className="mb-1 map-radio-row",
                                            ),
                                            dbc.Select(
                                                id=IDs.MAP_PAGE_FREQ_SEL,
                                                options=[],
                                                placeholder="Load data first…",
                                                size="sm",
                                                className="mb-2",
                                            ),
                                            _lbl("Component"),
                                            dbc.Select(
                                                id=IDs.MAP_PAGE_COMP,
                                                options=[
                                                    {
                                                        "label": "XY (TE mode)",
                                                        "value": "xy",
                                                    },
                                                    {
                                                        "label": "YX (TM mode)",
                                                        "value": "yx",
                                                    },
                                                    {
                                                        "label": "Determinant",
                                                        "value": "det",
                                                    },
                                                ],
                                                value="xy",
                                                size="sm",
                                                className="mb-2",
                                            ),
                                        ],
                                        style={"display": "none"},
                                    ),
                                    html.Div(
                                        id=IDs.MAP_PAGE_DEPTH_SEC,
                                        children=[
                                            _lbl("Target depth (m)"),
                                            dbc.Input(
                                                id=IDs.MAP_PAGE_DEPTH_TARGET,
                                                type="number",
                                                value=500,
                                                min=10,
                                                max=100000,
                                                step=50,
                                                size="sm",
                                                className="mb-1",
                                            ),
                                            html.Div(
                                                id=IDs.MAP_PAGE_DEPTH_INFO,
                                                className="map-crs-info",
                                                children="",
                                            ),
                                        ],
                                        style={"display": "none"},
                                    ),
                                    html.Div(
                                        id=IDs.MAP_PAGE_INV_SEC,
                                        children=[
                                            html.Hr(className="tool-sep"),
                                            _lbl("Inversion line"),
                                            dbc.Select(
                                                id=IDs.MAP_PAGE_INV_LINE_SEL,
                                                options=[],
                                                placeholder="Run inversion first…",
                                                size="sm",
                                                className="mb-2",
                                            ),
                                            html.Div(
                                                id=IDs.MAP_PAGE_INV_DEPTHS,
                                                children=[
                                                    _lbl("Depth level (m)"),
                                                    dbc.Input(
                                                        id=IDs.MAP_PAGE_INV_DEPTH_M,
                                                        type="number",
                                                        value=500,
                                                        min=0,
                                                        max=50000,
                                                        step=50,
                                                        size="sm",
                                                        className="mb-1",
                                                    ),
                                                    html.Div(
                                                        id="map-inv-depth-info",
                                                        className="map-crs-info",
                                                        children="",
                                                    ),
                                                ],
                                            ),
                                        ],
                                        style={"display": "none"},
                                    ),
                                ],
                                title=html.Span(
                                    [
                                        html.I(className="bi bi-layers me-2"),
                                        "Map Type",
                                    ],
                                    className="map3d-settings-title",
                                ),
                                item_id="map-sec-type",
                            ),
                            # ── Display ───────────────────────────────────────────
                            dbc.AccordionItem(
                                [
                                    _lbl("Colormap"),
                                    dbc.Select(
                                        id=IDs.MAP_PAGE_CMAP,
                                        options=[
                                            {"label": c, "value": c}
                                            for c in _MAP_CMAPS
                                        ],
                                        value="plasma",
                                        size="sm",
                                        className="mb-2",
                                    ),
                                    _lbl("Marker size"),
                                    dcc.Slider(
                                        id=IDs.MAP_PAGE_MARKER_SIZE,
                                        min=6,
                                        max=22,
                                        step=2,
                                        value=10,
                                        marks={6: "6", 14: "14", 22: "22"},
                                        className="mb-2",
                                    ),
                                    _lbl("Opacity"),
                                    dcc.Slider(
                                        id=IDs.MAP_PAGE_OPACITY,
                                        min=20,
                                        max=100,
                                        step=5,
                                        value=90,
                                        marks={
                                            20: "20%",
                                            60: "60%",
                                            100: "100%",
                                        },
                                        className="mb-2",
                                    ),
                                    _lbl("Overlay"),
                                    dbc.Checklist(
                                        id=IDs.MAP_PAGE_SHOW_LABELS,
                                        options=[
                                            {
                                                "label": "Station labels",
                                                "value": "labels",
                                            },
                                            {
                                                "label": "Profile lines",
                                                "value": "profiles",
                                            },
                                        ],
                                        value=["labels", "profiles"],
                                        className="mb-2",
                                    ),
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Span(
                                                        [
                                                            html.I(
                                                                className="bi bi-compass me-1"
                                                            ),
                                                            "Rotation",
                                                        ],
                                                        className="ctrl-sublabel",
                                                    ),
                                                    dbc.Button(
                                                        [
                                                            html.I(
                                                                className="bi bi-arrow-up me-1"
                                                            ),
                                                            "N↑",
                                                        ],
                                                        id="map-bearing-reset",
                                                        size="sm",
                                                        color="secondary",
                                                        outline=True,
                                                        className="map-bearing-reset-btn",
                                                        n_clicks=0,
                                                    ),
                                                ],
                                                className="d-flex align-items-center justify-content-between mb-1",
                                            ),
                                            dcc.Slider(
                                                id=IDs.MAP_PAGE_BEARING,
                                                min=0,
                                                max=359,
                                                step=1,
                                                value=0,
                                                marks={
                                                    0: "0°",
                                                    90: "90°",
                                                    180: "180°",
                                                    270: "270°",
                                                },
                                                tooltip={
                                                    "placement": "bottom",
                                                    "always_visible": False,
                                                },
                                                className="mb-1",
                                            ),
                                        ],
                                        className="mb-1",
                                    ),
                                ],
                                title=html.Span(
                                    [
                                        html.I(
                                            className="bi bi-sliders me-2"
                                        ),
                                        "Display",
                                    ],
                                    className="map3d-settings-title",
                                ),
                                item_id="map-sec-display",
                            ),
                            # ── Coordinate System ─────────────────────────────────
                            dbc.AccordionItem(
                                [
                                    dbc.Select(
                                        id=IDs.MAP_PAGE_CRS_MODE,
                                        options=[
                                            {
                                                "label": "Geographic (lat / lon)  EPSG:4326",
                                                "value": "geo",
                                            },
                                            {
                                                "label": "UTM Zone",
                                                "value": "utm",
                                            },
                                            {
                                                "label": "Custom EPSG",
                                                "value": "custom",
                                            },
                                        ],
                                        value="geo",
                                        size="sm",
                                        className="mb-1",
                                    ),
                                    html.Div(
                                        id=IDs.MAP_PAGE_UTM_SEC,
                                        children=html.Div(
                                            [
                                                dbc.InputGroup(
                                                    [
                                                        dbc.InputGroupText(
                                                            "Zone"
                                                        ),
                                                        dbc.Input(
                                                            id=IDs.MAP_PAGE_UTM_ZONE,
                                                            type="number",
                                                            value=50,
                                                            min=1,
                                                            max=60,
                                                            step=1,
                                                            size="sm",
                                                        ),
                                                        dbc.Select(
                                                            id=IDs.MAP_PAGE_UTM_HEM,
                                                            options=[
                                                                {
                                                                    "label": "N",
                                                                    "value": "N",
                                                                },
                                                                {
                                                                    "label": "S",
                                                                    "value": "S",
                                                                },
                                                            ],
                                                            value="N",
                                                            size="sm",
                                                        ),
                                                    ],
                                                    size="sm",
                                                    className="mb-1",
                                                ),
                                            ]
                                        ),
                                        style={"display": "none"},
                                    ),
                                    html.Div(
                                        id=IDs.MAP_PAGE_EPSG_SEC,
                                        children=dbc.InputGroup(
                                            [
                                                dbc.InputGroupText("EPSG"),
                                                dbc.Input(
                                                    id=IDs.MAP_PAGE_EPSG,
                                                    type="text",
                                                    value="4326",
                                                    placeholder="e.g. 32650",
                                                    size="sm",
                                                ),
                                            ],
                                            size="sm",
                                            className="mb-1",
                                        ),
                                        style={"display": "none"},
                                    ),
                                    html.Div(
                                        id=IDs.MAP_PAGE_CRS_INFO,
                                        className="map-crs-info mb-1",
                                    ),
                                ],
                                title=html.Span(
                                    [
                                        html.I(className="bi bi-globe me-2"),
                                        "Coordinate System",
                                    ],
                                    className="map3d-settings-title",
                                ),
                                item_id="map-sec-crs",
                            ),
                            # ── Basemap ───────────────────────────────────────────
                            dbc.AccordionItem(
                                [
                                    dbc.Select(
                                        id=IDs.MAP_PAGE_BASEMAP,
                                        options=_BASEMAP_OPTIONS,
                                        value="carto-darkmatter",
                                        size="sm",
                                    ),
                                ],
                                title=html.Span(
                                    [
                                        html.I(className="bi bi-map me-2"),
                                        "Basemap",
                                    ],
                                    className="map3d-settings-title",
                                ),
                                item_id="map-sec-basemap",
                            ),
                            # ── Contour Overlay ───────────────────────────────────
                            dbc.AccordionItem(
                                [
                                    html.Div(
                                        [
                                            dbc.Checklist(
                                                id=IDs.MAP_PAGE_CONTOUR_ENABLE,
                                                options=[
                                                    {
                                                        "label": "Enable contour overlay",
                                                        "value": "on",
                                                    }
                                                ],
                                                value=[],
                                                switch=True,
                                                className="contour-toggle-check",
                                            ),
                                        ],
                                        className="d-flex align-items-center justify-content-between mb-1",
                                    ),
                                    html.Div(
                                        [
                                            html.Span(
                                                "Interpolates station values onto a grid and overlays "
                                                "filled / line contours on the basemap — similar to Surfer.",
                                                className="text-muted",
                                                style={
                                                    "fontSize": "10px",
                                                    "lineHeight": "1.45",
                                                },
                                            ),
                                        ],
                                        className="mb-1",
                                    ),
                                    html.Div(
                                        id=IDs.MAP_PAGE_CONTOUR_SEC,
                                        children=[
                                            _lbl("Display mode"),
                                            dbc.Select(
                                                id=IDs.MAP_PAGE_CONTOUR_MODE,
                                                options=[
                                                    {
                                                        "label": "Filled + Lines  (Surfer-style)",
                                                        "value": "filled+lines",
                                                    },
                                                    {
                                                        "label": "Filled only",
                                                        "value": "filled",
                                                    },
                                                    {
                                                        "label": "Lines only",
                                                        "value": "lines",
                                                    },
                                                ],
                                                value="filled+lines",
                                                size="sm",
                                                className="mb-2",
                                            ),
                                            _lbl("Colormap"),
                                            dbc.Select(
                                                id=IDs.MAP_PAGE_CONTOUR_CMAP,
                                                options=[
                                                    {"label": c, "value": c}
                                                    for c in _CONTOUR_CMAPS
                                                ],
                                                value="jet",
                                                size="sm",
                                                className="mb-2",
                                            ),
                                            html.Div(
                                                [
                                                    html.Div(
                                                        [
                                                            _lbl("Levels"),
                                                            dbc.Input(
                                                                id=IDs.MAP_PAGE_CONTOUR_LEVELS,
                                                                type="number",
                                                                value=12,
                                                                min=4,
                                                                max=30,
                                                                step=1,
                                                                size="sm",
                                                            ),
                                                        ],
                                                        style={"flex": "1"},
                                                    ),
                                                    html.Div(
                                                        [
                                                            _lbl("Opacity %"),
                                                            dbc.Input(
                                                                id=IDs.MAP_PAGE_CONTOUR_OPACITY,
                                                                type="number",
                                                                value=60,
                                                                min=10,
                                                                max=90,
                                                                step=5,
                                                                size="sm",
                                                            ),
                                                        ],
                                                        style={"flex": "1"},
                                                    ),
                                                ],
                                                className="d-flex gap-2 mb-2",
                                            ),
                                            html.Div(
                                                [
                                                    html.Div(
                                                        [
                                                            _lbl(
                                                                "Interpolation"
                                                            ),
                                                            dbc.Select(
                                                                id=IDs.MAP_PAGE_CONTOUR_INTERP,
                                                                options=[
                                                                    {
                                                                        "label": "Cubic",
                                                                        "value": "cubic",
                                                                    },
                                                                    {
                                                                        "label": "Linear",
                                                                        "value": "linear",
                                                                    },
                                                                    {
                                                                        "label": "Nearest",
                                                                        "value": "nearest",
                                                                    },
                                                                ],
                                                                value="cubic",
                                                                size="sm",
                                                            ),
                                                        ],
                                                        style={"flex": "1"},
                                                    ),
                                                    html.Div(
                                                        [
                                                            _lbl(
                                                                "Extrapolate"
                                                            ),
                                                            dbc.Select(
                                                                id=IDs.MAP_PAGE_CONTOUR_EXTRA,
                                                                options=[
                                                                    {
                                                                        "label": "None",
                                                                        "value": "0.0",
                                                                    },
                                                                    {
                                                                        "label": "Light 12%",
                                                                        "value": "0.12",
                                                                    },
                                                                    {
                                                                        "label": "Medium 28%",
                                                                        "value": "0.28",
                                                                    },
                                                                ],
                                                                value="0.12",
                                                                size="sm",
                                                            ),
                                                        ],
                                                        style={"flex": "1"},
                                                    ),
                                                ],
                                                className="d-flex gap-2 mb-2",
                                            ),
                                            html.Div(
                                                [
                                                    html.Div(
                                                        [
                                                            _lbl(
                                                                "Smoothing σ"
                                                            ),
                                                            dbc.Input(
                                                                id=IDs.MAP_PAGE_CONTOUR_SMOOTH,
                                                                type="number",
                                                                value=1.0,
                                                                min=0.0,
                                                                max=5.0,
                                                                step=0.5,
                                                                size="sm",
                                                            ),
                                                        ],
                                                        style={"flex": "1"},
                                                    ),
                                                    html.Div(
                                                        [
                                                            _lbl(
                                                                "Line width"
                                                            ),
                                                            dbc.Input(
                                                                id=IDs.MAP_PAGE_CONTOUR_LWIDTH,
                                                                type="number",
                                                                value=1.0,
                                                                min=0.3,
                                                                max=4.0,
                                                                step=0.3,
                                                                size="sm",
                                                            ),
                                                        ],
                                                        style={"flex": "1"},
                                                    ),
                                                ],
                                                className="d-flex gap-2 mb-2",
                                            ),
                                            _lbl("Grid resolution"),
                                            dbc.Select(
                                                id=IDs.MAP_PAGE_CONTOUR_RES,
                                                options=[
                                                    {
                                                        "label": "Coarse (100 × 100)",
                                                        "value": "100",
                                                    },
                                                    {
                                                        "label": "Normal (150 × 150)",
                                                        "value": "150",
                                                    },
                                                    {
                                                        "label": "Fine   (200 × 200)",
                                                        "value": "200",
                                                    },
                                                    {
                                                        "label": "Ultra  (300 × 300)",
                                                        "value": "300",
                                                    },
                                                ],
                                                value="150",
                                                size="sm",
                                                className="mb-2",
                                            ),
                                            dbc.Checklist(
                                                id=IDs.MAP_PAGE_CONTOUR_OPTS,
                                                options=[
                                                    {
                                                        "label": "Show level labels",
                                                        "value": "labels",
                                                    },
                                                    {
                                                        "label": "Log scale (ρ / depth)",
                                                        "value": "log",
                                                    },
                                                ],
                                                value=[],
                                                className="mb-1 contour-opts-check",
                                            ),
                                        ],
                                        style={"display": "none"},
                                    ),
                                ],
                                title=html.Span(
                                    [
                                        html.I(
                                            className="bi bi-bezier2 me-2"
                                        ),
                                        "Contour Overlay",
                                    ],
                                    className="map3d-settings-title",
                                ),
                                item_id="map-sec-contour",
                            ),
                            # ── Survey Lines ──────────────────────────────────────
                            dbc.AccordionItem(
                                [
                                    html.Div(
                                        [
                                            html.Span(
                                                "Filter visible profiles",
                                                className="text-muted",
                                                style={"fontSize": "10px"},
                                            ),
                                            html.Span(
                                                "(all)",
                                                id="map-line-filter-count",
                                                className="text-muted",
                                                style={"fontSize": "10px"},
                                            ),
                                        ],
                                        className="d-flex justify-content-between align-items-center mb-1",
                                    ),
                                    dbc.Checklist(
                                        id=IDs.MAP_PAGE_LINE_FILTER,
                                        options=[],
                                        value=[],
                                        className="mb-1 map-line-filter",
                                    ),
                                ],
                                title=html.Span(
                                    [
                                        html.I(
                                            className="bi bi-list-ul me-2"
                                        ),
                                        "Survey Lines",
                                    ],
                                    className="map3d-settings-title",
                                ),
                                item_id="map-sec-lines",
                            ),
                        ],
                        active_item=["map-sec-type"],
                        always_open=True,
                        className="map3d-settings-accordion",
                    ),
                ],
                className="map-ctrl-scroll",
            ),
        ],
        className="page-sidebar map-sidebar",
    )

    def _tb_sep():
        return html.Div(className="map-tb-sep")

    def _tb_btn(icon, label, btn_id, title="", active=False):
        cls = "map-tb-btn active" if active else "map-tb-btn"
        inner = [html.I(className=f"bi {icon}")]
        if label:
            inner += [html.Span(label, className="map-tb-label")]
        return html.Button(
            inner, id=btn_id, n_clicks=0, className=cls, title=title
        )

    map_area = html.Div(
        [
            # ── Fixed toolbar ──────────────────────────────────────────────
            html.Div(
                [
                    # ── Left: data-summary chip ────────────────────────────
                    html.Div(
                        id="map-tb-info",
                        className="map-tb-info",
                        children=html.Span(
                            "No data loaded",
                            className="text-muted",
                            style={"fontSize": "11px"},
                        ),
                    ),
                    _tb_sep(),
                    # ── View controls ──────────────────────────────────────
                    _tb_btn(
                        "bi-arrows-fullscreen",
                        "Fit",
                        "map-tb-fit",
                        "Zoom to fit all stations",
                    ),
                    _tb_btn(
                        "bi-crosshair",
                        "",
                        "map-tb-center",
                        "Re-centre on data (preserve zoom)",
                    ),
                    _tb_sep(),
                    # ── Layer toggles ──────────────────────────────────────
                    _tb_btn(
                        "bi-card-text",
                        "Labels",
                        "map-tb-labels",
                        "Toggle station labels",
                    ),
                    _tb_btn(
                        "bi-bezier",
                        "Profiles",
                        "map-tb-profiles",
                        "Toggle survey-line polylines",
                    ),
                    _tb_btn(
                        "bi-layers-half",
                        "Contour",
                        "map-tb-contour",
                        "Toggle contour overlay",
                    ),
                    _tb_sep(),
                    # ── Basemap quick-switch ───────────────────────────────
                    html.Span("Basemap", className="map-tb-grp-label"),
                    _tb_btn(
                        "bi-moon-stars-fill",
                        "",
                        "map-tb-bm-dark",
                        "Dark  (Carto)",
                    ),
                    _tb_btn(
                        "bi-sun-fill", "", "map-tb-bm-light", "Light (Carto)"
                    ),
                    _tb_btn(
                        "bi-globe-americas",
                        "",
                        "map-tb-bm-sat",
                        "Satellite (ESRI)",
                    ),
                    _tb_btn(
                        "bi-map-fill",
                        "",
                        "map-tb-bm-street",
                        "Street map (OSM)",
                    ),
                    _tb_btn(
                        "bi-mountain",
                        "",
                        "map-tb-bm-topo",
                        "Topographic (ESRI)",
                    ),
                    _tb_sep(),
                    # ── Marker-size stepper ────────────────────────────────
                    html.Span("Markers", className="map-tb-grp-label"),
                    html.Button(
                        html.I(className="bi bi-dash"),
                        id="map-tb-sz-down",
                        n_clicks=0,
                        className="map-tb-btn",
                        title="Decrease marker size",
                    ),
                    html.Div(
                        id="map-tb-sz-val",
                        className="map-tb-sz-val",
                        children="10",
                    ),
                    html.Button(
                        html.I(className="bi bi-plus"),
                        id="map-tb-sz-up",
                        n_clicks=0,
                        className="map-tb-btn",
                        title="Increase marker size",
                    ),
                    # ── Spacer pushes right group to the end ───────────────
                    html.Div(style={"flex": "1"}),
                    # ── Right: export + refresh ────────────────────────────
                    _tb_btn(
                        "bi-image",
                        "PNG",
                        "map-tb-export",
                        "Export map as PNG",
                    ),
                    dbc.Button(
                        [
                            html.I(className="bi bi-arrow-clockwise"),
                            html.Span("Refresh", className="map-tb-label"),
                        ],
                        id=IDs.MAP_PAGE_REFRESH,
                        color="primary",
                        size="sm",
                        n_clicks=0,
                        title="Force full redraw (auto-updates on any control change)",
                        className="map-tb-btn map-tb-refresh",
                    ),
                ],
                className="map-topbar",
            ),
            dcc.Graph(
                id=IDs.MAP_GRAPH,
                figure=_empty_map(dark=True),
                responsive=True,
                config={
                    "displayModeBar": True,
                    "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                    "scrollZoom": True,
                    "displaylogo": False,
                    "toImageButtonOptions": {
                        "format": "png",
                        "filename": "pycsamt_map",
                        "scale": 2,
                    },
                },
                style={"flex": "1", "minHeight": "0", "width": "100%"},
                clear_on_unhover=True,
            ),
            html.Div(id=IDs.MAP_PAGE_STATS, className="map-stats-strip"),
            # Hidden anchor — clientside resize trigger writes here
            html.Div(id="map-resize-anchor", style={"display": "none"}),
        ],
        className="page-main-area",
        style={"display": "flex", "flexDirection": "column"},
    )

    return html.Div(
        [sidebar, map_area],
        className="page-two-zone",
    )


def _profile_page_content() -> html.Div:
    """Profile View — fixed sub-tab bar at top + fixed view panel below."""
    _img_style = {"width": "100%", "height": "100%", "objectFit": "contain"}

    # ── Tab definitions ────────────────────────────────────────────────────────
    _TABS = [
        ("rho-phi", "ρₐ / φ", "bi-activity"),
        ("rho-ps", "ρₐ Section", "bi-grid-1x2"),
        ("phi-ps", "φ Section", "bi-grid-3x2-gap"),
        ("pt", "Phase Tensor", "bi-hexagon"),
        ("tipper", "Tipper", "bi-arrows-expand"),
        ("section", "2D Section", "bi-layers"),
        ("pub", "→ Publication", "bi-journal-richtext"),
    ]
    _DEFAULT = "rho-ps"

    _freq_slider = dcc.RangeSlider(
        id=IDs.FREQ_SLIDER,
        min=-4,
        max=4,
        step=0.25,
        value=[-4, 4],
        marks=None,
        tooltip={"placement": "bottom", "always_visible": False},
        disabled=True,
        className="px-1",
    )

    # ── Left sidebar ──────────────────────────────────────────────────────────
    sidebar = html.Div(
        [
            # ── Run bar — FIXED at top ─────────────────────────────────────
            html.Div(
                [
                    dbc.Button(
                        [
                            html.I(className="bi bi-arrow-repeat me-1"),
                            "Refresh",
                        ],
                        id="profile-page-refresh",
                        color="primary",
                        size="sm",
                        className="w-100 mb-1",
                        n_clicks=0,
                    ),
                    html.Button(
                        [
                            html.I(className="bi bi-house me-2"),
                            "Back to Dashboard",
                        ],
                        id="profile-page-home-btn",
                        className="btn btn-sm btn-outline-secondary w-100",
                        n_clicks=0,
                    ),
                ],
                className="prof-run-bar",
            ),
            # ── Scrollable controls ────────────────────────────────────────
            html.Div(
                [
                    html.Div("Station", className="ctrl-label mt-1"),
                    html.Div(
                        [
                            html.Button(
                                html.I(className="bi bi-chevron-left"),
                                id=IDs.PROFILE_PAGE_PREV,
                                className="btn btn-sm btn-outline-secondary profile-nav-btn",
                                n_clicks=0,
                                title="Previous station",
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
                                n_clicks=0,
                                title="Next station",
                            ),
                        ],
                        className="d-flex align-items-center gap-1 mb-3",
                    ),
                    html.Hr(className="tool-sep"),
                    html.Div(
                        [
                            html.Span("Period range", className="ctrl-label"),
                            html.Span(
                                "log₁₀(T)",
                                className="text-muted",
                                style={"fontSize": "10px"},
                            ),
                        ],
                        className="d-flex align-items-center justify-content-between mb-1",
                    ),
                    _freq_slider,
                    html.Hr(className="tool-sep"),
                    html.Div("Components", className="ctrl-label"),
                    dbc.Checklist(
                        id="profile-page-comps",
                        options=[
                            {"label": "Zxy (TE)", "value": "xy"},
                            {"label": "Zyx (TM)", "value": "yx"},
                            {"label": "Zxx", "value": "xx"},
                            {"label": "Zyy", "value": "yy"},
                        ],
                        value=["xy", "yx"],
                        className="mb-2",
                    ),
                    html.Div("Phase range", className="ctrl-label mt-1"),
                    dbc.Select(
                        id=IDs.PROFILE_PAGE_PHASE_RANGE,
                        options=[
                            {"label": "Auto (data range)", "value": "auto"},
                            {"label": "0° to 90°", "value": "0_90"},
                            {"label": "−90° to 90°", "value": "-90_90"},
                            {"label": "−180° to 180°", "value": "-180_180"},
                            {"label": "−360° to 360°", "value": "-360_360"},
                        ],
                        value="auto",
                        size="sm",
                        className="mb-2",
                    ),
                    dbc.Switch(
                        id=IDs.PROFILE_PAGE_ERRBAR,
                        label="Error bars",
                        value=False,
                        className="mb-2",
                    ),
                    html.Hr(className="tool-sep"),
                    dbc.Button(
                        [
                            html.I(className="bi bi-download me-1"),
                            "Export Tab",
                        ],
                        id="profile-page-export",
                        color="secondary",
                        size="sm",
                        className="w-100",
                        n_clicks=0,
                    ),
                    dcc.Download(id="profile-page-download"),
                ],
                className="prof-ctrl-scroll",
            ),
        ],
        className="page-sidebar prof-sidebar",
    )

    # ── Fixed sub-tab bar ─────────────────────────────────────────────────────
    tab_bar = html.Div(
        [
            html.Button(
                [html.I(className=f"bi {icon} me-1"), label],
                id=f"prof-tab-{tid}",
                className=f"prof-tab-btn{' active' if tid == _DEFAULT else ''}",
                n_clicks=0,
                title=label,
            )
            for tid, label, icon in _TABS
        ],
        className="prof-tab-bar",
    )

    # ── Station info header bar ───────────────────────────────────────────────
    info_bar = html.Div(
        id=IDs.PROF_INFO_BAR,
        children=[
            html.Span(
                [
                    html.I(
                        className="bi bi-broadcast me-1",
                        style={"fontSize": "10px"},
                    ),
                    "No station selected — load data and choose a station",
                ],
                className="prof-info-station",
                style={"fontWeight": "400", "color": "var(--sub0)"},
            ),
        ],
        className="prof-info-bar",
    )

    # ── Shared graph config ───────────────────────────────────────────────────
    _GCfg = {
        "displayModeBar": True,
        "modeBarButtonsToRemove": ["lasso2d", "select2d"],
        "scrollZoom": True,
        "displaylogo": False,
        "toImageButtonOptions": {
            "format": "png",
            "filename": "pycsamt_section",
            "scale": 2,
        },
    }

    def _img_panel(panel_id, img_id, visible=False):
        return html.Div(
            html.Div(
                html.Img(
                    id=img_id, src=empty_src(dark=True), style=_img_style
                ),
                className="fig-img-wrap profile-page-fig",
            ),
            id=panel_id,
            className="prof-panel",
            style={"display": "flex"} if visible else {"display": "none"},
        )

    def _ps_panel(panel_id, graph_id, header_id, visible=False):
        return html.Div(
            [
                html.Div(id=header_id, className="prof-ps-comp-bar"),
                dcc.Graph(
                    id=graph_id,
                    figure=_empty_plotly_heatmap(dark=True),
                    config=_GCfg,
                    style={
                        "flex": "1",
                        "minHeight": "0",
                        "width": "100%",
                        "height": "100%",
                    },
                ),
            ],
            id=panel_id,
            className="prof-panel",
            style={"display": "flex"} if visible else {"display": "none"},
        )

    # ── All content panels (always in DOM, show/hide via CSS) ─────────────────
    panels = html.Div(
        [
            # 0  ρₐ / φ curves
            _img_panel("prof-panel-rho-phi", IDs.IMG_RHO_PHI),
            # 1  ρₐ pseudosection  (default visible)
            _ps_panel(
                "prof-panel-rho-ps",
                IDs.IMG_RHO_PS,
                "prof-ps-rho-header",
                visible=True,
            ),
            # 2  φ pseudosection
            _ps_panel(
                "prof-panel-phi-ps", IDs.IMG_PHASE_PS, "prof-ps-phi-header"
            ),
            # 3  Phase Tensor
            _img_panel("prof-panel-pt", IDs.IMG_PT),
            # 4  Tipper
            _img_panel("prof-panel-tipper", IDs.IMG_TIPPER),
            # 5  2D Section
            _img_panel("prof-panel-section", IDs.IMG_SECTION),
            # 6  Publication view
            html.Div(
                [
                    html.Div(
                        html.Button(
                            [
                                html.I(className="bi bi-arrow-left me-1"),
                                "Exit Publication",
                            ],
                            id="prof-pub-back-btn",
                            className="btn btn-sm btn-outline-secondary",
                            n_clicks=0,
                        ),
                        className="prof-pub-back-bar",
                    ),
                    html.Div(
                        html.Img(
                            id=IDs.IMG_PUB,
                            src=empty_src(dark=True),
                            style=_img_style,
                        ),
                        className="fig-img-wrap profile-page-fig",
                        style={"flex": "1", "minHeight": "0"},
                    ),
                ],
                id="prof-panel-pub",
                className="prof-panel",
                style={"display": "none"},
            ),
        ],
        className="prof-view-panel",
    )

    main_area = html.Div(
        [
            dcc.Store(id=IDs.PROF_ACTIVE_TAB, data=_DEFAULT),
            tab_bar,
            info_bar,
            panels,
        ],
        className="page-main-area profile-page-main",
    )

    return html.Div(
        [sidebar, main_area],
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
                            [
                                html.I(
                                    className="bi bi-arrows-fullscreen me-1"
                                ),
                                "Esc  to redock",
                            ],
                            className="lb-hint",
                        ),
                        html.A(
                            [
                                html.I(className="bi bi-download me-1"),
                                "Save PNG",
                            ],
                            id="lightbox-download",
                            className="lb-btn",
                            href="#",
                            download="pycsamt_figure.png",
                        ),
                        html.Button(
                            [
                                html.I(
                                    className="bi bi-arrows-angle-contract me-1"
                                ),
                                "Redock",
                            ],
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
        (8, 8, 15, 0.0, 8),
        (5, 18, 72, 1.2, 6),
        (12, 28, 35, 0.4, 10),
        (4, 40, 80, 2.1, 7),
        (7, 52, 25, 0.8, 9),
        (10, 62, 60, 1.6, 8),
        (5, 73, 88, 0.2, 6),
        (8, 82, 42, 1.9, 9),
        (4, 90, 18, 0.6, 7),
        (6, 95, 65, 2.5, 8),
        (9, 35, 92, 1.1, 10),
        (3, 48, 8, 0.9, 6),
        (6, 68, 78, 1.7, 8),
        (4, 15, 50, 2.3, 7),
    ]
    particles = html.Div(
        [
            html.Div(
                className="wlc-particle",
                style={
                    "width": f"{sz}px",
                    "height": f"{sz}px",
                    "left": f"{lft}%",
                    "top": f"{top}%",
                    "animationDelay": f"{dly}s",
                    "animationDuration": f"{dur}s",
                },
            )
            for sz, lft, top, dly, dur in _DOTS
        ],
        className="wlc-particles",
    )

    feature_cards = html.Div(
        [
            html.Button(
                [
                    html.I(
                        className="bi bi-cloud-upload-fill wlc-card-icon",
                        style={"color": "var(--blue)"},
                    ),
                    html.Div("Load Data", className="wlc-card-title"),
                    html.Div(
                        "EDI · AVG · J · multi-line",
                        className="wlc-card-desc",
                    ),
                ],
                className="wlc-card",
                id=IDs.WELCOME_CARD_LOAD,
                n_clicks=0,
            ),
            html.Button(
                [
                    html.I(
                        className="bi bi-bezier2 wlc-card-icon",
                        style={"color": "var(--green)"},
                    ),
                    html.Div("Pipeline", className="wlc-card-title"),
                    html.Div(
                        "8-step processing workflow",
                        className="wlc-card-desc",
                    ),
                ],
                className="wlc-card",
                id=IDs.WELCOME_CARD_PIPE,
                n_clicks=0,
            ),
            html.Button(
                [
                    html.I(
                        className="bi bi-robot wlc-card-icon",
                        style={"color": "var(--mauve)"},
                    ),
                    html.Div("AI Agents", className="wlc-card-title"),
                    html.Div(
                        "ML · DL · LLM solvers", className="wlc-card-desc"
                    ),
                ],
                className="wlc-card",
                id=IDs.WELCOME_CARD_AGENTS,
                n_clicks=0,
            ),
            html.Button(
                [
                    html.I(
                        className="bi bi-bar-chart-line-fill wlc-card-icon",
                        style={"color": "var(--teal)"},
                    ),
                    html.Div("Visualize", className="wlc-card-title"),
                    html.Div(
                        "25 QC · 16 advanced plots", className="wlc-card-desc"
                    ),
                ],
                className="wlc-card",
                id=IDs.WELCOME_CARD_VIZ,
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
                        src="/icons/pycsamt-v2-symbol.svg",
                        className="wlc-logo",
                    ),
                    html.H1(
                        [
                            "Welcome to ",
                            html.Span("py", className="wlc-py"),
                            html.Span("CSAMT", className="wlc-csamt"),
                        ],
                        className="wlc-title",
                    ),
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
                        [
                            html.I(className="bi bi-cloud-upload-fill me-2"),
                            "Load EDI Data  —  Start",
                        ],
                        id=IDs.WELCOME_CTA,
                        className="wlc-cta",
                        n_clicks=0,
                    ),
                    html.P(
                        "Point to a parent folder to load all survey lines at once",
                        className="wlc-hint",
                    ),
                ],
                className="wlc-hero",
            ),
        ],
        id=IDs.WELCOME_OVERLAY,
        className="wlc-overlay",
    )


# ──────────────────────────────────────────────────────────────────────────────
# Tools off-canvas panel
# ──────────────────────────────────────────────────────────────────────────────

# Registry drives both the menu item list above and the panel content below
_TOOL_REGISTRY = [
    # (tool_id, icon_name, label, description)
    (
        "strike",
        "strike-analyzer",
        "Strike Analyzer",
        "Compute regional-strike angles and dimensionality from the impedance tensor.",
    ),
    (
        "validator",
        "EDI-validator",
        "EDI Validator",
        "Per-station data-quality checklist — flag bad impedances, missing frequencies.",
    ),
    (
        "converter",
        "format-converter",
        "Format Converter",
        "Export the loaded survey to EDI, CSV, or JSON.",
    ),
    (
        "batch",
        "batch-export",
        "Batch Export Plots",
        "Save every profile / pseudosection figure to a folder in one shot.",
    ),
    (
        "coords",
        "coordinate-transformer",
        "Coordinate Transformer",
        "Convert station coordinates between UTM and geographic Lat / Lon (pyproj).",
    ),
    (
        "elevation",
        "elevation",
        "Elevation Enrichment",
        "Fetch SRTM elevation for every loaded station via open-elevation API.",
    ),
    (
        "response",
        "station-response",
        "Station Response Inspector",
        "Full impedance-tensor response curves for a selected station.",
    ),
    (
        "strike-profile",
        "strike-profile",
        "Strike Profile Viewer",
        "Strike angle vs. station-position line chart with IQR ribbon.",
    ),
    (
        "phase-tensor",
        "phase-tensor",
        "Phase Tensor Map",
        "Geographic map of phase-tensor ellipses at a chosen period.",
    ),
    (
        "dimensionality",
        "dimensionnality",
        "Dimensionality Classifier",
        "Label each station × frequency as 1-D / 2-D / 3-D automatically.",
    ),
    (
        "freq-editor",
        "frequency-editor",
        "Frequency Editor",
        "Confidence-based QC: drop, mask, or recover frequency bands per station.",
    ),
    (
        "layered-model",
        "layered-model",
        "Layered Model Builder",
        "Build and preview a 1-D layered-earth resistivity model.",
    ),
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
        return html.Div(
            "Select a tool from the menu.", className="tool-empty"
        )

    tid, icon, label, desc = _TOOL_BY_ID[tool_id]

    return html.Div(
        [
            html.Div(
                [
                    _icon(icon, size=22, cls="tool-panel-icon"),
                    html.Div(
                        [
                            html.H6(label, className="tool-panel-title"),
                            html.P(desc, className="tool-panel-desc"),
                        ]
                    ),
                ],
                className="tool-panel-head",
            ),
            html.Hr(className="tool-sep"),
            html.Div(body, className="tool-panel-body"),
        ],
        className="tool-panel",
    )


def _tools_offcanvas() -> dbc.Offcanvas:
    return dbc.Offcanvas(
        [
            dcc.Store(id=IDs.ACTIVE_TOOL, data="strike"),
            dcc.Store(id="tools-run-store", storage_type="memory", data={}),
            # Tool picker: compact icon buttons for quick switching
            html.Div(
                [
                    html.Span(
                        "Available tools", className="offcanvas-section-lbl"
                    ),
                    html.Div(
                        [
                            html.Button(
                                _icon(icon, size=14, cls="menu-item-icon"),
                                id={"type": "tool-picker-btn", "index": tid},
                                className="tool-picker-btn",
                                title=label,
                                n_clicks=0,
                            )
                            for tid, icon, label, _ in _TOOL_REGISTRY
                        ],
                        className="tool-picker-rail",
                    ),
                ],
                className="offcanvas-section",
            ),
            # Active tool content — updated by callback
            html.Div(id="tools-panel-content", className="tool-panel-wrap"),
        ],
        id=IDs.TOOLS_CANVAS,
        title=html.Span(
            [_icon("tools", size=16, cls="menu-item-icon"), " Tools"],
            className="d-flex align-items-center gap-2",
        ),
        placement="end",
        scrollable=True,
        is_open=False,
        className="pycsamt-offcanvas",
        style={"width": "480px"},
    )


def _tools_path_browser_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle(
                    [
                        html.I(className="bi bi-folder2-open me-2"),
                        html.Span(
                            "Browse Server Path", id="tools-path-title"
                        ),
                    ]
                ),
                close_button=True,
            ),
            dbc.ModalBody(
                [
                    dcc.Store(id="tools-path-current", data="~"),
                    dcc.Store(id="tools-path-target", data={}),
                    html.Div(
                        [
                            html.Span(
                                [
                                    html.I(className="bi bi-hdd me-1"),
                                    "Current:",
                                ],
                                className="tools-path-kicker",
                            ),
                            html.Code(
                                id="tools-path-cwd",
                                children="~",
                                className="tools-path-cwd",
                            ),
                        ],
                        className="tools-path-head mb-2",
                    ),
                    dbc.Button(
                        [
                            html.I(className="bi bi-arrow-up me-1"),
                            "Parent folder",
                        ],
                        id="tools-path-up",
                        color="secondary",
                        outline=True,
                        size="sm",
                        className="mb-2",
                    ),
                    html.Div(
                        id="tools-path-list", className="tools-path-list"
                    ),
                    html.Hr(className="tool-sep"),
                    html.Div(
                        [
                            html.Div(
                                "Create folder in current location",
                                className="tools-path-kicker mb-1",
                            ),
                            dbc.InputGroup(
                                [
                                    dbc.Input(
                                        id="tools-path-new",
                                        placeholder="new_folder_name",
                                        type="text",
                                        size="sm",
                                    ),
                                    dbc.Button(
                                        [
                                            html.I(
                                                className="bi bi-folder-plus me-1"
                                            ),
                                            "Create",
                                        ],
                                        id="tools-path-mkdir",
                                        color="secondary",
                                        size="sm",
                                    ),
                                ],
                                size="sm",
                            ),
                        ],
                    ),
                ]
            ),
            dbc.ModalFooter(
                [
                    dbc.Button(
                        "Cancel",
                        id="tools-path-cancel",
                        color="secondary",
                        size="sm",
                        className="me-auto",
                    ),
                    dbc.Button(
                        [
                            html.I(className="bi bi-check2-circle me-1"),
                            "Select folder",
                        ],
                        id="tools-path-select-folder",
                        color="success",
                        size="sm",
                    ),
                ]
            ),
        ],
        id="tools-path-modal",
        is_open=False,
        size="lg",
        backdrop="static",
        className="tools-path-modal",
    )


# ──────────────────────────────────────────────────────────────────────────────
# Settings off-canvas panel
# ──────────────────────────────────────────────────────────────────────────────


def _settings_offcanvas() -> dbc.Offcanvas:

    def _sec_hdr(icon_cls: str, title: str) -> html.Div:
        return html.Div(
            [
                html.I(className=f"bi {icon_cls} me-2 settings-sec-icon"),
                html.Span(title, className="settings-sec-title"),
            ],
            className="settings-sec-hdr",
        )

    def _row_label(text: str) -> html.Div:
        return html.Div(text, className="ctrl-label mt-1 mb-1")

    # ── AI & Agents ───────────────────────────────────────────────────────────
    ai_section = html.Div(
        [
            _sec_hdr("bi-stars", "AI & Agents"),
            # Provider pill buttons
            _row_label("Provider"),
            html.Div(
                [
                    html.Button(
                        [html.I(className="bi bi-cloud me-1"), "Claude"],
                        id={
                            "type": "settings-provider-btn",
                            "index": "claude",
                        },
                        className="settings-provider-btn active",
                        n_clicks=0,
                    ),
                    html.Button(
                        [html.I(className="bi bi-stars me-1"), "OpenAI"],
                        id={
                            "type": "settings-provider-btn",
                            "index": "openai",
                        },
                        className="settings-provider-btn",
                        n_clicks=0,
                    ),
                    html.Button(
                        [
                            html.I(className="bi bi-gem me-1"),
                            "Gemini",
                        ],
                        id={
                            "type": "settings-provider-btn",
                            "index": "gemini",
                        },
                        className="settings-provider-btn",
                        n_clicks=0,
                    ),
                    html.Button(
                        [
                            html.I(className=("bi bi-lightning-charge me-1")),
                            "DeepSeek",
                        ],
                        id={
                            "type": "settings-provider-btn",
                            "index": "deepseek",
                        },
                        className="settings-provider-btn",
                        n_clicks=0,
                    ),
                ],
                className="settings-provider-row mb-3",
            ),
            # API Key
            _row_label("API Key"),
            dbc.InputGroup(
                [
                    dbc.Input(
                        id=IDs.SETTINGS_AI_KEY,
                        type="password",
                        placeholder="Paste your API key…",
                        size="sm",
                        autoComplete="new-password",
                    ),
                    dbc.Button(
                        html.I(
                            className="bi bi-eye", id="settings-ai-eye-icon"
                        ),
                        id=IDs.SETTINGS_AI_KEY_SHOW,
                        color="secondary",
                        size="sm",
                        outline=True,
                        n_clicks=0,
                        title="Show / hide key",
                    ),
                ],
                size="sm",
                className="mb-1",
            ),
            html.Div(
                [
                    html.I(className="bi bi-shield-lock me-1"),
                    "Stored in browser localStorage only",
                ],
                className="settings-key-note mb-3",
            ),
            # Model
            _row_label("Model"),
            dbc.Select(
                id=IDs.SETTINGS_AI_MODEL,
                options=[
                    {
                        "label": "claude-sonnet-4-6  (default)",
                        "value": "claude-sonnet-4-6",
                    },
                    {"label": "claude-opus-4-8", "value": "claude-opus-4-8"},
                    {
                        "label": "claude-haiku-4-5",
                        "value": "claude-haiku-4-5-20251001",
                    },
                ],
                value="claude-sonnet-4-6",
                size="sm",
                className="mb-3",
            ),
            # Test connection
            dbc.Button(
                [html.I(className="bi bi-plug-fill me-1"), "Test Connection"],
                id=IDs.SETTINGS_AI_TEST,
                color="primary",
                size="sm",
                outline=True,
                className="w-100 mb-2",
                n_clicks=0,
            ),
            html.Div(
                id=IDs.SETTINGS_AI_STATUS, className="settings-ai-status"
            ),
            # Hidden store for selected provider
            dcc.Store(id=IDs.SETTINGS_AI_PROVIDER, data="claude"),
        ],
        className="settings-section settings-ai",
    )

    # ── Computation ───────────────────────────────────────────────────────────
    comp_section = html.Div(
        [
            _sec_hdr("bi-cpu", "Computation"),
            _row_label("Pseudosection interpolation"),
            dbc.Select(
                id="settings-interp",
                options=[
                    {"label": "Bilinear", "value": "bilinear"},
                    {"label": "Nearest", "value": "nearest"},
                    {"label": "Cubic", "value": "cubic"},
                ],
                value="bilinear",
                size="sm",
                className="mb-2",
            ),
            _row_label("Skin-depth formula"),
            dbc.Select(
                id="settings-skindepth",
                options=[
                    {"label": "Standard (Cagniard)", "value": "cagniard"},
                    {"label": "True skin depth", "value": "true"},
                ],
                value="cagniard",
                size="sm",
                className="mb-2",
            ),
            _row_label("Topography DEM source"),
            dbc.Select(
                id="settings-dem",
                options=[
                    {
                        "label": "Open-Elevation API",
                        "value": "open-elevation",
                    },
                    {"label": "SRTM local file", "value": "srtm"},
                    {"label": "None (disable)", "value": "none"},
                ],
                value="open-elevation",
                size="sm",
                className="mb-2",
            ),
            _row_label("Default colormap"),
            dbc.Select(
                id="settings-cmap",
                options=[
                    {"label": c, "value": c}
                    for c in [
                        "viridis",
                        "plasma",
                        "seismic",
                        "jet",
                        "RdBu_r",
                        "rainbow",
                    ]
                ],
                value="viridis",
                size="sm",
            ),
        ],
        className="settings-section settings-comp",
    )

    # ── Profile ───────────────────────────────────────────────────────────────
    profile_section = html.Div(
        [
            _sec_hdr("bi-person-gear", "Profile"),
            html.Div(
                [
                    dbc.Button(
                        [html.I(className="bi bi-download me-1"), "Save"],
                        id=IDs.BTN_SETTINGS_SAVE,
                        size="sm",
                        color="secondary",
                        className="flex-grow-1",
                        n_clicks=0,
                    ),
                    dcc.Upload(
                        id=IDs.SETTINGS_UPLOAD,
                        children=dbc.Button(
                            [html.I(className="bi bi-upload me-1"), "Load"],
                            size="sm",
                            color="secondary",
                            className="w-100",
                        ),
                        accept=".json",
                        multiple=False,
                        className="flex-grow-1",
                    ),
                ],
                className="d-flex gap-2 mb-2",
            ),
            dbc.Button(
                [
                    html.I(className="bi bi-arrow-counterclockwise me-1"),
                    "Reset to Defaults",
                ],
                id=IDs.BTN_SETTINGS_RESET,
                size="sm",
                color="danger",
                outline=True,
                className="w-100",
                n_clicks=0,
            ),
        ],
        className="settings-section settings-profile",
    )

    return dbc.Offcanvas(
        [
            dcc.Store(id=IDs.SETTINGS_STORE, storage_type="local", data={}),
            ai_section,
            comp_section,
            profile_section,
            html.Div(id=IDs.SETTINGS_FEEDBACK, className="settings-feedback"),
        ],
        id=IDs.SETTINGS_CANVAS,
        title=html.Span(
            [html.I(className="bi bi-gear-fill me-2"), "Settings"],
            className="d-flex align-items-center",
        ),
        placement="end",
        scrollable=True,
        is_open=False,
        className="pycsamt-offcanvas",
        style={"width": "400px"},
    )


# ──────────────────────────────────────────────────────────────────────────────
# Recompute off-canvas panel
# ──────────────────────────────────────────────────────────────────────────────


def _lines_offcanvas() -> dbc.Offcanvas:
    """Offcanvas panel for managing active / muted survey lines."""
    _mode_hint = {
        "folder": "Lines come from subfolder names (works with multi-folder uploads).",
        "auto": "Prefix before the first '-' or '_' in station IDs becomes the line name. "
        "Purely numeric prefixes get an 'L' prepended (e.g. '18-001A' → L18).",
        "edit": "Type new names in the inputs below, then click Apply.",
    }
    return dbc.Offcanvas(
        [
            # ── Header ──────────────────────────────────────────────────
            html.Div(
                [
                    html.Div(
                        [
                            html.I(
                                className="bi bi-layers-half me-2",
                                style={
                                    "color": "var(--teal)",
                                    "fontSize": "16px",
                                },
                            ),
                            html.Span(
                                "Survey Lines", className="rcp-canvas-title"
                            ),
                        ],
                        className="d-flex align-items-center",
                    ),
                    html.Div(
                        "Active lines flow through QC, Correction and Profile. "
                        "Muted lines remain on the map.",
                        className="lines-canvas-hint",
                    ),
                ],
                className="mb-3",
            ),
            # ── Line assignment mode ─────────────────────────────────────
            html.Div(
                [
                    html.Div(
                        "Line assignment",
                        className="offcanvas-section-lbl mb-1",
                    ),
                    dbc.RadioItems(
                        id=IDs.LINES_DETECT_MODE,
                        options=[
                            {
                                "label": html.Span(
                                    [
                                        html.I(
                                            className="bi bi-folder2-open me-2",
                                            style={"color": "#f9e2af"},
                                        ),
                                        "Folder names",
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
                                        "Auto-detect from ID",
                                    ],
                                    className="d-flex align-items-center",
                                ),
                                "value": "auto",
                            },
                            {
                                "label": html.Span(
                                    [
                                        html.I(
                                            className="bi bi-pencil-square me-2",
                                            style={"color": "#a6e3a1"},
                                        ),
                                        "Edit / Rename",
                                    ],
                                    className="d-flex align-items-center",
                                ),
                                "value": "edit",
                            },
                        ],
                        value="folder",
                        className="lines-mode-radio mb-2",
                        inputStyle={"cursor": "pointer"},
                        labelStyle={
                            "cursor": "pointer",
                            "fontSize": "12px",
                            "padding": "3px 0",
                        },
                    ),
                    # Action buttons — visibility toggled by callback
                    dbc.Button(
                        [
                            html.I(className="bi bi-stars me-2"),
                            "Detect Lines",
                        ],
                        id=IDs.BTN_LINES_DETECT,
                        color="info",
                        size="sm",
                        outline=True,
                        className="w-100 mb-1",
                        style={"display": "none"},
                        n_clicks=0,
                    ),
                    dbc.Button(
                        [
                            html.I(className="bi bi-check2-circle me-2"),
                            "Apply Renames",
                        ],
                        id=IDs.BTN_LINES_RENAME_APPLY,
                        color="success",
                        size="sm",
                        outline=True,
                        className="w-100 mb-1",
                        style={"display": "none"},
                        n_clicks=0,
                    ),
                    # Inline hint for current mode
                    html.Div(
                        id="lines-mode-hint",
                        children=_mode_hint["folder"],
                        style={
                            "fontSize": "10px",
                            "color": "#6c7086",
                            "marginTop": "4px",
                        },
                    ),
                ],
                className="lines-mode-section mb-3",
            ),
            html.Hr(className="my-2"),
            # ── All / None quick select ──────────────────────────────────
            html.Div(
                [
                    dbc.Button(
                        "All",
                        id=IDs.BTN_LINES_ALL,
                        color="success",
                        size="sm",
                        outline=True,
                        className="me-2",
                    ),
                    dbc.Button(
                        "None",
                        id=IDs.BTN_LINES_NONE,
                        color="danger",
                        size="sm",
                        outline=True,
                    ),
                ],
                className="d-flex mb-3",
            ),
            html.Hr(className="my-2"),
            # ── Dynamic line list (populated by callback) ────────────────
            html.Div(
                id=IDs.LINES_PANEL_BODY,
                children=[
                    html.Div(
                        "Load survey data to see lines.",
                        className="text-muted small text-center mt-3",
                    )
                ],
            ),
        ],
        id=IDs.LINES_OFFCANVAS,
        title="",
        placement="end",
        scrollable=True,
        is_open=False,
        style={"width": "360px"},
        className="rcp-offcanvas",
    )


def _recompute_offcanvas() -> dbc.Offcanvas:
    def _sec(title: str, children: list) -> html.Div:
        return html.Div(
            [
                html.Span(title, className="offcanvas-section-lbl"),
                *children,
            ],
            className="offcanvas-section rcp-section",
        )

    # ── Data source ──────────────────────────────────────────────────────
    source_sec = _sec(
        "Data Source",
        [
            dbc.RadioItems(
                id=IDs.RECOMPUTE_SOURCE,
                options=[
                    {
                        "label": html.Span(
                            [
                                html.I(className="bi bi-layers me-1"),
                                "Use currently loaded data",
                            ]
                        ),
                        "value": "current",
                    },
                    {
                        "label": html.Span(
                            [
                                html.I(className="bi bi-folder2-open me-1"),
                                "Load from folder path",
                            ]
                        ),
                        "value": "path",
                    },
                ],
                value="current",
                className="rcp-source-radio",
            ),
            # Folder path input — shown only when source == "path"
            html.Div(
                id=IDs.RECOMPUTE_PATH_SEC,
                style={"display": "none"},
                children=[
                    dbc.Input(
                        id=IDs.RECOMPUTE_PATH,
                        type="text",
                        placeholder="/path/to/edi/folder",
                        size="sm",
                        className="mt-2",
                    ),
                    html.Small(
                        "Folder is scanned recursively — sub-folders become profile lines.",
                        className="text-muted rcp-hint",
                    ),
                ],
            ),
        ],
    )

    # ── Options ──────────────────────────────────────────────────────────
    options_sec = _sec(
        "Transfer Function Options",
        [
            dbc.Checklist(
                id=IDs.RECOMPUTE_RESPHASE,
                options=[
                    {
                        "label": "Recompute ρₐ / φ from impedance Z  (recommended)",
                        "value": "on",
                    }
                ],
                value=["on"],
                switch=True,
                className="recompute-switch mb-2 mt-2",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label(
                                "Rotation angle (°)",
                                className="recompute-lbl",
                            ),
                            dbc.Input(
                                id=IDs.RECOMPUTE_ROTATE_ANG,
                                type="number",
                                placeholder="0 = none",
                                min=-180,
                                max=180,
                                step=0.1,
                                size="sm",
                            ),
                        ],
                        style={"flex": "1"},
                    ),
                    html.Div(
                        [
                            html.Label(
                                "Rotate components", className="recompute-lbl"
                            ),
                            dbc.Checklist(
                                id=IDs.RECOMPUTE_COMPS,
                                options=[
                                    {"label": "Z", "value": "Z"},
                                    {"label": "Tipper", "value": "Tip"},
                                ],
                                value=["Z", "Tip"],
                                inline=True,
                                className="recompute-comp-check",
                            ),
                        ],
                        style={"flex": "1.5"},
                    ),
                ],
                className="d-flex gap-2 mb-2 align-items-start",
            ),
            html.Div(
                [
                    html.Label(
                        "Frequency band (Hz)", className="recompute-lbl"
                    ),
                    html.Div(
                        [
                            dbc.Input(
                                id=IDs.RECOMPUTE_FMIN,
                                type="number",
                                placeholder="f min",
                                size="sm",
                                style={"flex": "1"},
                            ),
                            html.Span(
                                "–",
                                className="text-muted align-self-center px-1",
                                style={"fontSize": "13px"},
                            ),
                            dbc.Input(
                                id=IDs.RECOMPUTE_FMAX,
                                type="number",
                                placeholder="f max",
                                size="sm",
                                style={"flex": "1"},
                            ),
                        ],
                        className="d-flex gap-1 align-items-center",
                    ),
                    html.Small(
                        "Leave blank to keep all frequencies.",
                        className="text-muted rcp-hint",
                    ),
                ],
                className="mb-2",
            ),
            html.Div(
                [
                    html.Label(
                        "Fill missing values", className="recompute-lbl"
                    ),
                    dbc.Select(
                        id=IDs.RECOMPUTE_FILL,
                        options=[
                            {"label": "None — leave as NaN", "value": ""},
                            {"label": "Fill with 0", "value": "zero"},
                            {"label": "Fill with NaN", "value": "nan"},
                        ],
                        value="",
                        size="sm",
                    ),
                ],
                className="mb-1",
            ),
        ],
    )

    # ── Run button ────────────────────────────────────────────────────────
    run_sec = html.Div(
        [
            dbc.Button(
                [
                    html.I(className="bi bi-arrow-repeat me-2"),
                    "Recompute & Load",
                ],
                id=IDs.BTN_RECOMPUTE,
                size="md",
                className="rcp-run-btn w-100",
            ),
            html.Div(
                id=IDs.RECOMPUTE_FEEDBACK,
                className="small mt-2 text-center rcp-feedback",
            ),
        ],
        className="px-3 pb-2 pt-1",
    )

    # ── Confirm re-run card (shown only when data already recomputed) ─────
    confirm_sec = html.Div(
        id=IDs.RECOMPUTE_CONFIRM_SEC,
        style={"display": "none"},
        className="px-3 pb-3",
        children=[
            html.Div(
                [
                    html.Div(
                        [
                            html.I(
                                className="bi bi-exclamation-triangle-fill me-2 rcp-confirm-icon"
                            ),
                            html.Span(
                                "Sites already recomputed",
                                className="rcp-confirm-title",
                            ),
                        ],
                        className="d-flex align-items-center mb-1",
                    ),
                    html.P(
                        "The active dataset was already processed by the recompute engine. "
                        "Run again with the current settings?",
                        className="rcp-confirm-body",
                    ),
                    html.Div(
                        [
                            dbc.Button(
                                "Cancel",
                                id=IDs.BTN_RECOMPUTE_CANCEL,
                                size="sm",
                                color="secondary",
                                outline=True,
                            ),
                            dbc.Button(
                                [
                                    html.I(
                                        className="bi bi-arrow-repeat me-1"
                                    ),
                                    "Yes, recompute again",
                                ],
                                id=IDs.BTN_RECOMPUTE_YES,
                                size="sm",
                                className="rcp-confirm-yes-btn",
                            ),
                        ],
                        className="d-flex gap-2 justify-content-end mt-2",
                    ),
                ],
                className="rcp-confirm-card",
            ),
        ],
    )

    # ── Progress panel ────────────────────────────────────────────────────
    progress_sec = html.Div(
        id=IDs.RECOMPUTE_PROG_SEC,
        style={"display": "none"},
        className="px-3 pb-3",
        children=[
            html.Div(
                [
                    html.Div(
                        [
                            html.I(
                                className="bi bi-arrow-repeat recompute-spin me-1",
                                id=IDs.RECOMPUTE_SPIN_ICON,
                            ),
                            html.Span(
                                id=IDs.RECOMPUTE_TICKER,
                                children="Initialising…",
                                className="recompute-ticker-text",
                            ),
                        ],
                        style={"flex": "1"},
                    ),
                    dbc.Badge(
                        "Waiting",
                        id=IDs.RECOMPUTE_BADGE,
                        color="secondary",
                        className="ms-2 align-self-center",
                    ),
                ],
                className="d-flex align-items-center mb-2",
            ),
            dbc.Progress(
                id=IDs.RECOMPUTE_PROGRESS,
                value=0,
                max=100,
                striped=True,
                animated=True,
                color="teal",
                style={"height": "8px", "borderRadius": "6px"},
                className="mb-2",
            ),
            html.Div(
                id=IDs.RECOMPUTE_LOG,
                className="recompute-log",
                children=[],
            ),
        ],
    )

    return dbc.Offcanvas(
        [
            source_sec,
            html.Hr(className="my-1"),
            options_sec,
            html.Hr(className="my-1"),
            run_sec,
            confirm_sec,
            progress_sec,
        ],
        id=IDs.RECOMPUTE_CANVAS,
        title=html.Span(
            [
                html.I(
                    className="bi bi-arrow-repeat me-2",
                    style={"color": "var(--teal)"},
                ),
                "Recompute Transfer Functions",
            ],
            className="d-flex align-items-center rcp-canvas-title",
        ),
        placement="end",
        scrollable=True,
        is_open=False,
        className="pycsamt-offcanvas rcp-offcanvas",
        style={"width": "440px"},
    )


# ──────────────────────────────────────────────────────────────────────────────
# About modal
# ──────────────────────────────────────────────────────────────────────────────


def _about_feature(icon: str, color: str, label: str) -> html.Div:
    """One feature chip in the About modal highlights grid."""
    return html.Div(
        [
            html.I(
                className=f"bi {icon}",
                style={"color": f"var(--{color})"},
            ),
            html.Span(label),
        ],
        className="about-feature",
    )


def _about_modal() -> dbc.Modal:
    return dbc.Modal(
        [
            dbc.ModalHeader(
                dbc.ModalTitle([_icon("help", size=16), " About pyCSAMT"]),
                close_button=True,
            ),
            dbc.ModalBody(
                [
                    html.Div(
                        [
                            html.Img(
                                src="/icons/pycsamt-v2-symbol.svg",
                                className="about-logo-anim",
                                style={
                                    "height": "72px",
                                    "marginBottom": "12px",
                                },
                            ),
                            html.H4(
                                [
                                    html.Span(
                                        "py",
                                        style={
                                            "color": "var(--peach)",
                                            "fontStyle": "italic",
                                        },
                                    ),
                                    html.Span(
                                        "CSAMT",
                                        style={"color": "var(--blue)"},
                                    ),
                                    html.Span(
                                        " v2",
                                        style={
                                            "color": "var(--ov0)",
                                            "fontSize": "0.8em",
                                        },
                                    ),
                                ],
                                className="mb-1",
                            ),
                            html.P(
                                "Controlled-Source Audio Magnetotellurics · "
                                "Python Processing Suite",
                                className="text-muted small mb-3",
                            ),
                        ],
                        className="text-center",
                    ),
                    html.P(
                        "An end-to-end toolkit for CSAMT / AMT / MT data — "
                        "from raw EDI files to quality control, distortion "
                        "correction, AI-assisted inversion, and publication "
                        "figures, with an interactive web studio and a "
                        "conversational agent.",
                        className="about-lead",
                    ),
                    html.Div("Highlights", className="about-sec-label"),
                    html.Div(
                        [
                            _about_feature(
                                "bi-clipboard-check",
                                "blue",
                                "Quality control",
                            ),
                            _about_feature(
                                "bi-sliders",
                                "teal",
                                "Static-shift correction",
                            ),
                            _about_feature(
                                "bi-compass", "mauve", "Phase-tensor & strike"
                            ),
                            _about_feature("bi-funnel", "green", "Denoising"),
                            _about_feature(
                                "bi-cpu", "peach", "AI inversion (1D/2D/3D)"
                            ),
                            _about_feature(
                                "bi-map",
                                "yellow",
                                "Interactive maps & sections",
                            ),
                        ],
                        className="about-features",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Span(
                                        "Author", className="about-key"
                                    ),
                                    html.Span(
                                        "L. Kouadio  ·  etanoyau@gmail.com",
                                        className="about-val",
                                    ),
                                ],
                                className="about-row",
                            ),
                            html.Div(
                                [
                                    html.Span(
                                        "License", className="about-key"
                                    ),
                                    html.Span(
                                        "LGPL-3.0", className="about-val"
                                    ),
                                ],
                                className="about-row",
                            ),
                            html.Div(
                                [
                                    html.Span(
                                        "Framework", className="about-key"
                                    ),
                                    html.Span(
                                        "Dash 4  ·  DBC 2  ·  Plotly 6",
                                        className="about-val",
                                    ),
                                ],
                                className="about-row",
                            ),
                            html.Div(
                                [
                                    html.Span("Stack", className="about-key"),
                                    html.Span(
                                        "NumPy · SciPy · Matplotlib · pyproj",
                                        className="about-val",
                                    ),
                                ],
                                className="about-row",
                            ),
                        ],
                        className="about-table",
                    ),
                ]
            ),
            dbc.ModalFooter(
                [
                    html.A(
                        [
                            _icon("docs", size=13, cls="menu-item-icon"),
                            " Documentation",
                        ],
                        href="https://pycsamt.readthedocs.io",
                        target="_blank",
                        className="btn btn-outline-secondary btn-sm me-2",
                    ),
                    html.A(
                        [
                            _icon("github", size=13, cls="menu-item-icon"),
                            " GitHub",
                        ],
                        href="https://github.com/earthai-tech/pycsamt",
                        target="_blank",
                        className="btn btn-outline-secondary btn-sm me-2",
                    ),
                    dbc.Button(
                        "Close",
                        id="about-close-btn",
                        color="secondary",
                        size="sm",
                    ),
                ]
            ),
        ],
        id=IDs.ABOUT_MODAL,
        is_open=False,
        size="md",
        centered=True,
    )


def _session_offcanvas() -> dbc.Offcanvas:
    """Session save / restore offcanvas panel."""

    def _sec(title, icon, children, border="blue"):
        return html.Div(
            [
                html.Div(
                    [
                        html.I(
                            className=f"bi bi-{icon} me-2",
                            style={"color": f"var(--{border})"},
                        ),
                        html.Span(title, className="settings-sec-title"),
                    ],
                    className="settings-sec-hdr",
                ),
                *children,
            ],
            className=f"settings-section settings-{border}",
        )

    what_is_saved = html.Ul(
        [
            html.Li("Station records & metadata (STORE_DATA)"),
            html.Li("Active survey lines"),
            html.Li("Correction pipeline steps"),
            html.Li("Elevation fetch / correction results"),
            html.Li("Computation preferences (Settings)"),
        ],
        className="session-what-list",
    )

    body = [
        # ── auto-save indicator ───────────────────────────────────────────
        html.Div(
            [
                html.I(
                    className="bi bi-info-circle me-1",
                    style={"color": "var(--blue)"},
                ),
                "Sessions are auto-saved to your browser on every workflow "
                "action.  Use Download / Upload to share across machines.",
            ],
            className="session-info-note mb-3",
        ),
        # ── Survey note ───────────────────────────────────────────────────
        _sec(
            "Survey Note",
            "pencil-square",
            [
                html.Div(
                    "Optional label for this snapshot:",
                    className="ctrl-label",
                ),
                dbc.Input(
                    id=IDs.SESSION_NOTE,
                    placeholder="e.g. WILLY AMT 2024 — pre-inversion",
                    size="sm",
                    className="mb-1",
                ),
            ],
            border="blue",
        ),
        # ── Save section ──────────────────────────────────────────────────
        _sec(
            "Save Session",
            "floppy",
            [
                html.Div(
                    [
                        dbc.Button(
                            [
                                html.I(className="bi bi-download me-1"),
                                "Download JSON",
                            ],
                            id=IDs.BTN_SESSION_SAVE,
                            color="primary",
                            size="sm",
                            className="me-2",
                            n_clicks=0,
                        ),
                        html.Div(
                            id=IDs.SESSION_AUTOSAVE,
                            className="session-autosave-chip",
                        ),
                    ],
                    className="d-flex align-items-center gap-2 mb-2",
                ),
                html.Details(
                    [
                        html.Summary(
                            "What is included?",
                            className="session-detail-summary",
                        ),
                        what_is_saved,
                    ],
                    className="mb-1",
                ),
            ],
            border="green",
        ),
        # ── Restore / Load ────────────────────────────────────────────────
        _sec(
            "Restore Session",
            "arrow-counterclockwise",
            [
                dcc.Upload(
                    id=IDs.SESSION_UL,
                    children=html.Div(
                        [
                            html.I(
                                className="bi bi-file-earmark-arrow-up me-2"
                            ),
                            "Drop session JSON here or ",
                            html.A("browse"),
                        ]
                    ),
                    accept=".json",
                    className="session-upload-zone mb-2",
                ),
                dbc.Button(
                    [
                        html.I(className="bi bi-arrow-counterclockwise me-1"),
                        "Restore browser snapshot",
                    ],
                    id=IDs.BTN_SESSION_RESTORE,
                    color="secondary",
                    size="sm",
                    outline=True,
                    className="w-100 mb-1",
                    n_clicks=0,
                ),
                dbc.Button(
                    [html.I(className="bi bi-trash me-1"), "Clear snapshot"],
                    id=IDs.BTN_SESSION_CLEAR,
                    color="danger",
                    size="sm",
                    outline=True,
                    className="w-100",
                    n_clicks=0,
                ),
            ],
            border="mauve",
        ),
        # ── Feedback ──────────────────────────────────────────────────────
        html.Div(id=IDs.SESSION_FEEDBACK, className="session-feedback mt-2"),
        dcc.Download(id=IDs.SESSION_DL),
    ]

    return dbc.Offcanvas(
        body,
        id=IDs.SESSION_CANVAS,
        title="Workflow Session",
        placement="end",
        is_open=False,
        style={"width": "370px"},
        className="settings-canvas",
    )


def layout() -> html.Div:
    """Return the complete Dash application layout."""
    from pycsamt.app.web.pages import (  # noqa: PLC0415
        advanced,
        agents_page,
        correction,
        forward,
        interpretation,
        inv_results,
        inversion,
        map3d,
        pipeline,
        qc_page,
        tdem,
    )

    _extra_pages = [
        ("qc", qc_page),
        ("correction", correction),
        ("advanced", advanced),
        ("tdem", tdem),
        ("pipeline", pipeline),
        ("forward", forward),
        ("inversion", inversion),
        ("interpretation", interpretation),
        ("map3d", map3d),
        ("inv-results", inv_results),
        ("agents", agents_page),
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
    extra_divs.insert(
        0,
        html.Div(
            _map_page_content(),
            id="page-map",
            className="page-content",
            style={"display": "none"},
        ),
    )
    extra_divs.insert(
        1,
        html.Div(
            _profile_page_content(),
            id="page-profile",
            className="page-content",
            style={"display": "none"},
        ),
    )

    return html.Div(
        [
            # Client-side state stores
            dcc.Store(id=IDs.SESSION_ID, storage_type="session"),
            dcc.Store(id=IDs.STORE_DATA, storage_type="memory"),
            dcc.Store(id=IDs.STORE_SELECTION, storage_type="memory"),
            dcc.Store(id=IDs.STORE_THEME, storage_type="local", data="dark"),
            dcc.Store(
                id=IDs.STORE_SIDEBAR, storage_type="local", data="expanded"
            ),
            dcc.Store(id=IDs.NAV_SECTION, storage_type="memory", data="home"),
            # Folder-browse store — written by JS via set_props
            dcc.Store(id=IDs.FOLDER_UPLOAD_STORE, storage_type="memory"),
            # Recompute — background-job polling state
            dcc.Store(id=IDs.RECOMPUTE_STORE, storage_type="memory"),
            # Correction page — applied-step list (JSON-serialisable)
            dcc.Store(
                id=IDs.CORR_STORE, data={"steps": []}, storage_type="memory"
            ),
            # Active lines filter — {"active": [...], "all": [...]}
            dcc.Store(id=IDs.STORE_ACTIVE_LINES, storage_type="memory"),
            # Station-panel local line filter — None means "show all"
            dcc.Store(
                id=IDs.STATION_LINE_FILTER, storage_type="memory", data=None
            ),
            dcc.Store(
                id=IDs.STORE_TDEM_LOADED, storage_type="memory", data=False
            ),
            dcc.Interval(
                id=IDs.RECOMPUTE_INTERVAL,
                interval=450,
                n_intervals=0,
                disabled=True,
            ),
            # Polling interval (agent result streaming)
            dcc.Interval(
                id=IDs.INTERVAL, interval=2000, n_intervals=0, disabled=True
            ),
            # Download triggers
            dcc.Download(id="settings-download"),
            dcc.Download(id=IDs.INTERP_DOWNLOAD),
            dcc.Download(id=IDs.CORR_DOWNLOAD),
            dcc.Download(id=IDs.MAP3D_DOWNLOAD),
            # ── Persistent tool result stores ─────────────────────────────
            # These live here (not inside tool bodies) so they survive tool
            # tab switches.  storage_type="session" keeps them within the
            # browser tab but clears on full page refresh.
            dcc.Store(
                id="tool-elev-raw-store", storage_type="session", data=None
            ),
            dcc.Store(
                id="tool-elev-corrected-store",
                storage_type="session",
                data=None,
            ),
            dcc.Store(
                id="tool-elev-mode-store",
                storage_type="session",
                data="fetch",
            ),
            dcc.Store(
                id="tool-elev-scope-store", storage_type="session", data="all"
            ),
            dcc.Store(
                id="tool-elev-upload-store", storage_type="session", data=None
            ),
            dcc.Store(
                id="tool-coord-mode-store",
                storage_type="session",
                data="single",
            ),
            dcc.Store(
                id="tool-coord-upload-store",
                storage_type="session",
                data=None,
            ),
            dcc.Store(
                id="tool-coord-result-store",
                storage_type="session",
                data=None,
            ),
            # Generic per-tool last-run output cache (all tools write here after a run)
            dcc.Store(
                id="tool-outputs-store", storage_type="session", data={}
            ),
            # Topography file upload (map 3D)
            dcc.Store(
                id=IDs.MAP3D_TOPO_UPLOAD_STORE,
                storage_type="session",
                data=None,
            ),
            dcc.Download(id="tool-elev-dl-csv"),
            dcc.Download(id="tool-elev-dl-h5"),
            dcc.Download(id="tool-coord-download"),
            # ── Session snapshot (localStorage — survives page refresh) ───
            dcc.Store(
                id=IDs.SESSION_SNAPSHOT, storage_type="local", data=None
            ),
            # Navbar
            _navbar(),
            # Load data modal
            _load_modal(),
            # Toasts (floating, top-right)
            _error_toast(),
            _success_toast(),
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
            # Tools, Settings, Lines, Recompute, About panels
            _tools_path_browser_modal(),
            _tools_offcanvas(),
            _settings_offcanvas(),
            _lines_offcanvas(),
            _recompute_offcanvas(),
            _session_offcanvas(),
            _about_modal(),
            # Welcome / landing overlay (hidden once data is loaded)
            _welcome_overlay(),
        ],
        id="app-root",
        style={"minHeight": "100vh"},
    )
