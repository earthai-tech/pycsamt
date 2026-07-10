# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Component ID constants for the map-view platform."""

from __future__ import annotations


class IDs:
    # ── root / stores ─────────────────────────
    ROOT             = "mv-root"
    SESSION_ID       = "mv-session-id"
    STORE_DATA       = "mv-store-data"        # station_records + meta (JSON)
    STORE_THEME      = "mv-store-theme"
    STORE_VIEW       = "mv-store-view"        # active view name
    STORE_CONTROLS   = "mv-store-controls"    # per-view control state
    STORE_LINES      = "mv-store-lines"       # {all, active, colors}
    STORE_LINE_FILTER = "mv-store-line-filter"
    STORE_SELECTION  = "mv-store-selection"   # clicked station
    STORE_MASKED     = "mv-store-masked"      # excluded station ids
    LOAD_MODE_STORE  = "mv-load-mode-store"   # replace | append
    FOLDER_STORE     = "mv-folder-store"      # JS folder upload payload
    UPLOAD_SELECTION = "mv-upload-selection"  # staged file entries
    SOURCE_SELECTION = "mv-source-selection"
    INV_FOLDER_STORE = "mv-inv-folder-store"  # JS ModEM folder upload payload

    # ── top bar ───────────────────────────────
    BTN_SIDEBAR      = "mv-btn-sidebar"
    BTN_INSPECTOR    = "mv-btn-inspector"
    BTN_LOAD         = "mv-btn-load"
    BTN_ADD_LINE     = "mv-btn-add-line"
    DATA_BADGE       = "mv-data-badge"
    DATA_BADGE_TEXT  = "mv-data-badge-text"
    BTN_THEME        = "mv-btn-theme"
    BTN_SESSION      = "mv-btn-session"
    BTN_EXPORT       = "mv-btn-export"
    BTN_HELP         = "mv-btn-help"

    # ── session save / restore ────────────────
    CANVAS_SESSION      = "mv-canvas-session"
    SESSION_NOTE        = "mv-session-note"
    BTN_SESSION_SAVE    = "mv-btn-session-save"
    SESSION_DL          = "mv-session-dl"
    SESSION_SNAPSHOT    = "mv-session-snapshot"      # localStorage auto-save
    SESSION_AUTOSAVE    = "mv-session-autosave"
    BTN_SESSION_RESTORE = "mv-btn-session-restore"
    SESSION_UL          = "mv-session-upload"
    BTN_SESSION_CLEAR   = "mv-btn-session-clear"
    SESSION_FEEDBACK    = "mv-session-feedback"

    # ── view rail ─────────────────────────────
    RAIL_MAP         = "mv-rail-map"
    RAIL_3D          = "mv-rail-3d"

    # ── data / lines panel ────────────────────
    LINE_PILLS       = "mv-line-pills"
    LINE_PANEL       = "mv-line-panel"
    STATION_LIST     = "mv-station-list"
    STA_LOAD_BAR     = "mv-sta-load-bar"

    # ── canvas ────────────────────────────────
    CANVAS_GRAPH     = "mv-canvas-graph"
    CANVAS_TITLE     = "mv-canvas-title"
    CANVAS_EMPTY     = "mv-canvas-empty"
    WELCOME          = "mv-welcome"
    WELCOME_CTA      = "mv-welcome-cta"

    # ── canvas toolbar (2-D map view) ─────────
    TOOLBAR_2D       = "mv-toolbar-2d"
    TB_INFO          = "mv-tb-info"
    TB_FIT           = "mv-tb-fit"
    TB_LABELS        = "mv-tb-labels"
    TB_PROFILES      = "mv-tb-profiles"
    TB_CONTOUR       = "mv-tb-contour"
    TB_BM_DARK       = "mv-tb-bm-dark"
    TB_BM_LIGHT      = "mv-tb-bm-light"
    TB_BM_SAT        = "mv-tb-bm-sat"
    TB_BM_STREET     = "mv-tb-bm-street"
    TB_BM_TOPO       = "mv-tb-bm-topo"
    TB_MARK_DEC      = "mv-tb-mark-dec"
    TB_MARK_INC      = "mv-tb-mark-inc"
    TB_MARK_VAL      = "mv-tb-mark-val"
    STORE_FIT        = "mv-store-fit"

    # ── canvas toolbar (3-D map view) ─────────
    TOOLBAR_3D       = "mv-toolbar-3d"
    TB3D_RESET       = "mv-tb3d-reset"
    TB3D_MODE_FENCE  = "mv-tb3d-mode-fence"
    TB3D_MODE_BLOCK  = "mv-tb3d-mode-block"
    TB3D_MODE_DEPTH  = "mv-tb3d-mode-depth"
    TB3D_MODE_SURFACE = "mv-tb3d-mode-surface"
    TB3D_DEPTH_FULL  = "mv-tb3d-depth-full"
    TB3D_DEPTH_500   = "mv-tb3d-depth-500"
    TB3D_DEPTH_1K    = "mv-tb3d-depth-1k"
    TB3D_DEPTH_2K    = "mv-tb3d-depth-2k"
    TB3D_TOPO        = "mv-tb3d-topo"

    # ── coordinate system ─────────────────────
    CTL_CRS_MODE     = "mv-ctl-crs-mode"
    CTL_UTM_ZONE     = "mv-ctl-utm-zone"
    CTL_UTM_HEM      = "mv-ctl-utm-hem"
    CTL_EPSG         = "mv-ctl-epsg"
    CRS_INFO         = "mv-crs-info"
    GRP_UTM          = "mv-grp-utm"
    GRP_EPSG         = "mv-grp-epsg"

    # ── right inspector ───────────────────────
    INSPECTOR        = "mv-inspector"
    CTL_COMPONENT    = "mv-ctl-component"
    CTL_QUANTITY     = "mv-ctl-quantity"
    CTL_OVERLAY      = "mv-ctl-overlay"
    CTL_FREQUENCY    = "mv-ctl-frequency"
    CTL_FREQ_LABEL   = "mv-ctl-freq-label"
    CTL_CMAP         = "mv-ctl-cmap"
    CTL_LOG          = "mv-ctl-log"
    CTL_MODE3D       = "mv-ctl-mode3d"
    CTL_LABELS       = "mv-ctl-labels"
    CTL_BASEMAP      = "mv-ctl-basemap"
    CTL_CONTOUR_LEVELS = "mv-ctl-contour-levels"
    CTL_CONTOUR_MODE = "mv-ctl-contour-mode"
    CTL_CONTOUR_ENABLE = "mv-ctl-contour-enable"
    CTL_CONTOUR_INTERP = "mv-ctl-contour-interp"
    CTL_CONTOUR_SMOOTH = "mv-ctl-contour-smooth"
    CTL_CONTOUR_RES  = "mv-ctl-contour-res"
    CTL_MARKER_SIZE  = "mv-ctl-marker-size"
    CTL_OPACITY_MAP  = "mv-ctl-opacity-map"
    CTL_PROFILES     = "mv-ctl-profiles"
    GRP_MAP          = "mv-grp-map"
    CTL_OPACITY      = "mv-ctl-opacity"
    CTL_AZIMUTH      = "mv-ctl-azimuth"
    CTL_SPACING      = "mv-ctl-spacing"
    CTL_DEPTH_LO     = "mv-ctl-depth-lo"
    CTL_DEPTH_HI     = "mv-ctl-depth-hi"
    CTL_NSLICES      = "mv-ctl-nslices"
    CTL_SURFACES     = "mv-ctl-surfaces"
    CTL_CONTOURS     = "mv-ctl-contours"
    CTL_SCALE        = "mv-ctl-scale"
    CTL_ASPECT       = "mv-ctl-aspect"
    CTL_X_UNIT       = "mv-ctl-x-unit"
    CTL_DEPTH_UNIT   = "mv-ctl-depth-unit"
    CTL_SMOOTH       = "mv-ctl-smooth"
    CTL_SECTION_RES  = "mv-ctl-section-res"
    BTN_DEPTH_FULL   = "mv-btn-depth-full"
    BTN_DEPTH_500    = "mv-btn-depth-500"
    BTN_DEPTH_1K     = "mv-btn-depth-1k"
    BTN_DEPTH_2K     = "mv-btn-depth-2k"
    CTL_RHO_LO       = "mv-ctl-rho-lo"
    CTL_RHO_HI       = "mv-ctl-rho-hi"
    BTN_RHO_ALL      = "mv-btn-rho-all"
    BTN_RHO_COND     = "mv-btn-rho-cond"
    BTN_RHO_MID      = "mv-btn-rho-mid"
    BTN_RHO_RES      = "mv-btn-rho-res"
    CTL_VMIN         = "mv-ctl-vmin"
    CTL_VMAX         = "mv-ctl-vmax"
    CTL_TOPO         = "mv-ctl-topo"
    CTL_TERRAIN      = "mv-ctl-terrain"
    CTL_SHOW_STA     = "mv-ctl-show-sta"
    CTL_STA_LABELS   = "mv-ctl-sta-labels"
    CTL_STA_SYMBOL   = "mv-ctl-sta-symbol"
    CTL_STA_SIZE     = "mv-ctl-sta-size"
    CTL_STA_COLOR    = "mv-ctl-sta-color"
    TOPO_SRC         = "mv-topo-src"
    TOPO_UPLOAD      = "mv-topo-upload"
    TOPO_UPLOAD_WRAP = "mv-topo-upload-wrap"
    TOPO_UPLOAD_STORE = "mv-topo-upload-store"
    TOPO_FETCH_WRAP  = "mv-topo-fetch-wrap"
    TOPO_API         = "mv-topo-api"
    BTN_TOPO_APPLY   = "mv-btn-topo-apply"
    TOPO_STATUS      = "mv-topo-status"
    TOPO_EXPORT_FMT  = "mv-topo-export-fmt"
    BTN_TOPO_EXPORT  = "mv-btn-topo-export"
    TOPO_EXPORT_DL   = "mv-topo-export-dl"
    GRP_STATION      = "mv-grp-station"
    GRP_3D           = "mv-grp-3d"
    STATION_INSPECT  = "mv-station-inspect"

    # ── bottom dock ───────────────────────────
    DOCK_TABLE       = "mv-dock-table"
    DOCK_TOGGLE      = "mv-dock-toggle"
    DOCK_CLOSE       = "mv-dock-close"
    DOCK_BODY        = "mv-dock-body"

    # ── load modal ────────────────────────────
    MODAL_LOAD       = "mv-modal-load"
    MODAL_TITLE      = "mv-modal-title"
    MODAL_CLOSE      = "mv-modal-close"
    MODE_BTN_REPLACE = "mv-mode-btn-replace"
    MODE_BTN_APPEND  = "mv-mode-btn-append"
    MODE_HINT        = "mv-mode-hint"
    BTN_BROWSE       = "mv-btn-browse"
    UPLOAD           = "mv-upload"
    DROP_WRAP        = "mv-drop-wrap"
    DROP_TITLE       = "mv-drop-title"
    FILE_MANAGER     = "mv-file-manager"
    PREFLIGHT        = "mv-preflight"
    DETECTED_SUMMARY = "mv-detected-summary"
    LOAD_LINE_FILTER_WRAP = "mv-load-line-filter-wrap"
    LOAD_LINE_FILTER = "mv-load-line-filter"
    LOADER_OVERLAY   = "mv-loader-overlay"
    LOADER_MSG       = "mv-loader-msg"
    FILE_COUNT       = "mv-file-count"
    BROWSE_STATUS    = "mv-browse-status"
    PROGRESS_WRAP    = "mv-progress-wrap"
    PROG_FILL        = "mv-prog-fill"
    PROG_LABEL       = "mv-prog-label"
    PROG_SUBLABEL    = "mv-prog-sublabel"
    BTN_LOAD_CONFIRM = "mv-btn-load-confirm"
    LOAD_FEEDBACK    = "mv-load-feedback"

    # ── load modal: inversion-results tab ─────
    LOAD_TABS        = "mv-load-tabs"
    BTN_INV_BROWSE   = "mv-btn-inv-browse"
    INV_FILE_COUNT   = "mv-inv-file-count"
    INV_BROWSE_STATUS = "mv-inv-browse-status"
    INV_LOADER_OVERLAY = "mv-inv-loader-overlay"
    INV_LOADER_MSG   = "mv-inv-loader-msg"
    CK_INV_KNOWN_STA = "mv-ck-inv-known-sta"
    INV_STATUS       = "mv-inv-status"
    BTN_INV_CONFIRM  = "mv-btn-inv-confirm"

    # ── export modal / download ───────────────
    EXPORT_DL        = "mv-export-dl"

    # ── sites / station settings ──────────────
    BTN_SETTINGS     = "mv-btn-settings"
    CANVAS_SETTINGS  = "mv-canvas-settings"
    SET_LINES        = "mv-set-lines"
    SET_SUMMARY      = "mv-set-summary"
    BTN_CLEAR_MASKS  = "mv-btn-clear-masks"
    BTN_MASK_HIDDEN  = "mv-btn-mask-hidden"
    SET_MASKED_LIST  = "mv-set-masked-list"

    # ── help modal ────────────────────────────
    MODAL_HELP       = "mv-modal-help"
    BTN_HELP_CLOSE   = "mv-btn-help-close"
