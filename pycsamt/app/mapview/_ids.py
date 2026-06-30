# -*- coding: utf-8 -*-
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
    LOAD_MODE_STORE  = "mv-load-mode-store"   # replace | append
    FOLDER_STORE     = "mv-folder-store"      # JS folder upload payload
    UPLOAD_SELECTION = "mv-upload-selection"  # staged file entries
    SOURCE_SELECTION = "mv-source-selection"

    # ── top bar ───────────────────────────────
    BTN_SIDEBAR      = "mv-btn-sidebar"
    BTN_INSPECTOR    = "mv-btn-inspector"
    BTN_LOAD         = "mv-btn-load"
    BTN_ADD_LINE     = "mv-btn-add-line"
    DATA_BADGE       = "mv-data-badge"
    DATA_BADGE_TEXT  = "mv-data-badge-text"
    BTN_THEME        = "mv-btn-theme"
    BTN_EXPORT       = "mv-btn-export"
    BTN_HELP         = "mv-btn-help"

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

    # ── canvas toolbar ────────────────────────
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
    CTL_RHO_LO       = "mv-ctl-rho-lo"
    CTL_RHO_HI       = "mv-ctl-rho-hi"
    CTL_VMIN         = "mv-ctl-vmin"
    CTL_VMAX         = "mv-ctl-vmax"
    CTL_TOPO         = "mv-ctl-topo"
    CTL_TERRAIN      = "mv-ctl-terrain"
    CTL_SHOW_STA     = "mv-ctl-show-sta"
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
    GRP_STATION      = "mv-grp-station"
    GRP_3D           = "mv-grp-3d"
    STATION_INSPECT  = "mv-station-inspect"

    # ── bottom dock ───────────────────────────
    DOCK_TABLE       = "mv-dock-table"
    DOCK_TOGGLE      = "mv-dock-toggle"
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

    # ── export modal / download ───────────────
    EXPORT_DL        = "mv-export-dl"

    # ── help modal ────────────────────────────
    MODAL_HELP       = "mv-modal-help"
    BTN_HELP_CLOSE   = "mv-btn-help-close"
