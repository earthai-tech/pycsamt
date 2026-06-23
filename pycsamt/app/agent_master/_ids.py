# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Component ID constants for agent_master."""

from __future__ import annotations


class IDs:
    # ── stores ────────────────────────────────
    STORE_MESSAGES   = "am-store-messages"
    STORE_EDI        = "am-store-edi"
    STORE_SESSION    = "am-store-session"
    STORE_SETTINGS   = "am-store-settings"
    STORE_JOB        = "am-store-job"
    STORE_THEME      = "am-store-theme"
    STORE_FIGS       = "am-store-figs"
    STORE_POSTPROC   = "am-store-postproc"

    # ── top-bar ───────────────────────────────
    BTN_LOAD_EDI     = "am-btn-load-edi"
    BTN_SAVE_SESSION = "am-btn-save-session"
    BTN_SETTINGS     = "am-btn-settings"
    BTN_HELP         = "am-btn-help"
    MODAL_HELP       = "am-modal-help"
    BTN_HELP_CLOSE   = "am-btn-help-close"
    HELP_CHIP        = "am-help-chip"
    BTN_THEME        = "am-theme-toggle"
    EDI_BADGE        = "am-edi-badge"
    EDI_BADGE_TEXT   = "am-edi-badge-text"

    # ── chat ──────────────────────────────────
    CHAT_WINDOW      = "am-chat-window"
    INPUT            = "am-input"
    BTN_SEND         = "am-btn-send"
    BTN_ATTACH       = "am-btn-attach"
    INTERVAL_POLL    = "am-interval-poll"
    WELCOME          = "am-welcome"

    # ── EDI offcanvas ─────────────────────────
    CANVAS_EDI       = "am-canvas-edi"
    BTN_CANVAS_CLOSE = "am-canvas-edi-close"
    UPLOAD_EDI       = "am-upload-edi"
    EDI_PATH_INPUT   = "am-edi-path"
    BTN_BROWSE       = "am-btn-browse"
    BTN_LOAD_CONFIRM = "am-btn-load-confirm"
    LINES_MODE       = "am-lines-mode"
    BTN_DETECT_LINES = "am-btn-detect-lines"
    BTN_APPLY_RENAME = "am-btn-apply-rename"
    LINES_PANEL      = "am-lines-panel"
    EDI_LOAD_STATUS  = "am-edi-load-status"
    LOAD_MODE_STORE  = "am-load-mode-store"

    # ── settings offcanvas ────────────────────
    CANVAS_SETTINGS  = "am-canvas-settings"
    PROVIDER_TABS    = "am-provider-tabs"
    KEY_CLAUDE       = "am-key-claude"
    KEY_OPENAI       = "am-key-openai"
    KEY_GEMINI       = "am-key-gemini"
    KEY_DEEPSEEK     = "am-key-deepseek"
    MODEL_CLAUDE     = "am-model-claude"
    MODEL_OPENAI     = "am-model-openai"
    MODEL_GEMINI     = "am-model-gemini"
    MODEL_DEEPSEEK   = "am-model-deepseek"
    EXPORT_FORMAT    = "am-export-format"
    ACTIVE_PROVIDER  = "am-active-provider"
    BTN_SAVE_KEYS    = "am-btn-save-keys"
    KEYS_STATUS      = "am-keys-status"
    OUTPUT_DIR       = "am-output-dir"
    LINE_REGISTRY    = "am-line-registry"

    # ── figure modal ──────────────────────────
    MODAL_FIG        = "am-modal-fig"
    MODAL_FIG_IMG    = "am-modal-fig-img"
    MODAL_FIG_TITLE  = "am-modal-fig-title"
    MODAL_FIG_KEY    = "am-modal-fig-key"
    BTN_EXPORT_FIG   = "am-btn-export-fig"
    EXPORT_DL        = "am-export-dl"

    # ── session ───────────────────────────────
    SAVE_STATUS      = "am-save-status"
    BTN_LOAD_SESSION = "am-btn-load-session"

    # ── folder upload store (JS → Python) ─────
    FOLDER_STORE     = "am-folder-store"

    # ── splash overlay ────────────────────────
    SPLASH_OVERLAY   = "am-splash-overlay"
    SPLASH_CTA       = "am-splash-cta"
    SPLASH_CARD_LOAD = "am-splash-card-load"
    SPLASH_CARD_CHAT = "am-splash-card-chat"
    SPLASH_CARD_AI   = "am-splash-card-ai"
    SPLASH_CARD_REPORT = (
        "am-splash-card-report"
    )

    # ── PINN / Hybrid parameter offcanvas ──────
    CANVAS_INV_PARAMS  = "am-canvas-inv-params"
    STORE_INV_CONFIG   = "am-store-inv-config"
    BTN_INV_PARAMS     = "am-btn-inv-params"
    INV_MODE           = "am-inv-mode"
    INV_DIM            = "am-inv-dim"
    INV_SOLVER         = "am-inv-solver"
    INV_N_LAYERS       = "am-inv-n-layers"
    INV_DEPTH_MAX      = "am-inv-depth-max"
    INV_EPOCHS         = "am-inv-epochs"
    INV_LR             = "am-inv-lr"
    INV_SMOOTH_W       = "am-inv-smooth-w"
    INV_LAT_W          = "am-inv-lat-w"
    INV_GRAPH_W        = "am-inv-graph-w"
    INV_RADIUS         = "am-inv-radius"
    INV_CHECKPOINT     = "am-inv-checkpoint"
    INV_PANEL_2D       = "am-inv-panel-2d"
    INV_PANEL_3D       = "am-inv-panel-3d"
    INV_PANEL_HYBRID   = "am-inv-panel-hybrid"
    BTN_INV_CONFIRM    = "am-btn-inv-confirm"
    INV_STATUS         = "am-inv-status"

    # ── "+" resource button ───────────────────
    BTN_PLUS      = "am-btn-plus"
    PLUS_POPOVER  = "am-plus-popover"
    PLUS_LOAD     = "am-plus-load"
    PLUS_SETTINGS = "am-plus-settings"
    PLUS_INV      = "am-plus-inv"
    PLUS_PATH     = "am-plus-path"
    PLUS_WEB      = "am-plus-web"

    # ── sidebar ───────────────────────────────
    AM_BODY          = "am-body"
    AM_MAIN          = "am-main"
    SIDEBAR          = "am-sidebar"
    BTN_SIDEBAR      = "am-btn-sidebar"
    BTN_NEW_CHAT     = "am-btn-new-chat"
    SIDEBAR_HISTORY  = "am-sidebar-history"
    SIDEBAR_FIGS     = "am-sidebar-figs"
    SIDEBAR_PINS     = "am-sidebar-pins"
    SIDEBAR_RUNS     = "am-sidebar-runs"
    STORE_HISTORY    = "am-store-history"
    STORE_PINS       = "am-store-pins"
    PIN_SCROLL_DUMMY = "am-pin-scroll-dummy"
    # segmented sidebar switcher (Chat / Session)
    SEG_CHAT         = "am-seg-chat"
    SEG_SESSION      = "am-seg-session"
    PANEL_CHAT       = "am-panel-chat"
    PANEL_SESSION    = "am-panel-session"
    BTN_NEW_SESSION  = "am-btn-new-session"
    STORE_SB_TAB     = "am-store-sb-tab"

    # ── smart param modal ─────────────────────
    MODAL_PARAMS       = "am-modal-params"
    STORE_PENDING      = "am-store-pending"
    PARAM_FORM_BODY    = "am-param-form-body"
    PARAM_MODAL_TITLE  = "am-param-modal-title"
    PARAM_MODAL_DESC   = "am-param-modal-desc"
    BTN_PARAM_RUN      = "am-btn-param-run"
    BTN_PARAM_CANCEL   = "am-btn-param-cancel"

    # ── line selector modal ───────────────
    MODAL_LINE_SEL      = "am-modal-line-sel"
    LINE_SEL_LIST       = "am-line-sel-list"
    BTN_LINE_RUN_SEL    = "am-btn-line-run-sel"
    BTN_LINE_RUN_ALL    = "am-btn-line-run-all"
    LINE_SEL_STATUS     = "am-line-sel-status"

    # ── output directory browser modal ───
    BTN_OUTPUT_BROWSE    = "am-btn-output-browse"
    MODAL_OUTPUT_BROWSE  = "am-modal-output-browse"
    OUTPUT_BROWSE_STORE  = "am-output-browse-store"
    OUTPUT_BROWSE_PATH   = "am-output-browse-path"
    OUTPUT_BROWSE_LIST   = "am-output-browse-list"
    OUTPUT_MKDIR_INPUT   = "am-output-mkdir-input"
    BTN_OUTPUT_MKDIR     = "am-btn-output-mkdir"
    BTN_OUTPUT_UP        = "am-btn-output-up"
    BTN_OUTPUT_CONFIRM   = "am-btn-output-confirm"
    OUTPUT_BROWSE_STATUS = "am-output-browse-status"

    # ── post-correction action modal ──────
    MODAL_POSTPROC      = "am-modal-postproc"
    POSTPROC_MSG        = "am-postproc-msg"
    BTN_POSTPROC_APPLY  = (
        "am-btn-postproc-apply"
    )
    BTN_POSTPROC_EXPORT = (
        "am-btn-postproc-export"
    )
    POSTPROC_PATH       = "am-postproc-path"
    BTN_POSTPROC_OK     = "am-btn-postproc-ok"
    POSTPROC_STATUS     = "am-postproc-status"
    POSTPROC_COLLAPSE   = (
        "am-postproc-collapse"
    )
