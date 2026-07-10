# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Callbacks for the TDEM analysis page.

Tab system
----------
decay     → Decay / Rho      (DECAY_PLOTS)
section   → Survey Section   (SECTION_PLOTS)
map       → Map & Overview   (MAP_PLOTS)
dashboard → Dashboard        (DASHBOARD_PLOTS)
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dash import (
    Input,
    Output,
    State,
    clientside_callback,
    ctx,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.tdem_controller import (
    DASHBOARD_PLOTS,
    DECAY_PLOTS,
    MAP_PLOTS,
    SECTION_PLOTS,
    TDEM_GROUPS,
    TDEMController,
)
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.pages.tdem import _SAMPLE
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    apply_web_light_theme,
    empty_src,
    fig_to_src,
)

# ── Singleton controller ──────────────────────────────────────────────────────
_CTRL = TDEMController()

# ── Tab / plot maps ───────────────────────────────────────────────────────────
_TDEM_TAB_IDS = ["decay", "section", "map", "dashboard"]
_TDEM_TAB_IDX = {t: i for i, t in enumerate(_TDEM_TAB_IDS)}

_TAB_TO_PLOTS = {
    "decay":     DECAY_PLOTS,
    "section":   SECTION_PLOTS,
    "map":       MAP_PLOTS,
    "dashboard": DASHBOARD_PLOTS,
}
_TAB_TO_GROUP = {
    "decay":     "Decay / Rho",
    "section":   "Survey Section",
    "map":       "Map & Overview",
    "dashboard": "Dashboard",
}
_TAB_ICONS = {
    "decay":     "bi-activity",
    "section":   "bi-layers",
    "map":       "bi-map",
    "dashboard": "bi-grid-3x3",
}

# Unified meta: class_name → (has_ax, data_key, label)
_TDEM_META: dict[str, tuple] = {}
for _grp, _plots in TDEM_GROUPS:
    for _lbl, _cls, _has_ax, _dk in _plots:
        _TDEM_META[_cls] = (_has_ax, _dk, _lbl)

# Conditional parameter panel sets
_SOUNDING_PLOTS     = {"PlotDecayCurve", "PlotTransformedRho"}
_SECTION_CMAP_PLOTS = {
    "PlotTEMAVGSection", "PlotTEMZSection", "PlotGateProfile",
    "PlotSurveyMap", "PlotSurveyOverview",
}
_GATE_PLOTS         = {"PlotGateProfile"}

_FIGSIZE_MAP = {
    "8x4":  (8,  4),
    "10x5": (10, 5),
    "13x5": (13, 5),
    "10x7": (10, 7),
}


def register_tdem(app) -> None:

    # ── 1. Browse folder (native OS dialog, properly sized) ───────────────────
    @app.callback(
        Output(IDs.TDEM_FOLDER,    "value"),
        Input(IDs.BTN_TDEM_BROWSE, "n_clicks"),
        prevent_initial_call=True,
    )
    def browse_tdem_folder(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        try:
            import tkinter as tk
            from tkinter import filedialog
            from tkinter import font as tkfont

            root = tk.Tk()
            root.withdraw()

            # Scale the Tk UI so the dialog is legible on HiDPI / Linux
            try:
                root.tk.call("tk", "scaling", 2.0)
            except Exception:
                pass

            # Increase font sizes inside the dialog
            for name in ("TkDefaultFont", "TkTextFont", "TkHeadingFont",
                         "TkMenuFont", "TkFixedFont"):
                try:
                    f = tkfont.nametofont(name)
                    f.configure(size=max(12, f.cget("size")))
                except Exception:
                    pass

            root.attributes("-topmost", True)
            root.update()

            folder = filedialog.askdirectory(
                parent=root,
                title="Select TDEM data folder  (.AVG / .Z / .LOG files)",
                mustexist=True,
            )
            root.destroy()
            return folder if folder else no_update
        except Exception:
            raise PreventUpdate

    # ── 2. Tab button clicks → TDEM_ACTIVE_TAB store ─────────────────────────
    @app.callback(
        Output(IDs.TDEM_ACTIVE_TAB, "data"),
        [Input(f"tdem-tab-{t}", "n_clicks") for t in _TDEM_TAB_IDS],
        prevent_initial_call=True,
    )
    def switch_tdem_tab(*_clicks):
        tid = ctx.triggered_id
        if not tid:
            return no_update
        return tid.replace("tdem-tab-", "")

    # ── 3. Store → panel visibility + button classes (clientside) ─────────────
    clientside_callback(
        """
        function(active_tab) {
            if (!active_tab)
                return Array(8).fill(window.dash_clientside.no_update);
            var tabs    = ["decay","section","map","dashboard"];
            var flex    = {display:"flex",flexDirection:"column",
                           height:"100%",minHeight:"0",flex:"1"};
            var none    = {display:"none"};
            var styles  = tabs.map(function(t){ return t===active_tab ? flex : none; });
            var classes = tabs.map(function(t){
                return "prof-tab-btn" + (t===active_tab ? " active" : "");
            });
            return styles.concat(classes);
        }
        """,
        Output("tdem-panel-decay",     "style"),
        Output("tdem-panel-section",   "style"),
        Output("tdem-panel-map",       "style"),
        Output("tdem-panel-dashboard", "style"),
        Output("tdem-tab-decay",       "className"),
        Output("tdem-tab-section",     "className"),
        Output("tdem-tab-map",         "className"),
        Output("tdem-tab-dashboard",   "className"),
        Input(IDs.TDEM_ACTIVE_TAB, "data"),
        prevent_initial_call=False,
    )

    # ── 4. Tab change → plot dropdown ─────────────────────────────────────────
    @app.callback(
        Output(IDs.TDEM_PLOT, "options"),
        Output(IDs.TDEM_PLOT, "value"),
        Input(IDs.TDEM_ACTIVE_TAB, "data"),
    )
    def sync_tdem_plots(active_tab):
        if not active_tab:
            raise PreventUpdate
        plots = _TAB_TO_PLOTS.get(active_tab, DECAY_PLOTS)
        opts  = [{"label": lbl, "value": cls} for lbl, cls, *_ in plots]
        return opts, (opts[0]["value"] if opts else None)

    # ── 5. Plot → parameter panel visibility ──────────────────────────────────
    @app.callback(
        Output("tdem-sounding-section", "style"),
        Output("tdem-cmap-section",     "style"),
        Output("tdem-gate-section",     "style"),
        Input(IDs.TDEM_PLOT, "value"),
    )
    def update_param_panels(class_name):
        none = {"display": "none"}
        blk  = {"display": "block"}
        if not class_name:
            return none, none, none
        return (
            blk if class_name in _SOUNDING_PLOTS     else none,
            blk if class_name in _SECTION_CMAP_PLOTS else none,
            blk if class_name in _GATE_PLOTS         else none,
        )

    # ── 6. Context info bar ───────────────────────────────────────────────────
    @app.callback(
        Output(IDs.TDEM_CTX_BAR,    "children"),
        Input(IDs.TDEM_ACTIVE_TAB,   "data"),
        Input(IDs.TDEM_PLOT,         "value"),
        Input(IDs.STORE_TDEM_LOADED, "data"),
        prevent_initial_call=False,
    )
    def update_ctx_bar(active_tab, class_name, tdem_loaded):
        group = _TAB_TO_GROUP.get(active_tab or "decay", "Decay / Rho")
        icon  = _TAB_ICONS.get(active_tab or "decay", "bi-activity")
        label = (_TDEM_META.get(class_name, (None, None, ""))[2]
                 if class_name else "")
        parts = [
            html.Span(
                [html.I(className=f"bi {icon} me-1"), group],
                className="adv-info-group",
            ),
        ]
        if label:
            parts += [html.Span("·", className="adv-info-sep"),
                      html.Span(label, className="adv-info-plot")]
        if tdem_loaded:
            parts += [html.Span("·", className="adv-info-sep"),
                      html.Span("data loaded", className="adv-info-stat")]
        else:
            parts.append(html.Span(
                " · Browse and load a TDEM folder to begin",
                style={"color": "var(--sub0)", "fontSize": "10.5px"},
            ))
        return parts

    # ── 7. Load folder (path or sample button) ────────────────────────────────
    @app.callback(
        Output(IDs.TDEM_INFO,         "children"),
        Output(IDs.STORE_TDEM_LOADED, "data"),
        Output(IDs.TDEM_FOLDER,       "value",   allow_duplicate=True),
        Output(IDs.TOAST_ERROR,       "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY,        "children",allow_duplicate=True),
        Input(IDs.BTN_TDEM_LOAD,   "n_clicks"),
        Input(IDs.BTN_TDEM_SAMPLE, "n_clicks"),
        State(IDs.TDEM_FOLDER,     "value"),
        prevent_initial_call=True,
    )
    def load_tdem(n_load, n_sample, folder):
        trigger = ctx.triggered_id
        if trigger == IDs.BTN_TDEM_SAMPLE:
            folder = _SAMPLE
        if not folder:
            raise PreventUpdate
        folder = folder.strip()
        try:
            ok = _CTRL.load_folder(folder)
            if ok:
                return _fmt_info(_CTRL.summary), True, folder, False, ""
            return ("No TDEM files found in that folder.",
                    False, folder, True, "TDEM: folder empty")
        except Exception as exc:
            return f"Error: {exc}", False, folder, True, f"TDEM load: {exc}"

    # ── 8. Generate → active tab's image ─────────────────────────────────────
    @app.callback(
        Output(IDs.IMG_TDEM_DECAY,   "src"),
        Output(IDs.IMG_TDEM_SECTION, "src"),
        Output(IDs.IMG_TDEM_MAP,     "src"),
        Output(IDs.IMG_TDEM_DASH,    "src"),
        Output(IDs.TDEM_SPINNER,     "children"),
        Output(IDs.TOAST_ERROR,      "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,       "children", allow_duplicate=True),
        Input(IDs.BTN_TDEM_RUN,      "n_clicks"),
        Input(IDs.TDEM_PLOT,         "value"),
        Input(IDs.NAV_SECTION,       "data"),
        State(IDs.TDEM_ACTIVE_TAB,   "data"),
        State(IDs.STORE_TDEM_LOADED, "data"),
        State(IDs.STORE_THEME,       "data"),
        State(IDs.TDEM_FIGSIZE,      "value"),
        State(IDs.TDEM_CMAP,         "value"),
        State(IDs.TDEM_SOUNDING_IDX, "value"),
        State(IDs.TDEM_GATE_IDX,     "value"),
        prevent_initial_call=True,
    )
    def run_tdem(n_clicks, class_name, nav_section,
                 active_tab, tdem_loaded, theme,
                 figsize, cmap, sounding_idx, gate_idx):
        trigger = ctx.triggered_id
        if trigger in (IDs.TDEM_PLOT, IDs.NAV_SECTION):
            if nav_section != "tdem" or not tdem_loaded:
                return (no_update,) * 7
        if not class_name:
            raise PreventUpdate
        dark = (theme or "dark") == "dark"
        if not _CTRL.is_loaded():
            _e = empty_src()
            return _e, _e, _e, _e, "", True, "TDEM: load a folder first"
        has_ax, data_key, _ = _TDEM_META.get(class_name, (True, "survey", class_name))
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()
        _CTRL.dark = dark
        try:
            w, h   = _FIGSIZE_MAP.get(figsize or "10x5", (10, 5))
            fig    = plt.figure(figsize=(w, h))
            _extra = {}
            if class_name in _SECTION_CMAP_PLOTS and cmap:
                _extra["cmap"] = cmap
            if class_name in _SOUNDING_PLOTS and sounding_idx is not None:
                _extra["sounding_index"] = int(sounding_idx)
            if class_name in _GATE_PLOTS and gate_idx is not None:
                _extra["gate_index"] = int(gate_idx)
            result = _ctrl_draw(class_name, has_ax, data_key, fig, **_extra)
            src    = fig_to_src(result if result is not None else fig)
            plt.close("all")
            outputs = [no_update] * 4
            outputs[_TDEM_TAB_IDX.get(active_tab or "decay", 0)] = src
            return (*outputs, "", False, "")
        except Exception as exc:
            plt.close("all")
            return (no_update, no_update, no_update, no_update,
                    "", True, f"{class_name}: {exc}")



# ── Private helpers ───────────────────────────────────────────────────────────

def _ctrl_draw(class_name: str, has_ax: bool, data_key: str, fig, **kwargs):
    """Try draw_ext(**kwargs) if available, fall back to draw()."""
    if kwargs and hasattr(_CTRL, "draw_ext"):
        try:
            return _CTRL.draw_ext(class_name, has_ax, data_key, fig, **kwargs)
        except TypeError:
            pass
    return _CTRL.draw(class_name, has_ax, data_key, fig)


def _fmt_info(summary: str) -> list:
    lines = (summary or "").strip().split("\n")
    return [html.Div(ln, className="tdem-info-line") for ln in lines if ln]
