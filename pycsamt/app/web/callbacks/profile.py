# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/profile.py — MT-curve / pseudosection profile tab renders.

Tab breakdown
-------------
tab-rho-phi  (idx 0) — ρₐ / φ response curves     → matplotlib PNG   (IMG_RHO_PHI.src)
tab-rho-ps   (idx 1) — ρₐ pseudosection heatmap    → Plotly figure    (IMG_RHO_PS.figure)
tab-phase-ps (idx 2) — φ  pseudosection heatmap    → Plotly figure    (IMG_PHASE_PS.figure)
tab-tipper   (idx 3) — Tipper components           → matplotlib PNG   (IMG_TIPPER.src)
tab-pt       (idx 4) — Phase tensor pseudosection  → matplotlib PNG   (IMG_PT.src)
tab-qc       (idx 5) — QC dashboard                → handled by qc.py
"""

from __future__ import annotations

import traceback

from dash import ctx, Input, Output, State, no_update

from pycsamt.app.desktop.controllers.plot_controller import PlotController
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    apply_web_light_theme,
    build_plotly_pseudosection,
    fig_to_src,
)

_TAB_IDX = {
    "tab-rho-phi":  0,
    "tab-rho-ps":   1,
    "tab-phase-ps": 2,
    "tab-tipper":   3,
    "tab-pt":       4,
    "tab-section":  5,
    # tab-qc is excluded — handled by callbacks/qc.py
}

_NO_SRC = no_update
_NO_FIG = no_update
_SKIP   = (_NO_SRC, _NO_SRC, _NO_SRC, _NO_FIG, _NO_FIG, _NO_SRC, False, "")


def register_profile(app) -> None:
    @app.callback(
        # matplotlib PNG outputs
        Output(IDs.IMG_RHO_PHI,  "src"),
        Output(IDs.IMG_TIPPER,   "src"),
        Output(IDs.IMG_PT,       "src"),
        # Plotly heatmap outputs
        Output(IDs.IMG_RHO_PS,   "figure"),
        Output(IDs.IMG_PHASE_PS, "figure"),
        # 2D section PNG
        Output(IDs.IMG_SECTION,  "src"),
        # Error toast
        Output(IDs.TOAST_ERROR,  "is_open",   allow_duplicate=True),
        Output(IDs.TOAST_BODY,   "children",  allow_duplicate=True),
        # Triggers
        Input(IDs.PROFILE_TABS,         "active_tab"),
        Input(IDs.STORE_SELECTION,      "data"),
        Input(IDs.STORE_DATA,           "data"),
        Input(IDs.STORE_THEME,          "data"),
        # Sidebar controls
        State(IDs.PROFILE_PAGE_ERRBAR,  "value"),
        State(IDs.FREQ_SLIDER,          "value"),
        State(IDs.SESSION_ID,           "data"),
        prevent_initial_call=True,
    )
    def update_profile(active_tab, selection, store_data, theme,
                       show_errbar, freq_range, session_id):
        if not store_data or not session_id:
            return _SKIP

        # Skip if the QC tab is active (qc.py handles that)
        if active_tab == "tab-qc" and ctx.triggered_id == IDs.PROFILE_TABS:
            return _SKIP

        sites = cache_get(session_id)
        if sites is None:
            return _SKIP

        dark = (theme or "dark") == "dark"
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        # Selection change always re-renders the primary ρₐ/φ tab (idx 0)
        idx = 0 if ctx.triggered_id == IDs.STORE_SELECTION else _TAB_IDX.get(active_tab, 1)

        # Period range from slider (log10 values → actual period limits)
        T_min, T_max = None, None
        if freq_range and len(freq_range) == 2:
            import math
            T_min = 10 ** freq_range[0]
            T_max = 10 ** freq_range[1]

        try:
            # ── Plotly pseudosection tabs (1 & 2) ────────────────────────
            if idx == 1:
                fig = build_plotly_pseudosection(sites, "rho_xy", dark=dark)
                return _NO_SRC, _NO_SRC, _NO_SRC, fig, _NO_FIG, _NO_SRC, False, ""

            if idx == 2:
                fig = build_plotly_pseudosection(sites, "phi_xy", dark=dark)
                return _NO_SRC, _NO_SRC, _NO_SRC, _NO_FIG, fig, _NO_SRC, False, ""

            # ── Matplotlib PNG tabs (0, 3, 4, 5) ─────────────────────────
            import matplotlib.pyplot as plt

            ctrl = PlotController()
            ctrl.set_sites(sites)
            ctrl.dark = dark
            ctrl.set_station((selection or {}).get("station_id"))
            if T_min is not None:
                ctrl.set_period_range(T_min, T_max)

            if idx == 0:
                fig_mpl = plt.figure(figsize=(9, 5))
                _draw_rho_phi(ctrl, fig_mpl, show_errbar=bool(show_errbar))
                src = fig_to_src(fig_mpl)
                return src, _NO_SRC, _NO_SRC, _NO_FIG, _NO_FIG, _NO_SRC, False, ""

            if idx == 3:
                fig_mpl, ax = plt.subplots(figsize=(9, 4))
                ctrl.draw_tipper(ax)
                return _NO_SRC, fig_to_src(fig_mpl), _NO_SRC, _NO_FIG, _NO_FIG, _NO_SRC, False, ""

            if idx == 4:
                fig_mpl, ax = plt.subplots(figsize=(9, 4))
                ctrl.draw_phase_tensor(ax)
                return _NO_SRC, _NO_SRC, fig_to_src(fig_mpl), _NO_FIG, _NO_FIG, _NO_SRC, False, ""

            if idx == 5:
                src = _draw_section(ctrl, dark)
                return _NO_SRC, _NO_SRC, _NO_SRC, _NO_FIG, _NO_FIG, src, False, ""

        except Exception as exc:
            msg = f"Profile render failed (tab {active_tab}):\n{exc}"
            return _NO_SRC, _NO_SRC, _NO_SRC, _NO_FIG, _NO_FIG, _NO_SRC, True, msg

        return _SKIP

    # ── Station prev / next navigation ───────────────────────────────────────
    @app.callback(
        Output(IDs.STORE_SELECTION, "data", allow_duplicate=True),
        Input(IDs.PROFILE_PAGE_PREV, "n_clicks"),
        Input(IDs.PROFILE_PAGE_NEXT, "n_clicks"),
        State(IDs.PROFILE_PAGE_STATION, "value"),
        State(IDs.PROFILE_PAGE_STATION, "options"),
        prevent_initial_call=True,
    )
    def navigate_station(prev_clicks, next_clicks, current_id, options):
        if not options or not ctx.triggered_id:
            return no_update
        ids = [o["value"] for o in options]
        try:
            i = ids.index(current_id)
        except ValueError:
            return no_update
        if ctx.triggered_id == IDs.PROFILE_PAGE_PREV:
            new_idx = (i - 1) % len(ids)
        else:
            new_idx = (i + 1) % len(ids)
        return {"station_id": ids[new_idx]}

    # ── Export current tab figure ─────────────────────────────────────────────
    @app.callback(
        Output("profile-page-download", "data"),
        Input("profile-page-export", "n_clicks"),
        State(IDs.PROFILE_TABS, "active_tab"),
        State(IDs.IMG_RHO_PHI,  "src"),
        State(IDs.IMG_TIPPER,   "src"),
        State(IDs.IMG_PT,       "src"),
        State(IDs.IMG_SECTION,  "src"),
        prevent_initial_call=True,
    )
    def export_profile_tab(n, active_tab, src_rho_phi, src_tipper, src_pt, src_section):
        if not n:
            return no_update
        src_map = {
            "tab-rho-phi":  (src_rho_phi,  "pycsamt_rho_phi.png"),
            "tab-tipper":   (src_tipper,   "pycsamt_tipper.png"),
            "tab-pt":       (src_pt,       "pycsamt_pt.png"),
            "tab-section":  (src_section,  "pycsamt_section.png"),
        }
        if active_tab not in src_map:
            return no_update
        src, filename = src_map[active_tab]
        if not src or not src.startswith("data:image"):
            return no_update
        import base64
        _, b64 = src.split(",", 1)
        return dict(content=b64, filename=filename, base64=True, type="image/png")


# ── Private helpers ───────────────────────────────────────────────────────────

def _draw_rho_phi(ctrl: "PlotController", fig_mpl, show_errbar: bool = False):
    """Call ctrl.draw_rho_phi, passing error-bar flag when supported."""
    try:
        ctrl.draw_rho_phi(fig_mpl, show_errbar=show_errbar)
    except TypeError:
        ctrl.draw_rho_phi(fig_mpl)


def _draw_section(ctrl: "PlotController", dark: bool) -> str:
    """Render the 2D resistivity section; return a data-URI PNG."""
    import matplotlib.pyplot as plt
    from pycsamt.app.web.utils import fig_to_src, empty_src

    fig_mpl, ax = plt.subplots(figsize=(11, 4))
    try:
        ctrl.draw_2d_section(ax)
    except AttributeError:
        ax.text(0.5, 0.5, "2D section not available\n(run inversion first)",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_axis_off()
    except Exception as exc:
        ax.text(0.5, 0.5, f"Section error:\n{exc}",
                ha="center", va="center", transform=ax.transAxes, fontsize=10)
        ax.set_axis_off()
    return fig_to_src(fig_mpl)
