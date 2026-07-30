# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/profile.py — Profile View page callbacks.

Tab layout (short IDs stored in PROF_ACTIVE_TAB store)
-------------------------------------------------------
rho-phi   — ρₐ / φ response curves         → matplotlib PNG  (IMG_RHO_PHI.src)
rho-ps    — ρₐ pseudosection (per-comps)   → Plotly subplots (IMG_RHO_PS.figure)
phi-ps    — φ  pseudosection (per-comps)   → Plotly subplots (IMG_PHASE_PS.figure)
pt        — Phase tensor pseudosection      → matplotlib PNG  (IMG_PT.src)
tipper    — Tipper components               → matplotlib PNG  (IMG_TIPPER.src)
section   — 2-D inversion section           → matplotlib PNG  (IMG_SECTION.src)
pub       — Publication-ready MT response   → matplotlib PNG  (IMG_PUB.src)
"""

from __future__ import annotations

from dash import (
    Input,
    Output,
    State,
    clientside_callback,
    ctx,
    no_update,
)

from pycsamt.app.desktop.controllers.plot_controller import (
    PlotController,
)
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    apply_web_light_theme,
    build_multi_pseudosection,
    fig_to_src,
    filter_sites_by_lines,
    no_active_lines_src,
)

# Ordered list matching layout.py _TABS
_TAB_IDS = ["rho-phi", "rho-ps", "phi-ps", "pt", "tipper", "section", "pub"]
_TAB_IDX = {t: i for i, t in enumerate(_TAB_IDS)}

_NO_SRC = no_update
_NO_FIG = no_update
_SKIP = (
    _NO_SRC,
    _NO_SRC,
    _NO_SRC,  # rho-phi, tipper, pt  src
    _NO_FIG,
    _NO_FIG,  # rho-ps, phi-ps  figure
    _NO_SRC,
    _NO_SRC,  # section, pub  src
    no_update,
    no_update,  # ps-comp bar headers
    False,
    "",  # toast
)

_COMP_CLS = {"xy": "te", "yx": "tm", "xx": "diag", "yy": "diag"}
_COMP_LBL = {"xy": "Z_XY", "yx": "Z_YX", "xx": "Z_XX", "yy": "Z_YY"}
_PHASE_MAP = {
    "0_90": (0.0, 90.0),
    "-90_90": (-90.0, 90.0),
    "-180_180": (-180.0, 180.0),
    "-360_360": (-360.0, 360.0),
}


def _comp_bar(comps, plot_type: str):
    """Component-info bar children for a pseudosection panel."""
    from dash import html

    prefix = "ρₐ" if plot_type == "rho" else "φ"
    pills = [
        html.Span(
            _COMP_LBL.get(c, c.upper()),
            className=f"prof-ps-pill {_COMP_CLS.get(c, '')}".strip(),
        )
        for c in comps
    ]
    return [
        html.Span(
            f"{prefix} sections:",
            style={
                "color": "var(--sub0)",
                "marginRight": "8px",
                "fontSize": "10.5px",
            },
        ),
        *pills,
    ]


def register_profile(app) -> None:

    # ── 1. Tab button clicks → PROF_ACTIVE_TAB store ─────────────────────────
    @app.callback(
        Output(IDs.PROF_ACTIVE_TAB, "data"),
        [Input(f"prof-tab-{t}", "n_clicks") for t in _TAB_IDS],
        prevent_initial_call=True,
    )
    def switch_tab(*_clicks):
        tid = ctx.triggered_id
        if not tid:
            return no_update
        return tid.replace("prof-tab-", "")

    # ── 2. Store → panel visibility + button classes (clientside, instant) ────
    clientside_callback(
        """
        function(active_tab) {
            if (!active_tab) return Array(14).fill(window.dash_clientside.no_update);
            var tabs = ["rho-phi","rho-ps","phi-ps","pt","tipper","section","pub"];
            var flex = {display:"flex",flexDirection:"column",
                        height:"100%",minHeight:"0",flex:"1"};
            var none = {display:"none"};
            var styles  = tabs.map(function(t){ return t===active_tab ? flex : none; });
            var classes = tabs.map(function(t){
                return "prof-tab-btn" + (t===active_tab ? " active" : "");
            });
            return styles.concat(classes);
        }
        """,
        Output("prof-panel-rho-phi", "style"),
        Output("prof-panel-rho-ps", "style"),
        Output("prof-panel-phi-ps", "style"),
        Output("prof-panel-pt", "style"),
        Output("prof-panel-tipper", "style"),
        Output("prof-panel-section", "style"),
        Output("prof-panel-pub", "style"),
        Output("prof-tab-rho-phi", "className"),
        Output("prof-tab-rho-ps", "className"),
        Output("prof-tab-phi-ps", "className"),
        Output("prof-tab-pt", "className"),
        Output("prof-tab-tipper", "className"),
        Output("prof-tab-section", "className"),
        Output("prof-tab-pub", "className"),
        Input(IDs.PROF_ACTIVE_TAB, "data"),
        prevent_initial_call=False,
    )

    # ── 3. Publication back button → rho-phi ─────────────────────────────────
    @app.callback(
        Output(IDs.PROF_ACTIVE_TAB, "data", allow_duplicate=True),
        Input("prof-pub-back-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def pub_back(n):
        return "rho-phi" if n else no_update

    # ── 4. Main render callback ───────────────────────────────────────────────
    @app.callback(
        Output(IDs.IMG_RHO_PHI, "src"),
        Output(IDs.IMG_TIPPER, "src"),
        Output(IDs.IMG_PT, "src"),
        Output(IDs.IMG_RHO_PS, "figure"),
        Output(IDs.IMG_PHASE_PS, "figure"),
        Output(IDs.IMG_SECTION, "src"),
        Output(IDs.IMG_PUB, "src"),
        Output("prof-ps-rho-header", "children"),
        Output("prof-ps-phi-header", "children"),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Input(IDs.PROF_ACTIVE_TAB, "data"),
        Input(IDs.STORE_SELECTION, "data"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.STORE_THEME, "data"),
        Input(IDs.PROFILE_PAGE_ERRBAR, "value"),
        Input("profile-page-comps", "value"),
        Input("profile-page-refresh", "n_clicks"),
        Input(IDs.PROFILE_PAGE_PHASE_RANGE, "value"),
        State(IDs.FREQ_SLIDER, "value"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def update_profile(
        active_tab,
        selection,
        store_data,
        theme,
        show_errbar,
        components,
        _refresh,
        phase_range,
        freq_range,
        session_id,
        active_lines_store,
    ):
        if not store_data or not session_id:
            return _SKIP

        sites = cache_get(session_id)
        if sites is None:
            return _SKIP

        dark = (theme or "dark") == "dark"

        # Active-lines filter
        _als = active_lines_store or {}
        active = _als.get("active", _als.get("all"))
        if active is not None:
            records = (store_data or {}).get("station_records", [])
            filtered = filter_sites_by_lines(sites, records, active)
            if filtered is None:
                warn = no_active_lines_src(dark)
                return (
                    warn,
                    warn,
                    warn,
                    no_update,
                    no_update,
                    warn,
                    warn,
                    no_update,
                    no_update,
                    False,
                    "",
                )
            sites = filtered

        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        idx = _TAB_IDX.get(active_tab, 0)
        comps = tuple(components) if components else ("xy", "yx")

        T_min, T_max = None, None
        if freq_range and len(freq_range) == 2:
            T_min = 10 ** freq_range[0]
            T_max = 10 ** freq_range[1]

        phase_ylim = _PHASE_MAP.get(phase_range or "auto")

        try:
            # ── Plotly pseudosection tabs ──────────────────────────────────────
            if idx == 1:
                fig = build_multi_pseudosection(
                    sites, list(comps), "rho", dark
                )
                hdr = _comp_bar(comps, "rho")
                return (
                    _NO_SRC,
                    _NO_SRC,
                    _NO_SRC,
                    fig,
                    _NO_FIG,
                    _NO_SRC,
                    _NO_SRC,
                    hdr,
                    no_update,
                    False,
                    "",
                )

            if idx == 2:
                fig = build_multi_pseudosection(
                    sites, list(comps), "phi", dark
                )
                hdr = _comp_bar(comps, "phi")
                return (
                    _NO_SRC,
                    _NO_SRC,
                    _NO_SRC,
                    _NO_FIG,
                    fig,
                    _NO_SRC,
                    _NO_SRC,
                    no_update,
                    hdr,
                    False,
                    "",
                )

            # ── matplotlib PNG tabs ────────────────────────────────────────────
            import matplotlib.pyplot as plt

            ctrl = PlotController()
            ctrl.set_sites(sites)
            ctrl.dark = dark
            ctrl.set_station((selection or {}).get("station_id"))
            ctrl.set_show_errbar(bool(show_errbar))
            ctrl.set_components(comps)
            if T_min is not None:
                ctrl.set_period_range(T_min, T_max)
            ctrl.set_phase_ylim(*(phase_ylim or (None, None)))

            if idx == 0:
                fig_mpl = plt.figure(figsize=(9, 5))
                ctrl.draw_rho_phi(fig_mpl)
                src = fig_to_src(fig_mpl)
                plt.close("all")
                return (
                    src,
                    _NO_SRC,
                    _NO_SRC,
                    _NO_FIG,
                    _NO_FIG,
                    _NO_SRC,
                    _NO_SRC,
                    no_update,
                    no_update,
                    False,
                    "",
                )

            if idx == 3:
                fig_mpl, ax = plt.subplots(figsize=(9, 4))
                ctrl.draw_phase_tensor(ax)
                src = fig_to_src(fig_mpl)
                plt.close("all")
                return (
                    _NO_SRC,
                    _NO_SRC,
                    src,
                    _NO_FIG,
                    _NO_FIG,
                    _NO_SRC,
                    _NO_SRC,
                    no_update,
                    no_update,
                    False,
                    "",
                )

            if idx == 4:
                fig_mpl, ax = plt.subplots(figsize=(9, 4))
                ctrl.draw_tipper(ax)
                src = fig_to_src(fig_mpl)
                plt.close("all")
                return (
                    _NO_SRC,
                    src,
                    _NO_SRC,
                    _NO_FIG,
                    _NO_FIG,
                    _NO_SRC,
                    _NO_SRC,
                    no_update,
                    no_update,
                    False,
                    "",
                )

            if idx == 5:
                src = _draw_section(ctrl)
                return (
                    _NO_SRC,
                    _NO_SRC,
                    _NO_SRC,
                    _NO_FIG,
                    _NO_FIG,
                    src,
                    _NO_SRC,
                    no_update,
                    no_update,
                    False,
                    "",
                )

            if idx == 6:
                fig_mpl = plt.figure(figsize=(14, 8))
                ctrl.draw_publication_view(fig_mpl)
                src = fig_to_src(fig_mpl)
                plt.close("all")
                return (
                    _NO_SRC,
                    _NO_SRC,
                    _NO_SRC,
                    _NO_FIG,
                    _NO_FIG,
                    _NO_SRC,
                    src,
                    no_update,
                    no_update,
                    False,
                    "",
                )

        except Exception as exc:
            import traceback

            msg = f"Profile render failed ({active_tab}):\n{exc}\n{traceback.format_exc()}"
            return (
                _NO_SRC,
                _NO_SRC,
                _NO_SRC,
                _NO_FIG,
                _NO_FIG,
                _NO_SRC,
                _NO_SRC,
                no_update,
                no_update,
                True,
                msg,
            )

        return _SKIP

    # ── 5. Station prev / next navigation ─────────────────────────────────────
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
        new_idx = (
            i - 1 if ctx.triggered_id == IDs.PROFILE_PAGE_PREV else i + 1
        ) % len(ids)
        return {"station_id": ids[new_idx]}

    # ── 6. Export current tab figure ──────────────────────────────────────────
    @app.callback(
        Output("profile-page-download", "data"),
        Input("profile-page-export", "n_clicks"),
        State(IDs.PROF_ACTIVE_TAB, "data"),
        State(IDs.IMG_RHO_PHI, "src"),
        State(IDs.IMG_TIPPER, "src"),
        State(IDs.IMG_PT, "src"),
        State(IDs.IMG_SECTION, "src"),
        State(IDs.IMG_PUB, "src"),
        prevent_initial_call=True,
    )
    def export_profile_tab(
        n, active_tab, src_rho_phi, src_tipper, src_pt, src_section, src_pub
    ):
        if not n:
            return no_update
        src_map = {
            "rho-phi": (src_rho_phi, "pycsamt_rho_phi.png"),
            "tipper": (src_tipper, "pycsamt_tipper.png"),
            "pt": (src_pt, "pycsamt_phase_tensor.png"),
            "section": (src_section, "pycsamt_2d_section.png"),
            "pub": (src_pub, "pycsamt_publication.png"),
        }
        if active_tab not in src_map:
            return no_update
        src, filename = src_map[active_tab]
        if not src or not src.startswith("data:image"):
            return no_update
        _, b64 = src.split(",", 1)
        return dict(
            content=b64, filename=filename, base64=True, type="image/png"
        )


# ── Private helpers ───────────────────────────────────────────────────────────


def _draw_section(ctrl: PlotController) -> str:
    import matplotlib.pyplot as plt

    from pycsamt.app.web.utils import fig_to_src

    fig_mpl, ax = plt.subplots(figsize=(11, 4))
    try:
        ctrl.draw_2d_section(ax)
    except AttributeError:
        ax.text(
            0.5,
            0.5,
            "2D section not available\n(run inversion first)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )
        ax.set_axis_off()
    except Exception as exc:
        ax.text(
            0.5,
            0.5,
            f"Section error:\n{exc}",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=10,
        )
        ax.set_axis_off()
    src = fig_to_src(fig_mpl)
    plt.close("all")
    return src
