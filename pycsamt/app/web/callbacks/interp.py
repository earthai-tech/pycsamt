# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Geological Interpretation page."""
from __future__ import annotations

import base64
import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dash import (
    Input,
    Output,
    State,
    clientside_callback,
    ctx,
    dcc,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.interp_controller import (
    WORKFLOW_CATALOGUE,
    InterpController,
)
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    empty_src,
    fig_to_src,
)

_CTRL = InterpController()

# ── Category metadata (must mirror pages/interpretation.py) ──────────────────
# (slug, full-name, short-label, bootstrap-icon)
_INTERP_CATS = [
    ("setup",  "Setup & Model",        "Setup",    "bi-diagram-3"),
    ("geo",    "Geological",           "Geology",  "bi-layers"),
    ("hydro",  "Hydrology",            "Hydrology","bi-droplet-half"),
    ("field",  "Field Constraints",    "Field",    "bi-pin-map"),
    ("emdiag", "EM Diagnostics",       "EM Diag",  "bi-activity"),
    ("unc",    "Uncertainty (MC)",     "Uncert.",  "bi-bar-chart-steps"),
    ("adv",    "Advanced Plots",       "Advanced", "bi-graph-up-arrow"),
    ("fus",    "Fusion & Time-Lapse",  "Fusion",   "bi-clock-history"),
    ("exp",    "Export",               "Export",   "bi-box-arrow-up"),
]

_DEFAULT_CAT_SLUG = "setup"
_CAT_BY_SLUG      = {slug: full for slug, full, _, _ in _INTERP_CATS}
_CAT_SLUG_LIST    = [slug for slug, _, _, _ in _INTERP_CATS]
_CAT_SLUG_IDX     = {slug: i for i, slug in enumerate(_CAT_SLUG_LIST)}
_CAT_ICON         = {slug: icon for slug, _, _, icon in _INTERP_CATS}

# ── Figsize presets ───────────────────────────────────────────────────────────
_FIGSIZES = {
    "auto":   (11, 5),
    "wide":   (15, 4.5),
    "square": (7,  7),
    "tall":   (7,  10),
    "large":  (14, 7),
}

_CAT_PREREQS = {
    "Geological":        ("run_geological",  "Run Geological Classification",
                          "Classifies lithology from resistivity and rock DB. "
                          "Required for: Stratigraphic log, Fence diagram, Calibrated model."),
    "Hydrology":         ("run_hydro",       "Run Hydrology Estimation",
                          "Fits an EMHydroModel to the 2-D resistivity section. "
                          "Required for: K map, Water table, Transmissivity, Aquifer char."),
    "Field Constraints": (None,              "Load borehole / constraint data",
                          "Field constraint plots require borehole or constraint data "
                          "loaded via the Desktop app. Not available in the web interface."),
    "Uncertainty (MC)":  ("run_monte_carlo", "Run Monte Carlo Uncertainty",
                          "Runs ensemble sampling on the hydro model (needs Hydrology done first). "
                          "Required for: K uncertainty section, WT uncertainty, Histograms."),
    "Fusion & Time-Lapse": (None,            "Fusion not yet available",
                            "Fusion & Time-Lapse plots require a fused TDEM+AMT model "
                            "built in the Desktop app."),
}

# ── JS snippet for clientside switcher ───────────────────────────────────────
_SLUGS_JS = '["' + '","'.join(_CAT_SLUG_LIST) + '"]'
_FLEX_STY = (
    '{"display":"flex","flexDirection":"column",'
    '"flex":"1","height":"100%","minHeight":"0","padding":"8px"}'
)


def _get_corrected_sites(session_id: str, corr_store: dict | None, dark: bool):
    from pycsamt.app.desktop.controllers.correction_controller import (
        CorrectionController,
    )
    raw = cache_get(session_id)
    if raw is None:
        return None
    steps = (corr_store or {}).get("steps", [])
    if not steps:
        return raw
    ctrl = CorrectionController()
    ctrl.dark = dark
    ctrl.set_raw_sites(raw)
    for step in steps:
        try:
            ctrl.apply(step["fn_name"], step.get("kwargs", {}),
                       label=step.get("label", step["fn_name"]))
        except Exception:
            pass
    return ctrl.current_sites


def _build_figure(fn_name, cat_slug, data_src, figsize_key, cmap,
                  depth_max, freq_lo, freq_hi,
                  session_id, active_lines_store, store_data, corr_store, theme):
    """Resolve sites, run the plot method, resize and return (fig, plot_label, cat).
    Raises ValueError on any failure."""
    slug = cat_slug or _DEFAULT_CAT_SLUG
    cat  = _CAT_BY_SLUG.get(slug, _INTERP_CATS[0][1])
    dark = (theme or "dark") == "dark"

    if dark:
        apply_web_dark_theme()

    if (data_src or "raw") == "corrected":
        sites = _get_corrected_sites(session_id, corr_store, dark)
    else:
        sites = cache_get(session_id)

    if sites is None:
        raise ValueError("no-data")

    from pycsamt.app.web.utils import filter_sites_by_lines
    _als   = active_lines_store or {}
    active = _als.get("active", _als.get("all"))
    if active is not None:
        records  = (store_data or {}).get("station_records", [])
        filtered = filter_sites_by_lines(sites, records, active)
        if filtered is not None:
            sites = filtered

    _CTRL.set_sites(sites)

    if getattr(_CTRL, fn_name, None) is None:
        raise ValueError(f"unknown-method:{fn_name}")

    extra_kw: dict = {}
    if cmap:           extra_kw["cmap"]      = cmap
    if depth_max:      extra_kw["depth_max"] = float(depth_max)
    if freq_lo is not None: extra_kw["freq_lo"] = float(freq_lo)
    if freq_hi is not None: extra_kw["freq_hi"] = float(freq_hi)

    fig = _CTRL.generate(fn_name, **extra_kw)
    if fig is None:
        raise ValueError("no-figure")

    if (figsize_key or "auto") != "auto":
        fw, fh = _FIGSIZES.get(figsize_key, (11, 5))
        fig.set_size_inches(fw, fh)
    fig.tight_layout()

    plot_label = fn_name
    for lbl, fn, _ in WORKFLOW_CATALOGUE.get(cat, []):
        if fn == fn_name:
            plot_label = lbl
            break

    return fig, plot_label, cat


def register_interp(app) -> None:

    # ── T1. Category pill → store ─────────────────────────────────────────
    @app.callback(
        Output(IDs.INTERP_ACTIVE_CAT, "data"),
        [Input(f"interp-cat-btn-{slug}", "n_clicks") for slug, *_ in _INTERP_CATS],
        prevent_initial_call=True,
    )
    def switch_interp_cat(*_clicks):
        if not any(_clicks):
            raise PreventUpdate
        triggered = ctx.triggered_id or ""
        for slug, *_ in _INTERP_CATS:
            if triggered == f"interp-cat-btn-{slug}":
                return slug
        raise PreventUpdate

    # ── T2. Clientside: panel visibility + button active class ───────────
    clientside_callback(
        f"""
        function(active) {{
            var slugs = {_SLUGS_JS};
            var flex  = {_FLEX_STY};
            var none  = {{display:"none"}};
            var panelStyles = slugs.map(function(s) {{
                return s === active ? flex : none;
            }});
            var btnClasses = slugs.map(function(s) {{
                return "interp-cat-btn" + (s === active ? " active" : "");
            }});
            return panelStyles.concat(btnClasses);
        }}
        """,
        [Output(f"interp-panel-{slug}", "style")      for slug, *_ in _INTERP_CATS] +
        [Output(f"interp-cat-btn-{slug}", "className") for slug, *_ in _INTERP_CATS],
        Input(IDs.INTERP_ACTIVE_CAT, "data"),
    )

    # ── T3. Update context bar when category switches ─────────────────────
    @app.callback(
        Output(IDs.INTERP_INFO_STRIP, "children", allow_duplicate=True),
        Input(IDs.INTERP_ACTIVE_CAT, "data"),
        prevent_initial_call=True,
    )
    def update_interp_ctx_bar(slug):
        if not slug:
            raise PreventUpdate
        cat   = _CAT_BY_SLUG.get(slug, "")
        total = len(WORKFLOW_CATALOGUE.get(cat, []))
        icon  = _CAT_ICON.get(slug, "bi-grid")
        return html.Div([
            html.I(className=f"bi {icon} me-2",
                   style={"fontSize": "12px", "color": "var(--blue)"}),
            html.Span(cat,
                      style={"fontSize": "11px", "fontWeight": "600",
                             "color": "var(--txt)", "marginRight": "8px"}),
            html.Span(f"{total} plot{'s' if total != 1 else ''}",
                      style={"fontSize": "10px", "color": "var(--sub0)"}),
        ], className="d-flex align-items-center")

    # ── 1. Sync plot options when category changes ────────────────────────
    @app.callback(
        Output(IDs.INTERP_PLOT, "options"),
        Output(IDs.INTERP_PLOT, "value"),
        Output(IDs.INTERP_DESC, "children"),
        Input(IDs.INTERP_ACTIVE_CAT, "data"),
    )
    def sync_plots(slug):
        cat   = _CAT_BY_SLUG.get(slug or _DEFAULT_CAT_SLUG,
                                  _INTERP_CATS[0][1])
        plots = WORKFLOW_CATALOGUE.get(cat, [])
        opts  = [{"label": lbl, "value": fn} for lbl, fn, _ in plots]
        desc  = plots[0][2] if plots else ""
        return opts, (opts[0]["value"] if opts else None), desc

    # ── 2. Update description when plot selection changes ─────────────────
    @app.callback(
        Output(IDs.INTERP_DESC, "children", allow_duplicate=True),
        Input(IDs.INTERP_PLOT, "value"),
        State(IDs.INTERP_ACTIVE_CAT, "data"),
        prevent_initial_call=True,
    )
    def update_desc(fn_name, slug):
        if not fn_name:
            raise PreventUpdate
        cat = _CAT_BY_SLUG.get(slug or _DEFAULT_CAT_SLUG,
                                _INTERP_CATS[0][1])
        for _lbl, fn, desc in WORKFLOW_CATALOGUE.get(cat, []):
            if fn == fn_name:
                return desc
        return no_update

    # ── 3. Category → prereq card visibility + button label ──────────────
    @app.callback(
        Output("interp-prereq-card",      "style"),
        Output("interp-prereq-hint",      "children"),
        Output("interp-prereq-btn-label", "children"),
        Input(IDs.INTERP_ACTIVE_CAT, "data"),
    )
    def sync_prereq_card(slug):
        cat  = _CAT_BY_SLUG.get(slug or _DEFAULT_CAT_SLUG,
                                 _INTERP_CATS[0][1])
        info = _CAT_PREREQS.get(cat)
        if info is None:
            return {"display": "none"}, "", "Run"
        fn_name, btn_label, hint_text = info
        style = {"display": "block"}
        if fn_name is None:
            return style, hint_text, "Not available"
        return style, hint_text, btn_label

    @app.callback(
        Output(IDs.INTERP_PREREQ_OUT, "children"),
        Input(IDs.BTN_INTERP_PREREQ,  "n_clicks"),
        State(IDs.INTERP_ACTIVE_CAT,  "data"),
        State(IDs.SESSION_ID,         "data"),
        State(IDs.CORR_STORE,         "data"),
        State(IDs.STORE_THEME,        "data"),
        prevent_initial_call=True,
    )
    def run_prereq(n_clicks, slug, session_id, corr_store, theme):
        if not n_clicks:
            raise PreventUpdate
        cat  = _CAT_BY_SLUG.get(slug or _DEFAULT_CAT_SLUG,
                                 _INTERP_CATS[0][1])
        info = _CAT_PREREQS.get(cat)
        if info is None:
            raise PreventUpdate
        fn_name = info[0]
        if fn_name is None:
            return info[2]
        dark = (theme or "dark") == "dark"
        raw  = cache_get(session_id)
        if raw is None:
            return "No data loaded — load a survey first."
        sites = _get_corrected_sites(session_id, corr_store, dark)
        _CTRL.set_sites(sites)
        try:
            fn = getattr(_CTRL, fn_name, None)
            if fn is None:
                return f"Method '{fn_name}' not found in controller."
            msg = fn()
            return msg or "Done."
        except Exception as exc:
            return f"Error: {exc}"

    # ── 4. Data-source hint text ──────────────────────────────────────────
    @app.callback(
        Output("interp-src-hint", "children"),
        Input(IDs.INTERP_DATA_SRC, "value"),
        State(IDs.CORR_STORE,      "data"),
    )
    def update_src_hint(src, corr_store):
        if src == "corrected":
            n = len((corr_store or {}).get("steps", []))
            if n == 0:
                return ("No corrections applied yet — will use raw data. "
                        "Go to Correction page to build a chain.")
            return f"Using corrected data ({n} step{'s' if n != 1 else ''})."
        return "Using raw loaded data."

    # ── 5. Generate plot ──────────────────────────────────────────────────
    _n_cats = len(_INTERP_CATS)

    @app.callback(
        [Output(f"img-interp-{slug}", "src") for slug, *_ in _INTERP_CATS],
        Output("interp-last-src",         "data"),
        Output(IDs.INTERP_INFO_STRIP,     "children",   allow_duplicate=True),
        Output(IDs.INTERP_SPINNER,        "children"),
        Output(IDs.TOAST_ERROR,           "is_open",    allow_duplicate=True),
        Output(IDs.TOAST_BODY,            "children",   allow_duplicate=True),
        Input(IDs.BTN_INTERP_RUN,         "n_clicks"),
        State(IDs.INTERP_PLOT,            "value"),
        State(IDs.INTERP_ACTIVE_CAT,      "data"),
        State(IDs.INTERP_DATA_SRC,        "value"),
        State(IDs.INTERP_FIGSIZE,         "value"),
        State(IDs.INTERP_CMAP,            "value"),
        State(IDs.INTERP_DEPTH_MAX,       "value"),
        State(IDs.INTERP_FREQ_LO,         "value"),
        State(IDs.INTERP_FREQ_HI,         "value"),
        State(IDs.SESSION_ID,             "data"),
        State(IDs.STORE_ACTIVE_LINES,     "data"),
        State(IDs.STORE_DATA,             "data"),
        State(IDs.CORR_STORE,             "data"),
        State(IDs.STORE_THEME,            "data"),
        prevent_initial_call=True,
    )
    def run_interp(n_clicks, fn_name, cat_slug, data_src, figsize_key, cmap,
                   depth_max, freq_lo, freq_hi,
                   session_id, active_lines_store, store_data, corr_store, theme):
        if not n_clicks or not fn_name:
            raise PreventUpdate

        slug = cat_slug or _DEFAULT_CAT_SLUG
        idx  = _CAT_SLUG_IDX.get(slug, 0)
        imgs = [no_update] * _n_cats
        dark = (theme or "dark") == "dark"

        def _err(msg, exc_str=""):
            imgs[idx] = empty_src(dark=dark)
            return (*imgs, no_update, _info_chip(msg, "danger"),
                    "", True, exc_str or msg)

        if not cache_get(session_id):
            return (*imgs, no_update,
                    _info_chip("No data loaded.", "warning"),
                    "", True, "Load survey data first.")

        try:
            fig, plot_label, cat = _build_figure(
                fn_name, cat_slug, data_src, figsize_key, cmap,
                depth_max, freq_lo, freq_hi,
                session_id, active_lines_store, store_data, corr_store, theme,
            )
            src = fig_to_src(fig)
            plt.close(fig)

            n_corr  = len((corr_store or {}).get("steps", []))
            src_tag = ("corrected" if data_src == "corrected" and n_corr > 0
                       else "raw")
            imgs[idx] = src
            strip = _info_strip(plot_label, cat, src_tag, n_corr,
                                figsize_key or "auto", cmap or "viridis")
            return (*imgs, src, strip, "", False, "")

        except ValueError as exc:
            msg = str(exc)
            if msg == "no-data":
                return (*imgs, no_update,
                        _info_chip("No data loaded.", "warning"),
                        "", True, "Load survey data first.")
            if msg.startswith("unknown-method:"):
                return _err(f"Unknown method: {msg.split(':', 1)[1]}")
            if msg == "no-figure":
                return _err("Method returned no figure.")
            return _err(f"Error: {msg}", msg)
        except Exception as exc:
            plt.close("all")
            return _err(f"Error: {exc}", str(exc))

    # ── 6. Export figure (PNG fast-path; SVG / PDF / EPS re-generate) ────
    @app.callback(
        Output(IDs.INTERP_DOWNLOAD,   "data"),
        Input(IDs.BTN_INTERP_EXPORT,  "n_clicks"),
        State(IDs.INTERP_EXPORT_FMT,  "value"),
        State("interp-last-src",      "data"),
        State(IDs.INTERP_PLOT,        "value"),
        State(IDs.INTERP_ACTIVE_CAT,  "data"),
        State(IDs.INTERP_DATA_SRC,    "value"),
        State(IDs.INTERP_FIGSIZE,     "value"),
        State(IDs.INTERP_CMAP,        "value"),
        State(IDs.INTERP_DEPTH_MAX,   "value"),
        State(IDs.INTERP_FREQ_LO,     "value"),
        State(IDs.INTERP_FREQ_HI,     "value"),
        State(IDs.SESSION_ID,         "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.STORE_DATA,         "data"),
        State(IDs.CORR_STORE,         "data"),
        State(IDs.STORE_THEME,        "data"),
        prevent_initial_call=True,
    )
    def export_plot(n_clicks, fmt, last_src, fn_name, cat_slug,
                    data_src, figsize_key, cmap, depth_max, freq_lo, freq_hi,
                    session_id, active_lines_store, store_data, corr_store, theme):
        if not n_clicks or not fn_name:
            raise PreventUpdate

        fmt   = (fmt or "png").lower()
        fname = f"{fn_name}.{fmt}"

        # PNG fast path — decode the stored data URI, no re-run needed
        if fmt == "png" and last_src:
            try:
                _, b64 = last_src.split(",", 1)
                return dcc.send_bytes(base64.b64decode(b64), fname)
            except Exception:
                pass  # fall through to re-generation

        # Vector / PDF formats — re-generate the figure
        try:
            fig, _label, _cat = _build_figure(
                fn_name, cat_slug, data_src, figsize_key, cmap,
                depth_max, freq_lo, freq_hi,
                session_id, active_lines_store, store_data, corr_store, theme,
            )
            buf = io.BytesIO()
            dpi = 150 if fmt == "png" else None
            fig.savefig(buf, format=fmt, bbox_inches="tight",
                        **({"dpi": dpi} if dpi else {}))
            plt.close(fig)
            buf.seek(0)
            return dcc.send_bytes(buf.read(), fname)
        except Exception:
            plt.close("all")
            raise PreventUpdate


# ── UI helpers ────────────────────────────────────────────────────────────────

def _info_chip(msg: str, color: str = "secondary"):
    colors = {
        "warning":   "#f9e2af",
        "danger":    "#f38ba8",
        "success":   "#a6e3a1",
        "secondary": "#6c7086",
    }
    return html.Span(msg, style={"fontSize": "11px",
                                  "color": colors.get(color, "#cdd6f4")})


def _info_strip(plot_label: str, cat: str, src_tag: str,
                n_corr: int, figsize: str, cmap: str):
    chips = [
        _chip(cat,           "#89b4fa"),
        _chip(plot_label,    "#cdd6f4"),
        _chip(f"src: {src_tag}",
              "#a6e3a1" if src_tag == "corrected" else "#6c7086"),
        _chip(f"fig: {figsize}", "#6c7086"),
        _chip(f"cmap: {cmap}",   "#6c7086"),
    ]
    if src_tag == "corrected" and n_corr > 0:
        chips.append(_chip(f"{n_corr} correction{'s' if n_corr != 1 else ''}",
                           "#cba6f7"))
    return html.Div(chips, className="interp-strip-chips")


def _chip(text: str, color: str):
    return html.Span(
        text,
        style={
            "fontSize":     "10px",
            "padding":      "2px 7px",
            "borderRadius": "10px",
            "background":   "rgba(255,255,255,0.06)",
            "border":       f"1px solid {color}44",
            "color":        color,
            "marginRight":  "4px",
            "display":      "inline-block",
        },
    )
