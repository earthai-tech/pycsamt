# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for the Forward Modelling page (1D, 2D MT, 3D MT)."""
from __future__ import annotations

import matplotlib
import numpy as np

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

_FWD_TAB_IDS  = ["1d", "2d", "3d"]
_FWD_TAB_KEYS = ["fwd-tab-1d", "fwd-tab-2d", "fwd-tab-3d"]
_FWD_TAB_IDX  = {k: i for i, k in enumerate(_FWD_TAB_KEYS)}
_FWD_TAB_ICON = {
    "fwd-tab-1d": "bi-reception-4",
    "fwd-tab-2d": "bi-grid",
    "fwd-tab-3d": "bi-cube",
}
_FWD_TAB_LABEL = {"fwd-tab-1d": "1-D", "fwd-tab-2d": "2-D MT", "fwd-tab-3d": "3-D MT"}

from pycsamt.app.desktop.controllers.forward_controller import (
    ForwardController,
)
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    apply_web_light_theme,
    fig_to_src,
)

_CTRL = ForwardController()

_SHOW = {"display": "block"}
_HIDE = {"display": "none"}


# ── Helpers ───────────────────────────────────────────────────────────────────

def _build_model(table_data, halfspace_rho):
    """Build LayeredModel from finite-layer table rows + halfspace ρ."""
    from pycsamt.forward import LayeredModel

    if not table_data:
        raise ValueError("Add at least one finite layer before running.")
    halfspace_rho = float(halfspace_rho or 1000.0)

    rho_finite = np.array([float(r["resistivity"]) for r in table_data])
    thk_finite = np.array([float(r["thickness"])   for r in table_data])

    if np.any(thk_finite <= 0):
        raise ValueError("All layer thicknesses must be > 0 m.")
    if np.any(rho_finite <= 0) or halfspace_rho <= 0:
        raise ValueError("All resistivities must be > 0 Ω·m.")

    rho = np.append(rho_finite, halfspace_rho)
    thk = thk_finite
    return LayeredModel(resistivity=rho, thickness=thk)


def _plot_model_1d(ax, model, dark: bool):
    rho   = model.resistivity
    thk   = model.thickness
    tops  = np.concatenate([[0.0], np.cumsum(thk)])

    depths_plot, rho_plot = [], []
    for i in range(len(rho) - 1):
        depths_plot += [tops[i], tops[i + 1]]
        rho_plot    += [rho[i],  rho[i]]
    bottom = tops[-1] * 1.2 if tops[-1] > 0 else 500.0
    depths_plot += [tops[-1], bottom]
    rho_plot    += [rho[-1],  rho[-1]]

    colour = "#89b4fa" if dark else "#1e66f5"
    ax.plot(rho_plot, depths_plot, color=colour, lw=2)
    ax.fill_betweenx(depths_plot, rho_plot, alpha=0.12, color=colour)
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Resistivity (Ω·m)")
    ax.set_ylabel("Depth (m)")
    ax.set_title("1-D Model")
    ax.set_xlim(left=0.1)


def _plot_mt_response(axes, resp, dark: bool):
    c_rho   = "#a6e3a1" if dark else "#40a02b"
    c_phase = "#f38ba8" if dark else "#d20f39"
    T = 1.0 / resp.freqs

    axes[0].loglog(T, resp.rho_a, color=c_rho, lw=2)
    axes[0].set_xlabel("Period (s)")
    axes[0].set_ylabel("ρₐ (Ω·m)")
    axes[0].set_title("Apparent Resistivity")

    phase = np.where(resp.phase < 0, resp.phase + 180, resp.phase)
    axes[1].semilogx(T, phase, color=c_phase, lw=2)
    axes[1].set_xlabel("Period (s)")
    axes[1].set_ylabel("Phase (°)")
    axes[1].set_title("Impedance Phase")
    axes[1].set_ylim(0, 90)


def _plot_tem_response(ax, resp, dark: bool):
    colour = "#f9e2af" if dark else "#df8e1d"
    vals   = np.abs(resp.dBz_dt)
    mask   = np.isfinite(vals) & (vals > 0)
    if mask.sum() < 2:
        ax.text(0.5, 0.5, "No valid TEM response",
                ha="center", va="center", transform=ax.transAxes)
        return
    ax.loglog(resp.times[mask], vals[mask], color=colour, lw=2)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("|dBz/dt| (V/m²)")
    ax.set_title("TEM Step-off Response")


# ── 1D runner ─────────────────────────────────────────────────────────────────

def _run_1d(table_data, halfspace_rho, method, f_min, f_max, n_pts,
            offset, loop_r, dark):
    model  = _build_model(table_data, halfspace_rho)
    method = (method or "MT1D").upper()

    if method == "TEM1D":
        from pycsamt.forward.em1d import TEM1DForward
        times = np.logspace(f_min, f_max, n_pts)
        resp  = TEM1DForward(times, loop_radius=float(loop_r or 50)).run(model)

        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        _plot_model_1d(axes[0], model, dark)
        _plot_tem_response(axes[1], resp, dark)

    elif method == "CSAMT1D":
        from pycsamt.forward.em1d import CSAMT1DForward
        freqs = np.logspace(f_min, f_max, n_pts)
        r     = float(offset or 5000)
        resp  = CSAMT1DForward(freqs, source_offset=r).run(model)

        fig, axes = plt.subplots(1, 3, figsize=(14, 5))
        _plot_model_1d(axes[0], model, dark)
        _plot_mt_response(axes[1:], resp, dark)

    else:
        from pycsamt.forward.em1d import MT1DForward
        freqs = np.logspace(f_min, f_max, n_pts)
        resp  = MT1DForward(freqs).run(model)

        fig, axes = plt.subplots(1, 3, figsize=(14, 5))
        _plot_model_1d(axes[0], model, dark)
        _plot_mt_response(axes[1:], resp, dark)

    fig.tight_layout()
    src = fig_to_src(fig)
    plt.close(fig)
    n_layers = model.n_layers
    gate_label = "time gates" if method == "TEM1D" else "frequencies"
    msg = (f"{method} done — {n_layers} layers "
           f"({n_layers - 1} finite + halfspace), {n_pts} {gate_label}")
    return src, msg


# ── 2D runner ─────────────────────────────────────────────────────────────────

def _run_2d(m2d_type, bg2, nx2, nz2, x2max, z2max, n2st,
            a2rho, ax2lo, ax2hi, az2lo, az2hi, n2lay, seed2, plot2,
            table_data, halfspace_rho, f_min, f_max, n_pts, dark):
    from pycsamt.forward import Grid2D, MT2DForward
    from pycsamt.forward.plot import (
        plot_model_2d,
        plot_pseudosection_2d,
        plot_response_profiles,
    )

    m2d_type = m2d_type or "halfspace"
    bg2      = float(bg2 or 100.0)
    nx2      = int(nx2 or 25)
    nz2      = int(nz2 or 20)
    x2max    = float(x2max or 8000)
    z2max    = float(z2max or 5000)
    n2st     = int(n2st or 10)
    plot2    = plot2 or "model"

    if m2d_type == "anomaly":
        bounds = (float(ax2lo or 2000), float(ax2hi or 4000),
                  float(az2lo or 200),  float(az2hi or 800))
        grid = Grid2D.with_anomaly(
            bg_rho=bg2, anomaly_rho=float(a2rho or 5),
            anomaly_bounds=bounds,
            nx=nx2, nz=nz2, x_max=x2max, z_max=z2max, n_stations=n2st,
        )

    elif m2d_type == "from_layers":
        model_1d = _build_model(table_data, halfspace_rho)
        grid = Grid2D.from_1d_layers(
            model_1d, nx=nx2, x_max=x2max, n_stations=n2st,
        )

    elif m2d_type == "random":
        grid = Grid2D.random(
            nx=nx2, nz=nz2, x_max=x2max, z_max=z2max,
            n_layers=int(n2lay or 4), n_stations=n2st,
            lateral_variation=True, seed=int(seed2) if seed2 is not None else 42,
        )

    else:  # halfspace
        grid = Grid2D.halfspace(
            rho=bg2, nx=nx2, nz=nz2, x_max=x2max, z_max=z2max, n_stations=n2st,
        )

    freqs  = np.logspace(f_min, f_max, n_pts)
    resp2d = MT2DForward(freqs, grid, verbose=False).run()

    # Choose plot
    if plot2 == "model":
        fig, ax = plt.subplots(figsize=(11, 4))
        plot_model_2d(grid, ax=ax, show_stations=True)

    elif plot2 == "psd_te_rho":
        fig, ax = plt.subplots(figsize=(11, 5))
        plot_pseudosection_2d(resp2d, mode="te", quantity="rho_a", ax=ax)

    elif plot2 == "psd_tm_rho":
        fig, ax = plt.subplots(figsize=(11, 5))
        plot_pseudosection_2d(resp2d, mode="tm", quantity="rho_a", ax=ax)

    elif plot2 == "psd_te_phi":
        fig, ax = plt.subplots(figsize=(11, 5))
        plot_pseudosection_2d(resp2d, mode="te", quantity="phase", ax=ax)

    elif plot2 == "psd_tm_phi":
        fig, ax = plt.subplots(figsize=(11, 5))
        plot_pseudosection_2d(resp2d, mode="tm", quantity="phase", ax=ax)

    elif plot2 == "profiles_te":
        fig, ax = plt.subplots(figsize=(10, 4))
        plot_response_profiles(resp2d, mode="te", quantity="rho_a", ax=ax)

    elif plot2 == "profiles_tm":
        fig, ax = plt.subplots(figsize=(10, 4))
        plot_response_profiles(resp2d, mode="tm", quantity="rho_a", ax=ax)

    else:
        fig, ax = plt.subplots(figsize=(11, 4))
        plot_model_2d(grid, ax=ax)

    fig.tight_layout()
    src = fig_to_src(fig)
    plt.close(fig)

    n_st = len(resp2d.stations_x)
    msg = (f"MT2D done — {m2d_type} model, "
           f"{n_pts} freqs, {n_st} stations, "
           f"grid {grid.nx}×{grid.nz}")
    return src, msg, resp2d


# ── 3D runner ─────────────────────────────────────────────────────────────────

def _run_3d(m3d_type, bg3, nx3, ny3, nz3, x3max, y3max, z3max, nx3st, ny3st,
            a3rho, ax3lo, ax3hi, ay3lo, ay3hi, az3lo, az3hi, n3lay, seed3,
            plot3, comp3, freq3idx, f_min, f_max, n_pts, dark):
    from pycsamt.forward import Grid3D, MT3DForward
    from pycsamt.forward.plot import (
        plot_model_3d,
        plot_response_map_3d,
        plot_response_section_3d,
        plot_tensor_components_3d,
    )

    m3d_type = m3d_type or "halfspace"
    bg3      = float(bg3 or 100.0)
    nx3      = int(nx3 or 15)
    ny3      = int(ny3 or 15)
    nz3      = int(nz3 or 12)
    x3max    = float(x3max or 8000)
    y3max    = float(y3max or 8000)
    z3max    = float(z3max or 5000)
    nx3st    = int(nx3st or 4)
    ny3st    = int(ny3st or 4)
    plot3    = plot3 or "model"
    comp3    = comp3 or "xy"
    freq3idx = int(freq3idx or 0)

    if m3d_type == "block":
        bounds = (
            float(ax3lo or 2000), float(ax3hi or 6000),
            float(ay3lo or 2000), float(ay3hi or 6000),
            float(az3lo or 300),  float(az3hi or 1500),
        )
        grid = Grid3D.block_anomaly(
            bg_rho=bg3, anomaly_rho=float(a3rho or 5),
            bounds=bounds,
            nx=nx3, ny=ny3, nz=nz3,
            x_max=x3max, y_max=y3max, z_max=z3max,
            nx_stations=nx3st, ny_stations=ny3st,
        )

    elif m3d_type == "random":
        grid = Grid3D.random_layered(
            nx=nx3, ny=ny3, nz=nz3,
            x_max=x3max, y_max=y3max, z_max=z3max,
            n_layers=int(n3lay or 4),
            nx_stations=nx3st, ny_stations=ny3st,
            lateral_variation=True,
            seed=int(seed3) if seed3 is not None else 42,
        )

    else:  # halfspace
        grid = Grid3D.halfspace(
            rho=bg3, nx=nx3, ny=ny3, nz=nz3,
            x_max=x3max, y_max=y3max, z_max=z3max,
            nx_stations=nx3st, ny_stations=ny3st,
        )

    freqs  = np.logspace(f_min, f_max, n_pts)
    resp3d = MT3DForward(freqs, grid, method="quasi3d", verbose=False).run()

    # Clamp freq index
    n_freqs    = len(resp3d.freqs)
    freq3idx   = min(freq3idx, n_freqs - 1)

    if plot3 == "model":
        axes = plot_model_3d(grid, show_stations=True)
        fig  = np.asarray(axes).flat[0].figure

    elif plot3 == "map":
        fig, ax = plt.subplots(figsize=(7, 6))
        plot_response_map_3d(resp3d, freq_idx=freq3idx, component=comp3,
                             quantity="rho_a", ax=ax)

    elif plot3 == "section":
        fig, ax = plt.subplots(figsize=(11, 5))
        plot_response_section_3d(resp3d, component=comp3, quantity="rho_a", ax=ax)

    elif plot3 == "tensor":
        axes = plot_tensor_components_3d(resp3d, freq_idx=freq3idx, quantity="rho_a")
        fig  = np.asarray(axes).flat[0].figure

    else:
        axes = plot_model_3d(grid, show_stations=True)
        fig  = np.asarray(axes).flat[0].figure

    try:
        fig.tight_layout()
    except Exception:
        pass
    src = fig_to_src(fig)
    plt.close(fig)

    n_st = len(resp3d.stations_xy)
    msg = (f"MT3D quasi-3D done — {m3d_type} model, "
           f"{n_pts} freqs, {n_st} stations, "
           f"grid {grid.nx}×{grid.ny}×{grid.nz}")
    return src, msg, resp3d


# ── Callbacks ─────────────────────────────────────────────────────────────────

def register_forward(app) -> None:

    # 1a. Dim-bar buttons → FWD_ACTIVE_TAB store
    @app.callback(
        Output(IDs.FWD_ACTIVE_TAB, "data"),
        [Input(f"fwd-dim-btn-{t}", "n_clicks") for t in _FWD_TAB_IDS],
        prevent_initial_call=True,
    )
    def switch_fwd_dim(*_):
        tid = ctx.triggered_id
        if not tid:
            return no_update
        return "fwd-tab-" + tid.replace("fwd-dim-btn-", "")

    # 1b. View-tab buttons → FWD_ACTIVE_TAB store (allow_duplicate)
    @app.callback(
        Output(IDs.FWD_ACTIVE_TAB, "data", allow_duplicate=True),
        [Input(f"fwd-view-tab-{t}", "n_clicks") for t in _FWD_TAB_IDS],
        prevent_initial_call=True,
    )
    def switch_fwd_view_tab(*_):
        tid = ctx.triggered_id
        if not tid:
            return no_update
        return "fwd-tab-" + tid.replace("fwd-view-tab-", "")

    # 1c. Store → all panel/button visibility (clientside, instant)
    clientside_callback(
        """
        function(active_tab) {
            if (!active_tab)
                return Array(12).fill(window.dash_clientside.no_update);
            var keys  = ["1d", "2d", "3d"];
            var flex  = {display:"flex", flexDirection:"column",
                         height:"100%", minHeight:"0", flex:"1"};
            var blk   = {display:"block"};
            var none  = {display:"none"};
            var ctrl  = keys.map(function(k){
                return active_tab === "fwd-tab-" + k ? blk : none;
            });
            var dCls  = keys.map(function(k){
                return "fwd-dim-btn" + (active_tab === "fwd-tab-" + k ? " active" : "");
            });
            var view  = keys.map(function(k){
                return active_tab === "fwd-tab-" + k ? flex : none;
            });
            var vCls  = keys.map(function(k){
                return "prof-tab-btn" + (active_tab === "fwd-tab-" + k ? " active" : "");
            });
            return ctrl.concat(dCls).concat(view).concat(vCls);
        }
        """,
        # sidebar ctrl panels
        Output("fwd-ctrl-1d",     "style"),
        Output("fwd-ctrl-2d",     "style"),
        Output("fwd-ctrl-3d",     "style"),
        # dim-bar btn classes
        Output("fwd-dim-btn-1d",  "className"),
        Output("fwd-dim-btn-2d",  "className"),
        Output("fwd-dim-btn-3d",  "className"),
        # view panels
        Output(IDs.FWD_1D_PANEL,  "style"),
        Output(IDs.FWD_2D_PANEL,  "style"),
        Output(IDs.FWD_3D_PANEL,  "style"),
        # view tab classes
        Output("fwd-view-tab-1d", "className"),
        Output("fwd-view-tab-2d", "className"),
        Output("fwd-view-tab-3d", "className"),
        Input(IDs.FWD_ACTIVE_TAB, "data"),
        prevent_initial_call=False,
    )

    # 1d. Context info bar update
    @app.callback(
        Output(IDs.FWD_CTX_BAR, "children"),
        Input(IDs.FWD_ACTIVE_TAB, "data"),
        Input(IDs.FWD_METHOD,     "value"),
    )
    def update_fwd_ctx_bar(active_tab, method):
        active_tab = active_tab or "fwd-tab-1d"
        icon  = _FWD_TAB_ICON.get(active_tab, "bi-reception-4")
        label = _FWD_TAB_LABEL.get(active_tab, "1-D")
        parts = [
            html.Span(
                [html.I(className=f"bi {icon} me-1"), label],
                className="adv-info-group",
            ),
        ]
        if active_tab == "fwd-tab-1d" and method:
            parts += [
                html.Span("·", className="adv-info-sep"),
                html.Span(method, className="adv-info-plot"),
            ]
        parts += [
            html.Span("·", className="adv-info-sep"),
            html.Span(
                "Configure parameters then click Run Forward",
                style={"color": "var(--sub0)", "fontSize": "10.5px"},
            ),
        ]
        return parts

    # 2. Show/hide 1D method-specific params + relabel freq axis
    @app.callback(
        Output(IDs.FWD_CSAMT_CARD,  "style"),
        Output(IDs.FWD_TEM_CARD,    "style"),
        Output(IDs.FWD_FREQ_LABEL,  "children"),
        Output(IDs.FWD_FREQ_MIN,    "value"),
        Output(IDs.FWD_FREQ_MAX,    "value"),
        Output(IDs.FWD_DIM,         "data"),
        Input(IDs.FWD_METHOD,       "value"),
    )
    def sync_method_ui(method):
        method = method or "MT1D"
        csamt_style = _SHOW if method == "CSAMT1D" else _HIDE
        tem_style   = _SHOW if method == "TEM1D"   else _HIDE
        if method == "TEM1D":
            label    = "Time range [log₁₀ s]"
            f_min, f_max = -6, -2
        else:
            label    = "Frequency range [log₁₀ Hz]"
            f_min, f_max = -3, 4
        return csamt_style, tem_style, label, f_min, f_max, method

    # 3. Show/hide 2D anomaly / random cards
    @app.callback(
        Output(IDs.FWD2_ANOMALY_CARD, "style"),
        Output(IDs.FWD2_RANDOM_CARD,  "style"),
        Input(IDs.FWD2_MODEL_TYPE,    "value"),
    )
    def sync_2d_model(m_type):
        return (
            _SHOW if m_type == "anomaly" else _HIDE,
            _SHOW if m_type == "random"  else _HIDE,
        )

    # 4. Show/hide 3D block / random cards
    @app.callback(
        Output(IDs.FWD3_BLOCK_CARD,  "style"),
        Output(IDs.FWD3_RANDOM_CARD, "style"),
        Input(IDs.FWD3_MODEL_TYPE,   "value"),
    )
    def sync_3d_model(m_type):
        return (
            _SHOW if m_type == "block"  else _HIDE,
            _SHOW if m_type == "random" else _HIDE,
        )

    # 5. Load geological preset (1D only)
    @app.callback(
        Output(IDs.FWD_LAYER_TABLE,   "data"),
        Output(IDs.FWD_HALFSPACE_RHO, "value"),
        Output(IDs.FWD_FEEDBACK,      "children"),
        Input(IDs.BTN_FWD_PRESET,     "n_clicks"),
        State(IDs.FWD_PRESET,         "value"),
        prevent_initial_call=True,
    )
    def load_preset(n_clicks, preset_name):
        if not n_clicks or not preset_name:
            raise PreventUpdate
        try:
            record = ForwardController.build_preset_1d(preset_name)
            rho    = record["resistivity"]
            thk    = record["thickness"]
            rows = [
                {"layer": i + 1, "resistivity": float(rho[i]), "thickness": float(thk[i])}
                for i in range(len(thk))
            ]
            return rows, float(rho[-1]), f"Preset '{preset_name}' loaded — {len(rho)} layers"
        except Exception as exc:
            return no_update, no_update, f"Error: {exc}"

    # 6. Add a finite layer
    @app.callback(
        Output(IDs.FWD_LAYER_TABLE, "data", allow_duplicate=True),
        Input(IDs.BTN_FWD_ADD_LAYER, "n_clicks"),
        State(IDs.FWD_LAYER_TABLE,   "data"),
        prevent_initial_call=True,
    )
    def add_layer(n_clicks, current_data):
        if not n_clicks:
            raise PreventUpdate
        data     = list(current_data or [])
        last_rho = float(data[-1]["resistivity"]) if data else 100.0
        data.append({"layer": len(data) + 1, "resistivity": last_rho, "thickness": 500})
        return data

    # 7. Renumber layers on row delete
    @app.callback(
        Output(IDs.FWD_LAYER_TABLE, "data", allow_duplicate=True),
        Input(IDs.FWD_LAYER_TABLE,  "data"),
        prevent_initial_call=True,
    )
    def renumber_layers(data):
        if not data:
            raise PreventUpdate
        for i, row in enumerate(data):
            row["layer"] = i + 1
        return data

    # 8. Run forward (dispatches to 1D / 2D / 3D, writes only to active tab image)
    @app.callback(
        Output(IDs.IMG_FWD_1D,        "src"),
        Output(IDs.IMG_FWD_2D,        "src"),
        Output(IDs.IMG_FWD_3D,        "src"),
        Output(IDs.FWD_SPINNER,       "children"),
        Output(IDs.FWD_FEEDBACK,      "children",  allow_duplicate=True),
        Output(IDs.TOAST_ERROR,       "is_open",   allow_duplicate=True),
        Output(IDs.TOAST_BODY,        "children",  allow_duplicate=True),
        Input(IDs.BTN_FWD_RUN,        "n_clicks"),
        # Dimension
        State(IDs.FWD_ACTIVE_TAB,     "data"),
        # 1D States
        State(IDs.FWD_LAYER_TABLE,    "data"),
        State(IDs.FWD_HALFSPACE_RHO,  "value"),
        State(IDs.FWD_METHOD,         "value"),
        # Shared freq
        State(IDs.FWD_FREQ_MIN,       "value"),
        State(IDs.FWD_FREQ_MAX,       "value"),
        State(IDs.FWD_N_FREQ,         "value"),
        # 1D method-specific
        State(IDs.FWD_OFFSET,         "value"),
        State(IDs.FWD_LOOP_R,         "value"),
        # 2D States
        State(IDs.FWD2_MODEL_TYPE,    "value"),
        State(IDs.FWD2_BG_RHO,        "value"),
        State(IDs.FWD2_NX,            "value"),
        State(IDs.FWD2_NZ,            "value"),
        State(IDs.FWD2_X_MAX,         "value"),
        State(IDs.FWD2_Z_MAX,         "value"),
        State(IDs.FWD2_N_STATIONS,    "value"),
        State(IDs.FWD2_ANOMALY_RHO,   "value"),
        State(IDs.FWD2_AX_LO,         "value"),
        State(IDs.FWD2_AX_HI,         "value"),
        State(IDs.FWD2_AZ_LO,         "value"),
        State(IDs.FWD2_AZ_HI,         "value"),
        State(IDs.FWD2_N_LAYERS,      "value"),
        State(IDs.FWD2_SEED,          "value"),
        State(IDs.FWD2_PLOT_TYPE,     "value"),
        # 3D States
        State(IDs.FWD3_MODEL_TYPE,    "value"),
        State(IDs.FWD3_BG_RHO,        "value"),
        State(IDs.FWD3_NX,            "value"),
        State(IDs.FWD3_NY,            "value"),
        State(IDs.FWD3_NZ,            "value"),
        State(IDs.FWD3_X_MAX,         "value"),
        State(IDs.FWD3_Y_MAX,         "value"),
        State(IDs.FWD3_Z_MAX,         "value"),
        State(IDs.FWD3_NX_ST,         "value"),
        State(IDs.FWD3_NY_ST,         "value"),
        State(IDs.FWD3_ANOMALY_RHO,   "value"),
        State(IDs.FWD3_AX_LO,         "value"),
        State(IDs.FWD3_AX_HI,         "value"),
        State(IDs.FWD3_AY_LO,         "value"),
        State(IDs.FWD3_AY_HI,         "value"),
        State(IDs.FWD3_AZ_LO,         "value"),
        State(IDs.FWD3_AZ_HI,         "value"),
        State(IDs.FWD3_N_LAYERS,      "value"),
        State(IDs.FWD3_SEED,          "value"),
        State(IDs.FWD3_PLOT_TYPE,     "value"),
        State(IDs.FWD3_COMPONENT,     "value"),
        State(IDs.FWD3_FREQ_IDX,      "value"),
        # Theme + session
        State(IDs.STORE_THEME,        "data"),
        State(IDs.SESSION_ID,         "data"),
        prevent_initial_call=True,
    )
    def run_forward(
        n_clicks, dimension,
        # 1D
        table_data, halfspace_rho, method,
        freq_min, freq_max, n_freq, offset, loop_r,
        # 2D
        m2d_type, bg2, nx2, nz2, x2max, z2max, n2st,
        a2rho, ax2lo, ax2hi, az2lo, az2hi, n2lay, seed2, plot2,
        # 3D
        m3d_type, bg3, nx3, ny3, nz3, x3max, y3max, z3max, nx3st, ny3st,
        a3rho, ax3lo, ax3hi, ay3lo, ay3hi, az3lo, az3hi, n3lay, seed3,
        plot3, comp3, freq3idx,
        # theme + session
        theme, session_id,
    ):
        if not n_clicks:
            raise PreventUpdate

        dark = (theme or "dark") == "dark"
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        dimension = dimension or "fwd-tab-1d"
        n_pts = int(n_freq or 20)
        f_min = float(freq_min if freq_min is not None else -3)
        f_max = float(freq_max if freq_max is not None else 4)

        # Outputs: [src_1d, src_2d, src_3d, spinner, feedback, toast_open, toast_body]
        imgs = [no_update, no_update, no_update]

        try:
            from pycsamt.app.web.cache import cache_set_fwd
            if dimension == "fwd-tab-2d":
                src, msg, resp = _run_2d(
                    m2d_type, bg2, nx2, nz2, x2max, z2max, n2st,
                    a2rho, ax2lo, ax2hi, az2lo, az2hi, n2lay, seed2, plot2,
                    table_data, halfspace_rho, f_min, f_max, n_pts, dark,
                )
                cache_set_fwd(session_id, "2d", resp)
                imgs[1] = src

            elif dimension == "fwd-tab-3d":
                src, msg, resp = _run_3d(
                    m3d_type, bg3, nx3, ny3, nz3, x3max, y3max, z3max, nx3st, ny3st,
                    a3rho, ax3lo, ax3hi, ay3lo, ay3hi, az3lo, az3hi, n3lay, seed3,
                    plot3, comp3, freq3idx, f_min, f_max, n_pts, dark,
                )
                cache_set_fwd(session_id, "3d", resp)
                imgs[2] = src

            else:  # 1D
                src, msg = _run_1d(
                    table_data, halfspace_rho, method,
                    f_min, f_max, n_pts, offset, loop_r, dark,
                )
                imgs[0] = src

            return (*imgs, "", msg, False, "")

        except Exception as exc:
            plt.close("all")
            return (*imgs, "", f"Error: {exc}", True, str(exc))

    # 9. Save model to library (1D only)
    @app.callback(
        Output(IDs.FWD_FEEDBACK,      "children",  allow_duplicate=True),
        Input(IDs.BTN_FWD_SAVE,       "n_clicks"),
        State(IDs.FWD_MODEL_NAME,     "value"),
        State(IDs.FWD_LAYER_TABLE,    "data"),
        State(IDs.FWD_HALFSPACE_RHO,  "value"),
        State(IDs.FWD_METHOD,         "value"),
        prevent_initial_call=True,
    )
    def save_model(n_clicks, name, table_data, halfspace_rho, method):
        if not n_clicks:
            raise PreventUpdate
        if not name:
            return "Enter a model name first."
        try:
            model  = _build_model(table_data, halfspace_rho)
            record = ForwardController.model_to_record(
                name, method or "MT1D", model.resistivity, model.thickness,
            )
            _CTRL.save_model(name, record)
            return f"Saved '{name}' ({model.n_layers} layers, {method})."
        except Exception as exc:
            return f"Save error: {exc}"
