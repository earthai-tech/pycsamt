# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callbacks for PINN / Hybrid inversion parameter
offcanvas."""

from __future__ import annotations

from dash import Input, Output, State, no_update
from dash.exceptions import PreventUpdate

from .._ids import IDs

_PINN_KWS = {
    "pinn", "physics-informed",
    "physics informed", "no training data",
    "gradient descent",
}
_HYBRID_KWS = {
    "hybrid", "two-stage", "two stage",
    "ai + physics", "warm start", "warmstart",
}


def _is_pinn_or_hybrid(text: str) -> bool:
    t = text.lower()
    return (
        any(kw in t for kw in _PINN_KWS)
        or any(kw in t for kw in _HYBRID_KWS)
    )


def register_inv_params(app) -> None:
    """Register all PINN/Hybrid parameter callbacks."""

    # ── Dim radio → show/hide 2-D / 3-D panels ───────
    @app.callback(
        Output(IDs.INV_PANEL_2D, "style"),
        Output(IDs.INV_PANEL_3D, "style"),
        Input(IDs.INV_DIM, "value"),
        prevent_initial_call=False,
    )
    def toggle_dim_panels(dim):
        show = {"display": "block"}
        hide = {"display": "none"}
        d2 = show if dim in (2, 3) else hide
        d3 = show if dim == 3 else hide
        return d2, d3

    # ── Mode radio → show/hide hybrid checkpoint ─────
    @app.callback(
        Output(IDs.INV_PANEL_HYBRID, "style"),
        Input(IDs.INV_MODE, "value"),
        prevent_initial_call=False,
    )
    def toggle_hybrid_panel(mode):
        if mode == "hybrid":
            return {"display": "block"}
        return {"display": "none"}

    # ── Confirm → save to store, close canvas ────────
    @app.callback(
        Output(IDs.STORE_INV_CONFIG, "data"),
        Output(IDs.INV_STATUS, "children"),
        Output(IDs.CANVAS_INV_PARAMS, "is_open",
               allow_duplicate=True),
        Input(IDs.BTN_INV_CONFIRM, "n_clicks"),
        State(IDs.INV_MODE, "value"),
        State(IDs.INV_DIM, "value"),
        State(IDs.INV_SOLVER, "value"),
        State(IDs.INV_N_LAYERS, "value"),
        State(IDs.INV_DEPTH_MAX, "value"),
        State(IDs.INV_EPOCHS, "value"),
        State(IDs.INV_LR, "value"),
        State(IDs.INV_SMOOTH_W, "value"),
        State(IDs.INV_LAT_W, "value"),
        State(IDs.INV_GRAPH_W, "value"),
        State(IDs.INV_RADIUS, "value"),
        State(IDs.INV_CHECKPOINT, "value"),
        prevent_initial_call=True,
    )
    def confirm_inv_params(
        n_clicks,
        mode, dim, solver,
        n_layers, depth_max,
        epochs, lr,
        smooth_w, lat_w,
        graph_w, radius,
        checkpoint,
    ):
        if not n_clicks:
            raise PreventUpdate

        # --- basic validation -------------------------
        errors = []
        if not n_layers or int(n_layers) < 2:
            errors.append(
                "n_layers must be >= 2."
            )
        if not depth_max or float(depth_max) <= 0:
            errors.append(
                "depth_max must be positive."
            )
        if not epochs or int(epochs) < 1:
            errors.append(
                "epochs must be >= 1."
            )
        if not lr or float(lr) <= 0:
            errors.append(
                "Learning rate must be > 0."
            )
        if errors:
            from dash import html
            msg = html.Span(
                "; ".join(errors),
                style={"color": "var(--red)"},
            )
            return no_update, msg, no_update

        cfg = {
            "mode": mode or "pinn",
            "dim": int(dim or 1),
            "solver": solver or "mt1d",
            "n_layers": int(n_layers),
            "depth_max": float(depth_max),
            "epochs": int(epochs),
            "lr": float(lr),
            "smoothness_weight": float(
                smooth_w or 0.01
            ),
            "lateral_weight": float(lat_w or 0.005),
            "graph_weight": float(graph_w or 0.005),
            "radius": float(radius or 5000.0),
            "checkpoint": checkpoint or "",
        }
        from dash import html
        ok_msg = html.Span(
            "Parameters saved.",
            style={"color": "var(--green)"},
        )
        return cfg, ok_msg, False

    # ── Load store values into form on canvas open ────
    @app.callback(
        Output(IDs.INV_MODE, "value"),
        Output(IDs.INV_DIM, "value"),
        Output(IDs.INV_SOLVER, "value"),
        Output(IDs.INV_N_LAYERS, "value"),
        Output(IDs.INV_DEPTH_MAX, "value"),
        Output(IDs.INV_EPOCHS, "value"),
        Output(IDs.INV_LR, "value"),
        Output(IDs.INV_SMOOTH_W, "value"),
        Output(IDs.INV_LAT_W, "value"),
        Output(IDs.INV_GRAPH_W, "value"),
        Output(IDs.INV_RADIUS, "value"),
        Output(IDs.INV_CHECKPOINT, "value"),
        Input(IDs.CANVAS_INV_PARAMS, "is_open"),
        State(IDs.STORE_INV_CONFIG, "data"),
        prevent_initial_call=True,
    )
    def load_inv_form(is_open, cfg):
        if not is_open or not cfg:
            raise PreventUpdate
        return (
            cfg.get("mode", "pinn"),
            cfg.get("dim", 1),
            cfg.get("solver", "mt1d"),
            cfg.get("n_layers", 10),
            cfg.get("depth_max", 2000.0),
            cfg.get("epochs", 500),
            cfg.get("lr", 0.01),
            cfg.get("smoothness_weight", 0.01),
            cfg.get("lateral_weight", 0.005),
            cfg.get("graph_weight", 0.005),
            cfg.get("radius", 5000.0),
            cfg.get("checkpoint", ""),
        )
