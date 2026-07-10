"""
Bridge helpers for PINN / Hybrid inverters in
the web app.

Public helpers
--------------
session_to_obs_1d
session_to_obs_2d
decode_npz_checkpoint
build_pinn_inv
build_hybrid_inv
plot_pinn_section
plot_pinn_convergence
plot_pinn_data_fit
pinn_stats_div
"""

from __future__ import annotations

import base64
import tempfile
from typing import Any

import numpy as np

# ---

_DARK_LINE   = "#89b4fa"   # Catppuccin Sapphire
_DARK_OBS    = "#f38ba8"   # Red
_DARK_PRED   = "#a6e3a1"   # Green
_DARK_LOSS   = "#f38ba8"

_LIGHT_LINE  = "#1e66f5"
_LIGHT_OBS   = "#d20f39"
_LIGHT_PRED  = "#40a02b"
_LIGHT_LOSS  = "#d20f39"


# ---
# Session -> observation list
# ---

def session_to_obs_1d(
    session_id: str,
    stations: list[str] | None = None,
    *,
    comp: str = "xy",
) -> list:
    r"""
    Return a ``List[SiteObs1D]`` from the session
    cache, optionally filtered to *stations*.

    Parameters
    ----------
    session_id : str
        Active browser session key.
    stations : list of str, optional
        Station names to keep. ``None`` keeps all.
    comp : str, default ``'xy'``
        Impedance component for apparent resistivity.

    Returns
    -------
    list of SiteObs1D
    """
    from pycsamt.ai.inversion._sites_bridge import (
        sites_to_obs_1d,
    )
    from pycsamt.app.web.cache import cache_get

    sites = cache_get(session_id)
    if not sites:
        raise ValueError(
            "No data loaded for this session."
        )
    obs = sites_to_obs_1d(
        sites, comp=comp, verbose=0
    )
    if stations:
        keep = set(stations)
        obs = [o for o in obs if o.name in keep]
        if not obs:
            raise ValueError(
                "None of the selected stations were"
                " found in the session data."
            )
    return obs


def session_to_obs_2d(
    session_id: str,
    stations: list[str] | None = None,
    *,
    comp_te: str = "xy",
    comp_tm: str = "yx",
) -> list:
    r"""
    Return a ``List[SiteObs2D]`` from the session
    cache, optionally filtered to *stations*.

    Parameters
    ----------
    session_id : str
        Active browser session key.
    stations : list of str, optional
        Station names to keep. ``None`` keeps all.
    comp_te : str, default ``'xy'``
        TE-mode impedance component.
    comp_tm : str, default ``'yx'``
        TM-mode impedance component.

    Returns
    -------
    list of SiteObs2D
    """
    from pycsamt.ai.inversion._sites_bridge import (
        sites_to_obs_2d,
    )
    from pycsamt.app.web.cache import cache_get

    sites = cache_get(session_id)
    if not sites:
        raise ValueError(
            "No data loaded for this session."
        )
    obs = sites_to_obs_2d(
        sites,
        comp_te=comp_te,
        comp_tm=comp_tm,
        verbose=0,
    )
    if stations:
        keep = set(stations)
        obs = [o for o in obs if o.name in keep]
        if not obs:
            raise ValueError(
                "None of the selected stations were"
                " found in the session data."
            )
    return obs


# ---
# Checkpoint decode
# ---

def decode_npz_checkpoint(
    content_b64: str,
) -> str:
    r"""
    Decode a base64 ``.npz`` Dash upload and write
    it to a temporary file.

    Parameters
    ----------
    content_b64 : str
        ``data:application/...;base64,<payload>``
        string as returned by ``dcc.Upload``.

    Returns
    -------
    str
        Path to the temp ``.npz`` file (caller is
        responsible for cleanup when done).
    """
    _, payload = content_b64.split(",", 1)
    raw = base64.b64decode(payload)
    fd, path = tempfile.mkstemp(suffix=".npz")
    with open(fd, "wb") as fh:
        fh.write(raw)
    return path


# ---
# PINN inverter factory
# ---

def build_pinn_inv(
    obs: list,
    dim: str,
    *,
    n_layers: int = 5,
    depth_max: float = 2000.0,
    n_freqs: int = 32,
    solver: str = "mt1d",
    mode: str = "te",
    smoothness_weight: float = 0.01,
    lateral_weight: float = 0.005,
    graph_weight: float = 0.005,
    radius: float = 5000.0,
    station_spacing: float = 500.0,
    epochs: int = 300,
    lr: float = 0.01,
    comp: str = "xy",
    comp_te: str = "xy",
    comp_tm: str = "yx",
    device: str | None = None,
) -> Any:
    r"""
    Instantiate a ``PINNInverter{1D,2D,3D}`` from a
    pre-built observation list.

    The short-circuit guards in
    ``_sites_bridge.sites_to_obs_1d/2d`` ensure the
    constructors accept the obs list directly without
    re-loading EDI files.

    Parameters
    ----------
    obs : list
        ``List[SiteObs1D]`` for ``dim='1d'``;
        ``List[SiteObs2D]`` for ``dim='2d'`` or
        ``dim='3d'``.
    dim : {'1d', '2d', '3d'}
        Inversion dimensionality.
    n_layers : int, default 5
    depth_max : float, default 2000.0
    n_freqs : int, default 32
    solver : str, default ``'mt1d'``
        Forward solver (1-D only).
    mode : str, default ``'te'``
        EM polarisation (2-D / 3-D).
    smoothness_weight : float, default 0.01
    lateral_weight : float, default 0.005
    graph_weight : float, default 0.005
    radius : float, default 5000.0
        Neighbour radius in metres (3-D).
    station_spacing : float, default 500.0
        Fallback uniform grid step (3-D).
    lr : float, default 0.01
    comp : str, default ``'xy'``
        Impedance component (1-D).
    comp_te : str, default ``'xy'``
    comp_tm : str, default ``'yx'``
    device : str or None
        ``'cpu'``, ``'cuda'``, ``'mps'``, or
        ``None`` for auto-detect.

    Returns
    -------
    PINNInverter1D | PINNInverter2D | PINNInverter3D
    """
    if device == "auto":
        device = None

    if dim == "1d":
        from pycsamt.ai.inversion import (
            PINNInverter1D,
        )
        return PINNInverter1D(
            obs,
            solver=solver,
            n_layers=n_layers,
            depth_max=depth_max,
            smoothness_weight=smoothness_weight,
            lr=lr,
            comp=comp,
            device=device,
        )

    if dim == "2d":
        from pycsamt.ai.inversion import (
            PINNInverter2D,
        )
        return PINNInverter2D(
            obs,
            n_layers=n_layers,
            depth_max=depth_max,
            n_freqs=n_freqs,
            mode=mode,
            smoothness_weight=smoothness_weight,
            lateral_weight=lateral_weight,
            epochs=epochs,
            lr=lr,
            comp_te=comp_te,
            comp_tm=comp_tm,
            device=device,
        )

    # dim == "3d"
    from pycsamt.ai.inversion import PINNInverter3D

    n_sta = len(obs)
    step  = float(station_spacing)
    xy = np.column_stack([
        np.arange(n_sta, dtype=float) * step,
        np.zeros(n_sta),
    ])
    return PINNInverter3D(
        obs,
        n_layers=n_layers,
        depth_max=depth_max,
        n_freqs=n_freqs,
        mode=mode,
        smoothness_weight=smoothness_weight,
        graph_weight=graph_weight,
        radius=radius,
        station_coords=xy,
        epochs=epochs,
        lr=lr,
        comp_te=comp_te,
        comp_tm=comp_tm,
        device=device,
    )


# ---
# Hybrid inverter factory
# ---

def build_hybrid_inv(
    obs: list,
    dim: str,
    checkpoint: Any,
    *,
    n_layers: int | None = None,
    n_freqs: int = 32,
    solver: str = "mt1d",
    mode: str = "te",
    smoothness_weight: float = 0.005,
    lateral_weight: float = 0.003,
    graph_weight: float = 0.003,
    radius: float = 5000.0,
    station_spacing: float = 500.0,
    max_iter: int = 150,
    lr: float = 0.005,
    comp: str = "xy",
    comp_te: str = "xy",
    comp_tm: str = "yx",
    device: str | None = None,
) -> Any:
    r"""
    Instantiate a ``HybridInverter{1D,2D,3D}`` from a
    pre-built observation list.

    Parameters
    ----------
    obs : list
        ``List[SiteObs1D]`` for ``dim='1d'``;
        ``List[SiteObs2D]`` otherwise.
    dim : {'1d', '2d', '3d'}
    checkpoint : str or fitted inverter
        Path to ``.npz`` checkpoint, or a fitted
        ``EMInverter*`` instance.
    n_layers : int or None
        Layer count; ``None`` uses the AI model's
        ``n_depth`` attribute.
    n_freqs : int, default 32
    solver : str, default ``'mt1d'``
    mode : str, default ``'te'``
    smoothness_weight : float, default 0.005
    lateral_weight : float, default 0.003
    graph_weight : float, default 0.003
    radius : float, default 5000.0
    station_spacing : float, default 500.0
    max_iter : int, default 150
        Stage-2 Adam epochs.
    lr : float, default 0.005
    comp : str, default ``'xy'``
    comp_te : str, default ``'xy'``
    comp_tm : str, default ``'yx'``
    device : str or None

    Returns
    -------
    HybridInverter1D | HybridInverter2D |
    HybridInverter3D
    """
    if device == "auto":
        device = None

    if dim == "1d":
        from pycsamt.ai.inversion import (
            HybridInverter1D,
        )
        return HybridInverter1D(
            obs,
            checkpoint,
            solver=solver,
            max_iter=max_iter,
            smoothness_weight=smoothness_weight,
            lr=lr,
            comp=comp,
            n_freqs=n_freqs,
            device=device,
        )

    if dim == "2d":
        from pycsamt.ai.inversion import (
            HybridInverter2D,
        )
        return HybridInverter2D(
            obs,
            checkpoint,
            n_layers=n_layers,
            n_freqs=n_freqs,
            mode=mode,
            smoothness_weight=smoothness_weight,
            lateral_weight=lateral_weight,
            epochs=max_iter,
            lr=lr,
            comp_te=comp_te,
            comp_tm=comp_tm,
            device=device,
        )

    # dim == "3d"
    from pycsamt.ai.inversion import HybridInverter3D

    n_sta = len(obs)
    step  = float(station_spacing)
    xy = np.column_stack([
        np.arange(n_sta, dtype=float) * step,
        np.zeros(n_sta),
    ])
    return HybridInverter3D(
        obs,
        checkpoint,
        n_layers=n_layers,
        n_freqs=n_freqs,
        mode=mode,
        smoothness_weight=smoothness_weight,
        graph_weight=graph_weight,
        radius=radius,
        station_coords=xy,
        epochs=max_iter,
        lr=lr,
        comp_te=comp_te,
        comp_tm=comp_tm,
        device=device,
    )


# ---
# Plot: resistivity section / volume
# ---

def plot_pinn_section(
    inv: Any,
    dim: str,
    *,
    dark: bool = True,
    cmap: str = "viridis_r",
) -> str:
    r"""
    Render the PINN resistivity result as a PNG and
    return a base64 ``data:image/png`` URI.

    For ``dim='1d'`` a panel of step plots is drawn
    (one per station).  For ``dim='2d'`` and
    ``dim='3d'`` the section / pseudo-volume is shown
    as a colour image.

    Parameters
    ----------
    inv : fitted inverter
    dim : {'1d', '2d', '3d'}
    dark : bool, default True
    cmap : str, default ``'viridis_r'``

    Returns
    -------
    str
        ``data:image/png;base64,...``
    """
    import matplotlib.pyplot as plt

    from pycsamt.app.web.utils import (
        apply_web_dark_theme,
        apply_web_light_theme,
        fig_to_src,
    )

    if dark:
        apply_web_dark_theme()
    else:
        apply_web_light_theme()

    if dim == "1d":
        return _plot_1d_models(
            inv, dark=dark, fig_to_src=fig_to_src
        )

    # 2-D / 3-D colour section
    if dim == "3d":
        mat = inv.resistivity_volume(as_log10=True)
    else:
        mat = inv.resistivity_section(as_log10=True)

    stations = list(getattr(inv, "stations", []))
    n_sta = mat.shape[1]

    fig, ax = plt.subplots(
        figsize=(max(8, n_sta * 0.45 + 2), 5)
    )
    im = ax.imshow(
        mat,
        aspect="auto",
        cmap=cmap,
        origin="upper",
        interpolation="bilinear",
    )
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(
        "log10(resistivity / Ohm.m)", fontsize=9
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("Depth layer")
    ax.set_title(
        f"PINN {dim.upper()} Resistivity Section"
    )
    if stations:
        step = max(1, len(stations) // 15)
        ax.set_xticks(
            range(0, len(stations), step)
        )
        ax.set_xticklabels(
            stations[::step],
            rotation=45,
            ha="right",
            fontsize=7,
        )
    fig.tight_layout()
    return fig_to_src(fig)


def _plot_1d_models(
    inv: Any,
    *,
    dark: bool,
    fig_to_src,
) -> str:
    """Step plots for all 1-D results."""
    import matplotlib.pyplot as plt

    c = _DARK_LINE if dark else _LIGHT_LINE
    stations = list(getattr(inv, "stations", []))
    models = inv.predict()
    n = max(len(models), 1)
    ncols = min(n, 5)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(ncols * 3, nrows * 4),
        squeeze=False,
    )
    flat = axes.ravel()

    for i, m in enumerate(models):
        ax = flat[i]
        rho  = np.asarray(m.resistivity)
        thk  = np.asarray(m.thickness)
        tops = np.concatenate([[0.0], np.cumsum(thk)])
        d, r = [], []
        for j in range(len(rho) - 1):
            d += [tops[j], tops[j + 1]]
            r += [rho[j], rho[j]]
        bot = max(tops[-1] * 1.3, 500.0)
        d += [tops[-1], bot]
        r += [rho[-1], rho[-1]]
        ax.plot(r, d, color=c, lw=1.5)
        ax.set_xscale("log")
        ax.invert_yaxis()
        ax.set_xlabel("rho (Ohm.m)", fontsize=7)
        ax.set_ylabel("Depth (m)", fontsize=7)
        name = (
            stations[i]
            if i < len(stations)
            else f"S{i}"
        )
        ax.set_title(name, fontsize=8)

    for i in range(len(models), len(flat)):
        flat[i].set_visible(False)

    fig.tight_layout()
    return fig_to_src(fig)


# ---
# Plot: Adam loss convergence
# ---

def plot_pinn_convergence(
    inv: Any,
    *,
    dark: bool = True,
    label: str = "PINN",
) -> str:
    r"""
    Plot the Adam loss history and return a base64 PNG.

    Tries ``convergence_curve()`` (2-D / 3-D) first,
    then ``loss_curves()`` (1-D, averaged across
    stations), then falls back to ``_history``.

    Parameters
    ----------
    inv : fitted inverter
    dark : bool, default True
    label : str, default ``'PINN'``

    Returns
    -------
    str
        ``data:image/png;base64,...``
    """
    import matplotlib.pyplot as plt

    from pycsamt.app.web.utils import (
        apply_web_dark_theme,
        apply_web_light_theme,
        fig_to_src,
    )

    if dark:
        apply_web_dark_theme()
    else:
        apply_web_light_theme()

    c = _DARK_LOSS if dark else _LIGHT_LOSS
    history = _extract_loss_history(inv)

    if not history:
        fig, ax = plt.subplots(figsize=(7, 3))
        ax.text(
            0.5, 0.5,
            "No loss history available.",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        fig.tight_layout()
        return fig_to_src(fig)

    epochs = np.arange(1, len(history) + 1)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogy(epochs, history, color=c, lw=1.8)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.set_title(f"{label} Adam Convergence")
    fig.tight_layout()
    return fig_to_src(fig)


def _extract_loss_history(inv: Any) -> list[float]:
    """Pull the mean loss-per-epoch from any inverter."""
    # PINNInverter2D / PINNInverter3D
    fn = getattr(inv, "convergence_curve", None)
    if callable(fn):
        try:
            df = fn()
            if hasattr(df, "columns"):
                if "loss" in df.columns:
                    return df["loss"].tolist()
        except Exception:
            pass

    # PINNInverter1D / HybridInverter1D
    fn = getattr(inv, "loss_curves", None)
    if callable(fn):
        try:
            df = fn()
            if "loss" in df.columns:
                return (
                    df.groupby("epoch")["loss"]
                    .mean()
                    .tolist()
                )
        except Exception:
            pass

    # Hybrid 2-D / 3-D stage-2 history stored
    # in a shared list attribute
    hist = getattr(inv, "_history", None)
    if hist:
        return list(hist)

    return []


# ---
# Plot: observed vs predicted data fit
# ---

def plot_pinn_data_fit(
    inv: Any,
    station: str,
    *,
    dark: bool = True,
    df=None,
) -> str:
    r"""
    Plot observed vs predicted apparent resistivity
    and phase for one station.

    Parameters
    ----------
    inv : fitted inverter or None
        Must support ``.residuals()`` returning a
        DataFrame with columns ``station, freq,
        rho_obs, rho_pred, phase_obs, phase_pred``.
        Ignored when *df* is provided.
    station : str
        Station name to plot.
    dark : bool, default True
    df : DataFrame or None
        Pre-computed residuals DataFrame.  When
        provided, *inv* is not called.

    Returns
    -------
    str
        ``data:image/png;base64,...``
    """
    import matplotlib.pyplot as plt

    from pycsamt.app.web.utils import (
        apply_web_dark_theme,
        apply_web_light_theme,
        fig_to_src,
    )

    if dark:
        apply_web_dark_theme()
    else:
        apply_web_light_theme()

    c_obs  = _DARK_OBS  if dark else _LIGHT_OBS
    c_pred = _DARK_PRED if dark else _LIGHT_PRED

    if df is None:
        try:
            df = inv.residuals()
        except Exception as exc:
            fig, ax = plt.subplots(figsize=(7, 3))
            ax.text(
                0.5, 0.5,
                f"Could not compute"
                f" residuals:\n{exc}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            fig.tight_layout()
            return fig_to_src(fig)

    df = df[df["station"] == station]

    if df.empty:
        fig, ax = plt.subplots(figsize=(7, 3))
        ax.text(
            0.5, 0.5,
            f"No data for station {station!r}.",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        fig.tight_layout()
        return fig_to_src(fig)

    freq     = df["freq"].values
    period   = 1.0 / freq
    rho_obs  = df["rho_obs"].values
    rho_pred = df["rho_pred"].values
    ph_obs   = df["phase_obs"].values
    ph_pred  = df["phase_pred"].values

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10, 4)
    )

    ax1.loglog(
        period, rho_obs,
        "o", color=c_obs, ms=4, label="Observed",
    )
    ax1.loglog(
        period, rho_pred,
        "-", color=c_pred, lw=2, label="Predicted",
    )
    ax1.set_xlabel("Period (s)")
    ax1.set_ylabel("App. resistivity (Ohm.m)")
    ax1.set_title(f"{station}  rho_a")
    ax1.legend(fontsize=8)

    ax2.semilogx(
        period, ph_obs,
        "o", color=c_obs, ms=4, label="Observed",
    )
    ax2.semilogx(
        period, ph_pred,
        "-", color=c_pred, lw=2, label="Predicted",
    )
    ax2.set_xlabel("Period (s)")
    ax2.set_ylabel("Phase (deg)")
    ax2.set_title(f"{station}  phase")
    ax2.legend(fontsize=8)

    fig.suptitle(
        f"Data Fit  --  {station}", fontsize=11
    )
    fig.tight_layout()
    return fig_to_src(fig)


# ---
# Statistics HTML div
# ---

def pinn_stats_div(
    inv: Any,
    dim: str,
    *,
    elapsed_s: float = 0.0,
    label: str = "PINN",
) -> Any:
    r"""
    Return an ``html.Div`` with a concise run summary.

    Parameters
    ----------
    inv : fitted inverter
    dim : {'1d', '2d', '3d'}
    elapsed_s : float, default 0.0
        Wall-clock seconds the inversion took.
    label : str, default ``'PINN'``

    Returns
    -------
    dash.html.Div
    """
    from dash import html

    n_sites   = getattr(inv, "n_sites", "?")
    n_layers  = getattr(inv, "n_layers", "?")
    depth_max = getattr(inv, "depth_max", "?")
    epochs    = getattr(
        inv, "epochs",
        getattr(inv, "max_iter", "?")
    )
    lr        = getattr(inv, "lr", "?")

    rows: list = [
        ("Track",          f"{label} {dim.upper()}"),
        ("Stations",       str(n_sites)),
        ("Layers",         str(n_layers)),
        ("Depth max (m)",  str(depth_max)),
        ("Epochs",         str(epochs)),
        ("Learning rate",  str(lr)),
        ("Elapsed (s)",    f"{elapsed_s:.1f}"),
    ]

    loss_history = _extract_loss_history(inv)
    if loss_history:
        rows.append(
            ("Final loss",
             f"{loss_history[-1]:.4e}")
        )

    _td_key = {
        "color":         "#a6adc8",
        "paddingRight":  "16px",
        "fontSize":      "11px",
        "whiteSpace":    "nowrap",
    }
    _td_val = {
        "color":    "#cdd6f4",
        "fontSize": "11px",
    }

    trs = [
        html.Tr([
            html.Td(k, style=_td_key),
            html.Td(v, style=_td_val),
        ])
        for k, v in rows
    ]
    return html.Div(
        html.Table(
            trs,
            style={
                "width":           "100%",
                "borderCollapse":  "collapse",
            },
        ),
        className="inv-stats-pinn",
    )
