# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Side-by-side true vs predicted 1-D resistivity profile panels.

:func:`plot_compare` is the primary entry point.  Each panel shows:

* Predicted model (solid, purple, ``EM_COLORS['pred']``)
* True model (dashed, green, ``EM_COLORS['true']``)
* Shaded difference band (alpha = 0.15)
* Per-site RMSE annotation

Usage
-----
>>> from pycsamt.ai.plot.compare import plot_compare
>>> fig = plot_compare(true_models, pred_models, n_cols=5)
>>> fig.savefig("comparison.png", dpi=300)
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

__all__ = ["plot_compare", "plot_profile_pair"]


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────


def plot_compare(
    true_models,
    pred_models,
    *,
    n_cols: int = 5,
    max_sites: int = 20,
    depth_max: float | None = None,
    log_scale: bool = True,
    show_rmse: bool = True,
    site_labels: Sequence[str] | None = None,
    figsize_per_panel: tuple = (2.2, 4.0),
    title: str | None = None,
    style: bool = True,
):
    """
    Multi-panel comparison of true and predicted 1-D resistivity profiles.

    Parameters
    ----------
    true_models : list of LayeredModel or ndarray (n_sites, n_params)
        Ground-truth models.
    pred_models : list of LayeredModel or ndarray (n_sites, n_params)
        Network-predicted models (same order as *true_models*).
    n_cols : int
        Number of columns in the panel grid.
    max_sites : int
        Maximum number of sites to plot (first *max_sites* are used).
    depth_max : float or None
        Maximum depth axis value [m].
    log_scale : bool
        Log₁₀ x-axis for resistivity.
    show_rmse : bool
        Annotate each panel with per-site log₁₀(ρ) RMSE.
    site_labels : list of str or None
        Panel titles (default: ``'Site 0'``, ``'Site 1'``, …).
    figsize_per_panel : (w, h)
        Size of each individual panel in inches.
    title : str or None
        Figure suptitle.
    style : bool
        Apply :class:`~pycsamt.ai.plot._style.EMStyle`.

    Returns
    -------
    fig : Figure
    """
    import matplotlib.pyplot as plt

    from ._style import EM_COLORS, EMStyle

    n_sites = min(len(true_models), len(pred_models), max_sites)
    n_rows = int(np.ceil(n_sites / n_cols))
    fw = figsize_per_panel[0] * n_cols
    fh = figsize_per_panel[1] * n_rows

    ctx = EMStyle() if style else _NullContext()
    with ctx:
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(fw, fh),
            sharey=False,
        )
        axes = np.asarray(axes).ravel()

        for i in range(n_sites):
            ax = axes[i]
            label = site_labels[i] if site_labels else f"Site {i}"
            tm = _to_model(true_models[i])
            pm = _to_model(pred_models[i])
            _draw_pair(
                ax,
                tm,
                pm,
                depth_max=depth_max,
                log_scale=log_scale,
                show_rmse=show_rmse,
                label=label,
                colors=EM_COLORS,
            )

        # Hide empty panels
        for j in range(n_sites, len(axes)):
            axes[j].set_visible(False)

        if title:
            fig.suptitle(title, fontsize=13, y=1.01)

        fig.tight_layout()
    return fig


def plot_profile_pair(
    true_model,
    pred_model,
    *,
    ax=None,
    depth_max: float | None = None,
    log_scale: bool = True,
    show_rmse: bool = True,
    legend: bool = True,
    style: bool = True,
):
    """
    Plot a single true/predicted resistivity–depth pair on *ax*.

    Parameters
    ----------
    true_model, pred_model : LayeredModel or ndarray (n_params,)
    ax : Axes or None
    depth_max : float or None
    log_scale : bool
    show_rmse : bool
    legend : bool
    style : bool

    Returns
    -------
    ax : Axes
    """
    import matplotlib.pyplot as plt

    from ._style import EM_COLORS, EMStyle

    if ax is None:
        ctx = EMStyle() if style else _NullContext()
        with ctx:
            _, ax = plt.subplots(figsize=(3.5, 5))

    tm = _to_model(true_model)
    pm = _to_model(pred_model)
    _draw_pair(
        ax,
        tm,
        pm,
        depth_max=depth_max,
        log_scale=log_scale,
        show_rmse=show_rmse,
        label="",
        colors=EM_COLORS,
        legend=legend,
    )
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────


def _draw_pair(
    ax,
    true_m,
    pred_m,
    *,
    depth_max,
    log_scale,
    show_rmse,
    label,
    colors,
    legend=False,
):
    """Draw one profile pair onto *ax*."""
    # True model
    _plot_profile(
        ax, true_m, color=colors["true"], linestyle="--", lw=1.4, label="True"
    )
    # Predicted model
    _plot_profile(
        ax,
        pred_m,
        color=colors["pred"],
        linestyle="-",
        lw=1.8,
        label="Predicted",
    )

    # Shade the difference
    _shade_diff(ax, true_m, pred_m)

    ax.invert_yaxis()
    if log_scale:
        ax.set_xscale("log")
    dmax = depth_max or _max_depth(true_m, pred_m)
    ax.set_ylim(dmax, 0.0)
    ax.set_xlabel(r"ρ (Ω·m)", fontsize=9)
    ax.set_ylabel("Depth (m)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_title(label, fontsize=9, pad=3)

    if show_rmse:
        err = _rmse_log_rho(true_m, pred_m)
        if np.isfinite(err):
            ax.text(
                0.97,
                0.97,
                f"RMSE={err:.3f}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=7.5,
                color=colors["error"],
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1),
            )
    if legend:
        ax.legend(fontsize=8, loc="lower left")


def _plot_profile(ax, model, **kw):
    """Staircase resistivity–depth plot."""
    rho, depth = _staircase(model)
    ax.plot(rho, depth, **kw)


def _shade_diff(ax, true_m, pred_m):
    """Shade area between true and predicted profiles."""
    rho_t, depth = _staircase(true_m)
    rho_p, _ = _staircase(pred_m)
    n = min(len(rho_t), len(rho_p))
    ax.fill_betweenx(depth[:n], rho_t[:n], rho_p[:n], alpha=0.12, color="#762a83")


def _staircase(model):
    """Convert LayeredModel or ndarray to (rho_steps, depth_steps)."""
    rho, depths = _model_arrays(model)
    # Build staircase: each layer contributes top and bottom value
    n = len(rho)
    rho_step = np.repeat(rho, 2)
    d_step = np.empty(2 * n)
    d_step[0::2] = depths
    d_step[1::2] = np.concatenate([depths[1:], [depths[-1] * 1.5]])
    return rho_step, d_step


def _model_arrays(model):
    """Return (rho_array, depth_array) from LayeredModel or vector."""
    try:
        from pycsamt.forward.synthetic import LayeredModel

        if isinstance(model, LayeredModel):
            return model.resistivity, model.depth
    except ImportError:
        pass
    # Assume flat vector: first half = log10(rho), second half = thick
    v = np.asarray(model, dtype=float).ravel()
    n = (len(v) + 1) // 2
    rho = 10.0 ** v[:n]
    thick = np.maximum(v[n:], 1.0)
    depths = np.concatenate([[0.0], np.cumsum(thick)])
    return rho, depths


def _max_depth(true_m, pred_m):
    _, dt = _model_arrays(true_m)
    _, dp = _model_arrays(pred_m)
    return float(max(dt[-1], dp[-1])) * 1.1


def _rmse_log_rho(true_m, pred_m):
    rho_t, _ = _model_arrays(true_m)
    rho_p, _ = _model_arrays(pred_m)
    n = min(len(rho_t), len(rho_p))
    diff = np.log10(np.maximum(rho_t[:n], 1e-6)) - np.log10(np.maximum(rho_p[:n], 1e-6))
    return float(np.sqrt(np.mean(diff**2)))


def _to_model(m):
    """Pass-through — already a model or array."""
    return m


class _NullContext:
    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass
