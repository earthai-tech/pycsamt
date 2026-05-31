# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Diagnostic plots for AI/ML model evaluation.

All functions follow the :class:`~pycsamt.ai.plot._style.EMStyle`
publication conventions and accept an optional ``ax`` parameter for
embedding in composite figures.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from ._style import EMStyle, EM_COLORS, EM_FIGSIZE, add_colorbar

__all__ = [
    "plot_confusion_matrix",
    "plot_residuals",
    "plot_layer_errors",
    "plot_uncertainty_bands",
    "plot_feature_importance",
]

_CLASS_NAMES_DEFAULT = ["1-D", "2-D", "3-D"]


# ─────────────────────────────────────────────────────────────────────────────
# plot_confusion_matrix
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_confusion_matrix(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    class_names: Optional[List[str]] = None,
    normalise: bool = True,
    cmap: str = "Blues",
    title: str = "Confusion Matrix",
    figsize: Optional[Tuple[float, float]] = None,
    ax: Optional[Axes] = None,
    style: bool = True,
) -> Figure:
    """
    Plot a confusion matrix for a classification model.

    Parameters
    ----------
    y_true : int ndarray (n_samples,)
        Ground-truth class labels.
    y_pred : int ndarray (n_samples,)
        Predicted class labels.
    class_names : list of str or None
        Display labels for each class.  Defaults to ``['1-D','2-D','3-D']``
        for three-class problems.
    normalise : bool
        Show row-normalised (recall) fractions.
    cmap : str
        Matplotlib colormap.
    title, figsize, ax, style : see :func:`plot_section`.

    Returns
    -------
    fig : Figure
    """
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    classes = np.unique(np.concatenate([y_true, y_pred]))
    n = len(classes)

    if class_names is None:
        class_names = (
            _CLASS_NAMES_DEFAULT[:n] if n <= 3
            else [str(c) for c in classes]
        )

    # Build confusion matrix
    cm = np.zeros((n, n), dtype=int)
    for t, p in zip(y_true, y_pred):
        ti = np.where(classes == t)[0]
        pi = np.where(classes == p)[0]
        if len(ti) and len(pi):
            cm[ti[0], pi[0]] += 1

    if normalise:
        row_sums = cm.sum(axis=1, keepdims=True)
        cm_plot = cm.astype(float) / (row_sums + 1e-12)
        fmt = ".2f"
        vmax = 1.0
    else:
        cm_plot = cm.astype(float)
        fmt = "d"
        vmax = cm.max() or 1

    if figsize is None:
        s = max(3.0, n * 1.2)
        figsize = (s, s)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    im = ax.imshow(cm_plot, cmap=cmap, vmin=0, vmax=vmax, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(class_names, fontsize=9)
    ax.set_yticklabels(class_names, fontsize=9)
    ax.set_xlabel("Predicted", fontsize=9)
    ax.set_ylabel("True", fontsize=9)
    ax.set_title(title, fontsize=10)

    thresh = vmax / 2.0
    for i in range(n):
        for j in range(n):
            val = cm_plot[i, j]
            color = "white" if val > thresh else EM_COLORS["text"]
            text = f"{val:{fmt}}" if fmt != "d" else str(int(cm[i, j]))
            ax.text(j, i, text, ha="center", va="center",
                    color=color, fontsize=9)

    add_colorbar(im, ax, label="Recall" if normalise else "Count", pad=0.03)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# plot_residuals
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_residuals(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    param_names: Optional[List[str]] = None,
    n_cols: int = 4,
    figsize_per_panel: Tuple[float, float] = (2.5, 2.5),
    style: bool = True,
) -> Figure:
    """
    Scatter plots of predicted vs. true for each model parameter.

    Each panel shows the 1:1 line, coloured scatter, and per-parameter
    R² annotation.

    Parameters
    ----------
    y_true : ndarray (n_samples, n_params)
    y_pred : ndarray (n_samples, n_params)
    param_names : list of str or None
    n_cols : int
    figsize_per_panel : (width, height)
    style : bool

    Returns
    -------
    fig : Figure
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    if y_true.ndim == 1:
        y_true = y_true[:, np.newaxis]
        y_pred = y_pred[:, np.newaxis]

    n_params = y_true.shape[1]
    if param_names is None:
        param_names = [f"param {i}" for i in range(n_params)]

    n_rows = int(np.ceil(n_params / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_per_panel[0] * n_cols, figsize_per_panel[1] * n_rows),
    )
    axes = np.array(axes).reshape(-1)

    color = EM_COLORS["primary"]

    for pi in range(n_params):
        ax = axes[pi]
        yt = y_true[:, pi]
        yp = y_pred[:, pi]
        mask = np.isfinite(yt) & np.isfinite(yp)
        yt_m, yp_m = yt[mask], yp[mask]

        lo = min(yt_m.min(), yp_m.min()) if len(yt_m) else 0
        hi = max(yt_m.max(), yp_m.max()) if len(yt_m) else 1
        ax.plot([lo, hi], [lo, hi], "--", color=EM_COLORS["error"],
                lw=0.8, zorder=1)
        ax.scatter(yt_m, yp_m, s=6, alpha=0.5, color=color,
                   linewidths=0, zorder=2)

        if len(yt_m) > 1:
            ss_res = np.sum((yt_m - yp_m) ** 2)
            ss_tot = np.sum((yt_m - yt_m.mean()) ** 2)
            r2 = 1.0 - ss_res / (ss_tot + 1e-12)
            ax.text(0.05, 0.93, f"R²={r2:.3f}", transform=ax.transAxes,
                    fontsize=7, va="top")

        ax.set_title(param_names[pi], fontsize=8)
        ax.set_xlabel("True", fontsize=7)
        ax.set_ylabel("Predicted", fontsize=7)
        ax.tick_params(labelsize=6)

    for pi in range(n_params, len(axes)):
        axes[pi].set_visible(False)

    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# plot_layer_errors
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_layer_errors(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    n_layers: int,
    *,
    log_rho: bool = True,
    ax: Optional[Axes] = None,
    figsize: Optional[Tuple[float, float]] = None,
    style: bool = True,
) -> Figure:
    """
    Per-layer mean absolute error bar chart.

    Parameters
    ----------
    y_true : ndarray (n_samples, 2*n_layers-1)
    y_pred : ndarray (n_samples, 2*n_layers-1)
    n_layers : int
    log_rho : bool
        Label ρ columns as log₁₀(ρ).
    ax : Axes or None
    figsize : (width, height) or None
    style : bool

    Returns
    -------
    fig : Figure
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    mae = np.nanmean(np.abs(y_true - y_pred), axis=0)
    n_rho = n_layers
    n_thick = n_layers - 1
    rho_mae = mae[:n_rho]
    thick_mae = mae[n_rho: n_rho + n_thick]

    rho_lbl = [
        (r"$\log_{10}\rho_{%d}$" if log_rho else r"$\rho_{%d}$") % (i + 1)
        for i in range(n_rho)
    ]
    thick_lbl = [r"$h_{%d}$" % (i + 1) for i in range(n_thick)]

    labels = rho_lbl + thick_lbl
    values = np.concatenate([rho_mae, thick_mae])
    colors = (
        [EM_COLORS["primary"]] * n_rho
        + [EM_COLORS["secondary"]] * n_thick
    )

    if figsize is None:
        figsize = (max(5.0, len(labels) * 0.5), 3.5)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    x = np.arange(len(labels))
    bars = ax.bar(x, values, color=colors, width=0.6, edgecolor="white", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, rotation=45, ha="right")
    ax.set_ylabel("Mean absolute error", fontsize=9)
    ax.set_title("Per-parameter MAE", fontsize=10)

    # Legend
    from matplotlib.patches import Patch
    ax.legend(
        handles=[
            Patch(color=EM_COLORS["primary"], label="Resistivity"),
            Patch(color=EM_COLORS["secondary"], label="Thickness"),
        ],
        fontsize=8,
        frameon=False,
    )

    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# plot_uncertainty_bands
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_uncertainty_bands(
    x: np.ndarray,
    y_pred: np.ndarray,
    y_upper: np.ndarray,
    y_lower: np.ndarray,
    y_true: Optional[np.ndarray] = None,
    *,
    ax: Optional[Axes] = None,
    xlabel: str = "",
    ylabel: str = "",
    title: str = "Prediction with Uncertainty",
    figsize: Optional[Tuple[float, float]] = None,
    style: bool = True,
) -> Figure:
    """
    1-D prediction curve with uncertainty bands.

    Suitable for showing per-site model parameter predictions (e.g.
    resistivity vs. depth) with ±1σ or 10/90 percentile bands.

    Parameters
    ----------
    x : ndarray (n_points,)
        X-axis values (e.g. depth or frequency).
    y_pred : ndarray (n_points,)
        Central prediction.
    y_upper, y_lower : ndarray (n_points,)
        Upper and lower uncertainty bounds.
    y_true : ndarray (n_points,) or None
        Ground-truth values to overlay.
    ax : Axes or None
    xlabel, ylabel, title : str
    figsize : (width, height) or None
    style : bool

    Returns
    -------
    fig : Figure
    """
    x = np.asarray(x, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    y_upper = np.asarray(y_upper, dtype=float)
    y_lower = np.asarray(y_lower, dtype=float)

    if figsize is None:
        figsize = EM_FIGSIZE["single"]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    ax.fill_betweenx(
        x, y_lower, y_upper,
        color=EM_COLORS["primary"], alpha=0.20, linewidth=0,
        label="Uncertainty band",
    )
    ax.plot(y_pred, x, color=EM_COLORS["primary"], lw=1.5,
            label="Predicted", zorder=3)

    if y_true is not None:
        y_true = np.asarray(y_true, dtype=float)
        ax.plot(y_true, x, color=EM_COLORS["secondary"], lw=1.5,
                ls="--", label="True", zorder=4)

    ax.invert_yaxis()
    ax.set_xlabel(xlabel or r"$\log_{10}(\rho)$ (Ω·m)", fontsize=9)
    ax.set_ylabel(ylabel or "Depth (m)", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# plot_feature_importance
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_feature_importance(
    importances: np.ndarray,
    feature_names: Optional[List[str]] = None,
    *,
    top_n: int = 20,
    horizontal: bool = True,
    ax: Optional[Axes] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: str = "Feature Importance",
    style: bool = True,
) -> Figure:
    """
    Horizontal bar chart of feature importances.

    Compatible with scikit-learn ``feature_importances_`` arrays and
    any non-negative importance measure.

    Parameters
    ----------
    importances : ndarray (n_features,)
    feature_names : list of str or None
    top_n : int
        Show only the top *n* features by importance.
    horizontal : bool
        Use horizontal bars (default) for long feature names.
    ax : Axes or None
    figsize : (width, height) or None
    title : str
    style : bool

    Returns
    -------
    fig : Figure
    """
    importances = np.asarray(importances, dtype=float)
    n_feats = len(importances)

    if feature_names is None:
        feature_names = [f"f{i}" for i in range(n_feats)]

    # Select top N
    order = np.argsort(importances)[::-1][:top_n]
    vals = importances[order]
    names = [feature_names[i] for i in order]

    if figsize is None:
        figsize = (5.0, max(3.0, 0.35 * len(vals)))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    y_pos = np.arange(len(vals))
    color = [EM_COLORS["primary"]] * len(vals)
    color[0] = EM_COLORS["secondary"]  # highlight top feature

    if horizontal:
        ax.barh(y_pos, vals[::-1], color=color[::-1],
                edgecolor="white", lw=0.4)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names[::-1], fontsize=8)
        ax.set_xlabel("Importance", fontsize=9)
    else:
        ax.bar(y_pos, vals, color=color, edgecolor="white", lw=0.4)
        ax.set_xticks(y_pos)
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Importance", fontsize=9)

    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    return fig
