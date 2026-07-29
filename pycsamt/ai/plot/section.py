# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
2-D resistivity section and pseudo-section visualisation.

All functions respect the :class:`~pycsamt.ai.plot._style.EMStyle`
publication conventions and return :class:`~matplotlib.figure.Figure`
or :class:`~matplotlib.axes.Axes` objects for downstream composition.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ._style import (
    EM_CMAPS,
    EM_COLORS,
    EM_FIGSIZE,
    EMStyle,
    add_colorbar,
)

__all__ = [
    "plot_section",
    "plot_section_pair",
    "plot_pseudo_section",
]

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────


def _default_depths(n_depth: int, depth_max: float = 2000.0) -> np.ndarray:
    return np.linspace(0.0, depth_max, n_depth + 1)


def _default_stations(n_stations: int, spacing: float = 1.0) -> np.ndarray:
    return np.arange(n_stations) * spacing


def _section_im(
    ax: Axes,
    rho_2d: np.ndarray,
    depths: np.ndarray,
    stations: np.ndarray,
    *,
    log_scale: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str = "RdYlBu_r",
    aspect: str = "auto",
) -> Any:
    """Draw a filled 2-D resistivity section on *ax*."""
    data = rho_2d.copy()
    if log_scale:
        data = np.log10(np.maximum(data, 1e-6))
        cbar_label = r"$\log_{10}(\rho)$ (Ω·m)"
    else:
        cbar_label = r"$\rho$ (Ω·m)"

    if vmin is None:
        vmin = np.nanpercentile(data, 2)
    if vmax is None:
        vmax = np.nanpercentile(data, 98)

    # pcolormesh: (n_depth+1, n_sta+1) edges, (n_depth, n_sta) values
    X, Y = np.meshgrid(
        stations,
        depths[:-1] if len(depths) == rho_2d.shape[0] + 1 else depths,
    )
    im = ax.pcolormesh(
        X,
        Y,
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="auto",
    )
    ax.set_aspect(aspect)
    return im, cbar_label


# ─────────────────────────────────────────────────────────────────────────────
# plot_section
# ─────────────────────────────────────────────────────────────────────────────


@EMStyle()
def plot_section(
    rho_2d: np.ndarray,
    *,
    depths: np.ndarray | None = None,
    stations: np.ndarray | None = None,
    depth_max: float = 2000.0,
    station_spacing: float = 1.0,
    log_scale: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | None = None,
    title: str = "",
    xlabel: str = "Station",
    ylabel: str = "Depth (m)",
    show_sites: bool = True,
    figsize: tuple[float, float] | None = None,
    ax: Axes | None = None,
    style: bool = True,
) -> Figure:
    """
    Plot a 2-D resistivity section.

    Parameters
    ----------
    rho_2d : ndarray, shape (n_depth, n_stations)
        Resistivity or log₁₀(ρ) 2-D model.
    depths : ndarray (n_depth,) or None
        Depth values in metres.  Default: linear 0 → ``depth_max``.
    stations : ndarray (n_stations,) or None
        Station positions (arbitrary units).
    depth_max : float
        Maximum depth used when ``depths`` is ``None``.
    station_spacing : float
        Station interval when ``stations`` is ``None``.
    log_scale : bool
        Apply log₁₀ transform before plotting.  Set ``False`` if
        ``rho_2d`` already contains log₁₀(ρ).
    vmin, vmax : float or None
        Colour scale limits.
    cmap : str or None
        Matplotlib colormap name.  Defaults to ``'RdYlBu_r'``.
    title : str
    xlabel, ylabel : str
    show_sites : bool
        Draw site-marker triangles at the surface.
    figsize : (width, height) or None
    ax : Axes or None
    style : bool
        Apply :class:`~pycsamt.ai.plot._style.EMStyle` context.

    Returns
    -------
    fig : Figure
    """
    n_depth, n_stations = rho_2d.shape

    if depths is None:
        depths = _default_depths(n_depth, depth_max)
    if stations is None:
        stations = _default_stations(n_stations, station_spacing)
    if cmap is None:
        cmap = EM_CMAPS["resistivity"]
    if figsize is None:
        figsize = EM_FIGSIZE["wide"]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    im, cbar_label = _section_im(
        ax,
        rho_2d,
        depths,
        stations,
        log_scale=log_scale,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
    )

    if show_sites:
        ax.plot(
            stations,
            np.full_like(stations, depths.min()),
            "v",
            ms=4,
            color=EM_COLORS["text"],
            zorder=5,
        )

    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=10)

    add_colorbar(im, ax, label=cbar_label)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# plot_section_pair
# ─────────────────────────────────────────────────────────────────────────────


@EMStyle()
def plot_section_pair(
    true_2d: np.ndarray,
    pred_2d: np.ndarray,
    *,
    depths: np.ndarray | None = None,
    stations: np.ndarray | None = None,
    depth_max: float = 2000.0,
    station_spacing: float = 1.0,
    log_scale: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | None = None,
    show_difference: bool = True,
    figsize: tuple[float, float] | None = None,
    style: bool = True,
) -> Figure:
    """
    Side-by-side comparison of true and predicted 2-D resistivity sections.

    Optionally shows the absolute difference (``show_difference=True``).

    Parameters
    ----------
    true_2d : ndarray (n_depth, n_stations)
    pred_2d : ndarray (n_depth, n_stations)
    depths, stations, depth_max, station_spacing, log_scale : see :func:`plot_section`
    vmin, vmax : float or None
        Shared colour scale.  When ``None``, computed from *true_2d*.
    cmap : str or None
    show_difference : bool
        Add a third panel with the signed difference.
    figsize : (width, height) or None
    style : bool

    Returns
    -------
    fig : Figure
    """
    n_depth, n_stations = true_2d.shape

    if depths is None:
        depths = _default_depths(n_depth, depth_max)
    if stations is None:
        stations = _default_stations(n_stations, station_spacing)
    if cmap is None:
        cmap = EM_CMAPS["resistivity"]

    # Compute in log space for consistent comparison
    true_log = np.log10(np.maximum(true_2d, 1e-6)) if log_scale else true_2d
    pred_log = np.log10(np.maximum(pred_2d, 1e-6)) if log_scale else pred_2d
    cbar_label = r"$\log_{10}(\rho)$ (Ω·m)" if log_scale else r"$\rho$ (Ω·m)"

    if vmin is None:
        vmin = np.nanpercentile(true_log, 2)
    if vmax is None:
        vmax = np.nanpercentile(true_log, 98)

    n_cols = 3 if show_difference else 2
    if figsize is None:
        figsize = (n_cols * 5.0, 4.0)

    fig, axes = plt.subplots(1, n_cols, figsize=figsize, sharey=True)

    for ax, data, ttl in zip(
        axes[:2],
        [true_log, pred_log],
        ["True model", "Predicted model"],
    ):
        X, Y = np.meshgrid(stations, depths[:n_depth])
        im = ax.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        ax.invert_yaxis()
        ax.set_title(ttl, fontsize=9)
        ax.set_xlabel("Station")

    axes[0].set_ylabel("Depth (m)")
    add_colorbar(im, axes[1], label=cbar_label)

    if show_difference:
        diff = pred_log - true_log
        dmax = np.nanpercentile(np.abs(diff), 95)
        X, Y = np.meshgrid(stations, depths[:n_depth])
        im_d = axes[2].pcolormesh(
            X,
            Y,
            diff,
            cmap=EM_CMAPS["misfit"],
            vmin=-dmax,
            vmax=dmax,
            shading="auto",
        )
        axes[2].invert_yaxis()
        axes[2].set_title("Difference", fontsize=9)
        axes[2].set_xlabel("Station")
        add_colorbar(im_d, axes[2], label=r"$\Delta\log_{10}(\rho)$")

    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# plot_pseudo_section
# ─────────────────────────────────────────────────────────────────────────────


@EMStyle()
def plot_pseudo_section(
    rho_a_2d: np.ndarray,
    *,
    freqs: np.ndarray | None = None,
    stations: np.ndarray | None = None,
    station_spacing: float = 1.0,
    log_freq: bool = True,
    log_rho: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | None = None,
    component: str = "xy",
    title: str = "",
    figsize: tuple[float, float] | None = None,
    ax: Axes | None = None,
    style: bool = True,
) -> Figure:
    """
    Apparent-resistivity or phase pseudo-section plot.

    Parameters
    ----------
    rho_a_2d : ndarray (n_freqs, n_stations)
        Apparent resistivity or phase values.
    freqs : ndarray (n_freqs,) or None
        Frequencies in Hz.
    stations : ndarray (n_stations,) or None
        Station positions.
    station_spacing : float
    log_freq : bool
        Use log₁₀(period) for the y-axis.
    log_rho : bool
        Apply log₁₀ transform to the data.
    vmin, vmax : float or None
    cmap : str or None
    component : str
        Label for the colour bar (e.g. ``'xy'``, ``'yx'``, ``'phase_xy'``).
    title : str
    figsize : (width, height) or None
    ax : Axes or None
    style : bool

    Returns
    -------
    fig : Figure
    """
    n_freqs, n_stations = rho_a_2d.shape

    if stations is None:
        stations = _default_stations(n_stations, station_spacing)
    if freqs is None:
        freqs = np.logspace(-3, 3, n_freqs)
    if cmap is None:
        cmap = EM_CMAPS["resistivity"]
    if figsize is None:
        figsize = EM_FIGSIZE["wide"]

    y_axis = np.log10(1.0 / freqs) if log_freq else freqs

    data = rho_a_2d.copy()
    if log_rho:
        data = np.log10(np.maximum(data, 1e-6))
        cbar_label = rf"$\log_{{10}}(\rho_{{a,{component}}})$ (Ω·m)"
    else:
        cbar_label = rf"$\rho_{{a,{component}}}$ (Ω·m)"

    if vmin is None:
        vmin = np.nanpercentile(data, 2)
    if vmax is None:
        vmax = np.nanpercentile(data, 98)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    X, Y = np.meshgrid(stations, y_axis)
    im = ax.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")

    ax.set_xlabel("Station")
    ax.set_ylabel(r"$\log_{10}(T)$ (s)" if log_freq else "Frequency (Hz)")
    if title:
        ax.set_title(title, fontsize=10)

    add_colorbar(im, ax, label=cbar_label)
    fig.tight_layout()
    return fig
