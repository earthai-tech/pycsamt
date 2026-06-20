# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.plot.utils
==================
Generic 2-D pseudosection plotter used by emtools functions
(plot_skew_2d, plot_sensitivity_depth_section, …).
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm

from ..utils.plot import _get_xticks_formatage


def plot2d(
    ar,
    *,
    y=None,
    x=None,
    distance: float = 50.0,
    stnlist: Optional[List[str]] = None,
    prefix: str = "S",
    how: str = "py",
    to_log10: bool = False,
    plot_contours: bool = False,
    top_label: str = "",
    cb_label: str = "",
    rotate_xlabel: float = 0.0,
    cmap: str = "jet_r",
    plt_style: str = "pcolormesh",
    fig_size: Tuple[float, float] = (8.0, 4.0),
    fig_dpi: int = 150,
    font_size: float = 7.0,
    ax: Optional[plt.Axes] = None,
    **_extra,          # absorb legacy BasePlot kwargs without error
) -> plt.Axes:
    """Plot a 2-D pseudosection from a matrix.

    Parameters
    ----------
    ar : array-like, shape (n_rows, n_cols)
        Data matrix.  Rows → y-axis (frequency/period), cols → stations.
    y : array-like, optional
        Y-axis coordinates (length == n_rows).
    x : array-like, optional
        X-axis coordinates (length == n_cols).  If omitted, derived from
        *distance* and column count.
    distance : float
        Station spacing in metres used when *x* is None.
    stnlist : list of str, optional
        Station labels for the top twin-axis.
    prefix : str
        Auto-label prefix when *stnlist* is None (e.g. ``'S'`` → S00, S01 …).
    how : str
        ``'py'`` starts counting from 0, anything else from 1.
    to_log10 : bool
        Convert *ar* (and *y*) to log10 before plotting.
    plot_contours : bool
        Overlay filled contours (pcolormesh mode only).
    top_label : str
        Label for the twin top-axis.
    cb_label : str
        Colorbar label.
    rotate_xlabel : float
        Rotation angle for the station tick labels.
    cmap : str
        Matplotlib colormap name.
    plt_style : {'pcolormesh', 'imshow'}
        Rendering backend.
    fig_size : tuple of float
        Figure size ``(width, height)`` in inches.
    fig_dpi : int
        Figure resolution.
    font_size : float
        Base font size in points.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw into.  A new figure is created when *None*.

    Returns
    -------
    matplotlib.axes.Axes
    """
    ar = np.asarray(ar, dtype=float)
    if ar.ndim != 2:
        raise ValueError(f"plot2d expects a 2-D array; got shape {ar.shape}")

    n_rows, n_cols = ar.shape

    # Optional log10 transform
    if to_log10:
        ar = np.log10(ar)
        if y is not None:
            y = np.log10(np.asarray(y, dtype=float))

    # Axis coordinates
    y_vec = np.arange(n_rows, dtype=float) if y is None else np.asarray(y, dtype=float)
    x_vec = (np.arange(n_cols) * float(distance)) if x is None else np.asarray(x, dtype=float)

    # Station labels
    if stnlist is None:
        start = 0 if how == "py" else 1
        stnlist = [f"{prefix}{i:02d}" for i in range(start, start + n_cols)]

    # --- Create axes if needed ---
    if ax is None:
        fig, axe = plt.subplots(1, figsize=fig_size, dpi=fig_dpi)
    else:
        axe = ax
        fig = axe.get_figure()

    # --- Plot ---
    if plt_style == "imshow":
        mappable = axe.imshow(
            ar,
            interpolation="nearest",
            cmap=cmap,
            aspect="auto",
            origin="lower",
            extent=(x_vec.min(), x_vec.max(), y_vec.min(), y_vec.max()),
        )
    else:  # pcolormesh (default)
        X, Y = np.meshgrid(x_vec, y_vec)
        vmin, vmax = np.nanmin(ar), np.nanmax(ar)

        if plot_contours:
            levels = mticker.MaxNLocator(nbins=15).tick_values(vmin, vmax)
            norm = BoundaryNorm(levels, ncolors=plt.get_cmap(cmap).N, clip=True)
            mappable = axe.pcolormesh(X, Y, np.flipud(ar), cmap=cmap, norm=norm,
                                      shading="auto")
            dx = dy = 0.05
            axe.contourf(X + dx / 2, Y + dy / 2, np.flipud(ar),
                         levels=levels, cmap=cmap)
        else:
            mappable = axe.pcolormesh(X, Y, np.flipud(ar), cmap=cmap,
                                      vmin=vmin, vmax=vmax, shading="auto")

    axe.set_xlim(x_vec.min(), x_vec.max())
    axe.set_ylim(y_vec.min(), y_vec.max())
    axe.set_xlabel("Distance (m)", fontsize=1.5 * font_size)
    axe.tick_params(labelsize=font_size)

    # Colorbar
    cbar = fig.colorbar(mappable, ax=axe)
    cbar.ax.tick_params(axis="y", direction="in", labelsize=font_size)
    if cb_label:
        cbar.set_label(cb_label, fontsize=1.2 * font_size)

    # Twin top-axis with station labels aligned to real x positions.
    # Subsample when there are many stations so labels don't overlap.
    axe2 = axe.twiny()
    axe2.set_xlim(axe.get_xlim())
    _max_lbl = 20
    if n_cols > _max_lbl:
        _step = max(1, n_cols // _max_lbl)
        _tick_x   = x_vec[::_step]
        _tick_lbl = [stnlist[i] for i in range(0, n_cols, _step)]
    else:
        _tick_x   = x_vec
        _tick_lbl = list(stnlist)
    axe2.set_xticks(_tick_x)
    axe2.set_xticklabels(_tick_lbl,
                         rotation=rotate_xlabel if rotate_xlabel else 90,
                         fontsize=font_size, ha="left")
    axe2.tick_params(labelsize=font_size)
    if top_label:
        axe2.set_xlabel(top_label, fontsize=1.5 * font_size)

    fig.tight_layout()
    return axe
