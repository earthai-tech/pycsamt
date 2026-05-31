# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Shared matplotlib style contract for all pycsamt.ai plots.

:class:`EMStyle` is a context manager (and decorator) that applies
a consistent set of ``rcParams`` for every plot in this package.
All plotting functions in :mod:`pycsamt.ai.plot` wrap their content
with this context so figures look identical whether generated
interactively or in batch.

Style principles
----------------
* Resistivity sections — ``RdYlBu_r`` (diverging) with log-norm
* Uncertainty overlays — ``YlOrRd`` (sequential)
* Misfit / residuals — ``RdBu_r``
* Primary line colour — ``#2166ac`` (steel blue)
* Secondary line colour — ``#d6604d`` (terracotta)
* Figure sizes:
    * Single-column  3.5 × 4.5 in
    * Double-column  7.0 × 4.5 in
    * Wide section   10.0 × 4.0 in
* Typography: axis labels 11 pt, title 13 pt, ticks 9 pt
* White background with light-grey grid (#edededed, lw 0.7)
* DPI 120 for screen; pass ``dpi=300`` to ``savefig`` for print

Usage
-----
>>> from pycsamt.ai.plot._style import EMStyle
>>> with EMStyle():
...     fig, ax = plt.subplots()
...     ax.plot(...)

Alternatively use as a decorator:

>>> @EMStyle()
... def make_figure():
...     ...
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Optional

__all__ = [
    "EMStyle",
    "EM_COLORS",
    "EM_CMAPS",
    "EM_FIGSIZE",
    "em_context",
]

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

EM_COLORS = {
    "primary": "#2166ac",      # steel blue — predicted model, main lines
    "secondary": "#d6604d",    # terracotta — true model, reference lines
    "true": "#1b7837",         # dark green — ground truth
    "pred": "#762a83",         # purple — predicted
    "error": "#b2182b",        # red — large error / warning
    "background": "white",
    "grid": "#ededed",
    "text": "#1a1a1a",
}

EM_CMAPS = {
    "resistivity": "RdYlBu_r",    # diverging, log-scaled
    "uncertainty": "YlOrRd",       # sequential, low-to-high spread
    "misfit": "RdBu_r",            # diverging, centred on zero
    "conductivity": "RdYlBu",      # inverted for conductivity maps
}

EM_FIGSIZE = {
    "single": (3.5, 4.5),    # single-column journal figure
    "double": (7.0, 4.5),    # double-column / two-panel
    "wide": (10.0, 4.0),     # wide section / profile plot
    "square": (5.0, 5.0),    # scatter / calibration plots
    "tall": (3.5, 6.0),      # deep 1-D profile
}

_EM_RC = {
    # Typography
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 13,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.title_fontsize": 10,
    # Lines
    "lines.linewidth": 1.5,
    "lines.markersize": 5,
    "patch.linewidth": 0.8,
    # Axes
    "axes.facecolor": "white",
    "axes.edgecolor": "#555555",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#ededed",
    "grid.linewidth": 0.7,
    "grid.alpha": 1.0,
    # Figure
    "figure.facecolor": "white",
    "figure.dpi": 120,
    "figure.autolayout": False,
    "savefig.bbox": "tight",
    "savefig.dpi": 200,
    # Colourbar
    "image.cmap": "RdYlBu_r",
    # Math
    "mathtext.fontset": "dejavusans",
}


# ─────────────────────────────────────────────────────────────────────────────
# Context manager / decorator
# ─────────────────────────────────────────────────────────────────────────────

class EMStyle:
    """
    Context manager that applies the shared EM plot style.

    Parameters
    ----------
    overrides : dict or None
        Additional ``rcParams`` to merge on top of the base style.

    Examples
    --------
    >>> with EMStyle():
    ...     fig, ax = plt.subplots(figsize=EM_FIGSIZE['double'])
    ...     ax.semilogy(period, rho_a)
    """

    def __init__(self, overrides: Optional[dict] = None):
        self._overrides = overrides or {}

    def __enter__(self):
        import matplotlib as mpl
        self._rc_backup = dict(mpl.rcParams)
        rc = dict(_EM_RC)
        rc.update(self._overrides)
        mpl.rcParams.update(rc)
        return self

    def __exit__(self, *_):
        import matplotlib as mpl
        mpl.rcParams.update(self._rc_backup)

    def __call__(self, func):
        """Use as a decorator."""
        from functools import wraps

        @wraps(func)
        def wrapper(*args, **kwargs):
            with self.__class__(self._overrides):
                return func(*args, **kwargs)
        return wrapper


@contextmanager
def em_context(**overrides):
    """
    Lightweight context manager shorthand.

    >>> with em_context(figure_dpi=150):
    ...     make_plot()
    """
    with EMStyle(overrides):
        yield


# ─────────────────────────────────────────────────────────────────────────────
# Shared colorbar helper
# ─────────────────────────────────────────────────────────────────────────────

def add_colorbar(
    mappable,
    ax,
    label: str = r"$\log_{10}(\rho)$ (Ω·m)",
    *,
    size: str = "4%",
    pad: float = 0.08,
):
    """
    Attach a right-side colorbar to *ax* using ``make_axes_locatable``.

    Parameters
    ----------
    mappable
        The ``ScalarMappable`` (output of ``imshow``, ``pcolormesh``, etc.)
    ax
        Host axes.
    label : str
        Colorbar label.
    size : str
        Width as a percentage of the parent axes.
    pad : float
        Padding between axes and colorbar.

    Returns
    -------
    cbar : Colorbar
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=size, pad=pad)
    import matplotlib.pyplot as plt
    cbar = plt.colorbar(mappable, cax=cax)
    cbar.set_label(label, fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    return cbar
