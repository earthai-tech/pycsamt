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

__all__ = [
    "EMStyle",
    "EM_COLORS",
    "EM_CMAPS",
    "EM_FIGSIZE",
    "em_context",
    "StationTickConfig",
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

    def __init__(self, overrides: dict | None = None):
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
# Station-axis tick configuration
# ─────────────────────────────────────────────────────────────────────────────

class StationTickConfig:
    """
    Reusable configuration for station-indexed x-axis tick spacing.

    Prevents label overlap automatically when many stations are displayed.
    All plot functions that draw a station axis accept either individual
    ``tick_label_rotation`` / ``tick_fontsize`` / ``tick_every`` keyword
    arguments *or* a pre-built ``StationTickConfig`` instance.

    Parameters
    ----------
    every : int or ``"auto"``, default ``"auto"``
        Show every *n*-th tick label.  ``"auto"`` picks the smallest
        'nice' step (1, 2, 5, 10, 20, 25, 50, …) that prevents label
        overlap given the figure width and label character length.
    rotation : float, default ``45.0``
        Tick-label rotation in degrees.
    fontsize : int, default ``7``
        Tick-label font size.
    fmt : str, default ``"{}"``
        ``str.format`` template applied to each *visible* label.
        Use ``"{:02d}"`` to zero-pad integer indices, for example.

    Examples
    --------
    >>> from pycsamt.ai.plot._style import StationTickConfig
    >>> cfg = StationTickConfig(every=5, rotation=30, fontsize=8)
    >>> cfg.apply(ax, x_positions, station_labels)

    Automatic mode (default) chooses the step at render time:

    >>> cfg_auto = StationTickConfig()          # every="auto"
    >>> cfg_auto.compute_every(128, figwidth_in=11.0)
    5
    """

    #: Nice step sequence for automatic spacing
    _NICE_STEPS = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000]

    def __init__(
        self,
        every: Union[int, str] = "auto",
        rotation: float = 45.0,
        fontsize: int = 7,
        fmt: str = "{}",
    ) -> None:
        self.every    = every
        self.rotation = rotation
        self.fontsize = fontsize
        self.fmt      = fmt

    # ── public ────────────────────────────────────────────────────────────

    def compute_every(
        self,
        n_stations: int,
        figwidth_in: float = 10.0,
        max_label_len: int = 4,
    ) -> int:
        """
        Return the minimum tick step that prevents label overlap.

        Parameters
        ----------
        n_stations : int
        figwidth_in : float
            Total figure width in inches.
        max_label_len : int
            Character length of the widest label string.

        Returns
        -------
        step : int
        """
        if isinstance(self.every, int):
            return max(1, self.every)

        # Available horizontal space (leave ≈15 % for margins / colorbars)
        avail_in = figwidth_in * 0.82
        space_per_station = avail_in / max(n_stations, 1)

        # Projected label footprint at given rotation
        # char_width ≈ fontsize × 0.006 in (empirical for 7-8 pt sans-serif)
        char_w   = self.fontsize * 0.006
        r        = abs(float(self.rotation))
        r_rad    = r * 3.14159 / 180.0
        lbl_w    = char_w * max(max_label_len, 2)
        lbl_h    = self.fontsize * 0.014           # line height
        eff_w    = lbl_w * abs(__import__("math").sin(r_rad)) \
                 + lbl_h * abs(__import__("math").cos(r_rad)) \
                 + 0.025   # small padding

        if space_per_station >= eff_w:
            return 1

        needed = eff_w / space_per_station
        for step in self._NICE_STEPS:
            if step >= needed:
                return step
        return int(needed) + 1

    def apply(
        self,
        ax: plt.Axes,
        positions: np.ndarray,
        labels: List[str],
        xlabel: str = "",
        xlim: Tuple[float, float] | None = None,
    ) -> None:
        """
        Set x-ticks on *ax*, showing only every ``compute_every``-th label.

        Invisible ticks are set to an empty string so the tick mark is
        still drawn; only the label is hidden.

        Parameters
        ----------
        ax : Axes
        positions : ndarray
            x-positions of the ticks (usually ``np.arange(n_st)``).
        labels : list of str
            Full label for every tick position.
        xlabel : str
            Optional x-axis label text.
        xlim : (lo, hi) or None
            Passed to ``ax.set_xlim`` when provided.
        """
        import matplotlib.pyplot as _plt  # noqa: F401 — ensure plt available

        n  = len(labels)
        fw = ax.figure.get_figwidth() if ax.figure is not None else 10.0
        ml = max((len(str(lbl)) for lbl in labels), default=4)
        every = self.compute_every(n, fw, ml)

        tick_labels = [
            self.fmt.format(lbl) if (i % every == 0) else ""
            for i, lbl in enumerate(labels)
        ]
        ha = "right" if self.rotation > 20 else "center"
        ax.set_xticks(positions)
        ax.set_xticklabels(
            tick_labels,
            rotation=self.rotation,
            ha=ha,
            fontsize=self.fontsize,
        )
        if xlim is not None:
            ax.set_xlim(*xlim)
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=8)

    def __repr__(self) -> str:
        return (
            f"StationTickConfig(every={self.every!r}, "
            f"rotation={self.rotation}, fontsize={self.fontsize})"
        )


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
