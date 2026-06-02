# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Visualisation for pyCSAMT 1-D and 2-D forward models and responses.

All plot functions read visual parameters from the package-wide API singletons
at call time, so a single :func:`~pycsamt.api.style.configure_style` or
:func:`~pycsamt.api.style.use_style` call propagates to every figure.

Functions
---------
**1-D forward**

.. autosummary::

   plot_response_1d
   plot_model_1d
   plot_response_and_model_1d

**2-D forward**

.. autosummary::

   plot_model_2d
   plot_pseudosection_2d
   plot_response_profiles

Quick start
-----------
Single 1-D sounding::

    from pycsamt.forward import MT1DForward, LayeredModel
    from pycsamt.forward.plot import plot_response_and_model_1d
    import numpy as np

    model = LayeredModel([100, 10, 500], [300, 800])
    resp  = MT1DForward(np.logspace(-3, 3, 30)).run(model)

    fig = plot_response_and_model_1d(resp, model, title="My 1-D model")
    fig.savefig("model_1d.png", dpi=150, bbox_inches="tight")

2-D pseudo-section::

    from pycsamt.forward.plot import plot_pseudosection_2d, plot_model_2d

    fig1 = plot_model_2d(grid)
    fig2 = plot_pseudosection_2d(resp2d, mode="both")
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from ..api._rose_style import _UNSET
from ..api.style import PYCSAMT_STYLE
from ..api.control import PYCSAMT_CONTROL
from ..api.station import PYCSAMT_STATION_RENDERING

__all__ = [
    "plot_response_1d",
    "plot_model_1d",
    "plot_response_and_model_1d",
    "plot_model_2d",
    "plot_pseudosection_2d",
    "plot_response_profiles",
]


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _x_vals(freqs: np.ndarray) -> np.ndarray:
    """Return x-axis values from frequency using the active control."""
    return PYCSAMT_CONTROL.x.transform(freqs)


def _x_label() -> str:
    return PYCSAMT_CONTROL.x.label()


def _x_log() -> bool:
    return PYCSAMT_CONTROL.x.use_log_scale()


def _rho_vals(rho: np.ndarray) -> np.ndarray:
    return PYCSAMT_CONTROL.rho.transform(rho)


def _rho_label() -> str:
    return PYCSAMT_CONTROL.rho.label()


def _phase_vals(phase: np.ndarray) -> np.ndarray:
    return PYCSAMT_CONTROL.phase.transform(phase)


def _phase_label() -> str:
    return PYCSAMT_CONTROL.phase.label()


def _spine_style(ax: Axes) -> None:
    """Apply a clean spine / grid style consistent with pyCSAMT figures."""
    ax.grid(True, which="both", ls=":", lw=0.4, color="0.75", zorder=0)
    ax.set_axisbelow(True)


def _add_colorbar(
    fig: Figure,
    ax: Axes,
    mappable,
    label: str,
    fontsize: float = 8,
    pad: float = 0.03,
    shrink: float = 0.95,
) -> None:
    cb = fig.colorbar(mappable, ax=ax, pad=pad, shrink=shrink, aspect=25)
    cb.set_label(label, fontsize=fontsize)
    cb.ax.tick_params(labelsize=fontsize - 0.5)


# ─────────────────────────────────────────────────────────────────────────────
# 1-D response plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_response_1d(
    response,
    *,
    modes: Union[str, Sequence[str]] = "both",
    show_te: bool = _UNSET,
    show_tm: bool = _UNSET,
    color_te=_UNSET,
    color_tm=_UNSET,
    lw: float = _UNSET,
    marker_te=_UNSET,
    marker_tm=_UNSET,
    ms: float = _UNSET,
    label_te: str = _UNSET,
    label_tm: str = _UNSET,
    title: str = "",
    figsize: Tuple[float, float] = (7, 5.5),
    axes=None,
) -> np.ndarray:
    """Plot 1-D MT/CSAMT apparent resistivity and phase vs period (or frequency).

    Reads visual defaults from :data:`~pycsamt.api.style.PYCSAMT_STYLE` and
    axis behaviour from :data:`~pycsamt.api.control.PYCSAMT_CONTROL`.

    Parameters
    ----------
    response : ForwardResponse
        Output of :class:`~pycsamt.forward.em1d.MT1DForward` or
        :class:`~pycsamt.forward.em1d.CSAMT1DForward`.
    modes : {'te', 'tm', 'both'}
        Which modes to plot.  For a 1-D response ``response.rho_a`` and
        ``response.phase`` hold a single polarisation; pass ``'both'`` to
        show both ρ_a curves with TE/TM styling on the same axes.
    show_te, show_tm : bool
        Override visibility of each mode independently.
    color_te, color_tm : colour spec
        Line colours.  Default: ``PYCSAMT_STYLE.mt.te.color`` /
        ``PYCSAMT_STYLE.mt.tm.color``.
    lw : float
        Line width.
    marker_te, marker_tm : str
        Marker styles.
    ms : float
        Marker size.
    label_te, label_tm : str
        Legend labels.
    title : str
        Figure suptitle.
    figsize : (float, float)
    axes : array-like of Axes or None
        Two pre-existing axes (rho, phase).  Created when not given.

    Returns
    -------
    axes : ndarray of Axes, shape (2,)
        ``[ax_rho, ax_phase]``
    """
    _st = PYCSAMT_STYLE.mt

    # Resolve defaults
    if color_te  is _UNSET: color_te  = _st.te.color
    if color_tm  is _UNSET: color_tm  = _st.tm.color
    if lw        is _UNSET: lw        = _st.te.lw
    if marker_te is _UNSET: marker_te = _st.te.marker
    if marker_tm is _UNSET: marker_tm = _st.tm.marker
    if ms        is _UNSET: ms        = _st.te.ms
    if label_te  is _UNSET: label_te  = _st.te.label or "TE"
    if label_tm  is _UNSET: label_tm  = _st.tm.label or "TM"

    freqs = response.freqs
    x     = _x_vals(freqs)

    # Handle both 1D (rho_a shape (nf,)) and 2D responses
    rho_a = np.asarray(response.rho_a)
    phase = np.asarray(response.phase)
    if rho_a.ndim > 1:
        rho_a = rho_a[:, 0]
        phase = phase[:, 0]

    if axes is None:
        fig, axs = plt.subplots(2, 1, figsize=figsize, sharex=True,
                                constrained_layout=True)
    else:
        axs = np.asarray(axes).ravel()
        fig = axs[0].get_figure()

    ax_r, ax_p = axs[0], axs[1]

    kw_te = dict(color=color_te, lw=lw, marker=marker_te, ms=ms,
                 mfc="white", mew=1.0, alpha=_st.te.alpha)
    kw_tm = dict(color=color_tm, lw=lw, marker=marker_tm, ms=ms,
                 mfc="white", mew=1.0, alpha=_st.tm.alpha)

    modes_str = modes if isinstance(modes, str) else ",".join(modes)
    plot_te = (modes_str in ("te", "both") and show_te is not False)
    plot_tm = (modes_str in ("tm", "both") and show_tm is not False)

    if plot_te:
        ax_r.plot(x, _rho_vals(rho_a), label=label_te, **kw_te)
        ax_p.plot(x, _phase_vals(phase), label=label_tm, **kw_te)
    if plot_tm and hasattr(response, "rho_a_tm"):
        rho_tm  = np.asarray(response.rho_a_tm)
        phi_tm  = np.asarray(response.phase_tm)
        if rho_tm.ndim > 1:
            rho_tm = rho_tm[:, 0]; phi_tm = phi_tm[:, 0]
        ax_r.plot(x, _rho_vals(rho_tm), label=label_tm, **kw_tm)
        ax_p.plot(x, _phase_vals(phi_tm), **kw_tm)

    for ax in (ax_r, ax_p):
        _spine_style(ax)
        if _x_log():
            ax.set_xscale("log")

    ax_r.set_ylabel(_rho_label(), fontsize=9)
    ax_p.set_ylabel(_phase_label(), fontsize=9)
    ax_p.set_xlabel(_x_label(), fontsize=9)

    if plot_te or plot_tm:
        ax_r.legend(fontsize=8, framealpha=0.8)

    if title:
        fig.suptitle(title, fontsize=10, y=1.01)

    return np.array([ax_r, ax_p])


# ─────────────────────────────────────────────────────────────────────────────
# 1-D model plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_model_1d(
    models,
    labels: Optional[Sequence[str]] = None,
    *,
    log_rho: bool = True,
    depth_max: Optional[float] = None,
    lw: float = _UNSET,
    alpha: float = _UNSET,
    title: str = "",
    figsize: Tuple[float, float] = (3.8, 5.5),
    ax: Optional[Axes] = None,
) -> Axes:
    """Plot one or more 1-D layered earth models as resistivity-depth profiles.

    Multiple models are coloured with :attr:`PYCSAMT_STYLE.multiline` so the
    gradient or cycle palette matches all other multi-model plots.

    Parameters
    ----------
    models : LayeredModel or list of LayeredModel
        One model or a list to overlay.
    labels : list of str, optional
        Legend labels; defaults to model name attributes.
    log_rho : bool
        Use log₁₀ scale on the resistivity axis.
    depth_max : float or None
        Maximum depth shown [m].  Defaults to the deepest interface × 1.2.
    lw : float
        Line width.  Default: ``PYCSAMT_STYLE.multiline.lw``.
    alpha : float
        Line alpha.  Default: ``PYCSAMT_STYLE.multiline.alpha``.
    title : str
    figsize : (float, float)
    ax : Axes or None

    Returns
    -------
    ax : Axes
    """
    from .synthetic import LayeredModel

    _ml = PYCSAMT_STYLE.multiline
    if lw    is _UNSET: lw    = _ml.lw
    if alpha is _UNSET: alpha = _ml.alpha

    if isinstance(models, LayeredModel):
        models = [models]
    n = len(models)
    colors = _ml.colors(n) if n > 1 else [PYCSAMT_STYLE.mt.xy.color]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    for k, model in enumerate(models):
        rho   = model.resistivity
        thick = model.thickness
        depth_top = np.concatenate([[0.0], np.cumsum(thick)])
        d_max = depth_max or (float(depth_top[-1] + thick[-1]) * 1.2
                               if len(thick) else float(rho[0]) * 2)

        depth_bot = np.concatenate([depth_top[1:], [d_max]])
        xs = np.repeat(rho, 2)
        ys = np.empty(2 * len(rho))
        ys[0::2] = depth_top
        ys[1::2] = depth_bot

        lab = (labels[k] if labels and k < len(labels) else
               model.name or f"model {k + 1}")
        ax.plot(xs, ys, color=colors[k], lw=lw, alpha=alpha,
                label=lab if n > 1 else None)

    ax.invert_yaxis()
    if log_rho:
        ax.set_xscale("log")
    ax.set_xlabel(r"Resistivity  ($\Omega\cdot$m)", fontsize=9)
    ax.set_ylabel("Depth  (m)", fontsize=9)
    if n > 1:
        ax.legend(fontsize=8, framealpha=0.8)
    if title:
        ax.set_title(title, fontsize=10)
    _spine_style(ax)
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# 1-D composite (model + response side-by-side)
# ─────────────────────────────────────────────────────────────────────────────

def plot_response_and_model_1d(
    response,
    model=None,
    *,
    title: str = "",
    figsize: Tuple[float, float] = (10, 5),
    gridspec_kw: Optional[dict] = None,
) -> Figure:
    """3-panel figure combining model depth profile, ρ_a, and phase.

    This is the canonical "validate and save" view for a 1-D forward run.

    Parameters
    ----------
    response : ForwardResponse
    model : LayeredModel or None
        If given, plotted in the left panel.  When ``None``, only the
        two response panels are shown.
    title : str
        Figure suptitle.
    figsize : (float, float)
    gridspec_kw : dict or None
        Passed to :func:`~matplotlib.pyplot.subplots`.

    Returns
    -------
    Figure
    """
    has_model = model is not None
    ncols = 3 if has_model else 2

    gs_kw = gridspec_kw or (
        {"width_ratios": [1, 2, 2]} if has_model else {"width_ratios": [1, 1]}
    )
    fig, axs = plt.subplots(
        1, ncols, figsize=figsize,
        gridspec_kw=gs_kw, constrained_layout=True,
    )

    if has_model:
        ax_m, ax_r, ax_p = axs
        plot_model_1d(model, ax=ax_m)
        ax_m.set_title("Earth model", fontsize=9, pad=6)
    else:
        ax_r, ax_p = axs

    freqs  = response.freqs
    x      = _x_vals(freqs)
    rho_a  = np.asarray(response.rho_a)
    phase  = np.asarray(response.phase)
    if rho_a.ndim > 1:
        rho_a = rho_a[:, 0]; phase = phase[:, 0]

    _st = PYCSAMT_STYLE.mt

    kw = _st.te.plot_kwargs(label="TE")
    ax_r.plot(x, _rho_vals(rho_a), **kw)
    ax_p.plot(x, _phase_vals(phase), **kw)

    for ax in (ax_r, ax_p):
        _spine_style(ax)
        if _x_log():
            ax.set_xscale("log")
        ax.set_xlabel(_x_label(), fontsize=9)

    ax_r.set_ylabel(_rho_label(), fontsize=9)
    ax_r.set_title(r"Apparent resistivity  $\rho_a$", fontsize=9, pad=6)
    ax_p.set_ylabel(_phase_label(), fontsize=9)
    ax_p.set_title("Impedance phase", fontsize=9, pad=6)

    if title:
        fig.suptitle(title, fontsize=11, y=1.02)

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 2-D model view
# ─────────────────────────────────────────────────────────────────────────────

def plot_model_2d(
    grid,
    *,
    cmap: str = "jet_r",
    log_scale: bool = True,
    clip_core: bool = True,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_stations: bool = True,
    station_preset: str = "inversion",
    title: str = "",
    figsize: Tuple[float, float] = (11, 4),
    ax: Optional[Axes] = None,
) -> Axes:
    """Plot the 2-D resistivity model on a colour map.

    Station markers and labels are rendered using
    :data:`~pycsamt.api.station.PYCSAMT_STATION_RENDERING`.

    Parameters
    ----------
    grid : Grid2D
    cmap : str
        Colourmap name.  ``"jet_r"`` is the geophysical convention
        (blue = conductive, red = resistive).
    log_scale : bool
        Display log₁₀(ρ) rather than ρ.
    clip_core : bool
        Clip to the non-padding core region.
    vmin, vmax : float or None
        Colour limits in log₁₀(Ω·m) when *log_scale* is True.
    show_stations : bool
        Render station markers and labels on the surface.
    station_preset : str
        ``"inversion"`` or ``"pseudosection"`` or ``"survey"``.
    title : str
    figsize : (float, float)
    ax : Axes or None

    Returns
    -------
    ax : Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.get_figure()

    p = grid.n_pad
    xs = grid.core_x_slice if clip_core else slice(None)
    zs = grid.core_z_slice if clip_core else slice(None)
    rho = grid.resistivity[zs, xs]
    xn  = (grid.x_nodes[p: grid.nx + 1 - p]
           if clip_core and p > 0 else grid.x_nodes)
    zn  = (grid.z_nodes[: grid.nz + 1 - p]
           if clip_core and p > 0 else grid.z_nodes)

    if log_scale:
        data  = np.log10(np.maximum(rho, 1e-12))
        clabel = r"$\log_{10}\rho$  ($\Omega\cdot$m)"
    else:
        data  = rho
        clabel = r"$\rho$  ($\Omega\cdot$m)"

    pc = ax.pcolormesh(
        xn, zn, data,
        cmap=cmap, shading="flat",
        vmin=vmin, vmax=vmax,
    )
    _add_colorbar(fig, ax, pc, clabel)

    ax.invert_yaxis()
    ax.set_xlabel("Distance  (m)", fontsize=9)
    ax.set_ylabel("Depth  (m)", fontsize=9)

    if show_stations and grid.n_stations > 0:
        # Map station x to the core-clipped coordinate system
        st_x = grid.x_stations
        xlim = (float(xn[0]), float(xn[-1]))
        labels = [f"{i+1}" for i in range(grid.n_stations)]
        sty = PYCSAMT_STATION_RENDERING.style_for(station_preset)
        sty.apply(ax, st_x, labels, xlim=xlim)

    ax.set_title(title or "2-D resistivity model", fontsize=10, pad=6)
    _spine_style(ax)
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# 2-D pseudo-section
# ─────────────────────────────────────────────────────────────────────────────

def plot_pseudosection_2d(
    response,
    *,
    mode: str = "te",
    quantity: str = "rho_a",
    cmap: str = _UNSET,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    n_contours: int = 0,
    show_stations: bool = True,
    station_preset: str = "pseudosection",
    title: str = "",
    figsize: Tuple[float, float] = (11, 5),
    ax: Optional[Axes] = None,
) -> Axes:
    """Plot a 2-D MT pseudo-section (period × station distance).

    The colour shows log₁₀(ρ_a) or phase at each (period, station) point.
    Station markers are rendered via
    :data:`~pycsamt.api.station.PYCSAMT_STATION_RENDERING`.

    Parameters
    ----------
    response : ForwardResponse2D
    mode : {'te', 'tm'}
        Which polarisation mode to display.
    quantity : {'rho_a', 'phase'}
        Quantity to colour.
    cmap : str
        Colour map.  Defaults to ``"jet_r"`` for ρ_a and ``"RdBu_r"``
        for phase.
    vmin, vmax : float or None
        Colour limits.
    n_contours : int
        Number of contour lines overlaid on the colour map (0 = none).
    show_stations : bool
        Add station axis with markers.
    station_preset : str
    title : str
    figsize : (float, float)
    ax : Axes or None

    Returns
    -------
    ax : Axes
    """
    # Resolve the data
    mode = mode.lower()
    quantity = quantity.lower()
    attr = f"rho_a_{mode}" if quantity == "rho_a" else f"phase_{mode}"
    data_raw = getattr(response, attr)   # shape (n_freqs, n_stations)

    freqs   = response.freqs
    st_x    = response.stations_x       # (n_stations,)

    # Build x-axis (station distance), y-axis (log10 period)
    y_vals = np.log10(1.0 / freqs)     # log10(period)
    x_vals = st_x

    # Prepare colour data
    if quantity == "rho_a":
        data_c = np.log10(np.maximum(data_raw, 1e-12))
        default_cmap = "jet_r"
        cb_label     = (r"$\log_{10}\rho_a$  ($\Omega\cdot$m)  "
                        + f"[{mode.upper()}]")
    else:
        data_c = _phase_vals(data_raw)
        default_cmap = "RdBu_r"
        cb_label = _phase_label() + f"  [{mode.upper()}]"

    if cmap is _UNSET:
        cmap = default_cmap

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.get_figure()

    # pcolormesh needs coordinate edges (n+1 values)
    dx_half = np.diff(x_vals) / 2.0 if len(x_vals) > 1 else np.array([50.0])
    x_edges = np.concatenate([
        [x_vals[0] - dx_half[0]],
        x_vals[:-1] + dx_half,
        [x_vals[-1] + dx_half[-1]],
    ])
    dy_half = np.abs(np.diff(y_vals)) / 2.0 if len(y_vals) > 1 else np.array([0.2])
    y_edges = np.concatenate([
        [y_vals[0] - dy_half[0]],
        y_vals[:-1] + (np.sign(np.diff(y_vals)) * dy_half),
        [y_vals[-1] + np.sign(y_vals[-1] - y_vals[-2]) * dy_half[-1]],
    ])

    pc = ax.pcolormesh(
        x_edges, y_edges, data_c,
        cmap=cmap, shading="flat",
        vmin=vmin, vmax=vmax,
    )
    _add_colorbar(fig, ax, pc, cb_label)

    if n_contours > 0:
        ax.contour(
            x_vals, y_vals, data_c,
            levels=n_contours,
            colors="k", linewidths=0.5, alpha=0.6,
        )

    # Y axis
    ax.set_ylabel(r"$\log_{10}T$  (s)", fontsize=9)
    # Make sure higher periods (deeper) are at top (reversed y-axis is NOT
    # standard for pseudo-sections — low period = shallow = top)
    if y_vals[0] > y_vals[-1]:   # frequencies sorted high → low → periods low→high
        ax.invert_yaxis()

    # Station axis via PYCSAMT_STATION_RENDERING
    if show_stations and len(st_x) > 0:
        labels = [f"{i+1}" for i in range(len(st_x))]
        xlim   = (float(x_edges[0]), float(x_edges[-1]))
        sty    = PYCSAMT_STATION_RENDERING.style_for(station_preset)
        sty.apply(ax, st_x, labels, xlim=xlim)
    else:
        ax.set_xlabel("Distance  (m)", fontsize=9)

    mode_label = mode.upper()
    qty_label  = "Apparent resistivity" if quantity == "rho_a" else "Phase"
    ax.set_title(
        title or f"2-D MT pseudo-section — {qty_label} ({mode_label})",
        fontsize=10, pad=6,
    )
    _spine_style(ax)
    return ax


# ─────────────────────────────────────────────────────────────────────────────
# Horizontal response profiles
# ─────────────────────────────────────────────────────────────────────────────

def plot_response_profiles(
    response,
    *,
    mode: str = "te",
    quantity: str = "rho_a",
    freq_indices: Optional[Sequence[int]] = None,
    n_freqs_shown: int = 5,
    lw: float = _UNSET,
    alpha: float = _UNSET,
    title: str = "",
    figsize: Tuple[float, float] = (9, 4),
    ax: Optional[Axes] = None,
) -> Axes:
    """Plot ρ_a (or phase) vs station distance at selected frequencies.

    Each frequency is a separate line coloured by
    :attr:`PYCSAMT_STYLE.multiline`, making it easy to see how the lateral
    anomaly signature changes with depth (period).

    Parameters
    ----------
    response : ForwardResponse2D
    mode : {'te', 'tm'}
    quantity : {'rho_a', 'phase'}
    freq_indices : sequence of int or None
        Indices into ``response.freqs`` to display.  When ``None``,
        *n_freqs_shown* equally spaced indices are chosen automatically.
    n_freqs_shown : int
        Number of frequency curves when *freq_indices* is ``None``.
    lw : float
        Line width.  Default: ``PYCSAMT_STYLE.multiline.lw``.
    alpha : float
        Line alpha.  Default: ``PYCSAMT_STYLE.multiline.alpha``.
    title : str
    figsize : (float, float)
    ax : Axes or None

    Returns
    -------
    ax : Axes
    """
    _ml = PYCSAMT_STYLE.multiline
    if lw    is _UNSET: lw    = _ml.lw
    if alpha is _UNSET: alpha = _ml.alpha

    attr   = f"rho_a_{mode.lower()}" if quantity == "rho_a" else f"phase_{mode.lower()}"
    data   = getattr(response, attr)   # (n_freqs, n_stations)
    freqs  = response.freqs
    st_x   = response.stations_x

    nf = len(freqs)
    if freq_indices is None:
        freq_indices = np.round(
            np.linspace(0, nf - 1, min(n_freqs_shown, nf))
        ).astype(int)
    freq_indices = list(freq_indices)

    colors = _ml.colors(len(freq_indices))

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    for ki, fi in enumerate(freq_indices):
        row = data[fi, :]
        if quantity == "rho_a":
            y = _rho_vals(row)
        else:
            y = _phase_vals(row)
        per  = 1.0 / freqs[fi]
        lab  = f"T = {per:.3g} s"
        ax.plot(st_x, y, color=colors[ki], lw=lw, alpha=alpha,
                marker=".", ms=4, label=lab)

    ax.set_xlabel("Station distance  (m)", fontsize=9)
    ylabel = (_rho_label() if quantity == "rho_a" else _phase_label())
    ax.set_ylabel(f"{ylabel}  [{mode.upper()}]", fontsize=9)
    ax.legend(fontsize=7.5, framealpha=0.8, ncol=2)

    ax.set_title(
        title or f"Lateral profiles — {quantity} ({mode.upper()})",
        fontsize=10, pad=6,
    )
    _spine_style(ax)
    return ax
