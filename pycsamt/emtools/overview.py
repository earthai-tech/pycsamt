# pycsamt/emtools/overview.py
"""Combined single-station MT/CSAMT response overview.

:func:`plot_response_overview` reproduces the classic multi-panel MT
"full response" figure (as produced by tools such as MTpy) -- apparent
resistivity and phase for the off-diagonal and diagonal impedance
components, induction arrows, and phase-tensor ellipses, all sharing one
period axis -- while wiring the whole thing through pyCSAMT's own style
and view-control system (:data:`~pycsamt.api.style.PYCSAMT_STYLE`,
:data:`~pycsamt.api.control.PYCSAMT_CONTROL`) so every colour, error bar,
and axis convention matches the rest of ``emtools``.
"""

from __future__ import annotations

import copy
from dataclasses import replace
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from matplotlib.transforms import Affine2D

from ..api._rose_style import _UNSET
from ..api.style import PYCSAMT_STYLE
from ._core import _axes_list, _get_t_block, ensure_sites
from .plot import (
    _component_style,
    _comp_slice,
    _control_or_default,
    _err_rhoa,
    _errorbar_from_style,
    _phase_display,
    _rho_display,
    _rhoa_from,
    _x_from_freq,
    _zblk_flex,
)
from .tensor import _attach_cbar, _resolve_cvals, build_phase_tensor_table
from .tf import _pick_station

__all__ = ["plot_response_overview"]


# ------------------------------- helpers --------------------------------- #


def _x_offset(x0: float, delta: float, log_x: bool) -> float:
    """Shift *x0* by *delta* "decades" on a shared period/frequency axis.

    When the shared axis is drawn with a log scale (``log_x=True``, the
    ``"period"``/``"frequency"`` views), *delta* is applied
    multiplicatively (``x0 * 10**delta``) so a fixed *delta* always reads
    as the same visual fraction of a decade regardless of *x0*.  When the
    axis already carries ``log10(...)`` values on a linear scale
    (``"log10_period"``/``"log10_frequency"``), *delta* is simply added.
    """
    return x0 * (10.0**delta) if log_x else x0 + delta


def _legend_handles(components, style_source, colors):
    handles, labels = [], []
    for comp in components:
        style = style_source(comp)
        color = colors.get(comp) if colors else None
        kw = style.plot_kwargs()
        if color is not None:
            kw["color"] = color
        kw.pop("label", None)
        handles.append(Line2D([], [], **kw))
        labels.append(f"$Z_{{{comp}}}$")
    return handles, labels


def _draw_component_column(
    ax_rho,
    ax_phase,
    *,
    z,
    ze,
    fr,
    x,
    components,
    ctl,
    phase_range,
    colors,
    show_error_bars,
    show_phase_error_bars,
    raw,
    force_style,
    log_log_rho,
):
    """Draw rho/phase curves for one impedance-component pair.

    When *log_log_rho* is true (the default, matching ``ctl.rho.view ==
    "log10"``), apparent resistivity is drawn from *raw* Ohm-metre values
    on a true Matplotlib log y-axis (``ax_rho.set_yscale("log")``) instead
    of pre-transforming values to :math:`\\log_{10}\\rho_a` on a linear
    axis. Both display the same information, but only the true log axis
    draws Matplotlib's own per-decade minor-tick grid -- the classic
    "log-log" look -- rather than sparse integer gridlines over already
    logged numbers.
    """
    rho_label = phase_label = ""
    if log_log_rho:
        ax_rho.set_yscale("log")
    for comp in components:
        zz = _comp_slice(z, comp)
        ee = None
        if isinstance(ze, np.ndarray) and ze.shape == z.shape:
            ee = _comp_slice(ze, comp)
        if log_log_rho:
            rho_y = _rhoa_from(zz, fr)
            rho_err = _err_rhoa(zz, ee, fr)
            rho_label = r"$\rho_a$ ($\Omega\,\mathrm{m}$)"
        else:
            rho_y, rho_err, rho_label = _rho_display(zz, ee, fr, ctl)
        phase_y, phase_err, phase_label, display_range = _phase_display(
            zz, ee, ctl, phase_range
        )
        style = _component_style(comp, raw, force_style)
        color = colors.get(comp) if colors else None
        _errorbar_from_style(
            ax_rho,
            x,
            rho_y,
            rho_err,
            style,
            color=color,
            show_error_bars=show_error_bars,
        )
        _errorbar_from_style(
            ax_phase,
            x,
            phase_y,
            phase_err,
            style,
            color=color,
            show_error_bars=show_phase_error_bars,
        )
    if display_range is not None:
        ax_phase.set_ylim(*display_range)
    return rho_label, phase_label


def _draw_induction_arrow_row(
    ax,
    x,
    u,
    v,
    *,
    color,
    tilt_decades,
    dy_scale,
    lw,
    mutation_scale,
    log_x,
    zorder=3,
):
    """Draw one real- or imaginary-part induction-arrow strip.

    Each arrow tail sits exactly at its sample's position on the shared
    period/frequency axis. Its 2-D direction encodes the true
    ``(u, v) = (Tx, Ty)`` vector: the vertical head offset is the
    physically meaningful part (``dy_scale * |T| * sin(angle)``, plotted
    in real data units on the row's own linear y-axis); the horizontal
    tilt (``tilt_decades * |T| * cos(angle)``) is a schematic fan
    reserved for visual separation between neighbouring periods and does
    not represent a period shift.
    """
    tips = []
    for xi, ui, vi in zip(x, u, v):
        if not np.all(np.isfinite([xi, ui, vi])):
            continue
        mag = float(np.hypot(ui, vi))
        if mag < 1e-9:
            continue
        ang = float(np.arctan2(vi, ui))
        x1 = _x_offset(float(xi), tilt_decades * mag * np.cos(ang), log_x)
        y1 = dy_scale * mag * np.sin(ang)
        tips.append(y1)
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(float(xi), 0.0),
            arrowprops=dict(
                arrowstyle="-|>",
                color=color,
                lw=lw,
                mutation_scale=mutation_scale,
                shrinkA=0,
                shrinkB=0,
            ),
            zorder=zorder,
        )
    return tips


def _draw_pt_ellipse_row(
    ax,
    x,
    s1,
    s2,
    theta,
    cvals,
    skew,
    *,
    cmap,
    norm,
    scale,
    min_aspect,
    edgecolor,
    linewidth,
    alpha,
    skew_threshold,
    mark_3d,
    cells_per_decade,
    log_x,
):
    """Draw a phase-tensor ellipse strip anchored on a shared period axis.

    Ellipse geometry is built in *display pixels* (the same technique
    used by ``normalise_by="shape"`` in
    :func:`pycsamt.emtools.tensor.plot_phase_tensor_psection`) rather
    than in axis data units.  This is required here, not merely
    cosmetic: this row's x-axis is the same log-scale period axis shared
    with the rho/phase panels above it, and a data-unit ellipse width
    would be nonlinearly warped by that log transform (a constant "width
    in seconds" balloons at long period and vanishes at short period).
    Pixel sizing keeps every ellipse a true, undistorted shape wherever
    it falls on the shared axis.
    """
    ax.figure.canvas.draw()
    finite_x = x[np.isfinite(x)]
    x0 = float(np.nanmedian(finite_x)) if finite_x.size else 0.0
    x1 = _x_offset(x0, 1.0 / max(cells_per_decade, 1e-6), log_x)
    px0 = ax.transData.transform((x0, 0.0))[0]
    px1 = ax.transData.transform((x1, 0.0))[0]
    cell_px = max(abs(px1 - px0), 1e-6)
    width_px = scale * cell_px

    bbox = ax.bbox
    for xi, s1i, s2i, thi, ci, beti in zip(x, s1, s2, theta, cvals, skew):
        if not np.all(np.isfinite([xi, s1i, s2i, thi, ci])):
            continue
        ratio = abs(s2i) / max(abs(s1i), 1e-12)
        ratio = min(1.0, max(float(min_aspect), ratio))
        disp = ax.transData.transform((xi, 0.0))
        x_frac, y_frac = ax.transAxes.inverted().transform(disp)
        is_3d = (
            skew_threshold is not None
            and mark_3d
            and np.isfinite(beti)
            and abs(beti) > skew_threshold
        )
        lw_i = linewidth * 3.0 if is_3d else linewidth
        shape_transform = (
            Affine2D()
            .scale(0.5 * width_px, 0.5 * width_px * ratio)
            .rotate_deg(thi)
            .scale(1.0 / bbox.width, 1.0 / bbox.height)
            .translate(x_frac, y_frac)
            + ax.transAxes
        )
        ell = Ellipse(
            (0.0, 0.0),
            width=2.0,
            height=2.0,
            angle=0.0,
            transform=shape_transform,
            facecolor=cmap(norm(ci)),
            edgecolor=edgecolor,
            linewidth=lw_i,
            alpha=alpha,
            zorder=3,
        )
        ell.set_clip_path(ax.patch)
        ax.add_patch(ell)


# ------------------------------ main plot --------------------------------- #


def plot_response_overview(
    sites: Any,
    *,
    station: str | None = None,
    control: Any | None = None,
    x_view: str | None = "period",
    log_log_rho=_UNSET,  # default: True when the resolved control.rho.view == "log10"
    # ── which components ──────────────────────────────────────────────────
    offdiag_components: tuple[str, str] = ("xy", "yx"),
    diag_components: tuple[str, str] = ("xx", "yy"),
    show_diag: bool = True,
    phase_range: tuple[float, float] | None = None,
    # ── layout ────────────────────────────────────────────────────────────
    height_ratios: tuple[float, float, float, float] = (2.0, 1.0, 1.4, 1.7),
    cbar_orientation: str = "horizontal",
    cbar_height_ratio: float = 0.12,
    cbar_pad_ratio: float = 0.55,
    cbar_width_ratio: float = 0.10,
    figsize: tuple[float, float] = (11.0, 9.8),
    wspace: float = 0.24,
    hspace: float = 0.10,
    axes=None,
    # ── impedance style ──────────────────────────────────────────────────
    colors: dict[str, str] | None = None,
    raw: bool = False,
    force_style: bool = False,
    show_error_bars: bool = True,
    show_phase_error_bars: bool = False,
    show_component_legend: bool = True,
    title: str | None = None,
    # ── induction-arrow row ──────────────────────────────────────────────
    show_arrows: bool = True,
    arrow_colors: tuple[str | None, str | None] = (None, None),
    arrow_tilt_decades: float = 0.4,
    arrow_dy_scale: float = 1.0,
    arrow_lw: float = 1.3,
    arrow_mutation_scale: float = 7.0,
    ylim_arrows: tuple[float, float] | None = None,
    show_arrow_legend: bool = True,
    # ── phase-tensor ellipse row ─────────────────────────────────────────
    show_ellipses: bool = True,
    c_by=_UNSET,  # default: PYCSAMT_STYLE.pt_ellipse.c_by
    cmap=_UNSET,  # default: auto from c_by via resolve_cmap()
    clim: tuple[float, float] | None = None,
    clim_pct=_UNSET,  # default: PYCSAMT_STYLE.pt_ellipse.clim_pct
    symmetric_clim=_UNSET,  # default: PYCSAMT_STYLE.pt_ellipse.resolve_symmetric_clim()
    ellipse_scale: float = 1.5,  # bigger than PYCSAMT_STYLE.pt_ellipse.scale
    # (0.85): this row's pixel/"shape"-based sizing (see
    # _draw_pt_ellipse_row) does not tile a fixed station grid the way
    # plot_phase_tensor_psection does, so headroom for a larger scale
    # exists without ellipses colliding into their neighbours' data.
    min_aspect=_UNSET,  # default: PYCSAMT_STYLE.pt_ellipse.min_aspect
    cells_per_decade: float = 6.0,
    edgecolor=_UNSET,  # default: PYCSAMT_STYLE.pt_ellipse.edgecolor
    linewidth=_UNSET,  # default: PYCSAMT_STYLE.pt_ellipse.linewidth
    ellipse_alpha=_UNSET,  # default: PYCSAMT_STYLE.pt_ellipse.alpha
    skew_threshold=_UNSET,  # default: PYCSAMT_STYLE.pt_ellipse.skew_threshold
    mark_3d=_UNSET,  # default: PYCSAMT_STYLE.pt_ellipse.mark_3d
    show_ellipse_colorbar: bool = True,
    ellipse_colorbar_label: str | None = None,
    # ── ticks / grid ─────────────────────────────────────────────────────
    tick_fontsize: int = 8,
    grid: bool = True,
    # ── standard emtools ─────────────────────────────────────────────────
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    r"""Single-station MT/CSAMT "full response" overview figure.

    Reproduces the classic multi-panel layout used to QC a wideband MT
    sounding -- apparent resistivity and phase for the off-diagonal
    (:math:`Z_{xy}, Z_{yx}`) and diagonal (:math:`Z_{xx}, Z_{yy}`)
    impedance components side by side, with induction arrows and
    phase-tensor ellipses spanning both columns underneath, all sharing
    one period/frequency axis -- while going entirely through pyCSAMT's
    own :data:`~pycsamt.api.style.PYCSAMT_STYLE` and
    :data:`~pycsamt.api.control.PYCSAMT_CONTROL` systems instead of
    hardcoded colours or axis conventions.

    Layout (rows, top to bottom; all four data rows share the x-axis)::

        ┌─────────────────┬─────────────────┐
        │  App. Res. xy/yx │ App. Res. xx/yy │        height_ratios[0]
        ├─────────────────┼─────────────────┤
        │  Phase xy/yx     │ Phase xx/yy     │        height_ratios[1]
        ├─────────────────┴─────────────────┤
        │      Induction arrows (real/imag)  │        height_ratios[2]
        ├─────────────────────────────────────┤
        │      Phase-tensor ellipse strip    │        height_ratios[3]
        ├─────────────────────────────────────┤
        │           skew β (°) colourbar     │        cbar_height_ratio
        └─────────────────────────────────────┘

    The rho/phase/arrow/ellipse rows all span the full data width; the
    ellipse row's colourbar lives in its own reserved row below it by
    default (see *cbar_orientation*) rather than a right-hand gutter, so
    it never narrows the ellipse row relative to the arrow row above it
    and the right margin stays free.

    Every row shares one x-axis. By default (*x_view* = ``"period"``) that
    axis is a true Matplotlib log scale over raw period in seconds --
    the classic MT "log-log" quicklook, complete with Matplotlib's own
    per-decade minor-tick grid. Apparent resistivity likewise defaults to
    a true log y-axis over raw :math:`\Omega\cdot\mathrm{m}` values
    (*log_log_rho*), so the top row is genuinely log-log and the phase
    row below it (linear degrees over the same log x-axis) is genuinely
    semilog -- both draw the dense reference-figure-style grid rather
    than sparse integer gridlines over pre-logged numbers. Pass
    ``x_view="log10_period"`` for the alternative
    :math:`\log_{10}T\,(\mathrm{s})` linear-axis convention used
    elsewhere in ``emtools`` (e.g. :func:`plot_raw_sites_1d`), or
    ``x_view=None`` to defer entirely to ``control.x.view``. Apparent
    resistivity similarly can be forced to the pre-logged
    :math:`\log_{10}\rho_a` linear-axis convention with
    ``log_log_rho=False`` (then following ``control.rho.view``), and
    phase follows ``control.phase`` throughout.

    Parameters
    ----------
    sites : Sites-like
        EDI path, glob pattern, :class:`~pycsamt.site.base.Sites`, or any
        input accepted by :func:`pycsamt.emtools.ensure_sites`.
    station : str or None
        Station to plot. Defaults to the first station (sorted by name)
        when *sites* resolves to more than one.
    control : object, optional
        Plot view control. Defaults to
        :data:`pycsamt.api.control.PYCSAMT_CONTROL`. Only ``control.phase``
        (and, when *log_log_rho* is ``False``, ``control.rho``) is read
        directly from this object -- the x-axis view is controlled
        separately by *x_view* so this function's default log-log look
        does not depend on, or silently change, the shared global
        control's ``x.view`` setting.
    x_view : str or None, default "period"
        X-axis convention for all four rows, applied as a local override
        on top of *control* (the shared control object is never mutated).
        One of ``"period"``, ``"log10_period"``, ``"frequency"``,
        ``"log10_frequency"`` (see
        :class:`pycsamt.api.control.FrequencyAxisControl`), or ``None`` to
        use whatever ``control.x.view`` is already set to.
    log_log_rho : bool, optional
        Force the apparent-resistivity row onto a true log y-axis over
        raw values (``True``) or the pre-logged
        :math:`\log_{10}\rho_a`-on-linear-axis convention (``False``).
        Defaults to ``True`` when the resolved ``control.rho.view`` is
        ``"log10"`` (the package default) and ``False`` otherwise, so
        passing ``rho__view="linear"`` through *control* still works as
        expected without needing to also touch this flag.
    offdiag_components, diag_components : (str, str)
        Impedance components drawn (overlaid on the same axes pair) in
        the left and right column respectively. Any two-letter component
        keys accepted by :attr:`PYCSAMT_STYLE.mt` work here, so a TE/TM
        pair (``("te", "tm")``) is equally valid if that is how a survey
        was processed.
    show_diag : bool, default True
        Draw the right-hand (*diag_components*) column. When ``False``
        the figure narrows to a single column and the arrow/ellipse rows
        span that one column's width.
    phase_range : (lo, hi) or None, optional
        Explicit phase display range shared by both columns. If omitted,
        the active ``control.phase`` policy is used (default
        :math:`\pm 180^\circ`) -- one consistent range for both columns,
        unlike tools that let the off-diagonal and diagonal phase axes
        drift to different ranges.
    height_ratios : (rho, phase, arrows, ellipses), default (2.0, 1.0, 1.4, 1.7)
        Relative row heights. Apparent resistivity gets twice the phase
        row's height by default (a 2:1 ratio, i.e. rho is 2/3 and phase
        1/3 of the combined rho+phase height) -- the conventional "big
        rho, small phase" look; pass e.g. ``(1.0, 1.0, ...)`` for equal
        rows instead.
    cbar_orientation : {"horizontal", "vertical"}, default "horizontal"
        Placement of the ellipse-row colourbar. ``"horizontal"`` (the
        default) adds it as its own thin row *below* the ellipse strip,
        spanning the full data width and freeing the right-hand margin
        entirely. ``"vertical"`` instead reserves a narrow gutter column
        to the right of every row (see *cbar_width_ratio*).  Either way
        the space is reserved in the gridspec *up front* rather than
        carved out of the ellipse axes after the fact -- the latter is
        what :func:`matplotlib.figure.Figure.colorbar`-via-divider
        approaches normally do, and it silently shrinks the ellipse row
        relative to the arrow row above it even though both still report
        "the same" x-limits, so a given period ends up at a different
        pixel column in each row. Reserving the space in advance keeps
        the arrow and ellipse rows pixel-aligned.
    cbar_height_ratio : float, default 0.12
        Height of the horizontal colourbar row (``cbar_orientation=
        "horizontal"``), relative to the same units as *height_ratios*
        (e.g. relative to the ellipse row's own ``height_ratios[3]``).
    cbar_pad_ratio : float, default 0.55
        Height of a blank spacer row inserted between the ellipse row and
        the horizontal colourbar row, in the same units as
        *height_ratios*. Needed because the ellipse row's own x tick
        labels and "Period (s)" label are drawn outside its axes box, in
        the margin below it -- with too little pad they collide with the
        colourbar's own tick labels. Increase if your *tick_fontsize* or
        a custom x label still overlaps the colourbar.
    cbar_width_ratio : float, default 0.10
        Width of the vertical colourbar column (``cbar_orientation=
        "vertical"``), as a fraction of one data column.
        Both *cbar_height_ratio* and *cbar_width_ratio* -- like the rest
        of the reserved-space guarantee described under
        *cbar_orientation* -- are ignored when *axes* is provided (the
        colourbar then falls back to narrowing the ellipse axes itself,
        since new gridspec rows/columns cannot be inserted into an
        already-built external layout).
    figsize : (float, float), default (11.0, 9.8)
        Figure size (ignored when *axes* is provided).
    wspace, hspace : float
        Column/row spacing (ignored when *axes* is provided).
    axes : sequence of :class:`~matplotlib.axes.Axes` or None
        Pre-built axes to draw into instead of creating a new figure and
        gridspec. Must supply, in order, the rho and phase axes for each
        column (``2 * ncols`` axes, ``ncols=2`` unless *show_diag* is
        ``False``), followed by the arrow-row axes and the ellipse-row
        axes -- e.g. for the default two-column layout: ``[ax_rho_off,
        ax_rho_diag, ax_phase_off, ax_phase_diag, ax_arrow, ax_ellipse]``.
        Useful for composing this overview into a larger multi-station
        figure built with your own gridspec. Note the arrow/ellipse
        pixel-alignment guarantee described under *cbar_width_ratio*
        does not apply here -- align them yourself if you supply axes.
    colors : dict, optional
        Optional per-component colour overrides (keys are component
        letters, e.g. ``{"xy": "black"}``). Omitted components keep
        their :data:`PYCSAMT_STYLE.mt` colour -- the package default is
        already the reference-figure convention (blue circles for
        :math:`Z_{xy}`, red squares for :math:`Z_{yx}`).
    raw, force_style : bool
        When ``raw=True``, curves use :data:`PYCSAMT_STYLE.raw` (a
        neutral diagnostic style) instead of per-component colours,
        unless *force_style* is also ``True``. See
        :func:`plot_raw_sites_1d` for the same convention.
    show_error_bars, show_phase_error_bars : bool, default True, False
        Toggle apparent-resistivity and phase error bars independently.
        Both use the *same colour* as their component's curve (through
        ``style.errorbar_kwargs()``), which is the "error bars based on
        that colour" behaviour.
    show_component_legend : bool, default True
        Draw a small legend (component colour + marker) inside each
        top-row (apparent-resistivity) panel.
    title : str or None
        Figure title. Defaults to the station name.
    show_arrows : bool, default True
        Draw the induction-arrow row. Skipped (with a "no tipper"
        placeholder) when the station has no tipper data.
    arrow_colors : (real, imag), default (None, None)
        Explicit colour override for the real- and imaginary-part
        arrows. ``None`` falls back to :data:`PYCSAMT_STYLE.mt.xy` and
        :data:`PYCSAMT_STYLE.mt.yx` respectively -- the same real/imag
        colour convention already used by
        :func:`pycsamt.emtools.plot.plot_response_tipper`.
    arrow_tilt_decades : float, default 0.4
        Schematic horizontal fan applied to each arrow (see
        :func:`_draw_induction_arrow_row`); purely cosmetic separation,
        not a period shift.
    arrow_dy_scale : float, default 1.0
        Scales the arrow's vertical (physically meaningful) extent.
        Tipper magnitudes are typically O(0.01-1), so the default already
        matches a sensible row height; increase for a very "quiet"
        station or decrease if arrows overrun neighbouring rows.
    arrow_lw, arrow_mutation_scale : float
        Arrow line width and matplotlib ``mutation_scale`` (head size).
    ylim_arrows : (lo, hi) or None
        Explicit y-limits for the arrow row. Auto-scaled from the drawn
        arrow tips (with margin) when omitted.
    show_arrow_legend : bool, default True
        Draw a "real"/"imag" legend below the arrow row.
    show_ellipses : bool, default True
        Draw the phase-tensor ellipse row. Skipped (with a placeholder)
        when phase-tensor invariants cannot be computed for the station.
    c_by, cmap, clim, clim_pct, symmetric_clim, ellipse_scale, min_aspect,
    edgecolor, linewidth, ellipse_alpha, skew_threshold, mark_3d :
        Phase-tensor ellipse controls -- identical semantics to
        :func:`pycsamt.emtools.tensor.plot_phase_tensor_strip`, defaulting
        to :data:`PYCSAMT_STYLE.pt_ellipse`. To reproduce the classic
        MTpy "phimin" colouring (blue-to-orange-to-navy over
        :math:`0^\circ`-:math:`90^\circ`), pass
        ``c_by="phimin_deg", clim=(0.0, 90.0)`` together with a diverging
        colormap of your choice -- note this differs from ``c_by="phi_min"``,
        which colours by the raw (tan-units) phase-tensor singular value
        instead of its arctan in degrees and saturates near one end of a
        0-90 scale; the package default (``c_by="skew"``) is the more
        diagnostic choice for flagging 3-D structure.
    cells_per_decade : float, default 8.0
        Visual ellipse pitch along the period axis (ellipse-widths per
        decade); see :func:`pycsamt.emtools.tensor.plot_phase_tensor_strip`
        for why this is independent of the actual sample spacing.
    show_ellipse_colorbar : bool, default True
        Attach a colourbar to the ellipse row.
    ellipse_colorbar_label : str or None
        Override the automatic colourbar label derived from *c_by*.
    tick_fontsize : int, default 8
        Tick-label size shared by all rows.
    grid : bool, default True
        Draw light panel grids on the rho/phase/arrow rows.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`pycsamt.emtools.ensure_sites`.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> from pycsamt.emtools import plot_response_overview
    >>> fig = plot_response_overview("data/gv_data/gv_final_edi", station="gv100")

    Raw resistivity (not log10) and a "phimin"-style ellipse colouring:

    >>> from pycsamt.api.control import PYCSAMT_CONTROL
    >>> with PYCSAMT_CONTROL.context(rho__view="linear"):
    ...     fig = plot_response_overview(
    ...         "data/gv_data/gv_final_edi", station="gv100",
    ...         c_by="phimin_deg", clim=(0.0, 90.0), cmap="turbo",
    ...     )

    See Also
    --------
    plot_response_tipper : Per-component grid across many stations.
    plot_phase_tensor_strip : Standalone ellipse strip for one station.
    plot_induction_arrows : Map-view induction arrows across a survey.
    """
    ctl = _control_or_default(control)
    if x_view is not None and str(x_view).lower() != ctl.x.view.lower():
        ctl = copy.copy(ctl)
        ctl.x = replace(ctl.x, view=str(x_view))
    if log_log_rho is _UNSET:
        log_log_rho = ctl.rho.view.lower() == "log10"

    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup, strict=strict, verbose=verbose
    )
    station, ed = _pick_station(S, station)

    ncols = 2 if show_diag else 1
    columns = [offdiag_components] + ([diag_components] if show_diag else [])

    cbar_orientation = str(cbar_orientation).lower()
    if cbar_orientation not in {"horizontal", "vertical"}:
        raise ValueError("cbar_orientation must be 'horizontal' or 'vertical'.")
    cbar_horizontal = cbar_orientation == "horizontal"

    axes_given = _axes_list(axes, 2 * ncols + 2) if axes is not None else None
    ax_cbar = None
    if axes_given is None:
        fig = plt.figure(figsize=figsize, constrained_layout=False)
        # A dedicated colourbar row/column is reserved up front rather
        # than carved out of ax_ell after the fact -- see
        # *cbar_orientation* in the docstring for why: the latter
        # silently shrinks the ellipse row relative to the arrow row
        # above it even though both still share the same x-limits.
        if cbar_horizontal:
            # A blank spacer row sits between the ellipse row and the
            # colourbar row. `hspace` alone is not enough here: it pads
            # between *axes boxes*, but the ellipse row's own x tick
            # labels and "Period (s)" label are drawn *outside* that box,
            # in the margin -- with no spacer they collide with the
            # colourbar's tick labels directly below.
            outer = fig.add_gridspec(
                6,
                ncols,
                height_ratios=[*height_ratios, cbar_pad_ratio, cbar_height_ratio],
                wspace=wspace,
                hspace=hspace,
            )
        else:
            outer = fig.add_gridspec(
                4,
                ncols + 1,
                width_ratios=[*([1.0] * ncols), cbar_width_ratio],
                height_ratios=height_ratios,
                wspace=wspace,
                hspace=hspace,
            )

        ax_rho = [fig.add_subplot(outer[0, j]) for j in range(ncols)]
        ax_phase = [
            fig.add_subplot(outer[1, j], sharex=ax_rho[j]) for j in range(ncols)
        ]
        for j in range(1, ncols):
            ax_rho[j].sharex(ax_rho[0])
        ax_arrow = fig.add_subplot(outer[2, :ncols], sharex=ax_rho[0])
        ax_ell = fig.add_subplot(outer[3, :ncols], sharex=ax_rho[0])
        ax_cbar = fig.add_subplot(outer[5, :] if cbar_horizontal else outer[3, ncols])
        ax_cbar.axis("off")
    else:
        ax_rho = axes_given[:ncols]
        ax_phase = axes_given[ncols : 2 * ncols]
        ax_arrow = axes_given[2 * ncols]
        ax_ell = axes_given[2 * ncols + 1]
        fig = ax_rho[0].figure

    fig.suptitle(title or str(station), fontsize=11, fontweight="bold")

    out = _zblk_flex(ed)
    if len(out) == 4:
        _, z, fr, ze = out
    else:
        _, z, fr = out[:3]
        ze = None

    if z is None or fr is None:
        for ax in ax_rho + ax_phase + [ax_arrow, ax_ell]:
            ax.text(
                0.5, 0.5, "no impedance data", ha="center", va="center",
                transform=ax.transAxes,
            )
        return fig

    fr = np.asarray(fr, dtype=float)
    x, x_label, log_x = _x_from_freq(fr, ctl)

    for j, comps in enumerate(columns):
        a_rho, a_phase = ax_rho[j], ax_phase[j]
        if log_x:
            a_rho.set_xscale("log")
            a_phase.set_xscale("log")
        if grid:
            a_rho.grid(True, alpha=0.25, which="both")
            a_phase.grid(True, alpha=0.25, which="both")
        rho_label, phase_label = _draw_component_column(
            a_rho,
            a_phase,
            z=z,
            ze=ze,
            fr=fr,
            x=x,
            components=comps,
            ctl=ctl,
            phase_range=phase_range,
            colors=colors,
            show_error_bars=show_error_bars,
            show_phase_error_bars=show_phase_error_bars,
            raw=raw,
            force_style=force_style,
            log_log_rho=log_log_rho,
        )
        a_rho.set_ylabel(rho_label if j == 0 else "")
        a_phase.set_ylabel(phase_label if j == 0 else "")
        a_rho.tick_params(axis="x", labelbottom=False)
        a_phase.tick_params(axis="x", labelbottom=False)
        a_rho.tick_params(axis="both", labelsize=tick_fontsize)
        a_phase.tick_params(axis="both", labelsize=tick_fontsize)
        if show_component_legend:
            handles, labels = _legend_handles(
                comps, lambda c: _component_style(c, raw, force_style), colors
            )
            a_rho.legend(
                handles, labels, loc="best", fontsize=8, framealpha=0.85
            )

    # ── induction-arrow row ───────────────────────────────────────────────
    if show_arrows:
        t_out = _get_t_block(ed, with_errors=False)
        _, tip, tfr = t_out[:3]
        if tip is None or tfr is None:
            ax_arrow.text(
                0.5, 0.5, "no tipper", ha="center", va="center",
                transform=ax_arrow.transAxes, color="0.4",
            )
        else:
            tfr = np.asarray(tfr, dtype=float)
            tx_x, _, _ = _x_from_freq(tfr, ctl)
            tip = np.asarray(tip, dtype=complex)
            mt = PYCSAMT_STYLE.mt
            real_color = arrow_colors[0] or mt.xy.color
            imag_color = arrow_colors[1] or mt.yx.color
            tips_r = _draw_induction_arrow_row(
                ax_arrow,
                tx_x,
                np.real(tip[:, 0]),
                np.real(tip[:, 1]),
                color=real_color,
                tilt_decades=arrow_tilt_decades,
                dy_scale=arrow_dy_scale,
                lw=arrow_lw,
                mutation_scale=arrow_mutation_scale,
                log_x=log_x,
            )
            tips_i = _draw_induction_arrow_row(
                ax_arrow,
                tx_x,
                np.imag(tip[:, 0]),
                np.imag(tip[:, 1]),
                color=imag_color,
                tilt_decades=arrow_tilt_decades,
                dy_scale=arrow_dy_scale,
                lw=arrow_lw,
                mutation_scale=arrow_mutation_scale,
                log_x=log_x,
            )
            if ylim_arrows is not None:
                ax_arrow.set_ylim(*ylim_arrows)
            else:
                all_tips = np.asarray(tips_r + tips_i, dtype=float)
                lim = (
                    float(np.nanmax(np.abs(all_tips))) * 1.25
                    if all_tips.size
                    else 0.1
                )
                lim = max(lim, 1e-3)
                ax_arrow.set_ylim(-lim, lim)
            ax_arrow.axhline(0.0, color="0.75", lw=0.7, zorder=0)
            if show_arrow_legend:
                handles = [
                    Line2D([], [], color=real_color, lw=arrow_lw, label="real"),
                    Line2D([], [], color=imag_color, lw=arrow_lw, label="imag"),
                ]
                ax_arrow.legend(
                    handles=handles, loc="upper right", fontsize=8,
                    framealpha=0.85,
                )
    else:
        ax_arrow.axis("off")
    if grid and show_arrows:
        ax_arrow.grid(True, alpha=0.25, which="both")
    ax_arrow.set_ylabel("Induction\narrows", fontsize=9)
    ax_arrow.tick_params(axis="both", labelsize=tick_fontsize)
    # The ellipse row directly below shares this same x-axis (sharex=
    # ax_rho[0]) and already draws the period/frequency tick labels, so
    # repeating them here would be redundant -- only tick *marks* stay,
    # not the numeric labels.
    ax_arrow.tick_params(axis="x", labelbottom=False)

    # ── phase-tensor ellipse row ─────────────────────────────────────────
    if show_ellipses:
        # `ed` is a Site-wrapped object here (from `_pick_station`); passing
        # it straight through in a list makes `ensure_sites`/`to_sites`
        # silently resolve to zero items (it expects raw EDI-like items, not
        # an already-Site-wrapped one). Unwrap to the underlying EDI first,
        # mirroring `_core._wrap_one`'s convention for the same situation.
        edi = getattr(ed, "edi", None)
        pt_item = edi if edi is not None else ed
        df = build_phase_tensor_table(
            [pt_item],
            recursive=False,
            on_dup=on_dup,
            strict=False,
            verbose=verbose,
        )
        if df.empty:
            ax_ell.text(
                0.5, 0.5, "no phase tensor data", ha="center", va="center",
                transform=ax_ell.transAxes, color="0.4",
            )
        else:
            _es = PYCSAMT_STYLE.pt_ellipse
            c_by_ = _es.c_by if c_by is _UNSET else c_by
            cmap_name = (
                _es.copy(c_by=c_by_).resolve_cmap() if cmap is _UNSET else cmap
            )
            clim_pct_ = _es.clim_pct if clim_pct is _UNSET else clim_pct
            symmetric_ = (
                _es.copy(c_by=c_by_).resolve_symmetric_clim()
                if symmetric_clim is _UNSET
                else symmetric_clim
            )
            scale_ = _es.scale if ellipse_scale is _UNSET else ellipse_scale
            min_aspect_ = _es.min_aspect if min_aspect is _UNSET else min_aspect
            edgecolor_ = _es.edgecolor if edgecolor is _UNSET else edgecolor
            linewidth_ = _es.linewidth if linewidth is _UNSET else linewidth
            alpha_ = _es.alpha if ellipse_alpha is _UNSET else ellipse_alpha
            skew_threshold_ = (
                _es.skew_threshold if skew_threshold is _UNSET else skew_threshold
            )
            mark_3d_ = _es.mark_3d if mark_3d is _UNSET else mark_3d

            df = df.sort_values("freq")
            x_ell = ctl.x.transform(df["freq"].to_numpy(float))
            cvals, cbar_label = _resolve_cvals(df, c_by_)
            if ellipse_colorbar_label:
                cbar_label = ellipse_colorbar_label
            finite_c = cvals[np.isfinite(cvals)]
            if clim is not None:
                vmin, vmax = float(clim[0]), float(clim[1])
            else:
                lo_pct, hi_pct = clim_pct_
                vmin = (
                    float(np.nanpercentile(finite_c, lo_pct))
                    if len(finite_c)
                    else -1.0
                )
                vmax = (
                    float(np.nanpercentile(finite_c, hi_pct))
                    if len(finite_c)
                    else 1.0
                )
                if symmetric_:
                    vlim = max(abs(vmin), abs(vmax))
                    vmin, vmax = -vlim, vlim
            if not np.isfinite(vmin) or not np.isfinite(vmax):
                vmin, vmax = -1.0, 1.0
            elif np.isclose(vmin, vmax):
                pad = max(abs(vmin) * 0.05, 1e-6)
                vmin, vmax = vmin - pad, vmax + pad
            norm = Normalize(vmin=vmin, vmax=vmax)
            cm = plt.get_cmap(cmap_name)

            ax_ell.set_xlim(ax_rho[0].get_xlim())
            ax_ell.set_ylim(-0.62, 0.62)
            _draw_pt_ellipse_row(
                ax_ell,
                x_ell,
                df["s1"].to_numpy(float),
                df["s2"].to_numpy(float),
                df["theta"].to_numpy(float),
                cvals,
                df["skew"].to_numpy(float),
                cmap=cm,
                norm=norm,
                scale=scale_,
                min_aspect=min_aspect_,
                edgecolor=edgecolor_,
                linewidth=linewidth_,
                alpha=alpha_,
                skew_threshold=skew_threshold_,
                mark_3d=mark_3d_,
                cells_per_decade=cells_per_decade,
                log_x=log_x,
            )
            ax_ell.set_yticks([])
            if show_ellipse_colorbar:
                sm = ScalarMappable(cmap=cm, norm=norm)
                sm.set_array([])
                if ax_cbar is not None:
                    # Dedicated gridspec row/column (see cbar_orientation):
                    # keeps ax_arrow and ax_ell pixel-aligned since neither
                    # loses space to the colourbar after the fact.
                    ax_cbar.axis("on")
                    cb = fig.colorbar(
                        sm,
                        cax=ax_cbar,
                        orientation="horizontal" if cbar_horizontal else "vertical",
                    )
                    cb.set_label(cbar_label, fontsize=9)
                    cb.ax.tick_params(labelsize=8)
                else:
                    _attach_cbar(ax_ell, sm, cbar_label)
    else:
        ax_ell.axis("off")
    ax_ell.set_ylabel("Phase\ntensor", fontsize=9)
    ax_ell.tick_params(axis="both", labelsize=tick_fontsize)
    ax_ell.tick_params(axis="x", labelbottom=True)
    ax_ell.set_xlabel(x_label, fontsize=10)
    if grid and show_ellipses:
        ax_ell.grid(True, ls=":", lw=0.4, color="#cccccc", alpha=0.7, zorder=0)

    return fig
