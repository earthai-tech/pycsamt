# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Topography rendering helpers for 2-D section and pseudosection plots.

Two rendering modes are provided:

**Depth-section mode** (``draw_topo_section``)
    For inversion / interpretation pcolormesh plots whose x-axis is
    along-profile distance (km) and y-axis is *absolute elevation* (km)
    in a terrain-following coordinate system.  The function fills the
    above-surface polygon, draws the terrain polyline, and places station
    marker pins at the real terrain elevation.

**Pseudosection strip mode** (``draw_topo_strip``)
    For period-vs-station imshow pseudosections.  A compact elevation
    profile strip is inserted *above* the main image axes so that the
    station positions are linked to real terrain.

Both functions accept an optional :class:`~pycsamt.topo.config.TopoConfig`
argument; when omitted they fall back to the global singleton
:data:`~pycsamt.topo.config.PYCSAMT_TOPO`.

Typical usage::

    from pycsamt.topo.overlay import draw_topo_section, draw_topo_strip

    # After pcolormesh in a depth-section axes:
    draw_topo_section(ax, chainage_km, elev_m, station_names)

    # Above a pseudosection image:
    ax_strip = draw_topo_strip(fig, ax_ps, chainage_km, elev_m, names)
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

__all__ = [
    "draw_topo_section",
    "draw_topo_strip",
    "add_station_labels",
]

# Catppuccin-compatible neutral palette (overridden by TopoConfig)
_DEFAULT_FILL = "#a89070"
_DEFAULT_LINE = "#6b4e2a"
_DEFAULT_PIN = "#f5c542"  # station marker colour
_DEFAULT_LABEL = "#cdd6f4"  # station label colour (light on dark)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def draw_topo_section(
    ax,
    chainage_km: np.ndarray,
    elev_m: np.ndarray,
    station_names: Sequence[str] | None = None,
    *,
    station_x_km: np.ndarray | None = None,
    cfg=None,
    dark: bool = True,
    marker_style=None,
) -> None:
    """Overlay terrain on a depth-section axes (terrain-following frame).

    Call this **after** :func:`~matplotlib.axes.Axes.pcolormesh` has
    already drawn the 2-D section.  The axes y-axis must represent
    *absolute elevation* (km a.s.l.) in a terrain-following coordinate
    frame (positive upward or positive downward — the function reads the
    current y-limits to decide direction).

    Draws:
    1. A filled polygon masking the above-surface space.
    2. A terrain surface polyline.
    3. Station marker pins at the real terrain elevation.
    4. Station name labels above the pins.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes that already has the pcolormesh drawn.
    chainage_km : array_like (n_stations,)
        Along-profile distances of stations (km).
    elev_m : array_like (n_stations,)
        Terrain elevation at each station (m a.s.l.).
    station_names : sequence of str, optional
        Station labels.  Omit for no labels.
    station_x_km : array_like (n_stations,), optional
        X positions of station markers.  Defaults to ``chainage_km``.
    cfg : TopoConfig, optional
        Configuration override.  Defaults to :data:`PYCSAMT_TOPO`.
    dark : bool
        Use dark-palette label colours when True.
    marker_style : pycsamt.api.station.StationMarkerStyle, optional
        Station-pin style override.  Defaults to
        :data:`pycsamt.api.station.PYCSAMT_STATION_RENDERING`'s
        ``inversion`` marker when omitted, so existing callers keep
        their current appearance; pass this to use a different marker
        for this call only, without touching the global rendering
        config.
    """
    cfg = _get_cfg(cfg)
    chain = np.asarray(chainage_km, dtype=float)
    elev = np.asarray(elev_m, dtype=float) / 1000.0  # → km
    sx = np.asarray(
        station_x_km if station_x_km is not None else chain, dtype=float
    )

    if len(chain) == 0:
        return

    elev_ex = elev * cfg.exaggeration

    # ── 1. Terrain fill polygon ────────────────────────────────────────────
    # Determine which direction is "above surface" from current y-limits
    ylim = ax.get_ylim()
    y_top = max(ylim)  # the upper edge of the current view

    if cfg.clip_below_surface:
        # Build a closed polygon: surface line → top-right → top-left
        x_fill = np.concatenate([[chain[0]], chain, [chain[-1]]])
        y_fill = np.concatenate([[y_top], elev_ex, [y_top]])
        ax.fill(
            x_fill,
            y_fill,
            color=cfg.fill_color,
            alpha=cfg.fill_alpha,
            zorder=4,
            linewidth=0,
        )

    # ── 2. Surface polyline ────────────────────────────────────────────────
    if cfg.show_surface_line:
        ax.plot(
            chain,
            elev_ex,
            color=cfg.line_color,
            linewidth=cfg.line_width,
            solid_capstyle="round",
            zorder=5,
        )

    # ── 3. Station pins ────────────────────────────────────────────────────
    if cfg.station_pins_at_surface:
        from pycsamt.topo.drape import interp_elev

        pin_elev = interp_elev(chain, elev, sx) * cfg.exaggeration
    else:
        pin_elev = np.zeros_like(sx)

    # Pad markers above the topo line toward the viewer.
    y_range = abs(y_top - min(ylim))
    toward_top = 1.0 if y_top >= min(ylim) else -1.0
    pad = toward_top * cfg.marker_pad_fraction * y_range
    marker_y = pin_elev + pad

    # Use the caller's marker style, else the global rendering config.
    if marker_style is not None:
        _mstyle = marker_style
    else:
        from pycsamt.api.station import PYCSAMT_STATION_RENDERING

        _mstyle = PYCSAMT_STATION_RENDERING.inversion.marker
    ax.scatter(
        sx,
        marker_y,
        marker=_mstyle.marker,
        s=_mstyle.size,
        facecolors=_mstyle.facecolor,
        edgecolors=_mstyle.edgecolor,
        linewidths=_mstyle.linewidth,
        alpha=_mstyle.alpha,
        zorder=_mstyle.zorder,
        clip_on=False,
    )

    # ── 4. Station labels (smart thinning) ────────────────────────────────
    if station_names:
        label_color = _DEFAULT_LABEL if dark else "#4c4f69"
        # Thin labels when crowded; markers always stay visible.
        from pycsamt.api.station import StationAxisStyle

        _dummy_style = StationAxisStyle()
        figwidth = ax.figure.get_figwidth() if ax.figure is not None else 10.0
        visible_idx = _dummy_style.label_indices(station_names, figwidth)
        label_offset = toward_top * cfg.marker_pad_fraction * 3.5 * y_range
        for i in visible_idx:
            ax.text(
                sx[i],
                marker_y[i] + label_offset,
                station_names[i],
                fontsize=7,
                rotation=90,
                va="bottom" if toward_top > 0 else "top",
                ha="center",
                color=label_color,
                # Labels deliberately sit above the terrain pins.  Clipping
                # them to the data axes hides the names whenever relief is
                # close to the upper y-limit, which defeats the shared API
                # convention that station names appear at the top.
                clip_on=False,
                zorder=_mstyle.zorder + 1,
            )

    # ── 5. Drop the redundant top frame line ───────────────────────────────
    # The terrain polyline/fill already marks "the surface" on the side
    # the markers and labels sit on; a plain straight box-frame line
    # there besides is not just redundant, it is what the labels/markers
    # visually collide with. Left, right, and bottom keep their normal
    # depth-axis / profile-distance frame -- only the terrain side loses
    # its spine.
    ax.spines["top" if toward_top > 0 else "bottom"].set_visible(False)

    # ── 6. Clear the figure of every marker/label ──────────────────────────
    # Markers and labels are drawn with clip_on=False so the plot boundary
    # never visually clips them, but without enough headroom in ylim they
    # can still hang past the axes into a neighbouring subplot's title or
    # off the edge of the figure. A fixed-fraction pad has no way to know
    # how tall a rotated label actually renders, so this measures the real
    # rendered extent and grows ylim past it, for any caller -- not just
    # this one section style.
    _expand_ylim_for_topside_artists(ax, toward_top)


def _expand_ylim_for_topside_artists(
    ax, toward_top: float, pad_frac: float = 0.02, max_iter: int = 4
) -> None:
    """Grow *ax*'s y-limit on the terrain side so its own frame (the
    spine at the current ``ylim`` edge) clears every station marker and
    label rather than cutting through them.

    Expanding ``ylim`` changes the axes' pixel-to-data scale, which
    changes how much *data* height a text artist's fixed *point* size
    represents -- so a single measure-then-expand pass systematically
    under-shoots (the label was measured against the old, tighter
    scale). This re-measures after each expansion and stops once the
    frame already clears everything, which converges in a couple of
    iterations since each correction is smaller than the last.

    Only text artists (the rotated station-name labels) are measured.
    The marker itself is small and already kept clear by
    ``marker_pad_fraction``; a scatter ``PathCollection``'s
    ``get_window_extent()`` can also return a degenerate, inverted box
    in some matplotlib states, which is not worth working around here.

    Requires a live canvas/renderer to measure real text extents; if
    none is available (e.g. a backend that has not drawn yet and
    cannot be forced to), this silently does nothing rather than
    raising -- the section still renders, just with the older
    fixed-fraction margin.
    """
    fig = ax.figure
    if fig is None or fig.canvas is None or not ax.texts:
        return

    for _ in range(max_iter):
        try:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
        except Exception:
            return

        inv = ax.transData.inverted()
        extreme_y = None
        for artist in ax.texts:
            try:
                bbox = artist.get_window_extent(renderer=renderer)
            except Exception:
                continue
            if not np.isfinite([bbox.x0, bbox.y0, bbox.x1, bbox.y1]).all():
                continue
            _, y0 = inv.transform((bbox.x0, bbox.y0))
            _, y1 = inv.transform((bbox.x1, bbox.y1))
            y_edge = max(y0, y1) if toward_top > 0 else min(y0, y1)
            if extreme_y is None:
                extreme_y = y_edge
            elif toward_top > 0:
                extreme_y = max(extreme_y, y_edge)
            else:
                extreme_y = min(extreme_y, y_edge)

        if extreme_y is None or not np.isfinite(extreme_y):
            return

        ylim = ax.get_ylim()
        lo, hi = min(ylim), max(ylim)
        pad = pad_frac * max(hi - lo, 1e-6)
        ascending = ylim[0] <= ylim[1]

        if toward_top > 0:
            if extreme_y <= hi:
                return  # frame already clears every label/marker
            new_hi = extreme_y + pad
            ax.set_ylim((lo, new_hi) if ascending else (new_hi, lo))
        else:
            if extreme_y >= lo:
                return
            new_lo = extreme_y - pad
            ax.set_ylim((new_lo, hi) if ascending else (hi, new_lo))


def draw_topo_strip(
    fig,
    main_ax,
    chainage_km: np.ndarray,
    elev_m: np.ndarray,
    station_names: Sequence[str] | None = None,
    *,
    cfg=None,
    dark: bool = True,
    marker_style=None,
    facecolor: str | None = None,
) -> Any:
    """Add an elevation-profile strip above a pseudosection image axes.

    For period-vs-station pseudosections the x-axis is station *index*,
    not a real distance.  This function maps station index ↔ chainage so
    the strip correctly displays the terrain shape.

    The strip is inserted by shrinking the main axes and placing a new
    ``Axes`` in the freed space above it.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    main_ax : matplotlib.axes.Axes
        The existing pseudosection image axes.
    chainage_km : array_like (n_stations,)
        Along-profile distances (km) — used as x in the strip.
    elev_m : array_like (n_stations,)
        Terrain elevation at each station (m).
    station_names : sequence of str, optional
        Labels for the strip tick marks.
    cfg : TopoConfig, optional
    dark : bool
    marker_style : pycsamt.api.station.StationMarkerStyle, optional
        Station-pin style override for the strip, analogous to
        :func:`draw_topo_section`'s *marker_style*.  Defaults to
        :data:`pycsamt.api.station.PYCSAMT_STATION_RENDERING`'s
        ``pseudosection`` marker when omitted.
    facecolor : str, optional
        Strip axes background colour.  Defaults to a dark slate in
        ``dark=True`` mode and ``"none"`` (transparent, so the strip
        blends into the figure background rather than sitting inside a
        visibly distinct box) in ``dark=False`` mode.

    Returns
    -------
    ax_strip : matplotlib.axes.Axes
        The newly created elevation-strip axes.
    """
    cfg = _get_cfg(cfg)
    chain = np.asarray(chainage_km, dtype=float)
    elev = np.asarray(elev_m, dtype=float)
    n = len(chain)
    if n == 0:
        return None

    h_ratio = cfg.strip_height_ratio

    # Shrink the main axes to make room at the top
    pos = main_ax.get_position()  # Bbox in figure fraction
    new_bottom = pos.y0
    new_height = pos.height * (1.0 - h_ratio)
    strip_bottom = pos.y0 + new_height
    strip_height = pos.height * h_ratio

    main_ax.set_position([pos.x0, new_bottom, pos.width, new_height])

    ax_s = fig.add_axes(
        [pos.x0, strip_bottom, pos.width, strip_height],
        sharex=None,
    )

    # ── Draw elevation profile in the strip ───────────────────────────────
    x_idx = np.arange(n)  # station index = x in pseudosection
    elev_ex = elev * cfg.exaggeration

    # Filled area
    ax_s.fill_between(
        x_idx,
        elev_ex,
        elev_ex.min() - 1,
        color=cfg.fill_color,
        alpha=cfg.fill_alpha + 0.1,
        linewidth=0,
        zorder=2,
    )
    # Surface line
    if cfg.show_surface_line:
        ax_s.plot(
            x_idx,
            elev_ex,
            color=cfg.line_color,
            linewidth=cfg.line_width,
            zorder=3,
        )
    # Station pins — use the caller's marker style, else the global config
    from pycsamt.api.station import (
        PYCSAMT_STATION_RENDERING,
        StationAxisStyle,
    )

    if marker_style is not None:
        _mstyle = marker_style
    else:
        _mstyle = PYCSAMT_STATION_RENDERING.pseudosection.marker
    ax_s.scatter(
        x_idx,
        elev_ex,
        marker=_mstyle.marker,
        s=_mstyle.size * 0.75,
        facecolors=_mstyle.facecolor,
        edgecolors=_mstyle.edgecolor,
        linewidths=_mstyle.linewidth,
        alpha=_mstyle.alpha,
        zorder=4,
        clip_on=True,
    )

    # ── Style the strip ────────────────────────────────────────────────────
    fg = _DEFAULT_LABEL if dark else "#4c4f69"
    if facecolor is not None:
        bg = facecolor
    else:
        bg = "#181825" if dark else "none"

    ax_s.set_facecolor(bg)
    ax_s.set_xlim(
        main_ax.get_xlim() if hasattr(main_ax, "get_xlim") else (-0.5, n - 0.5)
    )
    ax_s.margins(x=0)
    ax_s.tick_params(
        bottom=False,
        top=False,
        labelbottom=False,
        left=True,
        labelsize=7,
        colors=fg,
    )
    ax_s.set_ylabel("Elev (m)", fontsize=7, color=fg)
    ax_s.spines[:].set_color("#45475a" if dark else "#bcc0cc")
    ax_s.yaxis.label.set_color(fg)
    ax_s.tick_params(axis="y", colors=fg)

    # Station name labels — smart thinning: markers always visible, labels subset
    if station_names and len(station_names) == n:
        _dummy = StationAxisStyle()
        figwidth = (
            ax_s.figure.get_figwidth() if ax_s.figure is not None else 10.0
        )
        visible_idx = _dummy.label_indices(station_names, figwidth)
        for i in visible_idx:
            ax_s.text(
                x_idx[i],
                elev_ex[i],
                f"  {station_names[i]}",
                fontsize=6,
                rotation=90,
                va="bottom",
                ha="center",
                color=fg,
                clip_on=True,
            )

    return ax_s


def add_station_labels(
    ax,
    x_km: np.ndarray,
    y_km: np.ndarray,
    names: Sequence[str],
    *,
    color: str = _DEFAULT_LABEL,
    fontsize: int = 7,
    offset_km: float = 0.05,
    rotation: int = 90,
) -> None:
    """Draw rotated station name labels above marker positions.

    Parameters
    ----------
    ax : Axes
    x_km, y_km : positions of each station marker
    names : station name strings
    color : text color
    fontsize : int
    offset_km : vertical offset in axes units
    rotation : label rotation in degrees
    """
    for xi, yi, name in zip(x_km, y_km, names):
        ax.text(
            xi,
            yi + offset_km,
            name,
            fontsize=fontsize,
            rotation=rotation,
            va="bottom",
            ha="center",
            color=color,
            clip_on=True,
            zorder=7,
        )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _get_cfg(cfg):
    """Return cfg if given, else fall back to global singleton."""
    if cfg is not None:
        return cfg
    from pycsamt.topo.config import PYCSAMT_TOPO

    return PYCSAMT_TOPO
