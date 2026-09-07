# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Style helpers for pyCSAMT map renderers."""

from __future__ import annotations

import math
from collections.abc import Sequence

from .config import MapTheme

# "Unclassified" band: shown wherever a geology legend leaves a gap in its
# resistivity coverage (see :func:`geology_colorscale`) rather than letting
# Plotly interpolate a meaningless gradient between two unrelated
# lithology colours.
GEOLOGY_GAP_COLOR = "#9aa0a6"

PLOTLY_CMAP_REMAP: dict[str, str] = {
    "coolwarm": "balance",
    "seismic": "rdbu",
    "terrain": "earth",
    "YlOrRd": "ylorrd",
    "RdBu_r": "rdbu_r",
    "copper": "thermal",
    "gnuplot2": "turbid",
    "tab10": "plotly3",
    "tab20": "plasma",
}

THEME_COLORS: dict[MapTheme, dict[str, str]] = {
    "dark": {
        "paper": "#1e1e2e",
        "plot": "#181825",
        "text": "#cdd6f4",
        "grid": "#313244",
        "accent": "#89b4fa",
    },
    "light": {
        "paper": "#eff1f5",
        "plot": "#ffffff",
        "text": "#4c4f69",
        "grid": "#ccd0da",
        "accent": "#1e66f5",
    },
    "publication": {
        "paper": "#ffffff",
        "plot": "#ffffff",
        "text": "#111111",
        "grid": "#cccccc",
        "accent": "#1565c0",
    },
}


def to_plotly_cmap(
    cmap: str | None,
    fallback: str = "plasma",
) -> str:
    """Return a Plotly-compatible colorscale name."""
    if not cmap:
        return fallback
    return PLOTLY_CMAP_REMAP.get(cmap, cmap)


# ---------------------------------------------------------------------------
# Geology-legend (Interpretation overlay) discrete colour banding
# ---------------------------------------------------------------------------
# A "band" is a plain ``(rho_min, rho_max, hex_color)`` tuple, Ω·m linear.


def geology_crange(
    bands: Sequence[tuple[float, float, str]],
) -> tuple[float, float]:
    """Return the ``(log10_lo, log10_hi)`` span covered by *bands*.

    *bands* is a plain, format-agnostic ``(rho_min, rho_max, color)`` list
    in linear Ω·m -- the app layer builds it from a
    :class:`pycsamt.geology.lithology.RockDatabase` /
    :class:`pycsamt.format.geology.GeologyLegend`; this module stays free
    of any format/geology import. Non-positive or degenerate bands are
    ignored. Falls back to ``(0.0, 1.0)`` when no band is usable.
    """
    usable = [
        (float(b[0]), float(b[1]))
        for b in bands
        if float(b[0]) > 0 and float(b[1]) > float(b[0])
    ]
    if not usable:
        return (0.0, 1.0)
    los = [lo for lo, _ in usable]
    his = [hi for _, hi in usable]
    return (math.log10(min(los)), math.log10(max(his)))


def geology_colorscale(
    bands: Sequence[tuple[float, float, str]],
    *,
    log_lo: float | None = None,
    log_hi: float | None = None,
    gap_color: str = GEOLOGY_GAP_COLOR,
    textured: dict[str, int] | None = None,
) -> list[list]:
    """Build a hard-stepped Plotly colorscale from a geology legend.

    Returns ``[[t, color], ...]`` stops with ``t`` in ``[0, 1]`` normalized
    against ``(log_lo, log_hi)`` (defaults to :func:`geology_crange` of
    *bands*), using repeated ``t`` values at each band boundary so Plotly
    renders flat colour bands instead of a smooth gradient -- pass the
    same ``(log_lo, log_hi)`` as ``cmin``/``cmax`` on the trace, in
    log10(Ω·m), so the bands line up with the data. Any part of the
    ``[log_lo, log_hi]`` span not covered by a band (including gaps
    between bands) renders as *gap_color* rather than an interpolated
    blend of neighbouring lithology colours.

    *textured* optionally maps a band's name (its 4th tuple element) to
    a stop count -- that band's sub-range is filled with an N-stop
    :func:`pattern_band_stops` gradient (light tint of the band colour
    -> full colour) instead of one flat colour. On its own this still
    just renders a gradient; the actual *pattern* only appears once the
    caller also overrides each cell's colour-axis value within that
    sub-range by a real spatial pattern-stencil sample (see
    :func:`pycsamt.map.volume._inject_pattern_values`) instead of its
    raw resistivity -- this function only prepares the palette that
    sampling then reads from.
    """
    lo, hi = (log_lo, log_hi) if log_lo is not None and log_hi is not None \
        else geology_crange(bands)
    if not (hi > lo):
        hi = lo + 1.0

    def _t(value: float) -> float:
        return min(1.0, max(0.0, (value - lo) / (hi - lo)))

    textured = textured or {}
    segments: list[tuple[float, float, str, str]] = []
    for band in bands:
        rho_min, rho_max, color = float(band[0]), float(band[1]), band[2]
        name = str(band[3]) if len(band) >= 4 else ""
        if not (rho_min > 0 and rho_max > rho_min):
            continue
        t0, t1 = _t(math.log10(rho_min)), _t(math.log10(rho_max))
        if t1 <= t0:
            continue
        segments.append((t0, t1, str(color), name))
    segments.sort(key=lambda s: s[0])

    stops: list[list] = []
    cursor = 0.0
    for t0, t1, color, name in segments:
        if t0 > cursor:
            stops += [[cursor, gap_color], [t0, gap_color]]
        n_stops = textured.get(name) if name else None
        if n_stops and n_stops > 1:
            palette = pattern_band_stops(color, int(n_stops))
            for i, stop_color in enumerate(palette):
                tt = t0 + (t1 - t0) * (i / (len(palette) - 1))
                stops.append([tt, stop_color])
        else:
            stops += [[t0, color], [t1, color]]
        cursor = max(cursor, t1)
    if cursor < 1.0:
        stops += [[cursor, gap_color], [1.0, gap_color]]
    if not stops:
        stops = [[0.0, gap_color], [1.0, gap_color]]
    stops[0][0] = 0.0
    stops[-1][0] = 1.0
    return stops


def pattern_band_stops(base_color: str, n: int = 16) -> list[str]:
    """Return *n* hex colours from a light tint of *base_color* ("no
    ink" -- a pattern stencil's background) to the full *base_color*
    ("full ink"), for texturing one geology band with a pattern
    stencil sampled to a per-cell density in ``[0, 1]``.

    A tile is read purely as a grayscale/alpha *stencil* (see
    :func:`pycsamt.geology.patterns` tile handling), not as a literal
    colour source -- the legend's own assigned colour is what actually
    tints the rendered pattern, so a hatch/stipple/brick swatch reads
    in the same colour its legend chip already shows, whatever the
    source pack's original ink colour was.
    """
    n = max(2, int(n))
    r, g, b = _hex_to_rgb(base_color)
    lr = round(r + (255 - r) * 0.85)
    lg = round(g + (255 - g) * 0.85)
    lb = round(b + (255 - b) * 0.85)
    out = []
    for i in range(n):
        t = i / (n - 1)
        rr = round(lr + (r - lr) * t)
        gg = round(lg + (g - lg) * t)
        bb = round(lb + (b - lb) * t)
        out.append(f"#{rr:02X}{gg:02X}{bb:02X}")
    return out


def _hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    text = str(hex_color).lstrip("#")
    if len(text) != 6:
        return (128, 128, 128)
    try:
        r, g, b = (int(text[i : i + 2], 16) for i in (0, 2, 4))
    except ValueError:
        return (128, 128, 128)
    return (r, g, b)


def geology_colorbar_ticks(
    bands: Sequence[tuple[float, float, str] | tuple[float, float, str, str]],
) -> tuple[list[float], list[str]]:
    """Return ``(tickvals, ticktext)`` labelling each band by name.

    *bands* is the same list :func:`geology_colorscale` takes, with an
    optional 4th element, the band's display name (e.g. a rock name).
    Each tick sits at the band's own log10 geometric-mean resistivity so
    it lands inside its colour, not on a boundary shared with its
    neighbour. Bands with no name (a plain 3-tuple) are skipped -- pass
    this straight to a trace's ``colorbar=dict(tickvals=..., ticktext=
    ...)`` in place of the default numeric Ω·m ticks, "masking" the
    resistivity reading with the geology legend itself.
    """
    tickvals: list[float] = []
    ticktext: list[str] = []
    for band in bands:
        if len(band) < 4 or not str(band[3]).strip():
            continue
        rho_min, rho_max = float(band[0]), float(band[1])
        if not (rho_min > 0 and rho_max > rho_min):
            continue
        mid = math.sqrt(rho_min * rho_max)
        tickvals.append(math.log10(mid))
        ticktext.append(str(band[3]))
    return tickvals, ticktext


def geology_legend_shapes_annotations(
    bands: Sequence[tuple[float, float, str] | tuple[float, float, str, str]],
    *,
    x: float = 1.02,
    y_top: float = 0.95,
    row_height: float = 0.05,
    min_row_height: float = 0.026,
    font_size: int = 11,
    min_font_size: int = 8,
) -> tuple[list[dict], list[dict]]:
    """Build a non-overlapping on-canvas geology legend: one coloured
    "chip" annotation per band, stacked at a fixed row height in the
    figure's right margin (``xref="paper"``, so the caller must reserve
    real margin space -- e.g. ``margin=dict(r=220)``).

    Unlike labelling each band at its own resistivity position on a
    colourbar (:func:`geology_colorbar_ticks`), row position here is by
    *index*, not by value -- two bands whose resistivity ranges sit
    close together (or a legend with many entries on a short colourbar)
    can never overlap, only ever crowd, and rows shrink toward
    *min_row_height*/*min_font_size* before that happens. If even the
    minimum still cannot fit every band in ``[0, y_top]``, the remainder
    collapses into one final "+N more" row rather than overlapping or
    running off the canvas.

    Returns ``(shapes, annotations)`` -- *shapes* is always empty (the
    chip's own background covers the swatch), kept for API symmetry with
    a future non-text legend element.
    """
    entries: list[tuple[str, str]] = []
    for band in bands:
        name = (
            str(band[3]).strip()
            if len(band) >= 4 and str(band[3]).strip()
            else f"{float(band[0]):g}–{float(band[1]):g} Ω·m"
        )
        entries.append((name, str(band[2])))
    if not entries:
        return [], []

    n = len(entries)
    available = max(y_top - 0.02, min_row_height)
    # Row height shrinks toward min_row_height as more bands need to
    # fit; once even that minimum can't seat all of them, cap the row
    # *count* instead (one slot reserved for a "+N more" summary) so
    # rows never have to overlap to make room.
    rh = max(min_row_height, min(row_height, available / n))
    max_rows = max(1, int(available / rh))
    fs = font_size if rh >= row_height else max(
        min_font_size, round(font_size * (rh / row_height))
    )

    truncated = n > max_rows
    n_shown = (max_rows - 1) if truncated else n
    n_shown = max(0, n_shown)

    annotations: list[dict] = []
    for i in range(n_shown):
        name, color = entries[i]
        annotations.append(
            dict(
                x=x, y=y_top - i * rh, xref="paper", yref="paper",
                xanchor="left", yanchor="middle",
                text=f"  {name}  ", showarrow=False, align="left",
                font=dict(size=fs, color=_contrast_text_color(color)),
                bgcolor=color, bordercolor="rgba(0,0,0,0.25)",
                borderwidth=1, borderpad=3,
            )
        )
    if truncated:
        annotations.append(
            dict(
                x=x, y=y_top - n_shown * rh, xref="paper", yref="paper",
                xanchor="left", yanchor="middle",
                text=f"+{n - n_shown} more…", showarrow=False,
                font=dict(size=fs), align="left",
            )
        )
    return [], annotations


def _contrast_text_color(hex_color: str) -> str:
    """Black or white, whichever reads better on *hex_color*."""
    text = str(hex_color).lstrip("#")
    if len(text) != 6:
        return "#111111"
    try:
        r, g, b = (int(text[i : i + 2], 16) for i in (0, 2, 4))
    except ValueError:
        return "#111111"
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "#111111" if luminance > 140 else "#f5f5f5"


def theme_colors(theme: MapTheme = "light") -> dict[str, str]:
    """Return a copied color palette for *theme*."""
    colors = THEME_COLORS.get(theme, THEME_COLORS["light"])
    return dict(colors)


def is_dark_theme(theme: MapTheme = "light") -> bool:
    """Return True for dark map styling."""
    return theme == "dark"
