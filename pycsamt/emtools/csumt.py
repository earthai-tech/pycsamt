"""Controlled-source ultra-audio MT depth and survey-planning tools.

The module converts apparent resistivity and transmitter frequency into
Bostick depth estimates, evaluates vertical resolution and depth coverage,
designs CSUMT frequency schedules, and plots depth sections for survey sites.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..api.station import PYCSAMT_STATION_RENDERING
from ._core import (
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
)

__all__ = [
    # pure survey-planning functions (no sites required)
    "bostick_depth_from_rho",
    "vertical_resolution_pair",
    "frequency_for_depth",
    "frequency_schedule",
    # sites-based analysis
    "bostick_depth",
    "vertical_resolution",
    "depth_coverage_table",
    "plot_depth_section",
]

# ----------------------------- constants ---------------------------------- #

BOSTICK_CONST: float = 356.0
"""
Bostick depth constant in metres:
    D(f) = 356 × √(ρ_a / f)

Derived from the skin depth δ = 503√(ρ/f) as D_B = δ / √2 ≈ 356√(ρ/f).
"""

F_MIN_CSUMT: float = 9.6e3
"""Lower bound of the CSUMT frequency range: 9.6 kHz (zhang2025)."""

F_MAX_CSUMT: float = 614.4e3
"""Upper bound of the CSUMT frequency range: 614.4 kHz (zhang2025)."""


# ------------------------------ helpers ----------------------------------- #


def _unwrap(ed: Any) -> Any:
    """Unwrap a Sites-level Site wrapper to the underlying EDI-like object."""
    edi = getattr(ed, "edi", None)
    if edi is not None and hasattr(edi, "Z"):
        return edi
    return ed


def _rho_a_det(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    """Geometric-mean apparent resistivity from off-diagonal Z (Ω·m)."""
    rxy = 0.2 * np.abs(z[:, 0, 1]) ** 2 / np.maximum(fr, 1e-24)
    ryx = 0.2 * np.abs(z[:, 1, 0]) ** 2 / np.maximum(fr, 1e-24)
    return np.sqrt(np.maximum(rxy * ryx, 1e-12))


# ====================== pure survey-planning API ========================== #


def bostick_depth_from_rho(
    rho: float | np.ndarray,
    freq: float | np.ndarray,
) -> np.ndarray:
    """
    Bostick depth estimate D(f) from apparent resistivity.

        D(f) = 356 × √(ρ_a / f)     [metres]

    Parameters
    ----------
    rho : float or array
        Apparent resistivity ρ_a in Ω·m.  Broadcastable with *freq*.
    freq : float or array
        Frequency in Hz.

    Returns
    -------
    numpy.ndarray
        Depth in metres, same shape as broadcast of *rho* and *freq*.

    References
    ----------
    Zhang et al. (2025), Eq. (1), *Measurement*.
    """
    rho = np.asarray(rho, dtype=float)
    freq = np.asarray(freq, dtype=float)
    return BOSTICK_CONST * np.sqrt(
        np.maximum(rho, 0.0) / np.maximum(freq, 1e-24)
    )


def vertical_resolution_pair(
    rho: float,
    f_lo: float,
    f_hi: float,
) -> float:
    """
    Vertical resolution ΔD between two adjacent frequencies.

        ΔD = 356 × √ρ_c × (1/√f_lo − 1/√f_hi)     [metres; f_lo < f_hi]

    Parameters
    ----------
    rho : float
        Characteristic (apparent) resistivity ρ_c in Ω·m.
    f_lo : float
        Lower frequency in Hz  (deeper penetration).
    f_hi : float
        Higher frequency in Hz (shallower penetration).

    Returns
    -------
    float
        Vertical resolution in metres.  Positive when f_lo < f_hi.

    References
    ----------
    Zhang et al. (2025), Eq. (2), *Measurement*.
    """
    rho = max(float(rho), 0.0)
    f_lo = max(float(f_lo), 1e-24)
    f_hi = max(float(f_hi), 1e-24)
    return (
        BOSTICK_CONST
        * np.sqrt(rho)
        * (1.0 / np.sqrt(f_lo) - 1.0 / np.sqrt(f_hi))
    )


def frequency_for_depth(
    depth_m: float | np.ndarray,
    rho: float,
) -> np.ndarray:
    """
    Invert the Bostick formula: return the frequency (Hz) that maps to
    a given depth for a background resistivity *rho*.

        f = ρ × (356 / D)²     [Hz]

    Parameters
    ----------
    depth_m : float or array
        Target depth(s) in metres.
    rho : float
        Background apparent resistivity in Ω·m.

    Returns
    -------
    numpy.ndarray
        Frequency in Hz, same shape as *depth_m*.
    """
    d = np.asarray(depth_m, dtype=float)
    return float(rho) * (BOSTICK_CONST / np.maximum(d, 1e-6)) ** 2


def frequency_schedule(
    target_depths: float | np.ndarray,
    rho_estimate: float,
    *,
    f_min: float = F_MIN_CSUMT,
    f_max: float = F_MAX_CSUMT,
    min_resolution_m: float | None = None,
    fill_decades: bool = False,
    per_decade: int = 3,
    as_khz: bool = False,
) -> np.ndarray:
    """
    Design a CSUMT frequency schedule that samples a set of target depths.

    Each target depth is converted to a frequency via
    :func:`frequency_for_depth`, then clipped to [*f_min*, *f_max*].
    Optionally, intermediate frequencies can be inserted to guarantee a
    minimum vertical resolution between consecutive depth levels.

    Parameters
    ----------
    target_depths : float or array
        Target depths in metres (deepest first or any order — sorted
        internally).
    rho_estimate : float
        Background apparent resistivity ρ (Ω·m) used for the conversion.
    f_min : float, default=9.6e3
        Minimum transmitter frequency in Hz (lower bound of CSUMT range).
    f_max : float, default=614.4e3
        Maximum transmitter frequency in Hz (upper bound of CSUMT range).
    min_resolution_m : float or None
        If given, insert additional frequencies between adjacent target
        depths whenever their vertical resolution would exceed this value.
    fill_decades : bool, default=False
        If True, add *per_decade* log-spaced frequencies within each
        decade of the schedule to smooth coverage.
    per_decade : int, default=3
        Number of extra frequencies to insert per decade when
        *fill_decades* is True.
    as_khz : bool, default=False
        If True, return frequencies in kHz instead of Hz.

    Returns
    -------
    numpy.ndarray
        Sorted frequencies in Hz (or kHz if *as_khz* is True).

    References
    ----------
    Zhang et al. (2025), "Controlled source ultra-audio frequency
    magnetotellurics (CSUMT) transmitter", *Measurement*.
    """
    depths = np.unique(np.asarray(target_depths, dtype=float).ravel())
    depths = depths[depths > 0]
    if depths.size == 0:
        return np.array([], dtype=float)

    freqs = frequency_for_depth(depths, rho_estimate)
    freqs = freqs[(freqs >= f_min) & (freqs <= f_max)]
    freqs = np.unique(freqs)

    if min_resolution_m is not None and freqs.size >= 2:
        extra: list[float] = []
        fs = np.sort(freqs)
        for i in range(len(fs) - 1):
            f_lo, f_hi = fs[i], fs[i + 1]
            rho_c = rho_estimate
            delta = vertical_resolution_pair(rho_c, f_lo, f_hi)
            if delta > min_resolution_m:
                n_insert = int(np.ceil(delta / min_resolution_m)) - 1
                extra.extend(
                    np.logspace(
                        np.log10(f_lo),
                        np.log10(f_hi),
                        n_insert + 2,
                    )[1:-1].tolist()
                )
        if extra:
            freqs = np.unique(np.concatenate([freqs, extra]))
            freqs = freqs[(freqs >= f_min) & (freqs <= f_max)]

    if fill_decades and freqs.size >= 2:
        lo_d = np.log10(freqs.min())
        hi_d = np.log10(freqs.max())
        n = max(2, int(np.ceil((hi_d - lo_d) * per_decade)))
        fill = np.logspace(lo_d, hi_d, n)
        fill = fill[(fill >= f_min) & (fill <= f_max)]
        freqs = np.unique(np.concatenate([freqs, fill]))

    freqs = np.sort(freqs)
    if as_khz:
        return freqs / 1e3
    return freqs


# ======================= sites-based analysis API ========================= #


def bostick_depth(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Bostick depth estimate per station per frequency from measured data.

    Uses the apparent resistivity derived from the off-diagonal impedance
    tensor components (geometric mean):

        D(f) = 356 × √(ρ_a(f) / f)     [metres]

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
        Any input accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        One row per (station, frequency) with columns:
        ``station``, ``freq_hz``, ``period_s``,
        ``rho_a_ohmm``, ``depth_m``.

    References
    ----------
    Zhang et al. (2025), Eq. (1).
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    _COLS = ["station", "freq_hz", "period_s", "rho_a_ohmm", "depth_m"]
    rows = []
    for i, ed in enumerate(_iter_items(S)):
        ed = _unwrap(ed)
        station = _name(ed, i)
        _, z, fr = _get_z_block(ed)
        if z is None or fr is None:
            continue
        rho_a = _rho_a_det(z, fr)
        depth = bostick_depth_from_rho(rho_a, fr)
        for j in range(fr.size):
            rows.append(
                {
                    "station": station,
                    "freq_hz": float(fr[j]),
                    "period_s": 1.0 / max(float(fr[j]), 1e-24),
                    "rho_a_ohmm": float(rho_a[j]),
                    "depth_m": float(depth[j]),
                }
            )
    if not rows:
        return pd.DataFrame(columns=_COLS)
    return pd.DataFrame(rows, columns=_COLS)


def vertical_resolution(
    sites: Any,
    *,
    rho_override: float | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Vertical resolution ΔD between adjacent frequencies per station.

    For each consecutive pair (f_lo, f_hi) in the station's frequency
    list (sorted ascending), computes:

        ΔD = D(f_lo) − D(f_hi)     [metres]

    using the Bostick depths derived from the measured ρ_a.
    Alternatively, supply *rho_override* to use a fixed background
    resistivity with the analytical formula
    ``356 × √ρ × (1/√f_lo − 1/√f_hi)``.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
    rho_override : float or None
        If given, use this constant resistivity for all ΔD calculations
        (analytical formula) instead of the per-frequency ρ_a.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        Columns: ``station``, ``freq_lo_hz``, ``freq_hi_hz``,
        ``depth_lo_m``, ``depth_hi_m``, ``delta_depth_m``,
        ``rho_a_ohmm``.
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    _COLS = [
        "station",
        "freq_lo_hz",
        "freq_hi_hz",
        "depth_lo_m",
        "depth_hi_m",
        "delta_depth_m",
        "rho_a_ohmm",
    ]
    rows = []
    for i, ed in enumerate(_iter_items(S)):
        ed = _unwrap(ed)
        station = _name(ed, i)
        _, z, fr = _get_z_block(ed)
        if z is None or fr is None or fr.size < 2:
            continue
        # sort ascending by frequency
        order = np.argsort(fr)
        fr_s = fr[order]
        z_s = z[order]
        rho_a = _rho_a_det(z_s, fr_s)

        for j in range(fr_s.size - 1):
            f_lo = float(fr_s[j])
            f_hi = float(fr_s[j + 1])
            if rho_override is not None:
                rho_c = float(rho_override)
                d_lo = bostick_depth_from_rho(rho_c, f_lo)
                d_hi = bostick_depth_from_rho(rho_c, f_hi)
            else:
                rho_c = float(np.sqrt(rho_a[j] * rho_a[j + 1]))
                d_lo = bostick_depth_from_rho(float(rho_a[j]), f_lo)
                d_hi = bostick_depth_from_rho(float(rho_a[j + 1]), f_hi)
            rows.append(
                {
                    "station": station,
                    "freq_lo_hz": f_lo,
                    "freq_hi_hz": f_hi,
                    "depth_lo_m": float(d_lo),
                    "depth_hi_m": float(d_hi),
                    "delta_depth_m": float(d_lo - d_hi),
                    "rho_a_ohmm": rho_c,
                }
            )
    if not rows:
        return pd.DataFrame(columns=_COLS)
    return pd.DataFrame(rows, columns=_COLS)


def depth_coverage_table(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Summary depth-coverage statistics per station.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        One row per station with columns:
        ``station``, ``n_freq``, ``freq_min_hz``, ``freq_max_hz``,
        ``depth_min_m``, ``depth_max_m``,
        ``mean_resolution_m``, ``median_resolution_m``.
    """
    bd = bostick_depth(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    vr = vertical_resolution(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    _COLS = [
        "station",
        "n_freq",
        "freq_min_hz",
        "freq_max_hz",
        "depth_min_m",
        "depth_max_m",
        "mean_resolution_m",
        "median_resolution_m",
    ]
    if bd.empty:
        return pd.DataFrame(columns=_COLS)

    rows = []
    for st in bd["station"].unique():
        sub_bd = bd[bd["station"] == st]
        sub_vr = vr[vr["station"] == st] if not vr.empty else pd.DataFrame()

        dres = (
            sub_vr["delta_depth_m"].dropna().values
            if not sub_vr.empty
            else np.array([])
        )
        rows.append(
            {
                "station": st,
                "n_freq": len(sub_bd),
                "freq_min_hz": float(sub_bd["freq_hz"].min()),
                "freq_max_hz": float(sub_bd["freq_hz"].max()),
                "depth_min_m": float(sub_bd["depth_m"].min()),
                "depth_max_m": float(sub_bd["depth_m"].max()),
                "mean_resolution_m": float(np.nanmean(dres))
                if dres.size
                else float("nan"),
                "median_resolution_m": float(np.nanmedian(dres))
                if dres.size
                else float("nan"),
            }
        )
    return pd.DataFrame(rows, columns=_COLS)


def plot_depth_section(
    sites: Any,
    *,
    log_color: bool = True,
    sort_by: str = "name",
    cmap: str = "viridis_r",
    figsize: tuple[float, float] = (10.0, 5.0),
    period_axis: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Any | None = None,
) -> Any:
    """
    Pseudosection of Bostick depth across stations and periods/frequencies.

    Each cell (station × period) is coloured by the Bostick depth
    D(f) = 356 √(ρ_a / f).

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
    log_color : bool, default=True
        Color by log10(depth) instead of depth.
    sort_by : {"name", "lon", "lat"}
        Station ordering along the x-axis.
    cmap : str, default="viridis_r"
        Matplotlib colormap name.
    figsize : (float, float), default=(10, 5)
    period_axis : bool, default=True
        If True y-axis is period (s); otherwise frequency (Hz).
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    df = bostick_depth(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if df.empty:
        ax.text(
            0.5,
            0.5,
            "no data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    # station order
    stations = df["station"].unique().tolist()
    if sort_by in ("lon", "lat"):
        S = ensure_sites(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
        coords: dict[str, float] = {}
        for ii, ed in enumerate(_iter_items(S)):
            ed = _unwrap(ed)
            nm = _name(ed, ii)
            v = getattr(ed, sort_by, None) or getattr(
                ed, "longitude" if sort_by == "lon" else "latitude", None
            )
            coords[nm] = float(v) if v is not None else float("inf")
        stations = sorted(stations, key=lambda s: coords.get(s, float("inf")))
    else:
        stations = sorted(stations)

    y_key = "period_s" if period_axis else "freq_hz"
    all_y = np.sort(df[y_key].unique())

    grid = np.full((len(all_y), len(stations)), np.nan)
    y_idx = {v: k for k, v in enumerate(all_y)}
    x_idx = {s: k for k, s in enumerate(stations)}

    for row in df.itertuples(index=False):
        yi = y_idx.get(getattr(row, y_key))
        xi = x_idx.get(row.station)
        if yi is not None and xi is not None:
            grid[yi, xi] = row.depth_m

    plot_data = np.log10(np.maximum(grid, 1e-3)) if log_color else grid

    xs = np.arange(len(stations) + 1) - 0.5
    ys = np.arange(len(all_y) + 1) - 0.5
    im = ax.pcolormesh(xs, ys, plot_data, cmap=cmap, shading="auto")

    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(len(stations), dtype=float),
        stations,
        preset="pseudosection",
        xlim=(-0.5, len(stations) - 0.5),
    )

    n_ytick = min(8, len(all_y))
    step = max(1, len(all_y) // n_ytick)
    tick_idx = np.arange(0, len(all_y), step)
    ax.set_yticks(tick_idx)
    ax.set_yticklabels([f"{all_y[k]:.3g}" for k in tick_idx], fontsize=8)
    ax.set_ylabel("Period (s)" if period_axis else "Frequency (Hz)")
    if period_axis and not ax.yaxis_inverted():
        ax.invert_yaxis()
    ax.set_title("Bostick Depth Section  D = 356√(ρ_a / f)")

    cb = plt.colorbar(im, ax=ax)
    cb.set_label("log10 depth (m)" if log_color else "depth (m)")
    return ax
