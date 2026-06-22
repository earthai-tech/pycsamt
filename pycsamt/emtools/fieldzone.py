from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ._core import (
    ensure_sites,
    _iter_items,
    _get_z_block,
    _name,
)

__all__ = [
    "classify_field_zones",
    "near_field_factor",
    "plot_field_zones",
]

# ----------------------------- constants ---------------------------------- #

_MU0 = 4.0 * np.pi * 1e-7        # H/m
_FAR = "far"
_TRANS = "transition"
_NEAR = "near"

_ZONE_COLORS = {
    _FAR:   "#2ca02c",   # green
    _TRANS: "#ff7f0e",   # orange
    _NEAR:  "#d62728",   # red
}

_ZONE_INT = {_FAR: 0, _TRANS: 1, _NEAR: 2}
_INT_ZONE = {0: _FAR, 1: _TRANS, 2: _NEAR}


# ------------------------------ helpers ----------------------------------- #

def _rho_a_det(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    """Geometric-mean apparent resistivity from off-diagonal Z (Ω·m)."""
    rxy = 0.2 * np.abs(z[:, 0, 1]) ** 2 / np.maximum(fr, 1e-24)
    ryx = 0.2 * np.abs(z[:, 1, 0]) ** 2 / np.maximum(fr, 1e-24)
    return np.sqrt(np.maximum(rxy * ryx, 1e-12))


def _bostick_depth(rho: np.ndarray, freq: np.ndarray) -> np.ndarray:
    """δ_B = 356 √(ρ_a / f)  (metres, Bostick approximation)."""
    return 356.0 * np.sqrt(
        np.maximum(rho, 1e-6) / np.maximum(freq, 1e-12)
    )


def _kr_abs(rho: np.ndarray, freq: np.ndarray, offset: float) -> np.ndarray:
    """
    |k·r| = r / δ_B = r / (356 √(ρ_a / f)).

    This is the dimensionless field-zone parameter for CSAMT.
    The plane-wave approximation is valid when |k·r| >> 1.
    """
    db = _bostick_depth(rho, freq)
    return offset / np.maximum(db, 1e-6)


def _zone_labels(
    kr: np.ndarray,
    far_threshold: float,
    near_threshold: float,
) -> np.ndarray:
    labels = np.where(
        kr >= far_threshold, _FAR,
        np.where(kr >= near_threshold, _TRANS, _NEAR),
    )
    return labels


def _nf_correction(
    rho: np.ndarray,
    freq: np.ndarray,
    offset: float,
) -> np.ndarray:
    """
    Near-field correction factor |F(p)| for E_y (equatorial, HED).

    For a horizontal electric dipole over a uniform half-space the ratio
    of the actual E_y to the far-field (plane-wave) approximation is:

        F(p) = 1 − 3/p² + 3/p³   where  p = k·r (complex),
        k = |k|·e^{iπ/4} = (1+i)/√2 · √(ωμ₀/ρ_a)

    F → 1 in the far field (|p| >> 1), deviates strongly in the near
    field.  |F|² gives the bias factor for apparent resistivity.

    References
    ----------
    Chen & Yan (2005) eqs. (8)–(10).
    """
    omega = 2.0 * np.pi * np.maximum(freq, 1e-12)
    k_abs = np.sqrt(omega * _MU0 / np.maximum(rho, 1e-6))
    # complex wave number: k = |k| * exp(i*pi/4)
    p = k_abs * offset * (1.0 + 1j) / np.sqrt(2.0)
    p = np.where(np.abs(p) < 1e-12, 1e-12 * (1.0 + 1j), p)
    F = 1.0 - 3.0 / p ** 2 + 3.0 / p ** 3
    return np.abs(F)


def _resolve_offset(
    ed: Any,
    source_offset: Union[float, Dict[str, float], None],
    station: str,
) -> Optional[float]:
    if isinstance(source_offset, (int, float)):
        return float(source_offset)
    if isinstance(source_offset, dict):
        r = source_offset.get(station)
        if r is not None:
            return float(r)
        low = {k.lower(): v for k, v in source_offset.items()}
        r = low.get(station.lower())
        if r is not None:
            return float(r)
    for attr in ("source_offset", "offset", "x", "east"):
        v = getattr(ed, attr, None)
        if isinstance(v, (int, float)) and float(v) > 0:
            return float(v)
    return None


def _unwrap(ed: Any) -> Any:
    """Unwrap a Sites-level Site object to its underlying EDI-like object."""
    edi = getattr(ed, "edi", None)
    if edi is not None and hasattr(edi, "Z"):
        return edi
    return ed


def _sorted_stations(stations: List[str], S: Any, sort_by: str) -> List[str]:
    if sort_by not in ("lon", "lat"):
        return sorted(stations)
    coords: Dict[str, float] = {}
    for i, ed in enumerate(_iter_items(S)):
        ed = _unwrap(ed)
        nm = _name(ed, i)
        v = getattr(ed, sort_by, None) or getattr(
            ed, "longitude" if sort_by == "lon" else "latitude", None
        )
        coords[nm] = float(v) if v is not None else np.inf
    return sorted(stations, key=lambda s: coords.get(s, np.inf))


# ----------------------------- public API --------------------------------- #

def classify_field_zones(
    sites: Any,
    source_offset: Union[float, Dict[str, float], None] = None,
    *,
    far_threshold: float = 3.0,
    near_threshold: float = 0.3,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Classify CSAMT measurement zones per station per frequency.

    For each (station, frequency) pair the dimensionless field-zone
    parameter is computed as:

        |k·r| = r / δ_B     where  δ_B = 356 √(ρ_a / f)  (metres)

    and the zone is assigned as:

    =========  =========================  =====================
    Zone       Condition                  Meaning
    =========  =========================  =====================
    ``"far"``          |k·r| ≥ far_threshold  Plane-wave approx. valid
    ``"transition"``   near ≤ |k·r| < far    Mixed zone
    ``"near"``         |k·r| < near_threshold Geometric near-field
    =========  =========================  =====================

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
        Any input accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    source_offset : float, dict {station: float}, or None
        Source–receiver separation **r** in metres.  If a dict is
        given, missing stations are skipped (or read from ``ed.offset``
        / ``ed.source_offset``).
    far_threshold : float, default=3.0
        |k·r| threshold for the far-field (plane-wave) zone.
    near_threshold : float, default=0.3
        |k·r| threshold below which the near-field zone is declared.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        Tidy table, one row per (station, frequency):

        ``station``, ``freq_hz``, ``period_s``, ``offset_m``,
        ``rho_a_ohmm``, ``delta_bostick_m``, ``kr``, ``zone``

    References
    ----------
    Chen & Yan (2005), *J. Geophysics and Engineering* 2, 105–120.
    Yan & Fu (2004), analytical shadow/overprint estimation.
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    rows = []
    for i, ed in enumerate(_iter_items(S)):
        ed = _unwrap(ed)
        station = _name(ed, i)
        _, z, fr = _get_z_block(ed)
        if z is None or fr is None:
            continue
        offset = _resolve_offset(ed, source_offset, station)
        if offset is None:
            if verbose > 0:
                import warnings
                warnings.warn(
                    f"classify_field_zones: no source_offset for "
                    f"'{station}' — station skipped.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            continue
        rho_a = _rho_a_det(z, fr)
        kr    = _kr_abs(rho_a, fr, offset)
        db    = _bostick_depth(rho_a, fr)
        zones = _zone_labels(kr, far_threshold, near_threshold)
        for j in range(fr.size):
            rows.append({
                "station":        station,
                "freq_hz":        float(fr[j]),
                "period_s":       1.0 / max(float(fr[j]), 1e-12),
                "offset_m":       offset,
                "rho_a_ohmm":     float(rho_a[j]),
                "delta_bostick_m": float(db[j]),
                "kr":             float(kr[j]),
                "zone":           str(zones[j]),
            })

    _COLS = [
        "station", "freq_hz", "period_s", "offset_m",
        "rho_a_ohmm", "delta_bostick_m", "kr", "zone",
    ]
    if not rows:
        return pd.DataFrame(columns=_COLS)
    return pd.DataFrame(rows, columns=_COLS)


def near_field_factor(
    sites: Any,
    source_offset: Union[float, Dict[str, float], None] = None,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Near-field correction factor for apparent resistivity (equatorial HED).

    For the equatorial E_y component from a horizontal electric dipole
    over a homogeneous half-space, the ratio of the measured E_y to the
    plane-wave (far-field) E_y is:

        F(p) = 1 − 3/p² + 3/p³     p = k·r  (complex)

    so the apparent resistivity is biased by factor |F(p)|².  When
    |F(p)| ≈ 1 the data are in the plane-wave zone; strong departures
    indicate near-field contamination.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
    source_offset : float, dict {station: float}, or None
        Source–receiver separation **r** in metres.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        Columns: ``station``, ``freq_hz``, ``period_s``, ``offset_m``,
        ``rho_a_ohmm``, ``kr``, ``nf_factor``.

        * ``nf_factor`` = |F(p)|; close to 1.0 → far-field (safe).
        * ``nf_factor`` far from 1.0 → near-field bias present.

    References
    ----------
    Chen & Yan (2005), eqs. (8)–(10).
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    rows = []
    for i, ed in enumerate(_iter_items(S)):
        ed = _unwrap(ed)
        station = _name(ed, i)
        _, z, fr = _get_z_block(ed)
        if z is None or fr is None:
            continue
        offset = _resolve_offset(ed, source_offset, station)
        if offset is None:
            continue
        rho_a = _rho_a_det(z, fr)
        kr    = _kr_abs(rho_a, fr, offset)
        nff   = _nf_correction(rho_a, fr, offset)
        for j in range(fr.size):
            rows.append({
                "station":    station,
                "freq_hz":    float(fr[j]),
                "period_s":   1.0 / max(float(fr[j]), 1e-12),
                "offset_m":   offset,
                "rho_a_ohmm": float(rho_a[j]),
                "kr":         float(kr[j]),
                "nf_factor":  float(nff[j]),
            })

    _COLS = [
        "station", "freq_hz", "period_s", "offset_m",
        "rho_a_ohmm", "kr", "nf_factor",
    ]
    if not rows:
        return pd.DataFrame(columns=_COLS)
    return pd.DataFrame(rows, columns=_COLS)


def plot_field_zones(
    sites: Any,
    source_offset: Union[float, Dict[str, float], None] = None,
    *,
    far_threshold: float = 3.0,
    near_threshold: float = 0.3,
    contour_kr: bool = True,
    kr_levels: Tuple[float, ...] = (0.1, 0.3, 1.0, 3.0, 10.0),
    sort_by: str = "name",
    period_axis: bool = True,
    log_y: bool = True,
    figsize: Tuple[float, float] = (10.0, 5.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[Any] = None,
) -> Any:
    """
    Pseudosection of CSAMT field zones across stations and frequencies.

    Each cell (station × period/frequency) is filled by zone colour:

    * **green**  → far field  (|k·r| ≥ *far_threshold*)
    * **orange** → transition
    * **red**    → near field (|k·r| < *near_threshold*)

    Dashed white contours of constant |k·r| can be overlaid.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
    source_offset : float, dict {station: float}, or None
        Source–receiver separation **r** in metres.
    far_threshold, near_threshold : float
        Zone boundaries in |k·r|.
    contour_kr : bool, default=True
        Draw |k·r| contours.
    kr_levels : tuple of float
        |k·r| values to contour.
    sort_by : {"name", "lon", "lat"}
        Station ordering along the x-axis.
    period_axis : bool, default=True
        If True y-axis shows period (s), else frequency (Hz).
    log_y : bool, default=True
        Use a quasi-log y-axis (log-spaced tick labels).
    figsize : (float, float), default=(10, 5)
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.
    ax : matplotlib.axes.Axes, optional
        Draw on existing axes.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.patches as mpatches

    df = classify_field_zones(
        sites,
        source_offset,
        far_threshold=far_threshold,
        near_threshold=near_threshold,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if df.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes)
        return ax

    # station order
    stations = df["station"].unique().tolist()
    if sort_by in ("lon", "lat"):
        S = ensure_sites(
            sites, recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )
        stations = _sorted_stations(stations, S, sort_by)
    else:
        stations = sorted(stations)

    y_key = "period_s" if period_axis else "freq_hz"
    all_y = np.sort(df[y_key].unique())

    # fill grid
    grid_zone = np.full((len(all_y), len(stations)), np.nan)
    grid_kr   = np.full((len(all_y), len(stations)), np.nan)
    y_idx = {v: k for k, v in enumerate(all_y)}
    x_idx = {s: k for k, s in enumerate(stations)}

    for row in df.itertuples(index=False):
        yi = y_idx.get(getattr(row, y_key))
        xi = x_idx.get(row.station)
        if yi is not None and xi is not None:
            grid_zone[yi, xi] = _ZONE_INT[row.zone]
            grid_kr[yi, xi]   = row.kr

    cmap = mcolors.ListedColormap(
        [_ZONE_COLORS[_FAR], _ZONE_COLORS[_TRANS], _ZONE_COLORS[_NEAR]]
    )
    xs = np.arange(len(stations) + 1) - 0.5
    ys = np.arange(len(all_y) + 1) - 0.5
    ax.pcolormesh(xs, ys, grid_zone, cmap=cmap, vmin=-0.5, vmax=2.5,
                  shading="auto")

    if (contour_kr
            and not np.all(np.isnan(grid_kr))
            and grid_kr.shape[0] >= 2
            and grid_kr.shape[1] >= 2):
        XX, YY = np.meshgrid(np.arange(len(stations)), np.arange(len(all_y)))
        valid_levels = [lv for lv in sorted(kr_levels)
                        if np.nanmin(grid_kr) < lv < np.nanmax(grid_kr)]
        if valid_levels:
            cs = ax.contour(XX, YY, grid_kr, levels=valid_levels,
                            colors="white", linewidths=0.8, linestyles="--")
            ax.clabel(cs, fmt="|kr|=%.2g", fontsize=7, inline=True)

    ax.set_xticks(range(len(stations)))
    ax.set_xticklabels(stations, rotation=45, ha="right", fontsize=8)

    # y-axis ticks
    n_ytick = min(8, len(all_y))
    step    = max(1, len(all_y) // n_ytick)
    tick_idx = np.arange(0, len(all_y), step)
    ax.set_yticks(tick_idx)
    ax.set_yticklabels([f"{all_y[k]:.3g}" for k in tick_idx], fontsize=8)
    ax.set_ylabel("Period (s)" if period_axis else "Frequency (Hz)")
    ax.set_xlabel("Station")
    ax.set_title("CSAMT Field Zone Classification (|k·r|)")

    patches = [
        mpatches.Patch(
            color=_ZONE_COLORS[_FAR],
            label=f"Far field  |kr|≥{far_threshold}",
        ),
        mpatches.Patch(
            color=_ZONE_COLORS[_TRANS],
            label=f"Transition  {near_threshold}≤|kr|<{far_threshold}",
        ),
        mpatches.Patch(
            color=_ZONE_COLORS[_NEAR],
            label=f"Near field  |kr|<{near_threshold}",
        ),
    ]
    ax.legend(handles=patches, loc="upper right", fontsize=8, framealpha=0.85)
    return ax
