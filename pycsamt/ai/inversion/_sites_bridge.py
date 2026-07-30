# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Convert arbitrary site inputs to numpy arrays
for AI inversion (PINN and Hybrid) modules.

Accepted inputs (passed through ``ensure_sites``):
    - filesystem path / glob / directory
    - ``EDIFile``, ``EDICollection``
    - ``Site``, ``Sites``
    - ``APISurvey``
    - iterables of any of the above
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "SiteObs1D",
    "SiteObs2D",
    "sites_to_obs_1d",
    "sites_to_obs_2d",
    "sites_to_features_1d",
    "sites_to_panel_2d",
    "sites_to_coords_3d",
    "obs_to_features_1d",
]

_COMP_IDX: dict = {
    "xy": (0, 1),
    "yx": (1, 0),
    "xx": (0, 0),
    "yy": (1, 1),
}


@dataclass
class SiteObs1D:
    r"""
    Observed 1-D MT/CSAMT data for one station.

    Parameters
    ----------
    name : str
        Station identifier.
    freq : ndarray (n_f,)
        Frequencies in Hz, sorted high to low.
    rho_obs : ndarray (n_f,)
        Apparent resistivity in Ohm.m.
    phase_obs : ndarray (n_f,)
        Impedance phase in degrees.
    """

    name: str
    freq: np.ndarray
    rho_obs: np.ndarray
    phase_obs: np.ndarray


def _normalize_sites(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    """Return a Sites object from any input."""
    from pycsamt.emtools._core import ensure_sites

    return ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )


def _extract_rho_phase(
    site: Any,
    comp: str,
) -> tuple[
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
]:
    """
    Return (freq, rho_c, phase_c) for one site.

    Returns three ``None`` values when any array is
    missing or when no finite values survive masking.
    """
    freq = site.freq
    rho = site.rho
    ph = site.phase
    if freq is None or rho is None or ph is None:
        return None, None, None

    freq = np.asarray(freq, dtype=float)
    rho = np.asarray(rho, dtype=float)
    ph = np.asarray(ph, dtype=float)

    if rho.ndim == 3:  # (n, 2, 2)
        i, j = _COMP_IDX[comp]
        rho_c = rho[:, i, j]
        ph_c = ph[:, i, j]
    elif rho.ndim == 1:
        rho_c = rho
        ph_c = ph
    else:
        name = getattr(site, "name", "?")
        warnings.warn(
            f"Unexpected rho shape {rho.shape} for site {name!r}. Skipping.",
            UserWarning,
            stacklevel=5,
        )
        return None, None, None

    # Sort high-to-low frequency
    order = np.argsort(freq)[::-1]
    freq = freq[order]
    rho_c = rho_c[order]
    ph_c = ph_c[order]

    valid = np.isfinite(rho_c) & np.isfinite(ph_c) & (rho_c > 0)
    if not valid.any():
        return None, None, None

    return freq[valid], rho_c[valid], ph_c[valid]


def _make_common_grid(
    all_freq: list[np.ndarray],
    n_freqs: int = 32,
    freq_min: float | None = None,
    freq_max: float | None = None,
) -> np.ndarray:
    """Build a log-spaced common frequency grid."""
    all_f = np.concatenate(all_freq)
    all_f = all_f[np.isfinite(all_f) & (all_f > 0)]
    fmin = freq_min if freq_min is not None else float(all_f.min())
    fmax = freq_max if freq_max is not None else float(all_f.max())
    return np.logspace(np.log10(fmin), np.log10(fmax), n_freqs)


def _interp_to_grid(
    freq: np.ndarray,
    rho_c: np.ndarray,
    ph_c: np.ndarray,
    freqs_grid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate rho/phase to a common log-freq grid."""
    from scipy.interpolate import interp1d

    n = len(freqs_grid)
    if len(freq) < 2:
        return np.full(n, np.nan), np.full(n, np.nan)

    log_f_src = np.log10(freq)
    log_f_dst = np.log10(freqs_grid)
    log_rho_src = np.log10(np.maximum(rho_c, 1e-12))

    fn_rho = interp1d(
        log_f_src,
        log_rho_src,
        kind="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    fn_ph = interp1d(
        log_f_src,
        ph_c,
        kind="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    log_rho_g = fn_rho(log_f_dst)
    ph_g = fn_ph(log_f_dst)

    # Mask frequencies outside station range
    f_lo, f_hi = freq.min(), freq.max()
    oob = (freqs_grid < f_lo) | (freqs_grid > f_hi)
    log_rho_g[oob] = np.nan
    ph_g[oob] = np.nan
    return log_rho_g, ph_g


def sites_to_obs_1d(
    sites: Any,
    *,
    comp: str = "xy",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> list[SiteObs1D]:
    """
    Convert any site input to observed-data containers.

    Parameters
    ----------
    sites : Any
        Path, glob, ``EDIFile``, ``EDICollection``,
        ``Site``, ``Sites``, ``APISurvey``, or
        iterable of any of the above.
    comp : {'xy', 'yx', 'xx', 'yy'}, default 'xy'
        Impedance tensor component to extract.
    recursive, on_dup, strict, verbose
        Forwarded to ``ensure_sites``.

    Returns
    -------
    list of SiteObs1D
        One entry per valid station.

    Raises
    ------
    ValueError
        If *comp* is invalid or no valid site data
        survives extraction.
    """
    if isinstance(sites, list) and sites and isinstance(sites[0], SiteObs1D):
        return sites
    if comp not in _COMP_IDX:
        raise ValueError(
            f"comp must be one of {list(_COMP_IDX)}; got {comp!r}."
        )
    S = _normalize_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    out: list[SiteObs1D] = []
    for site in S:
        name = getattr(site, "name", str(id(site)))
        freq, rho_c, ph_c = _extract_rho_phase(site, comp)
        if freq is None:
            if verbose > 0:
                warnings.warn(
                    f"No valid data for {name!r}. Skipping.",
                    UserWarning,
                    stacklevel=2,
                )
            continue
        out.append(
            SiteObs1D(
                name=name,
                freq=freq,
                rho_obs=rho_c,
                phase_obs=ph_c,
            )
        )
    if not out:
        raise ValueError(
            "sites_to_obs_1d: no valid sites after "
            "extraction.  Check comp= and data quality."
        )
    return out


def sites_to_features_1d(
    sites: Any,
    *,
    comp: str = "xy",
    n_freqs: int = 32,
    freq_min: float | None = None,
    freq_max: float | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Convert any site input to an ML feature matrix.

    The feature layout matches
    :meth:`~pycsamt.forward.em1d.ForwardResponse.to_array`
    so the result is directly compatible with a fitted
    :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`::

        X = [log10(rho_f1), ..., log10(rho_fn), phase_f1, ..., phase_fn]

    Parameters
    ----------
    sites : Any
        Same accepted types as :func:`sites_to_obs_1d`.
    comp : str, default ``'xy'``
        Tensor component.
    n_freqs : int, default 32
        Frequencies in the common interpolation grid.
    freq_min, freq_max : float or None
        Override the auto-detected frequency range.

    Returns
    -------
    X : ndarray, shape (n_stations, 2*n_freqs)
        Block-format feature matrix (rho block first,
        then phase block).
    freqs : ndarray, shape (n_freqs,)
        Common frequency grid in Hz.
    names : list of str
        Station names in the same row order as ``X``.
    """
    obs = sites_to_obs_1d(
        sites,
        comp=comp,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    freqs_grid = _make_common_grid(
        [o.freq for o in obs],
        n_freqs=n_freqs,
        freq_min=freq_min,
        freq_max=freq_max,
    )
    rows: list[np.ndarray] = []
    names: list[str] = []
    for o in obs:
        log_rho_g, ph_g = _interp_to_grid(
            o.freq, o.rho_obs, o.phase_obs, freqs_grid
        )
        # Block layout: rho block then phase block
        feat = np.empty(2 * n_freqs, dtype=np.float32)
        feat[:n_freqs] = log_rho_g.astype(np.float32)
        feat[n_freqs:] = ph_g.astype(np.float32)
        rows.append(feat)
        names.append(o.name)

    X = np.vstack(rows)  # (n_stations, 2*n_freqs)
    return X, freqs_grid, names


def obs_to_features_1d(
    obs: list,
    n_freqs: int = 32,
    *,
    freq_min: float | None = None,
    freq_max: float | None = None,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Build a feature matrix from a list of obs objects.

    Works with both :class:`SiteObs1D` (uses
    ``rho_obs``, ``phase_obs``) and
    :class:`SiteObs2D` (uses ``rho_te``,
    ``phase_te``).  Skips ``ensure_sites``, so it
    is safe to call with pre-extracted obs lists.

    Parameters
    ----------
    obs : list of SiteObs1D or SiteObs2D
    n_freqs : int, default 32

    Returns
    -------
    X : ndarray (n_stations, 2*n_freqs)
    freqs : ndarray (n_freqs,)
    names : list of str
    """
    freqs_grid = _make_common_grid(
        [o.freq for o in obs],
        n_freqs=n_freqs,
        freq_min=freq_min,
        freq_max=freq_max,
    )
    rows: list[np.ndarray] = []
    names: list[str] = []
    for o in obs:
        # Support SiteObs1D (rho_obs) and SiteObs2D
        # (rho_te).
        rho_src = getattr(
            o,
            "rho_obs",
            getattr(o, "rho_te", None),
        )
        ph_src = getattr(
            o,
            "phase_obs",
            getattr(o, "phase_te", None),
        )
        if rho_src is None or ph_src is None:
            continue
        lr_g, ph_g = _interp_to_grid(o.freq, rho_src, ph_src, freqs_grid)
        # Fill NaN with fallback constants
        lr_g = np.where(np.isfinite(lr_g), lr_g, 2.0)
        ph_g = np.where(np.isfinite(ph_g), ph_g, 45.0)
        feat = np.empty(2 * n_freqs, dtype=np.float32)
        feat[:n_freqs] = lr_g.astype(np.float32)
        feat[n_freqs:] = ph_g.astype(np.float32)
        rows.append(feat)
        names.append(o.name)

    if not rows:
        raise ValueError("obs_to_features_1d: no valid obs found.")
    return np.vstack(rows), freqs_grid, names


# ── 2-D helpers ──


@dataclass
class SiteObs2D:
    r"""
    Observed 2-D MT profile data for one station.

    Parameters
    ----------
    name : str
        Station identifier.
    freq : ndarray (n_f,)
        Frequencies in Hz, sorted high to low.
    rho_te : ndarray (n_f,)
        TE apparent resistivity in Ohm.m
        (Zxy convention).
    phase_te : ndarray (n_f,)
        TE phase in degrees.
    rho_tm : ndarray (n_f,)
        TM apparent resistivity in Ohm.m
        (Zyx convention).
    phase_tm : ndarray (n_f,)
        TM phase magnitude in degrees
        (absolute value stored).
    """

    name: str
    freq: np.ndarray
    rho_te: np.ndarray
    phase_te: np.ndarray
    rho_tm: np.ndarray
    phase_tm: np.ndarray


def _extract_rho_phase_2d(
    site: Any,
    comp_te: str,
    comp_tm: str,
) -> tuple[
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
]:
    """
    Return (freq, rho_te, ph_te, rho_tm, ph_tm)
    from a site object.

    TM phase is returned as absolute value.
    Returns five ``None`` values on any failure.
    """
    freq = site.freq
    rho = site.rho
    ph = site.phase
    if freq is None or rho is None or ph is None:
        return None, None, None, None, None

    freq = np.asarray(freq, dtype=float)
    rho = np.asarray(rho, dtype=float)
    ph = np.asarray(ph, dtype=float)

    order = np.argsort(freq)[::-1]
    freq = freq[order]

    if rho.ndim == 3:
        i_te, j_te = _COMP_IDX[comp_te]
        i_tm, j_tm = _COMP_IDX[comp_tm]
        rho_te = rho[:, i_te, j_te][order]
        ph_te = ph[:, i_te, j_te][order]
        rho_tm = rho[:, i_tm, j_tm][order]
        ph_tm = ph[:, i_tm, j_tm][order]
    elif rho.ndim == 1:
        rho_te = rho[order]
        ph_te = ph[order]
        rho_tm = rho_te.copy()
        ph_tm = np.abs(ph_te)
    else:
        return None, None, None, None, None

    valid = np.isfinite(rho_te) & np.isfinite(ph_te) & (rho_te > 0)
    valid_tm = np.isfinite(rho_tm) & np.isfinite(ph_tm) & (rho_tm > 0)
    # Use TE-valid mask; fill TM where TM missing
    if not valid.any():
        valid = valid_tm
    if not valid.any():
        return None, None, None, None, None

    rho_tm = np.where(valid_tm, rho_tm, rho_te)
    ph_tm = np.where(valid_tm, np.abs(ph_tm), np.abs(ph_te))

    return (
        freq[valid],
        rho_te[valid],
        ph_te[valid],
        rho_tm[valid],
        ph_tm[valid],
    )


def sites_to_obs_2d(
    sites: Any,
    *,
    comp_te: str = "xy",
    comp_tm: str = "yx",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> list[SiteObs2D]:
    """
    Convert any site input to 2-D observation containers.

    Parameters
    ----------
    sites : Any
        Path, glob, ``EDIFile``, ``EDICollection``,
        ``Site``, ``Sites``, ``APISurvey``, or
        iterable of any of the above.
    comp_te : str, default ``'xy'``
        Tensor component treated as TE mode.
    comp_tm : str, default ``'yx'``
        Tensor component treated as TM mode.
    recursive, on_dup, strict, verbose
        Forwarded to ``ensure_sites``.

    Returns
    -------
    list of SiteObs2D
        One entry per valid station, same order as
        input.

    Raises
    ------
    ValueError
        If no valid station data is found.
    """
    if isinstance(sites, list) and sites and isinstance(sites[0], SiteObs2D):
        return sites
    for comp in (comp_te, comp_tm):
        if comp not in _COMP_IDX:
            raise ValueError(
                f"comp must be one of {list(_COMP_IDX)}; got {comp!r}."
            )
    S = _normalize_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    out: list[SiteObs2D] = []
    for site in S:
        name = getattr(site, "name", str(id(site)))
        result = _extract_rho_phase_2d(site, comp_te, comp_tm)
        freq, rho_te, ph_te, rho_tm, ph_tm = result
        if freq is None:
            if verbose > 0:
                warnings.warn(
                    f"No valid data for {name!r}. Skipping.",
                    UserWarning,
                    stacklevel=2,
                )
            continue
        out.append(
            SiteObs2D(
                name=name,
                freq=freq,
                rho_te=rho_te,
                phase_te=ph_te,
                rho_tm=rho_tm,
                phase_tm=ph_tm,
            )
        )
    if not out:
        raise ValueError(
            "sites_to_obs_2d: no valid stations after "
            "extraction.  Check comp_te/comp_tm and "
            "data quality."
        )
    return out


def sites_to_panel_2d(
    sites: Any,
    *,
    n_freqs: int = 32,
    n_components: int = 4,
    comp_te: str = "xy",
    comp_tm: str = "yx",
    freq_min: float | None = None,
    freq_max: float | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Build an ``EMInverter2D`` input panel from sites.

    Returns
    -------
    X_panel : ndarray
        Shape ``(1, n_components, n_freqs, n_stations)``.
        Channel layout (for ``n_components=4``):
        ``[log10(rho_te), phase_te,
           log10(rho_tm), phase_tm]``.
        For ``n_components=2``: TE channels only.
    freqs : ndarray (n_freqs,)
        Common frequency grid in Hz (high to low).
    names : list of str
        Station names in the same column order as
        ``X_panel``.
    """
    # Accept pre-built SiteObs2D list directly
    if isinstance(sites, list) and sites and isinstance(sites[0], SiteObs2D):
        obs: list[SiteObs2D] = sites
    else:
        obs = sites_to_obs_2d(
            sites,
            comp_te=comp_te,
            comp_tm=comp_tm,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    freqs_grid = _make_common_grid(
        [o.freq for o in obs],
        n_freqs=n_freqs,
        freq_min=freq_min,
        freq_max=freq_max,
    )
    n_sta = len(obs)
    n_ch = min(n_components, 4)
    panel = np.full(
        (1, n_ch, n_freqs, n_sta),
        np.nan,
        dtype=np.float32,
    )
    names: list[str] = []
    for k, o in enumerate(obs):
        log_rho_te, ph_te = _interp_to_grid(
            o.freq, o.rho_te, o.phase_te, freqs_grid
        )
        log_rho_tm, ph_tm = _interp_to_grid(
            o.freq, o.rho_tm, o.phase_tm, freqs_grid
        )
        # freqs_grid is low-to-high from logspace;
        # panel freq axis = high-to-low
        lr_te = log_rho_te[::-1].astype(np.float32)
        p_te = ph_te[::-1].astype(np.float32)
        lr_tm = log_rho_tm[::-1].astype(np.float32)
        p_tm = ph_tm[::-1].astype(np.float32)
        channels = [lr_te, p_te, lr_tm, p_tm]
        for c in range(n_ch):
            panel[0, c, :, k] = channels[c]
        names.append(o.name)
    return panel, freqs_grid[::-1], names


def sites_to_coords_3d(
    sites: Any,
    *,
    station_spacing: float = 500.0,
    recursive: bool = True,
    on_dup: str = "replace",
    verbose: int = 0,
) -> np.ndarray:
    """
    Extract station (x, y) positions in metres.

    Tries ``site.coords`` (lat, lon); falls back to
    a uniform grid when coordinates are missing or
    all zero.

    Parameters
    ----------
    sites : Any
        Same accepted types as :func:`sites_to_obs_1d`.
    station_spacing : float, default 500.0
        Spacing in metres used when coordinates are
        unavailable (uniform-grid fallback).
    recursive, on_dup, verbose
        Forwarded to ``ensure_sites``.

    Returns
    -------
    xy : ndarray (n_stations, 2)
        Station (x, y) positions in metres.
    """
    S = _normalize_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        verbose=verbose,
    )
    station_list = list(S)
    n = len(station_list)

    lats = np.full(n, np.nan)
    lons = np.full(n, np.nan)

    for i, site in enumerate(station_list):
        c = getattr(site, "coords", None)
        if c is not None and len(c) >= 2:
            try:
                lats[i] = float(c[0])
                lons[i] = float(c[1])
            except (TypeError, ValueError):
                pass

    # Use flat-earth projection when coords span
    # more than 1 m in either direction (~1e-5 deg)
    lat_ok = np.isfinite(lats)
    lon_ok = np.isfinite(lons)
    all_ok = lat_ok & lon_ok
    lat_rng = (
        float(lats[all_ok].max() - lats[all_ok].min()) if all_ok.any() else 0.0
    )
    lon_rng = (
        float(lons[all_ok].max() - lons[all_ok].min()) if all_ok.any() else 0.0
    )
    _DEG_THRESH = 1e-5  # ~1 m

    if all_ok.sum() >= 2 and (lat_rng > _DEG_THRESH or lon_rng > _DEG_THRESH):
        lat_ref = float(np.nanmean(lats))
        lon_ref = float(np.nanmean(lons))
        cos_lat = float(np.cos(np.radians(lat_ref)))
        x = (lons - lon_ref) * cos_lat * 111320.0
        y = (lats - lat_ref) * 111320.0
        # Fill any remaining NaN with 0.0
        x = np.where(np.isfinite(x), x, 0.0)
        y = np.where(np.isfinite(y), y, 0.0)
        return np.column_stack([x, y])

    # Fallback: row-major uniform grid
    side = int(np.ceil(np.sqrt(n)))
    sp = float(station_spacing)
    xs = np.array([(i % side) * sp for i in range(n)])
    ys = np.array([(i // side) * sp for i in range(n)])
    return np.column_stack([xs, ys])
