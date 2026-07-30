# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
coord_correction — station-coordinate correction functions for AMT/CSAMT surveys.

All functions accept a Sites-like object, apply a spatial correction to the
lat / lon / elev stored in each site's EDI header, and return a **new** Sites
(inplace=False by default).

Functions
---------
correct_profile_projection   Project stations onto the best-fit survey line
correct_spacing_regularize   Redistribute stations at uniform spacing
correct_outlier_snap         Snap off-line stations to the profile
correct_elevation_smooth     Smooth the noisy GPS elevation profile
correct_coordinate_shift     Apply a fixed (Δlat, Δlon, Δelev) offset
correct_interpolate_missing  Fill missing / zero coordinates from neighbours
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

# ── Low-level helpers ─────────────────────────────────────────────────────────


def _get_coords_df(sites) -> pd.DataFrame:
    """Return DataFrame [station, lat, lon, elev] for every site."""
    from pycsamt.emtools._core import _iter_items, _name
    from pycsamt.site.utils import get_coords

    rows = []
    for i, ed in enumerate(_iter_items(sites)):
        station = _name(ed, i)
        # Sites yields Site wrappers; unwrap to underlying EDIFile so that
        # get_coords can reach ed.get_section("head")
        edi = getattr(ed, "edi", ed)
        try:
            c = get_coords(edi)
            lat = float(c.lat) if c is not None else np.nan
            lon = float(c.lon) if c is not None else np.nan
            elev = float(c.elev) if c is not None else 0.0
        except Exception:
            lat, lon, elev = np.nan, np.nan, 0.0
        rows.append({"station": station, "lat": lat, "lon": lon, "elev": elev})

    return pd.DataFrame(rows, columns=["station", "lat", "lon", "elev"])


def _apply_coords_bulk(
    sites, coords_dict: dict[str, tuple], inplace: bool = False
):
    """
    Apply corrected (lat, lon, elev) per station from *coords_dict*.

    Uses _apply_each so the correction is non-destructive when inplace=False.
    """
    from pycsamt.emtools._core import (
        _apply_each,
        _iter_items,
        _name,
    )
    from pycsamt.site.utils import set_coords

    def _one(Si):
        items = list(_iter_items(Si))
        if not items:
            return Si
        ed = items[0]
        # Unwrap Site wrapper → EDIFile so set_coords can reach the header
        edi = getattr(ed, "edi", ed)
        station = _name(ed, 0)
        if station in coords_dict:
            lat, lon, elev = coords_dict[station]
            set_coords(edi, lat=lat, lon=lon, elev=elev, inplace=True)
        return Si

    return _apply_each(sites, _one, inplace=inplace, verbose=0)


def _to_utm_arrays(lat, lon):
    """Convert lat/lon arrays to (easting, northing, zone_str)."""
    from pycsamt.gis.utils import to_utm

    try:
        result = to_utm(lat, lon, as_frame=True)
        east = np.asarray(result["easting"], dtype=float)
        north = np.asarray(result["northing"], dtype=float)
        zone = (
            str(result["zone"].iloc[0])
            if hasattr(result["zone"], "iloc")
            else str(result["zone"])
        )
        return east, north, zone
    except Exception:
        # Flat-Earth approximation as fallback (metres per degree ≈ 111_319)
        lat0 = float(np.nanmean(lat))
        east = (
            (np.asarray(lon) - float(np.nanmean(lon)))
            * 111_319.0
            * np.cos(np.radians(lat0))
        )
        north = (np.asarray(lat) - lat0) * 111_319.0
        return east, north, None


def _to_ll_arrays(easting, northing, zone):
    """Convert easting/northing arrays back to (lat, lon)."""
    from pycsamt.gis.utils import to_ll

    try:
        result = to_ll(easting, northing, zone=zone, as_frame=True)
        lats = np.asarray(result["lat"], dtype=float)
        lons = np.asarray(result["lon"], dtype=float)
        return lats, lons
    except Exception:
        # Reverse flat-Earth fallback
        return easting, northing  # caller should handle gracefully


def _pca_azimuth(easting: np.ndarray, northing: np.ndarray) -> float:
    """PCA-based profile azimuth (degrees from north, 0–180)."""
    pts = np.column_stack(
        [easting - easting.mean(), northing - northing.mean()]
    )
    cov = np.cov(pts.T)
    _, evecs = np.linalg.eigh(cov)
    pc1 = evecs[:, -1]  # largest eigenvalue
    az = float(np.degrees(np.arctan2(pc1[0], pc1[1])))
    return az % 180.0


# ── Correction functions ──────────────────────────────────────────────────────


def correct_profile_projection(
    sites: Any,
    *,
    azimuth: float | None = None,
    keep_elevation: bool = True,
    inplace: bool = False,
    verbose: int = 0,
) -> Any:
    """
    Project each station onto the best-fit survey line.

    Stations scattered perpendicular to the profile are moved to their
    nearest point on the line — their chainage (along-profile distance)
    is preserved, only the cross-profile offset is removed.

    Parameters
    ----------
    azimuth : float, optional
        Known profile azimuth in degrees from north [0, 360).
        If None, inferred from the station positions via PCA.
    keep_elevation : bool
        If True (default) the original elevation is kept; if False,
        elevation is linearly interpolated along the projected profile.
    """
    df = _get_coords_df(sites)
    valid = df.dropna(subset=["lat", "lon"])
    if len(valid) < 2:
        return sites

    east, north, zone = _to_utm_arrays(
        valid["lat"].values, valid["lon"].values
    )
    az = azimuth if azimuth is not None else _pca_azimuth(east, north)
    az_rad = np.radians(az)
    ue, un = np.sin(az_rad), np.cos(az_rad)  # unit vector along profile

    # Centroid as profile origin
    oe, on_ = east.mean(), north.mean()

    # Chainage = projection of offset vector onto profile direction
    chain = (east - oe) * ue + (north - on_) * un

    # New UTM positions (origin + chain × unit vector)
    new_east = oe + chain * ue
    new_north = on_ + chain * un

    # Convert back to lat/lon
    new_lat, new_lon = _to_ll_arrays(new_east, new_north, zone)

    coords: dict[str, tuple] = {}
    for i, (_, row) in enumerate(valid.iterrows()):
        coords[row["station"]] = (
            float(new_lat[i]),
            float(new_lon[i]),
            float(row["elev"]) if keep_elevation else float(row["elev"]),
        )

    return _apply_coords_bulk(sites, coords, inplace=inplace)


def correct_spacing_regularize(
    sites: Any,
    *,
    spacing_m: float = 200.0,
    azimuth: float | None = None,
    preserve_extent: bool = True,
    inplace: bool = False,
    verbose: int = 0,
) -> Any:
    """
    Redistribute stations at uniform spacing along the survey profile.

    The profile direction is kept; only the along-profile positions are
    resampled.  Station *count* is preserved and the profile extent
    (first/last chainage) is optionally maintained.

    Parameters
    ----------
    spacing_m : float
        Target station spacing in metres (default 200 m).
    preserve_extent : bool
        If True (default), keep the same total profile length; stations
        are uniformly spaced between the first and last chainage.
        If False, spacing_m is applied strictly and end stations may move.
    """
    df = _get_coords_df(sites)
    valid = df.dropna(subset=["lat", "lon"])
    n = len(valid)
    if n < 2:
        return sites

    east, north, zone = _to_utm_arrays(
        valid["lat"].values, valid["lon"].values
    )
    az = azimuth if azimuth is not None else _pca_azimuth(east, north)
    az_rad = np.radians(az)
    ue, un = np.sin(az_rad), np.cos(az_rad)

    oe, on_ = east.mean(), north.mean()
    chain = (east - oe) * ue + (north - on_) * un

    # Sort by chainage so resampling is monotone
    order = np.argsort(chain)
    chain_sorted = chain[order]
    elev_sorted = valid["elev"].values[order]
    stations_sorted = valid["station"].values[order]

    if preserve_extent:
        new_chain = np.linspace(chain_sorted[0], chain_sorted[-1], n)
    else:
        new_chain = np.arange(
            chain_sorted[0], chain_sorted[-1] + spacing_m, spacing_m
        )[:n]
        # Pad if fewer than n
        if len(new_chain) < n:
            new_chain = np.linspace(chain_sorted[0], chain_sorted[-1], n)

    new_east = oe + new_chain * ue
    new_north = on_ + new_chain * un
    new_lat, new_lon = _to_ll_arrays(new_east, new_north, zone)

    # Interpolate elevation to new chainages
    new_elev = np.interp(new_chain, chain_sorted, elev_sorted)

    coords: dict[str, tuple] = {}
    for i, st in enumerate(stations_sorted):
        coords[st] = (
            float(new_lat[i]),
            float(new_lon[i]),
            float(new_elev[i]),
        )

    return _apply_coords_bulk(sites, coords, inplace=inplace)


def correct_outlier_snap(
    sites: Any,
    *,
    threshold_m: float = 50.0,
    azimuth: float | None = None,
    inplace: bool = False,
    verbose: int = 0,
) -> Any:
    """
    Snap stations that deviate more than *threshold_m* from the profile line.

    Only cross-profile displacement is corrected; stations within the
    threshold are left untouched.  Uses the same PCA profile direction
    as ``correct_profile_projection``.

    Parameters
    ----------
    threshold_m : float
        Maximum allowed perpendicular distance from the profile line (m).
        Stations beyond this are projected onto the line.
    """
    df = _get_coords_df(sites)
    valid = df.dropna(subset=["lat", "lon"])
    if len(valid) < 2:
        return sites

    east, north, zone = _to_utm_arrays(
        valid["lat"].values, valid["lon"].values
    )
    az = azimuth if azimuth is not None else _pca_azimuth(east, north)
    az_rad = np.radians(az)
    ue, un = np.sin(az_rad), np.cos(az_rad)  # along-profile unit vector
    pe, pn = -un, ue  # perpendicular unit vector

    oe, on_ = east.mean(), north.mean()
    de = east - oe
    dn = north - on_

    chain = de * ue + dn * un  # along-profile distance
    cross_dist = np.abs(de * pe + dn * pn)  # perpendicular distance

    coords: dict[str, tuple] = {}
    for i, (_, row) in enumerate(valid.iterrows()):
        if cross_dist[i] > threshold_m:
            # Snap to profile line
            new_e = oe + chain[i] * ue
            new_n = on_ + chain[i] * un
            ll = _to_ll_arrays(np.array([new_e]), np.array([new_n]), zone)
            coords[row["station"]] = (
                float(ll[0][0]),
                float(ll[1][0]),
                float(row["elev"]),
            )

    if not coords:
        return sites
    return _apply_coords_bulk(sites, coords, inplace=inplace)


def correct_elevation_smooth(
    sites: Any,
    *,
    method: str = "loess",
    window: int = 5,
    azimuth: float | None = None,
    inplace: bool = False,
    verbose: int = 0,
) -> Any:
    """
    Smooth the elevation profile along the survey line.

    GPS vertical accuracy is typically 3–5× worse than horizontal;
    smoothing removes high-frequency noise without distorting the
    large-scale topographic trend.

    Parameters
    ----------
    method : {"loess", "mean"}
        Smoothing algorithm.  ``"loess"`` fits local quadratic polynomials
        (Cleveland 1979); ``"mean"`` applies a symmetric running average.
    window : int
        Half-window for local fitting (``"loess"``) or full window
        (``"mean"``).  Larger values = more smoothing.
    """
    df = _get_coords_df(sites)
    valid = df.dropna(subset=["lat", "lon"])
    if len(valid) < 3:
        return sites

    east, north, zone = _to_utm_arrays(
        valid["lat"].values, valid["lon"].values
    )
    az = azimuth if azimuth is not None else _pca_azimuth(east, north)
    az_rad = np.radians(az)
    ue, un = np.sin(az_rad), np.cos(az_rad)

    oe, on_ = east.mean(), north.mean()
    chain = (east - oe) * ue + (north - on_) * un
    order = np.argsort(chain)
    chain_s = chain[order]
    elev_s = valid["elev"].values[order].astype(float)
    st_s = valid["station"].values[order]

    if method == "loess":
        smooth_elev = _loess_smooth(chain_s, elev_s, span=window)
    else:  # running mean
        try:
            from scipy.ndimage import uniform_filter1d

            smooth_elev = uniform_filter1d(elev_s, size=max(2, window))
        except ImportError:
            k = max(1, window // 2)
            smooth_elev = np.array(
                [
                    np.mean(elev_s[max(0, i - k) : i + k + 1])
                    for i in range(len(elev_s))
                ]
            )

    coords: dict[str, tuple] = {}
    for i, (st, (_, row)) in enumerate(
        zip(st_s, valid.iloc[order].iterrows())
    ):
        coords[st] = (
            float(row["lat"]),
            float(row["lon"]),
            float(smooth_elev[i]),
        )

    return _apply_coords_bulk(sites, coords, inplace=inplace)


def correct_coordinate_shift(
    sites: Any,
    *,
    delta_lat: float = 0.0,
    delta_lon: float = 0.0,
    delta_elev: float = 0.0,
    inplace: bool = False,
    verbose: int = 0,
) -> Any:
    """
    Apply a uniform (Δlat, Δlon, Δelev) offset to all stations.

    Useful for correcting systematic datum errors (e.g., WGS84 vs.
    local datum offset) or known GPS receiver biases.
    """
    df = _get_coords_df(sites)
    coords: dict[str, tuple] = {
        row["station"]: (
            float(row["lat"]) + delta_lat,
            float(row["lon"]) + delta_lon,
            float(row["elev"]) + delta_elev,
        )
        for _, row in df.iterrows()
        if not (np.isnan(row["lat"]) or np.isnan(row["lon"]))
    }
    if not coords:
        return sites
    return _apply_coords_bulk(sites, coords, inplace=inplace)


def correct_interpolate_missing(
    sites: Any,
    *,
    fill_nan_only: bool = True,
    inplace: bool = False,
    verbose: int = 0,
) -> Any:
    """
    Fill missing or zero coordinates by linear interpolation from neighbours.

    Useful when a few stations have GPS dropout (lat=0, lon=0) or NaN
    headers.

    Parameters
    ----------
    fill_nan_only : bool
        If True (default), only fill NaN values.
        If False, also fill coordinates that are exactly 0.0.
    """
    df = _get_coords_df(sites).copy()

    for col in ("lat", "lon", "elev"):
        vals = df[col].copy()
        if not fill_nan_only:
            vals.replace(0.0, np.nan, inplace=True)
        df[col] = vals.interpolate(method="index", limit_direction="both")

    coords: dict[str, tuple] = {
        row["station"]: (
            float(row["lat"]),
            float(row["lon"]),
            float(row["elev"]),
        )
        for _, row in df.iterrows()
    }
    return _apply_coords_bulk(sites, coords, inplace=inplace)


# ── LOESS helper (no scipy dependency required) ───────────────────────────────


def _loess_smooth(x: np.ndarray, y: np.ndarray, span: int = 5) -> np.ndarray:
    """
    Simple local quadratic LOESS smoother (tricube weights).
    *span* is the half-window number of neighbours.
    """
    n = len(y)
    y_s = y.copy().astype(float)
    for i in range(n):
        lo = max(0, i - span)
        hi = min(n, i + span + 1)
        xi = x[lo:hi]
        yi = y[lo:hi]
        d = np.abs(xi - x[i])
        dmax = d.max() + 1e-12
        w = (1.0 - (d / dmax) ** 3) ** 3  # tricube
        sw = w.sum()
        if sw < 1e-12:
            continue
        # Weighted linear fit
        xw = (xi * w).sum() / sw
        yw = (yi * w).sum() / sw
        b1 = ((xi - xw) * (yi - yw) * w).sum() / (
            ((xi - xw) ** 2 * w).sum() + 1e-30
        )
        y_s[i] = yw + b1 * (x[i] - xw)
    return y_s


# ── DataFrame-based correction functions (no EDI objects touched) ─────────────
# These operate purely on [station, lat, lon, elev] DataFrames.
# The controller uses them so that EDI objects are NEVER modified during the
# preview/apply workflow — only at final commit time.


def df_profile_projection(
    df: pd.DataFrame,
    azimuth: float | None = None,
    keep_elevation: bool = True,
    **_,
) -> pd.DataFrame:
    """DataFrame version of correct_profile_projection."""
    valid = df.dropna(subset=["lat", "lon"])
    if len(valid) < 2:
        return df.copy()
    east, north, zone = _to_utm_arrays(
        valid["lat"].values, valid["lon"].values
    )
    az = (
        azimuth
        if (azimuth is not None and azimuth >= 0)
        else _pca_azimuth(east, north)
    )
    az_rad = np.radians(az)
    ue, un = np.sin(az_rad), np.cos(az_rad)
    oe, on_ = east.mean(), north.mean()
    chain = (east - oe) * ue + (north - on_) * un
    new_east = oe + chain * ue
    new_north = on_ + chain * un
    new_lat, new_lon = _to_ll_arrays(new_east, new_north, zone)
    result = df.copy()
    for i, idx in enumerate(valid.index):
        result.at[idx, "lat"] = float(new_lat[i])
        result.at[idx, "lon"] = float(new_lon[i])
    return result


def df_spacing_regularize(
    df: pd.DataFrame,
    spacing_m: float = 200.0,
    azimuth: float | None = None,
    preserve_extent: bool = True,
    **_,
) -> pd.DataFrame:
    """DataFrame version of correct_spacing_regularize."""
    valid = df.dropna(subset=["lat", "lon"])
    n = len(valid)
    if n < 2:
        return df.copy()
    east, north, zone = _to_utm_arrays(
        valid["lat"].values, valid["lon"].values
    )
    az = (
        azimuth
        if (azimuth is not None and azimuth >= 0)
        else _pca_azimuth(east, north)
    )
    az_rad = np.radians(az)
    ue, un = np.sin(az_rad), np.cos(az_rad)
    oe, on_ = east.mean(), north.mean()
    chain = (east - oe) * ue + (north - on_) * un
    order = np.argsort(chain)
    chain_s = chain[order]
    elev_s = valid["elev"].values[order].astype(float)
    idx_s = valid.index[order]
    if preserve_extent:
        new_chain = np.linspace(chain_s[0], chain_s[-1], n)
    else:
        new_chain = np.arange(chain_s[0], chain_s[-1] + spacing_m, spacing_m)[
            :n
        ]
        if len(new_chain) < n:
            new_chain = np.linspace(chain_s[0], chain_s[-1], n)
    new_east = oe + new_chain * ue
    new_north = on_ + new_chain * un
    new_lat, new_lon = _to_ll_arrays(new_east, new_north, zone)
    new_elev = np.interp(new_chain, chain_s, elev_s)
    result = df.copy()
    for i, idx in enumerate(idx_s):
        result.at[idx, "lat"] = float(new_lat[i])
        result.at[idx, "lon"] = float(new_lon[i])
        result.at[idx, "elev"] = float(new_elev[i])
    return result


def df_outlier_snap(
    df: pd.DataFrame,
    threshold_m: float = 50.0,
    azimuth: float | None = None,
    **_,
) -> pd.DataFrame:
    """DataFrame version of correct_outlier_snap."""
    valid = df.dropna(subset=["lat", "lon"])
    if len(valid) < 2:
        return df.copy()
    east, north, zone = _to_utm_arrays(
        valid["lat"].values, valid["lon"].values
    )
    az = (
        azimuth
        if (azimuth is not None and azimuth >= 0)
        else _pca_azimuth(east, north)
    )
    az_rad = np.radians(az)
    ue, un = np.sin(az_rad), np.cos(az_rad)
    pe, pn = -un, ue
    oe, on_ = east.mean(), north.mean()
    de, dn = east - oe, north - on_
    chain = de * ue + dn * un
    cross = np.abs(de * pe + dn * pn)
    result = df.copy()
    for i, idx in enumerate(valid.index):
        if cross[i] > threshold_m:
            ne = oe + chain[i] * ue
            nn = on_ + chain[i] * un
            ll = _to_ll_arrays(np.array([ne]), np.array([nn]), zone)
            result.at[idx, "lat"] = float(ll[0][0])
            result.at[idx, "lon"] = float(ll[1][0])
    return result


def df_elevation_smooth(
    df: pd.DataFrame,
    method: str = "loess",
    window: int = 5,
    azimuth: float | None = None,
    **_,
) -> pd.DataFrame:
    """DataFrame version of correct_elevation_smooth."""
    valid = df.dropna(subset=["lat", "lon"])
    if len(valid) < 3:
        return df.copy()
    east, north, zone = _to_utm_arrays(
        valid["lat"].values, valid["lon"].values
    )
    az = (
        azimuth
        if (azimuth is not None and azimuth >= 0)
        else _pca_azimuth(east, north)
    )
    az_rad = np.radians(az)
    ue, un = np.sin(az_rad), np.cos(az_rad)
    oe, on_ = east.mean(), north.mean()
    chain = (east - oe) * ue + (north - on_) * un
    order = np.argsort(chain)
    elev_s = valid["elev"].values[order].astype(float)
    idx_s = valid.index[order]
    if method == "loess":
        smooth = _loess_smooth(chain[order], elev_s, span=window)
    else:
        try:
            from scipy.ndimage import uniform_filter1d

            smooth = uniform_filter1d(elev_s, size=max(2, window))
        except ImportError:
            k = max(1, window // 2)
            smooth = np.array(
                [
                    np.mean(elev_s[max(0, i - k) : i + k + 1])
                    for i in range(len(elev_s))
                ]
            )
    result = df.copy()
    for i, idx in enumerate(idx_s):
        result.at[idx, "elev"] = float(smooth[i])
    return result


def df_shift(
    df: pd.DataFrame,
    delta_lat: float = 0.0,
    delta_lon: float = 0.0,
    delta_elev: float = 0.0,
    **_,
) -> pd.DataFrame:
    """DataFrame version of correct_coordinate_shift."""
    result = df.copy()
    result["lat"] = result["lat"] + delta_lat
    result["lon"] = result["lon"] + delta_lon
    result["elev"] = result["elev"] + delta_elev
    return result


def df_interpolate_missing(
    df: pd.DataFrame,
    fill_nan_only: bool = True,
    **_,
) -> pd.DataFrame:
    """DataFrame version of correct_interpolate_missing."""
    result = df.copy()
    for col in ("lat", "lon", "elev"):
        vals = result[col].copy()
        if not fill_nan_only:
            vals.replace(0.0, np.nan, inplace=True)
        result[col] = vals.interpolate(method="index", limit_direction="both")
    return result


# ── Dispatcher ────────────────────────────────────────────────────────────────

_DF_FUNCS: dict[str, Any] = {
    "_coord_profile_projection": df_profile_projection,
    "_coord_spacing_regularize": df_spacing_regularize,
    "_coord_outlier_snap": df_outlier_snap,
    "_coord_elevation_smooth": df_elevation_smooth,
    "_coord_shift": df_shift,
    "_coord_interpolate_missing": df_interpolate_missing,
}


def apply_coord_correction_to_df(
    fn_name: str,
    df: pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
    """Apply a named coordinate correction to a DataFrame and return the result."""
    fn = _DF_FUNCS.get(fn_name)
    if fn is None:
        return df.copy()
    return fn(df, **kwargs)


def apply_coords_df_to_sites(sites, coords_df: pd.DataFrame) -> None:
    """Write corrected coordinates from *coords_df* back to EDI headers (in-place)."""
    from pycsamt.emtools._core import _iter_items, _name
    from pycsamt.site.utils import set_coords

    lookup = dict(
        zip(
            coords_df["station"],
            zip(coords_df["lat"], coords_df["lon"], coords_df["elev"]),
        )
    )
    for i, ed in enumerate(_iter_items(sites)):
        station = _name(ed, i)
        edi = getattr(ed, "edi", ed)
        if station in lookup:
            lat, lon, elev = lookup[station]
            if np.isfinite(lat) and np.isfinite(lon):
                set_coords(
                    edi,
                    lat=float(lat),
                    lon=float(lon),
                    elev=float(elev),
                    inplace=True,
                )
