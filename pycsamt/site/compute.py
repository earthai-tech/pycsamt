# -*- coding: utf-8 -*-

from __future__ import annotations

from typing import Any, Iterable, Tuple, Optional
import math
import numpy as np
import pandas as pd 

from ..constants import MU_0 
from .utils import (
    iter_edifiles, station_name, get_freq, get_z,
)


def _as_sites_iter(obj: Any) -> Iterable[tuple[str, Any]]:
    for ed in iter_edifiles(obj):
        yield station_name(ed), ed


def _rho_from_z(z: complex, f: float) -> float:
    if not np.isfinite(f) or f <= 0:
        return np.nan
    return (abs(z) ** 2) / (MU_0 * 2.0 * math.pi * f) # # H/m


def _rho_xy_yx(Z: np.ndarray, f: np.ndarray) -> tuple:
    zxy = Z[..., 0, 1]
    zyx = Z[..., 1, 0]
    rho_xy = np.array([_rho_from_z(z, ff) for z, ff
                       in zip(zxy, f)], float)
    rho_yx = np.array([_rho_from_z(z, ff) for z, ff
                       in zip(zyx, f)], float)
    return rho_xy, rho_yx


def _phase_deg(a: np.ndarray) -> np.ndarray:
    return np.degrees(np.angle(a))


def _nearest_idx(f: np.ndarray, fx: float) -> int:
    return int(np.argmin(np.abs(f - fx)))


def _interp(x: np.ndarray, y: np.ndarray, xq: float) -> float:
    try:
        return float(np.interp(xq, x, y))
    except Exception:
        return float("nan")


def _rotmat(theta_deg: float) -> np.ndarray:
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    return np.array([[c, s], [-s, c]], float)


def _rotate_tensor(Z: np.ndarray, theta: float) -> np.ndarray:
    R = _rotmat(theta)
    RT = R.T
    # Z shape (nf, 2, 2)
    return np.einsum("ab,nbc,cd->nad", R, Z, RT)


def _swift_cost(Z: np.ndarray, theta: float) -> float:
    Zr = _rotate_tensor(Z, theta)
    d = np.abs(Zr[:, 0, 0]) ** 2 + np.abs(Zr[:, 1, 1]) ** 2
    if d.size == 0:
        return float("inf")
    return float(np.nanmedian(d))


def _swift_theta(Z: np.ndarray) -> float:
    # grid search 0..179 deg (Swift-style)
    thetas = np.arange(0.0, 180.0, 1.0, float)
    costs = np.array([_swift_cost(Z, th) for th in thetas],
                     float)
    if not np.isfinite(costs).any():
        return float("nan")
    k = int(np.nanargmin(costs))
    return float(thetas[k])


def _phase_diff_theta(Z: np.ndarray) -> float:
    # fallback: pick quadrant via |Zxy| vs |Zyx|
    zxy = np.nanmedian(np.abs(Z[:, 0, 1]))
    zyx = np.nanmedian(np.abs(Z[:, 1, 0]))
    if not np.isfinite(zxy + zyx):
        return float("nan")
    return 0.0 if zxy >= zyx else 90.0


def strike_estimate(
    obj: Any,
    *,
    method: str = "swift",
) -> Any:
    """
    Single site -> float angle (deg).
    Multi-site -> pandas.DataFrame with columns:
    station, method, theta_deg.
    """
    multi = False
    rows = []
    for st, ed in _as_sites_iter(obj):
        Z = get_z(ed)
        f = get_freq(ed)
        if Z is None or f is None or Z.ndim != 3 or \
           Z.shape[1:] != (2, 2):
            ang = float("nan")
        else:
            m = (method or "swift").lower()
            if m == "swift" or m == "groom":
                ang = _swift_theta(Z)
            elif m == "phase_diff":
                ang = _phase_diff_theta(Z)
            else:
                ang = _swift_theta(Z)
        rows.append((st, method, ang))
        multi = True
    if len(rows) == 1 and not isinstance(obj, Iterable):
        return rows[0][2]
    if len(rows) == 1 and not multi:
        return rows[0][2]

    return pd.DataFrame(
        rows, columns=["station", "method", "theta_deg"]
    )


def res_at_freq(
    obj: Any,
    freq: float,
    *,
    how: str = "nearest",
) -> Any:
    """
    Single site -> dict with res_xy, res_yx, f_used.
    Multi-site -> DataFrame with station, res_xy, res_yx, f_used.
    """
    rows = []
    for st, ed in _as_sites_iter(obj):
        f = get_freq(ed)
        Z = get_z(ed)
        if f is None or Z is None:
            rows.append((st, np.nan, np.nan, np.nan))
            continue
        if how.lower().startswith("near"):
            k = _nearest_idx(f, float(freq))
            fx = float(f[k])
            zxy = Z[k, 0, 1]
            zyx = Z[k, 1, 0]
            rxy = _rho_from_z(zxy, fx)
            ryx = _rho_from_z(zyx, fx)
            rows.append((st, rxy, ryx, fx))
        else:
            rxy, ryx = _rho_xy_yx(Z, f)
            rxy = _interp(f, rxy, float(freq))
            ryx = _interp(f, ryx, float(freq))
            rows.append((st, rxy, ryx, float(freq)))
    if len(rows) == 1 and not isinstance(obj, Iterable):
        _, rxy, ryx, fx = rows[0]
        return {"res_xy": rxy, "res_yx": ryx, "f_used": fx}

    return pd.DataFrame(
        rows, columns=["station", "res_xy", "res_yx", "f_used"]
    )


def phase_slope(
    obj: Any,
    band: Tuple[float, float],
) -> Any:
    """
    Slope d(phase)/d(log10 f) per component in band.
    Single site -> dict with slope_xy, slope_yx.
    Multi-site -> DataFrame with station, slope_xy, slope_yx.
    """
    fmin, fmax = float(band[0]), float(band[1])
    rows = []
    for st, ed in _as_sites_iter(obj):
        f = get_freq(ed)
        Z = get_z(ed)
        if f is None or Z is None:
            rows.append((st, np.nan, np.nan))
            continue
        mask = (f >= min(fmin, fmax)) & (f <= max(fmin, fmax))
        if not np.any(mask):
            rows.append((st, np.nan, np.nan))
            continue
        xf = np.log10(f[mask])
        ph_xy = _phase_deg(Z[mask, 0, 1])
        ph_yx = _phase_deg(Z[mask, 1, 0])
        try:
            sx = float(np.polyfit(xf, ph_xy, 1)[0])
        except Exception:
            sx = float("nan")
        try:
            sy = float(np.polyfit(xf, ph_yx, 1)[0])
        except Exception:
            sy = float("nan")
        rows.append((st, sx, sy))
    if len(rows) == 1 and not isinstance(obj, Iterable):
        _, sx, sy = rows[0]
        return {"slope_xy": sx, "slope_yx": sy}

    return pd.DataFrame(
        rows, columns=["station", "slope_xy", "slope_yx"]
    )


def tipper_magnitude(
    obj: Any,
    *,
    per_freq: bool = False,
) -> Any:
    """
    If single site:
      - per_freq=False -> dict(mean, median, max)
      - per_freq=True  -> dict(freq, mag)
    If multi-site:
      - per_freq=False -> DataFrame(station, mean, median, max)
      - per_freq=True  -> DataFrame(station, freq, mag)
    """
    def _tip_arr(ed: Any) -> Optional[np.ndarray]:
        t = getattr(ed, "Tip", None)
        arr = getattr(t, "_tipper", None)
        if arr is None:
            # sometimes public name is 'tipper'
            arr = getattr(t, "tipper", None)
        if arr is None:
            return None
        a = np.asarray(arr)
        if a.ndim == 1:
            a = a[..., None]
        return a

    def _freq_arr(ed: Any) -> Optional[np.ndarray]:
        f = get_freq(ed)
        return None if f is None else np.asarray(f, float)

    rows = []
    long_rows = []
    for st, ed in _as_sites_iter(obj):
        f = _freq_arr(ed)
        T = _tip_arr(ed)
        if f is None or T is None:
            if per_freq:
                continue
            rows.append((st, np.nan, np.nan, np.nan))
            continue
        # norm across horizontal components
        mag = np.sqrt(np.sum(np.abs(T) ** 2, axis=-1))
        if per_freq:
            for ff, mm in zip(f, mag):
                long_rows.append((st, float(ff), float(mm)))
        else:
            rows.append(
                (st,
                 float(np.nanmean(mag)),
                 float(np.nanmedian(mag)),
                 float(np.nanmax(mag)))
            )

    # single site return dict
    if not isinstance(obj, Iterable):
        if per_freq:
            if not long_rows:
                return {"freq": np.array([]), "mag": np.array([])}
            _, f0, _ = long_rows[0]
            # rebuild arrays in order
            fa = np.array([r[1] for r in long_rows], float)
            ma = np.array([r[2] for r in long_rows], float)
            return {"freq": fa, "mag": ma}
        if not rows:
            return {"mean": np.nan, "median": np.nan,
                    "max": np.nan}
        _, m, md, mx = rows[0]
        return {"mean": m, "median": md, "max": mx}

    if per_freq:
        return pd.DataFrame(
            long_rows, columns=["station", "freq", "mag"]
        )
    return pd.DataFrame(
        rows, columns=["station", "mean", "median", "max"]
    )
