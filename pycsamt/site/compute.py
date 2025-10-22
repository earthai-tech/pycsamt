# -*- coding: utf-8 -*- 
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from typing import Any, Iterable, Tuple, Optional
import math
import numpy as np
import pandas as pd

from ..constants import MU_0
from .utils import (
    iter_edifiles,
    station_name,
    get_freq,
)

__all__= [ 
    
    "strike_estimate", 
    "res_at_freq", 
    "phase_slope", 
    "tipper_magnitude", 
    ]

def strike_estimate(
    obj: Any,
    *,
    method: str = "swift",
) -> Any:
    r"""
    Estimate a strike angle from impedance tensors.
    
    This computes a 2D geoelectric strike angle in degrees.  The
    routine supports either a single site or many sites.  For a
    single site a scalar angle is returned.  For multiple sites a
    ``pandas.DataFrame`` is returned with one row per site.
    
    Parameters
    ----------
    obj : Any
        A single EDI-like object (e.g. ``EDIFile``) or an iterable
        of EDI-like objects.  Each item must expose a ``.Z`` section
        of shape ``(n_freq, 2, 2)`` or be convertible to such.
    method : str, optional
        Strike method.  Allowed values are:
    
        - ``"swift"`` (default): grid search over 0..179 degrees that
          minimizes the diagonal power after rotation.
        - ``"groom"``: alias of ``"swift"`` in this lightweight mode.
        - ``"phase_diff"``: heuristic that returns ``0`` or ``90``
          degrees based on the relative magnitude of off-diagonals.
    
    Returns
    -------
    float or pandas.DataFrame
        If ``obj`` is a single site, returns a float angle in
        degrees within ``[0, 180)``.  If ``obj`` is iterable,
        returns a ``DataFrame`` with columns:
    
        ``station``, ``method``, ``theta_deg``.
    
    Notes
    -----
    The Swift-style criterion rotates the impedance tensor
    :math:`Z` by a test angle :math:`\\theta` and minimizes
    
    .. math::
    
       J(\\theta) = \lvert Z'_{xx} \rvert^2 +
                    \lvert Z'_{yy} \rvert^2 ,
    
    where :math:`Z'` is the rotated tensor.  The returned
    angle is the argmin over a 1 degree grid in 0..179.
    
    The ``"phase_diff"`` fallback returns ``0`` if
    median :math:`\\lvert Z_{xy} \rvert \ge
    \lvert Z_{yx} \rvert`, else ``90``.  It is intended for
    degraded or sparse data.
    
    This function does not alter data.  If you need
    deterministic behavior for incomplete arrays, consider
    preparing tensors with
    :func:`pycsamt.site.edit.fill_missing`.
    
    Examples
    --------
    Single site, Swift estimate:
    
    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site import compute as cmp, edit as ed
    >>> edf = EDIFile("S01.edi")  # doctest: +SKIP
    >>> edf = ed.fill_missing(edf, how="zero",
    ...                       components=("Z",), inplace=False)
    >>> ang = cmp.strike_estimate(edf, method="swift")  # doctest: +SKIP
    >>> 0.0 <= ang < 180.0  # doctest: +SKIP
    True
    
    Many sites, returning a DataFrame:
    
    >>> e1 = EDIFile("S01.edi")  # doctest: +SKIP
    >>> e2 = EDIFile("S02.edi")  # doctest: +SKIP
    >>> df = cmp.strike_estimate([e1, e2], method="phase_diff")
    ... # doctest: +SKIP
    >>> list(df.columns)  # doctest: +SKIP
    ['station', 'method', 'theta_deg']
    
    See Also
    --------
    pycsamt.site.edit.rotate
        Rotate site tensors by a user angle.
    pycsamt.site.compute.phase_slope
        Phase slope diagnostic over a frequency band.
    
    References
    ----------
    .. [1] Swift, C. M., 1967. A magnetotelluric investigation of an
           electrical conductivity anomaly in the southwestern
           United States. PhD thesis, MIT.
    .. [2] Groom, R. W., and R. C. Bailey, 1989. Decomposition of
           magnetotelluric impedance tensors in the presence of
           local three dimensional galvanic distortion. JGR.
    """

    rows = []
    for st, ed in _as_sites_iter(obj):
        Z = _get_z(ed)
        f = get_freq(ed)
        if Z is None or f is None or Z.ndim != 3 or \
           Z.shape[1:] != (2, 2):
            ang = float("nan")
        else:
            m = (method or "swift").lower()
            if m in {"swift", "groom"}:
                ang = _swift_theta(Z)
            elif m == "phase_diff":
                ang = _phase_diff_theta(Z)
            else:
                ang = _swift_theta(Z)
        rows.append((st, method, ang))

    if len(rows) == 1 and not isinstance(obj, Iterable):
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
    r"""
    Evaluate apparent resistivity at a target frequency.
    
    Computes apparent resistivity for the :math:`Z_{xy}` and
    :math:`Z_{yx}` components at a requested frequency.  Works with a
    single site or a collection.  For a single site, a dict is
    returned.  For multiple sites, a ``pandas.DataFrame`` is returned.
    
    Parameters
    ----------
    obj : Any
        A single EDI-like object (e.g. ``EDIFile``) or an iterable of
        such objects.  Each item must expose a ``.Z`` section of shape
        ``(n_freq, 2, 2)`` and a frequency vector.
    freq : float
        Query frequency in Hz.
    how : str, optional
        Selection mode:
    
        - ``"nearest"`` (default): choose the nearest available
          frequency in the site data and report that value.
        - ``"interp"``: linearly interpolate resistivity versus
          frequency using ``numpy.interp``.  Interpolation occurs on
          linear frequency, not log frequency.
    
    Returns
    -------
    dict or pandas.DataFrame
        If ``obj`` is a single site, returns a dictionary with keys
        ``"res_xy"``, ``"res_yx"``, ``"f_used"``.  If ``obj`` is
        iterable, returns a ``DataFrame`` with columns
        ``station``, ``res_xy``, ``res_yx``, ``f_used``.
    
    Notes
    -----
    Apparent resistivity :math:`\\rho_a` is computed as
    
    .. math::
    
       \rho_a = \frac{\lvert Z \rvert^2}
                       {\mu_0\,2\pi\\,f} ,
    
    where :math:`Z` is the complex impedance for the selected
    component, :math:`\\mu_0` is the magnetic permeability of free
    space, and :math:`f` is frequency in Hz.
    
    When ``how="interp"``, the function first computes
    :math:`\rho_a` at all native frequencies, then interpolates the
    result to the query frequency using linear interpolation in
    frequency.  If the frequency vector or impedance is missing,
    ``NaN`` values are returned.
    
    Examples
    --------
    Single site, nearest selection:
    
    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site import compute as cmp, edit as ed
    >>> edf = EDIFile("S01.edi")  # doctest: +SKIP
    >>> edf = ed.fill_missing(edf, how="zero",
    ...                       components=("Z",), inplace=False)
    >>> out = cmp.res_at_freq(edf, 150.0, how="nearest")
    ... # doctest: +SKIP
    >>> set(out.keys()) == {"res_xy", "res_yx", "f_used"}
    ... # doctest: +SKIP
    True
    
    Single site, interpolated:
    
    >>> out = cmp.res_at_freq(edf, 150.0, how="interp")
    ... # doctest: +SKIP
    >>> out["f_used"]  # doctest: +SKIP
    150.0
    
    Many sites, DataFrame:
    
    >>> e1 = EDIFile("S01.edi")  # doctest: +SKIP
    >>> e2 = EDIFile("S02.edi")  # doctest: +SKIP
    >>> df = cmp.res_at_freq([e1, e2], 1.0, how="interp")
    ... # doctest: +SKIP
    >>> list(df.columns)  # doctest: +SKIP
    ['station', 'res_xy', 'res_yx', 'f_used']
    
    See Also
    --------
    pycsamt.site.compute.strike_estimate
        Estimate 2D strike angle from Z.
    pycsamt.site.edit.select_freq
        Subset site data by frequency criteria.
    
    References
    ----------
    .. [1] Vozoff, K., 1991. The magnetotelluric method. In
           Electromagnetic methods in applied geophysics.
    .. [2] Simpson, F., and K. Bahr, 2005. Practical Magnetotellurics.
           Cambridge University Press.
    """

    rows = []
    fq = float(freq)
    for st, ed in _as_sites_iter(obj):
        f = get_freq(ed)
        Z = _get_z(ed)
        if f is None or Z is None:
            rows.append((st, np.nan, np.nan, np.nan))
            continue

        if how.lower().startswith("near"):
            k = _nearest_idx(f, fq)
            fx = float(f[k])
            zxy = Z[k, 0, 1]
            zyx = Z[k, 1, 0]
            rxy = _rho_from_z(zxy, fx)
            ryx = _rho_from_z(zyx, fx)
            rows.append((st, rxy, ryx, fx))
        else:
            rxy, ryx = _rho_xy_yx(Z, f)
            rxyv = _interp(f, rxy, fq)
            ryxv = _interp(f, ryx, fq)
            rows.append((st, rxyv, ryxv, fq))

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
    r"""
    Estimate phase slopes within a frequency band.
    
    For each site, this computes the least-squares slope of phase
    (degrees) versus :math:`\log_{10}(f)` over the requested band.
    Two slopes are reported, one for :math:`Z_{xy}` and one for
    :math:`Z_{yx}`.
    
    If a single site is provided, a dictionary is returned.  If an
    iterable of sites is provided, a ``pandas.DataFrame`` is
    returned with one row per station.
    
    Parameters
    ----------
    obj : Any
        A single EDI-like object (e.g. ``EDIFile``) or an iterable
        of such objects.
    band : tuple of float
        Inclusive frequency band as ``(fmin, fmax)`` in Hz.  The
        order does not matter; the function uses the numeric min and
        max.
    
    Returns
    -------
    dict or pandas.DataFrame
        Single site -> ``{"slope_xy": float, "slope_yx": float}``.
        Multi-site  -> DataFrame with columns
        ``["station", "slope_xy", "slope_yx"]``.
    
    Notes
    -----
    The phase series for each off-diagonal component is computed as
    
    .. math::
    
       \phi(f) = \operatorname{angle}(Z(f)) \times 180/\\pi ,
    
    then a straight line is fit
    
    .. math::
    
       \phi(f) \approx a\\,\\log_{10}(f) + b
    
    using ``numpy.polyfit(x, y, 1)`` where
    :math:`x=\\log_{10}(f)`.  The reported slope is :math:`a` in
    units of degrees per decade.
    
    Rows or sites with missing data in the band are reported as
    ``NaN``.  The function does not unwrap phase.
    
    Examples
    --------
    Single site:
    
    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site import compute as cmp, edit as ed
    >>> edf = EDIFile("S01.edi")  # doctest: +SKIP
    >>> edf = ed.fill_missing(edf, how="zero",
    ...                       components=("Z",), inplace=False)
    >>> out = cmp.phase_slope(edf, band=(1.0, 1000.0))
    ... # doctest: +SKIP
    >>> set(out.keys()) == {"slope_xy", "slope_yx"}
    ... # doctest: +SKIP
    True
    
    Many sites:
    
    >>> e1 = EDIFile("S01.edi")  # doctest: +SKIP
    >>> e2 = EDIFile("S02.edi")  # doctest: +SKIP
    >>> df = cmp.phase_slope([e1, e2], band=(0.1, 10.0))
    ... # doctest: +SKIP
    >>> list(df.columns)  # doctest: +SKIP
    ['station', 'slope_xy', 'slope_yx']
    
    See Also
    --------
    pycsamt.site.compute.strike_estimate
        Strike angle by Swift-style criterion.
    pycsamt.site.compute.res_at_freq
        Apparent resistivity at a target frequency.
    
    References
    ----------
    .. [1] Simpson, F., and K. Bahr, 2005. Practical
           Magnetotellurics. Cambridge University Press.
    .. [2] Vozoff, K., 1991. The magnetotelluric method. In
           Electromagnetic methods in applied geophysics.
    """

    fmin, fmax = float(band[0]), float(band[1])
    rows = []
    for st, ed in _as_sites_iter(obj):
        f = get_freq(ed)
        Z = _get_z(ed)
        if f is None or Z is None:
            rows.append((st, np.nan, np.nan))
            continue
        m = (f >= min(fmin, fmax)) & (f <= max(fmin, fmax))
        if not np.any(m):
            rows.append((st, np.nan, np.nan))
            continue
        xf = np.log10(f[m])
        ph_xy = _phase_deg(Z[m, 0, 1])
        ph_yx = _phase_deg(Z[m, 1, 0])
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
    r"""
    Summarize or tabulate tipper magnitudes.
    
    Computes the magnitude of the tipper vector per frequency as
    
    .. math::
    
       \lVert \mathbf{T} \rVert =
       \sqrt{\lvert T_x \rvert^2 + \lvert T_y \rvert^2} ,
    
    where :math:`T_x, T_y` are the complex tipper components.  The
    result can be returned as per-frequency values or summarized
    statistics.
    
    For a single site, returns a dict.  For an iterable of sites,
    returns a ``pandas.DataFrame``.
    
    Parameters
    ----------
    obj : Any
        A single EDI-like object (e.g. ``EDIFile``) or an iterable
        of such objects.  The tipper may be attached as ``ed.Tip``,
        ``ed.T``, or ``ed.TIP`` and must expose a 2-column array
        shaped ``(n_freq, 2)``.
    per_freq : bool, optional
        If ``False`` (default), return summary statistics
        (mean, median, max).  If ``True``, return per-frequency
        values.
    
    Returns
    -------
    dict or pandas.DataFrame
        Single site:
          * ``per_freq=False`` -> ``{"mean", "median", "max"}``
          * ``per_freq=True``  -> ``{"freq", "mag"}``
        Multi-site:
          * ``per_freq=False`` -> DataFrame with columns
            ``["station", "mean", "median", "max"]``
          * ``per_freq=True``  -> DataFrame with columns
            ``["station", "freq", "mag"]``
    
    Notes
    -----
    If the site has no tipper section, summary statistics are
    ``NaN`` and per-frequency mode yields an empty result for that
    site.  To initialize missing arrays, consider
    :func:`pycsamt.site.edit.fill_missing` with ``components=("Tip",)``.
    
    Frequencies are reported from the site frequency vector.  The
    function assumes the tipper array and frequency vector are
    aligned along their first dimension.
    
    Examples
    --------
    Single site, summary stats:
    
    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site import compute as cmp, edit as ed
    >>> edf = EDIFile("S01.edi")  # doctest: +SKIP
    >>> edf = ed.fill_missing(edf, how="zero",
    ...                       components=("Tip",), inplace=False)
    >>> s = cmp.tipper_magnitude(edf, per_freq=False)
    ... # doctest: +SKIP
    >>> set(s.keys()) == {"mean", "median", "max"}
    ... # doctest: +SKIP
    True
    
    Single site, per-frequency:
    
    >>> out = cmp.tipper_magnitude(edf, per_freq=True)
    ... # doctest: +SKIP
    >>> list(out.keys())  # doctest: +SKIP
    ['freq', 'mag']
    
    Many sites, summary:
    
    >>> e1 = EDIFile("S01.edi")  # doctest: +SKIP
    >>> e2 = EDIFile("S02.edi")  # doctest: +SKIP
    >>> df = cmp.tipper_magnitude([e1, e2], per_freq=False)
    ... # doctest: +SKIP
    >>> list(df.columns)  # doctest: +SKIP
    ['station', 'mean', 'median', 'max']
    
    See Also
    --------
    pycsamt.site.edit.fill_missing
        Initialize or sanitize Z/Tip arrays in a site.
    pycsamt.site.compute.res_at_freq
        Apparent resistivity at a target frequency.
    
    References
    ----------
    .. [1] Simpson, F., and K. Bahr, 2005. Practical
           Magnetotellurics. Cambridge University Press.
    """

    rows = []
    long_rows = []

    for st, ed in _as_sites_iter(obj):
        f = get_freq(ed)
        T = _tip_arr(ed)
        if f is None or T is None:
            if per_freq:
                continue
            rows.append((st, np.nan, np.nan, np.nan))
            continue

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

    if not isinstance(obj, Iterable):
        if per_freq:
            if not long_rows:
                return {
                    "freq": np.array([]),
                    "mag": np.array([]),
                }
            fa = np.array([r[1] for r in long_rows], float)
            ma = np.array([r[2] for r in long_rows], float)
            return {"freq": fa, "mag": ma}
        if not rows:
            return {
                "mean": np.nan,
                "median": np.nan,
                "max": np.nan,
            }
        _, m, md, mx = rows[0]
        return {"mean": m, "median": md, "max": mx}

    if per_freq:
        return pd.DataFrame(
            long_rows, columns=["station", "freq", "mag"]
        )
    return pd.DataFrame(
        rows, columns=["station", "mean", "median", "max"]
    )



# Local helpers

def _get_z(ed: Any) -> Optional[np.ndarray]:
    """
    Return Z as ndarray of shape (n, 2, 2) or None.
    Looks under ed.Z.{z, impedance, _z}. Tolerates missing.
    """
    Z = getattr(ed, "Z", None)
    if Z is None:
        return None

    for nm in ("z", "impedance", "_z"):
        try:
            a = getattr(Z, nm)
        except Exception:
            a = None
        if a is None:
            continue
        arr = np.asarray(a)
        if arr.ndim == 3 and arr.shape[-2:] == (2, 2):
            return arr
    return None


def _as_sites_iter(obj: Any) -> Iterable[tuple[str, Any]]:
    for ed in iter_edifiles(obj):
        yield station_name(ed), ed


def _rho_from_z(z: complex, f: float) -> float:
    if not np.isfinite(f) or f <= 0:
        return np.nan
    return (abs(z) ** 2) / (MU_0 * 2.0 * math.pi * f)


def _rho_xy_yx(Z: np.ndarray, f: np.ndarray) -> tuple:
    zxy = Z[..., 0, 1]
    zyx = Z[..., 1, 0]
    rho_xy = np.array(
        [_rho_from_z(z, ff) for z, ff in zip(zxy, f)],
        float,
    )
    rho_yx = np.array(
        [_rho_from_z(z, ff) for z, ff in zip(zyx, f)],
        float,
    )
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
    return np.einsum("ab,nbc,cd->nad", R, Z, RT)


def _swift_cost(Z: np.ndarray, theta: float) -> float:
    Zr = _rotate_tensor(Z, theta)
    d = np.abs(Zr[:, 0, 0]) ** 2 + np.abs(Zr[:, 1, 1]) ** 2
    if d.size == 0:
        return float("inf")
    return float(np.nanmedian(d))


def _swift_theta(Z: np.ndarray) -> float:
    thetas = np.arange(0.0, 180.0, 1.0, float)
    costs = np.array([_swift_cost(Z, th) for th in thetas],
                     float)
    if not np.isfinite(costs).any():
        return float("nan")
    k = int(np.nanargmin(costs))
    return float(thetas[k])


def _phase_diff_theta(Z: np.ndarray) -> float:
    zxy = np.nanmedian(np.abs(Z[:, 0, 1]))
    zyx = np.nanmedian(np.abs(Z[:, 1, 0]))
    if not np.isfinite(zxy + zyx):
        return float("nan")
    return 0.0 if zxy >= zyx else 90.0


def _tip_arr(ed: Any) -> Optional[np.ndarray]:
    """
    Return tipper array with shape (n, 2) or None.
    Checks ed.T, ed.TIP, ed.Tip, and ed.Z.tipper.
    """
    tip_obj = None
    for cand in ("T", "TIP", "Tip"):
        tip_obj = getattr(ed, cand, None)
        if tip_obj is not None:
            break

    if tip_obj is not None:
        for nm in ("tipper", "_tipper"):
            try:
                a = getattr(tip_obj, nm)
            except Exception:
                a = None
            if a is not None:
                arr = np.asarray(a)
                if arr.ndim == 2 and arr.shape[1] == 2:
                    return arr
                if arr.ndim == 1:
                    return arr[..., None]
        return None

    Z = getattr(ed, "Z", None)
    if Z is not None:
        for nm in ("tipper", "tip", "_tipper"):
            try:
                a = getattr(Z, nm)
            except Exception:
                a = None
            if a is None:
                continue
            arr = np.asarray(a)
            if arr.ndim == 2 and arr.shape[1] == 2:
                return arr
            if arr.ndim == 1:
                return arr[..., None]
    return None
