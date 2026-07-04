
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ._core import (
    ensure_sites,
    _axes_list,
    _iter_items,
    _apply_each,
    _get_z_block,
    _name,
    hide_polar_radius_labels,
)
from .tensor import build_phase_tensor_table
from ..api._rose_style import _UNSET
from ..api.style import PYCSAMT_STYLE
from ..api.station import PYCSAMT_STATION_RENDERING
from ..api.labels import LOG10_PERIOD_LABEL, PERIOD_LABEL
from ..api.view import maybe_wrap_frame


def _rho_det_from_z(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    # ρa_det ≈ sqrt(ρ_xy ρ_yx) ; ρ = 0.2|Z|^2/f
    zx = z[:, 0, 1]
    zy = z[:, 1, 0]
    rx = 0.2 * (np.abs(zx) ** 2) / (fr + 1e-24)
    ry = 0.2 * (np.abs(zy) ** 2) / (fr + 1e-24)
    return np.sqrt(rx * ry)


def _site_coords(ed: Any) -> Optional[Tuple[float, float]]:
    """Return (lat, lon) from a Site's ``.coords``, if available.

    Real ``Site`` objects expose coordinates only via ``.coords``
    (returning ``(lat, lon, elev)``), not flat ``.lat``/``.lon``
    attributes — checking the latter alone silently falls back to
    "no coordinates" for every real station.
    """
    coords = getattr(ed, "coords", None)
    if coords is None or len(coords) < 2:
        return None
    try:
        return float(coords[0]), float(coords[1])
    except (TypeError, ValueError):
        return None


def _site_order_key(ed: Any, key: str) -> Tuple[int, float, str]:
    # sorting key for along-line order
    st = _name(ed, 0)
    if key == "lon":
        x = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
        if x is None:
            coords = _site_coords(ed)
            x = coords[1] if coords is not None else None
        return (0, float(x)) if x is not None else (1, np.inf, st)
    if key == "lat":
        y = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
        if y is None:
            coords = _site_coords(ed)
            y = coords[0] if coords is not None else None
        return (0, float(y)) if y is not None else (1, np.inf, st)
    return (0, 0.0, st)  # by name later

def _order_sites(S, sort_by: str) -> List[Any]:
    items = list(_iter_items(S))
    if sort_by in ("lon", "lat"):
        items = sorted(items, key=lambda e: _site_order_key(e, sort_by))
    else:
        items = sorted(items, key=lambda e: _name(e, 0))
    return items


def _neighbors(i: int, n: int, k: int) -> List[int]:
    lo = max(0, i - k)
    hi = min(n - 1, i + k)
    ids = list(range(lo, hi + 1))
    if i in ids:
        ids.remove(i)
    return ids


def _w_of_dist(d: np.ndarray, scheme: str, k: int) -> np.ndarray:
    d = np.abs(d).astype(float)
    if scheme == "tri":
        w = np.maximum(0.0, (k + 1) - d)
    elif scheme == "gauss":
        sig = max(1.0, 0.5 * k)
        w = np.exp(-0.5 * (d / sig) ** 2)
    else:
        w = np.ones_like(d)
    return w / (np.sum(w) + 1e-12)


def _nearest_idx(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    # nearest index in x for each y (in log-freq space)
    lx, ly = np.log10(x), np.log10(y)
    idx = np.searchsorted(lx, ly)
    idx = np.clip(idx, 1, lx.size - 1)
    left = np.abs(ly - lx[idx - 1])
    right = np.abs(ly - lx[idx])
    pick_left = left <= right
    idx[pick_left] -= 1
    return idx

def estimate_ss_ama(
    sites: Any,
    *,
    sort_by: str = "lon",    # lon|lat|name
    half_window: int = 3,     # k neighbors each side
    weights: str = "tri",     # tri|gauss|uniform
    pband: Optional[Tuple[float, float]] = None,  # (s,s)
    max_skew: Optional[float] = 6.0,  # ignore |β|>th
    robust_freq: str = "median",      # median|mean
    robust_overall: str = "median",   # median|mean
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    r"""Estimate AMA static-shift correction factors.

    Computes the Adaptive Moving-Average (AMA)
    spatial log10-resistivity trend across
    *half_window* neighbours, then returns the
    per-station deviation from that trend as
    a correction-factor table.

    Parameters
    ----------
    sites : Sites, str, Path, list, EDICollection
        EDI data source accepted by
        :func:`ensure_sites`.
    sort_by : str, default ``'lon'``
        Along-line order axis.
        ``'lon'``, ``'lat'``, or ``'name'``.
    half_window : int, default 3
        Neighbours on each side of the target.
    weights : str, default ``'tri'``
        Spatial weight scheme:
        ``'tri'`` (triangular), ``'gauss'``,
        or ``'uniform'``.
    pband : tuple of float or None
        Period band ``(p_min_s, p_max_s)``
        in seconds.  ``None`` uses all periods.
    max_skew : float or None, default 6.0
        Phase-tensor skew threshold.  Points
        where ``|beta| > max_skew`` are excluded.
    robust_freq : str, default ``'median'``
        Neighbour aggregation per frequency.
    robust_overall : str, default ``'median'``
        Reduce per-frequency deltas to a scalar.
    recursive : bool, default True
        Recursive EDI directory search.
    on_dup : str, default ``'replace'``
        Duplicate-station resolution.
    strict : bool, default False
        Raise on EDI parse errors.
    verbose : int, default 0
        Verbosity level.
    api : bool or None
        Return an APIFrame when True.

    Returns
    -------
    pandas.DataFrame
        One row per station with columns:

        ``station``
            Station identifier.
        ``delta_log10_rho``
            Estimated log10 shift.
            Positive = rho above spatial trend.
        ``fac_rho``
            Resistivity correction factor
            :math:`10^{-\delta}`.
        ``fac_z``
            Impedance correction factor
            :math:`10^{-0.5\delta}`.
        ``n_used``
            Frequencies used in the estimate.

    See Also
    --------
    correct_ss_ama : estimate + apply in one call.
    apply_ss_factors : apply a pre-built table.

    Examples
    --------
    ::

        from pycsamt.api import read_edis
        from pycsamt.emtools.ss import (
            estimate_ss_ama,
        )
        survey = read_edis("L22PLT/")
        sites  = survey.collection
        tbl = estimate_ss_ama(
            sites,
            half_window=3,
            sort_by="lon",
        )
        print(
            tbl[[
                "station",
                "delta_log10_rho",
                "fac_z",
            ]]
        )
    """
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    items = _order_sites(S, sort_by=sort_by)
    n = len(items)
    if n == 0:
        df = pd.DataFrame(
            columns=["station", "delta_log10_rho", "fac_rho",
                     "fac_z", "n_used"]
        )

        return maybe_wrap_frame(
            df, api=api, name="estimate_ss_ama",
            kind="emtools.ss.ama", source=sites,
        )

    # optional skew mask via phase-tensor table
    pt = build_phase_tensor_table(
        S, recursive=False, on_dup=on_dup,
        strict=False, verbose=verbose,
    )

    # precompute per site arrays
    ST, FR, LR = [], [], []  # name, freq, log10 ρ_det
    for i, ed in enumerate(items):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        rho = _rho_det_from_z(z, fr)
        per = 1.0 / fr
        m = np.isfinite(rho)
        if pband is not None:
            lo, hi = pband
            m &= (per >= lo) & (per <= hi)
        if max_skew is not None and not pt.empty:
            sdf = pt[pt["station"] == st]
            if not sdf.empty:
                # align skew by nearest period
                p_ref = sdf["period"].to_numpy()
                sk = np.abs(sdf["skew"].to_numpy())
                idx = _nearest_idx(1.0 / fr, p_ref)
                m &= (sk[idx] <= float(max_skew))
        fr1 = fr[m]
        lr1 = np.log10(np.maximum(rho[m], 1e-24))
        if fr1.size == 0:
            continue
        ST.append(st)
        FR.append(fr1)
        LR.append(lr1)

    if not ST:
        df = pd.DataFrame(
            columns=["station", "delta_log10_rho", "fac_rho",
                     "fac_z", "n_used"]
        )

        return maybe_wrap_frame(
            df, api=api, name="estimate_ss_ama",
            kind="emtools.ss.ama", source=sites,
        )

    # compute AMA trend and deltas
    rows = []
    for i, st in enumerate(ST):
        fr = FR[i]
        lr = LR[i]
        nbr_ids = _neighbors(i, len(ST), half_window)
        if not nbr_ids:
            continue
        # spatial weights by index distance
        dist = np.array([abs(j - i) for j in nbr_ids], dtype=float)
        w = _w_of_dist(dist, weights, half_window)
        # trend at each freq: combine neighbors nearest values
        t = np.full(fr.size, np.nan, dtype=float)
        for kf, f in enumerate(fr):
            vals = []
            for jj, j in enumerate(nbr_ids):
                frj = FR[j]
                lrj = LR[j]
                ij = _nearest_idx(frj, np.array([f]))[0]
                vals.append(lrj[ij])
            vals = np.array(vals, dtype=float)
            if robust_freq == "mean":
                t[kf] = np.nansum(w * vals)
            else:
                # weighted median (approx via repeat)
                rr = np.repeat(vals, np.maximum(1, (w * 100).astype(int)))
                t[kf] = np.nanmedian(rr)
        d = lr - t  # ≈ log10(s_i) per freq
        # A station with no overlapping finite data vs its neighbours
        # yields an all-NaN d → NaN delta → NaN factors that would
        # destroy the impedance when applied. Skip it (it stays
        # uncorrected, factor 1.0) rather than emitting a NaN row.
        if not np.any(np.isfinite(d)):
            continue
        if robust_overall == "mean":
            delta = float(np.nanmean(d))
        else:
            delta = float(np.nanmedian(d))
        if not np.isfinite(delta):
            continue
        fac_rho = 10.0 ** (-delta)
        fac_z = 10.0 ** (-0.5 * delta)
        rows.append(
            dict(
                station=st,
                delta_log10_rho=delta,
                fac_rho=fac_rho,
                fac_z=fac_z,
                n_used=int(np.isfinite(d).sum()),
            )
        )

    tbl = pd.DataFrame.from_records(rows)
    # ``rows`` is empty when no station had usable neighbours (e.g. a
    # single-station survey: _neighbors returns []). from_records([])
    # yields a column-less frame, so guard before sorting to avoid a
    # KeyError on "station" — return an empty, correctly-typed table so
    # correction degrades to a no-op instead of crashing.
    if tbl.empty or "station" not in tbl.columns:
        df = pd.DataFrame(
            columns=["station", "delta_log10_rho", "fac_rho",
                     "fac_z", "n_used"]
        )
    else:
        df = tbl.sort_values("station").reset_index(drop=True)

    return maybe_wrap_frame(
        df,
        api=api,
        name="estimate_ss_ama",
        kind="emtools.ss.ama",
        source=sites,
        description="AMA static-shift correction factors by station.",
    )


def _scale_site_Z(ed: Any, s: float) -> None:
    Z, z, fr = _get_z_block(ed)
    if Z is None:
        return
    # Guard against non-finite / non-positive factors. A NaN factor
    # (e.g. a station whose AMA delta was all-NaN) or a 0/negative one
    # would otherwise turn the impedance into NaN/zeros, silently
    # destroying the data — corrupted EDIs then fail to load with
    # "No stations with valid impedance data found". Leave Z unchanged.
    try:
        s = float(s)
    except (TypeError, ValueError):
        return
    if not np.isfinite(s) or s <= 0:
        return
    try:
        Z.z = z * s
    except Exception:
        pass
    # scale errors if present
    try:
        ze = getattr(Z, "z_err", None)
        if isinstance(ze, np.ndarray) and ze.shape == z.shape:
            setattr(Z, "z_err", ze * s)
    except Exception:
        pass


def apply_ss_factors(
    sites: Any,
    factors: Dict[str, float] | pd.DataFrame,
    *,
    key: str = "fac_z",  # fac_z: multiply Z by this
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    r"""Apply pre-computed static-shift correction factors to sites.

    Scales each site's impedance tensor Z by a per-station
    correction factor from a table (e.g. from
    :func:`estimate_ss_ama`, :func:`estimate_ss_loess`,
    etc.) or dictionary.

    Parameters
    ----------
    sites : any
        EDI data source accepted by
        :func:`ensure_sites`.
    factors : dict or pandas.DataFrame
        If DataFrame, must contain ``'station'`` and
        *key* columns. If dict, maps station names to
        correction factors.
    key : str, default ``'fac_z'``
        Column name or dict key holding the impedance
        scaling factors. Common choices are
        ``'fac_z'`` (impedance) or
        ``'fac_rho'`` (resistivity).
    inplace : bool, default False
        Modify the input Sites object.  When False,
        a corrected copy is returned.
    recursive : bool, default True
        Recursive EDI directory search.
    on_dup : str, default ``'replace'``
        Duplicate-station resolution.
    strict : bool, default False
        Raise on EDI parse errors.
    verbose : int, default 0
        Verbosity level.

    Returns
    -------
    Sites
        Corrected Sites object (same type as input).
        When *inplace* is True the original is modified
        and returned.

    See Also
    --------
    estimate_ss_ama : Estimate factors via AMA.
    estimate_ss_loess : Estimate factors via LOESS.
    """
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    # Unwrap APIFrame so the isinstance check below works
    _inner = getattr(factors, "df", None)
    if isinstance(_inner, pd.DataFrame):
        factors = _inner

    if isinstance(factors, pd.DataFrame):
        if "station" in factors.columns and key in factors.columns:
            fmap = {r.station: float(r[key]) for _, r in
                    factors.iterrows()}
        else:
            raise ValueError("bad factors table")
    else:
        fmap = {str(k): float(v) for k, v in factors.items()}

    def _one(Si):
        ed = next(_iter_items(Si))
        st = _name(ed, 0)
        s = fmap.get(st, 1.0)
        _scale_site_Z(ed, s)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def correct_ss_ama(
    sites: Any,
    *,
    sort_by: str = "lon",
    half_window: int = 3,
    weights: str = "tri",
    pband: Optional[Tuple[float, float]] = None,
    max_skew: Optional[float] = 6.0,
    robust_freq: str = "median",
    robust_overall: str = "median",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    r"""Correct static shift by the AMA method.

    Estimates per-station log10-resistivity
    shift factors with :func:`estimate_ss_ama`,
    then scales each site's impedance tensor Z
    by the corresponding ``fac_z`` column.

    Parameters
    ----------
    sites : Sites, str, Path, list, EDICollection
        EDI data source.
    sort_by : str, default ``'lon'``
        Along-line order axis for AMA estimation.
    half_window : int, default 3
        Neighbours on each side of the target.
    weights : str, default ``'tri'``
        Spatial weight scheme (``'tri'``,
        ``'gauss'``, or ``'uniform'``).
    pband : tuple of float or None
        Period band ``(p_min_s, p_max_s)``
        in seconds.
    max_skew : float or None, default 6.0
        Phase-tensor skew exclusion threshold.
    robust_freq : str, default ``'median'``
        Neighbour aggregation per frequency.
    robust_overall : str, default ``'median'``
        Reduce per-frequency deltas to a scalar.
    inplace : bool, default False
        Modify the input Sites object in place.
        When False, returns a corrected copy.
    recursive : bool, default True
        Recursive EDI directory search.
    on_dup : str, default ``'replace'``
        Duplicate-station resolution.
    strict : bool, default False
        Raise on EDI parse errors.
    verbose : int, default 0
        Verbosity level.

    Returns
    -------
    Sites
        Corrected Sites object (same type as
        input).  When *inplace* is True the
        original object is modified and returned.

    See Also
    --------
    estimate_ss_ama : inspect factors before apply.
    apply_ss_factors : apply a custom factor table.

    Examples
    --------
    ::

        from pycsamt.api import read_edis
        from pycsamt.emtools.ss import (
            correct_ss_ama,
        )
        survey     = read_edis("L22PLT/")
        sites      = survey.collection
        sites_corr = correct_ss_ama(
            sites,
            half_window=3,
            sort_by="lon",
        )
    """
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    tbl = estimate_ss_ama(
        S,
        sort_by=sort_by,
        half_window=half_window,
        weights=weights,
        pband=pband,
        max_skew=max_skew,
        robust_freq=robust_freq,
        robust_overall=robust_overall,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    return apply_ss_factors(
        S,
        tbl,
        key="fac_z",
        inplace=inplace,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )


def _prep_lr_curves(
    sites: Any,
    *,
    pband: Optional[Tuple[float, float]],
    max_skew: Optional[float],
    recursive: bool,
    on_dup: str,
    strict: bool,
    verbose: int,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    pt = build_phase_tensor_table(
        S, recursive=False, on_dup=on_dup,
        strict=False, verbose=verbose,
    )
    ST, FR, LR = [], [], []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        rho = _rho_det_from_z(z, fr)
        per = 1.0 / fr
        m = np.isfinite(rho)
        if pband is not None:
            lo, hi = pband
            m &= (per >= lo) & (per <= hi)
        if max_skew is not None and not pt.empty:
            sdf = pt[pt["station"] == st]
            if not sdf.empty:
                p_ref = sdf["period"].to_numpy()
                sk = np.abs(sdf["skew"].to_numpy())
                idx = _nearest_idx(1.0 / fr, p_ref)
                m &= (sk[idx] <= float(max_skew))
        fr1 = fr[m]
        lr1 = np.log10(np.maximum(rho[m], 1e-24))
        if fr1.size == 0:
            continue
        ST.append(st); FR.append(fr1); LR.append(lr1)
    return ST, FR, LR


def _tricube(u: np.ndarray) -> np.ndarray:
    v = np.clip(1.0 - np.abs(u) ** 3, 0.0, 1.0)
    return v ** 3


def _loess_at_center(
    x: np.ndarray, y: np.ndarray, w: np.ndarray, poly: int
) -> float:
    # eval at x=0 (center); poly 0 or 1
    if y.size == 0:
        return np.nan
    if poly <= 0 or np.allclose(x, 0.0):
        num = np.sum(w * y)
        den = np.sum(w) + 1e-12
        return float(num / den)
    X = np.vstack([np.ones_like(x), x]).T
    W = np.diag(w)
    A = X.T @ W @ X
    b = X.T @ W @ y
    try:
        beta = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        beta = np.linalg.pinv(A) @ b
    return float(beta[0])


def _loess_trend_for_site(
    i: int,
    FR: List[np.ndarray],
    LR: List[np.ndarray],
    *,
    k: int,
    poly: int,
    it: int,
) -> Tuple[np.ndarray, np.ndarray]:
    fr = FR[i]; lr = LR[i]
    n = len(FR)
    ids = list(range(max(0, i - k), min(n - 1, i + k) + 1))
    if i in ids:
        ids.remove(i)
    # station coord as index distance
    xi = np.array([j - i for j in ids], dtype=float)
    t = np.full(fr.size, np.nan, dtype=float)
    for m, f in enumerate(fr):
        ys = []
        for j in ids:
            ij = _nearest_idx(FR[j], np.array([f]))[0]
            ys.append(LR[j][ij])
        ys = np.asarray(ys, dtype=float)
        if ys.size == 0:
            continue
        # tricube weights on normalized |x|
        u = np.abs(xi) / (k + 1e-12)
        w = _tricube(u)
        # robust bisquare iterations
        for _ in range(max(1, it)):
            mu = _loess_at_center(xi, ys, w, poly)
            r = ys - mu
            s = np.nanmedian(np.abs(r)) + 1e-12
            u = np.clip(r / (6.0 * s), -1.0, 1.0)
            wb = (1.0 - u * u) ** 2
            w = w * wb
        t[m] = _loess_at_center(xi, ys, w, poly)
    return fr, t


def estimate_ss_loess(
    sites: Any,
    *,
    half_window: int = 3,
    poly: int = 1,
    it: int = 2,
    pband: Optional[Tuple[float, float]] = None,
    max_skew: Optional[float] = 6.0,
    summary: str = "median",  # median|mean
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    r"""Estimate static-shift factors via locally-weighted regression (LOESS).

    Fits a local polynomial trend across neighbouring
    stations in the along-line direction, then returns
    the per-station deviation from that trend as
    correction factors.

    Parameters
    ----------
    sites : Sites, str, Path, list, EDICollection
        EDI data source accepted by
        :func:`ensure_sites`.
    half_window : int, default 3
        Neighbours on each side of the target.
    poly : int, default 1
        Polynomial degree (0=constant, 1=linear).
    it : int, default 2
        Robust iteration count.
    pband : tuple of float or None
        Period band :math:`(p_{min}, p_{max})` in seconds.
        ``None`` uses all periods.
    max_skew : float or None, default 6.0
        Phase-tensor skew threshold.  Points where
        :math:`|\\beta| > ` *max_skew* are excluded.
    summary : str, default ``'median'``
        Per-station aggregation: ``'median'`` or ``'mean'``.
    recursive : bool, default True
        Recursive EDI directory search.
    on_dup : str, default ``'replace'``
        Duplicate-station resolution.
    strict : bool, default False
        Raise on EDI parse errors.
    verbose : int, default 0
        Verbosity level.
    api : bool or None
        Return an APIFrame when True.

    Returns
    -------
    pandas.DataFrame
        One row per station with columns:
        ``station``, ``delta_log10_rho``,
        ``fac_rho``, ``fac_z``, ``n_used``.

    See Also
    --------
    estimate_ss_ama : AMA (moving average) method.
    estimate_ss_bilateral : Bilateral filtering method.
    """
    ST, FR, LR = _prep_lr_curves(
        sites, pband=pband, max_skew=max_skew,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if not ST:
        df = pd.DataFrame(
            columns=["station", "delta_log10_rho",
                     "fac_rho", "fac_z", "n_used"]
        )

        return maybe_wrap_frame(
            df, api=api, name="estimate_ss_loess",
            kind="emtools.ss.loess", source=sites,
        )
    rows = []
    for i, st in enumerate(ST):
        fr, tr = _loess_trend_for_site(
            i, FR, LR, k=int(half_window),
            poly=int(poly), it=int(it),
        )
        lr = LR[i]
        d = lr - tr
        if summary == "mean":
            delta = float(np.nanmean(d))
        else:
            delta = float(np.nanmedian(d))
        rows.append(
            dict(
                station=st,
                delta_log10_rho=delta,
                fac_rho=10.0 ** (-delta),
                fac_z=10.0 ** (-0.5 * delta),
                n_used=int(np.isfinite(d).sum()),
            )
        )
    df = pd.DataFrame.from_records(rows).sort_values(
        "station"
    ).reset_index(drop=True)

    return maybe_wrap_frame(
        df,
        api=api,
        name="estimate_ss_loess",
        kind="emtools.ss.loess",
        source=sites,
        description="LOESS static-shift correction factors by station.",
    )


def _bilateral_trend_for_site(
    i: int,
    FR: List[np.ndarray],
    LR: List[np.ndarray],
    *,
    k: int,
    sig_dist: Optional[float],
    sig_val: Optional[float],
) -> Tuple[np.ndarray, np.ndarray]:
    fr = FR[i]; lr = LR[i]
    n = len(FR)
    ids = list(range(max(0, i - k), min(n - 1, i + k) + 1))
    if i in ids:
        ids.remove(i)
    xi = np.array([j - i for j in ids], dtype=float)
    sd = float(sig_dist) if sig_dist else max(1.0, 0.5 * k)
    t = np.full(fr.size, np.nan, dtype=float)
    for m, f in enumerate(fr):
        ys = []
        for j in ids:
            ij = _nearest_idx(FR[j], np.array([f]))[0]
            ys.append(LR[j][ij])
        ys = np.asarray(ys, dtype=float)
        if ys.size == 0:
            continue
        # spatial kernel
        ws = np.exp(-0.5 * (xi / sd) ** 2)
        # range kernel (value similarity)
        sv = (np.nanmedian(np.abs(ys - np.nanmedian(ys)))
              + 1e-12) if sig_val is None else float(sig_val)
        wr = np.exp(-0.5 * ((ys - lr[m]) / sv) ** 2)
        w = ws * wr
        t[m] = float(np.sum(w * ys) / (np.sum(w) + 1e-12))
    return fr, t


def estimate_ss_bilateral(
    sites: Any,
    *,
    half_window: int = 4,
    sig_dist: Optional[float] = None,
    sig_val: Optional[float] = None,
    pband: Optional[Tuple[float, float]] = None,
    max_skew: Optional[float] = 6.0,
    summary: str = "median",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    r"""Estimate static-shift factors via bilateral filtering.

    Applies a combined spatial and range-based Gaussian
    filter (bilateral filter) to compute a local trend,
    then returns per-station deviations as correction
    factors.

    Parameters
    ----------
    sites : Sites, str, Path, list, EDICollection
        EDI data source accepted by
        :func:`ensure_sites`.
    half_window : int, default 4
        Spatial window (neighbours each side).
    sig_dist : float or None
        Spatial Gaussian width (in index units).  When
        ``None``, defaults to :math:`0.5 \\times
        \\texttt{half\\_window}`.
    sig_val : float or None
        Range (value) Gaussian width.  When ``None``,
        estimated from data.
    pband : tuple of float or None
        Period band :math:`(p_{min}, p_{max})` in seconds.
    max_skew : float or None, default 6.0
        Phase-tensor skew threshold.
    summary : str, default ``'median'``
        Aggregation: ``'median'`` or ``'mean'``.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`ensure_sites`.
    api : bool or None
        Return an APIFrame when True.

    Returns
    -------
    pandas.DataFrame
        One row per station with columns:
        ``station``, ``delta_log10_rho``,
        ``fac_rho``, ``fac_z``, ``n_used``.

    See Also
    --------
    estimate_ss_ama : Moving-average method.
    estimate_ss_loess : Local polynomial method.
    """
    ST, FR, LR = _prep_lr_curves(
        sites, pband=pband, max_skew=max_skew,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if not ST:
        df = pd.DataFrame(
            columns=["station", "delta_log10_rho",
                     "fac_rho", "fac_z", "n_used"]
        )

        return maybe_wrap_frame(
            df, api=api, name="estimate_ss_bilateral",
            kind="emtools.ss.bilateral", source=sites,
        )
    rows = []
    for i, st in enumerate(ST):
        fr, tr = _bilateral_trend_for_site(
            i, FR, LR, k=int(half_window),
            sig_dist=sig_dist, sig_val=sig_val,
        )
        lr = LR[i]; d = lr - tr
        delta = float(np.nanmedian(d)) if summary == "median" \
            else float(np.nanmean(d))
        rows.append(
            dict(
                station=st,
                delta_log10_rho=delta,
                fac_rho=10.0 ** (-delta),
                fac_z=10.0 ** (-0.5 * delta),
                n_used=int(np.isfinite(d).sum()),
            )
        )
    df = pd.DataFrame.from_records(rows).sort_values(
        "station"
    ).reset_index(drop=True)

    return maybe_wrap_frame(
        df,
        api=api,
        name="estimate_ss_bilateral",
        kind="emtools.ss.bilateral",
        source=sites,
        description="Bilateral static-shift correction factors by station.",
    )


def estimate_ss_refmedian(
    sites: Any,
    *,
    pband: Optional[Tuple[float, float]] = None,
    max_skew: Optional[float] = 6.0,
    smooth_sites: int = 0,  # optional along-site median
    summary: str = "median",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    r"""Estimate static-shift factors via reference-median method.

    Computes a global frequency-wise median resistivity
    across all stations, then estimates per-station shifts
    as deviations from this reference curve.

    Parameters
    ----------
    sites : Sites, str, Path, list, EDICollection
        EDI data source.
    pband : tuple of float or None
        Period band :math:`(p_{min}, p_{max})` in seconds.
    max_skew : float or None, default 6.0
        Phase-tensor skew threshold.
    smooth_sites : int, default 0
        Optional smoothing window (reserved for future use).
    summary : str, default ``'median'``
        Aggregation: ``'median'`` or ``'mean'``.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`ensure_sites`.
    api : bool or None
        Return an APIFrame when True.

    Returns
    -------
    pandas.DataFrame
        One row per station with columns:
        ``station``, ``delta_log10_rho``,
        ``fac_rho``, ``fac_z``, ``n_used``.

    See Also
    --------
    estimate_ss_ama : Moving-average method.
    estimate_ss_loess : Local polynomial method.
    """
    ST, FR, LR = _prep_lr_curves(
        sites, pband=pband, max_skew=max_skew,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if not ST:
        df = pd.DataFrame(
            columns=["station", "delta_log10_rho",
                     "fac_rho", "fac_z", "n_used"]
        )

        return maybe_wrap_frame(
            df, api=api, name="estimate_ss_refmedian",
            kind="emtools.ss.refmedian", source=sites,
        )
    # build a global ref via frequency-wise median
    # grid = union of all frequencies
    G = np.unique(np.concatenate(FR))
    Ref = np.full(G.size, np.nan, dtype=float)
    for k, f in enumerate(G):
        vals = []
        for fr, lr in zip(FR, LR):
            j = _nearest_idx(fr, np.array([f]))[0]
            vals.append(lr[j])
        vals = np.asarray(vals, dtype=float)
        if smooth_sites > 0:
            # along-site median in a window around each station
            # here we fallback to global median to stay simple
            pass
        Ref[k] = np.nanmedian(vals)
    rows = []
    for st, fr, lr in zip(ST, FR, LR):
        idx = _nearest_idx(G, fr)
        d = lr - Ref[idx]
        delta = float(np.nanmedian(d)) if summary == "median" \
            else float(np.nanmean(d))
        rows.append(
            dict(
                station=st,
                delta_log10_rho=delta,
                fac_rho=10.0 ** (-delta),
                fac_z=10.0 ** (-0.5 * delta),
                n_used=int(np.isfinite(d).sum()),
            )
        )
    df = pd.DataFrame.from_records(rows).sort_values(
        "station"
    ).reset_index(drop=True)

    return maybe_wrap_frame(
        df,
        api=api,
        name="estimate_ss_refmedian",
        kind="emtools.ss.refmedian",
        source=sites,
        description="Reference-median static-shift factors by station.",
    )

# ----------------------- SS visualization (QC) --------------------------- #


def _pair_sites(
    before: Any, after: Any, *, verbose: int = 0
) -> Dict[str, Tuple[Any, Any]]:
    B = ensure_sites(before, recursive=False, strict=False)
    A = ensure_sites(after, recursive=False, strict=False)
    bm = {}
    for i, ed in enumerate(_iter_items(B)):
        bm[_name(ed, i)] = ed
    am = {}
    for i, ed in enumerate(_iter_items(A)):
        am[_name(ed, i)] = ed
    common = {}
    for st, edb in bm.items():
        eda = am.get(st, None)
        if eda is not None:
            common[st] = (edb, eda)
    return common


def plot_ss_delta_psection(
    before: Any,
    after: Any,
    *,
    axis_y: str = "logperiod",
    vlim: Optional[float] = None,
    pband: Optional[Tuple[float, float]] = None,
    figsize: Tuple[float, float] = (9.0, 4.8),
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    r"""Plot pseudosection of static-shift change (corrected minus original).

    Displays a heatmap showing the pointwise difference
    :math:`\Delta\log_{10}\rho = \rho_{after} - 
    \rho_{before}` across all stations and frequencies
    on a log-period y-axis.

    Parameters
    ----------
    before : any
        EDI data source (uncorrected sites).
    after : any
        EDI data source (corrected sites).
    axis_y : str, default ``'logperiod'``
        Y-axis scale: ``'logperiod'`` or ``'period'``.
    vlim : float or None
        Symmetric colour range :math:`\pm \texttt{vlim}`.
        When ``None``, auto-scales from data.
    pband : tuple of float or None
        Period band :math:`(p_{min}, p_{max})` in seconds.
    figsize : (float, float), default (9, 4.8)
        Figure size.
    verbose : int, default 0
        Verbosity level.
    ax : matplotlib.axes.Axes or None
        Draw on existing axes.

    Returns
    -------
    matplotlib.axes.Axes
    """
    pairs = _pair_sites(before, after, verbose=verbose)
    rows = []
    yvals = []
    labs = []
    for k, (edb, eda) in enumerate(pairs.values()):
        Zb, zb, frb = _get_z_block(edb)
        Za, za, fra = _get_z_block(eda)
        if Zb is None or Za is None:
            continue
        rb = _rho_det_from_z(zb, frb)
        ra = _rho_det_from_z(za, fra)
        perb = 1.0 / frb
        if pband is not None:
            lo, hi = pband
            m = (perb >= lo) & (perb <= hi)
        else:
            m = np.isfinite(rb)
        if not np.any(m):
            continue
        j = _nearest_idx(fra, frb[m])
        dlog = np.log10(np.maximum(ra[j], 1e-24)) - np.log10(
            np.maximum(rb[m], 1e-24)
        )
        rows.append(dlog)
        yvals.append(np.log10(perb[m]) if axis_y == "logperiod"
                     else perb[m])
        labs.append(list(pairs.keys())[k])
    if not rows:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no overlap", ha="center", va="center")
        return ax
    # union y grid and assemble image
    yg = np.unique(np.concatenate([y for y in yvals]))
    M = np.full((len(labs), yg.size), np.nan, dtype=float)
    for i, (yy, dl) in enumerate(zip(yvals, rows)):
        ii = np.searchsorted(yg, yy)
        ii = np.clip(ii, 0, yg.size - 1)
        M[i, ii] = dl
    Z = M.T  # (y, station)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    v = Z[np.isfinite(Z)]
    if vlim is None and v.size:
        a = np.nanpercentile(v, 95)
        vlim = float(max(a, 0.1))
    im = ax.imshow(
        Z,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-(vlim or 0.5),
        vmax=(vlim or 0.5),
    )
    ax.set_ylabel(LOG10_PERIOD_LABEL if axis_y == "logperiod"
                  else PERIOD_LABEL)
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(len(labs), dtype=float),
        labs,
        preset="pseudosection",
        xlim=(-0.5, len(labs) - 0.5),
    )
    yt = np.linspace(0, Z.shape[0] - 1, num=min(8, Z.shape[0]))
    yv = np.linspace(yg.min(), yg.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("Δ log10 ρ_det (after − before)")
    return ax


def plot_ss_station_curves(
    before: Any,
    after: Any,
    *,
    station: Optional[str] = None,
    pband: Optional[Tuple[float, float]] = None,
    figsize: Tuple[float, float] = (7.8, 4.2),
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    r"""Plot before-and-after apparent-resistivity curves for a single station.

    Overlays two 1-D sounding curves (before correction
    and after correction) on a period x-axis to visualize
    the magnitude and frequency-dependence of the
    correction at one location.

    Parameters
    ----------
    before : any
        Uncorrected EDI data.
    after : any
        Corrected EDI data.
    station : str or None
        Station identifier. When ``None``, the first
        common station is used.
    pband : tuple of float or None
        Period band :math:`(p_{min}, p_{max})` in seconds.
    figsize : (float, float), default (7.8, 4.2)
        Figure size.
    verbose : int, default 0
        Verbosity level.
    ax : matplotlib.axes.Axes or None
        Draw on existing axes.

    Returns
    -------
    matplotlib.axes.Axes
    """
    pairs = _pair_sites(before, after, verbose=verbose)
    if not pairs:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no common stations",
                ha="center", va="center")
        return ax
    if station is None:
        station = list(pairs.keys())[0]
    edb, eda = pairs.get(station, (None, None))
    if edb is None or eda is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "station not found",
                ha="center", va="center")
        return ax
    Zb, zb, frb = _get_z_block(edb)
    Za, za, fra = _get_z_block(eda)
    if Zb is None or Za is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no Z blocks", ha="center", va="center")
        return ax
    rb = _rho_det_from_z(zb, frb)
    ra = _rho_det_from_z(za, fra)
    pb = 1.0 / frb
    pa = 1.0 / fra
    if pband is not None:
        lo, hi = pband
        mb = (pb >= lo) & (pb <= hi)
        ma = (pa >= lo) & (pa <= hi)
    else:
        mb = np.isfinite(rb)
        ma = np.isfinite(ra)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    _cs = PYCSAMT_STYLE.correction
    ax.set_xscale("log")
    ax.plot(pb[mb], rb[mb], **_cs.before.plot_kwargs(ms=3.5))
    ax.plot(pa[ma], ra[ma], **_cs.after.plot_kwargs(ms=3.5))
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("ρ_det (Ω·m)")
    ax.set_title(str(station))
    ax.grid(True, alpha=0.25, which="both")
    ax.legend()
    return ax


def plot_ss_delta_profile(
    before: Any,
    after: Any,
    *,
    pband: Optional[Tuple[float, float]] = None,
    robust: str = "median",  # median|mean
    figsize: Tuple[float, float] = (8.6, 3.6),
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    r"""Plot per-station static-shift correction amplitudes as a bar chart.

    Shows the median (or mean) of the frequency-dependent
    correction :math:`\Delta\log_{10}\rho` at each station,
    making it easy to identify spatial patterns in the
    applied corrections.

    Parameters
    ----------
    before : any
        Uncorrected EDI data.
    after : any
        Corrected EDI data.
    pband : tuple of float or None
        Period band :math:`(p_{min}, p_{max})` in seconds.
    robust : str, default ``'median'``
        Aggregation method: ``'median'`` or ``'mean'``.
    figsize : (float, float), default (8.6, 3.6)
        Figure size.
    verbose : int, default 0
        Verbosity level.
    ax : matplotlib.axes.Axes or None
        Draw on existing axes.

    Returns
    -------
    matplotlib.axes.Axes
    """
    pairs = _pair_sites(before, after, verbose=verbose)
    labs = []
    deltas = []
    for st, (edb, eda) in pairs.items():
        Zb, zb, frb = _get_z_block(edb)
        Za, za, fra = _get_z_block(eda)
        if Zb is None or Za is None:
            continue
        rb = _rho_det_from_z(zb, frb)
        ra = _rho_det_from_z(za, fra)
        pb = 1.0 / frb
        if pband is not None:
            lo, hi = pband
            m = (pb >= lo) & (pb <= hi)
        else:
            m = np.isfinite(rb)
        if not np.any(m):
            continue
        j = _nearest_idx(fra, frb[m])
        d = np.log10(np.maximum(ra[j], 1e-24)) - np.log10(
            np.maximum(rb[m], 1e-24)
        )
        val = float(np.nanmedian(d)) if robust == "median" \
            else float(np.nanmean(d))
        labs.append(st)
        deltas.append(val)
    if not labs:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no overlap", ha="center", va="center")
        return ax
    order = np.argsort(labs)
    labs = [labs[i] for i in order]
    deltas = [deltas[i] for i in order]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(labs))
    ax.axhline(0.0, color="0.7", lw=1.0)
    ax.bar(x, deltas, width=0.8)
    ax.set_ylabel("Δ log10 ρ_det (after − before)")
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        x.astype(float),
        labs,
        preset="pseudosection",
        xlim=(-0.5, len(labs) - 0.5),
    )
    return ax


# ═══════════════════════════════════════════════════════════════════════════════
# Publication-quality static-shift comparison plots
# ═══════════════════════════════════════════════════════════════════════════════

# ── internal rendering helpers ─────────────────────────────────────────────── #

def _ss_sort_freqs(
    freqs: np.ndarray, *arrays: np.ndarray
) -> Tuple[np.ndarray, ...]:
    """Return (sorted_freqs, sorted_arr1, …) all ascending in Hz."""
    order = np.argsort(freqs)
    return (freqs[order],) + tuple(a[:, order] for a in arrays)


def _logT_edges(freqs: np.ndarray) -> np.ndarray:
    """Cell boundary positions in log10-period space (sorted ascending Hz)."""
    lT = np.log10(1.0 / np.asarray(freqs, float))
    n  = lT.size
    if n > 1:
        d      = np.diff(lT)
        e      = np.empty(n + 1)
        e[0]   = lT[0]   - 0.5 * abs(d[0])
        e[1:-1] = lT[:-1] + 0.5 * d
        e[-1]  = lT[-1]  + 0.5 * abs(d[-1])
    else:
        e = np.array([lT[0] - 0.5, lT[0] + 0.5])
    return e


def _x_edges(x_centres: np.ndarray) -> np.ndarray:
    n = x_centres.size
    if n == 1:
        return np.array([x_centres[0] - 0.5, x_centres[0] + 0.5])
    e = np.empty(n + 1)
    e[0]    = x_centres[0]  - 0.5 * (x_centres[1]    - x_centres[0])
    e[1:-1] = 0.5 * (x_centres[:-1] + x_centres[1:])
    e[-1]   = x_centres[-1] + 0.5 * (x_centres[-1]   - x_centres[-2])
    return e


def _pcolor_lT(
    ax: plt.Axes,
    data: np.ndarray,       # (n_f, n_st)
    freqs: np.ndarray,
    x_centres: np.ndarray,
    vmin: float,
    vmax: float,
    *,
    cmap: str = "RdYlBu_r",
    period_up: bool = True,
) -> Any:
    """Render one panel with pcolormesh on a log10-period y-axis."""
    ye = _logT_edges(freqs)
    xe = _x_edges(x_centres)
    X, Y = np.meshgrid(xe, ye)
    qm = ax.pcolormesh(
        X, Y, data,
        cmap=cmap, vmin=vmin, vmax=vmax,
        shading="flat", rasterized=True,
    )
    if period_up:
        ax.invert_yaxis()
    return qm


def _set_lT_yticks(
    ax: plt.Axes,
    freqs: np.ndarray,
    *,
    n: int = 7,
    fontsize: int = 7,
    ylabel: str = "Period (s)",
) -> None:
    lT   = np.log10(1.0 / freqs)
    pos  = np.linspace(lT.min(), lT.max(), n)
    labs = []
    for v in pos:
        r = round(v)
        labs.append(
            f"$10^{{{r}}}$" if abs(r - v) < 0.04 else f"$10^{{{v:.1f}}}$"
        )
    ax.set_yticks(pos)
    ax.set_yticklabels(labs, fontsize=fontsize)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=8)


def _set_station_xticks(
    ax: plt.Axes,
    n_st: int,
    labels: List[str],
    *,
    rotation: float = 45.0,
    fontsize: int = 7,
    xlabel: str = "",
) -> None:
    x  = np.arange(n_st, dtype=float)
    ha = "right" if rotation > 20 else "center"
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=rotation, ha=ha, fontsize=fontsize)
    ax.set_xlim(-0.5, n_st - 0.5)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=8)


def _joint_clim(
    *arrays: np.ndarray,
    pct: Tuple[float, float] = (2.0, 98.0),
    hard_min: float = 0.5,
    hard_max: float = 4.5,
) -> Tuple[float, float]:
    flat = np.concatenate([a.ravel() for a in arrays])
    fin  = flat[np.isfinite(flat)]
    if not fin.size:
        return hard_min, hard_max
    return (
        max(float(np.percentile(fin, pct[0])), hard_min),
        min(float(np.percentile(fin, pct[1])), hard_max),
    )


# ─────────────────────────────────────────────────────────────────────────────

def plot_ss_comparison_psection(
    logRho_before: np.ndarray,
    logRho_after: np.ndarray,
    *,
    freqs: np.ndarray,
    station_labels: Optional[List[str]] = None,
    show_delta: bool = True,
    cmap: str = "RdYlBu_r",
    delta_cmap: str = "RdBu_r",
    clim: Optional[Tuple[float, float]] = None,
    clim_pct: Tuple[float, float] = (2.0, 98.0),
    delta_vlim: Optional[float] = None,
    delta_vlim_pct: float = 95.0,
    period_up: bool = True,
    title_before: str = "(a) Before static-shift correction",
    title_after: str = "(b) After static-shift correction",
    title_delta: str = r"(c) Correction amplitude $\Delta\log_{10}\rho$",
    suptitle: str = "",
    xlabel: str = "Station",
    ylabel: str = "Period (s)",
    n_yticks: int = 7,
    colorbar_label: str = r"$\log_{10}\,\rho_a$ (Ω·m)",
    delta_colorbar_label: str = r"$\Delta\log_{10}\rho$",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    figsize: Optional[Tuple[float, float]] = None,
    axes: Optional[Any] = None,
) -> plt.Figure:
    """
    Two- or three-panel pseudo-section comparison for static-shift correction.

    The *before* and *after* panels share a colour scale so that the
    station-dependent vertical offsets are directly visible.  The optional
    third panel shows the pointwise difference
    Δ log₁₀ ρ = after − before on a diverging scale, making the spatial
    pattern of the correction explicit.

    Parameters
    ----------
    logRho_before : ndarray, shape ``(n_st, n_f)``
        Log₁₀ apparent resistivity *before* static-shift correction (Ω·m).
    logRho_after : ndarray, shape ``(n_st, n_f)``
        Log₁₀ apparent resistivity *after* static-shift correction (Ω·m).
    freqs : ndarray, shape ``(n_f,)``
        Frequency array in Hz.  Need not be sorted.
    station_labels : list of str or None
        X-axis tick labels.  Defaults to ``"0", "1", …``.
    show_delta : bool, default ``True``
        Append a third panel showing Δ log₁₀ ρ.
    cmap : str, default ``"RdYlBu_r"``
        Colormap for the before/after panels.
    delta_cmap : str, default ``"RdBu_r"``
        Diverging colormap for the Δ panel.
    clim : (vmin, vmax) or None
        Explicit colour limits (log₁₀ Ω·m) shared by the before/after panels.
    clim_pct : (lo, hi), default ``(2.0, 98.0)``
        Percentile bounds for automatic *clim*.
    delta_vlim : float or None
        Symmetric limit ``(−δ, +δ)`` for the Δ panel.  When ``None``,
        derived from *delta_vlim_pct* of ``|Δ|``.
    delta_vlim_pct : float, default ``95.0``
    period_up : bool, default ``True``
        Long period at the top of each panel (MT convention).
    title_before, title_after, title_delta : str
        Per-panel titles.  Pass ``""`` to suppress.
    suptitle : str
        Figure-level title.
    xlabel, ylabel : str
        Axis labels.
    n_yticks : int, default ``7``
        Number of log-period y-ticks.
    colorbar_label, delta_colorbar_label : str
    tick_label_rotation : float, default ``45.0``
        Station tick rotation (degrees).
    tick_fontsize : int, default ``7``
    figsize : (w, h) or None
        Override automatic size.
    axes : sequence of Axes or None
        Pre-created axes (length 2 without delta, 3 with).

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    logRho_before = np.asarray(logRho_before, dtype=float)
    logRho_after  = np.asarray(logRho_after,  dtype=float)
    freqs         = np.asarray(freqs, dtype=float).ravel()

    n_st = logRho_before.shape[0]
    freqs, logRho_before, logRho_after = _ss_sort_freqs(
        freqs, logRho_before, logRho_after
    )

    if station_labels is None:
        station_labels = [str(i) for i in range(n_st)]
    x_centres = np.arange(n_st, dtype=float)

    n_panels = 3 if show_delta else 2
    if axes is None:
        if figsize is None:
            figsize = (11.0, 3.8 * n_panels)
        fig, axes = plt.subplots(
            n_panels, 1,
            figsize=figsize,
            sharex=True,
            gridspec_kw={"hspace": 0.42},
        )
    else:
        axes = list(axes)
        fig  = axes[0].get_figure()

    # ── shared colour limits ───────────────────────────────────────────────
    if clim is None:
        vmin, vmax = _joint_clim(logRho_before, logRho_after, pct=clim_pct)
    else:
        vmin, vmax = float(clim[0]), float(clim[1])

    # ── before panel ──────────────────────────────────────────────────────
    qm_b = _pcolor_lT(
        axes[0], logRho_before.T, freqs, x_centres, vmin, vmax,
        cmap=cmap, period_up=period_up,
    )
    if title_before:
        axes[0].set_title(title_before, fontsize=9, fontweight="bold", pad=3)
    _set_lT_yticks(axes[0], freqs, n=n_yticks,
                   fontsize=tick_fontsize, ylabel=ylabel)

    # ── after panel ───────────────────────────────────────────────────────
    _pcolor_lT(
        axes[1], logRho_after.T, freqs, x_centres, vmin, vmax,
        cmap=cmap, period_up=period_up,
    )
    if title_after:
        axes[1].set_title(title_after, fontsize=9, fontweight="bold", pad=3)
    _set_lT_yticks(axes[1], freqs, n=n_yticks,
                   fontsize=tick_fontsize, ylabel=ylabel)

    # shared colorbar spanning the two main panels
    cb_main = fig.colorbar(
        qm_b,
        ax=[axes[0], axes[1]],
        fraction=0.018, pad=0.01, aspect=35,
    )
    cb_main.set_label(colorbar_label, fontsize=8)
    cb_main.ax.tick_params(labelsize=7)

    # ── delta panel ───────────────────────────────────────────────────────
    if show_delta:
        delta = logRho_after - logRho_before
        fin_d = np.abs(delta)[np.isfinite(delta)]
        if delta_vlim is None:
            delta_vlim = (
                float(np.percentile(fin_d, delta_vlim_pct))
                if fin_d.size else 0.5
            )
        qm_d = _pcolor_lT(
            axes[2], delta.T, freqs, x_centres,
            -delta_vlim, delta_vlim,
            cmap=delta_cmap, period_up=period_up,
        )
        if title_delta:
            axes[2].set_title(title_delta, fontsize=9, fontweight="bold", pad=3)
        _set_lT_yticks(axes[2], freqs, n=n_yticks,
                       fontsize=tick_fontsize, ylabel=ylabel)
        cb_d = fig.colorbar(
            qm_d, ax=axes[2], fraction=0.018, pad=0.01, aspect=30,
        )
        cb_d.set_label(delta_colorbar_label, fontsize=8)
        cb_d.ax.tick_params(labelsize=7)

    # ── station axis (top panel; shared x across panels) ─────────────────
    PYCSAMT_STATION_RENDERING.apply(
        axes[0],
        x_centres,
        station_labels,
        preset="pseudosection",
        xlim=(x_centres[0] - 0.5, x_centres[-1] + 0.5),
    )
    for ax in axes[1:]:
        ax.tick_params(
            axis="x",
            which="both",
            top=False,
            bottom=False,
            labeltop=False,
            labelbottom=False,
        )

    if suptitle:
        fig.suptitle(suptitle, fontsize=10, fontweight="bold", y=1.005)

    return fig


# ─────────────────────────────────────────────────────────────────────────────

def plot_ss_1d_curves(
    logRho_before: np.ndarray,
    logRho_after: np.ndarray,
    *,
    freqs: np.ndarray,
    stations: Optional[Any] = None,
    station_labels: Optional[List[str]] = None,
    n_cols: int = 4,
    max_stations: int = 16,
    color_before   = _UNSET,   # default: PYCSAMT_STYLE.correction.before.color
    color_after    = _UNSET,   # default: PYCSAMT_STYLE.correction.after.color
    ls_before      = _UNSET,   # default: PYCSAMT_STYLE.correction.before.ls
    ls_after       = _UNSET,   # default: PYCSAMT_STYLE.correction.after.ls
    marker_before  = _UNSET,   # default: PYCSAMT_STYLE.correction.before.marker
    marker_after   = _UNSET,   # default: PYCSAMT_STYLE.correction.after.marker
    marker_size    = _UNSET,   # default: PYCSAMT_STYLE.correction.before.ms
    lw             = _UNSET,   # default: PYCSAMT_STYLE.correction.before.lw
    log_period: bool = True,
    show_shift_annotation: bool = True,
    annotation_fontsize: int = 7,
    ylabel: str = r"$\log_{10}\,\rho_a$ (Ω·m)",
    xlabel: str = "Period (s)",
    axes=None,
    figsize: Optional[Tuple[float, float]] = None,
    title: str = "",
    legend_loc: str = "best",
    show_grid: bool = True,
) -> plt.Figure:
    """
    Per-station 1-D apparent-resistivity curves: before and after correction.

    Lays out a grid of subplots (one per selected station) each showing
    the before/after sounding curves on a period x-axis.  A small annotation
    reports the mean correction amplitude Δ per station, making it easy to
    spot outliers.

    Parameters
    ----------
    logRho_before : ndarray, shape ``(n_st, n_f)``
    logRho_after  : ndarray, shape ``(n_st, n_f)``
    freqs : ndarray, shape ``(n_f,)``  Hz.
    stations : list of int, list of str, or None
        Stations to display.  Integers are row indices into *logRho_before*.
        Strings are matched against *station_labels*.  ``None`` → all
        stations, capped at *max_stations*.
    station_labels : list of str or None
        Label for each row.  Defaults to ``"0", "1", …``.
    n_cols : int, default ``4``
        Subplot grid columns.
    max_stations : int, default ``16``
        Cap when *stations* is ``None``.
    color_before : str, default ``"#2c7bb6"`` (blue)
    color_after  : str, default ``"#d7191c"`` (red)
    ls_before : str, default ``"--"``
    ls_after  : str, default ``"-"``
    marker_before, marker_after : str
    marker_size : float, default ``3.0``
    lw : float, default ``1.2``
    log_period : bool, default ``True``
        Log-scale period x-axis.
    show_shift_annotation : bool, default ``True``
        Print mean Δ log₁₀ ρ in the lower-right corner of each subplot.
    annotation_fontsize : int, default ``7``
    ylabel, xlabel : str
    figsize : (w, h) or None
    title : str
        Figure-level title.
    legend_loc : str, default ``"best"``
        Legend location (first subplot only).
    show_grid : bool, default ``True``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    # ── resolve visual style from PYCSAMT_STYLE.correction ───────────────
    _cs = PYCSAMT_STYLE.correction
    if color_before  is _UNSET: color_before  = _cs.before.color
    if color_after   is _UNSET: color_after   = _cs.after.color
    if ls_before     is _UNSET: ls_before     = _cs.before.ls
    if ls_after      is _UNSET: ls_after      = _cs.after.ls
    if marker_before is _UNSET: marker_before = _cs.before.marker
    if marker_after  is _UNSET: marker_after  = _cs.after.marker
    if marker_size   is _UNSET: marker_size   = _cs.before.ms
    if lw            is _UNSET: lw            = _cs.before.lw

    logRho_before = np.asarray(logRho_before, dtype=float)
    logRho_after  = np.asarray(logRho_after,  dtype=float)
    freqs         = np.asarray(freqs, dtype=float).ravel()
    n_st_total    = logRho_before.shape[0]

    if station_labels is None:
        station_labels = [str(i) for i in range(n_st_total)]

    # ── resolve station selection ─────────────────────────────────────────
    if stations is None:
        idx = list(range(min(n_st_total, max_stations)))
    else:
        idx = []
        for s in stations:
            if isinstance(s, (int, np.integer)):
                if 0 <= int(s) < n_st_total:
                    idx.append(int(s))
            else:
                try:
                    idx.append(station_labels.index(str(s)))
                except ValueError:
                    pass
        if not idx:
            idx = list(range(min(n_st_total, max_stations)))

    n_shown = len(idx)
    n_rows  = max(1, int(np.ceil(n_shown / n_cols)))
    if figsize is None:
        figsize = (n_cols * 3.2, n_rows * 2.8)

    axes_given = _axes_list(axes, n_shown) if axes is not None else None
    if axes_given is None:
        fig, axes_grid = plt.subplots(
            n_rows, n_cols,
            figsize=figsize,
            squeeze=False,
        )
        axes_flat = axes_grid.ravel()
    else:
        axes_flat = np.asarray(axes_given, dtype=object)
        fig = axes_flat[0].figure

    # sort periods ascending for clean curves
    order  = np.argsort(1.0 / freqs)
    per_s  = (1.0 / freqs)[order]

    for k, si in enumerate(idx):
        ax = axes_flat[k]
        rb = logRho_before[si][order]
        ra = logRho_after[si][order]
        fin = np.isfinite(rb) & np.isfinite(ra)

        ax.plot(
            per_s[fin], rb[fin],
            color=color_before, ls=ls_before, lw=lw,
            marker=marker_before, ms=marker_size,
            label="before",
        )
        ax.plot(
            per_s[fin], ra[fin],
            color=color_after, ls=ls_after, lw=lw,
            marker=marker_after, ms=marker_size,
            label="after",
        )

        if log_period:
            ax.set_xscale("log")

        ax.set_title(station_labels[si], fontsize=8, fontweight="bold", pad=2)

        if show_shift_annotation and np.any(fin):
            delta_mean = float(np.nanmean(ra[fin] - rb[fin]))
            sign = "+" if delta_mean >= 0 else ""
            ax.text(
                0.97, 0.04,
                f"Δ={sign}{delta_mean:.2f}",
                transform=ax.transAxes,
                ha="right", va="bottom",
                fontsize=annotation_fontsize,
                color="#555555",
                bbox=dict(fc="white", ec="none", alpha=0.7, pad=1.5),
            )

        if show_grid:
            ax.grid(True, which="both", alpha=0.2, lw=0.5)

        ax.tick_params(labelsize=7)
        ax.set_xlabel(xlabel, fontsize=7)
        ax.set_ylabel(ylabel, fontsize=7)

        if k == 0:
            ax.legend(loc=legend_loc, fontsize=7, framealpha=0.8)

    for k in range(n_shown, len(axes_flat)):
        axes_flat[k].set_visible(False)

    if title:
        fig.suptitle(title, fontsize=10, fontweight="bold", y=1.01)

    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────

def plot_ss_summary(
    logRho_before: np.ndarray,
    logRho_after: np.ndarray,
    *,
    freqs: np.ndarray,
    station_labels: Optional[List[str]] = None,
    cmap: str = "RdYlBu_r",
    delta_cmap: str = "RdBu_r",
    clim: Optional[Tuple[float, float]] = None,
    clim_pct: Tuple[float, float] = (2.0, 98.0),
    delta_vlim: Optional[float] = None,
    delta_vlim_pct: float = 95.0,
    period_up: bool = True,
    n_yticks: int = 7,
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    colorbar_label: str = r"$\log_{10}\,\rho_a$ (Ω·m)",
    shift_bar_color: str = "#4c72b0",
    shift_bar_neg_color: str = "#c44e52",
    shift_robust: str = "median",
    suptitle: str = "",
    axes=None,
    figsize: Optional[Tuple[float, float]] = None,
) -> plt.Figure:
    """
    Four-panel summary figure for static-shift correction.

    Layout::

        ┌──────────────┬──────────────┐
        │  (a) Before  │  (b) After   │  shared y-axis · shared colorbar
        ├──────────────┴──────────────┤
        │  (c) Δ log₁₀ ρ section     │  diverging colorbar
        ├──────────────────────────── ┤
        │  (d) Per-station shift bar  │  positive/negative coloured bars
        └─────────────────────────────┘

    Parameters
    ----------
    logRho_before : ndarray, shape ``(n_st, n_f)``
    logRho_after  : ndarray, shape ``(n_st, n_f)``
    freqs : ndarray, shape ``(n_f,)``  Hz.
    station_labels : list of str or None
        X-axis tick labels for all panels.
    cmap : str, default ``"RdYlBu_r"``
    delta_cmap : str, default ``"RdBu_r"``
    clim, clim_pct : see :func:`plot_ss_comparison_psection`.
    delta_vlim, delta_vlim_pct : see :func:`plot_ss_comparison_psection`.
    period_up : bool, default ``True``
    n_yticks : int, default ``7``
    tick_label_rotation : float, default ``45.0``
    tick_fontsize : int, default ``7``
    colorbar_label : str
    shift_bar_color : str
        Bar colour for positive per-station shifts (default blue).
    shift_bar_neg_color : str
        Bar colour for negative shifts (default red).
    shift_robust : ``"median"`` | ``"mean"``
        Aggregation used to reduce per-frequency shifts to a scalar per
        station for panel (d).
    suptitle : str
        Figure-level title.
    figsize : (w, h) or None

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    logRho_before = np.asarray(logRho_before, dtype=float)
    logRho_after  = np.asarray(logRho_after,  dtype=float)
    freqs         = np.asarray(freqs, dtype=float).ravel()

    n_st = logRho_before.shape[0]
    freqs, logRho_before, logRho_after = _ss_sort_freqs(
        freqs, logRho_before, logRho_after
    )

    if station_labels is None:
        station_labels = [str(i) for i in range(n_st)]
    x_centres = np.arange(n_st, dtype=float)

    axes_given = _axes_list(axes, 4)
    if axes_given is None:
        if figsize is None:
            figsize = (13.0, 14.0)

        fig = plt.figure(figsize=figsize)
        gs  = fig.add_gridspec(
            3, 2,
            height_ratios=[1, 1, 0.55],
            hspace=0.46,
            wspace=0.20,
        )
        ax_before = fig.add_subplot(gs[0, 0])
        ax_after  = fig.add_subplot(gs[0, 1], sharey=ax_before)
        ax_delta  = fig.add_subplot(gs[1, :])
        ax_bar    = fig.add_subplot(gs[2, :])
    else:
        ax_before, ax_after, ax_delta, ax_bar = axes_given
        fig = ax_before.figure

    # ── colour limits ──────────────────────────────────────────────────────
    if clim is None:
        vmin, vmax = _joint_clim(logRho_before, logRho_after, pct=clim_pct)
    else:
        vmin, vmax = float(clim[0]), float(clim[1])

    # ── (a) before ────────────────────────────────────────────────────────
    qm_b = _pcolor_lT(
        ax_before, logRho_before.T, freqs, x_centres, vmin, vmax,
        cmap=cmap, period_up=period_up,
    )
    ax_before.set_title("(a) Before correction",
                        fontsize=9, fontweight="bold", pad=3)
    _set_lT_yticks(ax_before, freqs, n=n_yticks,
                   fontsize=tick_fontsize, ylabel="Period (s)")
    _set_station_xticks(
        ax_before, n_st, station_labels,
        rotation=tick_label_rotation, fontsize=tick_fontsize,
        xlabel="Station",
    )

    # ── (b) after ─────────────────────────────────────────────────────────
    _pcolor_lT(
        ax_after, logRho_after.T, freqs, x_centres, vmin, vmax,
        cmap=cmap, period_up=period_up,
    )
    ax_after.set_title("(b) After correction",
                       fontsize=9, fontweight="bold", pad=3)
    _set_lT_yticks(ax_after, freqs, n=n_yticks, fontsize=tick_fontsize, ylabel="")
    ax_after.tick_params(axis="y", labelleft=False)
    _set_station_xticks(
        ax_after, n_st, station_labels,
        rotation=tick_label_rotation, fontsize=tick_fontsize,
        xlabel="Station",
    )

    # shared colorbar for (a)+(b)
    cb_main = fig.colorbar(
        qm_b, ax=[ax_before, ax_after], fraction=0.015, pad=0.01, aspect=35,
    )
    cb_main.set_label(colorbar_label, fontsize=8)
    cb_main.ax.tick_params(labelsize=7)

    # ── (c) delta section (full-width) ────────────────────────────────────
    delta = logRho_after - logRho_before
    fin_d = np.abs(delta)[np.isfinite(delta)]
    if delta_vlim is None:
        delta_vlim = float(np.percentile(fin_d, delta_vlim_pct)) if fin_d.size else 0.5
    qm_d = _pcolor_lT(
        ax_delta, delta.T, freqs, x_centres,
        -delta_vlim, delta_vlim,
        cmap=delta_cmap, period_up=period_up,
    )
    ax_delta.set_title(
        r"(c) Correction amplitude  $\Delta\log_{10}\rho$ (after − before)",
        fontsize=9, fontweight="bold", pad=3,
    )
    _set_lT_yticks(ax_delta, freqs, n=n_yticks,
                   fontsize=tick_fontsize, ylabel="Period (s)")
    _set_station_xticks(
        ax_delta, n_st, station_labels,
        rotation=tick_label_rotation, fontsize=tick_fontsize,
        xlabel="Station",
    )
    cb_d = fig.colorbar(qm_d, ax=ax_delta, fraction=0.015, pad=0.01, aspect=30)
    cb_d.set_label(r"$\Delta\log_{10}\rho$", fontsize=8)
    cb_d.ax.tick_params(labelsize=7)

    # ── (d) per-station shift bar chart ───────────────────────────────────
    delta_per_st = np.where(np.isfinite(delta), delta, np.nan)
    if shift_robust == "median":
        shift_vals = np.nanmedian(delta_per_st, axis=1)
    else:
        shift_vals = np.nanmean(delta_per_st, axis=1)

    bar_colors = [
        shift_bar_color if v >= 0 else shift_bar_neg_color
        for v in shift_vals
    ]
    ax_bar.bar(x_centres, shift_vals, color=bar_colors, width=0.75, alpha=0.85)
    ax_bar.axhline(0.0, color="0.4", lw=0.8, ls="--")
    ax_bar.set_title(
        r"(d) Per-station shift  $\langle\Delta\log_{10}\rho\rangle$",
        fontsize=9, fontweight="bold", pad=3,
    )
    ax_bar.set_ylabel(r"$\Delta\log_{10}\rho$ (dex)", fontsize=8)
    ax_bar.grid(True, axis="y", alpha=0.25, lw=0.5)
    ax_bar.tick_params(labelsize=tick_fontsize)
    _set_station_xticks(
        ax_bar, n_st, station_labels,
        rotation=tick_label_rotation, fontsize=tick_fontsize,
        xlabel="Station",
    )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold", y=1.005)

    return fig


# ------------------- one-shot QC wrappers (sites in) -------------------- #

def _select_kwargs(kws: Dict[str, Any], allowed: set) -> Dict[str, Any]:
    return {k: v for k, v in kws.items() if k in allowed}


_ALLOWED = {
    "ama": {
        "sort_by", "half_window", "weights", "pband",
        "max_skew", "robust_freq", "robust_overall",
        "recursive", "on_dup", "strict", "verbose",
    },
    "loess": {
        "half_window", "poly", "it", "pband", "max_skew",
        "summary", "recursive", "on_dup", "strict", "verbose",
    },
    "bilateral": {
        "half_window", "sig_dist", "sig_val", "pband",
        "max_skew", "summary", "recursive", "on_dup",
        "strict", "verbose",
    },
    "refmedian": {
        "pband", "max_skew", "smooth_sites", "summary",
        "recursive", "on_dup", "strict", "verbose",
    },
}


def _correct_sites(
    sites: Any,
    method: str,
    **corr: Any,
):
    S = ensure_sites(
        sites, recursive=corr.get("recursive", True),
        on_dup=corr.get("on_dup", "replace"),
        strict=corr.get("strict", False),
        verbose=corr.get("verbose", 0),
    )
    m = method.lower()
    if m == "ama":
        kw = _select_kwargs(corr, _ALLOWED["ama"])
        return correct_ss_ama(S, inplace=False, **kw)
    if m == "loess":
        kw = _select_kwargs(corr, _ALLOWED["loess"])
        tbl = estimate_ss_loess(S, **kw)
        return apply_ss_factors(S, tbl, inplace=False)
    if m == "bilateral":
        kw = _select_kwargs(corr, _ALLOWED["bilateral"])
        tbl = estimate_ss_bilateral(S, **kw)
        return apply_ss_factors(S, tbl, inplace=False)
    if m == "refmedian":
        kw = _select_kwargs(corr, _ALLOWED["refmedian"])
        tbl = estimate_ss_refmedian(S, **kw)
        return apply_ss_factors(S, tbl, inplace=False)
    raise ValueError(f"unknown method: {method}")


def ss_qc_psection(
    sites: Any,
    *,
    method: str = "ama",
    return_sites: bool = False,
    # plot opts
    axis_y: str = "logperiod",
    vlim: Optional[float] = None,
    pband: Optional[Tuple[float, float]] = None,
    figsize: Tuple[float, float] = (9.0, 4.8),
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
    # correction kwargs (forwarded)
    **corr: Any,
):
    r"""Estimate static-shift correction and plot delta pseudosection.

    Combines automatic static-shift estimation with
    a heatmap visualization in one call.  A convenience
    wrapper around a correction estimator and
    :func:`plot_ss_delta_psection`.

    Parameters
    ----------
    sites : any
        EDI paths or :class:`~pycsamt.site.base.Sites`.
    method : str, default ``'ama'``
        Correction method: ``'ama'``, ``'loess'``,
        ``'bilateral'``, or ``'refmedian'``.
    return_sites : bool, default False
        When ``True``, return ``(ax, corrected_sites)``.
    axis_y, vlim, pband, figsize
        Forwarded to :func:`plot_ss_delta_psection`.
    verbose : int, default 0
        Verbosity level.
    ax : matplotlib.axes.Axes or None
        Draw on existing axes.
    **corr :
        Forwarded to the correction estimator.

    Returns
    -------
    matplotlib.axes.Axes or (Axes, Sites)
    """
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _correct_sites(S0, method, **corr)
    ax = plot_ss_delta_psection(
        S0,
        S1,
        axis_y=axis_y,
        vlim=vlim,
        pband=pband,
        figsize=figsize,
        verbose=verbose,
        ax=ax,
    )
    return (ax, S1) if return_sites else ax


def ss_qc_station_curves(
    sites: Any,
    *,
    method: str = "ama",
    station: Optional[str] = None,
    return_sites: bool = False,
    # plot opts
    pband: Optional[Tuple[float, float]] = None,
    figsize: Tuple[float, float] = (7.8, 4.2),
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
    # correction kwargs
    **corr: Any,
):
    r"""Estimate correction and plot before/after curves for one station.

    A convenience wrapper combining automatic
    static-shift estimation with 1-D curve visualization.

    Parameters
    ----------
    sites : any
        EDI paths or Sites object.
    method : str, default ``'ama'``
        Correction method.
    station : str or None
        Station identifier. When ``None``, uses the first.
    return_sites : bool, default False
        When ``True``, return ``(ax, corrected_sites)``.
    pband, figsize, verbose, ax
        Forwarded to :func:`plot_ss_station_curves`.
    **corr :
        Forwarded to the correction estimator.

    Returns
    -------
    matplotlib.axes.Axes or (Axes, Sites)
    """
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _correct_sites(S0, method, **corr)
    ax = plot_ss_station_curves(
        S0,
        S1,
        station=station,
        pband=pband,
        figsize=figsize,
        verbose=verbose,
        ax=ax,
    )
    return (ax, S1) if return_sites else ax


def ss_qc_profile(
    sites: Any,
    *,
    method: str = "ama",
    return_sites: bool = False,
    # plot opts
    pband: Optional[Tuple[float, float]] = None,
    robust: str = "median",
    figsize: Tuple[float, float] = (8.6, 3.6),
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
    # correction kwargs
    **corr: Any,
):
    r"""Estimate correction and plot per-station shift profile.

    A convenience wrapper for automatic static-shift
    estimation with bar-chart visualization of the
    per-station amplitudes.

    Parameters
    ----------
    sites : any
        EDI paths or Sites object.
    method : str, default ``'ama'``
        Correction method.
    return_sites : bool, default False
        When ``True``, return ``(ax, corrected_sites)``.
    pband, robust, figsize, verbose, ax
        Forwarded to :func:`plot_ss_delta_profile`.
    **corr :
        Forwarded to the correction estimator.

    Returns
    -------
    matplotlib.axes.Axes or (Axes, Sites)
    """
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _correct_sites(S0, method, **corr)
    ax = plot_ss_delta_profile(
        S0,
        S1,
        pband=pband,
        robust=robust,
        figsize=figsize,
        verbose=verbose,
        ax=ax,
    )
    return (ax, S1) if return_sites else ax


def ss_comparison_psection(
    sites: Any,
    *,
    method: str = "ama",
    return_sites: bool = False,
    station_labels: Optional[List[str]] = None,
    show_delta: bool = True,
    cmap: str = "RdYlBu_r",
    delta_cmap: str = "RdBu_r",
    clim: Optional[Tuple[float, float]] = None,
    clim_pct: Tuple[float, float] = (2.0, 98.0),
    delta_vlim: Optional[float] = None,
    delta_vlim_pct: float = 95.0,
    period_up: bool = True,
    suptitle: str = "",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    figsize: Optional[Tuple[float, float]] = None,
    verbose: int = 0,
    **corr: Any,
) -> Any:
    """
    Correct *sites* for static shift and plot a comparison pseudo-section.

    A convenience wrapper that combines :func:`~pycsamt.emtools.correct_ss_ama`
    (or the chosen *method*) with :func:`plot_ss_comparison_psection`.

    Parameters
    ----------
    sites : any
        EDI paths, glob pattern, or :class:`~pycsamt.site.base.Sites`
        accepted by :func:`~pycsamt.emtools.ensure_sites`.
    method : ``"ama"`` | ``"loess"`` | ``"bilateral"`` | ``"refmedian"``
        Static-shift estimator.
    return_sites : bool, default ``False``
        When ``True``, return ``(fig, corrected_sites)`` instead of *fig*.
    **corr :
        Forwarded to the correction estimator.

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
        Or ``(fig, corrected_sites)`` when *return_sites* is ``True``.

    See Also
    --------
    plot_ss_comparison_psection : Lower-level function that accepts pre-built
        arrays directly.
    """
    S0 = ensure_sites(sites, recursive=True, strict=False, verbose=verbose)
    S1 = _correct_sites(S0, method, **corr)

    items0 = list(_iter_items(S0))
    items1 = list(_iter_items(S1))

    labels = [_name(e, k) for k, e in enumerate(items0)]
    n_st   = len(labels)

    # collect rho_det arrays from each site pair
    all_f: set = set()
    rho0_map: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    rho1_map: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}

    for k, (e0, e1) in enumerate(zip(items0, items1)):
        _, z0, fr0 = _get_z_block(e0)
        _, z1, fr1 = _get_z_block(e1)
        if z0 is None:
            continue
        st = labels[k]
        rho0_map[st] = (_rho_det_from_z(z0, fr0), fr0)
        rho1_map[st] = (_rho_det_from_z(z1, fr1), fr1)
        all_f.update(fr0.tolist())

    if not all_f:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No Z-tensor data found in sites",
                ha="center", va="center")
        return (fig, S1) if return_sites else fig

    freqs_union = np.array(sorted(all_f))
    n_f         = freqs_union.size
    logRho_b    = np.full((n_st, n_f), np.nan)
    logRho_a    = np.full((n_st, n_f), np.nan)

    for k, st in enumerate(labels):
        if st not in rho0_map:
            continue
        rho0, fr0 = rho0_map[st]
        rho1, fr1 = rho1_map[st]
        j0 = _nearest_idx(freqs_union, fr0)
        j1 = _nearest_idx(freqs_union, fr1)
        logRho_b[k, j0] = np.log10(np.maximum(rho0, 1e-24))
        logRho_a[k, j1] = np.log10(np.maximum(rho1, 1e-24))

    fig = plot_ss_comparison_psection(
        logRho_b, logRho_a,
        freqs=freqs_union,
        station_labels=station_labels if station_labels is not None else labels,
        show_delta=show_delta,
        cmap=cmap,
        delta_cmap=delta_cmap,
        clim=clim,
        clim_pct=clim_pct,
        delta_vlim=delta_vlim,
        delta_vlim_pct=delta_vlim_pct,
        period_up=period_up,
        suptitle=suptitle,
        tick_label_rotation=tick_label_rotation,
        tick_fontsize=tick_fontsize,
        figsize=figsize,
    )
    return (fig, S1) if return_sites else fig


# ---------------------- Static-shift radar (polar) ---------------------- #

def _rot_mat(th: np.ndarray) -> np.ndarray:
    c = np.cos(th); s = np.sin(th)
    R = np.zeros((th.size, 2, 2), dtype=float)
    R[:, 0, 0] = c; R[:, 0, 1] = s
    R[:, 1, 0] = -s; R[:, 1, 1] = c
    return R


def _rotate_z(z: np.ndarray, ang_deg: np.ndarray) -> np.ndarray:
    th = np.radians(ang_deg.astype(float))
    R = _rot_mat(th)
    Rt = np.transpose(R, (0, 2, 1))
    return R @ z @ Rt


def _pt_phi_for_station(
    S, station: str, fr: np.ndarray, stat: str
) -> np.ndarray:
    tb = build_phase_tensor_table(
        S, recursive=False, on_dup="replace",
        strict=False, verbose=0,
    )
    if tb is None or getattr(tb, "empty", False):
        return np.zeros(fr.size, dtype=float)
    sdf = tb[tb["station"] == station]
    if sdf.empty:
        return np.zeros(fr.size, dtype=float)
    # try common strike/azimuth column names
    for col in ("azimuth", "strike", "phi", "theta"):
        if col in sdf.columns:
            p_ref = sdf["period"].to_numpy(dtype=float)
            phi = sdf[col].to_numpy(dtype=float)
            j = _nearest_idx(1.0 / fr, p_ref)
            if stat == "median":
                return phi[j]
            if stat == "mean":
                return phi[j]
            # fallback to per-frequency direct mapping
            return phi[j]
    return np.zeros(fr.size, dtype=float)


def plot_ss_radar(
    sites: Any,
    *,
    station: Optional[str] = None,
    pband: Optional[Tuple[float, float]] = None,
    rotate: str = "pt",        # pt|none|deg
    rotate_stat: str = "median",
    rotate_deg: float = 0.0,   # used when rotate="deg"
    radial: str = "log10rho",  # log10rho|rho
    theta_axis: str = "logperiod",  # logperiod|period|freq
    fill_between: bool = False,
    colors         = _UNSET,   # default: (mt.xy.color, mt.yx.color)
    marker         = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.marker
    ms             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.ms
    lw             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.lw
    ls             = _UNSET,   # default: PYCSAMT_STYLE.mt.xy.ls
    figsize: Tuple[float, float] = (4.8, 4.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    eps: float = 1e-24, 
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    r"""Plot apparent resistivity against period on a polar grid.

    Displays the off-diagonal impedance components (xy and
    yx) as radial curves on a polar coordinate system,
    where the azimuthal angle encodes frequency (or period)
    and the radius encodes resistivity magnitude. Useful for
    detecting anisotropy and strike angles across the full
    frequency spectrum.

    Parameters
    ----------
    sites : any
        EDI data source.
    station : str or None
        Station identifier. When ``None``, uses the first.
    pband : tuple of float or None
        Period band :math:`(p_{min}, p_{max})` in seconds.
    rotate : str, default ``'pt'``
        Rotation mode: ``'pt'`` (phase-tensor strike),
        ``'deg'`` (fixed angle), or ``'none'``
        (no rotation).
    rotate_stat : str, default ``'median'``
        Per-frequency aggregation for phase-tensor rotation.
    rotate_deg : float, default 0.0
        Fixed rotation angle (degrees) when rotate='deg'.
    radial : str, default ``'log10rho'``
        Radial scale: ``'log10rho'`` (log base 10
        of apparent resistivity) or ``'rho'``
        (linear resistivity).
    theta_axis : str, default ``'logperiod'``
        Angular axis: ``'logperiod'``, ``'period'``, or
        ``'freq'`` (Hz).
    fill_between : bool, default False
        Shade the region between xy and yx curves.
    colors : tuple or _UNSET
        (color_xy, color_yx).  Defaults from style.
    marker, ms, lw, ls : _UNSET or values
        Line and marker style. Defaults from style.
    figsize : (float, float), default (4.8, 4.8)
        Figure size.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`ensure_sites`.
    eps : float, default 1e-24
        Numerical floor to avoid division by zero.
    ax : matplotlib.axes.Axes or None
        Draw on existing axes (auto-creates polar if needed).

    Returns
    -------
    matplotlib.axes.Axes
        Polar axes object.
    """
    # ── resolve visual style from PYCSAMT_STYLE.mt ───────────────────────
    _mt = PYCSAMT_STYLE.mt
    if colors is _UNSET: colors = (_mt.xy.color, _mt.yx.color)
    if marker is _UNSET: marker = _mt.xy.marker
    if ms     is _UNSET: ms     = _mt.xy.ms
    if lw     is _UNSET: lw     = _mt.xy.lw
    if ls     is _UNSET: ls     = _mt.xy.ls

    def _ensure_polar_axis(axis: Optional[plt.Axes]) -> plt.Axes:
        if axis is None:
            _, new_ax = plt.subplots(
                figsize=figsize, subplot_kw={"polar": True}
            )
            return new_ax
        if getattr(axis, "name", "") == "polar":
            return axis
        fig = axis.figure
        pos = axis.get_position()
        axis.remove()
        return fig.add_axes(pos, projection="polar")

    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    # pick station
    sel = {}
    for i, ed in enumerate(_iter_items(S)):
        sel[_name(ed, i)] = ed
    if not sel:
        ax = _ensure_polar_axis(ax)
        ax.text(0.5, 0.5, "no sites", ha="center", va="center")
        return ax
    if station is None:
        station = sorted(sel.keys())[0]
    ed = sel.get(station, None)
    if ed is None:
        ax = _ensure_polar_axis(ax)
        ax.text(0.5, 0.5, "station not found",
                ha="center", va="center")
        return ax
    Z, z, fr = _get_z_block(ed)
    if Z is None:
        ax = _ensure_polar_axis(ax)
        ax.text(0.5, 0.5, "no Z", ha="center", va="center")
        return ax

    # rotation angles per frequency
    if rotate == "pt":
        ang = _pt_phi_for_station(S, station, fr, rotate_stat)
    elif rotate == "deg":
        ang = np.full(fr.size, float(rotate_deg), dtype=float)
    else:
        ang = np.zeros(fr.size, dtype=float)

    zr = _rotate_z(z, ang)

    # select band + compute radii
    per = 1.0 / fr
    m = np.ones(fr.size, dtype=bool)
    if pband is not None:
        lo, hi = pband; m &= (per >= lo) & (per <= hi)
    xy = zr[:, 0, 1]; yx = zr[:, 1, 0]
    if radial == "rho":
        r_xy = 0.2 * (np.abs(xy) ** 2) / (fr + eps)
        r_yx = 0.2 * (np.abs(yx) ** 2) / (fr + eps)
    else:
        r_xy = np.log10(0.2 * (np.abs(xy) ** 2) / (fr + eps))
        r_yx = np.log10(0.2 * (np.abs(yx) ** 2) / (fr + eps))

    # theta mapping
    if theta_axis == "freq":
        x = fr
    elif theta_axis == "period":
        x = per
    else:
        x = np.log10(np.maximum(per, 1e-24))
        # normalize to [0, 2π)
        x = (x - x.min()) / (x.max() - x.min() + eps)
        x = 2.0 * np.pi * x
    th = x if theta_axis == "logperiod" else \
         (2.0 * np.pi * (x - x.min()) /
          (x.max() - x.min() + 1e-24))

    th = th[m]; r1 = r_xy[m]; r2 = r_yx[m]

    ax = _ensure_polar_axis(ax)
    # set polar style: 0 at north, CW
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    hide_polar_radius_labels(ax)

    # plot
    suffix = " (rot)" if rotate != "none" else ""
    ax.plot(th, r1, ls=ls, lw=lw, marker=marker,
            ms=ms, color=colors[0], label=f"ρa_xy{suffix}")
    ax.plot(th, r2, ls=ls, lw=lw, marker=marker,
            ms=ms, color=colors[1], label=f"ρa_yx{suffix}")
    if fill_between:
        lo = np.minimum(r1, r2); hi = np.maximum(r1, r2)
        ax.fill_between(th, lo, hi, color="0.5", alpha=0.10)

    ax.grid(True, alpha=0.25)
    hide_polar_radius_labels(ax)
    ax.set_title(str(station), pad=10)
    ax.set_ylabel("")
    # Outside the polar axes entirely: bbox_to_anchor=(0.02, 0.02) used
    # to sit right on top of the 225-degree angular tick label, which
    # always lands near that same lower-left corner of the bounding box.
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.05),
              frameon=False, fontsize=8)
    return ax


# ========= Near-surface effect detection (lei2017) ======================== #

_TYPE_COLORS: Dict[str, str] = {
    "clean":        "#2ca02c",   # green
    "static":       "#1f77b4",   # blue
    "near_surface": "#ff7f0e",   # orange
    "mixed":        "#d62728",   # red
}

_NS_COLS = [
    "station", "n_hf", "n_lf",
    "sigma_hf", "sigma_lf", "ns_index",
    "slope_hf", "slope_lf", "gradient_delta",
    "ss_delta_log10", "ns_flag", "ss_flag", "distortion_type",
]


def _unwrap_ns(ed: Any) -> Any:
    """Unwrap a Sites-level Site wrapper to its underlying EDI-like object."""
    edi = getattr(ed, "edi", None)
    if edi is not None and hasattr(edi, "Z"):
        return edi
    return ed


def _log_slope(log_f: np.ndarray, log_rho: np.ndarray) -> float:
    """Median d(log10 ρ)/d(log10 f) via finite differences."""
    if log_f.size < 2:
        return float("nan")
    dlr = np.diff(log_rho)
    dlf = np.diff(log_f)
    valid = np.abs(dlf) > 1e-10
    if not valid.any():
        return float("nan")
    return float(np.nanmedian(dlr[valid] / dlf[valid]))


def _ama_residuals_ns(
    FR: List[np.ndarray],
    LR: List[np.ndarray],
    *,
    half_window: int,
    weights: str,
) -> List[np.ndarray]:
    """Per-frequency log10ρ residuals vs AMA spatial trend for every site."""
    n = len(FR)
    out = []
    for i in range(n):
        fr = FR[i]
        lr = LR[i]
        nbr_ids = _neighbors(i, n, half_window)
        t = np.full(fr.size, np.nan, dtype=float)
        if nbr_ids:
            dist = np.array([abs(j - i) for j in nbr_ids], dtype=float)
            w = _w_of_dist(dist, weights, half_window)
            for kf, f in enumerate(fr):
                vals = np.array(
                    [LR[j][_nearest_idx(FR[j], np.array([f]))[0]]
                     for j in nbr_ids],
                    dtype=float,
                )
                rr = np.repeat(vals, np.maximum(1, (w * 100).astype(int)))
                t[kf] = np.nanmedian(rr)
        out.append(lr - t)
    return out


def detect_near_surface(
    sites: Any,
    *,
    f_split: float = 1.0,
    pband: Optional[Tuple[float, float]] = None,
    ns_threshold: float = 2.0,
    ss_threshold: float = 0.1,
    sort_by: str = "lon",
    half_window: int = 3,
    weights: str = "tri",
    max_skew: Optional[float] = 6.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    """
    Detect and classify near-surface distortion in CSAMT/MT apparent
    resistivity curves.

    Distinguishes between two types of distortion:

    * **Static effect** — frequency-independent multiplicative shift of the
      whole ρ_a curve.  Addressable by AMA/LOESS static-shift correction.
    * **Near-surface effect** — frequency-dependent distortion concentrated
      at high frequencies (f ≥ *f_split*), caused by shallow inhomogeneities.
      *Not* correctable by conventional static-shift methods; 2-D inversion
      is recommended.

    Three per-station diagnostics are computed from the residuals of the
    ρ_a curve relative to an AMA spatial trend:

    1. **NS index**  η = σ_HF / σ_LF — spread ratio between the
       high-frequency (f ≥ *f_split*) and low-frequency bands.  η >> 1 is
       the hallmark of near-surface contamination.

    2. **Gradient delta**  Δγ = |slope_HF − slope_LF| — absolute difference
       of the log-log slope d(log ρ_a)/d(log f) between the two bands.

    3. **Static shift**  δ = median(log10 ρ_a − AMA trend) — classic AMA
       shift estimate over the full frequency range.

    Classification:

    ===================  ===========================
    ``"clean"``          η ≤ ns_threshold, |δ| ≤ ss_threshold
    ``"static"``         η ≤ ns_threshold, |δ| > ss_threshold
    ``"near_surface"``   η > ns_threshold, |δ| ≤ ss_threshold
    ``"mixed"``          η > ns_threshold, |δ| > ss_threshold
    ===================  ===========================

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
        Any input accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    f_split : float, default=1.0
        Frequency boundary in Hz separating the HF (f ≥ f_split) from
        the LF (f < f_split) band.
    pband : (float, float) or None
        Period band ``(lo_s, hi_s)`` pre-filter applied before
        all computations.
    ns_threshold : float, default=2.0
        η > this → near-surface flag.
    ss_threshold : float, default=0.1
        |δ| > this (log10 units) → static-shift flag.
    sort_by : {"lon", "lat", "name"}, default="lon"
        Station ordering for the AMA spatial trend.
    half_window : int, default=3
        Number of neighbouring stations each side in the AMA trend.
    weights : {"tri", "gauss", "uniform"}, default="tri"
        Spatial weighting for the AMA trend.
    max_skew : float or None, default=6.0
        Phase-tensor skew ceiling; data above this are excluded.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    pandas.DataFrame
        One row per station with columns:
        ``station``, ``n_hf``, ``n_lf``,
        ``sigma_hf``, ``sigma_lf``, ``ns_index``,
        ``slope_hf``, ``slope_lf``, ``gradient_delta``,
        ``ss_delta_log10``, ``ns_flag``, ``ss_flag``, ``distortion_type``.

    References
    ----------
    Lei et al. (2017), "The non-static effect of near-surface
    inhomogeneity on CSAMT data", *Geophysics*.
    """
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    pt = build_phase_tensor_table(
        S, recursive=False, on_dup=on_dup,
        strict=False, verbose=verbose,
    )

    items = _order_sites(S, sort_by=sort_by)
    ST: List[str] = []
    FR: List[np.ndarray] = []
    LR: List[np.ndarray] = []

    for i, ed in enumerate(items):
        ed = _unwrap_ns(ed)
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        rho = _rho_det_from_z(z, fr)
        per = 1.0 / fr
        m = np.isfinite(rho)
        if pband is not None:
            lo, hi = pband
            m &= (per >= lo) & (per <= hi)
        if max_skew is not None and not pt.empty:
            sdf = pt[pt["station"] == st]
            if not sdf.empty:
                p_ref = sdf["period"].to_numpy()
                sk = np.abs(sdf["skew"].to_numpy())
                idx = _nearest_idx(1.0 / fr, p_ref)
                m &= sk[idx] <= float(max_skew)
        fr1 = fr[m]
        lr1 = np.log10(np.maximum(rho[m], 1e-24))
        if fr1.size == 0:
            continue
        ST.append(st)
        FR.append(fr1)
        LR.append(lr1)

    if not ST:
        df = pd.DataFrame(columns=_NS_COLS)

        return maybe_wrap_frame(
            df,
            api=api,
            name="near_surface_detection",
            kind="emtools.ss.near_surface",
            source=sites,
        )

    residuals = _ama_residuals_ns(FR, LR,
                                  half_window=half_window, weights=weights)
    rows = []
    for i, (st, fr, lr, delta) in enumerate(zip(ST, FR, LR, residuals)):
        hf = fr >= f_split
        lf = ~hf

        d_hf = delta[hf]
        d_lf = delta[lf]

        σ_hf = float(np.nanstd(d_hf)) if d_hf.size >= 2 else float("nan")
        σ_lf = float(np.nanstd(d_lf)) if d_lf.size >= 2 else float("nan")
        η = (σ_hf / (σ_lf + 1e-12)
             if np.isfinite(σ_hf) and np.isfinite(σ_lf)
             else float("nan"))

        slope_hf = _log_slope(
            np.log10(np.maximum(fr[hf], 1e-24)), lr[hf]
        )
        slope_lf = _log_slope(
            np.log10(np.maximum(fr[lf], 1e-24)), lr[lf]
        )
        grad_delta = (abs(slope_hf - slope_lf)
                      if np.isfinite(slope_hf) and np.isfinite(slope_lf)
                      else float("nan"))

        fin = delta[np.isfinite(delta)]
        ss_delta = float(np.nanmedian(fin)) if fin.size else float("nan")

        ns_flag = bool(np.isfinite(η) and η > ns_threshold)
        ss_flag = bool(np.isfinite(ss_delta) and abs(ss_delta) > ss_threshold)

        dtype = (
            "mixed"        if ns_flag and ss_flag else
            "near_surface" if ns_flag else
            "static"       if ss_flag else
            "clean"
        )

        rows.append({
            "station":        st,
            "n_hf":           int(hf.sum()),
            "n_lf":           int(lf.sum()),
            "sigma_hf":       σ_hf,
            "sigma_lf":       σ_lf,
            "ns_index":       float(η)          if np.isfinite(η)          else float("nan"),
            "slope_hf":       float(slope_hf)   if np.isfinite(slope_hf)   else float("nan"),
            "slope_lf":       float(slope_lf)   if np.isfinite(slope_lf)   else float("nan"),
            "gradient_delta": float(grad_delta) if np.isfinite(grad_delta) else float("nan"),
            "ss_delta_log10": ss_delta,
            "ns_flag":        ns_flag,
            "ss_flag":        ss_flag,
            "distortion_type": dtype,
        })

    df = pd.DataFrame(rows, columns=_NS_COLS)

    return maybe_wrap_frame(
        df,
        api=api,
        name="near_surface_detection",
        kind="emtools.ss.near_surface",
        source=sites,
        description="Near-surface and static-shift distortion diagnostics.",
    )


def plot_ns_detection(
    sites: Any,
    *,
    f_split: float = 1.0,
    pband: Optional[Tuple[float, float]] = None,
    ns_threshold: float = 2.0,
    ss_threshold: float = 0.1,
    sort_by: str = "lon",
    half_window: int = 3,
    weights: str = "tri",
    max_skew: Optional[float] = 6.0,
    show_ss: bool = True,
    figsize: Tuple[float, float] = (9.0, 4.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Bar chart of the NS index per station, colored by distortion type.

    Each bar height is η = σ_HF / σ_LF.  A dashed line marks
    *ns_threshold*.  An optional secondary y-axis shows the static-shift
    estimate δ (log10 units) as a stem plot.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
    f_split : float, default=1.0
        HF/LF split frequency in Hz.
    pband : (float, float) or None
    ns_threshold, ss_threshold : float
    sort_by : {"lon", "lat", "name"}
    half_window, weights, max_skew
        Forwarded to :func:`detect_near_surface`.
    show_ss : bool, default=True
        If True and ax has room, overlay static-shift δ as a grey
        stem plot on a secondary y-axis.
    figsize : (float, float), default=(9, 4.5)
    recursive, on_dup, strict, verbose
        Forwarded to :func:`~pycsamt.emtools._core.ensure_sites`.
    ax : matplotlib.axes.Axes, optional
        Draw on existing axes.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.patches as mpatches

    df = detect_near_surface(
        sites,
        f_split=f_split,
        pband=pband,
        ns_threshold=ns_threshold,
        ss_threshold=ss_threshold,
        sort_by=sort_by,
        half_window=half_window,
        weights=weights,
        max_skew=max_skew,
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

    x = np.arange(len(df))
    colors = [_TYPE_COLORS[t] for t in df["distortion_type"]]
    ax.bar(x, df["ns_index"].fillna(0).values,
           color=colors, width=0.7, edgecolor="0.3", linewidth=0.5)
    ax.axhline(ns_threshold, color="k", lw=1.2, ls="--")

    # secondary axis for static-shift δ
    if show_ss and "ss_delta_log10" in df.columns:
        ax2 = ax.twinx()
        ax2.stem(
            x,
            df["ss_delta_log10"].fillna(0).values,
            linefmt="0.55",
            markerfmt="D",
            basefmt="none",
        )
        ax2.axhline(0, color="0.6", lw=0.7, ls=":")
        ax2.set_ylabel("δ (log10 ρ_a shift)", fontsize=8, color="0.4")
        ax2.tick_params(axis="y", labelcolor="0.4")

    ax.set_xticks(x)
    ax.set_xticklabels(df["station"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("NS index  η = σ_HF / σ_LF")
    ax.set_xlabel("Station")
    ax.set_title(
        f"Near-surface effect detection  "
        f"(f_split = {f_split} Hz,  η threshold = {ns_threshold})"
    )
    ax.grid(axis="y", alpha=0.25)

    present = df["distortion_type"].unique()
    patches = [
        mpatches.Patch(
            color=_TYPE_COLORS[k],
            label=k.replace("_", " ").title(),
        )
        for k in ("clean", "static", "near_surface", "mixed")
        if k in present
    ]
    patches.append(
        plt.Line2D([0], [0], color="k", ls="--", lw=1.2,
                   label=f"η threshold ({ns_threshold})")
    )
    ax.legend(handles=patches, fontsize=8, loc="upper right", framealpha=0.85)
    return ax
