
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ._core import (
    ensure_sites,
    _iter_items,
    _apply_each,
    _get_z_block,
    _name,
)
from .tensor import build_phase_tensor_table


def _rho_det_from_z(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    # ρa_det ≈ sqrt(ρ_xy ρ_yx) ; ρ = 0.2|Z|^2/f
    zx = z[:, 0, 1]
    zy = z[:, 1, 0]
    rx = 0.2 * (np.abs(zx) ** 2) / (fr + 1e-24)
    ry = 0.2 * (np.abs(zy) ** 2) / (fr + 1e-24)
    return np.sqrt(rx * ry)


def _site_order_key(ed: Any, key: str) -> Tuple[int, float, str]:
    # sorting key for along-line order
    st = _name(ed, 0)
    if key == "lon":
        x = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
        return (0, float(x)) if x is not None else (1, np.inf, st)
    if key == "lat":
        y = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
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
) -> pd.DataFrame:
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    items = _order_sites(S, sort_by=sort_by)
    n = len(items)
    if n == 0:
        return pd.DataFrame(
            columns=["station", "delta_log10_rho", "fac_rho",
                     "fac_z", "n_used"]
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
        return pd.DataFrame(
            columns=["station", "delta_log10_rho", "fac_rho",
                     "fac_z", "n_used"]
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
        if robust_overall == "mean":
            delta = float(np.nanmean(d))
        else:
            delta = float(np.nanmedian(d))
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
    return tbl.sort_values("station").reset_index(drop=True)


def _scale_site_Z(ed: Any, s: float) -> None:
    Z, z, fr = _get_z_block(ed)
    if Z is None:
        return
    try:
        Z.z = z * float(s)
    except Exception:
        pass
    # scale errors if present
    try:
        ze = getattr(Z, "z_err", None)
        if isinstance(ze, np.ndarray) and ze.shape == z.shape:
            setattr(Z, "z_err", ze * float(s))
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
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
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
) -> pd.DataFrame:
    ST, FR, LR = _prep_lr_curves(
        sites, pband=pband, max_skew=max_skew,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if not ST:
        return pd.DataFrame(
            columns=["station", "delta_log10_rho",
                     "fac_rho", "fac_z", "n_used"]
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
    return pd.DataFrame.from_records(rows).sort_values(
        "station"
    ).reset_index(drop=True)


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
) -> pd.DataFrame:
    ST, FR, LR = _prep_lr_curves(
        sites, pband=pband, max_skew=max_skew,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if not ST:
        return pd.DataFrame(
            columns=["station", "delta_log10_rho",
                     "fac_rho", "fac_z", "n_used"]
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
    return pd.DataFrame.from_records(rows).sort_values(
        "station"
    ).reset_index(drop=True)


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
) -> pd.DataFrame:
    ST, FR, LR = _prep_lr_curves(
        sites, pband=pband, max_skew=max_skew,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if not ST:
        return pd.DataFrame(
            columns=["station", "delta_log10_rho",
                     "fac_rho", "fac_z", "n_used"]
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
    return pd.DataFrame.from_records(rows).sort_values(
        "station"
    ).reset_index(drop=True)

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
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)" if axis_y == "logperiod"
                  else "Period (s)")
    ax.set_xticks(np.arange(len(labs)))
    ax.set_xticklabels(labs, rotation=90)
    yt = np.linspace(0, Z.shape[0] - 1, num=min(8, Z.shape[0]))
    yv = np.linspace(yg.min(), yg.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
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
    ax.set_xscale("log")
    ax.plot(pb[mb], rb[mb], "o-", lw=1.4, label="before")
    ax.plot(pa[ma], ra[ma], "s-", lw=1.4, label="after")
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
    ax.set_xlabel("Station")
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=90)
    return ax

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
    colors: Tuple[str, str] = ("#1f77b4", "#d62728"),
    marker: str = "o",
    ms: float = 2.5,
    lw: float = 1.2,
    ls: str = "-",
    figsize: Tuple[float, float] = (4.8, 4.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    eps: float = 1e-24, 
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    # pick station
    sel = {}
    for i, ed in enumerate(_iter_items(S)):
        sel[_name(ed, i)] = ed
    if not sel:
        if ax is None: _, ax = plt.subplots(subplot_kw={"polar": True})
        ax.text(0.5, 0.5, "no sites", ha="center", va="center")
        return ax
    if station is None:
        station = sorted(sel.keys())[0]
    ed = sel.get(station, None)
    if ed is None:
        if ax is None: _, ax = plt.subplots(subplot_kw={"polar": True})
        ax.text(0.5, 0.5, "station not found",
                ha="center", va="center")
        return ax
    Z, z, fr = _get_z_block(ed)
    if Z is None:
        if ax is None: _, ax = plt.subplots(subplot_kw={"polar": True})
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

    if ax is None:
        _, ax = plt.subplots(
            figsize=figsize, subplot_kw={"polar": True}
        )
    # set polar style: 0 at north, CW
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    # plot
    ax.plot(th, r1, ls=ls, lw=lw, marker=marker,
            ms=ms, color=colors[0], label="ρa_xy (rot)")
    ax.plot(th, r2, ls=ls, lw=lw, marker=marker,
            ms=ms, color=colors[1], label="ρa_yx (rot)")
    if fill_between:
        lo = np.minimum(r1, r2); hi = np.maximum(r1, r2)
        ax.fill_between(th, lo, hi, color="0.5", alpha=0.10)

    ax.grid(True, alpha=0.25)
    ax.set_title(str(station), pad=10)
    if radial == "rho":
        ax.set_rlabel_position(135)
        ax.set_ylabel("ρa (Ω·m)")
    else:
        ax.set_rlabel_position(135)
        ax.set_ylabel("log10 ρa")
    ax.legend(loc="lower left", bbox_to_anchor=(0.02, 0.02),
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
) -> pd.DataFrame:
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
        return pd.DataFrame(columns=_NS_COLS)

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

    return pd.DataFrame(rows, columns=_NS_COLS)


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
