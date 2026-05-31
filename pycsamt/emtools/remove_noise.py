from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import (
    ensure_sites,
    _apply_each,
    _iter_items,
    _get_z_block,
    _get_t_block,
    _name,
    _station_positions,
)

_MU0: float = 4.0 * np.pi * 1e-7   # H/m

# ----------------------------- SNR table -------------------------------- #

def snr_table(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    rows = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        out = _get_z_block(ed, with_errors=True)
        if len(out) == 4:
            Z, z, fr, ze = out
        else:
            Z, z, fr, ze = (out + (None,))[:4]
        if Z is None:
            continue
        a = np.sqrt(np.nanmean(np.abs(z) ** 2, axis=(1, 2)))
        if isinstance(ze, np.ndarray) and ze.shape == z.shape:
            e = np.sqrt(np.nanmean(np.abs(ze) ** 2, axis=(1, 2)))
        else:
            e = np.full_like(a, np.nan)
        snr = a / (e + 1e-12)
        for f, s in zip(fr, snr):
            rows.append(dict(station=st, freq=float(f), snr=float(s)))
    return pd.DataFrame.from_records(rows)


# --------------------------- power-line notching ------------------------- #

def _harm_mask(fr: np.ndarray, mains: float, n_harm: int,
               tol_hz: float) -> np.ndarray:
    kk = np.arange(1, int(n_harm) + 1, dtype=float)
    fH = kk * float(mains)
    m = np.zeros(fr.size, dtype=bool)
    for fh in fH:
        m |= np.abs(fr - fh) <= float(tol_hz)
    return m


def _interp_rows(
    y: np.ndarray,
    good: np.ndarray,
) -> np.ndarray:
    # y: (n, ...) complex; interp along axis 0
    if y.ndim == 1:
        y = y[:, None]
    x = np.arange(y.shape[0])
    yi = y.copy()
    g = good & np.all(np.isfinite(y), axis=tuple(range(1, y.ndim)))
    if g.sum() < 2:
        yi[~g] = np.nan
        return yi.squeeze()
    xi = x[g]
    for idx in np.argwhere(~g).ravel():
        # nearest two good neighbors (linear)
        j = np.searchsorted(xi, idx)
        j0 = max(0, min(j - 1, xi.size - 1))
        j1 = max(0, min(j, xi.size - 1))
        a, b = xi[j0], xi[j1]
        if a == b:
            yi[idx] = y[a]
        else:
            t = (idx - a) / (b - a)
            yi[idx] = (1 - t) * y[a] + t * y[b]
    return yi.squeeze()


def notch_powerline(
    sites: Any,
    *,
    mains_hz: float = 50.0,   # 50 or 60
    n_harm: int = 30,
    tol_hz: float = 0.08,     # Hz window around each harmonic
    mode: str = "interp",     # mask|interp
    also: str = "both",       # z|tipper|both
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        mH = _harm_mask(fr, mains_hz, n_harm, tol_hz)
        z2 = z.copy()
        if mode == "mask":
            z2[mH] = np.nan
        else:
            z2[mH] = np.nan  # mark, then fill
            good = ~np.isnan(z2).any(axis=(1, 2))
            z2 = _interp_rows(z2, good)
        Z.z = z2
        if also in ("tipper", "both"):
            T, t, ft = _get_t_block(ed)
            if T is not None and t is not None:
                mt = _harm_mask(ft, mains_hz, n_harm, tol_hz)
                t2 = t.copy()
                if mode == "mask":
                    t2[mt] = np.nan
                else:
                    t2[mt] = np.nan
                    good = ~np.isnan(t2).any(axis=1)
                    t2 = _interp_rows(t2, good)
                T.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ------------------------- log-frequency smoothing ----------------------- #

def _smooth1d(y: np.ndarray, win: int, kind: str) -> np.ndarray:
    if win <= 1:
        return y
    w = np.ones(int(win), dtype=float)
    if kind == "tri":
        w = np.convolve(w, w, mode="full")
    w = w / np.sum(w)
    if y.ndim == 1:
        yy = np.convolve(y, w, mode="same")
        return yy
    # apply per column
    out = np.empty_like(y)
    for j in range(y.shape[1]):
        out[:, j] = np.convolve(y[:, j], w, mode="same")
    return out


def smooth_logfreq(
    sites: Any,
    *,
    win: int = 5,
    kind: str = "tri",         # box|tri
    also: str = "both",        # z|tipper|both
    gate_snr: Optional[float] = None,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    ST = snr_table(S) if gate_snr is not None else None

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        y = z.reshape(z.shape[0], -1)
        ok = np.isfinite(y).all(axis=1)
        if gate_snr is not None and not ST.empty:
            st = _name(ed, 0)
            sdf = ST[ST["station"] == st]
            if not sdf.empty:
                idx = np.searchsorted(sdf["freq"].to_numpy(), fr)
                idx = np.clip(idx, 0, len(sdf) - 1)
                ok &= (sdf["snr"].to_numpy()[idx] >= gate_snr)
        y2 = y.copy()
        y2[ok] = _smooth1d(y[ok], win, kind)
        Z.z = y2.reshape(z.shape)
        if also in ("tipper", "both"):
            T, t, ft = _get_t_block(ed)
            if T is not None and t is not None:
                yt = t
                okt = np.isfinite(yt).all(axis=1)
                T.tipper = _smooth1d(yt, win, kind) if okt.any() else yt
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --------------------- group-trend shrinkage (robust) -------------------- #

def _build_group_trend(
    sites: Any,
    groups: Dict[str, List[str]],
) -> Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    # trend[grp][station_ref] = (fr_union, z_med_on_union)
    trend: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}
    # union frequency per group
    for g, sts in groups.items():
        G: List[float] = []
        pool = []
        for i, ed in enumerate(_iter_items(sites)):
            st = _name(ed, i)
            if st not in sts:
                continue
            Z, z, fr = _get_z_block(ed)
            if Z is None:
                continue
            G.append(fr)
            pool.append((st, fr, z))
        if not pool:
            continue
        Gu = np.unique(np.concatenate(G))
        Zm = []
        for k in range(Gu.size):
            vals = []
            for st, fr, z in pool:
                j = np.searchsorted(fr, Gu[k])
                j = np.clip(j, 0, fr.size - 1)
                vals.append(z[j])
            vals = np.asarray(vals)
            Zm.append(np.nanmedian(vals, axis=0))
        trend[g] = {"_union": (Gu, np.asarray(Zm))}
    return trend


def shrink_to_group_trend(
    sites: Any,
    *,
    groups: Optional[Dict[str, List[str]]] = None,
    group_key: Optional[str] = None,
    lam: float = 0.25,          # 0..1; 0=no change, 1=trend
    gate_harm: bool = True,
    mains_hz: float = 50.0,
    n_harm: int = 30,
    tol_hz: float = 0.08,
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if groups is None:
        # auto: single group with all stations
        names = [_name(ed, i) for i, ed in enumerate(_iter_items(S))]
        groups = {"ALL": names}
    Tref = _build_group_trend(S, groups)

    def _one(Si):
        ed = next(_iter_items(Si))
        st = _name(ed, 0)
        # find group containing st
        g = None
        for k, v in groups.items():
            if st in v:
                g = k; break
        if g is None or g not in Tref:
            return Si
        (Gu, Zu) = Tref[g]["_union"]
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        idx = np.searchsorted(Gu, fr)
        idx = np.clip(idx, 0, Gu.size - 1)
        tr = Zu[idx]
        if gate_harm:
            mH = _harm_mask(fr, mains_hz, n_harm, tol_hz)
        else:
            mH = np.ones(fr.size, dtype=bool)
        z2 = z.copy()
        z2[mH] = (1.0 - lam) * z[mH] + lam * tr[mH]
        Z.z = z2
        if also in ("tipper", "both"):
            T, t, ft = _get_t_block(ed)
            if T is not None and t is not None:
                # simple shrink on tipper to its group median
                tt = t.copy()
                T.tipper = (1.0 - lam) * t + lam * tt
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ------------------------------ pipeline --------------------------------- #

def remove_noise_pipeline(
    sites: Any,
    *,
    mains_hz: float = 50.0,
    n_harm: int = 30,
    tol_hz: float = 0.08,
    notch_mode: str = "interp",
    smooth_win: int = 5,
    smooth_kind: str = "tri",
    gate_snr: Optional[float] = 2.5,
    group_shrink: bool = False,
    shrink_lam: float = 0.25,
    groups: Optional[Dict[str, List[str]]] = None,
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = notch_powerline(
        sites,
        mains_hz=mains_hz,
        n_harm=n_harm,
        tol_hz=tol_hz,
        mode=notch_mode,
        also=also,
        inplace=False,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    S = smooth_logfreq(
        S,
        win=smooth_win,
        kind=smooth_kind,
        also=also,
        gate_snr=gate_snr,
        inplace=False,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if group_shrink:
        S = shrink_to_group_trend(
            S,
            groups=groups,
            lam=shrink_lam,
            mains_hz=mains_hz,
            n_harm=n_harm,
            tol_hz=tol_hz,
            also=also,
            inplace=False,
            recursive=False,
            on_dup=on_dup,
            strict=False,
            verbose=verbose,
        )
    if inplace:
        # copy back into input sites
        _ = _apply_each(
            sites, lambda Si: Si, inplace=True, verbose=verbose
        )
        return sites
    return S

# ------------------------- 1) Hampel outlier filter --------------------- #

def _hampel_1d(y: np.ndarray, k: int, nsig: float) -> np.ndarray:
    if y.size == 0 or k <= 0:
        return y
    y2 = y.copy()
    n = y.size
    for i in range(n):
        lo = max(0, i - k); hi = min(n, i + k + 1)
        win = y[lo:hi]
        med = np.nanmedian(win)
        mad = np.nanmedian(np.abs(win - med)) + 1e-12
        if np.abs(y[i] - med) > nsig * 1.4826 * mad:
            y2[i] = med
    return y2


def hampel_filter_freq(
    sites: Any,
    *,
    win: int = 3,
    nsig: float = 3.0,
    on: str = "both",          # z|tipper|both
    domain: str = "reim",      # reim|magphase
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is not None:
            Y = z.reshape(z.shape[0], -1)
            if domain == "magphase":
                m = np.abs(Y); p = np.angle(Y)
                m2 = np.vstack([
                    _hampel_1d(m[:, j], win, nsig)
                    for j in range(m.shape[1])
                ]).T
                Y2 = m2 * np.exp(1j * p)
            else:
                R = np.real(Y); I = np.imag(Y)
                R2 = np.vstack([
                    _hampel_1d(R[:, j], win, nsig)
                    for j in range(R.shape[1])
                ]).T
                I2 = np.vstack([
                    _hampel_1d(I[:, j], win, nsig)
                    for j in range(I.shape[1])
                ]).T
                Y2 = R2 + 1j * I2
            Z.z = Y2.reshape(z.shape)
        if on in ("tipper", "both"):
            from ._core import _get_t_block
            T, t, ft = _get_t_block(ed)
            if T is not None:
                r = np.real(t); im = np.imag(t)
                r2 = np.vstack([
                    _hampel_1d(r[:, j], win, nsig)
                    for j in range(r.shape[1])
                ]).T
                i2 = np.vstack([
                    _hampel_1d(im[:, j], win, nsig)
                    for j in range(im.shape[1])
                ]).T
                T.tipper = r2 + 1j * i2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# -------- 2) Spatial median smooth (across stations at fixed freq) ------ #

def spatial_median_filter(
    sites: Any,
    *,
    half_window: int = 2,
    lam: float = 0.25,          # 0..1 shrink toward local median
    on: str = "z",              # z|tipper|both
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    items = list(_iter_items(S))
    n = len(items)

    def _nbrs(i: int) -> List[int]:
        lo = max(0, i - half_window)
        hi = min(n - 1, i + half_window)
        ids = list(range(lo, hi + 1))
        if i in ids:
            ids.remove(i)
        return ids

    # prefetch Z-blocks
    ZB = []
    for ed in items:
        ZB.append(_get_z_block(ed))

    def _one(Si):
        ed = next(_iter_items(Si))
        i = items.index(ed)
        Z, z, fr = ZB[i]
        if Z is not None:
            z2 = z.copy()
            for k in range(z.shape[0]):
                pool = []
                for j in _nbrs(i):
                    Zj, zj, frj = ZB[j]
                    if Zj is None:
                        continue
                    jj = np.clip(
                        np.searchsorted(frj, fr[k]), 0, frj.size - 1
                    )
                    pool.append(zj[jj])
                if pool:
                    med = np.nanmedian(np.asarray(pool), axis=0)
                    z2[k] = (1.0 - lam) * z[k] + lam * med
            Z.z = z2
        if on in ("tipper", "both"):
            T, t, ft = _get_t_block(ed)
            if T is not None:
                # simple neighbor median on tipper too
                t2 = t.copy()
                for k in range(t.shape[0]):
                    pool = []
                    for j in _nbrs(i):
                        Tj, tj, ftj = None, None, None
                        try:
                            Tj, tj, ftj = _get_t_block(items[j])
                        except Exception:
                            pass
                        if Tj is None:
                            continue
                        jj = np.clip(
                            np.searchsorted(ftj, ft[k]),
                            0, ftj.size - 1
                        )
                        pool.append(tj[jj])
                    if pool:
                        med = np.nanmedian(np.asarray(pool), axis=0)
                        t2[k] = (1.0 - lam) * t[k] + lam * med
                T.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --------- 3) Low-rank denoise (profile RPCA-ish on |Z_xy|,|Z_yx|) ------ #

def _svd_rank_k(M: np.ndarray, r: int) -> np.ndarray:
    try:
        U, s, Vh = np.linalg.svd(M, full_matrices=False)
    except np.linalg.LinAlgError:
        U, s, Vh = np.linalg.svd(
            M + 1e-12 * np.random.randn(*M.shape),
            full_matrices=False,
        )
    r = int(max(1, min(r, min(M.shape))))
    return (U[:, :r] * s[:r]) @ Vh[:r, :]


def rpca_offdiag_denoise(
    sites: Any,
    *,
    rank: int = 2,
    keep_phase: bool = True,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    items = list(_iter_items(S))
    # build matrix M(station, freq) from log |Z_xy| median(|Z_xy|,|Z_yx|)
    ST, FR, MAG = [], [], []
    for i, ed in enumerate(items):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        m = np.nanmedian(
            np.stack([np.abs(z[:, 0, 1]), np.abs(z[:, 1, 0])], axis=1),
            axis=1,
        )
        ST.append(i); FR.append(fr); MAG.append(np.log10(m + 1e-24))
    if not ST:
        return S
    G = np.unique(np.concatenate(FR))
    M = np.full((len(ST), G.size), np.nan, dtype=float)
    for row, (fr, lg) in enumerate(zip(FR, MAG)):
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        M[row, idx] = lg
    # simple fill of gaps by nearest valid along freq
    for i in range(M.shape[0]):
        r = M[i]
        g = np.isfinite(r)
        if g.sum() >= 2:
            xi = np.where(g)[0]
            for j in np.where(~g)[0]:
                k = np.searchsorted(xi, j)
                k0 = max(0, min(k - 1, xi.size - 1))
                k1 = max(0, min(k, xi.size - 1))
                a, b = xi[k0], xi[k1]
                t = 0.0 if a == b else (j - a) / (b - a)
                r[j] = (1 - t) * r[a] + t * r[b]
            M[i] = r
    L = _svd_rank_k(M, rank=rank)
    # push back magnitudes; keep phase if requested
    def _one(Si):
        ed = next(_iter_items(Si))
        i = items.index(ed)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        mag_clean = 10 ** L[i, idx]
        # per off-diagonal component scale
        for (a, b) in [(0, 1), (1, 0)]:
            comp = z[:, a, b]
            if keep_phase:
                ph = np.angle(comp)
                Z.z[:, a, b] = mag_clean * np.exp(1j * ph)
            else:
                # scale original by ratio on magnitude
                r = mag_clean / (np.abs(comp) + 1e-24)
                Z.z[:, a, b] = comp * r
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# -------- 4) Off-diagonal consistency (sym/anti-sym blend) -------------- #

def enforce_offdiag_consistency(
    sites: Any,
    *,
    mode: str = "anti",        # anti|sym
    lam: float = 0.5,          # blend toward constraint
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            return Si
        xy = z[:, 0, 1]; yx = z[:, 1, 0]
        if mode == "sym":
            target_xy = 0.5 * (xy + yx)
            target_yx = target_xy
        else:
            target_xy = 0.5 * (xy - yx)
            target_yx = -target_xy
        z[:, 0, 1] = (1.0 - lam) * xy + lam * target_xy
        z[:, 1, 0] = (1.0 - lam) * yx + lam * target_yx
        Z.z = z
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ------------- 5) Incoherent-frequency vote mask (SNR-based) ------------ #

def mask_incoherent_freqs(
    sites: Any,
    *,
    snr_thresh: float = 2.5,
    min_frac: float = 0.4,     # keep if ≥ this fraction pass
    also: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    T = snr_table(S)
    if T.empty:
        return S
    # union freq grid
    G = np.unique(T["freq"].to_numpy(dtype=float))
    # for each freq, fraction of stations with SNR ≥ thresh
    frac = np.zeros(G.size, dtype=float)
    for k, f in enumerate(G):
        m = T["freq"] == f
        sn = T.loc[m, "snr"].to_numpy(dtype=float)
        frac[k] = np.nanmean(sn >= snr_thresh)
    keepG = frac >= float(min_frac)

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is not None:
            idx = np.searchsorted(G, fr)
            idx = np.clip(idx, 0, G.size - 1)
            keep = keepG[idx]
            z2 = z.copy(); z2[~keep] = np.nan
            Z.z = z2
        if also in ("tipper", "both"):
            Tt, t, ft = _get_t_block(ed)
            if Tt is not None:
                idx = np.searchsorted(G, ft)
                idx = np.clip(idx, 0, G.size - 1)
                keep = keepG[idx]
                t2 = t.copy(); t2[~keep] = np.nan
                Tt.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --------- 6) Static-shift removal — Torres-Verdín & Bostick (1992) ----- #

def _hanning_weights(dx: np.ndarray, w_H: float) -> np.ndarray:
    """Hanning (von Hann) spatial filter weights (Torres-Verdín & Bostick 1992).

    h(x) = (1 + cos(2π x / W_H)) / W_H  for |x| ≤ W_H/2, else 0.
    W_H is the full window width [m].
    """
    half = w_H / 2.0
    return np.where(
        np.abs(dx) <= half,
        np.maximum((1.0 + np.cos(2.0 * np.pi * dx / w_H)) / w_H, 0.0),
        0.0,
    )


def correct_static_shift(
    sites: Any,
    *,
    window_m: float = 1500.0,
    spacing_m: float = 200.0,
    comp: str = "det",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    """
    Remove static shift via Hanning adaptive moving-average (AMA) spatial filter.

    Implements the Torres-Verdín & Bostick (1992) approach used in
    Kouadio et al. (2024):

    1. For each frequency, build the spatial log(ρ_a) profile across stations.
    2. Apply a Hanning low-pass spatial filter with full-width ``window_m``.
    3. The static-shift correction factor at station *i* is
       ``C_i = sqrt(ρ_smooth_i / ρ_obs_i)`` (in log space:
       ``log C = 0.5 (log ρ_smooth − log ρ_obs)``).
    4. Update every Z component: ``Z_corrected = Z × C``.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
        Input sites.
    window_m : float
        Full Hanning window width [m] (Torres-Verdín W_H).  Stations
        further than ``window_m/2`` contribute zero weight.
    spacing_m : float
        Fallback station spacing [m] used when EDI metadata carries no
        coordinate information.
    comp : {"det", "xy", "yx"}
        Apparent-resistivity component used to estimate the static shift.
        ``"det"`` uses the arithmetic mean of |Z_xy|² and |Z_yx|².
    inplace : bool
        If *True*, modify the input sites in place.  If *False* (default),
        return a new Sites object with corrected Z tensors.
    recursive, on_dup, strict, verbose
        Passed to :func:`ensure_sites`.

    Returns
    -------
    Sites
        Sites with static-shift-corrected Z tensors (when ``inplace=False``).
    """
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )

    items = list(_iter_items(S))
    if not items:
        return S

    names     = [_name(ed, i) for i, ed in enumerate(items)]
    positions = _station_positions(items, spacing_m)
    N         = len(items)

    # --- collect ρ_a on native frequency grids ---------------------------
    rho_data: List[Tuple[Optional[np.ndarray], Optional[np.ndarray]]] = []
    all_freqs: List[np.ndarray] = []

    for ed in items:
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            rho_data.append((None, None))
            continue
        omega = 2.0 * np.pi * np.maximum(fr, 1e-24)
        if comp == "xy":
            mag2 = np.abs(z[:, 0, 1]) ** 2
        elif comp == "yx":
            mag2 = np.abs(z[:, 1, 0]) ** 2
        else:   # det
            mag2 = 0.5 * (np.abs(z[:, 0, 1]) ** 2 + np.abs(z[:, 1, 0]) ** 2)
        rho = mag2 / (omega * _MU0)
        lr  = np.where(rho > 0.0, np.log(rho), np.nan)
        rho_data.append((fr, lr))
        all_freqs.append(fr)

    if not all_freqs:
        return S

    G = np.unique(np.concatenate(all_freqs))
    F = G.size

    # --- map each station's log(ρ_a) onto the union frequency grid --------
    log_rho_mat = np.full((N, F), np.nan, dtype=float)
    for i, (fr_i, lr_i) in enumerate(rho_data):
        if fr_i is None:
            continue
        idx = np.searchsorted(G, fr_i)
        idx = np.clip(idx, 0, F - 1)
        log_rho_mat[i, idx] = lr_i

    # --- vectorised Hanning spatial smooth --------------------------------
    dx_mat = positions[:, None] - positions[None, :]   # [N, N]
    w_mat  = _hanning_weights(dx_mat, window_m)         # [N, N]

    valid = np.isfinite(log_rho_mat).astype(float)     # [N, F]
    # broadcast: w_mat[N,N,1] × valid[1,N,F] × log_rho[1,N,F]
    w3    = w_mat[:, :, None]
    lr3   = log_rho_mat[None, :, :]
    v3    = valid[None, :, :]
    num   = np.sum(w3 * v3 * lr3, axis=1)              # [N, F]
    denom = np.sum(w3 * v3,       axis=1)               # [N, F]
    log_rho_smooth = np.where(denom > 1e-30, num / denom, log_rho_mat)

    # --- log correction factor C = sqrt(ρ_smooth / ρ_obs) ----------------
    log_C = np.where(
        np.isfinite(log_rho_mat) & np.isfinite(log_rho_smooth),
        0.5 * (log_rho_smooth - log_rho_mat),
        0.0,
    )
    corr_map = {names[i]: (G, np.exp(log_C[i])) for i in range(N)}

    # --- apply per station -----------------------------------------------
    def _one(Si):
        ed = next(_iter_items(Si))
        st = _name(ed, 0)
        Z, z, fr = _get_z_block(ed)
        if Z is None or st not in corr_map:
            return Si
        G_c, C_arr = corr_map[st]
        idx = np.searchsorted(G_c, fr)
        idx = np.clip(idx, 0, G_c.size - 1)
        C   = C_arr[idx]
        Z.z = z * C[:, None, None]
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# -------------------- NOISE REMOVAL QC PLOTS ---------------------------- #

def _offdiag_logmag(z: np.ndarray) -> np.ndarray:
    m = np.nanmedian(
        np.stack([np.abs(z[:, 0, 1]), np.abs(z[:, 1, 0])], axis=1),
        axis=1,
    )
    return np.log10(np.maximum(m, 1e-24))


def _union_offdiag_matrix(
    sites: Any,
) -> Tuple[List[str], np.ndarray, np.ndarray]:
    S = ensure_sites(sites, recursive=False, strict=False)
    sts, Gs, Ms = [], [], []
    for i, ed in enumerate(_iter_items(S)):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        sts.append(_name(ed, i))
        Gs.append(fr)
        Ms.append(_offdiag_logmag(z))
    if not sts:
        return [], np.array([]), np.empty((0, 0))
    G = np.unique(np.concatenate(Gs))
    M = np.full((len(sts), G.size), np.nan, dtype=float)
    for row, (fr, lm) in enumerate(zip(Gs, Ms)):
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        M[row, idx] = lm
        # fill small gaps by nearest neighbors
        r = M[row]
        g = np.isfinite(r)
        if g.sum() >= 2:
            xi = np.where(g)[0]
            for j in np.where(~g)[0]:
                k = np.searchsorted(xi, j)
                k0 = max(0, min(k - 1, xi.size - 1))
                k1 = max(0, min(k, xi.size - 1))
                a, b = xi[k0], xi[k1]
                t = 0.0 if a == b else (j - a) / (b - a)
                r[j] = (1 - t) * r[a] + t * r[b]
            M[row] = r
    return sts, G, M


def _denoise_sites(
    sites: Any,
    method: str,
    **kws: Any,
):
    m = method.lower()
    if m == "pipeline":
        return remove_noise_pipeline(sites, **kws)
    if m == "notch":
        return notch_powerline(sites, **kws)
    if m == "smooth":
        return smooth_logfreq(sites, **kws)
    if m == "rpca":
        return rpca_offdiag_denoise(sites, **kws)
    if m == "spatial":
        return spatial_median_filter(sites, **kws)
    if m == "hampel":
        return hampel_filter_freq(sites, **kws)
    raise ValueError(f"unknown method: {method}")


# 1) Δ log|Z_off| pseudosection (station × log-period) ------------------- #

def nr_qc_delta_offdiag_psection(
    sites: Any,
    *,
    method: str = "pipeline",
    vlim: Optional[float] = None,
    figsize: Tuple[float, float] = (9.0, 4.8),
    ax: Optional[plt.Axes] = None,
    **denoise: Any,
) -> plt.Axes:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _denoise_sites(S0, method, inplace=False, **denoise)
    st0, G0, M0 = _union_offdiag_matrix(S0)
    st1, G1, M1 = _union_offdiag_matrix(S1)
    if not st0 or not st1:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    # align stations intersection, and G union
    labs = sorted(list(set(st0) & set(st1)))
    I0 = [st0.index(s) for s in labs]
    I1 = [st1.index(s) for s in labs]
    G = np.unique(np.concatenate([G0, G1]))
    def _resample(M, Gs):
        out = np.full((M.shape[0], G.size), np.nan)
        for i in range(M.shape[0]):
            idx = np.searchsorted(G, Gs)
            idx = np.clip(idx, 0, G.size - 1)
            out[i, idx] = M[i]
            # fill as above
            r = out[i]; g = np.isfinite(r)
            if g.sum() >= 2:
                xi = np.where(g)[0]
                for j in np.where(~g)[0]:
                    k = np.searchsorted(xi, j)
                    k0 = max(0, min(k - 1, xi.size - 1))
                    k1 = max(0, min(k, xi.size - 1))
                    a, b = xi[k0], xi[k1]
                    t = 0 if a == b else (j - a) / (b - a)
                    r[j] = (1 - t) * r[a] + t * r[b]
        return out
    A = _resample(M0[I0], G0)
    B = _resample(M1[I1], G1)
    D = (B - A).T  # (freq, station)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    v = D[np.isfinite(D)]
    if vlim is None and v.size:
        vlim = float(max(0.1, np.nanpercentile(np.abs(v), 95)))
    lp = np.log10(1.0 / G)
    im = ax.imshow(
        D,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-(vlim or 0.5),
        vmax=(vlim or 0.5),
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(len(labs)))
    ax.set_xticklabels(labs, rotation=90)
    yt = np.linspace(0, len(lp) - 1, num=min(8, len(lp)))
    yv = np.linspace(lp.min(), lp.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("Δ log10 |Z_off| (after−before)")
    return ax


# 2) SNR gain profile (per-station dB improvement) ----------------------- #

def nr_qc_snr_gain_profile(
    sites: Any,
    *,
    method: str = "pipeline",
    pband: Optional[Tuple[float, float]] = None,
    figsize: Tuple[float, float] = (8.6, 3.6),
    ax: Optional[plt.Axes] = None,
    **denoise: Any,
) -> plt.Axes:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _denoise_sites(S0, method, inplace=False, **denoise)
    T0 = snr_table(S0)
    T1 = snr_table(S1)
    if T0.empty or T1.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no snr", ha="center", va="center")
        return ax
    def _band(v, f):
        if pband is None:
            return np.nanmedian(v)
        lo, hi = pband
        m = (1.0 / f >= lo) & (1.0 / f <= hi)
        return np.nanmedian(v[m]) if np.any(m) else np.nan
    labs, gain = [], []
    for st in sorted(set(T0["station"]) & set(T1["station"])):
        s0 = T0[T0["station"] == st]
        s1 = T1[T1["station"] == st]
        if s0.empty or s1.empty:
            continue
        f0 = s0["freq"].to_numpy(dtype=float)
        f1 = s1["freq"].to_numpy(dtype=float)
        sn0 = s0["snr"].to_numpy(dtype=float)
        sn1 = s1["snr"].to_numpy(dtype=float)
        g0 = _band(sn0, f0)
        g1 = _band(sn1, f1)
        if np.isfinite(g0) and np.isfinite(g1) and g0 > 0:
            labs.append(st)
            gain.append(20.0 * np.log10(g1 / g0))
    if not labs:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no gain", ha="center", va="center")
        return ax
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(labs))
    ax.axhline(0.0, color="0.7", lw=1.0)
    ax.bar(x, gain, width=0.8)
    ax.set_ylabel("SNR gain (dB)")
    ax.set_xlabel("Station")
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=90)
    return ax


# 3) Harmonic waterfall (Δ dB at mains harmonics) ------------------------ #

def nr_qc_harmonic_waterfall(
    sites: Any,
    *,
    method: str = "notch",
    mains_hz: float = 50.0,
    n_harm: int = 30,
    tol_hz: float = 0.08,
    figsize: Tuple[float, float] = (9.0, 4.6),
    ax: Optional[plt.Axes] = None,
    **denoise: Any,
) -> plt.Axes:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _denoise_sites(
        S0, method, inplace=False,
        mains_hz=mains_hz, n_harm=n_harm, tol_hz=tol_hz,
        **denoise,
    )
    sts, H = [], []
    for i, (S, tag) in enumerate([(S0, "b"), (S1, "a")]):
        rows = []
        labs = []
        for j, ed in enumerate(_iter_items(S)):
            Z, z, fr = _get_z_block(ed)
            if Z is None:
                continue
            lm = _offdiag_logmag(z)
            kk = np.arange(1, int(n_harm) + 1, dtype=float)
            hv = []
            for k in kk:
                m = np.abs(fr - k * mains_hz) <= float(tol_hz)
                hv.append(
                    float(np.nanmedian(lm[m])) if np.any(m) else np.nan
                )
            rows.append(hv)
            if i == 0:
                labs.append(_name(ed, j))
        if i == 0:
            sts = labs
            Hb = np.array(rows, float)
        else:
            Ha = np.array(rows, float)
    if not sts:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    D = (Hb - Ha).T  # reduction (dB-like log units)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        D,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="viridis",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("Harmonic index k (k·mains)")
    ax.set_xticks(np.arange(len(sts)))
    ax.set_xticklabels(sts, rotation=90)
    yt = np.linspace(0, D.shape[0] - 1, num=min(8, D.shape[0]))
    yv = np.linspace(1, D.shape[0], num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{int(v)}" for v in yv])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("Δ log10 |Z_off| (reduction)")
    return ax


# 4) Station curves (off-diag mag with mains guides) --------------------- #

def nr_qc_station_offdiag_curves(
    sites: Any,
    *,
    method: str = "pipeline",
    station: Optional[str] = None,
    mains_hz: float = 50.0,
    n_harm: int = 12,
    tol_hz: float = 0.08,
    figsize: Tuple[float, float] = (8.0, 4.2),
    ax: Optional[plt.Axes] = None,
    **denoise: Any,
) -> plt.Axes:
    S0 = ensure_sites(sites, recursive=False, strict=False)
    S1 = _denoise_sites(S0, method, inplace=False, **denoise)
    # pick station
    fm = {}
    for i, ed in enumerate(_iter_items(S0)):
        fm[_name(ed, i)] = ed
    if station is None and fm:
        station = sorted(fm.keys())[0]
    ed0 = fm.get(station, None)
    if ed0 is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "station not found",
                ha="center", va="center")
        return ax
    # match from S1
    for ed in _iter_items(S1):
        if _name(ed, 0) == station:
            ed1 = ed; break
    else:
        ed1 = None
    if ed1 is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no after", ha="center", va="center")
        return ax
    Z0, z0, f0 = _get_z_block(ed0)
    Z1, z1, f1 = _get_z_block(ed1)
    if Z0 is None or Z1 is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no Z", ha="center", va="center")
        return ax
    m0 = 10 ** _offdiag_logmag(z0)
    m1 = 10 ** _offdiag_logmag(z1)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.set_xscale("log")
    ax.plot(1.0 / f0, m0, "o-", lw=1.4, label="before")
    ax.plot(1.0 / f1, m1, "s-", lw=1.4, label="after")
    # mains guides
    kk = np.arange(1, int(n_harm) + 1, dtype=float)
    for k in kk:
        f = k * mains_hz
        ax.axvspan(
            1.0 / (f + tol_hz),
            1.0 / max(1e-12, (f - tol_hz)),
            alpha=0.08, color="r",
        )
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("|Z_off|")
    ax.set_title(str(station))
    ax.grid(True, alpha=0.25, which="both")
    ax.legend()
    return ax
