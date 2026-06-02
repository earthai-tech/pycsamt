
from __future__ import annotations 

from typing import List, Any, Optional, Sequence, Tuple 
import numpy as np

from ._core import (
    ensure_sites,
    _apply_each,
    _iter_items,
    _get_z_block,
    _get_t_block,
)

# re-use package editors when it saves code
from pycsamt.site import edit as _edit

_BACKWARD_SINCE = "2.0.0"
_BACKWARD_REMOVE = "2.17.0"


# ------------------------------- helpers -------------------------------- #

def _sorted_unique(fr: np.ndarray) -> np.ndarray:
    fr = np.asarray(fr, dtype=float)
    fr = fr[np.isfinite(fr)]
    if fr.size == 0:
        return fr
    fr = np.unique(fr)
    fr = fr[fr > 0.0]
    return np.sort(fr)


def _nearest_idx(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    # for each y, nearest index in x
    idx = np.searchsorted(x, y)
    idx = np.clip(idx, 1, x.size - 1)
    left = x[idx - 1]
    right = x[idx]
    pick_left = (y - left) <= (right - y)
    idx[pick_left] -= 1
    return idx


def _interp_complex(
    x: np.ndarray,
    y: np.ndarray,
    xnew: np.ndarray,
    *,
    method: str = "nearest",
) -> np.ndarray:
    if y.ndim == 1:
        y = y[:, None]
    r = y.real
    im = y.imag
    if method == "nearest":
        idx = _nearest_idx(x, xnew)
        rr = r[idx]
        ii = im[idx]
    else:
        rr = np.vstack(
            [np.interp(xnew, x, r[:, j]) for j in range(r.shape[1])]
        ).T
        ii = np.vstack(
            [np.interp(xnew, x, im[:, j]) for j in range(im.shape[1])]
        ).T
    out = rr + 1j * ii
    return out.squeeze()




def _set_block_freq(obj: Any, fr: np.ndarray) -> None:
    try:
        obj.freq = fr
    except Exception:
        pass


def _regrid_z(
    z: np.ndarray,
    fr: np.ndarray,
    fnew: np.ndarray,
    *,
    method: str,
) -> np.ndarray:
    # z: (n,2,2) complex
    n = fnew.size
    out = np.empty((n, 2, 2), dtype=complex)
    x = np.log10(fr)
    xn = np.log10(fnew)
    for a in range(2):
        for b in range(2):
            y = z[:, a, b]
            out[:, a, b] = _interp_complex(
                x, y, xn, method=method
            )
    return out


def _regrid_t(
    t: np.ndarray,
    fr: np.ndarray,
    fnew: np.ndarray,
    *,
    method: str,
) -> np.ndarray:
    # t: (n,2) complex
    n = fnew.size
    out = np.empty((n, 2), dtype=complex)
    x = np.log10(fr)
    xn = np.log10(fnew)
    for k in range(2):
        y = t[:, k]
        out[:, k] = _interp_complex(x, y, xn, method=method)
    return out


def _apply_grid_one(
    ed: Any,
    fnew: np.ndarray,
    *,
    method: str,
) -> None:
    Z, z, frz = _get_z_block(ed)
    if Z is not None:
        z2 = _regrid_z(z, frz, fnew, method=method)
        try:
            Z.z = z2
        except Exception:
            pass
        _set_block_freq(Z, fnew)
    T, t, frt = _get_t_block(ed)
    if T is not None:
        t2 = _regrid_t(t, frt, fnew, method=method)
        try:
            T.tipper = t2
        except Exception:
            pass
        _set_block_freq(T, fnew)


def _valid_freq_of(ed: Any) -> Optional[np.ndarray]:
    Z, z, fr = _get_z_block(ed)
    if isinstance(fr, np.ndarray) and fr.size:
        return _sorted_unique(fr)
    T, t, fr = _get_t_block(ed)
    if isinstance(fr, np.ndarray) and fr.size:
        return _sorted_unique(fr)
    return None


def _union_grid(SitesObj) -> np.ndarray:
    frs: List[np.ndarray] = []
    for ed in _iter_items(SitesObj):
        fr = _valid_freq_of(ed)
        if isinstance(fr, np.ndarray) and fr.size:
            frs.append(fr)
    if not frs:
        return np.array([], dtype=float)
    g = np.unique(np.concatenate(frs))
    return _sorted_unique(g)


def _intersect_grid(SitesObj) -> np.ndarray:
    frs: List[np.ndarray] = []
    for ed in _iter_items(SitesObj):
        fr = _valid_freq_of(ed)
        if isinstance(fr, np.ndarray) and fr.size:
            frs.append(fr)
    if not frs:
        return np.array([], dtype=float)
    g = frs[0]
    for fr in frs[1:]:
        g = np.intersect1d(g, fr)
    return _sorted_unique(g)


# ------------------------------ public API ------------------------------- #

def select_band(
    sites: Any,
    *,
    fmin: Optional[float] = None,
    fmax: Optional[float] = None,
    keep: Optional[Sequence[float]] = None,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # delegate to existing editor for robustness
    return _edit.select_freq(
        S,
        fmin=fmin,
        fmax=fmax,
        keep=keep,
        inplace=inplace,
    )


def drop_duplicates(
    sites: Any,
    *,
    tol: float = 1e-10,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        for objname in ("Z", "z", "Tipper", "tipper"):
            O = getattr(ed, objname, None)
            if O is None:
                continue
            fr = getattr(O, "freq", None)
            if not isinstance(fr, np.ndarray):
                continue
            fr = np.asarray(fr, dtype=float)
            if fr.size == 0:
                continue
            frs = _sorted_unique(fr)
            if frs.size == fr.size and np.allclose(fr, frs):
                continue
            # nearest index map
            idx = _nearest_idx(frs, fr)
            # keep first occurrence only
            seen = set()
            keep_idx = []
            for i, j in enumerate(idx):
                if abs(fr[i] - frs[j]) <= tol and j not in seen:
                    keep_idx.append(i)
                    seen.add(j)
            keep_idx = np.array(keep_idx, dtype=int)
            for fld in ("z", "z_err", "tipper", "tipper_err"):
                val = getattr(O, fld, None)
                if isinstance(val, np.ndarray) and val.shape[0] == fr.size:
                    setattr(O, fld, val[keep_idx])
            _set_block_freq(O, fr[keep_idx])
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def regrid_to(
    sites: Any,
    target_freq: Sequence[float],
    *,
    method: str = "nearest",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    fnew = _sorted_unique(np.asarray(target_freq, dtype=float))

    def _one(Si):
        ed = next(_iter_items(Si))
        if fnew.size == 0:
            return Si
        _apply_grid_one(ed, fnew, method=method)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def regrid_logspace(
    sites: Any,
    *,
    fmin: Optional[float] = None,
    fmax: Optional[float] = None,
    per_decade: int = 6,
    method: str = "nearest",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # infer band from union
    U = _union_grid(S)
    if U.size == 0:
        return S
    lo = fmin if fmin is not None else U.min()
    hi = fmax if fmax is not None else U.max()
    lo, hi = float(lo), float(hi)
    if lo <= 0 or hi <= 0 or hi <= lo:
        return S
    ndec = np.log10(hi) - np.log10(lo)
    npts = max(2, int(np.ceil(ndec * per_decade)))
    fnew = np.logspace(np.log10(lo), np.log10(hi), npts)
    return regrid_to(
        S,
        fnew,
        method=method,
        inplace=inplace,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )


def decimate_step(
    sites: Any,
    *,
    step: int = 2,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, frz = _get_z_block(ed)
        if Z is not None:
            idx = np.arange(0, frz.size, step, dtype=int)
            Z.z = z[idx]
            _set_block_freq(Z, frz[idx])
        T, t, frt = _get_t_block(ed)
        if T is not None:
            idx = np.arange(0, frt.size, step, dtype=int)
            T.tipper = t[idx]
            _set_block_freq(T, frt[idx])
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def smooth_mavg(
    sites: Any,
    *,
    k: int = 3,
    on: str = "both",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    k = max(1, int(k))
    if k == 1:
        return S

    def _roll(y: np.ndarray) -> np.ndarray:
        if y.ndim == 1:
            y = y[:, None]
        out = np.empty_like(y)
        for j in range(y.shape[1]):
            v = y[:, j].astype(complex)
            re = np.convolve(v.real, np.ones(k) / k, "same")
            im = np.convolve(v.imag, np.ones(k) / k, "same")
            out[:, j] = re + 1j * im
        return out.squeeze()

    def _one(Si):
        ed = next(_iter_items(Si))
        if on in ("z", "both"):
            Z, z, fr = _get_z_block(ed)
            if Z is not None:
                z2 = z.copy()
                for a in range(2):
                    for b in range(2):
                        z2[:, a, b] = _roll(z[:, a, b])
                Z.z = z2
        if on in ("tipper", "both"):
            T, t, fr = _get_t_block(ed)
            if T is not None:
                t2 = t.copy()
                for j in range(2):
                    t2[:, j] = _roll(t[:, j])
                T.tipper = t2
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def align_grid(
    sites: Any,
    *,
    mode: str = "union",
    ref_station: Optional[str] = None,
    method: str = "nearest",
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if mode == "union":
        fnew = _union_grid(S)
    elif mode == "intersection":
        fnew = _intersect_grid(S)
    elif mode == "ref":
        fnew = None
        for ed in _iter_items(S):
            nm = getattr(ed, "station", None) or getattr(ed, "name", None)
            if nm == ref_station:
                fr = _valid_freq_of(ed)
                if isinstance(fr, np.ndarray):
                    fnew = fr
                break
        if fnew is None:
            return S
    else:
        raise ValueError("mode must be union|intersection|ref")
    if fnew.size == 0:
        return S
    return regrid_to(
        S,
        fnew,
        method=method,
        inplace=inplace,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )


# ------------------------- coverage + quality ---------------------------- #

import matplotlib.pyplot as plt


def _rel_err_block(z: np.ndarray, ze: Optional[np.ndarray]) -> np.ndarray:
    if ze is None:
        return np.zeros(z.shape[0], dtype=float)
    # avg relative err across off-diagonals (fallback: all)
    zx = z[:, 0, 1]
    zy = z[:, 1, 0]
    ex = ze[:, 0, 1] if ze.ndim == 3 else None
    ey = ze[:, 1, 0] if ze.ndim == 3 else None
    if ex is not None and ey is not None:
        rx = np.abs(ex) / (np.abs(zx) + 1e-12)
        ry = np.abs(ey) / (np.abs(zy) + 1e-12)
        r = 0.5 * (rx + ry)
        return np.asarray(r, dtype=float)
    # fallback to full-tensor mean
    num = np.linalg.norm(ze.reshape(ze.shape[0], -1), axis=1)
    den = np.linalg.norm(z.reshape(z.shape[0], -1), axis=1) + 1e-12
    return num / den


def plot_coverage_quality_heatmap(
    sites: Any,
    *,
    axis: str = "period",
    figsize: Tuple[float, float] = (7.5, 4.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # assemble presence & quality per site on union grid
    labs: List[str] = []
    frs: List[np.ndarray] = []
    quals: List[np.ndarray] = []
    for i, ed in enumerate(_iter_items(S)):
        nm = getattr(ed, "station", None) or getattr(ed, "name", None)
        nm = nm or f"site_{i}"
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        ze = getattr(Z, "z_err", None)
        r = _rel_err_block(z, ze)
        q = 1.0 / (1.0 + r)
        frs.append(_sorted_unique(fr))
        # re-align q to unique set by averaging duplicates
        uf = frs[-1]
        if uf.size != fr.size:
            idx = _nearest_idx(uf, fr)
            acc = np.zeros(uf.size, dtype=float)
            cnt = np.zeros(uf.size, dtype=float)
            for k, j in enumerate(idx):
                acc[j] += q[k]
                cnt[j] += 1.0
            q = acc / (cnt + 1e-12)
        else:
            q = q
        quals.append(q)
        labs.append(nm)
    grid = _union_grid(S)
    if grid.size == 0 or not labs:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    # project site vectors onto union grid
    M = np.zeros((len(labs), grid.size), dtype=float)
    for i, (fr, q) in enumerate(zip(frs, quals)):
        idx = _nearest_idx(grid, fr)
        M[i, idx] = q
    # 0..1 quality; 0 means absent/low confidence
    if axis == "period":
        y = 1.0 / grid
        ylab = "period (s)"
    else:
        y = grid
        ylab = "freq (Hz)"
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        M.T,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        vmin=0.0,
        vmax=1.0,
    )
    ax.set_xlabel("station")
    ax.set_ylabel(ylab)
    ax.set_xticks(np.arange(len(labs)))
    ax.set_xticklabels(labs, rotation=90)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("quality (1/(1+rel_err))")
    return ax


# ------------------------- apparent depth image -------------------------- #

def _rho_det_from_z(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    # ρa_xy/yx = 0.2 * |Z|^2 / f (ohm·m); det via geom. mean
    zx = z[:, 0, 1]
    zy = z[:, 1, 0]
    rx = 0.2 * (np.abs(zx) ** 2) / (fr + 1e-24)
    ry = 0.2 * (np.abs(zy) ** 2) / (fr + 1e-24)
    rdet = np.sqrt(rx * ry)
    return rdet


def _depth_from_rho_f(rho: np.ndarray, fr: np.ndarray) -> np.ndarray:
    # δ ≈ 503 * sqrt(ρ / f)  (meters)
    with np.errstate(invalid="ignore"):
        d = 503.0 * np.sqrt(np.maximum(rho, 0.0) / (fr + 1e-24))
    return d


def plot_apparent_depth_psection(
    sites: Any,
    *,
    axis_y: str = "period",
    agg: str = "median",
    figsize: Tuple[float, float] = (7.5, 4.5),
    log_color: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # gather rows: station, freq, period, depth
    rows: List[Tuple[str, float, float, float]] = []
    for i, ed in enumerate(_iter_items(S)):
        nm = getattr(ed, "station", None) or getattr(ed, "name", None)
        nm = nm or f"site_{i}"
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        rho = _rho_det_from_z(z, fr)
        dep = _depth_from_rho_f(rho, fr)
        per = 1.0 / fr
        for f, p, d in zip(fr, per, dep):
            if not np.isfinite(d):
                continue
            rows.append((nm, f, p, d))
    if not rows:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    # pivot to station × axis_y
    dtype = [("station", "U64"), ("freq", "f8"),
             ("period", "f8"), ("depth", "f8")]
    arr = np.array(rows, dtype=dtype)
    import pandas as pd
    df = pd.DataFrame(arr)
    ykey = "period" if axis_y == "period" else "freq"
    piv = df.pivot_table(
        index=ykey,
        columns="station",
        values="depth",
        aggfunc=agg,
    )
    piv = piv.sort_index()
    Zm = piv.to_numpy(dtype=float).T
    if log_color:
        Zplot = np.log10(np.maximum(Zm, 1e-3))
        cblab = "log10 depth (m)"
    else:
        Zplot = Zm
        cblab = "depth (m)"
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        Zplot,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("station")
    ax.set_ylabel(ykey)
    ax.set_xticks(np.arange(Zplot.shape[0]))
    ax.set_xticklabels(piv.columns, rotation=90)
    yt = np.linspace(0, Zplot.shape[1] - 1,
                     num=min(8, Zplot.shape[1]))
    yv = np.linspace(piv.index.min(), piv.index.max(),
                     num=min(8, len(piv.index)))
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.3g}" for v in yv])
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(cblab)
    return ax


# ------------------------- band-collapsed microplots --------------------- #

def _det_phase_from_z(z: np.ndarray) -> np.ndarray:
    detz = z[:, 0, 0] * z[:, 1, 1] - z[:, 0, 1] * z[:, 1, 0]
    return np.degrees(np.angle(detz))


def _tip_amp(t: Optional[np.ndarray]) -> Optional[np.ndarray]:
    if t is None:
        return None
    return np.sqrt(np.abs(t[:, 0]) ** 2 + np.abs(t[:, 1]) ** 2)


def _bands_auto(all_period: np.ndarray, n: int) -> np.ndarray:
    lo = float(np.nanmin(all_period))
    hi = float(np.nanmax(all_period))
    lo = max(lo, 1e-6)
    edges = np.logspace(np.log10(lo), np.log10(hi), n + 1)
    return edges


def plot_band_microstrips(
    sites: Any,
    *,
    bands: Optional[Sequence[Tuple[float, float]]] = None,
    n_bands: int = 6,
    figsize: Tuple[float, float] = (9.0, 6.0),
    marker_size: float = 16.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # collect per-site arrays
    sites_data = []
    allP = []
    for i, ed in enumerate(_iter_items(S)):
        nm = getattr(ed, "station", None) or getattr(ed, "name", None)
        nm = nm or f"site_{i}"
        Z, z, fr = _get_z_block(ed)
        T, t, frt = _get_t_block(ed)
        if Z is None:
            continue
        per = 1.0 / fr
        rho = _rho_det_from_z(z, fr)
        lgr = np.log10(np.maximum(rho, 1e-12))
        ph = _det_phase_from_z(z)
        ta = _tip_amp(t) if T is not None else None
        sites_data.append((nm, per, lgr, ph, ta))
        allP.append(per)
    if not sites_data:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    allP = np.concatenate(allP)
    if bands is None:
        edges = _bands_auto(allP, n_bands)
        bands = list(zip(edges[:-1], edges[1:]))
    # aggregate per band
    rows = []
    for (nm, per, lgr, ph, ta) in sites_data:
        for bi, (lo, hi) in enumerate(bands):
            m = (per >= lo) & (per < hi)
            if not np.any(m):
                rows.append((nm, bi, np.nan, np.nan, np.nan))
                continue
            r = float(np.nanmedian(lgr[m]))
            p = float(np.nanmedian(ph[m]))
            if ta is not None and ta.size == per.size:
                tt = float(np.nanmedian(ta[m]))
            else:
                tt = np.nan
            rows.append((nm, bi, r, p, tt))
    # normalize per metric for color
    import pandas as pd
    df = pd.DataFrame(
        rows, columns=["station", "band", "logrho", "phi", "tamp"]
    )
    # global min/max (robust) for color normalization
    def _rng(s: pd.Series) -> Tuple[float, float]:
        v = s.to_numpy()
        v = v[np.isfinite(v)]
        if v.size == 0:
            return 0.0, 1.0
        return (np.nanpercentile(v, 5), np.nanpercentile(v, 95))
    r0, r1 = _rng(df["logrho"])
    p0, p1 = _rng(df["phi"])
    t0, t1 = _rng(df["tamp"])
    def _norm(x, a, b):
        return np.clip((x - a) / (b - a + 1e-12), 0.0, 1.0)
    df["cr"] = _norm(df["logrho"], r0, r1)
    df["cp"] = _norm(df["phi"], p0, p1)
    df["ct"] = _norm(df["tamp"], t0, t1)
    # plot: 3 dots per band (circle, square, triangle)
    st = list(df["station"].unique())
    st_map = {s: i for i, s in enumerate(st)}
    y0 = np.arange(len(st), dtype=float)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    for _, row in df.iterrows():
        yi = st_map[row["station"]]
        xi = float(row["band"])
        # offsets for 3 metrics
        yR = yi - 0.22
        yP = yi + 0.00
        yT = yi + 0.22
        ax.scatter(
            [xi], [yR],
            s=marker_size,
            c=[[row["cr"]]],
            cmap="viridis",
            marker="o",
            vmin=0.0, vmax=1.0,
        )
        ax.scatter(
            [xi], [yP],
            s=marker_size,
            c=[[row["cp"]]],
            cmap="viridis",
            marker="s",
            vmin=0.0, vmax=1.0,
        )
        ax.scatter(
            [xi], [yT],
            s=marker_size,
            c=[[row["ct"]]],
            cmap="viridis",
            marker="^",
            vmin=0.0, vmax=1.0,
        )
    ax.set_xlim(-0.5, len(bands) - 0.5)
    ax.set_ylim(-0.6, len(st) - 0.4)
    ax.set_xlabel("band")
    ax.set_ylabel("station")
    ax.set_yticks(y0)
    ax.set_yticklabels(st)
    # simple mini-legend using text markers
    ax.text(1.01, 0.85, "o log10 ρ_det",
            transform=ax.transAxes)
    ax.text(1.01, 0.70, "s φ_det",
            transform=ax.transAxes)
    ax.text(1.01, 0.55, "^ |T|",
            transform=ax.transAxes)
    return ax
