# pycsamt/emtools/strike.py
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple 

import re
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection


import numpy as np
import pandas as pd

from ._core import (
    ensure_sites,
    _iter_items,
    _apply_each,
    _get_z_block,
    _name,
)
from ..site import edit as _edit
from .tensor import build_phase_tensor_table

# -------------------------- small helpers ------------------------------- #

def _rotmat(deg: float) -> np.ndarray:
    th = np.radians(float(deg))
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, s], [-s, c]], dtype=float)


def _rotate_tensor(z: np.ndarray, deg: float) -> np.ndarray:
    R = _rotmat(deg)
    Rt = R.T
    out = np.empty_like(z)
    for i in range(z.shape[0]):
        out[i] = R @ z[i] @ Rt
    return out


def _diag_offdiag_norm(z: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    d = np.sqrt(np.abs(z[:, 0, 0]) ** 2 + np.abs(z[:, 1, 1]) ** 2)
    o = np.sqrt(np.abs(z[:, 0, 1]) ** 2 + np.abs(z[:, 1, 0]) ** 2)
    return d, o


def _score(z: np.ndarray, kind: str) -> np.ndarray:
    d, o = _diag_offdiag_norm(z)
    if kind == "diag_ratio":
        return d / (o + 1e-12)
    if kind == "offdiag_neg":
        return -o
    if kind == "det_diag":
        det = np.abs(z[:, 0, 0] * z[:, 1, 1])
        return det / (o + 1e-12)
    return d / (o + 1e-12)


def _wrap90(a: np.ndarray) -> np.ndarray:
    x = ((a + 90.0) % 180.0) - 90.0
    return x


def _unwrap_deg_180(a: np.ndarray) -> np.ndarray:
    x = a.copy().astype(float)
    for i in range(1, x.size):
        d = x[i] - x[i - 1]
        if d > 90.0:
            x[i:] -= 180.0
        elif d < -90.0:
            x[i:] += 180.0
    return x


def _band_edges(p: np.ndarray, band: Optional[Tuple[float, float]]):
    if band is None:
        lo = max(1e-6, float(np.nanmin(p)))
        hi = float(np.nanmax(p))
        return lo, hi
    return float(band[0]), float(band[1])


# ----------------------- strike by sweep (Z) ----------------------------- #

def estimate_strike_sweep(
    sites: Any,
    *,
    angles: np.ndarray = np.arange(-90.0, 91.0, 1.0),
    metric: str = "diag_ratio",
    band: Optional[Tuple[float, float]] = None,
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
    rows: List[Dict[str, float]] = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        p = 1.0 / fr
        lo, hi = _band_edges(p, band)
        m = (p >= lo) & (p <= hi)
        if not np.any(m):
            continue
        zB = z[m]
        best = np.zeros(zB.shape[0], dtype=float)
        val = np.full(zB.shape[0], np.inf, dtype=float)
        for a in angles:
            zr = _rotate_tensor(zB, a)
            sc = _score(zr, metric)
            upd = sc < val
            val[upd] = sc[upd]
            best[upd] = a
        best_u = _unwrap_deg_180(best)
        ang = _wrap90(np.nanmedian(best_u))
        iqr = np.nanpercentile(best_u, 75) - np.nanpercentile(
            best_u, 25
        )
        rows.append(
            dict(
                station=st, ang=ang, iqr=iqr,
                lo=lo, hi=hi, n=int(zB.shape[0]),
            )
        )
    cols = ["station", "ang", "iqr", "lo", "hi", "n"]
    return pd.DataFrame.from_records(rows, columns=cols)


# ------------------ strike from phase tensor (theta) --------------------- #

def estimate_strike_phase_tensor(
    sites: Any,
    *,
    band: Optional[Tuple[float, float]] = None,
    robust: bool = True,
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
    df = build_phase_tensor_table(
        S, recursive=False, on_dup=on_dup,
        strict=False, verbose=verbose,
    )
    if df.empty:
        return pd.DataFrame(
            columns=["station", "ang", "iqr", "lo", "hi", "n"]
        )
    rows: List[Dict[str, float]] = []
    for st, sdf in df.groupby("station"):
        p = sdf["period"].to_numpy()
        lo, hi = _band_edges(p, band)
        m = (p >= lo) & (p <= hi)
        if not np.any(m):
            continue
        th = sdf.loc[m, "theta"].to_numpy(dtype=float)
        th = _unwrap_deg_180(th)
        if robust:
            ang = _wrap90(np.nanmedian(th))
            iqr = np.nanpercentile(th, 75) - np.nanpercentile(
                th, 25
            )
        else:
            ang = _wrap90(float(np.nanmean(th)))
            iqr = float(np.nanstd(th))
        rows.append(
            dict(station=st, ang=ang, iqr=iqr, lo=lo, hi=hi,
                 n=int(np.sum(m)))
        )
    return pd.DataFrame.from_records(
        rows, columns=["station", "ang", "iqr", "lo", "hi", "n"]
    )


# ---------------------- consensus strike (blend) ------------------------- #

def estimate_strike_consensus(
    sites: Any,
    *,
    band: Optional[Tuple[float, float]] = None,
    w_sweep: float = 0.5,
    w_pt: float = 0.5,
    angles: np.ndarray = np.arange(-90.0, 91.0, 1.0),
    metric: str = "diag_ratio",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    t1 = estimate_strike_sweep(
        sites,
        angles=angles, metric=metric,
        band=band, recursive=recursive,
        on_dup=on_dup, strict=strict, verbose=verbose,
    )
    t2 = estimate_strike_phase_tensor(
        sites,
        band=band, robust=True,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if t1.empty and t2.empty:
        return pd.DataFrame(
            columns=["station", "ang", "iqr", "lo", "hi", "n"]
        )
    df = pd.merge(
        t1, t2, on="station", how="outer",
        suffixes=("_sweep", "_pt"),
    )
    def pick(a, b, ws, wp):
        if np.isnan(a) and np.isnan(b):
            return np.nan
        a = 0.0 if np.isnan(a) else a
        b = 0.0 if np.isnan(b) else b
        u = _unwrap_deg_180(np.array([a, b]))
        u = _wrap90(u)
        return float(ws * u[0] + wp * u[1])
    ang = []
    iqr = []
    lo = []
    hi = []
    n = []
    for _, r in df.iterrows():
        ang.append(
            pick(r.get("ang_sweep", np.nan),
                 r.get("ang_pt", np.nan),
                 w_sweep, w_pt)
        )
        i1 = r.get("iqr_sweep", np.nan)
        i2 = r.get("iqr_pt", np.nan)
        iqr.append(
            float(np.nanmedian([i1, i2]))
        )
        lo.append(
            float(np.nanmin([
                r.get("lo_sweep", np.nan), r.get("lo_pt", np.nan)
            ]))
        )
        hi.append(
            float(np.nanmax([
                r.get("hi_sweep", np.nan), r.get("hi_pt", np.nan)
            ]))
        )
        n.append(int(np.nansum([r.get("n_sweep", 0), r.get("n_pt", 0)])))
    out = pd.DataFrame(
        dict(station=df["station"], ang=ang, iqr=iqr,
             lo=lo, hi=hi, n=n)
    )
    return out


# ----------------------- rotation applicators ---------------------------- #

def rotate_to_strike(
    sites: Any,
    *,
    method: str = "consensus",  # consensus|sweep|pt
    band: Optional[Tuple[float, float]] = None,
    angles: np.ndarray = np.arange(-90.0, 91.0, 1.0),
    metric: str = "diag_ratio",
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
    if method == "sweep":
        TB = estimate_strike_sweep(
            S, angles=angles, metric=metric,
            band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    elif method == "pt":
        TB = estimate_strike_phase_tensor(
            S, band=band, robust=True,
            recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    else:
        TB = estimate_strike_consensus(
            S, band=band, angles=angles, metric=metric,
            recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    amap = {r.station: float(r.ang) for _, r in TB.iterrows()}

    def _one(Si):
        ed = next(_iter_items(Si))
        st = getattr(ed, "station", None) or getattr(ed, "name", None)
        ang = float(amap.get(st, 0.0))
        return _edit.rotate(Si, angle=ang, inplace=inplace)

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --------------------- per-frequency strike curve ------------------------ #

def strike_curve_sweep(
    sites: Any,
    *,
    angles: np.ndarray = np.arange(-90.0, 91.0, 1.0),
    metric: str = "diag_ratio",
    smooth: int = 5,
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
    rows: List[Dict[str, float]] = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        best = np.zeros(z.shape[0], dtype=float)
        val = np.full(z.shape[0], np.inf, dtype=float)
        for a in angles:
            zr = _rotate_tensor(z, a)
            sc = _score(zr, metric)
            upd = sc < val
            val[upd] = sc[upd]
            best[upd] = a
        best = _unwrap_deg_180(best)
        if smooth > 1 and best.size >= smooth:
            k = int(smooth)
            w = np.ones(k) / k
            u = np.convolve(best, w, mode="same")
            best = u
        for f, ang in zip(fr, _wrap90(best)):
            rows.append(dict(station=st, freq=float(f), ang=float(ang)))
    return pd.DataFrame.from_records(rows, columns=["station", "freq", "ang"])


def _auto_line(st: str) -> str:
    m = re.match(r"^([A-Za-z]+[0-9]+)", str(st))
    return m.group(1) if m else str(st)


def _axial_mean_deg(a: np.ndarray, w: np.ndarray) -> float:
    # axial mean: double angles, mean vector, halve back
    th = np.radians(2.0 * (a % 180.0))
    cw = np.cos(th) * w
    sw = np.sin(th) * w
    ang = 0.5 * np.degrees(np.arctan2(sw.sum(), cw.sum()))
    ang = (ang + 180.0) % 180.0
    return float(ang)

def plot_strike_rose_by_line(
    sites: Any,
    *,
    groups: dict[str, list[str]] | None = None,
    group_key: str | None = None,
    band: tuple[float, float] | None = None,
    method: str = "consensus",  # consensus|sweep|pt
    bins: int = 36,
    weight: str = "inv_iqr",  # inv_iqr|uniform
    figsize: tuple[float, float] = (8.6, 4.6),
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
    # 1) per-station strike + iqr (weight proxy)
    if method == "sweep":
        TB = estimate_strike_sweep(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    elif method == "pt":
        TB = estimate_strike_phase_tensor(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    else:
        TB = estimate_strike_consensus(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    if TB.empty:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, polar=True)
        ax.text(0.5, 0.5, "no strikes", ha="center", va="center")
        return fig
    TB = TB.copy()
    TB["ang"] = (TB["ang"] % 180.0).astype(float)
    TB["w"] = 1.0
    if weight == "inv_iqr":
        TB["w"] = 1.0 / (TB["iqr"].abs().to_numpy() + 1e-6)
    # 2) build group membership
    if groups is None:
        # from attribute on EDI, else from station prefix like E1
        lab = {}
        for i, ed in enumerate(_iter_items(S)):
            st = _name(ed, i)
            if group_key and hasattr(ed, group_key):
                lab[st] = str(getattr(ed, group_key))
            else:
                lab[st] = _auto_line(st)
        groups = {}
        for st, g in lab.items():
            groups.setdefault(g, []).append(st)
    # keep groups with at least 2 stations
    groups = {g: v for g, v in groups.items() if len(v) >= 2}
    if not groups:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, polar=True)
        ax.text(0.5, 0.5, "no groups", ha="center", va="center")
        return fig
    # 3) figure grid
    G = list(groups.keys())
    n = len(G)
    fig = plt.figure(figsize=figsize)
    for i, g in enumerate(G, 1):
        ax = fig.add_subplot(1, n, i, polar=True)
        subset = TB[TB["station"].isin(groups[g])]
        if subset.empty:
            ax.text(0.5, 0.5, "empty", ha="center", va="center")
            continue
        ang = subset["ang"].to_numpy(dtype=float)
        w = subset["w"].to_numpy(dtype=float)
        # histogram on 0..180 (axial), then mirror to 0..360
        bins = int(max(12, bins))
        edges = np.linspace(0.0, 180.0, bins + 1)
        h, _ = np.histogram(ang, bins=edges, weights=w)
        cen = 0.5 * (edges[1:] + edges[:-1])
        th = np.radians(np.concatenate([cen, cen + 180.0]))
        rr = np.concatenate([h, h])
        # color by height (nice gradient)
        vmax = rr.max() if rr.size else 1.0
        cols = plt.cm.YlOrRd(rr / (vmax + 1e-12))
        ax.bar(
            th, rr,
            width=np.radians(edges[1] - edges[0]),
            bottom=0.0,
            color=cols,
            edgecolor="none",
        )
        # axial mean + label box
        mu = _axial_mean_deg(ang, w)
        rmax = float(vmax) * 1.05
        for add in (0.0, 180.0):
            ax.plot(
                [np.radians(mu + add), np.radians(mu + add)],
                [0.0, rmax],
                color="crimson",
                lw=2.5,
            )
        ax.text(
            0.08, 0.90, f"{mu:.1f}°",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.2",
                      fc="white", ec="0.2", lw=0.6),
        )
        # polar cosmetics
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_title(str(g), pad=12.0)
    fig.tight_layout()
    return fig

# ---- STRIKE VIEWS: ribbon, profile, and map-sticks --------------------- #

def _hsv_rgb(h, s, v):
    hsv = np.stack([h, s, v], axis=-1)
    return mcolors.hsv_to_rgb(hsv)

def plot_strike_ribbon(
    sites: Any,
    *,
    method: str = "sweep",  # sweep|pt|consensus
    win: int = 5,
    figsize: tuple[float, float] = (9.0, 4.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    if method == "sweep":
        df = strike_curve_sweep(
            sites,
            recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )
    else:
        # consensus/pt → per-station single angle; expand flat
        tb = estimate_strike_consensus(
            sites,
            recursive=recursive, on_dup=on_dup,
            strict=strict, verbose=verbose,
        )
        rows = []
        for _, r in tb.iterrows():
            # fake a thin band so it still renders
            for f in (1e-3, 1e3):
                rows.append(
                    dict(station=r.station, freq=f, ang=r.ang)
                )
        df = pd.DataFrame.from_records(rows)
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no strike", ha="center", va="center")
        return ax
    df = df.copy()
    df["lp"] = np.log10(1.0 / df["freq"].to_numpy())
    sts = list(df["station"].unique())
    X = []
    H = []
    for st in sts:
        s = df[df["station"] == st].sort_values("lp")
        th = s["ang"].to_numpy(dtype=float) % 180.0
        lp = s["lp"].to_numpy(dtype=float)
        h = th / 180.0
        k = max(3, int(win))
        if th.size >= k:
            vv = np.convolve(
                ((th - np.nanmean(th)) ** 2),
                np.ones(k) / k,
                mode="same",
            )
        else:
            vv = np.full_like(th, np.nan)
        v0 = (np.nanpercentile(vv, 5)
              if np.isfinite(vv).any() else 0.0)
        v1 = (np.nanpercentile(vv, 95)
              if np.isfinite(vv).any() else 1.0)
        s_sat = 1.0 - np.clip((vv - v0) / (v1 - v0 + 1e-12),
                              0.0, 1.0)
        H.append(np.vstack([h, s_sat, np.ones_like(h)]))
        X.append(lp)
    ygrid = np.unique(np.concatenate(X))
    img = np.zeros((ygrid.size, len(sts), 3))
    for j, (lp, hs) in enumerate(zip(X, H)):
        i = np.searchsorted(ygrid, lp)
        i = np.clip(i, 0, ygrid.size - 1)
        h, s, v = hs
        rgb = _hsv_rgb(h, s, v)
        for r, k in enumerate(i):
            img[k, j, :] = rgb[r]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.imshow(
        img,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(len(sts)))
    ax.set_xticklabels(sts, rotation=90)
    yt = np.linspace(0, len(ygrid) - 1, num=min(8, len(ygrid)))
    yv = np.linspace(ygrid.min(), ygrid.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    return ax


def plot_strike_profile(
    sites: Any,
    *,
    method: str = "consensus",  # consensus|sweep|pt
    band: tuple[float, float] | None = None,
    sort_by: str = "auto",  # auto|lon|lat|name
    figsize: tuple[float, float] = (8.6, 3.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if method == "sweep":
        tb = estimate_strike_sweep(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    elif method == "pt":
        tb = estimate_strike_phase_tensor(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    else:
        tb = estimate_strike_consensus(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    if tb.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no strikes", ha="center", va="center")
        return ax
    def _key(st, ed):
        if sort_by == "lon":
            x = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
            return (1, st) if x is None else (0, float(x))
        if sort_by == "lat":
            y = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
            return (1, st) if y is None else (0, float(y))
        if sort_by == "name":
            return (0, st)
        # auto: lon then name
        x = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
        return (0, float(x)) if x is not None else (1, st)
    order = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        order.append((st, _key(st, ed)))
    order = [st for st, _ in sorted(order, key=lambda t: t[1])]
    tb = tb.set_index("station").reindex(order).reset_index()
    ang = tb["ang"].to_numpy(dtype=float) % 180.0
    iq = tb["iqr"].to_numpy(dtype=float)
    x = np.arange(tb.shape[0])
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.plot(x, ang, "-", lw=1.5)
    # IQR ribbon
    lo = (ang - 0.5 * iq)
    hi = (ang + 0.5 * iq)
    ax.fill_between(x, lo, hi, alpha=0.25)
    ax.set_ylim(-5.0, 185.0)
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylabel("Strike (deg)")
    ax.set_xlabel("Station")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=90)
    ax.grid(True, alpha=0.2, which="both")
    return ax


def plot_strike_mapsticks(
    sites: Any,
    *,
    method: str = "consensus",  # consensus|sweep|pt
    band: tuple[float, float] | None = None,
    len_deg: float = 0.02,
    alpha_scale: float = 0.9,  # from confidence
    figsize: tuple[float, float] = (7.8, 6.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    if method == "sweep":
        tb = estimate_strike_sweep(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    elif method == "pt":
        tb = estimate_strike_phase_tensor(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    else:
        tb = estimate_strike_consensus(
            S, band=band, recursive=False, on_dup=on_dup,
            strict=False, verbose=verbose,
        )
    if tb.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no strikes", ha="center", va="center")
        return ax
    segs = []
    alphas = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        row = tb[tb["station"] == st]
        if row.empty:
            continue
        lat = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
        lon = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
        if lat is None or lon is None:
            continue
        a = float(row["ang"].iloc[0]) % 180.0
        c = 1.0 / (float(row["iqr"].iloc[0]) + 1e-6)
        # line segment centered at (lon,lat), axial symmetry
        th = np.radians(a)
        dx = 0.5 * len_deg * np.sin(th)
        dy = 0.5 * len_deg * np.cos(th)
        segs.append([(lon - dx, lat - dy), (lon + dx, lat + dy)])
        alphas.append(
            alpha_scale * np.clip(c / np.nanmax([c, 1.0]), 0.1, 1.0)
        )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    lc = LineCollection(
        segs, colors=[(0.1, 0.1, 0.1, a) for a in alphas],
        linewidths=2.0,
    )
    ax.add_collection(lc)
    xs = [s[0][0] for s in segs] + [s[1][0] for s in segs]
    ys = [s[0][1] for s in segs] + [s[1][1] for s in segs]
    if xs and ys:
        ax.set_xlim(min(xs) - len_deg, max(xs) + len_deg)
        ax.set_ylim(min(ys) - len_deg, max(ys) + len_deg)
    ax.set_xlabel("Lon")
    ax.set_ylabel("Lat")
    ax.set_aspect("equal", adjustable="box")
    return ax
