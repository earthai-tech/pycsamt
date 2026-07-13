# pycsamt/emtools/skew.py
from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..api.labels import LOG10_PERIOD_LABEL, PERIOD_LABEL
from ..api.station import PYCSAMT_STATION_RENDERING
from ._core import (
    _apply_each,
    _get_t_block,
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
)
from .tensor import build_phase_tensor_table


# ---------- Bahr skewness ----------
def _z_to_2x2(z):
    z = np.asarray(z)
    if z.ndim == 3 and z.shape[-2:] == (2, 2):
        return z
    if z.ndim == 2 and z.shape[1] == 4:
        return z.reshape(-1, 2, 2)
    raise ValueError("Z must be (n,2,2) or (n,4) with Zxx,Zxy,Zyx,Zyy.")


def bahr_skewness(Z):
    Z = _z_to_2x2(Z)
    Zxx, Zxy = Z[:, 0, 0], Z[:, 0, 1]
    Zyx, Zyy = Z[:, 1, 0], Z[:, 1, 1]
    s1, s2 = Zxx + Zyy, Zxy - Zyx
    d1, d2 = Zxx - Zyy, Zxy + Zyx
    num = np.abs(s1) ** 2 + np.abs(s2) ** 2
    den = np.abs(d1) ** 2 + np.abs(d2) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = np.sqrt(num / den)
    eta[~np.isfinite(eta)] = np.nan
    return eta


def _skew_track_for(
    ed: Any, pt: pd.DataFrame
) -> tuple[np.ndarray | None, np.ndarray | None]:
    st = _name(ed, 0)
    Z, z, fr = _get_z_block(ed)
    if Z is None:
        return None, None
    sdf = pt[pt["station"] == st]
    if sdf.empty:
        return fr, np.full(fr.size, np.nan, dtype=float)
    per = 1.0 / fr
    p_ref = sdf["period"].to_numpy(dtype=float)
    sk_ref = sdf["skew"].to_numpy(dtype=float)
    idx = np.searchsorted(p_ref, per)
    idx = np.clip(idx, 0, p_ref.size - 1)
    sk = sk_ref[idx]
    return fr, sk


def _mask_apply(
    ed: Any,
    keep: np.ndarray,
    *,
    also: str = "both",  # z|tipper|both
) -> None:
    Z, z, fr = _get_z_block(ed)
    if Z is not None and z is not None:
        z2 = z.copy()
        z2[~keep] = np.nan
        try:
            Z.z = z2
        except:
            pass
    if also in ("tipper", "both"):
        T, t, ft = _get_t_block(ed)
        if T is not None and t is not None:
            t2 = t.copy()
            t2[~keep] = np.nan
            try:
                T.tipper = t2
            except:
                pass


def _runs_bool(m: np.ndarray) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    if m.size == 0:
        return out
    s = None
    for i, v in enumerate(m):
        if v and s is None:
            s = i
        if (not v) and s is not None:
            out.append((s, i - 1))
            s = None
    if s is not None:
        out.append((s, m.size - 1))
    return out


def _fill_small_gaps(m: np.ndarray, max_gap: int) -> np.ndarray:
    if max_gap <= 0 or m.size == 0:
        return m
    mi = m.astype(int)
    d = np.diff(np.pad(mi, (1, 1)))
    # rising/falling edges
    on = np.where(d == 1)[0]
    off = np.where(d == -1)[0]
    if on.size == 0 or off.size == 0:
        return m
    # ensure paired
    if off[0] < on[0]:
        off = off[1:]
    n = min(on.size, off.size)
    on, off = on[:n], off[:n]
    out = mi.copy()
    for a, b in zip(off[:-1], on[1:]):
        gap = b - a
        if 0 < gap <= max_gap:
            out[a:b] = 1
    return out.astype(bool)


def skew_table(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    pt = build_phase_tensor_table(
        S,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    return pt


def mask_by_skew(
    sites: Any,
    *,
    thresh: float = 6.0,
    mode: str = "abs_gt",  # abs_gt|gt|lt|abs_lt
    also: str = "both",
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
    pt = build_phase_tensor_table(
        S,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        fr, sk = _skew_track_for(ed, pt)
        if fr is None or sk is None:
            return Si
        a = np.abs(sk)
        if mode == "gt":
            keep = sk <= thresh
        elif mode == "lt":
            keep = sk >= thresh
        elif mode == "abs_lt":
            keep = a <= thresh
        else:
            keep = (
                a <= thresh
                if np.isfinite(thresh)
                else np.ones(sk.size, dtype=bool)
            )
        _mask_apply(ed, keep, also=also)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def keep_longest_low_skew(
    sites: Any,
    *,
    thresh: float = 3.0,
    min_len: int = 3,
    pad: int = 0,
    also: str = "both",
    fallback: str = "keep_all",  # keep_all|drop_all
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    inplace: bool = False,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    pt = build_phase_tensor_table(
        S,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        fr, sk = _skew_track_for(ed, pt)
        if fr is None or sk is None:
            return Si
        good = np.isfinite(sk) & (np.abs(sk) <= thresh)
        runs = _runs_bool(good)
        if runs:
            # longest run
            lens = [j - i + 1 for (i, j) in runs]
            k = int(np.argmax(lens))
            i0, i1 = runs[k]
            if lens[k] < min_len:
                if fallback == "drop_all":
                    keep = np.zeros_like(good, dtype=bool)
                else:
                    keep = np.ones_like(good, dtype=bool)
            else:
                i0 = max(0, i0 - pad)
                i1 = min(good.size - 1, i1 + pad)
                keep = np.zeros_like(good, dtype=bool)
                keep[i0 : i1 + 1] = True
        else:
            keep = (
                np.ones_like(good, dtype=bool)
                if fallback == "keep_all"
                else np.zeros_like(good, dtype=bool)
            )
        _mask_apply(ed, keep, also=also)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def close_skew_gaps(
    sites: Any,
    *,
    thresh: float = 3.0,
    max_gap: int = 1,
    also: str = "both",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    inplace: bool = False,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    pt = build_phase_tensor_table(
        S,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )

    def _one(Si):
        ed = next(_iter_items(Si))
        fr, sk = _skew_track_for(ed, pt)
        if fr is None or sk is None:
            return Si
        good = np.isfinite(sk) & (np.abs(sk) <= thresh)
        keep = _fill_small_gaps(good, max_gap=max_gap)
        _mask_apply(ed, keep, also=also)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


def select_low_skew_band(
    sites: Any,
    *,
    thresh: float = 3.0,
    frac: float = 0.6,
    min_len: int = 3,
    pad: int = 0,
    also: str = "both",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    inplace: bool = False,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # build per-site keep masks; then intersect by fraction
    pt = build_phase_tensor_table(
        S,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    # gather masks on each site's own grid
    masks = []
    grids = []
    items = []
    for ed in _iter_items(S):
        fr, sk = _skew_track_for(ed, pt)
        if fr is None or sk is None:
            continue
        good = np.isfinite(sk) & (np.abs(sk) <= thresh)
        # longest run (optional min_len/pad)
        runs = _runs_bool(good)
        if runs:
            lens = [j - i + 1 for (i, j) in runs]
            k = int(np.argmax(lens))
            i0, i1 = runs[k]
            if lens[k] >= min_len:
                i0 = max(0, i0 - pad)
                i1 = min(good.size - 1, i1 + pad)
                keep = np.zeros_like(good, dtype=bool)
                keep[i0 : i1 + 1] = True
            else:
                keep = good
        else:
            keep = good
        masks.append(keep)
        grids.append(fr)
        items.append(ed)
    if not masks:
        return S
    # vote on union grid
    G = np.unique(np.concatenate(grids))
    vote = np.zeros(G.size, dtype=float)
    for fr, m in zip(grids, masks):
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        vote[idx] += m.astype(float)
    keep_union = vote >= (frac * len(masks))

    # apply per site by nearest
    def _one(Si):
        ed = next(_iter_items(Si))
        Z, z, fr = _get_z_block(ed)
        if Z is None or fr is None:
            return Si
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        keep = keep_union[idx]
        _mask_apply(ed, keep, also=also)
        return Si

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# --- BEAUTIFUL SKEW VIEWS ------------------------------------------------ #

# 1) Traffic-light pseudosection (green/amber/red + alpha=confidence)


def plot_skew_traffic_psection(
    sites: Any,
    *,
    t1: float = 3.0,
    t2: float = 6.0,
    figsize: tuple[float, float] = (9.0, 4.8),
    axis_y: str = "logperiod",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    df = build_phase_tensor_table(
        S,
        recursive=False,
        on_dup=on_dup,
        strict=False,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor", ha="center", va="center")
        return ax
    df = df.copy()
    b = np.abs(df["beta"].to_numpy(dtype=float))
    cls = np.full(b.size, 1, dtype=int)
    cls[b <= t1] = 0
    cls[b > t2] = 2
    # confidence: margin from nearest boundary, 0..1
    m0 = np.maximum(t1 - b, 0.0)
    m1 = np.maximum(np.minimum(b - t1, t2 - b), 0.0)
    m2 = np.maximum(b - t2, 0.0)
    conf = np.zeros_like(b, dtype=float)
    conf[cls == 0] = m0[cls == 0]
    conf[cls == 1] = m1[cls == 1]
    conf[cls == 2] = m2[cls == 2]
    if np.isfinite(conf).any():
        v0 = np.nanpercentile(conf, 5)
        v1 = np.nanpercentile(conf, 95)
        conf = np.clip((conf - v0) / (v1 - v0 + 1e-12), 0.0, 1.0)
    else:
        conf[:] = 1.0
    df["cls"] = cls
    df["conf"] = conf
    if axis_y == "logperiod":
        df["yy"] = np.log10(df["period"].to_numpy())
        ylab = LOG10_PERIOD_LABEL
    else:
        df["yy"] = df["period"].to_numpy()
        ylab = PERIOD_LABEL
    sts = list(df["station"].unique())
    sidx = {s: i for i, s in enumerate(sts)}
    X = df["station"].map(sidx).to_numpy(dtype=int)
    Y = df["yy"].to_numpy(dtype=float)
    # grid in y
    yall = np.unique(Y)
    H = np.zeros((yall.size, len(sts), 4))
    pal = {
        0: (0.20, 0.60, 0.20),  # green
        1: (0.95, 0.70, 0.20),  # amber
        2: (0.85, 0.25, 0.20),  # red
    }
    for x, y, c, a in zip(X, Y, cls, conf):
        yi = int(np.searchsorted(yall, y))
        r, g, b = pal[int(c)]
        if a >= H[yi, x, 3]:
            H[yi, x, :] = (r, g, b, a)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.imshow(
        H,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_ylabel(ylab)
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(len(sts), dtype=float),
        sts,
        preset="pseudosection",
        xlim=(-0.5, len(sts) - 0.5),
    )
    yt = np.linspace(0, yall.size - 1, num=min(8, yall.size))
    yv = np.linspace(yall.min(), yall.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    return ax


# 2) Percentile ribbon (line-level summary of |beta| vs period)


def plot_skew_percentile_ribbon(
    sites: Any,
    *,
    n_bins: int = 30,
    q_lo: float = 25.0,
    q_hi: float = 75.0,
    extra: tuple[float, float] | None = (10.0, 90.0),
    figsize: tuple[float, float] = (8.6, 3.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor", ha="center", va="center")
        return ax
    p = df["period"].to_numpy(dtype=float)
    b = np.abs(df["beta"].to_numpy(dtype=float))
    lp = np.log10(np.maximum(p, 1e-9))
    lo, hi = float(np.nanmin(lp)), float(np.nanmax(lp))
    edges = np.linspace(lo, hi, int(max(8, n_bins)) + 1)
    cen = 0.5 * (edges[1:] + edges[:-1])
    Q1 = np.full(cen.size, np.nan)
    Q2 = np.full(cen.size, np.nan)
    Q3 = np.full(cen.size, np.nan)
    E1 = np.full(cen.size, np.nan)
    E2 = np.full(cen.size, np.nan)
    for i in range(cen.size):
        m = (lp >= edges[i]) & (lp < edges[i + 1])
        if not np.any(m):
            continue
        Q1[i] = np.nanpercentile(b[m], q_lo)
        Q2[i] = np.nanpercentile(b[m], 50.0)
        Q3[i] = np.nanpercentile(b[m], q_hi)
        if extra is not None:
            E1[i] = np.nanpercentile(b[m], extra[0])
            E2[i] = np.nanpercentile(b[m], extra[1])
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    x = 10**cen
    ax.set_xscale("log")
    if extra is not None:
        ax.fill_between(x, E1, E2, alpha=0.15)
    ax.fill_between(x, Q1, Q3, alpha=0.35)
    ax.plot(x, Q2, "-", lw=2.0)
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("|beta| (deg)")
    ax.grid(True, alpha=0.25, which="both")
    return ax


# 3) Vote-band curve (fraction of sites with |beta| <= t)


def plot_skew_vote_band(
    sites: Any,
    *,
    thresh: float = 3.0,
    n_bins: int = 40,
    figsize: tuple[float, float] = (8.6, 3.4),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    df = build_phase_tensor_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if df.empty:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no phase tensor", ha="center", va="center")
        return ax
    # bin in log-period; vote per station in each bin
    p = df["period"].to_numpy(dtype=float)
    b = np.abs(df["beta"].to_numpy(dtype=float))
    lp = np.log10(np.maximum(p, 1e-9))
    lo, hi = float(np.nanmin(lp)), float(np.nanmax(lp))
    edges = np.linspace(lo, hi, int(max(8, n_bins)) + 1)
    cen = 0.5 * (edges[1:] + edges[:-1])
    frac = np.zeros(cen.size, dtype=float)
    sts = df["station"].astype(str).unique().tolist()
    for i in range(cen.size):
        m = (lp >= edges[i]) & (lp < edges[i + 1])
        if not np.any(m):
            frac[i] = np.nan
            continue
        # count per-station pass/fail inside bin
        S = df.loc[m, ["station"]].copy()
        S["ok"] = b[m] <= thresh
        g = S.groupby("station")["ok"].mean() > 0.5
        frac[i] = g.sum() / max(1, len(sts))
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    x = 10**cen
    ax.set_xscale("log")
    ax.plot(x, frac, "-", lw=2.0)
    ax.fill_between(x, 0.0, frac, alpha=0.25)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("fraction |beta| ≤ thresh")
    ax.grid(True, alpha=0.25, which="both")
    return ax


def plot_skewness(f_hz, Z, *, threshold=0.4, ax=None, title=None):
    if ax is None:
        ax = plt.gca()

    f = np.asarray(f_hz, float)
    if f.size == 0:
        raise ValueError("Empty frequency array.")

    T = 1.0 / f
    x = np.log10(T)
    eta = bahr_skewness(Z)

    ax.plot(x, eta, "+", ms=3, mew=0.9, color="0.35", label=None)
    ax.axhline(
        threshold, color="red", lw=1.5, label=f"Threshold: η = {threshold:g}"
    )

    # Average skewness annotation
    avg_eta = float(np.nanmean(eta))
    txt = f"Aver. skewness: Bahr = {avg_eta:.3f}"
    ax.text(
        0.02,
        0.95,
        txt,
        transform=ax.transAxes,
        ha="left",
        va="top",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black"),
    )

    ax.set_ylabel("Skewness (η)")
    ax.set_xlabel(r"$\log_{10}$ Period (s)")
    if title:
        ax.set_title(title)

    # Simple “2D / 3D” guide on the right
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xr = xmax + 0.02 * (xmax - xmin)
    ax.plot([xmax, xmax], [0, threshold], color="tab:blue")
    ax.plot([xmax, xmax], [threshold, ymax], color="tab:orange")
    ax.text(xr, 0.5 * threshold, "2D", va="center", color="tab:blue")
    ax.text(
        xr, 0.5 * (threshold + ymax), "3D", va="center", color="tab:orange"
    )
    ax.legend(loc="upper right", frameon=True)
    ax.grid(True, alpha=0.25)
    return ax
