from __future__ import annotations

from typing import Any, List, Optional, Sequence, Tuple
import numpy as np
import matplotlib.pyplot as plt

from ._core import (
    ensure_sites,
    _iter_items,
    _get_t_block,
    _get_z_block,
    _name,
)

from .tensor import build_phase_tensor_table


# ------------------------------- helpers -------------------------------- #

def _pick_station(S, station: Optional[str]) -> Tuple[str, Any]:
    pool = {}
    for i, ed in enumerate(_iter_items(S)):
        pool[_name(ed, i)] = ed
    if not pool:
        raise RuntimeError("no sites")
    if station is None:
        station = sorted(pool.keys())[0]
    ed = pool.get(station, None)
    if ed is None:
        raise RuntimeError("station not found")
    return station, ed


def _bands_from_periods(
    per: np.ndarray,
    *,
    bands: Optional[Sequence[Tuple[float, float]]] = None,
    n_bands: int = 3,
) -> List[np.ndarray]:
    if bands:
        ms = []
        for (lo, hi) in bands:
            ms.append((per >= float(lo)) & (per <= float(hi)))
        return ms
    # quantile bands
    q = np.linspace(0, 1, num=int(n_bands) + 1)
    bb = np.quantile(per, q)
    ms = []
    for i in range(len(bb) - 1):
        lo, hi = bb[i], bb[i + 1]
        # include low edge, exclude high (except last)
        if i == len(bb) - 2:
            ms.append((per >= lo) & (per <= hi))
        else:
            ms.append((per >= lo) & (per < hi))
    return ms


def _station_xy(ed: Any, i: int) -> Tuple[float, float]:
    for kx, ky in [
        ("east", "north"), ("easting", "northing"),
        ("x", "y"), ("lon", "lat"),
    ]:
        x = getattr(ed, kx, None); y = getattr(ed, ky, None)
        if isinstance(x, (int, float)) and isinstance(y, (int, float)):
            return float(x), float(y)
    return float(i), 0.0  # fallback: index on a line


def _nearest_idx(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    j = np.searchsorted(x, y)
    j = np.clip(j, 1, x.size - 1)
    use_left = np.abs(y - x[j - 1]) <= np.abs(y - x[j])
    j[use_left] -= 1
    return j


def _pt_angle(
    S, station: str, per: np.ndarray
) -> Optional[np.ndarray]:
    if build_phase_tensor_table is None:
        return None
    tb = build_phase_tensor_table(
        S, recursive=False, on_dup="replace",
        strict=False, verbose=0,
    )
    if getattr(tb, "empty", False):
        return None
    sdf = tb[tb["station"] == station]
    if sdf.empty:
        return None
    p_ref = sdf["period"].to_numpy(dtype=float)
    for col in ("azimuth", "strike", "phi", "theta"):
        if col in sdf.columns:
            phi = sdf[col].to_numpy(dtype=float)
            j = _nearest_idx(p_ref, per)
            return phi[j]
    return None


# ------------------------- 11) Tipper hodograms ------------------------- #

def plot_tipper_hodograms(
    sites: Any,
    *,
    station: Optional[str] = None,
    bands: Optional[Sequence[Tuple[float, float]]] = None,
    n_bands: int = 3,
    normalize: bool = False,
    colors: Optional[Sequence[str]] = None,
    marker: str = "o",
    ms: float = 3.0,
    lw: float = 1.0,
    ls: str = "-",
    unit_circle: bool = True,
    figsize: Tuple[float, float] = (6.4, 3.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    st, ed = _pick_station(S, station)
    T, t, fr = _get_t_block(ed)
    if T is None or t is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "no tipper", ha="center", va="center")
        return fig
    per = 1.0 / fr
    Ms = _bands_from_periods(per, bands=bands, n_bands=n_bands)
    if colors is None:
        cols = [plt.cm.viridis(c) for c in np.linspace(0, 1, len(Ms))]
    else:
        cols = list(colors)[: len(Ms)]
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, wspace=0.25)
    axX = fig.add_subplot(gs[0]); axY = fig.add_subplot(gs[1])
    for k, m in enumerate(Ms):
        if not np.any(m):
            continue
        tx = t[m, 0]; ty = t[m, 1]
        X = np.real(tx); Y = np.imag(tx)
        U = np.real(ty); V = np.imag(ty)
        if normalize:
            s = np.nanpercentile(
                np.hypot(np.r_[X, U], np.r_[Y, V]), 95
            ) + 1e-24
            X, Y, U, V = X / s, Y / s, U / s, V / s
        col = cols[k]
        axX.plot(X, Y, ls=ls, lw=lw, color=col)
        axX.scatter(X, Y, s=12 * ms, color=col)
        axY.plot(U, V, ls=ls, lw=lw, color=col)
        axY.scatter(U, V, s=12 * ms, color=col)
    for ax, lab in [(axX, "Tx"), (axY, "Ty")]:
        ax.axhline(0, color="0.85", lw=0.8)
        ax.axvline(0, color="0.85", lw=0.8)
        if unit_circle:
            th = np.linspace(0, 2 * np.pi, 256)
            ax.plot(np.cos(th), np.sin(th), ":", color="0.7", lw=0.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Real")
        ax.set_ylabel("Imag")
        ax.set_title(f"{st} • {lab}")
        ax.grid(True, alpha=0.25)
    return fig


# --------- 12) Induction arrows + phase-tensor strike overlay ----------- #

def _arrow_from_tipper(
    t: np.ndarray,
    *,
    convention: str = "park",  # park|wiese|real|imag
) -> np.ndarray:
    tx = t[:, 0]; ty = t[:, 1]
    if convention == "real":
        vx = np.real(tx); vy = np.real(ty)
    else:
        # default: use imaginary (Parkinson-style)
        vx = -np.imag(tx); vy = -np.imag(ty)
        if convention == "wiese":
            # 90° rotation of Parkinson vector
            vx, vy = -vy, vx
        if convention == "imag":
            pass  # already imag-based
    return np.vstack([vx, vy]).T


def plot_induction_arrows(
    sites: Any,
    *,
    periods: Sequence[float] = (1.0,),
    convention: str = "park",   # park|wiese|real|imag
    scale: float = 1.0,         # quiver scale factor
    normalize: bool = True,
    strike_ticks: bool = True,
    tick_len: float = 0.25,
    figsize: Tuple[float, float] = (7.2, 3.4),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    # collect site positions and arrows per requested period
    sts, XY, AR = [], [], []
    per_layers: List[Tuple[float, np.ndarray]] = []
    for p in periods:
        xs, ys, u, v = [], [], [], []
        for i, ed in enumerate(_iter_items(S)):
            T, t, fr = _get_t_block(ed)
            if T is None or t is None:
                continue
            per = 1.0 / fr
            j = _nearest_idx(per, np.array([p], float))
            tx = t[j[0]]
            vec = _arrow_from_tipper(
                np.asarray([tx]), convention=convention
            )[0]
            x, y = _station_xy(ed, i)
            xs.append(x); ys.append(y)
            u.append(vec[0]); v.append(vec[1])
            if p == periods[0]:
                sts.append(_name(ed, i))
        per_layers.append(
            (p, np.vstack([np.array(xs), np.array(ys),
                           np.array(u), np.array(v)]))
        )
    if not per_layers:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no tipper", ha="center", va="center")
        return ax

    # normalization
    if normalize:
        all_mag = []
        for _, L in per_layers:
            mag = np.hypot(L[2], L[3])
            all_mag.append(mag)
        m95 = np.nanpercentile(np.hstack(all_mag), 95) + 1e-24
        for k in range(len(per_layers)):
            p, L = per_layers[k]
            L[2] /= m95; L[3] /= m95
            per_layers[k] = (p, L)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    # axis as map or profile depending on coords
    XY0 = per_layers[0][1]
    if np.unique(XY0[1]).size == 1:  # one line
        ax.set_xlabel("Station index / x")
        ax.set_ylabel("Arrow (arb.)")
    else:
        ax.set_xlabel("East / lon"); ax.set_ylabel("North / lat")

    # draw layers, one color per period
    cols = plt.cm.viridis(
        np.linspace(0.15, 0.85, len(per_layers))
    )
    for (c, (p, L)) in zip(cols, per_layers):
        ax.quiver(
            L[0], L[1], L[2] * scale, L[3] * scale,
            angles="xy", scale_units="xy", scale=1.0,
            color=c, width=0.003, headlength=4, headaxislength=3,
            minlength=0.0,
        )
    # optional strike ticks from PT
    if strike_ticks and build_phase_tensor_table is not None:
        # one average angle per site (first period layer’s sites)
        per = np.array(periods, dtype=float)
        # use the first period to sample strike
        p0 = float(per[0])
        labs = sts
        th_map = {}
        for i, ed in enumerate(_iter_items(S)):
            st = _name(ed, i)
            if st not in labs:
                continue
            Z, z, fr = _get_z_block(ed)[:3]
            if Z is None:
                continue
            th = _pt_angle(S, st, np.array([p0]))
            if th is None:
                continue
            th_map[st] = float(np.radians(th[0]))
        for i, ed in enumerate(_iter_items(S)):
            st = _name(ed, i)
            if st not in th_map:
                continue
            x, y = _station_xy(ed, i)
            th = th_map[st]
            dx = tick_len * np.cos(th)
            dy = tick_len * np.sin(th)
            ax.plot([x - dx, x + dx], [y - dy, y + dy],
                    "-", color="0.2", lw=1.2, alpha=0.9)

    # cosmetics
    if np.unique(XY0[1]).size == 1:
        # profile: center around zero vertically
        ax.axhline(0.0, color="0.8", lw=0.8)
    ax.grid(True, alpha=0.25)
    leg = [f"P={p:g}s" for p, _ in per_layers]
    fig = ax.figure
    fig.legend(
        [plt.Line2D([], [], color=c, lw=2.0) for c in cols],
        leg, ncol=min(len(leg), 4), frameon=False,
        loc="upper center", bbox_to_anchor=(0.5, 1.02), fontsize=9,
    )
    return ax
