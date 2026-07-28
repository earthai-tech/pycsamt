from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _as_1d(a: Any) -> np.ndarray:
    v = np.asarray(a, dtype=float).ravel()
    m = np.isfinite(v) & (v > 0.0)
    return v[m]


def _prep_curve(
    misfit: Any,
    rough: Any,
    lam: Any | None,
    *,
    sort: str = "auto",  # auto|x|lambda
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = _as_1d(rough)
    y = _as_1d(misfit)
    n = int(min(x.size, y.size))
    x, y = x[:n], y[:n]
    if lam is None:
        l = np.arange(n, dtype=float)
    else:
        l = np.asarray(lam, dtype=float).ravel()[:n]
    if sort == "x":
        idx = np.argsort(x)
    elif sort == "lambda":
        idx = np.argsort(l)
    else:
        mono = np.all(np.diff(l) >= 0) or np.all(np.diff(l) <= 0)
        idx = np.argsort(l) if mono else np.argsort(x)
    return x[idx], y[idx], l[idx]


def _movavg(v: np.ndarray, w: int) -> np.ndarray:
    w = int(max(1, w))
    if v.size == 0 or w == 1:
        return v
    k = np.ones(w, dtype=float) / float(w)
    return np.convolve(v, k, mode="same")


def _curvature_kappa(
    x: np.ndarray,
    y: np.ndarray,
    *,
    smooth: int = 3,
) -> np.ndarray:
    lx = np.log10(np.maximum(x, 1e-300))
    ly = np.log10(np.maximum(y, 1e-300))
    if smooth > 1:
        lx = _movavg(lx, smooth)
        ly = _movavg(ly, smooth)
    t = np.arange(lx.size, dtype=float)
    x1 = np.gradient(lx, t, edge_order=2)
    y1 = np.gradient(ly, t, edge_order=2)
    x2 = np.gradient(x1, t, edge_order=2)
    y2 = np.gradient(y1, t, edge_order=2)
    den = (x1 * x1 + y1 * y1 + 1e-24) ** 1.5
    k = np.abs(x1 * y2 - y1 * x2) / den
    return k


def _maxdist_knee(
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    lx = np.log10(np.maximum(x, 1e-300))
    ly = np.log10(np.maximum(y, 1e-300))
    p0 = np.array([lx[0], ly[0]], dtype=float)
    p1 = np.array([lx[-1], ly[-1]], dtype=float)
    v = p1 - p0
    den = np.hypot(v[0], v[1]) + 1e-24
    # signed distance to end-point line
    d = np.abs(v[1] * (lx - p0[0]) - v[0] * (ly - p0[1])) / den
    return d


def _pick_corner(
    x: np.ndarray,
    y: np.ndarray,
    *,
    method: str = "curvature",  # curvature|maxdist
    smooth: int = 3,
    skip: int = 1,
) -> tuple[int, np.ndarray]:
    n = int(min(x.size, y.size))
    if n <= 2:
        return 0, np.zeros(n)
    if method == "maxdist":
        s = _maxdist_knee(x, y)
    else:
        s = _curvature_kappa(x, y, smooth=smooth)
    lo = int(max(0, skip))
    hi = int(max(lo + 1, n - skip))
    j = int(lo + np.nanargmax(s[lo:hi]))
    return j, s


def lcurve_table(
    misfit: Any,
    rough: Any,
    lam: Any | None = None,
    *,
    sort: str = "auto",
    method: str = "curvature",
    smooth: int = 3,
    skip: int = 1,
    return_dict: bool = False,
):
    x, y, l = _prep_curve(misfit, rough, lam, sort=sort)
    if x.size == 0:
        return pd.DataFrame() if pd is not None else None
    j, s = _pick_corner(x, y, method=method, smooth=smooth, skip=skip)
    lx = np.log10(x)
    ly = np.log10(y)
    # local slope d log(misfit) / d log(rough)
    sl = np.gradient(ly) / (np.gradient(lx) + 1e-24)
    if return_dict:
        return dict(rough=x, misfit=y, lam=l, curv=s, slope=sl, corner=j)
    df = pd.DataFrame(dict(rough=x, misfit=y, lam=l, curv=s, slope=sl))
    df.attrs["corner_idx"] = j
    return df


# ------------------------------- plotting ------------------------------- #


def plot_lcurve(
    misfit: Any | Sequence[Any],
    rough: Any | Sequence[Any],
    lam: Any | Sequence[Any] | None = None,
    *,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    cmap: str = "viridis",
    marker: str = "o",
    ms: float = 3.0,
    lw: float = 1.4,
    alpha: float = 0.9,
    show_points: bool = True,
    show_path: bool = True,
    arrow_every: int = 0,  # 0 disables
    method: str = "curvature",
    smooth: int = 3,
    skip: int = 1,
    show_corner: bool = True,
    corner_style: dict[str, Any] | None = None,
    show_inset: bool = True,
    inset_loc: tuple[float, float, float, float] = (0.62, 0.12, 0.32, 0.32),
    figsize: tuple[float, float] = (6.0, 4.6),
    ax: plt.Axes | None = None,
):
    # normalize to list of curves
    if isinstance(misfit, (list, tuple)):
        Ms = list(misfit)
        Rs = list(rough)  # type: ignore
        Ls = list(lam) if isinstance(lam, (list, tuple)) else [None] * len(Ms)
        labs = list(labels) if labels is not None else [f"C{i}" for i in range(len(Ms))]
    else:
        Ms = [misfit]
        Rs = [rough]
        Ls = [lam]
        labs = [labels[0]] if labels else ["curve"]
    if colors is None:
        cols = [None] * len(Ms)
    else:
        cols = list(colors)

    # main axes
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # store corners for legend text
    leg = []
    # Shared inset axes: created once and reused across curves, rather
    # than once per curve at the identical `inset_loc` (which used to
    # silently stack N insets on top of each other, leaving only the
    # last curve's score visible).
    axins = None

    for i, (m, r, la) in enumerate(zip(Ms, Rs, Ls)):
        x, y, l = _prep_curve(m, r, la, sort="auto")
        if x.size == 0:
            continue
        j, score = _pick_corner(x, y, method=method, smooth=smooth, skip=skip)
        c = cols[i]
        if c is None:
            c = plt.get_cmap(cmap)(0.15 + 0.65 * i / max(1, len(Ms) - 1))
        # path and points
        ax.set_xscale("log")
        ax.set_yscale("log")
        if show_path:
            ax.plot(x, y, "-", color=c, lw=lw, alpha=alpha)
        if show_points:
            ax.scatter(
                x,
                y,
                c=np.linspace(0, 1, x.size),
                s=12 * ms,
                cmap=cmap,
                edgecolors="none",
                alpha=0.85,
            )
        # optional arrows showing λ direction
        if arrow_every and x.size > 2:
            step = int(max(1, arrow_every))
            for k in range(0, x.size - step, step):
                ax.annotate(
                    "",
                    xy=(x[k + step], y[k + step]),
                    xytext=(x[k], y[k]),
                    arrowprops=dict(
                        arrowstyle="->", lw=0.8, color=c, shrinkA=0, shrinkB=0
                    ),
                )
        # corner highlight
        if show_corner:
            sty = dict(marker="*", ms=9, mec="k", mfc=c, mew=0.8)
            if corner_style:
                sty.update(corner_style)
            ax.plot([x[j]], [y[j]], **sty)
            leg.append(f"{labs[i]}  λ*≈{l[j]:.3g}")
        else:
            leg.append(labs[i])

        # inset curvature / distance -- one shared inset for all curves
        if show_inset:
            if axins is None:
                axins = fig.add_axes(inset_loc)
                axins.set_xticks([])
                axins.set_yticks([])
                axins.set_title(
                    "curv" if method == "curvature" else "knee",
                    fontsize=8,
                )
            idx = np.arange(score.size)
            axins.plot(idx, score, "-", color=c, lw=1.0)
            axins.axvline(j, color=c, ls="--", lw=1.0)

    ax.grid(True, alpha=0.25, which="both")
    ax.set_xlabel("||Lm|| (model roughness)")
    ax.set_ylabel("||Gm−d|| (data misfit)")
    if leg:
        ax.legend(leg, loc="best", frameon=False, fontsize=9)
    return ax
