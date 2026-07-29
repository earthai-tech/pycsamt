from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from ..api.labels import LOG10_PERIOD_LABEL
from ..api.station import PYCSAMT_STATION_RENDERING
from ._core import (
    _axes_list,
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
    hide_polar_radius_labels,
)

# default colors consistent with plot.py
_COL = {
    "xy": "#1f77b4",  # blue
    "yx": "#d62728",  # red
    "xx": "#2ca02c",  # green
    "yy": "#9467bd",  # purple
}


# --------------------------- local helpers ------------------------------ #


def _zblk_flex(ed: Any, need_err: bool = False):
    try:
        return _get_z_block(ed, with_errors=need_err)
    except TypeError:
        return _get_z_block(ed)


def _comp(z: np.ndarray, comp: str) -> np.ndarray:
    i = 0 if comp[0] == "x" else 1
    j = 0 if comp[1] == "x" else 1
    return z[:, i, j]


def _logper(fr: np.ndarray) -> np.ndarray:
    return np.log10(np.maximum(1.0 / fr, 1e-24))


# ---------------------- 7) Impedance phasor wheel ----------------------- #


def plot_phasor_wheel(
    sites: Any,
    *,
    station: str | None = None,
    components: tuple[str, ...] = ("xy", "yx"),
    pband: tuple[float, float] | None = None,
    radius: str = "abs",  # abs|norm
    cmap: str = "viridis",
    colors: dict[str, str] | None = None,
    marker: str = "o",
    ms: float = 3.0,
    lw: float = 1.0,
    connect: bool = True,
    figsize: tuple[float, float] = (4.8, 4.8),
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
    cmap_cols = {**_COL, **(colors or {})}
    sel = {}
    for i, ed in enumerate(_iter_items(S)):
        sel[_name(ed, i)] = ed
    if not sel:
        if ax is None:
            _, ax = plt.subplots(subplot_kw={"polar": True})
        hide_polar_radius_labels(ax)
        ax.text(0.5, 0.5, "no sites", ha="center", va="center")
        return ax
    if station is None:
        station = sorted(sel.keys())[0]
    ed = sel.get(station, None)
    if ed is None:
        if ax is None:
            _, ax = plt.subplots(subplot_kw={"polar": True})
        hide_polar_radius_labels(ax)
        ax.text(0.5, 0.5, "station not found", ha="center", va="center")
        return ax
    Z, z, fr = _zblk_flex(ed)[:3]
    if Z is None:
        if ax is None:
            _, ax = plt.subplots(subplot_kw={"polar": True})
        hide_polar_radius_labels(ax)
        ax.text(0.5, 0.5, "no Z", ha="center", va="center")
        return ax

    per = 1.0 / fr
    m = np.ones(fr.size, dtype=bool)
    if pband is not None:
        lo, hi = pband
        m &= (per >= lo) & (per <= hi)

    angN = {}
    radN = {}
    for c in components:
        zz = _comp(z, c)
        ang = np.angle(zz)  # radians
        rad = np.abs(zz)
        if radius == "norm":
            s = np.nanpercentile(rad[m], 95) + 1e-24
            rad = rad / s
        angN[c] = ang[m]
        radN[c] = rad[m]

    th = _logper(fr[m])
    th = (th - th.min()) / (th.max() - th.min() + 1e-12)
    th = 2.0 * np.pi * th  # angle ~ log-period

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, subplot_kw={"polar": True})
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(-1)
    hide_polar_radius_labels(ax)
    for c in components:
        col = cmap_cols.get(c, "k")
        ax.scatter(
            angN[c],
            radN[c],
            c=th,
            cmap=cmap,
            s=10 * ms,
            alpha=0.9,
            edgecolors="none",
            label=f"Z{c.upper()}",
        )
        if connect:
            # sort by theta for a gentle path
            k = np.argsort(th)
            ax.plot(
                angN[c][k],
                radN[c][k],
                "-",
                lw=lw,
                color=col,
                alpha=0.6,
            )
    ax.grid(True, alpha=0.25)
    hide_polar_radius_labels(ax)
    ax.set_title(str(station), pad=8)
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0.02, 0.02),
        frameon=False,
        fontsize=8,
    )
    return ax


# --- 8) Off-diagonal antisymmetry residual pseudosection --------------- #


def plot_offdiag_antisym_residual(
    sites: Any,
    *,
    vlim: float | None = None,
    cmap: str = "magma",
    figsize: tuple[float, float] = (9.0, 4.8),
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
    labs: list[str] = []
    Gs: list[np.ndarray] = []
    Vs: list[np.ndarray] = []
    for i, ed in enumerate(_iter_items(S)):
        Z, z, fr = _zblk_flex(ed)[:3]
        if Z is None:
            continue
        xy = np.abs(z[:, 0, 1])
        yx = np.abs(z[:, 1, 0])
        num = np.abs(z[:, 0, 1] + z[:, 1, 0])
        den = xy + yx + 1e-24
        r = np.clip(num / den, 0.0, 1.0)
        labs.append(_name(ed, i))
        Gs.append(fr)
        Vs.append(r)
    if not labs:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    G = np.unique(np.concatenate(Gs))
    M = np.full((len(labs), G.size), np.nan, dtype=float)
    for i, (fr, v) in enumerate(zip(Gs, Vs)):
        idx = np.searchsorted(G, fr)
        idx = np.clip(idx, 0, G.size - 1)
        M[i, idx] = v
        # fill small gaps
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
    ZI = M.T  # (freq, station)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    lp = _logper(G)
    order = np.argsort(lp)
    lp = lp[order]
    ZI = ZI[order]
    v = ZI[np.isfinite(ZI)]
    if vlim is None and v.size:
        vlim = float(max(0.2, np.nanpercentile(v, 95)))
    im = ax.imshow(
        ZI,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        vmin=0.0,
        vmax=(vlim or 1.0),
    )
    ax.set_ylabel(LOG10_PERIOD_LABEL)
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(len(labs), dtype=float),
        labs,
        preset="pseudosection",
        xlim=(-0.5, len(labs) - 0.5),
    )
    yt = np.linspace(0, len(lp) - 1, num=min(8, len(lp)))
    yv = np.linspace(lp.min(), lp.max(), num=yt.size)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{v:.2g}" for v in yv])
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("|Zxy+Zyx|/(|Zxy|+|Zyx|)")
    return ax


# ----------------------- 9) Determinant track (per site) ---------------- #


def _det_ci(
    z: np.ndarray,
    fr: np.ndarray,
    ze: np.ndarray | None,
    *,
    pcts: tuple[float, float, float] = (10.0, 50.0, 90.0),
    n_draws: int = 200,
    seed: int | None = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a = z[:, 0, 0]
    b = z[:, 0, 1]
    c = z[:, 1, 0]
    d = z[:, 1, 1]
    det = a * d - b * c
    if ze is None:
        mag = np.abs(det)
        ph = np.degrees(np.angle(det))
        return mag, ph, np.vstack([mag, mag, mag]).T
    ea = ze[:, 0, 0]
    eb = ze[:, 0, 1]
    ec = ze[:, 1, 0]
    ed = ze[:, 1, 1]
    rng = np.random.default_rng(seed)
    n = int(max(16, n_draws))
    nf = det.size

    def E(s):
        return (
            (rng.standard_normal((n, nf)) + 1j * rng.standard_normal((n, nf)))
            / np.sqrt(2.0)
        ) * s[None, :]

    Ad = a[None, :] + E(ea)
    Bd = b[None, :] + E(eb)
    Cd = c[None, :] + E(ec)
    Dd = d[None, :] + E(ed)
    Dt = Ad * Dd - Bd * Cd
    mag = np.nanmedian(np.abs(Dt), axis=0)
    ph = np.nanmedian(np.degrees(np.angle(Dt)), axis=0)
    mag_lo = np.nanpercentile(np.abs(Dt), pcts[0], axis=0)
    mag_hi = np.nanpercentile(np.abs(Dt), pcts[2], axis=0)
    return mag, ph, np.vstack([mag_lo, mag_hi]).T


def plot_determinant_track(
    sites: Any,
    *,
    station: str | None = None,
    pband: tuple[float, float] | None = None,
    pcts: tuple[float, float, float] = (10.0, 50.0, 90.0),
    n_draws: int = 200,
    height_ratio: tuple[int, int] = (2, 1),
    axes=None,
    figsize: tuple[float, float] = (6.4, 3.8),
    color_mag: str = "C0",
    color_phase: str = "C3",
    fill_alpha: float = 0.20,
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
    sel = {}
    for i, ed in enumerate(_iter_items(S)):
        sel[_name(ed, i)] = ed
    axes_given = _axes_list(axes, 1) if axes is not None else None
    if not sel:
        if axes_given is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no sites", ha="center", va="center")
        return fig
    if station is None:
        station = sorted(sel.keys())[0]
    ed = sel.get(station, None)
    if ed is None:
        if axes_given is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "station not found", ha="center", va="center")
        return fig
    out = _zblk_flex(ed, need_err=True)
    if len(out) == 4:
        _, z, fr, ze = out
    else:
        _, z, fr = out[:3]
        ze = None
    if z is None or fr is None:
        if axes_given is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        else:
            ax = axes_given[0]
            fig = ax.figure
        ax.text(0.5, 0.5, "no Z", ha="center", va="center")
        return fig

    per = 1.0 / fr
    m = np.ones(fr.size, dtype=bool)
    if pband is not None:
        lo, hi = pband
        m &= (per >= lo) & (per <= hi)

    mag, ph, band = _det_ci(
        z[m],
        fr[m],
        ze[m] if isinstance(ze, np.ndarray) else None,
        pcts=pcts,
        n_draws=n_draws,
    )
    x = per[m]
    axes_given = _axes_list(axes, 2) if axes is not None else None
    if axes_given is None:
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 1, height_ratios=height_ratio, hspace=0.06)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
    else:
        ax1, ax2 = axes_given
        fig = ax1.figure
    ax1.set_xscale("log")
    ax2.set_xscale("log")
    ax1.fill_between(x, band[:, 0], band[:, 1], color=color_mag, alpha=fill_alpha)
    ax1.plot(x, mag, "-", color=color_mag, lw=1.6)
    ax2.plot(x, ph, "-", color=color_phase, lw=1.6)
    ax1.grid(True, alpha=0.25, which="both")
    ax2.grid(True, alpha=0.25, which="both")
    ax2.set_xlabel("Period (s)")
    ax1.set_ylabel("|det(Z)|")
    ax2.set_ylabel("Phase det(Z) (°)")
    ax1.set_title(str(station), pad=6)
    return fig
