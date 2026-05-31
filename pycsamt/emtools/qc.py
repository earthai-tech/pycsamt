
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as _Rect

from ._core import (
    ensure_sites,
    _iter_items,
    _get_z_block,
    _get_t_block,
    _name,
    _station_positions,
)
from .tensor import build_phase_tensor_table


# ------------------------------ helpers --------------------------------- #

def _row_ok_z(z: np.ndarray) -> np.ndarray:
    y = z.reshape(z.shape[0], -1)
    return np.isfinite(y).all(axis=1)


def _row_ok_t(t: np.ndarray) -> np.ndarray:
    y = t.reshape(t.shape[0], -1)
    return np.isfinite(y).all(axis=1)


def _snr_rows(z: np.ndarray, ze: Optional[np.ndarray]) -> np.ndarray:
    if ze is None:
        return np.full(z.shape[0], np.nan, dtype=float)
    a = np.sqrt(np.nanmean(np.abs(z) ** 2, axis=(1, 2)))
    e = np.sqrt(np.nanmean(np.abs(ze) ** 2, axis=(1, 2)))
    return a / (e + 1e-12)


def _offdiag_logmag(z: np.ndarray) -> np.ndarray:
    m = np.nanmedian(
        np.stack([np.abs(z[:, 0, 1]), np.abs(z[:, 1, 0])], axis=1),
        axis=1,
    )
    return np.log10(np.maximum(m, 1e-24))


def _y_ticks(yall: np.ndarray, ny: int) -> Tuple[np.ndarray, List[str]]:
    yt = np.linspace(0, yall.size - 1, num=min(ny, yall.size))
    yv = np.linspace(yall.min(), yall.max(), num=yt.size)
    lab = [f"{v:.2g}" for v in yv]
    return yt, lab


# ------------------------------ tables ---------------------------------- #

def build_qc_table(
    sites: Any,
    *,
    include_skew: bool = True,
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
    pt = None
    if include_skew:
        pt = build_phase_tensor_table(
            S,
            recursive=False,
            on_dup=on_dup,
            strict=False,
            verbose=verbose,
        )
    rows: List[Dict[str, Any]] = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        T, t, ft = _get_t_block(ed)
        if Z is None:
            continue
        n = z.shape[0]
        ko = _row_ok_z(z)
        ze = getattr(Z, "z_err", None)
        snr = _snr_rows(z, ze)
        med_snr = float(np.nanmedian(snr))
        n_ok = int(np.nansum(ko))
        n_t = 0
        n_t_ok = 0
        if T is not None and t is not None:
            n_t = t.shape[0]
            n_t_ok = int(np.nansum(_row_ok_t(t)))
        per = 1.0 / fr
        pmin = float(np.nanmin(per)) if per.size else np.nan
        pmax = float(np.nanmax(per)) if per.size else np.nan
        rec = dict(
            station=st,
            n_freq=int(n),
            n_ok=int(n_ok),
            frac_ok=float(n_ok / max(1, n)),
            n_tip=int(n_t),
            n_tip_ok=int(n_t_ok),
            snr_med=med_snr,
            pmin=pmin,
            pmax=pmax,
        )
        if include_skew and pt is not None:
            sdf = pt[pt["station"] == st]
            if not sdf.empty:
                sb = np.abs(sdf["beta"].to_numpy(dtype=float))
                rec["skew_med"] = float(np.nanmedian(sb))
                rec["skew_iqr"] = float(
                    np.nanpercentile(sb, 75)
                    - np.nanpercentile(sb, 25)
                )
            else:
                rec["skew_med"] = np.nan
                rec["skew_iqr"] = np.nan
        rows.append(rec)
    cols = [
        "station", "n_freq", "n_ok", "frac_ok", "n_tip",
        "n_tip_ok", "snr_med", "pmin", "pmax",
    ]
    if include_skew:
        cols += ["skew_med", "skew_iqr"]
    return pd.DataFrame.from_records(rows, columns=cols)


def qc_flags(
    sites: Any,
    *,
    min_frac_ok: float = 0.6,
    min_snr_med: float = 2.0,
    max_skew_med: float = 6.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    tb = build_qc_table(
        sites,
        include_skew=True,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if tb.empty:
        return tb
    flags = []
    for _, r in tb.iterrows():
        f = []
        if float(r["frac_ok"]) < float(min_frac_ok):
            f.append("low_coverage")
        if np.isfinite(r["snr_med"]) and r["snr_med"] < min_snr_med:
            f.append("low_snr")
        if np.isfinite(r.get("skew_med", np.nan)) and \
           r["skew_med"] > max_skew_med:
            f.append("high_skew")
        flags.append(",".join(f))
    out = tb.copy()
    out["flags"] = flags
    return out


# -------------------- confidence profile (Kouadio et al. 2024 Fig. 3) --- #

def plot_confidence_profile(
    sites: Any,
    *,
    ci_hi: float = 0.95,
    ci_lo: float = 0.50,
    spacing_m: float = 200.0,
    figsize: Tuple[float, float] = (9.0, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Profile confidence-index (CI) scatter plot along the survey line.

    Reproduces the Fig. 3 style from Kouadio et al. (2024): one dot per
    station coloured green (CI ≥ ``ci_hi``), pink (``ci_lo`` ≤ CI < ``ci_hi``),
    or red (CI < ``ci_lo``), with dashed threshold lines.

    CI is the fraction of frequencies with a valid (finite) Z tensor at
    each station — ``n_ok / n_freq`` from :func:`build_qc_table`.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
        Input sites.
    ci_hi : float
        Upper CI threshold (default 0.95 — "safe", green).
    ci_lo : float
        Lower CI threshold (default 0.50 — "recoverable", pink).
    spacing_m : float
        Fallback station spacing [m] used when no coordinate metadata is
        available on the EDI objects.
    figsize : tuple
        Figure size when a new figure is created.
    recursive, on_dup, strict, verbose
        Passed to :func:`ensure_sites`.
    ax : matplotlib.axes.Axes or None
        Axes to draw on; created if *None*.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    tb = build_qc_table(
        S, include_skew=False, recursive=False,
        on_dup=on_dup, strict=False, verbose=verbose,
    )

    items     = list(_iter_items(S))
    positions = _station_positions(items, spacing_m)
    names     = [_name(ed, i) for i, ed in enumerate(items)]

    xs, ys, cs = [], [], []
    for k, name in enumerate(names):
        row = tb[tb["station"] == name]
        if row.empty:
            continue
        ci = float(row["frac_ok"].values[0])
        xs.append(float(positions[k]))
        ys.append(ci)
        if ci >= ci_hi:
            cs.append("#2ca02c")      # green
        elif ci >= ci_lo:
            cs.append("#ff9999")      # pink
        else:
            cs.append("#d62728")      # red

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if xs:
        ax.scatter(xs, ys, c=cs, s=60, zorder=3, edgecolors="none")

    ax.axhline(
        ci_hi, ls="--", color="#2ca02c", lw=1.2,
        label=f"CI = {ci_hi:.2f} (safe)",
    )
    ax.axhline(
        ci_lo, ls="--", color="#ff9999", lw=1.2,
        label=f"CI = {ci_lo:.2f} (recoverable)",
    )
    ax.set_ylim(-0.05, 1.10)
    ax.set_xlabel("Distance along profile (m)")
    ax.set_ylabel("Confidence index (CI)")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.4)
    return ax


# ----------------------------- coverage plot ----------------------------- #

def plot_coverage_psection(
    sites: Any,
    *,
    metric: str = "presence",  # presence|snr|offdiag
    alpha_by: str = "none",    # none|snr
    figsize: Tuple[float, float] = (9.0, 4.8),
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
    sts: List[str] = []
    Ys: List[np.ndarray] = []
    Ms: List[np.ndarray] = []
    As: List[np.ndarray] = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed, return_errors=True)
        if Z is None:
            continue
        sts.append(st)
        lp = np.log10(np.maximum(1.0 / fr, 1e-9))
        Ys.append(lp)
        if metric == "offdiag":
            M = _offdiag_logmag(z)
        elif metric == "snr":
            ze = Z[3] if isinstance(Z, tuple) else None
            ze = getattr(ed.Z, "z_err", None) if ze is None else ze
            M = _snr_rows(z, ze)
        else:
            M = _row_ok_z(z).astype(float)
        Ms.append(M.astype(float))
        if alpha_by == "snr":
            ze = getattr(ed.Z, "z_err", None)
            A = _snr_rows(z, ze)
        else:
            A = np.ones_like(M, dtype=float)
        As.append(A)
    if not sts:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    yall = np.unique(np.concatenate(Ys))
    nx = len(sts)
    Zm = np.zeros((yall.size, nx, 4), dtype=float)
    v = []
    a = []
    for j, (lp, m, al) in enumerate(zip(Ys, Ms, As)):
        i = np.searchsorted(yall, lp)
        i = np.clip(i, 0, yall.size - 1)
        vv = np.nan_to_num(m, nan=np.nan)
        aa = np.nan_to_num(al, nan=0.0)
        v.append(vv); a.append(aa)
        # map metric to color
        if metric == "presence":
            col = (0.20, 0.60, 0.20)
            Zm[i, j, :3] = col
            Zm[i, j, 3] = aa
        else:
            # use viridis for metric
            pass
    if metric != "presence":
        V = np.concatenate(v)
        V = V[np.isfinite(V)]
        v0 = np.nanpercentile(V, 5) if V.size else 0.0
        v1 = np.nanpercentile(V, 95) if V.size else 1.0
        for j, (lp, m, al) in enumerate(zip(Ys, Ms, As)):
            i = np.searchsorted(yall, lp)
            i = np.clip(i, 0, yall.size - 1)
            sc = (m - v0) / (v1 - v0 + 1e-12)
            sc = np.clip(sc, 0.0, 1.0)
            rgb = plt.cm.viridis(sc)
            Zm[i, j, :3] = rgb[:, :3]
            Zm[i, j, 3] = al
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.imshow(
        Zm,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(nx))
    ax.set_xticklabels(sts, rotation=90)
    yt, yl = _y_ticks(yall, 8)
    ax.set_yticks(yt)
    ax.set_yticklabels(yl)
    return ax


# ----------------------------- SNR histogram ----------------------------- #

def plot_snr_hist(
    sites: Any,
    *,
    bins: int = 40,
    figsize: Tuple[float, float] = (7.2, 3.6),
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
    vals: List[float] = []
    for _, ed in enumerate(_iter_items(S)):
        Z, z, fr = _get_z_block(ed, return_errors=True)
        if Z is None:
            continue
        if isinstance(Z, tuple) and len(Z) == 4:
            _, z, fr, ze = Z
        else:
            ze = getattr(ed.Z, "z_err", None)
        snr = _snr_rows(z, ze)
        vals.extend(list(snr))
    v = np.array(vals, dtype=float)
    v = v[np.isfinite(v)]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.hist(v, bins=int(max(8, bins)))
    ax.set_xlabel("row SNR (|Z|/σ)")
    ax.set_ylabel("count")
    return ax


# ----------------------------- quicklook -------------------------------- #

def plot_qc_quicklook(
    sites: Any,
    *,
    figsize: Tuple[float, float] = (10.0, 8.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    plot_coverage_psection(
        sites,
        metric="presence",
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        ax=ax1,
    )
    plot_coverage_psection(
        sites,
        metric="snr",
        alpha_by="snr",
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        ax=ax2,
    )
    plot_snr_hist(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        ax=ax3,
    )
    return fig

# ---------------------- z-block helper (errors) ------------------------- #

def _zblk(ed: Any, need_err: bool = False):
    try:
        return _get_z_block(ed, return_errors=need_err)
    except TypeError:
        try:
            return _get_z_block(ed, with_errors=need_err)
        except TypeError:
            return _get_z_block(ed)


# ----------------------- rho_a + error propagation ---------------------- #

def _rhoa_xy_yx(z: np.ndarray, fr: np.ndarray) -> Tuple[
    np.ndarray, np.ndarray
]:
    c = 0.2 / (fr + 1e-24)
    rxy = c * (np.abs(z[:, 0, 1]) ** 2)
    ryx = c * (np.abs(z[:, 1, 0]) ** 2)
    return rxy, ryx


def _rhoa_ci(
    z: np.ndarray,
    ze: Optional[np.ndarray],
    fr: np.ndarray,
    *,
    comp: str = "xy",        # xy|yx
    pcts: Tuple[float, ...] = (10.0, 50.0, 90.0),
    n_draws: int = 200,
    seed: Optional[int] = 0,
) -> np.ndarray:
    a, b = (0, 1) if comp == "xy" else (1, 0)
    zz = z[:, a, b]
    if ze is None:
        c = 0.2 / (fr + 1e-24)
        m = c * (np.abs(zz) ** 2)
        P = [np.zeros_like(m) for _ in pcts]
        return np.vstack([m] + P).T
    ee = ze[:, a, b]
    g = np.isfinite(zz) & np.isfinite(ee)
    if not np.any(g):
        m = np.full(z.shape[0], np.nan, dtype=float)
        P = [np.full_like(m, np.nan) for _ in pcts]
        return np.vstack([m] + P).T
    rng = np.random.default_rng(seed)
    nf = zz.size
    n = int(max(16, n_draws))
    # complex Gaussian, σ equals |ze|
    E = (rng.standard_normal((n, nf))
         + 1j * rng.standard_normal((n, nf))) / np.sqrt(2.0)
    E = E * ee[None, :]
    Zs = zz[None, :] + E
    c = 0.2 / (fr + 1e-24)
    R = c[None, :] * (np.abs(Zs) ** 2)
    M = np.nanmedian(R, axis=0)
    Q = [np.nanpercentile(R, q, axis=0) for q in pcts]
    return np.vstack([M] + Q).T


def _shade_band(
    ax: plt.Axes,
    x: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    *,
    alpha: float = 0.25,
    color: str = "C0",
):
    ax.fill_between(x, lo, hi, alpha=alpha, color=color)


# ----------------------- 19) Consistency fan chart ---------------------- #

def plot_consistency_fan(
    sites: Any,
    *,
    station: Optional[str] = None,
    other: Any | None = None,      # optional comparison Sites
    comps: Tuple[str, str] = ("xy", "yx"),
    pcts: Tuple[float, float, float] = (10.0, 50.0, 90.0),
    n_draws: int = 200,
    figsize: Tuple[float, float] = (8.6, 4.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    ed_map = {}
    for i, ed in enumerate(_iter_items(S)):
        ed_map[_name(ed, i)] = ed
    if not ed_map:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no sites", ha="center", va="center")
        return ax
    if station is None:
        station = sorted(ed_map.keys())[0]
    ed = ed_map.get(station, None)
    if ed is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "station not found",
                ha="center", va="center")
        return ax
    out = _zblk(ed, need_err=True)
    if len(out) == 4:
        Z, z, fr, ze = out
    else:
        Z, z, fr = out[:3]
        ze = None
    if Z is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no Z", ha="center", va="center")
        return ax
    per = 1.0 / fr
    x = per
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.set_xscale("log")
    cols = {"xy": "C0", "yx": "C2"}
    for c in comps:
        CI = _rhoa_ci(
            z, ze, fr, comp=c, pcts=pcts, n_draws=n_draws
        )
        med = CI[:, 0]
        lo = np.minimum(CI[:, 1], CI[:, 2])
        hi = np.maximum(CI[:, 1], CI[:, 2])
        _shade_band(ax, x, lo, hi, color=cols[c], alpha=0.20)
        ax.plot(x, med, "-", lw=2.0, color=cols[c],
                label=f"ρa_{c}")
    if other is not None:
        So = ensure_sites(other, recursive=False, strict=False)
        # overlay only medians (dashed)
        for i2, edo in enumerate(_iter_items(So)):
            if _name(edo, i2) != station:
                continue
            Z2, z2, fr2 = _zblk(edo)[:3]
            if Z2 is None:
                break
            x2 = 1.0 / fr2
            rxy2, ryx2 = _rhoa_xy_yx(z2, fr2)
            if "xy" in comps:
                ax.plot(x2, rxy2, "--", lw=1.2, color=cols["xy"],
                        label="after xy")
            if "yx" in comps:
                ax.plot(x2, ryx2, "--", lw=1.2, color=cols["yx"],
                        label="after yx")
            break
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("ρa (Ω·m)")
    ax.grid(True, alpha=0.25, which="both")
    ax.set_title(str(station))
    ax.legend(ncol=2, fontsize=8)
    return ax


# ---------------------- 20) XY–YX crossover map ------------------------- #

def plot_xyyx_crossover_map(
    sites: Any,
    *,
    figsize: Tuple[float, float] = (9.0, 4.6),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    Xs, Ys, labels = [], [], []
    sts = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i); sts.append(st)
        Z, z, fr = _zblk(ed)[:3]
        if Z is None:
            continue
        rxy, ryx = _rhoa_xy_yx(z, fr)
        d = rxy - ryx
        p = 1.0 / fr
        if d.size < 2:
            continue
        s = np.sign(d)
        sc = s[:-1] * s[1:] <= 0.0
        idx = np.where(sc)[0]
        if idx.size == 0:
            continue
        lp = np.log10(np.maximum(p, 1e-9))
        for k in idx:
            w1 = np.abs(d[k]); w2 = np.abs(d[k + 1])
            t = w1 / (w1 + w2 + 1e-24)
            y = (1.0 - t) * lp[k] + t * lp[k + 1]
            Xs.append(i)
            Ys.append(y)
            labels.append(st)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if not Xs:
        ax.text(0.5, 0.5, "no crossovers",
                ha="center", va="center")
        return ax
    ax.scatter(Xs, Ys, s=16, c="crimson", alpha=0.8)
    ax.set_xlabel("Station")
    ax.set_ylabel("LogPeriod (s)")
    ax.set_xticks(np.arange(len(sts)))
    ax.set_xticklabels(sts, rotation=90)
    # y ticks from Ys
    yall = np.array(Ys, dtype=float)
    yt, yl = _y_ticks(yall, 8)
    ax.set_yticks(yt)
    lo, hi = float(np.nanmin(yall)), float(np.nanmax(yall))
    ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
    return ax


# ---------------------- 21) Noise cone overlay -------------------------- #

def overlay_noise_cone(
    ax: plt.Axes,
    period: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    *,
    color: str = "0.6",
    alpha: float = 0.18,
):
    x = period
    _shade_band(ax, x, lo, hi, color=color, alpha=alpha)


# ---------------------- 22) Spectral hole finder ------------------------ #

def overlay_spectral_holes(
    ax: plt.Axes,
    sites: Any,
    *,
    thresh_dec: float = 0.30,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    # assumes x = station index, y = log10(period)
    xmap = {}
    for i, ed in enumerate(_iter_items(S)):
        xmap[_name(ed, i)] = i
    for i, ed in enumerate(_iter_items(S)):
        Z, z, fr = _zblk(ed)[:3]
        if Z is None:
            continue
        p = 1.0 / fr
        lp = np.sort(np.log10(np.maximum(p, 1e-9)))
        if lp.size < 2:
            continue
        d = np.diff(lp)
        holes = np.where(d > float(thresh_dec))[0]
        for h in holes:
            y0, y1 = lp[h], lp[h + 1]
            r = _Rect(
                (i - 0.45, y0),
                0.90,
                (y1 - y0),
                facecolor=(0.8, 0.2, 0.2, 0.08),
                edgecolor="none",
                zorder=0.0,
            )
            ax.add_patch(r)
