"""
Polar uncertainty diagnostics for CSAMT impedance data.

Adapts the polar-based visualization framework from:
  kouadio2025 : Kouadio K.L. (2025), "k-diagram: Rethinking Forecasting
                Uncertainty via Polar-based Visualization",
                J. Open Source Softw. 10(116), 8661.
                DOI: 10.21105/joss.08661

The three core diagnostics from k-diagram are translated to CSAMT:

Coverage evaluation (kouadio2025 eq. 1):
    c_j = 1(L_j ≤ ρ_a,obs,j ≤ U_j)
    Binary coverage of per-frequency observed apparent resistivity within
    the predicted quantile interval [L_j, U_j].

Frequency-width drift (analogous to Horizon Drift, kouadio2025 Fig 2c):
    w̄_b = mean_{j ∈ band_b} (U_j − L_j) / ρ_a,obs,j × 100   [%]
    Mean relative interval width per frequency band — tracks how
    prediction uncertainty grows with period (a proxy for depth).

Relative-error polar histogram (analogous to polar violin, Fig 2b):
    ε_j = (ρ_a,pred,j − ρ_a,obs,j) / ρ_a,obs,j × 100   [%]
    Rose diagram of residuals binned by frequency decade.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ._core import (
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
    hide_polar_radius_labels,
)

__all__ = [
    "COVERAGE_THRESH",
    "coverage_score",
    "rho_coverage",
    "rho_error_stats",
    "coverage_table",
    "plot_polar_coverage",
    "plot_polar_errors",
    "plot_width_drift",
]

COVERAGE_THRESH: float = 0.9            # default nominal probability level


# ─────────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────────

def _unwrap(ed: Any) -> Any:
    edi = getattr(ed, "edi", None)
    if edi is not None and hasattr(edi, "Z"):
        return edi
    return ed


def _rho_a_from_z(
    z_block: np.ndarray,
    freqs: np.ndarray,
    comp: str,
) -> np.ndarray:
    """Cagniard ρ_a = 0.2 |Z_pq|² / f, for Z in practical units (mV/km per nT).

    Matches the convention used by :mod:`pycsamt.emtools.csumt`'s
    ``_rho_a_det`` — this used to divide by ``ωμ₀`` instead, which
    assumes Z is in SI ohms and is wrong by a ~10^5-10^6 factor for the
    practical-unit Z stored in EDI files.
    """
    if comp == "xy":
        Z = z_block[:, 0, 1]
    elif comp == "yx":
        Z = z_block[:, 1, 0]
    else:
        raise ValueError(f"rho_comp must be 'xy' or 'yx', got {comp!r}")
    return 0.2 * np.abs(Z) ** 2 / np.maximum(freqs, 1e-24)


def _fmt_hz(f: float) -> str:
    """Format a frequency compactly, avoiding scientific notation."""
    if f >= 1000.0:
        return f"{f / 1000.0:.3g} kHz"
    return f"{f:.3g} Hz"


def _resolve_bounds(
    bounds: dict | np.ndarray | float,
    station: str,
    n: int,
) -> np.ndarray | None:
    """Return per-frequency bound array for *station*, or None if missing."""
    if isinstance(bounds, dict):
        arr = bounds.get(station)
        if arr is None:
            return None
        return np.asarray(arr, dtype=float)
    arr = np.asarray(bounds, dtype=float)
    return np.full(n, float(arr)) if arr.ndim == 0 else arr


def _rho_dict_from_sites(sites_obj: Any, comp: str) -> dict[str, np.ndarray]:
    """Build {station: rho_a_array} from a sites-like object."""
    result: dict[str, np.ndarray] = {}
    try:
        for i, ed in enumerate(_iter_items(sites_obj)):
            ed = _unwrap(ed)
            _, z_block, freqs = _get_z_block(ed)
            if z_block is None or freqs is None or freqs.size == 0:
                continue
            station = _name(ed, i)
            result[station] = _rho_a_from_z(z_block, freqs, comp)
    except (TypeError, AttributeError):
        pass
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Pure-math helper — no sites dependency
# ─────────────────────────────────────────────────────────────────────────────

def coverage_score(
    y_true: np.ndarray | float,
    y_lo: np.ndarray | float,
    y_hi: np.ndarray | float,
) -> float:
    """
    Empirical coverage fraction of a prediction interval (kouadio2025 eq. 1).

    Parameters
    ----------
    y_true : array-like
        Observed values.
    y_lo, y_hi : array-like
        Lower and upper bounds of the predicted interval.

    Returns
    -------
    cov : float
        Fraction of observations that fall inside [y_lo, y_hi] ∈ [0, 1].
    """
    y  = np.asarray(y_true, dtype=float).ravel()
    lo = np.asarray(y_lo,   dtype=float).ravel()
    hi = np.asarray(y_hi,   dtype=float).ravel()
    return float(np.mean((y >= lo) & (y <= hi)))


# ─────────────────────────────────────────────────────────────────────────────
# Sites-based: per-frequency DataFrames
# ─────────────────────────────────────────────────────────────────────────────

def rho_coverage(
    sites: Any,
    q_lo: dict[str, np.ndarray] | np.ndarray | float,
    q_hi: dict[str, np.ndarray] | np.ndarray | float,
    *,
    rho_comp: str = "xy",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Per-frequency coverage of observed ρ_a within predicted quantile bounds.

    For each site and frequency, checks whether the Cagniard apparent
    resistivity extracted from the observed Z tensor falls inside the
    predicted interval [L_j, U_j] (kouadio2025 eq. 1):
        c_j = 1(L_j ≤ ρ_a,obs,j ≤ U_j)

    Parameters
    ----------
    sites : Sites | list
        Observed CSAMT sites.
    q_lo, q_hi : dict {station: array} or array or scalar
        Lower / upper quantile bounds aligned with each site's frequency
        array.  When a dict, keys must match station names; sites without
        a key are skipped.  A scalar broadcasts to every frequency of
        every site.
    rho_comp : {"xy", "yx"}
        Impedance component used to derive ρ_a (default ``"xy"``).
    recursive, on_dup, strict, verbose
        Forwarded to :func:`ensure_sites`.

    Returns
    -------
    pd.DataFrame
        Columns: station, freq_hz, period_s, rho_obs, q_lo, q_hi,
        covered, width_pct.
        ``covered`` is bool; ``width_pct`` = 100 × (q_hi−q_lo) / ρ_obs.
    """
    _COLS = ["station", "freq_hz", "period_s", "rho_obs",
             "q_lo", "q_hi", "covered", "width_pct"]

    sites = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                         strict=strict, verbose=verbose)
    rows: list[dict] = []

    for i, ed in enumerate(_iter_items(sites)):
        ed = _unwrap(ed)
        station = _name(ed, i)
        _, z_block, freqs = _get_z_block(ed)
        if z_block is None or freqs is None or freqs.size == 0:
            continue

        rho_obs = _rho_a_from_z(z_block, freqs, rho_comp)
        lo = _resolve_bounds(q_lo, station, freqs.size)
        hi = _resolve_bounds(q_hi, station, freqs.size)
        if lo is None or hi is None:
            continue

        for j in range(freqs.size):
            f   = float(freqs[j])
            r_o = float(rho_obs[j])
            l_j = float(lo[j])
            h_j = float(hi[j])
            w   = 100.0 * (h_j - l_j) / r_o if r_o > 0 else np.nan
            rows.append({
                "station":   station,
                "freq_hz":   f,
                "period_s":  1.0 / f if f > 0 else np.nan,
                "rho_obs":   r_o,
                "q_lo":      l_j,
                "q_hi":      h_j,
                "covered":   bool(l_j <= r_o <= h_j),
                "width_pct": w,
            })

    if not rows:
        return pd.DataFrame(columns=_COLS)
    return pd.DataFrame(rows, columns=_COLS)


def rho_error_stats(
    sites: Any,
    model_rho: dict[str, np.ndarray] | Any,
    *,
    rho_comp: str = "xy",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Per-frequency relative error between observed and predicted ρ_a.

    Computes ε_j = (ρ_a,pred,j − ρ_a,obs,j) / ρ_a,obs,j × 100 % for
    each station and frequency, analogous to the error distribution
    visualised in the k-diagram polar violin (kouadio2025 Fig 2b).

    Parameters
    ----------
    sites : Sites | list
        Observed CSAMT sites.
    model_rho : dict {station: array} or Sites-like
        Predicted apparent resistivity.  Either a dict mapping station
        names to 1-D arrays (same length as the corresponding site's
        frequency array), or a Sites-like container from which ρ_a is
        extracted with the same ``rho_comp`` setting.
    rho_comp : {"xy", "yx"}
        Impedance component.
    recursive, on_dup, strict, verbose
        Forwarded to :func:`ensure_sites`.

    Returns
    -------
    pd.DataFrame
        Columns: station, freq_hz, period_s, rho_obs, rho_pred,
        rel_err_pct, abs_err_pct.
    """
    _COLS = ["station", "freq_hz", "period_s", "rho_obs",
             "rho_pred", "rel_err_pct", "abs_err_pct"]

    sites = ensure_sites(sites, recursive=recursive, on_dup=on_dup,
                         strict=strict, verbose=verbose)

    if isinstance(model_rho, dict):
        pred_dict = {k: np.asarray(v, dtype=float) for k, v in model_rho.items()}
    else:
        pred_dict = _rho_dict_from_sites(model_rho, rho_comp)

    rows: list[dict] = []
    for i, ed in enumerate(_iter_items(sites)):
        ed = _unwrap(ed)
        station = _name(ed, i)
        _, z_block, freqs = _get_z_block(ed)
        if z_block is None or freqs is None or freqs.size == 0:
            continue

        rho_obs  = _rho_a_from_z(z_block, freqs, rho_comp)
        rho_pred = pred_dict.get(station)
        if rho_pred is None:
            continue

        for j in range(freqs.size):
            f   = float(freqs[j])
            r_o = float(rho_obs[j])
            r_p = float(rho_pred[j]) if j < len(rho_pred) else np.nan
            if r_o > 0 and np.isfinite(r_p):
                rel_err = 100.0 * (r_p - r_o) / r_o
                abs_err = abs(rel_err)
            else:
                rel_err = abs_err = np.nan
            rows.append({
                "station":     station,
                "freq_hz":     f,
                "period_s":    1.0 / f if f > 0 else np.nan,
                "rho_obs":     r_o,
                "rho_pred":    r_p,
                "rel_err_pct": rel_err,
                "abs_err_pct": abs_err,
            })

    if not rows:
        return pd.DataFrame(columns=_COLS)
    return pd.DataFrame(rows, columns=_COLS)


def coverage_table(
    sites: Any,
    q_lo: dict[str, np.ndarray] | np.ndarray | float,
    q_hi: dict[str, np.ndarray] | np.ndarray | float,
    *,
    rho_comp: str = "xy",
    nominal: float = COVERAGE_THRESH,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Per-station coverage summary.

    Parameters
    ----------
    nominal : float
        Target coverage probability (e.g. 0.9 for a 90 % interval).
        ``calibrated_flag`` is True when the empirical coverage ≥ nominal.

    Returns
    -------
    pd.DataFrame
        Columns: station, n_freq, empirical_cov, mean_width_pct,
        calibrated_flag.
    """
    _COLS = ["station", "n_freq", "empirical_cov",
             "mean_width_pct", "calibrated_flag"]

    detail = rho_coverage(sites, q_lo, q_hi, rho_comp=rho_comp,
                          recursive=recursive, on_dup=on_dup,
                          strict=strict, verbose=verbose)
    if detail.empty:
        return pd.DataFrame(columns=_COLS)

    rows: list[dict] = []
    for station, grp in detail.groupby("station", sort=False):
        emp = float(grp["covered"].mean())
        w   = float(grp["width_pct"].dropna().mean()) \
              if grp["width_pct"].notna().any() else np.nan
        rows.append({
            "station":         station,
            "n_freq":          len(grp),
            "empirical_cov":   emp,
            "mean_width_pct":  w,
            "calibrated_flag": bool(emp >= nominal),
        })
    return pd.DataFrame(rows, columns=_COLS)


# ─────────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────────

def plot_polar_coverage(
    sites: Any,
    q_lo: dict | np.ndarray | float,
    q_hi: dict | np.ndarray | float,
    *,
    rho_comp: str = "xy",
    log_radius: bool = True,
    n_freq_ticks: int = 8,
    figsize: tuple = (7, 7),
    title: str = "Coverage evaluation",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax=None,
):
    """
    Polar coverage plot: angle ∝ log₁₀(f), radius ∝ ρ_a,obs.

    Green markers = observed ρ_a within predicted interval (covered);
    red = outside.  Thin radial segments show each [q_lo, q_hi] range.
    Every station shares (almost) the same frequency grid, so each
    angular position is really one frequency shared by every station —
    the angle axis is labelled with that frequency directly (rather
    than the otherwise-meaningless default degree ticks) so a reader
    can tell which part of the band a cluster of misses falls in.

    Parameters
    ----------
    n_freq_ticks : int, default 8
        Number of evenly (log-)spaced frequency labels drawn around the
        ring. Set to 0 to fall back to the default degree ticks.

    Returns
    -------
    ax : matplotlib.axes.Axes (polar projection)
    """
    import matplotlib.pyplot as plt

    df = rho_coverage(sites, q_lo, q_hi, rho_comp=rho_comp,
                      recursive=recursive, on_dup=on_dup,
                      strict=strict, verbose=verbose)

    if ax is None:
        _, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=figsize)
    hide_polar_radius_labels(ax)

    if df.empty:
        ax.set_title(title + " (no data)")
        return ax

    f_vals = df["freq_hz"].values
    f_min, f_max = f_vals.min(), f_vals.max()
    if f_max <= f_min:
        f_max = f_min + 1.0
    log_f_min, log_f_max = np.log10(f_min), np.log10(f_max)
    log_f = np.log10(f_vals)
    theta = 2.0 * np.pi * (log_f - log_f_min) / (log_f_max - log_f_min)

    def _r(v):
        return np.log10(np.maximum(v, 1e-30)) if log_radius else v

    r_obs = _r(df["rho_obs"].values)
    r_lo  = _r(df["q_lo"].values)
    r_hi  = _r(df["q_hi"].values)
    covered = df["covered"].values

    ax.scatter(theta[covered],  r_obs[covered],  c="green", s=18, alpha=0.8,
               zorder=3, label="covered")
    ax.scatter(theta[~covered], r_obs[~covered], c="red",   s=18, alpha=0.8,
               zorder=3, label="not covered")

    for t, lo, hi in zip(theta, r_lo, r_hi):
        ax.plot([t, t], [lo, hi], color="gray", lw=0.5, alpha=0.35, zorder=1)

    if n_freq_ticks > 0:
        tick_theta = np.linspace(0.0, 2.0 * np.pi, n_freq_ticks, endpoint=False)
        tick_freq = 10.0 ** (
            log_f_min + (tick_theta / (2.0 * np.pi)) * (log_f_max - log_f_min)
        )
        ax.set_xticks(tick_theta)
        ax.set_xticklabels([_fmt_hz(f) for f in tick_freq], fontsize=8)

    emp_cov = float(covered.mean())
    hide_polar_radius_labels(ax)
    ax.set_title(f"{title}\ncoverage = {emp_cov:.3f}", pad=15)
    ax.legend(loc="lower right", fontsize=7)
    return ax


def plot_polar_errors(
    sites: Any,
    model_rho: dict[str, np.ndarray] | Any,
    *,
    rho_comp: str = "xy",
    n_bins: int = 18,
    figsize: tuple = (7, 7),
    title: str = "Error distribution",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax=None,
):
    """
    Polar rose diagram of relative residuals (ε = (ρ_pred − ρ_obs)/ρ_obs × 100 %).

    Each angular sector spans one frequency decade.  Bar length = mean |ε|
    within that sector; red = over-prediction (mean ε > 0), blue = under.
    Analogous to the polar violin in kouadio2025 Fig 2b.

    Returns
    -------
    ax : matplotlib.axes.Axes (polar projection)
    """
    import matplotlib.pyplot as plt

    df = rho_error_stats(sites, model_rho, rho_comp=rho_comp,
                         recursive=recursive, on_dup=on_dup,
                         strict=strict, verbose=verbose)

    if ax is None:
        _, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=figsize)
    hide_polar_radius_labels(ax)

    if df.empty or df["rel_err_pct"].isna().all():
        ax.set_title(title + " (no data)")
        return ax

    f_vals = df["freq_hz"].values
    f_min, f_max = f_vals.min(), f_vals.max()
    if f_max <= f_min:
        f_max = f_min + 1.0

    log_f = np.log10(f_vals)
    bins  = np.linspace(np.log10(f_min), np.log10(f_max), n_bins + 1)
    theta = np.deg2rad(np.linspace(0, 360, n_bins, endpoint=False))
    bw    = 2.0 * np.pi / n_bins

    errs     = df["rel_err_pct"].values
    heights  = np.zeros(n_bins)
    signs    = np.zeros(n_bins)

    for b in range(n_bins):
        mask = (log_f >= bins[b]) & (log_f < bins[b + 1])
        if mask.any():
            signs[b]   = float(np.sign(np.nanmean(errs[mask])))
            heights[b] = float(np.nanmean(np.abs(errs[mask])))

    colors = ["#e74c3c" if s >= 0 else "#2980b9" for s in signs]
    ax.bar(theta, heights, width=bw * 0.85, color=colors, alpha=0.72, bottom=0)
    hide_polar_radius_labels(ax)
    ax.set_title(f"{title}\nred=over-pred, blue=under-pred", pad=15)
    return ax


def plot_width_drift(
    sites: Any,
    q_lo: dict | np.ndarray | float,
    q_hi: dict | np.ndarray | float,
    *,
    rho_comp: str = "xy",
    n_bands: int = 8,
    polar: bool = False,
    figsize: tuple = (8, 4),
    title: str = "Frequency-width drift",
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax=None,
):
    """
    Mean relative interval width per frequency band (horizon-drift analogue).

    Visualises how the predicted interval width — relative to observed ρ_a —
    changes across the frequency spectrum, a proxy for how model uncertainty
    grows with probing depth (kouadio2025 § Forecast Horizon Drift).

    Parameters
    ----------
    n_bands : int
        Number of frequency bands.
    polar : bool
        Radial bar chart (k-diagram style) when True; Cartesian otherwise.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    df = rho_coverage(sites, q_lo, q_hi, rho_comp=rho_comp,
                      recursive=recursive, on_dup=on_dup,
                      strict=strict, verbose=verbose)

    if ax is None:
        if polar:
            _, ax = plt.subplots(subplot_kw={"projection": "polar"},
                                 figsize=figsize)
        else:
            _, ax = plt.subplots(figsize=figsize)
    if polar:
        hide_polar_radius_labels(ax)

    if df.empty:
        ax.set_title(title + " (no data)")
        return ax

    log_f  = np.log10(df["freq_hz"].values)
    f_min, f_max = log_f.min(), log_f.max()
    if f_max <= f_min:
        f_max = f_min + 1.0

    bins      = np.linspace(f_min, f_max, n_bands + 1)
    centers   = 10.0 ** (0.5 * (bins[:-1] + bins[1:]))
    widths_pct = df["width_pct"].values
    mean_w    = np.zeros(n_bands)

    for b in range(n_bands):
        mask = (log_f >= bins[b]) & (log_f < bins[b + 1])
        if mask.any():
            mean_w[b] = float(np.nanmean(widths_pct[mask]))

    labels = [f"{c:.2g}" for c in centers]

    if polar:
        theta = np.linspace(0, 2.0 * np.pi, n_bands, endpoint=False)
        bw    = 2.0 * np.pi / n_bands
        ax.bar(theta, mean_w, width=bw * 0.85, alpha=0.78)
        ax.set_xticks(theta)
        ax.set_xticklabels(labels, fontsize=7)
        hide_polar_radius_labels(ax)
    else:
        x = np.arange(n_bands)
        ax.bar(x, mean_w, color="#3498db", alpha=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_xlabel("Frequency band centre [Hz]")
        ax.set_ylabel("Mean interval width [% of ρ_a]")
        ax.grid(True, linestyle=":", alpha=0.5)

    ax.set_title(title)
    return ax
