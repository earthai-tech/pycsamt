"""
3-D axial anisotropy analysis for CSAMT impedance tensor data.

Implements the anisotropy metrics derived in:
  wang2017 : Wang & Tan (2017), "Research on the forward modeling of
             controlled-source audio-frequency magnetotellurics in
             three-dimensional axial anisotropic media",
             J. Appl. Geophys. 146, 27–36.

For axial anisotropy the conductivity tensor is diagonal:
    σ* = diag(σ_xx, σ_yy, σ_zz)

The CSAMT impedance tensor Z = [[Z_xx, Z_xy], [Z_yx, Z_yy]] yields
two independent Cagniard apparent resistivities (eqs 17–18):
    ρ_xy = (1/ωμ₀)|Z_xy|²   (equatorial configuration, sensitive to σ_xx)
    ρ_yx = (1/ωμ₀)|Z_yx|²   (axial configuration, sensitive to σ_yy, σ_zz)

Their ratio Λ = ρ_xy/ρ_yx, expressed in log₁₀ units, is the primary
anisotropy indicator (Λ = 0 for isotropic media).

The Swift (1967) skew S = |Z_xx − Z_yy| / |Z_xy + Z_yx| and the
corresponding rotation angle provide additional dimensionality and
strike information (non-zero for anisotropic / 3-D structures).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..api.station import PYCSAMT_STATION_RENDERING
from ._core import (
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
)

__all__ = [
    "ANISO_RATIO_THRESH",
    "SWIFT_SKEW_THRESH",
    "analyze_anisotropy",
    "anisotropy_table",
    "plot_anisotropy",
]

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

ANISO_RATIO_THRESH: float = 0.1  # |log10(Λ)| above which flag is raised
SWIFT_SKEW_THRESH: float = 0.2  # Swift skew above which 3-D is suspected

_DETAIL_COLS = [
    "station",
    "freq_hz",
    "period_s",
    "rho_xy_ohmm",
    "rho_yx_ohmm",
    "phi_xy_deg",
    "phi_yx_deg",
    "ratio_log10",
    "phase_diff_deg",
    "swift_skew",
    "strike_deg",
]
_TABLE_COLS = [
    "station",
    "n_freq",
    "mean_ratio_log10",
    "max_abs_ratio_log10",
    "mean_phase_diff_deg",
    "mean_swift_skew",
    "median_strike_deg",
    "anisotropy_flag",
]


# ─────────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────────


def _unwrap(ed: Any) -> Any:
    edi = getattr(ed, "edi", None)
    if edi is not None and hasattr(edi, "Z"):
        return edi
    return ed


def _rho_and_phase(
    z_block: np.ndarray,
    freq: np.ndarray,
) -> tuple:
    """
    Cagniard ρ_xy, ρ_yx, φ_xy, φ_yx from Z tensor (eqs 17–18, wang2017).

    ρ_pq = 0.2|Z_pq|²/f,  φ_pq = arctan(Im(Z_pq)/Re(Z_pq)), for Z in
    practical units (mV/km per nT) — matches the convention used by
    :mod:`pycsamt.emtools.csumt`'s ``_rho_a_det``. This used to divide
    by ``ωμ₀`` instead, which assumes Z is in SI ohms and is wrong by a
    ~10^5-10^6 factor for the practical-unit Z stored in EDI files;
    ``ratio_log10 = log10(ρ_xy/ρ_yx)`` was unaffected since the missing
    factor cancels in the ratio, but the absolute ``rho_xy_ohmm`` /
    ``rho_yx_ohmm`` columns were wrong.
    """
    Zxy = z_block[:, 0, 1]
    Zyx = z_block[:, 1, 0]

    rho_xy = 0.2 * np.abs(Zxy) ** 2 / np.maximum(freq, 1e-24)
    rho_yx = 0.2 * np.abs(Zyx) ** 2 / np.maximum(freq, 1e-24)
    phi_xy = np.angle(Zxy, deg=True)
    phi_yx = np.angle(Zyx, deg=True)
    return rho_xy, rho_yx, phi_xy, phi_yx


def _swift_skew(z_block: np.ndarray) -> np.ndarray:
    """
    Swift (1967) skew: S = |Z_xx − Z_yy| / |Z_xy + Z_yx|.

    S = 0 for a 1-D / 2-D isotropic earth; S > 0 indicates 3-D structure
    or anisotropy.
    """
    Zxx = z_block[:, 0, 0]
    Zxy = z_block[:, 0, 1]
    Zyx = z_block[:, 1, 0]
    Zyy = z_block[:, 1, 1]
    denom = np.abs(Zxy + Zyx)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denom < 1e-30, 0.0, np.abs(Zxx - Zyy) / denom)


def _swift_strike(z_block: np.ndarray) -> np.ndarray:
    """
    Per-frequency Swift strike angle [degrees] from Z tensor rotation.

    Minimises the diagonal elements of the rotated Z tensor:
        tan(2θ) = 2 Re[(Z_xy + Z_yx)(Z_xx − Z_yy)*]
                  / (|Z_xy + Z_yx|² − |Z_xx − Z_yy|²)
    """
    Zxx = z_block[:, 0, 0]
    Zxy = z_block[:, 0, 1]
    Zyx = z_block[:, 1, 0]
    Zyy = z_block[:, 1, 1]

    alpha = Zxy + Zyx  # sum of off-diagonal
    delta = Zxx - Zyy  # difference of diagonal

    num = 2.0 * np.real(alpha * np.conj(delta))
    denom = np.abs(alpha) ** 2 - np.abs(delta) ** 2

    theta = 0.5 * np.arctan2(num, denom)
    return np.rad2deg(theta)


# ─────────────────────────────────────────────────────────────────────────────
# Public: sites-based
# ─────────────────────────────────────────────────────────────────────────────


def analyze_anisotropy(
    sites: Any,
    *,
    ratio_threshold: float = ANISO_RATIO_THRESH,
    skew_threshold: float = SWIFT_SKEW_THRESH,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Per-frequency anisotropy metrics for a set of CSAMT sites.

    Computes the two Cagniard apparent resistivities ρ_xy and ρ_yx
    (wang2017 eqs 17–18), their log-ratio Λ = log₁₀(ρ_xy/ρ_yx), the
    phase difference, and the Swift skew from the full Z tensor.

    Parameters
    ----------
    sites : Sites | list
        EDI-like objects or a ``Sites`` container.
    ratio_threshold : float
        |log₁₀(Λ)| above which the anisotropy flag is raised (default 0.1).
    skew_threshold : float
        Swift skew above which 3-D / anisotropy is suspected (default 0.2).
    recursive, on_dup, strict, verbose
        Forwarded to :func:`ensure_sites`.

    Returns
    -------
    pd.DataFrame
        Columns: station, freq_hz, period_s, rho_xy_ohmm, rho_yx_ohmm,
        phi_xy_deg, phi_yx_deg, ratio_log10, phase_diff_deg,
        swift_skew, strike_deg.

    Notes
    -----
    ``ratio_log10 = 0`` and ``swift_skew = 0`` indicate a perfectly
    isotropic 1-D earth.  Non-zero diagonal Z elements (contributing to
    swift_skew > 0) suggest 3-D structure or electrical anisotropy
    (wang2017 §5.3).
    """
    sites = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    rows: list[dict] = []

    for i, ed in enumerate(_iter_items(sites)):
        ed = _unwrap(ed)
        station = _name(ed, i)

        _, z_block, freqs = _get_z_block(ed)
        if z_block is None or freqs is None or freqs.size == 0:
            continue

        rho_xy, rho_yx, phi_xy, phi_yx = _rho_and_phase(z_block, freqs)
        skew = _swift_skew(z_block)
        strike = _swift_strike(z_block)

        for j in range(freqs.size):
            f = float(freqs[j])
            rxy = float(rho_xy[j])
            ryx = float(rho_yx[j])

            if rxy > 0 and ryx > 0:
                ratio = float(np.log10(rxy / ryx))
            else:
                ratio = np.nan

            rows.append(
                {
                    "station": station,
                    "freq_hz": f,
                    "period_s": 1.0 / f if f > 0 else np.nan,
                    "rho_xy_ohmm": rxy,
                    "rho_yx_ohmm": ryx,
                    "phi_xy_deg": float(phi_xy[j]),
                    "phi_yx_deg": float(phi_yx[j]),
                    "ratio_log10": ratio,
                    "phase_diff_deg": float(phi_xy[j] - phi_yx[j]),
                    "swift_skew": float(skew[j]),
                    "strike_deg": float(strike[j]),
                }
            )

    if not rows:
        return pd.DataFrame(columns=_DETAIL_COLS)

    return pd.DataFrame(rows, columns=_DETAIL_COLS)


def anisotropy_table(
    sites: Any,
    *,
    ratio_threshold: float = ANISO_RATIO_THRESH,
    skew_threshold: float = SWIFT_SKEW_THRESH,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Per-station summary of anisotropy metrics.

    Parameters
    ----------
    sites : Sites | list
    ratio_threshold : float
        |log₁₀(Λ)| threshold for anisotropy flag (default 0.1).
    skew_threshold : float
        Swift skew threshold for anisotropy flag (default 0.2).

    Returns
    -------
    pd.DataFrame
        Columns: station, n_freq, mean_ratio_log10, max_abs_ratio_log10,
        mean_phase_diff_deg, mean_swift_skew, median_strike_deg,
        anisotropy_flag.

        ``anisotropy_flag`` is True when
        ``|mean_ratio_log10| > ratio_threshold`` OR
        ``mean_swift_skew > skew_threshold``.
    """
    detail = analyze_anisotropy(
        sites,
        ratio_threshold=ratio_threshold,
        skew_threshold=skew_threshold,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    if detail.empty:
        return pd.DataFrame(columns=_TABLE_COLS)

    rows: list[dict] = []
    for station, grp in detail.groupby("station", sort=False):
        ratio = grp["ratio_log10"].dropna()
        pdiff = grp["phase_diff_deg"].dropna()
        skew = grp["swift_skew"].dropna()
        strike = grp["strike_deg"].dropna()

        mean_ratio = float(ratio.mean()) if len(ratio) else np.nan
        max_abs_r = float(ratio.abs().max()) if len(ratio) else np.nan
        mean_pdiff = float(pdiff.mean()) if len(pdiff) else np.nan
        mean_skew = float(skew.mean()) if len(skew) else np.nan
        med_strike = float(strike.median()) if len(strike) else np.nan

        flag = bool(
            (np.isfinite(mean_ratio) and abs(mean_ratio) > ratio_threshold)
            or (np.isfinite(mean_skew) and mean_skew > skew_threshold)
        )

        rows.append(
            {
                "station": station,
                "n_freq": len(grp),
                "mean_ratio_log10": mean_ratio,
                "max_abs_ratio_log10": max_abs_r,
                "mean_phase_diff_deg": mean_pdiff,
                "mean_swift_skew": mean_skew,
                "median_strike_deg": med_strike,
                "anisotropy_flag": flag,
            }
        )

    return pd.DataFrame(rows, columns=_TABLE_COLS)


def plot_anisotropy(
    sites: Any,
    *,
    metric: str = "ratio_log10",
    ratio_threshold: float = ANISO_RATIO_THRESH,
    skew_threshold: float = SWIFT_SKEW_THRESH,
    cmap: str = "RdBu_r",
    figsize: tuple = (10, 5),
    period_axis: bool = True,
    log_y: bool = True,
    contour_zero: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax=None,
):
    """
    Plot anisotropy metric pseudo-section (station × frequency).

    Parameters
    ----------
    sites : Sites | list
    metric : str
        Column from :func:`analyze_anisotropy` to map to colour:
        ``"ratio_log10"`` (default), ``"swift_skew"``,
        ``"phase_diff_deg"``, or ``"strike_deg"``.
    ratio_threshold : float
    skew_threshold : float
    cmap : str
        Colormap (default ``"RdBu_r"`` — diverging, centred at 0).
    period_axis : bool
        Show period on y-axis (default) rather than frequency.
    log_y : bool
        Logarithmic y-axis.
    contour_zero : bool
        Draw a white contour at value = 0 (relevant for ratio_log10).
    ax : matplotlib.axes.Axes or None

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    df = analyze_anisotropy(
        sites,
        ratio_threshold=ratio_threshold,
        skew_threshold=skew_threshold,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if df.empty or metric not in df.columns:
        ax.set_xlabel("Station")
        ax.set_ylabel("Period (s)" if period_axis else "Frequency (Hz)")
        ax.set_title(f"Anisotropy — {metric} (no data)")
        return ax

    stations = list(dict.fromkeys(df["station"]))
    s_idx = {s: k for k, s in enumerate(stations)}
    freqs_all = np.sort(df["freq_hz"].unique())
    f_idx = {f: k for k, f in enumerate(freqs_all)}

    grid = np.full((len(freqs_all), len(stations)), np.nan)
    for _, row in df.iterrows():
        fi = f_idx.get(row["freq_hz"])
        si = s_idx.get(row["station"])
        if fi is not None and si is not None:
            grid[fi, si] = row[metric]

    y_vals = 1.0 / freqs_all if period_axis else freqs_all
    if period_axis:
        order = np.argsort(y_vals)
        y_vals = y_vals[order]
        grid = grid[order]
    x_vals = np.arange(len(stations))
    X, Y = np.meshgrid(x_vals, y_vals)

    valid = grid[np.isfinite(grid)]
    if len(valid) == 0:
        ax.set_title(f"Anisotropy — {metric} (all NaN)")
        return ax

    vmax = max(abs(valid.max()), abs(valid.min()), 1e-6)
    vmin = -vmax if metric == "ratio_log10" else valid.min()

    norm = Normalize(vmin=vmin, vmax=vmax)
    pcm = ax.pcolormesh(X, Y, grid, cmap=cmap, norm=norm, shading="nearest")
    cb = plt.colorbar(pcm, ax=ax)
    cb.set_label(_METRIC_LABELS.get(metric, metric))

    if contour_zero and metric == "ratio_log10" and np.isfinite(grid).any():
        if grid.shape[0] >= 2 and grid.shape[1] >= 2:
            ax.contour(
                X,
                Y,
                grid,
                levels=[0.0],
                colors="white",
                linewidths=1.2,
                linestyles="--",
            )

    PYCSAMT_STATION_RENDERING.apply(
        ax,
        x_vals,
        stations,
        preset="pseudosection",
        xlim=(-0.5, len(stations) - 0.5),
    )
    if log_y:
        ax.set_yscale("log")
    if period_axis:
        ax.set_ylabel("Period (s)")
        if not ax.yaxis_inverted():
            ax.invert_yaxis()
    else:
        ax.set_ylabel("Frequency (Hz)")

    ax.set_title(f"Anisotropy: {_METRIC_LABELS.get(metric, metric)} (wang2017)")
    return ax


_METRIC_LABELS: dict = {
    "ratio_log10": "log₁₀(ρ_xy/ρ_yx)",
    "swift_skew": "Swift skew S",
    "phase_diff_deg": "Δφ = φ_xy − φ_yx (°)",
    "strike_deg": "Strike angle θ (°)",
}
