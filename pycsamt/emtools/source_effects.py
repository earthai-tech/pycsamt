"""
CSAMT source overprint and shadow effect analysis.

Implements the analytical β-ratio method and spectral slope criterion:

  yan2004 : Yan & Fu (2004), "An analytical method to estimate shadow and
            source overprint effects in CSAMT sounding",
            Geophysics 69(1), 161–163.
  da2016  : Da et al. (2016), "Modeling and analysis of CSAMT field source
            effect and its characteristics", J. Geophys. Eng. 13, 49–58.

The β ratio is the ground-wave / surface-wave amplitude ratio at the
receiver location.  When β > 3 % the shadow / overprint effect may be
significant (yan2004 threshold).  Companion spectral slope analysis
(da2016) flags low-frequency ρ_a anomalies that are characteristic of
a resistivity contrast beneath the source dipole.
"""
from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from ._core import (
    ensure_sites,
    _iter_items,
    _get_z_block,
    _name,
)

__all__ = [
    "BETA_THRESH_PCT",
    "overprint_beta",
    "detect_source_overprint",
    "source_overprint_table",
    "plot_overprint_section",
]

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

_MU0: float = 4.0 * np.pi * 1e-7   # H/m
BETA_THRESH_PCT: float = 3.0        # yan2004: β > 3 % = significant effect

_DETAIL_COLS = [
    "station", "freq_hz", "period_s", "offset_m",
    "rho_a_ohmm", "kr", "beta_pct", "overprint_flag",
]
_TABLE_COLS = [
    "station", "n_freq", "offset_m",
    "beta_max_pct", "beta_mean_pct",
    "n_overprint", "overprint_frac",
    "lf_slope", "hf_slope", "slope_delta",
    "overprint_flag",
]

# ─────────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────────

def _unwrap(ed: Any) -> Any:
    edi = getattr(ed, "edi", None)
    if edi is not None and hasattr(edi, "Z"):
        return edi
    return ed


def _rho_a_det(z: np.ndarray, fr: np.ndarray) -> np.ndarray:
    rxy = 0.2 * np.abs(z[:, 0, 1]) ** 2 / np.maximum(fr, 1e-24)
    ryx = 0.2 * np.abs(z[:, 1, 0]) ** 2 / np.maximum(fr, 1e-24)
    return np.sqrt(np.maximum(rxy * ryx, 1e-12))


def _k1_scalar(rho: float, freq: float) -> complex:
    """Complex wavenumber k₁ = √(iωμ₀/ρ)."""
    omega = 2.0 * np.pi * freq
    return complex(np.sqrt(1j * omega * _MU0 / rho))


def _kr(rho: np.ndarray, freq: np.ndarray, offset: float) -> np.ndarray:
    """Dimensionless field-zone parameter |k₁|·r = r / δ_Bostick."""
    omega = 2.0 * np.pi * np.asarray(freq, dtype=float)
    k1_abs = np.sqrt(omega * _MU0 / np.maximum(rho, 1e-12))
    return k1_abs * abs(offset)


def _resolve_offset(ed: Any, source_offset: Any, station: str) -> Optional[float]:
    if isinstance(source_offset, dict):
        return source_offset.get(station, None)
    if source_offset is not None:
        return float(source_offset)
    for attr in ("source_offset", "offset", "dist"):
        val = getattr(ed, attr, None)
        if val is not None:
            try:
                return float(val)
            except (TypeError, ValueError):
                pass
    return None


def _bessel_I0K0(p: complex, q: complex) -> complex:
    """I₀(p) K₀(q) with complex arguments via scipy (AMOS)."""
    from scipy.special import iv, kv  # noqa: F401 — lazy import
    return complex(iv(0, p)) * complex(kv(0, q))


def _P_func(x: float, y: float, z: float, k1: complex) -> complex:
    """Sommerfeld (ground-wave) term  P = e^{-k₁r₃D} / r₃D."""
    r3 = float(np.sqrt(x * x + y * y + z * z))
    return np.exp(-k1 * r3) / r3


def _N_func(x: float, y: float, z: float, k1: complex) -> complex:
    """Foster (surface-wave) term  N = I₀(p) K₀(q)."""
    r3 = float(np.sqrt(x * x + y * y + z * z))
    p = k1 * (r3 + z) / 2.0
    q = k1 * (r3 - z) / 2.0
    return _bessel_I0K0(p, q)


def _beta_Ey_scalar(
    rho: float,
    freq: float,
    offset: float,
    dh_frac: float = 1e-3,
) -> float:
    """
    Ground-wave/surface-wave ratio β_Ey (dimensionless) at the surface
    receiver located at (x=offset, y=0, z=0), broadside to a y-directed
    horizontal electric dipole at the origin.

    Uses central finite differences to evaluate:
        β_Ey = |∂²P/∂z²| / |∂³N/∂x²∂z|      (yan2004 eq. 6)

    Returns the ratio as a fraction; multiply by 100 for %.
    """
    if offset <= 0.0 or rho <= 0.0 or freq <= 0.0:
        return np.nan

    k1 = _k1_scalar(rho, freq)
    x0, y0, z0 = float(offset), 0.0, 0.0
    h = max(abs(offset) * dh_frac, 0.5)   # step ≥ 0.5 m

    # ∂²P/∂z² at (x0, y0, z0) — central difference
    d2P_dz2 = (
        _P_func(x0, y0, z0 + h, k1)
        - 2.0 * _P_func(x0, y0, z0, k1)
        + _P_func(x0, y0, z0 - h, k1)
    ) / (h * h)

    # ∂³N/∂x²∂z = ∂/∂z[ (∂²N/∂x²) ]
    def _d2N_dx2(zz: float) -> complex:
        return (
            _N_func(x0 + h, y0, zz, k1)
            - 2.0 * _N_func(x0, y0, zz, k1)
            + _N_func(x0 - h, y0, zz, k1)
        ) / (h * h)

    d3N_dx2dz = (_d2N_dx2(z0 + h) - _d2N_dx2(z0 - h)) / (2.0 * h)

    denom = abs(d3N_dx2dz)
    if denom < 1e-200:
        return 0.0
    beta = float(abs(d2P_dz2) / denom)
    # guard against numerical blow-up at near-field (kr → 0)
    return min(beta, 1e3)


def _log_slope(log_f: np.ndarray, log_rho: np.ndarray) -> float:
    """Median d(log10 ρ) / d(log10 f) via 1-D linear regression."""
    if log_f.size < 2:
        return np.nan
    coef = np.polyfit(log_f, log_rho, 1)
    return float(coef[0])


# ─────────────────────────────────────────────────────────────────────────────
# Public: pure-math interface
# ─────────────────────────────────────────────────────────────────────────────

def overprint_beta(
    rho: Union[float, np.ndarray],
    freq: Union[float, np.ndarray],
    offset: Union[float, np.ndarray],
    *,
    dh_frac: float = 1e-3,
) -> np.ndarray:
    """
    Ground-wave / surface-wave amplitude ratio β_Ey (%).

    Evaluates equation (6) of Yan & Fu (2004) analytically at the
    surface receiver position broadside to the source dipole.

    Parameters
    ----------
    rho : float or ndarray
        Half-space apparent resistivity [Ω·m].
    freq : float or ndarray
        Frequency [Hz].
    offset : float or ndarray
        Source–receiver horizontal offset r [m].
    dh_frac : float
        Step size as a fraction of *offset* used for numerical
        differentiation (default 1e-3).

    Returns
    -------
    beta_pct : ndarray
        β × 100  [%].  Values above ``BETA_THRESH_PCT`` (3 %) indicate
        potential shadow / source overprint (yan2004).

    Notes
    -----
    The function uses central finite differences to evaluate the partial
    derivatives of the Sommerfeld term P = e^{−k₁r}/r and the Foster
    term N = I₀(p) K₀(q), where k₁ = √(iωμ₀/ρ) is the complex
    wavenumber and p, q are related to the 3-D distance and depth.
    """
    rho    = np.asarray(rho,    dtype=float)
    freq   = np.asarray(freq,   dtype=float)
    offset = np.asarray(offset, dtype=float)

    shape  = np.broadcast_shapes(rho.shape, freq.shape, offset.shape)
    rho_b  = np.broadcast_to(rho,    shape).ravel()
    freq_b = np.broadcast_to(freq,   shape).ravel()
    off_b  = np.broadcast_to(offset, shape).ravel()

    result = np.empty(rho_b.size, dtype=float)
    for i in range(rho_b.size):
        result[i] = _beta_Ey_scalar(
            float(rho_b[i]), float(freq_b[i]), float(off_b[i]), dh_frac
        ) * 100.0

    return result.reshape(shape) if shape else result.ravel()[0]


# ─────────────────────────────────────────────────────────────────────────────
# Public: sites-based interface
# ─────────────────────────────────────────────────────────────────────────────

def detect_source_overprint(
    sites: Any,
    source_offset: Any = None,
    *,
    beta_threshold: float = BETA_THRESH_PCT,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Per-frequency source overprint β index for a set of CSAMT sites.

    Computes the ground-wave / surface-wave ratio β_Ey (yan2004) for
    every measurement frequency at each site and returns a long-form
    DataFrame.

    Parameters
    ----------
    sites : Sites | list
        EDI-like objects or a ``Sites`` container.
    source_offset : float | dict | None
        Source–receiver offset [m].  A scalar applies to all sites;
        a dict maps ``{station: offset}``.  If *None* the function
        tries to read the offset from site attributes
        (``source_offset``, ``offset``, ``dist``).
    beta_threshold : float
        β [%] above which the overprint flag is raised (default 3.0).

    Returns
    -------
    pd.DataFrame
        Columns: station, freq_hz, period_s, offset_m, rho_a_ohmm,
        kr, beta_pct, overprint_flag.
        Rows with unknown offset have NaN in kr/beta_pct.
    """
    sites = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    rows: List[Dict] = []

    for i, ed in enumerate(_iter_items(sites)):
        ed = _unwrap(ed)
        station = _name(ed, i)
        off = _resolve_offset(ed, source_offset, station)

        if off is None:
            warnings.warn(
                f"detect_source_overprint: no source offset for '{station}'; "
                "beta_pct / kr will be NaN.",
                UserWarning,
                stacklevel=2,
            )

        _, z_block, freqs = _get_z_block(ed)
        if z_block is None or freqs is None or freqs.size == 0:
            continue

        rho = _rho_a_det(z_block, freqs)

        for j in range(freqs.size):
            f  = float(freqs[j])
            ra = float(rho[j])

            if off is not None and off > 0.0 and ra > 0.0 and f > 0.0:
                kr_val   = float(_kr(np.array([ra]), np.array([f]), off)[0])
                beta_val = _beta_Ey_scalar(ra, f, off) * 100.0
            else:
                kr_val   = np.nan
                beta_val = np.nan

            rows.append({
                "station":       station,
                "freq_hz":       f,
                "period_s":      1.0 / f if f > 0 else np.nan,
                "offset_m":      float(off) if off is not None else np.nan,
                "rho_a_ohmm":    ra,
                "kr":            kr_val,
                "beta_pct":      beta_val,
                "overprint_flag": bool(
                    np.isfinite(beta_val) and beta_val > beta_threshold
                ),
            })

    if not rows:
        return pd.DataFrame(columns=_DETAIL_COLS)

    df = pd.DataFrame(rows, columns=_DETAIL_COLS)
    return df


def source_overprint_table(
    sites: Any,
    source_offset: Any = None,
    *,
    beta_threshold: float = BETA_THRESH_PCT,
    f_split: float = 1.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Per-station summary of source overprint metrics.

    In addition to the maximum and mean β values (yan2004), the table
    includes the log-log ρ_a–frequency slope in the low-frequency (LF)
    and high-frequency (HF) bands and their difference (da2016).  A
    strongly negative ``slope_delta`` (LF slope << HF slope) indicates
    a resistivity contrast beneath the source (da2016 §2.2–2.3).

    Parameters
    ----------
    sites : Sites | list
    source_offset : float | dict | None
    beta_threshold : float
        β [%] threshold (default ``BETA_THRESH_PCT`` = 3.0).
    f_split : float
        Frequency [Hz] dividing LF from HF bands for slope analysis.

    Returns
    -------
    pd.DataFrame
        Columns: station, n_freq, offset_m, beta_max_pct, beta_mean_pct,
        n_overprint, overprint_frac, lf_slope, hf_slope, slope_delta,
        overprint_flag.
    """
    detail = detect_source_overprint(
        sites,
        source_offset,
        beta_threshold=beta_threshold,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    if detail.empty:
        return pd.DataFrame(columns=_TABLE_COLS)

    rows: List[Dict] = []
    for station, grp in detail.groupby("station", sort=False):
        grp = grp.sort_values("freq_hz")
        valid_beta = grp["beta_pct"].dropna()
        n_ov   = int(grp["overprint_flag"].sum())
        n_freq = len(grp)

        log_f   = np.log10(np.maximum(grp["freq_hz"].values, 1e-12))
        log_rho = np.log10(np.maximum(grp["rho_a_ohmm"].values, 1e-12))

        mask_lf = grp["freq_hz"].values < f_split
        mask_hf = grp["freq_hz"].values >= f_split

        lf_slope = _log_slope(log_f[mask_lf], log_rho[mask_lf])
        hf_slope = _log_slope(log_f[mask_hf], log_rho[mask_hf])
        slope_delta = (
            float(lf_slope - hf_slope)
            if np.isfinite(lf_slope) and np.isfinite(hf_slope)
            else np.nan
        )

        rows.append({
            "station":       station,
            "n_freq":        n_freq,
            "offset_m":      grp["offset_m"].iloc[0],
            "beta_max_pct":  float(valid_beta.max())  if len(valid_beta) else np.nan,
            "beta_mean_pct": float(valid_beta.mean()) if len(valid_beta) else np.nan,
            "n_overprint":   n_ov,
            "overprint_frac":n_ov / n_freq if n_freq else np.nan,
            "lf_slope":      lf_slope,
            "hf_slope":      hf_slope,
            "slope_delta":   slope_delta,
            "overprint_flag":n_ov > 0,
        })

    return pd.DataFrame(rows, columns=_TABLE_COLS)


def plot_overprint_section(
    sites: Any,
    source_offset: Any = None,
    *,
    beta_threshold: float = BETA_THRESH_PCT,
    log_color: bool = True,
    cmap: str = "hot_r",
    figsize: tuple = (10, 5),
    period_axis: bool = True,
    log_y: bool = True,
    contour_beta: bool = True,
    beta_levels: tuple = (1.0, 3.0, 10.0, 30.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax=None,
):
    """
    Plot source overprint β pseudo-section (station × frequency).

    A colour-coded pseudo-section of the ground-wave / surface-wave
    ratio β_Ey is drawn for each site.  Contour lines at selected β
    levels highlight the overprint-prone zones.

    Parameters
    ----------
    sites : Sites | list
    source_offset : float | dict | None
    beta_threshold : float
        Dashed contour drawn at this level [%] (default 3.0).
    log_color : bool
        Use log₁₀(β) colour scale.
    cmap : str
        Matplotlib colormap name.
    period_axis : bool
        Show periods on the right y-axis when *True*.
    log_y : bool
        Logarithmic frequency axis.
    contour_beta : bool
        Overlay β contour lines.
    beta_levels : tuple
        β [%] values for contour lines.
    ax : matplotlib.axes.Axes or None
        Axes to draw on; created if *None*.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    df = detect_source_overprint(
        sites,
        source_offset,
        beta_threshold=beta_threshold,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if df.empty:
        ax.set_xlabel("Station")
        ax.set_ylabel("Period (s)" if period_axis else "Frequency (Hz)")
        ax.set_title("Source overprint β (no data)")
        return ax

    stations = list(dict.fromkeys(df["station"]))
    s_idx    = {s: k for k, s in enumerate(stations)}
    df = df.copy()
    df["_sx"] = df["station"].map(s_idx)

    freqs_all = np.sort(df["freq_hz"].unique())
    grid_beta = np.full((len(freqs_all), len(stations)), np.nan)
    f_idx     = {f: k for k, f in enumerate(freqs_all)}
    for _, row in df.iterrows():
        fi = f_idx.get(row["freq_hz"])
        si = s_idx.get(row["station"])
        if fi is not None and si is not None and np.isfinite(row["beta_pct"]):
            grid_beta[fi, si] = row["beta_pct"]

    y_vals = 1.0 / freqs_all if period_axis else freqs_all
    x_vals = np.arange(len(stations))

    vmin = max(grid_beta[np.isfinite(grid_beta)].min(), 1e-3) \
        if np.isfinite(grid_beta).any() else 1e-3
    vmax = grid_beta[np.isfinite(grid_beta)].max() \
        if np.isfinite(grid_beta).any() else 100.0

    norm  = LogNorm(vmin=max(vmin, 1e-3), vmax=max(vmax, vmin * 10)) \
        if log_color else Normalize(vmin=vmin, vmax=vmax)
    X, Y  = np.meshgrid(x_vals, y_vals)

    pcm = ax.pcolormesh(X, Y, grid_beta, cmap=cmap, norm=norm, shading="nearest")
    plt.colorbar(pcm, ax=ax, label="β_Ey (%)")

    if (contour_beta
            and np.isfinite(grid_beta).any()
            and grid_beta.shape[0] >= 2
            and grid_beta.shape[1] >= 2):
        valid_levels = [lvl for lvl in beta_levels
                        if vmin < lvl < vmax]
        thr_label = [beta_threshold] if vmin < beta_threshold < vmax else []
        for lvl in valid_levels:
            ax.contour(X, Y, grid_beta, levels=[lvl],
                       colors="grey", linewidths=0.7, alpha=0.6)
        for lvl in thr_label:
            ax.contour(X, Y, grid_beta, levels=[lvl],
                       colors="white", linewidths=1.4,
                       linestyles="--")

    ax.set_xticks(x_vals)
    ax.set_xticklabels(stations, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Station")
    if log_y:
        ax.set_yscale("log")

    if period_axis:
        ax.set_ylabel("Period (s)")
        ax.invert_yaxis()
    else:
        ax.set_ylabel("Frequency (Hz)")

    ax.set_title(
        f"Source overprint β_Ey — threshold {beta_threshold:.1f} %"
        " (yan2004)"
    )
    return ax
