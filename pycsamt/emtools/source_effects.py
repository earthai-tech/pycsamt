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
from typing import Any

import numpy as np
import pandas as pd

from ..api.station import PYCSAMT_STATION_RENDERING
from ._core import (
    _apply_each,
    _get_z_block,
    _iter_items,
    _name,
    ensure_sites,
)

__all__ = [
    "BETA_THRESH_PCT",
    "overprint_beta",
    "detect_source_overprint",
    "source_overprint_table",
    "plot_overprint_section",
    # Wang & Lin (2023)
    "normalize_response",
    "correct_near_field",
    "plot_normalized_response",
]

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

_MU0: float = 4.0 * np.pi * 1e-7  # H/m
BETA_THRESH_PCT: float = 3.0  # yan2004: β > 3 % = significant effect

_DETAIL_COLS = [
    "station",
    "freq_hz",
    "period_s",
    "offset_m",
    "rho_a_ohmm",
    "kr",
    "beta_pct",
    "overprint_flag",
]
_TABLE_COLS = [
    "station",
    "n_freq",
    "offset_m",
    "beta_max_pct",
    "beta_mean_pct",
    "n_overprint",
    "overprint_frac",
    "lf_slope",
    "hf_slope",
    "slope_delta",
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
    """Geometric-mean apparent resistivity from off-diagonal Z (Ω·m).

    Falls back to whichever off-diagonal component is actually
    populated when the other is identically zero, e.g. scalar
    single-component CSAMT surveys that only record Zxy (or Zyx).
    Naively multiplying rxy and ryx in that case silently collapses
    to zero rather than signalling missing data.
    """
    zxy = np.abs(z[:, 0, 1])
    zyx = np.abs(z[:, 1, 0])
    rxy = 0.2 * zxy**2 / np.maximum(fr, 1e-24)
    ryx = 0.2 * zyx**2 / np.maximum(fr, 1e-24)
    xy_ok = zxy > 0
    yx_ok = zyx > 0
    both = xy_ok & yx_ok
    out = np.full(rxy.shape, np.nan)
    out[both] = np.sqrt(np.maximum(rxy[both] * ryx[both], 1e-12))
    only_xy = xy_ok & ~yx_ok
    out[only_xy] = rxy[only_xy]
    only_yx = yx_ok & ~xy_ok
    out[only_yx] = ryx[only_yx]
    return out


def _k1_scalar(rho: float, freq: float) -> complex:
    """Complex wavenumber k₁ = √(iωμ₀/ρ)."""
    omega = 2.0 * np.pi * freq
    return complex(np.sqrt(1j * omega * _MU0 / rho))


def _kr(rho: np.ndarray, freq: np.ndarray, offset: float) -> np.ndarray:
    """Dimensionless field-zone parameter |k₁|·r = r / δ_Bostick."""
    omega = 2.0 * np.pi * np.asarray(freq, dtype=float)
    k1_abs = np.sqrt(omega * _MU0 / np.maximum(rho, 1e-12))
    return k1_abs * abs(offset)


def _resolve_offset(ed: Any, source_offset: Any, station: str) -> float | None:
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
    from scipy.special import (  # noqa: F401 — lazy import
        iv,
        kv,
    )

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
    h = max(abs(offset) * dh_frac, 0.5)  # step ≥ 0.5 m

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
    rho: float | np.ndarray,
    freq: float | np.ndarray,
    offset: float | np.ndarray,
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
    rho = np.asarray(rho, dtype=float)
    freq = np.asarray(freq, dtype=float)
    offset = np.asarray(offset, dtype=float)

    shape = np.broadcast_shapes(rho.shape, freq.shape, offset.shape)
    rho_b = np.broadcast_to(rho, shape).ravel()
    freq_b = np.broadcast_to(freq, shape).ravel()
    off_b = np.broadcast_to(offset, shape).ravel()

    result = np.empty(rho_b.size, dtype=float)
    for i in range(rho_b.size):
        result[i] = (
            _beta_Ey_scalar(
                float(rho_b[i]), float(freq_b[i]), float(off_b[i]), dh_frac
            )
            * 100.0
        )

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
    rows: list[dict] = []

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
            f = float(freqs[j])
            ra = float(rho[j])

            if off is not None and off > 0.0 and ra > 0.0 and f > 0.0:
                kr_val = float(_kr(np.array([ra]), np.array([f]), off)[0])
                beta_val = _beta_Ey_scalar(ra, f, off) * 100.0
            else:
                kr_val = np.nan
                beta_val = np.nan

            rows.append(
                {
                    "station": station,
                    "freq_hz": f,
                    "period_s": 1.0 / f if f > 0 else np.nan,
                    "offset_m": float(off) if off is not None else np.nan,
                    "rho_a_ohmm": ra,
                    "kr": kr_val,
                    "beta_pct": beta_val,
                    "overprint_flag": bool(
                        np.isfinite(beta_val) and beta_val > beta_threshold
                    ),
                }
            )

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

    rows: list[dict] = []
    for station, grp in detail.groupby("station", sort=False):
        grp = grp.sort_values("freq_hz")
        valid_beta = grp["beta_pct"].dropna()
        n_ov = int(grp["overprint_flag"].sum())
        n_freq = len(grp)

        log_f = np.log10(np.maximum(grp["freq_hz"].values, 1e-12))
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

        rows.append(
            {
                "station": station,
                "n_freq": n_freq,
                "offset_m": grp["offset_m"].iloc[0],
                "beta_max_pct": float(valid_beta.max())
                if len(valid_beta)
                else np.nan,
                "beta_mean_pct": float(valid_beta.mean())
                if len(valid_beta)
                else np.nan,
                "n_overprint": n_ov,
                "overprint_frac": n_ov / n_freq if n_freq else np.nan,
                "lf_slope": lf_slope,
                "hf_slope": hf_slope,
                "slope_delta": slope_delta,
                "overprint_flag": n_ov > 0,
            }
        )

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
    s_idx = {s: k for k, s in enumerate(stations)}
    df = df.copy()
    df["_sx"] = df["station"].map(s_idx)

    freqs_all = np.sort(df["freq_hz"].unique())
    grid_beta = np.full((len(freqs_all), len(stations)), np.nan)
    f_idx = {f: k for k, f in enumerate(freqs_all)}
    for _, row in df.iterrows():
        fi = f_idx.get(row["freq_hz"])
        si = s_idx.get(row["station"])
        if fi is not None and si is not None and np.isfinite(row["beta_pct"]):
            grid_beta[fi, si] = row["beta_pct"]

    y_vals = 1.0 / freqs_all if period_axis else freqs_all
    if period_axis:
        order = np.argsort(y_vals)
        y_vals = y_vals[order]
        grid_beta = grid_beta[order]
    x_vals = np.arange(len(stations))

    vmin = (
        max(grid_beta[np.isfinite(grid_beta)].min(), 1e-3)
        if np.isfinite(grid_beta).any()
        else 1e-3
    )
    vmax = (
        grid_beta[np.isfinite(grid_beta)].max()
        if np.isfinite(grid_beta).any()
        else 100.0
    )

    norm = (
        LogNorm(vmin=max(vmin, 1e-3), vmax=max(vmax, vmin * 10))
        if log_color
        else Normalize(vmin=vmin, vmax=vmax)
    )
    X, Y = np.meshgrid(x_vals, y_vals)

    pcm = ax.pcolormesh(
        X, Y, grid_beta, cmap=cmap, norm=norm, shading="nearest"
    )
    plt.colorbar(pcm, ax=ax, label="β_Ey (%)")

    if (
        contour_beta
        and np.isfinite(grid_beta).any()
        and grid_beta.shape[0] >= 2
        and grid_beta.shape[1] >= 2
    ):
        valid_levels = [lvl for lvl in beta_levels if vmin < lvl < vmax]
        thr_label = [beta_threshold] if vmin < beta_threshold < vmax else []
        for lvl in valid_levels:
            ax.contour(
                X,
                Y,
                grid_beta,
                levels=[lvl],
                colors="grey",
                linewidths=0.7,
                alpha=0.6,
            )
        for lvl in thr_label:
            ax.contour(
                X,
                Y,
                grid_beta,
                levels=[lvl],
                colors="white",
                linewidths=1.4,
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
        ax.invert_yaxis()
    else:
        ax.set_ylabel("Frequency (Hz)")

    ax.set_title(
        f"Source overprint β_Ey — threshold {beta_threshold:.1f} % (yan2004)"
    )
    return ax


# =============================================================================
# Wang & Lin (2023) — normalized response analysis and near-field correction
# =============================================================================

_NORM_COLS = [
    "station",
    "freq_hz",
    "period_s",
    "offset_m",
    "rho_a_ohmm",
    "rho_n",
    "phi_obs_deg",
    "phi_ref_deg",
    "phi_diff_deg",
    "zone",
    "kr",
]

# ─────────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────────


def _rho_a_comp(z: np.ndarray, fr: np.ndarray, comp: str) -> np.ndarray:
    """Apparent resistivity for a named Z component."""
    if comp == "xy":
        return 0.2 * np.abs(z[:, 0, 1]) ** 2 / np.maximum(fr, 1e-24)
    if comp == "yx":
        return 0.2 * np.abs(z[:, 1, 0]) ** 2 / np.maximum(fr, 1e-24)
    return _rho_a_det(z, fr)


def _phase_comp_deg(z: np.ndarray, comp: str) -> np.ndarray:
    """Phase [°] for a named Z component."""
    if comp == "xy":
        return np.degrees(np.angle(z[:, 0, 1]))
    if comp == "yx":
        return np.degrees(np.angle(z[:, 1, 0]))
    phi_xy = np.degrees(np.angle(z[:, 0, 1]))
    phi_yx = np.degrees(np.angle(z[:, 1, 0]))
    return 0.5 * (phi_xy + phi_yx)


def _skin_depth_wang(rho: np.ndarray, freq: np.ndarray) -> np.ndarray:
    """δ = 503 √(ρ / f)  [m]  (Wang & Lin 2023, eq. 1)."""
    return 503.0 * np.sqrt(np.maximum(rho, 1e-6) / np.maximum(freq, 1e-12))


def _kr_wang(rho: np.ndarray, freq: np.ndarray, offset: float) -> np.ndarray:
    """r / δ using skin depth δ = 503 √(ρ/f)  (Wang & Lin 2023)."""
    return abs(offset) / np.maximum(_skin_depth_wang(rho, freq), 1e-6)


def _zone_wang(kr: np.ndarray) -> np.ndarray:
    """Field-zone labels from Wang & Lin (2023) thresholds (0.5δ, 4δ)."""
    return np.where(
        kr >= 4.0,
        "far",
        np.where(kr >= 0.5, "transition", "near"),
    )


def _F_complex(rho: np.ndarray, freq: np.ndarray, offset: float) -> np.ndarray:
    """
    Complex near-field factor F(p) = 1 − 3/p² + 3/p³  (equatorial HED).

    p = k · r,  k = √(i·ω·μ₀ / ρ_a)
    F → 1 in the far field; diverges in the geometric near field.
    """
    omega = 2.0 * np.pi * np.maximum(np.asarray(freq, dtype=float), 1e-12)
    k_abs = np.sqrt(
        omega * _MU0 / np.maximum(np.asarray(rho, dtype=float), 1e-6)
    )
    p = k_abs * abs(offset) * (1.0 + 1j) / np.sqrt(2.0)
    tiny = np.abs(p) < 1e-10
    if np.any(tiny):
        p = p.copy()
        p[tiny] = 1e-10 * (1.0 + 1j)
    return 1.0 - 3.0 / p**2 + 3.0 / p**3


def _label_pseudo_ax(
    ax: Any,
    stations: list[str],
    all_y: np.ndarray,
    period_axis: bool,
    title: str,
) -> None:
    """Apply station xticks, y-ticks, and title to a pseudosection axes."""
    ax.set_xticks(range(len(stations)))
    ax.set_xticklabels(stations, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Station")
    n_ytick = min(8, len(all_y))
    step = max(1, len(all_y) // n_ytick)
    tick_idx = np.arange(0, len(all_y), step)
    ax.set_yticks(tick_idx)
    ax.set_yticklabels([f"{all_y[k]:.3g}" for k in tick_idx], fontsize=8)
    ax.set_ylabel("Period (s)" if period_axis else "Frequency (Hz)")
    ax.set_title(title)


# ─────────────────────────────────────────────────────────────────────────────
# normalize_response
# ─────────────────────────────────────────────────────────────────────────────


def normalize_response(
    sites: Any,
    rho_ref: float = 100.0,
    source_offset: Any = None,
    *,
    comp: str = "det",
    phi_ref_deg: float = 45.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Normalized apparent resistivity and subtracted phase (Wang & Lin 2023).

    For each (station, frequency) pair computes::

        ρ_n    = ρ_obs / ρ_ref
        φ_diff = φ_obs − φ_ref

    and classifies the measurement zone using the skin-depth formula
    proposed by Wang & Lin (2023, eq. 1):

        δ = 503 √(ρ_a / f)  [m]

    with thresholds: near (r/δ < 0.5), transition (0.5–4), far (>4).

    Parameters
    ----------
    sites : Sites | list
        EDI-like objects or a ``Sites`` container.
    rho_ref : float
        Reference half-space resistivity [Ω·m] (default 100).
    source_offset : float | dict | None
        Source–receiver offset **r** [m].  A dict maps
        ``{station: r}``.  If *None*, zone and kr are NaN.
    comp : {"det", "xy", "yx"}
        Impedance component used for ρ_a and φ
        (``"det"`` = geometric-mean determinant).
    phi_ref_deg : float
        Reference half-space phase [°].  45° (default) is the far-field
        plane-wave value for a homogeneous 1-D half-space.

    Returns
    -------
    pandas.DataFrame
        Columns: station, freq_hz, period_s, offset_m, rho_a_ohmm, rho_n,
        phi_obs_deg, phi_ref_deg, phi_diff_deg, zone, kr.
        ``zone`` / ``kr`` are ``None`` / NaN when no offset is available.

    References
    ----------
    Wang & Lin (2023), *Geophysics* **88**(6), E215–E230.
    """
    if comp not in ("det", "xy", "yx"):
        raise ValueError(f"comp must be 'det', 'xy', or 'yx'; got {comp!r}")

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    rows: list[dict] = []

    for i, ed in enumerate(_iter_items(S)):
        ed = _unwrap(ed)
        station = _name(ed, i)
        off = _resolve_offset(ed, source_offset, station)

        _, z, fr = _get_z_block(ed)
        if z is None or fr is None or fr.size == 0:
            continue

        rho_a = _rho_a_comp(z, fr, comp)
        phi = _phase_comp_deg(z, comp)
        rho_n = rho_a / max(float(rho_ref), 1e-12)
        phi_diff = phi - float(phi_ref_deg)

        has_off = off is not None and off > 0.0
        kr_arr = (
            _kr_wang(rho_a, fr, off) if has_off else np.full(fr.size, np.nan)
        )
        zone_arr = (
            _zone_wang(kr_arr)
            if has_off
            else np.full(fr.size, None, dtype=object)
        )

        for j in range(fr.size):
            rows.append(
                {
                    "station": station,
                    "freq_hz": float(fr[j]),
                    "period_s": 1.0 / max(float(fr[j]), 1e-12),
                    "offset_m": float(off) if has_off else np.nan,
                    "rho_a_ohmm": float(rho_a[j]),
                    "rho_n": float(rho_n[j]),
                    "phi_obs_deg": float(phi[j]),
                    "phi_ref_deg": float(phi_ref_deg),
                    "phi_diff_deg": float(phi_diff[j]),
                    "zone": zone_arr[j],
                    "kr": float(kr_arr[j]) if has_off else np.nan,
                }
            )

    if not rows:
        return pd.DataFrame(columns=_NORM_COLS)
    return pd.DataFrame(rows, columns=_NORM_COLS)


# ─────────────────────────────────────────────────────────────────────────────
# correct_near_field
# ─────────────────────────────────────────────────────────────────────────────


def correct_near_field(
    sites: Any,
    source_offset: Any,
    *,
    inplace: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> Any:
    """
    Correct impedance tensor for CSAMT near-field contamination.

    Divides each element of **Z** by the complex near-field factor F(p)::

        Z_corrected = Z_obs / F(p)

    where F(p) = 1 − 3/p² + 3/p³ is the equatorial HED transfer-function
    ratio and p = k · r,  k = √(i·ω·μ₀ / ρ_a).  In the far field
    F(p) → 1 so no correction is applied; in the near/transition zone the
    correction restores the plane-wave equivalent impedance.

    Uses :func:`~pycsamt.emtools._core._apply_each` to apply the
    per-site correction and return a new ``Sites`` (or modify in-place).

    Parameters
    ----------
    sites : Sites | list
        EDI-like objects or a ``Sites`` container.
    source_offset : float | dict
        Source–receiver separation **r** [m].  Dict maps
        ``{station: r}``.
    inplace : bool, default False
        Modify **Z.z** in-place and return the original ``Sites``.

    Returns
    -------
    pycsamt.site.base.Sites
        Sites with corrected impedance tensors.

    References
    ----------
    Wang & Lin (2023), *Geophysics* **88**(6), E215–E230.
    Chen & Yan (2005), eqs. (8)–(10).
    """
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    def _one(Si: Any) -> None:
        for ii, edd in enumerate(_iter_items(Si)):
            edd_raw = _unwrap(edd)
            stn = _name(edd_raw, ii)
            Z_wrap, z, fr = _get_z_block(edd_raw)
            if Z_wrap is None or z is None or fr is None:
                return
            off = _resolve_offset(edd_raw, source_offset, stn)
            if off is None or off <= 0.0:
                if verbose > 0:
                    warnings.warn(
                        f"correct_near_field: no source offset for "
                        f"'{stn}'; station skipped.",
                        UserWarning,
                        stacklevel=4,
                    )
                return
            rho_a = _rho_a_det(z, fr)
            F = _F_complex(rho_a, fr, off)
            F_mag = np.abs(F)
            # guard: don't amplify Z by more than 1000× (|F| floor at 1e-3)
            F_eff = np.where(
                F_mag >= 1e-3, F, F * (1e-3 / np.maximum(F_mag, 1e-30))
            )
            Z_wrap.z = z / F_eff[:, None, None]

    return _apply_each(S, _one, inplace=inplace, verbose=verbose)


# ─────────────────────────────────────────────────────────────────────────────
# plot_normalized_response
# ─────────────────────────────────────────────────────────────────────────────


def plot_normalized_response(
    sites: Any,
    rho_ref: float = 100.0,
    source_offset: Any = None,
    *,
    comp: str = "det",
    phi_ref_deg: float = 45.0,
    period_axis: bool = True,
    figsize: tuple = (12.0, 5.0),
    cmap_rho: str = "RdBu_r",
    cmap_phi: str = "RdBu",
    rho_n_lim: tuple | None = None,
    phi_diff_lim: tuple | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    axes: Any = None,
) -> tuple:
    """
    Pseudosection of normalized ρ_a and subtracted phase (Wang & Lin 2023).

    Produces a two-panel figure analogous to Fig. 8(e–f) of Wang & Lin (2023):

    * **Left panel**: ρ_n = ρ_obs / ρ_ref  (centred at 1.0; red = high).
    * **Right panel**: φ_diff = φ_obs − φ_ref [°] (centred at 0°).

    Parameters
    ----------
    sites : Sites | list
    rho_ref : float
        Reference half-space resistivity [Ω·m].
    source_offset : float | dict | None
    comp : {"det", "xy", "yx"}
    phi_ref_deg : float
        Reference half-space phase [°] (default 45°).
    period_axis : bool
        Use period (s) on the y-axis when *True* (default).
    figsize : (float, float), default (12, 5)
    cmap_rho, cmap_phi : str
        Matplotlib colormap names for the two panels.
    rho_n_lim : (vmin, vmax) or None
        Colour limits for ρ_n.  Default: symmetric about 1.
    phi_diff_lim : (vmin, vmax) or None
        Colour limits for φ_diff.  Default: symmetric about 0.
    axes : (ax1, ax2) or None
        Draw on existing axes; created if *None*.

    Returns
    -------
    (ax1, ax2) : tuple of matplotlib.axes.Axes

    References
    ----------
    Wang & Lin (2023), *Geophysics* **88**(6), E215–E230 (Figs. 8e–f).
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm

    df = normalize_response(
        sites,
        rho_ref=rho_ref,
        source_offset=source_offset,
        comp=comp,
        phi_ref_deg=phi_ref_deg,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )

    if axes is None:
        _, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    else:
        ax1, ax2 = axes

    if df.empty:
        for ax in (ax1, ax2):
            ax.text(
                0.5,
                0.5,
                "no data",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
        return ax1, ax2

    stations = list(dict.fromkeys(df["station"]))
    y_col = "period_s" if period_axis else "freq_hz"
    all_y = np.sort(df[y_col].unique())
    s_idx = {s: k for k, s in enumerate(stations)}
    y_idx_map = {float(v): k for k, v in enumerate(all_y)}

    grid_rho_n = np.full((len(all_y), len(stations)), np.nan)
    grid_phi = np.full((len(all_y), len(stations)), np.nan)

    for _, row in df.iterrows():
        si = s_idx.get(row["station"])
        yi = y_idx_map.get(float(row[y_col]))
        if si is not None and yi is not None:
            grid_rho_n[yi, si] = row["rho_n"]
            grid_phi[yi, si] = row["phi_diff_deg"]

    X, Y = np.meshgrid(np.arange(len(stations)), np.arange(len(all_y)))

    # ── ρ_n panel ─────────────────────────────────────────────────────────
    valid_rn = grid_rho_n[np.isfinite(grid_rho_n)]
    if rho_n_lim is None:
        vdev = max(
            float(np.nanmax(np.abs(valid_rn - 1.0))) if len(valid_rn) else 1.0,
            0.01,
        )
        rn_min, rn_max = 1.0 - vdev, 1.0 + vdev
    else:
        rn_min, rn_max = rho_n_lim
    try:
        norm_rn = TwoSlopeNorm(
            vmin=rn_min, vcenter=1.0, vmax=max(rn_max, 1.0 + 1e-6)
        )
    except Exception:
        norm_rn = Normalize(vmin=rn_min, vmax=rn_max)
    pcm1 = ax1.pcolormesh(
        X, Y, grid_rho_n, cmap=cmap_rho, norm=norm_rn, shading="nearest"
    )
    plt.colorbar(pcm1, ax=ax1, label="ρ_n = ρ_obs / ρ_ref")
    _label_pseudo_ax(
        ax1,
        stations,
        all_y,
        period_axis,
        f"Normalized ρ_a  (ρ_ref = {rho_ref:.0f} Ω·m)",
    )

    # ── φ_diff panel ──────────────────────────────────────────────────────
    valid_pd = grid_phi[np.isfinite(grid_phi)]
    if phi_diff_lim is None:
        vdev_p = max(
            float(np.nanmax(np.abs(valid_pd))) if len(valid_pd) else 10.0, 1.0
        )
        pd_min, pd_max = -vdev_p, vdev_p
    else:
        pd_min, pd_max = phi_diff_lim
    try:
        norm_pd = TwoSlopeNorm(
            vmin=pd_min, vcenter=0.0, vmax=max(pd_max, 1e-6)
        )
    except Exception:
        norm_pd = Normalize(vmin=pd_min, vmax=pd_max)
    pcm2 = ax2.pcolormesh(
        X, Y, grid_phi, cmap=cmap_phi, norm=norm_pd, shading="nearest"
    )
    plt.colorbar(pcm2, ax=ax2, label="φ_diff = φ_obs − φ_ref  (°)")
    _label_pseudo_ax(
        ax2,
        stations,
        all_y,
        period_axis,
        f"Subtracted phase  (φ_ref = {phi_ref_deg:.0f}°)",
    )

    return ax1, ax2
