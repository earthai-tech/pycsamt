"""
Phased-array (PAS) transmitter design and radiation pattern analysis for CSAMT.

Implements the element and array-factor formulas from:
  gxac023 : Fan, Zhang & Wang (2022), "A novel phased-array transmitting
             source in controlled-source audio-frequency magnetotellurics",
             J. Geophys. Eng. 19, 595–614.

A traditional CSAMT transmitter is a single-dipole antenna source (SDAS).
The novel PAS consists of N co-linear SDASes with independent phase control.
The energy is focussed into a steerable beam, improving SNR and enlarging
the area of interest (AoI) without increasing transmitter power.

Key formulas
------------
Element pattern (eq. 7):
    F(θ) = [cos(kl cosθ/2) − cos(kl/2)] / sinθ
    θ measured from dipole/array axis (y).  θ=0 → null; θ=90° → maximum.

Array factor (eq. 19):
    AF_n = sin(N ψ/2) / [N sin(ψ/2)],  ψ = k d sin θ_b + β
    θ_b = broadside angle from perpendicular to array.

Beam-steering condition (eq. 23):
    β = −k d sin θ_m

Earth wavenumber for CSAMT (real part of complex k₁):
    k_eff = sqrt(π f μ₀ / ρ)   [m⁻¹]
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ..compat.numpy import trapz as _trapz
from ._core import hide_polar_radius_labels

__all__ = [
    "wavenumber",
    "sdas_element_pattern",
    "array_factor",
    "pas_pattern",
    "beam_steer",
    "steering_angles",
    "sdas_directivity",
    "snr_gain_db",
    "plot_radiation_pattern",
]

_MU0: float = 4.0 * np.pi * 1e-7  # H/m
_C0: float = 299_792_458.0  # m/s  (free-space speed of light)


# ─────────────────────────────────────────────────────────────────────────────
# Wavenumber helper
# ─────────────────────────────────────────────────────────────────────────────


def wavenumber(
    freq: float,
    rho: float | None = None,
) -> float:
    """
    Effective real wavenumber k [m⁻¹] for CSAMT or free-space propagation.

    Parameters
    ----------
    freq : float
        Frequency [Hz].
    rho : float or None
        Half-space resistivity [Ω·m].  If given, returns the earth (CSAMT)
        effective wavenumber ``Re(k₁) = sqrt(π f μ₀ / ρ)``.  If *None*,
        returns the free-space wavenumber ``2π f / c``.

    Returns
    -------
    k : float
        Wavenumber [m⁻¹].

    Notes
    -----
    The complex earth wavenumber is k₁ = √(i ω μ₀ / ρ).  Its real part
    equals |k₁| / √2 = √(π f μ₀ / ρ).  The corresponding wavelength is
    λ = 2π / k_eff ≈ 2π × 503 × √(ρ/f)  [m].
    """
    if rho is None:
        return float(2.0 * np.pi * freq / _C0)
    return float(np.sqrt(np.pi * freq * _MU0 / max(rho, 1e-12)))


# ─────────────────────────────────────────────────────────────────────────────
# Single-dipole (SDAS) element pattern
# ─────────────────────────────────────────────────────────────────────────────


def sdas_element_pattern(
    theta_deg: float | np.ndarray,
    l: float,
    k: float,
    *,
    normalize: bool = True,
) -> np.ndarray:
    """
    Far-field element pattern for a single finite-length SDAS (eq. 7).

    The dipole / array axis is the y-axis.  The angle θ is measured FROM
    the y-axis, so θ = 0° is along the dipole (null) and θ = 90° is
    broadside (maximum).

    Parameters
    ----------
    theta_deg : float or ndarray
        Angle(s) from the dipole axis [degrees], range [0, 180].
    l : float
        SDAS (dipole) physical length [m].
    k : float
        Wavenumber [m⁻¹].  Use :func:`wavenumber` to compute for given
        frequency and resistivity.
    normalize : bool
        Normalize the peak to 1.0 (default *True*).

    Returns
    -------
    F : ndarray
        |F(θ)| pattern values (≥ 0).

    Notes
    -----
    F(θ) = |[cos(kl cosθ/2) − cos(kl/2)]| / |sinθ|.
    The singularity at θ = 0° and 180° resolves to zero by L'Hôpital's rule.
    """
    theta = np.deg2rad(np.asarray(theta_deg, dtype=float))
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    num = np.cos(k * l * cos_t / 2.0) - np.cos(k * l / 2.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        F = np.where(np.abs(sin_t) < 1e-10, 0.0, np.abs(num / sin_t))

    if normalize:
        peak = F.max() if F.ndim > 0 else float(F)
        if peak > 0.0:
            F = F / peak

    return F


# ─────────────────────────────────────────────────────────────────────────────
# Array factor
# ─────────────────────────────────────────────────────────────────────────────


def array_factor(
    theta_b_deg: float | np.ndarray,
    N: int,
    d: float,
    k: float,
    beta: float = 0.0,
) -> np.ndarray:
    """
    Normalised array factor AF_n for an N-element linear PAS (eq. 19).

    The angle θ_b is measured FROM BROADSIDE (perpendicular to the array
    axis).  θ_b = 0° is the maximum direction when β = 0; θ_b = ±90° is
    along the array (end-fire direction).

    Parameters
    ----------
    theta_b_deg : float or ndarray
        Broadside angle(s) [degrees], range [−90, 90].
    N : int
        Number of SDAS elements.
    d : float
        Element-to-element spacing [m].
    k : float
        Wavenumber [m⁻¹].
    beta : float
        Inter-element phase shift [rad].  β = 0 → broadside array;
        use :func:`beam_steer` to compute β for a target angle.

    Returns
    -------
    AF : ndarray
        Normalised |AF_n(θ_b)| ∈ [0, 1].

    Notes
    -----
    AF_n = sin(N ψ/2) / [N sin(ψ/2)],  ψ = k d sinθ_b + β.
    """
    theta_b = np.deg2rad(np.asarray(theta_b_deg, dtype=float))
    psi = k * d * np.sin(theta_b) + beta
    half = psi / 2.0

    with np.errstate(divide="ignore", invalid="ignore"):
        AF = np.where(
            np.abs(half) < 1e-10,
            1.0,
            np.abs(np.sin(N * half)) / (N * np.abs(np.sin(half))),
        )

    return AF


# ─────────────────────────────────────────────────────────────────────────────
# Combined PAS pattern
# ─────────────────────────────────────────────────────────────────────────────


def pas_pattern(
    theta_b_deg: float | np.ndarray,
    N: int,
    d: float,
    k: float,
    beta: float = 0.0,
    l: float = 1000.0,
    *,
    normalize: bool = True,
) -> np.ndarray:
    """
    Total normalised far-field pattern of an N-element PAS.

    The combined pattern is the product of the SDAS element pattern and the
    array factor, evaluated at the same observation angle.

    Parameters
    ----------
    theta_b_deg : float or ndarray
        Broadside angle(s) [degrees], range [−90, 90].
    N : int
        Number of SDAS elements.
    d : float
        Element spacing [m].
    k : float
        Wavenumber [m⁻¹].
    beta : float
        Inter-element phase shift [rad].
    l : float
        SDAS length [m] (default 1000 m matching gxac023).
    normalize : bool
        Normalize peak to 1.0 (default *True*).

    Returns
    -------
    pattern : ndarray
        Combined |E_total(θ_b)| pattern (≥ 0).
    """
    theta_b_deg = np.asarray(theta_b_deg, dtype=float)

    # Element pattern: θ from y-axis = 90° − broadside angle
    theta_dipole_deg = 90.0 - np.abs(theta_b_deg)  # [0°, 90°]
    F = sdas_element_pattern(theta_dipole_deg, l, k, normalize=False)
    AF = array_factor(theta_b_deg, N, d, k, beta)

    pattern = F * AF
    if normalize:
        peak = pattern.max() if pattern.ndim > 0 else float(pattern)
        if peak > 0.0:
            pattern = pattern / peak

    return pattern


# ─────────────────────────────────────────────────────────────────────────────
# Beam-steering design
# ─────────────────────────────────────────────────────────────────────────────


def beam_steer(
    theta_m_deg: float,
    d: float,
    k: float,
) -> float:
    """
    Inter-element phase shift β [rad] to steer the main lobe to θ_m (eq. 23).

    Parameters
    ----------
    theta_m_deg : float
        Target main-lobe broadside angle [degrees].
    d : float
        Element spacing [m].
    k : float
        Wavenumber [m⁻¹].

    Returns
    -------
    beta : float
        Required phase shift [rad].  Apply the same β to each SDAS via the
        feed-network inductance delay.

    Notes
    -----
    Condition (eq. 23): β = −k d sinθ_m.
    """
    return float(-k * d * np.sin(np.deg2rad(theta_m_deg)))


def steering_angles(
    N: int,
    d: float,
    k: float,
    beta: float,
    *,
    n_range: int = 3,
) -> np.ndarray:
    """
    All main-lobe broadside angles [degrees] for the given PAS configuration.

    Solves k d sinθ + β = ±2nπ (eq. 21) for n = 0, ±1, ±2, ...

    Parameters
    ----------
    N, d, k : int / float
        Array parameters (as in :func:`array_factor`).
    beta : float
        Inter-element phase shift [rad].
    n_range : int
        Search over n = −n_range … +n_range (default 3).

    Returns
    -------
    angles : ndarray
        Sorted array of main-lobe broadside angles [degrees] inside [−90°, 90°].
    """
    angles = []
    for n in range(-n_range, n_range + 1):
        val = (-beta + 2.0 * np.pi * n) / (k * d)
        if abs(val) <= 1.0:
            angles.append(float(np.rad2deg(np.arcsin(val))))
    return np.sort(np.unique(np.round(angles, 6)))


# ─────────────────────────────────────────────────────────────────────────────
# Directivity and gain
# ─────────────────────────────────────────────────────────────────────────────


def sdas_directivity(
    l: float,
    k: float,
    *,
    n_theta: int = 2000,
) -> float:
    """
    2-D horizontal-plane directivity D₀ = 2π U_max / ∫ U(θ) dθ (eq. 12).

    Parameters
    ----------
    l : float
        SDAS length [m].
    k : float
        Wavenumber [m⁻¹].
    n_theta : int
        Number of angular samples for numerical integration.

    Returns
    -------
    D0 : float
        Directivity (dimensionless).  A perfect omnidirectional source has
        D₀ = 1 in 2-D.
    """
    # Integrate over full circle, using element pattern vs dipole axis angle
    theta = np.linspace(1e-6, np.pi - 1e-6, n_theta)
    F = sdas_element_pattern(np.rad2deg(theta), l, k, normalize=False)
    U = F**2
    P_rad = float(_trapz(U * np.sin(theta), theta) * 2.0 * np.pi)
    U_max = float(U.max())
    if P_rad < 1e-100:
        return np.nan
    return float(4.0 * np.pi * U_max / P_rad)


def snr_gain_db(N: int) -> float:
    """
    SNR improvement of an N-element PAS relative to a single SDAS [dB].

    For coherent beam forming, the gain scales as N²:
    G_PAS / G_SDAS = N²  →  10 log₁₀(N²) = 20 log₁₀(N) dB.

    Parameters
    ----------
    N : int
        Number of SDAS elements.

    Returns
    -------
    gain_dB : float
        SNR gain [dB].
    """
    return float(20.0 * np.log10(max(N, 1)))


# ─────────────────────────────────────────────────────────────────────────────
# Visualisation
# ─────────────────────────────────────────────────────────────────────────────


def plot_radiation_pattern(
    theta_b_deg: np.ndarray | Sequence,
    patterns: np.ndarray | Sequence,
    *,
    labels: Sequence[str] | None = None,
    polar: bool = True,
    normalize: bool = True,
    log_scale: bool = False,
    db_floor: float = -40.0,
    title: str = "Radiation pattern",
    figsize: tuple[float, float] = (7.0, 7.0),
    ax=None,
):
    """
    Plot one or more radiation patterns in polar or Cartesian format.

    Parameters
    ----------
    theta_b_deg : array-like
        Broadside angles [degrees], range [−90, 90].
    patterns : array-like or list of array-like
        Pattern amplitude(s).  A 2-D array is treated as multiple patterns
        with shape (n_patterns, n_angles).
    labels : list of str or None
        Legend labels for each pattern.
    polar : bool
        Polar (default) or Cartesian plot.
    normalize : bool
        Normalise each pattern to its peak before plotting.
    log_scale : bool
        Convert to dB (20 log₁₀) for the radial / y-axis.
    db_floor : float
        Minimum dB value when ``log_scale=True`` (default −40 dB).
    title : str
        Axes title.
    figsize : tuple
        Figure size if a new figure is created.
    ax : matplotlib.axes.Axes or None
        Axes to draw on; created if *None*.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    theta_b = np.asarray(theta_b_deg, dtype=float)
    pats = np.atleast_2d(np.asarray(patterns, dtype=float))

    if normalize:
        peak = pats.max(axis=1, keepdims=True)
        peak[peak == 0] = 1.0
        pats = pats / peak

    if log_scale:
        with np.errstate(divide="ignore"):
            pats = 20.0 * np.log10(np.maximum(pats, 1e-10))
        pats = np.maximum(pats, db_floor) - db_floor  # shift to ≥ 0

    n_pat = pats.shape[0]
    lbls = labels if labels is not None else [f"pattern {i + 1}" for i in range(n_pat)]

    if ax is None:
        if polar:
            _, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=figsize)
        else:
            _, ax = plt.subplots(figsize=figsize)

    for i in range(n_pat):
        if polar:
            # Polar: angle 0 at top (N), clockwise → use θ measured from top
            # Map broadside (θ_b) to polar angle: θ_b=0 → π/2 (right)
            phi = np.deg2rad(90.0 - theta_b)
            ax.plot(phi, pats[i], label=lbls[i])
        else:
            ax.plot(theta_b, pats[i], label=lbls[i])

    if polar:
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        hide_polar_radius_labels(ax)
        ax.set_title(title, pad=15)
    else:
        ax.set_xlabel("Broadside angle θ [deg]")
        y_label = "Amplitude [dB]" if log_scale else "Normalised amplitude"
        ax.set_ylabel(y_label)
        ax.set_title(title)
        ax.grid(True, linestyle=":")

    if n_pat > 1 or labels is not None:
        ax.legend(fontsize=8)

    return ax
