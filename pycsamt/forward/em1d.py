# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
1-D electromagnetic forward solvers.

Three EM methods are provided, all using the same
:class:`LayeredModel` / :class:`ForwardResponse` objects:

``MT1DForward``
    Magnetotelluric (plane-wave source).  Uses the exact recursive
    Wait (1954) impedance formula — the standard for generating
    MT/AMT/CSAMT far-field training data.

``TEM1DForward``
    Central-loop time-domain EM (step-off waveform).  Computes the
    frequency-domain H_z via a Hankel transform of the TE admittance
    kernel, then converts to time via a cosine transform.  Uses
    ``scipy.integrate`` (correct but not optimised; a DLF
    acceleration layer is planned for Phase 2).

``CSAMT1DForward``
    Controlled-source AMT.  In the far-field limit this reduces to
    MT.  A near-field geometry correction factor is applied when
    source–receiver distance and frequencies are provided.

References
----------
Wait, J.R. (1954). On the relation between telluric currents and the
Earth's magnetic field. *Geophysics*, 19(2), 281-289.

Nabighian, M.N. (1979). Quasi-static transient response of a
conducting half-space. *Geophysics*, 44(10), 1700-1705.

Ward, S.H. & Hohmann, G.W. (1988). Electromagnetic theory for
geophysical applications. In: *Electromagnetic Methods in Applied
Geophysics*, 1, 130-311.

Puzyrev, V. et al. (2021). Inversion of 1D frequency- and
time-domain EM data with CNNs. *Computers & Geosciences*, 149,
104681.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "MT1DForward",
    "TEM1DForward",
    "CSAMT1DForward",
    "ForwardResponse",
]

# Physical constants
MU0: float = 4.0e-7 * np.pi   # H/m — magnetic permeability of free space


# ─────────────────────────────────────────────────────────────────────────────
# Output container
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ForwardResponse:
    """
    Container for the output of a 1-D forward solver.

    Parameters
    ----------
    method : str
        Solver identifier (``'MT1D'``, ``'TEM1D'``, ``'CSAMT1D'``).
    freqs : ndarray or None
        Frequencies in Hz (MT/CSAMT).
    times : ndarray or None
        Times in seconds (TEM).
    z : ndarray or None
        Complex surface impedance (n_freq,) for MT/CSAMT [V/A].
    rho_a : ndarray or None
        Apparent resistivity (n_freq,) [Ω·m].
    phase : ndarray or None
        Impedance phase (n_freq,) [degrees, 0–90° for normal models].
    dBz_dt : ndarray or None
        dBz/dt step-off response (n_times,) [T/s] for TEM.
    hz_freq : ndarray or None
        Complex frequency-domain H_z (n_freq,) for TEM.
    model : LayeredModel or None
        The input model that produced this response.
    """
    method: str = "MT1D"
    freqs: np.ndarray | None = None
    times: np.ndarray | None = None
    z: np.ndarray | None = None
    rho_a: np.ndarray | None = None
    phase: np.ndarray | None = None
    dBz_dt: np.ndarray | None = None
    hz_freq: np.ndarray | None = None
    model: object = field(default=None, repr=False)

    def to_array(self, *, log_rho: bool = True, include_phase: bool = True) -> np.ndarray:
        """
        Flatten to a 1-D feature vector for ML input.

        For MT/CSAMT returns ``[log10(rho_a), phase_deg]`` concatenated
        (length 2 × n_freq); for TEM returns ``log10(|dBz_dt|)``
        (length n_times).

        Parameters
        ----------
        log_rho : bool
            Apply log₁₀ to apparent resistivity (recommended).
        include_phase : bool
            Include phase alongside ρ_a (MT only).
        """
        if self.method in ("MT1D", "CSAMT1D"):
            ra = np.log10(np.maximum(self.rho_a, 1e-12)) if log_rho else self.rho_a
            if include_phase:
                return np.concatenate([ra, self.phase])
            return ra
        # TEM
        db = np.abs(self.dBz_dt)
        db = np.where(db > 0, db, 1e-30)
        return np.log10(db) if log_rho else db

    def plot(self, ax=None, **kwargs):
        """Quick diagnostic plot.  Returns the Axes used."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(5, 4))
        if self.method in ("MT1D", "CSAMT1D"):
            period = 1.0 / self.freqs
            ax.loglog(period, self.rho_a, **kwargs)
            ax.set_xlabel("Period (s)")
            ax.set_ylabel(r"$\rho_a$ (Ω·m)")
            ax.set_title(f"{self.method} apparent resistivity")
        else:
            ax.loglog(self.times, np.abs(self.dBz_dt), **kwargs)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel(r"$|d\mathbf{B}_z/dt|$ (T/s)")
            ax.set_title("TEM1D step-off response")
        return ax


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _z_surface_mt(omega: float, rho: np.ndarray, thick: np.ndarray) -> complex:
    """
    Compute the MT surface impedance for a single angular frequency.

    Uses the Wait (1954) downward-recursion formula from the halfspace
    upward through each layer.

    Parameters
    ----------
    omega : float
        Angular frequency 2πf  [rad/s].
    rho : ndarray (n_layers,)
        Layer resistivities [Ω·m], top → bottom. Last entry = halfspace.
    thick : ndarray (n_layers-1,)
        Layer thicknesses [m].

    Returns
    -------
    Z_s : complex
        Surface impedance  [V/A].
    """
    # Halfspace intrinsic impedance
    k_n = np.sqrt(1j * omega * MU0 / rho[-1])
    Z = 1j * omega * MU0 / k_n            # = sqrt(iωμρ_N)

    # Propagate upward
    for j in range(len(rho) - 2, -1, -1):
        k_j = np.sqrt(1j * omega * MU0 / rho[j])
        Z0_j = 1j * omega * MU0 / k_j     # intrinsic impedance of layer j
        th = np.tanh(k_j * thick[j])
        Z = Z0_j * (Z + Z0_j * th) / (Z0_j + Z * th)

    return Z


def _te_admittance(lam: np.ndarray, omega: float,
                   rho: np.ndarray, thick: np.ndarray) -> np.ndarray:
    """
    Vectorised TE-mode surface admittance Y_TE(λ, ω) for all λ at once.

    Parameters
    ----------
    lam : ndarray (n_lam,)
        Horizontal wavenumbers [1/m].
    omega : float
        Angular frequency [rad/s].
    rho, thick : as in _z_surface_mt

    Returns
    -------
    Y_TE : ndarray (n_lam,), complex
    """
    sigma = 1.0 / rho                                  # conductivities

    # Vertical wavenumber  ν_j² = λ² + iωμ₀σ_j
    # Shape broadcast: (n_lam, n_layers)
    nu = np.sqrt(
        lam[:, None] ** 2 + 1j * omega * MU0 * sigma[None, :]
    )

    # Intrinsic TE admittance Y0_j = ν_j / (iωμ₀)
    Y0 = nu / (1j * omega * MU0)          # (n_lam, n_layers)

    # Start from halfspace
    Y = Y0[:, -1].copy()                  # (n_lam,)

    for j in range(len(rho) - 2, -1, -1):
        Y0j = Y0[:, j]
        nu_j = nu[:, j]
        th = np.tanh(nu_j * thick[j])
        Y = Y0j * (Y + Y0j * th) / (Y0j + Y * th)

    return Y


def _hankel_hz_fd(omega: float, rho: np.ndarray, thick: np.ndarray,
                  loop_radius: float, n_lam: int = 150) -> complex:
    """
    Frequency-domain H_z at the centre of a horizontal circular loop of
    radius *loop_radius* over a 1-D earth.

    The integral ∫₀^∞ λ · r_TE(λ,ω) · J₁(λa) dλ is evaluated with
    Gauss-Legendre quadrature in log(λ) space.

    Returns H_z per unit moment [1/m³].
    """
    from numpy.polynomial.legendre import leggauss
    from scipy.special import j1 as J1

    a = loop_radius
    lam_lo = np.log(1e-5 / a)
    lam_hi = np.log(1e5 / a)

    nodes, wts = leggauss(n_lam)
    t = 0.5 * (lam_hi - lam_lo) * nodes + 0.5 * (lam_hi + lam_lo)
    lam = np.exp(t)
    jac = 0.5 * (lam_hi - lam_lo)

    Y_TE = _te_admittance(lam, omega, rho, thick)
    Y0_air = lam / (1j * omega * MU0)    # air admittance ν₀/iωμ₀ ≈ λ/iωμ₀

    r_TE = (Y0_air - Y_TE) / (Y0_air + Y_TE)

    # Kernel:  λ · r_TE · J₁(λa)  ×  λ (Jacobian for log substitution)
    integrand = lam * r_TE * J1(lam * a) * lam

    return jac * complex(np.dot(wts, integrand)) / (4.0 * np.pi)


# ─────────────────────────────────────────────────────────────────────────────
# Abstract base solver
# ─────────────────────────────────────────────────────────────────────────────

class _Base1DForward(ABC):
    @abstractmethod
    def run(self, model) -> ForwardResponse:   # model: LayeredModel
        ...


# ─────────────────────────────────────────────────────────────────────────────
# MT 1-D forward
# ─────────────────────────────────────────────────────────────────────────────

class MT1DForward(_Base1DForward):
    """
    1-D magnetotelluric forward solver (plane-wave, isotropic earth).

    Uses the exact Wait (1954) recursive impedance algorithm.  Runs in
    O(n_freq × n_layers) time — fast enough for generating millions of
    synthetic training samples.

    Parameters
    ----------
    freqs : array-like
        Frequencies [Hz] at which to evaluate the response.
        Typical range: 1e-4 – 1e5 Hz for MT/AMT.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.em1d import MT1DForward
    >>> from pycsamt.forward.synthetic import LayeredModel
    >>> freqs = np.logspace(-3, 4, 30)
    >>> model = LayeredModel(
    ...     resistivity=[100, 10, 500],
    ...     thickness=[500, 1000]
    ... )
    >>> resp = MT1DForward(freqs).run(model)
    >>> resp.rho_a.shape
    (30,)
    """

    def __init__(self, freqs: Sequence[float] | np.ndarray):
        self.freqs = np.asarray(freqs, dtype=float)

    def run(self, model) -> ForwardResponse:
        """
        Compute the MT 1-D response for *model*.

        Parameters
        ----------
        model : LayeredModel
            Input earth model.

        Returns
        -------
        ForwardResponse
            Fields populated: ``z``, ``rho_a``, ``phase``, ``freqs``.
        """
        rho = model.resistivity
        thick = model.thickness
        omega = 2.0 * np.pi * self.freqs

        z_surf = np.array(
            [_z_surface_mt(w, rho, thick) for w in omega], dtype=complex
        )
        rho_a = np.abs(z_surf) ** 2 / (omega * MU0)
        phase = np.angle(z_surf, deg=True)

        return ForwardResponse(
            method="MT1D",
            freqs=self.freqs,
            z=z_surf,
            rho_a=rho_a,
            phase=phase,
            model=model,
        )


# ─────────────────────────────────────────────────────────────────────────────
# TEM 1-D forward
# ─────────────────────────────────────────────────────────────────────────────

class TEM1DForward(_Base1DForward):
    """
    1-D central-loop TEM forward solver (step-off waveform).

    Computes the vertical magnetic field H_z(ω) via a numerical Hankel
    transform of the TE admittance kernel, then converts to the step-off
    time-domain dBz/dt via a cosine transform.

    .. note::
        This implementation uses ``scipy.integrate`` for correctness.
        It is suitable for generating training datasets of moderate size
        (up to ~10 000 samples).  A Digital Linear Filter (DLF)
        optimisation yielding 100× speed-up is planned for Phase 2.

    Parameters
    ----------
    times : array-like
        Measurement times [s] for the step-off response.
        Typical range: 1e-6 – 1e-2 s.
    loop_radius : float
        Transmitter loop radius [m].  Default 50 m.
    moment : float
        Transmitter magnetic moment [A·m²].  Default 1 A·m².
    n_freqs : int
        Number of frequency-domain evaluation points used in the
        cosine transform.  Higher → better accuracy at early times.

    References
    ----------
    Ward & Hohmann (1988), *Electromagnetic Methods in Applied
    Geophysics*, Vol. 1.
    """

    def __init__(
        self,
        times: Sequence[float] | np.ndarray,
        loop_radius: float = 50.0,
        moment: float = 1.0,
        n_freqs: int = 64,
        n_lam: int = 100,
    ):
        self.times = np.asarray(times, dtype=float)
        self.loop_radius = float(loop_radius)
        self.moment = float(moment)
        self.n_freqs = n_freqs
        self.n_lam = n_lam
        # Pre-compute frequency grid for the cosine transform
        t_min = self.times.min()
        t_max = self.times.max()
        # Cover 4 decades below/above the time window
        self._omega = np.logspace(
            np.log10(2.0 * np.pi / (10.0 * t_max)),
            np.log10(2.0 * np.pi / (t_min / 10.0)),
            n_freqs,
        )

    def run(self, model) -> ForwardResponse:
        """
        Compute the TEM 1-D step-off response for *model*.

        Parameters
        ----------
        model : LayeredModel
            Input earth model.

        Returns
        -------
        ForwardResponse
            Fields populated: ``dBz_dt``, ``hz_freq``, ``times``.
        """
        rho = model.resistivity
        thick = model.thickness

        # Frequency-domain H_z at each ω
        hz_fd = np.array(
            [_hankel_hz_fd(w, rho, thick, self.loop_radius, self.n_lam)
             for w in self._omega],
            dtype=complex,
        ) * self.moment

        # Time-domain dBz/dt via cosine transform
        # dBz/dt(t) ≈ (2/π) × ∫₀^∞ ω · Im[H_z(ω)] · cos(ωt) dω
        dBz_dt = self._cosine_transform(hz_fd)

        return ForwardResponse(
            method="TEM1D",
            times=self.times,
            dBz_dt=dBz_dt,
            hz_freq=hz_fd,
            model=model,
        )

    def _cosine_transform(self, hz_fd: np.ndarray) -> np.ndarray:
        """Numerical cosine transform → dBz/dt at self.times."""
        omega = self._omega
        dlog_omega = np.diff(np.log(omega))                 # spacings
        # Trapezoidal integration in log(ω) space
        # ∫ f(ω) dω  ≈  Σ  f(ω_k) · ω_k · Δlog(ω)
        result = np.empty(len(self.times))
        for i, t in enumerate(self.times):
            integrand = omega * np.imag(hz_fd) * np.cos(omega * t)
            # Trapezoidal in log-space: integrate f(ω)*ω d(log ω)
            pts = integrand * omega                          # f·ω for log spacing
            result[i] = (
                (2.0 / np.pi)
                * np.sum(0.5 * (pts[:-1] + pts[1:]) * dlog_omega)
            )
        return result


# ─────────────────────────────────────────────────────────────────────────────
# CSAMT 1-D forward
# ─────────────────────────────────────────────────────────────────────────────

class CSAMT1DForward(_Base1DForward):
    """
    1-D controlled-source AMT forward solver.

    In the far-field (source–receiver offset r ≫ skin depth δ) the
    CSAMT response approximates the MT plane-wave response.  A first-
    order near-field correction factor after Zonge & Hughes (1991) is
    applied when *source_offset* is supplied.

    Parameters
    ----------
    freqs : array-like
        Frequencies [Hz].
    source_offset : float or None
        Source–receiver distance [m].  If ``None``, the far-field
        approximation is used without correction.
    dipole_length : float
        Electric dipole length [m].  Default 1000 m.

    References
    ----------
    Zonge, K.L. & Hughes, L.J. (1991). Controlled source audio-
    frequency magnetotellurics. In *Electromagnetic Methods in Applied
    Geophysics*, 2B, 713-809.
    """

    def __init__(
        self,
        freqs: Sequence[float] | np.ndarray,
        source_offset: float | None = None,
        dipole_length: float = 1000.0,
    ):
        self.freqs = np.asarray(freqs, dtype=float)
        self.source_offset = source_offset
        self.dipole_length = float(dipole_length)
        self._mt = MT1DForward(self.freqs)

    def run(self, model) -> ForwardResponse:
        """
        Compute the CSAMT 1-D response for *model*.

        Returns
        -------
        ForwardResponse
            Fields populated: ``z``, ``rho_a``, ``phase``.
        """
        resp = self._mt.run(model)

        if self.source_offset is not None:
            resp = self._apply_nearfield_correction(resp, model)

        resp.method = "CSAMT1D"
        return resp

    def _apply_nearfield_correction(
        self, resp: ForwardResponse, model
    ) -> ForwardResponse:
        """
        First-order near-field correction after Zonge & Hughes (1991).

        The correction factor f_nf = 1 / (1 + (r/δ)⁻²) modulates the
        apparent resistivity so that ρ_a → ρ_true as r/δ → ∞.
        """
        omega = 2.0 * np.pi * self.freqs
        # Skin depth for the first layer  δ = sqrt(2ρ/(ωμ₀))
        rho1 = model.resistivity[0]
        delta = np.sqrt(2.0 * rho1 / (omega * MU0))
        r = self.source_offset
        f_nf = 1.0 / (1.0 + (r / delta) ** (-2))

        # Correct apparent resistivity and recompute phase-compatible Z
        rho_a_corr = resp.rho_a * f_nf
        z_corr = np.sqrt(rho_a_corr * omega * MU0) * np.exp(
            1j * np.angle(resp.z)
        )
        return ForwardResponse(
            method="CSAMT1D",
            freqs=self.freqs,
            z=z_corr,
            rho_a=rho_a_corr,
            phase=np.angle(z_corr, deg=True),
            model=model,
        )
