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
    Central-loop time-domain EM (step-off waveform).  Delegates the
    Hankel/Fourier transform to :mod:`empymod`'s validated digital
    linear filters rather than a hand-rolled quadrature.

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
MU0: float = 4.0e-7 * np.pi  # H/m — magnetic permeability of free space


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

    def to_array(
        self, *, log_rho: bool = True, include_phase: bool = True
    ) -> np.ndarray:
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
            ra = (
                np.log10(np.maximum(self.rho_a, 1e-12))
                if log_rho
                else self.rho_a
            )
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
    Z = 1j * omega * MU0 / k_n  # = sqrt(iωμρ_N)

    # Propagate upward
    for j in range(len(rho) - 2, -1, -1):
        k_j = np.sqrt(1j * omega * MU0 / rho[j])
        Z0_j = 1j * omega * MU0 / k_j  # intrinsic impedance of layer j
        th = np.tanh(k_j * thick[j])
        Z = Z0_j * (Z + Z0_j * th) / (Z0_j + Z * th)

    return Z


# TEM1DForward's frequency-domain Hankel transform and time-domain Fourier
# transform are delegated to empymod (Werthmüller, 2017) rather than
# hand-rolled here -- see TEM1DForward's docstring for why.


# ─────────────────────────────────────────────────────────────────────────────
# Abstract base solver
# ─────────────────────────────────────────────────────────────────────────────


class _Base1DForward(ABC):
    @abstractmethod
    def run(self, model) -> ForwardResponse:  # model: LayeredModel
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
    >>> model = LayeredModel(resistivity=[100, 10, 500], thickness=[500, 1000])
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

    Computes the vertical magnetic field via :mod:`empymod`'s validated
    digital-linear-filter Hankel and Fourier transforms
    (Werthmüller, 2017), rather than a hand-rolled quadrature.

    .. note::
        An earlier from-scratch implementation (fixed-order Gauss-Legendre
        Hankel transform + naive trapezoidal cosine transform) did not
        converge at realistic time ranges -- its own frequency-domain
        kernel grew unboundedly at high frequency instead of decaying, and
        refining either quadrature's resolution changed the answer's
        *sign*, not just its accuracy, at several requested times. That
        implementation is gone; :mod:`empymod` is a peer-reviewed,
        independently validated EM modelling library built specifically
        for this class of problem, and is verified byte-for-byte against
        Ward & Hohmann (1988)'s own closed-form half-space dBz/dt formula
        (their eq. 4.70) as part of this module's test suite.

    The loop itself is represented as a small tangential electric bipole
    at the loop radius, scaled by the loop's circumference -- the same
    "loop as a scaled dipole segment" construction :mod:`empymod` uses to
    reproduce Ward & Hohmann's own central-loop figures (4.7-4.8), valid
    by the axisymmetry of a horizontal circular loop.

    Parameters
    ----------
    times : array-like
        Measurement times [s] for the step-off response.
        Typical range: 1e-6 – 1e-2 s.
    loop_radius : float
        Transmitter loop radius [m].  Default 50 m.
    moment : float
        Transmitter magnetic moment [A·m²].  Default 1 A·m².
    n_freqs, n_lam : int
        Unused. Accepted only so existing callers (e.g.
        :class:`pycsamt.inversion.backends.builtin.Builtin1DBackend`,
        which threads ``backend_options`` straight through) do not break;
        :mod:`empymod`'s own filters pick their resolution automatically.

    References
    ----------
    Ward & Hohmann (1988), *Electromagnetic Methods in Applied
    Geophysics*, Vol. 1.

    Werthmüller, D. (2017). An open-source full 3D electromagnetic
    modeler for 1D VTI media in Python: empymod. *Geophysics*, 82(6),
    WB9-WB19.
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
        # Frequency grid for the (informational) hz_freq field only --
        # not used by the dBz/dt computation itself.
        t_min = self.times.min()
        t_max = self.times.max()
        self._omega = np.logspace(
            np.log10(2.0 * np.pi / (10.0 * t_max)),
            np.log10(2.0 * np.pi / (t_min / 10.0)),
            n_freqs,
        )

    def _empymod_kwargs(self, model) -> dict:
        rho = np.asarray(model.resistivity, dtype=float)
        thick = np.asarray(model.thickness, dtype=float)
        a = self.loop_radius

        res = np.concatenate(([2e14], rho))  # air, then earth layers
        depth = np.concatenate(([0.0], np.cumsum(thick)))

        # A loop driven by current I has magnetic moment m = I * area;
        # empymod's loop-as-bipole trick wants strength = I * circumference
        # = (m / area) * (2*pi*a) = 2*m/a.
        strength = 2.0 * self.moment / a

        return dict(
            src=[a, 0.0, 0.0, 90.0, 0.0],
            rec=[0.0, 0.0, 0.0, 0.0, 90.0],
            depth=depth.tolist(),
            res=res.tolist(),
            strength=strength,
            mrec=True,
            epermH=[0.0] * res.size,  # reduces early-time numerical noise
            verb=0,
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
        import empymod

        kwargs = self._empymod_kwargs(model)

        # signal=0 (impulse response) is dBz/dt for a step-off source --
        # standard LTI result: d/dt of the switch-off step equals minus
        # the impulse response, and empymod returns that derivative
        # directly (verified against Ward & Hohmann 1988 eq. 4.70).
        dBz_dt = np.asarray(
            empymod.bipole(freqtime=self.times, signal=0, **kwargs).real,
            dtype=float,
        )

        hz_fd = np.asarray(
            empymod.bipole(
                freqtime=self._omega / (2.0 * np.pi), signal=None, **kwargs
            ),
            dtype=complex,
        )

        return ForwardResponse(
            method="TEM1D",
            times=self.times,
            dBz_dt=dBz_dt,
            hz_freq=hz_fd,
            model=model,
        )


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
