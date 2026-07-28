# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Noise models for synthetic EM training data.

Noise injection turns clean forward responses into realistic
training samples that improve neural-network generalisation.
Two families are provided:

``GaussianNoise``
    Simple additive relative Gaussian noise, identical at all
    frequencies/times.  Fast and widely used in the literature
    (e.g. Puzyrev 2021 applies 1–2 % relative noise).

``FieldRealisticNoise``
    Frequency-dependent noise model that mimics real-world MT/CSAMT
    data quality:

    * Cultural EM noise peaks at 50/60 Hz power-line frequencies.
    * Dead-band effect (elevated noise at ~1 000–3 000 s period) due to
      reduced natural signal in the MT spectrum.
    * Early-time TEM noise from transmitter turn-off and system
      bandwidth.

``MultiplicativeNoise``
    Log-space additive noise — useful when the dynamic range of the
    data spans several decades (log₁₀ space perturbation).

All noise classes follow the same interface::

    noisy_resp = noise_model.apply(clean_resp)

or the functional helpers ``add_gaussian_noise`` / ``add_noise`` can be
used directly.

References
----------
Puzyrev, V. & Swidinsky, A. (2021). Inversion of 1D frequency- and
time-domain EM data with CNNs. *Computers & Geosciences*, 149, 104681.

Egbert, G.D. (1997). Robust multiple station magnetotelluric data
processing. *Geophysical Journal International*, 130, 475–496.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from .em1d import ForwardResponse

__all__ = [
    "GaussianNoise",
    "FieldRealisticNoise",
    "MultiplicativeNoise",
    "add_gaussian_noise",
    "add_noise",
]


# ─────────────────────────────────────────────────────────────────────────────
# Abstract base
# ─────────────────────────────────────────────────────────────────────────────


class _BaseNoiseModel(ABC):
    """Abstract base for all noise models."""

    @abstractmethod
    def apply(
        self,
        response: ForwardResponse,
        *,
        seed: int | None = None,
    ) -> ForwardResponse:
        """
        Return a new :class:`~pycsamt.forward.em1d.ForwardResponse`
        with noise added.

        Parameters
        ----------
        response : ForwardResponse
            Clean forward response.
        seed : int or None
            Random seed for reproducibility.
        """

    def __call__(self, response: ForwardResponse, **kwargs) -> ForwardResponse:
        return self.apply(response, **kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# Gaussian noise
# ─────────────────────────────────────────────────────────────────────────────


class GaussianNoise(_BaseNoiseModel):
    """
    Add relative Gaussian noise to an EM forward response.

    For MT/CSAMT the noise is applied to log₁₀(ρ_a) and the phase
    independently.  For TEM it is applied to log₁₀(|dBz/dt|).

    Parameters
    ----------
    level : float
        Relative noise standard deviation (e.g. 0.05 = 5 %).
    apply_to : {'rho_phase', 'z', 'both', 'rho', 'phase'}
        Which quantities to perturb.  ``'rho_phase'``, ``'both'`` and
        ``'z'`` are equivalent: they jointly perturb ρₐ and phase, which
        together determine the perturbed complex impedance Z.  ``'rho'``
        and ``'phase'`` perturb only one of the two. Default
        ``'rho_phase'``.
    phase_level : float or None
        Separate noise level for phase in degrees.  If ``None``,
        *level* is used scaled to the typical 45° phase range.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.em1d import MT1DForward
    >>> from pycsamt.forward.synthetic import LayeredModel
    >>> from pycsamt.forward.noise import GaussianNoise
    >>> m = LayeredModel([100, 10, 500], [300, 800])
    >>> resp = MT1DForward(np.logspace(-3, 4, 20)).run(m)
    >>> noisy = GaussianNoise(level=0.05).apply(resp, seed=0)
    >>> noisy.rho_a.shape
    (20,)
    """

    def __init__(
        self,
        level: float = 0.05,
        apply_to: str = "rho_phase",
        phase_level: float | None = None,
    ):
        if not 0.0 < level < 1.0:
            raise ValueError(f"level must be in (0, 1), got {level}")
        self.level = level
        self.apply_to = apply_to
        self.phase_level = phase_level if phase_level is not None else level * 45.0

    def apply(
        self,
        response: ForwardResponse,
        *,
        seed: int | None = None,
    ) -> ForwardResponse:
        rng = np.random.default_rng(seed)

        if response.method in ("MT1D", "CSAMT1D"):
            return self._apply_mt(response, rng)
        return self._apply_tem(response, rng)

    def _apply_mt(self, resp: ForwardResponse, rng) -> ForwardResponse:
        rho_a = resp.rho_a.copy()
        phase = resp.phase.copy()

        if self.apply_to in ("rho_phase", "both", "z", "rho"):
            # Noise in log₁₀ space then exponentiate
            log_rho = np.log10(rho_a)
            log_rho += rng.normal(0.0, self.level, rho_a.shape)
            rho_a = 10.0**log_rho

        if self.apply_to in ("rho_phase", "both", "z", "phase"):
            phase += rng.normal(0.0, self.phase_level, phase.shape)

        # Recompute Z consistent with perturbed rho_a and phase
        omega = 2.0 * np.pi * resp.freqs
        from .em1d import MU0

        z = np.sqrt(rho_a * omega * MU0) * np.exp(1j * np.deg2rad(phase))

        import copy

        out = copy.copy(resp)
        out.rho_a = rho_a
        out.phase = phase
        out.z = z
        return out

    def _apply_tem(self, resp: ForwardResponse, rng) -> ForwardResponse:
        db = resp.dBz_dt.copy()
        log_db = np.log10(np.maximum(np.abs(db), 1e-30))
        log_db += rng.normal(0.0, self.level, db.shape)

        import copy

        out = copy.copy(resp)
        out.dBz_dt = np.sign(db) * 10.0**log_db
        return out


# ─────────────────────────────────────────────────────────────────────────────
# Multiplicative log-space noise
# ─────────────────────────────────────────────────────────────────────────────


class MultiplicativeNoise(_BaseNoiseModel):
    """
    Log-space (multiplicative) Gaussian noise.

    Equivalent to a log-normal perturbation on the data value.
    Appropriate when the data spans many orders of magnitude.

    Parameters
    ----------
    sigma_log10 : float
        Standard deviation of the log₁₀ perturbation.
        ``sigma_log10=0.05`` corresponds to ≈ ±12 % variation.
    """

    def __init__(self, sigma_log10: float = 0.05):
        self.sigma_log10 = sigma_log10

    def apply(
        self,
        response: ForwardResponse,
        *,
        seed: int | None = None,
    ) -> ForwardResponse:
        rng = np.random.default_rng(seed)
        import copy

        out = copy.copy(response)

        if response.rho_a is not None:
            log_rho = np.log10(np.maximum(response.rho_a, 1e-12))
            out.rho_a = 10.0 ** (
                log_rho + rng.normal(0.0, self.sigma_log10, log_rho.shape)
            )

        if response.dBz_dt is not None:
            db = response.dBz_dt
            log_db = np.log10(np.maximum(np.abs(db), 1e-30))
            out.dBz_dt = np.sign(db) * 10.0 ** (
                log_db + rng.normal(0.0, self.sigma_log10, db.shape)
            )
        return out


# ─────────────────────────────────────────────────────────────────────────────
# Field-realistic noise (frequency-dependent)
# ─────────────────────────────────────────────────────────────────────────────


class FieldRealisticNoise(_BaseNoiseModel):
    """
    Frequency-dependent noise model for MT/CSAMT training data.

    Simulates three dominant noise sources present in field data:

    1. **Power-line harmonics** — Gaussian spikes at 50 Hz, 100 Hz,
       150 Hz (or 60/120/180 Hz for North American data).
    2. **MT dead band** — elevated noise at ~1 000–3 000 s period (≈
       3×10⁻⁴ – 10⁻³ Hz) due to reduced natural source energy.
    3. **Background floor** — minimum noise level applied uniformly.

    Parameters
    ----------
    base_level : float
        Background relative noise floor.  Default 0.02 (2 %).
    powerline_freq : float
        Fundamental power-line frequency [Hz].  50 or 60.
    powerline_level : float
        Additional noise at power-line harmonics (relative).
    dead_band : bool
        If ``True``, inflate noise in the MT dead band.
    dead_band_freq_range : (low, high)
        Frequency range [Hz] of the dead band.  Default (3e-4, 1e-3).
    dead_band_level : float
        Additional noise level in the dead band.  Default 0.15.
    """

    def __init__(
        self,
        base_level: float = 0.02,
        powerline_freq: float = 50.0,
        powerline_level: float = 0.30,
        dead_band: bool = True,
        dead_band_freq_range: tuple = (3e-4, 1e-3),
        dead_band_level: float = 0.15,
    ):
        self.base_level = base_level
        self.powerline_freq = powerline_freq
        self.powerline_level = powerline_level
        self.dead_band = dead_band
        self.dead_band_freq_range = dead_band_freq_range
        self.dead_band_level = dead_band_level

    def noise_profile(self, freqs: np.ndarray) -> np.ndarray:
        """
        Return the frequency-dependent noise level array.

        Parameters
        ----------
        freqs : ndarray
            Frequencies [Hz].

        Returns
        -------
        sigma : ndarray, same shape as *freqs*
            Relative noise standard deviation at each frequency.
        """
        sigma = np.full(len(freqs), self.base_level)

        # Power-line harmonics (50, 100, 150 Hz etc.)
        for k in range(1, 5):
            f_harm = k * self.powerline_freq
            mask = np.abs(freqs - f_harm) / f_harm < 0.05  # within 5 %
            sigma[mask] = np.maximum(sigma[mask], self.powerline_level / k)

        # MT dead band
        if self.dead_band:
            f_lo, f_hi = self.dead_band_freq_range
            db_mask = (freqs >= f_lo) & (freqs <= f_hi)
            sigma[db_mask] = np.maximum(sigma[db_mask], self.dead_band_level)

        return sigma

    def apply(
        self,
        response: ForwardResponse,
        *,
        seed: int | None = None,
    ) -> ForwardResponse:
        if response.freqs is None:
            raise ValueError(
                "FieldRealisticNoise requires freqs; use GaussianNoise for TEM."
            )

        rng = np.random.default_rng(seed)
        sigma = self.noise_profile(response.freqs)

        import copy

        out = copy.copy(response)

        if response.rho_a is not None:
            log_rho = np.log10(np.maximum(response.rho_a, 1e-12))
            log_rho += rng.normal(0.0, 1.0, sigma.shape) * sigma
            out.rho_a = 10.0**log_rho

        if response.phase is not None:
            out.phase = (
                response.phase + rng.normal(0.0, 1.0, sigma.shape) * sigma * 45.0
            )

        if out.rho_a is not None and out.phase is not None:
            omega = 2.0 * np.pi * response.freqs
            from .em1d import MU0

            out.z = np.sqrt(out.rho_a * omega * MU0) * np.exp(
                1j * np.deg2rad(out.phase)
            )
        return out

    def plot_profile(self, freqs: np.ndarray, ax=None):
        """Visualise the noise level vs frequency.  Returns Axes."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(5, 3))
        sigma = self.noise_profile(freqs)
        ax.semilogx(freqs, sigma * 100.0, color="firebrick")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Noise level (%)")
        ax.set_title("Field-realistic noise profile")
        ax.grid(True, which="both", alpha=0.3)
        return ax


# ─────────────────────────────────────────────────────────────────────────────
# Functional helpers
# ─────────────────────────────────────────────────────────────────────────────


def add_gaussian_noise(
    response: ForwardResponse,
    level: float = 0.05,
    *,
    seed: int | None = None,
) -> ForwardResponse:
    """
    Convenience wrapper: add Gaussian noise at a given relative level.

    Parameters
    ----------
    response : ForwardResponse
    level : float
        Relative noise standard deviation.
    seed : int or None

    Returns
    -------
    ForwardResponse
        New response with noise applied.
    """
    return GaussianNoise(level).apply(response, seed=seed)


def add_noise(
    response: ForwardResponse,
    noise_model: _BaseNoiseModel | str = "gaussian",
    level: float = 0.05,
    *,
    seed: int | None = None,
    **kwargs,
) -> ForwardResponse:
    """
    Apply a named or pre-built noise model to a forward response.

    Parameters
    ----------
    response : ForwardResponse
    noise_model : str or noise instance
        ``'gaussian'``, ``'multiplicative'``, ``'field'``, or an
        already-instantiated noise object.
    level : float
        Base noise level (used when constructing named models).
    seed : int or None
    **kwargs
        Extra keyword arguments forwarded to the noise model constructor.

    Returns
    -------
    ForwardResponse
    """
    if isinstance(noise_model, str):
        nm = noise_model.lower()
        if nm in ("gaussian", "gauss"):
            model = GaussianNoise(level, **kwargs)
        elif nm in ("multiplicative", "mult", "log"):
            model = MultiplicativeNoise(level, **kwargs)
        elif nm in ("field", "realistic", "field_realistic"):
            model = FieldRealisticNoise(base_level=level, **kwargs)
        else:
            raise ValueError(
                f"Unknown noise_model {noise_model!r}. "
                "Use 'gaussian', 'multiplicative', or 'field'."
            )
    else:
        model = noise_model

    return model.apply(response, seed=seed)
