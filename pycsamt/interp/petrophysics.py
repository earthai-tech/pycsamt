# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Petrophysical transforms for EM hydrogeophysics.

Converts EM-derived resistivity (from AMT, MT, EMAP, TDEM inversions) into
quantitative hydrogeological parameters: saturation, porosity, hydraulic
conductivity, transmissivity, and water chemistry indicators.

All transforms are **bidirectional** — forward (hydro → ρ) and inverse (ρ →
hydro) — and accept numpy arrays so they can be applied cell-by-cell to a
full 2-D resistivity section from :class:`~pycsamt.interp.ResistivityModel`.

Background
----------
Two petrophysical frameworks are implemented:

* **Archie (1942)** — clean sandstones and crystalline rocks:

  .. math::

      \\rho = a \\cdot \\rho_w \\cdot \\phi^{-m} \\cdot S_w^{-n}

  where *a* is the tortuosity factor (≈ 1), *m* the cementation exponent
  (1.3–2.5), and *n* the saturation exponent (1.8–2.5).

* **Waxman-Smits (1968)** — clay-bearing sediments (EMAP, shallow AMT):

  .. math::

      \\sigma = \\frac{S_w^n}{F}(\\sigma_w + \\sigma_s / S_w)

  where σ_s is the surface (clay) conductivity.

Additional module-level functions cover:

* Kozeny-Carman hydraulic conductivity from porosity
* Transmissivity and storativity from layer geometry
* Water chemistry: ρ_w ↔ TDS ↔ EC
* Hashin-Shtrikman effective-medium bounds
* EM-specific helpers: Bostick depth (TDEM/AMT), skin depth, fracture-zone K

References
----------
.. [1] Archie, G. E. (1942). The electrical resistivity log as an aid in
   determining some reservoir characteristics. *Trans. AIME*, 146, 54–62.
.. [2] Waxman, M. H. & Smits, L. J. M. (1968). Electrical conductivities in
   oil-bearing shaly sands. *Soc. Pet. Eng. J.*, 8, 107–122.
.. [3] Keller, G. V. (1966). Electrical methods in geophysical prospecting.
   *Pergamon Press*, New York.
.. [4] Kozeny, J. (1927). Über kapillare Leitung des Wassers im Boden.
   *Sitzungsber. Akad. Wiss. Wien*, 136(2a), 271–306.
.. [5] Hashin, Z. & Shtrikman, S. (1962). A variational approach to the
   theory of the elastic behaviour of multiphase materials.
   *J. Mech. Phys. Solids*, 11, 127–140.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Tuple, Union

import numpy as np
from scipy.optimize import brentq

from ..api.property import PyCSAMTObject

__all__ = [
    # petrophysical models
    "ArchieModel",
    "WaxmanSmitsModel",
    "HashinShtrikmanBounds",
    # hydraulic functions
    "kozeny_carman_K",
    "rho_to_hydraulic_conductivity",
    "transmissivity",
    "storativity",
    # water chemistry
    "rho_w_to_tds",
    "tds_to_rho_w",
    "ec_mscm_to_rho",
    "rho_to_ec_mscm",
    # EM-specific helpers
    "skin_depth",
    "bostick_depth",
    "aquifer_top_from_profile",
    "water_table_from_profile",
    "fractured_zone_K",
]

# ── physical constants ─────────────────────────────────────────────────────────

_MU0 = 4.0 * np.pi * 1e-7          # H/m
_PI2 = 2.0 * np.pi
_BOSTICK_CONST = 503.3              # sqrt(1 / (μ₀ π))


# ── helpers ────────────────────────────────────────────────────────────────────

def _arr(x: Any) -> np.ndarray:
    """Cast x to a float64 ndarray."""
    return np.asarray(x, dtype=float)


def _broadcast(*arrays) -> tuple:
    """Broadcast all arrays to a common shape and return as float64."""
    arrs = [_arr(a) for a in arrays]
    shape = np.broadcast_shapes(*[a.shape for a in arrs])
    return tuple(np.broadcast_to(a, shape).copy() for a in arrs)


def _scalar_or_array(result: np.ndarray, was_scalar: bool) -> Any:
    if was_scalar and result.ndim == 0:
        return float(result)
    if was_scalar and result.size == 1:
        return float(result.flat[0])
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Archie model
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ArchieModel(PyCSAMTObject):
    r"""Archie (1942) petrophysical model.

    Relates formation resistivity to porosity, water saturation, and pore-water
    resistivity via:

    .. math::

        \rho = a \cdot \rho_w \cdot \phi^{-m} \cdot S_w^{-n}

    Typical parameter ranges
    ------------------------
    +------------------+--------------+---------+------------------------------+
    | Parameter        | Range        | Default | Notes                        |
    +==================+==============+=========+==============================+
    | *a* (tortuosity) | 0.62 – 1.0   | 1.0     | 0.62 for sandstone (Humble) |
    +------------------+--------------+---------+------------------------------+
    | *m* (cementation)| 1.3 – 2.5   | 1.8     | 1.3 unconsolidated; 2+ frac |
    +------------------+--------------+---------+------------------------------+
    | *n* (saturation) | 1.8 – 2.5   | 2.0     | 2.0 clean sand (standard)   |
    +------------------+--------------+---------+------------------------------+

    Parameters
    ----------
    m : float
        Cementation exponent.
    n : float
        Saturation exponent.
    a : float
        Tortuosity factor (default 1.0).

    Examples
    --------
    >>> archie = ArchieModel(m=1.8, n=2.0, a=1.0)
    >>> archie.forward(phi=0.30, Sw=1.0, rho_w=0.025)   # fully saturated
    0.309...
    >>> archie.saturation(rho=1.2, phi=0.30, rho_w=0.025)
    0.508...
    """

    m: float = 1.8
    n: float = 2.0
    a: float = 1.0

    def __post_init__(self):
        if self.m <= 0:
            raise ValueError("Cementation exponent m must be positive.")
        if self.n <= 0:
            raise ValueError("Saturation exponent n must be positive.")
        if self.a <= 0:
            raise ValueError("Tortuosity factor a must be positive.")

    # ── forward (hydro → ρ) ───────────────────────────────────────────────────

    def formation_factor(self, phi: Any) -> Any:
        r"""Formation factor :math:`F = a \phi^{-m}`.

        Parameters
        ----------
        phi : array-like
            Porosity (fraction, 0–1).

        Returns
        -------
        F : ndarray
        """
        phi = _arr(phi)
        phi = np.clip(phi, 1e-4, 0.99)
        return self.a * phi ** (-self.m)

    def forward(
        self,
        phi: Any,
        Sw: Any,
        rho_w: Any,
    ) -> Any:
        r"""Formation resistivity from Archie's law.

        Parameters
        ----------
        phi : array-like
            Porosity (fraction, 0–1).
        Sw : array-like
            Water saturation (fraction, 0–1).
        rho_w : float or array-like
            Pore-water resistivity (Ω·m).

        Returns
        -------
        rho : ndarray  — formation resistivity (Ω·m)
        """
        phi, Sw, rho_w = _broadcast(phi, Sw, rho_w)
        was_scalar = all(np.ndim(x) == 0 for x in (phi, Sw, rho_w))
        phi = np.clip(phi, 1e-4, 0.99)
        Sw  = np.clip(Sw,  1e-4, 1.00)
        rho = self.a * rho_w * phi ** (-self.m) * Sw ** (-self.n)
        rho = np.clip(rho, 1e-2, 1e7)
        return _scalar_or_array(rho, was_scalar)

    # ── inverse (ρ → hydro) ───────────────────────────────────────────────────

    def saturation(
        self,
        rho: Any,
        phi: Any,
        rho_w: Any,
    ) -> Any:
        r"""Water saturation from measured resistivity (Archie inverse).

        .. math::

            S_w = \left(\frac{a \rho_w}{\rho \phi^m}\right)^{1/n}

        Parameters
        ----------
        rho : array-like
            Formation resistivity (Ω·m).
        phi : array-like
            Porosity (fraction, 0–1).
        rho_w : float or array-like
            Pore-water resistivity (Ω·m).

        Returns
        -------
        Sw : ndarray — water saturation (0–1)
        """
        rho, phi, rho_w = _broadcast(rho, phi, rho_w)
        was_scalar = all(np.ndim(x) == 0 for x in (rho, phi, rho_w))
        phi = np.clip(phi, 1e-4, 0.99)
        rho = np.clip(rho, 1e-2, 1e7)
        F   = self.a * phi ** (-self.m)
        Sw  = (F * rho_w / rho) ** (1.0 / self.n)
        Sw  = np.clip(Sw, 0.0, 1.0)
        return _scalar_or_array(Sw, was_scalar)

    def porosity(
        self,
        rho: Any,
        Sw: Any,
        rho_w: Any,
    ) -> Any:
        r"""Porosity from resistivity and saturation (Archie inverse).

        .. math::

            \phi = \left(\frac{a \rho_w S_w^{-n}}{\rho}\right)^{1/m}

        Parameters
        ----------
        rho : array-like
            Formation resistivity (Ω·m).
        Sw : array-like
            Water saturation (fraction, 0–1).
        rho_w : float or array-like
            Pore-water resistivity (Ω·m).

        Returns
        -------
        phi : ndarray — porosity (0–1)
        """
        rho, Sw, rho_w = _broadcast(rho, Sw, rho_w)
        was_scalar = all(np.ndim(x) == 0 for x in (rho, Sw, rho_w))
        Sw  = np.clip(Sw,  1e-4, 1.0)
        rho = np.clip(rho, 1e-2, 1e7)
        phi = (self.a * rho_w * Sw ** (-self.n) / rho) ** (1.0 / self.m)
        phi = np.clip(phi, 1e-4, 0.99)
        return _scalar_or_array(phi, was_scalar)

    def fluid_resistivity(
        self,
        rho: Any,
        phi: Any,
        Sw: Any,
    ) -> Any:
        r"""Pore-water resistivity from formation resistivity (Archie inverse).

        .. math::

            \rho_w = \frac{\rho \phi^m S_w^n}{a}

        Parameters
        ----------
        rho : array-like
            Formation resistivity (Ω·m).
        phi : array-like
            Porosity (fraction, 0–1).
        Sw : array-like
            Water saturation (fraction, 0–1).

        Returns
        -------
        rho_w : ndarray — pore-water resistivity (Ω·m)
        """
        rho, phi, Sw = _broadcast(rho, phi, Sw)
        was_scalar = all(np.ndim(x) == 0 for x in (rho, phi, Sw))
        phi = np.clip(phi, 1e-4, 0.99)
        Sw  = np.clip(Sw,  1e-4, 1.0)
        rho_w = rho * phi ** self.m * Sw ** self.n / self.a
        return _scalar_or_array(rho_w, was_scalar)

    def water_content(self, phi: Any, Sw: Any) -> Any:
        r"""Volumetric water content :math:`\theta = \phi \cdot S_w`."""
        phi, Sw = _broadcast(phi, Sw)
        was_scalar = all(np.ndim(x) == 0 for x in (phi, Sw))
        theta = np.clip(phi, 0, 1) * np.clip(Sw, 0, 1)
        return _scalar_or_array(theta, was_scalar)

    def __repr__(self) -> str:
        return (f"ArchieModel(m={self.m}, n={self.n}, a={self.a})")


# ─────────────────────────────────────────────────────────────────────────────
# Waxman-Smits model
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class WaxmanSmitsModel(PyCSAMTObject):
    r"""Waxman-Smits (1968) model for clay-bearing formations.

    Extends Archie's law with a surface conductivity term σ_s that represents
    the excess conductance associated with clay minerals (cation exchange
    capacity).  Use this model when clay content > ~10 % or when Archie
    overestimates saturation in shaly formations.

    .. math::

        \sigma = \frac{S_w^n}{F}\bigl(\sigma_w + \sigma_s / S_w\bigr)

    where *F = a φ^{-m}* is the formation factor.

    Parameters
    ----------
    m : float
        Cementation exponent (default 1.8).
    n : float
        Saturation exponent (default 2.0).
    a : float
        Tortuosity factor (default 1.0).
    sigma_s : float
        Surface (clay) conductivity S/m (default 0 → reduces to Archie).

    Notes
    -----
    σ_s can be estimated from CEC (cation exchange capacity) measurements
    or calibrated against borehole resistivity logs.  Values for AMT/MT
    targets typically range 0.001–0.05 S/m for moderate clay content.

    Examples
    --------
    >>> ws = WaxmanSmitsModel(m=1.8, n=2.0, sigma_s=0.01)
    >>> ws.forward(phi=0.30, Sw=0.70, sigma_w=40.0)   # sigma_w in mS/m
    8.31...  # Ω·m
    """

    m: float = 1.8
    n: float = 2.0
    a: float = 1.0
    sigma_s: float = 0.0          # S/m surface conductivity

    def __post_init__(self):
        if self.m <= 0 or self.n <= 0 or self.a <= 0:
            raise ValueError("Exponents m, n, a must all be positive.")
        if self.sigma_s < 0:
            raise ValueError("Surface conductivity sigma_s must be >= 0.")

    def _formation_factor(self, phi: np.ndarray) -> np.ndarray:
        phi = np.clip(phi, 1e-4, 0.99)
        return self.a * phi ** (-self.m)

    # ── forward ───────────────────────────────────────────────────────────────

    def forward(
        self,
        phi: Any,
        Sw: Any,
        sigma_w: Any,
    ) -> Any:
        r"""Formation resistivity from Waxman-Smits equation.

        Parameters
        ----------
        phi : array-like
            Porosity (fraction, 0–1).
        Sw : array-like
            Water saturation (fraction, 0–1).
        sigma_w : float or array-like
            Pore-water conductivity in **mS/m** (≈ 40 mS/m for fresh water).
            Internally converted to S/m.

        Returns
        -------
        rho : ndarray — formation resistivity (Ω·m)
        """
        phi, Sw, sigma_w = _broadcast(phi, Sw, sigma_w)
        was_scalar = all(np.ndim(x) == 0 for x in (phi, Sw, sigma_w))
        phi = np.clip(phi, 1e-4, 0.99)
        Sw  = np.clip(Sw,  1e-4, 1.0)
        sigma_w_si = sigma_w * 1e-3          # mS/m → S/m
        F  = self._formation_factor(phi)
        s  = (Sw ** self.n / F) * (sigma_w_si + self.sigma_s / Sw)
        s  = np.maximum(s, 1e-10)
        rho = np.clip(1.0 / s, 1e-2, 1e7)
        return _scalar_or_array(rho, was_scalar)

    # ── inverse ───────────────────────────────────────────────────────────────

    def saturation(
        self,
        rho: Any,
        phi: Any,
        sigma_w: Any,
        *,
        tol: float = 1e-8,
        max_iter: int = 50,
    ) -> Any:
        r"""Water saturation by numerical inversion of Waxman-Smits equation.

        Parameters
        ----------
        rho : array-like
            Formation resistivity (Ω·m).
        phi : array-like
            Porosity (fraction, 0–1).
        sigma_w : float or array-like
            Pore-water conductivity in mS/m.
        tol : float
            Convergence tolerance for brentq (default 1e-8).
        max_iter : int
            Maximum root-finding iterations (default 50).

        Returns
        -------
        Sw : ndarray — water saturation (0–1)
        """
        rho, phi, sigma_w = _broadcast(rho, phi, sigma_w)
        was_scalar = all(np.ndim(x) == 0 for x in (rho, phi, sigma_w))
        phi = np.clip(phi, 1e-4, 0.99)
        rho = np.clip(rho, 1e-2, 1e7)
        sigma_w_si = sigma_w * 1e-3

        result = np.empty(rho.shape)
        F_arr = self._formation_factor(phi)

        it = np.nditer([rho, phi, sigma_w_si, F_arr, result],
                       op_flags=[['readonly']] * 4 + [['writeonly']])
        for rho_i, phi_i, sw_i, F_i, out_i in it:
            sigma_obs = 1.0 / float(rho_i)
            sw_si     = float(sw_i)
            F_v       = float(F_i)
            n, ss     = self.n, self.sigma_s

            def residual(S):
                return (S ** n / F_v) * (sw_si + ss / max(S, 1e-10)) - sigma_obs

            # Archie initial guess
            if sw_si > 1e-12:
                S0 = np.clip((F_v * sw_si / (sigma_obs * F_v)) ** (1.0 / n),
                             1e-3, 1.0)
            else:
                S0 = 0.5

            # Check bracket
            try:
                fa, fb = residual(1e-4), residual(1.0)
                if fa * fb < 0:
                    Sw_sol = brentq(residual, 1e-4, 1.0,
                                    xtol=tol, maxiter=max_iter)
                else:
                    Sw_sol = float(np.clip(S0, 0.0, 1.0))
            except Exception:
                Sw_sol = float(np.clip(S0, 0.0, 1.0))

            out_i[...] = np.clip(Sw_sol, 0.0, 1.0)

        return _scalar_or_array(result, was_scalar)

    def porosity(
        self,
        rho: Any,
        Sw: Any,
        sigma_w: Any,
    ) -> Any:
        r"""Porosity by direct inversion of Waxman-Smits (known saturation).

        Rearranging WS:

        .. math::

            F = \frac{S_w^n \,(\sigma_w + \sigma_s / S_w)}{\sigma_{obs}}
            \Rightarrow \phi = \left(F / a\right)^{-1/m}

        Parameters
        ----------
        rho : array-like
            Formation resistivity (Ω·m).
        Sw : array-like
            Water saturation (fraction, 0–1).
        sigma_w : float or array-like
            Pore-water conductivity in mS/m.

        Returns
        -------
        phi : ndarray — porosity (0–1)
        """
        rho, Sw, sigma_w = _broadcast(rho, Sw, sigma_w)
        was_scalar = all(np.ndim(x) == 0 for x in (rho, Sw, sigma_w))
        Sw  = np.clip(Sw,  1e-4, 1.0)
        rho = np.clip(rho, 1e-2, 1e7)
        sw_si    = sigma_w * 1e-3
        sigma_obs = 1.0 / rho
        F_needed = Sw ** self.n * (sw_si + self.sigma_s / Sw) / sigma_obs
        F_needed = np.maximum(F_needed, 1e-6)
        phi = (F_needed / self.a) ** (-1.0 / self.m)
        phi = np.clip(phi, 1e-4, 0.99)
        return _scalar_or_array(phi, was_scalar)

    def __repr__(self) -> str:
        return (f"WaxmanSmitsModel(m={self.m}, n={self.n}, "
                f"a={self.a}, sigma_s={self.sigma_s})")


# ─────────────────────────────────────────────────────────────────────────────
# Hashin-Shtrikman bounds
# ─────────────────────────────────────────────────────────────────────────────

class HashinShtrikmanBounds(PyCSAMTObject):
    r"""Hashin-Shtrikman (1962) bounds on effective resistivity.

    Provides theoretical upper and lower bounds for a two-phase medium
    (rock matrix + fluid) without assumptions about microstructure.  The
    bounds are tighter than Voigt-Reuss and are used to validate whether
    an Archie-derived porosity is physically plausible.

    For conductivity σ (= 1/ρ):

    .. math::

        \sigma^{HS-} \leq \sigma_{eff} \leq \sigma^{HS+}

    Parameters
    ----------
    rho_matrix : float
        Resistivity of the rock matrix (Ω·m).  Typical: granite 10⁴–10⁶,
        sandstone 100–1000, limestone 500–5000.
    rho_fluid : float
        Resistivity of the pore fluid (Ω·m).  Typical: fresh water 20–100,
        saline water 0.1–5.

    Examples
    --------
    >>> hs = HashinShtrikmanBounds(rho_matrix=1000.0, rho_fluid=25.0)
    >>> lower, upper = hs.bounds(phi=0.25)
    >>> print(f"{lower:.1f}  {upper:.1f}")
    """

    def __init__(
        self,
        rho_matrix: float,
        rho_fluid: float,
    ) -> None:
        if rho_matrix <= 0 or rho_fluid <= 0:
            raise ValueError("Resistivities must be positive.")
        self.rho_matrix = float(rho_matrix)
        self.rho_fluid  = float(rho_fluid)

    def bounds(
        self,
        phi: Any,
    ) -> Tuple[np.ndarray, np.ndarray]:
        r"""Compute Hashin-Shtrikman lower and upper conductivity bounds.

        Parameters
        ----------
        phi : array-like
            Porosity (volume fraction of fluid, 0–1).

        Returns
        -------
        rho_lower : ndarray — lower resistivity bound (HS+, fluid-connected)
        rho_upper : ndarray — upper resistivity bound (HS-, isolated fluid)
        """
        phi = _arr(phi)
        phi = np.clip(phi, 1e-6, 1.0 - 1e-6)
        f_phi   = phi              # fluid volume fraction
        f_mat   = 1.0 - phi        # matrix volume fraction
        s_f     = 1.0 / self.rho_fluid
        s_m     = 1.0 / self.rho_matrix

        # HS+ (fluid is the more conductive reference = lower ρ bound)
        denom_plus  = f_mat / (s_m - s_f + 3.0 * s_f)
        sigma_plus  = s_f + f_mat / (1.0 / (s_m - s_f) + f_phi / (3.0 * s_f))

        # HS- (matrix is reference = upper ρ bound)
        sigma_minus = s_m + f_phi / (1.0 / (s_f - s_m) + f_mat / (3.0 * s_m))

        rho_lower = np.clip(1.0 / sigma_plus,  1e-4, 1e8)
        rho_upper = np.clip(1.0 / sigma_minus, 1e-4, 1e8)
        return rho_lower, rho_upper

    def in_bounds(
        self,
        rho: Any,
        phi: Any,
        *,
        margin: float = 0.0,
    ) -> np.ndarray:
        """Return boolean mask: True where *rho* lies within HS bounds ± margin (log₁₀)."""
        rho = _arr(rho)
        lower, upper = self.bounds(phi)
        log_rho   = np.log10(np.clip(rho,   1e-4, 1e8))
        log_lower = np.log10(np.clip(lower, 1e-4, 1e8))
        log_upper = np.log10(np.clip(upper, 1e-4, 1e8))
        return (log_rho >= log_lower - margin) & (log_rho <= log_upper + margin)

    def __repr__(self) -> str:
        return (f"HashinShtrikmanBounds(rho_matrix={self.rho_matrix}, "
                f"rho_fluid={self.rho_fluid})")


# ─────────────────────────────────────────────────────────────────────────────
# Hydraulic property functions
# ─────────────────────────────────────────────────────────────────────────────

def kozeny_carman_K(
    phi: Any,
    *,
    d50_m: float = 2.5e-4,
    C: float = 180.0,
    T: float = 0.5,
    gravity: float = 9.81,
    kinematic_viscosity: float = 1e-6,
) -> Any:
    r"""Hydraulic conductivity from porosity via Kozeny-Carman equation.

    .. math::

        K = \frac{g}{\nu} \cdot \frac{d_{50}^2}{C}
            \cdot \frac{\phi^3}{(1-\phi)^2} \cdot T

    where :math:`g/\nu` converts intrinsic permeability (m²) to hydraulic
    conductivity (m/s): :math:`g = 9.81` m/s², :math:`\nu \approx 10^{-6}`
    m²/s at 20 °C.

    Parameters
    ----------
    phi : array-like
        Porosity (fraction, 0–1).
    d50_m : float
        Median grain size in **metres** (default 2.5×10⁻⁴ m = 0.25 mm,
        upper fine sand / lower medium sand).
    C : float
        Kozeny constant × shape factor (default 180 for spherical grains).
    T : float
        Tortuosity correction factor (default 0.5; range 0.3–0.7).
    gravity : float
        Gravitational acceleration (m/s²; default 9.81).
    kinematic_viscosity : float
        Kinematic viscosity of water (m²/s; default 10⁻⁶ at 20 °C).
        Use 1.31×10⁻⁶ at 10 °C or 0.80×10⁻⁶ at 30 °C if temperature matters.

    Returns
    -------
    K : ndarray — hydraulic conductivity (m/s)

    Notes
    -----
    Typical K values at φ ≈ 0.30 with matching d50:

    * Clean gravel (d50 = 5 mm)   : ~10⁻² – 10⁻¹ m/s
    * Coarse sand  (d50 = 1 mm)   : ~10⁻³ – 10⁻² m/s
    * Medium sand  (d50 = 0.5 mm) : ~10⁻⁴ – 10⁻³ m/s
    * Fine sand    (d50 = 0.1 mm) : ~10⁻⁶ – 10⁻⁵ m/s
    * Silty sand   (d50 = 0.05 mm): ~10⁻⁷ – 10⁻⁶ m/s
    * Clay / silt                  : < 10⁻⁹ m/s

    Examples
    --------
    >>> kozeny_carman_K(0.30, d50_m=1e-3)   # coarse sand, φ = 0.30
    1.16e-03  # m/s  (within typical coarse-sand range)
    """
    phi = _arr(phi)
    was_scalar = phi.ndim == 0
    phi = np.clip(phi, 1e-4, 0.99)
    K = (gravity / kinematic_viscosity) * T * d50_m**2 * phi**3 / (C * (1.0 - phi)**2)
    K = np.clip(K, 1e-15, 1e2)
    return _scalar_or_array(K, was_scalar)


def rho_to_hydraulic_conductivity(
    rho: Any,
    archie: ArchieModel,
    *,
    rho_w: float = 20.0,
    phi_prior: float = 0.25,
    Sw: float = 1.0,
    d50_m: float = 2.5e-4,
    C: float = 180.0,
    T: float = 0.5,
    gravity: float = 9.81,
    kinematic_viscosity: float = 1e-6,
) -> Any:
    r"""Estimate hydraulic conductivity K from EM resistivity.

    Chain: ρ →  φ (Archie inverse)  →  K (Kozeny-Carman).

    This is valid only in **fully or near-fully saturated** zones (S_w ≈ 1).
    For unsaturated zones use the forward Archie model to estimate S_w first,
    then re-run on the saturated resistivity.

    Parameters
    ----------
    rho : array-like
        Formation resistivity (Ω·m) from EM inversion.
    archie : ArchieModel
        Calibrated Archie model (m, n, a).
    rho_w : float
        Pore-water resistivity (Ω·m).  Default 20 Ω·m (fresh water at 25 °C,
        EC ≈ 0.5 mS/cm).  Use ~0.2 Ω·m for seawater, ~2–5 Ω·m for brackish.
    phi_prior : float
        Prior (reference) porosity used to clip unrealistic values.
    Sw : float
        Water saturation assumed for the inversion (default 1.0 = saturated).
    d50_m : float
        Median grain size (m) for Kozeny-Carman.
    C, T : float
        Kozeny constant and tortuosity (see :func:`kozeny_carman_K`).
    gravity : float
        Gravitational acceleration (m/s²; default 9.81).
    kinematic_viscosity : float
        Kinematic viscosity of water (m²/s; default 10⁻⁶ at 20 °C).

    Returns
    -------
    K : ndarray — hydraulic conductivity estimate (m/s)
    """
    phi = archie.porosity(rho, Sw, rho_w)
    # soft-clip by prior (outliers clamped to 3× prior)
    phi = np.clip(_arr(phi), 1e-4, min(0.99, 3.0 * phi_prior))
    return kozeny_carman_K(
        phi, d50_m=d50_m, C=C, T=T,
        gravity=gravity, kinematic_viscosity=kinematic_viscosity,
    )


def transmissivity(
    K: Any,
    thickness: Any,
) -> Any:
    r"""Aquifer transmissivity :math:`T = K \cdot b` (m²/s).

    Parameters
    ----------
    K : array-like
        Hydraulic conductivity (m/s).
    thickness : array-like
        Saturated thickness of the aquifer layer (m).

    Returns
    -------
    T : ndarray — transmissivity (m²/s)
    """
    K, b = _broadcast(K, thickness)
    was_scalar = all(np.ndim(x) == 0 for x in (K, thickness))
    result = np.clip(_arr(K), 0.0, None) * np.clip(_arr(b), 0.0, None)
    return _scalar_or_array(result, was_scalar)


def storativity(
    phi: Any,
    thickness: Any,
    *,
    specific_storage: float = 1e-4,
) -> Any:
    r"""Aquifer storativity :math:`S = S_s \cdot b` (confined) or φ (unconfined).

    For a **confined** aquifer: *S = Ss × b* (dimensionless).
    For an **unconfined** aquifer: *S ≈ φ* (specific yield).

    Parameters
    ----------
    phi : array-like
        Porosity (fraction, 0–1).  Used as specific yield for unconfined.
    thickness : array-like
        Saturated thickness (m).
    specific_storage : float
        Specific storage Ss (m⁻¹) for confined conditions (default 10⁻⁴).

    Returns
    -------
    S_confined : ndarray — confined storativity (dimensionless)
    S_unconfined : ndarray — unconfined storativity ≈ specific yield (fraction)
    """
    phi, b = _broadcast(phi, thickness)
    was_scalar = all(np.ndim(x) == 0 for x in (phi, thickness))
    S_confined   = np.clip(_arr(b), 0.0, None) * specific_storage
    S_unconfined = np.clip(_arr(phi), 0.0, 0.99)
    return (
        _scalar_or_array(S_confined,   was_scalar),
        _scalar_or_array(S_unconfined, was_scalar),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Water chemistry
# ─────────────────────────────────────────────────────────────────────────────

def rho_w_to_tds(rho_w: Any, *, temp_c: float = 25.0) -> Any:
    r"""Pore-water resistivity → Total Dissolved Solids (mg/L).

    Uses the empirical Keller (1966) relation:

    .. math::

        \mathrm{TDS} \approx \frac{640}{\mathrm{EC}_{(dS/m)}}
        = \frac{6400}{\sigma_{w(mS/cm)}}

    With a temperature correction from the reference 25 °C:

    .. math::

        \rho_{w,25} = \rho_{w,T} / [1 + 0.02(T - 25)]

    Parameters
    ----------
    rho_w : array-like
        Pore-water resistivity (Ω·m) at temperature *temp_c*.
    temp_c : float
        Measurement temperature (°C; default 25).

    Returns
    -------
    tds : ndarray — TDS in mg/L (drinking limit: 500 mg/L; saline: > 3 000)
    """
    rho_w = _arr(rho_w)
    was_scalar = rho_w.ndim == 0
    rho_w = np.clip(rho_w, 1e-3, 1e4)
    # temperature correction to 25 °C
    rho_25 = rho_w / (1.0 + 0.02 * (temp_c - 25.0))
    rho_25 = np.clip(rho_25, 1e-3, 1e4)
    # EC in mS/cm from Ω·m: EC_mScm = 100 / (10 * rho) = 10 / rho
    ec_mScm = 10.0 / rho_25
    tds = 640.0 * ec_mScm        # Keller (1966)
    return _scalar_or_array(tds, was_scalar)


def tds_to_rho_w(tds: Any, *, temp_c: float = 25.0) -> Any:
    r"""TDS (mg/L) → pore-water resistivity (Ω·m).

    Inverse of :func:`rho_w_to_tds`.

    Parameters
    ----------
    tds : array-like
        Total dissolved solids (mg/L).
    temp_c : float
        Temperature (°C; default 25).

    Returns
    -------
    rho_w : ndarray — pore-water resistivity (Ω·m)
    """
    tds = np.clip(_arr(tds), 1.0, 1e7)
    was_scalar = tds.ndim == 0
    ec_mScm = tds / 640.0
    rho_25  = 10.0 / np.maximum(ec_mScm, 1e-6)
    # undo temperature correction
    rho_w   = rho_25 * (1.0 + 0.02 * (temp_c - 25.0))
    return _scalar_or_array(np.clip(rho_w, 1e-3, 1e4), was_scalar)


def ec_mscm_to_rho(ec: Any, *, temp_c: float = 25.0) -> Any:
    """Electrical conductivity (mS/cm) → resistivity (Ω·m)."""
    ec  = np.clip(_arr(ec), 1e-6, 1e4)
    was_scalar = ec.ndim == 0
    rho_25 = 10.0 / ec
    rho    = rho_25 * (1.0 + 0.02 * (temp_c - 25.0))
    return _scalar_or_array(np.clip(rho, 1e-4, 1e6), was_scalar)


def rho_to_ec_mscm(rho: Any, *, temp_c: float = 25.0) -> Any:
    """Formation/water resistivity (Ω·m) → EC (mS/cm)."""
    rho = np.clip(_arr(rho), 1e-4, 1e6)
    was_scalar = rho.ndim == 0
    rho_25 = rho / (1.0 + 0.02 * (temp_c - 25.0))
    ec     = 10.0 / np.maximum(rho_25, 1e-6)
    return _scalar_or_array(ec, was_scalar)


# ─────────────────────────────────────────────────────────────────────────────
# EM-specific helpers
# ─────────────────────────────────────────────────────────────────────────────

def skin_depth(rho: Any, freq: Any) -> Any:
    r"""EM skin depth (penetration depth) in metres.

    .. math::

        \delta = 503 \sqrt{\rho / f}

    This gives the 1/e amplitude depth for a plane EM wave, useful for
    estimating the sensitivity depth of AMT/MT measurements.

    Parameters
    ----------
    rho : array-like
        Apparent (or formation) resistivity (Ω·m).
    freq : array-like
        Frequency (Hz).

    Returns
    -------
    delta : ndarray — skin depth (m)
    """
    rho, freq = _broadcast(rho, freq)
    was_scalar = all(np.ndim(x) == 0 for x in (rho, freq))
    rho  = np.clip(rho,  1e-2, 1e7)
    freq = np.clip(freq, 1e-6, 1e6)
    delta = _BOSTICK_CONST * np.sqrt(rho / freq)
    return _scalar_or_array(delta, was_scalar)


def bostick_depth(rho_a: Any, freq: Any) -> Any:
    r"""Bostick (1977) pseudo-depth for MT/AMT/TDEM data.

    .. math::

        d_B = 503 \sqrt{\rho_a / f}

    Identical numerically to skin depth but interpreted as the depth of
    maximum sensitivity for a given apparent resistivity and frequency.
    Use the full sequence {(ρ_a, f)} to build a pseudo-1D depth profile.

    Parameters
    ----------
    rho_a : array-like
        Apparent resistivity (Ω·m).
    freq : array-like
        Frequency (Hz).  For TDEM use equivalent frequency 1/(2πt).

    Returns
    -------
    d : ndarray — Bostick pseudo-depth (m)
    """
    return skin_depth(rho_a, freq)


def aquifer_top_from_profile(
    rho_log10: Any,
    z_centers: Any,
    *,
    rho_threshold_ohm_m: float = 300.0,
    direction: str = "low",
    min_depth: float = 0.0,
) -> Optional[float]:
    r"""Detect the top of an aquifer/conductor in a 1-D resistivity profile.

    Scans the profile from the surface down and returns the depth of the
    first cell that crosses the threshold.

    Parameters
    ----------
    rho_log10 : array-like (n_z,)
        Log₁₀ resistivity profile (one station column).
    z_centers : array-like (n_z,)
        Depth of each cell centre (m, positive downward).
    rho_threshold_ohm_m : float
        Resistivity threshold.  For aquifer detection (low resistivity):
        crossing downward indicates aquifer top.
    direction : {'low', 'high'}
        ``'low'`` — detect transition to low ρ (aquifer).
        ``'high'`` — detect transition to high ρ (basement / caprock).
    min_depth : float
        Ignore cells shallower than *min_depth* metres (default 0).

    Returns
    -------
    depth : float or None — top depth (m), or None if threshold not crossed.
    """
    rho  = 10.0 ** np.asarray(rho_log10, dtype=float)
    z    = np.asarray(z_centers, dtype=float)
    thr  = float(rho_threshold_ohm_m)
    mask = z >= min_depth

    for zi, ri in zip(z[mask], rho[mask]):
        if direction == "low" and ri <= thr:
            return float(zi)
        if direction == "high" and ri >= thr:
            return float(zi)
    return None


def water_table_from_profile(
    rho_log10: Any,
    z_centers: Any,
    archie: ArchieModel,
    *,
    rho_w: float = 0.025,
    Sw_threshold: float = 0.85,
    min_depth: float = 0.5,
) -> Optional[float]:
    r"""Estimate water table depth from a 1-D EM resistivity column.

    Converts each cell's log₁₀ρ to water saturation S_w via Archie's
    inverse, then finds the shallowest depth where S_w ≥ *Sw_threshold*.

    Parameters
    ----------
    rho_log10 : array-like (n_z,)
        Log₁₀ resistivity column from EM inversion.
    z_centers : array-like (n_z,)
        Depth of each cell (m).
    archie : ArchieModel
        Calibrated Archie model.
    rho_w : float
        Pore-water resistivity (Ω·m).
    Sw_threshold : float
        Saturation value defining the water table (default 0.85).
    min_depth : float
        Minimum search depth (m) to skip near-surface noise.

    Returns
    -------
    depth : float or None — estimated water table depth (m).
    """
    rho_log10 = np.asarray(rho_log10, dtype=float)
    z         = np.asarray(z_centers, dtype=float)
    rho       = 10.0 ** rho_log10

    # porosity prior from surface ρ (shallowest cell)
    phi_prior = 0.25
    for zi, rho_i in zip(z, rho):
        if zi < min_depth or not np.isfinite(rho_i) or rho_i <= 0:
            continue
        Sw_i = float(archie.saturation(rho_i, phi_prior, rho_w))
        if Sw_i >= Sw_threshold:
            return float(zi)
    return None


def fractured_zone_K(
    rho: Any,
    *,
    rho_matrix: float = 5000.0,
    aperture_m: float = 1e-3,
    cubic_factor: float = 1.0 / 12.0,
    kinematic_viscosity: float = 1e-6,
    gravity: float = 9.81,
) -> Any:
    r"""Hydraulic conductivity of a fractured zone from EM resistivity.

    Uses the parallel-plate cubic law for fracture flow:

    .. math::

        K = \frac{g \, b^3}{12 \nu} \cdot f_\mathrm{vol}

    where the volumetric fracture fraction *f_vol* is estimated from the
    resistivity contrast against the matrix:

    .. math::

        f_\mathrm{vol} \approx 1 - \frac{\rho}{\rho_\mathrm{matrix}}

    This is a rough first-pass estimate for fractured basement in AMT
    surveys.  For rigorous fracture characterisation use the Bahr or
    distortion framework in :mod:`pycsamt.emtools`.

    Parameters
    ----------
    rho : array-like
        Measured formation resistivity (Ω·m).
    rho_matrix : float
        Background (intact rock) resistivity (Ω·m).
    aperture_m : float
        Representative fracture aperture (m).
    cubic_factor : float
        Cubic-law pre-factor (1/12 for parallel plates).
    kinematic_viscosity : float
        Kinematic viscosity of water at ~20 °C (m²/s, default 10⁻⁶).
    gravity : float
        Gravitational acceleration (m/s²).

    Returns
    -------
    K : ndarray — effective hydraulic conductivity (m/s)
    """
    rho = np.clip(_arr(rho), 1e-2, rho_matrix)
    was_scalar = rho.ndim == 0
    # floor at 1e-4 instead of 0 so cells with rho ≈ rho_matrix yield a finite
    # (very negative) log₁₀K rather than K=0 → NaN in the section plot
    f_vol = np.clip(1.0 - rho / rho_matrix, 1e-4, 0.99)
    K_cubic = gravity * aperture_m**3 / (kinematic_viscosity / cubic_factor)
    K = K_cubic * f_vol
    return _scalar_or_array(np.clip(K, 0.0, 1e-1), was_scalar)
