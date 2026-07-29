# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MT utility functions — Python translation of MATLAB ioAscii/ helper scripts.

Translates:
  skindepth.m  → :func:`skin_depth`

Also provides unit and format conversion helpers used across the iotools
subpackage.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "skin_depth",
    "imp_units_factor",
    "loge_to_log10",
    "log10_to_loge",
    "loge_to_linear",
    "linear_to_loge",
]

# Permeability of free space (H/m)
_MU_0 = 4.0 * np.pi * 1e-7


# -------------------------------------------------------------------------
# Skin depth
# -------------------------------------------------------------------------


def skin_depth(
    period: float | np.ndarray,
    rho: float | np.ndarray = 100.0,
) -> float | np.ndarray:
    """Electromagnetic skin depth δ in metres.

    Equivalent to the ``skindepth.m`` script::

        δ = √(2ρ / (ω μ₀))

    where ω = 2π/T.

    Parameters
    ----------
    period : float or array-like
        Period in seconds.
    rho : float or array-like, optional
        Resistivity in Ω·m.  Defaults to 100 Ω·m.

    Returns
    -------
    float or ndarray
        Skin depth in **metres**.

    Examples
    --------
    >>> skin_depth(1.0)  # 1 s period, 100 Ω·m
    503.292...
    >>> skin_depth([1.0, 10.0, 100.0], rho=10.0)
    array([159.154..., 503.292..., 1591.549...])
    """
    T = np.asarray(period, dtype=float)
    rho = np.asarray(rho, dtype=float)
    omega = 2.0 * np.pi / T
    return np.sqrt(2.0 * rho / (omega * _MU_0))


# -------------------------------------------------------------------------
# Impedance unit conversion
# -------------------------------------------------------------------------

# Conversion factors to SI ([V/m]/[T])
_UNIT_TO_SI: dict[str, float] = {
    "[v/m]/[t]": 1.0,
    "v/m/t": 1.0,
    "[mv/km]/[nt]": 1.0,  # same numerical value as SI for MT
    "mv/km/nt": 1.0,
    "ohm": 1.0,
    "mt": 1.0,
}


def imp_units_factor(from_units: str, to_units: str) -> float:
    """Return the multiplicative conversion factor between impedance unit strings.

    Equivalent to ``ImpUnits(from, to)`` called implicitly in
    ``readZ_3D.m``.  In the ModEM context only SI and ``[mV/km]/[nT]``
    are common; both map to the same numerical value so the factor is
    always 1.0 for the supported units.

    Parameters
    ----------
    from_units, to_units : str
        Unit strings (case-insensitive).  Recognised values:
        ``'[V/m]/[T]'``, ``'[mV/km]/[nT]'``, ``'Ohm'``.

    Returns
    -------
    float

    Raises
    ------
    ValueError
        If either unit string is not recognised.
    """
    fu = from_units.strip().lower()
    tu = to_units.strip().lower()
    for key in (fu, tu):
        if key not in _UNIT_TO_SI:
            raise ValueError(
                f"Unknown impedance unit '{key}'. Known: {list(_UNIT_TO_SI)}"
            )
    return _UNIT_TO_SI[fu] / _UNIT_TO_SI[tu]


# -------------------------------------------------------------------------
# Resistivity encoding conversions
# -------------------------------------------------------------------------


def loge_to_log10(rho_loge: np.ndarray) -> np.ndarray:
    """Convert ``ln(ρ)`` to ``log₁₀(ρ)``."""
    return np.asarray(rho_loge, dtype=float) / np.log(10.0)


def log10_to_loge(rho_log10: np.ndarray) -> np.ndarray:
    """Convert ``log₁₀(ρ)`` to ``ln(ρ)``."""
    return np.asarray(rho_log10, dtype=float) * np.log(10.0)


def loge_to_linear(rho_loge: np.ndarray) -> np.ndarray:
    """Convert ``ln(ρ)`` to linear resistivity (Ω·m)."""
    return np.exp(np.asarray(rho_loge, dtype=float))


def linear_to_loge(rho: np.ndarray) -> np.ndarray:
    """Convert linear resistivity (Ω·m) to ``ln(ρ)``."""
    rho = np.asarray(rho, dtype=float)
    safe = np.where(rho > 0, rho, np.finfo(float).tiny)
    return np.log(safe)
