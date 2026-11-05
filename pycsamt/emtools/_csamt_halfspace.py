"""Exact CSAMT half-space response for near-field diagnosis and correction.

For a horizontal electric point dipole on a homogeneous half-space, the
broadside (equatorial) Cagniard impedance :math:`Z = E_x/H_y` departs from
the plane-wave impedance :math:`Z_{pw}=\\sqrt{i\\omega\\mu_0\\rho}` by a complex
ratio that depends only on the induction number :math:`u = r/\\delta`, with
:math:`\\delta=\\sqrt{2\\rho/(\\omega\\mu_0)}`:

.. math::  R(u) = Z_{hs}(u) / Z_{pw}.

``R`` is tabulated once with :mod:`empymod` (quasi-static, full
Hankel-transform solution, Werthmüller 2017) and interpolated afterwards.
Both the electric and the magnetic field depart from plane-wave behaviour in
the near field, so ``R`` is **not** the electric-field-only factor
:math:`1-3/p^2+3/p^3` used by earlier versions of pyCSAMT, which failed a
check against this exact solution.

The correction strategy is the standard *equivalent half-space* one: for
each datum, find the half-space resistivity whose exact CSAMT response has
the observed apparent resistivity at the same offset and frequency, then
divide the observed impedance by ``R`` at that half-space. On a homogeneous
half-space this recovers the plane-wave impedance exactly; on a layered earth
it is an approximation that degrades when resistivity varies strongly over
the source-receiver distance.

Limitations: broadside geometry, point dipole (grounded wires much shorter
than the offset), 1-D earth under the whole spread.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

MU0 = 4.0e-7 * np.pi
U_MIN = 1.0e-2
U_MAX = 5.0e2
_RHO0 = 100.0
_R0 = 1000.0


@lru_cache(maxsize=1)
def _table() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(log_u, log_abs_R, arg_R)`` on a dense induction-number grid."""
    import empymod

    u = np.logspace(np.log10(U_MIN), np.log10(U_MAX), 700)
    omega = 2.0 * _RHO0 * u**2 / (MU0 * _R0**2)
    freq = omega / (2.0 * np.pi)
    common = dict(
        src=[0.0, 0.0, 0.001],
        depth=[0.0],
        res=[2e14, _RHO0],
        freqtime=freq,
        epermH=[0.0, 0.0],
        epermV=[0.0, 0.0],
        verb=0,
    )
    ex = np.asarray(empymod.dipole(rec=[0.0, _R0, 0.001], ab=11, **common))
    hy = np.asarray(empymod.dipole(rec=[0.0, _R0, 0.001], ab=51, **common))
    z_hs = ex / hy
    z_pw = np.sqrt(1j * omega * MU0 * _RHO0)
    ratio = z_hs / z_pw
    return np.log(u), np.log(np.abs(ratio)), np.unwrap(np.angle(ratio))


def halfspace_ratio(u: np.ndarray) -> np.ndarray:
    """Complex ratio ``Z_hs/Z_pw`` at induction number(s) ``u = r/delta``.

    Values above :data:`U_MAX` return 1 (plane-wave). Values below
    :data:`U_MIN` are clamped to :data:`U_MIN`.
    """
    log_u, log_abs, arg = _table()
    lu = np.log(np.clip(np.asarray(u, dtype=float), U_MIN, None))
    mag = np.exp(np.interp(lu, log_u, log_abs, right=0.0))
    ph = np.interp(lu, log_u, arg, right=0.0)
    return mag * np.exp(1j * ph)


def equivalent_halfspace(
    rho_a: np.ndarray, freq: np.ndarray, offset: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Equivalent half-space for observed Cagniard apparent resistivities.

    Solves ``rho * |R(u(rho))|**2 = rho_a`` for each datum. Since
    ``rho = pi f mu0 r**2 / u**2``, this is ``|R(u)|**2 / u**2 = rho_a /
    (pi f mu0 r**2)``, whose left side decreases monotonically with ``u``.

    Returns
    -------
    rho_h : ndarray
        Equivalent half-space resistivity (ohm-m).
    u : ndarray
        Induction number ``r/delta`` of that half-space.
    in_range : ndarray of bool
        ``False`` where the datum lies beyond the tabulated near-field limit
        (``u < U_MIN``); such data are too deep in the near field to correct.
    """
    log_u, log_abs, _ = _table()
    rho_a = np.asarray(rho_a, dtype=float)
    freq = np.asarray(freq, dtype=float)
    target = np.log(np.maximum(rho_a, 1e-30) / (np.pi * freq * MU0 * offset**2))
    log_h = 2.0 * log_abs - 2.0 * log_u  # decreasing in u
    order = np.argsort(log_h)
    lu = np.interp(target, log_h[order], log_u[order])
    # Beyond the far end of the table the response is plane-wave: u from rho_a.
    far = target < log_h.min()
    lu = np.where(far, 0.5 * (np.log(np.pi * freq * MU0 * offset**2) - np.log(np.maximum(rho_a, 1e-30))), lu)
    u = np.exp(lu)
    in_range = target <= log_h.max()
    rho_h = np.pi * freq * MU0 * offset**2 / u**2
    return rho_h, u, in_range


def near_field_ratio(
    rho_a: np.ndarray, freq: np.ndarray, offset: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(R, u, in_range)`` evaluated at the equivalent half-space."""
    _, u, ok = equivalent_halfspace(rho_a, freq, offset)
    return halfspace_ratio(u), u, ok
