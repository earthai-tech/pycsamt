# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Optional Numba-accelerated kernel for the 1-D layered-earth recursion.

:mod:`.forward` evaluates the same impedance recursion once per frequency,
independently -- the frequencies never interact. The NumPy implementation in
:mod:`.forward` vectorizes that recursion *over* frequency and loops only
over layers, which is already the efficient loop order; the remaining cost is
per-layer NumPy dispatch overhead on small arrays (tens of frequencies), not
floating-point work. Rewriting the same recursion as a per-frequency scalar
loop and compiling it with Numba removes that dispatch overhead entirely.

This module is optional. When Numba is not installed, :func:`impedance_kernel`
is ``None`` and callers fall back to the pure-NumPy recursion in
:mod:`.forward`, which remains the source of truth: this kernel must stay
numerically equivalent to it, verified by
``pycsamt/models/occam1d/tests/test_occam1d_numba_kernel.py``.
"""

from __future__ import annotations

import cmath

import numpy as np

__all__ = ["HAS_NUMBA", "impedance_kernel"]

try:
    from numba import njit

    HAS_NUMBA = True
except ImportError:  # pragma: no cover - exercised only without numba
    HAS_NUMBA = False

if HAS_NUMBA:

    @njit(cache=True)
    def _ctanh(z, thick_limit):
        """Stable complex tanh for an argument with non-negative real part."""
        if z.real >= thick_limit:
            return complex(1.0, 0.0)
        e2z = cmath.exp(2.0 * z)
        return (e2z - 1.0) / (e2z + 1.0)

    @njit(cache=True)
    def _impedance_kernel(rho, thickness, omega, mu, thick_limit, eps, tiny):
        """Per-frequency scalar recursion, JIT-compiled.

        Returns
        -------
        impedance : complex128[:]
        bad_frequency, bad_layer : int64
            ``-1`` when the recursion never became singular; otherwise the
            zero-based frequency index and layer index of the first (in
            top-down recursion order, i.e. deepest first) singular step.
        """
        n_frequency = omega.shape[0]
        n_layers = rho.shape[0]
        impedance = np.empty(n_frequency, dtype=np.complex128)
        bad_frequency = -1
        bad_layer = -1
        for f in range(n_frequency):
            iwmu = complex(0.0, omega[f] * mu)
            z = cmath.sqrt(iwmu * rho[n_layers - 1])
            singular = False
            layer_at_fail = -1
            for j in range(n_layers - 2, -1, -1):
                intrinsic = cmath.sqrt(iwmu * rho[j])
                propagation = cmath.sqrt(iwmu / rho[j])
                transfer = _ctanh(propagation * thickness[j], thick_limit)
                numerator = z + intrinsic * transfer
                denominator = intrinsic + z * transfer
                scale = max(abs(intrinsic), abs(z * transfer))
                tolerance = eps * max(scale, tiny)
                if abs(denominator) <= tolerance:
                    singular = True
                    layer_at_fail = j
                    break
                z = intrinsic * numerator / denominator
            impedance[f] = z
            if singular and layer_at_fail > bad_layer:
                bad_layer = layer_at_fail
                bad_frequency = f
        return impedance, bad_frequency, bad_layer

    def impedance_kernel(rho, thickness, omega, mu, thick_limit):
        """Return ``(impedance, bad_frequency, bad_layer)`` via the JIT kernel."""
        finfo = np.finfo(float)
        return _impedance_kernel(
            np.ascontiguousarray(rho, dtype=float),
            np.ascontiguousarray(thickness, dtype=float),
            np.ascontiguousarray(omega, dtype=float),
            float(mu),
            float(thick_limit),
            float(finfo.eps),
            float(finfo.tiny),
        )

else:
    impedance_kernel = None
