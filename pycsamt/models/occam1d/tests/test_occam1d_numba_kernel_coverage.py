# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Line-coverage tests for the numba-JIT'd internals of
:mod:`pycsamt.models.occam1d._numba_kernels`.

``test_occam1d_numba_kernel.py`` (singular) already proves the compiled
kernel is numerically equivalent to the NumPy fallback, but coverage.py
cannot see inside compiled ``@njit`` bytecode -- it only sees the Python
wrapper. This file temporarily flips ``numba.config.DISABLE_JIT`` and
reloads the module so ``@njit`` becomes a plain-Python passthrough,
letting coverage.py trace ``_ctanh`` and ``_impedance_kernel`` directly.
The module is always reloaded back to its normal (compiled) state
afterwards so it does not leak into the rest of the test session.
"""

from __future__ import annotations

import importlib

import numpy as np
import pytest

from pycsamt.models.occam1d import _numba_kernels as _nk_module

pytestmark = pytest.mark.skipif(
    not _nk_module.HAS_NUMBA, reason="numba is not installed"
)


@pytest.fixture
def interpreted_kernels():
    """Reload ``_numba_kernels`` with JIT compilation disabled.

    Yields the module with ``_ctanh``/``_impedance_kernel``/
    ``impedance_kernel`` running as plain Python functions, then restores
    the compiled dispatchers on teardown.
    """
    import numba

    original = numba.config.DISABLE_JIT
    numba.config.DISABLE_JIT = True
    try:
        importlib.reload(_nk_module)
        assert not hasattr(_nk_module._impedance_kernel, "py_func"), (
            "expected a plain function, not a numba Dispatcher"
        )
        yield _nk_module
    finally:
        numba.config.DISABLE_JIT = original
        importlib.reload(_nk_module)


def _layered_rho_thickness(seed=0, n_layers=6):
    rng = np.random.default_rng(seed)
    thickness = rng.uniform(20.0, 200.0, n_layers - 1)
    resistivity = rng.uniform(1.0, 2000.0, n_layers)
    return resistivity, thickness


class TestInterpretedKernelMatchesCompiled:
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_matches_compiled_dispatcher_across_thin_and_thick_layers(
        self, interpreted_kernels, seed
    ):
        """Frequencies span a wide range so both the ``z.real < thick_limit``
        and ``z.real >= thick_limit`` branches of ``_ctanh`` are exercised
        (thick_limit=50.0, matching pycsamt.models.occam1d.forward)."""
        rho, thickness = _layered_rho_thickness(seed)
        frequency = np.logspace(-3, 4, 41)
        omega = 2.0 * np.pi * frequency
        mu = 4.0e-7 * np.pi
        thick_limit = 50.0

        interpreted, bad_f, bad_l = interpreted_kernels.impedance_kernel(
            rho, thickness, omega, mu, thick_limit
        )
        assert bad_f == -1
        assert bad_l == -1
        assert np.all(np.isfinite(interpreted))

        # cross-check against the true compiled path (JIT re-enabled).
        import numba

        saved = numba.config.DISABLE_JIT
        numba.config.DISABLE_JIT = False
        try:
            importlib.reload(_nk_module)
            compiled, bad_f2, bad_l2 = _nk_module.impedance_kernel(
                rho, thickness, omega, mu, thick_limit
            )
        finally:
            numba.config.DISABLE_JIT = saved
            importlib.reload(_nk_module)
            numba.config.DISABLE_JIT = True
            importlib.reload(_nk_module)

        np.testing.assert_allclose(interpreted, compiled, rtol=1e-10)
        assert (bad_f, bad_l) == (bad_f2, bad_l2)

    def test_two_layer_model_thin_layer_uses_tanh_branch(
        self, interpreted_kernels
    ):
        """A very thin top layer keeps ``propagation*thickness`` below
        ``thick_limit``, exercising the ``e2z`` tanh computation instead of
        the early ``complex(1.0, 0.0)`` return in ``_ctanh``."""
        rho = np.array([100.0, 500.0])
        thickness = np.array([0.05])
        omega = np.array([2.0 * np.pi * 1000.0])
        mu = 4.0e-7 * np.pi
        impedance, bad_f, bad_l = interpreted_kernels.impedance_kernel(
            rho, thickness, omega, mu, 50.0
        )
        assert bad_f == -1 and bad_l == -1
        assert np.isfinite(impedance[0])


class TestSingularDetection:
    def test_dc_frequency_triggers_singular_at_deepest_layer(
        self, interpreted_kernels
    ):
        """omega=0 collapses every intrinsic impedance and transfer factor
        to zero, making the denominator exactly zero on the first (deepest,
        i.e. j=n_layers-2) recursion step."""
        rho = np.array([100.0, 50.0, 200.0, 10.0])
        thickness = np.array([100.0, 300.0, 50.0])
        omega = np.array([1.0, 0.0])
        mu = 4.0e-7 * np.pi

        impedance, bad_frequency, bad_layer = (
            interpreted_kernels.impedance_kernel(
                rho, thickness, omega, mu, 50.0
            )
        )
        assert bad_frequency == 1
        assert bad_layer == len(rho) - 2
        # a value is still written for the singular frequency (last `z`
        # before the inner loop broke), not left uninitialized.
        assert impedance.shape == (2,)

    def test_worst_layer_tracking_does_not_regress_on_tie(
        self, interpreted_kernels
    ):
        """Two frequencies that both fail at the same layer keep the first
        one as bad_frequency (the `layer_at_fail > bad_layer` comparison is
        strict, so a tie must not overwrite it)."""
        rho = np.array([100.0, 50.0])
        thickness = np.array([100.0])
        omega = np.array([0.0, 0.0])
        mu = 4.0e-7 * np.pi

        _, bad_frequency, bad_layer = interpreted_kernels.impedance_kernel(
            rho, thickness, omega, mu, 50.0
        )
        assert bad_frequency == 0
        assert bad_layer == 0

    def test_only_second_frequency_singular_updates_tracking(
        self, interpreted_kernels
    ):
        rho = np.array([100.0, 50.0])
        thickness = np.array([100.0])
        omega = np.array([2.0 * np.pi * 10.0, 0.0])
        mu = 4.0e-7 * np.pi

        impedance, bad_frequency, bad_layer = (
            interpreted_kernels.impedance_kernel(
                rho, thickness, omega, mu, 50.0
            )
        )
        assert bad_frequency == 1
        assert bad_layer == 0
        assert np.isfinite(impedance[0])


class TestCtanhDirect:
    def test_thick_layer_early_return(self, interpreted_kernels):
        result = interpreted_kernels._ctanh(complex(100.0, 5.0), 50.0)
        assert result == complex(1.0, 0.0)

    def test_thin_layer_computes_tanh(self, interpreted_kernels):
        z = complex(0.1, 0.2)
        result = interpreted_kernels._ctanh(z, 50.0)
        import cmath

        expected = cmath.tanh(z)
        assert abs(result - expected) < 1e-10

    def test_boundary_real_part_equals_limit_uses_early_return(
        self, interpreted_kernels
    ):
        result = interpreted_kernels._ctanh(complex(50.0, 0.0), 50.0)
        assert result == complex(1.0, 0.0)


def test_module_reports_no_numba_flag_correctly_after_teardown():
    """Sanity check that the fixture above always restores the real
    compiled dispatcher rather than leaking the interpreted one into
    the rest of the test session."""
    assert _nk_module.HAS_NUMBA is True
    assert type(_nk_module._impedance_kernel).__name__ == "CPUDispatcher"
