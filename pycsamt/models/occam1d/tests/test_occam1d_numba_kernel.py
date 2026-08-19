"""Equivalence contract between the optional Numba kernel and NumPy path.

:mod:`pycsamt.models.occam1d._numba_kernels` is a soft dependency: when
Numba is installed, :meth:`Occam1DForwardModel.impedance` dispatches to a
compiled per-frequency kernel instead of the NumPy loop vectorized over
frequency. Both must produce the same physics; this file is the contract
that keeps them in sync as either implementation changes.
"""

import numpy as np
import pytest

from pycsamt.models.occam1d import _numba_kernels
from pycsamt.models.occam1d.forward import Occam1DForwardModel
from pycsamt.models.occam1d.model import Occam1DModel

pytestmark = pytest.mark.skipif(
    not _numba_kernels.HAS_NUMBA, reason="numba is not installed"
)


def _layered_model(seed=0, n_layers=12):
    rng = np.random.default_rng(seed)
    depth = np.concatenate([[0.0], np.cumsum(rng.uniform(20.0, 200.0, n_layers - 1))])
    resistivity = rng.uniform(1.0, 2000.0, n_layers)
    return Occam1DModel(depth, resistivity)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_numba_kernel_matches_numpy_recursion(seed):
    model = _layered_model(seed)
    forward = Occam1DForwardModel(model)
    frequency = np.logspace(-3, 4, 41)
    thickness = np.diff(model.depth)
    omega = 2.0 * np.pi * frequency
    numba_impedance = forward._recurse_numba(
        frequency, model.resistivity, thickness, omega
    )
    numpy_impedance = forward._recurse_numpy(
        frequency, model.resistivity, thickness, omega
    )
    np.testing.assert_allclose(numba_impedance, numpy_impedance, rtol=1e-10)


def test_public_impedance_call_uses_the_kernel_transparently():
    model = _layered_model(seed=3)
    forward = Occam1DForwardModel(model)
    frequency = np.logspace(-2, 3, 25)
    dispatched = forward.impedance(frequency)
    forced_numpy = forward._recurse_numpy(
        frequency,
        model.resistivity,
        np.diff(model.depth),
        2.0 * np.pi * frequency,
    )
    np.testing.assert_allclose(dispatched, forced_numpy, rtol=1e-10)
