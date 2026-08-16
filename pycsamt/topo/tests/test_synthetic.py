# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.topo.synthetic."""

from __future__ import annotations

import numpy as np

from pycsamt.topo.synthetic import synthetic_elevation_profile


def test_deterministic_and_shape():
    chainage = np.linspace(0, 2400, 25)
    a = synthetic_elevation_profile(chainage)
    b = synthetic_elevation_profile(chainage)
    assert a.shape == chainage.shape
    np.testing.assert_array_equal(a, b)


def test_base_and_amplitude_bound_the_range():
    chainage = np.linspace(0, 5000, 500)
    elev = synthetic_elevation_profile(chainage, base_m=100.0, amplitude_m=30.0)
    # dominant + 0.4x secondary sinusoid: total range is bounded by
    # base +/- amplitude*(1 + 0.4)
    assert elev.min() >= 100.0 - 30.0 * 1.4 - 1e-6
    assert elev.max() <= 100.0 + 30.0 * 1.4 + 1e-6


def test_phase_m_shifts_without_changing_shape():
    chainage = np.linspace(0, 2400, 25)
    a = synthetic_elevation_profile(chainage, phase_m=0.0)
    b = synthetic_elevation_profile(chainage + 500.0, phase_m=0.0)
    c = synthetic_elevation_profile(chainage, phase_m=500.0)
    np.testing.assert_allclose(b, c)
    assert not np.allclose(a, c)


def test_scalar_input():
    val = synthetic_elevation_profile(0.0)
    assert np.isscalar(val) or val.shape == ()
