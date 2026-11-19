# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regression: phase-tensor skew follows Caldwell et al. (2004)."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.tensor import _angles_deg


@pytest.mark.parametrize(
    "a, b, c, d",
    [(1.0, 0.3, -0.1, 0.8),      # normal tensor, trace > 0
     (-1.0, 0.3, -0.1, -0.8),    # trace < 0 (out-of-quadrant phases)
     (-0.2, 1.0, 0.4, 0.1)],
)
def test_beta_matches_arctan_definition(a, b, c, d):
    _, beta = _angles_deg(*(np.array([v]) for v in (a, b, c, d)))
    expect = 0.5 * np.degrees(np.arctan((b - c) / (a + d)))
    assert -45.0 < beta[0] <= 45.0
    np.testing.assert_allclose(beta, [expect], atol=1e-12)


def test_beta_zero_over_zero_is_zero():
    z = np.array([0.0])
    _, beta = _angles_deg(np.array([1.0]), z, z, np.array([-1.0]))
    assert beta[0] == 0.0


def test_pure_skew_is_plus_45():
    one = np.array([1.0])
    _, beta = _angles_deg(0 * one, one, -one, 0 * one)
    assert beta[0] == pytest.approx(45.0)
