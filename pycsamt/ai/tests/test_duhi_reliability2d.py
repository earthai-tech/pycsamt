"""Tests for DUHI observation-reliability factors."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.inversion import (
    combine_observation_reliability,
    dimensionality_reliability,
)


def test_dimensionality_reliability_matches_frozen_equation():
    actual = dimensionality_reliability([0.0, 5.0, 10.0])
    expected = np.exp(-np.square([0.0, 1.0, 2.0]))
    np.testing.assert_allclose(actual, expected)


def test_dimensionality_reliability_handles_nan_and_floor():
    actual = dimensionality_reliability(
        [np.nan, 100.0], beta_scale_deg=5.0, minimum=0.01
    )
    np.testing.assert_allclose(actual, [0.01, 0.01])


def test_combined_reliability_broadcasts_and_applies_floor():
    actual = combine_observation_reliability(
        [[0.8], [0.5]],
        [0.5, 0.2],
        minimum=0.11,
    )
    np.testing.assert_allclose(actual, [[0.4, 0.16], [0.25, 0.11]])


@pytest.mark.parametrize(
    "function,args,match",
    [
        (dimensionality_reliability, ([],), "at least one"),
        (dimensionality_reliability, ([0.0],), "beta_scale_deg"),
        (
            combine_observation_reliability,
            ([1.2], [1.0]),
            "measurement",
        ),
        (
            combine_observation_reliability,
            ([1.0, 0.5], [1.0, 0.5, 0.2]),
            "broadcastable",
        ),
    ],
)
def test_reliability_functions_reject_invalid_inputs(
    function, args, match
):
    kwargs = {"beta_scale_deg": 0.0} if match == "beta_scale_deg" else {}
    with pytest.raises(ValueError, match=match):
        function(*args, **kwargs)
