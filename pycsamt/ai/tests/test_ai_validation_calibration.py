"""Contracts for :mod:`pycsamt.ai.validation.calibration`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.validation import (
    ReliabilityCurve,
    empirical_coverage,
    predictive_sharpness,
    reliability_curve,
)

_TRUE = np.array([0.0, 0.0, 0.0, 0.0, 10.0])
_MEAN = np.zeros(5)
_STD = np.ones(5)


def test_empirical_coverage_known_case_across_default_levels():
    levels, coverage, n_valid = empirical_coverage(_TRUE, _MEAN, _STD)
    np.testing.assert_allclose(levels, [0.5, 0.8, 0.9, 0.95, 0.99])
    np.testing.assert_allclose(coverage, 0.8)
    assert n_valid == 5


def test_predictive_sharpness_basic():
    assert predictive_sharpness(np.array([1.0, 2.0, 3.0])) == (
        pytest.approx(2.0)
    )


def test_reliability_curve_matches_component_functions():
    curve = reliability_curve(_TRUE, _MEAN, _STD)
    levels, coverage, n_valid = empirical_coverage(_TRUE, _MEAN, _STD)
    assert curve.n_valid == n_valid == 5
    np.testing.assert_allclose(curve.levels, levels)
    np.testing.assert_allclose(curve.coverage, coverage)
    assert curve.sharpness == pytest.approx(predictive_sharpness(_STD))
    assert curve.calibration.kind == "calibration"
    assert curve.shape == (5,)


def test_levels_are_validated():
    with pytest.raises(ValueError, match=r"within \(0, 1\)"):
        empirical_coverage(
            np.zeros(2), np.zeros(2), np.ones(2), levels=[0.0]
        )
    with pytest.raises(ValueError, match=r"within \(0, 1\)"):
        empirical_coverage(
            np.zeros(2), np.zeros(2), np.ones(2), levels=[1.0]
        )
    with pytest.raises(ValueError, match="non-empty 1-D array"):
        empirical_coverage(
            np.zeros(2), np.zeros(2), np.ones(2), levels=[]
        )


def test_non_positive_and_nan_std_are_excluded():
    true = np.array([0.0, 0.0, 0.0])
    mean = np.array([0.0, 0.0, 0.0])
    std = np.array([1.0, 0.0, np.nan])
    _, _, n_valid = empirical_coverage(true, mean, std, levels=[0.5])
    assert n_valid == 1


def test_shape_and_empty_checks_raise():
    with pytest.raises(ValueError, match="share one shape"):
        empirical_coverage(np.zeros(2), np.zeros(3), np.ones(2))
    with pytest.raises(ValueError, match="not be empty"):
        empirical_coverage(np.array([]), np.array([]), np.array([]))
    with pytest.raises(ValueError, match="not be empty"):
        predictive_sharpness(np.array([]))


def test_valid_mask_shape_mismatch_raises():
    with pytest.raises(ValueError, match="same shape as y_true"):
        empirical_coverage(
            np.zeros(3),
            np.zeros(3),
            np.ones(3),
            valid=np.ones(2, dtype=bool),
        )
    with pytest.raises(ValueError, match="same shape as y_pred_std"):
        predictive_sharpness(np.ones(3), valid=np.ones(2, dtype=bool))


def test_no_valid_cell_raises_for_every_entry_point():
    empty_mask = np.zeros(3, dtype=bool)
    with pytest.raises(ValueError, match="no valid cell"):
        empirical_coverage(
            np.zeros(3), np.zeros(3), np.ones(3), valid=empty_mask
        )
    with pytest.raises(ValueError, match="no valid cell"):
        predictive_sharpness(np.ones(3), valid=empty_mask)
    with pytest.raises(ValueError, match="no valid cell"):
        reliability_curve(
            np.zeros(3), np.zeros(3), np.ones(3), valid=empty_mask
        )


def test_reliability_curve_validates_shape_and_type():
    good = reliability_curve(
        np.array([0.0, 10.0]), np.zeros(2), np.ones(2), levels=[0.5]
    )
    with pytest.raises(ValueError, match="share one shape"):
        ReliabilityCurve(
            levels=np.array([0.5, 0.8]),
            coverage=np.array([0.5]),
            calibration=good.calibration,
            sharpness=1.0,
            n_valid=2,
            shape=(2,),
        )
    with pytest.raises(TypeError, match="UncertaintyLossResult"):
        ReliabilityCurve(
            levels=np.array([0.5]),
            coverage=np.array([0.5]),
            calibration="not a result",
            sharpness=1.0,
            n_valid=2,
            shape=(2,),
        )
