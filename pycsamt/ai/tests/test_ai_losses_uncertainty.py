"""Contracts for :mod:`pycsamt.ai.losses.uncertainty`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.losses import (
    UncertaintyLoss,
    UncertaintyLossResult,
    calibration_loss,
    gaussian_nll_loss,
)


def test_gaussian_nll_matches_the_closed_form():
    pred = np.array([0.0, 1.0])
    true = np.array([0.0, 0.0])
    log_var = np.array([0.0, 0.0])
    result = gaussian_nll_loss(pred, true, log_var)
    per_cell = 0.5 * (
        np.square(pred - true) + log_var + np.log(2.0 * np.pi)
    )
    assert result.value == pytest.approx(float(np.mean(per_cell)))
    assert result.kind == "gaussian_nll"
    assert result.n_valid == 2
    assert result.weight_sum == pytest.approx(2.0)


def test_gaussian_nll_is_minimized_at_the_analytic_optimal_variance():
    diff = 3.0
    pred = np.array([diff])
    true = np.array([0.0])
    optimal_log_var = np.array([2.0 * np.log(abs(diff))])
    result = gaussian_nll_loss(pred, true, optimal_log_var)
    nearby = gaussian_nll_loss(pred, true, optimal_log_var + 0.5)
    expected = 0.5 * (1.0 + optimal_log_var[0] + np.log(2.0 * np.pi))
    assert result.value == pytest.approx(expected)
    assert result.value < nearby.value


def test_gaussian_nll_shape_and_empty_checks_raise():
    with pytest.raises(ValueError, match="share one shape"):
        gaussian_nll_loss(np.zeros(2), np.zeros(3), np.zeros(2))
    with pytest.raises(ValueError, match="not be empty"):
        gaussian_nll_loss(np.array([]), np.array([]), np.array([]))


def test_gaussian_nll_excludes_nan_and_respects_valid_mask():
    pred = np.array([0.0, np.nan])
    true = np.array([0.0, 0.0])
    log_var = np.array([0.0, 0.0])
    result = gaussian_nll_loss(pred, true, log_var)
    assert result.n_valid == 1

    masked = gaussian_nll_loss(
        np.array([0.0, 1.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        valid=np.array([True, False]),
    )
    assert masked.n_valid == 1
    assert masked.value == pytest.approx(0.5 * np.log(2.0 * np.pi))


def test_gaussian_nll_weights_scale_the_reduced_value():
    pred = np.array([0.0, 1.0])
    true = np.array([0.0, 0.0])
    log_var = np.array([0.0, 0.0])
    per_cell = 0.5 * (
        np.square(pred - true) + log_var + np.log(2.0 * np.pi)
    )
    unweighted = gaussian_nll_loss(pred, true, log_var, reduction="sum")
    weighted = gaussian_nll_loss(
        pred, true, log_var, weights=np.array([2.0, 1.0]), reduction="sum"
    )
    assert unweighted.value == pytest.approx(float(np.sum(per_cell)))
    expected = float(np.sum(per_cell * np.array([2.0, 1.0])))
    assert weighted.value == pytest.approx(expected)
    with pytest.raises(ValueError, match="non-negative"):
        gaussian_nll_loss(
            pred, true, log_var, weights=np.array([-1.0, 1.0])
        )


def test_calibration_l1_and_l2_penalties():
    coverage = np.array([0.4, 0.9])
    nominal = np.array([0.5, 0.8])
    l1 = calibration_loss(coverage, nominal, kind="l1")
    l2 = calibration_loss(coverage, nominal, kind="l2")
    assert l1.value == pytest.approx(0.1)
    assert l2.value == pytest.approx(0.01)
    assert l1.kind == "calibration"
    with pytest.raises(ValueError, match="kind"):
        calibration_loss(coverage, nominal, kind="huber")


def test_calibration_rejects_out_of_range_values():
    with pytest.raises(ValueError, match="coverage must be within"):
        calibration_loss(np.array([1.5, 0.5]), np.array([0.5, 0.5]))
    with pytest.raises(ValueError, match="nominal_levels must be within"):
        calibration_loss(np.array([0.5, 0.5]), np.array([-0.1, 0.5]))


def test_calibration_shape_and_empty_checks_raise():
    with pytest.raises(ValueError, match="share one shape"):
        calibration_loss(np.array([0.5, 0.5]), np.array([0.5]))
    with pytest.raises(ValueError, match="not be empty"):
        calibration_loss(np.array([]), np.array([]))


def test_calibration_excludes_nan_and_respects_valid_mask():
    result = calibration_loss(
        np.array([0.5, np.nan]), np.array([0.5, 0.5])
    )
    assert result.n_valid == 1
    assert result.value == pytest.approx(0.0)

    masked = calibration_loss(
        np.array([0.4, 0.9]),
        np.array([0.5, 0.8]),
        kind="l1",
        valid=np.array([True, False]),
    )
    assert masked.n_valid == 1
    assert masked.value == pytest.approx(0.1)


def test_calibration_weights_scale_the_reduced_value():
    coverage = np.array([0.4, 0.9])
    nominal = np.array([0.5, 0.7])
    result = calibration_loss(
        coverage,
        nominal,
        kind="l1",
        weights=np.array([2.0, 1.0]),
        reduction="sum",
    )
    assert result.value == pytest.approx(0.1 * 2.0 + 0.2 * 1.0)


def test_uncertainty_loss_result_rejects_invalid_state():
    with pytest.raises(ValueError, match="kind"):
        UncertaintyLossResult(
            value=0.0,
            kind="bogus",
            reduction="mean",
            n_valid=1,
            weight_sum=1.0,
        )
    with pytest.raises(ValueError, match="reduction"):
        UncertaintyLossResult(
            value=0.0,
            kind="gaussian_nll",
            reduction="bogus",
            n_valid=1,
            weight_sum=1.0,
        )


def test_uncertainty_loss_is_callable_and_reusable():
    loss = UncertaintyLoss(reduction="sum")
    pred = np.array([0.0, 1.0])
    true = np.array([0.0, 0.0])
    log_var = np.array([0.0, 0.0])
    result = loss(pred, true, log_var)
    assert isinstance(result, UncertaintyLossResult)
    direct = gaussian_nll_loss(pred, true, log_var, reduction="sum")
    assert result.value == pytest.approx(direct.value)
    with pytest.raises(ValueError, match="reduction"):
        UncertaintyLoss(reduction="bogus")
