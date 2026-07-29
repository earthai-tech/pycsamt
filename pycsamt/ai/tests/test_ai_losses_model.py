"""Contracts for :mod:`pycsamt.ai.losses.model`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.losses import (
    ModelLoss,
    ModelLossResult,
    depth_weights,
    model_huber_loss,
    model_l1_loss,
    model_l2_loss,
)


def test_l1_and_l2_losses_on_clean_inputs():
    pred = np.array([1.0, 3.0, -2.0])
    true = np.array([1.0, 1.0, -2.0])
    l1 = model_l1_loss(pred, true)
    l2 = model_l2_loss(pred, true)
    assert l1.value == pytest.approx(2.0 / 3.0)
    assert l2.value == pytest.approx(4.0 / 3.0)
    assert l1.kind == "l1"
    assert l2.kind == "l2"
    assert l1.n_valid == 3
    assert l1.weight_sum == pytest.approx(3.0)


def test_sum_reduction_matches_mean_times_weight():
    pred = np.array([1.0, 3.0])
    true = np.array([1.0, 1.0])
    mean_result = model_l1_loss(pred, true, reduction="mean")
    sum_result = model_l1_loss(pred, true, reduction="sum")
    assert sum_result.value == pytest.approx(
        mean_result.value * mean_result.weight_sum
    )


def test_huber_loss_matches_quadratic_and_linear_regimes():
    small = model_huber_loss(
        np.array([0.5]), np.array([0.0]), delta=1.0
    )
    large = model_huber_loss(
        np.array([5.0]), np.array([0.0]), delta=1.0
    )
    assert small.value == pytest.approx(0.125)
    assert large.value == pytest.approx(4.5)
    with pytest.raises(ValueError, match="delta"):
        model_huber_loss(np.array([1.0]), np.array([0.0]), delta=0.0)


def test_nan_and_explicit_mask_exclude_cells():
    pred = np.array([1.0, np.nan, 5.0])
    true = np.array([1.0, 2.0, 0.0])
    result = model_l1_loss(pred, true)
    assert result.n_valid == 2
    assert result.value == pytest.approx(2.5)

    masked = model_l1_loss(
        np.array([1.0, 5.0]),
        np.array([1.0, 0.0]),
        valid=np.array([True, False]),
    )
    assert masked.n_valid == 1
    assert masked.value == pytest.approx(0.0)


def test_weights_scale_contribution_and_reject_negative():
    pred = np.array([1.0, 3.0])
    true = np.array([0.0, 0.0])
    weighted = model_l1_loss(pred, true, weights=np.array([1.0, 0.0]))
    assert weighted.value == pytest.approx(1.0)
    assert weighted.weight_sum == pytest.approx(1.0)
    with pytest.raises(ValueError, match="non-negative"):
        model_l1_loss(pred, true, weights=np.array([-1.0, 0.0]))


def test_mismatched_shapes_and_empty_inputs_raise():
    with pytest.raises(ValueError, match="share one shape"):
        model_l1_loss(np.array([1.0, 2.0]), np.array([1.0]))
    with pytest.raises(ValueError, match="not be empty"):
        model_l1_loss(np.array([]), np.array([]))


def test_model_loss_result_rejects_invalid_state():
    with pytest.raises(ValueError, match="kind"):
        ModelLossResult(
            value=0.0,
            kind="bogus",
            reduction="mean",
            n_valid=1,
            weight_sum=1.0,
        )
    with pytest.raises(ValueError, match="reduction"):
        ModelLossResult(
            value=0.0,
            kind="l1",
            reduction="bogus",
            n_valid=1,
            weight_sum=1.0,
        )


def test_depth_weights_are_normalized_and_decreasing():
    weights = depth_weights(4)
    assert weights.sum() == pytest.approx(1.0)
    assert np.all(np.diff(weights) < 0)
    with pytest.raises(ValueError, match="positive integer"):
        depth_weights(0)


def test_model_loss_is_callable_and_reusable():
    loss = ModelLoss(kind="l2", reduction="sum")
    result = loss(np.array([1.0, 2.0]), np.array([0.0, 0.0]))
    assert isinstance(result, ModelLossResult)
    assert result.value == pytest.approx(5.0)
    with pytest.raises(ValueError, match="kind"):
        ModelLoss(kind="bogus")


def test_with_depth_weights_broadcasts_over_a_2d_grid():
    loss = ModelLoss.with_depth_weights(3, dimension=2, kind="l1")
    assert loss.weights.shape == (3, 1)
    grid_pred = np.zeros((3, 2))
    grid_true = np.ones((3, 2))
    result = loss(grid_pred, grid_true)
    assert result.n_valid == 6
    assert result.value == pytest.approx(1.0)
    shallow_only = loss(
        grid_pred, grid_true, weights=np.array([[1.0], [0.0], [0.0]])
    )
    assert shallow_only.value == pytest.approx(1.0)
    assert shallow_only.weight_sum == pytest.approx(2.0)


def test_with_depth_weights_rejects_bad_dimension():
    with pytest.raises(ValueError, match="dimension"):
        ModelLoss.with_depth_weights(3, dimension=4)
