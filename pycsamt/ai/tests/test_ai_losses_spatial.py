"""Contracts for :mod:`pycsamt.ai.losses.spatial`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.losses import (
    SpatialLoss,
    SpatialLossResult,
    gradient_smoothness_loss,
    total_variation_loss,
)


def test_gradient_smoothness_l1_and_l2_along_an_axis():
    grid = np.array([[0.0, 1.0, 3.0], [0.0, 0.0, 0.0]])
    l1 = gradient_smoothness_loss(grid, axis=1, kind="l1")
    l2 = gradient_smoothness_loss(grid, axis=1, kind="l2")
    assert l1.value == pytest.approx(0.75)
    assert l2.value == pytest.approx((1 + 4) / 4)
    assert l1.label == "grad_axis1"
    assert l1.n_valid == 4


def test_negative_axis_is_normalized_in_the_label():
    grid = np.zeros((3, 4))
    result = gradient_smoothness_loss(grid, axis=-1)
    assert result.label == "grad_axis1"


def test_axis_out_of_range_and_too_short_raise():
    grid = np.zeros((3, 4))
    with pytest.raises(ValueError, match="axis must be in"):
        gradient_smoothness_loss(grid, axis=5)
    with pytest.raises(ValueError, match="at least two cells"):
        gradient_smoothness_loss(np.zeros((1, 4)), axis=0)


def test_invalid_kind_and_reduction_raise():
    grid = np.zeros((2, 2))
    with pytest.raises(ValueError, match="kind"):
        gradient_smoothness_loss(grid, axis=0, kind="huber")
    with pytest.raises(ValueError, match="reduction"):
        SpatialLossResult(
            value=0.0,
            kind="l1",
            label="tv",
            reduction="bogus",
            n_valid=1,
            weight_sum=1.0,
        )


def test_nan_and_valid_mask_exclude_touching_differences():
    grid = np.array([1.0, np.nan, 3.0, 4.0])
    result = gradient_smoothness_loss(grid, axis=0, kind="l1")
    # pairs (0,1) and (1,2) are excluded by the NaN at index 1;
    # only pair (2, 3) survives.
    assert result.n_valid == 1
    assert result.value == pytest.approx(1.0)

    masked = gradient_smoothness_loss(
        np.array([1.0, 5.0, 3.0]),
        axis=0,
        kind="l1",
        valid=np.array([True, False, True]),
    )
    assert masked.n_valid == 0
    assert np.isnan(masked.value)


def test_weights_use_the_minimum_of_the_two_endpoints():
    grid = np.array([0.0, 2.0, 10.0])
    result = gradient_smoothness_loss(
        grid, axis=0, kind="l1", weights=np.array([1.0, 0.0, 1.0])
    )
    assert result.n_valid == 2
    assert result.weight_sum == pytest.approx(0.0)
    assert np.isnan(result.value)
    with pytest.raises(ValueError, match="non-negative"):
        gradient_smoothness_loss(
            grid, axis=0, weights=np.array([-1.0, 1.0, 1.0])
        )


def test_sum_reduction_matches_mean_times_weight():
    grid = np.array([0.0, 1.0, 5.0])
    mean_result = gradient_smoothness_loss(grid, axis=0, kind="l1")
    sum_result = gradient_smoothness_loss(
        grid, axis=0, kind="l1", reduction="sum"
    )
    assert sum_result.value == pytest.approx(
        mean_result.value * mean_result.weight_sum
    )


def test_total_variation_combines_every_axis():
    grid = np.array([[0.0, 1.0], [0.0, 3.0]])
    result = total_variation_loss(grid)
    assert result.label == "tv"
    assert result.n_valid == 4
    assert result.value == pytest.approx(1.5)
    summed = total_variation_loss(grid, reduction="sum")
    assert summed.value == pytest.approx(6.0)


def test_empty_and_scalar_inputs_raise():
    with pytest.raises(ValueError, match="non-empty"):
        gradient_smoothness_loss(np.array([]), axis=0)
    with pytest.raises(ValueError, match="non-empty"):
        total_variation_loss(np.array(1.0))


def test_spatial_loss_combines_weighted_terms():
    grid = np.array([[0.0, 1.0], [0.0, 3.0]])
    loss = SpatialLoss(
        lambda_x=1.0, lambda_z=0.0, lambda_tv=0.0, kind="l1"
    )
    assert loss(grid) == pytest.approx(2.0)

    both = SpatialLoss(lambda_x=1.0, lambda_z=1.0, kind="l1")
    grad_x = gradient_smoothness_loss(grid, axis=1, kind="l1").value
    grad_z = gradient_smoothness_loss(grid, axis=0, kind="l1").value
    assert both(grid) == pytest.approx(grad_x + grad_z)


def test_spatial_loss_rejects_bad_configuration_and_shape():
    with pytest.raises(ValueError, match="kind"):
        SpatialLoss(kind="bogus")
    with pytest.raises(ValueError, match="non-negative"):
        SpatialLoss(lambda_x=-1.0)
    with pytest.raises(ValueError, match="2-D"):
        SpatialLoss()(np.zeros((2, 2, 2)))
