"""Contracts for :mod:`pycsamt.ai.losses.boundary`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.losses import (
    BoundaryLoss,
    ModelLossResult,
    boundary_condition_loss,
)

_GRID = np.array([[1.0, 5.0], [0.0, 0.0]])
_AIR = np.array([[True, True], [False, False]])


def test_l1_l2_and_huber_penalties_on_boundary_cells():
    l1 = boundary_condition_loss(
        _GRID, boundary_mask=_AIR, target=0.0, kind="l1"
    )
    l2 = boundary_condition_loss(
        _GRID, boundary_mask=_AIR, target=0.0, kind="l2"
    )
    huber = boundary_condition_loss(
        _GRID, boundary_mask=_AIR, target=0.0, kind="huber", delta=1.0
    )
    assert isinstance(l1, ModelLossResult)
    assert l1.value == pytest.approx(3.0)
    assert l2.value == pytest.approx(13.0)
    assert huber.value == pytest.approx(2.5)
    assert l1.n_valid == 2
    assert l1.weight_sum == pytest.approx(2.0)


def test_scalar_and_array_targets_broadcast():
    scalar = boundary_condition_loss(
        _GRID, boundary_mask=_AIR, target=0.0, kind="l1"
    )
    array_target = boundary_condition_loss(
        _GRID,
        boundary_mask=_AIR,
        target=np.array([[0.0, 0.0], [0.0, 0.0]]),
        kind="l1",
    )
    assert scalar.value == pytest.approx(array_target.value)


def test_extra_valid_mask_further_restricts_boundary_cells():
    result = boundary_condition_loss(
        _GRID,
        boundary_mask=_AIR,
        target=0.0,
        kind="l1",
        valid=np.array([[True, False], [True, True]]),
    )
    assert result.n_valid == 1
    assert result.value == pytest.approx(1.0)


def test_weights_scale_boundary_cell_contribution():
    result = boundary_condition_loss(
        _GRID,
        boundary_mask=_AIR,
        target=0.0,
        kind="l1",
        weights=np.array([[2.0, 1.0], [1.0, 1.0]]),
    )
    assert result.weight_sum == pytest.approx(3.0)
    assert result.value == pytest.approx(7.0 / 3.0)


def test_boundary_mask_must_match_shape_and_select_a_cell():
    with pytest.raises(ValueError, match="same shape as y_pred"):
        boundary_condition_loss(
            _GRID, boundary_mask=np.zeros((3, 3), dtype=bool), target=0.0
        )
    with pytest.raises(ValueError, match="at least one cell"):
        boundary_condition_loss(
            _GRID, boundary_mask=np.zeros((2, 2), dtype=bool), target=0.0
        )


def test_bad_target_raises():
    with pytest.raises(ValueError, match="broadcastable"):
        boundary_condition_loss(
            _GRID, boundary_mask=_AIR, target=np.array([1.0, 2.0, 3.0])
        )
    with pytest.raises(ValueError, match="finite"):
        boundary_condition_loss(_GRID, boundary_mask=_AIR, target=np.nan)


def test_bad_valid_shape_raises():
    with pytest.raises(ValueError, match="same shape as y_pred"):
        boundary_condition_loss(
            _GRID,
            boundary_mask=_AIR,
            target=0.0,
            valid=np.zeros((3, 3), dtype=bool),
        )


def test_boundary_loss_is_callable_and_reusable():
    loss = BoundaryLoss(kind="l1")
    result = loss(_GRID, boundary_mask=_AIR, target=0.0)
    assert result.value == pytest.approx(3.0)
    with pytest.raises(ValueError, match="kind"):
        BoundaryLoss(kind="bogus")
    with pytest.raises(ValueError, match="reduction"):
        BoundaryLoss(reduction="bogus")
