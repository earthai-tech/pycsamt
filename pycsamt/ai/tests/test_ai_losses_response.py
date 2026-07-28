"""Contracts for :mod:`pycsamt.ai.losses.response`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.data import SurveyData
from pycsamt.ai.losses import (
    ResponseLoss,
    ResponseLossResult,
    response_loss_from_contracts,
    response_residual_loss,
)
from pycsamt.forward.maxwell import ForwardResult, SolverDiagnostics


def test_l1_and_l2_residuals_on_clean_complex_inputs():
    pred = np.array([1 + 1j, 2 + 2j])
    obs = np.array([1 + 1j, 0 + 0j])
    l2 = response_residual_loss(pred, obs, kind="l2")
    l1 = response_residual_loss(pred, obs, kind="l1")
    assert l2.value == pytest.approx(4.0)
    assert l1.value == pytest.approx(np.sqrt(8.0) / 2.0)
    assert l2.normalized is False
    assert l2.n_valid == 2
    assert l2.weight_sum == pytest.approx(2.0)


def test_errors_normalize_residuals_before_the_penalty():
    pred = np.array([1 + 1j, 2 + 2j])
    obs = np.array([1 + 1j, 0 + 0j])
    result = response_residual_loss(
        pred, obs, errors=np.array([1.0, 2.0]), kind="l2"
    )
    assert result.normalized is True
    assert result.value == pytest.approx(1.0)


def test_bad_errors_are_excluded_not_raised():
    pred = np.array([1 + 1j, 2 + 2j])
    obs = np.array([1 + 1j, 0 + 0j])
    result = response_residual_loss(
        pred, obs, errors=np.array([1.0, 0.0]), kind="l2"
    )
    assert result.n_valid == 1
    assert result.value == pytest.approx(0.0)

    with_nan_error = response_residual_loss(
        pred, obs, errors=np.array([1.0, np.nan]), kind="l2"
    )
    assert with_nan_error.n_valid == 1


def test_nan_and_valid_mask_exclude_cells():
    pred = np.array([1 + 1j, np.nan + 1j])
    obs = np.array([1 + 1j, 2 + 2j])
    result = response_residual_loss(pred, obs)
    assert result.n_valid == 1
    assert result.value == pytest.approx(0.0)

    masked = response_residual_loss(
        np.array([1 + 1j, 5 + 5j]),
        np.array([1 + 1j, 0 + 0j]),
        valid=np.array([True, False]),
    )
    assert masked.n_valid == 1
    assert masked.value == pytest.approx(0.0)


def test_shape_and_empty_mismatches_raise():
    with pytest.raises(ValueError, match="share one shape"):
        response_residual_loss(np.array([1 + 1j]), np.array([1 + 1j, 2j]))
    with pytest.raises(ValueError, match="not be empty"):
        response_residual_loss(np.array([]), np.array([]))
    with pytest.raises(ValueError, match="errors must have"):
        response_residual_loss(
            np.array([1 + 1j]), np.array([1 + 1j]), errors=np.array([1, 2])
        )


def test_response_loss_result_rejects_invalid_state():
    with pytest.raises(ValueError, match="kind"):
        ResponseLossResult(
            value=0.0,
            kind="bogus",
            reduction="mean",
            n_valid=1,
            weight_sum=1.0,
            normalized=False,
        )


def test_response_loss_is_callable_and_reusable():
    loss = ResponseLoss(kind="l2", reduction="sum")
    pred = np.array([1 + 1j, 2 + 2j])
    obs = np.array([0 + 0j, 0 + 0j])
    result = loss(pred, obs)
    assert isinstance(result, ResponseLossResult)
    assert result.value == pytest.approx(2.0 + 8.0)
    with pytest.raises(ValueError, match="kind"):
        ResponseLoss(kind="bogus")


def _forward_result(
    impedance, *, receiver_names=("S1",), components=("zxy",),
    frequencies=(10.0,),
):
    diagnostics = SolverDiagnostics(
        [[True]] * len(frequencies), [[1]] * len(frequencies),
        [[0.0]] * len(frequencies), 0.01,
    )
    return ForwardResult(
        "a" * 64,
        list(frequencies),
        list(receiver_names),
        list(components),
        impedance,
        None,
        "demo",
        "1",
        diagnostics,
    )


def test_response_loss_from_contracts_matches_raw_arrays():
    z = np.array([[[1 + 1j]]])
    observed = SurveyData(z, [10.0], ["S1"], ["zxy"], [[0.0, 0.0]])
    forward = _forward_result(z)
    result = response_loss_from_contracts(forward, observed)
    assert result.value == pytest.approx(0.0)
    assert result.n_valid == 1


def test_response_loss_from_contracts_uses_survey_errors():
    predicted = np.array([[[2 + 2j]]])
    observed = SurveyData(
        np.array([[[0 + 0j]]]),
        [10.0],
        ["S1"],
        ["zxy"],
        [[0.0, 0.0]],
        impedance_error=np.array([[[2.0]]]),
    )
    forward = _forward_result(predicted)
    result = response_loss_from_contracts(forward, observed)
    assert result.normalized is True
    assert result.value == pytest.approx(2.0)
    unnormalized = response_loss_from_contracts(
        forward, observed, use_errors=False
    )
    assert unnormalized.normalized is False
    assert unnormalized.value == pytest.approx(8.0)


def test_response_loss_from_contracts_rejects_type_and_misalignment():
    z = np.array([[[1 + 1j]]])
    observed = SurveyData(z, [10.0], ["S1"], ["zxy"], [[0.0, 0.0]])
    forward = _forward_result(z)
    with pytest.raises(TypeError, match="ForwardResult"):
        response_loss_from_contracts(observed, observed)
    with pytest.raises(TypeError, match="SurveyData"):
        response_loss_from_contracts(forward, forward)

    other_station = _forward_result(z, receiver_names=("S2",))
    with pytest.raises(ValueError, match="station_names"):
        response_loss_from_contracts(other_station, observed)

    other_component = _forward_result(z, components=("zyx",))
    with pytest.raises(ValueError, match="components"):
        response_loss_from_contracts(other_component, observed)

    other_frequency = _forward_result(z, frequencies=(20.0,))
    with pytest.raises(ValueError, match="frequencies_hz"):
        response_loss_from_contracts(other_frequency, observed)
