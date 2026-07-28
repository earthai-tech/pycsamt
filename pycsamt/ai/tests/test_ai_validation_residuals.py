"""Contracts for :mod:`pycsamt.ai.validation.residuals`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.data import SurveyData
from pycsamt.ai.losses import response_residual_loss
from pycsamt.ai.validation import (
    ResponseResidualReport,
    response_residual_report,
    response_residual_report_from_contracts,
)
from pycsamt.forward.maxwell import ForwardResult, SolverDiagnostics

_PRED = np.array([[[1 + 0j], [2 + 0j]], [[0j], [0j]]])
_OBS = np.array([[[1 + 0j], [0j]], [[0j], [3 + 0j]]])
_OVERALL = response_residual_loss(np.array([1 + 0j]), np.array([1 + 0j]))


def test_report_matches_response_residual_loss_overall():
    report = response_residual_report(_PRED, _OBS)
    direct = response_residual_loss(_PRED, _OBS)
    assert report.overall.value == pytest.approx(direct.value)
    assert report.overall.n_valid == 4
    assert report.shape == (2, 2, 1)


def test_report_l2_axis_breakdowns():
    report = response_residual_report(_PRED, _OBS, kind="l2")
    np.testing.assert_allclose(report.by_station, [2.0, 4.5])
    np.testing.assert_allclose(report.by_frequency, [0.0, 6.5])
    np.testing.assert_allclose(report.by_component, [3.25])
    assert report.overall.value == pytest.approx(3.25)


def test_report_l1_axis_breakdowns():
    report = response_residual_report(_PRED, _OBS, kind="l1")
    np.testing.assert_allclose(report.by_station, [1.0, 1.5])
    np.testing.assert_allclose(report.by_frequency, [0.0, 2.5])
    np.testing.assert_allclose(report.by_component, [1.25])


def test_report_valid_mask_can_exclude_a_whole_station():
    valid = np.ones((2, 2, 1), dtype=bool)
    valid[1, :, :] = False
    report = response_residual_report(_PRED, _OBS, valid=valid)
    assert report.by_station[0] == pytest.approx(2.0)
    assert np.isnan(report.by_station[1])
    np.testing.assert_allclose(report.by_frequency, [0.0, 4.0])
    assert report.overall.n_valid == 2
    assert report.overall.value == pytest.approx(2.0)


def test_report_errors_normalize_consistently_with_the_loss():
    errors = np.array([[[1.0], [2.0]], [[1.0], [1.0]]])
    report = response_residual_report(_PRED, _OBS, errors=errors)
    direct = response_residual_loss(_PRED, _OBS, errors=errors)
    assert report.overall.normalized is True
    assert report.overall.value == pytest.approx(direct.value)


def test_report_shape_ndim_and_empty_checks_raise():
    with pytest.raises(ValueError, match="share one shape"):
        response_residual_report(
            _PRED, np.zeros((2, 2, 2), dtype=complex)
        )
    with pytest.raises(ValueError, match="canonical shape"):
        response_residual_report(
            np.zeros((2, 2), dtype=complex),
            np.zeros((2, 2), dtype=complex),
        )
    with pytest.raises(ValueError, match="not be empty"):
        response_residual_report(
            np.zeros((0, 2, 1), dtype=complex),
            np.zeros((0, 2, 1), dtype=complex),
        )


def test_report_labels_are_attached_and_validated():
    report = response_residual_report(
        _PRED,
        _OBS,
        station_names=["S1", "S2"],
        frequencies_hz=[10.0, 1.0],
        components=["zxy"],
    )
    assert report.station_names == ("S1", "S2")
    np.testing.assert_allclose(report.frequencies_hz, [10.0, 1.0])
    assert report.components == ("zxy",)


def test_report_validates_shapes_and_types():
    with pytest.raises(TypeError, match="ResponseLossResult"):
        ResponseResidualReport(
            overall="not a result",
            by_station=np.zeros(2),
            by_frequency=np.zeros(2),
            by_component=np.zeros(1),
            station_names=None,
            frequencies_hz=None,
            components=None,
            shape=(2, 2, 1),
        )
    with pytest.raises(ValueError, match="by_station must have shape"):
        ResponseResidualReport(
            overall=_OVERALL,
            by_station=np.zeros(3),
            by_frequency=np.zeros(2),
            by_component=np.zeros(1),
            station_names=None,
            frequencies_hz=None,
            components=None,
            shape=(2, 2, 1),
        )
    with pytest.raises(
        ValueError, match="station_names must have length"
    ):
        ResponseResidualReport(
            overall=_OVERALL,
            by_station=np.zeros(2),
            by_frequency=np.zeros(2),
            by_component=np.zeros(1),
            station_names=("only-one",),
            frequencies_hz=None,
            components=None,
            shape=(2, 2, 1),
        )
    with pytest.raises(
        ValueError, match="frequencies_hz must have shape"
    ):
        ResponseResidualReport(
            overall=_OVERALL,
            by_station=np.zeros(2),
            by_frequency=np.zeros(2),
            by_component=np.zeros(1),
            station_names=None,
            frequencies_hz=np.zeros(3),
            components=None,
            shape=(2, 2, 1),
        )
    with pytest.raises(ValueError, match="components must have length"):
        ResponseResidualReport(
            overall=_OVERALL,
            by_station=np.zeros(2),
            by_frequency=np.zeros(2),
            by_component=np.zeros(1),
            station_names=None,
            frequencies_hz=None,
            components=("a", "b"),
            shape=(2, 2, 1),
        )


def _forward_result(impedance):
    diagnostics = SolverDiagnostics(
        [[True], [True]], [[1], [1]], [[0.0], [0.0]], 0.01
    )
    return ForwardResult(
        "a" * 64,
        [10.0, 1.0],
        ["S1", "S2"],
        ["zxy"],
        impedance,
        None,
        "demo",
        "1",
        diagnostics,
    )


def test_report_from_contracts_matches_raw_array_report_and_labels():
    observed = SurveyData(
        _OBS, [10.0, 1.0], ["S1", "S2"], ["zxy"], [[0.0, 0.0], [1.0, 0.0]]
    )
    forward = _forward_result(_PRED)
    report = response_residual_report_from_contracts(forward, observed)
    raw = response_residual_report(_PRED, _OBS)
    np.testing.assert_allclose(report.by_station, raw.by_station)
    np.testing.assert_allclose(report.by_frequency, raw.by_frequency)
    assert report.station_names == ("S1", "S2")
    np.testing.assert_allclose(report.frequencies_hz, [10.0, 1.0])
    assert report.components == ("zxy",)
