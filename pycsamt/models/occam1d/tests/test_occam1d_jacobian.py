"""Scientific contracts for analytic Occam1D response sensitivities."""

import json

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DData,
    Occam1DForwardModel,
    Occam1DJacobian,
    Occam1DJacobianMethod,
    Occam1DModel,
)


def _model():
    return Occam1DModel(
        [0.0, 10.0, 40.0, 150.0, 600.0],
        [30.0, 300.0, 15.0, 800.0, 100.0],
    )


def _data(mode="determinant"):
    frequency = np.logspace(3, -3, 13)
    return Occam1DData(
        frequency,
        np.full(frequency.size, 100.0),
        np.full(frequency.size, 45.0),
        np.full(frequency.size, 0.05),
        np.full(frequency.size, 2.0),
        mode=mode,
        station="S01",
    )


def test_analytic_jacobian_matches_central_reference():
    check = Occam1DJacobian(_model()).check(_data())
    assert check.passed
    assert check.max_absolute_error < 2e-6
    assert check.max_relative_error < 5e-4


@pytest.mark.parametrize("mode", ["xy", "yx", "te", "tm", "determinant"])
def test_analytic_derivative_supports_all_modes(mode):
    check = Occam1DJacobian(_model()).check(_data(mode))
    assert check.passed


def test_uniform_model_common_shift_has_exact_scale_invariant():
    model = Occam1DModel(
        [0.0, 10.0, 100.0, 1000.0],
        [100.0] * 4,
    )
    result = Occam1DJacobian(model).compute(_data())
    rho_rows = result.type_code == 103
    phase_rows = result.type_code == 104
    np.testing.assert_allclose(
        np.sum(result.matrix[rho_rows], axis=1),
        1.0,
        atol=3e-14,
    )
    np.testing.assert_allclose(
        np.sum(result.matrix[phase_rows], axis=1),
        0.0,
        atol=3e-13,
    )


def test_rows_follow_native_missing_observation_order():
    data = Occam1DData(
        [100.0, 10.0, 1.0],
        [100.0, np.nan, 100.0],
        [45.0, 45.0, np.nan],
        [0.05, np.nan, 0.05],
        [2.0, 2.0, np.nan],
        mode="xy",
    )
    result = Occam1DJacobian(_model()).compute(data)
    np.testing.assert_array_equal(result.frequency_index, [1, 1, 2, 3])
    np.testing.assert_array_equal(result.type_code, [103, 104, 104, 103])
    assert result.matrix.shape == (data.n_data, _model().n_layers)


def test_weighting_scales_rows_and_predictions_by_solver_error():
    data = _data()
    jacobian = Occam1DJacobian(_model())
    raw = jacobian.compute(data)
    weighted = jacobian.compute(data, weighted=True)
    errors = jacobian.solver_errors(data)
    np.testing.assert_allclose(
        weighted.matrix,
        raw.matrix / errors[:, None],
    )
    np.testing.assert_allclose(
        weighted.prediction,
        raw.prediction / errors,
    )


def test_central_result_records_per_layer_steps():
    steps = np.linspace(5e-5, 2e-4, _model().n_layers)
    result = Occam1DJacobian(_model()).compute(
        _data(),
        method=Occam1DJacobianMethod.CENTRAL,
        difference_step=steps,
    )
    np.testing.assert_allclose(result.step, steps)
    assert result.method == "central"


def test_central_n_jobs_matches_sequential_columns():
    pytest.importorskip("joblib")
    data = _data()
    jacobian = Occam1DJacobian(_model())
    sequential = jacobian.compute(data, method="central", n_jobs=1)
    parallel = jacobian.compute(data, method="central", n_jobs=2)
    np.testing.assert_array_equal(sequential.matrix, parallel.matrix)
    np.testing.assert_array_equal(sequential.prediction, parallel.prediction)


def test_result_is_immutable_and_exposes_linear_algebra_diagnostics():
    result = Occam1DJacobian(_model()).compute(_data())
    with pytest.raises(ValueError):
        result.matrix[0, 0] = 0.0
    assert 0 < result.rank <= min(result.shape)
    assert result.singular_values.size == min(result.shape)
    assert result.column_norm.shape == (_model().n_layers,)
    values = json.loads(json.dumps(result.diagnostics()))
    assert values["shape"] == list(result.shape)
    assert values["method"] == "analytic"


def test_supplied_forward_permeability_is_preserved_without_mutation():
    model = _model()
    forward = Occam1DForwardModel(model, permeability=1.1e-6)
    jacobian = Occam1DJacobian(model, forward=forward)
    assert jacobian.permeability == pytest.approx(1.1e-6)
    assert forward.is_configured


def test_parameter_and_method_validation_is_explicit():
    jacobian = Occam1DJacobian(_model())
    with pytest.raises(ValueError, match="shape"):
        jacobian.compute(_data(), parameters=[2.0])
    with pytest.raises(ValueError, match="method"):
        jacobian.compute(_data(), method="complex-step")
    with pytest.raises(TypeError, match="weighted"):
        jacobian.compute(_data(), weighted=1)


def test_difference_step_validation_is_explicit():
    jacobian = Occam1DJacobian(_model())
    with pytest.raises(ValueError, match="positive"):
        jacobian.compute(
            _data(),
            method="central",
            difference_step=0.0,
        )
    with pytest.raises(ValueError, match="shape"):
        jacobian.compute(
            _data(),
            method="central",
            difference_step=[1e-4],
        )


def test_check_diagnostics_are_json_safe():
    check = Occam1DJacobian(_model()).check(_data())
    values = json.loads(json.dumps(check.diagnostics()))
    assert values["passed"] is True
    assert values["shape"] == [26, 5]
