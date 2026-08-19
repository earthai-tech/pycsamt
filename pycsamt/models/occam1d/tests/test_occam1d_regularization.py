"""Operator, objective, and stable-solve tests for Occam1D regularization."""

import json

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DCandidate,
    Occam1DData,
    Occam1DJacobian,
    Occam1DLinearSolveError,
    Occam1DModel,
    Occam1DRegularization,
    Occam1DSolverPolicy,
)


def _model():
    return Occam1DModel(
        [0.0, 10.0, 30.0, 100.0, 400.0],
        [30.0, 200.0, 15.0, 500.0, 100.0],
    )


def _data():
    frequency = np.logspace(2, -2, 9)
    return Occam1DData(
        frequency,
        np.full(frequency.size, 100.0),
        np.full(frequency.size, 45.0),
        np.full(frequency.size, 0.05),
        np.full(frequency.size, 2.0),
    )


def test_first_difference_operator_matches_paper_definition():
    regularization = Occam1DRegularization(4, order=1)
    np.testing.assert_array_equal(
        regularization.operator,
        [
            [-1.0, 1.0, 0.0, 0.0],
            [0.0, -1.0, 1.0, 0.0],
            [0.0, 0.0, -1.0, 1.0],
        ],
    )
    assert regularization.evaluate([2.0] * 4).roughness == 0.0


def test_second_difference_operator_matches_paper_definition():
    regularization = Occam1DRegularization(5, order=2)
    np.testing.assert_array_equal(
        regularization.operator,
        [
            [1.0, -2.0, 1.0, 0.0, 0.0],
            [0.0, 1.0, -2.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, -2.0, 1.0],
        ],
    )
    assert regularization.evaluate([1, 2, 3, 4, 5]).roughness == 0.0


def test_depth_scaled_first_derivative_approximates_integral():
    depth = [0.0, 1.0, 5.0]
    regularization = Occam1DRegularization(
        3,
        order=1,
        depth=depth,
        scaling="depth",
    )
    values = np.asarray(depth)
    result = regularization.evaluate(values)
    assert result.roughness == pytest.approx(5.0)


def test_depth_scaled_second_derivative_annihilates_linear_depth_model():
    depth = np.array([0.0, 1.0, 4.0, 10.0])
    regularization = Occam1DRegularization(
        4,
        order=2,
        depth=depth,
        scaling="depth",
    )
    result = regularization.evaluate(2.0 + 0.3 * depth)
    np.testing.assert_allclose(result.smoothness, 0.0, atol=1e-15)


def test_row_weights_are_residual_amplitudes():
    regularization = Occam1DRegularization(
        3,
        row_weights=[2.0, 3.0],
    )
    result = regularization.evaluate([0.0, 1.0, 3.0])
    np.testing.assert_allclose(result.smoothness, [2.0, 6.0])
    assert result.roughness == pytest.approx(40.0)


def test_reference_penalty_is_separate_from_roughness():
    regularization = Occam1DRegularization(
        3,
        reference=[2.0, 2.0, 2.0],
        reference_weights=[1.0, 0.0, 2.0],
    )
    result = regularization.evaluate([3.0, 3.0, 3.0])
    assert result.roughness == pytest.approx(0.0)
    assert result.reference_penalty == pytest.approx(5.0)
    assert result.objective == pytest.approx(5.0)


def test_augmented_solve_matches_direct_ridge_solution():
    matrix = np.eye(3)
    rhs = np.array([1.0, 3.0, 5.0])
    regularization = Occam1DRegularization(3)
    multiplier = 4.0
    candidate = regularization.solve(matrix, rhs, multiplier)
    operator = regularization.operator
    expected = np.linalg.solve(
        matrix.T @ matrix + multiplier * operator.T @ operator,
        matrix.T @ rhs,
    )
    np.testing.assert_allclose(candidate.parameters, expected)
    assert candidate.objective == pytest.approx(
        candidate.data_objective
        + multiplier * candidate.structure_objective
    )


def test_larger_multiplier_produces_no_rougher_candidate():
    matrix = np.eye(5)
    rhs = np.array([1.0, 4.0, -2.0, 5.0, 0.0])
    regularization = Occam1DRegularization(5)
    weak = regularization.solve(matrix, rhs, 0.01)
    strong = regularization.solve(matrix, rhs, 100.0)
    assert strong.structure_objective <= weak.structure_objective
    assert strong.data_objective >= weak.data_objective


def test_rank_deficient_data_system_is_solved_without_normal_inverse():
    matrix = np.ones((2, 4))
    candidate = Occam1DRegularization(4).solve(matrix, [4.0, 4.0], 1.0)
    assert isinstance(candidate, Occam1DCandidate)
    assert candidate.rank > 0
    assert np.all(np.isfinite(candidate.parameters))
    np.testing.assert_allclose(candidate.parameters, 1.0)


def test_linearized_system_respects_jacobian_weighting():
    model = _model()
    data = _data()
    current = np.log10(model.resistivity)
    jacobian_engine = Occam1DJacobian(model)
    regularization = Occam1DRegularization.from_model(model)
    raw = jacobian_engine.compute(data)
    weighted = jacobian_engine.compute(data, weighted=True)
    raw_matrix, raw_rhs = regularization.linearized_system(
        data,
        raw,
        current,
    )
    weighted_matrix, weighted_rhs = regularization.linearized_system(
        data,
        weighted,
        current,
    )
    errors = jacobian_engine.solver_errors(data)
    np.testing.assert_allclose(weighted_matrix, raw_matrix / errors[:, None])
    np.testing.assert_allclose(weighted_rhs, raw_rhs / errors)


def test_native_preferences_are_opt_in():
    model = Occam1DModel(
        [0.0, 10.0, 30.0],
        [100.0] * 3,
        preference=[1.0, 2.0, 3.0],
        weight=[0.0, 2.0, 0.0],
    )
    default = Occam1DRegularization.from_model(model)
    native = Occam1DRegularization.from_model(
        model,
        use_native_preference=True,
    )
    assert default.n_reference == 0
    assert native.n_reference == 1
    assert native.evaluate([1.0, 3.0, 3.0]).reference_penalty == 4.0


def test_results_are_immutable_and_diagnostics_are_json_safe():
    regularization = Occam1DRegularization(3)
    result = regularization.evaluate([1.0, 2.0, 4.0])
    candidate = regularization.solve(np.eye(3), [1, 2, 3], 1.0)
    with pytest.raises(ValueError):
        result.parameters[0] = 0.0
    with pytest.raises(ValueError):
        candidate.parameters[0] = 0.0
    assert json.loads(json.dumps(result.diagnostics()))["roughness"] == 5.0
    values = json.loads(json.dumps(candidate.diagnostics()))
    assert values["solver"] == "svd_lstsq"


def test_configuration_and_system_validation_is_explicit():
    with pytest.raises(ValueError, match="order"):
        Occam1DRegularization(4, order=3)
    with pytest.raises(ValueError, match="depth is required"):
        Occam1DRegularization(4, scaling="depth")
    with pytest.raises(ValueError, match="reference is required"):
        Occam1DRegularization(3, reference_weights=[1, 0, 0])
    regularization = Occam1DRegularization(3)
    with pytest.raises(ValueError, match="matrix"):
        regularization.solve(np.eye(2), [1, 2], 1.0)
    with pytest.raises(ValueError, match="multiplier"):
        regularization.solve(np.eye(3), [1, 2, 3], -1.0)


def test_bounded_solver_minimizes_inside_log_resistivity_limits():
    regularization = Occam1DRegularization(3, order=1)
    candidate = regularization.solve(
        np.eye(3),
        np.array([-4.0, 2.0, 9.0]),
        0.0,
        bounds=(0.0, 4.0),
    )
    assert np.all(candidate.parameters >= 0.0)
    assert np.all(candidate.parameters <= 4.0)
    assert candidate.parameters == pytest.approx([0.0, 2.0, 4.0])
    assert candidate.solver == "scipy_lsq_linear_trf"
    assert candidate.n_active_bounds == 2
    assert candidate.active_lower.tolist() == [True, False, False]
    assert candidate.active_upper.tolist() == [False, False, True]


def test_bounds_accept_layer_vectors_and_reject_invalid_intervals():
    regularization = Occam1DRegularization(3)
    candidate = regularization.solve(
        np.eye(3),
        np.array([0.0, 3.0, 8.0]),
        0.0,
        bounds=([0.0, 1.0, 2.0], [1.0, 4.0, 5.0]),
    )
    assert candidate.parameters == pytest.approx(
        [0.0, 3.0, 5.0],
        abs=2e-5,
    )
    with pytest.raises(ValueError, match="below its upper"):
        regularization.solve(
            np.eye(3),
            np.ones(3),
            0.0,
            bounds=(2.0, 2.0),
        )
    with pytest.raises(ValueError, match="shape"):
        regularization.solve(
            np.eye(3),
            np.ones(3),
            0.0,
            bounds=([0.0, 1.0], 4.0),
        )


def test_numpy_svd_failure_uses_explicit_pivoted_qr_fallback(monkeypatch):
    regularization = Occam1DRegularization(3)

    def failed_svd(*args, **kwargs):
        raise np.linalg.LinAlgError("synthetic SVD nonconvergence")

    monkeypatch.setattr(np.linalg, "lstsq", failed_svd)
    candidate = regularization.solve(
        np.eye(3),
        np.array([1.0, 2.0, 3.0]),
        1.0,
        policy=Occam1DSolverPolicy(svd_failure_action="fallback"),
    )
    assert candidate.solver == "scipy_gelsy_fallback"
    assert np.all(np.isfinite(candidate.parameters))


def test_svd_raise_policy_exposes_exhausted_solver_error(monkeypatch):
    regularization = Occam1DRegularization(3)

    def failed_svd(*args, **kwargs):
        raise np.linalg.LinAlgError("synthetic SVD nonconvergence")

    monkeypatch.setattr(np.linalg, "lstsq", failed_svd)
    with pytest.raises(Occam1DLinearSolveError, match="SVD"):
        regularization.solve(
            np.eye(3),
            np.ones(3),
            1.0,
            policy=Occam1DSolverPolicy(svd_failure_action="raise"),
        )


def test_ill_conditioned_system_is_retried_with_diagonal_damping():
    regularization = Occam1DRegularization(2)
    candidate = regularization.solve(
        np.diag([1.0, 1e-10]),
        np.ones(2),
        0.0,
        policy=Occam1DSolverPolicy(
            condition_limit=1e4,
            ill_condition_action="damp",
            max_retries=2,
            damping_start=1e-4,
            damping_growth=10.0,
        ),
    )
    assert candidate.initial_condition_number > 1e4
    assert candidate.condition_number <= 1e4
    assert candidate.solve_attempts == 2
    assert candidate.damping == pytest.approx(1e-4)


def test_ill_condition_reject_policy_is_explicit():
    regularization = Occam1DRegularization(2)
    with pytest.raises(Occam1DLinearSolveError, match="condition number"):
        regularization.solve(
            np.diag([1.0, 1e-14]),
            np.ones(2),
            0.0,
            policy=Occam1DSolverPolicy(
                condition_limit=1e6,
                ill_condition_action="reject",
            ),
        )
