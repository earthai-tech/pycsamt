"""Deterministic equivalence contracts for native Occam1D restarts."""

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DConfig,
    Occam1DData,
    Occam1DForwardModel,
    Occam1DInversion,
    Occam1DModel,
    Occam1DRestart,
    Occam1DSolverPolicy,
)


def _problem(max_iterations):
    depth = np.array([0.0, 30.0, 100.0, 300.0, 800.0])
    truth = Occam1DModel(
        depth,
        [60.0, 250.0, 18.0, 500.0, 90.0],
    )
    frequency = np.logspace(-1.5, 2.5, 14)
    response = Occam1DForwardModel(truth).predict(frequency)
    data = Occam1DData(
        frequency,
        response.resistivity,
        response.phase,
        np.full(frequency.size, 0.025),
        np.full(frequency.size, 1.2),
        station="RESTART_EQ",
        mode="determinant",
    )
    model = Occam1DModel(depth, np.full(depth.size, 100.0))
    config = Occam1DConfig(
        n_layers=depth.size,
        first_thickness=30.0,
        depth_max=800.0,
        starting_resistivity=100.0,
        target_misfit=0.01,
        max_iterations=max_iterations,
        stepsize_cut_count=5,
        lagrange_start=3.0,
        lagrange_min=1e-8,
        lagrange_max=1e8,
    )
    return data, model, config


def _solver(data, model, config, **kwargs):
    return Occam1DInversion(
        data,
        model,
        config=config,
        multiplier_factors=[0.01, 0.1, 1.0, 10.0, 100.0],
        **kwargs,
    )


def _assert_iteration_equal(resumed, uninterrupted):
    assert resumed.number == uninterrupted.number
    np.testing.assert_allclose(
        resumed.parameters,
        uninterrupted.parameters,
        rtol=1e-12,
        atol=1e-13,
    )
    np.testing.assert_allclose(
        resumed.prediction,
        uninterrupted.prediction,
        rtol=1e-12,
        atol=1e-13,
    )
    np.testing.assert_allclose(
        resumed.residual,
        uninterrupted.residual,
        rtol=1e-12,
        atol=1e-13,
    )
    assert resumed.rms == pytest.approx(uninterrupted.rms, rel=1e-12)
    assert resumed.roughness == pytest.approx(
        uninterrupted.roughness,
        rel=1e-12,
    )
    assert resumed.multiplier == uninterrupted.multiplier
    assert resumed.step_scale == uninterrupted.step_scale
    assert resumed.rank == uninterrupted.rank
    assert resumed.condition_number == pytest.approx(
        uninterrupted.condition_number,
        rel=1e-12,
    )
    assert np.array_equal(resumed.active_lower, uninterrupted.active_lower)
    assert np.array_equal(resumed.active_upper, uninterrupted.active_upper)
    assert resumed.solver == uninterrupted.solver
    assert resumed.solve_attempts == uninterrupted.solve_attempts
    assert resumed.damping == uninterrupted.damping
    if uninterrupted.initial_condition_number is None:
        assert resumed.initial_condition_number is None
    else:
        assert resumed.initial_condition_number == pytest.approx(
            uninterrupted.initial_condition_number,
            rel=1e-12,
        )


def _assert_result_equal(resumed, uninterrupted):
    assert resumed.convergence is uninterrupted.convergence
    assert resumed.target_rms == uninterrupted.target_rms
    assert resumed.message == uninterrupted.message
    assert resumed.n_iterations == uninterrupted.n_iterations
    assert len(resumed.iterations) == len(uninterrupted.iterations)
    for left, right in zip(resumed.iterations, uninterrupted.iterations):
        _assert_iteration_equal(left, right)
    assert [item.to_dict() for item in resumed.rejected_candidates] == [
        item.to_dict() for item in uninterrupted.rejected_candidates
    ]
    assert [item.to_dict() for item in resumed.failed_iterations] == [
        item.to_dict() for item in uninterrupted.failed_iterations
    ]


@pytest.mark.parametrize("split_iteration", [1, 2, 3])
def test_json_restart_matches_uninterrupted_at_multiple_splits(
    tmp_path,
    split_iteration,
):
    data, model, full_config = _problem(max_iterations=4)
    uninterrupted = _solver(data, model, full_config).run()

    short_config = full_config.updated(max_iterations=split_iteration)
    partial_solver = _solver(data, model, short_config)
    partial_solver.run()
    path = partial_solver.restart().write(
        tmp_path / f"split-{split_iteration}.json"
    )
    checkpoint = Occam1DRestart.read(path)
    resumed = _solver(data, model, full_config).run(restart=checkpoint)

    _assert_result_equal(resumed, uninterrupted)


def test_bounded_restart_matches_uninterrupted_execution(tmp_path):
    data, model, full_config = _problem(max_iterations=4)
    bounds = (
        np.array([1.5, 1.7, 1.0, 1.8, 1.5]),
        np.array([2.5, 2.6, 2.3, 2.8, 2.5]),
    )
    uninterrupted = _solver(
        data,
        model,
        full_config,
        log_resistivity_bounds=bounds,
    ).run()
    partial_solver = _solver(
        data,
        model,
        full_config.updated(max_iterations=2),
        log_resistivity_bounds=bounds,
    )
    partial_solver.run()
    checkpoint = Occam1DRestart.read(
        partial_solver.restart().write(tmp_path / "bounded.json")
    )
    resumed = _solver(
        data,
        model,
        full_config,
        log_resistivity_bounds=bounds,
    ).run(restart=checkpoint)

    _assert_result_equal(resumed, uninterrupted)


def test_fallback_solver_restart_preserves_numerical_path(
    monkeypatch,
    tmp_path,
):
    data, model, full_config = _problem(max_iterations=3)
    policy = Occam1DSolverPolicy(
        condition_limit=1e14,
        ill_condition_action="accept",
        svd_failure_action="fallback",
    )

    def failed_numpy_svd(*args, **kwargs):
        raise np.linalg.LinAlgError("forced restart-equivalence fallback")

    monkeypatch.setattr(np.linalg, "lstsq", failed_numpy_svd)
    uninterrupted = _solver(
        data,
        model,
        full_config,
        solver_policy=policy,
    ).run()
    partial_solver = _solver(
        data,
        model,
        full_config.updated(max_iterations=1),
        solver_policy=policy,
    )
    partial_solver.run()
    checkpoint = Occam1DRestart.read(
        partial_solver.restart().write(tmp_path / "fallback.json")
    )
    resumed = _solver(
        data,
        model,
        full_config,
        solver_policy=policy,
    ).run(restart=checkpoint)

    _assert_result_equal(resumed, uninterrupted)
    assert all(
        item.solver in {"initial", "scipy_gelsy_fallback"}
        for item in resumed.iterations
    )
