"""Scientific and lifecycle tests for native Occam1D inversion."""

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DConfig,
    Occam1DConvergence,
    Occam1DData,
    Occam1DForwardModel,
    Occam1DInversion,
    Occam1DInversionResult,
    Occam1DIteration,
    Occam1DModel,
    Occam1DRestart,
    Occam1DRejectionReason,
    Occam1DSolverPolicy,
    Occam1DStartup,
)


def _problem(*, target=1.0, iterations=6, noise=0.02):
    depth = np.array([0.0, 40.0, 150.0, 500.0])
    truth = Occam1DModel(depth, [80.0, 20.0, 400.0, 100.0])
    frequency = np.logspace(-1.0, 2.0, 10)
    response = Occam1DForwardModel(truth).predict(frequency)
    data = Occam1DData(
        frequency,
        response.resistivity,
        response.phase,
        np.full(frequency.size, noise),
        np.full(frequency.size, 1.0),
        station="SYN01",
        mode="determinant",
    )
    start = Occam1DModel(depth, np.full(depth.size, 100.0))
    config = Occam1DConfig(
        n_layers=4,
        first_thickness=40.0,
        depth_max=500.0,
        starting_resistivity=100.0,
        target_misfit=target,
        max_iterations=iterations,
        stepsize_cut_count=6,
    )
    return data, start, config


def test_inversion_reduces_true_nonlinear_rms():
    data, model, config = _problem(target=0.8)
    result = Occam1DInversion(
        data,
        model,
        config=config,
        multiplier_factors=[0.01, 0.1, 1.0, 10.0],
    ).run()
    assert result.final.rms < result.initial.rms
    assert result.n_iterations >= 1
    assert np.all(np.isfinite(result.final.resistivity))


def test_native_result_export_and_main_images(tmp_path):
    data, model, config = _problem(target=0.8, iterations=1)
    inversion = Occam1DInversion(data, model, config=config)
    result = inversion.run()

    text_paths = inversion.export_result(tmp_path / "model-text", result)
    image_paths = inversion.save_main_images(
        tmp_path / "model-image",
        result,
        dpi=72,
        style_overrides={
            "observed__marker": "x",
            "model__color": "green",
            "target__linestyle": ":",
            "model_legend": False,
        },
    )

    assert set(text_paths) == {
        "model",
        "response",
        "iterations",
        "rejected",
        "failures",
        "summary",
        "metadata",
    }
    assert set(image_paths) == {
        "model", "response", "convergence", "summary"
    }
    assert all(
        path.is_file() and path.stat().st_size
        for path in text_paths.values()
    )
    assert all(
        path.is_file() and path.stat().st_size
        for path in image_paths.values()
    )
    assert "pycsamt.occam1d.native-result/v1" in text_paths[
        "metadata"
    ].read_text(encoding="utf8")


def test_exact_starting_model_stops_at_iteration_zero():
    data, _, config = _problem(target=1.0)
    model = Occam1DModel(
        [0.0, 40.0, 150.0, 500.0],
        [80.0, 20.0, 400.0, 100.0],
    )
    startup = Occam1DStartup.from_model(model, config)
    result = Occam1DInversion(
        data,
        model,
        config=config,
        startup=startup,
    ).run()
    assert result.convergence is Occam1DConvergence.TARGET
    assert result.n_iterations == 0


def test_callback_receives_initial_and_accepted_models():
    data, model, config = _problem(target=0.01, iterations=2)
    seen = []
    result = Occam1DInversion(
        data,
        model,
        config=config,
        multiplier_factors=[0.1, 1.0],
    ).run(callback=lambda item: seen.append(item.number))
    assert seen == list(range(result.n_iterations + 1))


def test_cooperative_cancellation_preserves_initial_model():
    data, model, config = _problem(target=0.01)
    result = Occam1DInversion(data, model, config=config).run(
        cancel=lambda: True
    )
    assert result.convergence is Occam1DConvergence.CANCELLED
    assert result.n_iterations == 0


def test_history_and_arrays_are_immutable():
    item = Occam1DIteration(
        0,
        np.array([2.0, 2.0]),
        np.array([2.0]),
        np.array([0.0]),
        0.0,
        0.0,
        1.0,
        0.0,
        0,
        float("inf"),
    )
    result = Occam1DInversionResult(
        (item,),
        Occam1DConvergence.TARGET,
        1.0,
        "done",
    )
    with pytest.raises(ValueError):
        item.parameters[0] = 3.0
    assert result.final is item
    assert result.diagnostics()["converged"] is True


def test_phase_residual_uses_nearest_180_degree_branch(monkeypatch):
    data, model, config = _problem()
    inversion = Occam1DInversion(data, model, config=config)
    observed = inversion.regularization.observed_solver_vector(data)
    phase_rows = np.tile([False, True], data.n_frequencies)
    prediction = observed.copy()
    prediction[phase_rows] += 179.0
    monkeypatch.setattr(inversion.forward, "solver_vector", lambda *a, **k: (
        prediction.copy()
    ))
    item = inversion._evaluate(
        inversion.startup.parameters,
        0,
        1.0,
        0.0,
        0,
        float("inf"),
    )
    assert np.allclose(item.residual[phase_rows], -1.0)


def test_constructor_rejects_layer_mismatch_and_bad_factors():
    data, model, config = _problem()
    bad = config.updated(n_layers=5, depth_max=540.0)
    with pytest.raises(ValueError, match="layer counts"):
        Occam1DInversion(data, model, config=bad)
    with pytest.raises(ValueError, match="strictly positive"):
        Occam1DInversion(
            data,
            model,
            config=config,
            multiplier_factors=[0.0, 1.0],
        )


def test_invalid_callbacks_are_rejected_at_public_boundary():
    data, model, config = _problem()
    inversion = Occam1DInversion(data, model, config=config)
    with pytest.raises(TypeError, match="callback"):
        inversion.run(callback=1)
    with pytest.raises(TypeError, match="cancel"):
        inversion.run(cancel=1)


def test_unset_model_resistivity_uses_startup_fill():
    data, model, config = _problem()
    unset = Occam1DModel(model.depth, np.full(model.n_layers, np.nan))
    inversion = Occam1DInversion(data, unset, config=config)
    assert np.allclose(
        inversion.startup.physical_resistivity,
        config.starting_resistivity,
    )
    assert np.isfinite(inversion._evaluate(
        inversion.startup.parameters,
        0,
        1.0,
        0.0,
        0,
        float("inf"),
    ).rms)


def test_inversion_enforces_layerwise_log_resistivity_bounds():
    data, model, config = _problem(target=0.01, iterations=3)
    lower = np.array([1.8, 1.7, 1.8, 1.9])
    upper = np.array([2.1, 2.1, 2.2, 2.1])
    inversion = Occam1DInversion(
        data,
        model,
        config=config,
        log_resistivity_bounds=(lower, upper),
        multiplier_factors=[0.01, 0.1, 1.0],
    )
    result = inversion.run()
    for iteration in result.iterations:
        assert np.all(iteration.parameters >= lower)
        assert np.all(iteration.parameters <= upper)
    returned_lower, returned_upper = inversion.log_resistivity_bounds
    returned_lower[0] = -99.0
    assert inversion.log_resistivity_bounds[0][0] == lower[0]
    assert np.array_equal(returned_upper, upper)


def test_inversion_rejects_infeasible_starting_model():
    data, model, config = _problem()
    with pytest.raises(ValueError, match="violate bounds"):
        Occam1DInversion(
            data,
            model,
            config=config,
            log_resistivity_bounds=(3.0, 4.0),
        )


def test_absolute_lagrange_bounds_saturate_and_deduplicate_search():
    data, model, config = _problem()
    inversion = Occam1DInversion(
        data,
        model,
        config=config,
        multiplier_factors=[1e-9, 0.1, 1.0, 10.0, 1e9],
        lagrange_bounds=(1.0, 10.0),
    )
    assert inversion.lagrange_bounds == (1.0, 10.0)
    assert inversion.candidate_multipliers() == pytest.approx(
        [1.0, 5.0, 10.0]
    )
    diagnostics = inversion.diagnostics()
    assert diagnostics["lagrange_min"] == 1.0
    assert diagnostics["lagrange_max"] == 10.0


def test_lagrange_bounds_reject_start_outside_interval():
    data, model, config = _problem()
    with pytest.raises(ValueError, match="lagrange_start"):
        Occam1DInversion(
            data,
            model,
            config=config,
            lagrange_bounds=(0.1, 2.0),
        )


def test_all_accepted_multipliers_remain_inside_absolute_bounds():
    data, model, config = _problem(target=0.01, iterations=3)
    inversion = Occam1DInversion(
        data,
        model,
        config=config,
        multiplier_factors=[1e-12, 0.1, 10.0, 1e12],
        lagrange_bounds=(0.5, 20.0),
    )
    result = inversion.run()
    for iteration in result.iterations[1:]:
        assert 0.5 <= iteration.multiplier <= 20.0


def test_restart_continuation_matches_uninterrupted_inversion():
    data, model, short_config = _problem(target=0.01, iterations=1)
    factors = [0.01, 0.1, 1.0, 10.0]
    short_solver = Occam1DInversion(
        data,
        model,
        config=short_config,
        multiplier_factors=factors,
    )
    partial = short_solver.run()
    checkpoint = short_solver.restart()

    long_config = short_config.updated(max_iterations=3)
    resumed = Occam1DInversion(
        data,
        model,
        config=long_config,
        multiplier_factors=factors,
    ).run(restart=checkpoint)
    uninterrupted = Occam1DInversion(
        data,
        model,
        config=long_config,
        multiplier_factors=factors,
    ).run()

    assert checkpoint.iteration_number == partial.final.number
    assert checkpoint.multiplier == partial.final.multiplier
    assert np.array_equal(checkpoint.parameters, partial.final.parameters)
    assert resumed.n_iterations == uninterrupted.n_iterations
    assert resumed.final.parameters == pytest.approx(
        uninterrupted.final.parameters
    )
    assert resumed.rms_history == pytest.approx(
        uninterrupted.rms_history
    )


def test_restart_json_roundtrip_preserves_complete_history(tmp_path):
    data, model, config = _problem(target=0.01, iterations=2)
    solver = Occam1DInversion(data, model, config=config)
    solver.run()
    checkpoint = solver.restart()
    path = checkpoint.write(tmp_path / "restart.json")
    restored = Occam1DRestart.read(path)

    assert restored.iteration_number == checkpoint.iteration_number
    assert restored.multiplier == checkpoint.multiplier
    assert restored.parameters == pytest.approx(checkpoint.parameters)
    assert restored.rms_history == pytest.approx(checkpoint.rms_history)
    assert restored.roughness_history == pytest.approx(
        checkpoint.roughness_history
    )
    assert restored.data_fingerprint == checkpoint.data_fingerprint
    assert restored.previous_convergence == checkpoint.previous_convergence
    assert restored.previous_message == checkpoint.previous_message
    assert restored.rejected_candidates == checkpoint.rejected_candidates
    assert restored.failed_iterations == checkpoint.failed_iterations
    assert restored.diagnostics()["rms_history"] == pytest.approx(
        checkpoint.rms_history
    )


def test_restart_rejects_different_sounding_and_geometry():
    data, model, config = _problem(target=0.01, iterations=1)
    solver = Occam1DInversion(data, model, config=config)
    solver.run()
    checkpoint = solver.restart()

    changed_data = Occam1DData(
        data.frequency,
        data.resistivity * 1.01,
        data.phase,
        data.resistivity_error,
        data.phase_error,
        station=data.station,
        mode=data.mode,
    )
    with pytest.raises(ValueError, match="fingerprint"):
        Occam1DInversion(
            changed_data,
            model,
            config=config,
        ).run(restart=checkpoint)

    changed_model = Occam1DModel(
        [0.0, 45.0, 150.0, 500.0],
        model.resistivity,
    )
    with pytest.raises(ValueError, match="geometry"):
        Occam1DInversion(
            data,
            changed_model,
            config=config,
        ).run(restart=checkpoint)


def test_restart_does_not_replay_callbacks_for_old_history():
    data, model, short_config = _problem(target=0.01, iterations=1)
    first = Occam1DInversion(data, model, config=short_config)
    first.run()
    checkpoint = first.restart()
    seen = []
    Occam1DInversion(
        data,
        model,
        config=short_config.updated(max_iterations=2),
    ).run(
        restart=checkpoint,
        callback=lambda item: seen.append(item.number),
    )
    assert all(number > checkpoint.iteration_number for number in seen)


def test_nonselected_candidates_are_structured_and_queryable():
    data, model, config = _problem(target=0.01, iterations=1)
    result = Occam1DInversion(
        data,
        model,
        config=config,
        multiplier_factors=[0.01, 0.1, 1.0, 10.0],
    ).run()
    assert result.rejected_candidates
    assert all(
        record.iteration == 1 for record in result.rejected_candidates
    )
    assert any(
        record.reason is Occam1DRejectionReason.NOT_SELECTED
        for record in result.rejected_candidates
    )
    diagnostics = result.diagnostics()
    assert diagnostics["n_rejected_candidates"] == len(
        result.rejected_candidates
    )


def test_one_failed_linear_candidate_does_not_abort_search(monkeypatch):
    data, model, config = _problem(target=0.01, iterations=1)
    solver = Occam1DInversion(
        data,
        model,
        config=config,
        multiplier_factors=[0.1, 1.0],
    )
    original = solver.regularization.solve_linearized
    calls = 0

    def partly_failing(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise RuntimeError("synthetic factorization failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(
        solver.regularization,
        "solve_linearized",
        partly_failing,
    )
    result = solver.run()
    failure = next(
        item
        for item in result.rejected_candidates
        if item.reason is Occam1DRejectionReason.LINEAR_SOLVE_FAILED
    )
    assert failure.exception_type == "RuntimeError"
    assert "synthetic" in failure.message
    assert result.n_iterations == 1
    assert not result.failed_iterations


def test_all_failed_candidates_return_structured_failed_result(
    monkeypatch,
    tmp_path,
):
    data, model, config = _problem(target=0.01, iterations=2)
    solver = Occam1DInversion(
        data,
        model,
        config=config,
        multiplier_factors=[0.1, 1.0, 10.0],
    )

    def always_failing(*args, **kwargs):
        raise RuntimeError("synthetic solver outage")

    monkeypatch.setattr(
        solver.regularization,
        "solve_linearized",
        always_failing,
    )
    result = solver.run()
    assert result.convergence is Occam1DConvergence.FAILED
    assert result.n_iterations == 0
    assert len(result.failed_iterations) == 1
    failure = result.failed_iterations[0]
    assert failure.reason == "no_valid_candidate"
    assert failure.attempted_candidates == 3
    assert failure.rejected_candidates == 3
    assert failure.recoverable is True
    assert all(
        record.reason is Occam1DRejectionReason.LINEAR_SOLVE_FAILED
        for record in result.rejected_candidates
    )
    restored = Occam1DRestart.read(
        solver.restart().write(tmp_path / "failed-restart.json")
    )
    assert restored.rejected_candidates == result.rejected_candidates
    assert restored.failed_iterations == result.failed_iterations


def test_condition_reject_policy_flows_into_failure_audit():
    data, model, config = _problem(target=0.01, iterations=1)
    policy = Occam1DSolverPolicy(
        condition_limit=1.0,
        ill_condition_action="reject",
    )
    solver = Occam1DInversion(
        data,
        model,
        config=config,
        multiplier_factors=[0.1, 1.0],
        solver_policy=policy,
    )
    result = solver.run()
    assert result.convergence is Occam1DConvergence.FAILED
    assert result.failed_iterations[0].reason == "no_valid_candidate"
    assert all(
        item.exception_type == "Occam1DLinearSolveError"
        for item in result.rejected_candidates
    )
    assert solver.diagnostics()["solver_policy"][
        "ill_condition_action"
    ] == "reject"
