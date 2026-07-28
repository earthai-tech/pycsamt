"""Tests for canonical Maxwell benchmark definitions."""

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    BackendCapabilities,
    BenchmarkThresholds,
    CallableMaxwellAdapter,
    ForwardResult,
    MaxwellMesh,
    ReceiverSet,
    SolverDiagnostics,
    half_space_benchmark,
    half_space_impedance,
    layered_earth_benchmark,
    layered_earth_impedance,
    run_benchmarks,
)


def _geometry():
    mesh = MaxwellMesh([0, 100, 200], [0, 100, 300])
    receivers = ReceiverSet([[50, 0], [150, 0]], ["A", "B"])
    return mesh, receivers


def _result(case, values=None, valid=None, converged=True):
    problem = case.problem
    count = len(problem.frequencies_hz)
    diagnostics = SolverDiagnostics(
        np.full((count, 1), converged),
        np.zeros((count, 1), int),
        np.zeros((count, 1)),
        0,
    )
    impedance = case.reference_impedance_v_a if values is None else values
    return ForwardResult(
        problem.problem_hash,
        problem.frequencies_hz,
        problem.receivers.names,
        problem.components,
        impedance,
        valid,
        "demo",
        "1",
        diagnostics,
    )


def _backend(cases):
    capabilities = BackendCapabilities(
        "demo",
        "1",
        (2,),
        ("zxy", "zyx"),
        verified_benchmarks=("half-space",),
    )
    references = {
        case.problem.problem_hash: case.reference_impedance_v_a for case in cases
    }

    def solve(problem):
        count = len(problem.frequencies_hz)
        diagnostics = SolverDiagnostics(
            np.ones((count, 1), bool),
            np.zeros((count, 1), int),
            np.zeros((count, 1)),
            0,
        )
        return ForwardResult(
            problem.problem_hash,
            problem.frequencies_hz,
            problem.receivers.names,
            problem.components,
            references[problem.problem_hash],
            None,
            "demo",
            "1",
            diagnostics,
        )

    return CallableMaxwellAdapter(capabilities, solve)


def test_half_space_impedance_conventions_and_scaling():
    positive = half_space_impedance(100, [1, 4])
    negative = half_space_impedance(
        100,
        [1, 4],
        time_dependence="exp(-iwt)",
    )
    np.testing.assert_allclose(negative, np.conjugate(positive))
    assert abs(positive[1] / positive[0]) == pytest.approx(2)
    assert np.angle(positive[0], deg=True) == pytest.approx(45)


def test_layered_single_halfspace_matches_closed_form():
    frequencies = [0.1, 1, 10]
    expected = half_space_impedance(50, frequencies)
    actual = layered_earth_impedance([50], [], frequencies)
    np.testing.assert_allclose(actual, expected)


def test_half_space_case_axes_signs_hash_and_exact_evaluation():
    mesh, receivers = _geometry()
    case = half_space_benchmark(mesh, receivers, [10, 1])
    np.testing.assert_allclose(
        case.reference_impedance_v_a[:, :, 1],
        -case.reference_impedance_v_a[:, :, 0],
    )
    assert len(case.benchmark_hash) == 64
    outcome = case.evaluate(_result(case))
    assert outcome.passed
    assert outcome.metrics.normalized_rms == 0
    assert outcome.to_dict()["passed"]


def test_evaluation_reports_amplitude_phase_validity_and_convergence():
    mesh, receivers = _geometry()
    limits = BenchmarkThresholds(
        maximum_normalized_rms=0.01,
        maximum_amplitude_relative_error=0.01,
        maximum_phase_error_deg=0.5,
        minimum_valid_fraction=1,
    )
    case = half_space_benchmark(
        mesh,
        receivers,
        [1],
        thresholds=limits,
    )
    changed = case.reference_impedance_v_a * 1.2 * np.exp(1j * 0.1)
    valid = np.ones(changed.shape, bool)
    valid[0, 0, 0] = False
    outcome = case.evaluate(_result(case, changed, valid=valid, converged=False))
    assert not outcome.passed
    assert len(outcome.failures) == 5
    no_values = case.evaluate(
        _result(
            case,
            case.reference_impedance_v_a,
            valid=np.zeros(case.reference_impedance_v_a.shape, bool),
        )
    )
    assert not no_values.passed
    assert np.isfinite(no_values.metrics.normalized_rms)


def test_layered_case_maps_lateral_uniform_conductivity():
    mesh = MaxwellMesh([0, 100, 200], [0, 100, 300])
    receivers = ReceiverSet([[50, 0]], ["A"])
    case = layered_earth_benchmark(
        mesh,
        receivers,
        [1],
        [10, 100],
        [100],
    )
    np.testing.assert_allclose(
        case.problem.conductivity_s_m,
        [[0.1, 0.1], [0.01, 0.01]],
    )
    with pytest.raises(ValueError, match="not a mesh z edge"):
        layered_earth_benchmark(
            mesh,
            receivers,
            [1],
            [10, 100],
            [50],
        )


def test_run_benchmarks_preserves_order_and_aggregates():
    mesh, receivers = _geometry()
    first = half_space_benchmark(mesh, receivers, [1])
    second = layered_earth_benchmark(
        mesh,
        receivers,
        [1],
        [10, 100],
        [100],
    )
    report = run_benchmarks(_backend((first, second)), (first, second))
    assert report.passed
    assert report.pass_fraction == 1
    assert [value.benchmark_name for value in report.outcomes] == [
        "half-space",
        "layered-earth",
    ]


def test_benchmark_input_validation():
    mesh, receivers = _geometry()
    with pytest.raises(ValueError, match="zxy and/or zyx"):
        half_space_benchmark(
            MaxwellMesh(
                [0, 1, 2],
                [0, 1, 2],
                [0, 1, 2],
            ),
            ReceiverSet([[0, 0, 0]], ["S"]),
            [1],
            components=("zxx",),
        )
    with pytest.raises(ValueError, match="omit only"):
        layered_earth_impedance([10, 100], [], [1])
    case = half_space_benchmark(mesh, receivers, [1])
    with pytest.raises(ValueError, match="unique"):
        run_benchmarks(_backend((case,)), (case, case))
