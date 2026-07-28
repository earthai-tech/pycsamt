"""Tests for validated Maxwell adapter execution."""

import warnings

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    AdapterPolicy,
    BackendCapabilities,
    BackendExecutionError,
    CallableMaxwellAdapter,
    ForwardResult,
    IncompatibleProblemError,
    InvalidBackendResultError,
    MaxwellBackend,
    MaxwellMesh,
    MaxwellProblem,
    ReceiverSet,
    SolverConvergenceError,
    SolverDiagnostics,
)


def _problem(frequencies=(10, 1)):
    mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
    return MaxwellProblem(
        mesh,
        np.ones(mesh.shape),
        frequencies,
        ReceiverSet([[0.5, 0]], ["S"]),
        ("zxy",),
    )


def _capabilities(**changes):
    values = dict(
        name="demo",
        version="1",
        dimensions=(2,),
        components=("zxy",),
        verified_benchmarks=("half-space",),
    )
    values.update(changes)
    return BackendCapabilities(**values)


def _result(
    problem,
    *,
    converged=True,
    residual=1e-9,
    valid=True,
    name="demo",
    version="1",
):
    count = len(problem.frequencies_hz)
    diagnostics = SolverDiagnostics(
        np.full((count, 1), converged),
        np.ones((count, 1), int),
        np.full((count, 1), residual),
        0.1,
    )
    return ForwardResult(
        problem.problem_hash,
        problem.frequencies_hz,
        problem.receivers.names,
        problem.components,
        np.full((1, count, 1), 1 + 2j),
        np.full((1, count, 1), valid),
        name,
        version,
        diagnostics,
    )


def test_callable_adapter_solves_and_matches_runtime_protocol():
    adapter = CallableMaxwellAdapter(_capabilities(), lambda problem: _result(problem))
    result = adapter.solve(_problem())
    assert isinstance(adapter, MaxwellBackend)
    assert result.success


def test_preflight_rejects_before_callback_is_called():
    called = []
    adapter = CallableMaxwellAdapter(
        _capabilities(dimensions=(3,)), lambda problem: called.append(problem)
    )
    with pytest.raises(IncompatibleProblemError, match="2-D problems"):
        adapter.solve(_problem())
    assert called == []


def test_capability_warnings_can_be_emitted_or_suppressed():
    problem = _problem()
    cap = _capabilities(verified_benchmarks=())
    with pytest.warns(RuntimeWarning, match="no verified benchmarks"):
        CallableMaxwellAdapter(cap, lambda value: _result(value)).solve(problem)
    policy = AdapterPolicy(emit_capability_warnings=False)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        CallableMaxwellAdapter(cap, lambda value: _result(value), policy).solve(problem)
    assert not captured


def test_backend_exceptions_are_chained_or_preserved_by_policy():
    def fail(problem):
        raise OSError("native failure")

    with pytest.raises(BackendExecutionError, match="native failure") as captured:
        CallableMaxwellAdapter(_capabilities(), fail).solve(_problem())
    assert isinstance(captured.value.__cause__, OSError)
    policy = AdapterPolicy(wrap_backend_exceptions=False)
    with pytest.raises(OSError, match="native failure"):
        CallableMaxwellAdapter(_capabilities(), fail, policy).solve(_problem())


@pytest.mark.parametrize(
    "callback, message",
    [
        (lambda problem: object(), "must return a ForwardResult"),
        (lambda problem: _result(problem, name="other"), "identity differs"),
        (lambda problem: _result(_problem((3, 2))), "problem_hash"),
    ],
)
def test_invalid_backend_results_are_rejected(callback, message):
    with pytest.raises(InvalidBackendResultError, match=message):
        CallableMaxwellAdapter(_capabilities(), callback).solve(_problem())


def test_convergence_residual_and_validity_policies():
    problem = _problem()
    with pytest.raises(SolverConvergenceError, match="unconverged"):
        CallableMaxwellAdapter(
            _capabilities(), lambda value: _result(value, converged=False)
        ).solve(problem)
    residual_policy = AdapterPolicy(maximum_relative_residual=1e-6)
    with pytest.raises(SolverConvergenceError, match="exceeds"):
        CallableMaxwellAdapter(
            _capabilities(),
            lambda value: _result(value, residual=1e-3),
            residual_policy,
        ).solve(problem)
    with pytest.raises(SolverConvergenceError, match="invalid impedance"):
        CallableMaxwellAdapter(
            _capabilities(), lambda value: _result(value, valid=False)
        ).solve(problem)
    relaxed = AdapterPolicy(require_convergence=False, require_all_valid=False)
    result = CallableMaxwellAdapter(
        _capabilities(),
        lambda value: _result(value, converged=False, valid=False),
        relaxed,
    ).solve(problem)
    assert not result.success


def test_solve_many_preserves_order_and_stops_on_failure():
    seen = []

    def solve(problem):
        seen.append(problem.problem_hash)
        return _result(problem)

    problems = [_problem((10,)), _problem((1,))]
    results = CallableMaxwellAdapter(_capabilities(), solve).solve_many(problems)
    assert [result.problem_hash for result in results] == seen


def test_adapter_policy_and_constructor_validation():
    with pytest.raises(ValueError, match="non-negative"):
        AdapterPolicy(maximum_relative_residual=-1)
    with pytest.raises(TypeError, match="callable"):
        CallableMaxwellAdapter(_capabilities(), None)
