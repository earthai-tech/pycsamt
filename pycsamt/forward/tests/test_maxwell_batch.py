"""Tests for robust batch solving, retries, and failure manifests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    AdapterPolicy,
    BackendCapabilities,
    BackendExecutionError,
    BatchAbortedError,
    BatchPolicy,
    BatchReport,
    CallableMaxwellAdapter,
    FailureManifest,
    ForwardResult,
    IncompatibleProblemError,
    MaxwellBackend,
    MaxwellMesh,
    MaxwellProblem,
    MaxwellResultCache,
    ProblemFailure,
    ReceiverSet,
    SolverDiagnostics,
    solve_batch,
)


def _problem(name="S", frequencies=(10, 1)):
    mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
    return MaxwellProblem(
        mesh,
        np.ones(mesh.shape),
        frequencies,
        ReceiverSet([[0.5, 0]], [name]),
        ("zxy",),
    )


def _capabilities(**changes):
    values = dict(
        name="demo", version="1", dimensions=(2,), components=("zxy",)
    )
    values.update(changes)
    return BackendCapabilities(**values)


def _result(problem, *, converged=True):
    count = len(problem.frequencies_hz)
    diagnostics = SolverDiagnostics(
        np.full((count, 1), converged),
        np.ones((count, 1), int),
        np.full((count, 1), 1e-9),
        0.01,
    )
    return ForwardResult(
        problem.problem_hash,
        problem.frequencies_hz,
        problem.receivers.names,
        problem.components,
        np.full((1, count, 1), 1 + 2j),
        None,
        "demo",
        "1",
        diagnostics,
    )


class _FlakyAdapter(CallableMaxwellAdapter):
    """Fails ``fail_times`` times per distinct problem, then succeeds."""

    def __init__(self, fail_times=0, exception=None, policy=None):
        self._fail_times = fail_times
        self._exception = exception or (
            lambda: BackendExecutionError("transient boom")
        )
        self.calls: dict[str, int] = {}
        super().__init__(_capabilities(), self._solve, policy)

    def _solve(self, problem):
        count = self.calls.get(problem.problem_hash, 0)
        self.calls[problem.problem_hash] = count + 1
        if count < self._fail_times:
            raise self._exception()
        return _result(problem)


# --------------------------------------------------------------------------- #
# BatchPolicy
# --------------------------------------------------------------------------- #


def test_batch_policy_defaults():
    policy = BatchPolicy()
    assert policy.max_attempts == 1
    assert policy.max_workers == 1
    assert BackendExecutionError in policy.retry_on


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_attempts": 0},
        {"retry_backoff_s": -1.0},
        {"retry_backoff_s": float("nan")},
        {"retry_on": ()},
        {"retry_on": (int,)},
        {"max_workers": 0},
    ],
)
def test_batch_policy_rejects_invalid_values(kwargs):
    with pytest.raises(ValueError):
        BatchPolicy(**kwargs)


# --------------------------------------------------------------------------- #
# ProblemFailure / FailureManifest
# --------------------------------------------------------------------------- #


def test_problem_failure_round_trips_through_dict():
    failure = ProblemFailure("a" * 64, 2, "BackendExecutionError", "boom")
    restored = ProblemFailure.from_dict(failure.to_dict())
    assert restored == failure


@pytest.mark.parametrize(
    "kwargs",
    [
        {"problem_hash": "not-a-hash"},
        {"attempts": 0},
        {"error_type": ""},
        {"message": "  "},
    ],
)
def test_problem_failure_rejects_invalid_fields(kwargs):
    values = dict(
        problem_hash="a" * 64, attempts=1, error_type="X", message="y"
    )
    values.update(kwargs)
    with pytest.raises(ValueError):
        ProblemFailure(**values)


def test_failure_manifest_rejects_duplicate_hashes():
    failure = ProblemFailure("a" * 64, 1, "X", "y")
    with pytest.raises(ValueError, match="unique"):
        FailureManifest((failure, failure))


def test_failure_manifest_len_bool_contains_hashes():
    failure = ProblemFailure("a" * 64, 1, "X", "y")
    manifest = FailureManifest((failure,))
    assert len(manifest) == 1
    assert bool(manifest) is True
    assert manifest.hashes == frozenset({"a" * 64})
    problem = _problem()
    assert problem not in manifest
    assert bool(FailureManifest()) is False


def test_failure_manifest_json_round_trip(tmp_path: Path):
    failure = ProblemFailure("a" * 64, 3, "SolverConvergenceError", "bad")
    manifest = FailureManifest((failure,))
    target = manifest.to_json_file(tmp_path / "failures.json")
    restored = FailureManifest.from_json_file(target)
    assert restored == manifest


def test_failure_manifest_from_dict_rejects_unknown_schema():
    with pytest.raises(ValueError, match="schema"):
        FailureManifest.from_dict({"schema_version": 99})


# --------------------------------------------------------------------------- #
# BatchReport
# --------------------------------------------------------------------------- #


def test_batch_report_success_fraction_and_empty_total():
    report = BatchReport(0, (), (), FailureManifest())
    assert report.success_fraction == 1.0
    report2 = BatchReport(2, ("a" * 64,), (), FailureManifest())
    assert report2.success_fraction == 0.5


def test_batch_report_rejects_inconsistent_counts():
    failure = ProblemFailure("b" * 64, 1, "X", "y")
    with pytest.raises(ValueError, match="cannot exceed total"):
        BatchReport(0, ("a" * 64,), (), FailureManifest((failure,)))


def test_batch_report_rejects_cache_hits_not_subset_of_solved():
    with pytest.raises(ValueError, match="subset"):
        BatchReport(1, ("a" * 64,), ("b" * 64,), FailureManifest())


def test_batch_aborted_error_carries_partial_report():
    report = BatchReport(1, (), (), FailureManifest())
    error = BatchAbortedError(report)
    assert error.report is report


# --------------------------------------------------------------------------- #
# solve_batch: type checks and the plain happy path
# --------------------------------------------------------------------------- #


def test_solve_batch_rejects_wrong_backend_type():
    with pytest.raises(TypeError, match="MaxwellBackend"):
        solve_batch([_problem()], object())


def test_solve_batch_rejects_wrong_policy_type():
    adapter = CallableMaxwellAdapter(_capabilities(), _result)
    with pytest.raises(TypeError, match="BatchPolicy"):
        solve_batch([_problem()], adapter, policy=object())


def test_solve_batch_solves_every_problem_without_a_cache():
    adapter = CallableMaxwellAdapter(_capabilities(), _result)
    problems = [_problem(name=f"S{i}") for i in range(3)]
    report = solve_batch(problems, adapter)
    assert report.total == 3
    assert len(report.solved) == 3
    assert report.cache_hits == ()
    assert not report.failed


def test_solve_batch_invokes_callbacks():
    adapter = CallableMaxwellAdapter(_capabilities(), _result)
    seen_results = []
    seen_failures = []
    solve_batch(
        [_problem()],
        adapter,
        on_result=lambda problem, result: seen_results.append(result),
        on_failure=lambda problem, failure: seen_failures.append(failure),
    )
    assert len(seen_results) == 1
    assert seen_failures == []


# --------------------------------------------------------------------------- #
# solve_batch: cache-based resumability
# --------------------------------------------------------------------------- #


def test_solve_batch_resumes_from_cache_without_recalling_backend(
    tmp_path: Path,
):
    adapter = _FlakyAdapter()
    cache = MaxwellResultCache(tmp_path / "cache")
    problem = _problem()

    first = solve_batch([problem], adapter, cache=cache)
    assert first.cache_hits == ()
    assert adapter.calls[problem.problem_hash] == 1

    second = solve_batch([problem], adapter, cache=cache)
    assert second.cache_hits == (problem.problem_hash,)
    assert adapter.calls[problem.problem_hash] == 1  # backend not re-called


# --------------------------------------------------------------------------- #
# solve_batch: retries
# --------------------------------------------------------------------------- #


def test_solve_batch_retries_transient_failures_until_success():
    adapter = _FlakyAdapter(fail_times=2)
    problem = _problem()
    report = solve_batch(
        [problem], adapter, policy=BatchPolicy(max_attempts=3, retry_backoff_s=0)
    )
    assert report.solved == (problem.problem_hash,)
    assert adapter.calls[problem.problem_hash] == 3


def test_solve_batch_records_failure_after_exhausting_attempts():
    adapter = _FlakyAdapter(fail_times=5)
    problem = _problem()
    report = solve_batch(
        [problem], adapter, policy=BatchPolicy(max_attempts=2, retry_backoff_s=0)
    )
    assert report.solved == ()
    assert len(report.failed) == 1
    failure = report.failed.failures[0]
    assert failure.attempts == 2
    assert failure.error_type == "BackendExecutionError"


def test_solve_batch_does_not_retry_non_transient_exceptions():
    # AdapterPolicy(wrap_backend_exceptions=False) lets the raw ValueError
    # through unwrapped; otherwise BaseMaxwellAdapter.solve() would coerce
    # it into a (retryable) BackendExecutionError before batch.py sees it.
    adapter = _FlakyAdapter(
        fail_times=10,
        exception=lambda: ValueError("permanent"),
        policy=AdapterPolicy(wrap_backend_exceptions=False),
    )
    problem = _problem()
    report = solve_batch(
        [problem], adapter, policy=BatchPolicy(max_attempts=5, retry_backoff_s=0)
    )
    assert len(report.failed) == 1
    assert report.failed.failures[0].attempts == 1
    assert report.failed.failures[0].error_type == "ValueError"
    assert adapter.calls[problem.problem_hash] == 1


def test_solve_batch_records_incompatible_problem_immediately():
    adapter = CallableMaxwellAdapter(_capabilities(dimensions=(3,)), _result)
    report = solve_batch(
        [_problem()], adapter, policy=BatchPolicy(max_attempts=4)
    )
    assert len(report.failed) == 1
    assert report.failed.failures[0].error_type == "IncompatibleProblemError"
    assert report.failed.failures[0].attempts == 1


# --------------------------------------------------------------------------- #
# solve_batch: partial batches, one bad problem does not sink the rest
# --------------------------------------------------------------------------- #


def test_solve_batch_continues_past_one_failing_problem_by_default():
    good = _problem(name="good")
    bad = _problem(name="bad", frequencies=(10, 1, 0.5))

    def solver(problem):
        if problem.receivers.names == ("bad",):
            raise BackendExecutionError("boom")
        return _result(problem)

    adapter = CallableMaxwellAdapter(_capabilities(), solver)
    report = solve_batch([good, bad], adapter, policy=BatchPolicy(max_attempts=1))
    assert report.total == 2
    assert report.solved == (good.problem_hash,)
    assert len(report.failed) == 1


def test_solve_batch_stop_on_first_failure_aborts_and_skips_remaining():
    calls = []

    def solver(problem):
        calls.append(problem.receivers.names)
        if problem.receivers.names == ("first",):
            raise BackendExecutionError("boom")
        return _result(problem)

    adapter = CallableMaxwellAdapter(_capabilities(), solver)
    problems = [_problem(name="first"), _problem(name="second")]
    policy = BatchPolicy(max_attempts=1, stop_on_first_failure=True)
    with pytest.raises(BatchAbortedError) as excinfo:
        solve_batch(problems, adapter, policy=policy)
    assert excinfo.value.report.total == 2
    assert len(excinfo.value.report.failed) == 1
    assert ("second",) not in calls  # never attempted


# --------------------------------------------------------------------------- #
# solve_batch: concurrency
# --------------------------------------------------------------------------- #


def test_solve_batch_with_multiple_workers_solves_everything():
    adapter = CallableMaxwellAdapter(_capabilities(), _result)
    problems = [_problem(name=f"S{i}") for i in range(8)]
    report = solve_batch(
        problems, adapter, policy=BatchPolicy(max_workers=4)
    )
    assert report.total == 8
    assert set(report.solved) == {p.problem_hash for p in problems}
    assert not report.failed


def test_solve_batch_with_multiple_workers_still_records_failures():
    def solver(problem):
        if problem.receivers.names == ("bad",):
            raise BackendExecutionError("boom")
        return _result(problem)

    adapter = CallableMaxwellAdapter(_capabilities(), solver)
    problems = [_problem(name=f"S{i}") for i in range(4)] + [
        _problem(name="bad")
    ]
    report = solve_batch(
        problems,
        adapter,
        policy=BatchPolicy(max_workers=3, max_attempts=1),
    )
    assert report.total == 5
    assert len(report.solved) == 4
    assert len(report.failed) == 1


# --------------------------------------------------------------------------- #
# Integration: filtering a rerun by a previous failure manifest
# --------------------------------------------------------------------------- #


def test_failure_manifest_can_filter_a_rerun():
    def solver(problem):
        if problem.receivers.names == ("bad",):
            raise BackendExecutionError("boom")
        return _result(problem)

    adapter = CallableMaxwellAdapter(_capabilities(), solver)
    problems = [_problem(name="good"), _problem(name="bad")]
    first = solve_batch(problems, adapter, policy=BatchPolicy(max_attempts=1))
    assert len(first.failed) == 1

    remaining = [p for p in problems if p not in first.failed]
    assert len(remaining) == 1
    assert remaining[0].receivers.names == ("good",)
