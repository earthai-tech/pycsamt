"""Tests for the content-addressed Maxwell result cache."""

import os
import time

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    BackendCapabilities,
    CacheCorruptionError,
    CacheLockTimeoutError,
    CallableMaxwellAdapter,
    ForwardResult,
    MaxwellMesh,
    MaxwellProblem,
    MaxwellResultCache,
    ReceiverSet,
    SolverDiagnostics,
)


def _problem(frequency=1.0):
    mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
    receivers = ReceiverSet([[0.5, 0]], ["S"])
    return MaxwellProblem(
        mesh,
        np.ones(mesh.shape),
        [frequency],
        receivers,
        ("zxy",),
    )


def _result(problem, value=1.0):
    diagnostics = SolverDiagnostics([[True]], [[0]], [[0]], 0)
    impedance = np.array([[[complex(value, 1)]]])
    return ForwardResult(
        problem.problem_hash,
        problem.frequencies_hz,
        problem.receivers.names,
        problem.components,
        impedance,
        None,
        "demo",
        "1",
        diagnostics,
    )


def _backend(calls):
    capabilities = BackendCapabilities(
        "demo",
        "1",
        (2,),
        ("zxy",),
        verified_benchmarks=("half-space",),
    )

    def solve(problem):
        calls.append(problem.problem_hash)
        return _result(problem)

    return CallableMaxwellAdapter(capabilities, solve)


def test_put_get_entry_and_statistics(tmp_path):
    cache = MaxwellResultCache(tmp_path / "cache")
    problem = _problem()
    assert cache.get(problem) is None
    entry = cache.put(problem, _result(problem))
    assert entry.key == problem.problem_hash
    assert cache.contains(problem)
    restored = cache.get(problem)
    np.testing.assert_array_equal(
        restored.impedance_v_a,
        _result(problem).impedance_v_a,
    )
    statistics = cache.statistics()
    assert statistics.entry_count == 1
    assert statistics.total_bytes == entry.size_bytes


def test_get_or_solve_uses_cache_after_first_call(tmp_path):
    calls = []
    cache = MaxwellResultCache(tmp_path)
    problem = _problem()
    backend = _backend(calls)
    first = cache.get_or_solve(problem, backend)
    second = cache.get_or_solve(problem, backend)
    assert first.problem_hash == second.problem_hash
    assert calls == [problem.problem_hash]


def test_put_retains_existing_unless_overwrite_requested(tmp_path):
    cache = MaxwellResultCache(tmp_path)
    problem = _problem()
    cache.put(problem, _result(problem, 1))
    cache.put(problem, _result(problem, 2))
    assert cache.get(problem).impedance_v_a.real.item() == 1
    cache.put(problem, _result(problem, 2), overwrite=True)
    assert cache.get(problem).impedance_v_a.real.item() == 2


def test_checksum_corruption_is_quarantined(tmp_path):
    cache = MaxwellResultCache(tmp_path)
    problem = _problem()
    entry = cache.put(problem, _result(problem))
    with entry.archive_path.open("ab") as stream:
        stream.write(b"damage")
    assert cache.get(problem) is None
    statistics = cache.statistics()
    assert statistics.entry_count == 0
    assert statistics.corrupt_count == 2


def test_corruption_can_be_reported_without_quarantine(tmp_path):
    cache = MaxwellResultCache(
        tmp_path,
        quarantine_corrupt=False,
    )
    problem = _problem()
    entry = cache.put(problem, _result(problem))
    entry.checksum_path.write_text("invalid\n", encoding="ascii")
    with pytest.raises(CacheCorruptionError, match="sidecar"):
        cache.get(problem)
    assert entry.archive_path.exists()


def test_remove_and_prune_are_scoped_to_entries(tmp_path):
    cache = MaxwellResultCache(tmp_path)
    first = _problem(1)
    second = _problem(2)
    cache.put(first, _result(first))
    time.sleep(0.01)
    cache.put(second, _result(second))
    removed = cache.prune(0)
    assert [value.key for value in removed] == [
        first.problem_hash,
        second.problem_hash,
    ]
    assert cache.statistics().entry_count == 0
    assert not cache.remove(first)


def test_orphan_files_are_counted(tmp_path):
    cache = MaxwellResultCache(tmp_path)
    problem = _problem()
    archive, _ = cache._paths(problem.problem_hash)
    archive.parent.mkdir(parents=True, exist_ok=True)
    archive.write_bytes(b"orphan")
    statistics = cache.statistics()
    assert statistics.entry_count == 0
    assert statistics.orphan_count == 1


def test_lock_timeout_and_stale_lock_recovery(tmp_path):
    cache = MaxwellResultCache(
        tmp_path,
        lock_timeout_s=0.02,
        poll_interval_s=0.005,
        stale_lock_s=10,
    )
    problem = _problem()
    lock = cache._locks / f"{problem.problem_hash}.lock"
    lock.write_text("busy", encoding="ascii")
    with pytest.raises(CacheLockTimeoutError, match="timed out"):
        cache.put(problem, _result(problem))
    old = time.time() - 20
    os.utime(lock, (old, old))
    cache.put(problem, _result(problem))
    assert cache.contains(problem)


def test_problem_result_and_configuration_validation(tmp_path):
    cache = MaxwellResultCache(tmp_path)
    first = _problem(1)
    second = _problem(2)
    with pytest.raises(ValueError, match="problem_hash"):
        cache.put(first, _result(second))
    with pytest.raises(ValueError, match="positive"):
        MaxwellResultCache(tmp_path / "bad", lock_timeout_s=0)
    with pytest.raises(ValueError, match="non-negative"):
        cache.prune(-1)
