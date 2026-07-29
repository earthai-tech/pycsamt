# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Robust, resumable batch solving with retries and failure manifests.

:func:`solve_batch` is the single entry point. It solves many
:class:`~pycsamt.forward.maxwell.contracts.MaxwellProblem` instances
against one :class:`~pycsamt.forward.maxwell.backends.MaxwellBackend`,
optionally through a
:class:`~pycsamt.forward.maxwell.cache.MaxwellResultCache` for
resumability, and never lets one bad problem silently corrupt or halt
the whole run: every terminal failure is recorded in a
:class:`FailureManifest` instead of entering a training dataset
unnoticed.

Three concerns are deliberately kept separate:

Resumability
    Pass a :class:`~pycsamt.forward.maxwell.cache.MaxwellResultCache`;
    a repeated run over the same problems and cache directory skips
    everything already solved, including across process restarts.

Retries
    :class:`BatchPolicy` retries only exceptions considered transient
    (by default
    :class:`~pycsamt.forward.maxwell.adapters.BackendExecutionError`
    and
    :class:`~pycsamt.forward.maxwell.adapters.SolverConvergenceError`).
    A deterministic failure such as
    :class:`~pycsamt.forward.maxwell.adapters.IncompatibleProblemError`
    is recorded immediately without wasting attempts.

Failure manifests
    Every problem that exhausts its attempts becomes a
    :class:`ProblemFailure` inside the returned :class:`BatchReport`,
    which can be persisted with
    :meth:`FailureManifest.to_json_file` and inspected or filtered
    back out of a future run's input problems.
"""

from __future__ import annotations

import json
import re
import time
from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .adapters import BackendExecutionError, SolverConvergenceError
from .backends import MaxwellBackend
from .cache import MaxwellResultCache
from .contracts import ForwardResult, MaxwellProblem

__all__ = [
    "BatchAbortedError",
    "BatchPolicy",
    "ProblemFailure",
    "FailureManifest",
    "BatchReport",
    "solve_batch",
]

_HEX_64 = re.compile(r"^[0-9a-f]{64}$")


def _problem_hash(value: str, name: str) -> str:
    hashed = str(value).strip().lower()
    if not _HEX_64.match(hashed):
        raise ValueError(f"{name} must be a 64-character hex digest.")
    return hashed


@dataclass(frozen=True)
class BatchPolicy:
    """Configure retries and concurrency for :func:`solve_batch`.

    Parameters
    ----------
    max_attempts : int, default=1
        Attempts per problem, including the first. Values above one
        retry a solve that raised an exception in ``retry_on``.
    retry_backoff_s : float, default=1.0
        Base delay before each retry; attempt *n* (n > 1) waits
        ``retry_backoff_s * n`` seconds before the next attempt.
    retry_on : tuple of exception types, optional
        Exception types treated as transient and worth retrying.
        Defaults to
        :class:`~pycsamt.forward.maxwell.adapters.BackendExecutionError`
        and
        :class:`~pycsamt.forward.maxwell.adapters.SolverConvergenceError`
        (which covers
        :class:`~pycsamt.forward.maxwell.external.ExternalProcessError`
        as a subclass). Anything else, including
        :class:`~pycsamt.forward.maxwell.adapters.IncompatibleProblemError`
        and
        :class:`~pycsamt.forward.maxwell.adapters.InvalidBackendResultError`,
        is recorded as a terminal failure after its first occurrence.
        Note that
        :class:`~pycsamt.forward.maxwell.adapters.BaseMaxwellAdapter`
        wraps ordinary exceptions raised inside a solve into
        ``BackendExecutionError`` by default (its own
        ``AdapterPolicy.wrap_backend_exceptions``), so a deterministic
        bug in a backend's own code is retried like any other
        ``BackendExecutionError`` unless that adapter was built with
        ``wrap_backend_exceptions=False``.
    stop_on_first_failure : bool, default=False
        Raise :class:`BatchAbortedError` as soon as any problem
        exhausts its attempts, instead of recording it and continuing.
        With ``max_workers > 1`` this only stops further submissions;
        futures already in flight are still allowed to finish.
    max_workers : int, default=1
        Number of solves run concurrently in a thread pool. Values
        above one only help when the backend releases the GIL during
        its work (true of scipy sparse solves and of any
        :class:`~pycsamt.forward.maxwell.external.BaseExternalMaxwellAdapter`,
        which spends most of its time waiting on a subprocess).

    Examples
    --------
    >>> policy = BatchPolicy(max_attempts=3, retry_backoff_s=0.5)
    >>> policy.max_attempts
    3
    """

    max_attempts: int = 1
    retry_backoff_s: float = 1.0
    retry_on: tuple[type[Exception], ...] = (
        BackendExecutionError,
        SolverConvergenceError,
    )
    stop_on_first_failure: bool = False
    max_workers: int = 1

    def __post_init__(self) -> None:
        if (
            not isinstance(self.max_attempts, int)
            or isinstance(self.max_attempts, bool)
            or self.max_attempts < 1
        ):
            raise ValueError("max_attempts must be a positive integer.")
        backoff = float(self.retry_backoff_s)
        if not backoff == backoff or backoff < 0:  # NaN-safe check
            raise ValueError("retry_backoff_s must be finite and >= 0.")
        retry_on = tuple(self.retry_on)
        if not retry_on or not all(
            isinstance(value, type) and issubclass(value, Exception)
            for value in retry_on
        ):
            raise ValueError(
                "retry_on must be a non-empty tuple of exception types."
            )
        if (
            not isinstance(self.max_workers, int)
            or isinstance(self.max_workers, bool)
            or self.max_workers < 1
        ):
            raise ValueError("max_workers must be a positive integer.")
        object.__setattr__(self, "retry_backoff_s", backoff)
        object.__setattr__(self, "retry_on", retry_on)


@dataclass(frozen=True)
class ProblemFailure:
    """Record one problem's terminal solve failure.

    Parameters
    ----------
    problem_hash : str
        Identity of the
        :class:`~pycsamt.forward.maxwell.contracts.MaxwellProblem`
        that failed.
    attempts : int
        Number of attempts made before this failure was recorded.
    error_type : str
        Exception class name, for quick triage without deserializing.
    message : str
        Human-readable exception message.

    Examples
    --------
    >>> failure = ProblemFailure("a" * 64, 2, "BackendExecutionError", "x")
    >>> failure.attempts
    2
    """

    problem_hash: str
    attempts: int
    error_type: str
    message: str

    def __post_init__(self) -> None:
        problem_hash = _problem_hash(self.problem_hash, "problem_hash")
        attempts = int(self.attempts)
        if attempts < 1:
            raise ValueError("attempts must be a positive integer.")
        error_type = str(self.error_type).strip()
        message = str(self.message).strip()
        if not error_type or not message:
            raise ValueError("error_type and message cannot be empty.")
        object.__setattr__(self, "problem_hash", problem_hash)
        object.__setattr__(self, "attempts", attempts)
        object.__setattr__(self, "error_type", error_type)
        object.__setattr__(self, "message", message)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation.

        Returns
        -------
        dict
            Problem identity, attempt count, and exception details.

        Examples
        --------
        >>> ProblemFailure("a" * 64, 1, "X", "y").to_dict()["attempts"]
        1
        """
        return {
            "problem_hash": self.problem_hash,
            "attempts": self.attempts,
            "error_type": self.error_type,
            "message": self.message,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ProblemFailure:
        """Restore a validated failure record.

        Parameters
        ----------
        data : mapping
            State previously returned by :meth:`to_dict`.

        Returns
        -------
        ProblemFailure
            Restored, validated record.

        Examples
        --------
        >>> state = ProblemFailure("a" * 64, 1, "X", "y").to_dict()
        >>> ProblemFailure.from_dict(state).error_type
        'X'
        """
        return cls(
            data["problem_hash"],
            data["attempts"],
            data["error_type"],
            data["message"],
        )


@dataclass(frozen=True)
class FailureManifest:
    """Ordered, JSON-persistable record of every terminal batch failure.

    Parameters
    ----------
    failures : sequence of ProblemFailure, optional
        Failures in the order they were recorded. Problem hashes must
        be unique within one manifest.

    Examples
    --------
    >>> failure = ProblemFailure("a" * 64, 1, "X", "boom")
    >>> manifest = FailureManifest((failure,))
    >>> len(manifest), bool(manifest)
    (1, True)
    """

    failures: tuple[ProblemFailure, ...] = ()

    def __post_init__(self) -> None:
        failures = tuple(self.failures)
        if any(not isinstance(value, ProblemFailure) for value in failures):
            raise TypeError("failures must contain ProblemFailure values.")
        hashes = [value.problem_hash for value in failures]
        if len(set(hashes)) != len(hashes):
            raise ValueError("problem hashes must be unique in one manifest.")
        object.__setattr__(self, "failures", failures)

    @property
    def hashes(self) -> frozenset[str]:
        """Return every failed problem's hash.

        Returns
        -------
        frozenset of str
            Set suitable for filtering a future run's input problems.

        Examples
        --------
        >>> failure = ProblemFailure("a" * 64, 1, "X", "boom")
        >>> "a" * 64 in FailureManifest((failure,)).hashes
        True
        """
        return frozenset(value.problem_hash for value in self.failures)

    def __len__(self) -> int:
        return len(self.failures)

    def __bool__(self) -> bool:
        return bool(self.failures)

    def __contains__(self, problem: MaxwellProblem) -> bool:
        return problem.problem_hash in self.hashes

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable, schema-versioned representation.

        Returns
        -------
        dict
            Every failure, in recorded order.

        Examples
        --------
        >>> FailureManifest().to_dict()["schema_version"]
        1
        """
        return {
            "schema_version": 1,
            "failures": [value.to_dict() for value in self.failures],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FailureManifest:
        """Restore a validated manifest.

        Parameters
        ----------
        data : mapping
            State previously returned by :meth:`to_dict`.

        Returns
        -------
        FailureManifest
            Restored manifest.

        Raises
        ------
        ValueError
            If the schema version is unsupported.

        Examples
        --------
        >>> state = FailureManifest().to_dict()
        >>> FailureManifest.from_dict(state).failures
        ()
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported FailureManifest schema version.")
        return cls(
            tuple(
                ProblemFailure.from_dict(value)
                for value in data.get("failures", ())
            )
        )

    def to_json_file(self, path: str | Path) -> Path:
        """Write this manifest as indented, deterministic JSON.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination file.

        Returns
        -------
        pathlib.Path
            The destination path.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> failure = ProblemFailure("a" * 64, 1, "X", "boom")
        >>> manifest = FailureManifest((failure,))
        >>> with TemporaryDirectory() as directory:
        ...     target = Path(directory) / "failures.json"
        ...     _ = manifest.to_json_file(target)
        ...     restored = FailureManifest.from_json_file(target)
        >>> restored == manifest
        True
        """
        target = Path(path)
        target.write_text(
            json.dumps(self.to_dict(), sort_keys=True, indent=2),
            encoding="utf-8",
        )
        return target

    @classmethod
    def from_json_file(cls, path: str | Path) -> FailureManifest:
        """Load a manifest written by :meth:`to_json_file`.

        Parameters
        ----------
        path : str or pathlib.Path
            Source file.

        Returns
        -------
        FailureManifest
            Restored manifest.

        Examples
        --------
        See :meth:`to_json_file` for a complete round trip.
        """
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_dict(data)


@dataclass(frozen=True)
class BatchReport:
    """Summarize one :func:`solve_batch` run.

    Parameters
    ----------
    total : int
        Number of problems submitted to the batch. Can exceed
        ``len(solved) + len(failed)`` when
        ``BatchPolicy.stop_on_first_failure`` ended the run before
        every problem was attempted.
    solved : tuple of str
        Problem hashes with a valid result, whether freshly computed
        or already cached.
    cache_hits : tuple of str
        Subset of ``solved`` that were already present in the cache
        (skipped re-solving). Empty when no cache was used.
    failed : FailureManifest
        Every problem that exhausted its attempts.

    Examples
    --------
    >>> report = BatchReport(1, ("a" * 64,), (), FailureManifest())
    >>> report.success_fraction
    1.0
    """

    total: int
    solved: tuple[str, ...]
    cache_hits: tuple[str, ...]
    failed: FailureManifest

    def __post_init__(self) -> None:
        total = int(self.total)
        if total < 0:
            raise ValueError("total must be non-negative.")
        solved = tuple(self.solved)
        cache_hits = tuple(self.cache_hits)
        if not isinstance(self.failed, FailureManifest):
            raise TypeError("failed must be a FailureManifest.")
        if len(solved) + len(self.failed) > total:
            raise ValueError("len(solved) + len(failed) cannot exceed total.")
        if not set(cache_hits) <= set(solved):
            raise ValueError("cache_hits must be a subset of solved.")
        object.__setattr__(self, "total", total)
        object.__setattr__(self, "solved", solved)
        object.__setattr__(self, "cache_hits", cache_hits)

    @property
    def success_fraction(self) -> float:
        """Return the fraction of submitted problems that were solved.

        Returns
        -------
        float
            ``len(solved) / total``; ``1.0`` when ``total`` is zero.

        Examples
        --------
        >>> BatchReport(0, (), (), FailureManifest()).success_fraction
        1.0
        """
        if self.total == 0:
            return 1.0
        return len(self.solved) / self.total

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation.

        Returns
        -------
        dict
            Totals, solved/cache-hit hashes, and the failure manifest.

        Examples
        --------
        >>> BatchReport(0, (), (), FailureManifest()).to_dict()["total"]
        0
        """
        return {
            "total": self.total,
            "solved": list(self.solved),
            "cache_hits": list(self.cache_hits),
            "success_fraction": self.success_fraction,
            "failed": self.failed.to_dict(),
        }


class BatchAbortedError(RuntimeError):
    """Raised when ``BatchPolicy.stop_on_first_failure`` aborts early.

    Parameters
    ----------
    report : BatchReport
        Partial report covering every problem resolved before the
        abort; every problem submitted after the triggering failure
        (in the sequential case) was never attempted.

    Examples
    --------
    >>> report = BatchReport(1, (), (), FailureManifest())
    >>> error = BatchAbortedError(report)
    >>> error.report is report
    True
    """

    def __init__(self, report: BatchReport) -> None:
        super().__init__(
            f"batch aborted after {len(report.failed)} failure(s); "
            f"{len(report.solved)}/{report.total} problems solved "
            "before stopping."
        )
        self.report = report


def _solve_one(
    problem: MaxwellProblem,
    backend: MaxwellBackend,
    cache: MaxwellResultCache | None,
    policy: BatchPolicy,
) -> tuple[ForwardResult, bool] | ProblemFailure:
    attempts = 0
    last_error: Exception | None = None
    for attempt in range(1, policy.max_attempts + 1):
        attempts = attempt
        if attempt > 1:
            time.sleep(policy.retry_backoff_s * attempt)
        cached = cache is not None and cache.contains(problem)
        try:
            if cache is not None:
                result = cache.get_or_solve(problem, backend)
            else:
                result = backend.solve(problem)
            return result, cached
        except policy.retry_on as exc:
            last_error = exc
            continue
        except Exception as exc:  # recorded as a failure, not swallowed
            return ProblemFailure(
                problem.problem_hash, attempt, type(exc).__name__, str(exc)
            )
    assert last_error is not None
    return ProblemFailure(
        problem.problem_hash,
        attempts,
        type(last_error).__name__,
        str(last_error),
    )


def solve_batch(
    problems: Iterable[MaxwellProblem],
    backend: MaxwellBackend,
    *,
    cache: MaxwellResultCache | None = None,
    policy: BatchPolicy | None = None,
    on_result: Callable[[MaxwellProblem, ForwardResult], None] | None = None,
    on_failure: Callable[[MaxwellProblem, ProblemFailure], None] | None = None,
) -> BatchReport:
    """Solve many problems robustly, with retries and a failure manifest.

    Parameters
    ----------
    problems : iterable of MaxwellProblem
        Problems to solve. Consumed fully before returning; order is
        preserved in the returned report only through per-problem
        hashes, not positionally.
    backend : MaxwellBackend
        Conforming backend used for every problem.
    cache : MaxwellResultCache or None, optional
        When given, problems already cached are skipped (resumability)
        and freshly computed results are written back to it.
    policy : BatchPolicy or None, optional
        Retry, concurrency, and abort configuration. Defaults to
        :class:`BatchPolicy`.
    on_result : callable, optional
        Invoked with ``(problem, result)`` for every problem that
        succeeds, in completion order. Useful for streaming results
        into a dataset without holding them all in memory.
    on_failure : callable, optional
        Invoked with ``(problem, failure)`` for every problem that
        exhausts its attempts, in completion order.

    Returns
    -------
    BatchReport
        Totals, solved/cache-hit problem hashes, and a
        :class:`FailureManifest` of every terminal failure.

    Raises
    ------
    TypeError
        If ``backend`` or ``policy`` has the wrong type.
    BatchAbortedError
        If ``policy.stop_on_first_failure`` is set and any problem
        exhausts its attempts.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.maxwell import (
    ...     CallableMaxwellAdapter,
    ...     BackendCapabilities,
    ...     ForwardResult,
    ...     MaxwellMesh,
    ...     MaxwellProblem,
    ...     ReceiverSet,
    ...     SolverDiagnostics,
    ... )
    >>> mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
    >>> problem = MaxwellProblem(
    ...     mesh, np.ones((2, 2)), [1], ReceiverSet([[0.5, 0]], ["S"])
    ... )
    >>> def solver(value):
    ...     diagnostics = SolverDiagnostics([[True]], [[0]], [[0]], 0)
    ...     return ForwardResult(
    ...         value.problem_hash,
    ...         value.frequencies_hz,
    ...         value.receivers.names,
    ...         value.components,
    ...         [[[1j, 1j]]],
    ...         None,
    ...         "demo",
    ...         "1",
    ...         diagnostics,
    ...     )
    >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy", "zyx"))
    >>> backend = CallableMaxwellAdapter(cap, solver)
    >>> report = solve_batch([problem], backend)
    >>> report.success_fraction
    1.0
    """
    if not isinstance(backend, MaxwellBackend):
        raise TypeError("backend must implement MaxwellBackend.")
    policy = BatchPolicy() if policy is None else policy
    if not isinstance(policy, BatchPolicy):
        raise TypeError("policy must be a BatchPolicy or None.")

    problem_list = list(problems)
    solved: list[str] = []
    cache_hits: list[str] = []
    failures: list[ProblemFailure] = []

    def _handle(problem: MaxwellProblem, outcome) -> None:
        if isinstance(outcome, ProblemFailure):
            failures.append(outcome)
            if on_failure is not None:
                on_failure(problem, outcome)
            return
        result, cached = outcome
        solved.append(problem.problem_hash)
        if cached:
            cache_hits.append(problem.problem_hash)
        if on_result is not None:
            on_result(problem, result)

    if policy.max_workers <= 1:
        for problem in problem_list:
            outcome = _solve_one(problem, backend, cache, policy)
            _handle(problem, outcome)
            if policy.stop_on_first_failure and failures:
                break
    else:
        with ThreadPoolExecutor(max_workers=policy.max_workers) as pool:
            pending = {}
            aborting = False
            for problem in problem_list:
                if aborting:
                    break
                future = pool.submit(
                    _solve_one, problem, backend, cache, policy
                )
                pending[future] = problem
            for future in as_completed(pending):
                _handle(pending[future], future.result())
                if policy.stop_on_first_failure and failures:
                    aborting = True

    report = BatchReport(
        len(problem_list),
        tuple(solved),
        tuple(cache_hits),
        FailureManifest(tuple(failures)),
    )
    if policy.stop_on_first_failure and failures:
        raise BatchAbortedError(report)
    return report
