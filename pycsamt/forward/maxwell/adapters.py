# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Validated execution layer for solver-specific Maxwell adapters.

Subclasses implement only the backend call and conversion to
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult`.  The public solve
path is fixed here so capability, axis, provenance, and convergence checks
cannot be accidentally skipped by individual integrations.
"""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .backends import BackendCapabilities, CompatibilityReport
from .contracts import ForwardResult, MaxwellProblem

__all__ = [
    "MaxwellAdapterError",
    "IncompatibleProblemError",
    "BackendExecutionError",
    "InvalidBackendResultError",
    "SolverConvergenceError",
    "AdapterPolicy",
    "BaseMaxwellAdapter",
    "CallableMaxwellAdapter",
]


class MaxwellAdapterError(RuntimeError):
    """Base exception raised by the validated adapter execution layer.

    Examples
    --------
    >>> error = MaxwellAdapterError("solver failed")
    >>> str(error)
    'solver failed'
    """


class IncompatibleProblemError(MaxwellAdapterError):
    """Indicate that declared backend capabilities reject a problem.

    Examples
    --------
    >>> isinstance(
    ...     IncompatibleProblemError("unsupported"), MaxwellAdapterError
    ... )
    True
    """


class BackendExecutionError(MaxwellAdapterError):
    """Wrap an exception raised inside a numerical backend.

    Examples
    --------
    >>> isinstance(BackendExecutionError("failed"), MaxwellAdapterError)
    True
    """


class InvalidBackendResultError(MaxwellAdapterError):
    """Indicate malformed, mislabeled, or mismatched backend output.

    Examples
    --------
    >>> isinstance(InvalidBackendResultError("bad axes"), MaxwellAdapterError)
    True
    """


class SolverConvergenceError(MaxwellAdapterError):
    """Indicate that a valid result violates the convergence policy.

    Examples
    --------
    >>> isinstance(
    ...     SolverConvergenceError("residual too large"), MaxwellAdapterError
    ... )
    True
    """


@dataclass(frozen=True)
class AdapterPolicy:
    """Configure solver-independent result acceptance rules.

    Parameters
    ----------
    require_convergence : bool, default=True
        Reject a result when any solve reports non-convergence.
    maximum_relative_residual : float or None, optional
        Reject a result whose worst reported residual exceeds this value.
        None delegates residual acceptance entirely to the backend.
    require_all_valid : bool, default=True
        Reject results containing masked or non-finite observations.
    emit_capability_warnings : bool, default=True
        Emit advisory messages from capability assessment.
    wrap_backend_exceptions : bool, default=True
        Wrap ordinary backend exceptions in :class:`BackendExecutionError`.

    Examples
    --------
    >>> policy = AdapterPolicy(maximum_relative_residual=1e-6)
    >>> policy.maximum_relative_residual
    1e-06
    """

    require_convergence: bool = True
    maximum_relative_residual: float | None = None
    require_all_valid: bool = True
    emit_capability_warnings: bool = True
    wrap_backend_exceptions: bool = True

    def __post_init__(self) -> None:
        threshold = self.maximum_relative_residual
        if threshold is not None:
            threshold = float(threshold)
            if not np.isfinite(threshold) or threshold < 0:
                raise ValueError(
                    "maximum_relative_residual must be finite and non-negative."
                )
        object.__setattr__(self, "maximum_relative_residual", threshold)


class BaseMaxwellAdapter(ABC):
    """Base class enforcing common preflight and postflight validation.

    Parameters
    ----------
    capabilities : BackendCapabilities
        Immutable declaration for the exact backend version.
    policy : AdapterPolicy or None, optional
        Result acceptance policy. Defaults to :class:`AdapterPolicy`.

    Notes
    -----
    Implementations override only :meth:`_solve_backend`. They must return a
    canonical :class:`ForwardResult`; all validation is performed by
    :meth:`solve`.

    Examples
    --------
    See :class:`CallableMaxwellAdapter` for a minimal concrete adapter.
    """

    def __init__(
        self,
        capabilities: BackendCapabilities,
        policy: AdapterPolicy | None = None,
    ) -> None:
        if not isinstance(capabilities, BackendCapabilities):
            raise TypeError("capabilities must be a BackendCapabilities object.")
        if policy is not None and not isinstance(policy, AdapterPolicy):
            raise TypeError("policy must be an AdapterPolicy or None.")
        self._capabilities = capabilities
        self._policy = AdapterPolicy() if policy is None else policy

    @property
    def capabilities(self) -> BackendCapabilities:
        """Return the immutable backend capability declaration.

        Returns
        -------
        BackendCapabilities
            Physical and numerical scope of this adapter.

        Examples
        --------
        >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
        >>> CallableMaxwellAdapter(
        ...     cap, lambda problem: None
        ... ).capabilities is cap
        True
        """
        return self._capabilities

    @property
    def policy(self) -> AdapterPolicy:
        """Return the immutable result acceptance policy.

        Returns
        -------
        AdapterPolicy
            Policy applied after every solve.

        Examples
        --------
        >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
        >>> CallableMaxwellAdapter(
        ...     cap, lambda problem: None
        ... ).policy.require_convergence
        True
        """
        return self._policy

    def assess(self, problem: MaxwellProblem) -> CompatibilityReport:
        """Assess a problem against declared backend capabilities.

        Parameters
        ----------
        problem : MaxwellProblem
            Candidate simulation problem.

        Returns
        -------
        CompatibilityReport
            Consolidated errors and warnings without invoking the solver.

        Examples
        --------
        >>> from .contracts import MaxwellMesh, ReceiverSet
        >>> mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
        >>> problem = MaxwellProblem(
        ...     mesh,
        ...     np.ones((2, 2)),
        ...     [1],
        ...     ReceiverSet([[0, 0]], ["S"]),
        ...     ("zxy",),
        ... )
        >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
        >>> CallableMaxwellAdapter(cap, lambda problem: None).assess(
        ...     problem
        ... ).compatible
        True
        """
        return self.capabilities.assess(problem)

    def solve(self, problem: MaxwellProblem) -> ForwardResult:
        """Validate, execute, and verify one Maxwell problem.

        Parameters
        ----------
        problem : MaxwellProblem
            Solver-neutral simulation input.

        Returns
        -------
        ForwardResult
            Canonical, problem-matched impedance result.

        Raises
        ------
        IncompatibleProblemError
            If preflight capability checks fail.
        BackendExecutionError
            If the numerical backend raises an ordinary exception.
        InvalidBackendResultError
            If returned output violates the result contract or problem axes.
        SolverConvergenceError
            If convergence, residual, or validity policy fails.

        Examples
        --------
        Concrete execution examples require a backend callback; see
        :class:`CallableMaxwellAdapter`.
        """
        report = self.assess(problem)
        if report.errors:
            raise IncompatibleProblemError(
                f"backend {self.capabilities.name!r} is incompatible: "
                + "; ".join(report.errors)
            )
        if self.policy.emit_capability_warnings:
            for message in report.warnings:
                warnings.warn(
                    f"backend {self.capabilities.name!r}: {message}",
                    RuntimeWarning,
                    stacklevel=2,
                )
        try:
            result = self._solve_backend(problem)
        except MaxwellAdapterError:
            raise
        except Exception as exc:
            if not self.policy.wrap_backend_exceptions:
                raise
            raise BackendExecutionError(
                f"backend {self.capabilities.name!r} failed: {exc}"
            ) from exc
        self._validate_result(problem, result)
        return result

    def solve_many(
        self, problems: Iterable[MaxwellProblem]
    ) -> tuple[ForwardResult, ...]:
        """Solve problems sequentially while preserving input order.

        Parameters
        ----------
        problems : iterable of MaxwellProblem
            Finite problem stream. Execution stops at the first failure.

        Returns
        -------
        tuple of ForwardResult
            Results in exactly the supplied order.

        Examples
        --------
        An empty collection performs no backend calls:

        >>> cap = BackendCapabilities("demo", "1", (2,), ("zxy",))
        >>> CallableMaxwellAdapter(cap, lambda problem: None).solve_many([])
        ()
        """
        return tuple(self.solve(problem) for problem in problems)

    @abstractmethod
    def _solve_backend(self, problem: MaxwellProblem) -> ForwardResult:
        raise NotImplementedError

    def _validate_result(self, problem: MaxwellProblem, result: ForwardResult) -> None:
        if not isinstance(result, ForwardResult):
            raise InvalidBackendResultError(
                "backend must return a ForwardResult; "
                f"received {type(result).__name__}."
            )
        try:
            result.validate_against(problem)
        except ValueError as exc:
            raise InvalidBackendResultError(str(exc)) from exc
        if (
            result.backend_name != self.capabilities.name
            or result.backend_version != self.capabilities.version
        ):
            raise InvalidBackendResultError(
                "result backend identity differs from adapter capabilities."
            )
        diagnostics = result.diagnostics
        if self.policy.require_convergence and not diagnostics.success:
            failed = int(
                np.size(diagnostics.converged) - np.count_nonzero(diagnostics.converged)
            )
            raise SolverConvergenceError(
                f"backend {self.capabilities.name!r} reported {failed} unconverged solve(s)."
            )
        threshold = self.policy.maximum_relative_residual
        if threshold is not None and diagnostics.maximum_relative_residual > threshold:
            raise SolverConvergenceError(
                "maximum relative residual "
                f"{diagnostics.maximum_relative_residual:.6g} exceeds {threshold:.6g}."
            )
        if self.policy.require_all_valid and not np.all(result.valid):
            invalid = int(result.valid.size - np.count_nonzero(result.valid))
            raise SolverConvergenceError(
                f"backend result contains {invalid} invalid impedance value(s)."
            )


class CallableMaxwellAdapter(BaseMaxwellAdapter):
    """Adapt a trusted callable to the validated Maxwell backend interface.

    Parameters
    ----------
    capabilities : BackendCapabilities
        Static declaration matching callback output identity.
    solver : callable
        Function accepting :class:`MaxwellProblem` and returning
        :class:`ForwardResult`.
    policy : AdapterPolicy or None, optional
        Solver-independent acceptance policy.

    Examples
    --------
    >>> from .contracts import MaxwellMesh, ReceiverSet, SolverDiagnostics
    >>> mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
    >>> problem = MaxwellProblem(
    ...     mesh, np.ones((2, 2)), [1], ReceiverSet([[0, 0]], ["S"]), ("zxy",)
    ... )
    >>> cap = BackendCapabilities(
    ...     "demo", "1", (2,), ("zxy",), verified_benchmarks=("half-space",)
    ... )
    >>> def solver(value):
    ...     diagnostics = SolverDiagnostics([[True]], [[0]], [[0]], 0)
    ...     return ForwardResult(
    ...         value.problem_hash,
    ...         value.frequencies_hz,
    ...         value.receivers.names,
    ...         value.components,
    ...         [[[1j]]],
    ...         None,
    ...         "demo",
    ...         "1",
    ...         diagnostics,
    ...     )
    >>> CallableMaxwellAdapter(cap, solver).solve(problem).success
    True
    """

    def __init__(
        self,
        capabilities: BackendCapabilities,
        solver: Callable[[MaxwellProblem], ForwardResult],
        policy: AdapterPolicy | None = None,
    ) -> None:
        if not callable(solver):
            raise TypeError("solver must be callable.")
        super().__init__(capabilities, policy)
        self._solver = solver

    def _solve_backend(self, problem: MaxwellProblem) -> ForwardResult:
        return self._solver(problem)
