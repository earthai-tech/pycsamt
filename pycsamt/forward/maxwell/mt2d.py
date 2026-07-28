# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Validated adapter for the in-repo 2-D MT finite-difference solver.

:class:`MT2DAdapter` bridges the solver-neutral
:class:`~pycsamt.forward.maxwell.contracts.MaxwellProblem` /
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult` contract to the
existing, tested finite-difference implementation in
:mod:`pycsamt.forward.em2d`
(:class:`~pycsamt.forward.em2d.MT2DForward`) over a
:class:`~pycsamt.forward.grid2d.Grid2D` model. No new physics is
introduced here: this module only translates problem/result shapes and
enforces the assumptions the wrapped solver actually makes.

Known solver-specific limits, declared through
:class:`~pycsamt.forward.maxwell.backends.BackendCapabilities` or checked
in :meth:`MT2DAdapter.assess`:

* only the ``zxy`` (TE) and ``zyx`` (TM) impedance components are
  produced, matching the physical 2-D TE/TM decomposition;
* every receiver must sit exactly at the surface (``z = 0``); the wrapped
  solver has no borehole or buried-receiver formulation;
* the earth model has no inactive (air) cells: the solver always treats
  the full mesh as earth, so a problem must use
  ``mark_air_inactive=False`` when built from
  :meth:`~pycsamt.forward.maxwell.mesh.SolverMeshModel.to_problem`;
* magnetic permeability is fixed at the vacuum value
  ``4*pi*1e-7`` H/m, matching ``pycsamt.forward.em2d.MU0``;
* only the ``exp(+iwt)`` phasor convention is supported, matching the
  sign convention hardcoded in :mod:`pycsamt.forward.em2d`.

A numerically degenerate solve (near-zero surface magnetic field at a
station) makes :mod:`pycsamt.forward.em2d` fall back to ``0 + 0j``
rather than raise or return ``NaN``; such an entry is currently
reported as valid because
:attr:`~pycsamt.forward.maxwell.contracts.ForwardResult.valid` defaults
to finiteness. This is a known, narrow limitation inherited from the
wrapped solver, not something this adapter can detect from its output
alone.

Scipy (a required pycsamt dependency, not optional) is imported only
through :mod:`pycsamt.forward.em2d`; this module performs no lazy or
optional-dependency handling of its own.
"""

from __future__ import annotations

import time

import numpy as np

from ..em2d import ForwardResponse2D, MT2DForward
from ..grid2d import Grid2D
from .adapters import AdapterPolicy, BaseMaxwellAdapter
from .backends import (
    BackendCapabilities,
    BackendRegistration,
    CompatibilityReport,
)
from .backends import register_backend as _register_backend
from .contracts import ForwardResult, MaxwellProblem, SolverDiagnostics

__all__ = [
    "MT2DAdapter",
    "register_mt2d_backend",
]

_MU0 = 4.0e-7 * np.pi
_SURFACE_TOLERANCE_M = 1e-6
_PERMEABILITY_RTOL = 1e-9

# Benchmarks this adapter is currently known to pass within the default
# pycsamt.forward.maxwell.benchmarks.BenchmarkThresholds(). Verified by
# pycsamt/forward/tests/test_maxwell_mt2d.py; update together with that
# test file, never in isolation.
_VERIFIED_BENCHMARKS = ("half-space", "layered-earth")


class MT2DAdapter(BaseMaxwellAdapter):
    """Validated 2-D MT adapter over the in-repo finite-difference solver.

    Parameters
    ----------
    version : str, default="1.0"
        Adapter version reported in every
        :class:`~pycsamt.forward.maxwell.contracts.ForwardResult`. Bump
        this when the translation logic in this module changes in a way
        that could alter numerical output.
    policy : AdapterPolicy or None, optional
        Solver-independent result acceptance policy. Defaults to
        :class:`~pycsamt.forward.maxwell.adapters.AdapterPolicy`.
    verbose : bool, default=False
        Forwarded to :class:`~pycsamt.forward.em2d.MT2DForward`; prints
        per-frequency progress when true.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.maxwell import (
    ...     MaxwellMesh,
    ...     MaxwellProblem,
    ...     ReceiverSet,
    ... )
    >>> mesh = MaxwellMesh(
    ...     np.linspace(0, 10_000, 41), np.linspace(0, 5_000, 31)
    ... )
    >>> problem = MaxwellProblem(
    ...     mesh,
    ...     np.full(mesh.shape, 1.0 / 100.0),
    ...     [10.0, 1.0],
    ...     ReceiverSet([[5_000.0, 0.0]], ["S00"]),
    ...     ("zxy", "zyx"),
    ... )
    >>> adapter = MT2DAdapter(verbose=False)
    >>> result = adapter.solve(problem)
    >>> result.shape
    (1, 2, 2)
    """

    def __init__(
        self,
        *,
        version: str = "1.0",
        policy: AdapterPolicy | None = None,
        verbose: bool = False,
    ) -> None:
        capabilities = BackendCapabilities(
            name="mt2d",
            version=version,
            dimensions=(2,),
            components=("zxy", "zyx"),
            time_conventions=("exp(+iwt)",),
            supports_nonuniform_mesh=True,
            supports_inactive_cells=False,
            supports_topography=False,
            supports_anisotropy=False,
            verified_benchmarks=_VERIFIED_BENCHMARKS,
        )
        super().__init__(capabilities, policy)
        self._verbose = bool(verbose)

    def assess(self, problem: MaxwellProblem) -> CompatibilityReport:
        """Assess a problem, adding this solver's surface-only checks.

        Parameters
        ----------
        problem : MaxwellProblem
            Candidate simulation problem.

        Returns
        -------
        CompatibilityReport
            The generic capability report from
            :meth:`~pycsamt.forward.maxwell.adapters.BaseMaxwellAdapter.assess`,
            extended with two solver-specific checks: every receiver must
            be at the surface, and the permeability must be the vacuum
            value the wrapped solver hardcodes.

        Examples
        --------
        See :class:`MT2DAdapter` for a complete solve example; a problem
        with a buried receiver is rejected before the solver runs.
        """
        base = super().assess(problem)
        errors = list(base.errors)
        depth = problem.receivers.coordinates_m[:, 1]
        if np.any(np.abs(depth) > _SURFACE_TOLERANCE_M):
            errors.append("mt2d only evaluates receivers at the surface (z=0)")
        x = problem.receivers.coordinates_m[:, 0]
        x_lo, x_hi = problem.mesh.x_edges_m[0], problem.mesh.x_edges_m[-1]
        if np.any((x < x_lo) | (x > x_hi)):
            errors.append(
                "receiver x coordinates must lie within the mesh "
                f"[{x_lo:g}, {x_hi:g}] m"
            )
        if not np.isclose(
            problem.magnetic_permeability_h_m,
            _MU0,
            rtol=_PERMEABILITY_RTOL,
        ):
            errors.append(
                f"mt2d assumes vacuum magnetic permeability ({_MU0:.6g} H/m)"
            )
        if len(errors) == len(base.errors):
            return base
        return CompatibilityReport(
            self.capabilities.name, False, tuple(errors), base.warnings
        )

    def _solve_backend(self, problem: MaxwellProblem) -> ForwardResult:
        grid = _grid_from_problem(problem)
        solver = MT2DForward(
            problem.frequencies_hz, grid, verbose=self._verbose
        )
        start = time.monotonic()
        response = solver.run()
        runtime_s = time.monotonic() - start
        return _result_from_response(
            problem, response, self.capabilities, runtime_s
        )


def register_mt2d_backend(*, replace: bool = False) -> None:
    """Register :class:`MT2DAdapter` in the process-wide backend registry.

    Parameters
    ----------
    replace : bool, default=False
        Explicitly replace an existing ``"mt2d"`` registration.

    Examples
    --------
    >>> from pycsamt.forward.maxwell import create_backend, list_backends
    >>> register_mt2d_backend(replace=True)
    >>> "mt2d" in list_backends()
    True
    >>> create_backend("mt2d").capabilities.name
    'mt2d'
    """
    capabilities = BackendCapabilities(
        name="mt2d",
        version="1.0",
        dimensions=(2,),
        components=("zxy", "zyx"),
        time_conventions=("exp(+iwt)",),
        supports_nonuniform_mesh=True,
        supports_inactive_cells=False,
        supports_topography=False,
        supports_anisotropy=False,
        verified_benchmarks=_VERIFIED_BENCHMARKS,
    )
    _register_backend(
        BackendRegistration(capabilities, MT2DAdapter), replace=replace
    )


def _grid_from_problem(problem: MaxwellProblem) -> Grid2D:
    mesh = problem.mesh
    resistivity = 1.0 / problem.conductivity_s_m
    return Grid2D(
        dx=mesh.cell_widths_m["x"],
        dz=mesh.cell_widths_m["z"],
        resistivity=resistivity,
        x_stations=problem.receivers.coordinates_m[:, 0],
        n_pad=0,
        name=problem.metadata.get("benchmark", "mt2d"),
    )


def _result_from_response(
    problem: MaxwellProblem,
    response: ForwardResponse2D,
    capabilities: BackendCapabilities,
    runtime_s: float,
) -> ForwardResult:
    n_station = response.n_stations
    n_freq = response.n_freqs
    n_component = len(problem.components)
    impedance = np.empty((n_station, n_freq, n_component), dtype=complex)
    for column, name in enumerate(problem.components):
        source = response.zxy if name == "zxy" else response.zyx
        impedance[:, :, column] = source.T

    converged = np.stack(
        [
            np.isfinite(response.residual_te),
            np.isfinite(response.residual_tm),
        ],
        axis=1,
    )
    residual = np.stack([response.residual_te, response.residual_tm], axis=1)
    diagnostics = SolverDiagnostics(
        converged,
        np.zeros(residual.shape, dtype=int),
        np.where(np.isfinite(residual), residual, 0.0),
        runtime_s,
    )
    return ForwardResult(
        problem.problem_hash,
        problem.frequencies_hz,
        problem.receivers.names,
        problem.components,
        impedance,
        None,
        capabilities.name,
        capabilities.version,
        diagnostics,
    )
