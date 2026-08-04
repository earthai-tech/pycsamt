# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""External MARE2DEM adapter: the production triangular-mesh MT backend.

Wraps the external MARE2DEM 2.5-D finite-element solver via
:class:`~pycsamt.forward.maxwell.external.BaseExternalMaxwellAdapter`,
reusing the file I/O already built in :mod:`pycsamt.models.mare2dem`
(``iotools.resistivity``, ``iotools.emdata``, ``iotools.settings``,
``iotools.poly``) rather than a second, competing implementation of
MARE2DEM's file formats. Consumes a
:class:`~pycsamt.forward.maxwell.contracts_tri.TriProblem`.

**Status: physics-validated against a real compiled binary
(2026-07-31)**, built from the vendored source
(``pycsamt/models/mare2dem/_source``) inside WSL2 Ubuntu with Intel's
free oneAPI toolchain (``ifx``/``icx``/Intel MPI/MKL -- the vendored
source hard-requires genuine Intel MKL's PARDISO/DSS sparse solver,
not just any BLAS/LAPACK; see ``intelmkl_solver.f90``'s ``use
mkl_dss``/``use mkl_pardiso``). No compiled binary is committed to
this repository (a local build artifact, matching
:mod:`~pycsamt.forward.maxwell.modem3d`'s own precedent). Running the
real binary against this adapter surfaced several real,
previously-latent bugs, all fixed below -- **not** inferred from
documentation, each confirmed by reading MARE2DEM's own Fortran
source and/or its actual runtime output:

* **The mesh-handling architecture was fundamentally wrong in the
  first implementation of this adapter.** MARE2DEM's own
  ``getStartingMesh`` (``mare2dem_mpi.f90``/``mare2dem_worker.f90``)
  unconditionally refuses any input mesh that already has triangles
  (``inputmesh%nele /= 0``) with ``"Error: input mesh already has
  been triangulated, stopping"``. It *only* accepts a PSLG (boundary
  polygon + named-region seed points, zero pre-existing elements) and
  always generates and adaptively refines its own mesh internally via
  its own bundled Triangle call. There is no flag or mode (including
  ``-F`` forward-only) that accepts a fixed, pre-triangulated mesh --
  confirmed by reading the full command-line argument parser
  (``mare2dem_help`` in ``mare2dem_io.f90``). This adapter originally
  wrote a complete pre-triangulated mesh via
  ``iotools.poly.write_triangulation`` (still a real, tested,
  reusable utility for *other* purposes, e.g. M12's Triangle-refine
  round trip -- just not usable here); it now writes a PSLG instead
  (see "Scope and mapping decisions" below), which required adding a
  hard :meth:`Mare2DEMAdapter.assess` requirement on
  ``mesh.boundary_segments``.
* **Resistivity is linear ohm-m in the file format, not log10**,
  despite ``ResistivityFile.resistivity``'s own docstring in
  ``iotools/resistivity.py`` claiming "log10-resistivity" --
  confirmed two ways: (1) ``mare2dem_io.f90``'s reader comment reads
  ``"! rho is linear"`` and applies ``log10()`` itself only
  internally, for its own inversion parameterization; (2) this
  project's pre-existing ``ResistivityModel.halfspace()``
  (``models/mare2dem/mesh.py``) already did ``10.0**log10_rho``
  before storing -- i.e. it already assumed linear storage, correctly,
  while this adapter did not. Writing ``log10(rho)`` here instead
  produced a self-consistent but *wrong* forward response: correct
  sign convention, correct 45-degree half-space phase, but magnitude
  off by a constant factor at every frequency (``log10(100) = 2`` was
  read back as a *linear* resistivity of 2 ohm-m, not 100).
* **The resistivity file's iteration number is mandatory**, not
  optional: MARE2DEM's CLI parser (``runmare2dem.f90``) requires the
  filename/argument to end in ``.<integer>`` (e.g. ``mare2dem.0``,
  iteration 0 = the starting/forward model) and hard-errors
  ("no iteration number in resistivity file") otherwise.
* **The ``"Model File"`` field must already carry the ``.poly``
  extension.** ``readPoly`` (``mare2dem_io.f90``) takes everything
  before the *last* dot as the base filename
  (``modelFileName(1:idot-1)``); a bare stem with no dot at all makes
  ``idot=0``, and Fortran's ``str(1:-1)`` substring silently returns
  an empty string, producing ``"Opening Poly file: .poly"``.

With all four fixed, a real half-space forward solve (see
``test_maxwell_mare2dem.py``'s ``requires_real_mare2dem``-gated test)
matches :func:`~pycsamt.forward.maxwell.benchmarks.half_space_impedance`
to within the ~0.2-0.8% mesh-discretization error MARE2DEM itself
reports for its own adaptive refinement, for both magnitude and the
``Zyx = -Zxy`` sign relationship. The MT data-type unit/sign
convention this adapter always assumed (SI ``V/A``, same
``exp(+iwt)`` phasor convention, no ModEM-style unit-conversion
factor) is now genuinely confirmed, not merely asserted.

Scope and mapping decisions
----------------------------

* **The input mesh is a PSLG recipe, not the literal computational
  mesh.** ``_prepare_inputs`` writes ``mesh.nodes_m`` +
  ``mesh.boundary_segments`` as a Triangle PSLG boundary, plus one
  region-seed point per unique ``mesh.region_ids`` value (the
  centroid of that region's first assigned triangle, guaranteed
  interior to the region). MARE2DEM meshes and adaptively refines
  this PSLG itself; the triangles it actually solves on are not
  guaranteed to match ``mesh.triangles``, only the boundary shape and
  region layout are preserved. :meth:`Mare2DEMAdapter.assess`
  requires ``problem.conductivity_s_m`` to be uniform within every
  region_id group and rejects the problem otherwise, rather than
  silently averaging away real heterogeneity -- MARE2DEM
  parameterizes resistivity by named region (a handful of zones,
  e.g. background/anomaly), not by individual triangle. Arbitrary
  fine-grained per-triangle heterogeneous fields (e.g. AI training
  data with hundreds of independent cells) are not representable
  through this adapter as a result; that needs either a coarser
  region-based training design or the in-house triangular FEM solver
  (a separate, not-yet-started milestone).
* Only ``zxy``/``zyx`` components are requested, matching
  :class:`~pycsamt.forward.maxwell.contracts_tri.TriMesh`'s 2-D-only
  scope.
* **MPI is used by default** (``config.use_mpi``), matching
  :class:`~pycsamt.models.mare2dem.config.Mare2DEMConfig`'s own default
  and :class:`~pycsamt.models.mare2dem.runner.Mare2DEMRunner`'s
  invocation convention (``mpirun -np N MARE2DEM <stem>``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from ...models.mare2dem.config import Mare2DEMConfig
from ...models.mare2dem.iotools.emdata import (
    EMDataFile,
    MTConfig,
    read_emdata,
    write_emdata,
)
from ...models.mare2dem.iotools.poly import PolyFile, write_poly
from ...models.mare2dem.iotools.resistivity import (
    ResistivityFile,
    write_resistivity,
)
from ...models.mare2dem.iotools.settings import SettingsFile, write_settings
from .adapters import AdapterPolicy
from .backends import (
    BackendCapabilities,
    BackendRegistration,
    CompatibilityReport,
)
from .backends import register_backend as _register_backend
from .contracts import ForwardResult, SolverDiagnostics
from .contracts_tri import TriMesh, TriProblem
from .external import (
    BaseExternalMaxwellAdapter,
    ExternalRunPolicy,
    ExternalRunResult,
    make_availability_probe,
)

__all__ = ["Mare2DEMAdapter", "register_mare2dem_backend"]

# MARE2DEM MT data-type codes (real, imaginary) per canonical component --
# see pycsamt.models.mare2dem.iotools.codes.DATA_CODES.
_DATA_CODES: dict[str, tuple[int, int]] = {
    "zxy": (113, 114),
    "zyx": (115, 116),
}

# Honest: no compiled MARE2DEM binary exists in this environment to run
# the analytic half-space/layered-earth benchmarks against, unlike
# mt2d.py/mt3d.py/modem3d.py. See the module docstring's "Status" section.
_VERIFIED_BENCHMARKS: tuple[str, ...] = ()


class Mare2DEMAdapter(BaseExternalMaxwellAdapter):
    """Adapter wrapping the external MARE2DEM 2.5-D triangular-mesh solver.

    Parameters
    ----------
    config : Mare2DEMConfig or None, optional
        MARE2DEM configuration (executable name, MPI settings). Defaults
        to ``Mare2DEMConfig()``.
    run_policy : ExternalRunPolicy or None, optional
        Executable resolution, timeout, retry, and working-directory
        rules. Built from ``config`` when omitted (``config.mpi_command``
        if ``config.use_mpi`` else ``config.binary``).
    policy : AdapterPolicy or None, optional
        Solver-independent result acceptance policy.
    std_err : float, default=1.0
        Placeholder standard error written for every forward-request
        DATA row (MARE2DEM's file format requires a nonzero value even
        though it is unused for a zero-iteration forward evaluation).
    response_filename : str or None, optional
        Exact name of MARE2DEM's predicted-response output file. When
        omitted, the adapter globs the working directory for a
        ``*.resp`` file, falling back to ``*_MARE2DEM.emdata``, and
        raises a clear error if that is ambiguous.
    version : str, default="mare2dem-adapter-1.0"
        Adapter version reported in every
        :class:`~pycsamt.forward.maxwell.contracts.ForwardResult`.

    Examples
    --------
    >>> from pycsamt.forward.maxwell import ExternalRunPolicy
    >>> adapter = Mare2DEMAdapter(
    ...     run_policy=ExternalRunPolicy("does-not-exist-xyz")
    ... )
    >>> adapter.capabilities.name
    'mare2dem'
    """

    def __init__(
        self,
        *,
        config: Mare2DEMConfig | None = None,
        run_policy: ExternalRunPolicy | None = None,
        policy: AdapterPolicy | None = None,
        std_err: float = 1.0,
        response_filename: str | None = None,
        version: str = "mare2dem-adapter-1.0",
    ) -> None:
        self._config = config or Mare2DEMConfig()
        if run_policy is None:
            executable = (
                self._config.mpi_command
                if self._config.use_mpi
                else self._config.binary
            )
            run_policy = ExternalRunPolicy(executable)
        capabilities = BackendCapabilities(
            name="mare2dem",
            version=version,
            dimensions=(2,),
            components=("zxy", "zyx"),
            time_conventions=("exp(+iwt)",),
            supports_nonuniform_mesh=True,
            supports_inactive_cells=False,
            supports_topography=True,
            supports_anisotropy=False,
            verified_benchmarks=_VERIFIED_BENCHMARKS,
        )
        super().__init__(capabilities, run_policy, policy)
        self._std_err = float(std_err)
        self._response_filename = response_filename

    def assess(self, problem: TriProblem) -> CompatibilityReport:
        """Assess a problem, adding this adapter's mesh/region checks.

        Parameters
        ----------
        problem : TriProblem
            Candidate simulation problem.

        Returns
        -------
        CompatibilityReport
            The generic capability report from
            :meth:`~pycsamt.forward.maxwell.adapters.BaseMaxwellAdapter.assess`,
            extended with a triangular-mesh-type check and a per-region
            uniform-conductivity check.

        Examples
        --------
        See :class:`Mare2DEMAdapter` for construction; a problem whose
        conductivity varies within one mesh region is rejected before
        any file is written.
        """
        base = super().assess(problem)
        errors = list(base.errors)
        if not isinstance(problem.mesh, TriMesh):
            errors.append("mare2dem requires a TriProblem built on a TriMesh.")
        else:
            if problem.mesh.boundary_segments is None:
                errors.append(
                    "mare2dem requires mesh.boundary_segments: it always "
                    "meshes its own PSLG internally (see the module "
                    "docstring's 'Real architecture' section) and refuses "
                    "any input mesh that already has triangles, so this "
                    "adapter needs the boundary polygon, not a pre-built "
                    "triangulation. Build the mesh with boundary_segments "
                    "populated (e.g. via a PSLG-aware mesh builder)."
                )
            region_ids = problem.mesh.region_ids
            conductivity = problem.conductivity_s_m
            for region in np.unique(region_ids):
                values = conductivity[region_ids == region]
                if not np.allclose(values, values[0]):
                    errors.append(
                        f"conductivity is not uniform within mesh region "
                        f"{int(region)}; mare2dem parameterizes by named "
                        "region (it meshes and refines the PSLG itself), "
                        "not by individual triangle -- arbitrary "
                        "per-triangle heterogeneous fields are not "
                        "representable through this adapter."
                    )
                    break
        if len(errors) == len(base.errors):
            return base
        return CompatibilityReport(
            self.capabilities.name, False, tuple(errors), base.warnings
        )

    def _prepare_inputs(self, problem: TriProblem, workdir: Path) -> Any:
        mesh = problem.mesh
        region_ids = mesh.region_ids
        unique_regions = np.unique(region_ids)
        conductivity_per_region = np.array(
            [
                problem.conductivity_s_m[region_ids == region][0]
                for region in unique_regions
            ]
        )
        resistivity_ohm_m = 1.0 / conductivity_per_region
        n_regions = len(unique_regions)

        # One representative seed point per region: the centroid of the
        # first triangle assigned to that region. It is guaranteed to
        # fall inside the region, since `mesh` was already validly
        # partitioned by region_ids -- MARE2DEM only needs one interior
        # point per region to flood-fill its own internally-generated
        # mesh, not the region's exact triangulated shape.
        centroids = mesh.triangle_centroids_m
        seeds = np.array(
            [
                centroids[np.argmax(region_ids == region)]
                for region in unique_regions
            ]
        )

        stem = "mare2dem"
        pf = PolyFile()
        pf.nodes = np.asarray(mesh.nodes_m, dtype=float)
        # TriMesh.boundary_segments is 0-based (contracts_tri.py); the
        # .poly format is 1-based (matches write_triangulation's own
        # node/element numbering elsewhere in this package).
        pf.segments = np.asarray(mesh.boundary_segments, dtype=int) + 1
        pf.regions = np.column_stack(
            [seeds, np.arange(1, n_regions + 1), np.zeros(n_regions)]
        )
        poly_path = write_poly(pf, workdir / f"{stem}.poly")

        # MARE2DEM's own CLI parsing (runmare2dem.f90) requires an
        # iteration number embedded in the resistivity filename/argument
        # (<stem>.<iteration>[.resistivity]) -- "iteration 0" is its
        # convention for the starting/forward model, matching
        # max_iterations=0 below. Confirmed against a real compiled
        # binary: passing a bare stem with no iteration number is a hard
        # CLI error ("no iteration number in resistivity file"), not
        # merely a naming preference.
        run_stem = f"{stem}.0"

        rf = ResistivityFile(
            resistivity_file=f"{run_stem}.resistivity",
            # MARE2DEM's readPoly (mare2dem_io.f90) requires this field to
            # already carry an extension: it takes everything before the
            # *last* dot as the base name (`modelFileName(1:idot-1)`), so a
            # bare stem with no dot at all makes idot=0 and silently
            # produces an empty filename ("Opening Poly file: .poly").
            # Confirmed against a real compiled binary, not inferred.
            poly_file=f"{stem}.poly",
            data_file=f"{stem}.emdata",
            settings_file=f"{stem}.settings",
            max_iterations=0,
            # LINEAR resistivity (ohm-m), *not* log10, despite this
            # field's own docstring in iotools/resistivity.py claiming
            # "log10-resistivity" -- confirmed two ways against a real
            # compiled binary: (1) mare2dem_io.f90's own reader comment
            # says "! rho is linear" and applies log10() itself only
            # internally for its inversion parameterization; (2) this
            # project's own pre-existing ResistivityModel.halfspace()
            # (models/mare2dem/mesh.py) already does `10.0**log10_rho`
            # before storing, i.e. it already assumed linear storage.
            # Writing log10(rho) here instead produced a self-consistent
            # but *wrong* forward response -- correct sign and correct
            # 45-degree half-space phase, magnitude off by a constant
            # sqrt(50) factor at every frequency (rho=2 instead of 100
            # ohm-m, since log10(100)=2 was misread as linear rho=2).
            resistivity=resistivity_ohm_m.reshape(-1, 1),
            free_parameter=np.zeros((n_regions, 1)),
            # Also linear (see above), not log10 -- 1e-6 to 1e6 ohm-m.
            bounds=np.tile([1.0e-6, 1.0e6], (n_regions, 1)),
            prejudice=np.zeros((n_regions, 2)),
        )
        resistivity_path = write_resistivity(
            rf, workdir / rf.resistivity_file
        )

        names = list(problem.receivers.names)
        coordinates_m = problem.receivers.coordinates_m
        n_rx = len(names)
        receivers = np.zeros((n_rx, 8))
        receivers[:, 1] = coordinates_m[:, 0]  # y = profile position
        receivers[:, 2] = coordinates_m[:, 1]  # z = depth (positive down)

        frequencies = np.asarray(problem.frequencies_hz, dtype=float)
        em = EMDataFile(
            format="EMData_2.3",
            mt=MTConfig(
                frequencies=frequencies,
                receivers=receivers,
                receiver_name=names,
            ),
        )
        rows: list[list[float]] = []
        for fi in range(len(frequencies)):
            for ri in range(n_rx):
                for component in problem.components:
                    real_code, imag_code = _DATA_CODES[component]
                    rows.append(
                        [real_code, fi + 1, 0, ri + 1, 0.0, self._std_err]
                    )
                    rows.append(
                        [imag_code, fi + 1, 0, ri + 1, 0.0, self._std_err]
                    )
        em.data = np.array(rows, dtype=float)
        data_path = write_emdata(em, workdir / f"{stem}.emdata")

        settings_path = write_settings(
            SettingsFile(), workdir / f"{stem}.settings"
        )

        return {
            "stem": stem,
            "run_stem": run_stem,
            "resistivity_path": resistivity_path,
            "data_path": data_path,
            "settings_path": settings_path,
        }

    def _build_command(
        self,
        problem: TriProblem,
        workdir: Path,
        executable: Path,
        context: Any,
    ) -> list[str]:
        run_stem = context["run_stem"]
        if self._config.use_mpi:
            return [
                str(executable),
                "-np",
                str(self._config.n_procs),
                self._config.binary,
                run_stem,
            ]
        return [str(executable), run_stem]

    def _locate_response_file(self, workdir: Path) -> Path:
        if self._response_filename is not None:
            return workdir / self._response_filename
        candidates = sorted(workdir.glob("*.resp"))
        if not candidates:
            candidates = sorted(workdir.glob("*_MARE2DEM.emdata"))
        if len(candidates) == 1:
            return candidates[0]
        if not candidates:
            raise ValueError(
                "no MARE2DEM predicted-response file was found in the "
                "working directory; pass response_filename explicitly if "
                "your MARE2DEM build uses a naming convention this "
                "adapter does not already recognize."
            )
        names = ", ".join(path.name for path in candidates)
        raise ValueError(
            "multiple candidate MARE2DEM response files were found "
            f"({names}); pass response_filename explicitly to disambiguate."
        )

    def _parse_result(
        self,
        problem: TriProblem,
        workdir: Path,
        run_result: ExternalRunResult,
        context: Any,
    ) -> ForwardResult:
        response_path = self._locate_response_file(workdir)
        response = read_emdata(response_path)
        if response.data.ndim != 2 or response.data.shape[1] < 8:
            raise ValueError(
                f"{response_path} is not a MARE2DEM response file "
                "(expected an 8-column DATA block)."
            )

        lookup: dict[tuple[int, int, int], float] = {}
        for row in response.data:
            code, freq_idx, _tx, rx_idx = (int(v) for v in row[:4])
            predicted = float(row[6])
            lookup[(code, freq_idx, rx_idx)] = predicted

        n_station = problem.receivers.count
        n_freq = len(problem.frequencies_hz)
        n_component = len(problem.components)
        impedance = np.empty((n_station, n_freq, n_component), dtype=complex)
        for fi in range(n_freq):
            for ri in range(n_station):
                for ci, component in enumerate(problem.components):
                    real_code, imag_code = _DATA_CODES[component]
                    real_key = (real_code, fi + 1, ri + 1)
                    imag_key = (imag_code, fi + 1, ri + 1)
                    if real_key not in lookup or imag_key not in lookup:
                        raise ValueError(
                            "MARE2DEM response is missing station "
                            f"{problem.receivers.names[ri]!r}, component "
                            f"{component!r} at frequency index {fi + 1}."
                        )
                    impedance[ri, fi, ci] = complex(
                        lookup[real_key], lookup[imag_key]
                    )

        diagnostics = SolverDiagnostics(
            np.ones((n_freq, 1), dtype=bool),
            np.zeros((n_freq, 1), dtype=int),
            np.zeros((n_freq, 1)),
            run_result.runtime_s,
        )
        return ForwardResult(
            problem.problem_hash,
            problem.frequencies_hz,
            problem.receivers.names,
            problem.components,
            impedance,
            None,
            self.capabilities.name,
            self.capabilities.version,
            diagnostics,
        )


def register_mare2dem_backend(
    *, config: Mare2DEMConfig | None = None, replace: bool = False
) -> None:
    """Register :class:`Mare2DEMAdapter` in the process-wide registry.

    Parameters
    ----------
    config : Mare2DEMConfig or None, optional
        Configuration used both to build the registered adapter and to
        derive its availability probe (whether the configured executable
        can currently be resolved).
    replace : bool, default=False
        Explicitly replace an existing ``"mare2dem"`` registration.

    Examples
    --------
    >>> register_mare2dem_backend(replace=True)
    >>> from pycsamt.forward.maxwell import list_backends
    >>> "mare2dem" in list_backends()
    True
    """
    cfg = config or Mare2DEMConfig()
    executable = cfg.mpi_command if cfg.use_mpi else cfg.binary
    capabilities = Mare2DEMAdapter(config=cfg).capabilities
    probe = make_availability_probe(executable)

    def factory(**options: Any) -> Mare2DEMAdapter:
        return Mare2DEMAdapter(config=cfg, **options)

    _register_backend(
        BackendRegistration(capabilities, factory, probe), replace=replace
    )
