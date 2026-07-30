# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""External ModEM adapter: the trusted production 3-D backend.

Per ``AI-INVERSION`` (Decision 1), ModEM — not an
in-house solver — is the recommended path for genuine 3-D production
forward modeling. :class:`ModEm3DAdapter` bridges the solver-neutral
:class:`~pycsamt.forward.maxwell.contracts.MaxwellProblem` /
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult` contract to
the vendored ModEM Fortran executable
(``pycsamt/models/modem/_source/3D/Mod3DMT.f90``, v6.2.6; Egbert,
Kelbert & Meqbel) via
:class:`~pycsamt.forward.maxwell.external.BaseExternalMaxwellAdapter`,
reusing the file I/O already built in :mod:`pycsamt.models.modem`
(:class:`~pycsamt.models.modem.data.ModEmData`,
:class:`~pycsamt.models.modem.model3d.ModEmModel3D`) rather than a
second, competing implementation of the ModEM file formats.

**Status: physics-validated against a real compiled binary (2026-07-29).**
No compiled ``Mod3DMT`` binary is committed to this repository (it is
a local build artifact, gitignored — see
``pycsamt/models/modem/_source/README`` for build instructions), but
one *was* built here with a MinGW-w64 gfortran/OpenBLAS toolchain and
used to run this project's own analytic benchmarks end-to-end. Doing
so surfaced and fixed several real, previously-latent bugs — this
adapter's generated files, and the vendored Fortran source itself,
had never actually been exercised against a live ModEM run before:

* the vendored ``Makefile`` was missing a build rule for
  ``sg_spherical.f90`` (`use`d directly by ``GridCalc.f90``) and never
  compiled ``Declaration_MPI``/``Sub_MPI``/``Main_MPI.f90`` (whose
  ``#ifdef MPI``-guarded bodies are unconditionally `use`d by several
  files regardless of MPI) nor passed ``-cpp`` to strip those guards
  for a serial build — all reproduce identically on Linux/Mac, not
  Windows-specific;
* :func:`~pycsamt.forward.maxwell.external.resolve_executable`'s
  ``search_paths`` fallback checked a literal
  ``Path(directory) / name`` rather than applying ``PATHEXT`` the way
  its own ``PATH`` lookup (via :func:`shutil.which`) already did, so a
  bare name like ``"Mod3DMT"`` never resolved to ``"Mod3DMT.exe"`` on
  Windows;
* :meth:`ModEm3DAdapter._build_command` omitted the third positional
  argument ModEM's own ``-F`` forward-mode requires (the predicted-data
  output filename), so it printed its own usage banner and exited 0
  having written nothing, indistinguishable from success until
  :meth:`ModEm3DAdapter._locate_predicted_file` found no output file;
* :class:`~pycsamt.models.modem.model3d.ModEmModel3D`'s WS-format
  reader/writer were missing a mandatory leading comment line ModEM's
  Fortran ``read_modelParam_ws`` (``WS.inc``) unconditionally reads
  and discards before the dimensions line — a real compiled binary
  immediately rejected the header-less file with a Fortran runtime
  error;
* ModEM's own ``setup_airlayers`` (``GridDef.f90``) hardcodes 10 air
  layers and its default "mirror" sizing method reads that many
  earth-layer widths with **no bounds check** against the actual
  earth cell count — fewer than 10 earth z-cells reads past the end
  of the array, producing garbage/NaN values that crash the solver
  with "b in QMR contains NaNs". :meth:`ModEm3DAdapter.assess` now
  rejects this before ever writing a file (see
  :data:`_MIN_EARTH_Z_CELLS`); this is a real latent bug in the
  vendored ModEM source itself, worked around here rather than patched
  there (patching unfamiliar, decades-old numerical Fortran to add a
  defensive bounds check was judged riskier than simply requiring
  enough earth cells).

With those fixed and at least :data:`_MIN_EARTH_Z_CELLS` earth
z-cells, this adapter passes both
:func:`~pycsamt.forward.maxwell.benchmarks.half_space_benchmark` and
:func:`~pycsamt.forward.maxwell.benchmarks.layered_earth_benchmark`
with real margin — see "Measured accuracy" below and
``pycsamt/forward/tests/test_maxwell_modem3d.py``'s
``requires_real_modem``-gated tests, which are skipped (not failed)
when no local binary is present, and were confirmed passing against
the real one built for this validation.

Scope and mapping decisions
----------------------------

* **No separate air layers** (``ModEmModel3D.n_air = 0``): the whole
  :class:`~pycsamt.forward.maxwell.contracts.MaxwellProblem` mesh is
  written as ModEM "earth" cells, matching the same "mesh top =
  physical surface" convention already used by
  :class:`~pycsamt.forward.maxwell.mt2d.MT2DAdapter` and
  :class:`~pycsamt.forward.maxwell.mt3d.MT3DAdapter` — not a ModEM
  limitation (ModEM supports real air/topography), a deliberate v1
  scope reduction for consistency and lower risk. Topography support
  is future work.
* **Non-uniform meshes are supported**: ModEM's own solver has no
  uniform-grid restriction, so ``supports_nonuniform_mesh=True`` (as
  of 2026-07-29, ``mt3d.py`` also supports non-uniform meshes, so this
  is no longer a point of difference between the two adapters).
* **Station and model coordinates are shifted** so the mesh's own
  minimum x/y edge maps to ModEM-local ``(0, 0)`` — the vendored
  :class:`~pycsamt.models.modem.data.ModEmData`/
  :class:`~pycsamt.models.modem.model3d.ModEmModel3D` pair has no
  shared origin/rotation field of its own (confirmed by reading both
  writers), so this adapter is the thing responsible for keeping
  station coordinates and model-grid coordinates in the same frame.
* **Units**: requests are written in ModEM's native
  ``[mV/km]/[nT]`` field-unit convention (the well-established
  default, unlike an unverified ``[V/A]`` request), and the response
  is converted to SI ``V/A`` on the way back using the standard
  ``4*pi*1e-4`` factor.
* **Full impedance tensor is always requested** from ModEM regardless
  of ``problem.components`` (no cost difference internally), then
  filtered down to the requested subset when building the result.
* Diagnostics are honest about what ModEM's plain predicted-data file
  does not expose: ``iterations`` is always 0 (forward mode does not
  iterate) and ``relative_residual`` is always 0.0 — not a measured
  quantity (the contract requires a finite value; ``0.0`` is a
  documented placeholder, not a claim of exact agreement). Only
  ``runtime_s`` (the external process's real wall-clock time) is a
  genuine measurement.

Measured accuracy
------------------
On a small uniform 8x8x10 grid (300 m cells, 10 earth z-cells — the
minimum :data:`_MIN_EARTH_Z_CELLS` requires), against a real compiled
``Mod3DMT``: :func:`~pycsamt.forward.maxwell.benchmarks.half_space_benchmark`
and :func:`~pycsamt.forward.maxwell.benchmarks.layered_earth_benchmark`
both pass the default
:class:`~pycsamt.forward.maxwell.benchmarks.BenchmarkThresholds` with
real margin (normalized RMS under 5%). Unlike
:mod:`~pycsamt.forward.maxwell.mt3d`'s research-only solver, no padded
mesh was needed here — ModEM synthesizes its own air layers and, as
an iterative-solver production code, does not carry this project's
own small-grid cell-budget restriction, so a plain uniform mesh
sufficed for this validation. This does not by itself validate every
mesh configuration (e.g. non-uniform meshes, receivers off-centre) --
only what the cited benchmarks actually exercise.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ...models.modem.config import ModEmConfig
from ...models.modem.data import ModEmData
from ...models.modem.model3d import ModEmModel3D
from .adapters import AdapterPolicy
from .backends import (
    BackendCapabilities,
    BackendRegistration,
    CompatibilityReport,
)
from .backends import register_backend as _register_backend
from .contracts import ForwardResult, MaxwellProblem, SolverDiagnostics
from .external import (
    BaseExternalMaxwellAdapter,
    ExternalRunPolicy,
    ExternalRunResult,
    make_availability_probe,
    resolve_executable,
)

__all__ = ["ModEm3DAdapter", "register_modem3d_backend"]

_MU0 = 4.0e-7 * np.pi
_SURFACE_TOLERANCE_M = 1e-6
_PERMEABILITY_RTOL = 1e-9
_FIELD_UNITS_TO_SI = 4.0e-4 * np.pi  # [mV/km]/[nT] -> V/A (Ohm)
_FULL_TENSOR = ("ZXX", "ZXY", "ZYX", "ZYY")

# ModEM's WS-format model reader hardcodes exactly 10 air layers
# (`nzAir = 10` in read_modelParam_ws, WS.inc) regardless of what
# n_air this adapter requests, and its "mirror" air-layer sizing
# (setup_airlayers, GridDef.f90) reads that many earth-layer widths
# without checking there are actually that many -- with fewer earth
# z-cells than this, it silently reads past the end of the array.
# Confirmed empirically (not just inferred from source) against a
# real compiled Mod3DMT binary: fewer than 10 earth z-cells produces
# a "b in QMR contains NaNs" crash from garbage memory; 10 or more
# solves correctly. See the module docstring's "Measured accuracy".
_MIN_EARTH_Z_CELLS = 10

# Confirmed 2026-07-29 against a real compiled Mod3DMT binary (built
# from the vendored source with a MinGW-w64 gfortran/OpenBLAS
# toolchain) -- not asserted from source reading alone. See the
# module docstring's "Measured accuracy" section for the actual
# numbers and the mesh constraint (_MIN_EARTH_Z_CELLS) required to
# reach them.
_VERIFIED_BENCHMARKS: tuple[str, ...] = ("half-space", "layered-earth")


def _default_search_paths() -> tuple[str, ...]:
    try:
        import pycsamt.models.modem as modem_package
    except ImportError:
        return ()
    source_root = Path(modem_package.__file__).resolve().parent / "_source"
    return (str(source_root / "3D"),)


def _unit_conversion_to_si(units: str) -> float:
    normalized = units.strip().lower().replace(" ", "")
    if normalized in ("[mv/km]/[nt]", "mv/km/nt"):
        return _FIELD_UNITS_TO_SI
    if normalized in ("[v/a]", "v/a", "ohm", "ohms", "si"):
        return 1.0
    raise ValueError(
        f"unrecognized ModEM data units {units!r}; expected "
        "'[mV/km]/[nT]' or '[V/A]'."
    )


class ModEm3DAdapter(BaseExternalMaxwellAdapter):
    """Adapter wrapping the external ModEM 3-D forward solver.

    Parameters
    ----------
    config : ModEmConfig or None, optional
        ModEM configuration (executable names, MPI settings, units,
        sign convention). Defaults to ``ModEmConfig(mode="3d")``.
    run_policy : ExternalRunPolicy or None, optional
        Executable resolution, timeout, retry, and working-directory
        rules. When omitted, a policy is built from ``config``
        (``config.mpi_command`` if ``config.use_mpi`` else
        ``config.binary_3d``), with the vendored
        ``pycsamt/models/modem/_source/3D`` directory as a search
        path fallback.
    policy : AdapterPolicy or None, optional
        Solver-independent result acceptance policy.
    predicted_data_filename : str or None, optional
        Exact name of ModEM's predicted-response output file within
        the run's working directory. When omitted, the adapter
        auto-detects it as the only ``*.dat`` file other than the
        one it wrote itself, and raises a clear error if that is
        ambiguous (zero or more than one candidate) — set this
        explicitly once you know your ModEM build's naming
        convention.
    version : str, default="modem-v6.2.6-adapter-1.0"
        Adapter version reported in every
        :class:`~pycsamt.forward.maxwell.contracts.ForwardResult`.
        Encodes both the vendored ModEM release this adapter targets
        and this adapter code's own revision.

    Examples
    --------
    >>> from pycsamt.forward.maxwell import ExternalRunPolicy
    >>> adapter = ModEm3DAdapter(
    ...     run_policy=ExternalRunPolicy("does-not-exist-xyz")
    ... )
    >>> adapter.capabilities.name
    'modem3d'
    """

    def __init__(
        self,
        *,
        config: ModEmConfig | None = None,
        run_policy: ExternalRunPolicy | None = None,
        policy: AdapterPolicy | None = None,
        predicted_data_filename: str | None = None,
        version: str = "modem-v6.2.6-adapter-1.0",
    ) -> None:
        self._config = config or ModEmConfig(mode="3d")
        if run_policy is None:
            executable = (
                self._config.mpi_command
                if self._config.use_mpi
                else self._config.binary_3d
            )
            run_policy = ExternalRunPolicy(
                executable, search_paths=_default_search_paths()
            )
        capabilities = BackendCapabilities(
            name="modem3d",
            version=version,
            dimensions=(3,),
            components=("zxx", "zxy", "zyx", "zyy"),
            time_conventions=("exp(+iwt)",),
            supports_nonuniform_mesh=True,
            supports_inactive_cells=False,
            supports_topography=False,
            supports_anisotropy=False,
            verified_benchmarks=_VERIFIED_BENCHMARKS,
        )
        super().__init__(capabilities, run_policy, policy)
        self._predicted_data_filename = predicted_data_filename

    def assess(self, problem: MaxwellProblem) -> CompatibilityReport:
        """Assess a problem, adding this adapter's mapping checks.

        Parameters
        ----------
        problem : MaxwellProblem
            Candidate simulation problem.

        Returns
        -------
        CompatibilityReport
            The generic capability report from
            :meth:`~pycsamt.forward.maxwell.adapters.BaseMaxwellAdapter.assess`,
            extended with surface-receiver, horizontal-bounds,
            surface-aligned-mesh, and vacuum-permeability checks.

        Examples
        --------
        See :class:`ModEm3DAdapter` for construction; a problem with
        a buried receiver is rejected before any file is written.
        """
        base = super().assess(problem)
        errors = list(base.errors)
        if problem.mesh.dimension == 3:
            coordinates = problem.receivers.coordinates_m
            depth = coordinates[:, 2]
            if np.any(np.abs(depth) > _SURFACE_TOLERANCE_M):
                errors.append(
                    "modem3d only evaluates receivers at the surface (z=0)"
                )
            x, y = coordinates[:, 0], coordinates[:, 1]
            x_lo, x_hi = (
                problem.mesh.x_edges_m[0],
                problem.mesh.x_edges_m[-1],
            )
            y_lo, y_hi = (
                problem.mesh.y_edges_m[0],
                problem.mesh.y_edges_m[-1],
            )
            if np.any((x < x_lo) | (x > x_hi) | (y < y_lo) | (y > y_hi)):
                errors.append(
                    "receiver x/y coordinates must lie within the "
                    f"mesh [{x_lo:g}, {x_hi:g}] x [{y_lo:g}, {y_hi:g}] m"
                )
            if abs(problem.mesh.z_edges_m[0]) > _SURFACE_TOLERANCE_M:
                errors.append(
                    "modem3d requires the mesh's top z edge at 0 "
                    "(no separate air layers are written; see the "
                    "module docstring)"
                )
            n_earth_z = problem.mesh.shape[0]
            if n_earth_z < _MIN_EARTH_Z_CELLS:
                errors.append(
                    f"modem3d requires at least {_MIN_EARTH_Z_CELLS} earth "
                    f"z-cells (got {n_earth_z}); ModEM's WS-format reader "
                    "hardcodes 10 air layers and its 'mirror' air-layer "
                    "sizing reads that many earth-layer widths without "
                    "bounds-checking against the actual earth cell count, "
                    "silently reading past the end of the array (undefined "
                    "values, observed as a 'b in QMR contains NaNs' crash) "
                    "when there are fewer earth cells than air layers -- "
                    "confirmed empirically with a real compiled Mod3DMT "
                    "binary, not merely inferred from source."
                )
        if not np.isclose(
            problem.magnetic_permeability_h_m,
            _MU0,
            rtol=_PERMEABILITY_RTOL,
        ):
            errors.append(
                "modem3d assumes vacuum magnetic permeability "
                f"({_MU0:.6g} H/m)"
            )
        if len(errors) == len(base.errors):
            return base
        return CompatibilityReport(
            self.capabilities.name, False, tuple(errors), base.warnings
        )

    def _prepare_inputs(self, problem: MaxwellProblem, workdir: Path):
        model = ModEmModel3D(config=self._config)
        widths = problem.mesh.cell_widths_m
        model.x_widths = np.array(widths["x"], dtype=float)
        model.y_widths = np.array(widths["y"], dtype=float)
        model.z_widths = np.array(widths["z"], dtype=float)
        model.n_air = 0
        model.log_type = "LOGE"
        model.rho_loge = np.log(1.0 / problem.conductivity_s_m)
        model_path = model.write(workdir / "model.ws")

        x0 = float(problem.mesh.x_edges_m[0])
        y0 = float(problem.mesh.y_edges_m[0])
        names = list(problem.receivers.names)
        coordinates_m = problem.receivers.coordinates_m

        data = ModEmData(config=self._config)
        data.site_names = names
        data.site_coords = {
            name: (
                float(coordinates_m[i, 0] - x0),
                float(coordinates_m[i, 1] - y0),
                0.0,
            )
            for i, name in enumerate(names)
        }
        periods = 1.0 / np.asarray(problem.frequencies_hz, dtype=float)
        data.periods = np.sort(periods)[::-1]
        rows = []
        for period in periods:
            for site_idx, name in enumerate(names):
                x, y, z = data.site_coords[name]
                for component in _FULL_TENSOR:
                    rows.append(
                        (
                            float(period),
                            site_idx,
                            x,
                            y,
                            z,
                            component,
                            0.0,
                            0.0,
                            1.0,
                        )
                    )
        data.blocks = [
            {
                "component_type": "Full_Impedance",
                "sign_convention": self._config.sign_convention,
                "units": self._config.units,
                "rotation_angle": 0.0,
                "origin": [0.0, 0.0],
                "n_periods": len(periods),
                "n_sites": len(names),
                "rows": rows,
            }
        ]
        data_path = data.write(workdir / "data.dat")
        predicted_path = workdir / "pred.dat"
        return {
            "model_path": model_path,
            "data_path": data_path,
            "predicted_path": predicted_path,
        }

    def _build_command(self, problem, workdir, executable, context):
        model_name = context["model_path"].name
        data_name = context["data_path"].name
        # ModEM's "-F" forward mode requires a third positional
        # argument naming the predicted-data output file
        # (usage: "-F rFile_Model rFile_Data wFile_Data [...]");
        # omitting it (as this adapter originally did) makes Mod3DMT
        # print its own usage banner (a "success", exit code 0) and
        # write nothing, which _locate_predicted_file then reports as
        # a generic "no predicted-data file was found" -- confirmed
        # against a real compiled Mod3DMT binary.
        predicted_name = context["predicted_path"].name
        if self._config.use_mpi:
            binary = resolve_executable(
                self._config.binary_3d,
                search_paths=self.run_policy.search_paths,
            )
            return [
                str(executable),
                "-np",
                str(self._config.n_procs),
                str(binary),
                "-F",
                model_name,
                data_name,
                predicted_name,
            ]
        return [str(executable), "-F", model_name, data_name, predicted_name]

    def _locate_predicted_file(self, workdir: Path, data_path: Path) -> Path:
        if self._predicted_data_filename is not None:
            return workdir / self._predicted_data_filename
        candidates = [
            path
            for path in workdir.glob("*.dat")
            if path.name != data_path.name
        ]
        if len(candidates) == 1:
            return candidates[0]
        if not candidates:
            raise ValueError(
                "no ModEM predicted-data file was found in the working "
                "directory; pass predicted_data_filename explicitly if "
                "your ModEM build uses a naming convention this adapter "
                "does not already recognize."
            )
        names = ", ".join(sorted(path.name for path in candidates))
        raise ValueError(
            "multiple candidate ModEM output files were found "
            f"({names}); pass predicted_data_filename explicitly to "
            "disambiguate."
        )

    def _parse_result(
        self,
        problem: MaxwellProblem,
        workdir: Path,
        run_result: ExternalRunResult,
        context,
    ) -> ForwardResult:
        predicted_path = self._locate_predicted_file(
            workdir, context["data_path"]
        )
        predicted = ModEmData.read(predicted_path)

        lookup: dict[tuple[str, str], complex] = {}
        periods_seen: list[float] = []
        for block in predicted.blocks:
            factor = _unit_conversion_to_si(
                block.get("units", self._config.units)
            )
            for row in block["rows"]:
                period, site_idx, _x, _y, _z, comp, real, imag, _err = row
                name = predicted.site_names[site_idx]
                lookup[(round(period, 6), name, comp.upper())] = (
                    factor * complex(real, imag)
                )
                periods_seen.append(period)
        period_array = np.array(sorted(set(periods_seen)))

        n_station = problem.receivers.count
        n_freq = len(problem.frequencies_hz)
        n_component = len(problem.components)
        impedance = np.empty((n_station, n_freq, n_component), dtype=complex)
        component_map = {
            "zxx": "ZXX",
            "zxy": "ZXY",
            "zyx": "ZYX",
            "zyy": "ZYY",
        }

        for fi, frequency in enumerate(problem.frequencies_hz):
            target_period = 1.0 / float(frequency)
            if period_array.size == 0:
                raise ValueError("ModEM predicted output contained no data.")
            nearest = period_array[
                np.argmin(np.abs(period_array - target_period))
            ]
            if abs(nearest - target_period) > 1e-6 * max(target_period, 1e-12):
                raise ValueError(
                    f"no ModEM predicted period matches {target_period:g} s "
                    f"(requested from {frequency:g} Hz) within tolerance."
                )
            key_period = round(float(nearest), 6)
            for si, name in enumerate(problem.receivers.names):
                for ci, component in enumerate(problem.components):
                    key = (key_period, name, component_map[component])
                    if key not in lookup:
                        raise ValueError(
                            f"ModEM predicted output is missing station "
                            f"{name!r}, component {component!r} at period "
                            f"{nearest:g} s."
                        )
                    impedance[si, fi, ci] = lookup[key]

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


def register_modem3d_backend(
    *, config: ModEmConfig | None = None, replace: bool = False
) -> None:
    """Register :class:`ModEm3DAdapter` in the process-wide registry.

    Parameters
    ----------
    config : ModEmConfig or None, optional
        Configuration used both to build the registered adapter and
        to derive its availability probe (whether the configured
        executable can currently be resolved).
    replace : bool, default=False
        Explicitly replace an existing ``"modem3d"`` registration.

    Examples
    --------
    >>> register_modem3d_backend(replace=True)
    >>> from pycsamt.forward.maxwell import list_backends
    >>> "modem3d" in list_backends()
    True
    >>> list_backends()["modem3d"]["available"]
    False
    """
    cfg = config or ModEmConfig(mode="3d")
    executable = cfg.mpi_command if cfg.use_mpi else cfg.binary_3d
    capabilities = ModEm3DAdapter(config=cfg).capabilities
    probe = make_availability_probe(
        executable, search_paths=_default_search_paths()
    )

    def factory(**options):
        return ModEm3DAdapter(config=cfg, **options)

    _register_backend(
        BackendRegistration(capabilities, factory, probe), replace=replace
    )
