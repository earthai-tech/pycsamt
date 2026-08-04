# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Generate versioned 2-D triangular-mesh training datasets (M14).

Triangular-mesh counterpart of :mod:`pycsamt.ai.training.dataset2d`:
correlated 2-D resistivity fields from :mod:`pycsamt.ai.geology` are
resampled onto triangle centroids, solved with
:class:`~pycsamt.forward.maxwell.mare2dem.Mare2DEMAdapter` through the
existing (solver-agnostic) :func:`~pycsamt.forward.maxwell.batch.solve_batch`,
and packaged into a :class:`MaxwellTri2DDataset` with a realization-level
split and :class:`~pycsamt.ai.data.manifest.DatasetManifest`, exactly as
:mod:`~pycsamt.ai.training.dataset2d` does for the rectilinear path.

Mesh generation and per-triangle resistivity
---------------------------------------------
``Mare2DEMAdapter`` parameterizes resistivity by named mesh region, not
by individual triangle (see its own module docstring). Genuinely
heterogeneous per-triangle training fields therefore need a mesh with
**one region per triangle**. This module builds that mesh with
:func:`~pycsamt.forward.maxwell.tri_mesh_gen.build_graded_tri_mesh` --
a real, quality-graded triangulation via Shewchuk's Triangle library
(fine near the receivers, growing smoothly with depth and lateral
distance), rather than through
:mod:`pycsamt.models.mare2dem.tri_mesh`'s Triangle-*subprocess* path
(:func:`~pycsamt.models.mare2dem.tri_mesh.build_survey_mesh`), which
needs a separately compiled MARE2DEM/Triangle binary. Using Triangle's
Python bindings in-process keeps that same usability property --
generating training data needs no compiled binary -- without giving up
real mesh quality: an earlier version of this module triangulated a
uniform rectangular lattice with plain :func:`scipy.spatial.Delaunay`,
which produced a mesh of uniform triangle size everywhere and no
grading at all; see :func:`~pycsamt.forward.maxwell.tri_mesh_gen.build_graded_tri_mesh`'s
own docstring for why that mattered beyond appearance.

**Status: unverified physics**, following exactly from
:mod:`~pycsamt.forward.maxwell.mare2dem`'s own honesty note. Generating
a real dataset (not just exercising this module's own data plumbing)
still requires a real compiled MARE2DEM binary.

**Known representability gap against the real solver.** Since the
adapter's 2026-07-31 validation, ``Mare2DEMAdapter`` writes ``mesh``
as a PSLG recipe and lets MARE2DEM mesh and adaptively refine it
internally; the triangles it actually solves on are not guaranteed to
match ``mesh.triangles``, and its ``.resistivity`` file only carries
one value per *named region*, not per triangle (see
``mare2dem.py``'s module docstring). This module's "one region per
triangle" design satisfies ``assess()``'s region-uniformity check only
because each region happens to contain exactly one triangle -- it does
*not* mean the real binary preserves this module's exact triangulation
or per-triangle field when it re-meshes the PSLG. A real (not
scripted-stand-in) generated dataset should therefore be treated as
approximate until cross-checked against the returned response, same
caution as any adaptively-refined external mesh. The in-house
triangular FEM solver (a separate, not-yet-started milestone) is the
eventual fix, since it would solve on this module's own mesh directly
instead of a re-meshed PSLG approximation of it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ...forward.maxwell import MaxwellResultCache, ReceiverSet, TriMesh, TriProblem
from ...forward.maxwell.batch import BatchPolicy, solve_batch
from ...forward.maxwell.mare2dem import Mare2DEMAdapter
from ...forward.maxwell.tri_mesh_gen import build_graded_tri_mesh
from ..data.contracts import SurveyData
from ..data.manifest import DatasetManifest
from ..data.splits import RealizationSplit, split_realizations
from ..experiments.config import SeedPlan
from ..geology import GaussianCorrelation, GeologyGrid, generate_gaussian_field

__all__ = [
    "MaxwellTri2DSample",
    "MaxwellTri2DDataset",
    "MaxwellTri2DDatasetConfig",
    "generate_2d_tri_maxwell_dataset",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _nearest_index(centres: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Return the nearest regular-grid index for each of *values*."""
    spacing = float(centres[1] - centres[0])
    index = np.round((values - centres[0]) / spacing).astype(int)
    return np.clip(index, 0, len(centres) - 1)


def _resistivity_per_triangle(
    mesh: TriMesh,
    grid: GeologyGrid,
    correlation: GaussianCorrelation,
    seed: int,
    *,
    log_resistivity_mean: float,
    log_resistivity_std: float,
) -> np.ndarray:
    """Sample a correlated field on *grid*, nearest-mapped to triangle centroids."""
    field = generate_gaussian_field(grid, correlation, seed=seed)
    centroids = mesh.triangle_centroids_m
    ix = _nearest_index(grid.x_m, centroids[:, 0])
    iz = _nearest_index(grid.z_m, centroids[:, 1])
    standardized = field.values[iz, ix]
    log_resistivity = log_resistivity_mean + log_resistivity_std * standardized
    return np.power(10.0, log_resistivity)


def _sample_correlation(
    config: "MaxwellTri2DDatasetConfig", seed: int
) -> GaussianCorrelation:
    rng = np.random.default_rng(seed)
    length_x = rng.uniform(*config.correlation_length_x_m)
    length_z = rng.uniform(*config.correlation_length_z_m)
    return GaussianCorrelation(length_x_m=length_x, length_z_m=length_z)


@dataclass(frozen=True)
class MaxwellTri2DDatasetConfig:
    """Immutable configuration for :func:`generate_2d_tri_maxwell_dataset`.

    Parameters
    ----------
    dataset_id : str
        Portable dataset identifier.
    x_range_m, z_range_m : (float, float)
        Shared mesh domain; ``z_range_m[0]`` must be 0 (surface).
    correlation_length_x_m, correlation_length_z_m : (float, float)
        Inclusive ``(low, high)`` ranges a correlation length is drawn
        from per realization.
    frequencies_hz : array-like
        Shared positive, unique simulated frequencies.
    station_x_m : array-like
        Shared receiver x-positions within ``x_range_m``.
    n_realizations : int
        Number of geological realizations to attempt.
    seed : int
        Root seed; every realization derives a labeled child seed via
        :class:`~pycsamt.ai.experiments.config.SeedPlan`.
    log_resistivity_mean, log_resistivity_std : float
        Affine map from the standardized correlated field to
        ``log10(resistivity_ohm_m)``.
    field_grid_cell_m : float
        Cell size of the regular :class:`~pycsamt.ai.geology.GeologyGrid`
        the correlated field is synthesized on before nearest-mapping
        onto triangle centroids. Should be finer than
        ``mesh_target_cell_m`` to resolve the field at triangle scale.
    mesh_target_cell_m : float
        Forwarded to
        :func:`~pycsamt.forward.maxwell.tri_mesh_gen.build_graded_tri_mesh`
        as ``surface_cell_m`` -- the target triangle edge length at a
        station, graded coarser with depth/lateral distance from there.
    topo_x_m, topo_z_m : array-like, optional
        Topography polyline (z positive down), forwarded to
        :func:`~pycsamt.forward.maxwell.tri_mesh_gen.build_graded_tri_mesh`.
        Both default to ``None`` (flat surface at ``z_range_m[0]=0``,
        unchanged from before this parameter existed). When given, every
        station's receiver z comes from interpolating this same
        polyline rather than ``0.0``, and ``z_range_m[0]``'s "must be
        depth 0" requirement no longer applies (the top is no longer
        flat).
    validation_fraction, test_fraction : float, default=0.1
        Forwarded to :func:`~pycsamt.ai.data.splits.split_realizations`.
    namespace : str, default="maxwell2dtri"
        :class:`SeedPlan` namespace.

    Examples
    --------
    >>> config = MaxwellTri2DDatasetConfig(
    ...     dataset_id="demo-2dtri-v1",
    ...     x_range_m=(0.0, 1000.0),
    ...     z_range_m=(0.0, 500.0),
    ...     correlation_length_x_m=(200.0, 400.0),
    ...     correlation_length_z_m=(100.0, 200.0),
    ...     frequencies_hz=[10.0, 1.0],
    ...     station_x_m=[200.0, 500.0, 800.0],
    ...     n_realizations=2,
    ...     seed=0,
    ... )
    >>> config.components
    ('zxy', 'zyx')
    """

    dataset_id: str
    x_range_m: tuple[float, float]
    z_range_m: tuple[float, float]
    correlation_length_x_m: tuple[float, float]
    correlation_length_z_m: tuple[float, float]
    frequencies_hz: Any
    station_x_m: Any
    n_realizations: int
    seed: int
    log_resistivity_mean: float = 2.0
    log_resistivity_std: float = 0.4
    components: tuple[str, ...] = ("zxy", "zyx")
    field_grid_cell_m: float = 50.0
    mesh_target_cell_m: float = 100.0
    topo_x_m: Any = None
    topo_z_m: Any = None
    validation_fraction: float = 0.1
    test_fraction: float = 0.1
    namespace: str = "maxwell2dtri"
    verbose: bool | int | str = False

    def __post_init__(self) -> None:
        x0, x1 = (float(self.x_range_m[0]), float(self.x_range_m[1]))
        z0, z1 = (float(self.z_range_m[0]), float(self.z_range_m[1]))
        if not (np.isfinite(x0) and np.isfinite(x1)) or x1 <= x0:
            raise ValueError("x_range_m must be finite and increasing.")
        if not (np.isfinite(z0) and np.isfinite(z1)) or z1 <= z0:
            raise ValueError("z_range_m must be finite and increasing.")
        has_topo = self.topo_x_m is not None or self.topo_z_m is not None
        if not has_topo and not np.isclose(z0, 0.0, atol=1e-6):
            raise ValueError(
                f"z_range_m[0] must be depth 0 (got {z0} m)."
            )
        object.__setattr__(self, "x_range_m", (x0, x1))
        object.__setattr__(self, "z_range_m", (z0, z1))
        if (self.topo_x_m is None) != (self.topo_z_m is None):
            raise ValueError("topo_x_m and topo_z_m must be given together.")
        if has_topo:
            topo_x = np.asarray(self.topo_x_m, dtype=float).ravel()
            topo_z = np.asarray(self.topo_z_m, dtype=float).ravel()
            if (
                topo_x.shape != topo_z.shape
                or topo_x.size < 2
                or not np.all(np.isfinite(topo_x))
                or not np.all(np.isfinite(topo_z))
            ):
                raise ValueError(
                    "topo_x_m/topo_z_m must be finite, equal-length, "
                    "and have at least two samples."
                )
            object.__setattr__(self, "topo_x_m", _readonly(topo_x))
            object.__setattr__(self, "topo_z_m", _readonly(topo_z))
        for name, value in (
            ("correlation_length_x_m", self.correlation_length_x_m),
            ("correlation_length_z_m", self.correlation_length_z_m),
        ):
            low, high = (float(value[0]), float(value[1]))
            if not (np.isfinite(low) and np.isfinite(high)):
                raise ValueError(f"{name} must be finite.")
            if low <= 0 or high < low:
                raise ValueError(f"{name} must satisfy 0 < low <= high.")
            object.__setattr__(self, name, (low, high))
        frequencies = np.asarray(self.frequencies_hz, dtype=float)
        if (
            frequencies.ndim != 1
            or frequencies.size < 1
            or not np.all(np.isfinite(frequencies))
            or np.any(frequencies <= 0)
            or len(set(frequencies.tolist())) != frequencies.size
        ):
            raise ValueError(
                "frequencies_hz must be a non-empty, positive, "
                "finite, unique vector."
            )
        station_x = np.asarray(self.station_x_m, dtype=float)
        if (
            station_x.ndim != 1
            or station_x.size < 1
            or not np.all(np.isfinite(station_x))
            or np.any(station_x < x0)
            or np.any(station_x > x1)
            or len(set(station_x.tolist())) != station_x.size
        ):
            raise ValueError(
                "station_x_m must be a non-empty, finite, unique "
                f"vector within x_range_m {(x0, x1)}."
            )
        components = tuple(str(c).strip().lower() for c in self.components)
        if not components or any(c not in ("zxy", "zyx") for c in components):
            raise ValueError(
                "components must be a non-empty subset of "
                "('zxy', 'zyx') for a 2-D triangular problem."
            )
        if (
            not isinstance(self.n_realizations, int)
            or isinstance(self.n_realizations, bool)
            or self.n_realizations < 1
        ):
            raise ValueError("n_realizations must be a positive integer.")
        if (
            not isinstance(self.seed, int)
            or isinstance(self.seed, bool)
            or self.seed < 0
        ):
            raise ValueError("seed must be a non-negative integer.")
        if not np.isfinite(self.log_resistivity_std) or (
            self.log_resistivity_std <= 0
        ):
            raise ValueError(
                "log_resistivity_std must be finite and positive."
            )
        if not np.isfinite(self.field_grid_cell_m) or self.field_grid_cell_m <= 0:
            raise ValueError("field_grid_cell_m must be finite and positive.")
        if not np.isfinite(self.mesh_target_cell_m) or self.mesh_target_cell_m <= 0:
            raise ValueError("mesh_target_cell_m must be finite and positive.")
        namespace = str(self.namespace).strip()
        if not namespace:
            raise ValueError("namespace cannot be empty.")
        object.__setattr__(self, "frequencies_hz", _readonly(frequencies))
        object.__setattr__(self, "station_x_m", _readonly(station_x))
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "namespace", namespace)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable configuration record.

        Examples
        --------
        >>> config = MaxwellTri2DDatasetConfig(
        ...     "d", (0.0, 100.0), (0.0, 50.0), (10.0, 20.0), (10.0, 20.0),
        ...     [1.0], [10.0, 50.0], 1, 0,
        ... )
        >>> config.to_dict()["dataset_id"]
        'd'
        """
        return {
            "schema_version": 1,
            "dataset_id": self.dataset_id,
            "x_range_m": list(self.x_range_m),
            "z_range_m": list(self.z_range_m),
            "correlation_length_x_m": list(self.correlation_length_x_m),
            "correlation_length_z_m": list(self.correlation_length_z_m),
            "frequencies_hz": self.frequencies_hz.tolist(),
            "station_x_m": self.station_x_m.tolist(),
            "n_realizations": self.n_realizations,
            "seed": self.seed,
            "log_resistivity_mean": self.log_resistivity_mean,
            "log_resistivity_std": self.log_resistivity_std,
            "components": list(self.components),
            "field_grid_cell_m": self.field_grid_cell_m,
            "mesh_target_cell_m": self.mesh_target_cell_m,
            "topo_x_m": None if self.topo_x_m is None else self.topo_x_m.tolist(),
            "topo_z_m": None if self.topo_z_m is None else self.topo_z_m.tolist(),
            "validation_fraction": self.validation_fraction,
            "test_fraction": self.test_fraction,
            "namespace": self.namespace,
        }


@dataclass(frozen=True)
class MaxwellTri2DSample:
    """One converged 2-D triangular-mesh realization.

    Parameters
    ----------
    realization_id : str
        Unique identifier, stable across regeneration with the same
        configuration.
    seed : int
        Seed used to generate this realization's correlated field.
    correlation : GaussianCorrelation
        Correlation lengths drawn for this realization.
    resistivity_ohm_m : ndarray, shape (n_triangles,)
        True per-triangle resistivity, the training target.
    survey : SurveyData
        Simulated, noiseless response at the configured stations and
        frequencies, the training input.
    relative_residual : float
        Worst-case relative solver residual over every frequency.
    """

    realization_id: str
    seed: int
    correlation: GaussianCorrelation
    resistivity_ohm_m: np.ndarray
    survey: SurveyData
    relative_residual: float

    def __post_init__(self) -> None:
        realization_id = str(self.realization_id).strip()
        if not realization_id:
            raise ValueError("realization_id cannot be empty.")
        if not isinstance(self.correlation, GaussianCorrelation):
            raise TypeError("correlation must be a GaussianCorrelation.")
        if not isinstance(self.survey, SurveyData):
            raise TypeError("survey must be a SurveyData.")
        resistivity = np.asarray(self.resistivity_ohm_m, dtype=float)
        if not np.all(np.isfinite(resistivity)) or np.any(resistivity <= 0):
            raise ValueError("resistivity_ohm_m must be positive and finite.")
        object.__setattr__(self, "realization_id", realization_id)
        object.__setattr__(self, "seed", int(self.seed))
        object.__setattr__(self, "resistivity_ohm_m", _readonly(resistivity))
        object.__setattr__(
            self, "relative_residual", float(self.relative_residual)
        )


@dataclass(frozen=True)
class MaxwellTri2DDataset:
    """Immutable collection of converged triangular-mesh realizations.

    Parameters
    ----------
    mesh : TriMesh
        Shared triangular geometry of every sample's
        ``resistivity_ohm_m`` (one region per triangle).
    samples : tuple of MaxwellTri2DSample
        Converged realizations, at least one.
    split : RealizationSplit
        Realization-level train/validation/test assignment.
    manifest : DatasetManifest
        Reproducibility record for this dataset.
    rejected : tuple of str
        Realization IDs excluded for not converging.
    """

    mesh: TriMesh
    samples: tuple[MaxwellTri2DSample, ...]
    split: RealizationSplit
    manifest: DatasetManifest
    rejected: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.mesh, TriMesh):
            raise TypeError("mesh must be a TriMesh.")
        samples = tuple(self.samples)
        if not samples:
            raise ValueError("samples must contain at least one entry.")
        if any(not isinstance(s, MaxwellTri2DSample) for s in samples):
            raise TypeError("samples entries must be MaxwellTri2DSample.")
        if not isinstance(self.split, RealizationSplit):
            raise TypeError("split must be a RealizationSplit.")
        if not isinstance(self.manifest, DatasetManifest):
            raise TypeError("manifest must be a DatasetManifest.")
        sample_ids = {s.realization_id for s in samples}
        split_ids = set(self.split.train) | set(self.split.validation)
        split_ids |= set(self.split.test)
        if sample_ids != split_ids:
            raise ValueError("split must cover exactly every sample id.")
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "rejected", tuple(self.rejected))

    def select(self, partition: str) -> tuple[MaxwellTri2DSample, ...]:
        """Return samples belonging to one split partition.

        Examples
        --------
        Samples in the training partition are disjoint from testing:

        >>> ids_train = {
        ...     s.realization_id for s in dataset.select("train")
        ... }  # doctest: +SKIP
        """
        return tuple(
            s
            for s in self.samples
            if self.split.partition(s.realization_id) == partition
        )


def generate_2d_tri_maxwell_dataset(
    config: MaxwellTri2DDatasetConfig,
    *,
    cache: MaxwellResultCache | None = None,
    policy: BatchPolicy | None = None,
    adapter: Mare2DEMAdapter | None = None,
) -> MaxwellTri2DDataset:
    """Generate a cached, split, versioned triangular-mesh 2-D dataset.

    Parameters
    ----------
    config : MaxwellTri2DDatasetConfig
        Complete, immutable dataset configuration.
    cache : MaxwellResultCache or None, optional
        When given, forward solves are cached by problem hash and
        resumable across calls.
    policy : BatchPolicy or None, optional
        Retry/concurrency policy forwarded to
        :func:`~pycsamt.forward.maxwell.batch.solve_batch`.
    adapter : Mare2DEMAdapter or None, optional
        Forward-solve backend. Defaults to ``Mare2DEMAdapter()``
        (resolved from the environment). Accepting a pre-built adapter
        rather than only executable-resolution options lets tests
        inject one pointed at a scripted stand-in executable, exactly
        as ``test_maxwell_mare2dem.py`` does for the adapter itself.

    Returns
    -------
    MaxwellTri2DDataset
        Every realization whose forward solve converged, split at the
        realization level, sharing one triangular mesh.

    Raises
    ------
    ValueError
        If no realization converges.

    Examples
    --------
    Building the mesh and geological realizations needs no external
    binary; only the forward solve does, so this example is
    illustrative only::

        >>> config = MaxwellTri2DDatasetConfig(
        ...     dataset_id="demo-2dtri-v1",
        ...     x_range_m=(0.0, 1000.0),
        ...     z_range_m=(0.0, 500.0),
        ...     correlation_length_x_m=(200.0, 400.0),
        ...     correlation_length_z_m=(100.0, 200.0),
        ...     frequencies_hz=[10.0, 1.0],
        ...     station_x_m=[200.0, 500.0, 800.0],
        ...     n_realizations=2,
        ...     seed=0,
        ... )
        >>> dataset = generate_2d_tri_maxwell_dataset(config)  # doctest: +SKIP
    """
    plan = SeedPlan(config.seed, namespace=config.namespace)
    has_topo = config.topo_x_m is not None
    mesh = build_graded_tri_mesh(
        config.x_range_m,
        config.z_range_m,
        config.station_x_m,
        surface_cell_m=config.mesh_target_cell_m,
        topo_x_m=config.topo_x_m,
        topo_z_m=config.topo_z_m,
    )
    # Computed once and reused for both the mesh builder (above) and the
    # receiver coordinates (below) so a receiver's z is bit-for-bit
    # identical to its mesh node's z -- required for
    # TriFEM2DAdapter's exact-match receiver/surface checks to agree.
    station_z = (
        np.interp(config.station_x_m, config.topo_x_m, config.topo_z_m)
        if has_topo
        else np.zeros_like(config.station_x_m)
    )
    field_nx = max(
        int(round((config.x_range_m[1] - config.x_range_m[0]) / config.field_grid_cell_m)),
        2,
    )
    field_nz = max(
        int(round((config.z_range_m[1] - config.z_range_m[0]) / config.field_grid_cell_m)),
        2,
    )
    field_grid = GeologyGrid.regular_2d(
        nx=field_nx,
        nz=field_nz,
        dx_m=config.field_grid_cell_m,
        dz_m=config.field_grid_cell_m,
        x_origin_m=config.x_range_m[0],
        z_origin_m=config.z_range_m[0],
    )

    station_names = tuple(f"S{i:04d}" for i in range(len(config.station_x_m)))
    receivers = ReceiverSet(
        [[float(x), float(z)] for x, z in zip(config.station_x_m, station_z)],
        station_names,
    )
    coordinates_m = [
        [float(x), 0.0, float(z)] for x, z in zip(config.station_x_m, station_z)
    ]

    from pycsamt.api.view.progress import get_progress_bar

    realizations: dict[str, dict[str, Any]] = {}
    problems: dict[str, TriProblem] = {}
    with get_progress_bar(
        total=config.n_realizations,
        desc="Generating 2-D triangular realizations",
        unit="realization",
        verbose=config.verbose,
    ) as bar:
        for index in range(config.n_realizations):
            realization_id = f"{config.dataset_id}-r{index:05d}"
            correlation_seed = plan.derive(f"{realization_id}/correlation")
            field_seed = plan.derive(f"{realization_id}/field")
            correlation = _sample_correlation(config, correlation_seed)
            resistivity = _resistivity_per_triangle(
                mesh,
                field_grid,
                correlation,
                field_seed,
                log_resistivity_mean=config.log_resistivity_mean,
                log_resistivity_std=config.log_resistivity_std,
            )
            problem = TriProblem(
                mesh,
                1.0 / resistivity,
                config.frequencies_hz,
                receivers,
                components=config.components,
                metadata={"realization_id": realization_id},
            )
            if problem.problem_hash in problems:
                raise RuntimeError(
                    "duplicate problem_hash across distinct realizations; "
                    "this should not happen."
                )
            realizations[problem.problem_hash] = {
                "realization_id": realization_id,
                "seed": field_seed,
                "correlation": correlation,
                "resistivity_ohm_m": resistivity,
            }
            problems[problem.problem_hash] = problem
            bar.update(1)

    samples: list[MaxwellTri2DSample] = []
    rejected: list[str] = []

    with get_progress_bar(
        total=len(problems),
        desc="Solving 2-D triangular Maxwell realizations",
        unit="problem",
        verbose=config.verbose,
    ) as solve_bar:

        def _on_result(problem: TriProblem, result: Any) -> None:
            meta = realizations[problem.problem_hash]
            if not result.success:
                rejected.append(meta["realization_id"])
                solve_bar.update(1)
                return
            survey = SurveyData(
                impedance=result.impedance_v_a,
                frequencies_hz=result.frequencies_hz,
                station_names=result.receiver_names,
                components=result.components,
                coordinates_m=coordinates_m,
                valid=result.valid,
            )
            samples.append(
                MaxwellTri2DSample(
                    realization_id=meta["realization_id"],
                    seed=meta["seed"],
                    correlation=meta["correlation"],
                    resistivity_ohm_m=meta["resistivity_ohm_m"],
                    survey=survey,
                    relative_residual=(
                        result.diagnostics.maximum_relative_residual
                    ),
                )
            )
            solve_bar.update(1)

        def _on_failure(problem: TriProblem, failure: Any) -> None:
            rejected.append(
                realizations[problem.problem_hash]["realization_id"]
            )
            solve_bar.update(1)

        solve_batch(
            problems.values(),
            adapter if adapter is not None else Mare2DEMAdapter(),
            cache=cache,
            policy=policy,
            on_result=_on_result,
            on_failure=_on_failure,
        )

    if not samples:
        raise ValueError(
            "no realization produced a converged forward solve; check "
            "that a MARE2DEM executable is available, or pass a "
            "cache pre-populated by a scripted stand-in for testing."
        )

    realization_ids = [sample.realization_id for sample in samples]
    split = split_realizations(
        realization_ids,
        validation_fraction=config.validation_fraction,
        test_fraction=config.test_fraction,
        seed=plan.derive("split"),
    )
    manifest = DatasetManifest(
        dataset_id=config.dataset_id,
        generator=(
            "pycsamt.ai.training.dataset2d_tri.generate_2d_tri_maxwell_dataset"
        ),
        generator_version="0.1.0",
        configuration=config.to_dict(),
        split=split,
        sample_count=len(samples),
    )
    return MaxwellTri2DDataset(
        mesh=mesh,
        samples=tuple(samples),
        split=split,
        manifest=manifest,
        rejected=tuple(rejected),
    )
