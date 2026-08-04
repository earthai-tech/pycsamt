# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Generate versioned 2-D Maxwell training datasets (M5 groundwork).

This module replaces the tiled station-wise 1-D approximation used by
the current ``Inv2DAgent`` with genuine 2-D physics, per the
AI-inversion plan's M4/M5 milestones: correlated 2-D resistivity
models from :mod:`pycsamt.ai.geology` are solved with the verified
:class:`~pycsamt.forward.maxwell.mt2d.MT2DAdapter`, cached and
retried through :func:`~pycsamt.forward.maxwell.batch.solve_batch`,
and packaged into a :class:`Maxwell2DDataset` with a realization-level
train/validation/test split and a ``DatasetManifest``
(:mod:`pycsamt.ai.data.manifest`).

A realization whose solve does not converge (see
:attr:`~pycsamt.forward.maxwell.contracts.ForwardResult.success`) is
excluded from the dataset rather than silently included, per the
plan's requirement that failed solves cannot silently enter training
data.

Solver mesh construction deliberately does *not* use
:func:`~pycsamt.forward.maxwell.mesh.build_solver_mesh`.  A mesh built
from a small uniform core plus separately geometrically-grown padding
corrupted ``MT2DAdapter``'s TM-mode (``zyx``) response for a laterally
uniform, vertically layered earth, even with ample lateral/vertical
extent by every standard mesh-quality check.  A single continuous
geometric z-progression with *uniform* lateral spacing, sized
generously from that solver's own calibrated test
(``forward/tests/test_maxwell_mt2d.py``), reproduces analytic
half-space and layered-earth solutions in both modes to within a few
percent.

``zyx`` (TM mode) is validated and included by default alongside
``zxy``. It was previously excluded here after appearing to carry a
30-65% error for laterally fine 2-D structure that survived every
mesh-parameter change tried; that investigation, however, always
varied lateral and vertical resolution independently rather than
together, and separately regenerated a new random field realization
at every resolution instead of resampling one fixed field. Both
turned out to matter: TM mode's governing equation carries the full
resistivity structure inside its diffusion coefficient (unlike TE
mode, where only the mild reaction term does), so it genuinely needs
finer *joint* x/z resolution than TE mode to reach comparable
accuracy, and comparing different realizations at different
resolutions looks exactly like non-convergence even when the solver
itself converges properly. Once resampling one field across
resolutions (a true refinement study) and refining x and z together,
TM mode converges cleanly through this module's own mesh
construction, changing by under 1-2% between realistic production
resolution and twice that resolution. A real, independent bug was
also found and fixed alongside this: ``pycsamt.forward.em2d``'s
TM-mode assembly (``_assemble_tm``) combined interface resistivities
with a plain arithmetic mean, inconsistent with the parallel-current,
harmonic-mean-equivalent combination its own surface-impedance
extraction (``_surface_impedance_tm``) already used; it now uses the
same thickness-weighted harmonic mean in both places, which changes
nothing for laterally or vertically uniform structure (averaging two
equal values is unaffected by the choice of mean) and improves
accuracy for genuinely heterogeneous 2-D structure.

This module only produces clean, noiseless forward-consistent survey
responses.  Realistic noise, dropout, static shift, and distortion
are a separate, already-built stage: see
:mod:`pycsamt.ai.domain_gap.simulator`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ...forward.maxwell import (
    MaxwellMesh,
    MaxwellProblem,
    MaxwellResultCache,
    ReceiverSet,
    skin_depth_m,
)
from ...forward.maxwell.batch import BatchPolicy, solve_batch
from ...forward.maxwell.mt2d import MT2DAdapter
from ..data.contracts import SurveyData
from ..data.manifest import DatasetManifest
from ..data.splits import RealizationSplit, split_realizations
from ..experiments.config import SeedPlan
from ..geology import GaussianCorrelation, GeologyGrid, generate_gaussian_field

__all__ = [
    "Maxwell2DSample",
    "Maxwell2DDataset",
    "Maxwell2DDatasetConfig",
    "generate_2d_maxwell_dataset",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    """Return a read-only copy of *value* as an ndarray."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


_Z_GROWTH = 1.22
"""Geometric z growth rate, matching the mesh already validated in
``forward/tests/test_maxwell_mt2d.py``.
"""


def _graded_z_edges(first_cell_m: float, target_m: float) -> np.ndarray:
    """Return z edges as one continuous geometric progression from
    zero, growing by :data:`_Z_GROWTH` per cell, reaching at least
    *target_m*.
    """
    widths = [first_cell_m]
    total = first_cell_m
    while total < target_m and len(widths) < 500:
        widths.append(widths[-1] * _Z_GROWTH)
        total += widths[-1]
    return np.concatenate([[0.0], np.cumsum(widths)])


def _uniform_x_edges(
    centre_m: float, cell_width_m: float, half_span_m: float
) -> np.ndarray:
    """Return one continuous *uniform* x progression centred on
    *centre_m*, reaching at least *half_span_m* on each side.
    """
    n = int(np.ceil(half_span_m / cell_width_m))
    half = np.arange(n + 1) * cell_width_m
    return centre_m + np.concatenate([-half[::-1][:-1], half])


def _centres_to_edges(centres: np.ndarray) -> np.ndarray:
    """Return edges from strictly increasing centres via midpoints,
    extrapolated at both ends.
    """
    middle = (centres[:-1] + centres[1:]) / 2.0
    first = centres[0] - (middle[0] - centres[0])
    last = centres[-1] + (centres[-1] - middle[-1])
    return np.concatenate([[first], middle, [last]])


def _resample_2d(
    values: np.ndarray,
    old_z_edges: np.ndarray,
    old_x_edges: np.ndarray,
    new_z_centres: np.ndarray,
    new_x_centres: np.ndarray,
) -> np.ndarray:
    """Piecewise-constant resample of a ``(z, x)`` model onto new
    cell centres; positions outside the original extent take the
    nearest edge row/column's value.
    """
    z_idx = np.searchsorted(old_z_edges, new_z_centres, side="right") - 1
    z_idx = np.clip(z_idx, 0, values.shape[0] - 1)
    x_idx = np.searchsorted(old_x_edges, new_x_centres, side="right") - 1
    x_idx = np.clip(x_idx, 0, values.shape[1] - 1)
    return values[np.ix_(z_idx, x_idx)]


def _solver_mesh_and_conductivity(
    grid: GeologyGrid,
    resistivity_ohm_m: np.ndarray,
    frequencies_hz: np.ndarray,
    *,
    safety_factor: float,
    max_mesh_cells: int,
) -> tuple[MaxwellMesh, np.ndarray]:
    """Build a 2-D solver mesh and matching conductivity, TE-mode
    validated.

    Uses a single continuous geometric z-progression and a single
    continuous *uniform* x-progression centred on the grid, rather
    than a small uniform core plus separately geometrically-grown
    padding. See the module docstring for why, and for the
    unresolved TM-mode (``zyx``) accuracy problem this construction
    does not fix.

    No air layer is added: ``MT2DAdapter`` always treats the full
    mesh as earth (its own module docstring).

    Raises
    ------
    ValueError
        If the required lateral extent at the grid's own resolution
        would need more than *max_mesh_cells* uniform x-cells.
    """
    dz, dx = grid.spacing_m
    shallow_skin_depth = float(
        skin_depth_m(
            float(np.min(resistivity_ohm_m)), float(np.max(frequencies_hz))
        )
    )
    deep_skin_depth = float(
        skin_depth_m(
            float(np.max(resistivity_ohm_m)), float(np.min(frequencies_hz))
        )
    )
    required = safety_factor * deep_skin_depth
    z_first_cell = min(dz, shallow_skin_depth / 10.0)

    z_target = max(required, grid.extent_m["z"][1])
    z_edges = _graded_z_edges(z_first_cell, z_target)
    z_centres = (z_edges[:-1] + z_edges[1:]) / 2.0

    x_min, x_max = grid.extent_m["x"]
    x_centre = (x_min + x_max) / 2.0
    x_half_target = max(required, (x_max - x_min) / 2.0)
    n_x_per_side = int(np.ceil(x_half_target / dx))
    if 2 * n_x_per_side > max_mesh_cells:
        raise ValueError(
            "requested configuration needs "
            f"{2 * n_x_per_side} uniform x-cells (lateral extent "
            f"{2 * x_half_target:.0f} m at {dx:.0f} m resolution), "
            f"exceeding max_mesh_cells={max_mesh_cells}. MT2DAdapter "
            "needs uniform lateral spacing (a graded x-mesh was "
            "found empirically to corrupt its TM-mode response), so "
            "cell count scales linearly with domain width. Reduce "
            "mesh_safety_factor, raise the survey frequency range, "
            "lower log_resistivity_mean/log_resistivity_std, coarsen "
            "the grid's dx_m, or raise max_mesh_cells if the runtime "
            "cost is acceptable."
        )
    x_edges = _uniform_x_edges(x_centre, dx, x_half_target)
    x_centres = (x_edges[:-1] + x_edges[1:]) / 2.0

    old_z_edges = _centres_to_edges(grid.z_m)
    old_x_edges = _centres_to_edges(grid.x_m)
    resistivity_on_mesh = _resample_2d(
        resistivity_ohm_m, old_z_edges, old_x_edges, z_centres, x_centres
    )
    return MaxwellMesh(x_edges, z_edges), 1.0 / resistivity_on_mesh


def _sample_correlation(
    config: Maxwell2DDatasetConfig, seed: int
) -> GaussianCorrelation:
    """Draw one seeded :class:`GaussianCorrelation` within config
    ranges.
    """
    rng = np.random.default_rng(seed)
    length_x = rng.uniform(*config.correlation_length_x_m)
    length_z = rng.uniform(*config.correlation_length_z_m)
    return GaussianCorrelation(length_x_m=length_x, length_z_m=length_z)


@dataclass(frozen=True)
class Maxwell2DDatasetConfig:
    """Immutable configuration for :func:`generate_2d_maxwell_dataset`.

    Parameters
    ----------
    dataset_id : str
        Portable dataset identifier.
    grid : GeologyGrid
        Shared 2-D ``(z, x)`` geometry for every realization. Its
        shallowest edge must be at depth zero.
    correlation_length_x_m, correlation_length_z_m : (float, float)
        Inclusive ``(low, high)`` ranges a horizontal/vertical
        correlation length is drawn from per realization.
    frequencies_hz : array-like
        Shared positive, unique simulated frequencies.
    station_x_m : array-like
        Shared receiver x-positions, within ``grid``'s x extent.
    n_realizations : int
        Number of geological realizations to attempt.
    seed : int
        Root seed; every realization derives a labeled child seed
        from it via :class:`~pycsamt.ai.experiments.config.SeedPlan`.
    log_resistivity_mean, log_resistivity_std : float
        Affine map from the standardized correlated field to
        ``log10(resistivity_ohm_m)``.
    components : sequence of {"zxy", "zyx"}, default=("zxy", "zyx")
        Requested 2-D impedance components. Both modes are validated
        by this module; see the module docstring for the convergence
        evidence behind that claim.
    mesh_safety_factor : float, default=8.0
        Skin-depth safety factor for the solver mesh's extent; see
        :func:`_solver_mesh_and_conductivity`.
    max_mesh_cells : int, default=200_000
        Upper bound on uniform lateral solver-mesh cells (both
        sides combined). Raises rather than silently attempting an
        intractable solve for configurations needing both fine
        lateral resolution and a wide domain.
    validation_fraction, test_fraction : float, default=0.1
        Forwarded to
        :func:`~pycsamt.ai.data.splits.split_realizations`.
    namespace : str, default="maxwell2d"
        :class:`SeedPlan` namespace.

    Examples
    --------
    >>> from pycsamt.ai.geology import GeologyGrid
    >>> grid = GeologyGrid.regular_2d(nx=8, nz=6, dx_m=200, dz_m=100)
    >>> config = Maxwell2DDatasetConfig(
    ...     dataset_id="demo-2d-v1",
    ...     grid=grid,
    ...     correlation_length_x_m=(400.0, 800.0),
    ...     correlation_length_z_m=(100.0, 200.0),
    ...     frequencies_hz=[10.0, 1.0],
    ...     station_x_m=[500.0, 900.0, 1300.0],
    ...     n_realizations=2,
    ...     seed=0,
    ... )
    >>> config.components
    ('zxy', 'zyx')
    """

    dataset_id: str
    grid: GeologyGrid
    correlation_length_x_m: tuple[float, float]
    correlation_length_z_m: tuple[float, float]
    frequencies_hz: Any
    station_x_m: Any
    n_realizations: int
    seed: int
    log_resistivity_mean: float = 2.0
    log_resistivity_std: float = 0.4
    components: tuple[str, ...] = ("zxy", "zyx")
    mesh_safety_factor: float = 8.0
    max_mesh_cells: int = 200_000
    validation_fraction: float = 0.1
    test_fraction: float = 0.1
    namespace: str = "maxwell2d"
    verbose: bool | int | str = False

    def __post_init__(self) -> None:
        if not isinstance(self.grid, GeologyGrid) or self.grid.dimension != 2:
            raise TypeError("grid must be a 2-D GeologyGrid.")
        shallow_edge = self.grid.extent_m["z"][0]
        if not np.isclose(shallow_edge, 0.0, atol=1e-6):
            raise ValueError(
                "grid's shallowest edge must be at depth 0 for "
                f"receivers placed at the surface (got {shallow_edge} m)."
            )
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
        x_min, x_max = self.grid.extent_m["x"]
        if (
            station_x.ndim != 1
            or station_x.size < 1
            or not np.all(np.isfinite(station_x))
            or np.any(station_x < x_min)
            or np.any(station_x > x_max)
            or len(set(station_x.tolist())) != station_x.size
        ):
            raise ValueError(
                "station_x_m must be a non-empty, finite, unique "
                f"vector within the grid's x extent {(x_min, x_max)}."
            )
        components = tuple(str(c).strip().lower() for c in self.components)
        if not components or any(c not in ("zxy", "zyx") for c in components):
            raise ValueError(
                "components must be a non-empty subset of "
                "('zxy', 'zyx') for a 2-D problem."
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
        if not np.isfinite(self.mesh_safety_factor) or (
            self.mesh_safety_factor <= 0
        ):
            raise ValueError("mesh_safety_factor must be finite and positive.")
        if (
            not isinstance(self.max_mesh_cells, int)
            or isinstance(self.max_mesh_cells, bool)
            or self.max_mesh_cells < 2
        ):
            raise ValueError(
                "max_mesh_cells must be an integer of at least 2."
            )
        namespace = str(self.namespace).strip()
        if not namespace:
            raise ValueError("namespace cannot be empty.")
        object.__setattr__(self, "frequencies_hz", _readonly(frequencies))
        object.__setattr__(self, "station_x_m", _readonly(station_x))
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "namespace", namespace)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable configuration record.

        Returns
        -------
        dict
            Every field of this configuration.

        Examples
        --------
        >>> from pycsamt.ai.geology import GeologyGrid
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=1, dz_m=1)
        >>> config = Maxwell2DDatasetConfig(
        ...     "d",
        ...     grid,
        ...     (1.0, 2.0),
        ...     (1.0, 2.0),
        ...     [1.0],
        ...     [1.0, 2.0],
        ...     1,
        ...     0,
        ... )
        >>> config.to_dict()["dataset_id"]
        'd'
        """
        return {
            "schema_version": 1,
            "dataset_id": self.dataset_id,
            "grid": self.grid.to_dict(),
            "correlation_length_x_m": list(self.correlation_length_x_m),
            "correlation_length_z_m": list(self.correlation_length_z_m),
            "frequencies_hz": self.frequencies_hz.tolist(),
            "station_x_m": self.station_x_m.tolist(),
            "n_realizations": self.n_realizations,
            "seed": self.seed,
            "log_resistivity_mean": self.log_resistivity_mean,
            "log_resistivity_std": self.log_resistivity_std,
            "components": list(self.components),
            "mesh_safety_factor": self.mesh_safety_factor,
            "max_mesh_cells": self.max_mesh_cells,
            "validation_fraction": self.validation_fraction,
            "test_fraction": self.test_fraction,
            "namespace": self.namespace,
        }


@dataclass(frozen=True)
class Maxwell2DSample:
    """One converged 2-D realization: model, response, and provenance.

    Parameters
    ----------
    realization_id : str
        Unique identifier, stable across regeneration with the same
        configuration.
    seed : int
        Seed used to generate this realization's correlated field.
    correlation : GaussianCorrelation
        Correlation lengths drawn for this realization.
    resistivity_ohm_m : ndarray
        True resistivity on the shared ``grid``, the training target.
    survey : SurveyData
        Simulated, noiseless response at the configured stations and
        frequencies, the training input.
    mesh_cells : int
        Total padded solver-mesh cell count.
    relative_residual : float
        Worst-case relative solver residual over every frequency.
    """

    realization_id: str
    seed: int
    correlation: GaussianCorrelation
    resistivity_ohm_m: np.ndarray
    survey: SurveyData
    mesh_cells: int
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
        object.__setattr__(self, "mesh_cells", int(self.mesh_cells))
        object.__setattr__(
            self, "relative_residual", float(self.relative_residual)
        )


@dataclass(frozen=True)
class Maxwell2DDataset:
    """Immutable collection of converged 2-D realizations with a split.

    Parameters
    ----------
    grid : GeologyGrid
        Shared geometry of every sample's ``resistivity_ohm_m``.
    samples : tuple of Maxwell2DSample
        Converged realizations, at least one.
    split : RealizationSplit
        Realization-level train/validation/test assignment covering
        every sample.
    manifest : DatasetManifest
        Reproducibility record for this dataset.
    rejected : tuple of str
        Realization IDs excluded for not converging.
    """

    grid: GeologyGrid
    samples: tuple[Maxwell2DSample, ...]
    split: RealizationSplit
    manifest: DatasetManifest
    rejected: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.grid, GeologyGrid):
            raise TypeError("grid must be a GeologyGrid.")
        samples = tuple(self.samples)
        if not samples:
            raise ValueError("samples must contain at least one entry.")
        if any(not isinstance(s, Maxwell2DSample) for s in samples):
            raise TypeError("samples entries must be Maxwell2DSample.")
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

    def select(self, partition: str) -> tuple[Maxwell2DSample, ...]:
        """Return samples belonging to one split partition.

        Parameters
        ----------
        partition : {"train", "validation", "test"}
            Partition name.

        Returns
        -------
        tuple of Maxwell2DSample
            Samples whose realization ID falls in that partition.

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


def generate_2d_maxwell_dataset(
    config: Maxwell2DDatasetConfig,
    *,
    cache: MaxwellResultCache | None = None,
    policy: BatchPolicy | None = None,
) -> Maxwell2DDataset:
    """Generate a cached, split, versioned 2-D Maxwell dataset.

    Parameters
    ----------
    config : Maxwell2DDatasetConfig
        Complete, immutable dataset configuration.
    cache : MaxwellResultCache or None, optional
        When given, forward solves are cached by problem hash and
        resumable across calls.
    policy : BatchPolicy or None, optional
        Retry/concurrency policy forwarded to
        :func:`~pycsamt.forward.maxwell.batch.solve_batch`.

    Returns
    -------
    Maxwell2DDataset
        Every realization whose forward solve converged, split at
        the realization level.

    Raises
    ------
    ValueError
        If no realization converges.

    Examples
    --------
    >>> from pycsamt.ai.geology import GeologyGrid
    >>> grid = GeologyGrid.regular_2d(nx=6, nz=4, dx_m=300, dz_m=150)
    >>> config = Maxwell2DDatasetConfig(
    ...     dataset_id="demo-2d-v1",
    ...     grid=grid,
    ...     correlation_length_x_m=(600.0, 900.0),
    ...     correlation_length_z_m=(150.0, 300.0),
    ...     frequencies_hz=[10.0, 3.0],
    ...     station_x_m=[600.0, 900.0, 1200.0],
    ...     n_realizations=2,
    ...     seed=0,
    ...     validation_fraction=0.0,
    ...     test_fraction=0.0,
    ... )
    >>> dataset = generate_2d_maxwell_dataset(config)
    >>> dataset.samples[0].survey.shape
    (3, 2, 2)
    >>> dataset.rejected
    ()
    """
    from pycsamt.api.view.progress import get_progress_bar

    plan = SeedPlan(config.seed, namespace=config.namespace)
    grid = config.grid
    station_names = tuple(f"S{i:04d}" for i in range(len(config.station_x_m)))
    receivers = ReceiverSet(
        [[float(x), 0.0] for x in config.station_x_m], station_names
    )
    coordinates_m = [[float(x), 0.0, 0.0] for x in config.station_x_m]

    realizations: dict[str, dict[str, Any]] = {}
    problems: dict[str, MaxwellProblem] = {}
    with get_progress_bar(
        total=config.n_realizations,
        desc="Generating 2-D realizations",
        unit="realization",
        verbose=config.verbose,
    ) as bar:
        for index in range(config.n_realizations):
            realization_id = f"{config.dataset_id}-r{index:05d}"
            correlation_seed = plan.derive(f"{realization_id}/correlation")
            field_seed = plan.derive(f"{realization_id}/field")
            correlation = _sample_correlation(config, correlation_seed)
            field = generate_gaussian_field(grid, correlation, seed=field_seed)
            log_resistivity = (
                config.log_resistivity_mean
                + config.log_resistivity_std * field.values
            )
            resistivity = np.power(10.0, log_resistivity)
            mesh, conductivity = _solver_mesh_and_conductivity(
                grid,
                resistivity,
                config.frequencies_hz,
                safety_factor=config.mesh_safety_factor,
                max_mesh_cells=config.max_mesh_cells,
            )
            problem = MaxwellProblem(
                mesh,
                conductivity,
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

    samples: list[Maxwell2DSample] = []
    rejected: list[str] = []

    with get_progress_bar(
        total=len(problems),
        desc="Solving 2-D Maxwell realizations",
        unit="problem",
        verbose=config.verbose,
    ) as solve_bar:

        def _on_result(problem: MaxwellProblem, result: Any) -> None:
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
                Maxwell2DSample(
                    realization_id=meta["realization_id"],
                    seed=meta["seed"],
                    correlation=meta["correlation"],
                    resistivity_ohm_m=meta["resistivity_ohm_m"],
                    survey=survey,
                    mesh_cells=int(np.prod(problem.mesh.shape)),
                    relative_residual=(
                        result.diagnostics.maximum_relative_residual
                    ),
                )
            )
            solve_bar.update(1)

        def _on_failure(problem: MaxwellProblem, failure: Any) -> None:
            rejected.append(
                realizations[problem.problem_hash]["realization_id"]
            )
            solve_bar.update(1)

        solve_batch(
            problems.values(),
            MT2DAdapter(),
            cache=cache,
            policy=policy,
            on_result=_on_result,
            on_failure=_on_failure,
        )

    if not samples:
        raise ValueError(
            "no realization produced a converged forward solve; "
            "check mesh_safety_factor and the requested frequency/"
            "resistivity range."
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
            "pycsamt.ai.training.dataset2d.generate_2d_maxwell_dataset"
        ),
        generator_version="0.1.0",
        configuration=config.to_dict(),
        split=split,
        sample_count=len(samples),
    )
    return Maxwell2DDataset(
        grid=grid,
        samples=tuple(samples),
        split=split,
        manifest=manifest,
        rejected=tuple(rejected),
    )
