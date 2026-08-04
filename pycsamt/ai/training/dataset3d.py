# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Generate versioned 3-D Maxwell training datasets (M8 groundwork).

This module replaces the tiled-1-D-plus-GCN-smoothing approximation
used by the current ``Inv3DAgent`` with genuine 3-D physics, mirroring
:mod:`pycsamt.ai.training.dataset2d`'s M4/M5 role but for M8:
correlated 3-D resistivity volumes from :mod:`pycsamt.ai.geology` are
solved with the research-only, small-grid
:class:`~pycsamt.forward.maxwell.mt3d.MT3DAdapter`, cached and retried
through :func:`~pycsamt.forward.maxwell.batch.solve_batch`, and
packaged into a :class:`Maxwell3DDataset` with a realization-level
train/validation/test split and a ``DatasetManifest``
(:mod:`pycsamt.ai.data.manifest`).

A realization whose solve does not converge (see
:attr:`~pycsamt.forward.maxwell.contracts.ForwardResult.success`) is
excluded from the dataset rather than silently included, per the
plan's requirement that failed solves cannot silently enter training
data.

Solver mesh construction: unlike :mod:`pycsamt.ai.training.dataset2d`
(which affords a large uniform-lateral mesh because 2-D cell count
scales linearly with domain width), :class:`MT3DAdapter` is capped at
``max_mesh_cells`` (small-grid research-only, per
``docs/source/development/adr/AI-INVERSION-M6-3D-ADR.md``) and 3-D
cell count scales as the *cube* of resolution, so a uniform mesh
cannot afford both fine resolution near the structure/stations and a
skin-depth-scale domain extent at the same time. This module instead
builds a **padded, non-uniform** mesh — a core region at the grid's
own native resolution (an exact, no-resample match to the geological
model there), geometrically padded outward in x/y and downward in z —
the same strategy :func:`~pycsamt.forward.maxwell.mt3d`'s own
calibrated benchmark mesh uses (see its module docstring and
``pycsamt/forward/tests/test_maxwell_mt3d.py``). ``mt3d.py`` only
gained non-uniform mesh support for this reason; see
:mod:`pycsamt.forward.maxwell.mt3d`'s module docstring for the
diagnosis this fixed.

This module only produces clean, noiseless forward-consistent survey
responses. Realistic noise, dropout, static shift, and distortion are
a separate, already-built stage: see
:mod:`pycsamt.ai.domain_gap.simulator`.

Frequency-dependent mesh accuracy (mitigated, not fully closed)
------------------------------------------------------------------
:class:`MT3DAdapter`'s own calibrated benchmark mesh
(``forward/tests/test_maxwell_mt3d.py``) is validated at <=2 Hz. An
earlier version of this module sized the solver's core resolution
from the geological grid's own spacing only (a *domain-extent* safety
factor, not a *resolution* one), which gave ~1-2% half-space error at
1-2 Hz but ~10% at 20-50 Hz -- a smaller skin depth needs finer
resolution than a fixed grid spacing affords, the same category of
finding :mod:`~pycsamt.forward.maxwell.mt2d`'s own mesh calibration
required (``test_maxwell_mt2d.py``'s module docstring).
``cells_per_skin_depth`` (:class:`Maxwell3DDatasetConfig`, default
8.0) mitigates this: the solver's core resolution is now chosen
independently of the geological grid's own spacing, targeting that
many cells across the shallowest skin depth (capped at the grid's own
spacing -- refining finer than the grid's own resolution adds solver
cost without adding information). At the default 8.0, the same
50/20 Hz half-space check improves from ~10.5%/6.3% error to
~8.8%/5.2% at an unchanged cell budget; raising ``cells_per_skin_depth``
(paired with ``max_mesh_cells``, since finer core cells cost more)
narrows it further (~6%/3.6% at ``cells_per_skin_depth=12`` and
roughly 10,000 cells) but has **not** been driven under 5% at 50 Hz
within a small-grid research budget -- this is a genuine cost/accuracy
trade-off of a fixed, non-adaptive cell count, not a bug still to fix.
Callers requesting frequencies much above a few hertz should
independently check accuracy for their own configuration (e.g. via a
low-``log_resistivity_std`` realization compared against
:func:`~pycsamt.forward.maxwell.benchmarks.half_space_impedance`, as
``pycsamt/ai/tests/test_ai_training_dataset3d.py`` does) rather than
assume the default holds.
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
from ...forward.maxwell.mt3d import MT3DAdapter
from ..data.contracts import SurveyData
from ..data.manifest import DatasetManifest
from ..data.splits import RealizationSplit, split_realizations
from ..experiments.config import SeedPlan
from ..geology import GaussianCorrelation, GeologyGrid, generate_gaussian_field

__all__ = [
    "Maxwell3DSample",
    "Maxwell3DDataset",
    "Maxwell3DDatasetConfig",
    "generate_3d_maxwell_dataset",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    """Return a read-only copy of *value* as an ndarray."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


_XY_GROWTH = 1.8
_Z_GROWTH = 1.6
"""Geometric padding growth rates, matching the mesh already validated
in ``forward/tests/test_maxwell_mt3d.py``'s calibrated benchmark mesh.
"""

_MAX_PAD_CELLS_PER_SIDE = 14
"""Safety ceiling on padding cells per side, independent of
``max_mesh_cells``: stops a pathological configuration (e.g. an
enormous safety factor) from looping indefinitely before the total
cell-count check below can even run.
"""


def _pad_widths(
    first_width: float, target_extra_m: float, growth: float
) -> np.ndarray:
    """Return geometrically growing pad widths reaching *target_extra_m*."""
    if target_extra_m <= 0:
        return np.zeros(0)
    widths = []
    total = 0.0
    width = first_width
    while total < target_extra_m and len(widths) < _MAX_PAD_CELLS_PER_SIDE:
        width *= growth
        widths.append(width)
        total += width
    return np.asarray(widths)


def _padded_symmetric_edges(
    core_min: float,
    core_max: float,
    core_width: float,
    half_extent_needed: float,
) -> np.ndarray:
    """Core edges at native resolution, padded symmetrically outward."""
    n_core = int(round((core_max - core_min) / core_width))
    core_edges = core_min + np.arange(n_core + 1) * core_width
    extra = max(0.0, half_extent_needed - (core_max - core_min) / 2.0)
    pad = _pad_widths(core_width, extra, _XY_GROWTH)
    if pad.size == 0:
        return core_edges
    left = core_min - np.cumsum(pad[::-1])[::-1]
    right = core_max + np.cumsum(pad)
    return np.concatenate([left, core_edges, right])


def _padded_z_edges(
    core_max: float, core_width: float, depth_needed: float
) -> np.ndarray:
    """Core edges from depth zero, padded downward only."""
    n_core = int(round(core_max / core_width))
    core_edges = np.arange(n_core + 1) * core_width
    extra = max(0.0, depth_needed - core_max)
    pad = _pad_widths(core_width, extra, _Z_GROWTH)
    if pad.size == 0:
        return core_edges
    return np.concatenate([core_edges, core_max + np.cumsum(pad)])


def _centres_to_edges(centres: np.ndarray) -> np.ndarray:
    """Return edges from strictly increasing centres via midpoints,
    extrapolated at both ends.
    """
    middle = (centres[:-1] + centres[1:]) / 2.0
    first = centres[0] - (middle[0] - centres[0])
    last = centres[-1] + (centres[-1] - middle[-1])
    return np.concatenate([[first], middle, [last]])


def _resample_3d(
    values: np.ndarray,
    old_z_edges: np.ndarray,
    old_y_edges: np.ndarray,
    old_x_edges: np.ndarray,
    new_z_centres: np.ndarray,
    new_y_centres: np.ndarray,
    new_x_centres: np.ndarray,
) -> np.ndarray:
    """Piecewise-constant resample of a ``(z, y, x)`` model onto new
    cell centres; positions outside the original extent take the
    nearest edge cell's value (the padding region's convention).
    """
    z_idx = np.clip(
        np.searchsorted(old_z_edges, new_z_centres, side="right") - 1,
        0,
        values.shape[0] - 1,
    )
    y_idx = np.clip(
        np.searchsorted(old_y_edges, new_y_centres, side="right") - 1,
        0,
        values.shape[1] - 1,
    )
    x_idx = np.clip(
        np.searchsorted(old_x_edges, new_x_centres, side="right") - 1,
        0,
        values.shape[2] - 1,
    )
    return values[np.ix_(z_idx, y_idx, x_idx)]


def _solver_mesh_and_conductivity(
    grid: GeologyGrid,
    resistivity_ohm_m: np.ndarray,
    station_xy_m: np.ndarray,
    frequencies_hz: np.ndarray,
    *,
    safety_factor: float,
    max_mesh_cells: int,
    cells_per_skin_depth: float | None,
) -> tuple[MaxwellMesh, np.ndarray]:
    """Build a padded, non-uniform 3-D solver mesh and matching
    conductivity for :class:`~pycsamt.forward.maxwell.mt3d.MT3DAdapter`.

    When *cells_per_skin_depth* is given, the solver's own core
    resolution is decoupled from the geological grid's resolution and
    chosen to resolve that many cells across the shallowest
    (highest-frequency, lowest-resistivity) skin depth, capped at the
    grid's own spacing (refining below the grid's own resolution would
    not add information, only cost). This is what makes the mesh
    frequency-aware -- see the module docstring's "Frequency-dependent
    mesh accuracy" section for the half-space evidence that motivated
    it, and why it defaults to ``None`` (the grid's own native
    resolution, this function's original behavior) rather than being
    on by default: refining core resolution can multiply cell count by
    an order of magnitude for a config that previously fit comfortably
    within ``max_mesh_cells``, so it must be an explicit choice.
    Whichever core resolution is chosen, the (possibly coarser)
    geological model is resampled onto it, the same nearest-cell
    convention already used for the padding region.

    Raises
    ------
    ValueError
        If the required extent at the chosen core resolution would
        need more than *max_mesh_cells* total cells.
    """
    dz, dy, dx = grid.spacing_m
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
    if cells_per_skin_depth is None:
        core_dx, core_dy, core_dz = dx, dy, dz
    else:
        target_core = shallow_skin_depth / cells_per_skin_depth
        core_dx = min(dx, target_core)
        core_dy = min(dy, target_core)
        core_dz = min(dz, target_core)

    x_min, x_max = grid.extent_m["x"]
    y_min, y_max = grid.extent_m["y"]
    station_half_x = float(
        max(
            np.max(station_xy_m[:, 0]) - x_min,
            x_max - np.min(station_xy_m[:, 0]),
        )
    )
    station_half_y = float(
        max(
            np.max(station_xy_m[:, 1]) - y_min,
            y_max - np.min(station_xy_m[:, 1]),
        )
    )
    x_half_needed = max(required, station_half_x, (x_max - x_min) / 2.0)
    y_half_needed = max(required, station_half_y, (y_max - y_min) / 2.0)

    x_edges = _padded_symmetric_edges(x_min, x_max, core_dx, x_half_needed)
    y_edges = _padded_symmetric_edges(y_min, y_max, core_dy, y_half_needed)
    z_max = grid.extent_m["z"][1]
    z_edges = _padded_z_edges(
        z_max, core_dz, max(required, shallow_skin_depth)
    )

    total_cells = (len(x_edges) - 1) * (len(y_edges) - 1) * (len(z_edges) - 1)
    if total_cells > max_mesh_cells:
        core_note = (
            f"{core_dx:.0f}x{core_dy:.0f}x{core_dz:.0f} m frequency-aware "
            f"core resolution (cells_per_skin_depth={cells_per_skin_depth}), "
            f"vs. {dx:.0f}x{dy:.0f}x{dz:.0f} m grid resolution"
            if cells_per_skin_depth is not None
            else f"{dx:.0f}x{dy:.0f}x{dz:.0f} m native grid resolution"
        )
        knob_note = (
            "lower cells_per_skin_depth, "
            if cells_per_skin_depth is not None
            else ""
        )
        raise ValueError(
            f"requested configuration needs {total_cells} solver cells "
            f"(padded extent {x_edges[-1] - x_edges[0]:.0f} x "
            f"{y_edges[-1] - y_edges[0]:.0f} x {z_edges[-1]:.0f} m at "
            f"{core_note}), exceeding max_mesh_cells={max_mesh_cells}. "
            "MT3DAdapter is a small-grid research-only solver (see its "
            f"module docstring). Reduce mesh_safety_factor, {knob_note}"
            "raise the survey frequency range, lower "
            "log_resistivity_mean/log_resistivity_std, coarsen the "
            "grid's spacing, shrink the station network's spatial "
            "extent, or raise max_mesh_cells if the runtime cost is "
            "acceptable."
        )

    old_z_edges = _centres_to_edges(grid.z_m)
    old_y_edges = _centres_to_edges(grid.y_m)
    old_x_edges = _centres_to_edges(grid.x_m)
    new_z_centres = (z_edges[:-1] + z_edges[1:]) / 2.0
    new_y_centres = (y_edges[:-1] + y_edges[1:]) / 2.0
    new_x_centres = (x_edges[:-1] + x_edges[1:]) / 2.0
    resistivity_on_mesh = _resample_3d(
        resistivity_ohm_m,
        old_z_edges,
        old_y_edges,
        old_x_edges,
        new_z_centres,
        new_y_centres,
        new_x_centres,
    )
    return MaxwellMesh(x_edges, z_edges, y_edges), 1.0 / resistivity_on_mesh


def _sample_correlation(
    config: Maxwell3DDatasetConfig, seed: int
) -> GaussianCorrelation:
    """Draw one seeded :class:`GaussianCorrelation` within config ranges."""
    rng = np.random.default_rng(seed)
    length_x = rng.uniform(*config.correlation_length_x_m)
    length_y = rng.uniform(*config.correlation_length_y_m)
    length_z = rng.uniform(*config.correlation_length_z_m)
    return GaussianCorrelation(
        length_x_m=length_x, length_z_m=length_z, length_y_m=length_y
    )


@dataclass(frozen=True)
class Maxwell3DDatasetConfig:
    """Immutable configuration for :func:`generate_3d_maxwell_dataset`.

    Parameters
    ----------
    dataset_id : str
        Portable dataset identifier.
    grid : GeologyGrid
        Shared 3-D ``(z, y, x)`` geometry for every realization. Its
        shallowest edge must be at depth zero.
    correlation_length_x_m, correlation_length_y_m, correlation_length_z_m
        : (float, float)
        Inclusive ``(low, high)`` ranges a horizontal/horizontal/
        vertical correlation length is drawn from per realization.
    frequencies_hz : array-like
        Shared positive, unique simulated frequencies.
    station_xy_m : array-like, shape (n_stations, 2)
        Shared receiver ``(x, y)`` positions, within ``grid``'s x/y
        extent.
    n_realizations : int
        Number of geological realizations to attempt.
    seed : int
        Root seed; every realization derives a labeled child seed
        from it via :class:`~pycsamt.ai.experiments.config.SeedPlan`.
    log_resistivity_mean, log_resistivity_std : float
        Affine map from the standardized correlated field to
        ``log10(resistivity_ohm_m)``.
    components : sequence of {"zxx", "zxy", "zyx", "zyy"},
        default=("zxy", "zyx")
        Requested 3-D impedance components.
    mesh_safety_factor : float, default=3.0
        Skin-depth safety factor for the solver mesh's padded extent;
        see :func:`_solver_mesh_and_conductivity`. Deliberately
        smaller than :mod:`~pycsamt.ai.training.dataset2d`'s default
        (8.0): :class:`MT3DAdapter`'s ``max_mesh_cells`` budget is
        far tighter, and a padded (rather than uniform) mesh already
        keeps far-field cells cheap.
    max_mesh_cells : int, default=6000
        Upper bound on total solver-mesh cells, matching
        :class:`~pycsamt.forward.maxwell.mt3d.MT3DAdapter`'s own
        default ``max_cells``; the same value is used to construct
        the adapter this module solves with, so the two stay
        consistent automatically. Raising this only helps if you have
        confirmed the resulting solve time is acceptable.
    cells_per_skin_depth : float or None, default=None
        Target solver-mesh core resolution, expressed as cells across
        the shallowest (highest-frequency, lowest-resistivity) skin
        depth, capped at the geological grid's own spacing (refining
        below the grid's own resolution costs cells without adding
        information). This is what makes the mesh frequency-aware; see
        the module docstring's "Frequency-dependent mesh accuracy"
        section for the half-space evidence a value of 8.0 gives
        (~1-2% error at 1-2 Hz, ~5-9% at 20-50 Hz, vs. ~10% at 20-50 Hz
        with the default). ``None`` (the default) preserves this
        function's original behavior -- core resolution equals the
        geological grid's own spacing exactly -- because refining it
        can multiply solver cell count (and therefore
        ``max_mesh_cells`` requirements) by an order of magnitude for
        a configuration that previously fit comfortably; this must be
        an explicit opt-in, not a silent default-behavior change.
    validation_fraction, test_fraction : float, default=0.1
        Forwarded to
        :func:`~pycsamt.ai.data.splits.split_realizations`.
    namespace : str, default="maxwell3d"
        :class:`SeedPlan` namespace.

    Examples
    --------
    >>> from pycsamt.ai.geology import GeologyGrid
    >>> grid = GeologyGrid.regular_3d(
    ...     nx=4, ny=4, nz=6, dx_m=200, dy_m=200, dz_m=100
    ... )
    >>> config = Maxwell3DDatasetConfig(
    ...     dataset_id="demo-3d-v1",
    ...     grid=grid,
    ...     correlation_length_x_m=(400.0, 800.0),
    ...     correlation_length_y_m=(400.0, 800.0),
    ...     correlation_length_z_m=(100.0, 200.0),
    ...     frequencies_hz=[50.0, 20.0],
    ...     station_xy_m=[[400.0, 400.0], [600.0, 600.0]],
    ...     n_realizations=2,
    ...     seed=0,
    ... )
    >>> config.components
    ('zxy', 'zyx')
    """

    dataset_id: str
    grid: GeologyGrid
    correlation_length_x_m: tuple[float, float]
    correlation_length_y_m: tuple[float, float]
    correlation_length_z_m: tuple[float, float]
    frequencies_hz: Any
    station_xy_m: Any
    n_realizations: int
    seed: int
    log_resistivity_mean: float = 2.0
    log_resistivity_std: float = 0.4
    components: tuple[str, ...] = ("zxy", "zyx")
    mesh_safety_factor: float = 3.0
    max_mesh_cells: int = 6_000
    cells_per_skin_depth: float | None = None
    validation_fraction: float = 0.1
    test_fraction: float = 0.1
    namespace: str = "maxwell3d"
    verbose: bool | int | str = False

    def __post_init__(self) -> None:
        if not isinstance(self.grid, GeologyGrid) or self.grid.dimension != 3:
            raise TypeError("grid must be a 3-D GeologyGrid.")
        shallow_edge = self.grid.extent_m["z"][0]
        if not np.isclose(shallow_edge, 0.0, atol=1e-6):
            raise ValueError(
                "grid's shallowest edge must be at depth 0 for "
                f"receivers placed at the surface (got {shallow_edge} m)."
            )
        for name, value in (
            ("correlation_length_x_m", self.correlation_length_x_m),
            ("correlation_length_y_m", self.correlation_length_y_m),
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
        station_xy = np.asarray(self.station_xy_m, dtype=float)
        x_min, x_max = self.grid.extent_m["x"]
        y_min, y_max = self.grid.extent_m["y"]
        if (
            station_xy.ndim != 2
            or station_xy.shape[1] != 2
            or station_xy.shape[0] < 1
            or not np.all(np.isfinite(station_xy))
            or np.any(station_xy[:, 0] < x_min)
            or np.any(station_xy[:, 0] > x_max)
            or np.any(station_xy[:, 1] < y_min)
            or np.any(station_xy[:, 1] > y_max)
        ):
            raise ValueError(
                "station_xy_m must be a non-empty, finite (n, 2) array "
                f"within the grid's x extent {(x_min, x_max)} and y "
                f"extent {(y_min, y_max)}."
            )
        if (
            len({tuple(row) for row in station_xy.tolist()})
            != station_xy.shape[0]
        ):
            raise ValueError("station_xy_m rows must be unique.")
        components = tuple(str(c).strip().lower() for c in self.components)
        if not components or any(
            c not in ("zxx", "zxy", "zyx", "zyy") for c in components
        ):
            raise ValueError(
                "components must be a non-empty subset of "
                "('zxx', 'zxy', 'zyx', 'zyy') for a 3-D problem."
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
            or self.max_mesh_cells < 8
        ):
            raise ValueError(
                "max_mesh_cells must be an integer of at least 8."
            )
        if self.cells_per_skin_depth is not None and (
            not np.isfinite(self.cells_per_skin_depth)
            or self.cells_per_skin_depth <= 0
        ):
            raise ValueError(
                "cells_per_skin_depth must be None or finite and positive."
            )
        namespace = str(self.namespace).strip()
        if not namespace:
            raise ValueError("namespace cannot be empty.")
        object.__setattr__(self, "frequencies_hz", _readonly(frequencies))
        object.__setattr__(self, "station_xy_m", _readonly(station_xy))
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
        >>> grid = GeologyGrid.regular_3d(
        ...     nx=4, ny=4, nz=4, dx_m=1, dy_m=1, dz_m=1
        ... )
        >>> config = Maxwell3DDatasetConfig(
        ...     "d",
        ...     grid,
        ...     (1.0, 2.0),
        ...     (1.0, 2.0),
        ...     (1.0, 2.0),
        ...     [1.0],
        ...     [[1.0, 1.0]],
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
            "correlation_length_y_m": list(self.correlation_length_y_m),
            "correlation_length_z_m": list(self.correlation_length_z_m),
            "frequencies_hz": self.frequencies_hz.tolist(),
            "station_xy_m": self.station_xy_m.tolist(),
            "n_realizations": self.n_realizations,
            "seed": self.seed,
            "log_resistivity_mean": self.log_resistivity_mean,
            "log_resistivity_std": self.log_resistivity_std,
            "components": list(self.components),
            "mesh_safety_factor": self.mesh_safety_factor,
            "max_mesh_cells": self.max_mesh_cells,
            "cells_per_skin_depth": self.cells_per_skin_depth,
            "validation_fraction": self.validation_fraction,
            "test_fraction": self.test_fraction,
            "namespace": self.namespace,
        }


@dataclass(frozen=True)
class Maxwell3DSample:
    """One converged 3-D realization: model, response, and provenance.

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
class Maxwell3DDataset:
    """Immutable collection of converged 3-D realizations with a split.

    Parameters
    ----------
    grid : GeologyGrid
        Shared geometry of every sample's ``resistivity_ohm_m``.
    samples : tuple of Maxwell3DSample
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
    samples: tuple[Maxwell3DSample, ...]
    split: RealizationSplit
    manifest: DatasetManifest
    rejected: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.grid, GeologyGrid):
            raise TypeError("grid must be a GeologyGrid.")
        samples = tuple(self.samples)
        if not samples:
            raise ValueError("samples must contain at least one entry.")
        if any(not isinstance(s, Maxwell3DSample) for s in samples):
            raise TypeError("samples entries must be Maxwell3DSample.")
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

    def select(self, partition: str) -> tuple[Maxwell3DSample, ...]:
        """Return samples belonging to one split partition.

        Parameters
        ----------
        partition : {"train", "validation", "test"}
            Partition name.

        Returns
        -------
        tuple of Maxwell3DSample
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


def generate_3d_maxwell_dataset(
    config: Maxwell3DDatasetConfig,
    *,
    cache: MaxwellResultCache | None = None,
    policy: BatchPolicy | None = None,
) -> Maxwell3DDataset:
    """Generate a cached, split, versioned 3-D Maxwell dataset.

    Parameters
    ----------
    config : Maxwell3DDatasetConfig
        Complete, immutable dataset configuration.
    cache : MaxwellResultCache or None, optional
        When given, forward solves are cached by problem hash and
        resumable across calls.
    policy : BatchPolicy or None, optional
        Retry/concurrency policy forwarded to
        :func:`~pycsamt.forward.maxwell.batch.solve_batch`.

    Returns
    -------
    Maxwell3DDataset
        Every realization whose forward solve converged, split at
        the realization level.

    Raises
    ------
    ValueError
        If no realization converges.

    Examples
    --------
    >>> from pycsamt.ai.geology import GeologyGrid
    >>> grid = GeologyGrid.regular_3d(
    ...     nx=4, ny=4, nz=6, dx_m=200, dy_m=200, dz_m=100
    ... )
    >>> config = Maxwell3DDatasetConfig(
    ...     dataset_id="demo-3d-v1",
    ...     grid=grid,
    ...     correlation_length_x_m=(400.0, 800.0),
    ...     correlation_length_y_m=(400.0, 800.0),
    ...     correlation_length_z_m=(100.0, 200.0),
    ...     frequencies_hz=[50.0, 20.0],
    ...     station_xy_m=[[400.0, 400.0], [600.0, 600.0]],
    ...     n_realizations=1,
    ...     seed=0,
    ...     validation_fraction=0.0,
    ...     test_fraction=0.0,
    ... )
    >>> dataset = generate_3d_maxwell_dataset(config)  # doctest: +SKIP
    >>> dataset.samples[0].survey.shape  # doctest: +SKIP
    (2, 2, 2)
    """
    plan = SeedPlan(config.seed, namespace=config.namespace)
    grid = config.grid
    station_names = tuple(
        f"S{i:04d}" for i in range(config.station_xy_m.shape[0])
    )
    receivers = ReceiverSet(
        [[float(x), float(y), 0.0] for x, y in config.station_xy_m],
        station_names,
    )
    coordinates_m = [[float(x), float(y), 0.0] for x, y in config.station_xy_m]

    from pycsamt.api.view.progress import get_progress_bar

    realizations: dict[str, dict[str, Any]] = {}
    problems: dict[str, MaxwellProblem] = {}
    with get_progress_bar(
        total=config.n_realizations,
        desc="Generating 3-D realizations",
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
                config.station_xy_m,
                config.frequencies_hz,
                safety_factor=config.mesh_safety_factor,
                max_mesh_cells=config.max_mesh_cells,
                cells_per_skin_depth=config.cells_per_skin_depth,
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

    samples: list[Maxwell3DSample] = []
    rejected: list[str] = []

    with get_progress_bar(
        total=len(problems),
        desc="Solving 3-D Maxwell realizations",
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
                Maxwell3DSample(
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
            MT3DAdapter(max_cells=config.max_mesh_cells),
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
            "pycsamt.ai.training.dataset3d.generate_3d_maxwell_dataset"
        ),
        generator_version="0.1.0",
        configuration=config.to_dict(),
        split=split,
        sample_count=len(samples),
    )
    return Maxwell3DDataset(
        grid=grid,
        samples=tuple(samples),
        split=split,
        manifest=manifest,
        rejected=tuple(rejected),
    )
