# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Construct solver meshes from geological models and topography.

The builder preserves the canonical ``(z, x)`` / ``(z, y, x)`` array order,
adds geometric padding and air layers, extends boundary conductivity by nearest
cells, and records earth/air regions explicitly.  It constructs meshes and
models only; discretization and boundary conditions remain backend concerns.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .contracts import MaxwellMesh, MaxwellProblem, ReceiverSet

__all__ = [
    "MeshDesign",
    "MeshQuality",
    "SolverMeshModel",
    "skin_depth_m",
    "build_solver_mesh",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def _count(value: int, name: str) -> int:
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError(f"{name} must be an integer.")
    result = int(value)
    if result < 0:
        raise ValueError(f"{name} cannot be negative.")
    return result


def _positive(value: float, name: str, *, minimum: float = 0.0) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= minimum:
        raise ValueError(f"{name} must be finite and greater than {minimum}.")
    return result


def _centres_to_edges(centres: np.ndarray, name: str) -> np.ndarray:
    values = np.asarray(centres, dtype=float)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError(f"{name} must contain at least two centres.")
    if not np.all(np.isfinite(values)) or not np.all(np.diff(values) > 0):
        raise ValueError(f"{name} must be finite and strictly increasing.")
    middle = (values[:-1] + values[1:]) / 2.0
    first = values[0] - (middle[0] - values[0])
    last = values[-1] + (values[-1] - middle[-1])
    return np.concatenate([[first], middle, [last]])


def _pad_axis(
    edges: np.ndarray, before: int, after: int, expansion: float
) -> np.ndarray:
    left_widths = np.diff(edges)[0] * expansion ** np.arange(before, 0, -1)
    right_widths = np.diff(edges)[-1] * expansion ** np.arange(1, after + 1)
    left = edges[0] - np.cumsum(left_widths[::-1])[::-1]
    right = edges[-1] + np.cumsum(right_widths)
    return np.concatenate([left, edges, right])


def _nearest_indices(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    insertion = np.searchsorted(source, target)
    insertion = np.clip(insertion, 1, len(source) - 1)
    left = insertion - 1
    choose_right = np.abs(target - source[insertion]) < np.abs(
        target - source[left]
    )
    return np.where(choose_right, insertion, left)


def _slice_state(value: slice) -> list[int]:
    return [int(value.start), int(value.stop)]


@dataclass(frozen=True)
class MeshDesign:
    """Configure geometric padding, air treatment, and quality targets.

    Parameters
    ----------
    horizontal_padding_cells : int or pair of int, default=6
        Number of padding cells before/after each horizontal core axis. A
        scalar is applied symmetrically.
    bottom_padding_cells : int, default=8
        Number of cells beneath the geological model.
    air_layers : int, default=8
        Number of cells above the geological reference surface.
    padding_expansion : float, default=1.35
        Geometric growth factor away from the core mesh.
    air_expansion : float, default=1.25
        Geometric growth factor upward through the air.
    air_conductivity_s_m : float, default=1e-8
        Positive numerical conductivity assigned to air cells.
    minimum_cells_per_skin_depth : float, default=4.0
        Advisory resolution target evaluated at the smallest skin depth.
    maximum_adjacent_ratio : float, default=1.5
        Advisory upper bound for adjacent cell-width ratios.
    maximum_aspect_ratio : float, default=20.0
        Advisory upper bound across cell widths.

    Examples
    --------
    >>> design = MeshDesign(horizontal_padding_cells=(3, 5), air_layers=4)
    >>> design.horizontal_padding
    (3, 5)
    """

    horizontal_padding_cells: int | tuple[int, int] = 6
    bottom_padding_cells: int = 8
    air_layers: int = 8
    padding_expansion: float = 1.35
    air_expansion: float = 1.25
    air_conductivity_s_m: float = 1e-8
    minimum_cells_per_skin_depth: float = 4.0
    maximum_adjacent_ratio: float = 1.5
    maximum_aspect_ratio: float = 20.0

    def __post_init__(self) -> None:
        padding = self.horizontal_padding_cells
        if isinstance(padding, (int, np.integer)) and not isinstance(
            padding, bool
        ):
            normalized = (_count(padding, "horizontal_padding_cells"),) * 2
        else:
            try:
                normalized = tuple(padding)
            except TypeError as exc:
                raise TypeError(
                    "horizontal_padding_cells must be an integer or pair."
                ) from exc
            if len(normalized) != 2:
                raise ValueError(
                    "horizontal_padding_cells must contain two values."
                )
            normalized = (
                _count(normalized[0], "horizontal_padding_cells[0]"),
                _count(normalized[1], "horizontal_padding_cells[1]"),
            )
        object.__setattr__(self, "horizontal_padding_cells", normalized)
        object.__setattr__(
            self,
            "bottom_padding_cells",
            _count(self.bottom_padding_cells, "bottom_padding_cells"),
        )
        object.__setattr__(
            self, "air_layers", _count(self.air_layers, "air_layers")
        )
        object.__setattr__(
            self,
            "padding_expansion",
            _positive(
                self.padding_expansion,
                "padding_expansion",
                minimum=1.0 - 1e-15,
            ),
        )
        object.__setattr__(
            self,
            "air_expansion",
            _positive(
                self.air_expansion, "air_expansion", minimum=1.0 - 1e-15
            ),
        )
        object.__setattr__(
            self,
            "air_conductivity_s_m",
            _positive(self.air_conductivity_s_m, "air_conductivity_s_m"),
        )
        object.__setattr__(
            self,
            "minimum_cells_per_skin_depth",
            _positive(
                self.minimum_cells_per_skin_depth,
                "minimum_cells_per_skin_depth",
            ),
        )
        object.__setattr__(
            self,
            "maximum_adjacent_ratio",
            _positive(
                self.maximum_adjacent_ratio,
                "maximum_adjacent_ratio",
                minimum=1.0 - 1e-15,
            ),
        )
        object.__setattr__(
            self,
            "maximum_aspect_ratio",
            _positive(
                self.maximum_aspect_ratio,
                "maximum_aspect_ratio",
                minimum=1.0 - 1e-15,
            ),
        )

    @property
    def horizontal_padding(self) -> tuple[int, int]:
        """Return normalized before/after horizontal padding counts.

        Returns
        -------
        tuple of int
            Padding cells on the low and high sides of each horizontal axis.

        Examples
        --------
        >>> MeshDesign(horizontal_padding_cells=3).horizontal_padding
        (3, 3)
        """
        return self.horizontal_padding_cells

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible design representation.

        Returns
        -------
        dict
            Versioned design state.

        Examples
        --------
        >>> MeshDesign(air_layers=2).to_dict()["air_layers"]
        2
        """
        return {
            "schema_version": 1,
            "horizontal_padding_cells": list(self.horizontal_padding),
            "bottom_padding_cells": self.bottom_padding_cells,
            "air_layers": self.air_layers,
            "padding_expansion": self.padding_expansion,
            "air_expansion": self.air_expansion,
            "air_conductivity_s_m": self.air_conductivity_s_m,
            "minimum_cells_per_skin_depth": self.minimum_cells_per_skin_depth,
            "maximum_adjacent_ratio": self.maximum_adjacent_ratio,
            "maximum_aspect_ratio": self.maximum_aspect_ratio,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> MeshDesign:
        """Restore a validated mesh design.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        MeshDesign
            Restored design.

        Examples
        --------
        >>> design = MeshDesign(horizontal_padding_cells=2)
        >>> MeshDesign.from_dict(design.to_dict()).horizontal_padding
        (2, 2)
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported MeshDesign schema version.")
        values = dict(data)
        values.pop("schema_version")
        values["horizontal_padding_cells"] = tuple(
            values["horizontal_padding_cells"]
        )
        return cls(**values)


@dataclass(frozen=True)
class MeshQuality:
    """Summarize numerical mesh quality and skin-depth resolution.

    Parameters
    ----------
    cell_count : int
        Total number of cells including air and padding.
    minimum_cell_width_m, maximum_cell_width_m : float
        Extreme cell widths over all axes.
    maximum_aspect_ratio : float
        Ratio of maximum to minimum cell width.
    maximum_adjacent_ratio : float
        Worst neighboring width expansion on any axis.
    minimum_skin_depth_m : float
        Smallest skin depth across the requested physics range.
    cells_per_minimum_skin_depth : float
        Skin depth divided by the largest core cell width.
    warnings : tuple of str
        Advisory quality violations.

    Examples
    --------
    >>> quality = MeshQuality(10, 1, 5, 5, 1.2, 100, 20, ())
    >>> quality.acceptable
    True
    """

    cell_count: int
    minimum_cell_width_m: float
    maximum_cell_width_m: float
    maximum_aspect_ratio: float
    maximum_adjacent_ratio: float
    minimum_skin_depth_m: float
    cells_per_minimum_skin_depth: float
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        cell_count = _count(self.cell_count, "cell_count")
        if cell_count < 1:
            raise ValueError("cell_count must be positive.")
        numeric_names = (
            "minimum_cell_width_m",
            "maximum_cell_width_m",
            "maximum_aspect_ratio",
            "maximum_adjacent_ratio",
            "minimum_skin_depth_m",
            "cells_per_minimum_skin_depth",
        )
        values = {
            name: _positive(getattr(self, name), name)
            for name in numeric_names
        }
        if values["maximum_cell_width_m"] < values["minimum_cell_width_m"]:
            raise ValueError(
                "maximum_cell_width_m cannot be smaller than minimum_cell_width_m."
            )
        messages = tuple(str(value).strip() for value in self.warnings)
        if any(not value for value in messages):
            raise ValueError("warnings cannot contain empty messages.")
        object.__setattr__(self, "cell_count", cell_count)
        for name, value in values.items():
            object.__setattr__(self, name, value)
        object.__setattr__(self, "warnings", messages)

    @property
    def acceptable(self) -> bool:
        """Return whether no advisory quality limits were violated.

        Returns
        -------
        bool
            True when :attr:`warnings` is empty.

        Examples
        --------
        >>> MeshQuality(1, 1, 1, 1, 1, 1, 1, ("coarse",)).acceptable
        False
        """
        return not self.warnings

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible mesh-quality diagnostics.

        Returns
        -------
        dict
            Numeric diagnostics and warnings.

        Examples
        --------
        >>> MeshQuality(1, 1, 1, 1, 1, 1, 1).to_dict()["cell_count"]
        1
        """
        return {
            "cell_count": self.cell_count,
            "minimum_cell_width_m": self.minimum_cell_width_m,
            "maximum_cell_width_m": self.maximum_cell_width_m,
            "maximum_aspect_ratio": self.maximum_aspect_ratio,
            "maximum_adjacent_ratio": self.maximum_adjacent_ratio,
            "minimum_skin_depth_m": self.minimum_skin_depth_m,
            "cells_per_minimum_skin_depth": self.cells_per_minimum_skin_depth,
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class SolverMeshModel:
    """Store a padded mesh, conductivity, regions, and construction record.

    Parameters
    ----------
    mesh : MaxwellMesh
        Solver-neutral padded mesh.
    conductivity_s_m : ndarray
        Positive conductivity shaped like ``mesh``. Air cells contain the
        configured small numerical conductivity.
    earth_mask : ndarray of bool
        True for cells on or below local terrain.
    core_slices : tuple of slice
        Geological core location in canonical array order.
    design : MeshDesign
        Construction settings.
    quality : MeshQuality
        Mesh-quality diagnostics for the requested frequency range.
    source_shape : tuple of int
        Original geological model shape.

    Examples
    --------
    Instances are normally created with :func:`build_solver_mesh`.
    """

    mesh: MaxwellMesh
    conductivity_s_m: np.ndarray
    earth_mask: np.ndarray
    core_slices: tuple[slice, ...]
    design: MeshDesign
    quality: MeshQuality
    source_shape: tuple[int, ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.mesh, MaxwellMesh)
            or not isinstance(self.design, MeshDesign)
            or not isinstance(self.quality, MeshQuality)
        ):
            raise TypeError(
                "mesh, design, and quality must use Maxwell mesh contract types."
            )
        conductivity = np.asarray(self.conductivity_s_m, dtype=float)
        earth = np.asarray(self.earth_mask, dtype=bool)
        if (
            conductivity.shape != self.mesh.shape
            or earth.shape != self.mesh.shape
        ):
            raise ValueError(
                "conductivity_s_m and earth_mask must have mesh shape."
            )
        if not np.all(np.isfinite(conductivity)) or np.any(conductivity <= 0):
            raise ValueError("conductivity_s_m must be positive and finite.")
        if len(self.core_slices) != self.mesh.dimension or tuple(
            conductivity[self.core_slices].shape
        ) != tuple(self.source_shape):
            raise ValueError("core_slices must select exactly source_shape.")
        object.__setattr__(self, "conductivity_s_m", _readonly(conductivity))
        object.__setattr__(self, "earth_mask", _readonly(earth, bool))
        object.__setattr__(self, "core_slices", tuple(self.core_slices))
        object.__setattr__(
            self,
            "source_shape",
            tuple(int(value) for value in self.source_shape),
        )

    @property
    def air_mask(self) -> np.ndarray:
        """Return the read-only complement of the earth mask.

        Returns
        -------
        ndarray of bool
            Air-region mask shaped like ``mesh``.

        Examples
        --------
        ``air_mask`` and ``earth_mask`` always partition the complete mesh.
        """
        return _readonly(~self.earth_mask, bool)

    @property
    def model_hash(self) -> str:
        """Return a deterministic digest of mesh, model, regions, and design.

        Returns
        -------
        str
            SHA-256 digest suitable for provenance checks.

        Examples
        --------
        A valid model hash always contains 64 hexadecimal characters.
        """
        digest = hashlib.sha256()
        digest.update(
            np.ascontiguousarray(self.conductivity_s_m, dtype="<f8").tobytes()
        )
        digest.update(
            np.ascontiguousarray(self.earth_mask, dtype=np.uint8).tobytes()
        )
        digest.update(
            json.dumps(
                self.provenance(), sort_keys=True, separators=(",", ":")
            ).encode("utf-8")
        )
        return digest.hexdigest()

    def assess_receivers(self, receivers: ReceiverSet) -> tuple[str, ...]:
        """Return receiver-placement errors without modifying coordinates.

        Parameters
        ----------
        receivers : ReceiverSet
            Candidate locations in the mesh coordinate system.

        Returns
        -------
        tuple of str
            Empty when all receivers lie inside mesh bounds and no receiver is
            below the discretized local terrain surface.

        Examples
        --------
        Use this check before :meth:`to_problem` when receiver coordinates are
        assembled independently from the mesh.
        """
        if not isinstance(receivers, ReceiverSet):
            raise TypeError("receivers must be a ReceiverSet.")
        if receivers.dimension != self.mesh.dimension:
            return ("receiver and mesh dimensions differ",)
        errors = []
        coordinates = receivers.coordinates_m
        x = coordinates[:, 0]
        if np.any(
            (x < self.mesh.x_edges_m[0]) | (x > self.mesh.x_edges_m[-1])
        ):
            errors.append("receiver x coordinates fall outside the mesh")
        if self.mesh.dimension == 3:
            y = coordinates[:, 1]
            z = coordinates[:, 2]
            if np.any(
                (y < self.mesh.y_edges_m[0]) | (y > self.mesh.y_edges_m[-1])
            ):
                errors.append("receiver y coordinates fall outside the mesh")
        else:
            z = coordinates[:, 1]
        if np.any(
            (z < self.mesh.z_edges_m[0]) | (z > self.mesh.z_edges_m[-1])
        ):
            errors.append("receiver z coordinates fall outside the mesh")
        if not errors:
            ix = np.clip(
                np.searchsorted(self.mesh.x_edges_m, x, side="right") - 1,
                0,
                self.mesh.shape[-1] - 1,
            )
            if self.mesh.dimension == 2:
                columns = self.earth_mask[:, ix]
            else:
                iy = np.clip(
                    np.searchsorted(self.mesh.y_edges_m, y, side="right") - 1,
                    0,
                    self.mesh.shape[1] - 1,
                )
                columns = self.earth_mask[:, iy, ix]
            first_earth = np.argmax(columns, axis=0)
            surface_depth = self.mesh.z_edges_m[first_earth]
            tolerance = (
                np.finfo(float).eps
                * np.maximum(1.0, np.abs(surface_depth))
                * 16
            )
            if np.any(z > surface_depth + tolerance):
                errors.append(
                    "receiver z coordinates fall below local terrain"
                )
        return tuple(errors)

    def to_problem(
        self,
        frequencies_hz: Sequence[float],
        receivers: ReceiverSet,
        *,
        components: Sequence[str] = ("zxy", "zyx"),
        mark_air_inactive: bool = False,
        time_dependence: str = "exp(+iwt)",
        magnetic_permeability_h_m: float = 4e-7 * np.pi,
        metadata: Mapping[str, Any] | None = None,
    ) -> MaxwellProblem:
        """Create a validated Maxwell problem from this mesh model.

        Parameters
        ----------
        frequencies_hz : sequence of float
            Positive simulation frequencies.
        receivers : ReceiverSet
            Receiver locations matching the mesh dimension.
        components : sequence of str, default=("zxy", "zyx")
            Requested canonical impedance components.
        mark_air_inactive : bool, default=False
            Use ``earth_mask`` as active cells. Keep false for formulations
            that solve conductive air explicitly.
        time_dependence : str, default="exp(+iwt)"
            Complex phasor convention.
        magnetic_permeability_h_m : float, default=4e-7*pi
            Uniform magnetic permeability.
        metadata : mapping or None, optional
            Additional problem provenance.

        Returns
        -------
        MaxwellProblem
            Solver-neutral problem ready for adapter assessment.

        Examples
        --------
        The generated problem includes ``mesh_model_hash`` in its metadata.
        """
        errors = self.assess_receivers(receivers)
        if errors:
            raise ValueError("; ".join(errors))
        provenance = {} if metadata is None else dict(metadata)
        provenance["mesh_model_hash"] = self.model_hash
        provenance["air_treatment"] = (
            "inactive" if mark_air_inactive else "conductive"
        )
        return MaxwellProblem(
            self.mesh,
            self.conductivity_s_m,
            frequencies_hz,
            receivers,
            tuple(components),
            self.earth_mask
            if mark_air_inactive
            else np.ones(self.mesh.shape, bool),
            time_dependence,
            magnetic_permeability_h_m,
            provenance,
        )

    def provenance(self) -> dict[str, Any]:
        """Return JSON-compatible mesh-construction provenance.

        Returns
        -------
        dict
            Mesh, design, quality, source shape, and core slices.

        Examples
        --------
        The returned schema version is currently one.
        """
        return {
            "schema_version": 1,
            "mesh": self.mesh.to_dict(),
            "design": self.design.to_dict(),
            "quality": self.quality.to_dict(),
            "source_shape": list(self.source_shape),
            "core_slices": [_slice_state(value) for value in self.core_slices],
        }

    def to_npz(self, path: str | Path) -> Path:
        """Persist a solver mesh model without enabling pickle.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination archive.

        Returns
        -------
        pathlib.Path
            Requested destination.

        Examples
        --------
        Archives can be restored with :meth:`from_npz`.
        """
        target = Path(path)
        np.savez_compressed(
            target,
            conductivity_s_m=self.conductivity_s_m,
            earth_mask=self.earth_mask,
            provenance_json=np.array(
                json.dumps(self.provenance(), sort_keys=True)
            ),
        )
        return target

    @classmethod
    def from_npz(cls, path: str | Path) -> SolverMeshModel:
        """Restore and validate a solver mesh archive without pickle.

        Parameters
        ----------
        path : str or pathlib.Path
            Archive written by :meth:`to_npz`.

        Returns
        -------
        SolverMeshModel
            Restored immutable mesh model.

        Examples
        --------
        Restored arrays remain read-only after construction.
        """
        with np.load(Path(path), allow_pickle=False) as archive:
            state = json.loads(str(archive["provenance_json"].item()))
            if state.get("schema_version") != 1:
                raise ValueError("unsupported SolverMeshModel schema version.")
            quality_state = dict(state["quality"])
            quality_state["warnings"] = tuple(
                quality_state.get("warnings", ())
            )
            return cls(
                MaxwellMesh.from_dict(state["mesh"]),
                archive["conductivity_s_m"],
                archive["earth_mask"],
                tuple(slice(*value) for value in state["core_slices"]),
                MeshDesign.from_dict(state["design"]),
                MeshQuality(**quality_state),
                tuple(state["source_shape"]),
            )


def skin_depth_m(resistivity_ohm_m: Any, frequency_hz: Any) -> np.ndarray:
    """Calculate electromagnetic skin depth for a non-magnetic conductor.

    Parameters
    ----------
    resistivity_ohm_m, frequency_hz : array-like
        Positive resistivity and frequency, broadcast using NumPy rules.

    Returns
    -------
    ndarray
        Skin depth in metres, ``sqrt(rho / (pi * mu0 * f))``.

    Examples
    --------
    >>> round(float(skin_depth_m(100, 1)))
    5033
    >>> skin_depth_m([100, 400], 1).shape
    (2,)
    """
    resistivity = np.asarray(resistivity_ohm_m, dtype=float)
    frequency = np.asarray(frequency_hz, dtype=float)
    if not np.all(np.isfinite(resistivity)) or np.any(resistivity <= 0):
        raise ValueError("resistivity_ohm_m must be positive and finite.")
    if not np.all(np.isfinite(frequency)) or np.any(frequency <= 0):
        raise ValueError("frequency_hz must be positive and finite.")
    return np.sqrt(resistivity / (np.pi * (4e-7 * np.pi) * frequency))


def build_solver_mesh(
    grid: Any,
    *,
    conductivity_s_m: Any | None = None,
    resistivity_ohm_m: Any | None = None,
    frequencies_hz: Sequence[float],
    topography: Any | None = None,
    design: MeshDesign | None = None,
) -> SolverMeshModel:
    """Build a padded Maxwell mesh from a geological cell-centre model.

    Parameters
    ----------
    grid : pycsamt.ai.geology.GeologyGrid
        Source grid in canonical geological order. Its upper cell edge must be
        depth zero, the reference used by topography and receiver coordinates.
    conductivity_s_m, resistivity_ohm_m : array-like or None
        Supply exactly one positive model shaped like ``grid``.
    frequencies_hz : sequence of float
        Frequencies used only for skin-depth quality diagnostics.
    topography : pycsamt.ai.geology.TopographicSurface or None, optional
        Terrain aligned to ``grid``. When omitted, the geological top edge is
        treated as a flat earth surface.
    design : MeshDesign or None, optional
        Padding and quality configuration.

    Returns
    -------
    SolverMeshModel
        Padded conductivity, earth/air regions, core mapping, and diagnostics.

    Examples
    --------
    >>> from pycsamt.ai.geology import GeologyGrid
    >>> grid = GeologyGrid.regular_2d(nx=3, nz=2, dx_m=100, dz_m=50)
    >>> model = build_solver_mesh(
    ...     grid,
    ...     resistivity_ohm_m=np.full(grid.shape, 100),
    ...     frequencies_hz=[10, 1],
    ...     design=MeshDesign(
    ...         horizontal_padding_cells=1,
    ...         bottom_padding_cells=1,
    ...         air_layers=1,
    ...     ),
    ... )
    >>> model.mesh.shape, model.core_slices
    ((4, 5), (slice(1, 3, None), slice(1, 4, None)))
    """
    from ...ai.geology import GeologyGrid, TopographicSurface

    if not isinstance(grid, GeologyGrid):
        raise TypeError("grid must be a GeologyGrid.")
    if (conductivity_s_m is None) == (resistivity_ohm_m is None):
        raise ValueError(
            "supply exactly one of conductivity_s_m or resistivity_ohm_m."
        )
    source = np.asarray(
        conductivity_s_m
        if conductivity_s_m is not None
        else resistivity_ohm_m,
        dtype=float,
    )
    if (
        source.shape != grid.shape
        or not np.all(np.isfinite(source))
        or np.any(source <= 0)
    ):
        raise ValueError(
            f"the source property must be positive, finite, and shaped {grid.shape}."
        )
    source_conductivity = (
        source if conductivity_s_m is not None else 1.0 / source
    )
    frequencies = np.asarray(frequencies_hz, dtype=float)
    if (
        frequencies.ndim != 1
        or len(frequencies) < 1
        or not np.all(np.isfinite(frequencies))
        or np.any(frequencies <= 0)
    ):
        raise ValueError(
            "frequencies_hz must be a non-empty positive finite vector."
        )
    configuration = MeshDesign() if design is None else design
    if not isinstance(configuration, MeshDesign):
        raise TypeError("design must be a MeshDesign or None.")
    aligned_topography = isinstance(topography, TopographicSurface) and (
        topography.grid.dimension == grid.dimension
        and topography.grid.crs == grid.crs
        and np.array_equal(topography.grid.x_m, grid.x_m)
        and np.array_equal(topography.grid.z_m, grid.z_m)
        and (
            (topography.grid.y_m is None and grid.y_m is None)
            or np.array_equal(topography.grid.y_m, grid.y_m)
        )
    )
    if topography is not None and not aligned_topography:
        raise ValueError(
            "topography must be a TopographicSurface aligned to grid."
        )

    before, after = configuration.horizontal_padding
    x_core_edges = _centres_to_edges(grid.x_m, "grid.x_m")
    x_edges = _pad_axis(
        x_core_edges, before, after, configuration.padding_expansion
    )
    z_core_edges = _centres_to_edges(grid.z_m, "grid.z_m")
    if not np.isclose(
        z_core_edges[0],
        0.0,
        rtol=0.0,
        atol=max(1e-10, abs(np.diff(z_core_edges)[0]) * 1e-10),
    ):
        raise ValueError("grid's upper cell edge must be depth zero.")
    z_earth_edges = _pad_axis(
        z_core_edges,
        0,
        configuration.bottom_padding_cells,
        configuration.padding_expansion,
    )
    if configuration.air_layers:
        first_width = np.diff(z_core_edges)[0]
        upward_widths = first_width * configuration.air_expansion ** np.arange(
            configuration.air_layers
        )
        air_edges = z_core_edges[0] - np.cumsum(upward_widths)[::-1]
        z_edges = np.concatenate([air_edges, z_earth_edges])
    else:
        z_edges = z_earth_edges
    y_edges = None
    if grid.dimension == 3:
        y_core_edges = _centres_to_edges(grid.y_m, "grid.y_m")
        y_edges = _pad_axis(
            y_core_edges, before, after, configuration.padding_expansion
        )
    mesh = MaxwellMesh(x_edges, z_edges, y_edges, grid.crs)

    centres = mesh.cell_centres_m
    ix = _nearest_indices(grid.x_m, centres["x"])
    iz = _nearest_indices(grid.z_m, centres["z"])
    if grid.dimension == 2:
        mapped = source_conductivity[np.ix_(iz, ix)]
        surface_core = (
            np.full(len(grid.x_m), z_core_edges[0])
            if topography is None
            else topography.surface_depth_m
        )
        surface = np.interp(
            centres["x"],
            grid.x_m,
            surface_core,
            left=surface_core[0],
            right=surface_core[-1],
        )
        earth = centres["z"][:, None] >= surface[None, :]
        core_slices = (
            slice(
                configuration.air_layers,
                configuration.air_layers + len(grid.z_m),
            ),
            slice(before, before + len(grid.x_m)),
        )
    else:
        iy = _nearest_indices(grid.y_m, centres["y"])
        mapped = source_conductivity[np.ix_(iz, iy, ix)]
        surface_core = (
            np.full((len(grid.y_m), len(grid.x_m)), z_core_edges[0])
            if topography is None
            else topography.surface_depth_m
        )
        surface_x = np.vstack(
            [
                np.interp(
                    centres["x"], grid.x_m, row, left=row[0], right=row[-1]
                )
                for row in surface_core
            ]
        )
        surface = np.vstack(
            [
                np.interp(
                    centres["y"],
                    grid.y_m,
                    surface_x[:, column],
                    left=surface_x[0, column],
                    right=surface_x[-1, column],
                )
                for column in range(surface_x.shape[1])
            ]
        ).T
        earth = centres["z"][:, None, None] >= surface[None, :, :]
        core_slices = (
            slice(
                configuration.air_layers,
                configuration.air_layers + len(grid.z_m),
            ),
            slice(before, before + len(grid.y_m)),
            slice(before, before + len(grid.x_m)),
        )
    mapped = np.where(earth, mapped, configuration.air_conductivity_s_m)

    widths = mesh.cell_widths_m
    flat_widths = np.concatenate(tuple(widths.values()))
    adjacent = []
    for values in widths.values():
        if len(values) > 1:
            ratios = values[1:] / values[:-1]
            adjacent.extend(np.maximum(ratios, 1.0 / ratios))
    maximum_adjacent = float(max(adjacent, default=1.0))
    minimum_skin = float(
        np.min(skin_depth_m(1.0 / source_conductivity, np.max(frequencies)))
    )
    core_widths = [
        np.max(np.diff(x_core_edges)),
        np.max(np.diff(z_core_edges)),
    ]
    if grid.dimension == 3:
        core_widths.append(np.max(np.diff(y_core_edges)))
    cells_per_skin = minimum_skin / float(max(core_widths))
    aspect = float(np.max(flat_widths) / np.min(flat_widths))
    quality_warnings = []
    if maximum_adjacent > configuration.maximum_adjacent_ratio:
        quality_warnings.append(
            f"adjacent cell ratio {maximum_adjacent:.3g} exceeds {configuration.maximum_adjacent_ratio:.3g}"
        )
    if aspect > configuration.maximum_aspect_ratio:
        quality_warnings.append(
            f"global cell-width ratio {aspect:.3g} exceeds {configuration.maximum_aspect_ratio:.3g}"
        )
    if cells_per_skin < configuration.minimum_cells_per_skin_depth:
        quality_warnings.append(
            f"minimum skin depth has {cells_per_skin:.3g} core cells; target is {configuration.minimum_cells_per_skin_depth:.3g}"
        )
    quality = MeshQuality(
        int(np.prod(mesh.shape)),
        float(np.min(flat_widths)),
        float(np.max(flat_widths)),
        aspect,
        maximum_adjacent,
        minimum_skin,
        float(cells_per_skin),
        tuple(quality_warnings),
    )
    return SolverMeshModel(
        mesh, mapped, earth, core_slices, configuration, quality, grid.shape
    )
