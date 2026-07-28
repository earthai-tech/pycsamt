# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Rotated ellipsoidal lenses embedded in layered electrical geology.

Lens geometry is evaluated at cell centres.  A lens may have a sharp boundary
or a finite transition shell blended in log10-resistivity space.  Overlap is
never implicit: callers choose a conflict policy and the resulting overlap
count is retained for auditing.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..data.manifest import canonical_hash
from .fields import GeologyGrid
from .layers import ElectricalLayer, LayeredGeology

__all__ = ["EllipsoidalLens", "LensGeology", "insert_lenses"]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class EllipsoidalLens:
    """Define a rotated 2-D elliptical or 3-D ellipsoidal electrical body.

    Parameters
    ----------
    name : str
        Unique non-empty body name.
    center_x_m, center_z_m : float
        Horizontal and depth coordinates of the body centre in metres.
    radius_x_m, radius_z_m : float
        Positive principal semi-axis lengths in metres.
    resistivity_ohm_m : float
        Positive target resistivity at the body core.
    center_y_m, radius_y_m : float or None, optional
        Second horizontal centre and radius. Both are required on a 3-D grid
        and both must be omitted on a 2-D grid.
    azimuth_deg : float, default=0.0
        Clockwise rotation of the horizontal x principal axis in 3-D.
    dip_deg : float, default=0.0
        Downward rotation of the x principal axis toward increasing depth.
        It rotates the x-z ellipse in 2-D and the azimuthal major-z plane in 3-D.
    transition_fraction : float, default=0.0
        Fraction of the normalized outer radius occupied by a smooth transition
        shell. It must lie in ``[0, 1)``. Zero gives a sharp boundary.

    Examples
    --------
    A dipping conductive lens in a 2-D section:

    >>> lens = EllipsoidalLens(
    ...     "conductor",
    ...     center_x_m=1000,
    ...     center_z_m=400,
    ...     radius_x_m=500,
    ...     radius_z_m=100,
    ...     resistivity_ohm_m=5,
    ...     dip_deg=20,
    ... )
    >>> lens.dimension
    2

    A 3-D ellipsoid additionally declares y geometry:

    >>> body = EllipsoidalLens(
    ...     "body", 0, 300, 400, 100, 10, center_y_m=0, radius_y_m=200
    ... )
    >>> body.dimension
    3
    """

    name: str
    center_x_m: float
    center_z_m: float
    radius_x_m: float
    radius_z_m: float
    resistivity_ohm_m: float
    center_y_m: float | None = None
    radius_y_m: float | None = None
    azimuth_deg: float = 0.0
    dip_deg: float = 0.0
    transition_fraction: float = 0.0

    def __post_init__(self) -> None:
        name = str(self.name).strip()
        if not name:
            raise ValueError("name cannot be empty.")
        finite_names = (
            "center_x_m",
            "center_z_m",
            "radius_x_m",
            "radius_z_m",
            "resistivity_ohm_m",
            "azimuth_deg",
            "dip_deg",
            "transition_fraction",
        )
        values = {field: float(getattr(self, field)) for field in finite_names}
        if not np.all(np.isfinite(tuple(values.values()))):
            raise ValueError("lens numerical parameters must be finite.")
        for field in ("radius_x_m", "radius_z_m", "resistivity_ohm_m"):
            if values[field] <= 0:
                raise ValueError(f"{field} must be positive.")
        if values["dip_deg"] < -90 or values["dip_deg"] > 90:
            raise ValueError("dip_deg must lie in [-90, 90].")
        if values["transition_fraction"] < 0 or values["transition_fraction"] >= 1:
            raise ValueError("transition_fraction must lie in [0, 1).")
        if (self.center_y_m is None) != (self.radius_y_m is None):
            raise ValueError("center_y_m and radius_y_m must be supplied together.")
        center_y = None
        radius_y = None
        if self.center_y_m is not None:
            center_y = float(self.center_y_m)
            radius_y = float(self.radius_y_m)
            if not np.isfinite(center_y) or not np.isfinite(radius_y) or radius_y <= 0:
                raise ValueError(
                    "center_y_m must be finite and radius_y_m finite and positive."
                )
        object.__setattr__(self, "name", name)
        for field, value in values.items():
            object.__setattr__(self, field, value)
        object.__setattr__(self, "center_y_m", center_y)
        object.__setattr__(self, "radius_y_m", radius_y)
        object.__setattr__(self, "azimuth_deg", values["azimuth_deg"] % 180.0)

    @property
    def dimension(self) -> int:
        """Return the dimensionality implied by the lens geometry.

        Returns
        -------
        {2, 3}
            Two when y geometry is absent, otherwise three.

        Examples
        --------
        >>> EllipsoidalLens("a", 0, 1, 2, 1, 10).dimension
        2
        """
        return 2 if self.center_y_m is None else 3

    def validate_grid(
        self, grid: GeologyGrid, *, require_intersection: bool = True
    ) -> None:
        """Validate dimensional and spatial compatibility with a grid.

        Parameters
        ----------
        grid : GeologyGrid
            Candidate 2-D section or 3-D volume.
        require_intersection : bool, default=True
            Require at least one cell centre inside the body.

        Returns
        -------
        None
            Successful return means the lens can be rasterized.

        Raises
        ------
        ValueError
            If dimensions differ or no cell centre intersects the lens.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
        >>> lens = EllipsoidalLens("a", 200, 100, 150, 75, 10)
        >>> lens.validate_grid(grid) is None
        True
        """
        if not isinstance(grid, GeologyGrid):
            raise TypeError("grid must be a GeologyGrid.")
        if grid.dimension != self.dimension:
            raise ValueError("lens and grid dimensions differ.")
        if grid.dimension == 2 and self.azimuth_deg != 0.0:
            raise ValueError("azimuth_deg must be zero for a 2-D lens.")
        if require_intersection and not np.any(self.normalized_radius(grid) <= 1.0):
            raise ValueError(
                f"lens {self.name!r} does not intersect any grid cell centre."
            )

    def normalized_radius(self, grid: GeologyGrid) -> np.ndarray:
        """Evaluate dimensionless ellipsoidal radius at every cell centre.

        Parameters
        ----------
        grid : GeologyGrid
            Grid with the same dimensionality as the lens.

        Returns
        -------
        ndarray
            Array shaped like ``grid``. Values at or below one are inside the
            lens; zero is the geometric centre.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
        >>> lens = EllipsoidalLens("a", 200, 100, 150, 75, 10)
        >>> lens.normalized_radius(grid).shape
        (4, 4)
        """
        if not isinstance(grid, GeologyGrid):
            raise TypeError("grid must be a GeologyGrid.")
        if grid.dimension != self.dimension:
            raise ValueError("lens and grid dimensions differ.")
        dip = np.deg2rad(self.dip_deg)
        if grid.dimension == 2:
            z, x = np.meshgrid(grid.z_m, grid.x_m, indexing="ij")
            dx = x - self.center_x_m
            dz = z - self.center_z_m
            major = dx * np.cos(dip) + dz * np.sin(dip)
            vertical = -dx * np.sin(dip) + dz * np.cos(dip)
            squared = np.square(major / self.radius_x_m) + np.square(
                vertical / self.radius_z_m
            )
        else:
            z, y, x = np.meshgrid(grid.z_m, grid.y_m, grid.x_m, indexing="ij")
            dx = x - self.center_x_m
            dy = y - self.center_y_m
            dz = z - self.center_z_m
            azimuth = np.deg2rad(self.azimuth_deg)
            horizontal_major = dx * np.cos(azimuth) + dy * np.sin(azimuth)
            horizontal_minor = -dx * np.sin(azimuth) + dy * np.cos(azimuth)
            major = horizontal_major * np.cos(dip) + dz * np.sin(dip)
            vertical = -horizontal_major * np.sin(dip) + dz * np.cos(dip)
            squared = (
                np.square(major / self.radius_x_m)
                + np.square(horizontal_minor / self.radius_y_m)
                + np.square(vertical / self.radius_z_m)
            )
        return np.sqrt(squared)

    def blend_weight(self, grid: GeologyGrid) -> np.ndarray:
        """Return the lens contribution weight at every cell centre.

        Parameters
        ----------
        grid : GeologyGrid
            Compatible geological grid.

        Returns
        -------
        ndarray
            Values in ``[0, 1]``. A sharp lens is one inside and zero outside;
            a transition shell uses cubic smoothstep interpolation.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
        >>> lens = EllipsoidalLens(
        ...     "a", 200, 100, 150, 75, 10, transition_fraction=0.2
        ... )
        >>> np.all(
        ...     (lens.blend_weight(grid) >= 0) & (lens.blend_weight(grid) <= 1)
        ... )
        True
        """
        radius = self.normalized_radius(grid)
        if self.transition_fraction == 0:
            return (radius <= 1.0).astype(float)
        core = 1.0 - self.transition_fraction
        position = np.clip((1.0 - radius) / self.transition_fraction, 0.0, 1.0)
        weight = position * position * (3.0 - 2.0 * position)
        weight[radius <= core] = 1.0
        return weight

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable lens definition.

        Returns
        -------
        dict
            Versioned geometry, orientation, electrical value, and transition.

        Examples
        --------
        >>> EllipsoidalLens("a", 0, 1, 2, 1, 10).to_dict()["name"]
        'a'
        """
        return {
            "schema_version": 1,
            "name": self.name,
            "center_x_m": self.center_x_m,
            "center_y_m": self.center_y_m,
            "center_z_m": self.center_z_m,
            "radius_x_m": self.radius_x_m,
            "radius_y_m": self.radius_y_m,
            "radius_z_m": self.radius_z_m,
            "resistivity_ohm_m": self.resistivity_ohm_m,
            "azimuth_deg": self.azimuth_deg,
            "dip_deg": self.dip_deg,
            "transition_fraction": self.transition_fraction,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EllipsoidalLens:
        """Restore a validated lens definition.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        EllipsoidalLens
            Immutable body definition.

        Examples
        --------
        >>> lens = EllipsoidalLens("a", 0, 1, 2, 1, 10)
        >>> EllipsoidalLens.from_dict(lens.to_dict()) == lens
        True
        """
        if data.get("schema_version", 1) != 1:
            raise ValueError("unsupported EllipsoidalLens schema version.")
        return cls(
            name=data["name"],
            center_x_m=data["center_x_m"],
            center_z_m=data["center_z_m"],
            radius_x_m=data["radius_x_m"],
            radius_z_m=data["radius_z_m"],
            resistivity_ohm_m=data["resistivity_ohm_m"],
            center_y_m=data.get("center_y_m"),
            radius_y_m=data.get("radius_y_m"),
            azimuth_deg=data.get("azimuth_deg", 0.0),
            dip_deg=data.get("dip_deg", 0.0),
            transition_fraction=data.get("transition_fraction", 0.0),
        )


def _compose_lenses(
    base: LayeredGeology,
    bodies: tuple[EllipsoidalLens, ...],
    conflict_policy: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    masks = []
    weights = []
    for lens in bodies:
        lens.validate_grid(base.grid)
        radius = lens.normalized_radius(base.grid)
        masks.append(radius <= 1.0)
        weights.append(lens.blend_weight(base.grid))
    overlap = np.sum(np.stack(masks), axis=0, dtype=np.int16)
    if conflict_policy == "error" and np.any(overlap > 1):
        raise ValueError("lens envelopes overlap under the 'error' conflict policy.")
    result = np.array(base.resistivity_ohm_m, copy=True)
    indices = np.full(base.grid.shape, -1, dtype=np.int16)
    for index, (lens, mask, weight) in enumerate(zip(bodies, masks, weights)):
        candidate = np.power(
            10.0,
            (1.0 - weight) * np.log10(base.resistivity_ohm_m)
            + weight * np.log10(lens.resistivity_ohm_m),
        )
        if conflict_policy in {"error", "last"}:
            selected = mask
        elif conflict_policy == "first":
            selected = mask & (indices < 0)
        elif conflict_policy == "most_conductive":
            selected = mask & (candidate < result)
        else:
            selected = mask & (candidate > result)
        result[selected] = candidate[selected]
        indices[selected] = index
    return result, indices, overlap


@dataclass(frozen=True)
class LensGeology:
    """Immutable layered model after one or more lens insertions.

    Parameters
    ----------
    base : LayeredGeology
        Self-contained stratigraphic model before lens insertion.
    lenses : sequence of EllipsoidalLens
        Bodies in declared precedence order.
    resistivity_ohm_m : ndarray
        Final positive resistivity model shaped like ``base.grid``.
    lens_index : ndarray of int
        Assigned body index per cell, or ``-1`` outside assigned bodies.
    overlap_count : ndarray of int
        Number of geometric lens envelopes covering each cell.
    conflict_policy : {"error", "first", "last", "most_conductive", "most_resistive"}
        Policy used where body envelopes overlap.

    Examples
    --------
    >>> grid = GeologyGrid.regular_2d(nx=8, nz=6, dx_m=100, dz_m=50)
    >>> base = __import__(
    ...     "pycsamt.ai.geology", fromlist=[""]
    ... ).generate_layered_geology(
    ...     grid, [ElectricalLayer("earth", 100)], [], seed=0
    ... )
    >>> model = insert_lenses(
    ...     base, [EllipsoidalLens("lens", 400, 150, 200, 75, 10)]
    ... )
    >>> model.resistivity_ohm_m.shape
    (6, 8)
    """

    base: LayeredGeology
    lenses: tuple[EllipsoidalLens, ...]
    resistivity_ohm_m: np.ndarray
    lens_index: np.ndarray
    overlap_count: np.ndarray
    conflict_policy: str = "error"

    def __post_init__(self) -> None:
        if not isinstance(self.base, LayeredGeology):
            raise TypeError("base must be a LayeredGeology.")
        lenses = tuple(self.lenses)
        if not lenses or any(not isinstance(lens, EllipsoidalLens) for lens in lenses):
            raise ValueError("lenses must contain at least one EllipsoidalLens.")
        if len({lens.name for lens in lenses}) != len(lenses):
            raise ValueError("lens names must be unique.")
        for lens in lenses:
            lens.validate_grid(self.base.grid)
        if self.conflict_policy not in {
            "error",
            "first",
            "last",
            "most_conductive",
            "most_resistive",
        }:
            raise ValueError("unsupported conflict_policy.")
        resistivity = np.asarray(self.resistivity_ohm_m, dtype=float)
        indices = np.asarray(self.lens_index)
        overlap = np.asarray(self.overlap_count)
        shape = self.base.grid.shape
        if (
            resistivity.shape != shape
            or not np.all(np.isfinite(resistivity))
            or np.any(resistivity <= 0)
        ):
            raise ValueError(
                f"resistivity_ohm_m must be finite, positive, and shaped {shape}."
            )
        if (
            indices.shape != shape
            or not np.issubdtype(indices.dtype, np.integer)
            or np.any(indices < -1)
            or np.any(indices >= len(lenses))
        ):
            raise ValueError(
                f"lens_index must be an integer array shaped {shape} with valid body indices."
            )
        if (
            overlap.shape != shape
            or not np.issubdtype(overlap.dtype, np.integer)
            or np.any(overlap < 0)
        ):
            raise ValueError(
                f"overlap_count must be a non-negative integer array shaped {shape}."
            )
        expected_overlap = np.zeros(shape, dtype=np.int16)
        for lens in lenses:
            expected_overlap += lens.normalized_radius(self.base.grid) <= 1.0
        if not np.array_equal(overlap, expected_overlap):
            raise ValueError("overlap_count is inconsistent with lens geometry.")
        if np.any((indices >= 0) & (overlap == 0)):
            raise ValueError("lens_index assigns cells outside all lens envelopes.")
        if self.conflict_policy == "error" and np.any(overlap > 1):
            raise ValueError("error conflict policy cannot contain overlapping lenses.")
        expected_resistivity, expected_indices, _ = _compose_lenses(
            self.base, lenses, self.conflict_policy
        )
        if not np.array_equal(indices, expected_indices):
            raise ValueError(
                "lens_index is inconsistent with geometry and conflict policy."
            )
        if not np.allclose(resistivity, expected_resistivity, rtol=1e-12, atol=0.0):
            raise ValueError(
                "resistivity_ohm_m is inconsistent with base, lenses, and conflict policy."
            )
        object.__setattr__(self, "lenses", lenses)
        object.__setattr__(self, "resistivity_ohm_m", _readonly(resistivity))
        object.__setattr__(self, "lens_index", _readonly(indices, np.int16))
        object.__setattr__(self, "overlap_count", _readonly(overlap, np.int16))

    @property
    def model_hash(self) -> str:
        """Return a digest of base, lenses, final values, and overlap policy.

        Returns
        -------
        str
            Platform-stable SHA-256 digest.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
        >>> base = __import__(
        ...     "pycsamt.ai.geology", fromlist=[""]
        ... ).generate_layered_geology(
        ...     grid, [ElectricalLayer("earth", 100)], [], seed=0
        ... )
        >>> model = insert_lenses(
        ...     base, [EllipsoidalLens("lens", 200, 100, 100, 50, 10)]
        ... )
        >>> len(model.model_hash)
        64
        """
        digest = hashlib.sha256()
        digest.update(
            np.ascontiguousarray(self.resistivity_ohm_m, dtype="<f8").tobytes()
        )
        digest.update(np.ascontiguousarray(self.lens_index, dtype="<i2").tobytes())
        digest.update(np.ascontiguousarray(self.overlap_count, dtype="<i2").tobytes())
        digest.update(canonical_hash(self.provenance()).encode("ascii"))
        return digest.hexdigest()

    def lens_mask(self, lens: int | str) -> np.ndarray:
        """Return cells assigned to one lens after conflict resolution.

        Parameters
        ----------
        lens : int or str
            Zero-based lens index or exact body name.

        Returns
        -------
        ndarray of bool
            Read-only assigned-cell mask.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
        >>> base = __import__(
        ...     "pycsamt.ai.geology", fromlist=[""]
        ... ).generate_layered_geology(
        ...     grid, [ElectricalLayer("earth", 100)], [], seed=0
        ... )
        >>> model = insert_lenses(
        ...     base, [EllipsoidalLens("lens", 200, 100, 100, 50, 10)]
        ... )
        >>> model.lens_mask("lens").any()
        True
        """
        if isinstance(lens, str):
            names = tuple(item.name for item in self.lenses)
            try:
                index = names.index(lens)
            except ValueError as exc:
                raise KeyError(f"unknown lens {lens!r}.") from exc
        else:
            if not isinstance(lens, (int, np.integer)) or isinstance(lens, bool):
                raise TypeError("lens must be an integer index or exact name.")
            index = int(lens)
            if index < 0 or index >= len(self.lenses):
                raise IndexError(
                    f"lens index {index} is outside [0, {len(self.lenses)})."
                )
        return _readonly(self.lens_index == index, bool)

    def summary(self) -> dict[str, Any]:
        """Return JSON-compatible lens occupancy and resistivity diagnostics.

        Returns
        -------
        dict
            Assigned fractions, overlap fraction, resistivity range, and hashes.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
        >>> base = __import__(
        ...     "pycsamt.ai.geology", fromlist=[""]
        ... ).generate_layered_geology(
        ...     grid, [ElectricalLayer("earth", 100)], [], seed=0
        ... )
        >>> model = insert_lenses(
        ...     base, [EllipsoidalLens("lens", 200, 100, 100, 50, 10)]
        ... )
        >>> model.summary()["n_lenses"]
        1
        """
        return {
            "shape": list(self.base.grid.shape),
            "n_lenses": len(self.lenses),
            "assigned_fractions": {
                lens.name: float(np.mean(self.lens_index == index))
                for index, lens in enumerate(self.lenses)
            },
            "overlap_fraction": float(np.mean(self.overlap_count > 1)),
            "maximum_overlap": int(np.max(self.overlap_count)),
            "resistivity_range_ohm_m": [
                float(np.min(self.resistivity_ohm_m)),
                float(np.max(self.resistivity_ohm_m)),
            ],
            "base_model_hash": self.base.model_hash,
            "conflict_policy": self.conflict_policy,
        }

    def provenance(self) -> dict[str, Any]:
        """Return generation provenance without final cell arrays.

        Returns
        -------
        dict
            Base provenance/hash, lens definitions, and conflict policy.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
        >>> base = __import__(
        ...     "pycsamt.ai.geology", fromlist=[""]
        ... ).generate_layered_geology(
        ...     grid, [ElectricalLayer("earth", 100)], [], seed=0
        ... )
        >>> model = insert_lenses(
        ...     base, [EllipsoidalLens("lens", 200, 100, 100, 50, 10)]
        ... )
        >>> model.provenance()["conflict_policy"]
        'error'
        """
        return {
            "schema_version": 1,
            "base": self.base.provenance(),
            "base_model_hash": self.base.model_hash,
            "lenses": [lens.to_dict() for lens in self.lenses],
            "conflict_policy": self.conflict_policy,
        }

    def to_npz(self, path: str | Path) -> Path:
        """Persist base geology and lenses in one pickle-free NPZ archive.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination archive.

        Returns
        -------
        pathlib.Path
            Requested destination path.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
        >>> base = __import__(
        ...     "pycsamt.ai.geology", fromlist=[""]
        ... ).generate_layered_geology(
        ...     grid, [ElectricalLayer("earth", 100)], [], seed=0
        ... )
        >>> model = insert_lenses(
        ...     base, [EllipsoidalLens("lens", 200, 100, 100, 50, 10)]
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     path = model.to_npz(Path(directory) / "lenses.npz")
        ...     restored = LensGeology.from_npz(path)
        >>> restored.model_hash == model.model_hash
        True
        """
        target = Path(path)
        np.savez_compressed(
            target,
            base_interface_depth_m=self.base.interface_depth_m,
            base_layer_index=self.base.layer_index,
            base_resistivity_ohm_m=self.base.resistivity_ohm_m,
            resistivity_ohm_m=self.resistivity_ohm_m,
            lens_index=self.lens_index,
            overlap_count=self.overlap_count,
            provenance_json=np.array(json.dumps(self.provenance(), sort_keys=True)),
        )
        return target

    @classmethod
    def from_npz(cls, path: str | Path) -> LensGeology:
        """Load and validate a self-contained lens-geology archive.

        Parameters
        ----------
        path : str or pathlib.Path
            Archive written by :meth:`to_npz`.

        Returns
        -------
        LensGeology
            Immutable restored base and body model.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=100, dz_m=50)
        >>> base = __import__(
        ...     "pycsamt.ai.geology", fromlist=[""]
        ... ).generate_layered_geology(
        ...     grid, [ElectricalLayer("earth", 100)], [], seed=0
        ... )
        >>> model = insert_lenses(
        ...     base, [EllipsoidalLens("lens", 200, 100, 100, 50, 10)]
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     restored = LensGeology.from_npz(
        ...         model.to_npz(Path(directory) / "m.npz")
        ...     )
        >>> np.array_equal(restored.resistivity_ohm_m, model.resistivity_ohm_m)
        True
        """
        with np.load(Path(path), allow_pickle=False) as archive:
            state = json.loads(str(archive["provenance_json"].item()))
            if state.get("schema_version") != 1:
                raise ValueError("unsupported LensGeology schema version.")
            base_state = state["base"]
            base = LayeredGeology(
                GeologyGrid.from_dict(base_state["grid"]),
                tuple(ElectricalLayer.from_dict(item) for item in base_state["layers"]),
                archive["base_interface_depth_m"],
                archive["base_layer_index"],
                archive["base_resistivity_ohm_m"],
                base_state["seed"],
                base_state["boundary"],
                base_state["interface_policy"],
                base_state["minimum_thickness_m"],
                base_state["adjusted_interface_fraction"],
                base_state.get("generation_config", {}),
            )
            if base.model_hash != state["base_model_hash"]:
                raise ValueError("restored base model hash does not match provenance.")
            return cls(
                base,
                tuple(EllipsoidalLens.from_dict(item) for item in state["lenses"]),
                archive["resistivity_ohm_m"],
                archive["lens_index"],
                archive["overlap_count"],
                state["conflict_policy"],
            )


def insert_lenses(
    base: LayeredGeology,
    lenses: Sequence[EllipsoidalLens],
    *,
    conflict_policy: str = "error",
) -> LensGeology:
    """Insert electrical lenses into a layered model with explicit overlap rules.

    Parameters
    ----------
    base : LayeredGeology
        Stratigraphic resistivity model to modify immutably.
    lenses : sequence of EllipsoidalLens
        One or more uniquely named compatible bodies.
    conflict_policy : {"error", "first", "last", "most_conductive", "most_resistive"}, default="error"
        ``error`` rejects geometric overlap. ``first`` or ``last`` gives
        precedence by declaration order. The remaining policies select the
        lowest or highest blended cell resistivity.

    Returns
    -------
    LensGeology
        New immutable model; ``base`` remains unchanged.

    Raises
    ------
    ValueError
        If bodies are empty, incompatible, duplicate-named, non-intersecting,
        or overlap under the ``error`` policy.

    Examples
    --------
    >>> from pycsamt.ai.geology import generate_layered_geology
    >>> grid = GeologyGrid.regular_2d(nx=12, nz=8, dx_m=100, dz_m=50)
    >>> base = generate_layered_geology(
    ...     grid, [ElectricalLayer("earth", 100)], [], seed=0
    ... )
    >>> lens = EllipsoidalLens(
    ...     "conductor", 600, 200, 250, 100, 5, transition_fraction=0.2
    ... )
    >>> result = insert_lenses(base, [lens])
    >>> np.min(result.resistivity_ohm_m) < 100
    True
    >>> np.array_equal(base.resistivity_ohm_m, np.full(grid.shape, 100.0))
    True
    """
    if not isinstance(base, LayeredGeology):
        raise TypeError("base must be a LayeredGeology.")
    bodies = tuple(lenses)
    if not bodies or any(not isinstance(lens, EllipsoidalLens) for lens in bodies):
        raise ValueError("lenses must contain at least one EllipsoidalLens.")
    if len({lens.name for lens in bodies}) != len(bodies):
        raise ValueError("lens names must be unique.")
    if conflict_policy not in {
        "error",
        "first",
        "last",
        "most_conductive",
        "most_resistive",
    }:
        raise ValueError("unsupported conflict_policy.")
    result, indices, overlap = _compose_lenses(base, bodies, conflict_policy)
    return LensGeology(base, bodies, result, indices, overlap, conflict_policy)
