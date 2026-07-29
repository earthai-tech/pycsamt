# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Correlated stratigraphic interfaces and layered resistivity priors.

Layer interfaces are stored as depth below the grid datum and have shape
``(n_interface, nx)`` in 2-D or ``(n_interface, ny, nx)`` in 3-D.  Resulting
resistivity arrays follow the geology convention ``(nz, nx)`` or
``(nz, ny, nx)``.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from ..data.manifest import canonical_hash
from .fields import (
    GaussianCorrelation,
    GeologyGrid,
    generate_gaussian_field,
)

__all__ = [
    "ElectricalLayer",
    "LayeredGeology",
    "generate_layered_geology",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _freeze_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _freeze_json(item) for key, item in value.items()}
        )
    if isinstance(value, list):
        return tuple(_freeze_json(item) for item in value)
    return value


def _thaw_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _thaw_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_thaw_json(item) for item in value]
    return value


def _seed(value: int) -> int:
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError("seed must be an integer.")
    result = int(value)
    if result < 0 or result >= 2**64:
        raise ValueError("seed must be in [0, 2**64).")
    return result


def _child_seed(root: int, label: str) -> int:
    payload = f"pycsamt.layers\0{root}\0{label}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def _sequence(value: float | Sequence[float], count: int, name: str) -> np.ndarray:
    if np.isscalar(value):
        result = np.full(count, float(value), dtype=float)
    else:
        result = np.asarray(value, dtype=float)
        if result.shape != (count,):
            raise ValueError(f"{name} must be scalar or contain {count} values.")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values.")
    return result


def _horizontal_shape(grid: GeologyGrid) -> tuple[int, ...]:
    return (len(grid.x_m),) if grid.dimension == 2 else (len(grid.y_m), len(grid.x_m))


def _horizontal_gaussian(
    grid: GeologyGrid,
    correlation: GaussianCorrelation,
    seed: int,
    boundary: str,
) -> np.ndarray:
    shape = _horizontal_shape(grid)
    spacing = (
        (grid.spacing_m[-1],)
        if grid.dimension == 2
        else (grid.spacing_m[1], grid.spacing_m[2])
    )
    rng = np.random.default_rng(seed)
    if boundary == "periodic":
        white = rng.standard_normal(shape)
        starts = None
    else:
        core = rng.standard_normal(shape)
        pads = tuple((size // 2, size - size // 2) for size in shape)
        white = np.pad(core, pads, mode="reflect")
        starts = tuple(pad[0] for pad in pads)
    axes = [
        2 * np.pi * np.fft.fftfreq(size, d=step)
        for size, step in zip(white.shape, spacing)
    ]
    if grid.dimension == 2:
        (kx,) = np.meshgrid(axes[0], indexing="ij")
        exponent = np.square(kx * correlation.length_x_m)
    else:
        ky, kx = np.meshgrid(axes[0], axes[1], indexing="ij")
        angle = np.deg2rad(correlation.azimuth_deg)
        major = kx * np.cos(angle) + ky * np.sin(angle)
        minor = -kx * np.sin(angle) + ky * np.cos(angle)
        exponent = np.square(major * correlation.length_x_m) + np.square(
            minor * correlation.length_y_m
        )
    values = np.fft.ifftn(np.fft.fftn(white) * np.exp(-0.25 * exponent)).real
    if starts is not None:
        values = values[
            tuple(slice(start, start + size) for start, size in zip(starts, shape))
        ]
    deviation = float(np.std(values))
    if deviation <= np.finfo(float).eps:
        raise ValueError("interface correlation produced a numerically constant field.")
    return (values - np.mean(values)) / deviation


def _project_interfaces(
    surfaces: np.ndarray,
    top: float,
    bottom: float,
    minimum: float,
) -> tuple[np.ndarray, float]:
    projected = np.array(surfaces, copy=True)
    count = projected.shape[0]
    if count == 0:
        return projected, 0.0
    for index in range(count):
        lower = top + minimum * (index + 1)
        upper = bottom - minimum * (count - index)
        projected[index] = np.clip(projected[index], lower, upper)
    for index in range(1, count):
        projected[index] = np.maximum(projected[index], projected[index - 1] + minimum)
    for index in range(count - 2, -1, -1):
        projected[index] = np.minimum(projected[index], projected[index + 1] - minimum)
    changed = ~np.isclose(projected, surfaces, rtol=0.0, atol=1e-10)
    return projected, float(np.mean(changed))


@dataclass(frozen=True)
class ElectricalLayer:
    """Define the electrical distribution of one stratigraphic unit.

    Parameters
    ----------
    name : str
        Unique non-empty unit name.
    resistivity_ohm_m : float
        Positive median resistivity in ohm metres.
    log10_std : float, default=0.0
        Non-negative standard deviation of within-unit log10 resistivity.
    heterogeneity : GaussianCorrelation or None, optional
        Spatial correlation model required when ``log10_std`` is positive.
    resistivity_bounds_ohm_m : tuple of float or None, optional
        Positive inclusive lower/upper clipping bounds.

    Examples
    --------
    A homogeneous conductive cover:

    >>> cover = ElectricalLayer("conductive cover", 10.0)
    >>> cover.log10_resistivity
    1.0

    A heterogeneous basement:

    >>> basement = ElectricalLayer(
    ...     "basement",
    ...     1000.0,
    ...     log10_std=0.2,
    ...     heterogeneity=GaussianCorrelation(1000, 200),
    ...     resistivity_bounds_ohm_m=(100, 5000),
    ... )
    """

    name: str
    resistivity_ohm_m: float
    log10_std: float = 0.0
    heterogeneity: GaussianCorrelation | None = None
    resistivity_bounds_ohm_m: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        name = str(self.name).strip()
        resistivity = float(self.resistivity_ohm_m)
        deviation = float(self.log10_std)
        if not name:
            raise ValueError("name cannot be empty.")
        if not np.isfinite(resistivity) or resistivity <= 0:
            raise ValueError("resistivity_ohm_m must be finite and positive.")
        if not np.isfinite(deviation) or deviation < 0:
            raise ValueError("log10_std must be finite and non-negative.")
        if deviation > 0 and not isinstance(self.heterogeneity, GaussianCorrelation):
            raise ValueError("positive log10_std requires a GaussianCorrelation.")
        if (
            deviation == 0
            and self.heterogeneity is not None
            and not isinstance(self.heterogeneity, GaussianCorrelation)
        ):
            raise TypeError("heterogeneity must be a GaussianCorrelation or None.")
        bounds = self.resistivity_bounds_ohm_m
        if bounds is not None:
            if len(bounds) != 2:
                raise ValueError("resistivity_bounds_ohm_m must contain two values.")
            bounds = (float(bounds[0]), float(bounds[1]))
            if (
                not np.all(np.isfinite(bounds))
                or bounds[0] <= 0
                or bounds[0] >= bounds[1]
            ):
                raise ValueError(
                    "resistivity bounds must be finite, positive, and increasing."
                )
            if not bounds[0] <= resistivity <= bounds[1]:
                raise ValueError(
                    "median resistivity must lie within resistivity bounds."
                )
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "resistivity_ohm_m", resistivity)
        object.__setattr__(self, "log10_std", deviation)
        object.__setattr__(self, "resistivity_bounds_ohm_m", bounds)

    @property
    def log10_resistivity(self) -> float:
        """Return the median resistivity in log10 ohm metres.

        Returns
        -------
        float
            ``log10(resistivity_ohm_m)``.

        Examples
        --------
        >>> ElectricalLayer("unit", 100).log10_resistivity
        2.0
        """
        return float(np.log10(self.resistivity_ohm_m))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable unit definition.

        Returns
        -------
        dict
            Electrical distribution and optional correlation state.

        Examples
        --------
        >>> ElectricalLayer("unit", 100).to_dict()["resistivity_ohm_m"]
        100.0
        """
        return {
            "schema_version": 1,
            "name": self.name,
            "resistivity_ohm_m": self.resistivity_ohm_m,
            "log10_std": self.log10_std,
            "heterogeneity": None
            if self.heterogeneity is None
            else self.heterogeneity.to_dict(),
            "resistivity_bounds_ohm_m": None
            if self.resistivity_bounds_ohm_m is None
            else list(self.resistivity_bounds_ohm_m),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ElectricalLayer:
        """Restore a validated electrical layer.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        ElectricalLayer
            Immutable unit definition.

        Examples
        --------
        >>> unit = ElectricalLayer("unit", 50)
        >>> ElectricalLayer.from_dict(unit.to_dict()) == unit
        True
        """
        if data.get("schema_version", 1) != 1:
            raise ValueError("unsupported ElectricalLayer schema version.")
        correlation = data.get("heterogeneity")
        return cls(
            data["name"],
            data["resistivity_ohm_m"],
            data.get("log10_std", 0.0),
            None if correlation is None else GaussianCorrelation.from_dict(correlation),
            data.get("resistivity_bounds_ohm_m"),
        )


@dataclass(frozen=True)
class LayeredGeology:
    """Immutable discretized layered electrical geology.

    Parameters
    ----------
    grid : GeologyGrid
        Target model grid.
    layers : sequence of ElectricalLayer
        Units ordered shallowest to deepest.
    interface_depth_m : ndarray
        Interface surfaces shaped ``(n_layer - 1, *horizontal_shape)``.
    layer_index : ndarray of int
        Zero-based unit index for every model cell, shaped like ``grid``.
    resistivity_ohm_m : ndarray
        Positive cell resistivities shaped like ``grid``.
    seed : int
        Root generation seed.
    boundary, interface_policy : str
        Recorded generation policies.
    minimum_thickness_m : float
        Enforced minimum vertical separation.
    adjusted_interface_fraction : float
        Fraction of interface samples changed by projection.
    generation_config : mapping, optional
        Finite JSON-compatible requested interface means, relief amplitudes,
        and correlation models needed to regenerate the surfaces.

    Examples
    --------
    >>> grid = GeologyGrid.regular_2d(nx=8, nz=6, dx_m=100, dz_m=50)
    >>> layers = [
    ...     ElectricalLayer("cover", 10),
    ...     ElectricalLayer("basement", 1000),
    ... ]
    >>> model = generate_layered_geology(grid, layers, [150], seed=1)
    >>> model.resistivity_ohm_m.shape
    (6, 8)
    >>> model.n_layers
    2
    """

    grid: GeologyGrid
    layers: tuple[ElectricalLayer, ...]
    interface_depth_m: np.ndarray
    layer_index: np.ndarray
    resistivity_ohm_m: np.ndarray
    seed: int
    boundary: str
    interface_policy: str
    minimum_thickness_m: float
    adjusted_interface_fraction: float = 0.0
    generation_config: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.grid, GeologyGrid):
            raise TypeError("grid must be a GeologyGrid.")
        layers = tuple(self.layers)
        if not layers or any(
            not isinstance(layer, ElectricalLayer) for layer in layers
        ):
            raise ValueError("layers must contain ElectricalLayer objects.")
        if len({layer.name for layer in layers}) != len(layers):
            raise ValueError("layer names must be unique.")
        horizontal = _horizontal_shape(self.grid)
        interfaces = np.asarray(self.interface_depth_m, dtype=float)
        expected_interfaces = (len(layers) - 1,) + horizontal
        if interfaces.shape != expected_interfaces or not np.all(
            np.isfinite(interfaces)
        ):
            raise ValueError(
                f"interface_depth_m must be finite and shaped {expected_interfaces}."
            )
        indices = np.asarray(self.layer_index)
        if indices.shape != self.grid.shape or not np.issubdtype(
            indices.dtype, np.integer
        ):
            raise ValueError(
                f"layer_index must be an integer array shaped {self.grid.shape}."
            )
        if np.any(indices < 0) or np.any(indices >= len(layers)):
            raise ValueError("layer_index contains an unavailable unit index.")
        resistivity = np.asarray(self.resistivity_ohm_m, dtype=float)
        if (
            resistivity.shape != self.grid.shape
            or not np.all(np.isfinite(resistivity))
            or np.any(resistivity <= 0)
        ):
            raise ValueError(
                f"resistivity_ohm_m must be finite, positive, and shaped {self.grid.shape}."
            )
        minimum = float(self.minimum_thickness_m)
        adjusted = float(self.adjusted_interface_fraction)
        if not np.isfinite(minimum) or minimum <= 0:
            raise ValueError("minimum_thickness_m must be finite and positive.")
        if not np.isfinite(adjusted) or adjusted < 0 or adjusted > 1:
            raise ValueError("adjusted_interface_fraction must be in [0, 1].")
        if self.boundary not in {"reflect", "periodic"}:
            raise ValueError("boundary must be 'reflect' or 'periodic'.")
        if self.interface_policy not in {"raise", "project"}:
            raise ValueError("interface_policy must be 'raise' or 'project'.")
        generation = dict(self.generation_config)
        try:
            generation = json.loads(
                json.dumps(generation, sort_keys=True, allow_nan=False)
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "generation_config must contain finite JSON-serializable values."
            ) from exc
        top, bottom = self.grid.extent_m["z"]
        boundaries = np.concatenate(
            [
                np.full((1,) + horizontal, top),
                interfaces,
                np.full((1,) + horizontal, bottom),
            ],
            axis=0,
        )
        if np.any(np.diff(boundaries, axis=0) < minimum - 1e-9):
            raise ValueError(
                "interfaces violate model boundaries or minimum thickness."
            )
        depth = self.grid.z_m.reshape((len(self.grid.z_m),) + (1,) * len(horizontal))
        expected_indices = np.zeros(self.grid.shape, dtype=np.int16)
        for surface in interfaces:
            expected_indices += depth >= surface
        if not np.array_equal(indices, expected_indices):
            raise ValueError("layer_index is inconsistent with interface depths.")
        for index, layer in enumerate(layers):
            values = resistivity[indices == index]
            if values.size == 0:
                continue
            if layer.log10_std == 0 and not np.allclose(
                values, layer.resistivity_ohm_m, rtol=1e-12, atol=0.0
            ):
                raise ValueError(
                    f"homogeneous layer {layer.name!r} has inconsistent resistivity values."
                )
            if layer.resistivity_bounds_ohm_m is not None and (
                np.any(values < layer.resistivity_bounds_ohm_m[0])
                or np.any(values > layer.resistivity_bounds_ohm_m[1])
            ):
                raise ValueError(
                    f"layer {layer.name!r} violates its resistivity bounds."
                )
        object.__setattr__(self, "layers", layers)
        object.__setattr__(self, "interface_depth_m", _readonly(interfaces))
        object.__setattr__(self, "layer_index", _readonly(indices, np.int16))
        object.__setattr__(self, "resistivity_ohm_m", _readonly(resistivity))
        object.__setattr__(self, "seed", _seed(self.seed))
        object.__setattr__(self, "minimum_thickness_m", minimum)
        object.__setattr__(self, "adjusted_interface_fraction", adjusted)
        object.__setattr__(self, "generation_config", _freeze_json(generation))

    @property
    def n_layers(self) -> int:
        """Return the number of stratigraphic units.

        Returns
        -------
        int
            Number of electrical layers.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> model = generate_layered_geology(
        ...     grid,
        ...     [ElectricalLayer("a", 1), ElectricalLayer("b", 2)],
        ...     [1],
        ...     seed=0,
        ... )
        >>> model.n_layers
        2
        """
        return len(self.layers)

    @property
    def model_hash(self) -> str:
        """Return a platform-stable digest of model values and provenance.

        Returns
        -------
        str
            SHA-256 digest.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> model = generate_layered_geology(
        ...     grid,
        ...     [ElectricalLayer("a", 1), ElectricalLayer("b", 2)],
        ...     [1],
        ...     seed=0,
        ... )
        >>> len(model.model_hash)
        64
        """
        digest = hashlib.sha256()
        digest.update(
            np.ascontiguousarray(self.resistivity_ohm_m, dtype="<f8").tobytes()
        )
        digest.update(np.ascontiguousarray(self.layer_index, dtype="<i2").tobytes())
        digest.update(
            np.ascontiguousarray(self.interface_depth_m, dtype="<f8").tobytes()
        )
        digest.update(canonical_hash(self.provenance()).encode("ascii"))
        return digest.hexdigest()

    def layer_mask(self, layer: int | str) -> np.ndarray:
        """Return a read-only Boolean mask for one unit.

        Parameters
        ----------
        layer : int or str
            Zero-based index or exact layer name.

        Returns
        -------
        ndarray of bool
            Mask shaped like the geological grid.

        Raises
        ------
        KeyError
            If a named unit is unavailable.
        IndexError
            If a numeric index is outside the layer range.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> model = generate_layered_geology(
        ...     grid,
        ...     [ElectricalLayer("a", 1), ElectricalLayer("b", 2)],
        ...     [1],
        ...     seed=0,
        ... )
        >>> model.layer_mask("a").shape
        (2, 2)
        """
        if isinstance(layer, str):
            names = tuple(item.name for item in self.layers)
            try:
                index = names.index(layer)
            except ValueError as exc:
                raise KeyError(f"unknown layer {layer!r}.") from exc
        else:
            if not isinstance(layer, (int, np.integer)) or isinstance(layer, bool):
                raise TypeError("layer must be an integer index or exact name.")
            index = int(layer)
            if index < 0 or index >= self.n_layers:
                raise IndexError(
                    f"layer index {index} is outside [0, {self.n_layers})."
                )
        return _readonly(self.layer_index == index, bool)

    def interface(self, index: int) -> np.ndarray:
        """Return one read-only interface-depth surface.

        Parameters
        ----------
        index : int
            Zero-based interface between layers ``index`` and ``index + 1``.

        Returns
        -------
        ndarray
            Shape ``(nx,)`` in 2-D or ``(ny, nx)`` in 3-D.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=3, nz=2, dx_m=1, dz_m=1)
        >>> model = generate_layered_geology(
        ...     grid,
        ...     [ElectricalLayer("a", 1), ElectricalLayer("b", 2)],
        ...     [1],
        ...     seed=0,
        ... )
        >>> model.interface(0).shape
        (3,)
        """
        if (
            not isinstance(index, int)
            or isinstance(index, bool)
            or index < 0
            or index >= self.n_layers - 1
        ):
            raise IndexError(f"interface index must be in [0, {self.n_layers - 1}).")
        return self.interface_depth_m[index]

    def summary(self) -> dict[str, Any]:
        """Return compact JSON-compatible geological diagnostics.

        Returns
        -------
        dict
            Shape, unit fractions, interface ranges, resistivity range, and
            interface-adjustment fraction.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> model = generate_layered_geology(
        ...     grid,
        ...     [ElectricalLayer("a", 1), ElectricalLayer("b", 2)],
        ...     [1],
        ...     seed=0,
        ... )
        >>> model.summary()["n_layers"]
        2
        """
        return {
            "shape": list(self.grid.shape),
            "n_layers": self.n_layers,
            "layer_fractions": {
                layer.name: float(np.mean(self.layer_index == index))
                for index, layer in enumerate(self.layers)
            },
            "interface_ranges_m": [
                [float(np.min(surface)), float(np.max(surface))]
                for surface in self.interface_depth_m
            ],
            "resistivity_range_ohm_m": [
                float(np.min(self.resistivity_ohm_m)),
                float(np.max(self.resistivity_ohm_m)),
            ],
            "adjusted_interface_fraction": self.adjusted_interface_fraction,
            "generation_config": _thaw_json(self.generation_config),
        }

    def provenance(self) -> dict[str, Any]:
        """Return generation provenance without cell arrays.

        Returns
        -------
        dict
            Grid, layer definitions, seed, and interface policies.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> model = generate_layered_geology(
        ...     grid,
        ...     [ElectricalLayer("a", 1), ElectricalLayer("b", 2)],
        ...     [1],
        ...     seed=0,
        ... )
        >>> model.provenance()["seed"]
        0
        """
        return {
            "schema_version": 1,
            "grid": self.grid.to_dict(),
            "layers": [layer.to_dict() for layer in self.layers],
            "seed": self.seed,
            "boundary": self.boundary,
            "interface_policy": self.interface_policy,
            "minimum_thickness_m": self.minimum_thickness_m,
            "adjusted_interface_fraction": self.adjusted_interface_fraction,
        }

    def to_npz(self, path: str | Path) -> Path:
        """Persist the layered model in a pickle-free compressed archive.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination NPZ path.

        Returns
        -------
        pathlib.Path
            Requested destination.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> model = generate_layered_geology(
        ...     grid,
        ...     [ElectricalLayer("a", 1), ElectricalLayer("b", 2)],
        ...     [1],
        ...     seed=0,
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     path = model.to_npz(Path(directory) / "layers.npz")
        ...     restored = LayeredGeology.from_npz(path)
        >>> restored.model_hash == model.model_hash
        True
        """
        target = Path(path)
        np.savez_compressed(
            target,
            interface_depth_m=self.interface_depth_m,
            layer_index=self.layer_index,
            resistivity_ohm_m=self.resistivity_ohm_m,
            provenance_json=np.array(json.dumps(self.provenance(), sort_keys=True)),
        )
        return target

    @classmethod
    def from_npz(cls, path: str | Path) -> LayeredGeology:
        """Load and validate a layered model without enabling pickle.

        Parameters
        ----------
        path : str or pathlib.Path
            Archive written by :meth:`to_npz`.

        Returns
        -------
        LayeredGeology
            Immutable restored model.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> model = generate_layered_geology(
        ...     GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1),
        ...     [ElectricalLayer("a", 1), ElectricalLayer("b", 2)],
        ...     [1],
        ...     seed=0,
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     path = model.to_npz(Path(directory) / "m.npz")
        ...     restored = LayeredGeology.from_npz(path)
        >>> np.array_equal(restored.layer_index, model.layer_index)
        True
        """
        with np.load(Path(path), allow_pickle=False) as archive:
            state = json.loads(str(archive["provenance_json"].item()))
            if state.get("schema_version") != 1:
                raise ValueError("unsupported LayeredGeology schema version.")
            return cls(
                GeologyGrid.from_dict(state["grid"]),
                tuple(ElectricalLayer.from_dict(item) for item in state["layers"]),
                archive["interface_depth_m"],
                archive["layer_index"],
                archive["resistivity_ohm_m"],
                state["seed"],
                state["boundary"],
                state["interface_policy"],
                state["minimum_thickness_m"],
                state["adjusted_interface_fraction"],
                state.get("generation_config", {}),
            )


def generate_layered_geology(
    grid: GeologyGrid,
    layers: Sequence[ElectricalLayer],
    mean_interface_depth_m: Sequence[float],
    *,
    seed: int,
    interface_relief_std_m: float | Sequence[float] = 0.0,
    interface_correlation: GaussianCorrelation
    | Sequence[GaussianCorrelation]
    | None = None,
    minimum_thickness_m: float | None = None,
    interface_policy: str = "project",
    boundary: str = "reflect",
) -> LayeredGeology:
    """Generate a correlated 2-D section or 3-D layered electrical volume.

    Parameters
    ----------
    grid : GeologyGrid
        Regular target grid.
    layers : sequence of ElectricalLayer
        Units ordered shallowest to deepest. At least one is required.
    mean_interface_depth_m : sequence of float
        One mean depth for every adjacent layer pair.
    seed : int
        Explicit root seed. Labeled child seeds independently drive interfaces
        and within-layer heterogeneity.
    interface_relief_std_m : float or sequence, default=0.0
        Interface depth standard deviations in metres.
    interface_correlation : GaussianCorrelation, sequence, or None, optional
        Horizontal correlation model(s). Required for every interface with
        positive relief. A single model is shared across interfaces.
    minimum_thickness_m : float or None, optional
        Minimum separation between top boundary, interfaces, and bottom
        boundary. Default is one vertical cell.
    interface_policy : {"project", "raise"}, default="project"
        Project invalid/crossing surfaces into the feasible domain or reject
        the realization. Projection fraction is recorded on the result.
    boundary : {"reflect", "periodic"}, default="reflect"
        Boundary policy for all correlated fields.

    Returns
    -------
    LayeredGeology
        Immutable interfaces, unit indices, resistivity, and provenance.

    Raises
    ------
    ValueError
        If layer/interface counts disagree, thickness is infeasible, relief
        lacks correlation, a 3-D correlation lacks ``length_y_m``, or strict
        interface validation fails.

    Examples
    --------
    Laterally varying three-layer 2-D model:

    >>> grid = GeologyGrid.regular_2d(nx=24, nz=16, dx_m=100, dz_m=50)
    >>> units = [
    ...     ElectricalLayer("cover", 20),
    ...     ElectricalLayer("sediments", 100),
    ...     ElectricalLayer("basement", 1000),
    ... ]
    >>> model = generate_layered_geology(
    ...     grid,
    ...     units,
    ...     [200, 500],
    ...     seed=5,
    ...     interface_relief_std_m=[30, 50],
    ...     interface_correlation=GaussianCorrelation(600, 100),
    ... )
    >>> model.interface_depth_m.shape
    (2, 24)
    >>> set(np.unique(model.layer_index)) <= {0, 1, 2}
    True
    """
    if not isinstance(grid, GeologyGrid):
        raise TypeError("grid must be a GeologyGrid.")
    units = tuple(layers)
    if not units or any(not isinstance(layer, ElectricalLayer) for layer in units):
        raise ValueError("layers must contain at least one ElectricalLayer.")
    if len({layer.name for layer in units}) != len(units):
        raise ValueError("layer names must be unique.")
    interface_count = len(units) - 1
    means = np.asarray(mean_interface_depth_m, dtype=float)
    if means.shape != (interface_count,) or not np.all(np.isfinite(means)):
        raise ValueError(
            f"mean_interface_depth_m must contain {interface_count} finite values."
        )
    if interface_count and not np.all(np.diff(means) > 0):
        raise ValueError("mean interface depths must be strictly increasing.")
    relief = _sequence(
        interface_relief_std_m, interface_count, "interface_relief_std_m"
    )
    if np.any(relief < 0):
        raise ValueError("interface relief standard deviations must be non-negative.")
    if boundary not in {"reflect", "periodic"}:
        raise ValueError("boundary must be 'reflect' or 'periodic'.")
    if interface_policy not in {"project", "raise"}:
        raise ValueError("interface_policy must be 'project' or 'raise'.")
    seed = _seed(seed)
    dz = grid.spacing_m[0]
    minimum = dz if minimum_thickness_m is None else float(minimum_thickness_m)
    if not np.isfinite(minimum) or minimum <= 0:
        raise ValueError("minimum_thickness_m must be finite and positive.")
    top, bottom = grid.extent_m["z"]
    if minimum * len(units) > bottom - top + 1e-12:
        raise ValueError(
            "minimum thickness is infeasible for the grid depth and layer count."
        )

    if interface_count == 0:
        correlations: tuple[GaussianCorrelation | None, ...] = ()
    elif isinstance(interface_correlation, GaussianCorrelation):
        correlations = (interface_correlation,) * interface_count
    elif interface_correlation is None:
        correlations = (None,) * interface_count
    else:
        correlations = tuple(interface_correlation)
        if len(correlations) != interface_count or any(
            not isinstance(item, GaussianCorrelation) for item in correlations
        ):
            raise ValueError(
                f"interface_correlation must contain {interface_count} GaussianCorrelation objects."
            )
    horizontal = _horizontal_shape(grid)
    surfaces = np.empty((interface_count,) + horizontal, dtype=float)
    for index in range(interface_count):
        if relief[index] == 0:
            surfaces[index].fill(means[index])
        else:
            correlation = correlations[index]
            if correlation is None:
                raise ValueError(
                    "positive interface relief requires interface_correlation."
                )
            correlation.validate_grid(grid)
            surfaces[index] = means[index] + relief[index] * _horizontal_gaussian(
                grid,
                correlation,
                _child_seed(seed, f"interface/{index}"),
                boundary,
            )
    projected, adjusted = _project_interfaces(surfaces, top, bottom, minimum)
    if interface_policy == "raise" and adjusted > 0:
        raise ValueError(
            "generated interfaces violate boundaries or minimum thickness."
        )
    if interface_policy == "project":
        surfaces = projected

    depth_shape = (len(grid.z_m),) + (1,) * len(horizontal)
    depth = grid.z_m.reshape(depth_shape)
    indices = np.zeros(grid.shape, dtype=np.int16)
    for surface in surfaces:
        indices += depth >= surface

    resistivity = np.empty(grid.shape, dtype=float)
    for index, layer in enumerate(units):
        if layer.log10_std > 0:
            layer.heterogeneity.validate_grid(grid)
            field = generate_gaussian_field(
                grid,
                layer.heterogeneity,
                seed=_child_seed(seed, f"heterogeneity/{index}"),
                boundary=boundary,
            ).values
            values = np.power(10.0, layer.log10_resistivity + layer.log10_std * field)
            if layer.resistivity_bounds_ohm_m is not None:
                values = np.clip(values, *layer.resistivity_bounds_ohm_m)
        else:
            values = layer.resistivity_ohm_m
        mask = indices == index
        resistivity[mask] = values if np.isscalar(values) else values[mask]
    return LayeredGeology(
        grid,
        units,
        surfaces,
        indices,
        resistivity,
        seed,
        boundary,
        interface_policy,
        minimum,
        adjusted,
        {
            "mean_interface_depth_m": means.tolist(),
            "interface_relief_std_m": relief.tolist(),
            "interface_correlation": [
                None if item is None else item.to_dict() for item in correlations
            ],
        },
    )
