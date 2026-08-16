# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Seeded geological families for two-dimensional EM benchmarks."""

from __future__ import annotations

import hashlib
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from ..data.manifest import canonical_hash
from .fields import GaussianCorrelation, GeologyGrid, generate_gaussian_field
from .layers import ElectricalLayer, generate_layered_geology
from .lenses import EllipsoidalLens, insert_lenses

__all__ = [
    "ID_BENCHMARK_FAMILIES",
    "OOD_BENCHMARK_FAMILIES",
    "BenchmarkGeology",
    "generate_benchmark_geology",
]

ID_BENCHMARK_FAMILIES = (
    "layered",
    "intrusion",
    "intrusion_halo",
    "dipping_fault",
    "multiple_body",
    "correlated_heterogeneous",
)

OOD_BENCHMARK_FAMILIES = (
    "extreme_dip",
    "deep_small_conductor",
    "overlapping_multibody",
    "high_contrast_halo",
    "anisotropic_correlation",
    "rugged_interfaces",
)


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    """Return a read-only array copy."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _seed(value: int) -> int:
    """Validate one unsigned seed."""
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError("seed must be an integer")
    result = int(value)
    if result < 0 or result >= 2**64:
        raise ValueError("seed must lie in [0, 2**64)")
    return result


def _range(config: Mapping[str, Any], name: str) -> tuple[float, float]:
    """Return a finite increasing numerical range."""
    value = config[name]
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"{name} must contain two values")
    low, high = float(value[0]), float(value[1])
    if not np.isfinite(low) or not np.isfinite(high) or high < low:
        raise ValueError(f"{name} must satisfy finite low <= high")
    return low, high


def _uniform(
    rng: np.random.Generator,
    config: Mapping[str, Any],
    name: str,
) -> float:
    """Draw one uniform value from a configured range."""
    low, high = _range(config, name)
    return float(rng.uniform(low, high))


def _log_uniform(
    rng: np.random.Generator,
    config: Mapping[str, Any],
    name: str,
) -> float:
    """Draw one positive value uniformly in log10 space."""
    low, high = _range(config, name)
    if low <= 0:
        raise ValueError(f"{name} must be positive for log-uniform sampling")
    return float(np.power(10.0, rng.uniform(np.log10(low), np.log10(high))))


def _integer(
    rng: np.random.Generator,
    config: Mapping[str, Any],
    name: str,
) -> int:
    """Draw one integer from an inclusive configured range."""
    low, high = _range(config, name)
    if not low.is_integer() or not high.is_integer():
        raise ValueError(f"{name} must contain integer bounds")
    return int(rng.integers(int(low), int(high) + 1))


def _child_seed(rng: np.random.Generator) -> int:
    """Draw a NumPy-compatible child seed."""
    return int(rng.integers(0, 2**32))


@dataclass(frozen=True)
class BenchmarkGeology:
    """One generated benchmark resistivity model with provenance.

    Parameters
    ----------
    family : str
        Frozen geological family name.
    regime : {"id", "ood"}
        Whether the realization belongs to the in-distribution or structural
        out-of-distribution family set.
    seed : int
        Root geology seed.
    grid : GeologyGrid
        Shared regular 2-D grid.
    resistivity_ohm_m : ndarray
        Positive finite resistivity shaped like ``grid``.
    parameters : mapping
        Concrete sampled geological parameters.
    """

    family: str
    regime: str
    seed: int
    grid: GeologyGrid
    resistivity_ohm_m: np.ndarray
    parameters: Mapping[str, Any]

    def __post_init__(self) -> None:
        family = str(self.family).strip()
        regime = str(self.regime).strip().lower()
        expected = (
            ID_BENCHMARK_FAMILIES
            if regime == "id"
            else OOD_BENCHMARK_FAMILIES
        )
        if regime not in {"id", "ood"} or family not in expected:
            raise ValueError("family and regime are incompatible")
        if not isinstance(self.grid, GeologyGrid) or self.grid.dimension != 2:
            raise TypeError("grid must be a 2-D GeologyGrid")
        resistivity = np.asarray(self.resistivity_ohm_m, dtype=float)
        if (
            resistivity.shape != self.grid.shape
            or not np.all(np.isfinite(resistivity))
            or np.any(resistivity <= 0)
        ):
            raise ValueError(
                "resistivity_ohm_m must be positive, finite, and shaped "
                f"{self.grid.shape}"
            )
        parameters = dict(self.parameters)
        canonical_hash(parameters)
        object.__setattr__(self, "family", family)
        object.__setattr__(self, "regime", regime)
        object.__setattr__(self, "seed", _seed(self.seed))
        object.__setattr__(
            self, "resistivity_ohm_m", _readonly(resistivity)
        )
        object.__setattr__(self, "parameters", parameters)

    @property
    def model_hash(self) -> str:
        """Return a SHA-256 covering model values and provenance.

        Returns
        -------
        str
            Platform-stable hexadecimal model digest.
        """
        digest = hashlib.sha256()
        digest.update(
            np.ascontiguousarray(self.resistivity_ohm_m, dtype="<f8")
        )
        digest.update(canonical_hash(self.to_dict()).encode("ascii"))
        return digest.hexdigest()

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible generation provenance.

        Returns
        -------
        dict
            Family, regime, seed, grid, and concrete sampled parameters.
        """
        return {
            "schema_version": 1,
            "family": self.family,
            "regime": self.regime,
            "seed": self.seed,
            "grid": self.grid.to_dict(),
            "parameters": dict(self.parameters),
        }


def _base_layered(
    grid: GeologyGrid,
    rng: np.random.Generator,
    config: Mapping[str, Any],
    *,
    relief_range: tuple[float, float],
    length_x_range: tuple[float, float] | None = None,
) -> tuple[Any, dict[str, Any]]:
    """Generate a randomized three-unit background."""
    units = (
        ElectricalLayer(
            "cover", _log_uniform(rng, config, "cover_resistivity_ohm_m")
        ),
        ElectricalLayer(
            "host", _log_uniform(rng, config, "host_resistivity_ohm_m")
        ),
        ElectricalLayer(
            "basement",
            _log_uniform(rng, config, "basement_resistivity_ohm_m"),
        ),
    )
    first = _uniform(rng, config, "first_interface_depth_m")
    second = _uniform(rng, config, "second_interface_depth_m")
    relief = [
        float(rng.uniform(*relief_range)),
        float(rng.uniform(*relief_range)),
    ]
    length_range = (
        _range(config, "interface_length_x_m")
        if length_x_range is None
        else length_x_range
    )
    correlation = GaussianCorrelation(
        float(rng.uniform(*length_range)),
        float(config["interface_length_z_m"]),
    )
    child_seed = _child_seed(rng)
    model = generate_layered_geology(
        grid,
        units,
        [first, second],
        seed=child_seed,
        interface_relief_std_m=relief,
        interface_correlation=correlation,
        minimum_thickness_m=float(config["minimum_thickness_m"]),
        interface_policy="project",
    )
    parameters = {
        "layer_resistivity_ohm_m": [
            layer.resistivity_ohm_m for layer in units
        ],
        "mean_interface_depth_m": [first, second],
        "interface_relief_std_m": relief,
        "interface_correlation": correlation.to_dict(),
        "layer_seed": child_seed,
        "adjusted_interface_fraction": model.adjusted_interface_fraction,
        "layered_model_hash": model.model_hash,
    }
    return model, parameters


def _fraction_x(
    grid: GeologyGrid,
    rng: np.random.Generator,
    config: Mapping[str, Any],
) -> float:
    """Draw an x coordinate as a fraction of grid width."""
    low, high = _range(config, "center_x_fraction")
    edge_low, edge_high = grid.extent_m["x"]
    return float(edge_low + rng.uniform(low, high) * (edge_high - edge_low))


def _fault_model(
    grid: GeologyGrid,
    base: Any,
    dip_deg: float,
    trace_x_m: float,
    throw_cells: int,
) -> np.ndarray:
    """Displace layered cells across one dipping fault plane."""
    source = np.asarray(base.resistivity_ohm_m)
    result = np.array(source, copy=True)
    tangent = np.tan(np.deg2rad(dip_deg))
    for z_index, depth_m in enumerate(grid.z_m):
        plane_x_m = trace_x_m + depth_m / tangent
        hanging_wall = grid.x_m >= plane_x_m
        source_index = min(z_index + throw_cells, source.shape[0] - 1)
        result[z_index, hanging_wall] = source[source_index, hanging_wall]
    return result


def _id_geology(
    grid: GeologyGrid,
    family: str,
    rng: np.random.Generator,
    base_config: Mapping[str, Any],
    family_config: Mapping[str, Any],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate one in-distribution family."""
    if family == "correlated_heterogeneous":
        correlation = GaussianCorrelation(
            _uniform(rng, family_config, "correlation_length_x_m"),
            _uniform(rng, family_config, "correlation_length_z_m"),
        )
        field_seed = _child_seed(rng)
        field = generate_gaussian_field(grid, correlation, seed=field_seed)
        mean = _uniform(rng, family_config, "log10_resistivity_mean")
        deviation = _uniform(rng, family_config, "log10_resistivity_std")
        values = np.power(10.0, mean + deviation * field.values)
        return values, {
            "field_seed": field_seed,
            "field_hash": field.field_hash,
            "correlation": correlation.to_dict(),
            "log10_resistivity_mean": mean,
            "log10_resistivity_std": deviation,
        }

    relief_range = (
        _range(family_config, "interface_relief_std_m")
        if family == "layered"
        else (15.0, 75.0)
    )
    base, parameters = _base_layered(
        grid, rng, base_config, relief_range=relief_range
    )
    if family == "layered":
        return base.resistivity_ohm_m, parameters

    if family in {"intrusion", "intrusion_halo"}:
        center_x = _fraction_x(grid, rng, family_config)
        center_z = _uniform(rng, family_config, "center_depth_m")
        dip = _uniform(rng, family_config, "dip_deg")
        if family == "intrusion":
            lens = EllipsoidalLens(
                "intrusion",
                center_x,
                center_z,
                _uniform(rng, family_config, "radius_x_m"),
                _uniform(rng, family_config, "radius_z_m"),
                _log_uniform(rng, family_config, "resistivity_ohm_m"),
                dip_deg=dip,
                transition_fraction=_uniform(
                    rng, family_config, "transition_fraction"
                ),
            )
            model = insert_lenses(base, [lens])
            parameters["lenses"] = [lens.to_dict()]
            return model.resistivity_ohm_m, parameters
        halo_x = _uniform(rng, family_config, "halo_radius_x_m")
        halo_z = _uniform(rng, family_config, "halo_radius_z_m")
        fraction = _uniform(rng, family_config, "core_radius_fraction")
        halo = EllipsoidalLens(
            "alteration_halo",
            center_x,
            center_z,
            halo_x,
            halo_z,
            _log_uniform(rng, family_config, "halo_resistivity_ohm_m"),
            dip_deg=dip,
            transition_fraction=0.15,
        )
        core = EllipsoidalLens(
            "intrusion_core",
            center_x,
            center_z,
            halo_x * fraction,
            halo_z * fraction,
            _log_uniform(rng, family_config, "core_resistivity_ohm_m"),
            dip_deg=dip,
            transition_fraction=0.08,
        )
        model = insert_lenses(base, [halo, core], conflict_policy="last")
        parameters["lenses"] = [halo.to_dict(), core.to_dict()]
        return model.resistivity_ohm_m, parameters

    if family == "dipping_fault":
        dip = _uniform(rng, family_config, "dip_deg")
        x_low, x_high = grid.extent_m["x"]
        fraction = rng.uniform(*_range(family_config, "trace_x_fraction"))
        trace = float(x_low + fraction * (x_high - x_low))
        throw = _integer(rng, family_config, "throw_cells")
        parameters["fault"] = {
            "dip_deg": dip,
            "surface_trace_x_m": trace,
            "throw_cells": throw,
            "throw_m": throw * grid.spacing_m[0],
        }
        return _fault_model(grid, base, dip, trace, throw), parameters

    count = _integer(rng, family_config, "body_count")
    x_edges = np.linspace(
        grid.extent_m["x"][0] + 600.0,
        grid.extent_m["x"][1] - 600.0,
        count,
    )
    lenses = []
    for index, centre in enumerate(x_edges):
        conductive = index % 3 != 1
        name = "conductive_resistivity_ohm_m"
        if not conductive:
            name = "resistive_resistivity_ohm_m"
        lenses.append(
            EllipsoidalLens(
                f"body_{index:02d}",
                float(centre + rng.uniform(-90.0, 90.0)),
                float(rng.uniform(450.0, 1200.0)),
                _uniform(rng, family_config, "radius_x_m"),
                _uniform(rng, family_config, "radius_z_m"),
                _log_uniform(rng, family_config, name),
                dip_deg=float(rng.uniform(-25.0, 25.0)),
            )
        )
    model = insert_lenses(base, lenses, conflict_policy="error")
    parameters["lenses"] = [lens.to_dict() for lens in lenses]
    return model.resistivity_ohm_m, parameters


def _ood_geology(
    grid: GeologyGrid,
    family: str,
    rng: np.random.Generator,
    base_config: Mapping[str, Any],
    family_config: Mapping[str, Any],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate one structural OOD family."""
    if family == "anisotropic_correlation":
        vertical = _uniform(rng, family_config, "correlation_length_z_m")
        ratio = _uniform(rng, family_config, "horizontal_vertical_ratio")
        correlation = GaussianCorrelation(vertical * ratio, vertical)
        field_seed = _child_seed(rng)
        field = generate_gaussian_field(grid, correlation, seed=field_seed)
        mean = _uniform(rng, family_config, "log10_resistivity_mean")
        deviation = _uniform(rng, family_config, "log10_resistivity_std")
        return np.power(10.0, mean + deviation * field.values), {
            "field_seed": field_seed,
            "field_hash": field.field_hash,
            "correlation": correlation.to_dict(),
            "horizontal_vertical_ratio": ratio,
            "log10_resistivity_mean": mean,
            "log10_resistivity_std": deviation,
        }

    relief_range = (20.0, 70.0)
    length_range = None
    if family == "rugged_interfaces":
        relief_range = _range(family_config, "interface_relief_std_m")
        length_range = _range(family_config, "interface_length_x_m")
    base, parameters = _base_layered(
        grid,
        rng,
        base_config,
        relief_range=relief_range,
        length_x_range=length_range,
    )
    if family == "rugged_interfaces":
        return base.resistivity_ohm_m, parameters

    if family == "extreme_dip":
        dip = _uniform(rng, family_config, "dip_deg")
        x_low, x_high = grid.extent_m["x"]
        fraction = rng.uniform(*_range(family_config, "trace_x_fraction"))
        trace = float(x_low + fraction * (x_high - x_low))
        throw = _integer(rng, family_config, "throw_cells")
        parameters["fault"] = {
            "dip_deg": dip,
            "surface_trace_x_m": trace,
            "throw_cells": throw,
        }
        return _fault_model(grid, base, dip, trace, throw), parameters

    if family == "deep_small_conductor":
        lens = EllipsoidalLens(
            "deep_small_conductor",
            _fraction_x(grid, rng, family_config),
            _uniform(rng, family_config, "center_depth_m"),
            _uniform(rng, family_config, "radius_x_m"),
            _uniform(rng, family_config, "radius_z_m"),
            _log_uniform(rng, family_config, "resistivity_ohm_m"),
        )
        snapped = False
        if not np.any(lens.normalized_radius(grid) <= 1.0):
            lens = replace(
                lens,
                center_x_m=float(
                    grid.x_m[
                        np.argmin(np.abs(grid.x_m - lens.center_x_m))
                    ]
                ),
                center_z_m=float(
                    grid.z_m[
                        np.argmin(np.abs(grid.z_m - lens.center_z_m))
                    ]
                ),
            )
            snapped = True
        model = insert_lenses(base, [lens])
        parameters["lenses"] = [lens.to_dict()]
        if snapped:
            parameters["center_snapped_to_grid"] = True
        return model.resistivity_ohm_m, parameters

    if family in {"high_contrast_halo", "overlapping_multibody"}:
        if family == "high_contrast_halo":
            center_x = _fraction_x(grid, rng, family_config)
            center_z = _uniform(rng, family_config, "center_depth_m")
            halo_x = _uniform(rng, family_config, "halo_radius_x_m")
            halo_z = _uniform(rng, family_config, "halo_radius_z_m")
            fraction = _uniform(rng, family_config, "core_radius_fraction")
            dip = _uniform(rng, family_config, "dip_deg")
            lenses = (
                EllipsoidalLens(
                    "ood_halo",
                    center_x,
                    center_z,
                    halo_x,
                    halo_z,
                    _log_uniform(
                        rng, family_config, "halo_resistivity_ohm_m"
                    ),
                    dip_deg=dip,
                    transition_fraction=0.15,
                ),
                EllipsoidalLens(
                    "ood_core",
                    center_x,
                    center_z,
                    halo_x * fraction,
                    halo_z * fraction,
                    _log_uniform(
                        rng, family_config, "core_resistivity_ohm_m"
                    ),
                    dip_deg=dip,
                ),
            )
        else:
            count = _integer(rng, family_config, "body_count")
            center_x = float(np.mean(grid.extent_m["x"]))
            center_z = float(np.mean(grid.extent_m["z"]))
            x_bounds = (float(grid.x_m[0]), float(grid.x_m[-1]))
            z_bounds = (float(grid.z_m[0]), float(grid.z_m[-1]))
            lenses = tuple(
                EllipsoidalLens(
                    f"overlap_{index:02d}",
                    float(
                        np.clip(
                            center_x
                            + rng.normal(
                                0,
                                family_config["centre_spread_x_m"],
                            ),
                            *x_bounds,
                        )
                    ),
                    float(
                        np.clip(
                            center_z
                            + rng.normal(
                                0,
                                family_config["centre_spread_z_m"],
                            ),
                            *z_bounds,
                        )
                    ),
                    _uniform(rng, family_config, "radius_x_m"),
                    _uniform(rng, family_config, "radius_z_m"),
                    _log_uniform(rng, family_config, "resistivity_ohm_m"),
                )
                for index in range(count)
            )
        model = insert_lenses(base, lenses, conflict_policy="last")
        parameters["lenses"] = [lens.to_dict() for lens in lenses]
        parameters["conflict_policy"] = "last"
        return model.resistivity_ohm_m, parameters
    raise ValueError(f"unsupported OOD family {family!r}")


def generate_benchmark_geology(
    grid: GeologyGrid,
    family: str,
    *,
    seed: int,
    configuration: Mapping[str, Any],
    regime: str = "id",
) -> BenchmarkGeology:
    """Generate one deterministic ID or structural-OOD geology.

    Parameters
    ----------
    grid : GeologyGrid
        Regular 2-D benchmark grid.
    family : str
        Family from :data:`ID_BENCHMARK_FAMILIES` or
        :data:`OOD_BENCHMARK_FAMILIES`.
    seed : int
        Explicit root seed.
    configuration : mapping
        Frozen mapping containing ``base_layered``, ``id_families``, and
        ``ood_families`` sections.
    regime : {"id", "ood"}, default="id"
        Distribution regime.

    Returns
    -------
    BenchmarkGeology
        Positive resistivity model and concrete sampled provenance.

    Examples
    --------
    >>> grid = GeologyGrid.regular_2d(nx=8, nz=6, dx_m=200, dz_m=100)
    >>> config = {  # doctest: +SKIP
    ...     "base_layered": {...},
    ...     "id_families": {...},
    ...     "ood_families": {...},
    ... }
    >>> model = generate_benchmark_geology(  # doctest: +SKIP
    ...     grid, "layered", seed=0, configuration=config
    ... )
    >>> model.resistivity_ohm_m.shape  # doctest: +SKIP
    (6, 8)
    """
    if not isinstance(grid, GeologyGrid) or grid.dimension != 2:
        raise TypeError("grid must be a 2-D GeologyGrid")
    regime = str(regime).strip().lower()
    family = str(family).strip()
    expected = (
        ID_BENCHMARK_FAMILIES
        if regime == "id"
        else OOD_BENCHMARK_FAMILIES
    )
    if regime not in {"id", "ood"} or family not in expected:
        raise ValueError("family and regime are incompatible")
    seed = _seed(seed)
    try:
        base_config = configuration["base_layered"]
        family_config = configuration[f"{regime}_families"][family]
    except KeyError as exc:
        raise ValueError(f"configuration is missing {exc.args[0]!r}") from exc
    rng = np.random.default_rng(seed)
    if regime == "id":
        values, parameters = _id_geology(
            grid, family, rng, base_config, family_config
        )
    else:
        values, parameters = _ood_geology(
            grid, family, rng, base_config, family_config
        )
    return BenchmarkGeology(
        family,
        regime,
        seed,
        grid,
        values,
        parameters,
    )
