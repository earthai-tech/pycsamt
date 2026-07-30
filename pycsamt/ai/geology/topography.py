# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Physics-facing topographic surfaces on geological model grids.

Survey extraction remains in :mod:`pycsamt.topo`.  This module converts those
station samples, or explicit projected samples, into immutable elevation
rasters and derives air/earth masks for geological and Maxwell meshes.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..data.manifest import canonical_hash
from .fields import GeologyGrid

__all__ = [
    "TopographicSurface",
    "interpolate_topography",
    "topography_from_sites",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _method(value: str) -> str:
    result = str(value).strip().lower()
    if result not in {"linear", "nearest", "cubic"}:
        raise ValueError(
            "interpolation_method must be 'linear', 'nearest', or 'cubic'."
        )
    return result


def _align_sites(
    names: list[str],
    elevation: np.ndarray,
    chainage: np.ndarray,
    requested: Sequence[str] | None,
) -> tuple[tuple[str, ...], np.ndarray, np.ndarray, np.ndarray]:
    if requested is None:
        indices = np.arange(len(names), dtype=int)
        return tuple(names), elevation, chainage, indices
    lookup = {}
    for index, name in enumerate(names):
        key = str(name).strip().casefold()
        if key in lookup:
            raise ValueError(
                f"duplicate station name after normalization: {name!r}."
            )
        lookup[key] = index
    requested_names = tuple(str(name).strip() for name in requested)
    missing = [
        name for name in requested_names if name.casefold() not in lookup
    ]
    if missing:
        raise ValueError(
            f"topography is missing requested stations: {missing}."
        )
    indices = np.asarray(
        [lookup[name.casefold()] for name in requested_names], dtype=int
    )
    selected_chainage = chainage[indices]
    selected_chainage = selected_chainage - selected_chainage[0]
    return requested_names, elevation[indices], selected_chainage, indices


@dataclass(frozen=True)
class TopographicSurface:
    """Immutable terrain elevation raster aligned to a geological grid.

    Parameters
    ----------
    grid : GeologyGrid
        Grid whose horizontal cell centres are sampled by the surface.
    elevation_m : ndarray
        Elevation above ``vertical_datum``. Shape is ``(nx,)`` in 2-D or
        ``(ny, nx)`` in 3-D.
    reference_elevation_m : float
        Elevation corresponding to geological depth zero. Using the maximum
        terrain elevation keeps all surface depths non-negative.
    vertical_datum : str, default="metres above sea level"
        Human-readable vertical datum and unit description.
    source : str, default="array"
        Provenance label such as ``"sites"``, ``"dem"``, or ``"array"``.
    interpolation_method : {"linear", "nearest", "cubic"}, default="linear"
        Method used to rasterize samples.
    sample_coordinates_m : ndarray or None, optional
        Projected sample x positions in 2-D or x/y positions in 3-D.
    sample_elevation_m : ndarray or None, optional
        Elevations corresponding to sample coordinates.
    station_names : sequence of str, optional
        Station identifiers corresponding to samples.

    Examples
    --------
    >>> grid = GeologyGrid.regular_2d(nx=4, nz=3, dx_m=100, dz_m=50)
    >>> surface = TopographicSurface(grid, [100, 110, 105, 95], 110)
    >>> surface.relief_m
    15.0
    >>> surface.earth_mask().shape
    (3, 4)
    """

    grid: GeologyGrid
    elevation_m: np.ndarray
    reference_elevation_m: float
    vertical_datum: str = "metres above sea level"
    source: str = "array"
    interpolation_method: str = "linear"
    sample_coordinates_m: np.ndarray | None = None
    sample_elevation_m: np.ndarray | None = None
    station_names: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.grid, GeologyGrid):
            raise TypeError("grid must be a GeologyGrid.")
        expected = (
            (len(self.grid.x_m),)
            if self.grid.dimension == 2
            else (len(self.grid.y_m), len(self.grid.x_m))
        )
        elevation = np.asarray(self.elevation_m, dtype=float)
        if elevation.shape != expected or not np.all(np.isfinite(elevation)):
            raise ValueError(
                f"elevation_m must be finite and shaped {expected}."
            )
        reference = float(self.reference_elevation_m)
        if not np.isfinite(reference):
            raise ValueError("reference_elevation_m must be finite.")
        datum = str(self.vertical_datum).strip()
        source = str(self.source).strip()
        if not datum or not source:
            raise ValueError("vertical_datum and source cannot be empty.")
        method = _method(self.interpolation_method)
        coordinates = self.sample_coordinates_m
        samples = self.sample_elevation_m
        names = tuple(str(name).strip() for name in self.station_names)
        if (coordinates is None) != (samples is None):
            raise ValueError(
                "sample_coordinates_m and sample_elevation_m must be supplied together."
            )
        if coordinates is not None:
            coordinates = np.asarray(coordinates, dtype=float)
            samples = np.asarray(samples, dtype=float)
            width = 1 if self.grid.dimension == 2 else 2
            if coordinates.ndim == 1 and width == 1:
                coordinates = coordinates[:, None]
            if coordinates.ndim != 2 or coordinates.shape[1] != width:
                raise ValueError(
                    f"sample_coordinates_m must have shape (n, {width})."
                )
            if (
                samples.shape != (len(coordinates),)
                or not np.all(np.isfinite(coordinates))
                or not np.all(np.isfinite(samples))
            ):
                raise ValueError(
                    "sample coordinates/elevations must be finite and have matching lengths."
                )
            if names and (
                len(names) != len(samples)
                or len(set(names)) != len(names)
                or any(not name for name in names)
            ):
                raise ValueError(
                    "station_names must uniquely label every sample."
                )
        elif names:
            raise ValueError(
                "station_names require sample coordinates and elevations."
            )
        object.__setattr__(self, "elevation_m", _readonly(elevation))
        object.__setattr__(self, "reference_elevation_m", reference)
        object.__setattr__(self, "vertical_datum", datum)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "interpolation_method", method)
        object.__setattr__(
            self,
            "sample_coordinates_m",
            None if coordinates is None else _readonly(coordinates),
        )
        object.__setattr__(
            self,
            "sample_elevation_m",
            None if samples is None else _readonly(samples),
        )
        object.__setattr__(self, "station_names", names)

    @property
    def relief_m(self) -> float:
        """Return maximum minus minimum terrain elevation.

        Returns
        -------
        float
            Relief in metres.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> TopographicSurface(grid, [10, 13], 13).relief_m
        3.0
        """
        return float(np.ptp(self.elevation_m))

    @property
    def surface_depth_m(self) -> np.ndarray:
        """Return terrain depth below the configured reference elevation.

        Returns
        -------
        ndarray
            Horizontal surface shaped like :attr:`elevation_m`. Positive values
            lie below reference depth zero; negative values lie above it.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> TopographicSurface(grid, [10, 8], 10).surface_depth_m.tolist()
        [0.0, 2.0]
        """
        return _readonly(self.reference_elevation_m - self.elevation_m)

    @property
    def surface_hash(self) -> str:
        """Return a platform-stable terrain and provenance digest.

        Returns
        -------
        str
            SHA-256 digest.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> len(TopographicSurface(grid, [10, 9], 10).surface_hash)
        64
        """
        digest = hashlib.sha256()
        digest.update(
            np.ascontiguousarray(self.elevation_m, dtype="<f8").tobytes()
        )
        digest.update(canonical_hash(self.provenance()).encode("ascii"))
        return digest.hexdigest()

    def local_depth_m(self) -> np.ndarray:
        """Return each cell centre's signed depth below local terrain.

        Returns
        -------
        ndarray
            Array shaped like ``grid``. Non-negative values are in the earth;
            negative values are above terrain and belong to the air region.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> surface = TopographicSurface(grid, [10, 9], 10)
        >>> surface.local_depth_m().shape
        (2, 2)
        """
        horizontal_rank = self.elevation_m.ndim
        depth = self.grid.z_m.reshape(
            (len(self.grid.z_m),) + (1,) * horizontal_rank
        )
        return _readonly(depth - self.surface_depth_m)

    def earth_mask(self) -> np.ndarray:
        """Return cells on or below the local terrain surface.

        Returns
        -------
        ndarray of bool
            Read-only physics-facing earth mask shaped like ``grid``.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> TopographicSurface(grid, [10, 9], 10).earth_mask().dtype == bool
        True
        """
        return _readonly(self.local_depth_m() >= 0, bool)

    def air_mask(self) -> np.ndarray:
        """Return cells strictly above the local terrain surface.

        Returns
        -------
        ndarray of bool
            Logical complement of :meth:`earth_mask`.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> surface = TopographicSurface(grid, [10, 9], 10)
        >>> np.array_equal(surface.air_mask(), ~surface.earth_mask())
        True
        """
        return _readonly(~self.earth_mask(), bool)

    def slope_degrees(self) -> np.ndarray:
        """Return terrain slope magnitude at horizontal cell centres.

        Returns
        -------
        ndarray
            Slope angle in degrees, shaped like :attr:`elevation_m`.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=3, nz=2, dx_m=10, dz_m=1)
        >>> slope = TopographicSurface(grid, [0, 10, 20], 20).slope_degrees()
        >>> np.allclose(slope, 45)
        True
        """
        if self.grid.dimension == 2:
            gradient = np.gradient(
                self.elevation_m, self.grid.spacing_m[-1], edge_order=1
            )
            magnitude = np.abs(gradient)
        else:
            gradient_y, gradient_x = np.gradient(
                self.elevation_m,
                self.grid.spacing_m[1],
                self.grid.spacing_m[2],
                edge_order=1,
            )
            magnitude = np.sqrt(np.square(gradient_x) + np.square(gradient_y))
        return _readonly(np.degrees(np.arctan(magnitude)))

    def summary(self) -> dict[str, Any]:
        """Return compact JSON-compatible terrain diagnostics.

        Returns
        -------
        dict
            Elevation range, relief, slopes, air fraction, source, and datum.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> TopographicSurface(grid, [10, 9], 10).summary()["relief_m"]
        1.0
        """
        slopes = self.slope_degrees()
        return {
            "shape": list(self.elevation_m.shape),
            "elevation_range_m": [
                float(np.min(self.elevation_m)),
                float(np.max(self.elevation_m)),
            ],
            "reference_elevation_m": self.reference_elevation_m,
            "relief_m": self.relief_m,
            "mean_slope_deg": float(np.mean(slopes)),
            "maximum_slope_deg": float(np.max(slopes)),
            "air_cell_fraction": float(np.mean(self.air_mask())),
            "source": self.source,
            "vertical_datum": self.vertical_datum,
            "sample_count": 0
            if self.sample_elevation_m is None
            else len(self.sample_elevation_m),
        }

    def provenance(self) -> dict[str, Any]:
        """Return complete JSON-compatible surface provenance.

        Returns
        -------
        dict
            Grid, datum, source, interpolation, samples, and station names.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> TopographicSurface(grid, [10, 9], 10).provenance()["source"]
        'array'
        """
        return {
            "schema_version": 1,
            "grid": self.grid.to_dict(),
            "reference_elevation_m": self.reference_elevation_m,
            "vertical_datum": self.vertical_datum,
            "source": self.source,
            "interpolation_method": self.interpolation_method,
            "sample_coordinates_m": None
            if self.sample_coordinates_m is None
            else self.sample_coordinates_m.tolist(),
            "sample_elevation_m": None
            if self.sample_elevation_m is None
            else self.sample_elevation_m.tolist(),
            "station_names": list(self.station_names),
        }

    def to_npz(self, path: str | Path) -> Path:
        """Persist terrain and provenance in a pickle-free NPZ archive.

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
        >>> from tempfile import TemporaryDirectory
        >>> surface = TopographicSurface(
        ...     GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1), [10, 9], 10
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     restored = TopographicSurface.from_npz(
        ...         surface.to_npz(Path(directory) / "topo.npz")
        ...     )
        >>> restored.surface_hash == surface.surface_hash
        True
        """
        target = Path(path)
        np.savez_compressed(
            target,
            elevation_m=self.elevation_m,
            provenance_json=np.array(
                json.dumps(self.provenance(), sort_keys=True)
            ),
        )
        return target

    @classmethod
    def from_npz(cls, path: str | Path) -> TopographicSurface:
        """Load and validate a topographic archive without enabling pickle.

        Parameters
        ----------
        path : str or pathlib.Path
            Archive written by :meth:`to_npz`.

        Returns
        -------
        TopographicSurface
            Immutable restored surface.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> surface = TopographicSurface(
        ...     GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1), [10, 9], 10
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     restored = TopographicSurface.from_npz(
        ...         surface.to_npz(Path(directory) / "t.npz")
        ...     )
        >>> np.array_equal(restored.elevation_m, surface.elevation_m)
        True
        """
        with np.load(Path(path), allow_pickle=False) as archive:
            state = json.loads(str(archive["provenance_json"].item()))
            if state.get("schema_version") != 1:
                raise ValueError(
                    "unsupported TopographicSurface schema version."
                )
            return cls(
                GeologyGrid.from_dict(state["grid"]),
                archive["elevation_m"],
                state["reference_elevation_m"],
                state["vertical_datum"],
                state["source"],
                state["interpolation_method"],
                state.get("sample_coordinates_m"),
                state.get("sample_elevation_m"),
                tuple(state.get("station_names", [])),
            )


def interpolate_topography(
    grid: GeologyGrid,
    sample_coordinates_m: np.ndarray,
    sample_elevation_m: np.ndarray,
    *,
    interpolation_method: str = "linear",
    reference_elevation_m: float | None = None,
    vertical_datum: str = "metres above sea level",
    source: str = "array",
    station_names: Sequence[str] = (),
) -> TopographicSurface:
    """Interpolate projected elevation samples onto a geological grid.

    Parameters
    ----------
    grid : GeologyGrid
        Target 2-D or 3-D grid.
    sample_coordinates_m : ndarray
        Shape ``(n,)`` or ``(n, 1)`` x positions in 2-D; shape ``(n, 2)`` x/y
        projected coordinates in 3-D.
    sample_elevation_m : ndarray, shape (n,)
        Finite elevations in metres above the declared datum.
    interpolation_method : {"linear", "nearest", "cubic"}, default="linear"
        Interpolation method. Values outside the convex hull are filled from
        nearest samples rather than extrapolated polynomials.
    reference_elevation_m : float or None, optional
        Geological depth-zero elevation. Default is maximum raster elevation.
    vertical_datum, source : str, optional
        Provenance labels.
    station_names : sequence of str, optional
        Unique labels for every sample.

    Returns
    -------
    TopographicSurface
        Immutable raster and original samples.

    Examples
    --------
    >>> grid = GeologyGrid.regular_2d(nx=4, nz=3, dx_m=100, dz_m=50)
    >>> surface = interpolate_topography(grid, [0, 400], [100, 120])
    >>> surface.elevation_m.shape
    (4,)
    """
    if not isinstance(grid, GeologyGrid):
        raise TypeError("grid must be a GeologyGrid.")
    method = _method(interpolation_method)
    coordinates = np.asarray(sample_coordinates_m, dtype=float)
    elevation = np.asarray(sample_elevation_m, dtype=float)
    width = 1 if grid.dimension == 2 else 2
    if coordinates.ndim == 1 and width == 1:
        coordinates = coordinates[:, None]
    if coordinates.ndim != 2 or coordinates.shape[1] != width:
        raise ValueError(f"sample_coordinates_m must have shape (n, {width}).")
    if (
        elevation.shape != (len(coordinates),)
        or len(elevation) < 1
        or not np.all(np.isfinite(coordinates))
        or not np.all(np.isfinite(elevation))
    ):
        raise ValueError(
            "sample coordinates/elevations must be finite, non-empty, and have matching lengths."
        )
    if grid.dimension == 2:
        order = np.argsort(coordinates[:, 0])
        x = coordinates[order, 0]
        values = elevation[order]
        if np.any(np.diff(x) <= 0):
            raise ValueError("2-D sample x positions must be unique.")
        if method == "nearest":
            index = np.argmin(np.abs(grid.x_m[:, None] - x[None, :]), axis=1)
            raster = values[index]
        elif method == "cubic":
            if len(x) < 4:
                raise ValueError(
                    "2-D cubic topography needs at least four samples."
                )
            from scipy.interpolate import interp1d

            raster = interp1d(
                x,
                values,
                kind="cubic",
                bounds_error=False,
                fill_value=(values[0], values[-1]),
            )(grid.x_m)
        else:
            raster = np.interp(
                grid.x_m, x, values, left=values[0], right=values[-1]
            )
    else:
        if len(np.unique(coordinates, axis=0)) != len(coordinates):
            raise ValueError(
                "3-D projected sample coordinates must be unique."
            )
        if len(coordinates) < 3 and method != "nearest":
            raise ValueError(
                "3-D linear/cubic topography needs at least three samples."
            )
        from scipy.interpolate import griddata

        y, x = np.meshgrid(grid.y_m, grid.x_m, indexing="ij")
        query = np.column_stack([x.ravel(), y.ravel()])
        kind = "cubic" if method == "cubic" else method
        raster = griddata(coordinates, elevation, query, method=kind)
        if np.any(~np.isfinite(raster)):
            nearest = griddata(coordinates, elevation, query, method="nearest")
            raster = np.where(np.isfinite(raster), raster, nearest)
        raster = raster.reshape(len(grid.y_m), len(grid.x_m))
    reference = (
        float(np.max(raster))
        if reference_elevation_m is None
        else float(reference_elevation_m)
    )
    return TopographicSurface(
        grid,
        raster,
        reference,
        vertical_datum,
        source,
        method,
        coordinates,
        elevation,
        tuple(station_names),
    )


def topography_from_sites(
    sites: Any,
    grid: GeologyGrid,
    *,
    station_names: Sequence[str] | None = None,
    coordinates_m: np.ndarray | None = None,
    interpolation_method: str = "linear",
    reference_elevation_m: float | None = None,
    vertical_datum: str = "metres above sea level",
    profile_origin_m: float | None = None,
    allow_all_zero: bool = False,
) -> TopographicSurface:
    """Extract station elevations with :mod:`pycsamt.topo` and rasterize them.

    Parameters
    ----------
    sites : Sites or EDI-like collection
        Any container accepted by :func:`pycsamt.topo.extract_elevation`.
    grid : GeologyGrid
        Target geology grid.
    station_names : sequence of str or None, optional
        Requested station order/subset. Matching is case-insensitive.
    coordinates_m : ndarray or None, optional
        Required projected x/y coordinates for a 3-D grid. For 2-D, when
        supplied as ``(n, 2)``, cumulative projected chainage replaces the
        latitude/longitude-derived chainage from :mod:`pycsamt.topo`.
    interpolation_method, reference_elevation_m, vertical_datum
        Forwarded to :func:`interpolate_topography`.
    profile_origin_m : float or None, optional
        X coordinate assigned to zero chainage in 2-D. Default is the grid's
        minimum outer x edge.
    allow_all_zero : bool, default=False
        Permit an all-zero elevation collection. False rejects the ambiguity
        between genuine sea-level terrain and missing EDI elevations.

    Returns
    -------
    TopographicSurface
        Physics-facing raster with original station samples retained.

    Raises
    ------
    ValueError
        If station alignment, elevations, chainage, or required 3-D projected
        coordinates are invalid.

    Examples
    --------
    >>> from types import SimpleNamespace
    >>> stations = [
    ...     SimpleNamespace(
    ...         Head=SimpleNamespace(dataid="S00", elev=100, lat=5, lon=-3)
    ...     ),
    ...     SimpleNamespace(
    ...         Head=SimpleNamespace(dataid="S01", elev=110, lat=5, lon=-2.999)
    ...     ),
    ... ]
    >>> grid = GeologyGrid.regular_2d(nx=3, nz=2, dx_m=50, dz_m=20)
    >>> surface = topography_from_sites(stations, grid)
    >>> surface.source, surface.station_names
    ('sites', ('S00', 'S01'))
    """
    if not isinstance(grid, GeologyGrid):
        raise TypeError("grid must be a GeologyGrid.")
    from ...topo import (
        extract_chainage,
        extract_elevation,
        extract_station_names,
    )

    names = extract_station_names(sites)
    original_count = len(names)
    elevation = np.asarray(extract_elevation(sites), dtype=float)
    chainage = np.asarray(extract_chainage(sites), dtype=float) * 1000.0
    if not names or elevation.shape != (len(names),):
        raise ValueError("Sites did not provide aligned station elevations.")
    names, elevation, chainage, indices = _align_sites(
        names, elevation, chainage, station_names
    )
    if not np.all(np.isfinite(elevation)):
        raise ValueError("Sites elevations must be finite.")
    if not allow_all_zero and not np.any(elevation != 0.0):
        raise ValueError(
            "Sites elevations are all zero and may represent missing data."
        )
    if coordinates_m is None:
        if grid.dimension == 3:
            raise ValueError(
                "3-D Sites topography requires projected coordinates_m."
            )
        if len(chainage) > 1 and np.any(np.diff(chainage) <= 0):
            raise ValueError(
                "Sites chainage must increase in requested station order."
            )
        origin = (
            grid.extent_m["x"][0]
            if profile_origin_m is None
            else float(profile_origin_m)
        )
        coordinates = (origin + chainage)[:, None]
    else:
        coordinates = np.asarray(coordinates_m, dtype=float)
        if (
            coordinates.ndim != 2
            or coordinates.shape[0] != original_count
            or coordinates.shape[1] < 2
        ):
            raise ValueError(
                "coordinates_m must provide x/y for every original Sites station."
            )
        coordinates = coordinates[indices, :2]
        if grid.dimension == 2:
            segments = np.linalg.norm(np.diff(coordinates, axis=0), axis=1)
            if np.any(segments <= 0):
                raise ValueError(
                    "projected station coordinates must advance along the profile."
                )
            origin = (
                grid.extent_m["x"][0]
                if profile_origin_m is None
                else float(profile_origin_m)
            )
            coordinates = (
                origin + np.concatenate([[0.0], np.cumsum(segments)])
            )[:, None]
    return interpolate_topography(
        grid,
        coordinates,
        elevation,
        interpolation_method=interpolation_method,
        reference_elevation_m=reference_elevation_m,
        vertical_datum=vertical_datum,
        source="sites",
        station_names=names,
    )
