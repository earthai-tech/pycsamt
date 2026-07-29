# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regular geological grids and anisotropic correlated Gaussian fields.

The canonical array order is ``(z, x)`` in 2-D and ``(z, y, x)`` in 3-D.
Coordinates represent cell centres in metres, with depth ``z`` increasing
downward.  Field generation uses spectral synthesis and never invokes an EM
forward solver.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np

from ..data.manifest import canonical_hash

__all__ = [
    "GeologyGrid",
    "GaussianCorrelation",
    "CorrelatedField",
    "DirectionalVariogram",
    "generate_gaussian_field",
    "directional_variogram",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _axis(value: Any, name: str) -> np.ndarray:
    axis = np.asarray(value, dtype=float)
    if axis.ndim != 1 or axis.size < 2:
        raise ValueError(f"{name} must be a 1-D array with at least two cells.")
    if not np.all(np.isfinite(axis)) or not np.all(np.diff(axis) > 0):
        raise ValueError(f"{name} must be finite and strictly increasing.")
    return _readonly(axis)


def _regular_spacing(axis: np.ndarray, name: str) -> float:
    differences = np.diff(axis)
    spacing = float(np.mean(differences))
    if not np.allclose(
        differences, spacing, rtol=1e-8, atol=max(1e-12, abs(spacing) * 1e-10)
    ):
        raise ValueError(f"{name} must be regularly spaced for spectral synthesis.")
    return spacing


def _seed(value: int) -> int:
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError("seed must be an integer.")
    result = int(value)
    if result < 0 or result >= 2**64:
        raise ValueError("seed must be in [0, 2**64).")
    return result


@dataclass(frozen=True)
class GeologyGrid:
    """Immutable regular cell-centre grid for 2-D or 3-D geological priors.

    Parameters
    ----------
    x_m, z_m : array-like
        Strictly increasing horizontal and depth cell centres in metres.
    y_m : array-like or None, optional
        Strictly increasing second horizontal axis. Omit it for a 2-D profile.
    crs : str or None, optional
        Coordinate reference system identifier for horizontal coordinates.

    Examples
    --------
    Construct a 2-D grid directly:

    >>> grid = GeologyGrid(x_m=[50, 150, 250], z_m=[25, 75])
    >>> grid.shape
    (2, 3)
    >>> grid.dimension
    2

    Use :meth:`regular_3d` for a volume:

    >>> volume = GeologyGrid.regular_3d(
    ...     nx=4, ny=3, nz=2, dx_m=100, dy_m=200, dz_m=50
    ... )
    >>> volume.shape
    (2, 3, 4)
    """

    x_m: np.ndarray
    z_m: np.ndarray
    y_m: np.ndarray | None = None
    crs: str | None = None

    def __post_init__(self) -> None:
        x = _axis(self.x_m, "x_m")
        z = _axis(self.z_m, "z_m")
        y = None if self.y_m is None else _axis(self.y_m, "y_m")
        _regular_spacing(x, "x_m")
        _regular_spacing(z, "z_m")
        if y is not None:
            _regular_spacing(y, "y_m")
        crs = None if self.crs is None else str(self.crs).strip()
        if self.crs is not None and not crs:
            raise ValueError("crs cannot be empty.")
        object.__setattr__(self, "x_m", x)
        object.__setattr__(self, "z_m", z)
        object.__setattr__(self, "y_m", y)
        object.__setattr__(self, "crs", crs)

    @classmethod
    def regular_2d(
        cls,
        *,
        nx: int,
        nz: int,
        dx_m: float,
        dz_m: float,
        x_origin_m: float = 0.0,
        z_origin_m: float = 0.0,
        crs: str | None = None,
    ) -> GeologyGrid:
        """Construct a uniform 2-D cell-centre grid.

        Parameters
        ----------
        nx, nz : int
            At least two cells in the horizontal and depth directions.
        dx_m, dz_m : float
            Positive cell widths in metres.
        x_origin_m, z_origin_m : float, default=0.0
            Coordinate of the outer edge before the first cell centre.
        crs : str or None, optional
            Coordinate reference system identifier.

        Returns
        -------
        GeologyGrid
            Grid whose centres begin one half-cell from each origin.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=3, nz=2, dx_m=100, dz_m=50)
        >>> grid.x_m.tolist(), grid.z_m.tolist()
        ([50.0, 150.0, 250.0], [25.0, 75.0])
        """
        x = _centres(nx, dx_m, x_origin_m, "nx", "dx_m")
        z = _centres(nz, dz_m, z_origin_m, "nz", "dz_m")
        return cls(x, z, crs=crs)

    @classmethod
    def regular_3d(
        cls,
        *,
        nx: int,
        ny: int,
        nz: int,
        dx_m: float,
        dy_m: float,
        dz_m: float,
        x_origin_m: float = 0.0,
        y_origin_m: float = 0.0,
        z_origin_m: float = 0.0,
        crs: str | None = None,
    ) -> GeologyGrid:
        """Construct a uniform 3-D cell-centre grid.

        Parameters
        ----------
        nx, ny, nz : int
            At least two cells along each axis.
        dx_m, dy_m, dz_m : float
            Positive cell widths in metres.
        x_origin_m, y_origin_m, z_origin_m : float, default=0.0
            Outer-edge coordinate before the first centre on each axis.
        crs : str or None, optional
            Coordinate reference system identifier.

        Returns
        -------
        GeologyGrid
            Uniform grid with canonical shape ``(nz, ny, nx)``.

        Examples
        --------
        >>> grid = GeologyGrid.regular_3d(
        ...     nx=2, ny=3, nz=4, dx_m=10, dy_m=20, dz_m=5
        ... )
        >>> grid.shape
        (4, 3, 2)
        """
        x = _centres(nx, dx_m, x_origin_m, "nx", "dx_m")
        y = _centres(ny, dy_m, y_origin_m, "ny", "dy_m")
        z = _centres(nz, dz_m, z_origin_m, "nz", "dz_m")
        return cls(x, z, y_m=y, crs=crs)

    @property
    def dimension(self) -> int:
        """Return whether the grid is two- or three-dimensional.

        Returns
        -------
        {2, 3}
            Spatial dimension.

        Examples
        --------
        >>> GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1).dimension
        2
        """
        return 2 if self.y_m is None else 3

    @property
    def shape(self) -> tuple[int, ...]:
        """Return the canonical geological array shape.

        Returns
        -------
        tuple of int
            ``(nz, nx)`` in 2-D or ``(nz, ny, nx)`` in 3-D.

        Examples
        --------
        >>> GeologyGrid.regular_2d(nx=5, nz=3, dx_m=1, dz_m=1).shape
        (3, 5)
        """
        if self.y_m is None:
            return (len(self.z_m), len(self.x_m))
        return (len(self.z_m), len(self.y_m), len(self.x_m))

    @property
    def spacing_m(self) -> tuple[float, ...]:
        """Return regular spacing in canonical array-axis order.

        Returns
        -------
        tuple of float
            ``(dz, dx)`` in 2-D or ``(dz, dy, dx)`` in 3-D.

        Raises
        ------
        ValueError
            If any coordinate axis is not regularly spaced.

        Examples
        --------
        >>> GeologyGrid.regular_2d(nx=2, nz=2, dx_m=100, dz_m=25).spacing_m
        (25.0, 100.0)
        """
        dz = _regular_spacing(self.z_m, "z_m")
        dx = _regular_spacing(self.x_m, "x_m")
        if self.y_m is None:
            return (dz, dx)
        return (dz, _regular_spacing(self.y_m, "y_m"), dx)

    @property
    def extent_m(self) -> dict[str, tuple[float, float]]:
        """Return outer cell-edge extents along all available axes.

        Returns
        -------
        dict
            Axis names mapped to ``(minimum_edge, maximum_edge)`` in metres.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=100, dz_m=50)
        >>> grid.extent_m["x"]
        (0.0, 200.0)
        """
        axes = {"x": self.x_m, "z": self.z_m}
        if self.y_m is not None:
            axes["y"] = self.y_m
        result = {}
        for name, axis in axes.items():
            spacing = _regular_spacing(axis, f"{name}_m")
            result[name] = (
                float(axis[0] - spacing / 2),
                float(axis[-1] + spacing / 2),
            )
        return result

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable grid definition.

        Returns
        -------
        dict
            Schema version, coordinates, and CRS.

        Examples
        --------
        >>> GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1).to_dict()[
        ...     "dimension"
        ... ]
        2
        """
        return {
            "schema_version": 1,
            "dimension": self.dimension,
            "x_m": self.x_m.tolist(),
            "y_m": None if self.y_m is None else self.y_m.tolist(),
            "z_m": self.z_m.tolist(),
            "crs": self.crs,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> GeologyGrid:
        """Restore a validated serialized grid.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        GeologyGrid
            Immutable grid.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> GeologyGrid.from_dict(grid.to_dict()).shape == grid.shape
        True
        """
        if data.get("schema_version", 1) != 1:
            raise ValueError("unsupported GeologyGrid schema version.")
        return cls(data["x_m"], data["z_m"], y_m=data.get("y_m"), crs=data.get("crs"))


def _centres(
    count: int,
    spacing: float,
    origin: float,
    count_name: str,
    spacing_name: str,
) -> np.ndarray:
    if not isinstance(count, int) or isinstance(count, bool) or count < 2:
        raise ValueError(f"{count_name} must be an integer of at least two.")
    spacing = float(spacing)
    origin = float(origin)
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError(f"{spacing_name} must be finite and positive.")
    if not np.isfinite(origin):
        raise ValueError("grid origins must be finite.")
    return origin + (np.arange(count, dtype=float) + 0.5) * spacing


@dataclass(frozen=True)
class GaussianCorrelation:
    """Anisotropic Gaussian spatial-correlation model.

    Parameters
    ----------
    length_x_m, length_z_m : float
        Positive horizontal and vertical correlation lengths in metres for
        ``C(h) = exp(-0.5 * (h / length)**2)``.
    length_y_m : float or None, optional
        Second horizontal length required for 3-D generation.
    azimuth_deg : float, default=0.0
        Clockwise rotation of horizontal principal axes. It affects 3-D fields
        and is normalized to ``[0, 180)`` because Gaussian axes are bidirectional.

    Examples
    --------
    >>> model = GaussianCorrelation(1000, 100, length_y_m=500, azimuth_deg=210)
    >>> model.azimuth_deg
    30.0
    >>> model.anisotropy_xz
    10.0
    """

    length_x_m: float
    length_z_m: float
    length_y_m: float | None = None
    azimuth_deg: float = 0.0

    def __post_init__(self) -> None:
        for name in ("length_x_m", "length_z_m"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive.")
            object.__setattr__(self, name, value)
        if self.length_y_m is not None:
            y = float(self.length_y_m)
            if not np.isfinite(y) or y <= 0:
                raise ValueError("length_y_m must be finite and positive or None.")
            object.__setattr__(self, "length_y_m", y)
        azimuth = float(self.azimuth_deg)
        if not np.isfinite(azimuth):
            raise ValueError("azimuth_deg must be finite.")
        object.__setattr__(self, "azimuth_deg", azimuth % 180.0)

    @property
    def anisotropy_xz(self) -> float:
        """Return horizontal-to-vertical correlation-length ratio.

        Returns
        -------
        float
            ``length_x_m / length_z_m``.

        Examples
        --------
        >>> GaussianCorrelation(500, 100).anisotropy_xz
        5.0
        """
        return self.length_x_m / self.length_z_m

    def validate_grid(self, grid: GeologyGrid) -> None:
        """Validate dimensional compatibility with a geological grid.

        Parameters
        ----------
        grid : GeologyGrid
            Grid on which the model will be sampled.

        Returns
        -------
        None
            Successful return means all required correlation lengths exist.

        Raises
        ------
        TypeError
            If ``grid`` is not :class:`GeologyGrid`.
        ValueError
            If a 3-D grid has no ``length_y_m``.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=2, nz=2, dx_m=1, dz_m=1)
        >>> GaussianCorrelation(2, 1).validate_grid(grid) is None
        True
        """
        if not isinstance(grid, GeologyGrid):
            raise TypeError("grid must be a GeologyGrid.")
        if grid.dimension == 3 and self.length_y_m is None:
            raise ValueError("length_y_m is required for a 3-D grid.")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable correlation model.

        Returns
        -------
        dict
            Versioned Gaussian lengths and azimuth.

        Examples
        --------
        >>> GaussianCorrelation(10, 2).to_dict()["kind"]
        'gaussian'
        """
        return {
            "schema_version": 1,
            "kind": "gaussian",
            "length_x_m": self.length_x_m,
            "length_y_m": self.length_y_m,
            "length_z_m": self.length_z_m,
            "azimuth_deg": self.azimuth_deg,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> GaussianCorrelation:
        """Restore a validated Gaussian correlation model.

        Parameters
        ----------
        data : mapping
            State returned by :meth:`to_dict`.

        Returns
        -------
        GaussianCorrelation
            Immutable model.

        Examples
        --------
        >>> model = GaussianCorrelation(10, 2)
        >>> GaussianCorrelation.from_dict(model.to_dict()) == model
        True
        """
        if (
            data.get("schema_version", 1) != 1
            or data.get("kind", "gaussian") != "gaussian"
        ):
            raise ValueError("unsupported GaussianCorrelation state.")
        return cls(
            data["length_x_m"],
            data["length_z_m"],
            data.get("length_y_m"),
            data.get("azimuth_deg", 0.0),
        )


@dataclass(frozen=True)
class CorrelatedField:
    """Immutable generated scalar field with complete generation provenance.

    Parameters
    ----------
    values : ndarray
        Finite values shaped exactly like ``grid.shape``.
    grid : GeologyGrid
        Spatial cell-centre grid.
    correlation : GaussianCorrelation
        Requested correlation model.
    seed : int
        Recorded random seed.
    boundary : {"periodic", "reflect"}, default="reflect"
        Spectral boundary policy used during generation.
    standardized : bool, default=True
        Whether the stored sample was standardized to zero sample mean and unit
        sample standard deviation.

    Examples
    --------
    >>> grid = GeologyGrid.regular_2d(nx=8, nz=6, dx_m=100, dz_m=50)
    >>> field = generate_gaussian_field(
    ...     grid, GaussianCorrelation(300, 100), seed=4
    ... )
    >>> field.values.shape
    (6, 8)
    >>> field.seed
    4
    """

    values: np.ndarray
    grid: GeologyGrid
    correlation: GaussianCorrelation
    seed: int
    boundary: str = "reflect"
    standardized: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.grid, GeologyGrid):
            raise TypeError("grid must be a GeologyGrid.")
        if not isinstance(self.correlation, GaussianCorrelation):
            raise TypeError("correlation must be a GaussianCorrelation.")
        self.correlation.validate_grid(self.grid)
        values = np.asarray(self.values, dtype=float)
        if values.shape != self.grid.shape or not np.all(np.isfinite(values)):
            raise ValueError(f"values must be finite and shaped {self.grid.shape}.")
        if self.boundary not in {"periodic", "reflect"}:
            raise ValueError("boundary must be 'periodic' or 'reflect'.")
        object.__setattr__(self, "values", _readonly(values))
        object.__setattr__(self, "seed", _seed(self.seed))
        object.__setattr__(self, "standardized", bool(self.standardized))

    @property
    def field_hash(self) -> str:
        """Return a digest covering values and generation configuration.

        Returns
        -------
        str
            Lowercase SHA-256 digest. Values are hashed as contiguous little-
            endian float64 bytes for platform-stable identity.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=1, dz_m=1)
        >>> field = generate_gaussian_field(
        ...     grid, GaussianCorrelation(2, 1), seed=0
        ... )
        >>> len(field.field_hash)
        64
        """
        digest = hashlib.sha256()
        digest.update(np.ascontiguousarray(self.values, dtype="<f8").tobytes())
        digest.update(canonical_hash(self.provenance()).encode("ascii"))
        return digest.hexdigest()

    def provenance(self) -> dict[str, Any]:
        """Return JSON-compatible generation provenance without field values.

        Returns
        -------
        dict
            Grid, correlation, seed, boundary, and standardization state.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=1, dz_m=1)
        >>> field = generate_gaussian_field(
        ...     grid, GaussianCorrelation(2, 1), seed=0
        ... )
        >>> field.provenance()["seed"]
        0
        """
        return {
            "schema_version": 1,
            "grid": self.grid.to_dict(),
            "correlation": self.correlation.to_dict(),
            "seed": self.seed,
            "boundary": self.boundary,
            "standardized": self.standardized,
        }

    def rescale(
        self, *, mean: float = 0.0, standard_deviation: float = 1.0
    ) -> CorrelatedField:
        """Return a copy with a requested sample mean and standard deviation.

        Parameters
        ----------
        mean : float, default=0.0
            Target arithmetic mean.
        standard_deviation : float, default=1.0
            Positive target population standard deviation.

        Returns
        -------
        CorrelatedField
            New field with unchanged spatial pattern and ``standardized=False``
            unless the requested values are exactly zero and one.

        Examples
        --------
        >>> grid = GeologyGrid.regular_2d(nx=8, nz=8, dx_m=1, dz_m=1)
        >>> field = generate_gaussian_field(
        ...     grid, GaussianCorrelation(2, 1), seed=1
        ... ).rescale(mean=2, standard_deviation=3)
        >>> (
        ...     round(float(np.mean(field.values)), 12),
        ...     round(float(np.std(field.values)), 12),
        ... )
        (2.0, 3.0)
        """
        target_mean = float(mean)
        target_std = float(standard_deviation)
        if (
            not np.isfinite(target_mean)
            or not np.isfinite(target_std)
            or target_std <= 0
        ):
            raise ValueError(
                "mean must be finite and standard_deviation finite and positive."
            )
        current_std = float(np.std(self.values))
        if current_std <= np.finfo(float).eps:
            raise ValueError("cannot rescale a numerically constant field.")
        values = (self.values - np.mean(self.values)) / current_std
        values = target_mean + target_std * values
        return replace(
            self,
            values=values,
            standardized=target_mean == 0.0 and target_std == 1.0,
        )

    def to_npz(self, path: str | Path) -> Path:
        """Persist values and provenance in a pickle-free NPZ archive.

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
        >>> grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=1, dz_m=1)
        >>> field = generate_gaussian_field(
        ...     grid, GaussianCorrelation(2, 1), seed=0
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     path = field.to_npz(Path(directory) / "field.npz")
        ...     restored = CorrelatedField.from_npz(path)
        >>> restored.field_hash == field.field_hash
        True
        """
        target = Path(path)
        np.savez_compressed(
            target,
            values=self.values,
            provenance_json=np.array(json.dumps(self.provenance(), sort_keys=True)),
        )
        return target

    @classmethod
    def from_npz(cls, path: str | Path) -> CorrelatedField:
        """Load and validate a field archive without enabling pickle.

        Parameters
        ----------
        path : str or pathlib.Path
            Archive written by :meth:`to_npz`.

        Returns
        -------
        CorrelatedField
            Immutable restored field.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> field = generate_gaussian_field(
        ...     GeologyGrid.regular_2d(nx=3, nz=3, dx_m=1, dz_m=1),
        ...     GaussianCorrelation(2, 1),
        ...     seed=3,
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     path = field.to_npz(Path(directory) / "f.npz")
        ...     restored = CorrelatedField.from_npz(path)
        >>> np.array_equal(restored.values, field.values)
        True
        """
        with np.load(Path(path), allow_pickle=False) as archive:
            state = json.loads(str(archive["provenance_json"].item()))
            if state.get("schema_version") != 1:
                raise ValueError("unsupported CorrelatedField schema version.")
            return cls(
                archive["values"],
                GeologyGrid.from_dict(state["grid"]),
                GaussianCorrelation.from_dict(state["correlation"]),
                state["seed"],
                state["boundary"],
                state["standardized"],
            )


@dataclass(frozen=True)
class DirectionalVariogram:
    """Empirical semivariance along one canonical geological axis.

    Parameters
    ----------
    axis : {"x", "y", "z"}
        Direction along which cell pairs were compared.
    lag_m, semivariance : ndarray
        Positive lag distances and corresponding finite semivariances.
    pair_count : ndarray of int
        Number of finite cell pairs supporting each lag.

    Examples
    --------
    >>> result = DirectionalVariogram("x", [1, 2], [0.2, 0.5], [10, 8])
    >>> result.n_lags
    2
    """

    axis: str
    lag_m: np.ndarray
    semivariance: np.ndarray
    pair_count: np.ndarray

    def __post_init__(self) -> None:
        if self.axis not in {"x", "y", "z"}:
            raise ValueError("axis must be 'x', 'y', or 'z'.")
        lag = np.asarray(self.lag_m, dtype=float)
        semivariance = np.asarray(self.semivariance, dtype=float)
        count = np.asarray(self.pair_count)
        if (
            lag.ndim != 1
            or lag.size == 0
            or semivariance.shape != lag.shape
            or count.shape != lag.shape
        ):
            raise ValueError(
                "lag_m, semivariance, and pair_count must be matching non-empty 1-D arrays."
            )
        if (
            not np.all(np.isfinite(lag))
            or np.any(lag <= 0)
            or not np.all(np.diff(lag) > 0)
        ):
            raise ValueError("lag_m must be finite, positive, and strictly increasing.")
        if not np.all(np.isfinite(semivariance)) or np.any(semivariance < 0):
            raise ValueError("semivariance must be finite and non-negative.")
        if not np.issubdtype(count.dtype, np.integer) or np.any(count <= 0):
            raise ValueError("pair_count must contain positive integers.")
        object.__setattr__(self, "lag_m", _readonly(lag))
        object.__setattr__(self, "semivariance", _readonly(semivariance))
        object.__setattr__(self, "pair_count", _readonly(count, np.int64))

    @property
    def n_lags(self) -> int:
        """Return the number of empirical lag bins.

        Returns
        -------
        int
            Length of the variogram arrays.

        Examples
        --------
        >>> DirectionalVariogram("z", [1], [0.1], [4]).n_lags
        1
        """
        return len(self.lag_m)


def generate_gaussian_field(
    grid: GeologyGrid,
    correlation: GaussianCorrelation,
    *,
    seed: int,
    boundary: str = "reflect",
    standardize: bool = True,
) -> CorrelatedField:
    """Generate a deterministic anisotropic Gaussian random field.

    Parameters
    ----------
    grid : GeologyGrid
        Regular 2-D or 3-D cell-centre grid.
    correlation : GaussianCorrelation
        Requested Gaussian covariance lengths and horizontal azimuth.
    seed : int
        Explicit seed recorded in the returned field.
    boundary : {"reflect", "periodic"}, default="reflect"
        ``periodic`` synthesizes directly on the requested grid. ``reflect``
        synthesizes on a grid doubled along each axis and crops the centre,
        reducing wrap-around correlation at opposite model edges.
    standardize : bool, default=True
        Shift and scale the realized sample to zero mean and unit population
        standard deviation.

    Returns
    -------
    CorrelatedField
        Immutable field with complete generation provenance.

    Raises
    ------
    ValueError
        If the grid is irregular, correlation model is incompatible, or the
        filtered realization is numerically constant.

    Examples
    --------
    >>> grid = GeologyGrid.regular_2d(nx=32, nz=16, dx_m=100, dz_m=50)
    >>> model = GaussianCorrelation(length_x_m=500, length_z_m=100)
    >>> first = generate_gaussian_field(grid, model, seed=12)
    >>> second = generate_gaussian_field(grid, model, seed=12)
    >>> np.array_equal(first.values, second.values)
    True
    >>> (
    ...     abs(round(float(np.mean(first.values)), 12)),
    ...     round(float(np.std(first.values)), 12),
    ... )
    (0.0, 1.0)
    """
    if not isinstance(grid, GeologyGrid):
        raise TypeError("grid must be a GeologyGrid.")
    if not isinstance(correlation, GaussianCorrelation):
        raise TypeError("correlation must be a GaussianCorrelation.")
    correlation.validate_grid(grid)
    seed = _seed(seed)
    if boundary not in {"reflect", "periodic"}:
        raise ValueError("boundary must be 'reflect' or 'periodic'.")
    spacing = grid.spacing_m
    rng = np.random.default_rng(seed)
    if boundary == "periodic":
        white = rng.standard_normal(grid.shape)
        crop_starts = None
    else:
        core = rng.standard_normal(grid.shape)
        pad_width = tuple((size // 2, size - size // 2) for size in grid.shape)
        white = np.pad(core, pad_width, mode="reflect")
        crop_starts = tuple(width[0] for width in pad_width)
    shape = white.shape
    spectral = np.fft.fftn(white)
    axes = [
        2 * np.pi * np.fft.fftfreq(size, d=step) for size, step in zip(shape, spacing)
    ]
    if grid.dimension == 2:
        kz, kx = np.meshgrid(axes[0], axes[1], indexing="ij")
        exponent = np.square(kz * correlation.length_z_m) + np.square(
            kx * correlation.length_x_m
        )
    else:
        kz, ky, kx = np.meshgrid(axes[0], axes[1], axes[2], indexing="ij")
        angle = np.deg2rad(correlation.azimuth_deg)
        k_major = kx * np.cos(angle) + ky * np.sin(angle)
        k_minor = -kx * np.sin(angle) + ky * np.cos(angle)
        exponent = (
            np.square(kz * correlation.length_z_m)
            + np.square(k_major * correlation.length_x_m)
            + np.square(k_minor * correlation.length_y_m)
        )
    amplitude = np.exp(-0.25 * exponent)
    values = np.fft.ifftn(spectral * amplitude).real
    if boundary == "reflect":
        values = values[
            tuple(
                slice(start, start + target)
                for start, target in zip(crop_starts, grid.shape)
            )
        ]
    sample_std = float(np.std(values))
    if not np.isfinite(sample_std) or sample_std <= np.finfo(float).eps:
        raise ValueError(
            "correlation lengths and grid produced a numerically constant field."
        )
    if standardize:
        values = (values - np.mean(values)) / sample_std
    return CorrelatedField(values, grid, correlation, seed, boundary, standardize)


def directional_variogram(
    field: CorrelatedField,
    axis: str,
    *,
    max_lag_cells: int | None = None,
) -> DirectionalVariogram:
    """Calculate an unbinned empirical directional semivariogram.

    Parameters
    ----------
    field : CorrelatedField
        Finite scalar field.
    axis : {"x", "y", "z"}
        Geological direction. ``"y"`` is unavailable for 2-D fields.
    max_lag_cells : int or None, optional
        Maximum positive integer cell offset. The default is half the selected
        axis length, with at least one lag.

    Returns
    -------
    DirectionalVariogram
        Lag distance, mean half-squared difference, and pair count.

    Examples
    --------
    >>> grid = GeologyGrid.regular_2d(nx=16, nz=8, dx_m=100, dz_m=50)
    >>> field = generate_gaussian_field(
    ...     grid, GaussianCorrelation(300, 100), seed=1
    ... )
    >>> variogram = directional_variogram(field, "x", max_lag_cells=3)
    >>> variogram.n_lags
    3
    >>> variogram.pair_count[0] > variogram.pair_count[-1]
    True
    """
    if not isinstance(field, CorrelatedField):
        raise TypeError("field must be a CorrelatedField.")
    axis_positions = {"z": 0, "x": field.values.ndim - 1}
    if field.grid.dimension == 3:
        axis_positions["y"] = 1
    if axis not in axis_positions:
        raise ValueError(f"axis must be one of {tuple(axis_positions)} for this field.")
    position = axis_positions[axis]
    size = field.values.shape[position]
    maximum = max(1, size // 2) if max_lag_cells is None else max_lag_cells
    if (
        not isinstance(maximum, int)
        or isinstance(maximum, bool)
        or maximum < 1
        or maximum >= size
    ):
        raise ValueError("max_lag_cells must be an integer in [1, axis_size).")
    spacing_index = {"z": 0, "y": 1, "x": field.values.ndim - 1}[axis]
    spacing = field.grid.spacing_m[spacing_index]
    semivariance = []
    pairs = []
    for lag in range(1, maximum + 1):
        left = [slice(None)] * field.values.ndim
        right = [slice(None)] * field.values.ndim
        left[position] = slice(None, -lag)
        right[position] = slice(lag, None)
        difference = field.values[tuple(right)] - field.values[tuple(left)]
        semivariance.append(0.5 * float(np.mean(np.square(difference))))
        pairs.append(difference.size)
    return DirectionalVariogram(
        axis,
        spacing * np.arange(1, maximum + 1, dtype=float),
        np.asarray(semivariance),
        np.asarray(pairs, dtype=np.int64),
    )
