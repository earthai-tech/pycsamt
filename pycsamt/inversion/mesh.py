# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Mesh descriptions and builders for EM inversion workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject

__all__ = [
    "InversionMesh",
    "build_1d_tensor_mesh",
    "build_3d_tensor_mesh",
    "build_fd2d_grid",
    "core_rho_from_start",
    "depth_widths",
]


@dataclass
class InversionMesh(PyCSAMTObject, MetadataMixin):
    """Lightweight mesh/grid descriptor.

    ``InversionMesh`` is the common mesh metadata object carried by
    :class:`pycsamt.inversion.results.InversionResult`. Numerical engines may
    keep their real discretization object in ``native`` while exposing common
    profile/depth centers through ``x_centers`` and ``z_centers``.

    Parameters
    ----------
    dimension : {"1d", "2d", "3d"}, default "1d"
        Mesh dimensionality represented by this descriptor.
    x_centers : array-like of float, optional
        Horizontal/profile cell centers in metres. For 1-D meshes this is
        usually ``[0.0]``.
    z_centers : array-like of float, optional
        Depth cell centers in metres, positive downward.
    native : object, optional
        Backend-native mesh/grid object, for example a SimPEG TensorMesh or
        built-in forward ``Grid2D``.
    metadata : dict, optional
        Free-form mesh provenance such as builder name, backend engine, padding,
        or core-cell shape.

    Examples
    --------
    >>> from pycsamt.inversion.mesh import InversionMesh
    >>> mesh = InversionMesh.for_1d([25.0, 125.0, 275.0])
    >>> mesh.dimension
    '1d'
    >>> mesh.z_centers.tolist()
    [25.0, 125.0, 275.0]

    References
    ----------
    .. [1] Oldenburg, D. W. and Li, Y. (2005). Inversion for applied
       geophysics: A tutorial. In *Near-Surface Geophysics*, SEG.
    """

    dimension: str = "1d"
    x_centers: Any = None
    z_centers: Any = None
    native: Any = field(default=None, repr=False)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.dimension = str(self.dimension).lower()
        self.x_centers = _array_or_none(self.x_centers)
        self.z_centers = _array_or_none(self.z_centers)
        self.validate()

    @classmethod
    def for_1d(cls, depths: Any) -> InversionMesh:
        """Build a 1-D mesh descriptor from depth centres.

        Parameters
        ----------
        depths : array-like of float
            Positive-downward depth cell centers in metres.

        Returns
        -------
        InversionMesh
            Descriptor with ``dimension="1d"``, ``x_centers=[0.0]``, and the
            supplied depth centers.

        Examples
        --------
        >>> from pycsamt.inversion.mesh import InversionMesh
        >>> InversionMesh.for_1d([10.0, 30.0]).z_centers.tolist()
        [10.0, 30.0]
        """
        return cls(dimension="1d", x_centers=np.array([0.0]), z_centers=depths)

    def validate(self) -> None:
        if self.dimension not in {"1d", "2d", "3d"}:
            raise ValueError("dimension must be '1d', '2d', or '3d'.")
        if self.z_centers is not None and np.any(self.z_centers < 0):
            raise ValueError("z_centers must be positive downward.")


def _array_or_none(value: Any) -> np.ndarray | None:
    if value is None:
        return None
    return np.asarray(value, dtype=float)


def depth_widths(
    depth_max: float,
    n_cells: int,
    options: dict[str, Any] | None = None,
) -> np.ndarray:
    """Return geometrically growing positive depth cell widths.

    Parameters
    ----------
    depth_max : float
        Total target depth in metres. The returned widths sum to this value.
    n_cells : int
        Number of depth cells.
    options : dict, optional
        Mesh options. Recognized keys are ``min_cell_size`` and
        ``growth_factor``. Defaults are chosen conservatively from
        ``depth_max`` and ``n_cells``.

    Returns
    -------
    ndarray
        Positive cell widths in metres.

    Raises
    ------
    ValueError
        If ``depth_max``, ``n_cells``, ``min_cell_size``, or ``growth_factor``
        are not positive.

    Examples
    --------
    >>> from pycsamt.inversion.mesh import depth_widths
    >>> depth_widths(100.0, 4, {"growth_factor": 1.0}).round(2).tolist()
    [25.0, 25.0, 25.0, 25.0]

    References
    ----------
    .. [1] Ward, S. H. and Hohmann, G. W. (1988). Electromagnetic theory for
       geophysical applications. In *Electromagnetic Methods in Applied
       Geophysics*, volume 1, SEG.
    """
    options = dict(options or {})
    depth_max = float(depth_max)
    n_cells = int(n_cells)
    if depth_max <= 0.0:
        raise ValueError("depth_max must be positive.")
    if n_cells <= 0:
        raise ValueError("n_cells must be positive.")
    min_cell = float(options.get("min_cell_size", max(depth_max / n_cells / 4.0, 1.0)))
    growth = float(options.get("growth_factor", 1.08))
    if min_cell <= 0.0:
        raise ValueError("min_cell_size must be positive.")
    if growth <= 0.0:
        raise ValueError("growth_factor must be positive.")
    widths = min_cell * growth ** np.arange(n_cells, dtype=float)
    widths *= depth_max / np.sum(widths)
    return widths


def build_1d_tensor_mesh(
    start: Any,
    options: dict[str, Any] | None,
    tensor_mesh_cls: Any,
) -> tuple[Any, np.ndarray]:
    """Build a 1-D TensorMesh-like object and positive-downward centres.

    This helper centralizes the depth discretization used by optional engines
    such as SimPEG. The caller supplies the actual TensorMesh class so this
    module stays free of optional dependencies.

    Parameters
    ----------
    start : object
        Starting model exposing ``n_layers``, ``resistivities``, and
        ``thicknesses``.
    options : dict, optional
        Mesh controls. Recognized keys include ``n_cells``, ``depth_max``,
        ``origin``, ``min_cell_size``, and ``growth_factor``.
    tensor_mesh_cls : type
        TensorMesh-like constructor accepting ``([widths], origin=...)``.

    Returns
    -------
    mesh : object
        Backend-native 1-D tensor mesh.
    z_centers : ndarray
        Positive-downward depth cell centers in metres.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.inversion.mesh import build_1d_tensor_mesh
    >>> from pycsamt.inversion.model import StartingModel
    >>> class TensorMesh:
    ...     def __init__(self, widths, origin="0"):
    ...         self.widths = widths
    ...         self.origin = origin
    >>> start = StartingModel([100.0, 300.0], [500.0])
    >>> mesh, z = build_1d_tensor_mesh(
    ...     start,
    ...     {"n_cells": 2, "depth_max": 100.0, "growth_factor": 1.0},
    ...     TensorMesh,
    ... )
    >>> np.round(z, 1).tolist()
    [25.0, 75.0]
    """
    options = dict(options or {})
    n_layers = int(getattr(start, "n_layers", len(start.resistivities)))
    thicknesses = np.asarray(start.thicknesses, dtype=float)
    n_cells = int(options.get("n_cells", max(32, n_layers * 12)))
    depth_max = float(
        options.get(
            "depth_max",
            max(float(np.sum(thicknesses)) * 3.0, float(thicknesses[-1]) * 4.0),
        )
    )
    widths = depth_widths(depth_max, n_cells, options)
    mesh = tensor_mesh_cls([widths], origin=options.get("origin", "0"))
    z_centers = np.cumsum(widths) - 0.5 * widths
    return mesh, z_centers


def build_3d_tensor_mesh(
    station_x: Any,
    station_y: Any,
    options: dict[str, Any] | None,
    tensor_mesh_cls: Any,
) -> tuple[Any, dict[str, np.ndarray]]:
    """Build a 3-D TensorMesh-like object around station coordinates.

    The helper constructs uniform horizontal cells around the station footprint
    and geometrically growing depth cells. It is used by optional 3-D physics
    paths while keeping the native mesh class injectable.

    Parameters
    ----------
    station_x, station_y : array-like of float
        Station coordinates in metres.
    options : dict, optional
        Mesh controls. Recognized keys include ``x_pad``, ``y_pad``,
        ``depth_max``, ``nx``, ``ny``, ``nz``, ``min_cell_size``, and
        ``growth_factor``.
    tensor_mesh_cls : type
        TensorMesh-like constructor accepting ``([hx, hy, hz], origin=...)``.

    Returns
    -------
    mesh : object
        Backend-native 3-D tensor mesh.
    centers : dict of ndarray
        Cell-center arrays with keys ``"x"``, ``"y"``, ``"z"``, and
        ``"z_depth"``. ``z`` follows the native mesh coordinate with depth
        negative; ``z_depth`` is positive downward.

    Examples
    --------
    >>> from pycsamt.inversion.mesh import build_3d_tensor_mesh
    >>> class TensorMesh:
    ...     def __init__(self, widths, origin=None):
    ...         self.widths = widths
    ...         self.origin = origin
    >>> mesh, centers = build_3d_tensor_mesh(
    ...     [0.0, 100.0],
    ...     [0.0, 50.0],
    ...     {"nx": 2, "ny": 2, "nz": 2, "depth_max": 100.0},
    ...     TensorMesh,
    ... )
    >>> sorted(centers)
    ['x', 'y', 'z', 'z_depth']
    """
    options = dict(options or {})
    station_x = np.asarray(station_x, dtype=float)
    station_y = np.asarray(station_y, dtype=float)
    if station_x.size == 0 or station_y.size == 0:
        raise ValueError("station coordinates are required for a 3-D mesh.")
    x_pad = float(options.get("x_pad", 1000.0))
    y_pad = float(options.get("y_pad", 1000.0))
    depth_max = float(options.get("depth_max", 3000.0))
    nx = int(options.get("nx", max(8, min(32, station_x.size * 4))))
    ny = int(options.get("ny", max(8, min(32, station_y.size * 4))))
    nz = int(options.get("nz", 16))
    if nx <= 0 or ny <= 0 or nz <= 0:
        raise ValueError("nx, ny, and nz must be positive.")
    x_min = float(np.min(station_x) - x_pad)
    x_max = float(np.max(station_x) + x_pad)
    y_min = float(np.min(station_y) - y_pad)
    y_max = float(np.max(station_y) + y_pad)
    hx = np.full(nx, (x_max - x_min) / nx, dtype=float)
    hy = np.full(ny, (y_max - y_min) / ny, dtype=float)
    hz = depth_widths(depth_max, nz, options)
    origin = np.array([x_min, y_min, -float(np.sum(hz))], dtype=float)
    mesh = tensor_mesh_cls([hx, hy, hz], origin=origin)
    centers = {
        "x": x_min + np.cumsum(hx) - 0.5 * hx,
        "y": y_min + np.cumsum(hy) - 0.5 * hy,
        "z": origin[2] + np.cumsum(hz) - 0.5 * hz,
    }
    centers["z_depth"] = -centers["z"]
    return mesh, centers


def build_fd2d_grid(
    start: Any,
    station_x: Any,
    options: dict[str, Any] | None,
    grid_cls: Any,
    *,
    make_padding_func: Any | None = None,
) -> tuple[Any, tuple[int, int]]:
    """Build the finite-difference 2-D grid used by built-in inversion.

    This builder creates the starting ``Grid2D``-like object for the built-in
    finite-difference MT/AMT/CSAMT profile inversion. It shifts station
    coordinates into the grid coordinate system, adds optional lateral/depth
    padding, and expands the layered starting model into a 2-D resistivity
    array.

    Parameters
    ----------
    start : object
        Starting model exposing ``resistivities`` and ``thicknesses``.
    station_x : array-like of float
        Station positions along profile in metres. At least two distinct
        stations are required.
    options : dict, optional
        Grid controls. Recognized keys include ``nx``/``fd2d_nx``,
        ``n_pad``/``fd2d_n_pad``, ``pad_factor``, ``x_margin``, ``x_max``, and
        ``halfspace_thickness``.
    grid_cls : type
        Grid2D-like constructor accepting ``dx``, ``dz``, ``resistivity``,
        ``x_stations``, ``n_pad``, and ``name``.
    make_padding_func : callable, optional
        Padding-width function. If omitted, :func:`pycsamt.forward.make_padding`
        is imported lazily.

    Returns
    -------
    grid : object
        Built finite-difference grid.
    core_shape : tuple of int
        ``(nz_core, nx_core)`` shape of the unpadded inversion core.

    Examples
    --------
    >>> from pycsamt.inversion.mesh import build_fd2d_grid
    >>> from pycsamt.inversion.model import StartingModel
    >>> class Grid:
    ...     def __init__(self, dx, dz, resistivity, x_stations, n_pad, name):
    ...         self.dx = dx
    ...         self.dz = dz
    ...         self.resistivity = resistivity
    ...         self.x_stations = x_stations
    >>> start = StartingModel([100.0, 300.0], [500.0])
    >>> grid, core_shape = build_fd2d_grid(
    ...     start, [0.0, 1000.0], {"nx": 2, "n_pad": 0}, Grid
    ... )
    >>> core_shape
    (2, 2)
    """
    options = dict(options or {})
    station_x = np.asarray(station_x, dtype=float)
    n_st = int(station_x.size)
    if n_st < 2:
        raise ValueError(
            "finite-difference 2-D inversion requires at least two stations."
        )
    x0 = float(np.nanmin(station_x))
    shifted = station_x - x0
    spread = float(np.nanmax(shifted))
    if spread <= 0.0:
        raise ValueError("station_x must span a non-zero profile length.")
    nx_core = int(options.get("nx", options.get("fd2d_nx", max(n_st - 1, 2))))
    n_pad = int(options.get("n_pad", options.get("fd2d_n_pad", 0)))
    pad_factor = float(options.get("pad_factor", options.get("fd2d_pad_factor", 1.3)))
    x_margin = float(
        options.get("x_margin", options.get("fd2d_x_margin", 0.05 * spread))
    )
    x_max = float(
        options.get("x_max", options.get("fd2d_x_max", spread + 2.0 * x_margin))
    )
    if nx_core <= 0 or n_pad < 0:
        raise ValueError("nx must be positive and n_pad must be non-negative.")
    dx_core = np.full(nx_core, x_max / nx_core, dtype=float)
    if n_pad > 0:
        make_padding_func = _resolve_make_padding(make_padding_func)
        dx_pad = make_padding_func(dx_core[0], n_pad, pad_factor)
    else:
        dx_pad = np.array([], dtype=float)
    dx = np.r_[dx_pad[::-1], dx_core, dx_pad]
    x_offset = float(dx_pad.sum() + x_margin)
    x_stations = shifted + x_offset

    thicknesses = np.asarray(start.thicknesses, dtype=float)
    dz_core = thicknesses.copy()
    last = dz_core[-1] if dz_core.size else float(options.get("dz_min", 100.0))
    dz_core = np.r_[dz_core, float(options.get("halfspace_thickness", 3.0 * last))]
    if n_pad > 0:
        make_padding_func = _resolve_make_padding(make_padding_func)
        dz_pad = make_padding_func(dz_core[-1], n_pad, pad_factor)
    else:
        dz_pad = np.array([], dtype=float)
    dz = np.r_[dz_core, dz_pad]
    nz_core = int(dz_core.size)
    rho_core = core_rho_from_start(start, nz_core, nx_core)
    rho = _pad_core_rho(rho_core, n_pad)
    grid = grid_cls(
        dx=dx,
        dz=dz,
        resistivity=rho,
        x_stations=x_stations,
        n_pad=n_pad,
        name="builtin_fd2d_start",
    )
    grid._pycsamt_x_offset = x_offset
    return grid, (nz_core, nx_core)


def core_rho_from_start(start: Any, nz: int, nx: int) -> np.ndarray:
    """Expand a layered starting model into a 2-D core resistivity grid.

    Parameters
    ----------
    start : object
        Starting model exposing ``resistivities``.
    nz, nx : int
        Number of vertical and horizontal core cells.

    Returns
    -------
    ndarray
        Resistivity grid with shape ``(nz, nx)``. If ``nz`` is larger than the
        number of supplied layers, the final layer resistivity is repeated.

    Examples
    --------
    >>> from pycsamt.inversion.mesh import core_rho_from_start
    >>> from pycsamt.inversion.model import StartingModel
    >>> start = StartingModel([100.0, 300.0], [500.0])
    >>> core_rho_from_start(start, 3, 2).tolist()
    [[100.0, 100.0], [300.0, 300.0], [300.0, 300.0]]
    """
    rho_col = np.asarray(start.resistivities, dtype=float)
    if rho_col.size < nz:
        rho_col = np.r_[rho_col, np.repeat(rho_col[-1], nz - rho_col.size)]
    rho_col = rho_col[:nz]
    return np.tile(rho_col[:, None], (1, nx))


def _pad_core_rho(rho_core: np.ndarray, n_pad: int) -> np.ndarray:
    if n_pad <= 0:
        return rho_core
    left = np.repeat(rho_core[:, :1], n_pad, axis=1)
    right = np.repeat(rho_core[:, -1:], n_pad, axis=1)
    rho_rows = np.c_[left, rho_core, right]
    bottom = np.repeat(rho_rows[-1:, :], n_pad, axis=0)
    return np.r_[rho_rows, bottom]


def _resolve_make_padding(make_padding_func: Any | None) -> Any:
    if make_padding_func is not None:
        return make_padding_func
    from ..forward import make_padding

    return make_padding
