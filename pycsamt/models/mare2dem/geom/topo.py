# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Topography interpolation utilities for MARE2DEM.

Port of ``m2d_parseTopo.m``.

MARE2DEM uses a right-handed coordinate system where the *y* axis runs
along the 2-D profile and *z* is positive downward.  Topography is
stored as a piecewise-linear ``(y, z)`` table.  This module provides
functions to interpolate bathymetry/topography depths and compute
seafloor/surface slope angles at arbitrary profile positions.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "parse_topo",
    "topo_depth",
    "topo_slope",
]


def parse_topo(
    topo: float | np.ndarray,
    y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate topography to profile positions *y*.

    Port of ``m2d_parseTopo.m``.

    Parameters
    ----------
    topo : float or array-like, shape (n_pts, 2)
        Topography input.  A single scalar gives a flat surface at
        constant depth.  An ``(n_pts, 2)`` array provides piecewise-
        linear topography with columns ``[y, z]``.
    y : array-like, shape (n,)
        Profile positions at which to interpolate.

    Returns
    -------
    z : numpy.ndarray, shape (n,)
        Depth of the topographic surface at each *y* position.
    slope_angle : numpy.ndarray, shape (n,)
        Surface slope angle in degrees, positive clockwise from the
        +y axis towards +z (down).  Zero for flat topography.
    on_node : numpy.ndarray of bool, shape (n,)
        ``True`` where a *y* position falls exactly on a topography
        node.  The slope is not well-defined at these points.

    Examples
    --------
    Flat topography at depth 1000 m (e.g. flat seafloor):

    >>> import numpy as np
    >>> from pycsamt.models.mare2dem.geom.topo import parse_topo
    >>> y = np.array([-5000.0, 0.0, 5000.0])
    >>> z, slope, on_node = parse_topo(1000.0, y)
    >>> z
    array([1000., 1000., 1000.])

    Variable topography:

    >>> topo_xy = np.array([[0.0, 500.0], [5000.0, 1000.0], [10000.0, 800.0]])
    >>> z, slope, on_node = parse_topo(topo_xy, np.array([2500.0, 5000.0]))
    """
    y = np.asarray(y, dtype=float).ravel()
    n = len(y)

    topo_arr = np.asarray(topo, dtype=float)

    # --- flat topography ---
    if topo_arr.ndim == 0 or (topo_arr.ndim == 1 and topo_arr.size == 1):
        depth = float(topo_arr.flat[0])
        return (
            np.full(n, depth),
            np.zeros(n),
            np.zeros(n, dtype=bool),
        )

    if topo_arr.ndim == 1 and topo_arr.size == 2:
        topo_arr = topo_arr.reshape(1, 2)

    # --- piecewise-linear topography ---
    # ensure unique y, sorted ascending
    _, idx = np.unique(topo_arr[:, 0], return_index=True)
    topo_arr = topo_arr[idx]

    ty = topo_arr[:, 0]
    tz = topo_arr[:, 1]

    # pad so topo extends well beyond all y positions
    y_min = min(y.min(), ty[0])
    y_max = max(y.max(), ty[-1])
    dy = y_max - y_min if y_max != y_min else 1e6

    ty_ext = np.concatenate(
        [
            [y_min - 2 * dy, y_min - dy],
            ty,
            [y_max + dy, y_max + 2 * dy],
        ]
    )
    tz_ext = np.concatenate(
        [
            [tz[0], tz[0]],
            tz,
            [tz[-1], tz[-1]],
        ]
    )

    # linear interpolation for depth
    z = np.interp(y, ty_ext, tz_ext)

    # piecewise-constant slope: atan2(dz, dy) at the left-node of each segment
    d_ty = np.diff(ty_ext)
    d_tz = np.diff(tz_ext)
    seg_angles = np.degrees(np.arctan2(d_tz, d_ty))
    seg_y = ty_ext[:-1]

    # "previous" interpolation: for each y, use slope of segment whose
    # left endpoint is just <= y
    indices = np.searchsorted(seg_y, y, side="right") - 1
    indices = np.clip(indices, 0, len(seg_angles) - 1)
    slope_angle = seg_angles[indices]

    on_node = np.isin(y, ty)

    return z, slope_angle, on_node


def topo_depth(
    topo: float | np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """Return interpolated topographic depth at positions *y*.

    Convenience wrapper around :func:`parse_topo`.
    """
    z, _, _ = parse_topo(topo, y)
    return z


def topo_slope(
    topo: float | np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """Return topographic slope angle (degrees) at positions *y*.

    Convenience wrapper around :func:`parse_topo`.
    """
    _, slope, _ = parse_topo(topo, y)
    return slope
