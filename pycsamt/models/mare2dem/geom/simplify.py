# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Douglas-Peucker polyline simplification.

Port of ``m2d_dpsimplify.m``.

Used to reduce the number of points in topography or geometry polygons
before passing them to the Triangle mesh generator.
"""

from __future__ import annotations

import numpy as np

__all__ = ["dp_simplify"]


def _perpendicular_distance(
    points: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
) -> np.ndarray:
    """Return the perpendicular distance from each point to the line p1→p2."""
    d = p2 - p1
    norm = np.linalg.norm(d)
    if norm < 1e-12:
        return np.linalg.norm(points - p1, axis=1)
    # signed area / base length = perpendicular distance
    cross = np.abs(
        (d[0]) * (p1[1] - points[:, 1]) - (p1[0] - points[:, 0]) * (d[1])
    )
    return cross / norm


def _dp_recursive(
    points: np.ndarray,
    indices: list[int],
    start: int,
    end: int,
    tolerance: float,
    keep: set[int],
) -> None:
    """Recursively mark points to keep."""
    if end <= start + 1:
        return
    seg = points[start : end + 1]
    dists = _perpendicular_distance(seg[1:-1], points[start], points[end])
    if len(dists) == 0:
        return
    max_dist = dists.max()
    max_idx = dists.argmax() + start + 1  # offset by start + 1

    if max_dist > tolerance:
        keep.add(max_idx)
        _dp_recursive(points, indices, start, max_idx, tolerance, keep)
        _dp_recursive(points, indices, max_idx, end, tolerance, keep)


def dp_simplify(
    points: np.ndarray,
    tolerance: float,
) -> np.ndarray:
    """Simplify a polyline using the Douglas-Peucker algorithm.

    Port of ``m2d_dpsimplify.m``.

    Parameters
    ----------
    points : array-like, shape (n, 2)
        Input polyline vertices.
    tolerance : float
        Maximum allowable perpendicular deviation. Points whose
        distance from the simplified segment exceeds this value
        are retained.

    Returns
    -------
    numpy.ndarray, shape (m, 2)
        Simplified polyline with ``m <= n`` vertices.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.models.mare2dem.geom.simplify import dp_simplify
    >>> pts = np.column_stack(
    ...     [np.linspace(0, 10, 100), np.sin(np.linspace(0, np.pi, 100))]
    ... )
    >>> simplified = dp_simplify(pts, tolerance=0.05)
    >>> len(simplified) < 100
    True
    """
    pts = np.asarray(points, dtype=float)
    if len(pts) <= 2:
        return pts.copy()

    keep: set[int] = {0, len(pts) - 1}
    _dp_recursive(pts, list(range(len(pts))), 0, len(pts) - 1, tolerance, keep)
    kept = sorted(keep)
    return pts[kept]
