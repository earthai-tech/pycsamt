# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Triangle-element centroid computation.

Port of ``m2d_getCentroids.m``.

Computes the area-weighted centroid of each region in a MARE2DEM
finite-element mesh, where each region is identified by a contiguous
integer index assigned to mesh triangles.
"""

from __future__ import annotations

import numpy as np

__all__ = ["get_centroids", "triangle_centroids", "triangle_areas"]


def triangle_centroids(
    nodes: np.ndarray,
    elements: np.ndarray,
) -> np.ndarray:
    """Return the centroid of each triangle element.

    Parameters
    ----------
    nodes : array-like, shape (n_nodes, 2)
        Node (y, z) coordinates.
    elements : array-like, shape (n_elements, 3)
        Element connectivity, 0-based or 1-based node indices.
        The function detects 1-based indices automatically.

    Returns
    -------
    numpy.ndarray, shape (n_elements, 2)
        Per-element centroid coordinates.
    """
    n = np.asarray(nodes, dtype=float)
    e = np.asarray(elements, dtype=int)
    # convert 1-based to 0-based if needed
    if e.min() >= 1:
        e = e - 1
    yt = n[e, 0]
    zt = n[e, 1]
    return np.column_stack([yt.mean(axis=1), zt.mean(axis=1)])


def triangle_areas(
    nodes: np.ndarray,
    elements: np.ndarray,
) -> np.ndarray:
    """Return the signed area of each triangle element.

    Parameters
    ----------
    nodes : array-like, shape (n_nodes, 2)
    elements : array-like, shape (n_elements, 3)

    Returns
    -------
    numpy.ndarray, shape (n_elements,)
        Area of each triangle (absolute value).
    """
    n = np.asarray(nodes, dtype=float)
    e = np.asarray(elements, dtype=int)
    if e.min() >= 1:
        e = e - 1
    y = n[:, 0]
    z = n[:, 1]
    yt = y[e]  # (n_el, 3)
    zt = z[e]
    # shoelace formula
    areas = 0.5 * np.abs(
        (yt[:, 1] - yt[:, 0]) * (zt[:, 2] - zt[:, 0])
        - (yt[:, 2] - yt[:, 0]) * (zt[:, 1] - zt[:, 0])
    )
    return areas


def get_centroids(
    nodes: np.ndarray,
    elements: np.ndarray,
    tri_index: np.ndarray,
) -> np.ndarray:
    """Compute area-weighted centroids of mesh regions.

    Port of ``m2d_getCentroids.m``.

    Parameters
    ----------
    nodes : array-like, shape (n_nodes, 2)
        Node (y, z) coordinates.
    elements : array-like, shape (n_elements, 3)
        Element connectivity (1-based indices accepted).
    tri_index : array-like, shape (n_elements,)
        Region index for each triangle (1-based).

    Returns
    -------
    numpy.ndarray, shape (n_regions, 3)
        Columns: y_centroid, z_centroid, total_area of each region.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.models.mare2dem.geom.centroids import get_centroids
    >>> nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    >>> elems = np.array([[1, 2, 3], [2, 4, 3]])
    >>> tri_idx = np.array([1, 1])
    >>> get_centroids(nodes, elems, tri_idx)
    array([[0.5, 0.5, 1. ]])
    """
    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    tri_index = np.asarray(tri_index, dtype=int)

    ctr = triangle_centroids(nodes, elements)
    areas = triangle_areas(nodes, elements)

    region_ids = np.unique(tri_index)
    n_regions = len(region_ids)
    result = np.zeros((n_regions, 3))

    for i, rid in enumerate(region_ids):
        mask = tri_index == rid
        a = areas[mask]
        c = ctr[mask]
        total_area = a.sum()
        if total_area > 0:
            yc = (c[:, 0] * a).sum() / total_area
            zc = (c[:, 1] * a).sum() / total_area
        else:
            yc, zc = c[:, 0].mean(), c[:, 1].mean()
        result[i] = [yc, zc, total_area]

    return result
