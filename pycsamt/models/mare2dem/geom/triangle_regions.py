# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Assign FEM mesh triangles to PSLG-bounded regions.

Port of ``m2d_getTriangleRegions.m``.

Given a Delaunay triangulation, a set of boundary segments that divide
the mesh into regions, and a list of seed points (one per region),
this function flood-fills from each seed to label all triangles in the
same connected component.  Triangles that straddle a boundary segment
are treated as walls that stop the flood fill.

Dependencies
------------
scipy.spatial.Delaunay (for the triangulation object)
scipy.sparse (for the adjacency matrix)
"""

from __future__ import annotations

import numpy as np

__all__ = ["get_triangle_regions"]


def get_triangle_regions(
    points: np.ndarray,
    triangles: np.ndarray,
    segments: np.ndarray,
    region_seeds: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Assign each triangle to a numbered region.

    Port of ``m2d_getTriangleRegions.m``.

    Parameters
    ----------
    points : array-like, shape (n_nodes, 2)
        Node (y, z) coordinates.
    triangles : array-like, shape (n_tri, 3)
        Triangle connectivity, 1-based node indices.
    segments : array-like, shape (n_segs, 2)
        Boundary segment node-index pairs (1-based) that separate
        regions.
    region_seeds : array-like, shape (n_seeds, 2) or None
        One seed point per pre-defined region.  Each seed is a (y, z)
        coordinate inside a region.  Pass ``None`` if no pre-defined
        seeds are given; all regions are then discovered automatically.

    Returns
    -------
    tri_index : numpy.ndarray, shape (n_tri,)
        Region index for each triangle (1-based).
    region_map : numpy.ndarray, shape (n_new_regions,)
        Maps new region indices back to input seed indices (1-based);
        0 for regions with no corresponding seed.

    Notes
    -----
    The algorithm builds a sparse adjacency matrix from the boundary
    segments, zeros out neighbour links that cross boundaries, then
    performs a BFS flood fill from each seed (or from untouched
    triangles after the seed pass).

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.models.mare2dem.geom.triangle_regions import get_triangle_regions
    >>> pts = np.array([[0.,0.],[1.,0.],[0.5,1.]])
    >>> tris = np.array([[1,2,3]])
    >>> segs = np.array([[1,2]])
    >>> tri_index, region_map = get_triangle_regions(pts, tris, segs)
    >>> tri_index
    array([1])
    """
    pts = np.asarray(points, dtype=float)
    tris = np.asarray(triangles, dtype=int)
    segs = np.asarray(segments, dtype=int)

    # convert to 0-based
    if tris.min() >= 1:
        tris = tris - 1
    if segs.min() >= 1:
        segs = segs - 1

    n_nodes = len(pts)
    n_tris = len(tris)

    # ---- build segment adjacency (boundary edges) ----
    try:
        from scipy.sparse import lil_matrix
        A = lil_matrix((n_nodes, n_nodes), dtype=np.int8)
        for s in segs:
            A[s[0], s[1]] = 1
            A[s[1], s[0]] = 1
        A = A.tocsr()
    except ImportError:
        # fallback: dense set of frozensets
        _segs_set: set[frozenset] = {frozenset((s[0], s[1])) for s in segs}
        A = None

    def _is_boundary(a: int, b: int) -> bool:
        if A is not None:
            return bool(A[a, b])
        return frozenset((a, b)) in _segs_set  # type: ignore[union-attr]

    # ---- build triangle-neighbour adjacency ----
    # edge → triangle index lookup
    edge_to_tri: dict[tuple[int, int], list[int]] = {}
    for ti, tri in enumerate(tris):
        for i in range(3):
            e = (tri[i], tri[(i + 1) % 3])
            e_sorted = (min(e), max(e))
            edge_to_tri.setdefault(e_sorted, []).append(ti)

    # neighbour list: for each tri, which tris share an edge that is not a boundary
    neighbours: list[list[int]] = [[] for _ in range(n_tris)]
    for (a, b), tlist in edge_to_tri.items():
        if len(tlist) == 2 and not _is_boundary(a, b):
            neighbours[tlist[0]].append(tlist[1])
            neighbours[tlist[1]].append(tlist[0])

    tri_index = np.zeros(n_tris, dtype=int)
    n_regions = 0

    # ---- find which triangle contains each seed ----
    seed_tris: list[int] = []
    if region_seeds is not None and len(region_seeds):
        seeds = np.asarray(region_seeds, dtype=float)
        try:
            from scipy.spatial import Delaunay as _DT
            # build Delaunay from pts to query point locations
            _dt = _DT(pts)
            for seed in seeds:
                ti = _dt.find_simplex(seed)
                seed_tris.append(int(ti))
        except Exception:
            seed_tris = [-1] * len(seeds)

    # ---- flood fill from seeds ----
    region_map: list[int] = []
    for seed_idx, start_tri in enumerate(seed_tris):
        if start_tri < 0 or tri_index[start_tri] != 0:
            continue
        n_regions += 1
        region_map.append(seed_idx + 1)
        queue = [start_tri]
        head = 0
        while head < len(queue):
            ti = queue[head]
            head += 1
            if tri_index[ti] != 0:
                continue
            tri_index[ti] = n_regions
            for nb in neighbours[ti]:
                if tri_index[nb] == 0:
                    queue.append(nb)

    # ---- flood fill remaining unassigned triangles ----
    for ti in range(n_tris):
        if tri_index[ti] != 0:
            continue
        n_regions += 1
        region_map.append(0)
        queue = [ti]
        head = 0
        while head < len(queue):
            t = queue[head]
            head += 1
            if tri_index[t] != 0:
                continue
            tri_index[t] = n_regions
            for nb in neighbours[t]:
                if tri_index[nb] == 0:
                    queue.append(nb)

    return tri_index, np.array(region_map, dtype=int)
