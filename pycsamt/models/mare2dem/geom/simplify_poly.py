# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""PSLG polygon simplification.

Port of ``m2d_simplify_poly.m``.

Cleans up a polygon defined by a node array and an adjacency matrix
of segments by:

1. Removing self-connected segments (diagonal entries).
2. Finding interior nodes that are shared by exactly two collinear
   segments and collapsing them.

This reduces unnecessary nodes that Mamba2D / Triangle would otherwise
mesh at high resolution without geometric benefit.
"""

from __future__ import annotations

import numpy as np

__all__ = ["simplify_poly"]


def _are_collinear(
    p0: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    tol: float = 1e-10,
) -> bool:
    """Return True when p0, p1, p2 are collinear to within *tol*."""
    v1 = p1 - p0
    v2 = p2 - p0
    cross = v1[0] * v2[1] - v1[1] * v2[0]
    norm = max(np.linalg.norm(v1), np.linalg.norm(v2), 1.0)
    return abs(cross) / norm < tol


def simplify_poly(
    nodes: np.ndarray,
    adjacency: np.ndarray,
    *,
    tol: float = 1e-10,
) -> tuple[np.ndarray, np.ndarray]:
    """Remove redundant interior nodes from a PSLG polygon.

    Port of ``m2d_simplify_poly.m``.

    Parameters
    ----------
    nodes : array-like, shape (n_nodes, 2)
        Node (y, z) coordinates.
    adjacency : array-like, shape (n_nodes, n_nodes)
        Dense or sparse symmetric adjacency matrix.  A non-zero
        ``adjacency[i, j]`` means node *i* and node *j* are connected
        by a segment.  Diagonal entries are ignored (self-loops).
    tol : float, default 1e-10
        Collinearity tolerance used to identify redundant nodes.

    Returns
    -------
    nodes_out : numpy.ndarray, shape (m_nodes, 2)
        Pruned node array with redundant interior nodes removed.
    adjacency_out : numpy.ndarray, shape (m_nodes, m_nodes)
        Updated adjacency matrix.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.models.mare2dem.geom.simplify_poly import simplify_poly
    >>> nodes = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    >>> adj = np.zeros((4, 4))
    >>> adj[0, 1] = adj[1, 0] = 1
    ... adj[1, 2] = adj[2, 1] = 1
    ... adj[2, 3] = adj[3, 2] = 1
    >>> n_out, a_out = simplify_poly(nodes, adj)
    >>> len(n_out)  # interior collinear nodes 1 and 2 removed
    2
    """
    nodes = np.asarray(nodes, dtype=float).copy()
    adj = np.asarray(adjacency, dtype=float).copy()

    # remove self-connections
    np.fill_diagonal(adj, 0.0)

    # record nodes that are isolated from the start (holes etc.) — keep them
    initial_degree = np.count_nonzero(adj, axis=1)
    originally_isolated = set(np.where(initial_degree == 0)[0].tolist())

    while True:
        degree = np.count_nonzero(adj, axis=1)
        candidates = np.where(degree == 2)[0]
        removed_any = False

        for idx in candidates:
            if np.count_nonzero(adj[idx]) != 2:
                continue
            j_list = np.where(adj[idx] != 0)[0]
            if len(j_list) != 2:
                continue
            j0, j1 = j_list
            p0, pm, p1 = nodes[j0], nodes[idx], nodes[j1]
            if _are_collinear(p0, pm, p1, tol=tol):
                val = min(adj[idx, j0], adj[idx, j1])
                adj[j0, idx] = 0
                adj[idx, j0] = 0
                adj[j1, idx] = 0
                adj[idx, j1] = 0
                adj[j0, j1] = val
                adj[j1, j0] = val
                removed_any = True

        if not removed_any:
            break

    # keep nodes that still have connections OR were originally isolated
    connected = set(np.where(np.count_nonzero(adj, axis=1) > 0)[0].tolist())
    keep = np.sort(list(connected | originally_isolated))

    # remap indices
    old_to_new = -np.ones(len(nodes), dtype=int)
    old_to_new[keep] = np.arange(len(keep))

    nodes_out = nodes[keep]
    adj_sub = adj[np.ix_(keep, keep)]

    return nodes_out, adj_sub
