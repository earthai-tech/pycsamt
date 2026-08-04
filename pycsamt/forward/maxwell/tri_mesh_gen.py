# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Real, graded, quality triangular FEM mesh generation.

:func:`build_graded_tri_mesh` is the mesh-generation module
:mod:`~pycsamt.forward.maxwell.contracts_tri`'s own docstring points to
(*"Building a TriMesh from real survey geometry belongs to a
mesh-generation module... this module only defines the validated data
shape that generator hands off"*). It replaces
:mod:`~pycsamt.ai.training.dataset2d_tri`'s original
``build_delaunay_mesh``, which triangulated a **uniform rectangular
lattice** of points with plain :func:`scipy.spatial.Delaunay` -- every
triangle the same size everywhere, no coarsening with depth or lateral
distance from the receivers. That produced a mesh that looked nothing
like a real FEM discretization (fine near the sources of interest,
growing smoothly away from them), and is not what
:class:`~pycsamt.forward.maxwell.tri_fem2d.TriFEM2DAdapter` actually
needs to solve accurately -- see that module's own note on vertical
resolution sensitivity near receivers.

This module instead uses `triangle
<https://rufat.be/triangle/>`_, the Python bindings for J. R.
Shewchuk's Triangle library (the same battle-tested mesh generator
already used, via subprocess, by
:mod:`pycsamt.models.mare2dem.triangle_exec` for the compiled-MARE2DEM
route) -- but called in-process, needing no external binary:

1. A PSLG (rectangle boundary + station positions as required
   vertices) is triangulated once with Triangle's own quality flag
   (``pq<min_angle>``), giving a real, proven min-angle-bounded
   (Ruppert's algorithm) conforming Delaunay triangulation of the
   domain.
2. A smooth **size function** grades a per-triangle target area from
   that first pass: the target edge length at a point grows
   geometrically with its distance to the *nearest station*,
   ``h(x,z) = surface_cell_m * growth_rate ** (d(x,z) / surface_cell_m)``,
   capped at ``max_cell_m``. This is a continuous generalization of the
   classic "N surrounding refinement layers around a receiver" scheme
   (see e.g. Usui (2015)'s hybrid-mesh figure) to an arbitrary number of
   receivers along a profile, rather than hand-building discrete
   nested regions.
3. Triangle's own adaptive refine mode (``pq<min_angle>ra``, with a
   ``triangle_max_area`` value attached per first-pass triangle) grows
   the mesh to match that target, still honoring the same min-angle
   bound.

``boundary_segments`` on the returned :class:`TriMesh` is derived
topologically -- any triangle edge belonging to exactly one triangle is
a boundary edge -- rather than from ``scipy.spatial.Delaunay``'s
``convex_hull`` (which only ever happened to work on the old uniform
lattice; once Triangle inserts Steiner points along a PSLG segment
during quality refinement, the convex hull of the node set is not
guaranteed to trace every boundary node). This matters beyond cosmetics:
:func:`~pycsamt.forward.maxwell.tri_fem2d._boundary_node_ids` reads
``mesh.boundary_segments`` directly to decide which nodes get the
Dirichlet boundary condition, so an incomplete boundary would silently
corrupt the FEM solve, not just the picture.

Every station is included as an explicit PSLG vertex, so it always sits
exactly on a mesh node at the surface -- the same convention the old
``build_delaunay_mesh`` established and
:class:`~pycsamt.forward.maxwell.tri_fem2d.TriFEM2DAdapter` /
:class:`~pycsamt.forward.maxwell.mare2dem.Mare2DEMAdapter` both already
depend on. By default that surface is flat at ``z_range_m[0]``; passing
``topo_x_m``/``topo_z_m`` builds the top boundary along a real
topography polyline instead (stations placed at their true interpolated
elevation, not ``z=0``) -- see :class:`~pycsamt.forward.maxwell.tri_fem2d.TriFEM2DAdapter`'s
own module docstring for the matching solver-side generalization
(``_local_surface_z``) this depends on.
"""

from __future__ import annotations

import numpy as np
import triangle as tr

from .contracts_tri import TriMesh

__all__ = ["build_graded_tri_mesh"]


def _boundary_edges(triangles: np.ndarray) -> np.ndarray:
    """Return every triangle edge shared by exactly one triangle.

    Purely topological, so it is correct for any manifold triangulation
    regardless of how many Steiner points quality refinement inserted
    along a boundary segment -- unlike relying on a convex hull, which
    is only correct by coincidence on an unrefined convex point set.

    Examples
    --------
    >>> _boundary_edges(np.array([[0, 1, 2]])).tolist()
    [[0, 1], [1, 2], [0, 2]]
    """
    edge_count: dict[tuple[int, int], int] = {}
    for triangle in triangles:
        a, b, c = int(triangle[0]), int(triangle[1]), int(triangle[2])
        for u, v in ((a, b), (b, c), (c, a)):
            key = (u, v) if u < v else (v, u)
            edge_count[key] = edge_count.get(key, 0) + 1
    return np.array(
        [edge for edge, count in edge_count.items() if count == 1],
        dtype=np.int64,
    )


def _graded_target_area(
    centroids: np.ndarray,
    station_x_m: np.ndarray,
    station_z_m: np.ndarray,
    *,
    surface_cell_m: float,
    growth_rate: float,
    max_cell_m: float,
) -> np.ndarray:
    """Return the graded max-area constraint at every centroid."""
    station_points = np.column_stack([station_x_m, station_z_m])
    distance = np.min(
        np.linalg.norm(
            centroids[:, None, :] - station_points[None, :, :], axis=2
        ),
        axis=1,
    )
    edge_length = np.minimum(
        surface_cell_m * np.power(growth_rate, distance / surface_cell_m),
        max_cell_m,
    )
    # Equilateral-triangle area for the graded target edge length.
    return (np.sqrt(3.0) / 4.0) * edge_length**2


def build_graded_tri_mesh(
    x_range_m: tuple[float, float],
    z_range_m: tuple[float, float],
    station_x_m,
    *,
    surface_cell_m: float,
    growth_rate: float = 1.3,
    max_cell_m: float | None = None,
    min_angle: float = 30.0,
    topo_x_m=None,
    topo_z_m=None,
) -> TriMesh:
    """Build a real, graded, quality FEM mesh via Shewchuk's Triangle.

    Parameters
    ----------
    x_range_m, z_range_m : (float, float)
        Domain extent. ``z_range_m[0]`` should be 0 (surface) when
        ``topo_z_m`` is omitted; ignored (the top follows the topography
        instead) when it is given. ``z_range_m[1]`` (the domain bottom)
        must always be deeper than every ``topo_z_m`` sample.
    station_x_m : array-like
        Receiver x-positions; included as explicit PSLG vertices so
        every receiver sits exactly on a mesh node. Placed at ``z=0``
        unless ``topo_z_m`` is given, in which case each station's z is
        interpolated from the topography polyline instead.
    topo_x_m, topo_z_m : array-like, optional
        Topography polyline (z positive down, so a ridge is negative)
        sampled at ``topo_x_m``. When given, the mesh's top boundary
        follows this polyline (merged with ``station_x_m``, sorted and
        deduplicated by x, each station's elevation interpolated from
        it) instead of a flat surface. Both must be given together, the
        same length, and finite. This only builds the *mesh*; using
        stations away from ``z=0`` for a real solve also needs
        :class:`~pycsamt.forward.maxwell.tri_fem2d.TriFEM2DAdapter`'s
        local-surface-aware boundary conditions (already the default
        there -- see that module's own docstring).
    surface_cell_m : float
        Target triangle edge length at a station (the finest part of
        the mesh). Smaller values give a finer, more expensive mesh
        near the receivers. Choose this relative to skin depth for FEM
        solver accuracy, not just picture quality: a first layer around
        ~0.03-0.05 skin depths (at the highest simulated frequency) is
        what
        :class:`~pycsamt.forward.maxwell.tri_fem2d.TriFEM2DAdapter`'s
        own analytic half-space benchmark needs to stay under 5%
        relative impedance error (see that module's "Station field
        extraction" docstring section for the same sensitivity on a
        structured mesh).
    growth_rate : float, default=1.3
        Geometric growth factor of the target edge length per
        ``surface_cell_m`` of distance from the nearest station. Must
        be greater than 1 (grading away from the receivers, never
        toward them).
    max_cell_m : float or None, optional
        Upper bound on the graded target edge length, reached far from
        every station / at depth. Defaults to ``25 * surface_cell_m``.
    min_angle : float, default=30.0
        Minimum triangle angle in degrees, Triangle's own quality
        constraint (same default as
        :func:`pycsamt.models.mare2dem.triangle_exec.run_triangle`).

    Returns
    -------
    TriMesh
        Graded mesh with one region per triangle
        (``region_ids = 1..n_triangles``, matching the old
        ``build_delaunay_mesh`` contract so
        :class:`~pycsamt.forward.maxwell.mare2dem.Mare2DEMAdapter`'s
        per-triangle-region convention still holds), ready for a
        per-triangle heterogeneous resistivity field.

    Raises
    ------
    ValueError
        If the ranges are not increasing, a size parameter is
        non-positive, or ``station_x_m`` is empty or outside
        ``x_range_m``.

    Examples
    --------
    >>> mesh = build_graded_tri_mesh(
    ...     (0.0, 1000.0), (0.0, 500.0), [200.0, 500.0, 800.0],
    ...     surface_cell_m=20.0,
    ... )
    >>> mesh.n_triangles > 0
    True
    >>> sorted(mesh.region_ids.tolist()) == list(
    ...     range(1, mesh.n_triangles + 1)
    ... )
    True
    """
    x0, x1 = float(x_range_m[0]), float(x_range_m[1])
    z0, z1 = float(z_range_m[0]), float(z_range_m[1])
    if x1 <= x0 or z1 <= z0:
        raise ValueError("x_range_m/z_range_m must be increasing.")

    if (topo_x_m is None) != (topo_z_m is None):
        raise ValueError("topo_x_m and topo_z_m must be given together.")
    topo_x = topo_z = None
    if topo_x_m is not None:
        topo_x = np.asarray(topo_x_m, dtype=float).ravel()
        topo_z = np.asarray(topo_z_m, dtype=float).ravel()
        if (
            topo_x.shape != topo_z.shape
            or topo_x.size < 2
            or not np.all(np.isfinite(topo_x))
            or not np.all(np.isfinite(topo_z))
        ):
            raise ValueError(
                "topo_x_m/topo_z_m must be finite, equal-length, and "
                "have at least two samples."
            )
        order = np.argsort(topo_x)
        topo_x, topo_z = topo_x[order], topo_z[order]
        if z1 <= float(topo_z.max()):
            raise ValueError(
                "z_range_m[1] (domain bottom) must be deeper than every "
                "topo_z_m sample."
            )

    surface_cell = float(surface_cell_m)
    if not np.isfinite(surface_cell) or surface_cell <= 0:
        raise ValueError("surface_cell_m must be finite and positive.")

    growth = float(growth_rate)
    if not np.isfinite(growth) or growth <= 1.0:
        raise ValueError("growth_rate must be finite and greater than 1.")

    max_cell = 25.0 * surface_cell if max_cell_m is None else float(max_cell_m)
    if not np.isfinite(max_cell) or max_cell < surface_cell:
        raise ValueError("max_cell_m must be finite and >= surface_cell_m.")

    angle = float(min_angle)
    if not np.isfinite(angle) or not (0.0 < angle <= 34.0):
        raise ValueError(
            "min_angle must be finite and in (0, 34] degrees "
            "(Triangle's own practical quality ceiling)."
        )

    station_x = np.asarray(station_x_m, dtype=float).ravel()
    if (
        station_x.size < 1
        or not np.all(np.isfinite(station_x))
        or np.any(station_x < x0)
        or np.any(station_x > x1)
    ):
        raise ValueError(
            f"station_x_m must be non-empty, finite, and within "
            f"x_range_m {(x0, x1)}."
        )

    # The top edge (where every station sits) is built as an explicit
    # chain of sub-segments through the stations, sorted by x, rather
    # than one long segment plus unconnected points that merely happen
    # to lie on it. An unconstrained vertex resting in a segment's
    # *interior* is legal PSLG input, but numerically fragile under
    # quality refinement -- confirmed empirically: it produced a
    # "ran out of precision" error (and, with more such points, a hard
    # segfault inside Triangle's C extension) while adapting this same
    # technique to a topography-following boundary in
    # docs/scripts/generate_ai_inversion_figures.py. Explicitly splitting
    # the edge at each station sidesteps the whole class of failure.
    if topo_x is None:
        interior_x = np.unique(station_x[(station_x > x0) & (station_x < x1)])
        chain_x = np.r_[x0, interior_x, x1]
        chain_z = np.full(chain_x.size, z0)
        station_z = np.zeros_like(station_x)
    else:
        chain_x = np.unique(
            np.concatenate([[x0, x1], topo_x, station_x])
        )
        chain_x = chain_x[(chain_x >= x0) & (chain_x <= x1)]
        chain_z = np.interp(chain_x, topo_x, topo_z)
        station_z = np.interp(station_x, topo_x, topo_z)
    top_chain = np.column_stack([chain_x, chain_z])
    other_corners = np.array([[x1, z1], [x0, z1]], dtype=float)
    vertices = np.vstack([top_chain, other_corners])

    n_total = len(vertices)
    segments = np.column_stack(
        [np.arange(n_total), np.r_[np.arange(1, n_total), 0]]
    )

    pslg = {"vertices": vertices, "segments": segments}
    base = tr.triangulate(pslg, f"pq{angle:g}")

    base_triangles = np.asarray(base["triangles"], dtype=np.int64)
    base_nodes = np.asarray(base["vertices"], dtype=float)
    centroids = base_nodes[base_triangles].mean(axis=1)
    base["triangle_max_area"] = _graded_target_area(
        centroids,
        station_x,
        station_z,
        surface_cell_m=surface_cell,
        growth_rate=growth,
        max_cell_m=max_cell,
    )

    refined = tr.triangulate(base, f"pq{angle:g}ra")

    nodes = np.asarray(refined["vertices"], dtype=float)
    triangles = np.asarray(refined["triangles"], dtype=np.int64)

    return TriMesh(
        nodes,
        triangles,
        region_ids=np.arange(1, len(triangles) + 1),
        boundary_segments=_boundary_edges(triangles),
    )
