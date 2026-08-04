# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.forward.maxwell.tri_mesh_gen`.

Exercises the real Shewchuk-Triangle-backed graded mesh builder that
replaced the old uniform-lattice ``scipy.spatial.Delaunay`` approach
(``pycsamt.ai.training.dataset2d_tri.build_delaunay_mesh``, removed).
No scripted stand-in is used here -- ``triangle`` is a pure Python
C-extension with prebuilt wheels, not an external binary, so the real
library is exercised directly.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.forward.maxwell.contracts_tri import TriMesh
from pycsamt.forward.maxwell.tri_mesh_gen import (
    _boundary_edges,
    build_graded_tri_mesh,
)

STATIONS = [50.0 * i for i in range(6)]


def _mesh(**overrides):
    kwargs = dict(
        x_range_m=(-25.0, 275.0),
        z_range_m=(0.0, 300.0),
        station_x_m=STATIONS,
        surface_cell_m=20.0,
    )
    kwargs.update(overrides)
    return build_graded_tri_mesh(**kwargs)


# ---------------------------------------------------------------------------
# Basic contract
# ---------------------------------------------------------------------------


def test_returns_a_valid_trimesh():
    mesh = _mesh()
    assert isinstance(mesh, TriMesh)
    assert mesh.n_triangles > 0
    assert mesh.dimension == 2


def test_stations_land_exactly_on_mesh_nodes():
    mesh = _mesh()
    for station_x in STATIONS:
        assert np.any(
            np.isclose(mesh.nodes_m[:, 0], station_x)
            & np.isclose(mesh.nodes_m[:, 1], 0.0)
        )


def test_one_region_per_triangle():
    mesh = _mesh()
    assert sorted(mesh.region_ids.tolist()) == list(
        range(1, mesh.n_triangles + 1)
    )


# ---------------------------------------------------------------------------
# Grading (the actual point of this module)
# ---------------------------------------------------------------------------


def test_mesh_grows_with_depth():
    mesh = _mesh(surface_cell_m=15.0, growth_rate=1.3)
    centroids = mesh.triangle_centroids_m
    areas = mesh.triangle_areas_m2
    shallow = areas[centroids[:, 1] < 20.0]
    deep = areas[centroids[:, 1] > 250.0]
    assert shallow.size > 0 and deep.size > 0
    assert deep.mean() > 5.0 * shallow.mean()


def test_higher_growth_rate_gives_fewer_total_triangles():
    slow = _mesh(growth_rate=1.1)
    fast = _mesh(growth_rate=1.6)
    assert fast.n_triangles < slow.n_triangles


def test_smaller_surface_cell_gives_more_triangles():
    coarse = _mesh(surface_cell_m=40.0)
    fine = _mesh(surface_cell_m=10.0)
    assert fine.n_triangles > coarse.n_triangles


# ---------------------------------------------------------------------------
# boundary_segments correctness (load-bearing for tri_fem2d's Dirichlet BC)
# ---------------------------------------------------------------------------


def test_boundary_segments_cover_the_full_rectangle_perimeter():
    x0, x1 = -25.0, 275.0
    z0, z1 = 0.0, 300.0
    mesh = _mesh(x_range_m=(x0, x1), z_range_m=(z0, z1))

    boundary_nodes = np.unique(mesh.boundary_segments.ravel())
    coords = mesh.nodes_m[boundary_nodes]
    on_perimeter = (
        np.isclose(coords[:, 0], x0)
        | np.isclose(coords[:, 0], x1)
        | np.isclose(coords[:, 1], z0)
        | np.isclose(coords[:, 1], z1)
    )
    assert on_perimeter.all()

    # Every corner must be present.
    for corner in ((x0, z0), (x1, z0), (x1, z1), (x0, z1)):
        assert np.any(
            np.isclose(mesh.nodes_m[:, 0], corner[0])
            & np.isclose(mesh.nodes_m[:, 1], corner[1])
        )


def test_boundary_segments_match_topological_ground_truth():
    mesh = _mesh()
    expected = _boundary_edges(mesh.triangles)
    got = np.asarray(mesh.boundary_segments)

    def _edge_set(edges):
        return {tuple(sorted((int(a), int(b)))) for a, b in edges}

    assert _edge_set(got) == _edge_set(expected)


def test_boundary_edges_topological_helper_direct():
    # Single triangle: every one of its 3 edges is a boundary edge.
    assert sorted(map(tuple, _boundary_edges(np.array([[0, 1, 2]])).tolist())) == [
        (0, 1),
        (0, 2),
        (1, 2),
    ]
    # Two triangles sharing edge (1, 2): only the 4 outer edges remain.
    shared = _boundary_edges(np.array([[0, 1, 2], [1, 3, 2]]))
    assert len(shared) == 4
    assert not any(set(e) == {1, 2} for e in shared.tolist())


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_rejects_bad_ranges():
    with pytest.raises(ValueError, match="increasing"):
        _mesh(x_range_m=(100.0, 0.0))
    with pytest.raises(ValueError, match="increasing"):
        _mesh(z_range_m=(100.0, 0.0))


def test_rejects_nonpositive_surface_cell():
    with pytest.raises(ValueError, match="surface_cell_m"):
        _mesh(surface_cell_m=0.0)


def test_rejects_growth_rate_not_greater_than_one():
    with pytest.raises(ValueError, match="growth_rate"):
        _mesh(growth_rate=1.0)


def test_rejects_max_cell_smaller_than_surface_cell():
    with pytest.raises(ValueError, match="max_cell_m"):
        _mesh(surface_cell_m=20.0, max_cell_m=10.0)


def test_rejects_stations_outside_domain():
    with pytest.raises(ValueError, match="station_x_m"):
        _mesh(station_x_m=[10.0, 9000.0])


def test_rejects_empty_stations():
    with pytest.raises(ValueError, match="station_x_m"):
        _mesh(station_x_m=[])


def test_rejects_bad_min_angle():
    with pytest.raises(ValueError, match="min_angle"):
        _mesh(min_angle=0.0)
    with pytest.raises(ValueError, match="min_angle"):
        _mesh(min_angle=40.0)


# ---------------------------------------------------------------------------
# Station on a domain corner (dedup path)
# ---------------------------------------------------------------------------


def test_station_coincident_with_domain_corner_does_not_error():
    mesh = build_graded_tri_mesh(
        (0.0, 200.0), (0.0, 100.0), [0.0, 100.0, 200.0], surface_cell_m=20.0
    )
    for station_x in (0.0, 100.0, 200.0):
        assert np.any(
            np.isclose(mesh.nodes_m[:, 0], station_x)
            & np.isclose(mesh.nodes_m[:, 1], 0.0)
        )


# ---------------------------------------------------------------------------
# Topography (topo_x_m/topo_z_m)
# ---------------------------------------------------------------------------


def _ridge(x):
    return -80.0 * np.exp(-(((x - 1100.0) / 500.0) ** 2)) - 20.0 * np.sin(x / 600.0)


_TOPO_STATION_X = np.linspace(200.0, 1800.0, 9)
_TOPO_X = np.linspace(0.0, 2000.0, 60)
_TOPO_Z = _ridge(_TOPO_X)


def _topo_mesh(**overrides):
    kwargs = dict(
        x_range_m=(0.0, 2000.0),
        z_range_m=(0.0, 700.0),
        station_x_m=_TOPO_STATION_X,
        surface_cell_m=40.0,
        topo_x_m=_TOPO_X,
        topo_z_m=_TOPO_Z,
    )
    kwargs.update(overrides)
    return kwargs["station_x_m"], build_graded_tri_mesh(**kwargs)


def test_topo_stations_land_on_the_topography_at_their_true_elevation():
    station_x, mesh = _topo_mesh()
    # The mesh places a station at its *piecewise-linear* interpolation of
    # the discrete topo_x_m/topo_z_m samples -- the same quantity
    # build_graded_tri_mesh itself computes -- not the continuous ridge()
    # function evaluated directly (which differs by curvature error
    # between sample points).
    station_z = np.interp(station_x, _TOPO_X, _TOPO_Z)
    for sx, sz in zip(station_x, station_z):
        assert np.any(
            np.isclose(mesh.nodes_m[:, 0], sx) & np.isclose(mesh.nodes_m[:, 1], sz)
        )
    # None of the stations sit at z=0 here -- proves this isn't silently
    # falling back to the flat-mode path.
    assert not np.any(np.isclose(station_z, 0.0))


def test_topo_mesh_still_grades_with_depth():
    _, mesh = _topo_mesh()
    centroids = mesh.triangle_centroids_m
    areas = mesh.triangle_areas_m2
    shallow = areas[centroids[:, 1] < 20.0]
    deep = areas[centroids[:, 1] > 400.0]
    assert shallow.size > 0 and deep.size > 0
    assert deep.mean() > 10.0 * shallow.mean()


def test_topo_boundary_segments_cover_the_full_perimeter():
    _, mesh = _topo_mesh()
    boundary_nodes = np.unique(mesh.boundary_segments.ravel())
    # Every boundary node must be on the topo line (its piecewise-linear
    # interpolation, matching what the mesh builder itself computes), the
    # two vertical sides (x=0 or x=2000), or the flat bottom (z=700) --
    # never a stray interior node.
    coords = mesh.nodes_m[boundary_nodes]
    on_topo = np.isclose(
        coords[:, 1], np.interp(coords[:, 0], _TOPO_X, _TOPO_Z), atol=1e-6
    )
    on_side = np.isclose(coords[:, 0], 0.0) | np.isclose(coords[:, 0], 2000.0)
    on_bottom = np.isclose(coords[:, 1], 700.0)
    assert np.all(on_topo | on_side | on_bottom)


def test_topo_flat_mode_is_unaffected_when_topo_omitted():
    # Same call as the very first test in this file (flat, no topo args)
    # must be byte-for-byte unaffected by the topo feature's addition.
    mesh = build_graded_tri_mesh(
        (-25.0, 275.0), (0.0, 300.0), STATIONS, surface_cell_m=20.0
    )
    assert mesh.n_triangles > 0
    for station_x in STATIONS:
        assert np.any(
            np.isclose(mesh.nodes_m[:, 0], station_x)
            & np.isclose(mesh.nodes_m[:, 1], 0.0)
        )


def test_rejects_topo_given_only_partially():
    with pytest.raises(ValueError, match="topo_x_m and topo_z_m"):
        _mesh(topo_x_m=[0.0, 100.0])


def test_rejects_topo_arrays_mismatched_or_too_short():
    with pytest.raises(ValueError, match="topo_x_m/topo_z_m"):
        _mesh(topo_x_m=[0.0, 100.0, 200.0], topo_z_m=[0.0, 0.0])
    with pytest.raises(ValueError, match="topo_x_m/topo_z_m"):
        _mesh(topo_x_m=[0.0], topo_z_m=[0.0])


def test_rejects_domain_bottom_shallower_than_topo():
    # default z_range_m=(0.0, 300.0); topo sitting at z=400 is deeper than
    # the domain bottom, which must be rejected rather than silently
    # producing a mesh with topography below its own floor.
    topo_x = np.array([-25.0, 275.0])
    topo_z = np.array([400.0, 400.0])
    with pytest.raises(ValueError, match="z_range_m\\[1\\]"):
        _mesh(topo_x_m=topo_x, topo_z_m=topo_z)
