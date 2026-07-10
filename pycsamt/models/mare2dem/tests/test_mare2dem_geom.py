# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MARE2DEM geometry utilities.

Covers :mod:`pycsamt.models.mare2dem.geom` — topo, intersections,
simplify, centroids, utm, area_of_interest, triangle_regions,
line_orientation, simplify_poly.
"""

from __future__ import annotations

import math

import numpy as np

from pycsamt.models.mare2dem import (
    do_rects_overlap,
    dp_simplify,
    estimate_area_of_interest,
    get_centroids,
    get_intersections,
    get_line_orientation,
    lonlat_to_utm,
    parse_topo,
    project_onto_line,
    simplify_poly,
    topo_depth,
    topo_slope,
    triangle_areas,
    triangle_centroids,
    utm_to_lonlat,
)
from pycsamt.models.mare2dem.iotools.emdata import (
    EMDataFile,
    MTConfig,
)

# ==========================================================================
# parse_topo
# ==========================================================================

class TestParseTopo:

    def test_flat_topo_constant_depth(self):
        y = np.linspace(-3000, 3000, 20)
        z, slope, on_node = parse_topo(-1500.0, y)
        assert z.shape == (20,)
        np.testing.assert_allclose(z, -1500.0)
        np.testing.assert_allclose(slope, 0.0)
        assert not on_node.any()

    def test_flat_topo_zero_slope(self):
        y = np.array([0.0, 1000.0, 2000.0])
        _, slope, _ = parse_topo(500.0, y)
        np.testing.assert_allclose(slope, 0.0, atol=1e-10)

    def test_linear_topo_slope(self):
        # 45-degree slope: dz/dy = 1 → slope_angle = 45°
        topo = np.column_stack([np.array([0., 1000.]), np.array([0., 1000.])])
        y = np.array([500.0])
        z, slope, on_node = parse_topo(topo, y)
        assert abs(z[0] - 500.0) < 1.0
        assert abs(slope[0] - 45.0) < 1.0

    def test_topo_interpolation_at_nodes(self):
        topo = np.array([[0., 0.], [1000., 100.], [2000., 50.]])
        z, _, on_node = parse_topo(topo, np.array([0., 1000., 2000.]))
        np.testing.assert_allclose(z, [0., 100., 50.], atol=1.0)
        assert on_node.sum() == 3

    def test_topo_depth_convenience(self):
        y = np.array([0.0, 500.0, 1000.0])
        z = topo_depth(-300.0, y)
        np.testing.assert_allclose(z, -300.0)

    def test_topo_slope_convenience(self):
        s = topo_slope(-300.0, np.array([0.0]))
        np.testing.assert_allclose(s, 0.0)

    def test_topo_extrapolation(self):
        topo = np.array([[0., 0.], [1000., 0.]])
        y = np.array([-500., 1500.])
        z, _, _ = parse_topo(topo, y)
        # Should extrapolate (nearest neighbour padding)
        assert np.isfinite(z).all()


# ==========================================================================
# get_intersections / do_rects_overlap
# ==========================================================================

class TestIntersections:

    def test_crossing_segments(self):
        segs_a = np.array([[0.0, 2.0, 1.0, 1.0]])  # horizontal y=1
        seg_b  = np.array([1.0, 1.0, 0.0, 2.0])    # vertical x=1
        inter, xi, yi, pa, pb = get_intersections(segs_a, seg_b)
        assert len(xi) == 1
        assert abs(xi[0] - 1.0) < 1e-9
        assert abs(yi[0] - 1.0) < 1e-9

    def test_parallel_no_intersection(self):
        segs_a = np.array([[0.0, 2.0, 0.0, 0.0]])  # y=0
        seg_b  = np.array([0.0, 2.0, 1.0, 1.0])    # y=1
        inter, xi, yi, _, _ = get_intersections(segs_a, seg_b)
        assert len(xi) == 0

    def test_touching_endpoints_excluded(self):
        # Segments touch at their endpoints → not classified as interior intersection
        segs_a = np.array([[0.0, 1.0, 0.0, 0.0]])
        seg_b  = np.array([1.0, 1.0, 0.0, 1.0])
        inter, xi, yi, _, _ = get_intersections(segs_a, seg_b)
        assert len(xi) == 0

    def test_multiple_segments(self):
        segs_a = np.array([
            [0.0, 2.0, 0.5, 0.5],   # intersects seg_b
            [0.0, 2.0, 2.5, 2.5],   # does not intersect
        ])
        seg_b = np.array([1.0, 1.0, 0.0, 2.0])
        inter, xi, yi, _, _ = get_intersections(segs_a, seg_b)
        assert len(xi) == 1

    def test_do_rects_overlap_true(self):
        a = np.array([0.0, 2.0, 0.0, 2.0])
        many = np.array([[1.0, 3.0, 1.0, 3.0]])
        idx = do_rects_overlap(a, many)
        assert 0 in idx

    def test_do_rects_overlap_false(self):
        a = np.array([0.0, 1.0, 0.0, 1.0])
        many = np.array([[5.0, 6.0, 5.0, 6.0]])
        idx = do_rects_overlap(a, many)
        assert len(idx) == 0


# ==========================================================================
# dp_simplify
# ==========================================================================

class TestDPSimplify:

    def test_straight_line_simplified(self):
        # Perfectly collinear points → only endpoints retained
        pts = np.column_stack([np.linspace(0, 10, 50), np.zeros(50)])
        s = dp_simplify(pts, tolerance=0.01)
        assert len(s) == 2

    def test_sine_wave_simplified(self):
        x = np.linspace(0, 2 * math.pi, 200)
        pts = np.column_stack([x, np.sin(x)])
        s = dp_simplify(pts, tolerance=0.1)
        assert len(s) < 200
        assert len(s) > 2

    def test_preserve_endpoints(self):
        pts = np.column_stack([np.linspace(0, 10, 30), np.random.default_rng(0).random(30)])
        s = dp_simplify(pts, tolerance=0.01)
        np.testing.assert_allclose(s[0], pts[0])
        np.testing.assert_allclose(s[-1], pts[-1])

    def test_zero_tolerance_keeps_all(self):
        pts = np.column_stack([np.linspace(0, 5, 10), np.random.default_rng(1).random(10)])
        s = dp_simplify(pts, tolerance=0.0)
        assert len(s) == len(pts)


# ==========================================================================
# Triangle centroids
# ==========================================================================

class TestCentroids:

    def test_single_triangle_centroid(self):
        nodes = np.array([[0., 0.], [3., 0.], [0., 3.]])
        elems = np.array([[1, 2, 3]])
        c = triangle_centroids(nodes, elems)
        assert c.shape == (1, 2)
        np.testing.assert_allclose(c[0], [1.0, 1.0])

    def test_triangle_area(self):
        nodes = np.array([[0., 0.], [1., 0.], [0., 1.]])
        elems = np.array([[1, 2, 3]])
        a = triangle_areas(nodes, elems)
        assert a.shape == (1,)
        np.testing.assert_allclose(a[0], 0.5)

    def test_get_centroids_one_region(self):
        nodes = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]])
        elems = np.array([[1, 2, 3], [2, 4, 3]])
        tri_idx = np.array([1, 1])
        c = get_centroids(nodes, elems, tri_idx)
        assert c.shape == (1, 3)  # 1 region, (yc, zc, area)
        np.testing.assert_allclose(c[0, 2], 1.0)  # total area = 0.5 + 0.5

    def test_get_centroids_two_regions(self):
        nodes = np.array([[0., 0.], [1., 0.], [2., 0.], [0., 1.], [2., 1.]])
        elems = np.array([[1, 2, 4], [2, 3, 5]])
        tri_idx = np.array([1, 2])
        c = get_centroids(nodes, elems, tri_idx)
        assert c.shape == (2, 3)

    def test_0based_vs_1based_consistent(self):
        nodes = np.array([[0., 0.], [1., 0.], [0., 1.]])
        e0 = np.array([[0, 1, 2]])
        e1 = np.array([[1, 2, 3]])
        c0 = triangle_centroids(nodes, e0)
        c1 = triangle_centroids(nodes, e1)
        np.testing.assert_allclose(c0, c1)


# ==========================================================================
# UTM conversion
# ==========================================================================

class TestUTM:

    def test_lonlat_to_utm_zone(self):
        _, _, zone, sh = lonlat_to_utm(np.array([-70.0]), np.array([42.0]))
        assert zone == 19
        assert not sh

    def test_utm_roundtrip_array(self):
        lon = np.array([-70.0, -71.0, -72.0])
        lat = np.array([42.0, 43.0, 44.0])
        e, n, zone, sh = lonlat_to_utm(lon, lat)
        lon2, lat2 = utm_to_lonlat(e, n, zone, sh)
        np.testing.assert_allclose(lon2, lon, atol=1e-4)
        np.testing.assert_allclose(lat2, lat, atol=1e-4)

    def test_southern_hemisphere(self):
        lon = np.array([150.0])
        lat = np.array([-35.0])
        e, n, zone, sh = lonlat_to_utm(lon, lat)
        assert sh

    def test_forced_zone(self):
        _, _, zone, _ = lonlat_to_utm(np.array([10.0]), np.array([50.0]), zone=32)
        assert zone == 32


# ==========================================================================
# estimate_area_of_interest
# ==========================================================================

class TestAreaOfInterest:

    def _make_mt_em(self, y_positions: np.ndarray) -> EMDataFile:
        em = EMDataFile()
        n = len(y_positions)
        em.mt = MTConfig(
            frequencies=np.array([1.0]),
            receivers=np.zeros((n, 8)),
            receiver_name=[f"S{i}" for i in range(n)],
        )
        em.mt.receivers[:, 1] = y_positions
        return em

    def test_elongated_profile(self):
        em = self._make_mt_em(np.linspace(-5000, 5000, 10))
        ylim, zlim = estimate_area_of_interest(em)
        assert ylim is not None and zlim is not None
        assert ylim[0] < -5000
        assert ylim[1] > 5000

    def test_single_station(self):
        em = self._make_mt_em(np.array([0.0]))
        ylim, zlim = estimate_area_of_interest(em)
        assert ylim is not None
        # 2 km box each side
        np.testing.assert_allclose(ylim, [-1000.0, 1000.0])

    def test_no_survey_returns_none(self):
        em = EMDataFile()
        ylim, zlim = estimate_area_of_interest(em)
        assert ylim is None
        assert zlim is None

    def test_ylim_always_smaller_than_zlim_in_profile(self):
        em = self._make_mt_em(np.linspace(-10000, 10000, 20))
        ylim, zlim = estimate_area_of_interest(em)
        dy = ylim[1] - ylim[0]
        assert dy > 19000  # at least the profile span

    def test_from_real_emdata(self, hill_dir):
        em = __import__(
            "pycsamt.models.mare2dem", fromlist=["read_emdata"]
        ).read_emdata(hill_dir / "hill.emdata")
        ylim, zlim = estimate_area_of_interest(em)
        assert ylim is not None
        assert ylim[0] < -1500
        assert ylim[1] > 1500


# ==========================================================================
# get_line_orientation / project_onto_line
# ==========================================================================

class TestLineOrientation:

    def test_north_south_profile(self):
        n = np.array([0.0, 1000.0, 2000.0])
        e = np.zeros(3)
        ori = get_line_orientation(n, e)
        assert abs(ori) < 1.0   # 0° → N-S

    def test_east_west_profile(self):
        n = np.zeros(3)
        e = np.array([0.0, 1000.0, 2000.0])
        ori = get_line_orientation(n, e)
        assert abs(ori - 90.0) < 1.0  # 90° → E-W

    def test_diagonal_profile(self):
        n = np.array([0.0, 1000.0, 2000.0])
        e = np.array([0.0, 1000.0, 2000.0])
        ori = get_line_orientation(n, e)
        assert abs(ori - 45.0) < 2.0

    def test_project_onto_ns_line(self):
        n = np.array([0.0, 1000.0, 2000.0])
        e = np.zeros(3)
        x, y = project_onto_line(n, e, utm0_north=0.0, utm0_east=0.0,
                                  line_orientation=0.0)
        # Along N-S profile, y should be northing
        np.testing.assert_allclose(y, n, atol=1.0)

    def test_project_preserves_distance(self):
        n = np.array([0.0, 1000.0])
        e = np.array([0.0, 0.0])
        x, y = project_onto_line(n, e, 0.0, 0.0, 0.0)
        dist = np.sqrt(x**2 + y**2)
        np.testing.assert_allclose(dist, [0.0, 1000.0], atol=1.0)


# ==========================================================================
# simplify_poly
# ==========================================================================

class TestSimplifyPoly:

    def test_collinear_chain_reduced(self):
        nodes = np.array([[0., 0.], [1., 0.], [2., 0.], [3., 0.]])
        adj = np.zeros((4, 4))
        adj[0, 1] = adj[1, 0] = 1
        adj[1, 2] = adj[2, 1] = 1
        adj[2, 3] = adj[3, 2] = 1
        n_out, a_out = simplify_poly(nodes, adj)
        assert len(n_out) == 2
        # Endpoints preserved
        np.testing.assert_allclose(n_out[0], [0.0, 0.0])
        np.testing.assert_allclose(n_out[1], [3.0, 0.0])

    def test_non_collinear_preserved(self):
        nodes = np.array([[0., 0.], [1., 0.], [1., 1.]])
        adj = np.zeros((3, 3))
        adj[0, 1] = adj[1, 0] = 1
        adj[1, 2] = adj[2, 1] = 1
        n_out, _ = simplify_poly(nodes, adj)
        assert len(n_out) == 3  # none removed: corner at node 1

    def test_self_loops_removed(self):
        nodes = np.array([[0., 0.], [1., 0.]])
        adj = np.diag([1.0, 1.0])  # self-connected
        adj[0, 1] = adj[1, 0] = 1
        n_out, a_out = simplify_poly(nodes, adj)
        # Self-loops removed; both nodes still connected
        assert np.count_nonzero(a_out) > 0

    def test_isolated_nodes_preserved(self):
        nodes = np.array([[0., 0.], [1., 0.], [5., 5.]])
        adj = np.zeros((3, 3))
        adj[0, 1] = adj[1, 0] = 1
        # Node 2 is isolated
        n_out, _ = simplify_poly(nodes, adj)
        assert len(n_out) == 3  # isolated node kept
