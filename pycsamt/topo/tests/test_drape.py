# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.topo.drape — terrain-following coordinate transforms.

Covers:
  * interp_elev  — linear / nearest / cubic interpolation and clamping
  * drape_section  — output shapes, flat / varied terrain, exaggeration
  * mask_above_topo  — NaN masking geometry
  * station_surface_z  — surface marker z-coordinates
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.topo.drape import (
    drape_section,
    interp_elev,
    mask_above_topo,
    station_surface_z,
)

# ---------------------------------------------------------------------------
# Helpers shared across test classes
# ---------------------------------------------------------------------------

def _flat_section(nx: int = 5, nz: int = 3):
    """Return (x_nodes, z_nodes, data, flat_elev) for a uniform grid."""
    x_nodes = np.linspace(0.0, 10.0, nx + 1)        # (nx+1,)
    z_nodes = np.linspace(0.0, 2.0, nz + 1)          # (nz+1,)
    data    = np.random.default_rng(0).random((nz, nx))
    elev    = np.full(nx, 1.0)                        # flat terrain at 1 km
    return x_nodes, z_nodes, data, elev


def _varied_section(nx: int = 6, nz: int = 4):
    """Return a section with linearly varying terrain elevation."""
    x_nodes = np.linspace(0.0, 12.0, nx + 1)
    z_nodes = np.linspace(0.0, 3.0, nz + 1)
    data    = np.ones((nz, nx))
    x_ctrs  = (x_nodes[:-1] + x_nodes[1:]) / 2.0
    # Terrain rises from 0.5 to 2.0 km across the profile
    elev    = 0.5 + (x_ctrs / x_ctrs[-1]) * 1.5
    return x_nodes, z_nodes, data, elev


# ---------------------------------------------------------------------------
# interp_elev
# ---------------------------------------------------------------------------

class TestInterpElev:
    def test_linear_midpoint(self):
        c = np.array([0.0, 10.0])
        e = np.array([0.0, 100.0])
        out = interp_elev(c, e, np.array([5.0]))
        assert out[0] == pytest.approx(50.0)

    def test_three_point_interpolation(self):
        c = np.array([0.0, 5.0, 10.0])
        e = np.array([100.0, 200.0, 150.0])
        q = np.array([2.5, 7.5])
        out = interp_elev(c, e, q)
        assert out[0] == pytest.approx(150.0)
        assert out[1] == pytest.approx(175.0)

    def test_extrapolation_clamped_left(self):
        c = np.array([5.0, 10.0])
        e = np.array([80.0, 120.0])
        out = interp_elev(c, e, np.array([0.0]))
        assert out[0] == pytest.approx(80.0)    # clamped to left boundary

    def test_extrapolation_clamped_right(self):
        c = np.array([0.0, 5.0])
        e = np.array([80.0, 120.0])
        out = interp_elev(c, e, np.array([20.0]))
        assert out[0] == pytest.approx(120.0)   # clamped to right boundary

    def test_single_station_constant(self):
        c = np.array([5.0])
        e = np.array([99.0])
        q = np.array([0.0, 2.5, 5.0, 7.5])
        out = interp_elev(c, e, q)
        np.testing.assert_allclose(out, 99.0)

    def test_empty_input_returns_zeros(self):
        out = interp_elev(np.array([]), np.array([]), np.array([1.0, 2.0]))
        np.testing.assert_array_equal(out, np.zeros(2))

    def test_nearest_neighbour(self):
        c = np.array([0.0, 10.0])
        e = np.array([50.0, 150.0])
        q = np.array([3.0])          # closer to c[0]
        out = interp_elev(c, e, q, method="nearest")
        assert out[0] == pytest.approx(50.0)

    def test_nearest_right(self):
        c = np.array([0.0, 10.0])
        e = np.array([50.0, 150.0])
        q = np.array([7.0])          # closer to c[1]
        out = interp_elev(c, e, q, method="nearest")
        assert out[0] == pytest.approx(150.0)

    def test_output_length_matches_query(self):
        c = np.array([0.0, 5.0, 10.0])
        e = np.array([1.0, 2.0, 3.0])
        q = np.linspace(0, 10, 50)
        assert len(interp_elev(c, e, q)) == 50

    def test_exact_knot_recovery(self):
        c = np.array([0.0, 5.0, 10.0])
        e = np.array([100.0, 200.0, 150.0])
        out = interp_elev(c, e, c)
        np.testing.assert_allclose(out, e, rtol=1e-10)

    def test_nan_elevation_ignored(self):
        """Stations with NaN elevation should be skipped in interpolation."""
        c = np.array([0.0, 5.0, 10.0])
        e = np.array([100.0, np.nan, 150.0])
        out = interp_elev(c, e, np.array([0.0, 10.0]))
        assert np.isfinite(out[0])
        assert np.isfinite(out[1])

    @pytest.mark.parametrize("method", ["linear", "nearest"])
    def test_returns_finite_values(self, method):
        c = np.linspace(0, 5, 10)
        e = np.sin(c) * 100 + 200
        out = interp_elev(c, e, np.linspace(0, 5, 30), method=method)
        assert np.all(np.isfinite(out))

    def test_real_willy_l18_chainage(self):
        """End-to-end: interpolate onto fine grid using real WILLY L18 data."""
        import glob
        import os

        from pycsamt.seg.edi import EDIFile
        from pycsamt.topo.extract import (
            extract_chainage,
            extract_elevation,
        )

        base = os.path.join(
            os.path.dirname(__file__),
            "..", "..", "..", "data", "AMT", "WILLY_DATA", "L18PLT",
        )
        paths = sorted(glob.glob(os.path.join(base, "*.edi")))
        if not paths:
            pytest.skip("WILLY L18PLT data not found")
        edis = [EDIFile(p) for p in paths]
        chain = extract_chainage(edis)
        elev  = extract_elevation(edis) / 1000.0    # m → km for drape
        fine  = np.linspace(chain[0], chain[-1], 200)
        out   = interp_elev(chain, elev, fine)
        assert len(out) == 200
        assert np.all(np.isfinite(out))
        assert out.min() > 0.0


# ---------------------------------------------------------------------------
# drape_section
# ---------------------------------------------------------------------------

class TestDrapeSection:
    def test_output_shapes_flat(self):
        xn, zn, d, el = _flat_section(nx=5, nz=3)
        xo, zd, do = drape_section(xn, zn, d, el)
        assert xo.shape == (6,)            # nx+1
        assert zd.shape == (4, 6)          # (nz+1, nx+1)
        assert do.shape == (3, 5)          # (nz, nx)

    def test_output_shapes_varied(self):
        xn, zn, d, el = _varied_section(nx=6, nz=4)
        xo, zd, do = drape_section(xn, zn, d, el)
        assert zd.shape == (5, 7)

    def test_x_nodes_unchanged(self):
        xn, zn, d, el = _flat_section()
        xo, _, _ = drape_section(xn, zn, d, el)
        np.testing.assert_array_equal(xo, xn)

    def test_flat_terrain_top_row_constant(self):
        """With uniform terrain=C and z[0]=0, top node row = C everywhere."""
        xn, zn, d, el = _flat_section()
        elev_km = 0.5
        el[:] = elev_km
        _, zd, _ = drape_section(xn, zn, d, el)
        np.testing.assert_allclose(zd[0, :], elev_km, rtol=1e-10)

    def test_flat_terrain_depth_rows(self):
        """With flat terrain, each row is shifted by the same depth."""
        xn, zn, d, _ = _flat_section(nx=4, nz=3)
        elev_km = 1.2
        el = np.full(4, elev_km)
        _, zd, _ = drape_section(xn, zn, d, el)
        for k in range(len(zn)):
            expected = elev_km - zn[k]
            np.testing.assert_allclose(zd[k, :], expected, rtol=1e-10)

    def test_varied_terrain_columns_differ(self):
        """With non-flat terrain, adjacent columns should have different z."""
        xn, zn, d, el = _varied_section()
        _, zd, _ = drape_section(xn, zn, d, el)
        # Top row (surface) should follow terrain
        assert not np.allclose(zd[0, :-1], zd[0, 1:])

    def test_data_passthrough_no_clip(self):
        xn, zn, d, el = _flat_section()
        _, _, do = drape_section(xn, zn, d, el, clip_above_surface=False)
        np.testing.assert_array_equal(do, d)

    def test_exaggeration_scales_depth(self):
        """Exaggeration of 2 should double the depth range in z_draped."""
        xn  = np.array([0.0, 1.0])
        zn  = np.array([0.0, 1.0, 2.0])
        d   = np.ones((2, 1))
        el  = np.array([0.0])
        _, zd1, _ = drape_section(xn, zn, d, el, exaggeration=1.0)
        _, zd2, _ = drape_section(xn, zn, d, el, exaggeration=2.0)
        # At z_nodes[2]=2.0, draped z_1 = 0 - 2*1 = -2; z_2 = 0 - 2*2 = -4
        assert zd2[-1, 0] == pytest.approx(2 * zd1[-1, 0])

    def test_clip_above_surface_nans_cells(self):
        """clip_above_surface=True must produce NaN in cells above terrain."""
        xn  = np.array([0.0, 1.0])
        zn  = np.array([-1.0, 0.0, 1.0])   # z starts negative — above surface
        d   = np.ones((2, 1))
        el  = np.array([0.0])
        _, _, do = drape_section(xn, zn, d, el, clip_above_surface=True)
        # At least one cell has z_centre < 0 and should be masked
        assert np.any(np.isnan(do))

    def test_output_dtype_float(self):
        xn, zn, d, el = _flat_section()
        _, zd, do = drape_section(xn, zn, d, el)
        assert zd.dtype.kind == "f"
        assert do.dtype.kind == "f"

    def test_single_cell_section(self):
        """Degenerate case: 1 × 1 grid should not crash."""
        xn = np.array([0.0, 1.0])
        zn = np.array([0.0, 1.0])
        d  = np.array([[42.0]])
        el = np.array([0.5])
        xo, zd, do = drape_section(xn, zn, d, el)
        assert zd.shape == (2, 2)
        assert do.shape == (1, 1)

    def test_real_willy_l18_section(self):
        """Build a synthetic inversion grid using L18 chainage/elevation.

        Stations are used as pcolormesh x-node boundaries (nx = n-1 cells).
        elev_at_centres is the elevation interpolated to midpoints between
        adjacent stations.
        """
        import glob
        import os

        from pycsamt.seg.edi import EDIFile
        from pycsamt.topo.extract import (
            extract_chainage,
            extract_elevation,
        )

        base = os.path.join(
            os.path.dirname(__file__),
            "..", "..", "..", "data", "AMT", "WILLY_DATA", "L18PLT",
        )
        paths = sorted(glob.glob(os.path.join(base, "*.edi")))
        if not paths:
            pytest.skip("WILLY L18PLT data not found")
        edis   = [EDIFile(p) for p in paths]
        chain  = extract_chainage(edis)           # (28,) km — station nodes
        elev   = extract_elevation(edis) / 1000.0 # (28,) km

        nx = len(chain) - 1                       # 27 cells
        nz = 20
        x_nodes = chain                           # stations ARE the nodes (nx+1)
        z_nodes = np.linspace(0.0, 1.5, nz + 1)
        rho     = np.random.default_rng(42).random((nz, nx)) * 3.0
        el_c    = (elev[:-1] + elev[1:]) / 2.0   # midpoint elevations (nx)

        xo, zd, do = drape_section(x_nodes, z_nodes, rho,
                                   elev_at_centres=el_c, exaggeration=1.0)
        assert xo.shape == (nx + 1,)
        assert zd.shape == (nz + 1, nx + 1)
        assert do.shape == (nz, nx)
        assert np.all(np.isfinite(zd))


# ---------------------------------------------------------------------------
# mask_above_topo
# ---------------------------------------------------------------------------

class TestMaskAboveTopo:
    def test_output_shape_preserved(self):
        xn = np.linspace(0, 5, 7)
        zn = np.linspace(0, 2, 5)
        d  = np.ones((4, 6))
        el = np.full(6, 1.0)
        out = mask_above_topo(xn, zn, d, el)
        assert out.shape == (4, 6)

    def test_no_masking_when_all_below_surface(self):
        """All cells below the surface — nothing should be NaN."""
        xn = np.array([0.0, 1.0, 2.0])
        zn = np.array([0.0, 0.5, 1.0])
        d  = np.ones((2, 2))
        el = np.array([1.0, 1.0])    # surface at 1 km
        out = mask_above_topo(xn, zn, d, el)
        assert not np.any(np.isnan(out))

    def test_cells_above_max_elev_are_nan(self):
        """Cells elevated above the terrain maximum should be NaN."""
        xn = np.array([0.0, 1.0])
        zn = np.array([-2.0, -1.0, 0.0])   # z_centres at -1.5 and -0.5
        d  = np.ones((2, 1))
        el = np.array([0.0])                # surface at 0
        out = mask_above_topo(xn, zn, d, el)
        # cells at z_centre < 0 → cell_elev = 0 - (-1.5) = +1.5 > max_elev=0
        assert np.isnan(out[0, 0])

    def test_input_data_not_mutated(self):
        xn = np.array([0.0, 1.0, 2.0])
        zn = np.array([0.0, 1.0, 2.0])
        d  = np.ones((2, 2))
        el = np.array([1.0, 1.0])
        d_copy = d.copy()
        mask_above_topo(xn, zn, d, el)
        np.testing.assert_array_equal(d, d_copy)

    def test_dtype_float(self):
        xn = np.linspace(0, 4, 5)
        zn = np.linspace(0, 2, 4)
        d  = np.ones((3, 4), dtype=np.float32)
        el = np.full(4, 0.5)
        out = mask_above_topo(xn, zn, d, el)
        assert out.dtype.kind == "f"


# ---------------------------------------------------------------------------
# station_surface_z
# ---------------------------------------------------------------------------

class TestStationSurfaceZ:
    def test_output_shape(self):
        c = np.linspace(0, 10, 5)
        e = np.array([1.0, 1.2, 1.5, 1.3, 1.1])
        q = np.linspace(0, 10, 20)
        out = station_surface_z(c, e, q)
        assert out.shape == (20,)

    def test_at_known_stations_no_exaggeration(self):
        """At station positions, draped z equals the known elevation."""
        c = np.array([0.0, 5.0, 10.0])
        e = np.array([0.5, 0.8, 0.6])
        out = station_surface_z(c, e, c, exaggeration=1.0)
        np.testing.assert_allclose(out, e, rtol=1e-10)

    def test_exaggeration_scales_output(self):
        c = np.array([0.0, 10.0])
        e = np.array([0.5, 1.0])
        q = np.array([5.0])
        out1 = station_surface_z(c, e, q, exaggeration=1.0)
        out2 = station_surface_z(c, e, q, exaggeration=2.0)
        assert out2[0] == pytest.approx(2.0 * out1[0])

    def test_values_in_elevation_range(self):
        c = np.linspace(0, 5, 10)
        e = np.linspace(0.5, 1.0, 10)
        q = np.linspace(0, 5, 50)
        out = station_surface_z(c, e, q)
        # Should stay within the elevation bounds (clamped extrapolation)
        assert out.min() >= e[0] - 1e-6
        assert out.max() <= e[-1] + 1e-6

    def test_real_willy_l18_markers(self):
        """Station marker z-coords for L18PLT should span the known elev range."""
        import glob
        import os

        from pycsamt.seg.edi import EDIFile
        from pycsamt.topo.extract import (
            extract_chainage,
            extract_elevation,
        )

        base = os.path.join(
            os.path.dirname(__file__),
            "..", "..", "..", "data", "AMT", "WILLY_DATA", "L18PLT",
        )
        paths = sorted(glob.glob(os.path.join(base, "*.edi")))
        if not paths:
            pytest.skip("WILLY L18PLT data not found")
        edis  = [EDIFile(p) for p in paths]
        chain = extract_chainage(edis)
        elev  = extract_elevation(edis) / 1000.0   # m → km
        out   = station_surface_z(chain, elev, chain, exaggeration=1.0)
        np.testing.assert_allclose(out, elev, rtol=1e-10)
        assert out.min() > 0.0
