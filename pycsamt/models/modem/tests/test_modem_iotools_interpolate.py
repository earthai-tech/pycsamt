"""Tests for pycsamt.models.modem.iotools.interpolate."""

import numpy as np
import pytest

pytest.importorskip("scipy")

from pycsamt.models.modem.iotools.impedance import (
    ImpedanceFile,
    ZBlock,
)
from pycsamt.models.modem.iotools.interpolate import (
    interp_model3d,
    interp_z3d,
)

# ---------------------------------------------------------------------------
# Synthetic model / data builders
# ---------------------------------------------------------------------------

def _uniform_model(nx=6, ny=5, nz=4, n_air=1, rho=100.0):
    """Small uniform-resistivity model."""
    from pycsamt.models.modem.model3d import ModEmModel3D
    m = ModEmModel3D()
    m.x_widths = np.full(nx, 2000.0)
    m.y_widths = np.full(ny, 2000.0)
    m.z_widths = np.concatenate([np.full(n_air, 50.0),
                                  np.geomspace(100.0, 1000.0, nz - n_air)])
    m.n_air = n_air
    m.rho_loge = np.full((nz, ny, nx), np.log(rho))
    return m


def _gradient_model(nx=8, ny=6, nz=5, n_air=1):
    """Model with a smooth resistivity gradient (easy to verify interpolation)."""
    from pycsamt.models.modem.model3d import ModEmModel3D
    m = ModEmModel3D()
    m.x_widths = np.full(nx, 1000.0)
    m.y_widths = np.full(ny, 1000.0)
    m.z_widths = np.concatenate([np.full(n_air, 50.0),
                                  np.full(nz - n_air, 500.0)])
    m.n_air = n_air
    # Linear rho gradient 100 → 1000 Ω·m along z-axis
    rho_vals = np.linspace(100.0, 1000.0, nz)
    m.rho_loge = np.log(rho_vals[:, None, None]) * np.ones((nz, ny, nx))
    return m


def _z3d_imp_on_grid(n_pts_per_dim=4, periods=(10.0, 100.0)):
    """ImpedanceFile with sites on a regular 2-D grid."""
    rng = np.random.default_rng(42)
    xs = np.linspace(-5000, 5000, n_pts_per_dim)
    ys = np.linspace(-5000, 5000, n_pts_per_dim)
    xx, yy = np.meshgrid(xs, ys)
    site_loc = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(n_pts_per_dim**2)])
    n_sites = len(site_loc)
    comps = ["ZXX", "ZXY", "ZYX", "ZYY"]
    blocks = []
    for T in periods:
        Z    = rng.standard_normal((n_sites, 4)) + 1j * rng.standard_normal((n_sites, 4))
        Zerr = np.abs(rng.standard_normal((n_sites, 4))) * 0.1 + 1e-3
        blocks.append(ZBlock(
            period=T, site_names=[f"S{k:03d}" for k in range(n_sites)],
            site_loc=site_loc, comp_names=comps, Z=Z, Zerr=Zerr,
        ))
    return ImpedanceFile(blocks=blocks)


# ===========================================================================
# interp_model3d
# ===========================================================================

class TestInterpModel3D:
    def test_returns_model3d(self):
        from pycsamt.models.modem.model3d import ModEmModel3D
        src = _uniform_model()
        tgt = _uniform_model(nx=4, ny=3, nz=3, n_air=1)
        result = interp_model3d(src, tgt, bg_rho=100.0)
        assert isinstance(result, ModEmModel3D)

    def test_output_shape_matches_target(self):
        src = _uniform_model(nx=6, ny=5, nz=4)
        tgt = _uniform_model(nx=3, ny=4, nz=5)
        result = interp_model3d(src, tgt)
        assert result.rho_loge.shape == (tgt.nz, tgt.ny, tgt.nx)

    def test_uniform_model_roundtrip(self):
        """Interpolating a uniform model preserves the value."""
        src = _uniform_model(rho=200.0)
        tgt = _uniform_model(nx=4, ny=4, nz=4, rho=1.0)  # different initial rho
        result = interp_model3d(src, tgt, bg_rho=200.0)
        # Interior points should be ≈ ln(200)
        np.testing.assert_allclose(
            result.rho_loge, np.log(200.0), atol=0.05
        )

    def test_no_nans_in_result(self):
        src = _uniform_model()
        tgt = _uniform_model(nx=4, ny=4, nz=4)
        result = interp_model3d(src, tgt, bg_rho=100.0)
        assert not np.any(np.isnan(result.rho_loge))

    def test_background_fills_outside_domain(self):
        """Target grid extends beyond source → outside cells get bg_rho."""
        from pycsamt.models.modem.model3d import ModEmModel3D
        src = _uniform_model(nx=4, ny=4, nz=4, rho=100.0)
        # Target is 10× wider → outer cells get background
        tgt = ModEmModel3D()
        tgt.x_widths = np.full(6, 5000.0)   # much wider than source
        tgt.y_widths = np.full(6, 5000.0)
        tgt.z_widths = np.full(4, 500.0)
        tgt.n_air = 1
        tgt.rho_loge = np.full((4, 6, 6), np.log(1.0))
        bg = 500.0
        result = interp_model3d(src, tgt, bg_rho=bg)
        assert not np.any(np.isnan(result.rho_loge))
        # Some cells at the extreme boundary should equal bg
        assert np.any(np.isclose(result.rho_loge, np.log(bg), atol=0.01))

    def test_scalar_bg_rho(self):
        src = _uniform_model()
        tgt = _uniform_model(nx=3, ny=3, nz=3)
        result = interp_model3d(src, tgt, bg_rho=50.0)
        assert result.rho_loge.shape == (3, 3, 3)

    def test_array_bg_rho(self):
        src = _uniform_model(nz=4)
        tgt = _uniform_model(nx=3, ny=3, nz=4)
        bg = np.array([10.0, 100.0, 1000.0, 10000.0])
        result = interp_model3d(src, tgt, bg_rho=bg)
        assert result.rho_loge.shape == (4, 3, 3)

    def test_target_grid_preserved(self):
        src = _gradient_model()
        tgt = _uniform_model(nx=3, ny=3, nz=4)
        result = interp_model3d(src, tgt)
        np.testing.assert_array_equal(result.x_widths, tgt.x_widths)
        np.testing.assert_array_equal(result.y_widths, tgt.y_widths)
        np.testing.assert_array_equal(result.z_widths, tgt.z_widths)

    def test_finer_target_more_cells(self):
        src = _uniform_model(nx=4, ny=4, nz=3)
        tgt = _uniform_model(nx=8, ny=8, nz=6)
        result = interp_model3d(src, tgt, bg_rho=100.0)
        assert result.nx == 8 and result.ny == 8 and result.nz == 6


# ===========================================================================
# interp_z3d
# ===========================================================================

class TestInterpZ3D:
    def test_returns_impedancefile(self):
        imp = _z3d_imp_on_grid()
        new_loc = np.column_stack([
            np.linspace(-3000, 3000, 5),
            np.linspace(-3000, 3000, 5),
            np.zeros(5),
        ])
        result = interp_z3d(imp, new_loc)
        assert isinstance(result, ImpedanceFile)

    def test_n_sites_updated(self):
        imp = _z3d_imp_on_grid(n_pts_per_dim=4)
        new_loc = np.column_stack([
            np.linspace(-2000, 2000, 7),
            np.linspace(-2000, 2000, 7),
            np.zeros(7),
        ])
        result = interp_z3d(imp, new_loc)
        assert result.blocks[0].n_sites == 7

    def test_site_loc_updated(self):
        imp = _z3d_imp_on_grid()
        new_loc = np.column_stack([np.arange(5, dtype=float) * 500,
                                   np.zeros(5), np.zeros(5)])
        result = interp_z3d(imp, new_loc)
        np.testing.assert_array_equal(result.blocks[0].site_loc, new_loc)

    def test_site_names_auto(self):
        imp = _z3d_imp_on_grid()
        new_loc = np.zeros((3, 3))
        result = interp_z3d(imp, new_loc)
        assert result.blocks[0].site_names == ["S000", "S001", "S002"]

    def test_site_names_explicit(self):
        imp = _z3d_imp_on_grid()
        new_loc = np.zeros((2, 3))
        result = interp_z3d(imp, new_loc, new_site_names=["AA", "BB"])
        assert result.blocks[0].site_names == ["AA", "BB"]

    def test_n_periods_preserved(self):
        imp = _z3d_imp_on_grid(periods=(1.0, 10.0, 100.0))
        new_loc = np.zeros((3, 3))
        result = interp_z3d(imp, new_loc)
        assert result.n_periods == 3

    def test_comp_names_preserved(self):
        imp = _z3d_imp_on_grid()
        new_loc = np.zeros((2, 3))
        result = interp_z3d(imp, new_loc)
        assert result.blocks[0].comp_names == ["ZXX", "ZXY", "ZYX", "ZYY"]

    def test_pct_error_applied(self):
        imp = _z3d_imp_on_grid(n_pts_per_dim=4)
        new_loc = np.column_stack([
            np.linspace(-3000, 3000, 5),
            np.linspace(-3000, 3000, 5),
            np.zeros(5),
        ])
        result = interp_z3d(imp, new_loc, pct_error=5.0)
        # Zerr should equal 5% of |Z|
        for blk in result.blocks:
            expected = np.abs(blk.Z) * 0.05
            np.testing.assert_allclose(blk.Zerr, expected, rtol=1e-10)

    def test_interpolation_at_original_site(self):
        """Querying at original locations should reproduce source values."""
        imp = _z3d_imp_on_grid(n_pts_per_dim=5, periods=(10.0,))
        blk = imp.blocks[0]
        # Query at the first site's location
        target_loc = blk.site_loc[:1, :]   # just one site
        result = interp_z3d(imp, target_loc)
        # Interpolated value should be close to original for point exactly on the grid
        np.testing.assert_allclose(
            result.blocks[0].Z[0].real, blk.Z[0].real, rtol=0.01
        )

    def test_z_shape_correct(self):
        imp = _z3d_imp_on_grid(n_pts_per_dim=4)
        new_loc = np.column_stack([
            np.linspace(-2000, 2000, 6),
            np.linspace(-2000, 2000, 6),
            np.zeros(6),
        ])
        result = interp_z3d(imp, new_loc)
        assert result.blocks[0].Z.shape == (6, 4)   # 6 sites, 4 components
