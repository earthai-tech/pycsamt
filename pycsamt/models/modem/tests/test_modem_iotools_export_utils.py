"""Tests for pycsamt.models.modem.iotools.export and .utils."""

import math

import numpy as np
import pytest

from pycsamt.models.modem.iotools.export import (
    write_meshtools3d,
)
from pycsamt.models.modem.iotools.utils import (
    imp_units_factor,
    linear_to_loge,
    log10_to_loge,
    loge_to_linear,
    loge_to_log10,
    skin_depth,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _model(nx=5, ny=4, nz=6, n_air=2, rho=100.0):
    from pycsamt.models.modem.model3d import ModEmModel3D
    m = ModEmModel3D()
    m.x_widths = np.full(nx, 2000.0)
    m.y_widths = np.full(ny, 2000.0)
    m.z_widths = np.concatenate([np.full(n_air, 50.0),
                                  np.full(nz - n_air, 500.0)])
    m.n_air    = n_air
    m.rho_loge = np.full((nz, ny, nx), math.log(rho))
    return m


# ===========================================================================
# export.write_meshtools3d
# ===========================================================================

class TestWriteMeshtools3D:
    def test_returns_two_paths(self, tmp_path):
        m = _model()
        msh, con = write_meshtools3d(m, tmp_path / "model")
        assert msh.exists() and con.exists()

    def test_msh_extension(self, tmp_path):
        m = _model()
        msh, _ = write_meshtools3d(m, tmp_path / "model.something")
        assert msh.suffix == ".msh"

    def test_con_extension(self, tmp_path):
        m = _model()
        _, con = write_meshtools3d(m, tmp_path / "model")
        assert con.suffix == ".con"

    def test_msh_first_line_dims(self, tmp_path):
        m = _model(nx=5, ny=4, nz=6, n_air=2)
        msh, _ = write_meshtools3d(m, tmp_path / "model")
        first = msh.read_text().splitlines()[0].split()
        nx, ny, nz_earth = int(first[0]), int(first[1]), int(first[2])
        assert nx == 5 and ny == 4 and nz_earth == 4   # 6 - 2 air

    def test_msh_origin_line(self, tmp_path):
        m = _model()
        msh, _ = write_meshtools3d(m, tmp_path / "model")
        lines = msh.read_text().splitlines()
        # second line = origin (should be 0 0 0 for default)
        origin = [float(v) for v in lines[1].split()]
        assert len(origin) == 3

    def test_con_line_count(self, tmp_path):
        m = _model(nx=5, ny=4, nz=6, n_air=2)
        _, con = write_meshtools3d(m, tmp_path / "model")
        lines = [ln for ln in con.read_text().splitlines() if ln.strip()]
        # ny * nx * nz_earth = 4 * 5 * 4 = 80
        assert len(lines) == 80

    def test_con_values_positive(self, tmp_path):
        m = _model(rho=100.0)
        _, con = write_meshtools3d(m, tmp_path / "model")
        vals = [float(ln) for ln in con.read_text().splitlines() if ln.strip()]
        assert all(v > 0 for v in vals)

    def test_con_uniform_rho(self, tmp_path):
        """Uniform 100 Ω·m model → σ = 0.01 S/m in every cell."""
        m = _model(rho=100.0)
        _, con = write_meshtools3d(m, tmp_path / "model")
        vals = np.array([float(ln) for ln in con.read_text().splitlines() if ln.strip()])
        np.testing.assert_allclose(vals, 0.01, rtol=1e-5)

    def test_air_layers_excluded(self, tmp_path):
        """Air cells (very high rho) should not appear in the .con file."""
        m = _model(nz=6, n_air=2, rho=100.0)
        m.rho_loge[:2] = math.log(1e12)   # set air cells explicitly
        _, con = write_meshtools3d(m, tmp_path / "model")
        vals = np.array([float(ln) for ln in con.read_text().splitlines() if ln.strip()])
        # All values should correspond to ρ=100, not 1e12
        np.testing.assert_allclose(vals, 0.01, rtol=1e-5)

    def test_msh_contains_x_widths(self, tmp_path):
        m = _model(nx=3)
        msh, _ = write_meshtools3d(m, tmp_path / "model")
        text = msh.read_text()
        assert "2000" in text   # x-widths should appear


# ===========================================================================
# utils.skin_depth
# ===========================================================================

class TestSkinDepth:
    def test_scalar(self):
        # δ = √(2·100 / (2π/1 · 4π×10⁻⁷)) ≈ 5033 m
        d = skin_depth(1.0, rho=100.0)
        assert abs(d - 5032.9) < 1.0

    def test_array_periods(self):
        d = skin_depth([1.0, 100.0], rho=100.0)
        assert d.shape == (2,)
        # longer period → deeper skin depth
        assert d[1] > d[0]

    def test_scales_with_sqrt_rho(self):
        d1 = skin_depth(1.0, rho=100.0)
        d4 = skin_depth(1.0, rho=400.0)
        np.testing.assert_allclose(d4 / d1, 2.0, rtol=1e-10)

    def test_scales_with_sqrt_period(self):
        d1 = skin_depth(1.0,  rho=100.0)
        d4 = skin_depth(4.0,  rho=100.0)
        np.testing.assert_allclose(d4 / d1, 2.0, rtol=1e-10)

    def test_known_value(self):
        # δ ≈ √(2·T·ρ / (4π²·10⁻⁷)) ≈ 503 √(T·ρ) metres
        T, rho = 10.0, 100.0
        expected = np.sqrt(2.0 * rho / (2.0 * np.pi / T * 4.0 * np.pi * 1e-7))
        np.testing.assert_allclose(skin_depth(T, rho), expected, rtol=1e-10)


# ===========================================================================
# utils.imp_units_factor
# ===========================================================================

class TestImpUnitsFactor:
    def test_same_units(self):
        assert imp_units_factor("[V/m]/[T]", "[V/m]/[T]") == 1.0

    def test_si_to_ohm(self):
        assert imp_units_factor("[V/m]/[T]", "Ohm") == 1.0

    def test_unknown_raises(self):
        with pytest.raises(ValueError):
            imp_units_factor("furlongs", "[V/m]/[T]")


# ===========================================================================
# utils encoding conversions
# ===========================================================================

class TestEncodingConversions:
    def test_loge_to_log10(self):
        loge = np.array([np.log(100.0)])
        log10 = loge_to_log10(loge)
        np.testing.assert_allclose(log10, [2.0], rtol=1e-10)

    def test_log10_to_loge(self):
        log10 = np.array([2.0])
        loge  = log10_to_loge(log10)
        np.testing.assert_allclose(loge, [np.log(100.0)], rtol=1e-10)

    def test_loge_to_linear(self):
        loge = np.array([np.log(100.0)])
        np.testing.assert_allclose(loge_to_linear(loge), [100.0], rtol=1e-10)

    def test_linear_to_loge(self):
        rho = np.array([100.0])
        np.testing.assert_allclose(linear_to_loge(rho), [np.log(100.0)], rtol=1e-10)

    def test_roundtrip_loge_log10(self):
        rho = np.array([1.0, 10.0, 100.0, 1000.0])
        loge = linear_to_loge(rho)
        np.testing.assert_allclose(loge_to_linear(loge), rho, rtol=1e-10)

    def test_linear_to_loge_nonpositive_safe(self):
        result = linear_to_loge(np.array([0.0, -1.0]))
        assert np.all(np.isfinite(result))
