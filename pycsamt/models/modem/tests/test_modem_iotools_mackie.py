"""Tests for pycsamt.models.modem.iotools.mackie."""

import math
from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Optional real-file fixtures (skipped when ModEMv626/ is absent)
# ---------------------------------------------------------------------------

_EX2D = (
    Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "2D_MT" / "BLOCK2"
)
_SKIP_REAL = pytest.mark.skipif(
    not _EX2D.exists(), reason="ModEMv626 example data not available"
)


# ---------------------------------------------------------------------------
# Helpers to build tiny synthetic models
# ---------------------------------------------------------------------------


def _make_model2d(nx=5, nz_earth=3, rho_val=100.0, n_air=0):
    """Return a ModEmModel2D with uniform resistivity (no air layers by default)."""
    from pycsamt.models.modem.model2d import ModEmModel2D

    obj = ModEmModel2D()
    obj.x_widths = np.full(nx, 1000.0)
    z_earth = np.geomspace(50.0, 500.0, nz_earth)
    obj.z_widths = z_earth
    obj.rho_loge = np.full((nz_earth, nx), math.log(rho_val))
    obj.n_air = 0
    return obj


def _make_model3d(nx=4, ny=3, nz=5, n_air=2, rho_val=100.0):
    """Return a ModEmModel3D with uniform resistivity."""
    from pycsamt.models.modem.model3d import ModEmModel3D

    obj = ModEmModel3D()
    obj.x_widths = np.full(nx, 2000.0)
    obj.y_widths = np.full(ny, 2000.0)
    obj.z_widths = np.concatenate(
        [
            np.full(n_air, 50.0),
            np.geomspace(100.0, 500.0, nz - n_air),
        ]
    )
    rho_grid = np.full((nz, ny, nx), math.log(rho_val))
    if n_air > 0:
        rho_grid[:n_air] = math.log(1e12)
    obj.rho_loge = rho_grid
    obj.n_air = n_air
    return obj


# ===========================================================================
# 2D — write / read roundtrip
# ===========================================================================


class TestMackie2DRoundtrip:
    def test_write_creates_file(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            write_mackie2d,
        )

        m = _make_model2d()
        out = write_mackie2d(m, tmp_path / "m.rho")
        assert out.exists()

    def test_header_line(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            write_mackie2d,
        )

        m = _make_model2d(nx=5, nz_earth=3)
        out = write_mackie2d(m, tmp_path / "m.rho")
        first = out.read_text().splitlines()[0]
        assert "5" in first and "3" in first

    def test_roundtrip_nx(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie2d,
            write_mackie2d,
        )

        m = _make_model2d(nx=6, nz_earth=4)
        write_mackie2d(m, tmp_path / "m.rho")
        back = read_mackie2d(tmp_path / "m.rho")
        assert back.nx == 6

    def test_roundtrip_nz_earth(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            _N_AIR_2D,
            read_mackie2d,
            write_mackie2d,
        )

        m = _make_model2d(nx=5, nz_earth=3)
        write_mackie2d(m, tmp_path / "m.rho")
        back = read_mackie2d(tmp_path / "m.rho")
        # read_mackie2d prepends 10 air layers
        assert back.nz == 3 + _N_AIR_2D

    def test_roundtrip_x_widths(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie2d,
            write_mackie2d,
        )

        m = _make_model2d(nx=5)
        write_mackie2d(m, tmp_path / "m.rho")
        back = read_mackie2d(tmp_path / "m.rho")
        np.testing.assert_allclose(back.x_widths, m.x_widths, rtol=1e-4)

    def test_roundtrip_z_earth_widths(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            _N_AIR_2D,
            read_mackie2d,
            write_mackie2d,
        )

        m = _make_model2d(nx=5, nz_earth=3)
        write_mackie2d(m, tmp_path / "m.rho")
        back = read_mackie2d(tmp_path / "m.rho")
        np.testing.assert_allclose(back.z_widths[_N_AIR_2D:], m.z_widths, rtol=1e-4)

    def test_roundtrip_rho_values(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            _N_AIR_2D,
            read_mackie2d,
            write_mackie2d,
        )

        m = _make_model2d(nx=5, nz_earth=3, rho_val=50.0)
        write_mackie2d(m, tmp_path / "m.rho")
        back = read_mackie2d(tmp_path / "m.rho")
        np.testing.assert_allclose(back.rho_loge[_N_AIR_2D:], m.rho_loge, rtol=1e-5)

    def test_air_layers_have_high_rho(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            _N_AIR_2D,
            read_mackie2d,
            write_mackie2d,
        )

        m = _make_model2d()
        write_mackie2d(m, tmp_path / "m.rho")
        back = read_mackie2d(tmp_path / "m.rho")
        assert np.all(back.rho_loge[:_N_AIR_2D] > 20.0)

    def test_sentinel_present(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            write_mackie2d,
        )

        m = _make_model2d()
        out = write_mackie2d(m, tmp_path / "m.rho")
        text = out.read_text()
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        # sentinel '1' appears as its own line somewhere after widths
        assert "1" in lines

    def test_loge_keyword_in_file(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            write_mackie2d,
        )

        m = _make_model2d()
        out = write_mackie2d(m, tmp_path / "m.rho", log_type="LOGE")
        assert "LOGE" in out.read_text().splitlines()[0]

    def test_n_air_explicit(self, tmp_path):
        """Explicit n_air overrides auto-detection."""
        from pycsamt.models.modem.iotools.mackie import (
            _N_AIR_2D,
            read_mackie2d,
            write_mackie2d,
        )

        m = _make_model2d(nx=5, nz_earth=4)
        # model has n_air=0, so n_air=0 means write all layers as earth
        write_mackie2d(m, tmp_path / "m.rho", n_air=0)
        back = read_mackie2d(tmp_path / "m.rho")
        assert back.nz == 4 + _N_AIR_2D


# ===========================================================================
# 2D — real file tests
# ===========================================================================


class TestMackie2DRealFile:
    @_SKIP_REAL
    def test_read_returns_model2d(self):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie2d,
        )
        from pycsamt.models.modem.model2d import ModEmModel2D

        m = read_mackie2d(_EX2D / "m0.rho")
        assert isinstance(m, ModEmModel2D)

    @_SKIP_REAL
    def test_real_nx(self):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie2d,
        )

        m = read_mackie2d(_EX2D / "m0.rho")
        assert m.nx == 100

    @_SKIP_REAL
    def test_real_nz(self):
        from pycsamt.models.modem.iotools.mackie import (
            _N_AIR_2D,
            read_mackie2d,
        )

        m = read_mackie2d(_EX2D / "m0.rho")
        # file has 31 earth layers + 10 air = 41 total
        assert m.nz == 31 + _N_AIR_2D

    @_SKIP_REAL
    def test_real_rho_loge_finite(self):
        from pycsamt.models.modem.iotools.mackie import (
            _N_AIR_2D,
            read_mackie2d,
        )

        m = read_mackie2d(_EX2D / "m0.rho")
        earth = m.rho_loge[_N_AIR_2D:]
        assert np.all(np.isfinite(earth))

    @_SKIP_REAL
    def test_real_roundtrip_rho(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            _N_AIR_2D,
            read_mackie2d,
            write_mackie2d,
        )

        m = read_mackie2d(_EX2D / "m0.rho")
        out = write_mackie2d(m, tmp_path / "m0_out.rho")
        back = read_mackie2d(out)
        np.testing.assert_allclose(
            back.rho_loge[_N_AIR_2D:], m.rho_loge[_N_AIR_2D:], rtol=1e-5
        )


# ===========================================================================
# 3D — write / read roundtrip
# ===========================================================================


class TestMackie3DRoundtrip:
    def test_write_creates_file(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            write_mackie3d,
        )

        m = _make_model3d()
        out = write_mackie3d(m, tmp_path / "m.mac")
        assert out.exists()

    def test_header_values(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            write_mackie3d,
        )

        m = _make_model3d(nx=4, ny=3, nz=5, n_air=2)
        out = write_mackie3d(m, tmp_path / "m.mac")
        first = out.read_text().splitlines()[0]
        assert "4" in first and "3" in first and "5" in first and "2" in first

    def test_values_keyword(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            write_mackie3d,
        )

        m = _make_model3d()
        out = write_mackie3d(m, tmp_path / "m.mac")
        assert "VALUES" in out.read_text().splitlines()[0]

    def test_roundtrip_shape(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
            write_mackie3d,
        )

        m = _make_model3d(nx=4, ny=3, nz=5, n_air=2)
        write_mackie3d(m, tmp_path / "m.mac")
        back = read_mackie3d(tmp_path / "m.mac")
        assert back.shape == m.shape

    def test_roundtrip_n_air(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
            write_mackie3d,
        )

        m = _make_model3d(n_air=2)
        write_mackie3d(m, tmp_path / "m.mac")
        back = read_mackie3d(tmp_path / "m.mac")
        assert back.n_air == 2

    def test_roundtrip_x_widths(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
            write_mackie3d,
        )

        m = _make_model3d(nx=4)
        write_mackie3d(m, tmp_path / "m.mac")
        back = read_mackie3d(tmp_path / "m.mac")
        np.testing.assert_allclose(back.x_widths, m.x_widths, rtol=1e-5)

    def test_roundtrip_y_widths(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
            write_mackie3d,
        )

        m = _make_model3d(ny=3)
        write_mackie3d(m, tmp_path / "m.mac")
        back = read_mackie3d(tmp_path / "m.mac")
        np.testing.assert_allclose(back.y_widths, m.y_widths, rtol=1e-5)

    def test_roundtrip_z_widths(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
            write_mackie3d,
        )

        m = _make_model3d(nz=5, n_air=2)
        write_mackie3d(m, tmp_path / "m.mac")
        back = read_mackie3d(tmp_path / "m.mac")
        np.testing.assert_allclose(back.z_widths, m.z_widths, rtol=1e-5)

    def test_roundtrip_rho_earth(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
            write_mackie3d,
        )

        m = _make_model3d(nx=4, ny=3, nz=5, n_air=2, rho_val=200.0)
        write_mackie3d(m, tmp_path / "m.mac")
        back = read_mackie3d(tmp_path / "m.mac")
        np.testing.assert_allclose(back.rho_loge[2:], m.rho_loge[2:], rtol=1e-5)

    def test_roundtrip_rho_air(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
            write_mackie3d,
        )

        m = _make_model3d(n_air=2)
        write_mackie3d(m, tmp_path / "m.mac")
        back = read_mackie3d(tmp_path / "m.mac")
        np.testing.assert_allclose(back.rho_loge[:2], m.rho_loge[:2], rtol=1e-5)

    def test_winglink_footer_present(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            write_mackie3d,
        )

        m = _make_model3d()
        out = write_mackie3d(m, tmp_path / "m.mac")
        assert "WINGLINK" in out.read_text()

    def test_origin_written(self, tmp_path):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
            write_mackie3d,
        )

        m = _make_model3d()
        write_mackie3d(m, tmp_path / "m.mac", origin=(1000.0, 2000.0, 0.0))
        back = read_mackie3d(tmp_path / "m.mac")
        np.testing.assert_allclose(back.origin[:2], [1000.0, 2000.0], atol=1.0)

    def test_non_uniform_rho_roundtrip(self, tmp_path):
        """Heterogeneous rho grid survives write–read cycle."""
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
            write_mackie3d,
        )
        from pycsamt.models.modem.model3d import ModEmModel3D

        rng = np.random.default_rng(42)
        nx, ny, nz = 3, 3, 4
        m = ModEmModel3D()
        m.x_widths = np.full(nx, 1000.0)
        m.y_widths = np.full(ny, 1000.0)
        m.z_widths = np.full(nz, 200.0)
        m.n_air = 0
        m.rho_loge = rng.uniform(3.0, 8.0, (nz, ny, nx))
        write_mackie3d(m, tmp_path / "m.mac")
        back = read_mackie3d(tmp_path / "m.mac")
        np.testing.assert_allclose(back.rho_loge, m.rho_loge, rtol=1e-5)

    def test_missing_file_raises(self):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie3d,
        )

        with pytest.raises(FileNotFoundError):
            read_mackie3d("/nonexistent/path/model.mac")

    def test_missing_file_2d_raises(self):
        from pycsamt.models.modem.iotools.mackie import (
            read_mackie2d,
        )

        with pytest.raises(FileNotFoundError):
            read_mackie2d("/nonexistent/path/model.rho")
