"""Self-contained tests for pycsamt.models.modem.model2d.ModEmModel2D.

Unlike test_modem_phase2.py, these tests build their own ModEM 2-D
model files and ``ModEmData`` objects in-memory, so they do not depend
on the (locally absent) ModEMv626 example data tree.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.models.modem.config import ModEmConfig
from pycsamt.models.modem.data import ModEmData
from pycsamt.models.modem.model2d import ModEmModel2D


def _write_model_file(path, nx=3, nz=2, log_type="LOGE", rho_rows=None):
    lines = [f"  {nx}  {nz}  {log_type}\n"]
    lines.append("  " + "  ".join(f"{100.0:.4f}" for _ in range(nx)) + "\n")
    lines.append("  " + "  ".join(f"{50.0:.4f}" for _ in range(nz)) + "\n")
    lines.append("  1\n")
    if rho_rows is None:
        rho_rows = [[4.6 for _ in range(nx)] for _ in range(nz)]
    for row in rho_rows:
        lines.append("  " + "  ".join(f"{v:.5E}" for v in row) + "\n")
    path.write_text("".join(lines))
    return path


class TestParseLogTypes:
    def test_read_loge(self, tmp_path):
        p = _write_model_file(tmp_path / "m.rho", log_type="LOGE")
        m = ModEmModel2D.read(p)
        assert m.log_type == "LOGE"
        np.testing.assert_allclose(m.rho_loge, 4.6)

    def test_read_log10_converted(self, tmp_path):
        rho_rows = [[2.0, 2.0, 2.0], [2.0, 2.0, 2.0]]
        p = _write_model_file(
            tmp_path / "m.rho", log_type="LOG10", rho_rows=rho_rows
        )
        m = ModEmModel2D.read(p)
        assert m.log_type == "LOG10"
        np.testing.assert_allclose(m.rho_loge, 2.0 * np.log(10.0))

    def test_read_linear_converted(self, tmp_path):
        rho_rows = [[100.0, 100.0, 100.0], [100.0, 100.0, 100.0]]
        p = _write_model_file(
            tmp_path / "m.rho", log_type="LINEAR", rho_rows=rho_rows
        )
        m = ModEmModel2D.read(p)
        assert m.log_type == "LINEAR"
        np.testing.assert_allclose(m.rho_loge, np.log(100.0))

    def test_read_linear_nonpositive_becomes_nan(self, tmp_path):
        rho_rows = [[0.0, -1.0, 100.0], [100.0, 100.0, 100.0]]
        p = _write_model_file(
            tmp_path / "m.rho", log_type="LINEAR", rho_rows=rho_rows
        )
        m = ModEmModel2D.read(p)
        assert np.isnan(m.rho_loge[0, 0])
        assert np.isnan(m.rho_loge[0, 1])
        assert np.isfinite(m.rho_loge[0, 2])

    def test_read_default_log_type_when_missing(self, tmp_path):
        p = tmp_path / "m.rho"
        p.write_text(
            "  2  2\n"
            "  100.0  100.0\n"
            "  50.0  50.0\n"
            "  1\n"
            "  4.6  4.6\n"
            "  4.6  4.6\n"
        )
        m = ModEmModel2D.read(p)
        assert m.log_type == "LOGE"

    def test_comments_and_blank_lines_skipped(self, tmp_path):
        p = tmp_path / "m.rho"
        p.write_text(
            "# comment header\n"
            "\n"
            "  3  2  LOGE\n"
            "  100.0  100.0  100.0\n"
            "  50.0  50.0\n"
            "# block count\n"
            "  1\n"
            "  4.6  4.6  4.6\n"
            "  4.6  4.6  4.6\n"
        )
        m = ModEmModel2D.read(p)
        assert m.nx == 3
        assert m.nz == 2


class TestModEmModel2DProperties:
    @pytest.fixture
    def model(self, tmp_path):
        p = _write_model_file(tmp_path / "m.rho", nx=4, nz=3)
        return ModEmModel2D.read(p)

    def test_nx_nz(self, model):
        assert model.nx == 4
        assert model.nz == 3

    def test_x_nodes_z_nodes(self, model):
        assert len(model.x_nodes) == model.nx + 1
        assert len(model.z_nodes) == model.nz + 1
        assert model.x_nodes[0] == 0.0
        assert model.x_nodes[-1] == pytest.approx(400.0)
        assert model.z_nodes[-1] == pytest.approx(150.0)

    def test_rho_linear(self, model):
        np.testing.assert_allclose(model.rho_linear, np.exp(model.rho_loge))

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            ModEmModel2D.read("/no/such/model2d.rho")

    def test_empty_model_defaults(self):
        m = ModEmModel2D()
        assert m.nx == 0
        assert m.nz == 0
        assert m.rho_loge.shape == (0, 0)
        assert isinstance(m.config, ModEmConfig)


class TestWriteRoundtrip:
    def test_roundtrip_preserves_grid(self, tmp_path):
        p = _write_model_file(tmp_path / "m.rho", nx=5, nz=4)
        m = ModEmModel2D.read(p)
        out = tmp_path / "sub" / "m_out.rho"
        result_path = m.write(out)
        assert result_path == out
        assert out.exists()

        back = ModEmModel2D.read(out)
        assert back.nx == m.nx
        assert back.nz == m.nz
        np.testing.assert_allclose(back.rho_loge, m.rho_loge, rtol=1e-4)

    def test_write_creates_parent_dirs(self, tmp_path):
        p = _write_model_file(tmp_path / "m.rho")
        m = ModEmModel2D.read(p)
        out = tmp_path / "a" / "b" / "c" / "out.rho"
        m.write(out)
        assert out.exists()

    def test_write_many_columns_wraps_rows(self, tmp_path):
        p = _write_model_file(tmp_path / "m.rho", nx=15, nz=2)
        m = ModEmModel2D.read(p)
        out = tmp_path / "wide.rho"
        m.write(out)
        back = ModEmModel2D.read(out)
        assert back.nx == 15
        np.testing.assert_allclose(back.rho_loge, m.rho_loge, rtol=1e-4)


def _make_data_from_offsets(offsets, mode="2d"):
    """Build a populated ModEmData with the given station easting offsets."""

    class _Site:
        pass

    n_freq = 4
    freqs = np.logspace(1, -1, n_freq)
    sites = []
    for i, off in enumerate(offsets):
        s = _Site()
        s.name = f"S{i:02d}"
        s.coords = (0.0, float(off), 0.0)
        omega = 2 * np.pi * freqs
        mu0 = 4 * np.pi * 1e-7
        z_mag = np.sqrt(omega * mu0 * 100.0)
        z_val = z_mag * (1.0 + 1.0j) / np.sqrt(2)
        z_arr = np.zeros((n_freq, 2, 2), dtype=complex)
        z_arr[:, 0, 1] = z_val
        z_arr[:, 1, 0] = -z_val
        s.freq = freqs
        s.z = z_arr
        s.z_err = np.abs(z_arr) * 0.05
        sites.append(s)
    cfg = ModEmConfig(mode=mode)
    return ModEmData.from_edi(sites, config=cfg), cfg


class TestHalfspace:
    def test_halfspace_basic_grid(self):
        data, cfg = _make_data_from_offsets(
            [0.0, 500.0, 1000.0, 1500.0, 2000.0]
        )
        m = ModEmModel2D.halfspace(data, config=cfg)
        assert m.nx > 0
        assert m.nz > 0
        assert m.rho_loge.shape == (m.nz, m.nx)

    def test_halfspace_earth_cells_match_initial_rho(self):
        data, cfg = _make_data_from_offsets([0.0, 300.0, 900.0])
        m = ModEmModel2D.halfspace(data, config=cfg)
        air = cfg.n_airlayers_2d
        np.testing.assert_allclose(
            m.rho_loge[air:, :], np.log(cfg.initial_rho), rtol=1e-6
        )

    def test_halfspace_air_layers_high_resistivity(self):
        cfg = ModEmConfig(mode="2d", n_airlayers_2d=3)
        data, _ = _make_data_from_offsets([0.0, 400.0, 800.0], mode="2d")
        m = ModEmModel2D.halfspace(data, config=cfg)
        np.testing.assert_allclose(m.rho_loge[:3, :], np.log(1e12))

    def test_halfspace_no_air_layers(self):
        cfg = ModEmConfig(mode="2d", n_airlayers_2d=0)
        data, _ = _make_data_from_offsets([0.0, 400.0], mode="2d")
        m = ModEmModel2D.halfspace(data, config=cfg)
        np.testing.assert_allclose(m.rho_loge, np.log(cfg.initial_rho))

    def test_halfspace_single_station(self):
        cfg = ModEmConfig(mode="2d", cell_size_h_2d=50.0)
        data, _ = _make_data_from_offsets([100.0], mode="2d")
        m = ModEmModel2D.halfspace(data, config=cfg)
        assert m.nx > 0

    def test_halfspace_odd_gap_count_padded_even(self):
        cfg = ModEmConfig(mode="2d", cell_size_h_2d=100.0, n_padding_x_2d=2)
        # 4 stations -> 3 gaps; whatever the per-gap cell count, exercise
        # the "append one more station-zone width" odd-length branch.
        data, _ = _make_data_from_offsets(
            [0.0, 250.0, 500.0, 750.0], mode="2d"
        )
        m = ModEmModel2D.halfspace(data, config=cfg)
        assert m.nx > 0

    def test_halfspace_default_config(self):
        data, _ = _make_data_from_offsets([0.0, 500.0])
        m = ModEmModel2D.halfspace(data)
        assert isinstance(m.config, ModEmConfig)

    def test_halfspace_verbose_logs(self, caplog):
        data, cfg = _make_data_from_offsets([0.0, 500.0])
        m = ModEmModel2D.halfspace(data, config=cfg, verbose=1)
        assert m.verbose == 1

    def test_read_verbose_logs(self, tmp_path):
        p = _write_model_file(tmp_path / "m.rho")
        m = ModEmModel2D.read(p, verbose=1)
        assert m.verbose == 1
