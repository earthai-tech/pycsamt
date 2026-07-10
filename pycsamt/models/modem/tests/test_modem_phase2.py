"""Phase 2 tests — ModEmData, ModEmModel2D, ModEmControl, ModEmLog."""

from pathlib import Path

import numpy as np
import pytest

_EX2D = Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "2D_MT" / "BLOCK2"
_EX3D = Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "3D_MT" / "BLOCK2"

pytestmark = pytest.mark.skipif(
    not _EX2D.exists(), reason="ModEMv626 example data not available"
)


# ---------------------------------------------------------------------------
# ModEmData — read
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def data_2d():
    from pycsamt.models.modem.data import ModEmData
    return ModEmData.read(_EX2D / "d0_JT.dat")


@pytest.fixture(scope="module")
def data_3d():
    from pycsamt.models.modem.data import ModEmData
    return ModEmData.read(_EX3D / "d0.dat")


def test_data_2d_n_sites(data_2d):
    assert data_2d.n_sites > 0


def test_data_2d_n_periods(data_2d):
    assert data_2d.n_periods > 0


def test_data_2d_periods_descending(data_2d):
    p = data_2d.periods
    assert np.all(p[:-1] >= p[1:])


def test_data_2d_blocks_have_rows(data_2d):
    assert len(data_2d.blocks) >= 1
    assert len(data_2d.blocks[0]["rows"]) > 0


def test_data_2d_component_type(data_2d):
    assert "TE_Impedance" in data_2d.component_types


def test_data_3d_n_sites(data_3d):
    assert data_3d.n_sites > 0


def test_data_3d_full_impedance(data_3d):
    assert "Full_Impedance" in data_3d.component_types


def test_data_roundtrip_2d(data_2d, tmp_path):
    out = tmp_path / "d0_out.dat"
    data_2d.write(out)
    from pycsamt.models.modem.data import ModEmData
    back = ModEmData.read(out)
    assert back.n_sites   == data_2d.n_sites
    assert back.n_periods == data_2d.n_periods


def test_data_missing_file_raises():
    from pycsamt.models.modem.data import ModEmData
    with pytest.raises(FileNotFoundError):
        ModEmData.read("/no/such/file.dat")


# ---------------------------------------------------------------------------
# ModEmData — from_edi with synthetic duck-typed sites
# ---------------------------------------------------------------------------

def _make_site(name, y_offset, n_freq=6):
    """Synthetic half-space site: Z_TE = Z_TM = (1+1j) * |Z|."""
    freqs = np.logspace(2, -2, n_freq)   # 100 Hz → 0.01 Hz
    rho   = 100.0
    omega = 2 * np.pi * freqs
    mu0   = 4 * np.pi * 1e-7
    z_mag = np.sqrt(omega * mu0 * rho)
    z_val = z_mag * (1.0 + 1.0j) / np.sqrt(2)

    z_arr = np.zeros((n_freq, 2, 2), dtype=complex)
    z_arr[:, 0, 1] =  z_val   # ZXY = Z_TE
    z_arr[:, 1, 0] = -z_val   # ZYX = Z_TM

    class _Site:
        pass
    s = _Site()
    s.name   = name
    s.coords = (0.0, y_offset, 0.0)
    s.freq   = freqs
    s.z      = z_arr
    s.z_err  = np.abs(z_arr) * 0.05
    return s


def test_from_edi_n_sites():
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.data import ModEmData
    sites = [_make_site(f"S{i:02d}", i * 500.0) for i in range(5)]
    cfg   = ModEmConfig(mode="3d", component_type="Full_Impedance")
    d     = ModEmData.from_edi(sites, config=cfg)
    assert d.n_sites == 5


def test_from_edi_periods_descending():
    from pycsamt.models.modem.data import ModEmData
    sites = [_make_site(f"S{i:02d}", i * 500.0) for i in range(3)]
    d     = ModEmData.from_edi(sites)
    p     = d.periods
    assert np.all(p[:-1] >= p[1:])


def test_from_edi_error_applied():
    from pycsamt.models.modem.data import ModEmData
    sites = [_make_site("S00", 0.0)]
    d     = ModEmData.from_edi(sites)
    rows  = d.blocks[0]["rows"]
    errs  = [r[8] for r in rows]
    assert all(e > 0 for e in errs)


def test_from_edi_empty_raises():
    from pycsamt.models.modem.data import ModEmData
    with pytest.raises(ValueError, match="empty"):
        ModEmData.from_edi([])


def test_from_edi_roundtrip(tmp_path):
    from pycsamt.models.modem.data import ModEmData
    sites = [_make_site(f"S{i:02d}", i * 500.0) for i in range(4)]
    d     = ModEmData.from_edi(sites)
    out   = tmp_path / "test.dat"
    d.write(out)
    back  = ModEmData.read(out)
    assert back.n_sites   == d.n_sites
    assert back.n_periods == d.n_periods


# ---------------------------------------------------------------------------
# ModEmModel2D
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def model2d():
    from pycsamt.models.modem.model2d import ModEmModel2D
    return ModEmModel2D.read(_EX2D / "m0.rho")


def test_model2d_nx(model2d):
    assert model2d.nx == 100


def test_model2d_nz(model2d):
    assert model2d.nz == 31


def test_model2d_rho_shape(model2d):
    assert model2d.rho_loge.shape == (31, 100)


def test_model2d_rho_finite(model2d):
    assert np.all(np.isfinite(model2d.rho_loge))


def test_model2d_nodes(model2d):
    assert len(model2d.x_nodes) == model2d.nx + 1
    assert len(model2d.z_nodes) == model2d.nz + 1


def test_model2d_roundtrip(model2d, tmp_path):
    out = tmp_path / "m0_out.rho"
    model2d.write(out)
    from pycsamt.models.modem.model2d import ModEmModel2D
    back = ModEmModel2D.read(out)
    assert back.nx == model2d.nx
    assert back.nz == model2d.nz
    assert np.allclose(back.rho_loge, model2d.rho_loge, rtol=1e-4)


def test_model2d_halfspace():
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.data import ModEmData
    from pycsamt.models.modem.model2d import ModEmModel2D
    sites = [_make_site(f"S{i:02d}", i * 500.0) for i in range(3)]
    cfg   = ModEmConfig(mode="2d")
    d     = ModEmData.from_edi(sites, config=cfg)
    m     = ModEmModel2D.halfspace(d, config=cfg)
    assert m.nx > 0 and m.nz > 0
    # uniform body should have same loge everywhere except air
    air = cfg.n_airlayers_2d
    np.testing.assert_allclose(
        m.rho_loge[air:, :],
        np.log(cfg.initial_rho),
        rtol=1e-6,
    )


def test_model2d_missing_file_raises():
    from pycsamt.models.modem.model2d import ModEmModel2D
    with pytest.raises(FileNotFoundError):
        ModEmModel2D.read("/no/such/file.rho")


# ---------------------------------------------------------------------------
# ModEmControl
# ---------------------------------------------------------------------------

def test_control_from_config():
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.control import ModEmControl
    cfg  = ModEmConfig(max_iterations=50, target_rms=1.03)
    ctrl = ModEmControl.from_config(cfg)
    assert ctrl.max_iterations == 50
    assert ctrl.target_rms     == pytest.approx(1.03)


def test_control_read():
    from pycsamt.models.modem.control import ModEmControl
    ctrl = ModEmControl.read(_EX3D / "block2.inv")
    assert ctrl.max_iterations == 120
    assert ctrl.target_rms     == pytest.approx(1.05)
    assert ctrl.initial_lambda == pytest.approx(1.0)


def test_control_roundtrip(tmp_path):
    from pycsamt.models.modem.control import (
        ModEmConfig,
        ModEmControl,
    )
    ctrl = ModEmControl.from_config(ModEmConfig(max_iterations=80))
    out  = tmp_path / "test.inv"
    ctrl.write(out)
    back = ModEmControl.read(out)
    assert back.max_iterations == 80
    assert back.initial_lambda == pytest.approx(ctrl.initial_lambda)


def test_control_missing_file_raises():
    from pycsamt.models.modem.control import ModEmControl
    with pytest.raises(FileNotFoundError):
        ModEmControl.read("/no/such/file.inv")


# ---------------------------------------------------------------------------
# ModEmLog
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def log_2d():
    from pycsamt.models.modem.log import ModEmLog
    return ModEmLog.read(_EX2D / "Modular_NLCG.log")


def test_log_n_iter_positive(log_2d):
    assert log_2d.n_iter > 0


def test_log_rms_finite(log_2d):
    assert np.all(np.isfinite(log_2d.rms))


def test_log_rms_decreasing_start(log_2d):
    assert log_2d.rms[0] > log_2d.rms[-1]


def test_log_iterations_length(log_2d):
    # log may have NLCG restarts (iteration numbers reset), so uniqueness
    # is not guaranteed — only verify array lengths are consistent
    assert len(log_2d.iterations) == log_2d.n_iter
    assert len(log_2d.rms)        == log_2d.n_iter


def test_log_final_rms(log_2d):
    assert log_2d.final_rms == pytest.approx(log_2d.rms[-1])


def test_log_3d():
    from pycsamt.models.modem.log import ModEmLog
    log = ModEmLog.read(Path(__file__).parents[4] /
                        "ModEMv626/ModEM/examples/3D_MT/Cascadia/inverse.log")
    assert log.n_iter > 0
    assert np.all(np.isfinite(log.rms))


def test_log_missing_file_raises():
    from pycsamt.models.modem.log import ModEmLog
    with pytest.raises(FileNotFoundError):
        ModEmLog.read("/no/such/file.log")
