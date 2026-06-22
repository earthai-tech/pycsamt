# -*- coding: utf-8 -*-
"""Phase 3 tests — ModEmModel3D, ModEmCovariance."""

import pytest
import numpy as np
from pathlib import Path

_EX3D = Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "3D_MT" / "BLOCK2"

pytestmark = pytest.mark.skipif(
    not _EX3D.exists(), reason="ModEMv626 example data not available"
)


# ---------------------------------------------------------------------------
# ModEmModel3D — read from file
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def model3d():
    from pycsamt.models.modem.model3d import ModEmModel3D
    return ModEmModel3D.read(_EX3D / "m0.ws")


def test_model3d_nx(model3d):
    assert model3d.nx == 21


def test_model3d_ny(model3d):
    assert model3d.ny == 28


def test_model3d_nz(model3d):
    assert model3d.nz == 11


def test_model3d_rho_shape(model3d):
    assert model3d.rho_loge.shape == (11, 28, 21)


def test_model3d_rho_finite(model3d):
    assert np.all(np.isfinite(model3d.rho_loge))


def test_model3d_nodes(model3d):
    assert len(model3d.x_nodes) == model3d.nx + 1
    assert len(model3d.y_nodes) == model3d.ny + 1
    assert len(model3d.z_nodes) == model3d.nz + 1


def test_model3d_rho_linear_positive(model3d):
    assert np.all(model3d.rho_linear > 0)


def test_model3d_shape_property(model3d):
    assert model3d.shape == (11, 28, 21)


def test_model3d_roundtrip(model3d, tmp_path):
    out = tmp_path / "m0_out.ws"
    model3d.write(out)
    from pycsamt.models.modem.model3d import ModEmModel3D
    back = ModEmModel3D.read(out)
    assert back.nx == model3d.nx
    assert back.ny == model3d.ny
    assert back.nz == model3d.nz
    assert np.allclose(back.rho_loge, model3d.rho_loge, rtol=1e-4)


def test_model3d_missing_file_raises():
    from pycsamt.models.modem.model3d import ModEmModel3D
    with pytest.raises(FileNotFoundError):
        ModEmModel3D.read("/no/such/file.ws")


# ---------------------------------------------------------------------------
# ModEmModel3D — halfspace construction
# ---------------------------------------------------------------------------

def _make_site(name, x_offset, y_offset, n_freq=5):
    """Synthetic 1-D site."""
    freqs = np.logspace(2, -1, n_freq)
    rho   = 100.0
    omega = 2 * np.pi * freqs
    mu0   = 4 * np.pi * 1e-7
    z_mag = np.sqrt(omega * mu0 * rho)
    z_val = z_mag * (1.0 + 1.0j) / np.sqrt(2)

    z_arr = np.zeros((n_freq, 2, 2), dtype=complex)
    z_arr[:, 0, 1] =  z_val
    z_arr[:, 1, 0] = -z_val

    class _Site:
        pass
    s = _Site()
    s.name   = name
    s.coords = (x_offset, y_offset, 0.0)
    s.freq   = freqs
    s.z      = z_arr
    s.z_err  = np.abs(z_arr) * 0.05
    return s


def test_model3d_halfspace_shape():
    from pycsamt.models.modem.data   import ModEmData
    from pycsamt.models.modem.model3d import ModEmModel3D
    from pycsamt.models.modem.config  import ModEmConfig
    sites = [_make_site(f"S{i:02d}", i * 1000.0, j * 1000.0)
             for i in range(3) for j in range(3)]
    cfg   = ModEmConfig(mode="3d")
    d     = ModEmData.from_edi(sites, config=cfg)
    m     = ModEmModel3D.halfspace(d, config=cfg)
    assert m.nx > 0 and m.ny > 0 and m.nz > 0
    assert m.rho_loge.shape == (m.nz, m.ny, m.nx)


def test_model3d_halfspace_air_layers():
    from pycsamt.models.modem.data    import ModEmData
    from pycsamt.models.modem.model3d import ModEmModel3D
    from pycsamt.models.modem.config   import ModEmConfig
    sites = [_make_site(f"S{i:02d}", i * 1000.0, 0.0) for i in range(3)]
    cfg   = ModEmConfig(mode="3d", n_airlayers=3)
    d     = ModEmData.from_edi(sites, config=cfg)
    m     = ModEmModel3D.halfspace(d, config=cfg)
    assert m.n_air == 3
    # air layers should have very high resistivity
    assert np.all(m.rho_loge[:3, :, :] > 20)


def test_model3d_halfspace_earth_uniform():
    from pycsamt.models.modem.data    import ModEmData
    from pycsamt.models.modem.model3d import ModEmModel3D
    from pycsamt.models.modem.config   import ModEmConfig
    sites = [_make_site(f"S{i:02d}", i * 1000.0, 0.0) for i in range(3)]
    cfg   = ModEmConfig(mode="3d", n_airlayers=3, initial_rho=50.0)
    d     = ModEmData.from_edi(sites, config=cfg)
    m     = ModEmModel3D.halfspace(d, config=cfg)
    n_air = cfg.n_airlayers
    np.testing.assert_allclose(
        m.rho_loge[n_air:, :, :],
        np.log(cfg.initial_rho),
        rtol=1e-6,
    )


def test_model3d_halfspace_roundtrip(tmp_path):
    from pycsamt.models.modem.data    import ModEmData
    from pycsamt.models.modem.model3d import ModEmModel3D
    from pycsamt.models.modem.config   import ModEmConfig
    sites = [_make_site(f"S{i:02d}", i * 1000.0, j * 1000.0)
             for i in range(2) for j in range(2)]
    cfg = ModEmConfig(mode="3d", nz=10, n_airlayers=2)
    d   = ModEmData.from_edi(sites, config=cfg)
    m   = ModEmModel3D.halfspace(d, config=cfg)
    out = tmp_path / "m_half.ws"
    m.write(out)
    back = ModEmModel3D.read(out)
    assert back.nx == m.nx and back.ny == m.ny and back.nz == m.nz
    assert np.allclose(back.rho_loge, m.rho_loge, rtol=1e-4)


# ---------------------------------------------------------------------------
# ModEmCovariance — read from file
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def cov():
    from pycsamt.models.modem.covariance import ModEmCovariance
    return ModEmCovariance.read(_EX3D / "example.cov")


def test_cov_dims(cov):
    assert cov.nx_earth == 21
    assert cov.ny_earth == 28
    assert cov.nz_earth == 11


def test_cov_smooth_x_len(cov):
    assert len(cov.smooth_x) == cov.nz_earth


def test_cov_smooth_y_len(cov):
    assert len(cov.smooth_y) == cov.nz_earth


def test_cov_smooth_z_positive(cov):
    assert cov.smooth_z > 0


def test_cov_n_smooth_iter(cov):
    assert cov.n_smooth_iter == 2


def test_cov_exceptions(cov):
    assert len(cov.exceptions) == 1
    assert cov.exceptions[0] == (2, 4, 0.0)


def test_cov_mask_blocks(cov):
    assert len(cov.mask_blocks) >= 1
    blk = cov.mask_blocks[0]
    assert blk["mask"].shape == (cov.nx_earth, cov.ny_earth)


def test_cov_mask_values(cov):
    for blk in cov.mask_blocks:
        assert np.all(blk["mask"] >= 0)
        assert np.all(blk["mask"] <= 9)


def test_cov_roundtrip(cov, tmp_path):
    out = tmp_path / "test.cov"
    cov.write(out)
    from pycsamt.models.modem.covariance import ModEmCovariance
    back = ModEmCovariance.read(out)
    assert back.nx_earth      == cov.nx_earth
    assert back.ny_earth      == cov.ny_earth
    assert back.nz_earth      == cov.nz_earth
    assert back.n_smooth_iter == cov.n_smooth_iter
    assert len(back.exceptions) == len(cov.exceptions)
    assert len(back.mask_blocks) == len(cov.mask_blocks)
    np.testing.assert_allclose(back.smooth_x, cov.smooth_x)
    np.testing.assert_allclose(back.smooth_y, cov.smooth_y)


def test_cov_missing_file_raises():
    from pycsamt.models.modem.covariance import ModEmCovariance
    with pytest.raises(FileNotFoundError):
        ModEmCovariance.read("/no/such/file.cov")


# ---------------------------------------------------------------------------
# ModEmCovariance — from_model construction
# ---------------------------------------------------------------------------

def test_cov_from_model_dims():
    from pycsamt.models.modem.data       import ModEmData
    from pycsamt.models.modem.model3d    import ModEmModel3D
    from pycsamt.models.modem.covariance import ModEmCovariance
    from pycsamt.models.modem.config      import ModEmConfig
    sites = [_make_site(f"S{i:02d}", i * 1000.0, j * 1000.0)
             for i in range(2) for j in range(2)]
    cfg = ModEmConfig(mode="3d", nz=8, n_airlayers=2)
    d   = ModEmData.from_edi(sites, config=cfg)
    m   = ModEmModel3D.halfspace(d, config=cfg)
    cov = ModEmCovariance.from_model(m, config=cfg)
    assert cov.nz_earth == m.nz - cfg.n_airlayers
    assert cov.nx_earth == m.nx
    assert cov.ny_earth == m.ny


def test_cov_from_model_uniform_mask():
    from pycsamt.models.modem.data       import ModEmData
    from pycsamt.models.modem.model3d    import ModEmModel3D
    from pycsamt.models.modem.covariance import ModEmCovariance
    from pycsamt.models.modem.config      import ModEmConfig
    sites = [_make_site(f"S{i:02d}", i * 1000.0, 0.0) for i in range(3)]
    cfg = ModEmConfig(mode="3d", nz=6, n_airlayers=2)
    d   = ModEmData.from_edi(sites, config=cfg)
    m   = ModEmModel3D.halfspace(d, config=cfg)
    cov = ModEmCovariance.from_model(m, config=cfg)
    assert len(cov.mask_blocks) == 1
    blk = cov.mask_blocks[0]
    assert blk["mask"].shape == (cov.nx_earth, cov.ny_earth)
    assert np.all(blk["mask"] == 1)
