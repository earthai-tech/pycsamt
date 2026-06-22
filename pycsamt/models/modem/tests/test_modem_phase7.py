# -*- coding: utf-8 -*-
"""Phase 7 tests — InputBuilder."""

import pytest
import numpy as np
from pathlib import Path


def _make_site(name, x_offset, y_offset=0.0, n_freq=5):
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


# ---------------------------------------------------------------------------
# 3D builder
# ---------------------------------------------------------------------------

@pytest.fixture
def sites_3d():
    return [_make_site(f"S{i:02d}", i * 1000.0, j * 1000.0)
            for i in range(3) for j in range(3)]


def test_builder_3d_creates_files(tmp_path, sites_3d):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.config   import ModEmConfig
    cfg = ModEmConfig(mode="3d", nz=6, n_airlayers=2)
    b   = InputBuilder(config=cfg)
    files = b.build(sites_3d, workdir=tmp_path)
    assert (tmp_path / "data.dat").exists()
    assert (tmp_path / "m0.ws").exists()
    assert (tmp_path / "covariance.cov").exists()
    assert (tmp_path / "control.inv").exists()


def test_builder_3d_files_key(tmp_path, sites_3d):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.config   import ModEmConfig
    cfg = ModEmConfig(mode="3d", nz=6, n_airlayers=2)
    b   = InputBuilder(config=cfg)
    files = b.build(sites_3d, workdir=tmp_path)
    assert set(files) >= {"data", "model", "covariance", "control"}


def test_builder_3d_model_attribute(tmp_path, sites_3d):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.model3d  import ModEmModel3D
    from pycsamt.models.modem.config   import ModEmConfig
    cfg = ModEmConfig(mode="3d", nz=6, n_airlayers=2)
    b   = InputBuilder(config=cfg)
    b.build(sites_3d, workdir=tmp_path)
    assert isinstance(b.model, ModEmModel3D)
    assert b.model.nx > 0


def test_builder_3d_data_attribute(tmp_path, sites_3d):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.data     import ModEmData
    from pycsamt.models.modem.config   import ModEmConfig
    cfg = ModEmConfig(mode="3d", nz=6, n_airlayers=2)
    b   = InputBuilder(config=cfg)
    b.build(sites_3d, workdir=tmp_path)
    assert isinstance(b.data, ModEmData)
    assert b.data.n_sites == len(sites_3d)


def test_builder_3d_covariance_dims(tmp_path, sites_3d):
    from pycsamt.models.modem.builder    import InputBuilder
    from pycsamt.models.modem.covariance import ModEmCovariance
    from pycsamt.models.modem.config      import ModEmConfig
    cfg = ModEmConfig(mode="3d", nz=6, n_airlayers=2)
    b   = InputBuilder(config=cfg)
    b.build(sites_3d, workdir=tmp_path)
    cov = ModEmCovariance.read(tmp_path / "covariance.cov")
    assert cov.nx_earth == b.model.nx
    assert cov.ny_earth == b.model.ny
    assert cov.nz_earth == b.model.nz - cfg.n_airlayers


def test_builder_3d_custom_filenames(tmp_path, sites_3d):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.config   import ModEmConfig
    cfg = ModEmConfig(mode="3d", nz=6, n_airlayers=2)
    b   = InputBuilder(config=cfg)
    b.build(sites_3d, workdir=tmp_path,
            data_filename="obs.dat",
            model_filename="start.ws",
            cov_filename="smooth.cov",
            ctrl_filename="run.inv")
    assert (tmp_path / "obs.dat").exists()
    assert (tmp_path / "start.ws").exists()
    assert (tmp_path / "smooth.cov").exists()
    assert (tmp_path / "run.inv").exists()


def test_builder_3d_empty_source_raises(tmp_path):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.config   import ModEmConfig
    cfg = ModEmConfig(mode="3d")
    b   = InputBuilder(config=cfg)
    with pytest.raises(ValueError, match="empty"):
        b.build([], workdir=tmp_path)


# ---------------------------------------------------------------------------
# 2D builder
# ---------------------------------------------------------------------------

@pytest.fixture
def sites_2d():
    return [_make_site(f"S{i:02d}", 0.0, i * 500.0) for i in range(4)]


def test_builder_2d_creates_files(tmp_path, sites_2d):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.config   import ModEmConfig
    cfg = ModEmConfig(mode="2d", nz_2d=10, n_airlayers_2d=3,
                      component_type="TE_Impedance")
    b   = InputBuilder(config=cfg)
    files = b.build(sites_2d, workdir=tmp_path)
    assert (tmp_path / "data.dat").exists()
    assert (tmp_path / "m0.rho").exists()
    assert (tmp_path / "control.inv").exists()


def test_builder_2d_no_covariance(tmp_path, sites_2d):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.config   import ModEmConfig
    cfg = ModEmConfig(mode="2d", nz_2d=10, n_airlayers_2d=3,
                      component_type="TE_Impedance")
    b   = InputBuilder(config=cfg)
    files = b.build(sites_2d, workdir=tmp_path)
    assert "covariance" not in files


def test_builder_2d_model_type(tmp_path, sites_2d):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.model2d  import ModEmModel2D
    from pycsamt.models.modem.config   import ModEmConfig
    cfg = ModEmConfig(mode="2d", nz_2d=10, n_airlayers_2d=3,
                      component_type="TE_Impedance")
    b   = InputBuilder(config=cfg)
    b.build(sites_2d, workdir=tmp_path)
    assert isinstance(b.model, ModEmModel2D)


# ---------------------------------------------------------------------------
# build_from_data
# ---------------------------------------------------------------------------

def test_builder_build_from_data(tmp_path, sites_3d):
    from pycsamt.models.modem.builder import InputBuilder
    from pycsamt.models.modem.data     import ModEmData
    from pycsamt.models.modem.config   import ModEmConfig
    cfg  = ModEmConfig(mode="3d", nz=6, n_airlayers=2)
    data = ModEmData.from_edi(sites_3d, config=cfg)
    b    = InputBuilder(config=cfg)
    files = b.build_from_data(data, workdir=tmp_path)
    assert "model" in files and "control" in files
    assert (tmp_path / "m0.ws").exists()
