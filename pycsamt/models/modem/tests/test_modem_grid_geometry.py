# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regressions: ModEM model files match Mod2DMT/Mod3DMT geometry.

* Mod2DMT reads 2-D (Mackie) models as earth layers only and adds 10
  air layers itself; receivers are located from the grid's left edge.
* Mod3DMT reads an explicit origin trailer as the grid corner.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from pycsamt.models.modem.builder import InputBuilder
from pycsamt.models.modem.config import ModEmConfig
from pycsamt.models.modem.model2d import ModEmModel2D
from pycsamt.models.modem.model3d import ModEmModel3D


class _Site:
    """Minimal site exposing what ModEmData.from_edi reads."""

    def __init__(self, name, lat, lon):
        self.name = name
        self.coords = (lat, lon, 0.0)
        self.freq = np.array([10.0, 1.0])
        z = np.zeros((2, 2, 2), complex)
        z[:, 0, 1] = 10 + 10j
        z[:, 1, 0] = -10 - 10j
        self.z = z
        self.z_err = np.abs(z) * 0.05


def _sites():
    # east-west profile, ~100 m spacing
    return [_Site(f"S{i}", 25.55, 113.91 + 0.001 * i)
            for i in range(5)]


def _read_y(path):
    ys = []
    for ln in path.read_text().splitlines():
        t = ln.split()
        if len(t) == 11 and not ln.startswith(("#", ">")):
            ys.append(float(t[5]))
    return np.array(ys)


def test_2d_builder_writes_earth_only_and_grid_relative_y(tmp_path):
    cfg = ModEmConfig(mode="2d", component_type="TE_Impedance",
                      n_airlayers_2d=4, nz_2d=12, cell_size_h_2d=25.0)
    b = InputBuilder(config=cfg)
    files = b.build(_sites(), workdir=tmp_path)
    # in-memory model keeps its air layers ...
    assert b.model.n_air == 4 and b.model.nz == 16
    # ... but the file has earth layers only (Mod2DMT adds air)
    back = ModEmModel2D.read(files["model"])
    assert back.nz == 12
    np.testing.assert_allclose(back.rho_loge,
                               np.log(cfg.initial_rho))
    # every receiver lies inside the grid, first one after padding
    y = _read_y(files["data"])
    width = float(np.sum(b.model.x_widths))
    pad = float(np.sum(b.model.x_widths[:cfg.n_padding_x_2d]))
    assert y.min() == pytest.approx(pad, abs=1e-3)
    assert 0 < y.min() and y.max() < width


def test_2d_write_include_air_keeps_old_behaviour(tmp_path):
    m = ModEmModel2D()
    m.x_widths = np.full(4, 10.0)
    m.z_widths = np.array([5.0, 5.0, 10.0, 20.0])
    m.rho_loge = np.zeros((4, 4))
    m.n_air = 2
    assert ModEmModel2D.read(m.write(tmp_path / "a.rho")).nz == 2
    assert ModEmModel2D.read(
        m.write(tmp_path / "b.rho", include_air=True)).nz == 4


def test_3d_from_widths_centres_origin(tmp_path):
    m = ModEmModel3D.from_widths([100.0] * 4, [50.0] * 6, [10.0, 30.0],
                                 rho=50.0)
    np.testing.assert_allclose(m.origin, [-200.0, -150.0, 0.0])
    np.testing.assert_allclose(m.rho_loge, np.log(50.0))
    assert m.rho_loge.shape == (2, 6, 4)
    back = ModEmModel3D.read(m.write(tmp_path / "m.ws"))
    np.testing.assert_allclose(back.origin, m.origin)


def test_3d_from_widths_validation():
    with pytest.raises(ValueError, match="x_widths"):
        ModEmModel3D.from_widths([0.0, 1.0], [1.0], [1.0])
    with pytest.raises(ValueError, match="shape"):
        ModEmModel3D.from_widths([1.0], [1.0], [1.0],
                                 rho=np.ones((2, 2, 2)))


def test_3d_zero_origin_warns(tmp_path, caplog):
    m = ModEmModel3D.from_widths([10.0] * 3, [10.0] * 3, [5.0],
                                 origin=[0.0, 0.0, 0.0])
    m.logger.propagate = True
    with caplog.at_level(logging.WARNING):
        m.write(tmp_path / "z.ws")
    assert "grid corner" in caplog.text
