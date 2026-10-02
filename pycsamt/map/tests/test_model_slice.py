"""Tests for pycsamt.map.model_slice (real 3-D model depth slices)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.map.model_slice import _affine, model_depth_slice

ROOT = Path(__file__).resolve().parents[3]
BH = ROOT / "data" / "MT" / "broken-hill" / "final-models"


@pytest.fixture(scope="module")
def bh_model(tmp_path_factory):
    # Only the .dat/.res files are tracked; the ~11 MB .rho model is not.
    if not BH.is_dir() or not list(BH.glob("*.rho")):
        pytest.skip("Broken Hill ModEM result not present")
    from pycsamt.format import convert_engine as ce
    from pycsamt.format.text import read_pcsf_or_pcsm

    out = tmp_path_factory.mktemp("pcsf") / "bh.pcsf"
    sk = ce.detect(BH, None)
    ce.write_model(ce.build_model(sk, ce.ConvertOptions()), out, "pcsf",
                   False)
    return read_pcsf_or_pcsm(out)


def test_slice_is_cropped_and_georeferenced(bh_model):
    sl = model_depth_slice(bh_model, 500.0)
    assert sl.geo and sl.mode == "at"
    ny, nx = sl.rho.shape
    assert sl.x.shape == (ny + 1, nx + 1)
    assert nx < bh_model.resistivity.shape[2]  # padding cropped
    lat = np.asarray(bh_model.stations.lat, float)
    lon = np.asarray(bh_model.stations.lon, float)
    # every station lies inside the sliced area
    assert sl.x.min() < lon.min() and sl.x.max() > lon.max()
    assert sl.y.min() < lat.min() and sl.y.max() > lat.max()
    assert np.nanmin(sl.rho) > 0


def test_slice_matches_model_layer_at_cell_centre(bh_model):
    z = np.asarray(bh_model.geometry.z)
    sl = model_depth_slice(bh_model, float(z[5]), georeference=False)
    full = model_depth_slice(bh_model, float(z[5]), margin=1e9,
                             georeference=False)
    np.testing.assert_allclose(full.rho, bh_model.resistivity[5],
                               rtol=1e-6)
    assert sl.rho.size < full.rho.size


def test_mean_mode_is_log_mean_down_to_depth(bh_model):
    z = np.asarray(bh_model.geometry.z)
    d = float(z[3])
    full = model_depth_slice(bh_model, d, mode="mean", margin=1e9,
                             georeference=False)
    expect = 10 ** np.mean(np.log10(bh_model.resistivity[:4]), axis=0)
    np.testing.assert_allclose(full.rho, expect, rtol=1e-6)


def test_outside_depth_raises(bh_model):
    with pytest.raises(ValueError, match="outside"):
        model_depth_slice(bh_model, 1e9)


def test_affine_needs_three_non_collinear_points():
    src = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    assert _affine(src, src) is None
    src = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    dst = src * 2 + 3
    coef = _affine(src, dst)
    np.testing.assert_allclose(np.column_stack([src, np.ones(4)]) @ coef,
                               dst, atol=1e-9)
