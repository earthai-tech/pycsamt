"""Tests for InversionResult._load() and iter2dat()."""

import math
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).parents[4] / "data" / "occam2D"

pytestmark = pytest.mark.skipif(
    not DATA_DIR.exists(), reason="sample occam2D data not available"
)


@pytest.fixture(scope="module")
def result():
    from pycsamt.models.occam2d.results import InversionResult
    return InversionResult(workdir=DATA_DIR)


# -----------------------------------------------------------------------
# _load() — file discovery
# -----------------------------------------------------------------------

def test_result_has_iter_files(result):
    assert len(result.iter_files) > 0


def test_result_has_resp_files(result):
    assert len(result.resp_files) > 0


def test_result_log_loaded(result):
    assert result.log is not None


def test_result_mesh_loaded(result):
    assert result.mesh is not None
    assert result.mesh.n_xcells == 576


def test_result_model_loaded(result):
    assert result.model is not None
    assert result.model.n_layers == 26


def test_result_best_iter_loaded(result):
    assert result.best_iter is not None
    assert result.best_iter.iteration == 17


def test_result_response_loaded(result):
    assert result.response is not None
    assert result.response.n_data == 1391


# -----------------------------------------------------------------------
# rho_2d
# -----------------------------------------------------------------------

def test_result_rho_2d_shape(result):
    assert result.rho_2d is not None
    assert result.rho_2d.shape == (31, 576)


def test_result_rho_2d_finite(result):
    finite = result.rho_2d[~np.isnan(result.rho_2d)]
    assert finite.size > 0


def test_result_rho_2d_range(result):
    finite = result.rho_2d[~np.isnan(result.rho_2d)]
    # log10(rho) values should be in a plausible range for crustal resistivity
    assert float(finite.min()) > -5.0
    assert float(finite.max()) <  6.0


# -----------------------------------------------------------------------
# Convenience properties
# -----------------------------------------------------------------------

def test_result_final_rms(result):
    assert math.isclose(result.final_rms, 0.9977012, rel_tol=1e-4)


def test_result_n_iterations(result):
    assert result.n_iterations >= 1


def test_result_summary(result):
    s = result.summary()
    assert "workdir" in s
    assert "RMS" in s


# -----------------------------------------------------------------------
# iter2dat
# -----------------------------------------------------------------------

def test_iter2dat_creates_file(result, tmp_path):
    out = tmp_path / "model.dat"
    result.iter2dat(out)
    assert out.exists()


def test_iter2dat_row_count(result, tmp_path):
    out = tmp_path / "model.dat"
    result.iter2dat(out)
    lines = out.read_text().splitlines()
    # Each non-NaN cell contributes one line; boundary NaN cells excluded
    assert len(lines) > 0


def test_iter2dat_column_count(result, tmp_path):
    out = tmp_path / "model.dat"
    result.iter2dat(out)
    first = out.read_text().splitlines()[0].split()
    assert len(first) == 3   # x  z  log10_rho


def test_iter2dat_x_centered(result, tmp_path):
    out = tmp_path / "model.dat"
    result.iter2dat(out)
    xs = [float(l.split()[0]) for l in out.read_text().splitlines()]
    # x-coordinates should straddle 0 (positive and negative)
    assert min(xs) < 0
    assert max(xs) > 0


def test_iter2dat_z_positive(result, tmp_path):
    out = tmp_path / "model.dat"
    result.iter2dat(out)
    zs = [float(l.split()[1]) for l in out.read_text().splitlines()]
    assert min(zs) >= 0


# -----------------------------------------------------------------------
# Error paths
# -----------------------------------------------------------------------

def test_bad_workdir_raises():
    from pycsamt.models.occam2d.results import InversionResult
    with pytest.raises(NotADirectoryError):
        InversionResult(workdir="/nonexistent/path")
