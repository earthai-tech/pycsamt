"""Tests for OccamResponse.read()."""

import math
from pathlib import Path

import numpy as np
import pytest

DATA_DIR  = Path(__file__).parents[4] / "data" / "occam2D"
RESP_FILE = DATA_DIR / "RESP17.resp"

pytestmark = pytest.mark.skipif(
    not RESP_FILE.exists(), reason="sample occam2D data not available"
)


@pytest.fixture(scope="module")
def resp():
    from pycsamt.models.occam2d.response import OccamResponse
    return OccamResponse.read(RESP_FILE)


def test_resp_n_data(resp):
    assert resp.n_data == 1391


def test_resp_data_shape(resp):
    assert resp.data.shape == (1391, 7)


def test_resp_observed_shape(resp):
    assert resp.observed.shape == (1391,)


def test_resp_modeled_shape(resp):
    assert resp.modeled.shape == (1391,)


def test_resp_residuals_shape(resp):
    assert resp.residuals.shape == (1391,)


def test_resp_site_indices_range(resp):
    assert int(resp.site_indices.min()) >= 1
    assert int(resp.site_indices.max()) <= 47


def test_resp_freq_indices_range(resp):
    assert int(resp.freq_indices.min()) >= 1
    assert int(resp.freq_indices.max()) <= 17


def test_resp_type_codes(resp):
    codes = set(resp.type_codes.tolist())
    assert 5 in codes
    assert 6 in codes


def test_resp_first_row_observed(resp):
    assert math.isclose(resp.observed[0], 4.8243, rel_tol=1e-3)


def test_resp_first_row_modeled(resp):
    assert math.isclose(resp.modeled[0], 4.9036, rel_tol=1e-3)


def test_resp_first_row_residual(resp):
    assert math.isclose(resp.residuals[0], -0.12094e-1, rel_tol=1e-2)


def test_resp_rms_positive(resp):
    assert resp.rms > 0.0


def test_resp_rms_value(resp):
    # RMS from residuals column
    expected = float(np.sqrt(np.mean(resp.residuals ** 2)))
    assert math.isclose(resp.rms, expected, rel_tol=1e-6)


def test_resp_misfit_per_site(resp):
    per_site = resp.misfit_per_site()
    assert len(per_site) > 0
    assert all(v > 0 for v in per_site.values())


def test_resp_misfit_per_frequency(resp):
    per_freq = resp.misfit_per_frequency()
    assert len(per_freq) > 0
    assert all(v > 0 for v in per_freq.values())


# -----------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------

def test_resp_defaults():
    from pycsamt.models.occam2d.response import OccamResponse
    r = OccamResponse()
    assert r.n_data == 0
    assert r.rms == 0.0
    assert r.type_codes.size == 0


# -----------------------------------------------------------------------
# Error handling
# -----------------------------------------------------------------------

def test_missing_file_raises():
    from pycsamt.models.occam2d.response import OccamResponse
    with pytest.raises(FileNotFoundError):
        OccamResponse.read("/nonexistent/RESP17.resp")


def test_empty_file_raises(tmp_path):
    from pycsamt.models.occam2d.response import OccamResponse
    empty = tmp_path / "empty.resp"
    empty.write_text("")
    with pytest.raises(ValueError, match="No valid data rows"):
        OccamResponse.read(empty)
