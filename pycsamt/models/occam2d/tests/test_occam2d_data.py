# -*- coding: utf-8 -*-
"""Tests for OccamData.read()."""

import math
import pytest
from pathlib import Path

import numpy as np

DATA_DIR  = Path(__file__).parents[4] / "data" / "occam2D"
DATA_FILE = DATA_DIR / "OccamDataFile.dat"

pytestmark = pytest.mark.skipif(
    not DATA_FILE.exists(), reason="sample occam2D data not available"
)


@pytest.fixture(scope="module")
def data():
    from pycsamt.models.occam2d.data import OccamData
    return OccamData.read(DATA_FILE)


def test_data_format(data):
    assert data.format_str == "OCCAM2MTDATA_1.0"


def test_data_title(data):
    assert "MTPY" in data.title.upper() or "OCCAM" in data.title.upper()


def test_data_n_sites(data):
    assert data.n_sites == 47


def test_data_sites_list(data):
    assert data.sites[0] == "S00"
    assert data.sites[-1] == "S46"


def test_data_n_offsets(data):
    assert len(data.offsets) == 47


def test_data_offset_first(data):
    assert math.isclose(data.offsets[0], 0.0, abs_tol=1e-6)


def test_data_offset_last(data):
    assert math.isclose(data.offsets[-1], 2312.2, rel_tol=1e-4)


def test_data_offsets_monotone(data):
    assert bool((np.diff(data.offsets) > 0).all())


def test_data_n_frequencies(data):
    assert data.n_frequencies == 17


def test_data_freq_first(data):
    assert math.isclose(data.frequencies[0], 10.0, rel_tol=1e-4)


def test_data_freq_last(data):
    assert math.isclose(data.frequencies[-1], 1e4, rel_tol=1e-4)


def test_data_n_data(data):
    assert data.n_data == 1391


def test_data_blocks_shape(data):
    assert data.data_blocks.shape == (1391, 5)


def test_data_blocks_site_range(data):
    sites = data.data_blocks[:, 0].astype(int)
    assert int(sites.min()) >= 1
    assert int(sites.max()) <= 47


def test_data_blocks_freq_range(data):
    freqs = data.data_blocks[:, 1].astype(int)
    assert int(freqs.min()) >= 1
    assert int(freqs.max()) <= 17


def test_data_type_codes_present(data):
    codes = set(data.type_codes.tolist())
    # File uses RhoTM (5) and PhsTM (6)
    assert 5 in codes
    assert 6 in codes


def test_data_first_row(data):
    row = data.data_blocks[0]
    assert int(row[0]) == 1   # site 1
    assert int(row[1]) == 1   # freq 1
    assert int(row[2]) == 5   # type RhoTM
    assert math.isclose(row[3], 4.8243, rel_tol=1e-3)
    assert math.isclose(row[4], 6.5578, rel_tol=1e-3)


# -----------------------------------------------------------------------
# Data type codes dict (independent of sample data)
# -----------------------------------------------------------------------

def test_data_file_type_code_keys():
    from pycsamt.models.occam2d.data import DATA_TYPE_CODES
    assert "RhoTE" in DATA_TYPE_CODES
    assert "PhsTM" in DATA_TYPE_CODES
    assert DATA_TYPE_CODES["RhoTE"] == 1
    assert DATA_TYPE_CODES["PhsTM"] == 6


# -----------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------

def test_data_defaults():
    from pycsamt.models.occam2d.data import OccamData
    d = OccamData()
    assert d.n_sites == 0
    assert d.n_frequencies == 0
    assert d.n_data == 0


# -----------------------------------------------------------------------
# Error handling
# -----------------------------------------------------------------------

def test_missing_file_raises():
    from pycsamt.models.occam2d.data import OccamData
    with pytest.raises(FileNotFoundError):
        OccamData.read("/nonexistent/OccamDataFile.dat")


def test_wrong_format_raises(tmp_path):
    from pycsamt.models.occam2d.data import OccamData
    bad = tmp_path / "bad.dat"
    bad.write_text("FORMAT: WRONG_FORMAT\nTITLE: test\n")
    with pytest.raises(ValueError, match="OCCAM2MTDATA"):
        OccamData.read(bad)
