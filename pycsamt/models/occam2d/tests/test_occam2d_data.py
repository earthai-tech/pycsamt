# -*- coding: utf-8 -*-
"""Tests for OccamData."""

import pytest
from pathlib import Path

DATA_DIR = Path(__file__).parents[4] / "data" / "occam2D"
DATA_FILE = DATA_DIR / "OccamDataFile.dat"


@pytest.mark.skipif(not DATA_FILE.exists(), reason="sample data not available")
def test_read_data_file():
    from pycsamt.models.occam2d.data import OccamData
    data = OccamData.read(DATA_FILE)
    assert data.n_sites > 0
    assert data.n_frequencies > 0
    assert data.n_data > 0


def test_data_file_type_code_keys():
    from pycsamt.models.occam2d.data import DATA_TYPE_CODES
    assert "RhoTE" in DATA_TYPE_CODES
    assert "PhsTM" in DATA_TYPE_CODES
    assert DATA_TYPE_CODES["RhoTE"] == 1
