# -*- coding: utf-8 -*-
"""Tests for file-type validators."""

import pytest
from pathlib import Path

DATA_DIR = Path(__file__).parents[5] / "data" / "occam2D"


@pytest.mark.skipif(not DATA_DIR.exists(), reason="sample data not available")
def test_detect_data_file():
    from pycsamt.models.occam2d.validation import is_data_file, detect_file_type, OccamFileType
    p = DATA_DIR / "OccamDataFile.dat"
    assert is_data_file(p)
    assert detect_file_type(p) == OccamFileType.DATA


@pytest.mark.skipif(not DATA_DIR.exists(), reason="sample data not available")
def test_detect_iter_file():
    from pycsamt.models.occam2d.validation import is_iter_file, detect_file_type, OccamFileType
    p = DATA_DIR / "ITER17.iter"
    assert is_iter_file(p)
    assert detect_file_type(p) == OccamFileType.ITER


@pytest.mark.skipif(not DATA_DIR.exists(), reason="sample data not available")
def test_detect_model_file():
    from pycsamt.models.occam2d.validation import is_model_file, OccamFileType, detect_file_type
    p = DATA_DIR / "Occam2DModel"
    assert is_model_file(p)
    assert detect_file_type(p) == OccamFileType.MODEL


def test_missing_file_raises():
    from pycsamt.models.occam2d.validation import is_data_file
    with pytest.raises(ValueError):
        is_data_file("/nonexistent/path/file.dat")
