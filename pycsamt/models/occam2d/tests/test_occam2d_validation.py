"""Tests for file-type validators."""

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parents[4] / "data" / "occam2D"


@pytest.mark.skipif(not DATA_DIR.exists(), reason="sample data not available")
def test_detect_data_file():
    from pycsamt.models.occam2d.validation import (
        OccamFileType,
        detect_file_type,
        is_data_file,
    )
    p = DATA_DIR / "OccamDataFile.dat"
    assert is_data_file(p)
    assert detect_file_type(p) == OccamFileType.DATA


@pytest.mark.skipif(not DATA_DIR.exists(), reason="sample data not available")
def test_detect_iter_file():
    from pycsamt.models.occam2d.validation import (
        OccamFileType,
        detect_file_type,
        is_iter_file,
    )
    p = DATA_DIR / "ITER17.iter"
    assert is_iter_file(p)
    assert detect_file_type(p) == OccamFileType.ITER


@pytest.mark.skipif(not DATA_DIR.exists(), reason="sample data not available")
def test_detect_model_file():
    from pycsamt.models.occam2d.validation import (
        OccamFileType,
        detect_file_type,
        is_model_file,
    )
    p = DATA_DIR / "Occam2DModel"
    assert is_model_file(p)
    assert detect_file_type(p) == OccamFileType.MODEL


def test_missing_file_returns_false():
    # Predicate functions return False for non-existent paths (no raise)
    from pycsamt.models.occam2d.validation import is_data_file
    assert is_data_file("/nonexistent/path/file.dat") is False
