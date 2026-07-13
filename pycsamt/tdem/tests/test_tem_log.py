"""Tests for Zonge TEMAVG processing log readers."""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.tdem import (
    TEMLog,
    is_tem_log_file,
    read_temavg_survey,
)
from pycsamt.tdem.io import read_tem_log

DATA_DIR = Path(__file__).parents[3] / "data" / "TEMAVG" / "JIANGSU"
LOG_FILE = DATA_DIR / "TEM100.LOG"

pytestmark = pytest.mark.skipif(
    not LOG_FILE.exists(),
    reason="TEMAVG sample data not available",
)


def test_is_tem_log_file():
    """LOG signature detection should recognize real TEMAVG logs."""
    assert is_tem_log_file(LOG_FILE) is True
    assert is_tem_log_file(DATA_DIR / "TEM100.AVG") is False


def test_read_log_metadata_and_shape():
    """Reader should parse stable provenance metadata."""
    log = TEMLog.read(LOG_FILE)

    assert log.version == "7.77"
    assert log.processed == "24 Oct 24"
    assert log.metadata["field_file"] == "TEM100.FLD"
    assert log.metadata["clock_type"] == "BINARY"
    assert log.metadata["clock_resolution"] == pytest.approx(1.9)
    assert log.metadata["avg_file"] == "TEM100.AVG"
    assert log.metadata["avg_dataset_count"] == 51
    assert log.metadata["closed"] is True
    assert log.n_records == 51
    assert len(log.stations) == 51


def test_first_log_record_values():
    """First acquisition summary row should match the sample log."""
    log = read_tem_log(LOG_FILE)
    rec = log.records[0]

    assert rec.station == pytest.approx(100.0)
    assert rec.frequency == pytest.approx(16.0)
    assert rec.loop == "INL"
    assert rec.component == "Hz"
    assert rec.duty == "50%"
    assert rec.first_window == pytest.approx(58.32)
    assert rec.rx_moment == pytest.approx(10000.0)
    assert rec.time_base == "31b"


def test_raw_processing_modes_are_preserved():
    """Legacy processing-mode tables should remain inspectable."""
    log = TEMLog.read(LOG_FILE)

    has_processing_modes = any(
        "PROCESSING MODES USED" in line for line in log.raw_modes
    )
    assert has_processing_modes
    assert len(log.raw_modes) > 10


def test_log_to_records_contains_expected_columns():
    """Record dictionaries should expose acquisition metadata."""
    log = TEMLog.read(LOG_FILE)
    row = log.to_records()[0]

    for key in (
        "source_file",
        "station",
        "frequency",
        "loop",
        "component",
        "first_window_us",
        "rx_moment",
        "field_file",
    ):
        assert key in row


def test_survey_parses_log_files():
    """Survey reader should expose parsed LOG files by stem."""
    survey = read_temavg_survey(DATA_DIR)

    assert survey.n_log_files == 55
    assert "TEM100" in survey.log_files
    assert survey.get_log("TEM100").n_records == 51
