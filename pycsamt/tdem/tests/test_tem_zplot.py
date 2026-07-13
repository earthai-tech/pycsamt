"""Tests for Zonge TEMAVG contour ``.Z`` readers."""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.tdem import (
    TEMAVG,
    TEMZPlot,
    is_tem_z_file,
    read_temavg_survey,
)
from pycsamt.tdem.io import read_tem_z

DATA_DIR = Path(__file__).parents[3] / "data" / "TEMAVG" / "JIANGSU"
Z_FILE = DATA_DIR / "TEM100.Z"

pytestmark = pytest.mark.skipif(
    not Z_FILE.exists(),
    reason="TEMAVG sample data not available",
)


def test_is_tem_z_file():
    """Z signature detection should recognize real contour files."""
    assert is_tem_z_file(Z_FILE) is True
    assert is_tem_z_file(DATA_DIR / "TEM100.AVG") is False


def test_read_zplot_metadata_and_shape():
    """Reader should parse header metadata and contour rows."""
    zplot = TEMZPlot.read(Z_FILE)

    assert zplot.version == "7.77"
    assert zplot.metadata["date"] == "21-07-04"
    assert zplot.metadata["window"] == "MAGNITUDE"
    assert zplot.metadata["component"] == "Hz"
    assert zplot.metadata["rx_area"] == pytest.approx(10000.0)
    assert zplot.n_records == 1275
    assert len(zplot.stations) == 51
    assert zplot.windows == list(range(1, 26))


def test_first_z_record_values():
    """First contour row should match the sample file."""
    zplot = read_tem_z(Z_FILE)
    rec = zplot.records[0]

    assert rec.line == 2
    assert rec.station == pytest.approx(100.0)
    assert rec.time == pytest.approx(0.05832)
    assert rec.time_s == pytest.approx(5.832e-5)
    assert rec.magnitude == pytest.approx(7.1650e4)
    assert rec.frequency == pytest.approx(16.0)
    assert rec.window == 1


def test_glued_negative_magnitude_parses():
    """Rows like ``12.2015-1.128e-01`` should parse cleanly."""
    zplot = TEMZPlot.read(Z_FILE)
    rows = zplot.rows_for_station(100.0)

    assert rows[-1].time == pytest.approx(12.2015)
    assert rows[-1].magnitude == pytest.approx(-1.128e-1)
    assert rows[-1].window == 25


def test_z_magnitude_matches_avg_magnitude():
    """Z contour magnitudes should agree with AVG magnitudes."""
    avg = TEMAVG.read(DATA_DIR / "TEM100.AVG")
    zplot = TEMZPlot.read(Z_FILE)

    avg_rows = avg.rows_for_station(160.0)
    z_rows = zplot.rows_for_station(160.0)

    assert len(avg_rows) == len(z_rows)
    for avg_row, z_row in zip(avg_rows, z_rows, strict=True):
        assert z_row.time == pytest.approx(avg_row.time, rel=5e-4)
        assert z_row.magnitude == pytest.approx(avg_row.magnitude)
        assert z_row.window == avg_row.window


def test_survey_parses_z_files():
    """Survey reader should expose parsed Z files by stem."""
    survey = read_temavg_survey(DATA_DIR)
    expected_z = len(list(DATA_DIR.glob("*.Z")))

    assert survey.n_z_files == expected_z
    assert "TEM100" in survey.z_files
    assert survey.get_z("TEM100").n_records == 1275
