"""Tests for Zonge TEMAVG processed-file readers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.tdem import (
    TEMAVG,
    LateTimeTransform,
    is_temavg_file,
    read_temavg_survey,
)
from pycsamt.tdem.io import read_temavg

DATA_DIR = Path(__file__).parents[3] / "data" / "TEMAVG" / "JIANGSU"
AVG_FILE = DATA_DIR / "TEM100.AVG"

pytestmark = pytest.mark.skipif(
    not AVG_FILE.exists(),
    reason="TEMAVG sample data not available",
)


def test_is_temavg_file():
    """AVG signature detection should recognize real TEMAVG files."""
    assert is_temavg_file(AVG_FILE) is True
    assert is_temavg_file(DATA_DIR / "TEM100.LOG") is False


def test_read_temavg_metadata_and_shape():
    """Reader should parse header metadata and gate rows."""
    avg = TEMAVG.read(AVG_FILE)

    assert avg.version == "7.77"
    assert avg.metadata["TXramp"] == pytest.approx(450.0)
    assert avg.metadata["TXarea"] == pytest.approx(129600.0)
    assert avg.metadata["RXarea"] == pytest.approx(10000.0)
    assert avg.n_records == 1275
    assert len(avg.stations) == 51
    assert avg.windows == list(range(1, 26))


def test_first_record_values():
    """The first processed gate should match the sample file."""
    avg = read_temavg(AVG_FILE)
    rec = avg.records[0]

    assert rec.station == pytest.approx(100.0)
    assert rec.frequency == pytest.approx(16.0)
    assert rec.component == "Hz"
    assert rec.current == pytest.approx(10.0)
    assert rec.window == 1
    assert rec.time == pytest.approx(0.05832)
    assert rec.time_s == pytest.approx(5.832e-5)
    assert rec.magnitude == pytest.approx(7.1650e4)
    assert rec.ramp_app_res == pytest.approx(2.6955e1)
    assert rec.depth == pytest.approx(3.5106e1)


def test_rows_for_station():
    """Station slicing should return one row per time gate."""
    avg = TEMAVG.read(AVG_FILE)
    rows = avg.rows_for_station(100.0)

    assert len(rows) == 25
    assert [row.window for row in rows] == list(range(1, 26))


def test_temavg_to_soundings_builds_station_decays():
    """AVG rows should convert to one TEMSounding per station."""
    avg = TEMAVG.read(AVG_FILE)

    soundings = avg.to_soundings(stations=[100.0])
    snd = soundings[0]

    assert len(soundings) == 1
    assert snd.station_name == "TEM100_100"
    assert snd.n_gates == 25
    assert snd.current == pytest.approx(10.0)
    assert snd.tx_area == pytest.approx(129600.0)
    assert snd.rx_area == pytest.approx(10000.0)
    assert snd.loop_shape == "square"
    assert snd.loop_dims == pytest.approx((360.0,))
    assert snd.time_gates[0] == pytest.approx(5.832e-5)
    assert snd.data[0] == pytest.approx(7.1650e4 * 10.0 * 1e-6)
    assert snd.error is None


def test_temavg_to_soundings_can_export_dbdt():
    """AVG magnitudes should convert to dB/dt when requested."""
    avg = TEMAVG.read(AVG_FILE)

    snd = avg.to_soundings(stations=[100.0], data_type="dBdt")[0]

    voltage = 7.1650e4 * 10.0 * 1e-6
    assert snd.data_type == "dBdt"
    assert snd.data[0] == pytest.approx(voltage / 10000.0)


def test_nan_quality_marker_is_preserved():
    """Asterisk quality markers should become nan, not errors."""
    avg = TEMAVG.read(DATA_DIR / "TEM1020.AVG")
    rows = avg.rows_for_station(400.0)

    assert np.isnan(rows[-1].percent_magnitude)


def test_read_temavg_survey_groups_files():
    """Folder reader should parse AVG files and group companions."""
    survey = read_temavg_survey(DATA_DIR)
    expected_avg = len(list(DATA_DIR.glob("*.AVG")))
    expected_records = sum(avg.n_records for avg in survey.avg_files.values())

    assert survey.n_avg_files == expected_avg
    assert len(survey.stations) == 51
    assert "TEM100" in survey.avg_files
    assert sorted(survey.companion_files["TEM100"]) == [
        "log",
        "mde",
        "z",
    ]
    assert len(survey.to_records()) == expected_records


def test_temavg_survey_to_soundings_exports_all_stations():
    """Survey conversion should flatten AVG files into soundings."""
    survey = read_temavg_survey(DATA_DIR)

    soundings = survey.to_soundings(stems=["TEM100"])

    assert len(soundings) == 51
    assert soundings[0].station_name == "TEM100_100"
    assert soundings[0].n_gates == 25


def test_temavg_sounding_can_run_late_time_transform():
    """Generated soundings should be valid transform inputs."""
    avg = TEMAVG.read(AVG_FILE)
    snd = avg.to_soundings(stations=[100.0])[0]

    result = LateTimeTransform().transform(snd)

    assert result["station_name"] == "TEM100_100"
    assert len(result["freq"]) == snd.n_gates
    assert len(result["rho_a"]) == snd.n_gates


def test_to_records_contains_expected_columns():
    """Record dictionaries should expose plotting-ready columns."""
    avg = TEMAVG.read(AVG_FILE)
    row = avg.to_records()[0]

    for key in (
        "source_file",
        "station",
        "time_s",
        "magnitude",
        "ramp_app_res",
        "depth",
        "tx_area",
        "rx_area",
    ):
        assert key in row
