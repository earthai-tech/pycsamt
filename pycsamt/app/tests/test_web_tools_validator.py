from __future__ import annotations

import numpy as np

from pycsamt.app.web.callbacks.tools import (
    _validate_sites_rows,
)


class _FakeZ:
    def __init__(self, freq, z, z_err=None):
        self.freq = np.asarray(freq, dtype=float)
        self.z = np.asarray(z, dtype=complex)
        self.z_err = z_err


class _FakeEDI:
    def __init__(self, station, lat, lon, z_obj):
        self.station = station
        self.lat = lat
        self.lon = lon
        self.Z = z_obj


def _z(freq):
    arr = np.zeros((len(freq), 2, 2), dtype=complex)
    arr[:, 0, 1] = 1.0 + 1.0j
    arr[:, 1, 0] = -1.0 - 1.0j
    return arr


def test_web_validator_rows_report_pass_warn_and_duplicates():
    freq = np.asarray([1.0, 10.0, 100.0])
    good = _FakeEDI(
        "S01",
        10.0,
        20.0,
        _FakeZ(freq, _z(freq), z_err=np.ones((3, 2, 2))),
    )
    warn = _FakeEDI(
        "S01",
        None,
        20.0,
        _FakeZ(freq[:1], _z(freq[:1]), z_err=None),
    )

    rows = _validate_sites_rows(
        [good, warn],
        store_data={"station_records": [{"ID": "S01", "Line": "L1"}]},
        checks={"coords", "z", "errors", "freq", "nan", "duplicates"},
        min_freq=3,
    )

    assert len(rows) == 2
    assert rows[0]["Status"] == "WARN"
    assert rows[0]["Z ok"] == "YES"
    assert rows[1]["Status"] == "WARN"
    assert rows[1]["N freq"] == 1
    assert "Duplicate station identifier" in rows[0]["Issues"]
    assert "Missing latitude/longitude" in rows[1]["Issues"]
    assert "Missing impedance error estimates" in rows[1]["Issues"]


def test_web_validator_rows_fail_missing_z():
    site = _FakeEDI("S02", 10.0, 20.0, None)

    rows = _validate_sites_rows(
        [site],
        store_data={"station_records": [{"ID": "S02", "Line": "L2"}]},
        checks={"coords", "z", "freq"},
        min_freq=3,
    )

    assert rows[0]["Status"] == "FAIL"
    assert rows[0]["Z ok"] == "NO"
    assert "Missing impedance tensor Z" in rows[0]["Issues"]
