"""Tests for :func:`pycsamt.iot.emtools_qc` -- routing IoT-acquired
impedance through emtools post-processing QC so field and downstream QC
agree. Requires the optional geospatial stack.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

pytest.importorskip("pyproj")

from pycsamt.iot import (  # noqa: E402
    FieldSession,
    StationConfig,
    emtools_qc,
    impedance_to_z,
)

FREQ = np.logspace(4, 0, 12)


def _session():
    session = FieldSession("SURV1", method="amt")
    session.add_station(
        StationConfig("S01", lat=6.5, lon=3.4, elevation=120.0)
    )
    session.add_station(
        StationConfig("S02", lat=6.6, lon=3.5, elevation=130.0)
    )
    return session


def _impedance():
    zxy = (1 + 1j) * np.sqrt(FREQ)
    err = np.full((FREQ.size, 2, 2), 0.05)
    return {
        "S01": impedance_to_z(zxy, FREQ, impedance_err=err, station="S01"),
        "S02": impedance_to_z(
            zxy * 1.5, FREQ, impedance_err=err, station="S02"
        ),
    }


@pytest.fixture(autouse=True)
def _silence_allnan():
    # a minimal EDI has no coherence spectra, so emtools nan-medians warn
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        yield


def test_qc_table_from_session():
    table = emtools_qc(_session(), _impedance())
    assert table.shape[0] == 2
    assert "station" in table.columns
    assert "skew_med" in table.columns


def test_qc_flags_from_session():
    flags = emtools_qc(_session(), _impedance(), flags=True)
    assert "flags" in flags.columns
    assert set(flags["station"]) == {"S01", "S02"}


def test_qc_from_sites_collection():
    sites = _session().to_sites_collection(_impedance())
    table = emtools_qc(sites)
    assert table.shape[0] == 2


def test_qc_from_edi_source(tmp_path):
    _session().to_edifiles(_impedance(), write=True, savepath=str(tmp_path))
    table = emtools_qc(str(tmp_path))
    assert table.shape[0] == 2
    assert set(table["station"]) == {"S01", "S02"}


def test_session_without_impedance_raises():
    with pytest.raises(ValueError):
        emtools_qc(_session())


def test_thresholds_forwarded():
    # a permissive threshold should not raise and returns per-station flags
    flags = emtools_qc(_session(), _impedance(), flags=True, min_frac_ok=0.1)
    assert flags.shape[0] == 2
