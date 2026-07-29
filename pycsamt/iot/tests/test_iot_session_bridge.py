"""Tests for the FieldSession -> data-model handoff methods
(:meth:`to_edifiles`, :meth:`to_sites_collection`).
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.iot import FieldSession, StationConfig, impedance_to_z

FREQ = np.logspace(4, 0, 10)


def _session():
    session = FieldSession("SURV1", method="amt")
    session.add_station(StationConfig("S01", lat=6.5, lon=3.4, elevation=120.0))
    session.add_station(StationConfig("S02", lat=6.6, lon=3.5, elevation=130.0))
    return session


def _impedance():
    zxy = (1 + 1j) * np.sqrt(FREQ)
    return {
        "S01": impedance_to_z(zxy, FREQ, station="S01"),
        "S02": impedance_to_z(zxy * 2, FREQ, station="S02"),
    }


def test_to_sites_still_returns_descriptors():
    # the original acquisition-descriptor behaviour is unchanged
    sites = _session().to_sites()
    assert all(isinstance(s, StationConfig) for s in sites)


def test_to_edifiles_objects_with_geometry():
    session = _session()
    edis = session.to_edifiles(_impedance())
    assert set(edis) == {"S01", "S02"}
    assert edis["S01"].station == "S01"
    assert edis["S01"].n_freq == FREQ.size
    # geometry recorded on the session flows into the EDI head
    head = edis["S02"].get_section("head")
    assert head.lat == pytest.approx(6.6)
    assert head.lon == pytest.approx(3.5)


def test_to_edifiles_write(tmp_path):
    from pycsamt.seg.edi import EDIFile

    paths = session_paths = _session().to_edifiles(
        _impedance(), write=True, savepath=str(tmp_path)
    )
    assert set(session_paths) == {"S01", "S02"}
    back = EDIFile(paths["S01"])
    assert back.station == "S01"


def test_to_edifiles_accepts_array_freq_pairs():
    session = _session()
    zxy = (1 + 1j) * np.sqrt(FREQ)
    edis = session.to_edifiles({"S01": (zxy, FREQ)})
    assert edis["S01"].n_freq == FREQ.size


def test_to_sites_collection_is_pipeline_ready():
    pytest.importorskip("pyproj")
    from pycsamt.site.base import Sites

    collection = _session().to_sites_collection(_impedance())
    assert isinstance(collection, Sites)
    assert len(collection) == 2
