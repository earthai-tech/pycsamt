from __future__ import annotations

import pytest

from pycsamt.airborne.ztem import ZTEMReferenceStation, ZTEMSystemSpec
from pycsamt.metadata import SiteMeta


def test_system_spec_defaults_build_instrument_meta():
    spec = ZTEMSystemSpec()
    meta = spec.to_instrument_meta(serial="SN1", software_version="1.0")
    assert meta.system == "ZTEM"
    assert meta.electric_sensor is None


def test_reference_station_rejects_non_sitemeta():
    with pytest.raises(TypeError):
        ZTEMReferenceStation(site="not-a-sitemeta")


def test_reference_station_rejects_wrong_magnetic_channels():
    with pytest.raises(ValueError):
        ZTEMReferenceStation(magnetic_channels=("Ex", "Ey"))


def test_reference_station_preferred_id_falls_back_to_site_name():
    ref = ZTEMReferenceStation(site=SiteMeta(site_id="BASE01"))
    assert ref.preferred_id == "BASE01"


def test_reference_station_preferred_id_none_when_nothing_available():
    ref = ZTEMReferenceStation()
    assert ref.preferred_id is None


def test_reference_station_explicit_id_takes_precedence():
    ref = ZTEMReferenceStation(
        station_id="EXPLICIT", site=SiteMeta(site_id="BASE01"),
    )
    assert ref.preferred_id == "EXPLICIT"
