from __future__ import annotations

import pytest

from pycsamt.airborne.afmag import (
    AFMAGReferenceStation,
    AirMtSystemSpec,
    OriginalAFMAGSystemSpec,
)
from pycsamt.metadata import SiteMeta


def test_original_afmag_system_spec_rejects_wrong_coil_count():
    with pytest.raises(ValueError):
        OriginalAFMAGSystemSpec(coil_count=3)


def test_original_afmag_system_spec_defaults_build_instrument_meta():
    spec = OriginalAFMAGSystemSpec()
    meta = spec.to_instrument_meta(serial="SN1")
    assert meta.system == "AFMAG (original comparator)"
    assert meta.electric_sensor is None


def test_airmt_system_spec_practical_frequency_mask_and_instrument_meta():
    spec = AirMtSystemSpec()
    mask = spec.practical_frequency_mask([1.0, 1e9, float("nan")])
    assert mask[-1] == False  # noqa: E712 - explicit numpy bool check
    meta = spec.to_instrument_meta(serial="SN2")
    assert meta.system == "AirMt / tensor AFMAG"


def test_afmag_reference_station_rejects_non_sitemeta():
    with pytest.raises(TypeError):
        AFMAGReferenceStation(site="not-a-sitemeta")


def test_afmag_reference_station_preferred_id_falls_back_to_site_name():
    ref = AFMAGReferenceStation(site=SiteMeta(site_id="REF01"))
    assert ref.preferred_id == "REF01"


def test_afmag_reference_station_preferred_id_none_when_nothing_available():
    ref = AFMAGReferenceStation()
    assert ref.preferred_id is None


def test_afmag_reference_station_explicit_id_takes_precedence():
    ref = AFMAGReferenceStation(
        station_id="EXPLICIT", site=SiteMeta(site_id="REF01"),
    )
    assert ref.preferred_id == "EXPLICIT"
