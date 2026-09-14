from __future__ import annotations

import pytest

from pycsamt.airborne.mobilemt import MobileMTReferenceStation, MobileMTSystemSpec
from pycsamt.metadata import SiteMeta


# ─────────────────────────────────────────────────────────────────────────
# MobileMTSystemSpec.validate
# ─────────────────────────────────────────────────────────────────────────


def test_system_spec_rejects_non_finite_frequency_range():
    with pytest.raises(ValueError):
        MobileMTSystemSpec(nominal_frequency_range_hz=(float("nan"), 100.0))
    with pytest.raises(ValueError):
        MobileMTSystemSpec(nominal_frequency_range_hz=(1.0, float("inf")))


def test_system_spec_rejects_bad_frequency_ordering():
    with pytest.raises(ValueError):
        MobileMTSystemSpec(nominal_frequency_range_hz=(0.0, 100.0))
    with pytest.raises(ValueError):
        MobileMTSystemSpec(nominal_frequency_range_hz=(100.0, 1.0))


def test_system_spec_rejects_non_positive_max_windows():
    with pytest.raises(ValueError):
        MobileMTSystemSpec(nominal_max_frequency_windows=0)


def test_system_spec_rejects_non_positive_sampling_rate():
    with pytest.raises(ValueError):
        MobileMTSystemSpec(nominal_sampling_rate_hz=0.0)
    with pytest.raises(ValueError):
        MobileMTSystemSpec(nominal_sampling_rate_hz=float("nan"))


def test_system_spec_rejects_wrong_input_channels():
    with pytest.raises(ValueError):
        MobileMTSystemSpec(input_channels=("Hx", "Hy"))


def test_system_spec_rejects_wrong_output_channels():
    with pytest.raises(ValueError):
        MobileMTSystemSpec(output_channels=("Ex", "Ey"))


def test_system_spec_nominal_frequency_mask_and_instrument_meta():
    spec = MobileMTSystemSpec()
    mask = spec.nominal_frequency_mask([1.0, 1e6, float("nan")])
    assert mask.tolist() == [
        spec.nominal_frequency_range_hz[0] <= 1.0
        <= spec.nominal_frequency_range_hz[1],
        1e6 <= spec.nominal_frequency_range_hz[1],
        False,
    ]
    meta = spec.to_instrument_meta(serial="SN1", software_version="1.0")
    assert meta.system == "MobileMT"
    assert meta.serial == "SN1"


# ─────────────────────────────────────────────────────────────────────────
# MobileMTReferenceStation.validate / preferred_id
# ─────────────────────────────────────────────────────────────────────────


def test_reference_station_rejects_non_sitemeta():
    with pytest.raises(TypeError):
        MobileMTReferenceStation(site="not-a-sitemeta")


def test_reference_station_rejects_wrong_electric_channels():
    with pytest.raises(ValueError):
        MobileMTReferenceStation(electric_channels=("Hx", "Hy"))


def test_reference_station_preferred_id_falls_back_to_site_name():
    ref = MobileMTReferenceStation(site=SiteMeta(site_id="BASE99"))
    assert ref.preferred_id == "BASE99"


def test_reference_station_preferred_id_none_when_nothing_available():
    ref = MobileMTReferenceStation()
    assert ref.preferred_id is None
