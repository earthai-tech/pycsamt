"""Tests for :mod:`pycsamt.iot.methods` — per-method acquisition profiles
and the method-aware :meth:`MonitoringConfig.for_method` factory.
"""

from __future__ import annotations

import pytest

from pycsamt.iot import (
    METHOD_PROFILES,
    EMMethod,
    MethodProfile,
    MonitoringConfig,
    assess_telemetry,
    method_profile,
    target_bands_for_method,
)


# ---------------------------------------------------------------------------
# profiles
# ---------------------------------------------------------------------------
def test_every_method_has_a_profile():
    for method in EMMethod:
        assert method in METHOD_PROFILES
        assert isinstance(METHOD_PROFILES[method], MethodProfile)


@pytest.mark.parametrize("name", ["amt", "mt", "csamt", "tdem", "tem"])
def test_method_profile_lookup_by_string(name):
    profile = method_profile(name)
    assert profile.method.value == name


def test_method_profile_accepts_enum():
    assert method_profile(EMMethod.CSAMT).method is EMMethod.CSAMT


def test_csamt_is_controlled_source():
    profile = method_profile("csamt")
    assert profile.controlled_source
    assert not profile.is_natural_source


def test_amt_and_mt_are_natural_source():
    assert method_profile("amt").is_natural_source
    assert method_profile("mt").is_natural_source
    assert method_profile("amt").powerline_sensitive


def test_mt_requires_vertical_field_channel():
    # long-period MT records Hz; AMT/CSAMT default to the 4-channel layout
    assert "hz" in method_profile("mt").required_channels
    assert "hz" not in method_profile("amt").required_channels


def test_unknown_method_has_empty_defaults():
    profile = method_profile("unknown")
    assert profile.frequency_band_hz is None
    assert profile.required_channels == ()


def test_unrecognised_method_raises():
    with pytest.raises(ValueError):
        method_profile("seismic")


def test_profile_as_dict_roundtrips_fields():
    d = method_profile("csamt").as_dict()
    assert d["method"] == "csamt"
    assert d["controlled_source"] is True
    assert d["frequency_band_hz"] == [0.125, 8192.0]


# ---------------------------------------------------------------------------
# target bands
# ---------------------------------------------------------------------------
def test_target_bands_span_and_clamp_to_method_band():
    bands = target_bands_for_method("amt")
    assert bands[0][0] == pytest.approx(1.0)          # clamped low edge
    assert bands[-1][1] == pytest.approx(10_000.0)    # clamped high edge
    # contiguous decade sub-bands
    for lo, hi in bands:
        assert lo < hi
    for i in range(len(bands) - 1):
        assert bands[i][1] == pytest.approx(bands[i + 1][0])


def test_target_bands_csamt_covers_subhertz():
    bands = target_bands_for_method("csamt")
    assert bands[0][0] == pytest.approx(0.125)
    assert bands[-1][1] == pytest.approx(8192.0)


def test_target_bands_empty_for_time_domain():
    assert target_bands_for_method("tdem") == []
    assert target_bands_for_method("unknown") == []


# ---------------------------------------------------------------------------
# MonitoringConfig.for_method
# ---------------------------------------------------------------------------
def test_for_method_seeds_band_and_channels():
    cfg = MonitoringConfig.for_method("csamt")
    assert cfg.method is EMMethod.CSAMT
    assert cfg.frequency_band_hz == (0.125, 8192.0)
    assert cfg.required_channels == ["ex", "ey", "hx", "hy"]


def test_for_method_overrides_win():
    cfg = MonitoringConfig.for_method(
        "amt", min_battery_v=11.5, required_channels=["ex", "ey"]
    )
    assert cfg.min_battery_v == pytest.approx(11.5)
    assert cfg.required_channels == ["ex", "ey"]      # override, not profile


def test_for_method_makes_assessment_method_aware():
    cfg = MonitoringConfig.for_method("csamt")
    # an MT packet with only two channels violates the CSAMT expectations
    packets = [{
        "device_id": "n1", "timestamp": 100.0, "topic": "t", "kind": "data",
        "payload": {"method": "mt", "channels": ["ex", "hy"],
                    "frequency_band_hz": [1.0, 500.0]},
    }]
    status = assess_telemetry(packets, config=cfg)
    assert "method_mismatch" in status.issues
    assert "required_channels_missing" in status.issues


def test_for_method_matching_stream_has_no_method_issues():
    cfg = MonitoringConfig.for_method("csamt")
    packets = [{
        "device_id": "n1", "timestamp": 100.0, "topic": "t", "kind": "data",
        "payload": {"method": "csamt", "channels": ["ex", "ey", "hx", "hy"],
                    "frequency_band_hz": [1.0, 5000.0]},
    }]
    status = assess_telemetry(packets, config=cfg)
    assert "method_mismatch" not in status.issues
    assert "required_channels_missing" not in status.issues
