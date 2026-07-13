"""Tests for the controlled-source transmitter telemetry schema
(:class:`~pycsamt.iot.schemas.SourcePayload`) and its ``PacketKind.SOURCE``
registration and flow through the telemetry path.
"""

from __future__ import annotations

import pytest

from pycsamt.iot import (
    DeviceConfig,
    FieldSession,
    PacketKind,
    SourcePayload,
    parse_payload,
    schema_for,
    validate_payload,
)


# ---------------------------------------------------------------------------
# registration
# ---------------------------------------------------------------------------
def test_source_packet_kind_exists():
    assert PacketKind.SOURCE.value == "source"


def test_schema_for_routes_source():
    assert schema_for("source") is SourcePayload
    assert schema_for(PacketKind.SOURCE) is SourcePayload


# ---------------------------------------------------------------------------
# alias folding and canonicalisation
# ---------------------------------------------------------------------------
def test_source_payload_folds_aliases():
    raw = {
        "tx_id": "TX1",
        "current": 9.8,
        "tx_voltage": 250.0,
        "frequency": 32.0,
        "ab_m": 100.0,
        "duty": 0.5,
        "transmitting": True,
        "tx_rx_offset": 5000.0,
        "bearing_deg": 45.0,
    }
    p = parse_payload("source", raw)
    assert isinstance(p, SourcePayload)
    assert p.source_id == "TX1"
    assert p.tx_current_a == pytest.approx(9.8)
    assert p.tx_voltage_v == pytest.approx(250.0)
    assert p.tx_frequency_hz == pytest.approx(32.0)
    assert p.dipole_length_m == pytest.approx(100.0)
    assert p.duty_cycle == pytest.approx(0.5)
    assert p.on is True
    assert p.offset_m == pytest.approx(5000.0)
    assert p.azimuth_deg == pytest.approx(45.0)


def test_source_payload_preserves_unknown_keys():
    p = parse_payload("source", {"tx_current_a": 5.0, "note": "keep me"})
    assert p.extra == {"note": "keep me"}


def test_validate_payload_canonical_keys():
    canon = validate_payload(
        "source", {"tx_id": "TX1", "current": 9.8, "extra_field": 1}
    )
    assert canon["source_id"] == "TX1"
    assert canon["tx_current_a"] == pytest.approx(9.8)
    assert canon["extra_field"] == 1
    # drop_none default removes unreported optional fields
    assert "tx_power_w" not in canon


def test_source_payload_as_dict_roundtrip():
    p = SourcePayload(source_id="TX1", tx_current_a=9.8, on=False)
    d = p.as_dict(drop_none=True)
    assert d["source_id"] == "TX1"
    assert d["on"] is False
    assert "azimuth_deg" not in d  # None dropped


def test_string_boolean_on_state():
    assert parse_payload("source", {"transmitting": "false"}).on is False
    assert parse_payload("source", {"tx_on": "yes"}).on is True


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------
def test_duty_cycle_range_enforced():
    with pytest.raises(ValueError):
        SourcePayload(duty_cycle=1.5)
    with pytest.raises(ValueError):
        SourcePayload(duty_cycle=-0.1)


def test_frequency_and_geometry_positivity():
    with pytest.raises(ValueError):
        SourcePayload(tx_frequency_hz=-1.0)
    with pytest.raises(ValueError):
        SourcePayload(dipole_length_m=0.0)
    with pytest.raises(ValueError):
        SourcePayload(offset_m=-10.0)


def test_current_and_voltage_may_be_zero_when_off():
    # a keyed-off transmitter legitimately reports zero current/voltage
    p = SourcePayload(tx_current_a=0.0, tx_voltage_v=0.0, on=False)
    assert p.tx_current_a == 0.0
    assert p.tx_voltage_v == 0.0


def test_non_finite_current_rejected():
    with pytest.raises(ValueError):
        SourcePayload(tx_current_a=float("nan"))


# ---------------------------------------------------------------------------
# flow through a session
# ---------------------------------------------------------------------------
def test_source_packet_through_field_session():
    session = FieldSession("CS1", method="csamt")
    session.add_device(DeviceConfig("tx-1", role="transmitter"))
    session.add_packet(
        {
            "device_id": "tx-1",
            "timestamp": 10.0,
            "topic": "t",
            "kind": "source",
            "payload": {
                "source_id": "TX1",
                "tx_current_a": 9.8,
                "tx_frequency_hz": 32.0,
            },
        }
    )
    assert session.n_packets == 1
    assert session.packets[0].kind is PacketKind.SOURCE
    # assessment must tolerate the new kind without error
    status = session.assess()
    assert status is not None


def test_source_topic_generation():
    device = DeviceConfig("tx-1", station="TX", role="transmitter")
    topic = device.topic("source", survey_id="CS1")
    assert topic == "pycsamt/CS1/TX/tx-1/source"
