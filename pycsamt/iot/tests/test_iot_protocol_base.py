"""Tests for :mod:`pycsamt.iot.protocols.base` (BaseTelemetryClient).

These exercise the transport-agnostic base behaviour directly: dry-run
recording, auto-connect, error wrapping, listen/subscribe bookkeeping,
and the default (un-overridden) transport hooks. ``FileTelemetryClient``
(dependency-free) doubles as the "real" concrete client where a live
connection is needed, following the convention used by the sibling
protocol test files (``dry_run=True`` for the offline contract, a real
concrete client for lifecycle/error-wrapping behaviour).
"""

from __future__ import annotations

import os

import pytest

from pycsamt.iot.core import TelemetryPacket
from pycsamt.iot.protocols.base import (
    BaseTelemetryClient,
    IoTProtocol,
    TelemetryAck,
    TelemetryClient,
    TelemetryError,
    _coerce_packet,
)
from pycsamt.iot.protocols.file import FileTelemetryClient


def _packet(topic="t/1", **kw):
    return TelemetryPacket(
        device_id="node-1", timestamp=1.0, topic=topic, payload={"a": 1}, **kw
    )


# ---------------------------------------------------------------------------
# _coerce_packet
# ---------------------------------------------------------------------------
def test_coerce_packet_accepts_packet_dict_and_rejects_other():
    pkt = _packet()
    assert _coerce_packet(pkt) is pkt

    from_dict = _coerce_packet(
        dict(device_id="n", timestamp=1.0, topic="t", payload={})
    )
    assert isinstance(from_dict, TelemetryPacket)

    with pytest.raises(TypeError):
        _coerce_packet(object())


# ---------------------------------------------------------------------------
# construction / protocol coercion
# ---------------------------------------------------------------------------
def test_protocol_defaults_to_file_and_accepts_string():
    c = BaseTelemetryClient()
    assert c.protocol == IoTProtocol.FILE

    c2 = BaseTelemetryClient(protocol="mqtt")
    assert c2.protocol == IoTProtocol.MQTT

    c3 = BaseTelemetryClient(protocol=IoTProtocol.HTTP)
    assert c3.protocol == IoTProtocol.HTTP


def test_options_are_captured():
    c = BaseTelemetryClient("ep", timeout=5, retries=2)
    assert c.options == {"timeout": 5, "retries": 2}
    assert c.sent == []
    assert c.subscriptions == []
    assert c.connected is False


# ---------------------------------------------------------------------------
# dry-run contract on the abstract base itself
# ---------------------------------------------------------------------------
def test_dry_run_send_records_and_does_not_touch_transport():
    c = BaseTelemetryClient(dry_run=True)
    ack = c.send(_packet())
    assert isinstance(ack, TelemetryAck)
    assert ack.ok is True
    assert ack.detail == "dry-run packet recorded"
    assert len(c.sent) == 1
    assert c.receive() is None
    assert c.listen(lambda p: None) == 0


def test_context_manager_connects_and_disconnects():
    c = BaseTelemetryClient(dry_run=True)
    with c as ctx:
        assert ctx is c
        assert c.connected is True
    assert c.connected is False


def test_connect_is_idempotent():
    c = BaseTelemetryClient(dry_run=True)
    c.connect()
    c.connect()  # no-op second time (already connected branch)
    assert c.connected is True


def test_disconnect_when_not_connected_is_noop():
    c = BaseTelemetryClient(dry_run=True)
    c.disconnect()
    assert c.connected is False


def test_healthcheck_true_in_dry_run():
    c = BaseTelemetryClient(dry_run=True)
    assert c.healthcheck() is True


def test_flush_default_returns_zero():
    c = BaseTelemetryClient(dry_run=True)
    assert c.flush() == 0


# ---------------------------------------------------------------------------
# default (un-overridden) transport hooks raise / return safe defaults
# ---------------------------------------------------------------------------
def test_transport_send_not_implemented_by_default():
    c = BaseTelemetryClient(dry_run=False)
    with pytest.raises(NotImplementedError):
        c._transport_send(_packet())


def test_transport_receive_raises_telemetry_error_by_default():
    c = BaseTelemetryClient(dry_run=False)
    with pytest.raises(TelemetryError, match="does not support receive"):
        c._transport_receive()


def test_transport_subscribe_raises_telemetry_error_by_default():
    c = BaseTelemetryClient(dry_run=False)
    with pytest.raises(TelemetryError, match="does not support subscribe"):
        c._transport_subscribe("topic")


def test_require_endpoint_raises_when_missing():
    c = BaseTelemetryClient(dry_run=False)
    with pytest.raises(TelemetryError, match="requires an endpoint"):
        c._require_endpoint()


def test_packet_id_and_payload_bytes_helpers():
    pkt = _packet(topic="pycsamt/x")
    pid = BaseTelemetryClient._packet_id(pkt)
    assert pkt.device_id in pid
    assert "pycsamt/x" in pid

    raw = BaseTelemetryClient._payload_bytes(pkt)
    assert isinstance(raw, bytes)
    assert b"pycsamt/x" in raw


# ---------------------------------------------------------------------------
# send(): real transport path, error wrapping, NotImplementedError passthrough
# ---------------------------------------------------------------------------
def test_send_real_transport_raises_notimplementederror_on_generic_client(tmp_path):
    client = TelemetryClient(str(tmp_path / "x.jsonl"), dry_run=False)
    with pytest.raises(NotImplementedError):
        client.send(_packet())


def test_send_auto_connects_when_not_connected(tmp_path):
    path = str(tmp_path / "out.jsonl")
    client = FileTelemetryClient(path, dry_run=False)
    assert client.connected is False
    ack = client.send(_packet())
    assert client.connected is True
    assert ack.ok is True
    assert os.path.exists(path)


def test_send_wraps_unexpected_exception_as_telemetry_error(monkeypatch, tmp_path):
    client = FileTelemetryClient(str(tmp_path / "out.jsonl"), dry_run=False)
    client.connect()

    def _boom(packet):
        raise RuntimeError("disk full")

    monkeypatch.setattr(client, "_transport_send", _boom)
    with pytest.raises(TelemetryError, match="send failed"):
        client.send(_packet())


def test_connect_wraps_unexpected_exception_as_telemetry_error(monkeypatch):
    client = BaseTelemetryClient(dry_run=False)

    def _boom():
        raise RuntimeError("cable unplugged")

    monkeypatch.setattr(client, "_connect", _boom)
    with pytest.raises(TelemetryError, match="connect failed"):
        client.connect()
    assert client.connected is False


def test_connect_reraises_telemetry_error_unwrapped(monkeypatch):
    client = BaseTelemetryClient(dry_run=False)

    def _boom():
        raise TelemetryError("already the right type")

    monkeypatch.setattr(client, "_connect", _boom)
    with pytest.raises(TelemetryError, match="already the right type"):
        client.connect()


def test_disconnect_is_best_effort_on_real_transport(monkeypatch):
    client = BaseTelemetryClient(dry_run=False)
    client.connect()

    def _boom():
        raise RuntimeError("already gone")

    monkeypatch.setattr(client, "_disconnect", _boom)
    client.disconnect()
    assert client.connected is False


def test_healthcheck_false_on_connect_failure(monkeypatch):
    client = BaseTelemetryClient(dry_run=False)

    def _boom():
        raise RuntimeError("unreachable")

    monkeypatch.setattr(client, "_connect", _boom)
    assert client.healthcheck() is False


# ---------------------------------------------------------------------------
# subscribe()
# ---------------------------------------------------------------------------
def test_subscribe_rejects_empty_topic():
    client = BaseTelemetryClient(dry_run=True)
    with pytest.raises(ValueError, match="cannot be empty"):
        client.subscribe("   ")


def test_subscribe_dry_run_records_without_transport_call():
    client = BaseTelemetryClient(dry_run=True)
    client.subscribe("a/b")
    client.subscribe("a/b")  # duplicate is not appended twice
    assert client.subscriptions == ["a/b"]


def test_subscribe_real_transport_auto_connects_and_delegates(monkeypatch):
    client = BaseTelemetryClient(dry_run=False)
    seen = []
    monkeypatch.setattr(client, "_transport_subscribe", seen.append)
    client.subscribe("topic/x")
    assert client.connected is True
    assert seen == ["topic/x"]
    assert client.subscriptions == ["topic/x"]


# ---------------------------------------------------------------------------
# receive() / listen()
# ---------------------------------------------------------------------------
def test_receive_auto_connects_on_real_transport(monkeypatch):
    client = BaseTelemetryClient(dry_run=False)
    monkeypatch.setattr(client, "_transport_receive", lambda timeout=None: {"v": 1})
    out = client.receive()
    assert client.connected is True
    assert out == {"v": 1}


def test_listen_dispatches_until_none_and_respects_max_messages(monkeypatch):
    client = BaseTelemetryClient(dry_run=False)
    client.connected = True
    queue = [{"v": 1}, {"v": 2}, {"v": 3}]

    def _recv(timeout=None):
        return queue.pop(0) if queue else None

    monkeypatch.setattr(client, "_transport_receive", _recv)
    seen = []
    n = client.listen(seen.append)
    assert n == 3
    assert [s["v"] for s in seen] == [1, 2, 3]


def test_listen_stops_at_max_messages(monkeypatch):
    client = BaseTelemetryClient(dry_run=False)
    client.connected = True
    monkeypatch.setattr(client, "_transport_receive", lambda timeout=None: {"v": 1})
    seen = []
    n = client.listen(seen.append, max_messages=2)
    assert n == 2
    assert len(seen) == 2


# ---------------------------------------------------------------------------
# TelemetryClient (generic dry-run recorder)
# ---------------------------------------------------------------------------
def test_telemetry_client_defaults_dry_run_true():
    client = TelemetryClient("anything")
    assert client.dry_run is True
    ack = client.send(_packet())
    assert ack.ok is True
    assert len(client.sent) == 1


def test_telemetry_client_accepts_string_protocol():
    client = TelemetryClient(protocol="serial", dry_run=True)
    assert client.protocol == IoTProtocol.SERIAL


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
