"""Tests for :mod:`pycsamt.iot.protocols.websocket` (WebSocket transport).

``websocket-client`` happens to be installed in this environment, but the
real-connection code paths are still exercised against a fake
``websocket`` module installed into ``sys.modules`` (rather than opening
a real socket), following the fake-optional-dependency convention used
in ``pycsamt/backends/tests``. The "package missing" path is simulated
via the documented ``sys.modules[name] = None`` trick, which forces
``import websocket`` to raise ``ImportError`` regardless of what is
actually installed.
"""

from __future__ import annotations

import sys
import types

import pytest

from pycsamt.iot.core import TelemetryPacket
from pycsamt.iot.protocols.base import IoTProtocol, TelemetryError
from pycsamt.iot.protocols.websocket import WebSocketTelemetryClient


class FakeWSConnection:
    """Stand-in for the object returned by ``websocket.create_connection``."""

    def __init__(self, url, timeout=None):
        self.url = url
        self.timeout = timeout
        self.connected = True
        self.sent = []
        self.queued_messages: list[str] = []
        self.raise_on_recv = False

    def send(self, data):
        self.sent.append(data)

    def recv(self):
        if self.raise_on_recv:
            raise OSError("socket closed")
        if self.queued_messages:
            return self.queued_messages.pop(0)
        return ""

    def close(self):
        self.connected = False


class FakeWebSocketModule(types.ModuleType):
    def __init__(self):
        super().__init__("websocket")
        self.raise_on_create = False
        self.created: list[FakeWSConnection] = []

    def create_connection(self, url, timeout=None):
        if self.raise_on_create:
            raise OSError("connection refused")
        conn = FakeWSConnection(url, timeout=timeout)
        self.created.append(conn)
        return conn


@pytest.fixture
def fake_websocket(monkeypatch):
    module = FakeWebSocketModule()
    monkeypatch.setitem(sys.modules, "websocket", module)
    return module


def _packet(topic="t/1", **kw):
    return TelemetryPacket(
        device_id="node-1", timestamp=1.0, topic=topic, payload={"a": 1}, **kw
    )


# ---------------------------------------------------------------------------
# dry-run contract
# ---------------------------------------------------------------------------
def test_dry_run_full_contract():
    client = WebSocketTelemetryClient("ws://example.org/telemetry", dry_run=True)
    assert client.protocol == IoTProtocol.WEBSOCKET

    with client:
        assert client.connected is True
        ack = client.send(_packet())
        assert ack.ok is True
        assert ack.protocol == "websocket"
        assert len(client.sent) == 1

        client.subscribe("ignored")
        assert client.receive() is None
        assert client.listen(lambda p: None) == 0
        assert client.flush() == 0
        assert client.healthcheck() is True
    assert client.connected is False


# ---------------------------------------------------------------------------
# missing dependency (forced via sys.modules[name] = None)
# ---------------------------------------------------------------------------
def test_real_connect_without_websocket_client_raises_telemetry_error(monkeypatch):
    monkeypatch.setitem(sys.modules, "websocket", None)
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    with pytest.raises(TelemetryError, match="websocket-client"):
        client.connect()


# ---------------------------------------------------------------------------
# URL validation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("endpoint", ["http://example.org", "example.org", "ftp://x"])
def test_connect_rejects_non_ws_scheme(endpoint):
    client = WebSocketTelemetryClient(endpoint, dry_run=False)
    with pytest.raises(TelemetryError, match="ws://"):
        client.connect()


def test_connect_accepts_wss_scheme(fake_websocket):
    client = WebSocketTelemetryClient("wss://example.org/telemetry", dry_run=False)
    client.connect()
    assert client.connected is True
    assert fake_websocket.created[0].url == "wss://example.org/telemetry"


# ---------------------------------------------------------------------------
# real connection lifecycle (fake websocket module)
# ---------------------------------------------------------------------------
def test_real_connect_passes_timeout(fake_websocket):
    client = WebSocketTelemetryClient(
        "ws://example.org/telemetry", dry_run=False, timeout=3.5
    )
    client.connect()
    assert fake_websocket.created[0].timeout == 3.5


def test_connect_failure_wrapped_as_telemetry_error(fake_websocket):
    fake_websocket.raise_on_create = True
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    with pytest.raises(TelemetryError, match="websocket connect failed"):
        client.connect()


def test_healthcheck_false_on_connect_failure(fake_websocket):
    fake_websocket.raise_on_create = True
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    assert client.healthcheck() is False


def test_send_serializes_packet_as_json(fake_websocket):
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    client.connect()
    ack = client.send(_packet(topic="pycsamt/node-1/data"))
    assert ack.ok is True
    assert "sent" in ack.detail
    sent_frame = client._handle.sent[0]
    assert "pycsamt/node-1/data" in sent_frame


def test_healthcheck_true_when_connected(fake_websocket):
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    client.connect()
    assert client.healthcheck() is True
    client.disconnect()
    assert client._handle is None


def test_disconnect_is_best_effort_when_close_raises(fake_websocket):
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    client.connect()

    def _boom():
        raise RuntimeError("close exploded")

    client._handle.close = _boom
    client.disconnect()
    assert client.connected is False
    assert client._handle is None


# ---------------------------------------------------------------------------
# receive parsing
# ---------------------------------------------------------------------------
def test_receive_returns_none_on_empty_message(fake_websocket):
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    client.connect()
    assert client.receive() is None


def test_receive_returns_none_when_recv_raises(fake_websocket):
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    client.connect()
    client._handle.raise_on_recv = True
    assert client.receive() is None


def test_receive_parses_json_message(fake_websocket):
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    client.connect()
    client._handle.queued_messages.append('{"topic": "a", "v": 1}')
    assert client.receive() == {"topic": "a", "v": 1}


def test_receive_wraps_non_json_message_as_raw(fake_websocket):
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    client.connect()
    client._handle.queued_messages.append("not-json")
    assert client.receive() == {"raw": "not-json"}


def test_listen_dispatches_queued_messages_then_stops(fake_websocket):
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    client.connect()
    client._handle.queued_messages.extend(['{"v": 1}', '{"v": 2}'])
    seen = []
    n = client.listen(seen.append, max_messages=10)
    assert n == 2
    assert [item["v"] for item in seen] == [1, 2]


# ---------------------------------------------------------------------------
# unreachable-via-public-API branches, exercised directly
# ---------------------------------------------------------------------------
def test_transport_send_requires_handle():
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    with pytest.raises(TelemetryError, match="not connected"):
        client._transport_send(_packet())


def test_transport_receive_requires_handle():
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    with pytest.raises(TelemetryError, match="not connected"):
        client._transport_receive()


def test_transport_healthcheck_false_without_handle():
    client = WebSocketTelemetryClient("ws://example.org", dry_run=False)
    assert client._transport_healthcheck() is False
