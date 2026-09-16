"""Tests for :mod:`pycsamt.iot.protocols.serial` (UART telemetry transport).

``pyserial`` is not installed in this environment, so the "package
missing" path is exercised for free via the real import failure. The
real-connection code paths are exercised against a fake ``serial`` module
installed into ``sys.modules``, following the fake-optional-dependency
convention used in ``pycsamt/backends/tests``.
"""

from __future__ import annotations

import sys
import types

import pytest

from pycsamt.iot.core import TelemetryPacket
from pycsamt.iot.protocols.base import IoTProtocol, TelemetryError
from pycsamt.iot.protocols.serial import SerialTelemetryClient


class FakeSerialPort:
    """Stand-in for ``serial.Serial``."""

    fail_open_message = None

    def __init__(self, port, baudrate=115200, timeout=1.0):
        if FakeSerialPort.fail_open_message:
            raise OSError(FakeSerialPort.fail_open_message)
        self.port = port
        self.baudrate = baudrate
        self.timeout = timeout
        self.is_open = True
        self.written = []
        self.flush_count = 0
        self.queued_lines: list[bytes] = []

    def write(self, data):
        self.written.append(data)

    def flush(self):
        self.flush_count += 1

    def close(self):
        self.is_open = False

    def readline(self):
        if self.queued_lines:
            return self.queued_lines.pop(0)
        return b""


@pytest.fixture
def fake_serial(monkeypatch):
    FakeSerialPort.fail_open_message = None
    module = types.ModuleType("serial")
    module.Serial = FakeSerialPort
    monkeypatch.setitem(sys.modules, "serial", module)
    return module


def _packet(topic="t/1", **kw):
    return TelemetryPacket(
        device_id="node-1", timestamp=1.0, topic=topic, payload={"a": 1}, **kw
    )


# ---------------------------------------------------------------------------
# dry-run contract
# ---------------------------------------------------------------------------
def test_dry_run_full_contract():
    client = SerialTelemetryClient("COM3", dry_run=True)
    assert client.protocol == IoTProtocol.SERIAL

    with client:
        assert client.connected is True
        ack = client.send(_packet())
        assert ack.ok is True
        assert ack.protocol == "serial"
        assert len(client.sent) == 1

        client.subscribe("ignored")
        assert client.receive() is None
        assert client.listen(lambda p: None) == 0
        assert client.flush() == 0
        assert client.healthcheck() is True
    assert client.connected is False


# ---------------------------------------------------------------------------
# missing dependency (pyserial is genuinely not installed here)
# ---------------------------------------------------------------------------
def test_real_connect_without_pyserial_raises_telemetry_error():
    client = SerialTelemetryClient("COM3", dry_run=False)
    with pytest.raises(TelemetryError, match="pyserial"):
        client.connect()


def test_real_connect_requires_endpoint(fake_serial):
    client = SerialTelemetryClient(dry_run=False)
    with pytest.raises(TelemetryError, match="requires an endpoint"):
        client.connect()


# ---------------------------------------------------------------------------
# real connection lifecycle (fake serial module)
# ---------------------------------------------------------------------------
def test_real_connect_opens_port_with_options(fake_serial):
    client = SerialTelemetryClient("/dev/ttyUSB0", dry_run=False, baudrate=9600, timeout=2.5)
    client.connect()
    assert client.connected is True
    assert client._handle.port == "/dev/ttyUSB0"
    assert client._handle.baudrate == 9600
    assert client._handle.timeout == 2.5


def test_connect_failure_wrapped_as_telemetry_error(fake_serial):
    FakeSerialPort.fail_open_message = "could not open port"
    client = SerialTelemetryClient("COM3", dry_run=False)
    with pytest.raises(TelemetryError, match="serial connect failed"):
        client.connect()


def test_healthcheck_false_on_connect_failure(fake_serial):
    FakeSerialPort.fail_open_message = "could not open port"
    client = SerialTelemetryClient("COM3", dry_run=False)
    assert client.healthcheck() is False


def test_send_writes_newline_delimited_json(fake_serial):
    client = SerialTelemetryClient("COM3", dry_run=False)
    client.connect()
    ack = client.send(_packet(topic="pycsamt/node-1/data"))
    assert ack.ok is True
    assert "wrote" in ack.detail
    written = client._handle.written[0]
    assert written.endswith(b"\n")
    assert b"pycsamt/node-1/data" in written
    assert client._handle.flush_count == 1


def test_healthcheck_true_when_open(fake_serial):
    client = SerialTelemetryClient("COM3", dry_run=False)
    client.connect()
    assert client.healthcheck() is True
    client.disconnect()
    assert client._handle is None


# ---------------------------------------------------------------------------
# receive parsing
# ---------------------------------------------------------------------------
def test_receive_returns_none_on_timeout_empty_line(fake_serial):
    client = SerialTelemetryClient("COM3", dry_run=False)
    client.connect()
    assert client.receive() is None


def test_receive_returns_none_on_whitespace_only_line(fake_serial):
    client = SerialTelemetryClient("COM3", dry_run=False)
    client.connect()
    client._handle.queued_lines.append(b"   \n")
    assert client.receive() is None


def test_receive_parses_json_line(fake_serial):
    client = SerialTelemetryClient("COM3", dry_run=False)
    client.connect()
    client._handle.queued_lines.append(b'{"topic": "a", "v": 1}\n')
    payload = client.receive()
    assert payload == {"topic": "a", "v": 1}


def test_receive_wraps_non_json_line_as_raw(fake_serial):
    client = SerialTelemetryClient("COM3", dry_run=False)
    client.connect()
    client._handle.queued_lines.append(b"not-json-garbage\n")
    payload = client.receive()
    assert payload == {"raw": "not-json-garbage"}


def test_listen_dispatches_queued_lines_then_stops(fake_serial):
    client = SerialTelemetryClient("COM3", dry_run=False)
    client.connect()
    client._handle.queued_lines.extend(
        [b'{"v": 1}\n', b'{"v": 2}\n']
    )
    seen = []
    n = client.listen(seen.append, max_messages=10)
    assert n == 2
    assert [item["v"] for item in seen] == [1, 2]


def test_disconnect_is_best_effort_when_close_raises(fake_serial):
    client = SerialTelemetryClient("COM3", dry_run=False)
    client.connect()

    def _boom():
        raise RuntimeError("close exploded")

    client._handle.close = _boom
    client.disconnect()
    assert client.connected is False
    assert client._handle is None


# ---------------------------------------------------------------------------
# unreachable-via-public-API branches, exercised directly
# ---------------------------------------------------------------------------
def test_transport_send_requires_handle():
    client = SerialTelemetryClient("COM3", dry_run=False)
    with pytest.raises(TelemetryError, match="not connected"):
        client._transport_send(_packet())


def test_transport_receive_requires_handle():
    client = SerialTelemetryClient("COM3", dry_run=False)
    with pytest.raises(TelemetryError, match="not connected"):
        client._transport_receive()


def test_transport_healthcheck_false_without_handle():
    client = SerialTelemetryClient("COM3", dry_run=False)
    assert client._transport_healthcheck() is False
