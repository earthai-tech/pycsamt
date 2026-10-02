"""Tests for :mod:`pycsamt.iot.protocols.mqtt` (MQTT telemetry transport).

``paho-mqtt`` is not installed in this environment, so the "package
missing" path is exercised for free via the real import failure. The
real-connection code paths (``_connect``/``_transport_send``/...) are
exercised against a fake ``paho.mqtt.client`` module installed into
``sys.modules``, following the fake-optional-dependency convention used
in ``pycsamt/backends/tests``.
"""

from __future__ import annotations

import sys
import types

import pytest

from pycsamt.iot.core import TelemetryPacket
from pycsamt.iot.protocols.base import IoTProtocol, TelemetryError
from pycsamt.iot.protocols.mqtt import MQTTTelemetryClient


class _FakePublishInfo:
    def __init__(self, raise_on_wait=False):
        self._raise = raise_on_wait

    def wait_for_publish(self, timeout=None):
        if self._raise:
            raise RuntimeError("this paho build has no wait_for_publish timeout")


class FakeMQTTClient:
    """Stand-in for ``paho.mqtt.client.Client``."""

    instances: list["FakeMQTTClient"] = []
    raise_on_connect = False
    raise_on_publish = False
    raise_on_wait_for_publish = False

    def __init__(self, client_id=""):
        self.client_id = client_id
        self.on_message = None
        self.connected = False
        self.published = []
        self.subscribed = []
        self.username = None
        self.password = None
        self.tls = False
        self.host = None
        self.port = None
        self.keepalive = None
        self.loop_started = False
        FakeMQTTClient.instances.append(self)

    def username_pw_set(self, username, password=None):
        self.username = username
        self.password = password

    def tls_set(self):
        self.tls = True

    def connect(self, host, port, keepalive=60):
        if FakeMQTTClient.raise_on_connect:
            raise OSError("connection refused")
        self.host, self.port, self.keepalive = host, port, keepalive
        self.connected = True

    def loop_start(self):
        self.loop_started = True

    def loop_stop(self):
        self.loop_started = False

    def disconnect(self):
        self.connected = False

    def publish(self, topic, payload=None, qos=0, retain=False):
        if FakeMQTTClient.raise_on_publish:
            raise RuntimeError("publish failed")
        self.published.append(
            {"topic": topic, "payload": payload, "qos": qos, "retain": retain}
        )
        return _FakePublishInfo(raise_on_wait=FakeMQTTClient.raise_on_wait_for_publish)

    def subscribe(self, topic):
        self.subscribed.append(topic)

    def is_connected(self):
        return self.connected


@pytest.fixture
def fake_paho(monkeypatch):
    FakeMQTTClient.instances = []
    FakeMQTTClient.raise_on_connect = False
    FakeMQTTClient.raise_on_publish = False
    FakeMQTTClient.raise_on_wait_for_publish = False
    paho = types.ModuleType("paho")
    paho_mqtt = types.ModuleType("paho.mqtt")
    paho_mqtt_client = types.ModuleType("paho.mqtt.client")
    paho_mqtt_client.Client = FakeMQTTClient
    paho.mqtt = paho_mqtt
    paho_mqtt.client = paho_mqtt_client
    monkeypatch.setitem(sys.modules, "paho", paho)
    monkeypatch.setitem(sys.modules, "paho.mqtt", paho_mqtt)
    monkeypatch.setitem(sys.modules, "paho.mqtt.client", paho_mqtt_client)
    return paho_mqtt_client


def _packet(topic="t/1", **kw):
    return TelemetryPacket(
        device_id="node-1", timestamp=1.0, topic=topic, payload={"a": 1}, **kw
    )


# ---------------------------------------------------------------------------
# dry-run contract
# ---------------------------------------------------------------------------
def test_dry_run_full_contract():
    client = MQTTTelemetryClient("mqtt://broker.example.com", dry_run=True)
    assert client.protocol == IoTProtocol.MQTT

    with client:
        assert client.connected is True
        ack = client.send(_packet())
        assert ack.ok is True
        assert ack.protocol == "mqtt"
        assert len(client.sent) == 1

        client.subscribe("topic/a")
        assert client.subscriptions == ["topic/a"]

        assert client.receive() is None
        assert client.listen(lambda p: None) == 0
        assert client.flush() == 0
        assert client.healthcheck() is True
    assert client.connected is False


# ---------------------------------------------------------------------------
# missing dependency (paho-mqtt is genuinely not installed here)
# ---------------------------------------------------------------------------
def test_real_connect_without_paho_raises_telemetry_error():
    client = MQTTTelemetryClient("mqtt://broker.example.com", dry_run=False)
    with pytest.raises(TelemetryError, match="paho-mqtt"):
        client.connect()


# ---------------------------------------------------------------------------
# host/port resolution
# ---------------------------------------------------------------------------
def test_resolve_host_port_prefers_explicit_host_option():
    client = MQTTTelemetryClient(dry_run=True, host="10.0.0.5", port=8883)
    assert client._resolve_host_port() == ("10.0.0.5", 8883)


def test_resolve_host_port_parses_endpoint_with_scheme():
    client = MQTTTelemetryClient("mqtt://broker.local:1884", dry_run=True)
    assert client._resolve_host_port() == ("broker.local", 1884)


def test_resolve_host_port_defaults_port_when_missing_from_url():
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=True)
    assert client._resolve_host_port() == ("broker.local", 1883)


def test_resolve_host_port_accepts_bare_host_without_scheme():
    client = MQTTTelemetryClient("broker.local:1885", dry_run=True)
    assert client._resolve_host_port() == ("broker.local", 1885)


def test_resolve_host_port_forces_tls_for_mqtts_scheme():
    client = MQTTTelemetryClient("mqtts://broker.local", dry_run=True)
    host, port = client._resolve_host_port()
    assert host == "broker.local"
    assert client.options["tls"] is True


def test_resolve_host_port_raises_when_hostname_unparseable():
    client = MQTTTelemetryClient("mqtt://:1883", dry_run=True)
    with pytest.raises(TelemetryError, match="Cannot parse MQTT host"):
        client._resolve_host_port()


def test_resolve_host_port_requires_endpoint_when_no_host_option():
    client = MQTTTelemetryClient(dry_run=True)
    with pytest.raises(TelemetryError, match="requires an endpoint"):
        client._resolve_host_port()


# ---------------------------------------------------------------------------
# real connection lifecycle (fake paho.mqtt.client)
# ---------------------------------------------------------------------------
def test_real_connect_configures_client_and_publishes(fake_paho):
    client = MQTTTelemetryClient(
        "mqtt://broker.local:1884",
        dry_run=False,
        username="alice",
        password="secret",
        client_id="node-1",
        keepalive=30,
    )
    client.connect()
    assert client.connected is True
    handle = FakeMQTTClient.instances[-1]
    assert handle.host == "broker.local"
    assert handle.port == 1884
    assert handle.keepalive == 30
    assert handle.username == "alice"
    assert handle.password == "secret"
    assert handle.loop_started is True

    ack = client.send(_packet(topic="pycsamt/node-1/data", qos=1, retained=True))
    assert ack.ok is True
    assert "published to" in ack.detail
    assert handle.published[0]["topic"] == "pycsamt/node-1/data"
    assert handle.published[0]["qos"] == 1
    assert handle.published[0]["retain"] is True

    client.subscribe("pycsamt/+/data")
    assert handle.subscribed == ["pycsamt/+/data"]

    assert client.healthcheck() is True

    client.disconnect()
    assert client.connected is False
    assert handle.connected is False


def test_tls_option_calls_tls_set(fake_paho):
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False, tls=True)
    client.connect()
    assert FakeMQTTClient.instances[-1].tls is True


def test_connect_failure_wrapped_as_telemetry_error(fake_paho):
    FakeMQTTClient.raise_on_connect = True
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    with pytest.raises(TelemetryError, match="mqtt connect failed"):
        client.connect()


def test_healthcheck_false_on_connect_failure(fake_paho):
    FakeMQTTClient.raise_on_connect = True
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    assert client.healthcheck() is False


def test_send_failure_wrapped_as_telemetry_error(fake_paho):
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    client.connect()
    FakeMQTTClient.raise_on_publish = True
    with pytest.raises(TelemetryError, match="mqtt send failed"):
        client.send(_packet())


def test_send_tolerates_missing_wait_for_publish_timeout_support(fake_paho):
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    client.connect()
    FakeMQTTClient.raise_on_wait_for_publish = True
    ack = client.send(_packet())
    assert ack.ok is True


def test_receive_pops_inbox_fifo(fake_paho):
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    client.connect()
    client._inbox.append({"topic": "a", "v": 1})
    client._inbox.append({"topic": "a", "v": 2})
    assert client.receive()["v"] == 1
    assert client.receive()["v"] == 2
    assert client.receive() is None


def test_listen_dispatches_until_inbox_drained(fake_paho):
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    client.connect()
    client._inbox.extend([{"v": 1}, {"v": 2}, {"v": 3}])
    seen = []
    n = client.listen(seen.append, max_messages=10)
    assert n == 3
    assert [item["v"] for item in seen] == [1, 2, 3]


def test_disconnect_is_best_effort_when_handle_misbehaves(fake_paho):
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    client.connect()
    handle = client._handle

    def _boom():
        raise RuntimeError("loop_stop exploded")

    handle.loop_stop = _boom
    client.disconnect()  # base.disconnect() swallows _disconnect() errors
    assert client.connected is False
    assert client._handle is None


# ---------------------------------------------------------------------------
# unreachable-via-public-API branches, exercised directly
# ---------------------------------------------------------------------------
def test_transport_send_requires_handle():
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    with pytest.raises(TelemetryError, match="not connected"):
        client._transport_send(_packet())


def test_transport_subscribe_requires_handle():
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    with pytest.raises(TelemetryError, match="not connected"):
        client._transport_subscribe("t")


def test_transport_healthcheck_false_without_handle():
    client = MQTTTelemetryClient("mqtt://broker.local", dry_run=False)
    assert client._transport_healthcheck() is False
