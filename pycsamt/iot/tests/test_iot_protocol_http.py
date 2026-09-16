"""Tests for :mod:`pycsamt.iot.protocols.http` (HTTP(S) telemetry transport).

This transport only relies on the standard library, so ``urllib.request``
is mocked directly at the module boundary rather than faking a
third-party package.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request

import pytest

from pycsamt.iot.core import TelemetryPacket
from pycsamt.iot.protocols.base import IoTProtocol, TelemetryError
from pycsamt.iot.protocols.http import HTTPTelemetryClient


class _FakeResponse:
    """Response with a ``status`` attribute (modern ``http.client``).

    Real ``http.client`` responses always implement both ``.status`` and
    ``.getcode()`` -- ``getattr(obj, "status", obj.getcode())`` evaluates
    the default eagerly, so ``getcode()`` must exist even though this
    branch is expected to read ``.status``.
    """

    def __init__(self, status):
        self.status = status

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def getcode(self):
        return self.status


class _FakeResponseLegacy:
    """Response without a ``status`` attribute, only ``getcode()``."""

    def __init__(self, status):
        self._status = status

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def getcode(self):
        return self._status


def _packet(topic="t/1", **kw):
    return TelemetryPacket(
        device_id="node-1", timestamp=1.0, topic=topic, payload={"a": 1}, **kw
    )


# ---------------------------------------------------------------------------
# dry-run contract
# ---------------------------------------------------------------------------
def test_dry_run_full_contract():
    client = HTTPTelemetryClient("http://example.org/telemetry", dry_run=True)
    assert client.protocol == IoTProtocol.HTTP

    with client:
        assert client.connected is True
        ack = client.send(_packet())
        assert ack.ok is True
        assert ack.protocol == "http"
        assert len(client.sent) == 1

        client.subscribe("ignored")
        assert client.receive() is None
        assert client.listen(lambda p: None) == 0
        assert client.flush() == 0
        assert client.healthcheck() is True
    assert client.connected is False


# ---------------------------------------------------------------------------
# endpoint / config validation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("endpoint", ["ftp://example.org", "example.org", "ws://x"])
def test_connect_rejects_non_http_scheme(endpoint):
    client = HTTPTelemetryClient(endpoint, dry_run=False)
    with pytest.raises(TelemetryError, match="http:// or https://"):
        client.connect()


def test_connect_requires_endpoint():
    client = HTTPTelemetryClient(dry_run=False)
    with pytest.raises(TelemetryError, match="requires an endpoint"):
        client.connect()


def test_connect_accepts_https_and_sets_connected():
    client = HTTPTelemetryClient("https://example.org/telemetry", dry_run=False)
    client.connect()
    assert client.connected is True


def test_healthcheck_false_without_endpoint():
    client = HTTPTelemetryClient(dry_run=False)
    assert client.healthcheck() is False


def test_healthcheck_true_with_valid_endpoint():
    client = HTTPTelemetryClient("http://example.org", dry_run=False)
    assert client.healthcheck() is True


# ---------------------------------------------------------------------------
# header construction
# ---------------------------------------------------------------------------
def test_headers_default_content_type_only():
    client = HTTPTelemetryClient("http://example.org", dry_run=False)
    assert client._headers() == {"Content-Type": "application/json"}


def test_headers_include_bearer_token():
    client = HTTPTelemetryClient(
        "http://example.org", dry_run=False, token="tok-123"
    )
    headers = client._headers()
    assert headers["Authorization"] == "Bearer tok-123"


def test_headers_merge_extra_headers():
    client = HTTPTelemetryClient(
        "http://example.org",
        dry_run=False,
        headers={"X-Api-Key": "k1"},
        token="tok",
    )
    headers = client._headers()
    assert headers["X-Api-Key"] == "k1"
    assert headers["Authorization"] == "Bearer tok"


def test_headers_extra_header_can_override_authorization():
    client = HTTPTelemetryClient(
        "http://example.org",
        dry_run=False,
        headers={"Authorization": "Basic xyz"},
        token="tok",
    )
    headers = client._headers()
    # setdefault: an explicit Authorization header wins over the token.
    assert headers["Authorization"] == "Basic xyz"


# ---------------------------------------------------------------------------
# real send (mocked urllib.request.urlopen)
# ---------------------------------------------------------------------------
def test_send_success_uses_status_attribute(monkeypatch):
    captured = {}

    def fake_urlopen(request, timeout=None):
        captured["request"] = request
        captured["timeout"] = timeout
        return _FakeResponse(status=200)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    client = HTTPTelemetryClient(
        "http://example.org/telemetry", dry_run=False, timeout=5.0, token="tok"
    )
    ack = client.send(_packet(topic="pycsamt/node-1/data"))
    assert ack.ok is True
    assert "200" in ack.detail

    request = captured["request"]
    assert request.full_url == "http://example.org/telemetry"
    assert request.get_method() == "POST"
    assert request.get_header("Authorization") == "Bearer tok"
    body = json.loads(request.data.decode("utf-8"))
    assert body["topic"] == "pycsamt/node-1/data"
    assert captured["timeout"] == 5.0


def test_send_success_falls_back_to_getcode(monkeypatch):
    monkeypatch.setattr(
        urllib.request, "urlopen", lambda request, timeout=None: _FakeResponseLegacy(201)
    )
    client = HTTPTelemetryClient("http://example.org", dry_run=False)
    ack = client.send(_packet())
    assert "201" in ack.detail


def test_send_uses_configured_method(monkeypatch):
    captured = {}

    def fake_urlopen(request, timeout=None):
        captured["method"] = request.get_method()
        return _FakeResponse(status=200)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    client = HTTPTelemetryClient("http://example.org", dry_run=False, method="PUT")
    client.send(_packet())
    assert captured["method"] == "PUT"


def test_send_raises_telemetry_error_on_http_error(monkeypatch):
    def fake_urlopen(request, timeout=None):
        raise urllib.error.HTTPError(
            request.full_url, 404, "Not Found", None, None
        )

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    client = HTTPTelemetryClient("http://example.org", dry_run=False)
    with pytest.raises(TelemetryError, match="HTTP 404"):
        client.send(_packet())


def test_send_raises_telemetry_error_on_url_error(monkeypatch):
    def fake_urlopen(request, timeout=None):
        raise urllib.error.URLError("connection refused")

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    client = HTTPTelemetryClient("http://example.org", dry_run=False)
    with pytest.raises(TelemetryError, match="HTTP request to"):
        client.send(_packet())


def test_send_raises_on_non_2xx_status_without_exception(monkeypatch):
    monkeypatch.setattr(
        urllib.request, "urlopen", lambda request, timeout=None: _FakeResponse(304)
    )
    client = HTTPTelemetryClient("http://example.org", dry_run=False)
    with pytest.raises(TelemetryError, match="status 304"):
        client.send(_packet())


def test_receive_is_unsupported_by_default():
    client = HTTPTelemetryClient("http://example.org", dry_run=False)
    client.connect()
    with pytest.raises(TelemetryError, match="does not support receive"):
        client.receive()
