"""MQTT telemetry transport (optional dependency: ``paho-mqtt``).

:class:`MQTTTelemetryClient` publishes packets to their topic and can
subscribe/listen for downlink messages. ``paho-mqtt`` is imported lazily
inside :meth:`_connect`, so the subpackage — and dry-run simulations —
work without it installed.
"""

from __future__ import annotations

import json
from typing import Any, Dict, Optional
from urllib.parse import urlparse

from ..core import TelemetryPacket
from .base import BaseTelemetryClient, IoTProtocol, TelemetryError

__all__ = ["MQTTTelemetryClient"]


def _import_paho():
    try:
        import paho.mqtt.client as mqtt  # type: ignore
    except Exception as exc:  # noqa: BLE001
        raise TelemetryError(
            "MQTT transport requires the 'paho-mqtt' package. Install it "
            "with `pip install paho-mqtt`."
        ) from exc
    return mqtt


class MQTTTelemetryClient(BaseTelemetryClient):
    """Publish/subscribe telemetry over MQTT."""

    protocol = IoTProtocol.MQTT

    def __init__(
        self,
        endpoint: Optional[str] = None,
        *,
        dry_run: bool = False,
        host: Optional[str] = None,
        port: int = 1883,
        username: Optional[str] = None,
        password: Optional[str] = None,
        tls: bool = False,
        keepalive: int = 60,
        client_id: Optional[str] = None,
        **options: Any,
    ) -> None:
        super().__init__(
            endpoint,
            dry_run=dry_run,
            host=host,
            port=port,
            username=username,
            password=password,
            tls=tls,
            keepalive=keepalive,
            client_id=client_id,
            **options,
        )
        self._inbox: list[Dict[str, Any]] = []

    def _resolve_host_port(self) -> tuple[str, int]:
        host = self.options.get("host")
        port = int(self.options.get("port", 1883))
        if host:
            return str(host), port
        endpoint = self._require_endpoint()
        parsed = urlparse(
            endpoint if "://" in endpoint else f"mqtt://{endpoint}"
        )
        if not parsed.hostname:
            raise TelemetryError(f"Cannot parse MQTT host from {endpoint!r}.")
        if parsed.scheme in {"mqtts", "ssl"}:
            self.options["tls"] = True
        return parsed.hostname, int(parsed.port or port)

    def _connect(self) -> None:
        mqtt = _import_paho()
        host, port = self._resolve_host_port()
        client = mqtt.Client(client_id=self.options.get("client_id") or "")
        username = self.options.get("username")
        if username:
            client.username_pw_set(username, self.options.get("password"))
        if self.options.get("tls"):
            client.tls_set()

        def _on_message(_client, _userdata, message):  # pragma: no cover
            try:
                payload = json.loads(message.payload.decode("utf-8"))
            except Exception:  # noqa: BLE001
                payload = {"raw": message.payload}
            payload.setdefault("topic", message.topic)
            self._inbox.append(payload)

        client.on_message = _on_message
        client.connect(host, port, keepalive=int(self.options.get("keepalive", 60)))
        client.loop_start()
        self._handle = client

    def _disconnect(self) -> None:
        if self._handle is not None:
            try:
                self._handle.loop_stop()
                self._handle.disconnect()
            finally:
                self._handle = None

    def _transport_send(self, packet: TelemetryPacket) -> str:
        if self._handle is None:
            raise TelemetryError("MQTT client is not connected.")
        info = self._handle.publish(
            packet.topic,
            payload=self._payload_bytes(packet),
            qos=int(packet.qos),
            retain=bool(packet.retained),
        )
        try:
            info.wait_for_publish(timeout=float(self.options.get("timeout", 10.0)))
        except Exception:  # noqa: BLE001 - older paho lacks timeout arg
            pass
        return f"published to {packet.topic}"

    def _transport_subscribe(self, topic: str) -> None:
        if self._handle is None:
            raise TelemetryError("MQTT client is not connected.")
        self._handle.subscribe(topic)

    def _transport_receive(
        self, *, timeout: Optional[float] = None
    ) -> Optional[Dict[str, Any]]:
        if self._inbox:
            return self._inbox.pop(0)
        return None

    def _transport_healthcheck(self) -> bool:
        return self._handle is not None and self._handle.is_connected()
