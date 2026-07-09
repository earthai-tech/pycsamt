"""Telemetry transports for IoT-enabled field acquisition.

This package supersedes the former single-module ``protocols.py``. It
keeps the original public names (:class:`IoTProtocol`,
:class:`TelemetryAck`, :class:`TelemetryClient`,
:func:`build_telemetry_client`) and adds concrete clients:

* :class:`FileTelemetryClient` — newline-delimited JSON (stdlib).
* :class:`HTTPTelemetryClient` — HTTP(S) POST (stdlib ``urllib``).
* :class:`MQTTTelemetryClient` — MQTT (needs ``paho-mqtt``).
* :class:`SerialTelemetryClient` — UART (needs ``pyserial``).
* :class:`WebSocketTelemetryClient` — WebSocket (needs ``websocket-client``).

Every client supports ``dry_run=True`` for offline recording.
"""

from __future__ import annotations

from typing import Any, Optional

from .base import (
    BaseTelemetryClient,
    IoTProtocol,
    TelemetryAck,
    TelemetryClient,
    TelemetryError,
)
from .file import FileTelemetryClient
from .http import HTTPTelemetryClient
from .mqtt import MQTTTelemetryClient
from .serial import SerialTelemetryClient
from .websocket import WebSocketTelemetryClient

__all__ = [
    "IoTProtocol",
    "TelemetryAck",
    "TelemetryError",
    "BaseTelemetryClient",
    "TelemetryClient",
    "FileTelemetryClient",
    "HTTPTelemetryClient",
    "MQTTTelemetryClient",
    "SerialTelemetryClient",
    "WebSocketTelemetryClient",
    "CLIENT_REGISTRY",
    "build_telemetry_client",
]


#: Concrete client for each protocol with a real transport implementation.
CLIENT_REGISTRY: dict[IoTProtocol, type[BaseTelemetryClient]] = {
    IoTProtocol.FILE: FileTelemetryClient,
    IoTProtocol.HTTP: HTTPTelemetryClient,
    IoTProtocol.MQTT: MQTTTelemetryClient,
    IoTProtocol.SERIAL: SerialTelemetryClient,
    IoTProtocol.WEBSOCKET: WebSocketTelemetryClient,
}


def build_telemetry_client(
    protocol: str | IoTProtocol,
    *,
    endpoint: Optional[str] = None,
    dry_run: bool = True,
    **options: Any,
) -> BaseTelemetryClient:
    """Construct a telemetry client for *protocol*.

    Parameters
    ----------
    protocol : str or IoTProtocol
        Transport identifier (``"mqtt"``, ``"http"``, ``"file"``,
        ``"serial"``, ``"websocket"``, or ``"lora"``).
    endpoint : str, optional
        Transport destination (broker URL, file path, HTTP URL, port).
    dry_run : bool
        When ``True`` (default), the returned client records packets
        without opening a real transport — safe for tests and demos.
    **options
        Transport-specific options forwarded to the client.

    Returns
    -------
    BaseTelemetryClient
        A concrete client when the protocol has a real transport, else a
        generic :class:`TelemetryClient` recorder (e.g. for ``lora``).
    """
    proto = (
        protocol if isinstance(protocol, IoTProtocol)
        else IoTProtocol(str(protocol).strip().lower())
    )
    cls = CLIENT_REGISTRY.get(proto)
    if cls is None:
        # No real transport (e.g. LoRa): fall back to the recorder so
        # dry-run simulations still work and real sends fail loudly.
        return TelemetryClient(
            endpoint, protocol=proto, dry_run=dry_run, **options
        )
    return cls(endpoint, dry_run=dry_run, **options)
