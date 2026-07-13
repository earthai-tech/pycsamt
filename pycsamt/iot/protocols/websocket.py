"""WebSocket telemetry transport (optional dependency: ``websocket-client``).

:class:`WebSocketTelemetryClient` sends JSON frames to a ``ws://`` or
``wss://`` endpoint and receives them back. The ``websocket-client``
package is imported lazily.
"""

from __future__ import annotations

import json
from typing import Any

from ..core import TelemetryPacket
from .base import (
    BaseTelemetryClient,
    IoTProtocol,
    TelemetryError,
)

__all__ = ["WebSocketTelemetryClient"]


def _import_websocket():
    try:
        import websocket  # type: ignore  # from the 'websocket-client' package
    except Exception as exc:  # noqa: BLE001
        raise TelemetryError(
            "WebSocket transport requires the 'websocket-client' package. "
            "Install it with `pip install websocket-client`."
        ) from exc
    return websocket


class WebSocketTelemetryClient(BaseTelemetryClient):
    """Send/receive telemetry over a WebSocket connection."""

    protocol = IoTProtocol.WEBSOCKET

    def __init__(
        self,
        endpoint: str | None = None,
        *,
        dry_run: bool = False,
        timeout: float = 10.0,
        **options: Any,
    ) -> None:
        super().__init__(
            endpoint, dry_run=dry_run, timeout=timeout, **options
        )

    def _connect(self) -> None:
        websocket = _import_websocket()
        url = self._require_endpoint()
        if not str(url).lower().startswith(("ws://", "wss://")):
            raise TelemetryError(
                f"WebSocket endpoint must start with ws:// or wss:// (got {url!r})."
            )
        self._handle = websocket.create_connection(
            url, timeout=float(self.options.get("timeout", 10.0))
        )

    def _disconnect(self) -> None:
        if self._handle is not None:
            try:
                self._handle.close()
            finally:
                self._handle = None

    def _transport_send(self, packet: TelemetryPacket) -> str:
        if self._handle is None:
            raise TelemetryError("WebSocket client is not connected.")
        self._handle.send(json.dumps(packet.as_dict(), default=str))
        return f"sent {packet.topic} over websocket"

    def _transport_receive(
        self, *, timeout: float | None = None
    ) -> dict[str, Any] | None:
        if self._handle is None:
            raise TelemetryError("WebSocket client is not connected.")
        try:
            message = self._handle.recv()
        except Exception:  # noqa: BLE001 - treat recv errors as no message
            return None
        if not message:
            return None
        try:
            return json.loads(message)
        except json.JSONDecodeError:
            return {"raw": message}

    def _transport_healthcheck(self) -> bool:
        return self._handle is not None and getattr(
            self._handle, "connected", False
        )
