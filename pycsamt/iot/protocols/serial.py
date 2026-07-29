"""Serial (UART) telemetry transport (optional dependency: ``pyserial``).

:class:`SerialTelemetryClient` writes newline-delimited JSON to a serial
port and reads it back. ``pyserial`` is imported lazily so the subpackage
imports without it.
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

__all__ = ["SerialTelemetryClient"]


def _import_pyserial():
    try:
        import serial  # type: ignore
    except Exception as exc:  # noqa: BLE001
        raise TelemetryError(
            "Serial transport requires the 'pyserial' package. Install it "
            "with `pip install pyserial`."
        ) from exc
    return serial


class SerialTelemetryClient(BaseTelemetryClient):
    """Send/receive telemetry over a serial port."""

    protocol = IoTProtocol.SERIAL

    def __init__(
        self,
        endpoint: str | None = None,
        *,
        dry_run: bool = False,
        baudrate: int = 115200,
        timeout: float = 1.0,
        **options: Any,
    ) -> None:
        super().__init__(
            endpoint,
            dry_run=dry_run,
            baudrate=baudrate,
            timeout=timeout,
            **options,
        )

    def _connect(self) -> None:
        serial = _import_pyserial()
        port = self._require_endpoint()
        self._handle = serial.Serial(
            port,
            baudrate=int(self.options.get("baudrate", 115200)),
            timeout=float(self.options.get("timeout", 1.0)),
        )

    def _disconnect(self) -> None:
        if self._handle is not None:
            try:
                self._handle.close()
            finally:
                self._handle = None

    def _transport_send(self, packet: TelemetryPacket) -> str:
        if self._handle is None:
            raise TelemetryError("Serial client is not connected.")
        self._handle.write(self._payload_bytes(packet) + b"\n")
        self._handle.flush()
        return f"wrote {packet.topic} to serial"

    def _transport_receive(
        self, *, timeout: float | None = None
    ) -> dict[str, Any] | None:
        if self._handle is None:
            raise TelemetryError("Serial client is not connected.")
        line = self._handle.readline()
        if not line:
            return None
        text = line.decode("utf-8", errors="replace").strip()
        if not text:
            return None
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            return {"raw": text}

    def _transport_healthcheck(self) -> bool:
        return self._handle is not None and getattr(self._handle, "is_open", False)
