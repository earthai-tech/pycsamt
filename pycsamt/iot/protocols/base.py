"""Telemetry transport base classes and the dry-run recorder.

:class:`BaseTelemetryClient` defines the transport interface shared by
every concrete client (MQTT, HTTP, file, serial, websocket): ``connect``,
``disconnect``, ``send``, ``receive``, ``subscribe``, ``listen``,
``flush``, and ``healthcheck``. Subclasses implement the ``_transport_*``
hooks; the base handles validation, dry-run recording, auto-connect, and
error wrapping.

Every client supports ``dry_run=True``, in which packets are recorded to
``client.sent`` and no network/hardware access occurs. This keeps demos,
tests, and documentation fully offline and makes optional dependencies
(paho-mqtt, pyserial, websocket-client) truly optional — they are only
imported when a real connection is opened.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable

from ..core import TelemetryPacket

__all__ = [
    "IoTProtocol",
    "TelemetryAck",
    "TelemetryError",
    "BaseTelemetryClient",
    "TelemetryClient",
]


class IoTProtocol(str, Enum):
    """Supported telemetry protocol identifiers."""

    MQTT = "mqtt"
    HTTP = "http"
    LORA = "lora"
    SERIAL = "serial"
    FILE = "file"
    WEBSOCKET = "websocket"


class TelemetryError(RuntimeError):
    """Raised when a telemetry transport operation fails."""


@dataclass
class TelemetryAck:
    """Acknowledgement returned by a telemetry client."""

    ok: bool
    protocol: str
    packet_id: str
    detail: str = ""


def _coerce_packet(packet: Any) -> TelemetryPacket:
    if isinstance(packet, TelemetryPacket):
        packet.validate()
        return packet
    if isinstance(packet, dict):
        return TelemetryPacket(**packet)
    raise TypeError("packet must be a TelemetryPacket or mapping.")


class BaseTelemetryClient:
    """Common behaviour for telemetry transports.

    Parameters
    ----------
    endpoint : str, optional
        Transport-specific destination (broker URL, file path, URL,
        serial port, ...).
    protocol : IoTProtocol
        The protocol this client implements.
    dry_run : bool
        When ``True``, packets are recorded to :attr:`sent` and no real
        transport is opened.
    **options
        Transport-specific options (credentials, timeouts, TLS, ...).
    """

    protocol: IoTProtocol = IoTProtocol.FILE

    def __init__(
        self,
        endpoint: str | None = None,
        *,
        protocol: IoTProtocol | None = None,
        dry_run: bool = False,
        **options: Any,
    ) -> None:
        if protocol is not None:
            self.protocol = (
                protocol
                if isinstance(protocol, IoTProtocol)
                else IoTProtocol(str(protocol))
            )
        self.endpoint = endpoint
        self.dry_run = bool(dry_run)
        self.options: dict[str, Any] = dict(options)
        self.sent: list[TelemetryPacket] = []
        self.subscriptions: list[str] = []
        self.connected = False
        self._handle: Any = None

    # -- context management ------------------------------------------------
    def __enter__(self) -> BaseTelemetryClient:
        self.connect()
        return self

    def __exit__(self, *exc: Any) -> None:
        self.disconnect()

    # -- lifecycle ---------------------------------------------------------
    def connect(self) -> None:
        """Open the transport (no-op in dry-run mode)."""
        if self.connected:
            return
        if not self.dry_run:
            try:
                self._connect()
            except TelemetryError:
                raise
            except Exception as exc:  # noqa: BLE001 - normalise transport errs
                raise TelemetryError(
                    f"{self.protocol.value} connect failed: {exc}"
                ) from exc
        self.connected = True

    def disconnect(self) -> None:
        """Close the transport (no-op in dry-run mode)."""
        if not self.connected:
            return
        if not self.dry_run:
            try:
                self._disconnect()
            except Exception:  # noqa: BLE001 - disconnect is best effort
                pass
        self.connected = False

    # -- messaging ---------------------------------------------------------
    def send(self, packet: Any) -> TelemetryAck:
        """Send *packet*, recording it in dry-run mode."""
        pkt = _coerce_packet(packet)
        packet_id = self._packet_id(pkt)
        if self.dry_run:
            self.sent.append(pkt)
            return TelemetryAck(
                ok=True,
                protocol=self.protocol.value,
                packet_id=packet_id,
                detail="dry-run packet recorded",
            )
        if not self.connected:
            self.connect()
        try:
            detail = self._transport_send(pkt)
        except (TelemetryError, NotImplementedError):
            # NotImplementedError marks an unimplemented transport (e.g. the
            # generic recorder); keep it distinct from real transport errors.
            raise
        except Exception as exc:  # noqa: BLE001 - normalise transport errs
            raise TelemetryError(f"{self.protocol.value} send failed: {exc}") from exc
        self.sent.append(pkt)
        return TelemetryAck(
            ok=True,
            protocol=self.protocol.value,
            packet_id=packet_id,
            detail=detail or "sent",
        )

    def receive(self, *, timeout: float | None = None) -> dict[str, Any] | None:
        """Receive one payload, or ``None`` in dry-run/unsupported mode."""
        if self.dry_run:
            return None
        if not self.connected:
            self.connect()
        return self._transport_receive(timeout=timeout)

    def subscribe(self, topic: str) -> None:
        """Register interest in *topic* (transport permitting)."""
        topic = str(topic).strip()
        if not topic:
            raise ValueError("topic cannot be empty.")
        if topic not in self.subscriptions:
            self.subscriptions.append(topic)
        if not self.dry_run:
            if not self.connected:
                self.connect()
            self._transport_subscribe(topic)

    def listen(
        self,
        callback: Callable[[dict[str, Any]], None],
        *,
        max_messages: int | None = None,
        timeout: float | None = None,
    ) -> int:
        """Poll :meth:`receive`, invoking *callback* per payload.

        Returns the number of payloads dispatched. In dry-run mode this
        returns ``0`` immediately.
        """
        if self.dry_run:
            return 0
        n = 0
        while max_messages is None or n < max_messages:
            payload = self.receive(timeout=timeout)
            if payload is None:
                break
            callback(payload)
            n += 1
        return n

    def flush(self) -> int:
        """Flush any buffered packets; returns count flushed (default 0)."""
        return 0

    def healthcheck(self) -> bool:
        """Return whether the transport is reachable."""
        if self.dry_run:
            return True
        try:
            self.connect()
            return self._transport_healthcheck()
        except Exception:  # noqa: BLE001 - health probe must not raise
            return False

    # -- helpers -----------------------------------------------------------
    @staticmethod
    def _packet_id(packet: TelemetryPacket) -> str:
        return f"{packet.device_id}:{packet.timestamp:g}:{packet.topic}"

    @staticmethod
    def _payload_bytes(packet: TelemetryPacket) -> bytes:
        import json

        return json.dumps(packet.as_dict(), default=str).encode("utf-8")

    def _require_endpoint(self) -> str:
        if not self.endpoint:
            raise TelemetryError(f"{self.protocol.value} client requires an endpoint.")
        return str(self.endpoint)

    # -- transport hooks (overridden by concrete clients) ------------------
    def _connect(self) -> None:  # pragma: no cover - trivial default
        return None

    def _disconnect(self) -> None:  # pragma: no cover - trivial default
        return None

    def _transport_send(self, packet: TelemetryPacket) -> str:
        raise NotImplementedError(
            f"{self.protocol.value} transport does not implement send()."
        )

    def _transport_receive(
        self, *, timeout: float | None = None
    ) -> dict[str, Any] | None:
        raise TelemetryError(
            f"{self.protocol.value} transport does not support receive()."
        )

    def _transport_subscribe(self, topic: str) -> None:
        raise TelemetryError(
            f"{self.protocol.value} transport does not support subscribe()."
        )

    def _transport_healthcheck(self) -> bool:  # pragma: no cover - default
        return self.connected


class TelemetryClient(BaseTelemetryClient):
    """Generic recorder used for dry-run simulation and unmapped protocols.

    Retains the historical dry-run behaviour: with ``dry_run=True``
    (default) packets are recorded to :attr:`sent`. With ``dry_run=False``
    it has no real transport and raises :class:`NotImplementedError` on
    send, since a concrete client should be used instead.
    """

    def __init__(
        self,
        endpoint: str | None = None,
        *,
        protocol: IoTProtocol | str = IoTProtocol.FILE,
        dry_run: bool = True,
        **options: Any,
    ) -> None:
        super().__init__(endpoint, protocol=protocol, dry_run=dry_run, **options)

    def _transport_send(self, packet: TelemetryPacket) -> str:
        raise NotImplementedError(
            f"{self.protocol.value} transport is not implemented; use a "
            "concrete client such as FileTelemetryClient or "
            "HTTPTelemetryClient."
        )
