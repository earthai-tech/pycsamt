"""Store-and-forward buffering for intermittent field links.

Remote AMT/CSAMT nodes routinely lose their uplink -- a gateway drops out,
a cellular backhaul flaps, a boat moves out of range. A bare telemetry
client raises :class:`TelemetryError` on every send in that window and the
packets are lost. :class:`StoreAndForwardClient` wraps any
:class:`~pycsamt.iot.protocols.base.BaseTelemetryClient` and, when a send
fails, queues the packet instead of dropping it. :meth:`flush` later drains
the queue in order, and an optional JSON-lines spool file lets a node
survive a restart with its backlog intact.

The semantics are at-least-once and order-preserving: once anything is
queued, later packets also queue (so nothing overtakes the backlog), and
:meth:`flush` stops at the first failure, leaving the rest for the next
attempt.
"""

from __future__ import annotations

import json
import os
from typing import Any

from ..core import TelemetryPacket
from .base import BaseTelemetryClient, TelemetryAck, TelemetryError, _coerce_packet

__all__ = ["StoreAndForwardClient"]


class StoreAndForwardClient:
    """Wrap a telemetry client with a persistent offline send buffer.

    Parameters
    ----------
    client : BaseTelemetryClient
        The underlying transport used for real sends.
    spool_path : str, optional
        Path to a JSON-lines spool file. When given, the queue is
        persisted on every change and reloaded on construction, so a node
        keeps its backlog across restarts.
    max_queue : int, optional
        Maximum number of buffered packets. When the buffer is full the
        oldest packet is dropped (the newest data is the most valuable for
        a live survey). ``None`` means unbounded.
    base_backoff_s : float
        Base delay for the exponential backoff hint (see
        :attr:`next_retry_delay_s`).
    max_backoff_s : float
        Ceiling for the backoff hint.
    """

    def __init__(
        self,
        client: BaseTelemetryClient,
        *,
        spool_path: str | None = None,
        max_queue: int | None = None,
        base_backoff_s: float = 1.0,
        max_backoff_s: float = 300.0,
    ) -> None:
        if not isinstance(client, BaseTelemetryClient):
            raise TypeError("client must be a BaseTelemetryClient.")
        if max_queue is not None and int(max_queue) <= 0:
            raise ValueError("max_queue must be positive or None.")
        self.client = client
        self.spool_path = spool_path
        self.max_queue = None if max_queue is None else int(max_queue)
        self.base_backoff_s = float(base_backoff_s)
        self.max_backoff_s = float(max_backoff_s)
        self.queue: list[TelemetryPacket] = []
        self.n_dropped = 0
        self._failures = 0
        if spool_path:
            self._load_spool()

    # -- context management ------------------------------------------------
    def __enter__(self) -> StoreAndForwardClient:
        self.client.connect()
        return self

    def __exit__(self, *exc: Any) -> None:
        self.client.disconnect()

    # -- properties --------------------------------------------------------
    @property
    def pending(self) -> int:
        """Number of packets currently buffered."""
        return len(self.queue)

    @property
    def next_retry_delay_s(self) -> float:
        """Exponential-backoff hint for when to next call :meth:`flush`.

        Returns ``0.0`` when the queue is empty. After consecutive flush
        failures the delay grows as ``base * 2**(failures-1)``, capped at
        ``max_backoff_s`` -- a caller's scheduling loop can honour this
        without this class ever sleeping itself.
        """
        if not self.queue:
            return 0.0
        if self._failures <= 0:
            return self.base_backoff_s
        delay = self.base_backoff_s * (2.0 ** (self._failures - 1))
        return min(delay, self.max_backoff_s)

    # -- messaging ---------------------------------------------------------
    def send(self, packet: Any) -> TelemetryAck:
        """Send *packet*, or queue it when the transport is unavailable.

        A direct send is attempted only when the buffer is empty, so a
        backlog is never overtaken by fresh packets. On transport failure
        the packet is queued and a ``queued`` acknowledgement is returned
        rather than raising.
        """
        pkt = _coerce_packet(packet)
        if not self.queue:
            try:
                return self.client.send(pkt)
            except TelemetryError as exc:
                self._enqueue(pkt)
                return self._queued_ack(pkt, f"queued: {exc}")
        self._enqueue(pkt)
        return self._queued_ack(pkt, "queued behind backlog")

    def flush(self) -> int:
        """Drain the buffer in order; return the number of packets sent.

        Sending stops at the first failure, leaving the remaining packets
        queued for a later attempt (at-least-once, order-preserving).
        """
        sent = 0
        while self.queue:
            pkt = self.queue[0]
            try:
                self.client.send(pkt)
            except TelemetryError:
                self._failures += 1
                break
            self.queue.pop(0)
            sent += 1
        else:
            # queue fully drained
            self._failures = 0
        if sent:
            self._persist()
        return sent

    # -- helpers -----------------------------------------------------------
    def _queued_ack(self, packet: TelemetryPacket, detail: str) -> TelemetryAck:
        return TelemetryAck(
            ok=False,
            protocol=self.client.protocol.value,
            packet_id=BaseTelemetryClient._packet_id(packet),
            detail=detail,
        )

    def _enqueue(self, packet: TelemetryPacket) -> None:
        self.queue.append(packet)
        if self.max_queue is not None and len(self.queue) > self.max_queue:
            self.queue.pop(0)  # drop the oldest
            self.n_dropped += 1
        self._persist()

    def _persist(self) -> None:
        if not self.spool_path:
            return
        parent = os.path.dirname(os.path.abspath(self.spool_path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        tmp = f"{self.spool_path}.tmp"
        with open(tmp, "w", encoding="utf-8") as handle:
            for pkt in self.queue:
                handle.write(json.dumps(pkt.as_dict(), default=str) + "\n")
        os.replace(tmp, self.spool_path)  # atomic swap

    def _load_spool(self) -> None:
        if not self.spool_path or not os.path.isfile(self.spool_path):
            return
        restored: list[TelemetryPacket] = []
        with open(self.spool_path, encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    restored.append(TelemetryPacket(**json.loads(line)))
                except Exception:  # noqa: BLE001 - skip a corrupt spool line
                    continue
        self.queue = restored
