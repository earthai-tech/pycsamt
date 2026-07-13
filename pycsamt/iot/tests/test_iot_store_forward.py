"""Tests for :class:`pycsamt.iot.StoreAndForwardClient` -- offline
buffering, in-order at-least-once delivery, backoff, and spool persistence.
"""

from __future__ import annotations

import pytest

from pycsamt.iot import StoreAndForwardClient, TelemetryError
from pycsamt.iot.protocols.base import (
    BaseTelemetryClient,
    IoTProtocol,
    TelemetryAck,
    _coerce_packet,
)


class FlakyClient(BaseTelemetryClient):
    """A test transport that can be toggled offline to raise on send."""

    protocol = IoTProtocol.HTTP

    def __init__(self):
        super().__init__("http://x", dry_run=False)
        self.online = True
        self.delivered = []

    def send(self, packet):
        pkt = _coerce_packet(packet)
        if not self.online:
            raise TelemetryError("offline")
        self.delivered.append(pkt)
        return TelemetryAck(ok=True, protocol="http", packet_id="x")


def _pkt(i):
    return {
        "device_id": f"n{i}",
        "timestamp": float(i),
        "topic": "t",
        "payload": {"i": i},
    }


def _order(client):
    return [p.payload["i"] for p in client.delivered]


# ---------------------------------------------------------------------------
# online / offline behaviour
# ---------------------------------------------------------------------------
def test_online_send_is_direct():
    inner = FlakyClient()
    saf = StoreAndForwardClient(inner)
    ack = saf.send(_pkt(1))
    assert ack.ok
    assert saf.pending == 0
    assert _order(inner) == [1]


def test_offline_send_queues_without_raising():
    inner = FlakyClient()
    inner.online = False
    saf = StoreAndForwardClient(inner)
    ack = saf.send(_pkt(1))
    assert not ack.ok  # queued, not delivered
    assert "queued" in ack.detail
    assert saf.pending == 1
    assert _order(inner) == []


def test_flush_drains_in_order_when_online():
    inner = FlakyClient()
    saf = StoreAndForwardClient(inner)
    inner.online = False
    for i in range(1, 4):
        saf.send(_pkt(i))
    assert saf.pending == 3
    inner.online = True
    assert saf.flush() == 3
    assert saf.pending == 0
    assert _order(inner) == [1, 2, 3]


def test_flush_offline_sends_nothing():
    inner = FlakyClient()
    inner.online = False
    saf = StoreAndForwardClient(inner)
    saf.send(_pkt(1))
    assert saf.flush() == 0
    assert saf.pending == 1


def test_new_packets_queue_behind_backlog():
    inner = FlakyClient()
    saf = StoreAndForwardClient(inner)
    inner.online = False
    saf.send(_pkt(1))  # queued (offline)
    inner.online = True
    saf.send(_pkt(2))  # queue non-empty -> also queued
    assert saf.pending == 2
    saf.flush()
    assert _order(inner) == [1, 2]  # order preserved, nothing overtakes


# ---------------------------------------------------------------------------
# backoff hint
# ---------------------------------------------------------------------------
def test_backoff_grows_on_repeated_failure():
    inner = FlakyClient()
    inner.online = False
    saf = StoreAndForwardClient(inner, base_backoff_s=1.0, max_backoff_s=10.0)
    saf.send(_pkt(1))
    assert saf.next_retry_delay_s == 1.0  # queued, not yet failed a flush
    saf.flush()
    assert saf.next_retry_delay_s == 1.0  # 1 failure: base
    saf.flush()
    assert saf.next_retry_delay_s == 2.0  # 2 failures: base * 2
    saf.flush()
    assert saf.next_retry_delay_s == 4.0


def test_backoff_zero_when_empty_and_resets_on_success():
    inner = FlakyClient()
    saf = StoreAndForwardClient(inner)
    assert saf.next_retry_delay_s == 0.0
    inner.online = False
    saf.send(_pkt(1))
    saf.flush()  # a failure
    inner.online = True
    saf.flush()  # drains
    assert saf.pending == 0
    assert saf.next_retry_delay_s == 0.0  # reset


# ---------------------------------------------------------------------------
# spool persistence
# ---------------------------------------------------------------------------
def test_spool_survives_restart(tmp_path):
    spool = str(tmp_path / "spool.jsonl")
    inner = FlakyClient()
    inner.online = False
    saf = StoreAndForwardClient(inner, spool_path=spool)
    for i in range(1, 4):
        saf.send(_pkt(i))
    # a fresh client on the same spool restores the backlog
    restored = StoreAndForwardClient(FlakyClient(), spool_path=spool)
    assert restored.pending == 3
    restored.client.online = True
    assert restored.flush() == 3
    assert _order(restored.client) == [1, 2, 3]


def test_spool_shrinks_after_flush(tmp_path):
    spool = str(tmp_path / "spool.jsonl")
    inner = FlakyClient()
    inner.online = False
    saf = StoreAndForwardClient(inner, spool_path=spool)
    for i in range(1, 4):
        saf.send(_pkt(i))
    inner.online = True
    saf.flush()
    reloaded = StoreAndForwardClient(FlakyClient(), spool_path=spool)
    assert reloaded.pending == 0  # spool rewritten empty


# ---------------------------------------------------------------------------
# bounded queue
# ---------------------------------------------------------------------------
def test_max_queue_drops_oldest():
    inner = FlakyClient()
    inner.online = False
    saf = StoreAndForwardClient(inner, max_queue=2)
    for i in range(5):
        saf.send(_pkt(i))
    assert saf.pending == 2
    assert saf.n_dropped == 3
    # the two newest survive
    inner.online = True
    saf.flush()
    assert _order(inner) == [3, 4]


def test_rejects_non_client():
    with pytest.raises(TypeError):
        StoreAndForwardClient(object())


def test_rejects_bad_max_queue():
    with pytest.raises(ValueError):
        StoreAndForwardClient(FlakyClient(), max_queue=0)
