# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for AgentWorker (Phase 4)."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

import logging

from pycsamt.app.desktop.workers.agent_worker import (
    AgentWorker,
    _SignalLogHandler,
)

# ── _SignalLogHandler ─────────────────────────────────────────────────────


def test_signal_log_handler_emits_line():
    received = []
    handler = _SignalLogHandler(received.append)
    handler.setFormatter(logging.Formatter("%(message)s"))
    record = logging.LogRecord("test", logging.INFO, "", 0, "hello", (), None)
    handler.emit(record)
    assert any("hello" in line for line in received)


def test_signal_log_handler_ignores_exceptions():
    def bad_fn(line):
        raise RuntimeError("sink error")

    handler = _SignalLogHandler(bad_fn)
    record = logging.LogRecord("test", logging.INFO, "", 0, "msg", (), None)
    handler.emit(record)  # must not propagate


# ── AgentWorker construction ──────────────────────────────────────────────


def test_agent_worker_creates(qapp):
    w = AgentWorker(
        agent_name="QC Quicklook",
        sites=None,
        params={},
    )
    assert w is not None


def test_agent_worker_unknown_name_emits_error(qapp):
    errors = []
    w = AgentWorker(agent_name="UNKNOWN_AGENT", sites=None, params={})
    w.error.connect(errors.append)
    w.run()  # run synchronously in test
    assert len(errors) == 1
    assert "UNKNOWN_AGENT" in errors[0]


def test_agent_worker_cancel_flag(qapp):
    w = AgentWorker(agent_name="QC Quicklook", sites=None, params={})
    assert not w._cancelled
    w.cancel()
    assert w._cancelled


# ── Processing mode — no-data graceful failure ────────────────────────────


def test_agent_worker_processing_no_sites_emits_error(qapp):
    errors = []
    w = AgentWorker(
        agent_name="QC Quicklook",
        sites=None,
        params={},
    )
    w.error.connect(errors.append)
    w.run()
    # Should emit an error because sites=None
    assert len(errors) == 1


# ── Signals exist ─────────────────────────────────────────────────────────


def test_agent_worker_has_required_signals(qapp):
    w = AgentWorker("QC Quicklook", sites=None, params={})
    assert hasattr(w, "log_line")
    assert hasattr(w, "progress")
    assert hasattr(w, "result")
    assert hasattr(w, "error")
