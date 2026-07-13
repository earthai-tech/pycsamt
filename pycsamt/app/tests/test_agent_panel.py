# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for AgentPanel (Phase 4)."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.panels.agent_panel import AgentPanel
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas


@pytest.fixture
def panel(qapp):
    p = AgentPanel()
    yield p
    p.close()


# ── Construction ──────────────────────────────────────────────────────────


def test_agent_panel_creates(qapp):
    p = AgentPanel()
    assert p is not None
    p.close()


def test_has_run_button(panel):
    assert panel._btn_run is not None


def test_has_export_button(panel):
    assert panel._btn_export is not None


def test_has_stop_button(panel):
    assert panel._btn_stop is not None
    assert not panel._btn_stop.isEnabled()


def test_has_progress_bar(panel):
    assert panel._progress is not None
    assert not panel._progress.isVisible()


def test_has_tab_widget(panel):
    assert panel._tabs is not None
    assert panel._tabs.count() == 3


def test_tab_labels(panel):
    labels = [panel._tabs.tabText(i) for i in range(3)]
    assert "Log" in labels
    assert "Result" in labels
    assert "Summary" in labels


def test_has_result_canvas(panel):
    assert isinstance(panel._result_canvas, MplCanvas)


def test_has_log_text(panel):
    assert panel._log_text is not None
    assert panel._log_text.isReadOnly()


def test_has_summary_browser(panel):
    assert panel._summary_browser is not None


# ── append_log ────────────────────────────────────────────────────────────


def test_append_log_adds_text(panel):
    panel.append_log("test message")
    assert "test message" in panel._log_text.toPlainText()


def test_append_log_adds_timestamp(panel):
    panel.append_log("ts check")
    text = panel._log_text.toPlainText()
    assert "[" in text and "]" in text


# ── AppController wiring ──────────────────────────────────────────────────


def test_set_app_controller(panel):
    class FakeCtrl:
        sites = None
        session = type("S", (), {"api_key": ""})()

    panel.set_app_controller(FakeCtrl())
    assert panel._ctrl is not None


def test_run_without_data_shows_message(panel):
    class FakeCtrl:
        sites = None
        session = type("S", (), {"api_key": ""})()

    panel.set_app_controller(FakeCtrl())
    panel._on_run()  # must not raise — shows status message instead
    assert "Load" in panel._status_lbl.text()


# ── Worker signal handlers ────────────────────────────────────────────────


def test_on_log_line_appends(panel):
    panel._on_log_line("signal line")
    assert "signal line" in panel._log_text.toPlainText()


def test_on_error_updates_status(panel):
    panel._on_error("something went wrong")
    assert (
        "Failed" in panel._status_lbl.text()
        or "ERROR" in panel._log_text.toPlainText()
    )


def test_on_worker_finished_restores_buttons(panel):
    panel._btn_run.setEnabled(False)
    panel._btn_stop.setEnabled(True)
    panel._status_lbl.setText("Running: QC…")
    panel._on_worker_finished()
    assert panel._btn_run.isEnabled()
    assert not panel._btn_stop.isEnabled()
