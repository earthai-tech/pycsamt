# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Extra tests for AgentPanel, on top of test_agent_panel.py (construction,
append_log, and the no-data _on_run guard). This file drives the full
run/stop/result/export flow and the module-level renderable-extraction
helpers, with RunAgentDialog and AgentWorker replaced by lightweight
stand-ins (no real LLM call, no real QThread).
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.app.desktop.panels.agent_panel import (
    AgentPanel,
    _extract_all_renderables,
    _extract_renderable,
    _is_renderable,
    _label,
)


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls():
    class _FakeWorker:
        instances = []

        def __init__(self, **kw):
            self.kw = kw
            self.log_line = _FakeSignal()
            self.progress = _FakeSignal()
            self.result = _FakeSignal()
            self.error = _FakeSignal()
            self.finished = _FakeSignal()
            self.started = False
            self.cancelled = False
            self.quit_called = False
            _FakeWorker.instances.append(self)

        def start(self):
            self.started = True

        def cancel(self):
            self.cancelled = True

        def quit(self):
            self.quit_called = True

    return _FakeWorker


class _FakeCtrl:
    def __init__(self, sites=None, api_key=""):
        self.sites = sites
        self.session = SimpleNamespace(api_key=api_key)


def _fake_dialog(accepted, agent="qc", params=None, api_key=None):
    _params = params or {}

    class _Fake:
        DialogCode = SimpleNamespace(Accepted=1, Rejected=0)
        selected_agent = agent
        params = _params

        def __init__(self, *a, **kw):
            self.api_key = api_key

        def exec(self):
            return 1 if accepted else 0

    return _Fake


@pytest.fixture
def panel(qapp):
    p = AgentPanel()
    yield p
    p.close()


# ── module helpers: _is_renderable / _label ──────────────────────────────


def test_is_renderable_figure():
    fig = plt.figure()
    assert _is_renderable(fig) is True
    plt.close(fig)


def test_is_renderable_axes():
    fig, ax = plt.subplots()
    assert _is_renderable(ax) is True
    plt.close(fig)


def test_is_renderable_plain_object_false():
    assert _is_renderable(object()) is False


def test_label_snake_case_to_title():
    assert _label("ss_summary_plot") == "Ss Summary Plot"


# ── _extract_renderable / _extract_all_renderables ───────────────────────


def test_extract_renderable_direct_figure():
    fig = plt.figure()
    assert _extract_renderable(fig) is fig
    plt.close(fig)


def test_extract_renderable_none_for_plain_result():
    assert _extract_renderable(SimpleNamespace(status="ok")) is None


def test_extract_all_renderables_from_data_dict():
    fig = plt.figure()
    result = SimpleNamespace(data={"main_plot": fig, "other": 42})
    found = _extract_all_renderables(result)
    assert found == [("Main Plot", fig)]
    plt.close(fig)


def test_extract_all_renderables_nested_figures_dict():
    fig1, fig2 = plt.figure(), plt.figure()
    result = SimpleNamespace(
        data={"figures": {"a_plot": fig1, "b_plot": fig2}}
    )
    found = _extract_all_renderables(result)
    assert found == [("A Plot", fig1), ("B Plot", fig2)]
    plt.close(fig1)
    plt.close(fig2)


def test_extract_all_renderables_figure_attribute():
    fig = plt.figure()
    result = SimpleNamespace(figure=fig)
    found = _extract_all_renderables(result)
    assert found == [("Figure", fig)]
    plt.close(fig)


def test_extract_all_renderables_no_duplicate_when_already_found():
    fig = plt.figure()
    result = SimpleNamespace(data={"x": fig}, figure=fig)
    found = _extract_all_renderables(result)
    assert found == [("X", fig)]
    plt.close(fig)


def test_extract_all_renderables_empty_for_plain_result():
    assert _extract_all_renderables(SimpleNamespace(status="ok")) == []


# ── _on_run full flow ──────────────────────────────────────────────────


def test_on_run_rejected_dialog_no_worker(panel, monkeypatch):
    panel.set_app_controller(_FakeCtrl(sites=["S1"]))
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.run_agent_dlg.RunAgentDialog",
        _fake_dialog(accepted=False),
    )
    panel._on_run()
    assert panel._worker is None


def test_on_run_accepted_starts_worker(panel, monkeypatch):
    panel.set_app_controller(_FakeCtrl(sites=["S1"], api_key="sk-existing"))
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.run_agent_dlg.RunAgentDialog",
        _fake_dialog(accepted=True, agent="strike", params={"k": 1}, api_key="sk-x"),
    )
    fake_worker_cls = _fake_worker_cls()
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.agent_worker.AgentWorker", fake_worker_cls
    )
    started_names = []
    panel.agent_started.connect(started_names.append)

    panel._on_run()

    assert started_names == ["strike"]
    assert panel._worker is not None
    assert panel._worker.started is True
    assert panel._worker.kw["agent_name"] == "strike"
    assert panel._worker.kw["sites"] == ["S1"]
    assert panel._worker.kw["api_key"] == "sk-x"
    assert not panel._btn_run.isEnabled()
    assert panel._btn_stop.isEnabled()
    assert "Running: strike" in panel._status_lbl.text()


def test_on_run_uses_session_api_key_when_dialog_has_none(panel, monkeypatch):
    panel.set_app_controller(_FakeCtrl(sites=["S1"], api_key="sk-session"))
    captured = {}

    class _Fake:
        DialogCode = SimpleNamespace(Accepted=1, Rejected=0)
        selected_agent = "qc"
        params = {}
        api_key = None

        def __init__(self, *a, **kw):
            captured["api_key_passed_in"] = kw.get("api_key")

        def exec(self):
            return 1

    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.run_agent_dlg.RunAgentDialog", _Fake
    )
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.agent_worker.AgentWorker",
        _fake_worker_cls(),
    )
    panel._on_run()
    assert captured["api_key_passed_in"] == "sk-session"


# ── _on_stop ──────────────────────────────────────────────────────────────


def test_on_stop_with_worker(panel, monkeypatch):
    panel.set_app_controller(_FakeCtrl(sites=["S1"]))
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.run_agent_dlg.RunAgentDialog",
        _fake_dialog(accepted=True),
    )
    fake_worker_cls = _fake_worker_cls()
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.agent_worker.AgentWorker", fake_worker_cls
    )
    panel._on_run()
    worker = panel._worker
    panel._on_stop()
    assert worker.cancelled is True
    assert worker.quit_called is True
    assert panel._status_lbl.text() == "Stopped."
    assert panel._btn_run.isEnabled()
    assert not panel._btn_stop.isEnabled()


def test_on_stop_without_worker_noop(panel):
    panel._on_stop()  # _worker is None -> guarded, must not raise


# ── _on_result / _render_result / _render_summary ────────────────────────


def test_on_result_with_figure_enables_export(panel):
    fig = plt.figure()
    panel._on_result(SimpleNamespace(figure=fig))
    assert panel._btn_export.isEnabled()
    assert panel._result_figure is fig
    plt.close(fig)


def test_on_result_without_figure_keeps_export_disabled(panel):
    panel._on_result(SimpleNamespace(status="ok"))
    assert not panel._btn_export.isEnabled()


def test_render_result_ends_on_summary_tab(panel):
    """_render_result sets the Result tab (1) active for a renderable, but
    the subsequent _render_summary call it always makes unconditionally
    re-switches to the Summary tab (2) — since _render_summary's `lines`
    list is never actually empty (a placeholder is appended when there's
    nothing else), `if lines: setCurrentIndex(2)` always fires. So the
    Result tab is never the one left visible, even though a figure was
    just rendered into it. Documented here as discovered (not fixed)."""
    fig = plt.figure()
    panel._render_result(SimpleNamespace(figure=fig))
    assert panel._tabs.currentIndex() == 2
    plt.close(fig)


def test_render_result_show_figure_exception_is_swallowed(panel, monkeypatch):
    fig = plt.figure()
    monkeypatch.setattr(
        panel._result_canvas,
        "show_figure",
        mock.Mock(side_effect=RuntimeError("boom")),
    )
    panel._render_result(SimpleNamespace(figure=fig))  # must not raise
    plt.close(fig)


def test_render_summary_all_fields(panel):
    result = SimpleNamespace(
        status="ok",
        summary="did the thing",
        elapsed_seconds=3.456,
        cost_estimate="$0.02",
        llm_interpretation="line1\nline2",
        warnings=["watch out", "careful"],
    )
    panel._render_summary(result)
    html = panel._summary_browser.toHtml()
    assert "did the thing" in html
    assert "3.5" in html
    assert "$0.02" in html
    assert "watch out" in html
    assert panel._tabs.currentIndex() == 2


def test_render_summary_no_fields_shows_placeholder(panel):
    panel._render_summary(SimpleNamespace())
    text = panel._summary_browser.toPlainText()
    assert "Processing result" in text


def test_render_summary_summary_empty_string_not_shown(panel):
    panel._render_summary(SimpleNamespace(summary=""))
    text = panel._summary_browser.toPlainText()
    assert "Summary:" not in text


# ── _on_export ────────────────────────────────────────────────────────────


def test_on_export_with_result_figure(panel, monkeypatch):
    fig = plt.figure()
    panel._result_figure = fig
    captured = []

    class _FakeExport:
        def __init__(self, *a, **kw):
            captured.append(kw)

        def exec(self):
            return None

    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog", _FakeExport
    )
    panel._on_export()
    assert captured[0]["figure"] is fig
    plt.close(fig)


def test_on_export_falls_back_to_canvas_figure(panel, monkeypatch):
    panel._result_figure = None
    captured = []

    class _FakeExport:
        def __init__(self, *a, **kw):
            captured.append(kw)

        def exec(self):
            return None

    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog", _FakeExport
    )
    panel._on_export()
    assert captured[0]["figure"] is panel._result_canvas.figure
