# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for agent_window.py: FigureViewerWindow, _HoverRevealButton,
_format_result_html, and AgentRunnerWindow.

Strategy
--------
* ``AgentWorker`` / ``ChatWorker`` are genuine ``QThread`` subclasses
  imported *locally* inside ``_on_run`` / ``_on_chat_send``; run-triggering
  tests monkeypatch the source module (``pycsamt.app.desktop.workers
  .agent_worker``) with a lightweight fake (the ``_FakeSignal`` idiom used
  throughout this test suite) instead of starting a real thread.
* ``QFileDialog`` / ``QMessageBox`` are monkeypatched to avoid real native
  dialogs under the offscreen platform.
"""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QKeyEvent
from PySide6.QtWidgets import QFileDialog, QMessageBox

from pycsamt.app.desktop.windows.agent_window import (
    AgentRunnerWindow,
    FigureViewerWindow,
    _format_result_html,
    _HoverRevealButton,
)


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_agent_worker_cls(on_start=None):
    class _FakeAgentWorker:
        instances = []

        def __init__(self, agent_name, sites, params, api_key, parent=None):
            self.agent_name = agent_name
            self.sites = sites
            self.params = params
            self.api_key = api_key
            self.log_line = _FakeSignal()
            self.progress = _FakeSignal()
            self.result = _FakeSignal()
            self.error = _FakeSignal()
            self.finished = _FakeSignal()
            self.agent_ready = _FakeSignal()
            self._running = False
            _FakeAgentWorker.instances.append(self)

        def isRunning(self):
            return self._running

        def cancel(self):
            pass

        def quit(self):
            pass

        def wait(self, ms):
            pass

        def start(self):
            self._running = True
            if on_start is not None:
                on_start(self)
            self._running = False

    return _FakeAgentWorker


def _fake_chat_worker_cls(on_start=None):
    class _FakeChatWorker:
        instances = []

        def __init__(self, agent, question, context, parent=None):
            self.agent = agent
            self.question = question
            self.context = context
            self.reply_done = _FakeSignal()
            self.error = _FakeSignal()
            _FakeChatWorker.instances.append(self)

        def start(self):
            if on_start is not None:
                on_start(self)

    return _FakeChatWorker


@pytest.fixture
def fig():
    f, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 4, 9])
    yield f
    plt.close(f)


@pytest.fixture
def win(qapp):
    w = AgentRunnerWindow(parent=None)
    w.show()
    yield w
    w.close()


# ═════════════════════════════════════════════════════════════════════════════
# FigureViewerWindow
# ═════════════════════════════════════════════════════════════════════════════


class TestFigureViewerWindow:
    def test_construction_renders_pixmap(self, qapp, fig):
        v = FigureViewerWindow(figure=fig, title="My Figure", parent=None)
        assert v.windowTitle() == "My Figure"
        assert not v._img_label.pixmap().isNull()
        v.close()

    def test_render_error_shows_message(self, qapp, fig, monkeypatch):
        def _boom(*a, **k):
            raise RuntimeError("savefig failed")

        monkeypatch.setattr(fig, "savefig", _boom)
        v = FigureViewerWindow(figure=fig, parent=None)
        assert "Render error" in v._img_label.text()
        v.close()

    def test_dpi_changed_rerenders(self, qapp, fig):
        v = FigureViewerWindow(figure=fig, parent=None)
        v._dpi_combo.setCurrentText("300")
        assert v._dpi == 300
        v.close()

    def test_dpi_changed_invalid_text_noop(self, qapp, fig):
        v = FigureViewerWindow(figure=fig, parent=None)
        v._on_dpi_changed("not-a-number")
        assert v._dpi == 150
        v.close()

    def test_wheel_event_ctrl_zooms_in(self, qapp, fig):
        from PySide6.QtCore import QPoint, QPointF
        from PySide6.QtGui import QWheelEvent

        v = FigureViewerWindow(figure=fig, parent=None)
        ev = QWheelEvent(
            QPointF(0, 0),
            QPointF(0, 0),
            QPoint(0, 0),
            QPoint(0, 120),
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.ControlModifier,
            Qt.ScrollPhase.NoScrollPhase,
            False,
        )
        v.wheelEvent(ev)
        assert v._dpi == 175
        v.close()

    def test_wheel_event_without_ctrl_passthrough(self, qapp, fig):
        from PySide6.QtCore import QPoint, QPointF
        from PySide6.QtGui import QWheelEvent

        v = FigureViewerWindow(figure=fig, parent=None)
        ev = QWheelEvent(
            QPointF(0, 0),
            QPointF(0, 0),
            QPoint(0, 0),
            QPoint(0, 120),
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
            Qt.ScrollPhase.NoScrollPhase,
            False,
        )
        v.wheelEvent(ev)  # delegates to super(); must not raise
        assert v._dpi == 150
        v.close()

    def test_export_cancelled(self, qapp, fig, monkeypatch):
        v = FigureViewerWindow(figure=fig, parent=None)
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        v._export("png")  # must not raise
        v.close()

    def test_export_success(self, qapp, fig, monkeypatch, tmp_path):
        v = FigureViewerWindow(figure=fig, parent=None)
        out = str(tmp_path / "out.png")
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (out, "")),
        )
        v._export("png")
        from pathlib import Path

        assert Path(out).exists()
        v.close()

    def test_export_failure_shows_warning(self, qapp, fig, monkeypatch):
        v = FigureViewerWindow(figure=fig, parent=None)
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("/no/such/dir/out.png", "")),
        )
        calls = []
        monkeypatch.setattr(
            QMessageBox,
            "warning",
            staticmethod(lambda *a, **k: calls.append(1)),
        )
        v._export("png")
        assert calls == [1]
        v.close()


# ═════════════════════════════════════════════════════════════════════════════
# _HoverRevealButton
# ═════════════════════════════════════════════════════════════════════════════


class TestHoverRevealButton:
    def test_enter_leave_toggles_style(self, qapp):
        from PySide6.QtWidgets import QWidget

        container = QWidget()
        container.resize(200, 200)
        btn = _HoverRevealButton(container, get_canvas_fn=lambda: None)
        enter_ev = QEvent(QEvent.Type.Enter)
        btn.eventFilter(container, enter_ev)
        assert btn._shown is True
        leave_ev = QEvent(QEvent.Type.Leave)
        btn.eventFilter(container, leave_ev)
        assert btn._shown is False

    def test_resize_repositions(self, qapp):
        from PySide6.QtWidgets import QWidget

        container = QWidget()
        container.resize(300, 300)
        btn = _HoverRevealButton(container, get_canvas_fn=lambda: None)
        resize_ev = QEvent(QEvent.Type.Resize)
        btn.eventFilter(container, resize_ev)  # must not raise

    def test_click_without_canvas_noop(self, qapp):
        from PySide6.QtWidgets import QWidget

        container = QWidget()
        btn = _HoverRevealButton(container, get_canvas_fn=lambda: None)
        btn._on_click()  # canvas is None -> early return

    def test_click_with_canvas_opens_viewer(self, qapp, fig, monkeypatch):
        from PySide6.QtWidgets import QWidget

        from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

        container = QWidget()
        canvas = MplCanvas(container, toolbar=False)
        canvas.figure.clear()
        ax = canvas.figure.add_subplot(111)
        ax.plot([1, 2], [3, 4])
        btn = _HoverRevealButton(container, get_canvas_fn=lambda: canvas)

        shown = []
        monkeypatch.setattr(
            "pycsamt.app.desktop.windows.agent_window.FigureViewerWindow.show",
            lambda self: shown.append(self),
        )
        btn._on_click()
        assert len(shown) == 1


# ═════════════════════════════════════════════════════════════════════════════
# _format_result_html
# ═════════════════════════════════════════════════════════════════════════════


class TestFormatResultHtml:
    def test_figure_only_result(self, fig):
        html = _format_result_html(fig)
        assert "Figure result" in html

    def test_full_metadata_result(self):
        result = SimpleNamespace(
            status="success",
            elapsed_seconds=1.5,
            cost_estimate_usd=0.002,
            summary="Everything looks fine.",
            llm_interpretation="The resistivity anomaly suggests <fluid>.",
            warnings=["low SNR at station A", "missing tipper"],
        )
        html = _format_result_html(result)
        assert "SUCCESS" in html
        assert "1.5 s" in html
        assert "Everything looks fine." in html
        assert "&lt;fluid&gt;" in html
        assert "Warnings (2)" in html

    def test_error_status_badge(self):
        result = SimpleNamespace(
            status="error", elapsed_seconds=None, cost_estimate_usd=None
        )
        html = _format_result_html(result)
        assert "ERROR" in html

    def test_empty_result_falls_back_to_repr(self):
        result = SimpleNamespace()
        html = _format_result_html(result)
        assert "namespace" in html.lower() or "SimpleNamespace" in html


# ═════════════════════════════════════════════════════════════════════════════
# AgentRunnerWindow
# ═════════════════════════════════════════════════════════════════════════════


class TestConstruction:
    def test_window_title(self, win):
        assert "Agent Runner" in win.windowTitle()

    def test_run_button_disabled_initially(self, win):
        assert not win._btn_run.isEnabled()
        assert not win._btn_stop.isEnabled()

    def test_param_groups_hidden_initially(self, win):
        assert not win._grp_params.isVisible()
        assert not win._grp_context.isVisible()
        assert not win._grp_apikey.isVisible()

    def test_chat_history_shows_welcome(self, win):
        assert "Chat with the agent" in win._chat_history.toHtml()

    def test_set_app_controller_populates_api_key(self, win):
        ctrl = SimpleNamespace(session=SimpleNamespace(api_key="sk-abc"))
        win.set_app_controller(ctrl)
        assert win._edit_apikey.text() == "sk-abc"

    def test_set_app_controller_no_session(self, win):
        ctrl = SimpleNamespace()
        win.set_app_controller(ctrl)
        assert win._edit_apikey.text() == ""

    def test_append_log_adds_timestamped_line(self, win):
        win.append_log("hello world")
        assert "hello world" in win._log_text.toPlainText()

    def test_result_figure_property_default_none(self, win):
        assert win.result_figure is None


class TestAgentSelection:
    def test_select_processing_agent(self, win):
        from pycsamt.app.desktop.agent_registry import AGENT_REGISTRY

        proc_name = next(
            n for n, e in AGENT_REGISTRY.items() if e["type"] == "processing"
        )
        win._on_agent_selected(proc_name)
        assert win._current_agent == proc_name
        assert win._detail_header.text() == proc_name
        assert win._btn_run.isEnabled()
        assert not win._grp_context.isVisible()
        assert not win._grp_apikey.isVisible()

    def test_select_llm_agent_shows_context_and_apikey(self, win):
        from pycsamt.app.desktop.agent_registry import AGENT_REGISTRY

        llm_name = next(n for n, e in AGENT_REGISTRY.items() if e["type"] == "llm")
        win._on_agent_selected(llm_name)
        assert win._grp_context.isVisible()
        assert win._grp_apikey.isVisible()

    def test_select_empty_name_disables_run(self, win):
        win._on_agent_selected("")
        assert not win._btn_run.isEnabled()

    def test_activate_agent_runs_immediately(self, win, monkeypatch):
        from pycsamt.app.desktop.agent_registry import AGENT_REGISTRY

        proc_name = next(
            n for n, e in AGENT_REGISTRY.items() if e["type"] == "processing"
        )
        calls = []
        monkeypatch.setattr(win, "_on_run", lambda: calls.append(1))
        win._on_agent_activated(proc_name)
        assert win._current_agent == proc_name
        assert calls == [1]


class TestParamForm:
    def test_rebuild_form_with_params(self, win):
        from pycsamt.app.desktop.agent_registry import AGENT_REGISTRY

        name = next(n for n, e in AGENT_REGISTRY.items() if e.get("params"))
        win._rebuild_params_form(name)
        assert win._grp_params.isVisible()
        assert len(win._params_widgets) == len(AGENT_REGISTRY[name]["params"])

    def test_rebuild_form_unknown_agent(self, win):
        win._rebuild_params_form("definitely-not-a-real-agent")
        assert not win._grp_params.isVisible()
        assert win._params_widgets == {}

    def test_make_widget_all_kinds(self, win):
        combo = win._make_param_widget(
            {"type": "combo", "options": ["a", "b"], "default": "b"}
        )
        assert combo.currentText() == "b"

        fl = win._make_param_widget(
            {"type": "float", "range": (0.0, 10.0), "default": 3.5}
        )
        assert fl.value() == pytest.approx(3.5)

        it = win._make_param_widget({"type": "int", "range": (0, 100), "default": 7})
        assert it.value() == 7

        bl = win._make_param_widget({"type": "bool", "default": True})
        assert bl.isChecked()

        le = win._make_param_widget({"type": "str", "default": "hi"})
        assert le.text() == "hi"

    def test_collect_params_includes_context(self, win):
        from pycsamt.app.desktop.agent_registry import AGENT_REGISTRY

        name = next(n for n, e in AGENT_REGISTRY.items() if e.get("params"))
        win._current_agent = name
        win._rebuild_params_form(name)
        win._ctx_edit.setPlainText("semi-arid terrain")
        vals = win._collect_params()
        assert vals["context"] == "semi-arid terrain"
        assert vals["user_prompt"] == "semi-arid terrain"

    def test_collect_params_no_context(self, win):
        from pycsamt.app.desktop.agent_registry import AGENT_REGISTRY

        name = next(n for n, e in AGENT_REGISTRY.items() if e.get("params"))
        win._current_agent = name
        win._rebuild_params_form(name)
        vals = win._collect_params()
        assert "context" not in vals


class TestRunStop:
    def test_run_no_controller(self, win):
        win._current_agent = "whatever"
        win._on_run()
        assert win._status_lbl.text() == "Load survey data first."

    def test_run_no_agent_selected(self, win):
        win._ctrl = SimpleNamespace(sites=object())
        win._current_agent = ""
        win._on_run()
        assert win._status_lbl.text() == "Select an agent first."

    def test_run_starts_worker(self, win, monkeypatch):
        win._ctrl = SimpleNamespace(sites=object())
        win._current_agent = "SomeAgent"
        fake_cls = _fake_agent_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.agent_worker.AgentWorker", fake_cls
        )
        win._on_run()
        assert len(fake_cls.instances) == 1
        assert not win._btn_run.isEnabled()
        assert win._btn_stop.isEnabled()

    def test_run_full_cycle_with_result(self, win, monkeypatch):
        win._ctrl = SimpleNamespace(sites=object())
        win._current_agent = "SomeAgent"

        result = SimpleNamespace(
            status="success",
            elapsed_seconds=0.5,
            cost_estimate_usd=None,
            summary="ok",
            llm_interpretation="",
            warnings=None,
            data={},
        )

        def _on_start(worker):
            worker.log_line.emit("started")
            worker.progress.emit(50)
            worker.agent_ready.emit(SimpleNamespace(api_key="sk-x"))
            worker.result.emit(result)
            worker.finished.emit()

        fake_cls = _fake_agent_worker_cls(_on_start)
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.agent_worker.AgentWorker", fake_cls
        )
        win._on_run()
        assert "started" in win._log_text.toPlainText()
        assert win._status_lbl.text() == "Finished."
        assert win._last_agent is not None
        assert "ok" in win._summary_browser.toHtml()

    def test_run_result_with_figure_shows_result_tab(self, win, monkeypatch, fig):
        win._ctrl = SimpleNamespace(sites=object())
        win._current_agent = "SomeAgent"
        result = SimpleNamespace(figure=fig, data={})

        def _on_start(worker):
            worker.result.emit(result)

        fake_cls = _fake_agent_worker_cls(_on_start)
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.agent_worker.AgentWorker", fake_cls
        )
        win._on_run()
        assert win._tabs.currentIndex() == 1  # Result tab
        assert win.result_figure is fig

    def test_run_error_and_stop(self, win, monkeypatch):
        win._ctrl = SimpleNamespace(sites=object())
        win._current_agent = "SomeAgent"

        def _on_start(worker):
            worker.error.emit("boom")

        fake_cls = _fake_agent_worker_cls(_on_start)
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.agent_worker.AgentWorker", fake_cls
        )
        win._on_run()
        assert "ERROR: boom" in win._status_lbl.text()
        assert "ERROR: boom" in win._log_text.toPlainText()

    def test_stop_cancels_running_worker(self, win, monkeypatch):
        fake_cls = _fake_agent_worker_cls()
        w = fake_cls("A", object(), {}, None)
        w._running = True
        win._worker = w
        win._current_agent = "A"
        win._on_stop()
        assert win._btn_run.isEnabled()
        assert not win._btn_stop.isEnabled()

    def test_stop_no_worker_noop(self, win):
        win._on_stop()  # must not raise


class TestCopyInterpretation:
    def test_copy_noop_when_empty(self, win):
        win._last_interpretation = ""
        win._on_copy_interpretation()  # must not raise

    def test_copy_sets_clipboard(self, win):
        win._last_interpretation = "some interpretation text"
        win._on_copy_interpretation()
        from PySide6.QtWidgets import QApplication

        assert QApplication.clipboard().text() == "some interpretation text"


class TestChat:
    def test_chat_clear_resets_welcome(self, win):
        win._chat_history.append("<p>something</p>")
        win._on_chat_clear()
        assert "Chat with the agent" in win._chat_history.toHtml()

    def test_send_empty_question_noop(self, win):
        win._chat_input.setPlainText("   ")
        win._on_chat_send()  # stripped to empty -> early return

    def test_send_without_agent_shows_status(self, win):
        win._last_agent = None
        win._chat_input.setPlainText("what does this mean?")
        win._on_chat_send()
        assert "Run an LLM agent first." in win._chat_status.text()

    def test_send_without_api_key_shows_status(self, win):
        win._last_agent = SimpleNamespace(api_key=None)
        win._chat_input.setPlainText("question")
        win._on_chat_send()
        assert "No API key" in win._chat_status.text()

    def test_send_starts_chat_worker(self, win, monkeypatch):
        win._last_agent = SimpleNamespace(api_key="sk-x")
        win._chat_input.setPlainText("what next?")
        fake_cls = _fake_chat_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.agent_worker.ChatWorker", fake_cls
        )
        win._on_chat_send()
        assert len(fake_cls.instances) == 1
        assert not win._btn_chat_send.isEnabled()
        assert "what next?" in win._chat_history.toHtml()

    def test_chat_reply_reenables_input(self, win, monkeypatch):
        win._last_agent = SimpleNamespace(api_key="sk-x")
        win._chat_input.setPlainText("q")

        def _on_start(worker):
            worker.reply_done.emit("here is the answer")

        fake_cls = _fake_chat_worker_cls(_on_start)
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.agent_worker.ChatWorker", fake_cls
        )
        win._on_chat_send()
        assert win._btn_chat_send.isEnabled()
        assert "here is the answer" in win._chat_history.toHtml()

    def test_chat_error_shows_status(self, win, monkeypatch):
        win._last_agent = SimpleNamespace(api_key="sk-x")
        win._chat_input.setPlainText("q")

        def _on_start(worker):
            worker.error.emit("llm down")

        fake_cls = _fake_chat_worker_cls(_on_start)
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.agent_worker.ChatWorker", fake_cls
        )
        win._on_chat_send()
        assert "Error: llm down" in win._chat_status.text()
        assert win._btn_chat_send.isEnabled()

    def test_append_chat_message_all_roles(self, win):
        win._append_chat_message("system", "sys msg")
        win._append_chat_message("user", "user msg")
        win._append_chat_message("agent", "agent msg")
        html = win._chat_history.toHtml()
        assert "sys msg" in html
        assert "user msg" in html
        assert "agent msg" in html

    def test_event_filter_enter_sends(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(win, "_on_chat_send", lambda: calls.append(1))
        ev = QKeyEvent(
            QEvent.Type.KeyPress,
            Qt.Key.Key_Return,
            Qt.KeyboardModifier.NoModifier,
        )
        handled = win.eventFilter(win._chat_input, ev)
        assert handled is True
        assert calls == [1]

    def test_event_filter_shift_enter_inserts_newline(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(win, "_on_chat_send", lambda: calls.append(1))
        ev = QKeyEvent(
            QEvent.Type.KeyPress,
            Qt.Key.Key_Return,
            Qt.KeyboardModifier.ShiftModifier,
        )
        handled = win.eventFilter(win._chat_input, ev)
        assert handled is False
        assert calls == []

    def test_event_filter_other_widget_passthrough(self, win):
        ev = QEvent(QEvent.Type.KeyPress)
        win.eventFilter(win, ev)  # not the chat input -> falls through


class TestSetDarkMode:
    def test_set_dark_mode_applies_theme_to_browser(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(win._browser, "apply_theme", lambda d: calls.append(d))
        win.set_dark_mode(False)
        assert calls == [False]
