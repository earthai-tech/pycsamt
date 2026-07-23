# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for PipelineWindow (pycsamt.app.desktop.windows.pipeline_window).

Strategy
--------
* ``PipelineController`` (a plain-Python dataclass model) is exercised for
  real through the window's own ``_ctrl`` — no mocking needed there.
* The real background worker, ``pycsamt.app.desktop.workers.pipeline_worker
  .PipelineWorker``, is a genuine ``QThread``. Rather than risk a hang by
  calling ``.start()`` under the offscreen platform, ``_start_worker``-
  triggering tests monkeypatch the module-level ``PipelineWorker`` name
  with a lightweight fake (mirroring the ``_FakeSignal``/fake-worker idiom
  used in test_recompute_dlg.py / test_frequency_editor_tool.py) whose
  ``.start()`` synchronously fires whichever signals a test wants. The
  worker *slot* methods (``_on_step_started`` etc.) are also exercised
  directly, independent of any worker.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QFileDialog

from pycsamt.app.desktop.controllers.pipeline_controller import StepStatus
from pycsamt.app.desktop.windows.pipeline_window import PipelineWindow


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls(on_start=None):
    class _FakeWorker:
        instances = []

        def __init__(self, ctrl, step_ids, parent=None):
            self.ctrl = ctrl
            self.step_ids = step_ids
            self.step_started = _FakeSignal()
            self.step_done = _FakeSignal()
            self.step_error = _FakeSignal()
            self.step_skipped = _FakeSignal()
            self.log_line = _FakeSignal()
            self.progress = _FakeSignal()
            self.all_done = _FakeSignal()
            self._running = False
            _FakeWorker.instances.append(self)

        def isRunning(self):
            return self._running

        def requestInterruption(self):
            pass

        def start(self):
            self._running = True
            if on_start is not None:
                on_start(self)
            self._running = False

    return _FakeWorker


@pytest.fixture
def win(qapp):
    w = PipelineWindow(parent=None)
    w.show()  # widgets only report real isVisible() once the top level is shown
    yield w
    w.close()


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "Processing Pipeline" in win.windowTitle()

    def test_step_list_populated(self, win):
        assert win._step_list.count() == len(win._ctrl.steps)

    def test_first_step_selected(self, win):
        assert win._selected_step == 0
        assert "Load Data" in win._step_title.text()

    def test_stop_button_disabled_initially(self, win):
        assert not win._btn_stop.isEnabled()

    def test_progress_label_initial(self, win):
        assert win._progress_lbl.text() == f"0 / {len(win._ctrl.steps)} steps done"

    def test_result_group_hidden_initially(self, win):
        assert not win._grp_result.isVisible()


# ── set_input_sites ──────────────────────────────────────────────────────────


class TestSetInputSites:
    def test_marks_step0_done_with_sites(self, win):
        win.set_input_sites(["s1", "s2"])
        assert win._ctrl.steps[0].status == StepStatus.DONE
        assert "stations" in win._ctrl.steps[0].result_info

    def test_none_sites_leaves_step0_pending(self, win):
        win.set_input_sites(None)
        assert win._ctrl.steps[0].status == StepStatus.PENDING

    def test_set_dark_mode_is_noop(self, win):
        win.set_dark_mode(True)  # must not raise
        win.set_dark_mode(False)


# ── Step selection / method combo / param form ──────────────────────────────


class TestStepSelectionAndParams:
    def test_select_step_out_of_range_is_noop(self, win):
        win._select_step(-1)
        win._select_step(999)

    def test_select_step_populates_method_combo(self, win):
        win._select_step(1)  # "QC Screening"
        assert win._combo_method.count() >= 1

    def test_method_changed_resets_params(self, win):
        win._select_step(1)
        step = win._ctrl.steps[1]
        label = win._combo_method.itemText(0)
        win._on_method_changed(label)
        assert step.active_method == step.get_method(step.active_method).name

    def test_on_method_changed_noop_when_no_step_selected(self, win):
        win._selected_step = -1
        win._on_method_changed("anything")  # must not raise

    def test_param_form_hidden_when_no_params(self, win):
        win._select_step(0)  # "current" method has no params
        assert not win._grp_params.isVisible()

    def test_param_form_shown_for_folder_method(self, win):
        win._select_step(0)
        step = win._ctrl.steps[0]
        for i in range(win._combo_method.count()):
            if win._combo_method.itemData(i) == "folder":
                win._combo_method.setCurrentIndex(i)
                win._on_method_changed(win._combo_method.itemText(i))
                break
        assert win._grp_params.isVisible()
        assert "folder" in win._param_widgets

    def test_path_widget_browse_button_updates_edit(self, win, monkeypatch):
        from PySide6.QtWidgets import QLineEdit

        win._select_step(0)
        for i in range(win._combo_method.count()):
            if win._combo_method.itemData(i) == "folder":
                win._combo_method.setCurrentIndex(i)
                win._on_method_changed(win._combo_method.itemText(i))
                break
        container = win._param_widgets["folder"]
        edit = container.findChild(QLineEdit)
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/chosen/dir"),
        )
        win._browse_path(edit)
        assert edit.text() == "/chosen/dir"

    def test_browse_path_cancelled_leaves_edit(self, win, monkeypatch):
        from PySide6.QtWidgets import QLineEdit

        edit = QLineEdit("orig")
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: ""),
        )
        win._browse_path(edit)
        assert edit.text() == "orig"

    def test_select_step_shows_result_for_done_step(self, win):
        win.set_input_sites(["s1"])
        win._select_step(0)
        assert win._grp_result.isVisible()
        assert "stations" in win._result_lbl.text()

    def test_select_step_shows_elapsed_time(self, win):
        step = win._ctrl.steps[0]
        step.status = StepStatus.DONE
        step.elapsed_s = 1.23
        win._select_step(0)
        assert "1.23" in win._elapsed_lbl.text()


# ── Combo widget kinds (via a synthetic ParamSpec-driven step) ──────────────


class TestParamWidgetKinds:
    def test_make_combo_widget(self, win):
        from pycsamt.app.desktop.controllers.pipeline_controller import (
            ParamSpec,
        )

        ps = ParamSpec("combo", "Mode", "b", options=["a", "b", "c"])
        w = win._make_param_widget(ps, "b")
        assert w.currentText() == "b"

    def test_make_float_widget(self, win):
        from pycsamt.app.desktop.controllers.pipeline_controller import (
            ParamSpec,
        )

        ps = ParamSpec("float", "Threshold", 0.5, min_val=0.0, max_val=1.0)
        w = win._make_param_widget(ps, None)
        assert w.value() == pytest.approx(0.5)

    def test_make_int_widget(self, win):
        from pycsamt.app.desktop.controllers.pipeline_controller import (
            ParamSpec,
        )

        ps = ParamSpec("int", "Window", 3, min_val=1, max_val=9)
        w = win._make_param_widget(ps, None)
        assert w.value() == 3

    def test_make_bool_widget(self, win):
        from pycsamt.app.desktop.controllers.pipeline_controller import (
            ParamSpec,
        )

        ps = ParamSpec("bool", "Enabled", True)
        w = win._make_param_widget(ps, True)
        assert w.isChecked()

    def test_make_default_lineedit_widget(self, win):
        from pycsamt.app.desktop.controllers.pipeline_controller import (
            ParamSpec,
        )

        ps = ParamSpec("string", "Label", "hi")
        w = win._make_param_widget(ps, "hi")
        assert w.text() == "hi"


# ── Run / skip / stop / reset (fake worker) ─────────────────────────────────


class TestPipelineExecution:
    def test_run_all_noop_when_all_done(self, win, monkeypatch):
        for s in win._ctrl.steps:
            s.status = StepStatus.DONE
        win._on_run_all()  # logs "All steps already complete." -- no worker

    def test_start_worker_requires_sites(self, win, monkeypatch):
        fake_cls = _fake_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.pipeline_worker.PipelineWorker",
            fake_cls,
        )
        win._on_run_step()  # is_ready() False (no sites loaded)
        assert "No survey data loaded" in win._log_text.toPlainText()
        assert fake_cls.instances == []

    def test_run_step_starts_worker_and_wires_signals(self, win, monkeypatch):
        win.set_input_sites(["s1"])
        progress_snapshot = {}

        def _on_start(worker):
            worker.step_started.emit(0)
            worker.log_line.emit("step 0 running")
            worker.progress.emit(50)
            # _refresh_stepper() (triggered by step_done/all_done below)
            # recomputes the bar from real step counts, so snapshot the
            # value the direct `progress` signal produced first.
            progress_snapshot["value"] = win._overall_progress.value()
            worker.step_done.emit(0, "result-0")
            worker.all_done.emit()

        fake_cls = _fake_worker_cls(_on_start)
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.pipeline_worker.PipelineWorker",
            fake_cls,
        )
        win._on_run_step()
        assert "step 0 running" in win._log_text.toPlainText()
        assert progress_snapshot["value"] == 50
        assert win._btn_run.isEnabled()  # re-enabled by all_done
        assert not win._btn_stop.isEnabled()

    def test_step_error_signal_handled(self, win, monkeypatch):
        win.set_input_sites(["s1"])

        def _on_start(worker):
            worker.step_started.emit(0)
            worker.step_error.emit(0, "boom")

        fake_cls = _fake_worker_cls(_on_start)
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.pipeline_worker.PipelineWorker",
            fake_cls,
        )
        win._on_run_step()  # must not raise

    def test_step_skipped_signal_handled(self, win, monkeypatch):
        win.set_input_sites(["s1"])

        def _on_start(worker):
            worker.step_skipped.emit(0)

        fake_cls = _fake_worker_cls(_on_start)
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.pipeline_worker.PipelineWorker",
            fake_cls,
        )
        win._on_run_step()  # must not raise

    def test_quick_pipeline_preloads_step0_and_starts_worker(
        self, win, monkeypatch
    ):
        win._ctrl._sites_input = ["s1", "s2", "s3"]
        started = []
        fake_cls = _fake_worker_cls(lambda w: started.append(w))
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.pipeline_worker.PipelineWorker",
            fake_cls,
        )
        win._on_quick_pipeline()
        assert len(started) == 1
        assert win._ctrl.steps[0].status == StepStatus.DONE

    def test_skip_step_marks_skipped(self, win):
        win.set_input_sites(["s1"])
        win._select_step(1)
        win._on_skip_step()
        assert win._ctrl.steps[1].status == StepStatus.SKIPPED
        assert "skipped" in win._log_text.toPlainText().lower()

    def test_stop_requests_interruption_when_running(self, win, monkeypatch):
        fake_cls = _fake_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.pipeline_worker.PipelineWorker",
            fake_cls,
        )
        win._worker = fake_cls(win._ctrl, [0])
        win._worker._running = True
        win._on_stop()
        assert "Stop requested" in win._log_text.toPlainText()

    def test_stop_noop_when_not_running(self, win):
        win._on_stop()  # no worker at all -> no-op

    def test_reset_clears_steps(self, win):
        win.set_input_sites(["s1"])
        win._ctrl.steps[1].status = StepStatus.DONE
        win._on_reset()
        assert win._ctrl.steps[1].status == StepStatus.PENDING
        # step 0 pre-done again since sites are still available
        assert win._ctrl.steps[0].status == StepStatus.DONE

    def test_reset_noop_while_worker_running(self, win, monkeypatch):
        fake_cls = _fake_worker_cls()
        win._worker = fake_cls(win._ctrl, [0])
        win._worker._running = True
        win._ctrl.steps[1].status = StepStatus.DONE
        win._on_reset()
        assert win._ctrl.steps[1].status == StepStatus.DONE  # untouched

    def test_reset_without_sites_just_refreshes(self, win):
        win._on_reset()  # no sites input at all


# ── Worker slot methods (direct) ─────────────────────────────────────────────


class TestWorkerSlots:
    def test_on_step_started_updates_status_and_progress_bar(self, win):
        win._on_step_started(0)
        assert win._ctrl.steps[0].status == StepStatus.RUNNING
        assert win._step_progress.isVisible()

    def test_on_step_done_advances_stepper_selection(self, win):
        win._select_step(0)
        win._on_step_done(0, "res")
        assert win._step_list.currentRow() == 1

    def test_on_step_done_no_advance_past_last_step(self, win):
        last = len(win._ctrl.steps) - 1
        win._select_step(last)
        win._on_step_done(last, "res")  # must not raise

    def test_on_step_error_refreshes(self, win):
        win._select_step(0)
        win._on_step_error(0, "boom")  # must not raise

    def test_on_step_skipped_refreshes(self, win):
        win._on_step_skipped(0)  # must not raise

    def test_on_log_line_appends(self, win):
        win._on_log_line("hello")
        assert "hello" in win._log_text.toPlainText()

    def test_on_all_done_emits_pipeline_finished_signal(self, win):
        win._ctrl._sites_chain[0] = ["out-sites"]
        received = []
        win.pipeline_finished.connect(lambda s: received.append(s))
        win._on_all_done()
        assert received == [["out-sites"]]
        assert win._btn_run.isEnabled()
        assert not win._btn_stop.isEnabled()

    def test_on_all_done_no_output_sites_no_signal(self, win):
        received = []
        win.pipeline_finished.connect(lambda s: received.append(s))
        win._on_all_done()
        assert received == []


# ── Preview / full view ──────────────────────────────────────────────────────


class TestPreview:
    def test_draw_step_preview_noop_without_diag_fn(self, win):
        step = win._ctrl.steps[0]
        step.diag_fn = None
        win._draw_step_preview(step)  # must not raise

    def test_draw_step_preview_noop_without_output_sites(self, win):
        step = win._ctrl.steps[1]
        step.output_sites = None
        win._draw_step_preview(step)  # must not raise

    def test_draw_step_preview_unknown_fn_noop(self, win):
        step = win._ctrl.steps[1]
        step.diag_fn = "definitely_not_a_real_emtools_fn"
        step.output_sites = object()
        win._draw_step_preview(step)  # must not raise

    def test_draw_step_preview_fn_raises_is_swallowed(self, win, monkeypatch):
        import pycsamt.emtools as et

        def _boom(*a, **k):
            raise RuntimeError("plot failure")

        monkeypatch.setattr(et, "fake_diag_fn", _boom, raising=False)
        step = win._ctrl.steps[1]
        step.diag_fn = "fake_diag_fn"
        step.output_sites = object()
        win._draw_step_preview(step)  # must not raise

    def test_on_full_view_no_diag_fn_closes_figure(self, win):
        step = win._ctrl.steps[0]
        step.diag_fn = None
        win._on_full_view()  # must not raise

    def test_on_full_view_fn_raises_closes_figure(self, win, monkeypatch):
        import pycsamt.emtools as et

        def _boom(*a, **k):
            raise RuntimeError("boom")

        monkeypatch.setattr(et, "fake_diag_fn2", _boom, raising=False)
        step = win._ctrl.steps[1]
        step.diag_fn = "fake_diag_fn2"
        step.output_sites = object()
        win._selected_step = 1
        win._on_full_view()  # must not raise

    def test_on_full_view_success_shows_plot(self, win, monkeypatch):
        import matplotlib.pyplot as plt

        import pycsamt.emtools as et

        calls = []
        monkeypatch.setattr(
            et, "fake_diag_fn3", lambda *a, **k: calls.append(1), raising=False
        )
        monkeypatch.setattr(plt, "show", lambda: calls.append("shown"))
        step = win._ctrl.steps[1]
        step.diag_fn = "fake_diag_fn3"
        step.output_sites = object()
        win._selected_step = 1
        win._on_full_view()
        assert "shown" in calls


# ── Save / load config ───────────────────────────────────────────────────────


class TestSaveLoadConfig:
    def test_save_config_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        win._on_save_config()  # must not raise

    def test_save_and_load_config_roundtrip(self, win, monkeypatch, tmp_path):
        cfg_path = str(tmp_path / "pipeline_cfg.json")
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (cfg_path, "")),
        )
        win._on_save_config()
        assert "saved" in win._log_text.toPlainText()

        with open(cfg_path, encoding="utf-8") as f:
            saved = json.load(f)
        assert isinstance(saved, dict)
        assert "steps" in saved

        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: (cfg_path, "")),
        )
        win._on_load_config()
        assert "loaded" in win._log_text.toPlainText()

    def test_load_config_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        win._on_load_config()  # must not raise


# ── Summary ──────────────────────────────────────────────────────────────────


class TestSummary:
    def test_update_summary_builds_html(self, win):
        win._ctrl.steps[0].status = StepStatus.DONE
        win._ctrl.steps[0].result_info = "5 stations"
        win._ctrl.steps[0].elapsed_s = 0.5
        win._update_summary()
        html = win._summary_browser.toHtml()
        assert "5 stations" in html


# ── Geometry persistence ─────────────────────────────────────────────────────


class TestGeometry:
    def test_save_and_restore_geometry(self, win):
        store = {}
        win.save_geometry_to(store)
        assert "pipeline_window" in store
        win.restore_geometry_from(store)  # must not raise

    def test_restore_geometry_missing_entry_noop(self, win):
        win.restore_geometry_from({})  # must not raise

    def test_restore_geometry_corrupt_base64_swallowed(self, win):
        win.restore_geometry_from(
            {"pipeline_window": {"geometry": "not-valid-base64!!"}}
        )  # must not raise

    def test_close_event_hides_instead_of_destroying(self, win):
        win.show()
        win.close()
        assert not win.isVisible()
