# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for PipelineWindow — Pipeline Studio.

Strategy
--------
* The real WorkflowController and library pipeline engine run on the small
  WILLY L18 survey.
* ``WorkflowWorker.start`` is patched to call ``run()`` synchronously, so no
  QThread is started under the offscreen platform (hang risk) and worker
  signals are delivered in order -- which also reproduces the ordering the
  status-race fix relies on.
* Run history is redirected to ``tmp_path``; nothing touches ~/.pycsamt.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QFileDialog, QMessageBox

import pycsamt.app.desktop.workers.workflow_worker as ww_mod
from pycsamt.app.desktop.controllers.workflow_controller import RunStatus
from pycsamt.app.desktop.windows.pipeline_window import (
    PipelineWindow,
    _category_title,
)

_ROOT = Path(__file__).parents[4]
_WILLY = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"


@pytest.fixture(scope="module")
def sites():
    if not (_WILLY.exists() and any(_WILLY.glob("*.edi"))):
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY))


@pytest.fixture(autouse=True)
def _sync_worker_and_tmp_history(monkeypatch, tmp_path):
    monkeypatch.setattr(ww_mod.WorkflowWorker, "start",
                        lambda self: self.run())
    from pycsamt.api.pipe import PYCSAMT_PIPE

    monkeypatch.setattr(PYCSAMT_PIPE, "history_path",
                        str(tmp_path / "history.jsonl"))


@pytest.fixture
def win(qapp):
    w = PipelineWindow()
    w.show()
    yield w
    w.close()


def _pills(win):
    return [r.pill.status for r in win._rows]


class TestConstruction:
    def test_opens_with_default_preset(self, win):
        assert win.windowTitle() == "pycsamt — Pipeline Studio"
        assert [s.code for s in win._ctrl.steps][0] == "NR001"
        assert len(win._rows) == len(win._ctrl.steps)

    def test_library_lists_registry_by_category(self, win):
        n = sum(win._lib_tree.topLevelItem(i).childCount()
                for i in range(win._lib_tree.topLevelItemCount()))
        assert n >= 55

    def test_category_titles(self):
        assert _category_title("noise_removal") == "Noise Removal"
        assert _category_title("qc") == "QC"

    def test_run_disabled_without_input(self, win):
        assert not win._btn_run_all.isEnabled()
        assert not win._btn_apply_main.isEnabled()


class TestEditing:
    def test_add_from_library_inserts_after_selection(self, win):
        win._step_list.setCurrentRow(0)
        parent = win._lib_tree.topLevelItem(0)
        win._lib_tree.setCurrentItem(parent.child(0))
        code = parent.child(0).data(0, Qt_UserRole())
        win._on_add_from_library()
        assert win._ctrl.steps[1].code == code
        assert win._step_list.currentRow() == 1

    def test_move_and_remove(self, win):
        first = win._ctrl.steps[0].code
        win._step_list.setCurrentRow(0)
        win._on_move(+1)
        assert win._ctrl.steps[1].code == first
        n = len(win._ctrl.steps)
        win._on_remove()
        assert len(win._ctrl.steps) == n - 1

    def test_library_search_filters(self, win):
        win._lib_search.setText("static")
        shown = [
            win._lib_tree.topLevelItem(i).text(0)
            for i in range(win._lib_tree.topLevelItemCount())
            if not win._lib_tree.topLevelItem(i).isHidden()
        ]
        # the Static Shift group, plus NR013 "Spatial Window Static Shift
        # Correction" filed under Noise Removal -- a genuine match
        assert "Static Shift" in shown and "Noise Removal" in shown
        assert "Tensor" not in shown and "Export" not in shown

    def test_library_drawer_toggle_and_session(self, win):
        win._toggle_library()
        assert not win.library_visible and not win._lib_panel.isVisible()
        store = {}
        win.save_geometry_to(store)
        win._set_library_visible(True)
        win.restore_geometry_from(store)
        assert not win.library_visible

    def test_apply_preset_confirms(self, win, monkeypatch):
        monkeypatch.setattr(QMessageBox, "question",
                            staticmethod(lambda *a, **k:
                                         QMessageBox.StandardButton.Yes))
        win._combo_preset.setCurrentIndex(
            win._combo_preset.findData("basic_qc"))
        win._on_apply_preset()
        assert [s.code for s in win._ctrl.steps][0] == "NR001"
        assert len(win._ctrl.steps) == 5


class TestRunning:
    def test_run_all_reaches_done_with_disabled_step(self, win, sites):
        win.set_input_sites(sites)
        win._rows[2].check.setChecked(False)
        win._start_run("all")
        statuses = _pills(win)
        assert statuses[2] is RunStatus.DISABLED
        assert all(s is RunStatus.DONE
                   for i, s in enumerate(statuses) if i != 2)
        assert win._btn_apply_main.isEnabled()
        assert "Run finished" in win._progress_lbl.text()

    def test_late_started_signal_cannot_overwrite_done(self, win, sites):
        """Regression: queued cross-thread 'started' arriving after the
        worker marked the step DONE used to leave it Running→Pending."""
        win.set_input_sites(sites)
        win._start_run("all")
        win._on_step_started(0)  # a stale, late-delivered signal
        assert win._ctrl.steps[0].status is RunStatus.DONE

    def test_apply_to_main_emits_output(self, win, sites):
        win.set_input_sites(sites)
        win._start_run("all")
        got = []
        win.pipeline_finished.connect(got.append)
        win._on_apply_main()
        assert got and len(got[0]) == len(sites)

    def test_run_step_needs_earlier_results(self, win, sites):
        win.set_input_sites(sites)
        win._step_list.setCurrentRow(3)
        win._start_run("one")
        assert win._ctrl.steps[3].status is RunStatus.PENDING
        assert "earlier steps" in win._log_text.toPlainText()

    def test_qc_before_after_figure_renders(self, win, sites):
        win.set_input_sites(sites)
        win._start_run("all")
        idx = next(i for i, s in enumerate(win._ctrl.steps)
                   if s.code == "SS001")
        win._step_list.setCurrentRow(idx)
        assert win._combo_qc.findText("plot_ss_delta_psection") >= 0
        win._on_show_qc()
        assert win._qc_view.showing_canvas

    def test_dashboard_drawn_natively(self, win, sites):
        win.set_input_sites(sites)
        win._start_run("all")
        fig = win._dash_view.canvas.figure
        assert win._dash_view.showing_canvas
        assert {ax.get_title() for ax in fig.axes} >= {
            "Step status", "Time per step", "Stations in → out"}

    def test_history_recorded_when_opted_in(self, win, sites, tmp_path):
        win.set_input_sites(sites)
        win._chk_history.setChecked(True)
        win._start_run("all")
        win._tabs.setCurrentIndex(3)
        assert win._hist_table.rowCount() == 1
        assert (tmp_path / "history.jsonl").exists()

    def test_reset_clears_results(self, win, sites):
        win.set_input_sites(sites)
        win._start_run("all")
        win._on_reset()
        assert all(s in (RunStatus.PENDING, RunStatus.DISABLED)
                   for s in _pills(win))
        assert not win._btn_apply_main.isEnabled()


class TestFiles:
    def test_save_and_open_workflow(self, win, tmp_path, monkeypatch):
        path = tmp_path / "wf.yaml"
        monkeypatch.setattr(QFileDialog, "getSaveFileName",
                            staticmethod(lambda *a, **k: (str(path), "")))
        win._on_save()
        assert path.exists()
        codes = [s.code for s in win._ctrl.steps]
        win._ctrl.clear()
        monkeypatch.setattr(QFileDialog, "getOpenFileName",
                            staticmethod(lambda *a, **k: (str(path), "")))
        win._on_open()
        assert [s.code for s in win._ctrl.steps] == codes

    def test_export_writes_outputs(self, win, sites, tmp_path, monkeypatch):
        win.set_input_sites(sites)
        win._start_run("all")
        out = tmp_path / "export"
        monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                            staticmethod(lambda *a, **k: str(out)))
        win._on_export()
        assert "Exported to" in win._log_text.toPlainText()
        assert any(out.rglob("*.edi"))


def Qt_UserRole():  # noqa: N802 -- tiny helper for readability above
    from PySide6.QtCore import Qt

    return Qt.ItemDataRole.UserRole
