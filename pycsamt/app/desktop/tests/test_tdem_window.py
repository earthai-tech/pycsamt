# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for TDEMWindow (pycsamt.app.desktop.windows.tdem_window).

Strategy mirrors test_advanced_window.py:
* QFileDialog is monkeypatched to avoid real native dialogs.
* TDEMController.draw() is monkeypatched directly for run/export tests so
  no real TDEM data or matplotlib rendering pipeline is required.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QFileDialog

from pycsamt.app.desktop.controllers.tdem_controller import TDEM_GROUPS
from pycsamt.app.desktop.windows.tdem_window import TDEMWindow


@pytest.fixture
def win(qapp):
    w = TDEMWindow(parent=None)
    w.show()
    yield w
    w.close()


def _select_category(win, label):
    for row, (grp_label, _plots) in enumerate(TDEM_GROUPS):
        if grp_label == label:
            win._combo_category.setCurrentIndex(row)
            return row
    raise AssertionError(f"category {label!r} not found")


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "TDEM Analysis" in win.windowTitle()

    def test_category_combo_populated(self, win):
        assert win._combo_category.count() == len(TDEM_GROUPS)

    def test_plot_combo_populated_for_first_category(self, win):
        _label, plots = TDEM_GROUPS[0]
        assert win._combo_plot.count() == len(plots)

    def test_load_progress_hidden_by_default(self, win):
        assert not win._load_progress.isVisible()

    def test_info_label_default_text(self, win):
        assert "No TDEM data loaded" in win._info_lbl.text()

    def test_desc_label_populated_on_init(self, win):
        assert win._desc_lbl.text() != ""


# ── Category switching ───────────────────────────────────────────────────────


class TestCategorySwitching:
    def test_switch_category_populates_plot_combo(self, win):
        row = _select_category(win, "Dashboard")
        _label, plots = TDEM_GROUPS[row]
        assert win._combo_plot.count() == len(plots)

    def test_on_category_changed_out_of_range_noop(self, win):
        win._on_category_changed(-1)
        win._on_category_changed(999)

    def test_switch_category_resets_plot_index(self, win):
        row = _select_category(win, "Survey Section")
        assert win._combo_plot.currentIndex() == 0

    def test_switch_category_updates_description(self, win):
        _select_category(win, "Map & Overview")
        assert win._desc_lbl.text() != ""

    def test_update_desc_bad_index_swallowed(self, win):
        win._update_desc(999, 999)
        assert win._desc_lbl.text() == ""


# ── Browse / load ─────────────────────────────────────────────────────────────


class TestBrowse:
    def test_browse_cancelled_noop(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: ""),
        )
        before = win._status_lbl.text()
        win._on_browse()
        assert win._status_lbl.text() == before

    def test_browse_success_updates_labels(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/fake/tdem_folder"),
        )
        monkeypatch.setattr(
            win._ctrl, "load_folder", lambda folder, progress_cb=None: True
        )
        win._ctrl.summary = "Folder: /fake/tdem_folder\nAVG files: 3   Z files: 3\nSoundings: 3"
        win._on_browse()
        assert win._info_lbl.text() == win._ctrl.summary
        assert "loaded" in win._status_lbl.text().lower()
        assert win._btn_run.isEnabled()
        assert not win._load_progress.isVisible()

    def test_browse_failure_updates_labels(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/fake/empty_folder"),
        )
        monkeypatch.setattr(
            win._ctrl, "load_folder", lambda folder, progress_cb=None: False
        )
        win._ctrl.summary = ""
        win._on_browse()
        assert win._info_lbl.text() == "No TDEM files found."
        assert win._status_lbl.text() == "Load failed."

    def test_browse_failure_uses_summary_when_present(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/fake/bad_folder"),
        )
        monkeypatch.setattr(
            win._ctrl, "load_folder", lambda folder, progress_cb=None: False
        )
        win._ctrl.summary = "Load error: boom"
        win._on_browse()
        assert win._info_lbl.text() == "Load error: boom"

    def test_browse_passes_progress_callback(self, win, monkeypatch):
        captured = {}

        def _fake_load(folder, progress_cb=None):
            captured["folder"] = folder
            captured["cb"] = progress_cb
            if progress_cb is not None:
                progress_cb(50)
            return True

        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/fake/folder"),
        )
        monkeypatch.setattr(win._ctrl, "load_folder", _fake_load)
        win._ctrl.summary = "ok"
        win._on_browse()
        assert captured["folder"] == "/fake/folder"
        assert captured["cb"] is not None


# ── Run / Export ──────────────────────────────────────────────────────────────


class TestRunExport:
    def test_run_no_category_or_plot_noop(self, win):
        win._combo_category.setCurrentIndex(-1)
        win._on_run()  # guarded early return, must not raise

    def test_run_success_draws_in_place(self, win, monkeypatch):
        _select_category(win, "Decay / Rho")
        monkeypatch.setattr(
            win._ctrl, "draw", lambda class_name, has_ax, data_key, fig: None
        )
        win._on_run()
        assert win._status_lbl.text() == "Done."
        assert win._btn_run.isEnabled()

    def test_run_success_replaces_figure(self, win, monkeypatch):
        import matplotlib.figure

        _select_category(win, "Decay / Rho")
        new_fig = matplotlib.figure.Figure()
        monkeypatch.setattr(
            win._ctrl,
            "draw",
            lambda class_name, has_ax, data_key, fig: new_fig,
        )
        win._on_run()
        assert win._status_lbl.text() == "Done."

    def test_run_exception_reported(self, win, monkeypatch):
        _select_category(win, "Decay / Rho")

        def _boom(class_name, has_ax, data_key, fig):
            raise RuntimeError("draw boom")

        monkeypatch.setattr(win._ctrl, "draw", _boom)
        win._on_run()
        assert "Error: draw boom" in win._status_lbl.text()
        assert win._btn_run.isEnabled()

    def test_run_plot_row_out_of_range_noop(self, win):
        row = _select_category(win, "Decay / Rho")
        _label, plots = TDEM_GROUPS[row]
        win._combo_plot.blockSignals(True)
        win._combo_plot.addItem("extra")
        win._combo_plot.setCurrentIndex(len(plots))  # stale index, out of range
        win._combo_plot.blockSignals(False)
        win._on_run()  # must not raise despite plot_row >= len(plots)

    def test_run_status_shows_class_name(self, win, monkeypatch):
        row = _select_category(win, "Decay / Rho")
        _label, plots = TDEM_GROUPS[row]
        _plot_label, class_name, _has_ax, _data_key = plots[0]
        seen = []

        def _record(cn, has_ax, data_key, fig):
            seen.append(cn)
            return None

        monkeypatch.setattr(win._ctrl, "draw", _record)
        win._on_run()
        assert seen == [class_name]

    def test_on_export_opens_dialog(self, win, monkeypatch):
        calls = []

        class _FakeExportDialog:
            def __init__(self, figure, parent):
                calls.append(figure)

            def exec(self):
                calls.append("exec")

        monkeypatch.setattr(
            "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog",
            _FakeExportDialog,
        )
        win._on_export()
        assert calls[-1] == "exec"


# ── Public API ────────────────────────────────────────────────────────────────


class TestPublicApi:
    def test_set_dark_mode_delegates(self, win):
        win.set_dark_mode(False)
        assert win._ctrl.dark is False
        win.set_dark_mode(True)
        assert win._ctrl.dark is True

    def test_set_sites_stores_on_base(self, win):
        sites = object()
        win.set_sites(sites)
        assert win._sites is sites
