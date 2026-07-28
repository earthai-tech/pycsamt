# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for AdvancedToolsWindow (pycsamt.app.desktop.windows.advanced_window).

Strategy
--------
* ``DimModelWorker`` / ``ConversionWorker`` are genuine ``QThread``
  subclasses; run-triggering tests monkeypatch their ``.start()`` method
  directly (they're imported as class attributes of the window module,
  not re-imported locally, so patching the class itself is sufficient
  here) instead of letting a real thread spin up under the offscreen
  platform.
* ``QFileDialog`` / ``QColorDialog`` are monkeypatched to avoid real
  native dialogs.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtGui import QColor
from PySide6.QtWidgets import QColorDialog, QFileDialog

from pycsamt.app.desktop.controllers.advanced_controller import (
    ADVANCED_GROUPS,
    CONV_INDEX,
    TOPO_INDEX,
    ConversionWorker,
    DimModelWorker,
)
from pycsamt.app.desktop.windows.advanced_window import AdvancedToolsWindow


@pytest.fixture
def win(qapp):
    w = AdvancedToolsWindow(parent=None)
    w.show()
    yield w
    w.close()


def _select_category(win, label):
    for row, (grp_label, _plots) in enumerate(ADVANCED_GROUPS):
        if grp_label == label:
            win._combo_category.setCurrentIndex(row)
            return row
    raise AssertionError(f"category {label!r} not found")


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "Advanced Tools" in win.windowTitle()

    def test_category_combo_populated(self, win):
        assert win._combo_category.count() == len(ADVANCED_GROUPS)

    def test_default_page_is_plot_page(self, win):
        assert win._params_stack.currentIndex() == 0
        assert win._content_stack.currentIndex() == 0

    def test_model_group_hidden_by_default(self, win):
        assert not win._grp_model.isVisible()

    def test_topo_file_row_hidden_by_default(self, win):
        assert not win._topo_file_row.isVisible()


# ── Category switching ───────────────────────────────────────────────────────


class TestCategorySwitching:
    def test_switch_to_topo_page(self, win):
        win._combo_category.setCurrentIndex(TOPO_INDEX)
        assert win._params_stack.currentIndex() == 1
        assert win._content_stack.currentIndex() == 1

    def test_switch_to_conversion_page(self, win):
        win._combo_category.setCurrentIndex(CONV_INDEX)
        assert win._params_stack.currentIndex() == 2
        assert win._content_stack.currentIndex() == 2

    def test_switch_to_plot_category_populates_combo(self, win):
        win._combo_category.setCurrentIndex(CONV_INDEX)
        row = _select_category(win, "Strike Analysis")
        _label, plots = ADVANCED_GROUPS[row]
        assert win._combo_plot.count() == len(plots)
        assert win._params_stack.currentIndex() == 0

    def test_on_category_changed_out_of_range_noop(self, win):
        win._on_category_changed(-1)
        win._on_category_changed(999)

    def test_atom_psection_shows_model_group(self, win):
        row = _select_category(win, "Depth Imaging")
        _label, plots = ADVANCED_GROUPS[row]
        for i, (_lbl, fn_name, _has_ax) in enumerate(plots):
            if fn_name == "plot_atom_psection":
                win._combo_plot.setCurrentIndex(i)
                assert win._grp_model.isVisible()
                return
        pytest.skip("plot_atom_psection not present in Depth Imaging group")

    def test_plot_changed_updates_description(self, win):
        row = _select_category(win, "Strike Analysis")
        _label, plots = ADVANCED_GROUPS[row]
        if len(plots) > 1:
            win._combo_plot.setCurrentIndex(1)
        assert win._desc_lbl.text() != ""

    def test_update_desc_empty_plots_list(self, win):
        win._update_desc(ADVANCED_GROUPS.index(("Topography", [])), 0)
        assert win._desc_lbl.text() == ""

    def test_update_desc_bad_index_swallowed(self, win):
        win._update_desc(999, 999)
        assert win._desc_lbl.text() == ""


# ── Run / Export ──────────────────────────────────────────────────────────────


class TestRunExport:
    def test_run_no_sites(self, win):
        _select_category(win, "Strike Analysis")
        win._on_run()
        assert win._status_lbl.text() == "Load survey data first."

    def test_run_no_category_or_plot_noop(self, win):
        win._combo_category.setCurrentIndex(-1)
        win._on_run()  # guarded early return

    def test_run_success_with_figure(self, win, monkeypatch):
        _select_category(win, "Strike Analysis")
        win._ctrl._sites = object()
        monkeypatch.setattr(win._ctrl, "draw", lambda fn, has_ax, fig: None)
        win._on_run()
        assert win._status_lbl.text() == "Done."
        assert win._btn_run.isEnabled()

    def test_run_success_replaces_figure(self, win, monkeypatch):
        import matplotlib.figure

        _select_category(win, "Strike Analysis")
        win._ctrl._sites = object()
        new_fig = matplotlib.figure.Figure()
        monkeypatch.setattr(win._ctrl, "draw", lambda fn, has_ax, fig: new_fig)
        win._on_run()
        assert win._status_lbl.text() == "Done."

    def test_run_exception_reported(self, win, monkeypatch):
        _select_category(win, "Strike Analysis")
        win._ctrl._sites = object()

        def _boom(fn, has_ax, fig):
            raise RuntimeError("draw boom")

        monkeypatch.setattr(win._ctrl, "draw", _boom)
        win._on_run()
        assert "Error: draw boom" in win._status_lbl.text()
        assert win._btn_run.isEnabled()

    def test_run_plot_row_out_of_range_noop(self, win):
        _select_category(win, "Strike Analysis")
        win._combo_plot.blockSignals(True)
        win._combo_plot.setCurrentIndex(-1)
        win._combo_plot.blockSignals(False)
        # currentIndex may clamp to 0; force plot_row invalid via direct call
        win._on_run()  # must not raise regardless

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


class TestAutoRender:
    def test_auto_render_noop_when_already_rendered(self, win):
        win._auto_rendered = True
        win._ctrl._sites = object()
        win._auto_render_if_ready()  # no crash, still marked rendered

    def test_auto_render_noop_without_sites(self, win):
        win._auto_rendered = False
        win._ctrl._sites = None
        win._auto_render_if_ready()
        assert win._auto_rendered is False

    def test_auto_render_noop_when_not_visible(self, win, monkeypatch):
        win._auto_rendered = False
        win._ctrl._sites = object()
        monkeypatch.setattr(win, "isVisible", lambda: False)
        win._auto_render_if_ready()
        assert win._auto_rendered is False

    def test_auto_render_noop_for_topo_category(self, win):
        win._combo_category.setCurrentIndex(TOPO_INDEX)
        win._auto_rendered = False
        win._ctrl._sites = object()
        win._auto_render_if_ready()  # ADVANCED_GROUPS[TOPO_INDEX] has no plots

    def test_auto_render_triggers_timer(self, win, monkeypatch):
        _select_category(win, "Strike Analysis")
        win._auto_rendered = False
        win._ctrl._sites = object()
        calls = []
        monkeypatch.setattr(
            "pycsamt.app.desktop.windows.advanced_window.QTimer.singleShot",
            staticmethod(lambda ms, fn: calls.append(fn)),
        )
        win._auto_render_if_ready()
        assert win._auto_rendered is True
        assert len(calls) == 1

    def test_show_event_triggers_auto_render(self, win, monkeypatch):
        from PySide6.QtGui import QShowEvent

        _select_category(win, "Strike Analysis")
        win._ctrl._sites = object()
        win._auto_rendered = False
        calls = []
        monkeypatch.setattr(win, "_auto_render_if_ready", lambda: calls.append(1))
        win.showEvent(QShowEvent())
        assert calls == [1]


# ── Train model ───────────────────────────────────────────────────────────────


class TestTrainModel:
    def test_train_no_sites(self, win):
        win._on_train_model()
        assert win._model_status_lbl.text() == "Load survey data first."

    def test_train_starts_worker(self, win, monkeypatch):
        win._ctrl._sites = object()
        started = []
        monkeypatch.setattr(DimModelWorker, "start", lambda self: started.append(self))
        win._on_train_model()
        assert len(started) == 1
        assert not win._btn_train_model.isEnabled()

    def test_on_model_trained_updates_labels(self, win):
        import numpy as np

        model = {"D": np.zeros((10, 6)), "meta": {"samples": 42}}
        win._on_model_trained(model)
        assert "6 atoms" in win._model_status_lbl.text()
        assert "42 samples" in win._model_status_lbl.text()
        assert win._btn_train_model.isEnabled()

    def test_on_model_trained_no_dictionary(self, win):
        win._on_model_trained({"meta": {}})
        assert "? atoms" in win._model_status_lbl.text()

    def test_on_model_train_error(self, win):
        win._on_model_train_error("training failed")
        assert "Error: training failed" in win._model_status_lbl.text()
        assert win._btn_train_model.isEnabled()


# ── Topo slots ────────────────────────────────────────────────────────────────


class TestTopoSlots:
    def test_source_changed_shows_file_row(self, win):
        win._combo_category.setCurrentIndex(TOPO_INDEX)  # realize the page
        win._combo_topo_source.setCurrentText("file")
        assert win._topo_file_row.isVisible()
        win._combo_topo_source.setCurrentText("sites")
        assert not win._topo_file_row.isVisible()

    def test_browse_topo_file_sets_text(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("/fake/elev.csv", "")),
        )
        win._browse_topo_file()
        assert win._edit_topo_file.text() == "/fake/elev.csv"

    def test_browse_topo_file_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        win._browse_topo_file()
        assert win._edit_topo_file.text() == ""

    def test_pick_fill_color_valid(self, win, monkeypatch):
        monkeypatch.setattr(
            QColorDialog,
            "getColor",
            staticmethod(lambda *a, **k: QColor("#ff0000")),
        )
        win._pick_color("fill")
        assert win._topo_fill_color == "#ff0000"

    def test_pick_line_color_valid(self, win, monkeypatch):
        monkeypatch.setattr(
            QColorDialog,
            "getColor",
            staticmethod(lambda *a, **k: QColor("#00ff00")),
        )
        win._pick_color("line")
        assert win._topo_line_color == "#00ff00"

    def test_pick_color_invalid_noop(self, win, monkeypatch):
        before = win._topo_fill_color
        monkeypatch.setattr(
            QColorDialog,
            "getColor",
            staticmethod(lambda *a, **k: QColor()),  # invalid
        )
        win._pick_color("fill")
        assert win._topo_fill_color == before

    def test_on_topo_apply_success(self, win, monkeypatch):
        import pycsamt.topo.config as topo_cfg

        monkeypatch.setattr(topo_cfg, "configure_topo", lambda **kw: None)

        class _FakeSummary:
            def summary(self):
                return "topo summary"

        monkeypatch.setattr(topo_cfg, "PYCSAMT_TOPO", _FakeSummary())
        win._on_topo_apply()
        assert win._topo_status.text() == "topo summary"

    def test_on_topo_apply_error(self, win, monkeypatch):
        import pycsamt.topo.config as topo_cfg

        def _boom(**kw):
            raise RuntimeError("bad config")

        monkeypatch.setattr(topo_cfg, "configure_topo", _boom)
        win._on_topo_apply()
        assert "Error: bad config" in win._topo_status.text()

    def test_on_topo_reset_success(self, win, monkeypatch):
        import pycsamt.topo.config as topo_cfg

        monkeypatch.setattr(topo_cfg, "reset_topo", lambda: None)

        class _FakeSummary:
            def summary(self):
                return "reset summary"

        monkeypatch.setattr(topo_cfg, "PYCSAMT_TOPO", _FakeSummary())
        win._on_topo_reset()
        assert "reset summary" in win._topo_status.text()

    def test_on_topo_reset_error(self, win, monkeypatch):
        import pycsamt.topo.config as topo_cfg

        def _boom():
            raise RuntimeError("reset failed")

        monkeypatch.setattr(topo_cfg, "reset_topo", _boom)
        win._on_topo_reset()
        assert "Error: reset failed" in win._topo_status.text()

    def test_sync_topo_widgets_from_config(self, win, monkeypatch):
        from types import SimpleNamespace

        import pycsamt.topo.config as topo_cfg

        fake_t = SimpleNamespace(
            enabled=True,
            source="file",
            elev_file="/some/file.csv",
            interp_method="cubic",
            exaggeration=2.0,
            fill_color="#111111",
            fill_alpha=0.6,
            line_color="#222222",
            line_width=2.0,
            show_surface_line=False,
            clip_below_surface=False,
            station_pins_at_surface=False,
            show_topo_strip=False,
            strip_height_ratio=0.25,
        )
        monkeypatch.setattr(topo_cfg, "PYCSAMT_TOPO", fake_t)
        win._sync_topo_widgets_from_config()
        assert win._chk_topo_enabled.isChecked()
        assert win._combo_topo_source.currentText() == "file"
        assert win._spin_topo_exag.value() == pytest.approx(2.0)
        assert win._topo_fill_color == "#111111"

    def test_sync_topo_widgets_exception_swallowed(self, win, monkeypatch):
        import pycsamt.topo.config as topo_cfg

        monkeypatch.setattr(
            topo_cfg,
            "PYCSAMT_TOPO",
            property(lambda self: (_ for _ in ()).throw(RuntimeError())),
        )
        win._sync_topo_widgets_from_config()  # must not raise

    def test_refresh_topo_preview_no_data(self, win):
        win._topo_ctrl._sites = None
        win._refresh_topo_preview()
        assert win._topo_stats_lbl.text() == "No data loaded"

    def test_refresh_topo_preview_with_stats(self, win, monkeypatch):
        win._topo_ctrl._sites = object()
        monkeypatch.setattr(
            win._topo_ctrl,
            "get_stats",
            lambda: {"n_stations": 5, "elev_min": 100.0, "elev_max": 500.0},
        )
        monkeypatch.setattr(win._topo_ctrl, "plot_elevation_profile", lambda fig: None)
        win._combo_topo_view.setCurrentIndex(0)
        win._refresh_topo_preview()
        assert "5 stations" in win._topo_stats_lbl.text()

    def test_refresh_topo_preview_stats_exception(self, win, monkeypatch):
        win._topo_ctrl._sites = object()

        def _boom():
            raise RuntimeError("stats failed")

        monkeypatch.setattr(win._topo_ctrl, "get_stats", _boom)
        monkeypatch.setattr(win._topo_ctrl, "plot_elevation_profile", lambda fig: None)
        win._refresh_topo_preview()
        assert win._topo_stats_lbl.text() == ""

    def test_refresh_topo_preview_plot_error_shows_message(self, win, monkeypatch):
        win._topo_ctrl._sites = object()
        monkeypatch.setattr(
            win._topo_ctrl,
            "get_stats",
            lambda: {"n_stations": 0, "elev_min": 0, "elev_max": 0},
        )

        def _boom(fig):
            raise RuntimeError("plot failed")

        monkeypatch.setattr(win._topo_ctrl, "plot_elevation_profile", _boom)
        win._combo_topo_view.setCurrentIndex(0)
        win._refresh_topo_preview()
        assert len(win._canvas_topo.figure.axes) >= 1

    def test_refresh_topo_preview_view_1_and_2(self, win, monkeypatch):
        win._topo_ctrl._sites = object()
        monkeypatch.setattr(
            win._topo_ctrl,
            "get_stats",
            lambda: {"n_stations": 1, "elev_min": 1, "elev_max": 2},
        )
        monkeypatch.setattr(win._topo_ctrl, "plot_fill_preview", lambda fig: None)
        monkeypatch.setattr(
            win._topo_ctrl, "plot_elevation_histogram", lambda fig: None
        )
        win._combo_topo_view.setCurrentIndex(1)
        win._refresh_topo_preview()
        win._combo_topo_view.setCurrentIndex(2)
        win._refresh_topo_preview()


# ── Conversion slots ──────────────────────────────────────────────────────────


class TestConversionSlots:
    def test_conv_type_changed_shows_avg_topo_group(self, win):
        win._combo_category.setCurrentIndex(CONV_INDEX)  # realize the page
        win._combo_conv_type.setCurrentIndex(1)  # J
        assert not win._avg_topo_group.isVisible()
        win._combo_conv_type.setCurrentIndex(0)  # AVG
        assert win._avg_topo_group.isVisible()

    def test_update_conv_run_state_disabled_without_path(self, win):
        win._edit_conv_path.setText("")
        assert not win._btn_conv_run.isEnabled()

    def test_update_conv_run_state_enabled_with_path(self, win):
        win._edit_conv_path.setText("/some/path")
        assert win._btn_conv_run.isEnabled()

    def test_browse_conv_path_file_selected(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("/some/file.avg", "")),
        )
        win._browse_conv_path()
        assert win._edit_conv_path.text() == "/some/file.avg"

    def test_browse_conv_path_falls_back_to_directory(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/some/dir"),
        )
        win._browse_conv_path()
        assert win._edit_conv_path.text() == "/some/dir"

    def test_browse_out_dir(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/out/dir"),
        )
        win._browse_out_dir()
        assert win._edit_out_dir.text() == "/out/dir"

    def test_browse_avg_stn_path(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("/some/k1.stn", "")),
        )
        win._browse_avg_stn_path()
        assert win._avg_stn_path.text() == "/some/k1.stn"

    def test_conv_run_no_path(self, win):
        win._edit_conv_path.setText("")
        win._on_conv_run()
        assert "provide an input" in win._conv_status.text()

    def test_conv_run_avg_options_and_worker_start(self, win, monkeypatch):
        win._edit_conv_path.setText("/some/path.avg")
        win._combo_conv_type.setCurrentIndex(0)
        win._avg_compute_z.setChecked(True)
        win._avg_stn_path.setText("/some/k1.stn")
        win._avg_epsg.setText("32650")
        win._chk_write_edis.setChecked(True)
        win._edit_out_dir.setText("/out")
        started = []
        monkeypatch.setattr(
            ConversionWorker, "start", lambda self: started.append(self)
        )
        win._on_conv_run()
        assert len(started) == 1
        assert win._conv_running is True

    def test_conv_run_j_options(self, win, monkeypatch):
        win._edit_conv_path.setText("/some/path.j")
        win._combo_conv_type.setCurrentIndex(1)
        win._j_station_name.setText("STA01")
        monkeypatch.setattr(ConversionWorker, "start", lambda self: None)
        win._on_conv_run()

    def test_conv_run_spectra_options(self, win, monkeypatch):
        win._edit_conv_path.setText("/some/path")
        win._combo_conv_type.setCurrentIndex(2)
        win._sp_estimate_errors.setChecked(True)
        win._sp_remote_ref.setChecked(True)
        win._sp_station_suffix.setText("_A")
        monkeypatch.setattr(ConversionWorker, "start", lambda self: None)
        win._on_conv_run()

    def test_on_conv_finished_populates_table_and_plots(self, win, monkeypatch):
        rows = [
            {
                "station": "S1",
                "n_freqs": 10,
                "f_min": 0.01,
                "f_max": 100.0,
                "lat": 1.234567,
                "lon": 2.345678,
                "elev": 500.0,
                "has_Z": True,
                "has_tipper": False,
            },
            {
                "station": "S2",
                "n_freqs": 5,
                "f_min": float("nan"),
                "f_max": float("nan"),
                "lat": None,
                "lon": None,
                "elev": None,
                "has_Z": False,
                "has_tipper": True,
            },
        ]
        monkeypatch.setattr(
            win._conv_ctrl,
            "build_stats",
            lambda collection, failures: {
                "rows": rows,
                "n_total": 2,
                "n_failures": 1,
            },
        )
        monkeypatch.setattr(win._conv_ctrl, "plot_impedance_curves", lambda fig: None)
        monkeypatch.setattr(win._conv_ctrl, "plot_station_map", lambda fig: None)
        win._on_conv_finished(object(), ["failure1"])
        assert win._conv_table.rowCount() == 2
        assert "2 stations" in win._conv_status.text()
        assert "1 failures" in win._conv_status.text()
        assert win._btn_conv_commit.isEnabled()
        assert win._btn_conv_export.isEnabled()

    def test_on_conv_error(self, win):
        win._on_conv_error("conversion failed")
        assert "Error: conversion failed" in win._conv_status.text()
        assert not win._conv_running
        assert not win._conv_progress.isVisible()

    def test_on_conv_commit_no_result_noop(self, win):
        received = []
        win.conversion_committed.connect(lambda r: received.append(r))
        win._on_conv_commit()
        assert received == []

    def test_on_conv_commit_emits_signal(self, win):
        win._conv_ctrl._result = "the-collection"
        received = []
        win.conversion_committed.connect(lambda r: received.append(r))
        win._on_conv_commit()
        assert received == ["the-collection"]

    def test_on_conv_export_no_result_noop(self, win):
        win._on_conv_export()  # has_result False -> early return

    def test_on_conv_export_cancelled(self, win, monkeypatch):
        win._conv_ctrl._result = ["fake"]
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: ""),
        )
        win._on_conv_export()  # must not raise

    def test_on_conv_export_success(self, win, monkeypatch):
        class _FakeEd:
            def write(self, save_dir):
                pass

        win._conv_ctrl._result = [_FakeEd(), _FakeEd()]
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/out"),
        )
        win._on_conv_export()
        assert "Exported 2 EDI files" in win._conv_status.text()

    def test_on_conv_export_per_item_failure_swallowed(self, win, monkeypatch):
        class _BadEd:
            def write(self, save_dir):
                raise RuntimeError("write failed")

        win._conv_ctrl._result = [_BadEd()]
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/out"),
        )
        win._on_conv_export()
        assert "Exported 0 EDI files" in win._conv_status.text()

    def test_on_conv_export_iteration_error_reported(self, win, monkeypatch):
        win._conv_ctrl._result = ["not-empty"]
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/out"),
        )
        import pycsamt.emtools._core as core_mod

        def _boom(_x):
            raise RuntimeError("iter failed")

        monkeypatch.setattr(core_mod, "_iter_items", _boom)
        win._on_conv_export()
        assert "Export error" in win._conv_status.text()

    def test_on_conv_clear_resets_ui(self, win):
        win._conv_ctrl._result = ["x"]
        win._conv_table.setColumnCount(1)
        win._conv_table.setRowCount(1)
        win._btn_conv_commit.setEnabled(True)
        win._btn_conv_export.setEnabled(True)
        win._on_conv_clear()
        assert win._conv_ctrl._result is None
        assert win._conv_table.rowCount() == 0
        assert not win._btn_conv_commit.isEnabled()
        assert not win._btn_conv_export.isEnabled()
        assert win._conv_status.text() == "Cleared."


# ── set_sites / set_dark_mode ─────────────────────────────────────────────────


class TestPublicApi:
    def test_set_sites_delegates_to_all_controllers(self, win):
        sites = object()
        win.set_sites(sites)
        assert win._ctrl._sites is sites
        assert win._topo_ctrl._sites is sites

    def test_set_sites_refreshes_topo_when_on_topo_page(self, win, monkeypatch):
        win._combo_category.setCurrentIndex(TOPO_INDEX)
        calls = []
        monkeypatch.setattr(win, "_refresh_topo_preview", lambda: calls.append(1))
        win.set_sites(object())
        assert calls == [1]

    def test_set_dark_mode_delegates(self, win):
        win.set_dark_mode(False)
        assert win._ctrl.dark is False
        assert win._topo_ctrl.dark is False
        assert win._conv_ctrl.dark is False
