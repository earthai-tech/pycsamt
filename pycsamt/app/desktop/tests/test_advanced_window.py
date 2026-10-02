# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the Advanced Tools studio
(pycsamt.app.desktop.windows.advanced_window).

Strategy
--------
* ``DimModelWorker`` / ``ConversionWorker`` are genuine ``QThread``
  subclasses; run-triggering tests monkeypatch their ``.start()`` method
  instead of letting a real thread spin up under the offscreen platform.
* ``QFileDialog`` / ``QColorDialog`` are monkeypatched to avoid real
  native dialogs.
* Plot tests draw real emtools figures from the bundled kap03 EDIs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtGui import QColor
from PySide6.QtWidgets import QColorDialog, QFileDialog, QLabel

from pycsamt.app.desktop.controllers.advanced_controller import (
    ConversionWorker,
    DimModelWorker,
)
from pycsamt.app.desktop.controllers.advanced_studio import (
    SECTIONS,
    line_groups,
    plot_inputs,
    plot_kwargs,
    plot_options,
)
from pycsamt.app.desktop.windows.advanced_window import AdvancedToolsWindow

TOPO = next(i for i, s in enumerate(SECTIONS) if s.key == "topo")
CONV = next(i for i, s in enumerate(SECTIONS) if s.key == "conv")
KAP03 = Path(__file__).resolve().parents[4] / "data" / "MT" / "kap03lmt_edis"


@pytest.fixture
def win(qapp):
    w = AdvancedToolsWindow(parent=None)
    w.show()
    yield w
    w.close()
    import matplotlib.pyplot as plt

    plt.close("all")


@pytest.fixture(scope="module")
def kap03():
    if not KAP03.is_dir():
        pytest.skip("kap03 sample EDIs missing")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(KAP03))


@pytest.fixture
def loaded(win, kap03):
    win.set_sites(kap03)
    return win


def _section(label):
    return next(i for i, s in enumerate(SECTIONS) if s.label == label)


def _shown(w) -> bool:
    return w._canvas_view._stack.currentWidget() is w._canvas


def _card(w) -> str:
    return " ".join(lbl.text() for lbl in
                    w._canvas_view._unavailable.findChildren(QLabel))


# ── studio catalogue (Qt-free) ───────────────────────────────────────────────


class TestStudioCatalogue:
    def test_sections_cover_groups_and_utilities(self):
        keys = [s.key for s in SECTIONS]
        assert keys.count("topo") == 1 and keys.count("conv") == 1
        assert all(s.plots for s in SECTIONS if s.key == "plots")

    def test_options_follow_the_signature(self):
        keys = [f.key for f in plot_options("plot_phase_tensor_psection")]
        assert "cmap" in keys and "pmin" in keys
        assert "period" not in keys
        assert "bins" in [f.key for f in plot_options("plot_strike_rose")]

    def test_kwargs_keep_defaults_when_auto(self):
        kw = plot_kwargs("plot_phase_tensor_psection",
                         {"cmap": "", "pmin": 0.0, "pmax": 0.0, "scale": 0})
        assert kw == {}
        kw = plot_kwargs("plot_phase_tensor_psection",
                         {"cmap": "magma", "pmin": 10.0, "pmax": 0.1})
        assert kw == {"cmap": "magma", "period_range": (0.1, 10.0)}

    def test_inputs(self):
        assert plot_inputs("plot_phase_tensor_strip") == ("station",)
        assert "lines" in plot_inputs("plot_phase_tensor_strip_grid")
        assert plot_inputs("plot_atom_psection") == ("model",)

    def test_line_groups(self):
        names = ["L1-01", "L1-02", "L2-01", "gv100", "gv101"]
        g = line_groups(names)
        assert g == {"L1": ["L1-01", "L1-02"], "L2": ["L2-01"],
                     "gv": ["gv100", "gv101"]}
        assert line_groups(names, "single") == {"all stations": names}
        assert line_groups([]) == {}


# ── construction & navigation ────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "Advanced Tools" in win.windowTitle()

    def test_nav_lists_every_section(self, win):
        labels = [win._nav.item(r).text() for r in range(win._nav.count())
                  if win._nav.item(r).data(0x0100) is not None]
        assert labels == [s.label for s in SECTIONS]

    def test_default_page_is_plot_page(self, win):
        assert win._page_stack.currentIndex() == 0
        assert win._plot_list.count() == len(SECTIONS[0].plots)

    def test_empty_state_card(self, win):
        assert not _shown(win)
        assert "EDI/XML" in _card(win)

    def test_model_group_hidden_by_default(self, win):
        assert not win._grp_model.isVisibleTo(win)

    def test_topo_file_row_hidden_by_default(self, win):
        win.select_section(TOPO)
        assert not win._topo_file_row.isVisibleTo(win)

    def test_figures_are_publication_white(self, win):
        win.set_dark_mode(True)
        assert win._ctrl.dark is False
        assert win._topo_ctrl.dark is False
        assert win._conv_ctrl.dark is False


class TestNavigation:
    def test_utility_pages(self, win):
        win.select_section(TOPO)
        assert win._page_stack.currentIndex() == 1
        win.select_section(CONV)
        assert win._page_stack.currentIndex() == 2

    def test_nav_click_selects_section(self, win):
        i = _section("Impedance / Z")
        row = next(r for r in range(win._nav.count())
                   if win._nav.item(r).data(0x0100) == i)
        win._nav.setCurrentRow(row)
        assert win._section == i
        assert win._plot_list.count() == len(SECTIONS[i].plots)

    def test_out_of_range_is_noop(self, win):
        win.select_section(99)
        assert win._section == 0

    def test_atom_psection_shows_model_group(self, win):
        win.select_plot("plot_atom_psection")
        assert win._grp_model.isVisibleTo(win)

    def test_station_input_only_for_strip(self, win):
        win.select_plot("plot_phase_tensor_strip")
        assert win._station_combo.isVisibleTo(win)
        assert not win._lines_combo.isVisibleTo(win)
        win.select_plot("plot_phase_tensor_strip_grid")
        assert win._lines_combo.isVisibleTo(win)
        assert not win._station_combo.isVisibleTo(win)

    def test_description_follows_plot(self, win):
        win.select_plot("plot_phasor_wheel")
        assert "Phasor wheel" in win._desc_lbl.text()
        assert "polar" in win._desc_lbl.text()

    def test_section_remembers_its_plot(self, win):
        win.select_plot("plot_strike_ribbon")
        win.select_section(TOPO)
        win.select_section(0)
        assert win.current_plot[1] == "plot_strike_ribbon"

    def test_options_form_only_when_plot_has_options(self, win):
        win.select_plot("plot_phase_tensor_psection")
        assert win._form_stack.currentWidget() is \
            win._forms["plot_phase_tensor_psection"]


# ── drawing (real data) ──────────────────────────────────────────────────────


class TestDrawing:
    def test_run_no_sites_explains(self, win):
        win._on_run()
        assert win._status_lbl.text() == "Load survey data first."
        assert not _shown(win)

    def test_set_sites_draws_current_plot(self, loaded):
        assert _shown(loaded)
        assert "stations" in loaded._data_lbl.text()
        assert loaded._station_combo.count() > 1

    @pytest.mark.parametrize("fn", ["plot_phase_tensor_strip",
                                    "plot_phase_tensor_strip_grid",
                                    "plot_strike_rose_by_line"])
    def test_plots_needing_survey_inputs_now_draw(self, loaded, fn):
        # these used to be "Cannot render from catalogue"
        loaded.select_plot(fn)
        assert _shown(loaded), _card(loaded)

    def test_station_choice_reaches_the_plot(self, loaded):
        loaded.select_plot("plot_phase_tensor_strip")
        name = loaded._station_combo.itemText(2)
        loaded._station_combo.setCurrentIndex(2)
        loaded._on_run()
        texts = " ".join(t.get_text() for ax in loaded._figure.axes
                         for t in ax.texts)
        assert name in texts

    def test_option_change_redraws(self, loaded):
        loaded.select_plot("plot_phase_tensor_psection")
        before = loaded._figure
        loaded._forms["plot_phase_tensor_psection"].set_values(
            {"cmap": "magma"})
        loaded._on_run()
        assert loaded._figure is not before and _shown(loaded)

    def test_failure_shows_card(self, loaded, monkeypatch):
        def boom(*a, **k):
            raise RuntimeError("kaboom")

        monkeypatch.setattr(loaded._ctrl, "draw", boom)
        loaded._on_run()
        assert not _shown(loaded)
        assert "kaboom" in _card(loaded)

    def test_superseded_figures_closed(self, loaded):
        import matplotlib.pyplot as plt

        for _ in range(4):
            loaded._on_run()
        assert len(plt.get_fignums()) <= 2

    def test_on_export_opens_dialog(self, loaded, monkeypatch):
        opened = []

        class _Dlg:
            def __init__(self, figure=None, parent=None):
                opened.append(figure)

            def exec(self):
                return 0

        monkeypatch.setattr(
            "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog", _Dlg)
        loaded._on_export()
        assert opened == [loaded._figure]


class TestLibrary:
    def test_pin_show_unpin_export(self, loaded, tmp_path):
        first = loaded._figure
        loaded._pin_current()
        loaded.select_plot("plot_strike_ribbon")
        loaded._pin_current()
        assert len(loaded.pinned_labels()) == 2
        loaded._show_pinned(loaded._gallery.item(0))
        assert loaded._figure is first
        paths = loaded.export_pinned(tmp_path, dpi=50)
        assert len(paths) == 2 and all(p.stat().st_size for p in paths)
        loaded._gallery.setCurrentRow(1)
        loaded._unpin()
        assert len(loaded.pinned_labels()) == 1


class TestAutoRender:
    def test_noop_when_already_rendered(self, win):
        win._auto_rendered = True
        win._ctrl._sites = object()
        win._auto_render_if_ready()

    def test_noop_without_sites(self, win):
        win._auto_render_if_ready()
        assert win._auto_rendered is False

    def test_noop_for_topo_section(self, win):
        win.select_section(TOPO)
        win._ctrl._sites = object()
        win._auto_render_if_ready()
        assert win._auto_rendered is False

    def test_triggers_timer(self, win, monkeypatch):
        win._ctrl._sites = object()
        calls = []
        monkeypatch.setattr(
            "pycsamt.app.desktop.windows.advanced_window.QTimer.singleShot",
            staticmethod(lambda ms, fn: calls.append(fn)))
        win._auto_render_if_ready()
        assert win._auto_rendered is True and len(calls) == 1


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
        win.select_section(TOPO)  # realize the page
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
        assert win._canvas_topo_view.showing_canvas is False
        assert "unavailable" in win._canvas_topo_view._unavailable._title.text().lower()

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
        win.select_section(CONV)  # realize the page
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

    def test_set_sites_refreshes_topo_when_on_topo_page(self, win,
                                                        monkeypatch):
        win.select_section(TOPO)
        calls = []
        monkeypatch.setattr(win, "_refresh_topo_preview",
                            lambda: calls.append(1))
        win.set_sites(object())
        assert calls == [1]

    def test_session_round_trip(self, win, qapp):
        win.select_section(CONV)
        win._chk_auto.setChecked(False)
        store = {}
        win.save_geometry_to(store)
        w2 = AdvancedToolsWindow()
        w2.restore_geometry_from(store)
        assert w2._section == CONV and not w2._chk_auto.isChecked()
        w2.close()

    def test_close_hides_and_signals(self, win):
        got = []
        win.panel_closed.connect(lambda: got.append(1))
        win.close()
        assert got and not win.isVisible()
