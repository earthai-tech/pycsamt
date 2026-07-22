# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for InterpretationWindow (pycsamt.app.desktop.windows.interp_window).

Strategy
--------
* No real EDI/model data is loaded for most tests: ``InterpController``
  degrades gracefully (placeholder figures / "No model loaded." strings)
  when there's no sites/model, which already exercises the vast majority
  of the window's wiring cheaply and deterministically.
* ``_RunWorker`` / ``_GeneratePlotWorker`` (QThread subclasses) are always
  exercised via a direct ``.run()`` call, never ``.start()`` — this runs
  the work synchronously on the test thread so signal results can be
  asserted immediately, matching the pattern used for the analogous
  workers in test_dimensionality_tool.py / test_frequency_editor_tool.py.
* QFileDialog / QInputDialog are monkeypatched to avoid real (blocking)
  native dialogs under the offscreen Qt platform.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QFileDialog, QInputDialog

from pycsamt.app.desktop.windows.interp_window import (
    InterpretationWindow,
    _GeneratePlotWorker,
    _RunWorker,
)


@pytest.fixture
def win(qapp):
    w = InterpretationWindow(parent=None)
    yield w
    w.close()


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "Interpretation" in win.windowTitle()

    def test_category_combo_populated(self, win):
        from pycsamt.app.desktop.controllers.interp_controller import (
            CATEGORIES,
        )

        assert win._combo_category.count() == len(CATEGORIES)

    def test_first_category_selected_by_default(self, win):
        assert win._combo_category.currentIndex() == 0
        assert win._combo_plot.count() > 0

    def test_status_card_shows_no_model_initially(self, win):
        assert "No model loaded" in win._status_card.text()

    def test_placeholder_tab_present(self, win):
        assert win._tab_widget.count() == 1
        assert win._tab_widget.widget(0) is win._placeholder


# ── Category / plot navigation ──────────────────────────────────────────────


class TestCategoryNavigation:
    def test_switch_to_geological_builds_station_combo(self, win):
        from pycsamt.app.desktop.controllers.interp_controller import (
            CATEGORIES,
        )

        idx = CATEGORIES.index("Geological")
        win._combo_category.setCurrentIndex(idx)
        assert "station" in win._param_widgets
        w = win._param_widgets["station"]
        assert w.itemText(0) == "(all)"

    def test_switch_to_hydrology_builds_spinboxes(self, win):
        from pycsamt.app.desktop.controllers.interp_controller import (
            CATEGORIES,
        )

        idx = CATEGORIES.index("Hydrology")
        win._combo_category.setCurrentIndex(idx)
        for key in ("m", "n", "rho_w", "phi"):
            assert key in win._param_widgets

    def test_switch_to_uncertainty_builds_n_samples_spin(self, win):
        from pycsamt.app.desktop.controllers.interp_controller import (
            CATEGORIES,
        )

        idx = CATEGORIES.index("Uncertainty (MC)")
        win._combo_category.setCurrentIndex(idx)
        assert "n_samples" in win._param_widgets
        assert win._param_widgets["n_samples"].value() == 300

    def test_switch_to_export_has_no_inline_params(self, win):
        from pycsamt.app.desktop.controllers.interp_controller import (
            CATEGORIES,
        )

        idx = CATEGORIES.index("Export")
        win._combo_category.setCurrentIndex(idx)
        assert win._param_widgets == {}

    def test_plot_changed_updates_description(self, win):
        from pycsamt.app.desktop.controllers.interp_controller import (
            CATEGORIES,
        )

        idx = CATEGORIES.index("Setup & Model")
        win._combo_category.setCurrentIndex(idx)
        if win._combo_plot.count() > 1:
            win._combo_plot.setCurrentIndex(1)
            assert win._plot_desc.text() != ""

    def test_update_plot_desc_out_of_range_clears_label(self, win):
        win._plot_desc.setText("stale")
        win._update_plot_desc(999, 999)
        assert win._plot_desc.text() == ""

    def test_on_category_changed_out_of_range_is_noop(self, win):
        win._on_category_changed(-1)  # must not raise
        win._on_category_changed(999)  # must not raise


# ── Generate plot (worker run synchronously) ────────────────────────────────


class TestGeneratePlot:
    def test_generate_worker_produces_tab_and_thumbnail(self, win):
        results = []
        worker = _GeneratePlotWorker(
            win._ctrl, "plot_model_summary", "Model Summary", {}
        )
        worker.finished.connect(lambda fig, label: results.append((fig, label)))
        worker.run()

        assert len(results) == 1
        fig, label = results[0]
        win._on_plot_ready(fig, label)
        # placeholder tab replaced by the new plot tab
        assert win._tab_widget.count() == 1
        assert win._tab_widget.widget(0) is not win._placeholder
        assert win._gallery.count() == 1
        assert win._btn_generate.isEnabled()

    def test_generate_worker_unknown_method_yields_placeholder(self, win):
        """
        ``InterpController.generate`` never raises -- an unknown method
        name falls into its own ``_not_implemented`` placeholder Figure,
        so ``_GeneratePlotWorker.failed`` only fires if constructing the
        worker's Figure itself throws, which normally can't happen; this
        exercises the ``finished`` (not ``failed``) path for a bogus name.
        """
        results = []
        worker = _GeneratePlotWorker(
            win._ctrl, "not_a_real_method", "Bogus", {}
        )
        worker.finished.connect(
            lambda fig, label: results.append((fig, label))
        )
        worker.run()
        assert len(results) == 1
        win._on_plot_ready(*results[0])
        assert win._run_status.text() == "✓  Bogus"

    def test_generate_plot_worker_failed_signal_on_thread_error(self, win):
        """Force the worker's own outer try/except (around ctrl.generate
        itself) by making generate() raise directly."""
        results = []

        class _BoomCtrl:
            def generate(self, *_a, **_k):
                raise RuntimeError("thread-level failure")

        worker = _GeneratePlotWorker(_BoomCtrl(), "anything", "X", {})
        worker.failed.connect(lambda msg: results.append(msg))
        worker.run()
        assert results == ["thread-level failure"]

    def test_on_generate_noop_when_no_category(self, win):
        win._combo_category.setCurrentIndex(-1)
        win._on_generate()  # must not raise; guarded by cat_row < 0

    def test_on_generate_starts_plot_worker(self, win, monkeypatch):
        """Exercises _on_generate itself (kwargs collection + worker wiring)
        without actually starting a background QThread."""
        started = []
        monkeypatch.setattr(
            _GeneratePlotWorker, "start", lambda self: started.append(self)
        )
        win._on_generate()
        assert len(started) == 1
        assert not win._btn_generate.isEnabled()


class TestTabManagement:
    def test_add_plot_tab_replaces_placeholder(self, win):
        fig = win._ctrl.generate("plot_model_summary")
        win._add_plot_tab(fig, "Test Plot")
        assert win._tab_widget.count() == 1
        assert "Test Plot" in win._tab_widget.tabText(0)

    def test_close_tab_ignores_placeholder(self, win):
        win._close_tab(0)  # widget(0) is the placeholder -> no-op
        assert win._tab_widget.count() == 1
        assert win._tab_widget.widget(0) is win._placeholder

    def test_close_tab_restores_placeholder_when_empty(self, win):
        fig = win._ctrl.generate("plot_model_summary")
        win._add_plot_tab(fig, "Test Plot")
        win._close_tab(0)
        assert win._tab_widget.count() == 1
        assert win._tab_widget.widget(0) is win._placeholder

    def test_rename_tab_applies_new_name(self, win, monkeypatch):
        fig = win._ctrl.generate("plot_model_summary")
        win._add_plot_tab(fig, "Old Name")
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("New Name", True)),
        )
        win._rename_tab(0)
        assert "New Name" in win._tab_widget.tabText(0)

    def test_rename_tab_cancelled_keeps_name(self, win, monkeypatch):
        fig = win._ctrl.generate("plot_model_summary")
        win._add_plot_tab(fig, "Keep Me")
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("Ignored", False)),
        )
        win._rename_tab(0)
        assert "Keep Me" in win._tab_widget.tabText(0)

    def test_rename_tab_on_placeholder_is_noop(self, win):
        win._rename_tab(0)  # widget(0) is placeholder -> returns early

    def test_gallery_click_switches_tab(self, win):
        fig = win._ctrl.generate("plot_model_summary")
        win._add_plot_tab(fig, "P1")
        item = win._gallery.item(0)
        win._tab_widget.setCurrentIndex(0)
        win._on_gallery_click(item)  # must not raise


# ── Run actions (worker run synchronously) ──────────────────────────────────


class TestRunActions:
    def test_run_geological_worker_no_model(self, win):
        results = []
        worker = _RunWorker(win._ctrl.run_geological)
        worker.finished.connect(lambda msg: results.append(msg))
        worker.run()
        assert results == ["No model loaded."]
        win._on_run_finished(results[0])
        assert win._run_status.text() == "No model loaded."

    def test_run_worker_catches_exception(self, win):
        def _boom():
            raise RuntimeError("synthetic")

        results = []
        worker = _RunWorker(_boom)
        worker.finished.connect(lambda msg: results.append(msg))
        worker.run()
        assert results == ["Error: synthetic"]

    def test_on_run_geological_wires_worker(self, win, monkeypatch):
        started = []
        monkeypatch.setattr(
            _RunWorker, "start", lambda self: started.append(self)
        )
        win._on_run_geological()
        assert len(started) == 1
        assert not win._btn_run_geo.isEnabled()

    def test_on_run_hydro_wires_worker_and_pushes_params(
        self, win, monkeypatch
    ):
        from pycsamt.app.desktop.controllers.interp_controller import (
            CATEGORIES,
        )

        idx = CATEGORIES.index("Hydrology")
        win._combo_category.setCurrentIndex(idx)
        started = []
        monkeypatch.setattr(
            _RunWorker, "start", lambda self: started.append(self)
        )
        win._on_run_hydro()
        assert len(started) == 1
        assert win._ctrl.state.petro_cfg is not None
        assert not win._btn_run_hydro.isEnabled()

    def test_on_run_mc_uses_default_when_no_param_widget(
        self, win, monkeypatch
    ):
        started = []
        monkeypatch.setattr(
            _RunWorker, "start", lambda self: started.append(self)
        )
        win._param_widgets = {}
        win._on_run_mc()
        assert len(started) == 1
        assert "200" in win._run_status.text()

    def test_on_run_mc_uses_spin_value(self, win, monkeypatch):
        from pycsamt.app.desktop.controllers.interp_controller import (
            CATEGORIES,
        )

        idx = CATEGORIES.index("Uncertainty (MC)")
        win._combo_category.setCurrentIndex(idx)
        win._param_widgets["n_samples"].setValue(500)
        started = []
        monkeypatch.setattr(
            _RunWorker, "start", lambda self: started.append(self)
        )
        win._on_run_mc()
        assert "500" in win._run_status.text()


# ── Model loading ────────────────────────────────────────────────────────────


class TestModelLoading:
    def test_load_from_inversion_no_parent(self, win):
        win._load_from_inversion()
        assert "No model available" in win._run_status.text()

    def test_load_from_inversion_parent_without_attr(self, win, qapp):
        from PySide6.QtWidgets import QWidget

        parent = QWidget()
        win.setParent(parent)
        win._load_from_inversion()
        assert "No model available" in win._run_status.text()
        win.setParent(None)
        parent.deleteLater()

    def test_load_from_inversion_pulls_model_from_parent(self, win, qapp):
        from PySide6.QtWidgets import QWidget

        class _FakeModel:
            n_x = 4
            n_z = 3
            depth_max = 100.0
            profile_length = 500.0
            method = "Occam2D"
            rms = 1.0

        class _FakeInvWin:
            _result_model = _FakeModel()

        parent = QWidget()
        parent._inversion_win = _FakeInvWin()
        win.setParent(parent)
        win._load_from_inversion()
        assert win._ctrl.has_model
        assert "loaded from Inversion" in win._run_status.text()
        win.setParent(None)
        parent.deleteLater()

    def test_load_from_inversion_model_missing_depth_max_crashes(
        self, win, qapp
    ):
        """
        Real bug: ``model_status`` deliberately uses
        ``getattr(m, "depth_max", None)`` so a model lacking that
        attribute degrades to None instead of raising, but
        ``_update_status_card`` then does ``f"{depth:.0f}"`` on it
        unconditionally -- unlike the adjacent ``rms`` field, which does
        guard with ``rms or '—'``. Any real model missing depth_max (or
        having it explicitly None) crashes the status-card update
        instead of degrading. Not fixed here, per instructions.
        """
        from PySide6.QtWidgets import QWidget

        class _FakeInvWin:
            _result_model = object()  # no depth_max/profile_length/etc.

        parent = QWidget()
        parent._inversion_win = _FakeInvWin()
        win.setParent(parent)
        with pytest.raises(TypeError):
            win._load_from_inversion()
        win.setParent(None)
        parent.deleteLater()

    def test_load_occam2d_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: ""),
        )
        win._load_occam2d()  # must not raise

    def test_load_occam2d_failure_reported(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/not/a/real/dir"),
        )
        win._load_occam2d()
        assert win._run_status.text().startswith("Load failed:")

    def test_load_occam2d_success(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/fake/dir"),
        )
        monkeypatch.setattr(
            win._ctrl, "set_model_from_occam2d", lambda path: None
        )
        win._load_occam2d()
        assert win._run_status.text() == "Occam2D model loaded."


# ── Borehole management ──────────────────────────────────────────────────────


class TestBoreholeManagement:
    def test_add_borehole_csv_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        win._add_borehole_csv()
        assert win._bh_list.count() == 0

    def test_add_borehole_csv_success(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("/fake/bh.csv", "")),
        )
        monkeypatch.setattr(
            win._ctrl, "add_borehole_csv", lambda path: "BH-1"
        )
        win._add_borehole_csv()
        assert win._bh_list.count() == 1
        assert win._bh_list.item(0).text() == "BH-1"

    def test_add_borehole_csv_failure_reported(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("/fake/bh.csv", "")),
        )

        def _boom(path):
            raise ValueError("bad csv")

        monkeypatch.setattr(win._ctrl, "add_borehole_csv", _boom)
        win._add_borehole_csv()
        assert "Borehole load failed" in win._run_status.text()

    def test_add_borehole_las_success(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("/fake/bh.las", "")),
        )
        monkeypatch.setattr(
            win._ctrl, "add_borehole_las", lambda path: "BH-LAS"
        )
        win._add_borehole_las()
        assert win._bh_list.item(0).text() == "BH-LAS"

    def test_add_borehole_las_failure_reported(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("/fake/bh.las", "")),
        )

        def _boom(path):
            raise ValueError("bad las")

        monkeypatch.setattr(win._ctrl, "add_borehole_las", _boom)
        win._add_borehole_las()
        assert "LAS load failed" in win._run_status.text()

    def test_remove_borehole(self, win, monkeypatch):
        monkeypatch.setattr(
            win._ctrl, "remove_borehole", lambda name: None
        )
        win._bh_list.addItem("BH-1")
        win._bh_list.setCurrentRow(0)
        win._remove_borehole()
        assert win._bh_list.count() == 0

    def test_remove_borehole_no_selection_is_noop(self, win):
        win._remove_borehole()  # must not raise


# ── Export ───────────────────────────────────────────────────────────────────


class TestExport:
    @pytest.mark.parametrize(
        "ui_method,ctrl_method",
        [
            ("_export_xyz", "export_xyz"),
            ("_export_las", "export_las"),
            ("_export_csv", "export_csv"),
            ("_export_vtk", "export_vtk"),
        ],
    )
    def test_export_cancelled_leaves_status_untouched(
        self, win, monkeypatch, ui_method, ctrl_method
    ):
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        win._run_status.setText("unchanged")
        getattr(win, ui_method)()
        assert win._run_status.text() == "unchanged"

    def test_export_xyz_success(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("/fake/out.xyz", "")),
        )
        monkeypatch.setattr(
            win._ctrl, "export_xyz", lambda path: "Exported XYZ."
        )
        win._export_xyz()
        assert win._run_status.text() == "Exported XYZ."

    def test_export_las_uses_selected_station(self, win, monkeypatch):
        from pycsamt.app.desktop.controllers.interp_controller import (
            CATEGORIES,
        )

        idx = CATEGORIES.index("Geological")
        win._combo_category.setCurrentIndex(idx)
        win._param_widgets["station"].addItem("STA01")
        win._param_widgets["station"].setCurrentText("STA01")

        captured = {}

        def _fake_export_las(path, station=""):
            captured["station"] = station
            return "Exported LAS."

        monkeypatch.setattr(win._ctrl, "export_las", _fake_export_las)
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("/fake/out.las", "")),
        )
        win._export_las()
        assert captured["station"] == "STA01"

    def test_export_current_figure_no_plot(self, win):
        win._export_current_figure()
        assert win._run_status.text() == "No plot to export."

    def test_export_current_figure_opens_dialog(self, win, monkeypatch):
        fig = win._ctrl.generate("plot_model_summary")
        win._add_plot_tab(fig, "P1")

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
        win._export_current_figure()
        assert calls[-1] == "exec"


# ── set_sites / set_dark_mode / status card / station names ─────────────────


class TestSetSitesAndStatus:
    def test_set_sites_updates_status_card(self, win):
        win.set_sites(object())
        # no model yet -> status card falls into "no model" branch, but the
        # sites-derived station count line is attempted (best-effort).
        assert "No resistivity model" in win._status_card.text()

    def test_set_sites_iteration_failure_is_swallowed(self, win):
        class _Bad:
            def __iter__(self):
                raise RuntimeError("boom")

        win.set_sites(_Bad())  # must not raise

    def test_set_dark_mode_delegates_to_controller(self, win):
        win.set_dark_mode(False)
        assert win._ctrl.dark is False

    def test_status_card_with_model(self, win):
        class _FakeModel:
            n_x = 10
            n_z = 5
            depth_max = 500.0
            profile_length = 2000.0
            method = "Occam2D"
            rms = 1.23

        win._ctrl.set_model(_FakeModel())
        win._update_status_card()
        txt = win._status_card.text()
        assert "Occam2D" in txt
        assert "10 × 5" in txt

    def test_station_names_from_model(self, win):
        class _FakeModel:
            station_names = ["A", "B", "C"]

        win._ctrl.set_model(_FakeModel())
        assert win._station_names() == ["A", "B", "C"]

    def test_station_names_empty_without_sites_or_model(self, win):
        assert win._station_names() == []

    def test_station_names_from_sites_fallback(self, win, monkeypatch):
        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "_iter_items", lambda s: ["ed1", "ed2"])
        monkeypatch.setattr(core_mod, "_name", lambda ed, i: f"ST{i}")
        win._ctrl.set_sites(object())
        assert win._station_names() == ["ST0", "ST1"]

    def test_station_names_swallows_exception(self, win, monkeypatch):
        import pycsamt.emtools._core as core_mod

        def _boom(_sites):
            raise RuntimeError("boom")

        monkeypatch.setattr(core_mod, "_iter_items", _boom)
        win._ctrl.set_sites(object())
        assert win._station_names() == []
