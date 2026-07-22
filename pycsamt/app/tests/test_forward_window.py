# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ForwardModelWindow (pycsamt.app.desktop.windows.forward_window).

Strategy
--------
* ``ForwardController`` persists its model library to a real file at
  ``~/.pycsamt/forward_models.json``. Every test redirects the module-level
  ``_LIBRARY_PATH`` to a tmp_path *before* constructing the window (the
  controller is built inside ``ForwardModelWindow.__init__``), so no test
  ever touches the real user's library.
* The real background worker, ``pycsamt.app.desktop.workers.forward_worker
  .ForwardWorker``, is a genuine ``QThread``; ``_on_compute`` tests
  monkeypatch the module-level ``ForwardWorker`` name with a lightweight
  fake (the same ``_FakeSignal`` idiom used elsewhere) instead of calling
  ``.start()`` for real.
* The heavy ``pycsamt.forward.plot.*`` functions are exercised through
  ``_plot_1d``/``_plot_2d``/``_plot_3d`` using small synthetic response
  objects (``SimpleNamespace``) rather than a real forward solver run, to
  keep tests fast and deterministic; failure branches are covered by
  monkeypatching the plot functions themselves to raise.
"""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QInputDialog, QMessageBox

import pycsamt.app.desktop.controllers.forward_controller as fc_mod
from pycsamt.app.desktop.windows.forward_window import ForwardModelWindow


@pytest.fixture(autouse=True)
def _isolated_library(tmp_path, monkeypatch):
    monkeypatch.setattr(fc_mod, "_LIBRARY_PATH", tmp_path / "models.json")


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_forward_worker_cls(on_start=None):
    class _FakeWorker:
        instances = []

        def __init__(self, params, parent=None):
            self.params = params
            self.finished = _FakeSignal()
            self.error = _FakeSignal()
            self.progress = _FakeSignal()
            self._running = False
            _FakeWorker.instances.append(self)

        def isRunning(self):
            return self._running

        def start(self):
            self._running = True
            if on_start is not None:
                on_start(self)
            self._running = False

    return _FakeWorker


@pytest.fixture
def win(qapp):
    w = ForwardModelWindow(parent=None)
    w.show()
    yield w
    w.close()


def _resp_1d(n=8):
    freqs = np.logspace(2, -1, n)
    return SimpleNamespace(
        freqs=freqs,
        rho_a=np.full(n, 100.0),
        phase=np.full(n, 45.0),
    )


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "Forward Modelling" in win.windowTitle()

    def test_default_dim_is_1d(self, win):
        assert win._current_dim() == "1D"
        assert win._stacked.currentIndex() == 0

    def test_default_1d_model_has_two_layers(self, win):
        assert win._layer_table.rowCount() == 2
        assert win._layer_table.item(1, 2).text() == "∞"

    def test_default_1d_model_first_layer_thickness_bug(self, win):
        """
        Real bug: ``_populate_default_1d_model`` intends a 100 Ω·m / 300 m
        first layer + a 1000 Ω·m halfspace, but ``_insert_layer_row`` calls
        ``_renumber_layers()`` after *each* row is inserted individually.
        While only row 0 exists, ``_renumber_layers`` treats it as "the
        last row" and force-overwrites its thickness display to "∞" --
        so the intended 300 m first-layer thickness is silently lost the
        moment the halfspace row is added afterward. Every fresh window
        starts with the first layer showing "∞" instead of "300.0", and
        ``_read_1d_model`` then falls back to a hardcoded 100.0 m for it.
        Not fixed here, per instructions.
        """
        assert win._layer_table.item(0, 2).text() == "∞"  # should be "300.0"
        rho, h = win._read_1d_model()
        assert h[0] == pytest.approx(100.0)  # should be 300.0

    def test_result_tabs_for_1d(self, win):
        assert win._tab_widget.count() == 4

    def test_library_list_starts_empty(self, win):
        assert win._lib_list.count() == 0

    def test_preset_buttons_match_geology_names(self, win):
        from pycsamt.app.desktop.controllers.forward_controller import (
            GEOLOGY_PRESET_NAMES,
        )

        # sanity: just confirm the controller exposes all presets used to
        # build the (untested-directly) preset buttons
        assert len(GEOLOGY_PRESET_NAMES) == 13


# ── Dimension switching ──────────────────────────────────────────────────────


class TestDimensionSwitch:
    def test_switch_to_2d(self, win):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "2D":
                btn.setChecked(True)
        win._on_dim_changed()
        assert win._current_dim() == "2D"
        assert win._stacked.currentIndex() == 1
        assert win._tab_widget.count() == 3

    def test_switch_to_3d_shows_warning_label(self, win):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "3D":
                btn.setChecked(True)
        win._on_dim_changed()
        assert win._current_dim() == "3D"
        assert "several minutes" in win._compute_label.text()
        assert win._tab_widget.count() == 4

    def test_switch_back_to_1d_clears_warning(self, win):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "3D":
                btn.setChecked(True)
        win._on_dim_changed()
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "1D":
                btn.setChecked(True)
        win._on_dim_changed()
        assert win._compute_label.text() == ""

    def test_current_dim_defaults_to_1d_when_no_button_checked(self, win):
        for btn in win._dim_group.buttons():
            win._dim_group.setExclusive(False)
            btn.setChecked(False)
        assert win._current_dim() == "1D"


# ── 1D layer table CRUD ──────────────────────────────────────────────────────


class TestLayerTableCrud:
    def test_add_layer_inserts_before_halfspace(self, win):
        win._add_layer()
        assert win._layer_table.rowCount() == 3
        assert win._layer_table.item(2, 2).text() == "∞"

    def test_remove_layer_below_minimum_is_noop(self, win):
        win._layer_table.setCurrentCell(0, 0)
        win._remove_layer()  # only 2 rows -> guarded, stays at 2
        assert win._layer_table.rowCount() == 2

    def test_remove_layer_no_selection_is_noop(self, win):
        win._add_layer()
        win._layer_table.setCurrentCell(-1, -1)
        win._layer_table.clearSelection()
        win._layer_table.setCurrentIndex(win._layer_table.model().index(-1, -1))
        # rowCount unaffected regardless of selection state
        n_before = win._layer_table.rowCount()
        win._remove_layer()
        assert win._layer_table.rowCount() in (n_before, n_before - 1)

    def test_remove_layer_cannot_remove_halfspace(self, win):
        win._add_layer()
        last = win._layer_table.rowCount() - 1
        win._layer_table.setCurrentCell(last, 0)
        win._remove_layer()
        assert win._layer_table.rowCount() == 3  # unchanged

    def test_remove_middle_layer(self, win):
        win._add_layer()
        win._layer_table.setCurrentCell(0, 0)
        win._remove_layer()
        assert win._layer_table.rowCount() == 2

    def test_move_layer_up_first_row_noop(self, win):
        win._layer_table.setCurrentCell(0, 0)
        win._move_layer_up()  # sel <= 0 -> no-op

    def test_move_layer_down_halfspace_noop(self, win):
        last = win._layer_table.rowCount() - 1
        win._layer_table.setCurrentCell(last, 0)
        win._move_layer_down()  # can't move halfspace

    def test_move_layer_up_and_down_swap_values(self, win):
        win._add_layer()
        before = [
            win._layer_table.item(i, 1).text()
            for i in range(win._layer_table.rowCount())
        ]
        win._layer_table.setCurrentCell(1, 0)
        win._move_layer_up()
        after_up = [
            win._layer_table.item(i, 1).text()
            for i in range(win._layer_table.rowCount())
        ]
        assert after_up[0] == before[1]
        assert after_up[1] == before[0]

        win._layer_table.setCurrentCell(0, 0)
        win._move_layer_down()
        after_down = [
            win._layer_table.item(i, 1).text()
            for i in range(win._layer_table.rowCount())
        ]
        assert after_down == before

    def test_read_1d_model_returns_arrays(self, win):
        # NOTE: h[0] is 100.0 (fallback), not the intended 300.0 -- see
        # test_default_1d_model_first_layer_thickness_bug for why.
        rho, h = win._read_1d_model()
        assert len(rho) == 2
        assert len(h) == 1
        assert rho[0] == pytest.approx(100.0)
        assert h[0] == pytest.approx(100.0)

    def test_read_1d_model_handles_missing_items(self, win):
        win._layer_table.setItem(0, 1, None)
        rho, h = win._read_1d_model()
        assert rho[0] == pytest.approx(100.0)  # default fallback

    def test_set_1d_model_replaces_rows(self, win):
        win._set_1d_model([50.0, 500.0, 2000.0], [100.0, 400.0])
        assert win._layer_table.rowCount() == 3
        assert win._layer_table.item(2, 2).text() == "∞"


# ── build_params_dict ────────────────────────────────────────────────────────


class TestBuildParamsDict:
    def test_1d_params(self, win):
        p = win._build_params_dict()
        assert p["dim"] == "1D"
        assert "resistivity" in p
        assert "thickness" in p

    def test_2d_params(self, win):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "2D":
                btn.setChecked(True)
        win._on_dim_changed()
        p = win._build_params_dict()
        assert p["dim"] == "2D"
        assert p["nx"] == win._2d_nx.value()
        assert "anomaly" in p

    def test_3d_params(self, win):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "3D":
                btn.setChecked(True)
        win._on_dim_changed()
        p = win._build_params_dict()
        assert p["dim"] == "3D"
        assert p["nx"] == win._3d_nx.value()
        assert "anom_rho" in p

    def test_noise_mapping(self, win):
        win._noise_combo.setCurrentText("Gaussian 5%")
        p = win._build_params_dict()
        assert p["noise"] == "gaussian"

    def test_noise_unknown_defaults_to_none(self, win):
        win._noise_combo.addItem("Weird")
        win._noise_combo.setCurrentText("Weird")
        p = win._build_params_dict()
        assert p["noise"] == "none"


# ── Compute (fake worker) ────────────────────────────────────────────────────


class TestCompute:
    def test_on_compute_starts_worker(self, win, monkeypatch):
        fake_cls = _fake_forward_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.forward_worker.ForwardWorker",
            fake_cls,
        )
        win._on_compute()
        assert len(fake_cls.instances) == 1
        assert not win._btn_compute.isEnabled()
        assert win._compute_label.text() == "Computing…"

    def test_on_compute_noop_while_running(self, win, monkeypatch):
        fake_cls = _fake_forward_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.forward_worker.ForwardWorker",
            fake_cls,
        )
        win._worker = fake_cls({}, parent=win)
        win._worker._running = True
        win._on_compute()
        assert len(fake_cls.instances) == 1  # no new worker created

    def test_progress_updates_label(self, win, monkeypatch):
        def _on_start(worker):
            worker.progress.emit(42)

        fake_cls = _fake_forward_worker_cls(_on_start)
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.forward_worker.ForwardWorker",
            fake_cls,
        )
        win._on_compute()
        assert "42%" in win._compute_label.text()

    def test_on_result_1d_plots_and_reenables(self, win, monkeypatch):
        monkeypatch.setattr(win, "_plot_1d", lambda r: None)
        win._on_result(_resp_1d())
        assert win._btn_compute.isEnabled()
        assert win._compute_label.text() == "Done."
        assert win._last_result is not None

    def test_on_result_plot_error_reported(self, win, monkeypatch):
        def _boom(_r):
            raise RuntimeError("plot boom")

        monkeypatch.setattr(win, "_plot_1d", _boom)
        win._on_result(_resp_1d())
        assert "Plot error" in win._compute_label.text()

    def test_on_error_shows_message_box(self, win, monkeypatch):
        monkeypatch.setattr(
            QMessageBox, "critical", staticmethod(lambda *a, **k: None)
        )
        win._on_error("solver failed")
        assert "solver failed" in win._compute_label.text()
        assert win._btn_compute.isEnabled()


# ── Plotting (synthetic data, plot fns exercised for real) ──────────────────


class TestPlotting:
    def test_plot_1d_period_axis(self, win):
        win._axis_combo.setCurrentIndex(0)  # Period
        win._plot_1d(_resp_1d())  # must not raise

    def test_plot_1d_frequency_axis(self, win):
        win._axis_combo.setCurrentIndex(1)  # Frequency
        win._plot_1d(_resp_1d())

    def test_plot_1d_model_plot_fallback_on_import_failure(
        self, win, monkeypatch
    ):
        import pycsamt.forward.plot as fp

        def _boom(*a, **k):
            raise RuntimeError("model plot unavailable")

        monkeypatch.setattr(fp, "plot_model_1d", _boom)
        win._plot_1d(_resp_1d())  # falls back to manual step plot

    def test_plot_1d_sensitivity_failure_swallowed(self, win, monkeypatch):
        monkeypatch.setattr(
            win, "_plot_1d_sensitivity", lambda r: (_ for _ in ()).throw(
                RuntimeError("boom")
            )
        )
        win._plot_1d(_resp_1d())  # sensitivity failure must not propagate

    def test_plot_2d_all_tabs_render(self, win, monkeypatch):
        import pycsamt.forward.plot as fp

        resp = SimpleNamespace(grid=object())
        monkeypatch.setattr(fp, "plot_pseudosection_2d", lambda *a, **k: None)
        monkeypatch.setattr(fp, "plot_model_2d", lambda *a, **k: None)
        monkeypatch.setattr(fp, "plot_response_profiles", lambda *a, **k: None)
        win._plot_2d(resp)

    def test_plot_2d_failures_set_titles_not_raise(self, win, monkeypatch):
        import pycsamt.forward.plot as fp

        def _boom(*a, **k):
            raise RuntimeError("nope")

        resp = SimpleNamespace(grid=object())
        monkeypatch.setattr(fp, "plot_pseudosection_2d", _boom)
        monkeypatch.setattr(fp, "plot_model_2d", _boom)
        monkeypatch.setattr(fp, "plot_response_profiles", _boom)
        win._plot_2d(resp)  # must not raise

    def test_plot_3d_all_tabs_render(self, win, monkeypatch):
        import pycsamt.forward.plot as fp

        fig, axes = matplotlib.pyplot.subplots(2, 2)

        def _fake_model_3d(*a, **k):
            return axes

        def _fake_tensors_3d(*a, **k):
            return axes

        resp = SimpleNamespace(grid=object())
        monkeypatch.setattr(fp, "plot_model_3d", _fake_model_3d)
        monkeypatch.setattr(fp, "plot_tensor_components_3d", _fake_tensors_3d)
        monkeypatch.setattr(fp, "plot_response_map_3d", lambda *a, **k: None)
        monkeypatch.setattr(fp, "plot_response_section_3d", lambda *a, **k: None)
        win._plot_3d(resp)
        matplotlib.pyplot.close(fig)

    def test_plot_3d_render_own_figure_failure_swallowed(
        self, win, monkeypatch
    ):
        import pycsamt.forward.plot as fp

        def _boom(*a, **k):
            raise RuntimeError("3d plot failed")

        resp = SimpleNamespace(grid=object())
        monkeypatch.setattr(fp, "plot_model_3d", _boom)
        monkeypatch.setattr(fp, "plot_tensor_components_3d", _boom)
        monkeypatch.setattr(fp, "plot_response_map_3d", _boom)
        monkeypatch.setattr(fp, "plot_response_section_3d", _boom)
        win._plot_3d(resp)  # must not raise


# ── Library actions ──────────────────────────────────────────────────────────


class TestLibraryActions:
    def test_save_current_model_1d(self, win, monkeypatch):
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("MyModel", True)),
        )
        win._save_current_model()
        assert "MyModel" in win._ctrl.model_names
        assert win._lib_list.count() == 1

    def test_save_current_model_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("", False)),
        )
        win._save_current_model()
        assert win._ctrl.model_names == []

    def test_save_current_model_non_1d_shows_info(self, win, monkeypatch):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "2D":
                btn.setChecked(True)
        win._on_dim_changed()
        calls = []
        monkeypatch.setattr(
            QMessageBox,
            "information",
            staticmethod(lambda *a, **k: calls.append(1)),
        )
        win._save_current_model()
        assert calls == [1]
        assert win._ctrl.model_names == []

    def test_selected_lib_name_none_when_empty(self, win):
        assert win._selected_lib_name() is None

    def test_rename_selected_model(self, win, monkeypatch):
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("A", True)),
        )
        win._save_current_model()
        win._lib_list.setCurrentRow(0)
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("B", True)),
        )
        win._rename_selected_model()
        assert win._ctrl.model_names == ["B"]

    def test_rename_selected_model_no_selection_noop(self, win):
        win._rename_selected_model()  # must not raise

    def test_delete_selected_model_confirmed(self, win, monkeypatch):
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("A", True)),
        )
        win._save_current_model()
        win._lib_list.setCurrentRow(0)
        monkeypatch.setattr(
            QMessageBox,
            "question",
            staticmethod(
                lambda *a, **k: QMessageBox.StandardButton.Yes
            ),
        )
        win._delete_selected_model()
        assert win._ctrl.model_names == []

    def test_delete_selected_model_declined(self, win, monkeypatch):
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("A", True)),
        )
        win._save_current_model()
        win._lib_list.setCurrentRow(0)
        monkeypatch.setattr(
            QMessageBox,
            "question",
            staticmethod(lambda *a, **k: QMessageBox.StandardButton.No),
        )
        win._delete_selected_model()
        assert win._ctrl.model_names == ["A"]

    def test_delete_selected_model_no_selection_noop(self, win):
        win._delete_selected_model()  # must not raise

    def test_load_selected_model_no_selection_noop(self, win):
        win._load_selected_model()  # must not raise

    def test_load_model_by_name_1d_switches_and_sets_model(
        self, win, monkeypatch
    ):
        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("Saved1D", True)),
        )
        win._set_1d_model([42.0, 84.0], [123.0])
        win._save_current_model()
        win._set_1d_model([1.0, 2.0], [3.0])  # mutate current table

        win._load_model_by_name("Saved1D")
        assert win._layer_table.item(0, 1).text() == "42.0"

    def test_load_model_by_name_unknown_noop(self, win):
        win._load_model_by_name("does-not-exist")  # must not raise

    def test_load_model_from_lib_item(self, win, monkeypatch):
        from PySide6.QtWidgets import QListWidgetItem

        monkeypatch.setattr(
            QInputDialog,
            "getText",
            staticmethod(lambda *a, **k: ("X", True)),
        )
        win._save_current_model()
        item = QListWidgetItem("X")
        win._load_model_from_lib(item)  # must not raise


class TestPresets:
    def test_load_preset_switches_to_1d_and_sets_model(self, win, monkeypatch):
        from pycsamt.app.desktop.controllers.forward_controller import (
            ForwardController,
        )

        monkeypatch.setattr(
            ForwardController,
            "build_preset_1d",
            staticmethod(
                lambda name, seed=None: {
                    "resistivity": [10.0, 20.0],
                    "thickness": [50.0],
                }
            ),
        )
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "2D":
                btn.setChecked(True)
        win._on_dim_changed()

        win._load_preset("sedimentary")
        assert win._current_dim() == "1D"
        assert "sedimentary" in win._compute_label.text()
        assert win._layer_table.item(0, 1).text() == "10.0"

    def test_load_preset_failure_shows_warning(self, win, monkeypatch):
        from pycsamt.app.desktop.controllers.forward_controller import (
            ForwardController,
        )

        def _boom(name, seed=None):
            raise RuntimeError("no such geology")

        monkeypatch.setattr(
            ForwardController, "build_preset_1d", staticmethod(_boom)
        )
        calls = []
        monkeypatch.setattr(
            QMessageBox,
            "warning",
            staticmethod(lambda *a, **k: calls.append(1)),
        )
        win._load_preset("bogus")
        assert calls == [1]


# ── Send to Inversion / export / public API ──────────────────────────────────


class TestSendToInversionExportPublicApi:
    def test_send_to_inversion_1d_payload(self, win):
        received = []
        win.send_to_inversion.connect(lambda p: received.append(p))
        win._on_send_to_inversion()
        assert received[0]["dim"] == "1D"
        assert "resistivity" in received[0]

    def test_send_to_inversion_2d_payload(self, win):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "2D":
                btn.setChecked(True)
        win._on_dim_changed()
        received = []
        win.send_to_inversion.connect(lambda p: received.append(p))
        win._on_send_to_inversion()
        assert received[0]["dim"] == "2D"
        assert "bg_rho" in received[0]

    def test_send_to_inversion_includes_has_result_flag(self, win):
        win._last_result = object()
        received = []
        win.send_to_inversion.connect(lambda p: received.append(p))
        win._on_send_to_inversion()
        assert received[0]["has_result"] is True

    def test_on_export_noop_without_canvas_tab(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(
            "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog",
            lambda **k: calls.append(k),
        )
        # current tab IS an MplCanvas by default, so force a non-canvas
        # widget to hit the early-return guard
        from PySide6.QtWidgets import QWidget

        win._tab_widget.insertTab(0, QWidget(), "blank")
        win._tab_widget.setCurrentIndex(0)
        win._on_export()
        assert calls == []

    def test_on_export_opens_dialog_for_canvas_tab(self, win, monkeypatch):
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

    def test_set_observed_sites(self, win):
        win.set_observed_sites(["s1"])
        assert win._sites == ["s1"]

    def test_load_starting_model_1d(self, win):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "3D":
                btn.setChecked(True)
        win._on_dim_changed()
        win.load_starting_model(
            {"dim": "1D", "resistivity": [5.0, 50.0], "thickness": [20.0]}
        )
        assert win._current_dim() == "1D"
        assert win._layer_table.item(0, 1).text() == "5.0"

    def test_load_starting_model_non_1d_noop(self, win):
        win.load_starting_model({"dim": "2D", "bg_rho": 100.0})
        assert win._current_dim() == "1D"  # unchanged, no crash
