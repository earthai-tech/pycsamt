# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for InversionWindow (pycsamt.app.desktop.windows.inversion_window).

Strategy
--------
* The three real background workers (``InversionWorker``,
  ``AIInversionWorker``) are genuine ``QThread`` subclasses; every
  run-triggering test monkeypatches the module they're imported from
  (inversion_worker.py / ai_inversion_worker.py) with a lightweight fake
  (the same ``_FakeSignal`` idiom used elsewhere) instead of calling
  ``.start()`` for real.
* ``InputBuilder``/``*Config`` classes for occam2d/modem/mare2dem are
  likewise monkeypatched at their source modules so no real solver input
  files are written to disk.
* Plotting helpers (``_plot_trad_result`` / ``_plot_ai_result``) are
  exercised with small synthetic result objects (``SimpleNamespace`` /
  plain dicts) rather than a real inversion run.
"""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QFileDialog, QMessageBox

from pycsamt.app.desktop.windows.inversion_window import InversionWindow


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


class _FakeWorkerBase:
    instances: list = []

    def __init__(self, *a, **kw):
        self.finished = _FakeSignal()
        self.error = _FakeSignal()
        self.progress = _FakeSignal()
        self._running = False
        type(self).instances.append(self)

    def isRunning(self):
        return self._running

    def quit(self):
        pass

    def terminate(self):
        pass

    def start(self):
        self._running = True


def _fake_inversion_worker_cls():
    class _FakeInversionWorker(_FakeWorkerBase):
        instances = []

        def __init__(self, *a, **kw):
            super().__init__(*a, **kw)
            self.stdout_line = _FakeSignal()
            self.kwargs = kw

        def cancel(self):
            pass

    return _FakeInversionWorker


def _fake_ai_worker_cls():
    class _FakeAIWorker(_FakeWorkerBase):
        instances = []

        def __init__(self, params, parent=None):
            super().__init__()
            self.params = params
            self.log_line = _FakeSignal()

    return _FakeAIWorker


@pytest.fixture
def win(qapp):
    w = InversionWindow(parent=None)
    w.show()
    yield w
    w.close()


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "Inversion" in win.windowTitle()

    def test_default_dim_is_2d(self, win):
        assert win._current_dim() == "2D"

    def test_default_engine_is_modem(self, win):
        assert win._current_engine() == "modem"

    def test_trad_stack_visible_by_default(self, win):
        assert win._trad_stack.isVisible()
        assert not win._ai_stack.isVisible()

    def test_stop_button_disabled_initially(self, win):
        assert not win._btn_stop.isEnabled()

    def test_starting_model_label_default(self, win):
        assert "no model" in win._fwd_model_label.text()

    def test_result_tabs_present(self, win):
        assert win._tabs.count() == 4


# ── Dimension / engine switching ─────────────────────────────────────────────


class TestDimEngineSwitching:
    def test_1d_only_allows_inv1d(self, win):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "1D":
                btn.setChecked(True)
        win._on_dim_changed()
        assert not win._rb_occam2d.isEnabled()
        assert not win._rb_modem.isEnabled()
        assert not win._rb_mare2dem.isEnabled()
        assert win._rb_inv1d.isEnabled()
        assert win._current_engine() == "inv1d"
        assert win._ai_stack.isVisible()

    def test_3d_only_allows_modem_and_inv3d(self, win):
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "3D":
                btn.setChecked(True)
        win._on_dim_changed()
        assert win._rb_modem.isEnabled()
        assert not win._rb_occam2d.isEnabled()
        assert not win._rb_mare2dem.isEnabled()
        assert win._rb_inv3d.isEnabled()

    def test_2d_allows_all_traditional(self, win):
        assert win._rb_occam2d.isEnabled()
        assert win._rb_modem.isEnabled()
        assert win._rb_mare2dem.isEnabled()
        assert win._rb_inv2d.isEnabled()

    def test_engine_switch_updates_stack_visibility(self, win):
        win._rb_occam2d.setChecked(True)
        win._on_engine_changed()
        assert win._trad_stack.currentIndex() == 0
        win._rb_inv2d.setChecked(True)
        win._on_engine_changed()
        assert win._ai_stack.isVisible()
        assert not win._trad_stack.isVisible()

    def test_on_engine_changed_noop_when_none_checked(self, win):
        win._engine_group.setExclusive(False)
        for btn in win._engine_group.buttons():
            btn.setChecked(False)
        win._on_engine_changed()  # must not raise
        win._engine_group.setExclusive(True)

    def test_current_engine_defaults_to_modem_when_none_checked(self, win):
        win._engine_group.setExclusive(False)
        for btn in win._engine_group.buttons():
            btn.setChecked(False)
        assert win._current_engine() == "modem"
        win._engine_group.setExclusive(True)

    def test_current_dim_defaults_to_2d_when_none_checked(self, win):
        win._dim_group.setExclusive(False)
        for btn in win._dim_group.buttons():
            btn.setChecked(False)
        assert win._current_dim() == "2D"
        win._dim_group.setExclusive(True)


# ── Starting model ────────────────────────────────────────────────────────────


class TestStartingModel:
    def test_load_starting_model_1d(self, win):
        win.load_starting_model(
            {"dim": "1D", "resistivity": [55.0, 550.0], "thickness": [10.0]}
        )
        assert win._current_dim() == "1D"
        assert win._init_rho.value() == pytest.approx(55.0)
        assert "1D" in win._fwd_model_label.text()
        assert "2 layers" in win._fwd_model_label.text()

    def test_load_starting_model_no_resistivity(self, win):
        win.load_starting_model({"dim": "2D"})
        assert win._init_rho.value() != 0  # untouched, still whatever default

    def test_clear_starting_model(self, win):
        win.load_starting_model({"dim": "1D", "resistivity": [1.0]})
        win._clear_starting_model()
        assert win._starting_model is None
        assert "no model" in win._fwd_model_label.text()


# ── set_sites / selected_sites ────────────────────────────────────────────────


class TestSetSitesAndSelection:
    def test_set_sites_none_clears_list(self, win):
        win.set_sites(None)
        assert win._station_list.count() == 0

    def test_set_sites_populates_list(self, win):
        sites = SimpleNamespace(
            as_list=lambda: [
                SimpleNamespace(name="A"),
                SimpleNamespace(name="B"),
            ]
        )
        win.set_sites(sites)
        assert win._station_list.count() == 2

    def test_set_sites_iteration_failure_swallowed(self, win):
        class _Bad:
            def as_list(self):
                raise RuntimeError("boom")

        win.set_sites(_Bad())  # must not raise

    def test_selected_sites_none_when_no_sites(self, win):
        assert win._selected_sites() is None

    def test_selected_sites_returns_all_when_none_selected(self, win):
        sites = SimpleNamespace(as_list=lambda: [SimpleNamespace(name="A")])
        win.set_sites(sites)
        assert win._selected_sites() is sites

    def test_selected_sites_filters_by_selection(self, win):
        selected = object()
        sites = SimpleNamespace(
            as_list=lambda: [
                SimpleNamespace(name="A"),
                SimpleNamespace(name="B"),
            ],
            select=lambda names: selected,
        )
        win.set_sites(sites)
        win._station_list.item(0).setSelected(True)
        assert win._selected_sites() is selected

    def test_selected_sites_select_raises_falls_back(self, win):
        sites = SimpleNamespace(
            as_list=lambda: [SimpleNamespace(name="A")],
            select=lambda names: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        win.set_sites(sites)
        win._station_list.item(0).setSelected(True)
        assert win._selected_sites() is sites


# ── Browse dialogs ────────────────────────────────────────────────────────────


class TestBrowseDialogs:
    def test_browse_workdir_sets_text(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/chosen/workdir"),
        )
        win._browse_workdir()
        assert win._workdir_edit.text() == "/chosen/workdir"

    def test_browse_workdir_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: ""),
        )
        win._workdir_edit.setText("unchanged")
        win._browse_workdir()
        assert win._workdir_edit.text() == "unchanged"

    def test_browse_binary_sets_text(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getOpenFileName",
            staticmethod(lambda *a, **k: ("/chosen/bin", "")),
        )
        win._browse_binary(win._occ_binary)
        assert win._occ_binary.text() == "/chosen/bin"


# ── Run / Stop ────────────────────────────────────────────────────────────────


class TestRunStop:
    def test_run_no_workdir_warns(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(
            QMessageBox,
            "warning",
            staticmethod(lambda *a, **k: calls.append(1)),
        )
        win._rb_modem.setChecked(True)
        win._on_engine_changed()
        win._on_run()
        assert calls == [1]

    def test_run_occam2d_builds_and_starts_worker(
        self, win, monkeypatch, tmp_path
    ):
        import pycsamt.models.occam2d as occam2d_mod

        class _FakeInputBuilder:
            def __init__(self, source, workdir, config):
                pass

            def build(self):
                return None

        monkeypatch.setattr(occam2d_mod, "InputBuilder", _FakeInputBuilder)
        fake_cls = _fake_inversion_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.inversion_worker.InversionWorker",
            fake_cls,
        )
        win._workdir_edit.setText(str(tmp_path))
        win._rb_occam2d.setChecked(True)
        win._on_engine_changed()
        win._on_run()
        assert len(fake_cls.instances) == 1
        assert not win._btn_run.isEnabled()
        assert win._btn_stop.isEnabled()

    def test_run_occam2d_input_builder_failure_logged(
        self, win, monkeypatch, tmp_path
    ):
        import pycsamt.models.occam2d as occam2d_mod

        class _BoomBuilder:
            def __init__(self, source, workdir, config):
                pass

            def build(self):
                raise RuntimeError("bad input")

        monkeypatch.setattr(occam2d_mod, "InputBuilder", _BoomBuilder)
        win._workdir_edit.setText(str(tmp_path))
        win._rb_occam2d.setChecked(True)
        win._on_engine_changed()
        win._on_run()
        assert "InputBuilder error" in win._log_edit.toPlainText()
        assert win._worker is None

    def test_run_modem_2d_builds_and_starts_worker(
        self, win, monkeypatch, tmp_path
    ):
        import pycsamt.models.modem as modem_mod

        class _FakeInputBuilder:
            def __init__(self, source, workdir, config):
                pass

            def build(self):
                return None

        monkeypatch.setattr(modem_mod, "InputBuilder", _FakeInputBuilder)
        fake_cls = _fake_inversion_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.inversion_worker.InversionWorker",
            fake_cls,
        )
        win._workdir_edit.setText(str(tmp_path))
        win._mod_mode.setCurrentText("2D")
        win._on_run()
        assert len(fake_cls.instances) == 1

    def test_run_mare2dem_builds_and_starts_worker(
        self, win, monkeypatch, tmp_path
    ):
        import pycsamt.models.mare2dem as m2d_mod

        class _FakeInputBuilder:
            def __init__(self, source, workdir, config):
                pass

            def build(self):
                return None

        monkeypatch.setattr(m2d_mod, "InputBuilder", _FakeInputBuilder)
        fake_cls = _fake_inversion_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.inversion_worker.InversionWorker",
            fake_cls,
        )
        win._workdir_edit.setText(str(tmp_path))
        win._rb_mare2dem.setChecked(True)
        win._on_engine_changed()
        win._on_run()
        assert len(fake_cls.instances) == 1

    def test_run_noop_while_worker_running(self, win, monkeypatch, tmp_path):
        fake_cls = _fake_inversion_worker_cls()
        w = fake_cls()
        w._running = True
        win._worker = w
        win._workdir_edit.setText(str(tmp_path))
        win._on_run()  # early return, no new worker attached
        assert win._worker is w

    def test_run_ai_starts_worker(self, win, monkeypatch):
        fake_cls = _fake_ai_worker_cls()
        monkeypatch.setattr(
            "pycsamt.app.desktop.workers.ai_inversion_worker.AIInversionWorker",
            fake_cls,
        )
        for btn in win._dim_group.buttons():
            if btn.property("dim") == "1D":
                btn.setChecked(True)
        win._on_dim_changed()
        win._on_run()
        assert len(fake_cls.instances) == 1
        assert win._tabs.currentIndex() == 0

    def test_on_stop_cancels_worker_with_cancel_method(self, win):
        fake_cls = _fake_inversion_worker_cls()
        w = fake_cls()
        w._running = True
        win._worker = w
        win._on_stop()
        assert win._btn_run.isEnabled()
        assert not win._btn_stop.isEnabled()
        assert "Stopped." in win._log_edit.toPlainText()

    def test_on_stop_falls_back_to_terminate(self, win):
        class _NoCancel(_FakeWorkerBase):
            pass

        w = _NoCancel()
        win._worker = w
        win._on_stop()  # AttributeError on cancel() -> terminate()
        assert win._run_status.text() == "Stopped."

    def test_on_stop_no_worker_noop(self, win):
        win._on_stop()  # must not raise


# ── build_ai_params ───────────────────────────────────────────────────────────


class TestBuildAiParams:
    def test_inv1d_params(self, win):
        p = win._build_ai_params("inv1d", "1D")
        assert p["arch"] == win._ai1_arch.currentText()
        assert p["geology"] is None  # "(none)" selected by default

    def test_inv1d_params_with_geology(self, win):
        win._ai1_geology.setCurrentText("volcanic")
        p = win._build_ai_params("inv1d", "1D")
        assert p["geology"] == "volcanic"

    def test_inv2d_params(self, win):
        p = win._build_ai_params("inv2d", "2D")
        assert p["n_components"] == win._ai2_n_comp.value()

    def test_inv3d_params_valid_hidden(self, win):
        win._ai3_hidden.setText("64, 32")
        p = win._build_ai_params("inv3d", "3D")
        assert p["hidden"] == [64, 32]

    def test_inv3d_params_invalid_hidden_falls_back(self, win):
        win._ai3_hidden.setText("not,numbers")
        p = win._build_ai_params("inv3d", "3D")
        assert p["hidden"] == [256, 128, 64]

    def test_build_x_obs_none_without_sites(self, win):
        p = {"f_min": 1e-3, "f_max": 1e3, "n_freq": 10}
        assert win._build_X_obs("inv1d", "1D", p) is None

    def test_build_x_obs_with_sites(self, win):
        class _FakeSite:
            def interpolate_rho_a(self, freqs):
                return np.full(len(freqs), 100.0)

            def interpolate_phase(self, freqs):
                return np.full(len(freqs), 45.0)

        sites = SimpleNamespace(as_list=lambda: [_FakeSite(), _FakeSite()])
        win.set_sites(sites)
        p = {"f_min": 1e-3, "f_max": 1e3, "n_freq": 10}
        X = win._build_X_obs("inv1d", "1D", p)
        assert X is not None
        assert X.shape == (2, 20)

    def test_build_x_obs_exception_returns_none(self, win):
        class _BadSites:
            def as_list(self):
                raise RuntimeError("boom")

        win._sites = _BadSites()
        p = {"f_min": 1e-3, "f_max": 1e3, "n_freq": 10}
        assert win._build_X_obs("inv1d", "1D", p) is None


# ── Result handlers / plotting ────────────────────────────────────────────────


class TestResultHandlersAndPlotting:
    def test_on_trad_finished_success(self, win):
        received = []
        win.result_ready.connect(lambda p: received.append(p))

        class _Result:
            def plot_model(self, ax):
                ax.set_title("model")

            def plot_response(self, ax):
                ax.set_title("resp")

            def plot_misfit(self, ax):
                ax.set_title("misfit")

        win._on_trad_finished(_Result())
        assert win._run_status.text() == "Done."
        assert received[0]["engine"] == win._current_engine()

    def test_on_trad_finished_plot_model_fallback(self, win, monkeypatch):
        class _Result:
            def plot_response(self, ax):
                pass

            def plot_misfit(self, ax):
                pass

        import pycsamt.models.occam2d as occam2d_mod

        class _FakePlotModel:
            def __init__(self, result):
                pass

            def plot(self, ax):
                ax.set_title("fallback")

        monkeypatch.setattr(occam2d_mod, "PlotModel", _FakePlotModel)
        win._on_trad_finished(_Result())  # plot_model raises AttributeError

    def test_on_trad_finished_all_plots_unavailable(self, win, monkeypatch):
        import pycsamt.models.occam2d as occam2d_mod

        def _boom(*a, **k):
            raise RuntimeError("nope")

        monkeypatch.setattr(occam2d_mod, "PlotModel", _boom)
        win._on_trad_finished(object())  # no plot_* methods at all

    def test_on_ai_finished_1d(self, win):
        received = []
        win.result_ready.connect(lambda p: received.append(p))
        y_pred = np.tile(
            np.concatenate([np.log10([100.0, 200.0]), [50.0]]), (3, 1)
        )
        result = {
            "dim": "1D",
            "y_pred": y_pred,
            "n_layers": 2,
            "freqs": np.logspace(2, -1, 8),
        }
        win._on_ai_finished(result)
        assert win._run_status.text() == "Done."
        assert received[0]["engine"] == win._current_engine()

    def test_on_ai_finished_2d(self, win):
        y_pred = np.random.rand(20, 40)
        result = {
            "dim": "2D",
            "y_pred": y_pred,
            "n_stations": 20,
            "n_depth": 40,
            "freqs": None,
        }
        win._on_ai_finished(result)

    def test_on_ai_finished_3d_with_coords(self, win):
        y_pred = np.random.rand(5, 3)
        coords = np.random.rand(5, 2)
        result = {
            "dim": "3D",
            "y_pred": y_pred,
            "coords": coords,
            "freqs": None,
        }
        win._on_ai_finished(result)

    def test_on_ai_finished_3d_without_coords(self, win):
        y_pred = np.random.rand(5)
        result = {"dim": "3D", "y_pred": y_pred, "freqs": None}
        win._on_ai_finished(result)

    def test_on_ai_finished_no_y_pred(self, win):
        win._on_ai_finished({"dim": "1D", "freqs": None})

    def test_on_ai_finished_with_loss_history(self, win):
        inv = SimpleNamespace(loss_history_=[1.0, 0.5, 0.1])
        result = {
            "dim": "1D",
            "y_pred": None,
            "freqs": None,
            "inverter": inv,
        }
        win._on_ai_finished(result)

    def test_ai_fit_plot_with_observed_data(self, win):
        n_freq = 8
        freqs = np.logspace(2, -1, n_freq)
        X_obs = np.tile(
            np.concatenate(
                [np.log10(np.full(n_freq, 50.0)), np.full(n_freq, 45.0)]
            ),
            (2, 1),
        )
        result = {"X_obs": X_obs}
        win._plot_ai_fit(win._tab_fit.figure.add_subplot(111), result, freqs)

    def test_ai_fit_plot_missing_data(self, win):
        ax = win._tab_fit.figure.add_subplot(111)
        win._plot_ai_fit(ax, {}, None)
        assert "No observed" in ax.get_title()

    def test_on_error_shows_message_and_resets_ui(self, win, monkeypatch):
        monkeypatch.setattr(
            QMessageBox, "critical", staticmethod(lambda *a, **k: None)
        )
        win._on_error("solver crashed")
        assert win._run_status.text() == "Error."
        assert "ERROR: solver crashed" in win._log_edit.toPlainText()
        assert win._btn_run.isEnabled()
