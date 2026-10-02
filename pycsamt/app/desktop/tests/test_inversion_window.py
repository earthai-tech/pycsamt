# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the Inversion Studio (pycsamt.app.desktop.windows.inversion).

* Workers run synchronously (``EngineRunWorker.start`` -> ``run()``), so a
  whole build → run → results cycle happens inside one test.
* Occam1D is pure Python and runs for real on the bundled CSAMT line.
* External solvers are replaced by a fake engine in the shared ``ENGINES``
  registry (build/run/cancel/error paths); their real result folders
  (Occam2D, ModEM 3-D Broken Hill, MARE2DEM demo) are opened in View
  Results.
* Binary discovery and the run/preset store are redirected (no WSL probe,
  no writes to the home folder).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtCore import Qt

import pycsamt.app.desktop.workers.inversion_worker as wmod
import pycsamt.models.solver_build as sb
from pycsamt.app.desktop.controllers import inversion_engines as ie
from pycsamt.app.desktop.controllers import inversion_runs as store
from pycsamt.app.desktop.windows.inversion_window import InversionWindow

ROOT = Path(__file__).resolve().parents[4]  # repository root
CSAMT = ROOT / Path("data/CSAMT")
OCCAM2D = ROOT / Path("data/occam2D")
MODEM3D = ROOT / Path("data/MT/broken-hill/final-models")
M2D_DEMO = ROOT / Path("data/mare2dem/demo_mt_inversion")


def _sites(n=None):
    if not CSAMT.is_dir():
        pytest.skip("bundled CSAMT EDIs not found")
    from pycsamt.io import read_transfer_function
    from pycsamt.site.base import Sites

    return Sites([read_transfer_function(p)
                  for p in sorted(CSAMT.glob("*.edi"))[:n]])


class FakeEngine(ie.Engine):
    """Stands in for an external solver (occam2d slot)."""

    key = "occam2d"
    label = "Fake2D"
    dim = "2D"
    binary_key = "occam2d"
    description = "fake external solver"

    def __init__(self):
        self.builds, self.binaries = [], []
        self.mode = "ok"
        self.hist = []

    def fields(self):
        return [
            ie.Field("n", "Cells", "int", 10, 1, 100, section="Mesh"),
            ie.Field("max_iterations", "Max", "int", 5, 1, 100),
            ie.Field("target_misfit", "Target", "float", 1.0, 0.1, 10, 0.1),
            ie.Field("use_mpi", "MPI", "bool", False, section="Run"),
        ]

    def check_sites(self, sites):
        return None

    def build(self, sites, workdir, values):
        self.builds.append(dict(values))
        if self.mode == "build_fail":
            raise ValueError("bad mesh")
        return ie.BuildInfo(self.key, Path(workdir), values=dict(values),
                            summary=[("Cells", str(values["n"]))],
                            stations=["A", "B"])

    def plot_mesh(self, fig, build):
        fig.add_subplot(111).plot([0, 1], [0, 1])

    def make_task(self, build, values, binary):
        def task(rep):
            self.binaries.append(binary)
            rep.stage("solving")
            rep.log("iter 1 rms 2.0")
            rep.iteration(1, 2.0)
            if self.mode == "fail":
                raise RuntimeError("solver crashed")
            if self.mode == "cancel":
                raise InterruptedError("stopped")
            return "RESULT"
        return task

    def history(self, workdir):
        return list(self.hist)

    def load(self, path, iteration=None, station=None):
        return ie.LoadedRun(self.key, Path(path), result=SimpleNamespace(),
                            rms=[(1, 2.0)], iterations=[1], iteration=1,
                            summary=[("Final RMS", "2.0")])

    def views(self, run):
        return [("model", "Model")]

    def render(self, key, fig, run):
        fig.add_subplot(111).plot([1, 2])
        return fig


@pytest.fixture
def env(monkeypatch, tmp_path):
    monkeypatch.setattr(wmod.EngineRunWorker, "start",
                        lambda self: self.run())
    found = {"occam2d": [], "modem2d": [],
             "modem3d": [sb.BinaryCandidate("C:/m/Mod3DMT.exe", "PATH")],
             "mare2dem": [sb.BinaryCandidate("wsl:/opt/MARE2DEM",
                                             "WSL build")]}
    monkeypatch.setattr(sb, "discover_binaries",
                        lambda key, source_dir=None, probe_wsl=True:
                        list(found.get(key, [])))
    monkeypatch.setattr(store, "runs_path", lambda: tmp_path / "runs.json")
    monkeypatch.setattr(store, "user_presets_path",
                        lambda: tmp_path / "presets.json")
    return SimpleNamespace(found=found, tmp=tmp_path)


@pytest.fixture
def win(qapp, env):
    w = InversionWindow()
    w.show()
    qapp.processEvents()
    yield w
    w._console_win.hide()
    w.hide()


@pytest.fixture
def fake(monkeypatch, env):
    eng = FakeEngine()
    monkeypatch.setitem(ie.ENGINES, "occam2d", eng)
    return eng


def _signals(win):
    got = {"state": [], "result": []}
    win.run_state.connect(lambda *a: got["state"].append(a))
    win.result_ready.connect(got["result"].append)
    return got


# ── construction & engines ──────────────────────────────────────────────


class TestConstruction:
    def test_title_and_engine_rows(self, win):
        assert win.windowTitle() == "pycsamt — Inversion Studio"
        keys = [win._engine_list.item(i).data(Qt.ItemDataRole.UserRole)
                for i in range(win._engine_list.count())]
        assert keys == ["occam1d", "occam2d", "modem2d", "modem3d",
                        "mare2dem"]
        assert win._engine_key == "occam2d"

    def test_engine_pills_follow_detected_binaries(self, win):
        st = {k: r.status for k, r in win._engine_rows.items()}
        assert st == {"occam1d": "builtin", "occam2d": "missing",
                      "modem2d": "missing", "modem3d": "ready",
                      "mare2dem": "wsl"}

    def test_stop_disabled_and_monitor_idle(self, win):
        assert not win._btn_stop.isEnabled()
        assert win._monitor.pill.state == "idle"

    def test_library_toggle(self, win):
        win._btn_library.setChecked(False)
        assert not win._lib_panel.isVisible()
        assert win._btn_library.text() == "◂ Library"
        win._btn_library.setChecked(True)
        assert win._lib_panel.isVisible()

    def test_close_hides_and_signals(self, win):
        got = []
        win.panel_closed.connect(lambda: got.append(1))
        win.close()
        assert not win.isVisible() and got == [1]


class TestEngineForms:
    def test_sections_follow_engine(self, win):
        win._select_engine("mare2dem")
        forms = win._engine_forms("mare2dem")
        assert "cell_y" in forms["Mesh"].widgets
        assert "n_procs" in forms["Run"].widgets
        assert win._mesh_host.currentWidget() is forms["Mesh"]
        assert not win._bin_group.isHidden()
        win._select_engine("occam1d")
        assert win._bin_group.isHidden()  # built in: no binary row
        assert win._run_group.isHidden()  # no MPI settings

    def test_values_persist_per_engine_and_band_is_optional(self, win):
        win._select_engine("occam1d")
        win._engine_forms("occam1d")["Mesh"].set_values({"n_layers": 55})
        win._select_engine("occam2d")
        win._select_engine("occam1d")
        v = win.values()
        assert v["n_layers"] == 55 and "freq_min" not in v
        win._chk_band.setChecked(True)
        win._f_min.setValue(0.5)
        assert win.values()["freq_min"] == pytest.approx(0.5)

    def test_advanced_fields_fold_away(self, win):
        form = win._engine_forms("occam2d")["Settings"]
        assert not form.btn_more.isHidden()
        assert form._adv_box.isHidden()
        form.btn_more.setChecked(True)
        assert not form._adv_box.isHidden()

    def test_presets_apply_save_delete(self, win):
        win._select_engine("occam2d")
        win._apply_preset("Fast preview")
        assert win.values()["max_iterations"] == 10
        store.save_preset("occam2d", "Mine", {"n_layers": 44})
        win._populate_presets()
        win._apply_preset("Mine")
        assert win.values()["n_layers"] == 44
        win._preset_list.setCurrentRow(win._preset_list.count() - 1)
        win._delete_preset()
        assert "Mine" not in [win._preset_list.item(i).data(
            Qt.ItemDataRole.UserRole) for i in range(
            win._preset_list.count())]
        win._preset_list.setCurrentRow(0)  # built-ins are protected
        win._delete_preset()
        assert win._preset_list.count() == 3

    def test_build_button_requests_solver_builder(self, win):
        got = []
        win.build_solver_requested.connect(got.append)
        for key in ("occam2d", "modem3d", "mare2dem"):
            win._select_engine(key)
            win._btn_build_solver.click()
        assert got == ["occam2d", "modem3d", "mare2dem"]


# ── data ────────────────────────────────────────────────────────────────


class TestData:
    def test_set_sites_lists_checkable_station_names(self, win):
        """Regression: the old window read ``site.name`` from EDIFile
        objects inside a bare ``except`` -- the list stayed empty and the
        station selection never worked."""
        win.set_sites(_sites())
        names = [win._station_list.item(i).text()
                 for i in range(win._station_list.count())]
        assert names[:3] == ["csa000", "csa050", "csa100"]
        assert win._station_count.text() == "10 of 10 stations used"
        assert win._data_lbl.text() == "10 stations loaded"

    def test_ticked_subset_is_selected(self, win):
        win.set_sites(_sites())
        win._tick_all(False)
        assert win._selected_sites() is None
        win._station_list.item(2).setCheckState(Qt.CheckState.Checked)
        win._station_list.item(5).setCheckState(Qt.CheckState.Checked)
        sub = win._selected_sites()
        assert [s.name for s in sub] == ["csa100", "csa250"]
        win._station_list.item(2).setSelected(True)
        win._tick_highlighted()
        assert win._ticked_names() == ["csa100"]

    def test_data_check_message(self, win):
        win.set_sites(_sites(1))
        win._select_engine("occam2d")
        assert "two stations" in win._data_check.text()
        win._select_engine("occam1d")
        assert win._data_check.text().startswith("✓")

    def test_starting_model_from_forward(self, win):
        win.load_starting_model({"dim": "1D", "resistivity": [42.0, 5.0]})
        assert win._engine_key == "occam1d"
        assert win.values()["starting_resistivity"] == pytest.approx(42.0)
        assert "2 layers" in win._fwd_model_label.text()
        win._clear_starting_model()
        assert win._fwd_model_label.text() == "(no model from Forward)"


# ── build & run: Occam1D for real ───────────────────────────────────────


class TestOccam1DRealRun:
    def test_build_run_and_open_results(self, win, tmp_path):
        got = _signals(win)
        win.set_sites(_sites(2))
        win._select_engine("occam1d")
        win._workdir_edit.setText(str(tmp_path / "o1d"))
        win._engine_forms("occam1d")["Settings"].set_values(
            {"max_iterations": 2})
        win._on_build()
        assert win._build is not None and not win._build_stale
        assert win._mesh_view.showing_canvas
        assert win._monitor.pill.state == "ready"
        assert win._summary_form.rowCount() == 4

        win._on_run()
        assert win._monitor.pill.state == "done"
        labels = {h[0] for h in win._monitor.history}
        assert labels == {"csa000", "csa050"}
        assert win._classical_stack.currentIndex() == 1  # View Results
        assert win._loaded.engine == "occam1d"
        assert win._view_list.count() == 5
        assert win._result_view.showing_canvas
        assert got["result"] and got["result"][-1]["engine"] == "occam1d"
        assert got["state"][0][2] is True and got["state"][-1][2] is False
        assert win._recent_list.count() == 1
        assert "iter" in win._console.text()
        # station / iteration pickers reload
        win._res_station.setCurrentText("csa050")
        assert win._loaded.station == "csa050"

    def test_settings_change_marks_build_stale(self, win, tmp_path):
        win.set_sites(_sites(2))
        win._select_engine("occam1d")
        win._workdir_edit.setText(str(tmp_path / "o1d"))
        win._on_build()
        win._engine_forms("occam1d")["Mesh"].set_values({"n_layers": 22})
        assert win._build_stale
        assert "rebuilt" in win._build_state_lbl.text()


# ── build & run: external solver (fake engine) ──────────────────────────


class TestExternalRun:
    def test_missing_binary_blocks_run(self, fake, win):
        win._select_engine("occam2d")
        win._binary_combo.setEditText("")
        win._on_run()
        assert fake.builds == [] and win._monitor.pill.state == "error"
        assert win._step_stack.currentIndex() == 3

    def test_run_builds_when_stale_and_uses_binary(self, fake, win,
                                                   tmp_path):
        got = _signals(win)
        win._select_engine("occam2d")
        win._workdir_edit.setText(str(tmp_path / "run"))
        win._binary_combo.setEditText("C:/solvers/Occam2D.exe")
        win._on_build()
        win._engine_forms("occam2d")["Mesh"].set_values({"n": 20})
        win._on_run()
        assert [b["n"] for b in fake.builds] == [10, 20]  # rebuilt
        assert fake.binaries == ["C:/solvers/Occam2D.exe"]
        assert win._monitor.pill.state == "done"
        assert win._loaded is not None and win._loaded.engine == "occam2d"
        assert got["result"][-1]["result"] == "RESULT"
        runs = store.recent_runs()
        assert runs[0]["status"] == "done" and runs[0]["final_rms"] == 2.0

    def test_binary_is_staged_into_run_folder(self, fake, win, tmp_path):
        exe = tmp_path / "Occam2D.exe"
        exe.write_bytes(b"exe")
        win._select_engine("occam2d")
        win._workdir_edit.setText(str(tmp_path / "run"))
        win._binary_combo.setEditText(str(exe))
        win._chk_stage.setChecked(True)
        win._on_run()
        staged = Path(fake.binaries[0])
        assert staged.parent == tmp_path / "run" and staged.is_file()
        assert exe.is_file()  # copy, not move

    def test_solver_error_and_cancel(self, fake, win, tmp_path):
        win._select_engine("occam2d")
        win._workdir_edit.setText(str(tmp_path / "run"))
        win._binary_combo.setEditText("C:/s/Occam2D.exe")
        fake.mode = "fail"
        win._on_run()
        assert win._monitor.pill.state == "error"
        assert "solver crashed" in win._console.text()
        assert store.recent_runs()[0]["status"] == "error"
        fake.mode = "cancel"
        win._on_run()
        assert win._monitor.pill.state == "stopped"
        assert win._btn_run.isEnabled() and not win._btn_stop.isEnabled()

    def test_build_error_shows_reason(self, fake, win, tmp_path):
        win._select_engine("occam2d")
        win._workdir_edit.setText(str(tmp_path / "run"))
        fake.mode = "build_fail"
        win._on_build()
        assert win._monitor.pill.state == "error"
        assert not win._mesh_view.showing_canvas
        assert win._build is None

    def test_log_history_is_polled(self, fake, win, tmp_path):
        win._select_engine("occam2d")
        win._run_engine, win._run_workdir = "occam2d", tmp_path
        win._monitor.start("x", max_iter=10, target=1.0)
        fake.hist = [(1, 5.0), (2, 3.0)]
        win._poll_history()
        assert [h[1:] for h in win._monitor.history] == [(1, 5.0), (2, 3.0)]
        win._poll_history()  # no duplicates
        assert len(win._monitor.history) == 2
        assert win._monitor.percent() == 20


class TestBinaryRow:
    def test_user_path_survives_refresh(self, win, env):
        win._select_engine("modem3d")
        assert win._binary_combo.currentText() == "C:/m/Mod3DMT.exe"
        assert "PATH" in win._binary_origin.text()
        win._binary_combo.setEditText("D:/mine/Mod3DMT.exe")
        env.found["modem3d"] = [sb.BinaryCandidate("C:/new/Mod3DMT.exe",
                                                   "Solver Builder")]
        win.refresh_solver_binaries()
        assert win._binary_combo.currentText() == "D:/mine/Mod3DMT.exe"

    def test_auto_path_follows_new_build(self, win, env):
        win._select_engine("occam2d")
        assert win._binary_combo.currentText() == ""
        env.found["occam2d"] = [sb.BinaryCandidate("C:/b/Occam2D.exe",
                                                   "Solver Builder")]
        win.refresh_solver_binaries()
        assert win._binary_combo.currentText() == "C:/b/Occam2D.exe"
        assert win._engine_rows["occam2d"].status == "ready"

    def test_each_solver_keeps_its_own_binary(self, win):
        """Regression (real run): the Occam2D path typed in the box was
        carried to ModEM and MARE2DEM as "your path" when switching
        engines -- ModEM then ran Occam2D.exe with its own arguments."""
        win._select_engine("occam2d")
        win._binary_combo.setEditText("C:/mine/Occam2D.exe")
        win._select_engine("modem3d")
        assert win._binary_combo.currentText() == "C:/m/Mod3DMT.exe"
        win._select_engine("mare2dem")
        assert win._binary_combo.currentText() == "wsl:/opt/MARE2DEM"
        win._select_engine("occam2d")
        assert win._binary_combo.currentText() == "C:/mine/Occam2D.exe"

    def test_wsl_binary_cannot_be_staged(self, win):
        win._select_engine("mare2dem")
        assert win._binary_combo.currentText().startswith("wsl:")
        assert not win._chk_stage.isEnabled()
        assert "WSL" in win._binary_origin.text()


# ── results ─────────────────────────────────────────────────────────────


class TestResults:
    def test_occam2d_folder_views_section_export_send(self, win, tmp_path):
        if not OCCAM2D.is_dir():
            pytest.skip("bundled Occam2D run not found")
        got = _signals(win)
        assert win.open_run(OCCAM2D)
        assert win._loaded.engine == "occam2d"
        keys = [win._view_list.item(i).data(Qt.ItemDataRole.UserRole)
                for i in range(win._view_list.count())]
        assert "section" in keys and "convergence" in keys
        win._view_list.setCurrentRow(keys.index("section"))
        assert win._res_stack.currentWidget() is win._tab_section
        assert win._tab_section._result is win._loaded.result
        win._view_list.setCurrentRow(keys.index("convergence"))
        assert win._result_view.showing_canvas
        assert win._btn_send.isEnabled()
        win._send_to_interpretation()
        assert got["result"][-1]["engine"] == "occam2d"
        out = win.export_pcsf(tmp_path / "occ.pcsf")
        assert out is not None and out.stat().st_size > 0

    def test_modem3d_folder(self, win):
        # The ~11 MB .rho model is not tracked; only .dat/.res are.
        if not MODEM3D.is_dir() or not list(MODEM3D.glob("*.rho")):
            pytest.skip("bundled ModEM result not found")
        assert win.open_run(MODEM3D)
        assert win._loaded.engine == "modem3d"
        assert win._result_view.showing_canvas  # model view
        keys = [win._view_list.item(i).data(Qt.ItemDataRole.UserRole)
                for i in range(win._view_list.count())]
        win._view_list.setCurrentRow(keys.index("convergence"))
        assert not win._result_view.showing_canvas  # no log: reason card
        assert not win._res_iter.isEnabled()

    def test_mare2dem_demo_folder(self, win):
        if not M2D_DEMO.is_dir():
            pytest.skip("bundled MARE2DEM demo not found")
        assert win.open_run(M2D_DEMO)
        for i in range(win._view_list.count()):
            win._view_list.setCurrentRow(i)
            assert win._result_view.showing_canvas

    def test_unknown_folder_explains(self, win, tmp_path):
        assert not win.open_run(tmp_path)
        assert not win._result_view.showing_canvas
        assert "Not an inversion" in win._result_view._unavailable._title \
            .text()


# ── console ─────────────────────────────────────────────────────────────


class TestConsole:
    def test_pop_out_and_dock(self, win):
        win._console.append("hello")
        win._pop_out_console()
        assert win.console_popped_out and win._console_win.isVisible()
        win._console.append("still streaming")
        win._console_win.close()  # closing docks it back
        assert not win.console_popped_out
        assert "still streaming" in win._console.text()

    def test_find_save_theme(self, win, tmp_path):
        for line in ("a", "iter 3 RMS 1.2", "b"):
            win._console.append(line)
        win._console.search.setText("rms 1.2")
        assert win._console.find_next()
        assert not win._console.btn_follow.isChecked()
        p = win._console.save_to(tmp_path / "run.log")
        assert "iter 3 RMS 1.2" in p.read_text()
        win._console.set_theme("amber")
        assert win._console.view.theme == "amber"


# ── session ─────────────────────────────────────────────────────────────


def test_session_round_trip(win, qapp, env):
    win._select_engine("mare2dem")
    win._btn_library.setChecked(False)
    win._console.set_theme("app")
    store_ = {}
    win.save_geometry_to(store_)
    entry = store_["inversion_window"]
    assert entry["engine"] == "mare2dem" and entry["console_theme"] == "app"
    other = InversionWindow()
    other.restore_geometry_from(store_)
    assert other._engine_key == "mare2dem"
    assert other._console.view.theme == "app"
    assert not other._btn_library.isChecked()


# ── AI tab ──────────────────────────────────────────────────────────────


class _FakeSignal:
    def __init__(self):
        self._fns = []

    def connect(self, fn):
        self._fns.append(fn)

    def disconnect(self, *a):
        self._fns = []

    def emit(self, *a):
        for fn in self._fns:
            fn(*a)


class _FakeAIWorker:
    instances: list = []

    def __init__(self, params, parent=None):
        self.params = params
        self.log_line, self.progress = _FakeSignal(), _FakeSignal()
        self.finished, self.error = _FakeSignal(), _FakeSignal()
        _FakeAIWorker.instances.append(self)

    def isRunning(self):
        return False

    def start(self):
        pass


class TestAIPage:
    def test_inv1d_uses_real_stations(self, win):
        """Regression: X_obs was built with a non-existent
        ``site.interpolate_rho_a`` inside a bare except -- every 1-D
        prediction silently ran on synthetic training samples."""
        win.set_sites(_sites(3))
        page = win._ai_page
        p = page._build_ai_params("inv1d")
        assert p["stations"] == ["csa000", "csa050", "csa100"]
        assert p["X_obs"].shape == (3, 60) and np.isfinite(p["X_obs"]).all()

    def test_inv2d_inv3d_attach_sites(self, win):
        win.set_sites(_sites(3))
        page = win._ai_page
        p2 = page._build_ai_params("inv2d")
        assert p2["sites"] is not None and "X_obs" not in p2
        assert p2["physics"] == "mt1d"
        page._ai2_physics.setCurrentIndex(1)
        assert page._build_ai_params("inv2d")["physics"] == "mt2d"
        page._ai3_hidden.setText("64, x")
        assert page._build_ai_params("inv3d")["hidden"] == [256, 128, 64]

    def test_run_and_plot_1d(self, win, monkeypatch):
        import pycsamt.app.desktop.workers.ai_inversion_worker as aim

        _FakeAIWorker.instances = []
        monkeypatch.setattr(aim, "AIInversionWorker", _FakeAIWorker)
        win.set_sites(_sites(2))
        win._mode_group.button(1).click()
        page = win._ai_page
        page.run()
        (worker,) = _FakeAIWorker.instances
        assert worker.params["stations"] == ["csa000", "csa050"]
        n = 5
        y = np.tile(np.r_[np.full(n, 2.0), np.full(n - 1, 50.0)], (2, 1))
        freqs = np.logspace(-3, 3, 30)
        worker.finished.emit({
            "dim": "1D", "y_pred": y, "n_layers": n,
            "X_obs": worker.params["X_obs"], "freqs": freqs,
            "stations": ["csa000", "csa050"],
            "inverter": SimpleNamespace(loss_history_=[1.0, 0.5, 0.2]),
        })
        assert page.monitor.pill.state == "done"
        assert page._tab_model_view.showing_canvas
        assert page._tab_fit_view.showing_canvas
        assert page._tab_convergence_view.showing_canvas
        assert not page._tab_compare_view.showing_canvas  # no classical

    def test_compare_with_occam1d(self, win, tmp_path):
        eng = ie.engine("occam1d")
        vals = eng.defaults()
        vals["max_iterations"] = 2
        build = eng.build(_sites(2), tmp_path, vals)
        eng.make_task(build, vals, None)(ie.RunReporter())
        win.open_run(tmp_path)
        page = win._ai_page
        n = 4
        page._result = {"dim": "1D", "n_layers": n,
                        "y_pred": np.tile(np.r_[np.full(n, 2.0),
                                                np.full(n - 1, 80.0)],
                                          (2, 1)),
                        "stations": ["csa000", "csa050"]}
        page._draw_compare()
        assert page._tab_compare_view.showing_canvas
        assert len(page._tab_compare_view.canvas.figure.axes) == 2

    def test_2d_needs_data(self, win):
        page = win._ai_page
        page._rb_inv2d.setChecked(True)
        page.run()
        assert page.monitor.pill.state == "error"
