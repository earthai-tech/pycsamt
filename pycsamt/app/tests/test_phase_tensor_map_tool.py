# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for PhaseTensorMapDialog and its ``_MapWorker`` (Tools menu).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — ~28 WILLY AMT EDIs with real lat/lon per
station, used here (rather than the KP TIPPER set) because the geographic
map needs real coordinates.

Strategy
--------
* ``_MapWorker.run()`` is called directly (not ``.start()``) with real
  small EDI data to exercise the actual phase-tensor-map computation,
  covering a normal in-range period, extreme/sparse periods at the edges
  of the data's frequency range, and error paths (bad ``sites`` input).
* ``PhaseTensorMapDialog`` orchestration (_build_ui/_on_plot/_on_done/
  _on_error) is exercised with a monkeypatched ``_MapWorker`` whose
  ``.start()`` synchronously emits through plain-callable fake signals,
  mirroring the idiom used in test_main_window_actions.py /
  test_inversion_dlg_extra.py / test_recompute_dlg.py.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.tools.phase_tensor_map_tool import (
    PhaseTensorMapDialog,
    _MapWorker,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"

_HAS_WILLY = _WILLY.exists() and any(_WILLY.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_sites():
    """Small real Sites collection (L18PLT profile) loaded once per session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY))


def _close():
    plt.close("all")


# ── _MapWorker — real computation ───────────────────────────────────────────


class TestMapWorkerRealComputation:
    def test_run_in_range_period_emits_done_with_real_figure(
        self, willy_sites
    ):
        done, err = [], []
        w = _MapWorker(
            willy_sites,
            period=0.05,
            c_by="skew",
            show_tipper=True,
            tipper_conv="parkinson",
            tipper_comp="real",
            station_labels=True,
        )
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()

        assert err == []
        assert len(done) == 1
        fig = done[0]
        assert isinstance(fig, plt.Figure)
        assert len(fig.axes) >= 1
        _close()

    def test_run_alternate_style_options(self, willy_sites):
        """Different c_by/tipper convention/component combo also renders."""
        done, err = [], []
        w = _MapWorker(
            willy_sites,
            period=0.01,
            c_by="ellipt",
            show_tipper=True,
            tipper_conv="wiese",
            tipper_comp="imag",
            station_labels=False,
        )
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()

        assert err == []
        assert len(done) == 1
        assert len(done[0].axes) >= 1
        _close()

    def test_run_no_tipper_overlay(self, willy_sites):
        done, err = [], []
        w = _MapWorker(
            willy_sites,
            period=0.1,
            c_by="theta",
            show_tipper=False,
            tipper_conv="parkinson",
            tipper_comp="amplitude",
            station_labels=True,
        )
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()

        assert err == []
        assert len(done) == 1
        _close()

    def test_run_period_far_outside_data_range_still_renders(
        self, willy_sites
    ):
        """A period way beyond the survey's frequency range (sparse/edge
        case: nearest available period is used per-station) must still
        complete and emit a figure, not raise/error."""
        done, err = [], []
        w = _MapWorker(
            willy_sites,
            period=1e5,
            c_by="skew",
            show_tipper=True,
            tipper_conv="parkinson",
            tipper_comp="real",
            station_labels=True,
        )
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()

        assert err == []
        assert len(done) == 1
        _close()

    def test_run_period_far_below_data_range_still_renders(
        self, willy_sites
    ):
        done, err = [], []
        w = _MapWorker(
            willy_sites,
            period=1e-5,
            c_by="skew",
            show_tipper=True,
            tipper_conv="parkinson",
            tipper_comp="real",
            station_labels=True,
        )
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()

        assert err == []
        assert len(done) == 1
        _close()


class TestMapWorkerErrors:
    def test_run_none_sites_emits_error_not_done(self):
        done, err = [], []
        w = _MapWorker(
            None,
            period=0.05,
            c_by="skew",
            show_tipper=True,
            tipper_conv="parkinson",
            tipper_comp="real",
            station_labels=True,
        )
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()

        assert done == []
        assert len(err) == 1
        assert isinstance(err[0], str) and err[0]

    def test_run_non_numeric_period_emits_error(self, willy_sites):
        """A non-numeric period breaks the internal period-distance
        arithmetic (numpy subtract on a str dtype) — a real error path
        surfaced through the real computation, not a mock."""
        done, err = [], []
        w = _MapWorker(
            willy_sites,
            period="not-a-number",
            c_by="skew",
            show_tipper=True,
            tipper_conv="parkinson",
            tipper_comp="real",
            station_labels=True,
        )
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()

        assert done == []
        assert len(err) == 1


# ── PhaseTensorMapDialog — construction ─────────────────────────────────────


class TestDialogConstruction:
    def test_creates_with_sites(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        assert dlg is not None
        assert dlg.windowTitle() == "Phase Tensor Map"
        dlg.close()

    def test_run_button_enabled_with_sites(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        assert dlg._run_btn.isEnabled()
        dlg.close()

    def test_no_sites_disables_run_button(self, qapp):
        dlg = PhaseTensorMapDialog(sites=None)
        assert not dlg._run_btn.isEnabled()
        assert "No survey" in dlg._status_lbl.text()
        dlg.close()

    def test_period_spin_default(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        assert dlg._period_spin.value() == pytest.approx(10.0)
        dlg.close()

    def test_cby_combo_populated(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        items = [
            dlg._cby_combo.itemText(i)
            for i in range(dlg._cby_combo.count())
        ]
        assert items == ["skew", "ellipt", "theta", "alpha", "s1", "s2"]
        dlg.close()

    def test_tipper_convention_combo_populated(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        items = [
            dlg._tipper_conv_combo.itemText(i)
            for i in range(dlg._tipper_conv_combo.count())
        ]
        assert items == ["parkinson", "wiese"]
        dlg.close()

    def test_tipper_component_combo_populated(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        items = [
            dlg._tipper_comp_combo.itemText(i)
            for i in range(dlg._tipper_comp_combo.count())
        ]
        assert items == ["real", "imag", "amplitude"]
        dlg.close()

    def test_show_tipper_checked_by_default(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        assert dlg._show_tipper_cb.isChecked()
        dlg.close()

    def test_station_labels_checked_by_default(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        assert dlg._labels_cb.isChecked()
        dlg.close()

    def test_worker_none_initially(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        assert dlg._worker is None
        dlg.close()

    def test_close_button_rejects(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        closed = []
        dlg.rejected.connect(lambda: closed.append(True))
        box = dlg.findChild(
            __import__(
                "PySide6.QtWidgets", fromlist=["QDialogButtonBox"]
            ).QDialogButtonBox
        )
        box.rejected.emit()
        assert closed == [True]
        dlg.close()


# ── PhaseTensorMapDialog — orchestration via fake worker ────────────────────


def _fake_worker_cls(*, fig=None, error=None, captured=None):
    class _FakeSignal:
        def __init__(self):
            self._fns = []

        def connect(self, fn):
            self._fns.append(fn)

        def emit(self, *a, **k):
            for fn in list(self._fns):
                fn(*a, **k)

    class _FakeWorker:
        def __init__(self, sites, period, c_by, show_tipper,
                     tipper_conv, tipper_comp, station_labels):
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            self._args = dict(
                sites=sites,
                period=period,
                c_by=c_by,
                show_tipper=show_tipper,
                tipper_conv=tipper_conv,
                tipper_comp=tipper_comp,
                station_labels=station_labels,
            )
            if captured is not None:
                captured.append(self)

        def start(self):
            if error is not None:
                self.error.emit(error)
            else:
                self.done.emit(fig)

    return _FakeWorker


class TestDialogOnPlotHappyPath:
    def test_on_plot_constructs_worker_with_current_ui_state(
        self, qapp, willy_sites, monkeypatch
    ):
        captured: list = []
        real_fig = plt.figure()
        fake_cls = _fake_worker_cls(fig=real_fig, captured=captured)
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.phase_tensor_map_tool._MapWorker",
            fake_cls,
        )

        dlg = PhaseTensorMapDialog(sites=willy_sites)
        dlg._period_spin.setValue(0.07)
        dlg._cby_combo.setCurrentText("ellipt")
        dlg._show_tipper_cb.setChecked(False)
        dlg._tipper_conv_combo.setCurrentText("wiese")
        dlg._tipper_comp_combo.setCurrentText("imag")
        dlg._labels_cb.setChecked(False)

        dlg._on_plot()

        assert len(captured) == 1
        args = captured[0]._args
        assert args["sites"] is willy_sites
        assert args["period"] == pytest.approx(0.07)
        assert args["c_by"] == "ellipt"
        assert args["show_tipper"] is False
        assert args["tipper_conv"] == "wiese"
        assert args["tipper_comp"] == "imag"
        assert args["station_labels"] is False
        dlg.close()
        _close()

    def test_on_plot_disables_run_button_and_updates_status(
        self, qapp, willy_sites, monkeypatch
    ):
        real_fig = plt.figure()
        fake_cls = _fake_worker_cls(fig=real_fig)
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.phase_tensor_map_tool._MapWorker",
            fake_cls,
        )

        dlg = PhaseTensorMapDialog(sites=willy_sites)
        dlg._on_plot()

        # The fake worker's start() ran synchronously and _on_done()
        # already re-enabled the button + set status to "Done."
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Done."
        dlg.close()
        _close()

    def test_on_done_shows_figure_on_canvas(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        real_fig = plt.figure()
        real_fig.add_subplot(111).plot([0, 1], [0, 1])

        dlg._run_btn.setEnabled(False)
        dlg._on_done(real_fig)

        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Done."
        assert dlg._canvas.figure is real_fig
        dlg.close()
        _close()


class TestDialogOnPlotErrorPath:
    def test_on_plot_error_path_updates_status_and_reenables_button(
        self, qapp, willy_sites, monkeypatch
    ):
        fake_cls = _fake_worker_cls(error="boom: bad period")
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.phase_tensor_map_tool._MapWorker",
            fake_cls,
        )

        dlg = PhaseTensorMapDialog(sites=willy_sites)
        dlg._on_plot()

        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: boom: bad period"
        dlg.close()
        _close()

    def test_on_error_directly(self, qapp, willy_sites):
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        dlg._run_btn.setEnabled(False)
        dlg._on_error("some failure")
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: some failure"
        dlg.close()


# ── End-to-end: real worker driven through the dialog (no monkeypatch) ──────


class TestDialogEndToEndRealWorker:
    def test_on_plot_with_real_worker_thread_runs_and_updates_canvas(
        self, qapp, willy_sites
    ):
        """Drive the real (unpatched) _MapWorker via QThread.start() and
        process events until it finishes, exercising the full plumbing
        end-to-end with real data."""
        dlg = PhaseTensorMapDialog(sites=willy_sites)
        dlg._period_spin.setValue(0.05)

        dlg._on_plot()
        worker = dlg._worker
        assert worker is not None
        # Real QThread: wait for it to actually finish.
        finished = worker.wait(15000)
        qapp.processEvents()

        assert finished
        assert dlg._status_lbl.text() == "Done."
        assert dlg._run_btn.isEnabled()
        assert dlg._canvas.figure is not None
        dlg.close()
        _close()
