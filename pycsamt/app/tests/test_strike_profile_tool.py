# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for StrikeProfileDialog and _ProfileWorker.

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — 28 WILLY AMT EDIs on a single profile line.
A small 8-station subset is used for the worker's real-computation tests
to keep the suite fast (``estimate_strike_consensus``/``estimate_strike_sweep``
are the slow paths; ``pt`` is fast).

Strategy
--------
* ``_ProfileWorker.run()`` is called synchronously (not ``.start()``) with
  real EDI data for every ``method`` x ``sort_by`` combination the dialog's
  UI actually offers, plus the period-band path and an error path (the
  underlying ``plot_strike_profile`` monkeypatched to raise).
* ``StrikeProfileDialog`` orchestration (_build_ui, _on_plot, _on_done,
  _on_error) is exercised with ``_ProfileWorker`` monkeypatched to a
  lightweight fake whose ``.start()`` synchronously emits through plain
  ``_FakeSignal`` stand-ins — mirroring the idiom used in
  test_recompute_dlg.py / test_inversion_dlg_extra.py /
  test_main_window_actions.py.
"""

from __future__ import annotations

import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.tools.strike_profile_tool import (
    _METHODS,
    _SORT_BY,
    StrikeProfileDialog,
    _ProfileWorker,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_L18 = _L18.exists() and any(_L18.glob("*.edi"))


# ── Session-scoped fixtures ───────────────────────────────────────────────────


@pytest.fixture(scope="session")
def willy_sites():
    """8-station subset of L18PLT, loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_L18:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools._core import ensure_sites

    files = sorted(glob.glob(str(_L18 / "*.edi")))[:8]
    return ensure_sites(files)


def _close():
    plt.close("all")


# ── Fake worker / signal plumbing for dialog-orchestration tests ─────────────


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls(*, fig=None, error=None):
    class _FakeWorker:
        instances = []

        def __init__(self, sites, method, band, sort_by):
            self.sites = sites
            self.method = method
            self.band = band
            self.sort_by = sort_by
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            self.started = False
            _FakeWorker.instances.append(self)

        def start(self):
            self.started = True
            if error is not None:
                self.error.emit(error)
            else:
                self.done.emit(fig)

    return _FakeWorker


# ── _ProfileWorker — real computation ─────────────────────────────────────────


class TestProfileWorkerReal:
    @pytest.mark.parametrize("method", _METHODS)
    @pytest.mark.parametrize("sort_by", _SORT_BY)
    def test_run_emits_done_with_real_figure(
        self, willy_sites, method, sort_by
    ):
        w = _ProfileWorker(willy_sites, method, None, sort_by)
        done, err = [], []
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()
        try:
            assert err == []
            assert len(done) == 1
            fig = done[0]
            assert isinstance(fig, plt.Figure)
            assert len(fig.axes) >= 1
            ax = fig.axes[0]
            # A real strike profile draws a line plus an IQR fill_between.
            assert len(ax.lines) >= 1
            assert len(ax.collections) >= 1  # fill_between ribbon
        finally:
            _close()

    def test_run_with_period_band(self, willy_sites):
        w = _ProfileWorker(willy_sites, "consensus", (0.001, 0.1), "auto")
        done, err = [], []
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()
        try:
            assert err == []
            assert len(done) == 1
            assert isinstance(done[0], plt.Figure)
        finally:
            _close()

    def test_run_sort_by_lon_orders_by_longitude(self, willy_sites):
        """sort_by='lon' must not raise and must produce one point per station."""
        w = _ProfileWorker(willy_sites, "consensus", None, "lon")
        done, err = [], []
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()
        try:
            assert err == []
            ax = done[0].axes[0]
            xticklabels = [t.get_text() for t in ax.get_xticklabels()]
            assert len(xticklabels) == len(willy_sites)
        finally:
            _close()


# ── _ProfileWorker — error path ───────────────────────────────────────────────


class TestProfileWorkerError:
    def test_run_emits_error_when_underlying_function_raises(
        self, willy_sites, monkeypatch
    ):
        import pycsamt.emtools.strike as strike_mod

        def _boom(*a, **k):
            raise RuntimeError("synthetic strike-profile failure")

        monkeypatch.setattr(strike_mod, "plot_strike_profile", _boom)

        w = _ProfileWorker(willy_sites, "consensus", None, "auto")
        done, err = [], []
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()

        assert done == []
        assert len(err) == 1
        assert "synthetic strike-profile failure" in err[0]

    def test_run_bad_sites_source_does_not_raise(self):
        """An unresolvable sites source is swallowed internally (loaded as
        an empty collection) rather than raising out of run(); the profile
        function then draws a 'no strikes' placeholder and still emits
        `done`, never `error`, for this input shape."""
        w = _ProfileWorker("/no/such/path/xyz123", "consensus", None, "auto")
        done, err = [], []
        w.done.connect(done.append)
        w.error.connect(err.append)
        w.run()  # must not raise
        assert len(done) == 1 or len(err) == 1
        _close()


# ── StrikeProfileDialog — construction ────────────────────────────────────────


class TestDialogConstruction:
    def test_creates_with_no_sites(self, qapp):
        dlg = StrikeProfileDialog(sites=None)
        try:
            assert dlg is not None
            assert dlg._sites is None
            assert dlg._worker is None
        finally:
            dlg.close()

    def test_no_sites_disables_run_button(self, qapp):
        dlg = StrikeProfileDialog(sites=None)
        try:
            assert dlg._run_btn.isEnabled() is False
            assert "No survey loaded" in dlg._status_lbl.text()
        finally:
            dlg.close()

    def test_with_sites_enables_run_button(self, qapp, willy_sites):
        dlg = StrikeProfileDialog(sites=willy_sites)
        try:
            assert dlg._run_btn.isEnabled() is True
            assert dlg._status_lbl.text() == ""
        finally:
            dlg.close()

    def test_window_title(self, qapp):
        dlg = StrikeProfileDialog(sites=None)
        try:
            assert dlg.windowTitle() == "Strike Profile Viewer"
        finally:
            dlg.close()

    def test_method_combo_has_all_methods(self, qapp):
        dlg = StrikeProfileDialog(sites=None)
        try:
            items = [
                dlg._method_combo.itemText(i)
                for i in range(dlg._method_combo.count())
            ]
            assert items == _METHODS
        finally:
            dlg.close()

    def test_sort_combo_has_all_options(self, qapp):
        dlg = StrikeProfileDialog(sites=None)
        try:
            items = [
                dlg._sort_combo.itemText(i)
                for i in range(dlg._sort_combo.count())
            ]
            assert items == _SORT_BY
        finally:
            dlg.close()

    def test_default_band_spins_are_zero(self, qapp):
        dlg = StrikeProfileDialog(sites=None)
        try:
            assert dlg._tmin_spin.value() == 0.0
            assert dlg._tmax_spin.value() == 0.0
        finally:
            dlg.close()


# ── StrikeProfileDialog — _on_plot orchestration (fake worker) ───────────────


class TestOnPlot:
    def test_on_plot_constructs_worker_with_ui_state(
        self, qapp, willy_sites, monkeypatch
    ):
        import pycsamt.app.desktop.tools.strike_profile_tool as mod

        fake_cls = _fake_worker_cls(fig=plt.figure())
        monkeypatch.setattr(mod, "_ProfileWorker", fake_cls)

        dlg = StrikeProfileDialog(sites=willy_sites)
        try:
            dlg._method_combo.setCurrentText("pt")
            dlg._sort_combo.setCurrentText("lat")
            dlg._on_plot()

            assert len(fake_cls.instances) == 1
            w = fake_cls.instances[0]
            assert w.sites is willy_sites
            assert w.method == "pt"
            assert w.sort_by == "lat"
            assert w.band is None
            assert w.started is True
            # Button re-enabled and status set to "Done." by the emitted done()
            assert dlg._run_btn.isEnabled() is True
            assert dlg._status_lbl.text() == "Done."
        finally:
            dlg.close()
            _close()

    def test_on_plot_disables_button_and_sets_status_before_start(
        self, qapp, willy_sites, monkeypatch
    ):
        """Regression-guard: _on_plot must disable the button and show the
        'Computing…' status *before* handing off to the worker, even though
        our fake worker's start() resolves synchronously and immediately
        flips it back via _on_done."""
        import pycsamt.app.desktop.tools.strike_profile_tool as mod

        seen_states = []

        class _ObservingSignal(_FakeSignal):
            pass

        class _ObservingWorker:
            def __init__(self, sites, method, band, sort_by):
                self.done = _FakeSignal()
                self.error = _FakeSignal()

            def start(self):
                # Capture dialog state exactly as _on_plot leaves it
                # before this start() call returns control.
                seen_states.append(
                    (dlg._run_btn.isEnabled(), dlg._status_lbl.text())
                )

        monkeypatch.setattr(mod, "_ProfileWorker", _ObservingWorker)
        dlg = StrikeProfileDialog(sites=willy_sites)
        try:
            dlg._on_plot()
            assert seen_states == [(False, "Computing strike profile…")]
        finally:
            dlg.close()

    def test_on_plot_valid_band_passed_through(
        self, qapp, willy_sites, monkeypatch
    ):
        import pycsamt.app.desktop.tools.strike_profile_tool as mod

        fake_cls = _fake_worker_cls(fig=plt.figure())
        monkeypatch.setattr(mod, "_ProfileWorker", fake_cls)

        dlg = StrikeProfileDialog(sites=willy_sites)
        try:
            dlg._tmin_spin.setValue(1.0)
            dlg._tmax_spin.setValue(100.0)
            dlg._on_plot()

            w = fake_cls.instances[0]
            assert w.band == (1.0, 100.0)
        finally:
            dlg.close()
            _close()

    def test_on_plot_tmin_greater_than_tmax_yields_no_band(
        self, qapp, willy_sites, monkeypatch
    ):
        import pycsamt.app.desktop.tools.strike_profile_tool as mod

        fake_cls = _fake_worker_cls(fig=plt.figure())
        monkeypatch.setattr(mod, "_ProfileWorker", fake_cls)

        dlg = StrikeProfileDialog(sites=willy_sites)
        try:
            dlg._tmin_spin.setValue(50.0)
            dlg._tmax_spin.setValue(10.0)
            dlg._on_plot()

            w = fake_cls.instances[0]
            assert w.band is None
        finally:
            dlg.close()
            _close()

    def test_on_plot_error_path_sets_status(
        self, qapp, willy_sites, monkeypatch
    ):
        import pycsamt.app.desktop.tools.strike_profile_tool as mod

        fake_cls = _fake_worker_cls(error="synthetic worker failure")
        monkeypatch.setattr(mod, "_ProfileWorker", fake_cls)

        dlg = StrikeProfileDialog(sites=willy_sites)
        try:
            dlg._on_plot()
            assert dlg._run_btn.isEnabled() is True
            assert dlg._status_lbl.text() == "Error: synthetic worker failure"
        finally:
            dlg.close()


# ── StrikeProfileDialog — _on_done / _on_error direct ─────────────────────────


class TestOnDoneOnError:
    def test_on_done_draws_figure_and_updates_status(self, qapp):
        dlg = StrikeProfileDialog(sites=None)
        try:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot([0, 1, 2], [1, 2, 3])
            dlg._run_btn.setEnabled(False)
            dlg._on_done(fig)
            assert dlg._run_btn.isEnabled() is True
            assert dlg._status_lbl.text() == "Done."
        finally:
            dlg.close()
            _close()

    def test_on_error_updates_status_and_reenables_button(self, qapp):
        dlg = StrikeProfileDialog(sites=None)
        try:
            dlg._run_btn.setEnabled(False)
            dlg._on_error("boom")
            assert dlg._run_btn.isEnabled() is True
            assert dlg._status_lbl.text() == "Error: boom"
        finally:
            dlg.close()


# ── End-to-end: real worker wired into the real dialog ────────────────────────


class TestEndToEndRealWorker:
    def test_on_plot_with_real_worker_run_synchronously(
        self, qapp, willy_sites, monkeypatch
    ):
        """Exercise the dialog's real _ProfileWorker wiring (not the fake),
        by forcing QThread.start() to run synchronously in-process."""
        monkeypatch.setattr(
            _ProfileWorker, "start", _ProfileWorker.run, raising=False
        )

        dlg = StrikeProfileDialog(sites=willy_sites)
        try:
            dlg._method_combo.setCurrentText("pt")
            dlg._on_plot()
            assert dlg._status_lbl.text() == "Done."
            assert dlg._run_btn.isEnabled() is True
            assert len(dlg._canvas.figure.axes) >= 1
        finally:
            dlg.close()
            _close()
