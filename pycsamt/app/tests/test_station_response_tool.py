# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for StationResponseDialog and _ResponseWorker.

Real data
---------
data/MT/kap03lmt_edis/       — 26 KP TIPPER EDIs (has real tipper data)
data/AMT/WILLY_DATA/L18PLT/  — ~28 WILLY AMT EDIs (one profile line)

Strategy
--------
* ``_ResponseWorker.run()`` is called directly (never ``.start()``) with
  real EDI data to exercise the real ``plot_station_response`` computation,
  per the ``InversionWorker`` / ``_StrikeWorker`` precedent in
  ``test_inversion_dlg.py`` / ``test_strike_tool.py``.
* Dialog-level orchestration (``_build_ui`` / ``_populate_stations`` /
  ``_on_plot`` / ``_on_done`` / ``_on_error``) is tested by monkeypatching
  ``_ResponseWorker`` with a lightweight fake whose ``.start()``
  synchronously emits through plain-callable "signals", mirroring the
  ``_FakeSignal`` idiom used in ``test_main_window_actions.py`` /
  ``test_recompute_dlg.py`` / ``test_strike_tool.py``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.tools.station_response_tool import (
    StationResponseDialog,
    _ResponseWorker,
)

# ── Paths ────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_TIPPER = _ROOT / "data" / "MT" / "kap03lmt_edis"
_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"

_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))
_HAS_L18 = _L18.exists() and any(_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def tipper_sites():
    """26-station TIPPER Sites loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_TIPPER:
        pytest.skip("TIPPER data not available")
    from pycsamt.emtools._core import ensure_sites

    return ensure_sites(str(_TIPPER))


@pytest.fixture(scope="session")
def l18_sites():
    """~28-station WILLY L18PLT Sites loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_L18:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools._core import ensure_sites

    return ensure_sites(str(_L18))


def _first_station_name(sites):
    from pycsamt.emtools._core import _iter_items, _name

    for i, ed in enumerate(_iter_items(sites)):
        return _name(ed, i)
    raise AssertionError("sites collection is empty")


def _close():
    plt.close("all")


# ── _FakeSignal / fake worker (dialog-level tests) ─────────────────────────


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls(*, fig=None, error=None):
    """Build a fake ``_ResponseWorker`` whose ``.start()`` fires synchronously."""

    class _FakeWorker:
        captured = []

        def __init__(self, sites, station, components, show_tipper, show_errors):
            self.sites = sites
            self.station = station
            self.components = components
            self.show_tipper = show_tipper
            self.show_errors = show_errors
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            _FakeWorker.captured.append(self)

        def start(self):
            if error is not None:
                self.error.emit(error)
            else:
                self.done.emit(fig)

    return _FakeWorker


def _sample_fig():
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 4, 9])
    return fig


# ── _ResponseWorker — real computation ──────────────────────────────────────


class TestResponseWorkerRealTipper:
    def test_run_emits_done_with_figure(self, qapp, tipper_sites):
        station = _first_station_name(tipper_sites)
        done_calls = []
        error_calls = []
        w = _ResponseWorker(
            tipper_sites,
            station,
            ("xx", "xy", "yx", "yy"),
            True,
            True,
        )
        w.done.connect(done_calls.append)
        w.error.connect(error_calls.append)
        w.run()

        assert error_calls == []
        assert len(done_calls) == 1
        fig = done_calls[0]
        assert isinstance(fig, plt.Figure)
        assert len(fig.axes) > 0
        # real content was drawn (lines and/or error-bar containers)
        n_lines = sum(len(ax.lines) for ax in fig.axes)
        n_containers = sum(len(ax.containers) for ax in fig.axes)
        assert n_lines > 0 or n_containers > 0
        _close()

    def test_run_subset_components_no_tipper(self, qapp, tipper_sites):
        station = _first_station_name(tipper_sites)
        done_calls = []
        error_calls = []
        w = _ResponseWorker(
            tipper_sites,
            station,
            ("xy", "yx"),
            False,
            True,
        )
        w.done.connect(done_calls.append)
        w.error.connect(error_calls.append)
        w.run()

        assert error_calls == []
        assert len(done_calls) == 1
        fig = done_calls[0]
        assert isinstance(fig, plt.Figure)
        assert len(fig.axes) > 0
        _close()

    def test_run_no_error_bars(self, qapp, tipper_sites):
        station = _first_station_name(tipper_sites)
        done_calls = []
        w = _ResponseWorker(
            tipper_sites,
            station,
            ("xx", "xy", "yx", "yy"),
            True,
            False,
        )
        w.done.connect(done_calls.append)
        w.error.connect(lambda *_: None)
        w.run()
        assert len(done_calls) == 1
        assert isinstance(done_calls[0], plt.Figure)
        _close()


class TestResponseWorkerRealWilly:
    def test_run_emits_done_with_figure(self, qapp, l18_sites):
        station = _first_station_name(l18_sites)
        done_calls = []
        error_calls = []
        w = _ResponseWorker(
            l18_sites,
            station,
            ("xx", "xy", "yx", "yy"),
            True,
            True,
        )
        w.done.connect(done_calls.append)
        w.error.connect(error_calls.append)
        w.run()

        assert error_calls == []
        assert len(done_calls) == 1
        fig = done_calls[0]
        assert isinstance(fig, plt.Figure)
        assert len(fig.axes) > 0
        _close()

    def test_run_single_component(self, qapp, l18_sites):
        station = _first_station_name(l18_sites)
        done_calls = []
        w = _ResponseWorker(
            l18_sites,
            station,
            ("xy",),
            True,
            True,
        )
        w.done.connect(done_calls.append)
        w.error.connect(lambda *_: None)
        w.run()
        assert len(done_calls) == 1
        _close()


class TestResponseWorkerError:
    def test_run_unknown_station_emits_error(self, qapp, l18_sites):
        done_calls = []
        error_calls = []
        w = _ResponseWorker(
            l18_sites,
            "NOT_A_REAL_STATION_XYZ",
            ("xx", "xy", "yx", "yy"),
            True,
            True,
        )
        w.done.connect(done_calls.append)
        w.error.connect(error_calls.append)
        w.run()

        assert done_calls == []
        assert len(error_calls) == 1
        assert "not found" in error_calls[0].lower()

    def test_run_generic_exception_emits_error(self, qapp, monkeypatch):
        import pycsamt.emtools.inspect as inspect_mod

        def _boom(*a, **k):
            raise RuntimeError("synthetic failure")

        monkeypatch.setattr(inspect_mod, "plot_station_response", _boom)
        done_calls = []
        error_calls = []
        w = _ResponseWorker(
            object(),
            "whatever",
            ("xy",),
            True,
            True,
        )
        w.done.connect(done_calls.append)
        w.error.connect(error_calls.append)
        w.run()

        assert done_calls == []
        assert len(error_calls) == 1
        assert "synthetic failure" in error_calls[0]


# ── StationResponseDialog — construction ────────────────────────────────────


@pytest.fixture
def dlg(qapp):
    d = StationResponseDialog(sites=None)
    yield d
    d.close()


class TestDialogBuildNoSites:
    def test_creates(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Station Response Inspector"

    def test_initial_state(self, dlg):
        assert dlg._worker is None
        assert dlg._station_names == []

    def test_no_sites_status_message(self, dlg):
        assert "No survey loaded" in dlg._status_lbl.text()

    def test_no_sites_run_button_disabled(self, dlg):
        assert not dlg._run_btn.isEnabled()

    def test_default_checkbox_state(self, dlg):
        assert not dlg._cb_xx.isChecked()
        assert dlg._cb_xy.isChecked()
        assert dlg._cb_yx.isChecked()
        assert not dlg._cb_yy.isChecked()
        assert dlg._cb_tipper.isChecked()
        assert dlg._cb_errors.isChecked()

    def test_canvas_exists(self, dlg):
        assert dlg._canvas is not None

    def test_station_combo_empty(self, dlg):
        assert dlg._station_combo.count() == 0


# ── _populate_stations — real sites ─────────────────────────────────────────


class TestPopulateStationsReal:
    def test_populate_from_tipper_sites(self, qapp, tipper_sites):
        d = StationResponseDialog(sites=tipper_sites)
        assert d._station_combo.count() == len(d._station_names)
        assert d._station_combo.count() > 0
        assert f"{len(d._station_names)} stations loaded." in d._status_lbl.text()
        assert d._run_btn.isEnabled()
        d.close()

    def test_populate_from_l18_sites(self, qapp, l18_sites):
        d = StationResponseDialog(sites=l18_sites)
        assert d._station_combo.count() == len(d._station_names)
        assert d._station_combo.count() > 0
        first_name = _first_station_name(l18_sites)
        assert d._station_combo.itemText(0) == first_name
        d.close()

    def test_populate_stations_raises_disables_button(self, qapp, monkeypatch):
        def _boom(_sites):
            raise RuntimeError("cannot iterate")

        # ``_populate_stations`` does a local ``from pycsamt.emtools._core
        # import _iter_items`` on every call, so patching the source module
        # attribute is enough to reach the except-branch.
        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "_iter_items", _boom)

        d = StationResponseDialog(sites=object())
        assert "Cannot read stations" in d._status_lbl.text()
        assert not d._run_btn.isEnabled()
        d.close()

    def test_populate_stations_unwrap_failure_falls_back_to_raw_item(
        self, qapp, monkeypatch
    ):
        """
        ``_unwrap`` itself has its own internal try/except and never raises
        in practice, so the per-item ``except Exception: pass`` around it
        is normally dead code; force it here to confirm the loop degrades
        to using the raw (un-unwrapped) item instead of aborting.
        """
        import pycsamt.emtools._core as core_mod

        monkeypatch.setattr(core_mod, "_iter_items", lambda _sites: ["ST-RAW"])

        def _boom(_ed):
            raise RuntimeError("cannot unwrap")

        monkeypatch.setattr(core_mod, "_unwrap", _boom)

        d = StationResponseDialog(sites=object())
        assert d._station_names == ["?"]  # "ST-RAW" str has no .station/.id
        d.close()


# ── _on_plot — no data / validation paths ───────────────────────────────────


class TestOnPlotValidation:
    def test_no_station_selected_shows_message(self, qapp):
        d = StationResponseDialog(sites=None)
        # combo is empty -> currentText() == ""
        d._on_plot()
        assert "Select a station first." in d._status_lbl.text()
        assert d._worker is None
        d.close()

    def test_no_components_selected_shows_message(self, qapp, l18_sites):
        d = StationResponseDialog(sites=l18_sites)
        for cb in (d._cb_xx, d._cb_xy, d._cb_yx, d._cb_yy):
            cb.setChecked(False)
        d._on_plot()
        assert "Select at least one component." in d._status_lbl.text()
        assert d._worker is None
        d.close()


# ── _on_plot — happy path via fake worker ───────────────────────────────────


class TestOnPlotWithFakeWorker:
    def test_happy_path_starts_worker_and_draws(self, qapp, monkeypatch):
        fig = _sample_fig()
        fake_cls = _fake_worker_cls(fig=fig)
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.station_response_tool._ResponseWorker",
            fake_cls,
        )
        d = StationResponseDialog(sites=object())
        d._station_names = ["S01", "S02"]
        d._station_combo.addItems(d._station_names)
        d._run_btn.setEnabled(True)

        d._on_plot()

        assert d._worker is not None
        assert d._run_btn.isEnabled()  # re-enabled by _on_done
        assert "Done." in d._status_lbl.text()
        _close()
        d.close()

    def test_worker_constructed_with_form_values(self, qapp, monkeypatch):
        fake_cls = _fake_worker_cls(fig=_sample_fig())
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.station_response_tool._ResponseWorker",
            fake_cls,
        )
        d = StationResponseDialog(sites=object())
        d._station_names = ["S01"]
        d._station_combo.addItems(d._station_names)
        d._run_btn.setEnabled(True)

        d._cb_xx.setChecked(True)
        d._cb_xy.setChecked(True)
        d._cb_yx.setChecked(False)
        d._cb_yy.setChecked(False)
        d._cb_tipper.setChecked(False)
        d._cb_errors.setChecked(False)

        d._on_plot()

        w = fake_cls.captured[-1]
        assert w.station == "S01"
        assert w.components == ("xx", "xy")
        assert w.show_tipper is False
        assert w.show_errors is False
        assert w.sites is d._sites
        _close()
        d.close()

    def test_error_path_updates_status_and_reenables(self, qapp, monkeypatch):
        fake_cls = _fake_worker_cls(error="boom happened")
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.station_response_tool._ResponseWorker",
            fake_cls,
        )
        d = StationResponseDialog(sites=object())
        d._station_names = ["S01"]
        d._station_combo.addItems(d._station_names)
        d._run_btn.setEnabled(True)

        d._on_plot()

        assert d._run_btn.isEnabled()
        assert "Error: boom happened" in d._status_lbl.text()
        d.close()

    def test_run_button_disabled_and_status_computing_during_run(
        self, qapp, monkeypatch
    ):
        """The 'Computing…' status is set before start(); fake worker
        resolves synchronously so by the time _on_plot returns the button
        is re-enabled again — verify the transient state via a spy worker."""

        class _SpyWorker(_fake_worker_cls(fig=_sample_fig())):
            def start(self):
                _SpyWorker.status_during_run = self._dlg_ref._status_lbl.text()
                _SpyWorker.btn_enabled_during_run = self._dlg_ref._run_btn.isEnabled()
                super().start()

        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.station_response_tool._ResponseWorker",
            _SpyWorker,
        )
        d = StationResponseDialog(sites=object())
        d._station_names = ["S01"]
        d._station_combo.addItems(d._station_names)
        d._run_btn.setEnabled(True)
        _SpyWorker._dlg_ref = d

        d._on_plot()

        assert "Computing" in _SpyWorker.status_during_run
        assert _SpyWorker.btn_enabled_during_run is False
        _close()
        d.close()


# ── _on_done / _on_error direct calls ───────────────────────────────────────


class TestOnDoneDirect:
    def test_on_done_draws_figure_and_reenables(self, dlg):
        dlg._run_btn.setEnabled(False)
        fig = _sample_fig()
        dlg._on_done(fig)
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Done."
        _close()

    def test_on_error_sets_message_and_reenables(self, dlg):
        dlg._run_btn.setEnabled(False)
        dlg._on_error("bad things happened")
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: bad things happened"


# ── End-to-end: real worker output feeds the dialog canvas ─────────────────


class TestEndToEndReal:
    def test_real_worker_output_feeds_on_done(self, dlg, l18_sites):
        station = _first_station_name(l18_sites)
        w = _ResponseWorker(
            l18_sites,
            station,
            ("xx", "xy", "yx", "yy"),
            True,
            True,
        )
        collected = []
        w.done.connect(collected.append)
        w.error.connect(lambda *_: None)
        w.run()

        assert len(collected) == 1
        dlg._on_done(collected[0])
        assert dlg._status_lbl.text() == "Done."
        _close()
