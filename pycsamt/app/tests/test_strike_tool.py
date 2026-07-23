# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for StrikeAnalyzerDialog and _StrikeWorker.

Real data
---------
data/AMT/WILLY_DATA/L18PLT/  — ~28 WILLY AMT EDIs (one profile line)

Strategy
--------
* ``_StrikeWorker.run()`` is called directly (never ``.start()``) with real
  EDI data to exercise the real ``estimate_strike_consensus`` computation
  and its DataFrame-shape contract, per the ``InversionWorker`` precedent
  in ``test_inversion_dlg.py``.
* Dialog-level orchestration (``_on_run`` / ``_on_done`` / ``_on_error`` /
  ``_fill_table`` / ``_draw_rose``) is tested by monkeypatching
  ``_StrikeWorker`` with a lightweight fake whose ``.start()`` synchronously
  emits through plain-callable "signals", mirroring the ``_FakeSignal``
  idiom used in ``test_main_window_actions.py`` / ``test_recompute_dlg.py``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.tools.strike_tool import (
    StrikeAnalyzerDialog,
    _StrikeWorker,
)

# ── Paths ────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"

_HAS_L18 = _L18.exists() and any(_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def l18_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_L18:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools._core import ensure_sites

    return ensure_sites(str(_L18))


# ── _FakeSignal / fake worker (dialog-level tests) ─────────────────────────


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls(*, df=None, error=None):
    """Build a fake ``_StrikeWorker`` whose ``.start()`` fires synchronously."""

    class _FakeWorker:
        captured = []

        def __init__(self, sites, band, w_sweep, w_pt):
            self.sites = sites
            self.band = band
            self.w_sweep = w_sweep
            self.w_pt = w_pt
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            _FakeWorker.captured.append(self)

        def start(self):
            if error is not None:
                self.error.emit(error)
            else:
                self.done.emit(df)

    return _FakeWorker


def _sample_df(n=3):
    return pd.DataFrame(
        {
            "station": [f"S{i}" for i in range(n)],
            "ang": [10.0 + i for i in range(n)],
            "iqr": [2.0 + i * 0.5 for i in range(n)],
            "lo": [0.001] * n,
            "hi": [1000.0] * n,
            "n": [5 + i for i in range(n)],
        }
    )


def _empty_df():
    return pd.DataFrame(columns=["station", "ang", "iqr", "lo", "hi", "n"])


# ── _StrikeWorker — real computation ────────────────────────────────────────


class TestStrikeWorkerReal:
    def test_run_emits_done_with_dataframe(self, qapp, l18_sites):
        done_calls = []
        error_calls = []
        w = _StrikeWorker(
            l18_sites, band=None, w_sweep=0.5, w_pt=0.5
        )
        w.done.connect(done_calls.append)
        w.error.connect(error_calls.append)
        w.run()

        assert error_calls == []
        assert len(done_calls) == 1
        df = done_calls[0]
        assert isinstance(df, pd.DataFrame)
        for col in ("station", "ang", "iqr", "lo", "hi", "n"):
            assert col in df.columns
        assert len(df) > 0

    def test_run_with_band_restricts_output(self, qapp, l18_sites):
        done_calls = []
        w = _StrikeWorker(
            l18_sites, band=(0.001, 1000.0), w_sweep=0.5, w_pt=0.5
        )
        w.done.connect(done_calls.append)
        w.error.connect(lambda *_: None)
        w.run()
        assert len(done_calls) == 1
        assert isinstance(done_calls[0], pd.DataFrame)

    def test_run_weights_zero_and_one(self, qapp, l18_sites):
        """w_sweep=1, w_pt=0 and vice versa should both run without error."""
        for w_sweep, w_pt in ((1.0, 0.0), (0.0, 1.0)):
            done_calls = []
            error_calls = []
            w = _StrikeWorker(
                l18_sites, band=None, w_sweep=w_sweep, w_pt=w_pt
            )
            w.done.connect(done_calls.append)
            w.error.connect(error_calls.append)
            w.run()
            assert error_calls == []
            assert len(done_calls) == 1


class TestStrikeWorkerError:
    def test_run_emits_error_on_exception(self, qapp, monkeypatch):
        import pycsamt.emtools.strike as strike_mod

        def _boom(*a, **k):
            raise RuntimeError("synthetic failure")

        monkeypatch.setattr(
            strike_mod, "estimate_strike_consensus", _boom
        )
        done_calls = []
        error_calls = []
        w = _StrikeWorker(object(), band=None, w_sweep=0.5, w_pt=0.5)
        w.done.connect(done_calls.append)
        w.error.connect(error_calls.append)
        w.run()

        assert done_calls == []
        assert len(error_calls) == 1
        assert "synthetic failure" in error_calls[0]

    def test_run_degenerate_input_emits_empty_dataframe(self, qapp):
        """A bare ``object()`` normalizes to a Sites with no items rather
        than raising — the worker emits ``done`` with an empty DataFrame,
        not ``error``."""
        done_calls = []
        error_calls = []
        w = _StrikeWorker(object(), band=None, w_sweep=0.5, w_pt=0.5)
        w.done.connect(done_calls.append)
        w.error.connect(error_calls.append)
        w.run()

        assert error_calls == []
        assert len(done_calls) == 1
        assert isinstance(done_calls[0], pd.DataFrame)
        assert done_calls[0].empty


# ── StrikeAnalyzerDialog — construction ─────────────────────────────────────


@pytest.fixture
def dlg(qapp):
    d = StrikeAnalyzerDialog(sites=None)
    yield d
    d.close()


class TestDialogBuild:
    def test_creates(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Strike Analyzer"

    def test_initial_state(self, dlg):
        assert dlg._worker is None
        assert dlg._df is None

    def test_default_band_values(self, dlg):
        assert dlg._band_min.value() == pytest.approx(0.001)
        assert dlg._band_max.value() == pytest.approx(1000.0)

    def test_default_weight_values(self, dlg):
        assert dlg._w_sweep.value() == pytest.approx(0.5)
        assert dlg._w_pt.value() == pytest.approx(0.5)

    def test_run_button_exists_and_enabled(self, dlg):
        assert dlg._run_btn is not None
        assert dlg._run_btn.isEnabled()

    def test_status_label_initially_empty(self, dlg):
        assert dlg._status_lbl.text() == ""

    def test_table_has_four_columns(self, dlg):
        assert dlg._table.columnCount() == 4

    def test_table_header_labels(self, dlg):
        labels = [
            dlg._table.horizontalHeaderItem(i).text() for i in range(4)
        ]
        assert labels == ["Station", "Strike (°)", "IQR (°)", "n obs"]

    def test_table_initially_empty(self, dlg):
        assert dlg._table.rowCount() == 0

    def test_canvas_exists(self, dlg):
        assert dlg._canvas is not None
        assert dlg._canvas.figure is not None


# ── _on_run ──────────────────────────────────────────────────────────────────


class TestOnRunNoData:
    def test_no_sites_warns_and_returns(self, dlg, monkeypatch):
        warned = []
        from PySide6.QtWidgets import QMessageBox

        monkeypatch.setattr(
            QMessageBox,
            "warning",
            staticmethod(lambda *a, **k: warned.append(a)),
        )
        dlg._on_run()
        assert len(warned) == 1
        assert dlg._worker is None


class TestOnRunWithFakeWorker:
    def test_happy_path_starts_worker_and_populates(
        self, qapp, monkeypatch
    ):
        df = _sample_df(3)
        fake_cls = _fake_worker_cls(df=df)
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.strike_tool._StrikeWorker",
            fake_cls,
        )
        d = StrikeAnalyzerDialog(sites=object())
        d._on_run()

        assert d._worker is not None
        assert d._run_btn.isEnabled()  # re-enabled by _on_done
        assert d._df is df
        assert d._table.rowCount() == 3
        assert "Done" in d._status_lbl.text()
        assert len(d._canvas.figure.axes) >= 1
        d.close()

    def test_worker_constructed_with_form_values(self, qapp, monkeypatch):
        fake_cls = _fake_worker_cls(df=_sample_df(1))
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.strike_tool._StrikeWorker",
            fake_cls,
        )
        d = StrikeAnalyzerDialog(sites=object())
        d._band_min.setValue(0.01)
        d._band_max.setValue(500.0)
        d._w_sweep.setValue(0.3)
        d._w_pt.setValue(0.7)
        d._on_run()

        w = fake_cls.captured[-1]
        assert w.band == (pytest.approx(0.01), pytest.approx(500.0))
        assert w.w_sweep == pytest.approx(0.3)
        assert w.w_pt == pytest.approx(0.7)
        assert w.sites is d._sites
        d.close()

    def test_error_path_updates_status_and_reenables(
        self, qapp, monkeypatch
    ):
        fake_cls = _fake_worker_cls(error="boom happened")
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.strike_tool._StrikeWorker",
            fake_cls,
        )
        d = StrikeAnalyzerDialog(sites=object())
        d._on_run()

        assert d._run_btn.isEnabled()
        assert "Error" in d._status_lbl.text()
        assert "boom happened" in d._status_lbl.text()
        d.close()

    def test_empty_dataframe_result_shows_message_no_table_fill(
        self, qapp, monkeypatch
    ):
        fake_cls = _fake_worker_cls(df=_empty_df())
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.strike_tool._StrikeWorker",
            fake_cls,
        )
        d = StrikeAnalyzerDialog(sites=object())
        d._on_run()

        assert d._df is not None
        assert d._df.empty
        assert "No strike estimates" in d._status_lbl.text()
        assert d._table.rowCount() == 0
        d.close()

    def test_run_button_disabled_during_run_status_running(
        self, qapp, monkeypatch
    ):
        """The 'Running…' status is set before start(); fake worker resolves
        synchronously so by the time _on_run returns the button is re-enabled
        again — verify the transient text was reachable via a spy worker."""

        class _SpyWorker(_fake_worker_cls(df=_sample_df(1))):
            def start(self):
                # Capture state right before the fake emits.
                _SpyWorker.status_during_run = self._dlg_ref._status_lbl.text()
                super().start()

        # Attach a back-reference so the spy can read dialog state.
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.strike_tool._StrikeWorker",
            _SpyWorker,
        )
        d = StrikeAnalyzerDialog(sites=object())
        _SpyWorker._dlg_ref = d
        d._on_run()
        assert _SpyWorker.status_during_run == "Running…"
        d.close()


# ── _on_done / _on_error direct calls ───────────────────────────────────────


class TestOnDoneDirect:
    def test_on_done_populates_table_and_draws_rose(self, dlg):
        df = _sample_df(4)
        dlg._on_done(df)
        assert dlg._df is df
        assert dlg._table.rowCount() == 4
        assert dlg._run_btn.isEnabled()
        assert "4 stations" in dlg._status_lbl.text()
        assert len(dlg._canvas.figure.axes) == 1

    def test_on_done_empty_df(self, dlg):
        df = _empty_df()
        dlg._on_done(df)
        assert dlg._table.rowCount() == 0
        assert "No strike estimates" in dlg._status_lbl.text()

    def test_on_error_sets_message(self, dlg):
        dlg._on_error("bad things happened")
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: bad things happened"


# ── _fill_table ──────────────────────────────────────────────────────────────


class TestFillTable:
    def test_fill_table_various_row_counts(self, dlg):
        for n in (0, 1, 5):
            dlg._fill_table(_sample_df(n))
            assert dlg._table.rowCount() == n

    def test_fill_table_resets_previous_rows(self, dlg):
        dlg._fill_table(_sample_df(5))
        assert dlg._table.rowCount() == 5
        dlg._fill_table(_sample_df(2))
        assert dlg._table.rowCount() == 2

    def test_fill_table_cell_values(self, dlg):
        df = _sample_df(1)
        dlg._fill_table(df)
        assert dlg._table.item(0, 0).text() == "S0"
        assert dlg._table.item(0, 1).text() == "10.0"
        assert dlg._table.item(0, 2).text() == "2.0"
        assert dlg._table.item(0, 3).text() == "5"

    def test_fill_table_nan_values_show_em_dash(self, dlg):
        df = pd.DataFrame(
            {
                "station": ["S0"],
                "ang": [float("nan")],
                "iqr": [float("nan")],
                "lo": [0.001],
                "hi": [1000.0],
                "n": [0],
            }
        )
        dlg._fill_table(df)
        assert dlg._table.item(0, 1).text() == "—"
        assert dlg._table.item(0, 2).text() == "—"

    def test_fill_table_missing_optional_columns_defaults(self, dlg):
        """station/ang/iqr use .get() with defaults when columns absent."""
        df = pd.DataFrame({"n": [3]})
        dlg._fill_table(df)
        assert dlg._table.rowCount() == 1
        assert dlg._table.item(0, 0).text() == ""
        assert dlg._table.item(0, 1).text() == "—"
        assert dlg._table.item(0, 2).text() == "—"
        assert dlg._table.item(0, 3).text() == "3"


# ── _draw_rose ───────────────────────────────────────────────────────────────


class TestDrawRose:
    def test_draw_rose_with_data_produces_polar_axes(self, dlg):
        df = _sample_df(6)
        dlg._draw_rose(df)
        assert len(dlg._canvas.figure.axes) == 1
        ax = dlg._canvas.figure.axes[0]
        assert "consensus" in ax.get_title().lower()

    def test_draw_rose_all_nan_angles_no_axes_added(self, dlg):
        df = pd.DataFrame(
            {
                "station": ["S0", "S1"],
                "ang": [float("nan"), float("nan")],
                "iqr": [1.0, 1.0],
                "lo": [0.001, 0.001],
                "hi": [1000.0, 1000.0],
                "n": [1, 1],
            }
        )
        dlg._canvas.figure.clear()
        dlg._draw_rose(df)
        assert len(dlg._canvas.figure.axes) == 0

    def test_draw_rose_clears_previous_figure(self, dlg):
        dlg._draw_rose(_sample_df(3))
        assert len(dlg._canvas.figure.axes) == 1
        dlg._draw_rose(_sample_df(5))
        # figure is cleared each call -> still exactly one polar axes
        assert len(dlg._canvas.figure.axes) == 1

    def test_draw_rose_real_worker_output(self, dlg, l18_sites):
        """End-to-end: real estimate_strike_consensus output feeds _draw_rose."""
        from pycsamt.emtools.strike import estimate_strike_consensus

        df = estimate_strike_consensus(
            l18_sites, band=None, w_sweep=0.5, w_pt=0.5, verbose=0
        )
        dlg._fill_table(df)
        dlg._draw_rose(df)
        assert dlg._table.rowCount() == len(df)
        assert len(dlg._canvas.figure.axes) == 1
