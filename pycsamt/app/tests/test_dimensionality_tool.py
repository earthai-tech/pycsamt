# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for DimensionalityDialog / _DimWorker (dimensionality_tool.py).

Real data
---------
data/MT/kap03lmt_edis/ — 26 KP TIPPER EDIs. Loaded once per session and
reused across the worker tests since ``classify_dimensionality`` on this
set completes in well under a second.

Strategy
--------
* ``_DimWorker.run()`` is called directly (never ``.start()``) against the
  real EDI data so the actual ``classify_dimensionality`` computation runs
  and its DataFrame shape/columns are asserted for real.
* ``DimensionalityDialog`` orchestration (_on_run/_on_done/_on_error) is
  tested with a lightweight fake worker (mirroring the ``_FakeSignal``
  idiom used in test_recompute_dlg.py / test_main_window_actions.py) whose
  ``.start()`` synchronously emits.
* ``_fill_table`` / ``_draw_bar`` are also exercised directly against
  synthetic DataFrames to cover column-fallback and empty/missing-column
  branches cheaply.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.app.desktop.tools.dimensionality_tool import (
    DimensionalityDialog,
    _DimWorker,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_TIPPER = _ROOT / "data" / "MT" / "kap03lmt_edis"
_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))


@pytest.fixture(scope="session")
def tipper_sites():
    """26-station TIPPER Sites loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_TIPPER:
        pytest.skip("TIPPER data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_TIPPER))


def _close():
    plt.close("all")


# ── _FakeSignal / fake worker idiom (mirrors test_recompute_dlg.py) ──────────


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls(*, done_result=None, error=None):
    """Build a fake _DimWorker replacement whose .start() emits synchronously."""

    class _FakeWorker:
        captured = []

        def __init__(self, sites, skew_th, ellipt_th, show_map):
            self.sites = sites
            self.skew_th = skew_th
            self.ellipt_th = ellipt_th
            self.show_map = show_map
            _FakeWorker.captured.append(
                (sites, skew_th, ellipt_th, show_map)
            )
            self.done = _FakeSignal()
            self.error = _FakeSignal()

        def start(self):
            if error is not None:
                self.error.emit(error)
            else:
                df, map_fig = done_result
                self.done.emit(df, map_fig)

    return _FakeWorker


# ── Synthetic DataFrames ──────────────────────────────────────────────────────


def _full_df(n_stations=3, per_station=2):
    rows = []
    for si in range(n_stations):
        for pi in range(per_station):
            rows.append(
                dict(
                    station=f"S{si}",
                    freq=1.0 / (pi + 1),
                    period=float(pi + 1),
                    beta_abs=float(pi),
                    ellipt_abs=0.1 * pi,
                    dim=(si + pi) % 3,
                )
            )
    return pd.DataFrame(rows)


# ── _DimWorker: real computation ──────────────────────────────────────────────


class TestDimWorkerConstruction:
    def test_creates(self, qapp):
        w = _DimWorker(None, 3.0, 0.2, False)
        assert w is not None
        assert w._skew_th == 3.0
        assert w._ellipt_th == 0.2
        assert w._show_map is False


class TestDimWorkerRealRun:
    def test_run_emits_done_with_real_dataframe(self, qapp, tipper_sites):
        done_calls = []
        error_calls = []
        w = _DimWorker(tipper_sites, 3.0, 0.2, show_map=False)
        w.done.connect(lambda df, fig: done_calls.append((df, fig)))
        w.error.connect(error_calls.append)

        w.run()

        assert error_calls == []
        assert len(done_calls) == 1
        df, map_fig = done_calls[0]
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        for col in ("station", "period", "beta_abs", "dim"):
            assert col in df.columns
        # show_map was False → no map figure requested
        assert map_fig is None

    def test_run_show_map_true_still_yields_no_map_figure(
        self, qapp, tipper_sites
    ):
        """
        Documents a real bug: the worker imports
        ``pycsamt.emtools.tensor.plot_dim_confidence_grid`` but that
        function actually lives in ``pycsamt.emtools.dimensionality``.
        The import fails, is swallowed by a bare ``except Exception: pass``,
        and map_fig silently stays None even when show_map=True and the
        classification itself succeeds. Not fixed here (tests only).
        """
        done_calls = []
        w = _DimWorker(tipper_sites, 3.0, 0.2, show_map=True)
        w.done.connect(lambda df, fig: done_calls.append((df, fig)))

        w.run()

        assert len(done_calls) == 1
        df, map_fig = done_calls[0]
        assert not df.empty
        assert map_fig is None  # bug: should be a Figure when show_map=True

    def test_run_show_map_true_produces_figure_with_import_patched(
        self, qapp, tipper_sites, monkeypatch
    ):
        """
        Exercises the success path (building `map_fig` from the returned
        Axes) that the import bug above normally makes unreachable, by
        patching the real ``plot_dim_confidence_grid`` (which actually
        lives in ``pycsamt.emtools.dimensionality``) onto
        ``pycsamt.emtools.tensor`` so the worker's buggy import resolves.
        """
        import pycsamt.emtools.tensor as tensor_mod
        from pycsamt.emtools.dimensionality import (
            plot_dim_confidence_grid,
        )

        monkeypatch.setattr(
            tensor_mod,
            "plot_dim_confidence_grid",
            plot_dim_confidence_grid,
            raising=False,
        )
        done_calls = []
        w = _DimWorker(tipper_sites, 3.0, 0.2, show_map=True)
        w.done.connect(lambda df, fig: done_calls.append((df, fig)))
        w.run()
        _close()

        assert len(done_calls) == 1
        df, map_fig = done_calls[0]
        assert not df.empty
        assert map_fig is not None

    def test_run_different_thresholds_change_labels(self, qapp, tipper_sites):
        done_calls = []
        w = _DimWorker(
            tipper_sites, skew_th=1000.0, ellipt_th=1000.0, show_map=False
        )
        w.done.connect(lambda df, fig: done_calls.append((df, fig)))
        w.run()
        df, _ = done_calls[0]
        # Extremely permissive thresholds → everything classifies as 1D (dim=0)
        assert (df["dim"] == 0).all()


class TestDimWorkerErrorPath:
    def test_run_none_sites_emits_error(self, qapp):
        done_calls = []
        error_calls = []
        w = _DimWorker(None, 3.0, 0.2, show_map=False)
        w.done.connect(lambda df, fig: done_calls.append((df, fig)))
        w.error.connect(error_calls.append)

        w.run()

        assert done_calls == []
        assert len(error_calls) == 1
        assert isinstance(error_calls[0], str)
        assert error_calls[0]  # non-empty message

    def test_run_classify_dimensionality_raises_generic_error(
        self, qapp, monkeypatch
    ):
        import pycsamt.emtools.dimensionality as dim_mod

        def _boom(*a, **k):
            raise RuntimeError("synthetic classification failure")

        monkeypatch.setattr(dim_mod, "classify_dimensionality", _boom)

        error_calls = []
        w = _DimWorker(object(), 3.0, 0.2, show_map=False)
        w.error.connect(error_calls.append)
        w.run()

        assert len(error_calls) == 1
        assert "synthetic classification failure" in error_calls[0]

    def test_run_unwraps_frame_attribute(self, qapp, monkeypatch):
        """Cover the `hasattr(df, "frame")` unwrap branch."""
        import pycsamt.emtools.dimensionality as dim_mod

        inner = _full_df()

        class _Wrapped:
            frame = inner

        monkeypatch.setattr(
            dim_mod, "classify_dimensionality", lambda *a, **k: _Wrapped()
        )

        done_calls = []
        w = _DimWorker(object(), 3.0, 0.2, show_map=False)
        w.done.connect(lambda df, fig: done_calls.append((df, fig)))
        w.run()

        assert len(done_calls) == 1
        df, _ = done_calls[0]
        assert df is inner

    def test_run_unwraps_underscore_frame_attribute(self, qapp, monkeypatch):
        """Cover the `hasattr(df, "_frame")` unwrap branch."""
        import pycsamt.emtools.dimensionality as dim_mod

        inner = _full_df()

        class _Wrapped:
            _frame = inner

        monkeypatch.setattr(
            dim_mod, "classify_dimensionality", lambda *a, **k: _Wrapped()
        )

        done_calls = []
        w = _DimWorker(object(), 3.0, 0.2, show_map=False)
        w.done.connect(lambda df, fig: done_calls.append((df, fig)))
        w.run()

        assert len(done_calls) == 1
        df, _ = done_calls[0]
        assert df is inner


# ── DimensionalityDialog construction ─────────────────────────────────────────


class TestDialogConstruction:
    def test_creates_with_sites(self, qapp, tipper_sites):
        dlg = DimensionalityDialog(sites=tipper_sites)
        assert dlg.windowTitle() == "Dimensionality Classifier"
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == ""
        dlg.close()

    def test_creates_without_sites_disables_run(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        assert not dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "No survey loaded."
        dlg.close()

    def test_default_threshold_values(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        assert dlg._skew_spin.value() == pytest.approx(3.0)
        assert dlg._ellipt_spin.value() == pytest.approx(0.20)
        assert dlg._map_cb.isChecked()
        dlg.close()

    def test_table_headers(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        assert dlg._table.columnCount() == 4
        labels = [
            dlg._table.horizontalHeaderItem(i).text() for i in range(4)
        ]
        assert labels == ["Station", "Period (s)", "Skew β (°)", "Dim"]
        dlg.close()

    def test_tabs_present(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        titles = [
            dlg._tabs.tabText(i) for i in range(dlg._tabs.count())
        ]
        assert titles == [
            "Classification Table",
            "Summary Chart",
            "Dimensionality Map",
        ]
        dlg.close()


# ── _on_run / _on_done / _on_error via fake worker ────────────────────────────


class TestOnRun:
    def test_on_run_disables_button_and_constructs_worker(
        self, qapp, monkeypatch
    ):
        df = _full_df()
        fake_cls = _fake_worker_cls(done_result=(df, None))
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.dimensionality_tool._DimWorker",
            fake_cls,
        )
        dlg = DimensionalityDialog(sites="dummy-sites")
        dlg._skew_spin.setValue(5.5)
        dlg._ellipt_spin.setValue(0.33)
        dlg._map_cb.setChecked(False)

        dlg._on_run()

        # start() ran synchronously and emitted done → button re-enabled
        assert dlg._run_btn.isEnabled()
        args = fake_cls.captured[-1]
        assert args == ("dummy-sites", 5.5, 0.33, False)
        dlg.close()

    def test_on_run_error_path_via_fake_worker(self, qapp, monkeypatch):
        fake_cls = _fake_worker_cls(error="boom")
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.dimensionality_tool._DimWorker",
            fake_cls,
        )
        dlg = DimensionalityDialog(sites="dummy-sites")

        dlg._on_run()

        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: boom"
        dlg.close()

    def test_on_run_clears_table_and_sets_status(self, qapp, monkeypatch):
        df = _full_df()
        fake_cls = _fake_worker_cls(done_result=(df, None))
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.dimensionality_tool._DimWorker",
            fake_cls,
        )
        dlg = DimensionalityDialog(sites="dummy-sites")
        # Pre-seed the table to verify _on_run clears it before running.
        dlg._table.setRowCount(2)

        dlg._on_run()

        # Table was repopulated by the (fake) successful run.
        assert dlg._table.rowCount() == len(df)
        dlg.close()


class TestOnDone:
    def test_on_done_none_df_shows_no_results(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        dlg._on_done(None, None)
        assert dlg._status_lbl.text() == "No results."
        assert dlg._run_btn.isEnabled()
        dlg.close()

    def test_on_done_empty_df_shows_no_results(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        dlg._on_done(pd.DataFrame(), None)
        assert dlg._status_lbl.text() == "No results."
        dlg.close()

    def test_on_done_normal_df_fills_table_and_status(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = _full_df(n_stations=3, per_station=2)
        dlg._on_done(df, None)

        assert dlg._table.rowCount() == len(df)
        assert dlg._df is df
        txt = dlg._status_lbl.text()
        assert f"{len(df)} rows" in txt
        assert "1D:" in txt and "2D:" in txt and "3D:" in txt
        dlg.close()
        _close()

    def test_on_done_missing_dim_column_counts_zero(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = pd.DataFrame({"station": ["S0", "S1"], "period": [1.0, 2.0]})
        dlg._on_done(df, None)
        txt = dlg._status_lbl.text()
        assert "1D: 0  2D: 0  3D: 0" in txt
        dlg.close()

    def test_on_done_with_map_fig_updates_map_canvas(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = _full_df()
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])

        dlg._on_done(df, fig)

        assert dlg._map_canvas.figure is fig
        dlg.close()
        _close()

    def test_on_done_without_map_fig_leaves_map_canvas_untouched(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        original_fig = dlg._map_canvas.figure
        df = _full_df()

        dlg._on_done(df, None)

        assert dlg._map_canvas.figure is original_fig
        dlg.close()
        _close()


class TestOnError:
    def test_on_error_sets_message_and_reenables_button(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        dlg._run_btn.setEnabled(False)
        dlg._on_error("something failed")
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: something failed"
        dlg.close()


# ── _fill_table ────────────────────────────────────────────────────────────────


class TestFillTable:
    def test_fill_table_full_columns(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = _full_df(n_stations=2, per_station=2)
        dlg._fill_table(df)

        assert dlg._table.rowCount() == len(df)
        first = df.iloc[0]
        assert dlg._table.item(0, 0).text() == str(first["station"])
        assert dlg._table.item(0, 3).text() in ("1D", "2D", "3D")
        dlg.close()

    def test_fill_table_empty_df_yields_zero_rows(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        empty = pd.DataFrame(
            columns=["station", "period", "beta_abs", "dim"]
        )
        dlg._fill_table(empty)
        assert dlg._table.rowCount() == 0
        dlg.close()

    def test_fill_table_alternate_column_names(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = pd.DataFrame(
            {
                "Station": ["A", "B"],
                "Period": [1.0, 2.0],
                "skew": [1.2, 4.5],
                "dim": [0, 2],
            }
        )
        dlg._fill_table(df)
        assert dlg._table.rowCount() == 2
        assert dlg._table.item(0, 0).text() == "A"
        assert dlg._table.item(0, 3).text() == "1D"
        assert dlg._table.item(1, 3).text() == "3D"
        dlg.close()

    def test_fill_table_id_column_fallback(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = pd.DataFrame({"id": ["Z1"], "T": [1.0], "beta": [0.5], "dim": [1]})
        dlg._fill_table(df)
        assert dlg._table.item(0, 0).text() == "Z1"
        assert dlg._table.item(0, 3).text() == "2D"
        dlg.close()

    def test_fill_table_missing_columns_uses_placeholder(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = pd.DataFrame({"dim": [0, 1, 2]})
        dlg._fill_table(df)
        assert dlg._table.rowCount() == 3
        assert dlg._table.item(0, 0).text() == "—"
        assert dlg._table.item(0, 1).text() == "—"
        assert dlg._table.item(0, 2).text() == "—"
        dlg.close()

    def test_fill_table_no_dim_column_defaults_to_3d(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = pd.DataFrame({"station": ["A"], "period": [1.0]})
        dlg._fill_table(df)
        assert dlg._table.item(0, 3).text() == "3D"
        dlg.close()

    def test_fill_table_clears_previous_rows(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        dlg._fill_table(_full_df(n_stations=3, per_station=2))
        assert dlg._table.rowCount() == 6
        dlg._fill_table(_full_df(n_stations=1, per_station=1))
        assert dlg._table.rowCount() == 1
        dlg.close()


# ── _draw_bar ──────────────────────────────────────────────────────────────────


class TestDrawBar:
    def test_draw_bar_creates_bars_for_each_station(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = _full_df(n_stations=4, per_station=3)

        dlg._draw_bar(df)

        fig = dlg._bar_canvas.figure
        ax = fig.axes[0]
        # 3 bar series (1D/2D/3D) x 4 stations = 12 patches
        assert len(ax.patches) == 3 * 4
        assert ax.get_legend() is not None
        xt_labels = [t.get_text() for t in ax.get_xticklabels()]
        assert xt_labels == ["S0", "S1", "S2", "S3"]
        dlg.close()
        _close()

    def test_draw_bar_missing_dim_column_no_op(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        original_fig = dlg._bar_canvas.figure
        df = pd.DataFrame({"station": ["A", "B"]})

        dlg._draw_bar(df)

        # returned early: bar canvas untouched
        assert dlg._bar_canvas.figure is original_fig
        dlg.close()

    def test_draw_bar_missing_station_column_no_op(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        original_fig = dlg._bar_canvas.figure
        df = pd.DataFrame({"dim": [0, 1, 2]})

        dlg._draw_bar(df)

        assert dlg._bar_canvas.figure is original_fig
        dlg.close()

    def test_draw_bar_single_station(self, qapp):
        dlg = DimensionalityDialog(sites=None)
        df = _full_df(n_stations=1, per_station=5)

        dlg._draw_bar(df)

        fig = dlg._bar_canvas.figure
        ax = fig.axes[0]
        assert len(ax.patches) == 3 * 1
        dlg.close()
        _close()


# ── End-to-end: real data through the fake-worker orchestration path ─────────


class TestEndToEndWithRealDataFakeWorker:
    def test_full_cycle_real_dataframe_through_dialog(
        self, qapp, tipper_sites, monkeypatch
    ):
        """Run the *real* worker computation once, then feed its real
        DataFrame through the dialog's fake-worker-driven _on_run/_on_done
        path to exercise table + chart rendering against real data."""
        real_worker = _DimWorker(tipper_sites, 3.0, 0.2, show_map=False)
        collected = []
        real_worker.done.connect(lambda df, fig: collected.append((df, fig)))
        real_worker.run()
        real_df, real_map_fig = collected[0]

        fake_cls = _fake_worker_cls(done_result=(real_df, real_map_fig))
        monkeypatch.setattr(
            "pycsamt.app.desktop.tools.dimensionality_tool._DimWorker",
            fake_cls,
        )
        dlg = DimensionalityDialog(sites=tipper_sites)

        dlg._on_run()

        assert dlg._table.rowCount() == len(real_df)
        assert f"{len(real_df)} rows" in dlg._status_lbl.text()
        fig = dlg._bar_canvas.figure
        assert len(fig.axes[0].patches) > 0
        dlg.close()
        _close()
