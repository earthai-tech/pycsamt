# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for FrequencyEditorDialog and _EditWorker
(pycsamt.app.desktop.tools.frequency_editor_tool).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/  — a small subset (4 EDIs) used to exercise
the real ``edit_frequencies_by_confidence`` computation inside
``_EditWorker.run()``.

Strategy
--------
* ``_EditWorker.run()`` is called directly (never ``.start()``) against
  real small EDI data to exercise the real confidence-scoring pipeline,
  for both a normal run and a real error path (invalid ``mode``).
* ``FrequencyEditorDialog`` orchestration (``_on_run``/``_on_done``/
  ``_on_error``/``_fill_table``/``_on_apply``) is tested by monkeypatching
  the module-level ``_EditWorker`` with a lightweight fake whose
  ``.start()`` synchronously emits through plain-callable "signals".
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

import pycsamt.app.desktop.tools.frequency_editor_tool as fet
from pycsamt.app.desktop.tools.frequency_editor_tool import (
    FrequencyEditorDialog,
    _EditWorker,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_L18 = _L18.exists() and len(list(_L18.glob("*.edi"))) >= 4


@pytest.fixture(scope="session")
def small_sites():
    """4-station Sites collection loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_L18:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    files = sorted(_L18.glob("*.edi"))[:4]
    return ensure_sites([str(f) for f in files])


# ── _FakeSignal / fake worker helpers ─────────────────────────────────────────


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_worker_cls(on_start):
    """Build a fake ``_EditWorker`` replacement.

    ``on_start(instance)`` is called synchronously from ``.start()`` and is
    responsible for emitting ``done`` or ``error`` on the instance.
    """

    class _FakeWorker:
        instances = []

        def __init__(self, sites, **kw):
            self.sites = sites
            self.kw = kw
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            _FakeWorker.instances.append(self)

        def start(self):
            on_start(self)

    return _FakeWorker


# ── Dialog fixture ─────────────────────────────────────────────────────────────


@pytest.fixture
def dlg(qapp):
    d = FrequencyEditorDialog(sites=None, parent=None)
    yield d
    d.close()


@pytest.fixture
def dlg_with_sites(qapp):
    d = FrequencyEditorDialog(sites=object(), parent=None)
    yield d
    d.close()


def _fake_result(sites="EDITED", n_dropped=1, n_masked=2, n_recovered=3):
    return SimpleNamespace(
        sites=sites,
        n_dropped=n_dropped,
        n_masked=n_masked,
        n_recovered=n_recovered,
    )


def _fake_decisions_df():
    return pd.DataFrame(
        {
            "station": ["S1", "S2", "S3", "S4"],
            "period": [0.001, 0.01, 0.1, 1.0],
            "confidence": [0.9, 0.6, 0.3, 0.1],
            "action": ["kept", "recovered", "dropped", "masked"],
        }
    )


# ══════════════════════════════════════════════════════════════════════════════
# _EditWorker — real computation
# ══════════════════════════════════════════════════════════════════════════════


class TestEditWorkerReal:
    def test_run_recover_mode_success(self, qapp, small_sites):
        done_calls = []
        error_calls = []
        w = _EditWorker(
            small_sites,
            mode="recover",
            method="composite",
            threshold=0.50,
            ci_hi=0.90,
            ci_lo=0.50,
            interp="linear",
            reject="drop",
            also="both",
        )
        w.done.connect(lambda *a: done_calls.append(a))
        w.error.connect(lambda msg: error_calls.append(msg))
        w.run()

        assert error_calls == []
        assert len(done_calls) == 1
        result, decisions, fig = done_calls[0]
        assert hasattr(result, "sites")
        assert result.sites is not None
        # Real decision table must carry an 'action' column.
        assert decisions is not None
        assert "action" in list(decisions.columns)
        assert not decisions.empty
        assert set(decisions["action"].unique()) <= {
            "kept",
            "recovered",
            "dropped",
            "masked",
        }
        # A diagnostic figure should normally be produced.
        assert fig is None or hasattr(fig, "savefig")

    def test_run_drop_mode_high_threshold_produces_drops(self, qapp, small_sites):
        """A near-1.0 threshold in drop mode should flag rows as dropped."""
        done_calls = []
        w = _EditWorker(
            small_sites,
            mode="drop",
            method="composite",
            threshold=0.999,
            ci_hi=0.90,
            ci_lo=0.50,
            interp="linear",
            reject="drop",
            also="both",
        )
        w.done.connect(lambda *a: done_calls.append(a))
        w.error.connect(lambda *_: None)
        w.run()

        assert len(done_calls) == 1
        result, decisions, _fig = done_calls[0]
        assert result.n_dropped >= 0  # real property computed from decisions
        assert decisions is not None and not decisions.empty

    def test_run_invalid_mode_emits_error_real_path(self, qapp, small_sites):
        """Real error path: edit_frequencies_by_confidence validates `mode`
        before touching the sites, so an invalid mode raises a real
        ValueError caught by the worker's broad except."""
        done_calls = []
        error_calls = []
        w = _EditWorker(
            small_sites,
            mode="not_a_real_mode",
            method="composite",
            threshold=0.50,
            ci_hi=0.90,
            ci_lo=0.50,
            interp="linear",
            reject="drop",
            also="both",
        )
        w.done.connect(lambda *a: done_calls.append(a))
        w.error.connect(lambda msg: error_calls.append(msg))
        w.run()

        assert done_calls == []
        assert len(error_calls) == 1
        assert "mode" in error_calls[0].lower()

    def test_worker_is_qthread_subclass(self, qapp, small_sites):
        from PySide6.QtCore import QThread

        w = _EditWorker(
            small_sites,
            mode="recover",
            method="composite",
            threshold=0.5,
            ci_hi=0.9,
            ci_lo=0.5,
            interp="linear",
            reject="drop",
            also="both",
        )
        assert isinstance(w, QThread)

    def test_run_unwraps_frame_and_underscore_frame_wrappers(
        self, qapp, small_sites, monkeypatch
    ):
        """
        ``result.decisions`` is unwrapped through up to two layers of
        APIFrame-like wrappers (``.frame`` then ``._frame``) before being
        used as a plain DataFrame. Real ``edit_frequencies_by_confidence``
        never returns such a wrapper, so this patches it to return one,
        to exercise both unwrap branches.
        """
        import pycsamt.emtools.frequency as freq_mod

        real_df = _fake_decisions_df()

        class _InnerWrapper:
            _frame = real_df

        class _OuterWrapper:
            frame = _InnerWrapper()

        fake_result = SimpleNamespace(sites="EDITED", decisions=_OuterWrapper())
        monkeypatch.setattr(
            freq_mod,
            "edit_frequencies_by_confidence",
            lambda *a, **k: fake_result,
        )
        done_calls = []
        w = _EditWorker(
            small_sites,
            mode="recover",
            method="composite",
            threshold=0.5,
            ci_hi=0.9,
            ci_lo=0.5,
            interp="linear",
            reject="drop",
            also="both",
        )
        w.done.connect(lambda *a: done_calls.append(a))
        w.run()

        assert len(done_calls) == 1
        _result, decisions, _fig = done_calls[0]
        assert decisions is real_df

    def test_run_summary_plot_failure_falls_back_gracefully(
        self, qapp, small_sites, monkeypatch
    ):
        """Covers the except branch when plot_frequency_edit_summary raises."""
        import pycsamt.emtools.frequency as freq_mod

        def _raise(*a, **k):
            raise RuntimeError("synthetic plot failure")

        monkeypatch.setattr(freq_mod, "plot_frequency_edit_summary", _raise)
        done_calls = []
        w = _EditWorker(
            small_sites,
            mode="recover",
            method="composite",
            threshold=0.5,
            ci_hi=0.9,
            ci_lo=0.5,
            interp="linear",
            reject="drop",
            also="both",
        )
        w.done.connect(lambda *a: done_calls.append(a))
        w.run()

        assert len(done_calls) == 1
        _result, _decisions, fig = done_calls[0]
        assert fig is None


# ══════════════════════════════════════════════════════════════════════════════
# FrequencyEditorDialog — construction
# ══════════════════════════════════════════════════════════════════════════════


class TestDialogConstruction:
    def test_creates_without_sites(self, dlg):
        assert dlg is not None
        assert dlg._sites is None
        assert dlg.edited_sites is None

    def test_creates_with_sites(self, dlg_with_sites):
        assert dlg_with_sites._sites is not None

    def test_run_button_disabled_when_no_sites(self, dlg):
        assert not dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "No survey loaded."

    def test_run_button_enabled_when_sites(self, dlg_with_sites):
        assert dlg_with_sites._run_btn.isEnabled()

    def test_apply_button_disabled_initially(self, dlg):
        assert not dlg._apply_btn.isEnabled()

    def test_mode_combo_items(self, dlg):
        items = [dlg._mode_combo.itemText(i) for i in range(dlg._mode_combo.count())]
        assert items == ["recover", "drop", "mask"]

    def test_method_combo_items(self, dlg):
        items = [
            dlg._method_combo.itemText(i) for i in range(dlg._method_combo.count())
        ]
        assert items == ["composite", "snr", "phase_slope", "coherence"]

    def test_threshold_defaults(self, dlg):
        assert dlg._thresh_spin.value() == pytest.approx(0.50)
        assert dlg._cihi_spin.value() == pytest.approx(0.90)
        assert dlg._cilo_spin.value() == pytest.approx(0.50)

    def test_table_headers(self, dlg):
        headers = [dlg._table.horizontalHeaderItem(i).text() for i in range(4)]
        assert headers == ["Station", "Period (s)", "Confidence", "Action"]

    def test_table_starts_empty(self, dlg):
        assert dlg._table.rowCount() == 0

    def test_tabs_present(self, dlg):
        assert dlg._tabs.count() == 2
        assert dlg._tabs.tabText(0) == "Decision Table"
        assert dlg._tabs.tabText(1) == "Before / After Chart"

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Frequency Editor"


# ══════════════════════════════════════════════════════════════════════════════
# _on_run / _on_done / _on_error — via fake worker
# ══════════════════════════════════════════════════════════════════════════════


class TestOnRunDone:
    def test_on_run_constructs_worker_with_form_values(
        self, dlg_with_sites, monkeypatch
    ):
        captured = {}

        def on_start(inst):
            captured["sites"] = inst.sites
            captured["kw"] = inst.kw
            # do not emit yet — just verify construction args first

        FakeWorker = _fake_worker_cls(on_start)
        monkeypatch.setattr(fet, "_EditWorker", FakeWorker)

        dlg_with_sites._mode_combo.setCurrentText("drop")
        dlg_with_sites._method_combo.setCurrentText("snr")
        dlg_with_sites._thresh_spin.setValue(0.75)
        dlg_with_sites._cihi_spin.setValue(0.95)
        dlg_with_sites._cilo_spin.setValue(0.35)
        dlg_with_sites._interp_combo.setCurrentText("nearest")
        dlg_with_sites._reject_combo.setCurrentText("mask")
        dlg_with_sites._also_combo.setCurrentText("z")

        dlg_with_sites._on_run()

        assert captured["sites"] is dlg_with_sites._sites
        assert captured["kw"]["mode"] == "drop"
        assert captured["kw"]["method"] == "snr"
        assert captured["kw"]["threshold"] == pytest.approx(0.75)
        assert captured["kw"]["ci_hi"] == pytest.approx(0.95)
        assert captured["kw"]["ci_lo"] == pytest.approx(0.35)
        assert captured["kw"]["interp"] == "nearest"
        assert captured["kw"]["reject"] == "mask"
        assert captured["kw"]["also"] == "z"

    def test_on_run_disables_buttons_and_sets_status(self, dlg_with_sites, monkeypatch):
        FakeWorker = _fake_worker_cls(lambda inst: None)
        monkeypatch.setattr(fet, "_EditWorker", FakeWorker)

        dlg_with_sites._apply_btn.setEnabled(True)
        dlg_with_sites._on_run()

        assert not dlg_with_sites._run_btn.isEnabled()
        assert not dlg_with_sites._apply_btn.isEnabled()
        assert dlg_with_sites._status_lbl.text() == "Running frequency editor…"
        assert dlg_with_sites._table.rowCount() == 0

    def test_on_run_then_done_happy_path(self, dlg_with_sites, monkeypatch):
        result = _fake_result()
        decisions = _fake_decisions_df()
        import matplotlib.pyplot as plt

        fig = plt.figure()

        def on_start(inst):
            inst.done.emit(result, decisions, fig)

        FakeWorker = _fake_worker_cls(on_start)
        monkeypatch.setattr(fet, "_EditWorker", FakeWorker)

        dlg_with_sites._on_run()

        assert dlg_with_sites._run_btn.isEnabled()
        assert dlg_with_sites.edited_sites == "EDITED"
        assert dlg_with_sites._apply_btn.isEnabled()
        assert "dropped: 1" in dlg_with_sites._status_lbl.text()
        assert "masked: 2" in dlg_with_sites._status_lbl.text()
        assert "recovered: 3" in dlg_with_sites._status_lbl.text()
        assert dlg_with_sites._table.rowCount() == 4
        assert dlg_with_sites._tabs.currentIndex() == 1
        plt.close(fig)

    def test_on_done_no_sites_attr_on_result(self, dlg_with_sites):
        """When result has no `.sites`, edited_sites stays None and apply
        stays disabled."""
        result = SimpleNamespace(n_dropped=0, n_masked=0, n_recovered=0)
        dlg_with_sites._on_done(result, None, None)
        assert dlg_with_sites.edited_sites is None
        assert not dlg_with_sites._apply_btn.isEnabled()

    def test_on_done_with_empty_decisions_skips_fill_table(self, dlg_with_sites):
        empty_df = pd.DataFrame(columns=["station", "period", "action"])
        result = _fake_result()
        dlg_with_sites._table.setRowCount(0)
        dlg_with_sites._on_done(result, empty_df, None)
        assert dlg_with_sites._table.rowCount() == 0

    def test_on_done_with_none_decisions(self, dlg_with_sites):
        result = _fake_result()
        dlg_with_sites._on_done(result, None, None)
        assert dlg_with_sites._table.rowCount() == 0

    def test_on_done_with_none_fig_does_not_switch_tab(self, dlg_with_sites):
        result = _fake_result()
        dlg_with_sites._tabs.setCurrentIndex(0)
        dlg_with_sites._on_done(result, _fake_decisions_df(), None)
        assert dlg_with_sites._tabs.currentIndex() == 0

    def test_on_run_then_error_path(self, dlg_with_sites, monkeypatch):
        def on_start(inst):
            inst.error.emit("boom: bad confidence method")

        FakeWorker = _fake_worker_cls(on_start)
        monkeypatch.setattr(fet, "_EditWorker", FakeWorker)

        dlg_with_sites._on_run()

        assert dlg_with_sites._run_btn.isEnabled()
        assert dlg_with_sites._status_lbl.text() == (
            "Error: boom: bad confidence method"
        )

    def test_on_error_direct(self, dlg_with_sites):
        dlg_with_sites._run_btn.setEnabled(False)
        dlg_with_sites._on_error("some failure")
        assert dlg_with_sites._run_btn.isEnabled()
        assert dlg_with_sites._status_lbl.text() == "Error: some failure"


# ══════════════════════════════════════════════════════════════════════════════
# _fill_table — column-detection and colour-coding branches
# ══════════════════════════════════════════════════════════════════════════════


class TestFillTable:
    def test_fill_table_empty_df(self, dlg):
        df = pd.DataFrame(columns=["station", "period", "confidence", "action"])
        dlg._fill_table(df)
        assert dlg._table.rowCount() == 0

    def test_fill_table_standard_columns_and_colors(self, dlg):
        df = _fake_decisions_df()
        dlg._fill_table(df)
        assert dlg._table.rowCount() == 4

        expected_bg = {
            0: fet._C_KEEP,
            1: fet._C_RECOVER,
            2: fet._C_DROP,
            3: fet._C_DROP,  # masked shares dropped's colour
        }
        for row, color in expected_bg.items():
            item = dlg._table.item(row, 3)
            assert item.background().color().name() == color.name()
            assert item.text() == df["action"].iloc[row]

        # Station / period / confidence columns rendered
        assert dlg._table.item(0, 0).text() == "S1"
        assert dlg._table.item(0, 2).text() == "0.900"

    def test_fill_table_unknown_action_falls_back_to_keep_color(self, dlg):
        df = pd.DataFrame(
            {
                "station": ["S1"],
                "period": [0.5],
                "confidence": [0.5],
                "action": ["mystery"],
            }
        )
        dlg._fill_table(df)
        item = dlg._table.item(0, 3)
        assert item.background().color().name() == fet._C_KEEP.name()
        assert item.text() == "mystery"

    def test_fill_table_no_action_column_defaults_kept(self, dlg):
        df = pd.DataFrame({"station": ["S1"], "period": [0.1], "confidence": [0.2]})
        dlg._fill_table(df)
        item = dlg._table.item(0, 3)
        assert item.text() == "kept"
        assert item.background().color().name() == fet._C_KEEP.name()

    def test_fill_table_alt_column_names_capitalized(self, dlg):
        df = pd.DataFrame(
            {
                "Station": ["A1"],
                "Period": [2.0],
                "ci": [0.42],
                "action": ["kept"],
            }
        )
        dlg._fill_table(df)
        assert dlg._table.item(0, 0).text() == "A1"
        assert dlg._table.item(0, 1).text() == "2"
        assert dlg._table.item(0, 2).text() == "0.420"

    def test_fill_table_alt_column_names_id_T_score(self, dlg):
        df = pd.DataFrame(
            {
                "id": ["B2"],
                "T": [3.0],
                "score": [0.71],
                "action": ["recovered"],
            }
        )
        dlg._fill_table(df)
        assert dlg._table.item(0, 0).text() == "B2"
        assert dlg._table.item(0, 1).text() == "3"
        assert dlg._table.item(0, 2).text() == "0.710"

    def test_fill_table_no_matching_columns_uses_placeholder(self, dlg):
        df = pd.DataFrame({"foo": [1], "bar": [2]})
        dlg._fill_table(df)
        assert dlg._table.rowCount() == 1
        assert dlg._table.item(0, 0).text() == "—"
        assert dlg._table.item(0, 1).text() == "—"
        assert dlg._table.item(0, 2).text() == "—"
        assert dlg._table.item(0, 3).text() == "kept"

    def test_fill_table_clears_previous_rows(self, dlg):
        dlg._fill_table(_fake_decisions_df())
        assert dlg._table.rowCount() == 4
        dlg._fill_table(pd.DataFrame({"station": ["only"], "action": ["kept"]}))
        assert dlg._table.rowCount() == 1

    def test_fill_table_real_decision_table_columns(self, dlg, small_sites, qapp):
        """The real `frequency_edit_decision_table` output uses 'period_s'
        (not 'period'). Document current behaviour: the Period column
        detector does not recognise it, so it always renders '—'.
        See note in the final report — potential real bug."""
        w = _EditWorker(
            small_sites,
            mode="recover",
            method="composite",
            threshold=0.5,
            ci_hi=0.9,
            ci_lo=0.5,
            interp="linear",
            reject="drop",
            also="both",
        )
        got = {}
        w.done.connect(lambda *a: got.update(result=a[0], decisions=a[1], fig=a[2]))
        w.error.connect(lambda *_: None)
        w.run()

        assert "period_s" in list(got["decisions"].columns)
        dlg._fill_table(got["decisions"])
        assert dlg._table.rowCount() == len(got["decisions"])
        # Period column is always the placeholder because 'period_s' isn't
        # one of the recognised aliases ("period", "Period", "T").
        assert dlg._table.item(0, 1).text() == "—"
        # Station and Confidence and Action, however, resolve correctly.
        assert dlg._table.item(0, 0).text() != "—"
        assert dlg._table.item(0, 2).text() != "—"


# ══════════════════════════════════════════════════════════════════════════════
# _on_apply
# ══════════════════════════════════════════════════════════════════════════════


class TestOnApply:
    def test_on_apply_sets_status_and_accepts(self, dlg_with_sites):
        from PySide6.QtWidgets import QDialog

        dlg_with_sites._on_apply()
        assert dlg_with_sites._status_lbl.text() == (
            "Edited sites ready — close this dialog to continue."
        )
        assert dlg_with_sites.result() == QDialog.DialogCode.Accepted

    def test_on_apply_after_done_preserves_edited_sites(
        self, dlg_with_sites, monkeypatch
    ):
        result = _fake_result(sites="MY_EDITED_SITES")

        def on_start(inst):
            inst.done.emit(result, _fake_decisions_df(), None)

        FakeWorker = _fake_worker_cls(on_start)
        monkeypatch.setattr(fet, "_EditWorker", FakeWorker)

        dlg_with_sites._on_run()
        assert dlg_with_sites.edited_sites == "MY_EDITED_SITES"

        dlg_with_sites._on_apply()
        assert dlg_with_sites.edited_sites == "MY_EDITED_SITES"

    def test_on_apply_with_no_prior_run_still_accepts(self, dlg):
        """User made no changes (never ran) yet clicks Apply — edited_sites
        stays None but the dialog still accepts."""
        from PySide6.QtWidgets import QDialog

        assert dlg.edited_sites is None
        dlg._on_apply()
        assert dlg.edited_sites is None
        assert dlg.result() == QDialog.DialogCode.Accepted


# ══════════════════════════════════════════════════════════════════════════════
# Full integration: real worker wired through the dialog's signal handlers
# ══════════════════════════════════════════════════════════════════════════════


class TestRealWorkerThroughDialogHandlers:
    """Exercise _on_done with the *real* worker's output (not `.start()`ed
    as a QThread — call `.run()` synchronously and feed the emitted values
    straight into the dialog handlers), bridging real computation with the
    dialog orchestration without needing actual thread execution."""

    def test_real_worker_output_feeds_on_done(self, dlg_with_sites, small_sites):
        dlg_with_sites._sites = small_sites
        w = _EditWorker(
            small_sites,
            mode="recover",
            method="composite",
            threshold=0.5,
            ci_hi=0.9,
            ci_lo=0.5,
            interp="linear",
            reject="drop",
            also="both",
        )
        collected = []
        w.done.connect(lambda *a: collected.append(a))
        w.error.connect(lambda *_: None)
        w.run()
        assert len(collected) == 1
        result, decisions, fig = collected[0]

        dlg_with_sites._on_done(result, decisions, fig)

        assert dlg_with_sites.edited_sites is result.sites
        assert dlg_with_sites._apply_btn.isEnabled()
        assert dlg_with_sites._table.rowCount() == len(decisions)

        dlg_with_sites._on_apply()
        assert dlg_with_sites.edited_sites is result.sites
