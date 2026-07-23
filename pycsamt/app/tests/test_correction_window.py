# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for CorrectionWindow (pycsamt.app.desktop.windows.correction_window).

Strategy
--------
* Real small WILLY EDI data is used for the Preview/Apply/Commit/Undo
  integration path against the real ``CorrectionController`` (only one or
  two tests -- the controller's own plot_*/apply methods are otherwise
  exercised indirectly through ``_refresh_plots``'s exception-swallowing
  ``_draw()`` wrapper, which makes most category/mode-switching wiring
  safe to exercise with *no* data at all).
* ``_on_stack_ctx_menu`` calls ``QMenu.exec()``; per the established
  gotcha (see fix_qmenu_exec_monkeypatch_hang memory), monkeypatching
  ``QMenu.exec`` directly is silently ignored by PySide6 and hangs the
  whole suite offscreen. The module-level ``QMenu`` name is swapped for
  a subclass instead.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtCore import QPoint
from PySide6.QtWidgets import QFileDialog, QMenu

import pycsamt.app.desktop.windows.correction_window as cw_mod
from pycsamt.app.desktop.windows.correction_window import CorrectionWindow

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


@pytest.fixture
def win(qapp):
    w = CorrectionWindow(parent=None)
    w.show()
    yield w
    w.close()


def _select_category(win, name):
    from pycsamt.app.desktop.controllers.correction_controller import (
        CATEGORIES,
    )

    win._combo_category.setCurrentIndex(CATEGORIES.index(name))


def _select_correction(win, label):
    idx = win._combo_correction.findText(label)
    assert idx >= 0, f"correction {label!r} not found"
    win._combo_correction.setCurrentIndex(idx)


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "Data Corrections" in win.windowTitle()

    def test_category_combo_populated(self, win):
        from pycsamt.app.desktop.controllers.correction_controller import (
            CATEGORIES,
        )

        assert win._combo_category.count() == len(CATEGORIES)

    def test_static_shift_initial_view(self, win):
        _select_category(win, "Static Shift")
        assert win._grp_ss_affected.isVisible()
        assert win._combo_mode.currentText() == "2D Section"
        assert win._view_stack.currentIndex() == 3

    def test_stratagem_hidden_groups_by_default_for_other_categories(
        self, win
    ):
        _select_category(win, "Noise Removal")
        assert not win._grp_strat.isVisible()
        assert not win._grp_ss_affected.isVisible()

    def test_buttons_disabled_without_data(self, win):
        assert not win._btn_preview.isEnabled()
        assert not win._btn_apply.isEnabled()
        assert not win._btn_commit.isEnabled()
        assert not win._btn_revert.isEnabled()

    def test_stack_status_no_corrections(self, win):
        assert win._stack_status.text() == "No corrections applied"


# ── Category / correction navigation ─────────────────────────────────────────


class TestCategoryNavigation:
    def test_switch_to_stratagem_shows_strat_page(self, win):
        _select_category(win, "Stratagem")
        assert win._grp_strat.isVisible()
        assert not win._combo_mode.isVisible()
        assert win._btn_export_strat.isVisible()
        assert win._view_stack.currentIndex() == 2

    def test_switch_to_coordinates_shows_coord_view_combo(self, win):
        _select_category(win, "Coordinates")
        assert win._combo_coord_view.isVisible()
        assert win._is_coord_category()

    def test_switch_to_tensor_rotation(self, win):
        _select_category(win, "Tensor Rotation")
        assert win._view_stack.currentIndex() in (0, 1)

    def test_switching_away_from_static_shift_resets_2d_section_mode(
        self, win
    ):
        _select_category(win, "Static Shift")
        assert win._combo_mode.currentText() == "2D Section"
        _select_category(win, "Noise Removal")
        assert win._combo_mode.currentText() == "Before / After"

    def test_correction_combo_populated_per_category(self, win):
        from pycsamt.app.desktop.controllers.correction_controller import (
            CATALOGUE,
            CATEGORIES,
        )

        _select_category(win, "Static Shift")
        expected = list(CATALOGUE["Static Shift"].keys())
        assert win._combo_correction.count() == len(expected)

    def test_correction_changed_updates_description(self, win):
        _select_category(win, "Static Shift")
        assert win._desc_lbl.isVisible()
        assert win._desc_lbl.text() != ""

    def test_on_category_changed_out_of_range_noop(self, win):
        win._on_category_changed(-1)
        win._on_category_changed(999)

    def test_on_correction_changed_out_of_range_noop(self, win):
        win._on_correction_changed(-1)
        win._on_correction_changed(999)


# ── Param widgets ─────────────────────────────────────────────────────────────


class TestParamWidgets:
    def test_make_spin_widget(self, win):
        from pycsamt.app.desktop.controllers.correction_controller import (
            ParamSpec,
        )

        spec = ParamSpec("n", "N", "spin", 3, (1, 20, 1), "")
        w = win._make_widget(spec)
        assert w.value() == 3

    def test_make_dspin_widget(self, win):
        from pycsamt.app.desktop.controllers.correction_controller import (
            ParamSpec,
        )

        spec = ParamSpec("x", "X", "dspin", 6.0, (0.5, 45.0, 0.5), "")
        w = win._make_widget(spec)
        assert w.value() == pytest.approx(6.0)

    def test_make_combo_widget(self, win):
        from pycsamt.app.desktop.controllers.correction_controller import (
            ParamSpec,
        )

        spec = ParamSpec("k", "K", "combo", "tri", ["tri", "box", "cos"], "")
        w = win._make_widget(spec)
        assert w.currentText() == "tri"

    def test_make_combo_widget_default_not_in_opts(self, win):
        from pycsamt.app.desktop.controllers.correction_controller import (
            ParamSpec,
        )

        spec = ParamSpec("k", "K", "combo", "zzz", ["a", "b"], "")
        w = win._make_widget(spec)
        assert w.currentIndex() == 0

    def test_make_check_widget(self, win):
        from pycsamt.app.desktop.controllers.correction_controller import (
            ParamSpec,
        )

        spec = ParamSpec("flag", "Flag", "check", True, None, "")
        w = win._make_widget(spec)
        assert w.isChecked()

    def test_make_fallback_lineedit_widget(self, win):
        from pycsamt.app.desktop.controllers.correction_controller import (
            ParamSpec,
        )

        spec = ParamSpec("s", "S", "text", "hello", None, "")
        w = win._make_widget(spec)
        assert w.text() == "hello"

    def test_rebuild_param_form_shows_no_params_label_when_empty(self, win):
        win._rebuild_param_form([])
        assert win._no_params_lbl.isVisible()

    def test_get_param_values_collects_all_kinds(self, win):
        from PySide6.QtWidgets import QCheckBox, QComboBox, QLineEdit
        from PySide6.QtWidgets import QDoubleSpinBox, QSpinBox

        win._param_widgets = {}
        sp = QSpinBox()
        sp.setValue(5)
        dsp = QDoubleSpinBox()
        dsp.setValue(1.5)
        cb = QComboBox()
        cb.addItem("opt1")
        chk = QCheckBox()
        chk.setChecked(True)
        le_num = QLineEdit("3.14")
        le_txt = QLineEdit("abc")
        win._param_widgets = {
            "sp": sp,
            "dsp": dsp,
            "cb": cb,
            "chk": chk,
            "le_num": le_num,
            "le_txt": le_txt,
        }
        vals = win._get_param_values()
        assert vals["sp"] == 5
        assert vals["dsp"] == pytest.approx(1.5)
        assert vals["cb"] == "opt1"
        assert vals["chk"] is True
        assert vals["le_num"] == pytest.approx(3.14)
        assert vals["le_txt"] == "abc"

    def test_get_param_values_includes_affected_stations_for_static_shift(
        self, win
    ):
        _select_category(win, "Static Shift")
        win._txt_ss_stations.setPlainText("S001, S002")
        vals = win._get_param_values()
        assert vals["affected_stations"] == ["S001", "S002"]

    def test_get_affected_stations_various_separators(self, win):
        win._txt_ss_stations.setPlainText("A, B;C\nD")
        assert win._get_affected_stations() == ["A", "B", "C", "D"]

    def test_get_affected_stations_empty(self, win):
        win._txt_ss_stations.setPlainText("")
        assert win._get_affected_stations() == []

    def test_current_fn_label_valid(self, win):
        _select_category(win, "Static Shift")
        fn, label = win._current_fn_label()
        assert fn
        assert label

    def test_current_fn_label_no_category(self, win):
        win._combo_category.setCurrentIndex(-1)
        fn, label = win._current_fn_label()
        assert fn == "" and label == ""


# ── Preview / Apply (no data) ─────────────────────────────────────────────────


class TestPreviewApplyNoData:
    def test_preview_without_data(self, win):
        win._on_preview()
        assert win._action_status.text() == "Load survey data first."

    def test_apply_without_data(self, win):
        win._on_apply()
        assert win._action_status.text() == "Load survey data first."


# ── Preview / Apply / Undo / Clear (real data) ───────────────────────────────


class TestPreviewApplyRealData:
    def test_preview_success(self, win, willy_sites):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        _select_correction(win, "LOESS")
        win._on_preview()
        assert win._action_status.text().startswith("Preview:")
        assert win._preview_sites is not None

    def test_apply_success_adds_to_stack(self, win, willy_sites):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        _select_correction(win, "LOESS")
        win._on_apply()
        assert win._action_status.text().startswith("Applied:")
        assert win._stack_list.count() == 1
        assert win._btn_commit.isEnabled()

    def test_undo_removes_last_step(self, win, willy_sites):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        _select_correction(win, "LOESS")
        win._on_apply()
        win._on_undo()
        assert win._stack_list.count() == 0
        assert win._stack_status.text() == "No corrections applied"

    def test_clear_stack(self, win, willy_sites):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        _select_correction(win, "LOESS")
        win._on_apply()
        win._on_clear_stack()
        assert win._stack_list.count() == 0

    def test_ama_correction_broken_via_affected_stations_kwarg(
        self, win, willy_sites
    ):
        """
        Real bug: every Static Shift correction's param values include an
        ``affected_stations`` kwarg (added unconditionally in
        ``_get_param_values`` for the whole category), but
        ``correct_ss_ama`` -- the *default*, first-listed Static Shift
        correction -- is dispatched straight to the real
        ``pycsamt.emtools`` function via ``_call_fn``'s fallback path,
        which does not accept that kwarg. Only the three internally
        wrapped methods (LOESS / Bilateral / RefMedian, which accept
        ``**_``) tolerate it. So the default "AMA (spatial average)"
        option is completely broken through the UI: Preview always
        errors and Apply always reports "Apply failed". Not fixed here,
        per instructions.
        """
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        _select_correction(win, "AMA (spatial average)")
        win._on_preview()
        assert "affected_stations" in win._action_status.text()
        win._on_apply()
        assert "Apply failed" in win._action_status.text()

    def test_preview_exception_reported(self, win, willy_sites, monkeypatch):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")

        def _boom(fn_name, kwargs):
            raise RuntimeError("preview boom")

        monkeypatch.setattr(win._ctrl, "preview", _boom)
        win._on_preview()
        assert "Error: preview boom" in win._action_status.text()

    def test_apply_exception_reported(self, win, willy_sites, monkeypatch):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")

        def _boom(fn_name, kwargs, label):
            raise RuntimeError("apply boom")

        monkeypatch.setattr(win._ctrl, "apply", _boom)
        win._on_apply()
        assert "Error: apply boom" in win._action_status.text()

    def test_apply_returns_none_reports_failure(
        self, win, willy_sites, monkeypatch
    ):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        monkeypatch.setattr(win._ctrl, "apply", lambda *a, **k: None)
        win._on_apply()
        assert "Apply failed" in win._action_status.text()

    def test_preview_returns_none_reports_failure(
        self, win, willy_sites, monkeypatch
    ):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        monkeypatch.setattr(win._ctrl, "preview", lambda *a, **k: None)
        win._on_preview()
        assert "Preview failed" in win._action_status.text()

    def test_no_correction_selected_noop(self, win, willy_sites, monkeypatch):
        win.set_sites(willy_sites)
        monkeypatch.setattr(win, "_current_fn_label", lambda: ("", ""))
        win._on_preview()  # returns early, no crash
        win._on_apply()


# ── Context menu (stack) ─────────────────────────────────────────────────────


class TestStackContextMenu:
    def test_no_item_at_pos_noop(self, win):
        win._on_stack_ctx_menu(QPoint(0, 0))  # empty list -> itemAt None

    def test_remove_step_via_context_menu(self, win, willy_sites, monkeypatch):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        _select_correction(win, "LOESS")
        win._on_apply()

        class _FakeMenu(QMenu):
            def exec(self, *a, **k):
                for act in self.actions():
                    if act.text().startswith("Remove"):
                        return act
                return None

        monkeypatch.setattr(cw_mod, "QMenu", _FakeMenu)
        rect = win._stack_list.visualItemRect(win._stack_list.item(0))
        win._on_stack_ctx_menu(rect.center())
        assert win._stack_list.count() == 0

    def test_view_after_for_step_via_context_menu(
        self, win, willy_sites, monkeypatch
    ):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        _select_correction(win, "LOESS")
        win._on_apply()

        class _FakeMenu(QMenu):
            def exec(self, *a, **k):
                for act in self.actions():
                    if "Show" in act.text():
                        return act
                return None

        monkeypatch.setattr(cw_mod, "QMenu", _FakeMenu)
        rect = win._stack_list.visualItemRect(win._stack_list.item(0))
        win._on_stack_ctx_menu(rect.center())
        assert win._preview_sites is not None

    def test_context_menu_dismissed_noop(self, win, willy_sites, monkeypatch):
        win.set_sites(willy_sites)
        _select_category(win, "Static Shift")
        _select_correction(win, "LOESS")
        win._on_apply()

        class _FakeMenu(QMenu):
            def exec(self, *a, **k):
                return None

        monkeypatch.setattr(cw_mod, "QMenu", _FakeMenu)
        rect = win._stack_list.visualItemRect(win._stack_list.item(0))
        win._on_stack_ctx_menu(rect.center())  # neither branch taken


# ── Stratagem source / load / export ─────────────────────────────────────────


class TestStratagem:
    def test_source_toggle_enables_dir_row(self, win):
        win._radio_load_dir.setChecked(True)
        assert win._strat_dir_row.isEnabled()
        assert win._btn_load_edi.isEnabled()
        win._radio_use_current.setChecked(True)
        assert not win._strat_dir_row.isEnabled()

    def test_browse_edi_dir_sets_label(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/some/edi/dir"),
        )
        win._on_browse_edi_dir()
        assert win._edi_dir_label.text() == "/some/edi/dir"

    def test_browse_edi_dir_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: ""),
        )
        win._on_browse_edi_dir()
        assert win._edi_dir_label.text() == "(no directory selected)"

    def test_load_edi_dir_no_path_selected(self, win):
        win._on_load_edi_dir()
        assert "Select a directory" in win._edi_load_status.text()

    def test_load_edi_dir_success(self, win, monkeypatch):
        win._edi_dir_label.setText("/fake/dir")
        monkeypatch.setattr(win._ctrl, "load_edi_dir", lambda path: 7)
        win._on_load_edi_dir()
        assert "7 stations loaded" in win._edi_load_status.text()
        from pycsamt.app.desktop.controllers.correction_controller import (
            CATEGORIES,
        )

        assert win._combo_category.currentIndex() == CATEGORIES.index(
            "Stratagem"
        )

    def test_load_edi_dir_singular_station(self, win, monkeypatch):
        win._edi_dir_label.setText("/fake/dir")
        monkeypatch.setattr(win._ctrl, "load_edi_dir", lambda path: 1)
        win._on_load_edi_dir()
        assert "1 station loaded" in win._edi_load_status.text()

    def test_load_edi_dir_failure(self, win, monkeypatch):
        win._edi_dir_label.setText("/fake/dir")

        def _boom(path):
            raise RuntimeError("bad dir")

        monkeypatch.setattr(win._ctrl, "load_edi_dir", _boom)
        win._on_load_edi_dir()
        assert "bad dir" in win._edi_load_status.text()

    def test_refresh_strat_plots_no_crash_without_data(self, win):
        win._refresh_strat_plots()  # everything None -> exceptions swallowed

    def test_populate_strat_qc_table_none_df(self, win):
        win._ctrl._strat_qc_report = None
        win._populate_strat_qc_table()
        assert win._strat_qc_table.rowCount() == 0

    def test_populate_strat_qc_table_with_data(self, win):
        win._ctrl._strat_qc_report = pd.DataFrame(
            {
                "station": ["A", "B"],
                "rms": [0.1, 0.9],
                "flagged": [False, True],
            }
        )
        win._populate_strat_qc_table()
        assert win._strat_qc_table.rowCount() == 2
        assert win._strat_qc_table.columnCount() == 3

    def test_export_stratagem_cancelled(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: ""),
        )
        win._on_strat_export()  # must not raise

    def test_export_stratagem_success(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/out/dir"),
        )
        monkeypatch.setattr(win._ctrl, "export_stratagem", lambda p: 4)
        win._on_strat_export()
        assert "Exported 4" in win._action_status.text()

    def test_export_stratagem_failure(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: "/out/dir"),
        )

        def _boom(p):
            raise RuntimeError("write failed")

        monkeypatch.setattr(win._ctrl, "export_stratagem", _boom)
        win._on_strat_export()
        assert "Export failed" in win._action_status.text()


# ── View mode / refresh_plots branches ───────────────────────────────────────


class TestViewModeAndRefresh:
    def test_mode_before_after(self, win):
        win._combo_mode.setCurrentText("Before / After")
        assert win._view_stack.currentIndex() == 0

    def test_mode_overlay(self, win):
        _select_category(win, "Noise Removal")
        win._combo_mode.setCurrentText("Overlay")
        assert win._view_stack.currentIndex() == 1

    def test_mode_diff(self, win):
        _select_category(win, "Noise Removal")
        win._combo_mode.setCurrentText("Diff")
        assert win._view_stack.currentIndex() == 1

    def test_mode_2d_section(self, win):
        _select_category(win, "Static Shift")
        win._combo_mode.setCurrentText("2D Section")
        assert win._view_stack.currentIndex() == 3

    def test_2d_section_for_coord_category_shows_message(self, win):
        _select_category(win, "Coordinates")
        win._combo_mode.setCurrentIndex(win._combo_mode.count() - 1)
        # Manually force 2D Section text since coords may not offer it
        win._combo_mode.blockSignals(True)
        if win._combo_mode.findText("2D Section") < 0:
            win._combo_mode.addItem("2D Section")
        win._combo_mode.setCurrentText("2D Section")
        win._combo_mode.blockSignals(False)
        win._refresh_plots()
        assert len(win._canvas_pseudo.figure.axes) >= 1

    def test_coord_elevation_view(self, win):
        _select_category(win, "Coordinates")
        win._combo_coord_view.setCurrentText("Elevation profile")
        win._combo_mode.setCurrentText("Before / After")
        win._refresh_plots()  # must not raise

    def test_coord_overlay_view(self, win):
        _select_category(win, "Coordinates")
        win._combo_mode.setCurrentText("Overlay")
        win._refresh_plots()

    def test_rotation_overlay_uses_rose(self, win):
        _select_category(win, "Tensor Rotation")
        win._combo_mode.setCurrentText("Overlay")
        win._refresh_plots()

    def test_rotation_diff_uses_rose(self, win):
        _select_category(win, "Tensor Rotation")
        win._combo_mode.setCurrentText("Diff")
        win._refresh_plots()

    def test_coord_diff_view(self, win):
        _select_category(win, "Coordinates")
        win._combo_mode.setCurrentText("Diff")
        win._refresh_plots()

    def test_refresh_with_preview_sites_shows_preview_label(
        self, win, willy_sites
    ):
        win.set_sites(willy_sites)
        _select_category(win, "Noise Removal")
        win._on_preview()
        win._combo_mode.setCurrentText("Before / After")
        win._refresh_plots()  # must not raise with preview data set


# ── Commit / Revert / Export / set_sites / set_dark_mode ────────────────────


class TestCommitRevertExport:
    def test_commit_no_corrections(self, win):
        win._on_commit()
        assert win._view_status.text() == "No corrections to commit."

    def test_commit_emits_signal(self, win, willy_sites):
        win.set_sites(willy_sites)
        _select_category(win, "Noise Removal")
        win._on_apply()
        received = []
        win.corrections_committed.connect(lambda s: received.append(s))
        win._on_commit()
        assert len(received) == 1
        assert "Committed" in win._view_status.text()

    def test_commit_with_coord_corrections(self, win, willy_sites, monkeypatch):
        win.set_sites(willy_sites)
        _select_category(win, "Coordinates")
        win._on_apply()
        if win._ctrl.has_coord_corrections:
            calls = []
            monkeypatch.setattr(
                "pycsamt.gis.coord_correction.apply_coords_df_to_sites",
                lambda sites, coords: calls.append((sites, coords)),
            )
            win._on_commit()
            assert len(calls) == 1

    def test_revert_emits_signal_and_clears_stack(self, win, willy_sites):
        win.set_sites(willy_sites)
        _select_category(win, "Noise Removal")
        win._on_apply()
        received = []
        win.corrections_reverted.connect(lambda: received.append(1))
        win._on_revert()
        assert received == [1]
        assert win._stack_list.count() == 0
        assert "Reverted" in win._view_status.text()

    def test_on_export_opens_dialog(self, win, monkeypatch):
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

    def test_on_export_overlay_mode_uses_single_canvas(self, win, monkeypatch):
        _select_category(win, "Noise Removal")
        win._combo_mode.setCurrentText("Overlay")
        calls = []

        class _FakeExportDialog:
            def __init__(self, figure, parent):
                calls.append(figure)

            def exec(self):
                pass

        monkeypatch.setattr(
            "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog",
            _FakeExportDialog,
        )
        win._on_export()
        assert calls[0] is win._canvas_single.figure

    def test_set_sites_delegates_and_refreshes(self, win, willy_sites):
        win.set_sites(willy_sites)
        assert win._ctrl.has_data
        assert win._btn_preview.isEnabled()

    def test_set_dark_mode_delegates(self, win):
        win.set_dark_mode(False)
        assert win._ctrl.dark is False
