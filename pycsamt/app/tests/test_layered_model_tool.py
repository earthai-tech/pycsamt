# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for LayeredModelDialog.

LayeredModelDialog is a fully synchronous, self-contained 1-D layered-earth
model builder — no survey data, no QThread worker. It is constructed with
just ``parent=None`` and every action (add/remove row, preset, preview, OK)
runs inline in the calling thread, so every method can be exercised
directly against real ``pycsamt.forward.synthetic.LayeredModel`` objects
with no fakes or monkeypatching required (except where noted).
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.tools.layered_model_tool import (
    _DEFAULT_LAYERS,
    _PRESETS,
    LayeredModelDialog,
)

# ── fixture ──────────────────────────────────────────────────────────────────


@pytest.fixture
def dlg(qapp):
    d = LayeredModelDialog(parent=None)
    yield d
    d.close()


# ── construction / _load_default ────────────────────────────────────────────


class TestDialogBuild:
    def test_creates(self, dlg):
        assert dlg is not None

    def test_window_title(self, dlg):
        assert dlg.windowTitle() == "Layered Model Builder"

    def test_model_initially_none(self, dlg):
        assert dlg.model is None

    def test_preset_combo_items(self, dlg):
        items = [
            dlg._preset_combo.itemText(i) for i in range(dlg._preset_combo.count())
        ]
        assert items == _PRESETS

    def test_n_spin_default_and_range(self, dlg):
        assert dlg._n_spin.value() == 3
        assert dlg._n_spin.minimum() == 2
        assert dlg._n_spin.maximum() == 20

    def test_table_has_three_columns(self, dlg):
        assert dlg._table.columnCount() == 3

    def test_table_header_labels(self, dlg):
        labels = [dlg._table.horizontalHeaderItem(i).text() for i in range(3)]
        assert labels == ["Layer", "ρ (Ω·m)", "Thickness (m)"]

    def test_default_layers_populate_table(self, dlg):
        assert dlg._table.rowCount() == len(_DEFAULT_LAYERS)

    def test_default_layer_cell_values(self, dlg):
        for i, (rho, thick) in enumerate(_DEFAULT_LAYERS):
            assert dlg._table.item(i, 0).text() == (
                "Halfspace" if thick is None else f"Layer {i + 1}"
            )
            assert dlg._table.item(i, 1).text() == f"{rho:.1f}"
            expected_thick = "" if thick is None else f"{thick:.1f}"
            assert dlg._table.item(i, 2).text() == expected_thick

    def test_last_default_row_is_halfspace(self, dlg):
        assert dlg._table.item(2, 0).text() == "Halfspace"
        assert dlg._table.item(2, 2).text() == ""

    def test_status_label_populated_after_load_default(self, dlg):
        # _load_default() calls _on_preview() which sets the status label.
        assert "layers" in dlg._status_lbl.text()

    def test_canvas_exists_and_drawn_after_init(self, dlg):
        assert dlg._canvas is not None
        assert len(dlg._canvas.figure.axes) >= 1


# ── _append_row ──────────────────────────────────────────────────────────────


class TestAppendRow:
    def test_append_row_with_thickness(self, dlg):
        n0 = dlg._table.rowCount()
        dlg._append_row(n0 + 1, 42.5, 123.4)
        r = dlg._table.rowCount() - 1
        assert dlg._table.item(r, 0).text() == f"Layer {n0 + 1}"
        assert dlg._table.item(r, 1).text() == "42.5"
        assert dlg._table.item(r, 2).text() == "123.4"

    def test_append_row_halfspace_thick_none(self, dlg):
        n0 = dlg._table.rowCount()
        dlg._append_row(n0 + 1, 7.0, None)
        r = dlg._table.rowCount() - 1
        assert dlg._table.item(r, 0).text() == "Halfspace"
        assert dlg._table.item(r, 2).text() == ""

    def test_append_row_increments_row_count(self, dlg):
        n0 = dlg._table.rowCount()
        dlg._append_row(n0 + 1, 1.0, 1.0)
        assert dlg._table.rowCount() == n0 + 1

    def test_append_row_rho_formatted_one_decimal(self, dlg):
        dlg._append_row(99, 1000.0, 50.0)
        r = dlg._table.rowCount() - 1
        assert dlg._table.item(r, 1).text() == "1000.0"


# ── _add_row ─────────────────────────────────────────────────────────────────


class TestAddRow:
    def test_add_row_increases_row_count(self, dlg):
        n0 = dlg._table.rowCount()
        dlg._add_row()
        assert dlg._table.rowCount() == n0 + 1

    def test_add_row_defaults(self, dlg):
        dlg._add_row()
        r = dlg._table.rowCount() - 1
        assert dlg._table.item(r, 1).text() == "100.0"
        assert dlg._table.item(r, 2).text() == "500.0"

    def test_add_row_after_default_halfspace_breaks_read(self, dlg):
        """Real behaviour (documented, not fixed): the default table's last
        row is the halfspace (empty thickness cell). ``_add_row`` appends a
        new *non*-halfspace row after it but never fills in a thickness for
        the row that is no longer last, so ``_read_model`` (invoked by the
        ``_on_preview`` this triggers) fails with a validation error instead
        of silently building a 4-layer model. See bug note in test module
        docstring / final report."""
        dlg._add_row()
        assert "Invalid thickness in row 3" in dlg._status_lbl.text()

    def test_add_row_after_filling_gap_thickness_succeeds(self, dlg):
        """Workaround for the bug above: manually filling in the thickness
        left blank by the former halfspace row lets the model build."""
        dlg._add_row()
        dlg._table.item(2, 2).setText("500.0")  # former halfspace row
        dlg._on_preview()
        assert "4 layers" in dlg._status_lbl.text()


# ── _remove_row ──────────────────────────────────────────────────────────────


class TestRemoveRow:
    def test_remove_row_decreases_count(self, dlg):
        n0 = dlg._table.rowCount()
        assert n0 == 3
        dlg._remove_row()
        assert dlg._table.rowCount() == n0 - 1

    def test_remove_row_blocked_at_two_rows(self, dlg):
        dlg._remove_row()  # 3 -> 2
        assert dlg._table.rowCount() == 2
        dlg._remove_row()  # blocked: n > 2 is False at n == 2
        assert dlg._table.rowCount() == 2

    def test_remove_row_never_reaches_zero(self, dlg):
        for _ in range(10):
            dlg._remove_row()
        assert dlg._table.rowCount() >= 2


# ── _on_preset / _apply_model_preset ────────────────────────────────────────


class TestOnPreset:
    @pytest.mark.parametrize("label", ["Random", "Blocky", "Smooth"])
    def test_each_preset_populates_table_with_n_rows(self, dlg, label):
        dlg._n_spin.setValue(5)
        idx = _PRESETS.index(label)
        dlg._preset_combo.setCurrentIndex(idx)
        assert dlg._table.rowCount() == 5

    def test_preset_last_row_is_halfspace(self, dlg):
        dlg._n_spin.setValue(4)
        dlg._preset_combo.setCurrentIndex(_PRESETS.index("Blocky"))
        last = dlg._table.rowCount() - 1
        assert dlg._table.item(last, 0).text() == "Halfspace"
        assert dlg._table.item(last, 2).text() == ""

    def test_placeholder_index_dispatches_nothing(self, dlg):
        n0 = dlg._table.rowCount()
        dlg._preset_combo.setCurrentIndex(0)  # "— choose preset —"
        assert dlg._table.rowCount() == n0

    def test_on_preset_direct_call_random(self, dlg):
        dlg._n_spin.setValue(6)
        dlg._on_preset(0)  # idx arg unused internally; reads combo text
        # combo still at default placeholder unless changed — verify no-op
        # then explicitly select Random and call again.
        dlg._preset_combo.setCurrentIndex(_PRESETS.index("Random"))
        dlg._on_preset(_PRESETS.index("Random"))
        assert dlg._table.rowCount() == 6


class TestApplyModelPresetDirect:
    def test_random_kind(self, dlg):
        dlg._apply_model_preset("random", 5)
        assert dlg._table.rowCount() == 5
        assert dlg._table.item(4, 0).text() == "Halfspace"

    def test_blocky_kind(self, dlg):
        dlg._apply_model_preset("blocky", 4)
        assert dlg._table.rowCount() == 4

    def test_smooth_kind(self, dlg):
        dlg._apply_model_preset("smooth", 7)
        assert dlg._table.rowCount() == 7

    def test_preset_error_sets_status_and_leaves_table(self, dlg, monkeypatch):
        import pycsamt.forward.synthetic as synth_mod

        def _boom(*a, **k):
            raise RuntimeError("synthetic preset failure")

        monkeypatch.setattr(synth_mod.LayeredModel, "random", staticmethod(_boom))
        n0 = dlg._table.rowCount()
        dlg._apply_model_preset("random", 5)
        assert "Preset error" in dlg._status_lbl.text()
        assert "synthetic preset failure" in dlg._status_lbl.text()
        # table untouched since the exception happens before setRowCount(0)
        assert dlg._table.rowCount() == n0


# ── _read_model ──────────────────────────────────────────────────────────────


class TestReadModel:
    def test_read_model_default_three_rows(self, dlg):
        m = dlg._read_model()
        assert m.n_layers == 3
        assert list(m.resistivity) == [100.0, 10.0, 500.0]
        assert list(m.thickness) == [300.0, 800.0]

    def test_read_model_single_row(self, dlg):
        dlg._table.setRowCount(0)
        dlg._append_row(1, 250.0, None)
        m = dlg._read_model()
        assert m.n_layers == 1
        assert list(m.resistivity) == [250.0]
        assert list(m.thickness) == []

    def test_read_model_many_rows(self, dlg):
        dlg._table.setRowCount(0)
        for i in range(6):
            dlg._append_row(i + 1, 10.0 * (i + 1), 100.0 if i < 5 else None)
        m = dlg._read_model()
        assert m.n_layers == 6
        assert len(m.thickness) == 5

    def test_read_model_last_row_thickness_ignored(self, dlg):
        """Even if the last row's thickness cell holds a value, it is not
        read — only rows before the last contribute a thickness."""
        dlg._table.setRowCount(0)
        dlg._append_row(1, 100.0, 300.0)
        dlg._append_row(2, 10.0, 800.0)  # give the "halfspace" a thickness
        m = dlg._read_model()
        assert m.n_layers == 2
        assert list(m.thickness) == [300.0]

    def test_read_model_invalid_resistivity_raises(self, dlg):
        dlg._table.item(0, 1).setText("not-a-number")
        with pytest.raises(ValueError, match="Invalid resistivity in row 1"):
            dlg._read_model()

    def test_read_model_invalid_thickness_raises(self, dlg):
        dlg._table.item(0, 2).setText("not-a-number")
        with pytest.raises(ValueError, match="Invalid thickness in row 1"):
            dlg._read_model()

    def test_read_model_empty_thickness_on_non_last_row_raises(self, dlg):
        dlg._table.item(0, 2).setText("")
        with pytest.raises(ValueError, match="Invalid thickness in row 1"):
            dlg._read_model()

    def test_read_model_missing_rho_item_raises(self, dlg):
        # Remove the QTableWidgetItem entirely (None) rather than blanking text.
        dlg._table.setItem(0, 1, None)
        with pytest.raises(ValueError, match="Invalid resistivity in row 1"):
            dlg._read_model()


# ── _on_preview ──────────────────────────────────────────────────────────────


class TestOnPreview:
    def test_on_preview_valid_model_updates_status(self, dlg):
        dlg._on_preview()
        text = dlg._status_lbl.text()
        assert "3 layers" in text
        assert "Ω·m" in text
        assert "depth" in text

    def test_on_preview_draws_into_canvas(self, dlg):
        dlg._on_preview()
        assert len(dlg._canvas.figure.axes) >= 1

    def test_on_preview_invalid_model_sets_error_status_no_raise(self, dlg):
        dlg._table.item(0, 1).setText("garbage")
        dlg._on_preview()  # must not raise
        assert "Invalid resistivity in row 1" in dlg._status_lbl.text()

    def test_on_preview_after_add_row_surfaces_gap_thickness_error(self, dlg):
        """See TestAddRow bug note: adding a row after the default
        halfspace leaves an invalid gap, surfaced here through the status
        label rather than a raised exception."""
        dlg._add_row()  # _add_row already calls _on_preview() internally
        assert "Invalid thickness in row 3" in dlg._status_lbl.text()

    def test_on_preview_after_remove_row(self, dlg):
        dlg._remove_row()
        dlg._on_preview()
        assert "2 layers" in dlg._status_lbl.text()


# ── _on_ok ───────────────────────────────────────────────────────────────────


class TestOnOk:
    def test_on_ok_valid_model_sets_model_and_accepts(self, dlg):
        result_holder = {}
        dlg.accept = lambda: result_holder.setdefault("accepted", True)
        dlg._on_ok()
        assert dlg.model is not None
        assert dlg.model.n_layers == 3
        assert result_holder.get("accepted") is True

    def test_on_ok_invalid_model_does_not_accept(self, dlg):
        dlg._table.item(0, 1).setText("garbage")
        accepted = {"called": False}
        dlg.accept = lambda: accepted.__setitem__("called", True)
        dlg._on_ok()
        assert accepted["called"] is False
        assert dlg.model is None
        assert "Cannot build model" in dlg._status_lbl.text()
        assert "Invalid resistivity in row 1" in dlg._status_lbl.text()

    def test_on_ok_model_matches_table_contents(self, dlg):
        dlg._table.setRowCount(0)
        dlg._append_row(1, 55.0, 200.0)
        dlg._append_row(2, 5.0, None)
        dlg.accept = lambda: None
        dlg._on_ok()
        assert list(dlg.model.resistivity) == [55.0, 5.0]
        assert list(dlg.model.thickness) == [200.0]
