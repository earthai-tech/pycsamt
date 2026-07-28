# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for OverlayCombo and CategoryComboBar widgets."""

from __future__ import annotations

from pycsamt.app.desktop.widgets.category_bar import CategoryComboBar
from pycsamt.app.desktop.widgets.overlay_combo import (
    OVERLAY_OPTIONS,
    OverlayCombo,
)

# ── OverlayCombo ─────────────────────────────────────────────────────────────


def test_overlay_combo_populates_options(qapp):
    combo = OverlayCombo()
    assert combo.count() == len(OVERLAY_OPTIONS)
    assert combo.current_overlay == OVERLAY_OPTIONS[0]


def test_overlay_combo_emits_signal_on_change(qapp):
    combo = OverlayCombo()
    received = []
    combo.overlay_changed.connect(received.append)
    combo.setCurrentText("Phase")
    assert received == ["Phase"]
    assert combo.current_overlay == "Phase"


# ── CategoryComboBar ─────────────────────────────────────────────────────────


def test_category_combo_bar_builds_with_labels(qapp):
    bar = CategoryComboBar(
        category_label="Category:", item_label="Item:", show_item=True
    )
    assert bar._lbl_item is not None
    assert not bar._combo_item.isHidden()


def test_category_combo_bar_builds_without_item_combo(qapp):
    bar = CategoryComboBar(show_item=False)
    assert bar._combo_item.isHidden()


def test_category_combo_bar_set_categories_and_items(qapp):
    bar = CategoryComboBar(category_label="Cat", item_label="Item")
    bar.set_categories(["A", "B", "C"])
    bar.set_items(["x", "y"])
    assert bar.current_category_text() == "A"
    assert bar.current_item_text() == "x"

    bar.set_current_category(2)
    bar.set_current_item(1)
    assert bar.current_category() == 2
    assert bar.current_item() == 1
    assert bar.current_category_text() == "C"
    assert bar.current_item_text() == "y"


def test_category_combo_bar_signals(qapp):
    bar = CategoryComboBar(category_label="Cat", item_label="Item")
    bar.set_categories(["A", "B"])
    bar.set_items(["x", "y"])

    cat_events = []
    item_events = []
    bar.category_changed.connect(cat_events.append)
    bar.item_changed.connect(item_events.append)

    bar.set_current_category(1)
    bar.set_current_item(1)

    assert cat_events == [1]
    assert item_events == [1]


def test_category_combo_bar_refresh_and_export_signals(qapp):
    bar = CategoryComboBar()
    refresh_events = []
    export_events = []
    bar.refresh_clicked.connect(lambda: refresh_events.append(True))
    bar.export_clicked.connect(lambda: export_events.append(True))

    bar._btn_refresh.click()
    bar._btn_export.click()

    assert refresh_events == [True]
    assert export_events == [True]


def test_category_combo_bar_visibility_toggles(qapp):
    bar = CategoryComboBar()
    bar.set_export_visible(False)
    assert bar._btn_export.isHidden()
    bar.set_refresh_visible(False)
    assert bar._btn_refresh.isHidden()
    bar.set_item_combo_visible(True)
    assert not bar._combo_item.isHidden()
