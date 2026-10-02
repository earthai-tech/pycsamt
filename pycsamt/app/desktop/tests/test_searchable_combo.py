# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for SearchableComboBox (station picker with a live-filter popup)."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QKeyEvent
from PySide6.QtTest import QTest

from pycsamt.app.desktop.widgets.searchable_combo import SearchableComboBox


def _names():
    return ["18-001A", "18-002U", "18-003A", "19-010B"]


def _key(key, text=""):
    return QKeyEvent(QKeyEvent.Type.KeyPress, key, Qt.KeyboardModifier.NoModifier, text)


# ── construction / placeholder ───────────────────────────────────────────


def test_starts_with_placeholder(qapp):
    combo = SearchableComboBox()
    assert combo.count() == 1
    assert combo.itemText(0) == SearchableComboBox._PLACEHOLDER
    assert combo.current_station() == ""


def test_set_names_populates_items_with_count_placeholder(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    assert combo.count() == len(_names()) + 1
    assert combo.itemText(0) == "— select station (4) —"
    assert combo.currentIndex() == 0
    assert combo.current_station() == ""


def test_set_names_empty_list_keeps_bare_placeholder(qapp):
    combo = SearchableComboBox()
    combo.set_names(["18-001A"])
    combo.set_names([])
    assert combo.itemText(0) == SearchableComboBox._PLACEHOLDER
    assert combo.count() == 1


def test_set_names_resets_current_selection(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo.select_station("18-002U")
    assert combo.current_station() == "18-002U"
    combo.set_names(_names())
    assert combo.current_station() == ""


def test_set_names_does_not_emit_signal(qapp):
    combo = SearchableComboBox()
    received = []
    combo.station_selected.connect(received.append)
    combo.set_names(_names())
    assert received == []


# ── select_station / current_station ─────────────────────────────────────


def test_select_station_sets_index_without_signal(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    received = []
    combo.station_selected.connect(received.append)
    combo.select_station("18-003A")
    assert combo.current_station() == "18-003A"
    assert received == []


def test_select_station_unknown_name_is_noop(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo.select_station("does-not-exist")
    assert combo.current_station() == ""


def test_current_station_empty_when_placeholder_selected(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo.setCurrentIndex(0)
    assert combo.current_station() == ""


# ── popup show/hide wiring ────────────────────────────────────────────────


def test_show_popup_resets_search_and_sizes_popup(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo.resize(250, 30)
    combo._popup._search.setText("stale query")
    combo.showPopup()
    try:
        assert combo._popup._search.text() == ""
        assert combo._popup.isVisible()
        assert combo._popup.width() >= 220
    finally:
        combo.hidePopup()


def test_show_popup_uses_minimum_width_when_combo_narrow(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo.resize(50, 30)
    combo.showPopup()
    try:
        assert combo._popup.width() == 220
    finally:
        combo.hidePopup()


def test_hide_popup_hides_internal_popup(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo.showPopup()
    assert combo._popup.isVisible()
    combo.hidePopup()
    assert not combo._popup.isVisible()


# ── popup selection -> combo update ───────────────────────────────────────


def test_popup_chosen_updates_combo_and_emits_signal(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    received = []
    combo.station_selected.connect(received.append)
    combo._on_popup_chosen("18-002U")
    assert combo.current_station() == "18-002U"
    assert received == ["18-002U"]
    assert not combo._popup.isVisible()


def test_popup_item_click_triggers_selection(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    received = []
    combo.station_selected.connect(received.append)
    combo.showPopup()
    item = combo._popup._list.item(1)
    combo._popup._list.itemClicked.emit(item)
    assert received == [item.text()]
    combo.hidePopup()


# ── _StationPopup: filtering ──────────────────────────────────────────────


def test_popup_filter_matches_substring_case_insensitive(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo._popup._search.setText("002")
    assert combo._popup._list.count() == 1
    assert combo._popup._list.item(0).text() == "18-002U"


def test_popup_filter_empty_text_restores_full_list(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo._popup._search.setText("003")
    combo._popup._search.setText("")
    assert combo._popup._list.count() == len(_names())


def test_popup_filter_no_match_empties_list(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo._popup._search.setText("zzz-nope")
    assert combo._popup._list.count() == 0


def test_popup_filter_auto_highlights_first_match(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo._popup._search.setText("18-")
    assert combo._popup._list.currentRow() == 0


def test_popup_reset_search_restores_all_names(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo._popup._search.setText("002")
    combo._popup.reset_search()
    assert combo._popup._search.text() == ""
    assert combo._popup._list.count() == len(_names())


def test_popup_preferred_height_scales_with_max_visible(qapp):
    popup = combo_popup = SearchableComboBox(max_visible=2)._popup
    combo_popup.set_names(_names())
    h_small = popup.preferred_height()

    popup2 = SearchableComboBox(max_visible=12)._popup
    popup2.set_names(_names())
    h_large = popup2.preferred_height()
    assert h_large >= h_small


def test_popup_preferred_height_with_no_names_uses_minimum_row(qapp):
    popup = SearchableComboBox()._popup
    popup.set_names([])
    assert popup.preferred_height() > 0


# ── _StationPopup: key navigation ─────────────────────────────────────────


def test_popup_enter_confirms_selected_row(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    popup = combo._popup
    popup.set_names(_names())
    popup._list.setCurrentRow(2)
    received = []
    popup.station_chosen.connect(received.append)
    popup.keyPressEvent(_key(Qt.Key.Key_Return))
    assert received == ["18-003A"]


def test_popup_enter_with_no_selection_confirms_first_row(qapp):
    combo = SearchableComboBox()
    popup = combo._popup
    popup.set_names(_names())
    popup._list.setCurrentRow(-1)
    received = []
    popup.station_chosen.connect(received.append)
    popup.keyPressEvent(_key(Qt.Key.Key_Enter))
    assert received == ["18-001A"]


def test_popup_enter_with_empty_list_emits_nothing(qapp):
    combo = SearchableComboBox()
    popup = combo._popup
    popup.set_names([])
    received = []
    popup.station_chosen.connect(received.append)
    popup.keyPressEvent(_key(Qt.Key.Key_Return))
    assert received == []


def test_popup_escape_hides_popup(qapp):
    combo = SearchableComboBox()
    combo.set_names(_names())
    combo.showPopup()
    assert combo._popup.isVisible()
    combo._popup.keyPressEvent(_key(Qt.Key.Key_Escape))
    assert not combo._popup.isVisible()


def test_popup_down_arrow_moves_selection_and_focuses_list(qapp):
    combo = SearchableComboBox()
    popup = combo._popup
    popup.set_names(_names())
    popup._list.setCurrentRow(0)
    popup.keyPressEvent(_key(Qt.Key.Key_Down))
    assert popup._list.currentRow() == 1


def test_popup_down_arrow_at_last_row_stays_put(qapp):
    combo = SearchableComboBox()
    popup = combo._popup
    popup.set_names(_names())
    last = popup._list.count() - 1
    popup._list.setCurrentRow(last)
    popup.keyPressEvent(_key(Qt.Key.Key_Down))
    assert popup._list.currentRow() == last


def test_popup_up_arrow_moves_selection_up(qapp):
    combo = SearchableComboBox()
    popup = combo._popup
    popup.set_names(_names())
    popup._list.setCurrentRow(2)
    popup.keyPressEvent(_key(Qt.Key.Key_Up))
    assert popup._list.currentRow() == 1


def test_popup_up_arrow_at_first_row_stays_put(qapp):
    combo = SearchableComboBox()
    popup = combo._popup
    popup.set_names(_names())
    popup._list.setCurrentRow(0)
    popup.keyPressEvent(_key(Qt.Key.Key_Up))
    assert popup._list.currentRow() == 0


def test_popup_other_key_routes_to_search_field(qapp):
    combo = SearchableComboBox()
    popup = combo._popup
    popup.set_names(_names())
    popup.show()
    try:
        popup._search.clearFocus()
        popup.keyPressEvent(_key(Qt.Key.Key_A, "a"))
        assert popup._search.hasFocus()
    finally:
        popup.hide()


def test_popup_other_key_when_search_already_focused_still_routes(qapp):
    combo = SearchableComboBox()
    popup = combo._popup
    popup.set_names(_names())
    popup.show()
    try:
        popup._search.setFocus()
        QTest.qWait(0)
        popup.keyPressEvent(_key(Qt.Key.Key_B, "b"))
        assert popup._search.hasFocus()
    finally:
        popup.hide()
