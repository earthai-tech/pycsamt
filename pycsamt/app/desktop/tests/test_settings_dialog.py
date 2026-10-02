# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for APIConfigDialog.

The dialog's own logic (tab dispatch, apply/cancel/reset orchestration) is
what's under test here — the underlying SettingsController and its
PYCSAMT_* singletons are covered separately in test_settings_controller.py.
To keep this file fast and decoupled from that singleton state, the real
per-tab pages are swapped for lightweight stand-ins after construction.
"""

from __future__ import annotations

from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")


@pytest.fixture
def ctrl():
    from pycsamt.app.desktop.controllers.settings_controller import (
        SettingsController,
    )

    return SettingsController()


@pytest.fixture
def dlg(qapp, ctrl):
    from pycsamt.app.desktop.dialogs.settings_dialog import APIConfigDialog

    d = APIConfigDialog(ctrl)
    yield d
    d.close()


class _StubPage:
    def __init__(self, collect_return=None, collect_raises=None, reset_raises=None):
        self._collect_return = collect_return if collect_return is not None else {}
        self._collect_raises = collect_raises
        self._reset_raises = reset_raises
        self.reset_called = False

    def collect(self):
        if self._collect_raises:
            raise self._collect_raises
        return self._collect_return

    def reset(self):
        self.reset_called = True
        if self._reset_raises:
            raise self._reset_raises


# ── Construction ──────────────────────────────────────────────────────────


def test_creates(dlg):
    assert dlg is not None


def test_window_title(dlg):
    assert dlg.windowTitle() == "API Configuration"


def test_all_pages_created(dlg):
    assert len(dlg._pages) == 9
    assert dlg._stack.count() == 9


def test_page_labels(dlg):
    assert dlg.page_labels() == [
        "Pseudosections",
        "View Controls",
        "Display",
        "Contours & Mesh",
        "Topography",
        "Interpretation",
        "Data Ordering",
        "Plot & Export",
        "Pipeline",
    ]


def test_sidebar_groups_pages(dlg):
    from pycsamt.app.desktop.dialogs.settings_dialog import _PAGE_ROLE

    headers = [
        dlg._nav.item(r).text()
        for r in range(dlg._nav.count())
        if dlg._nav.item(r).data(_PAGE_ROLE) is None
    ]
    assert headers == ["PLOTS", "DATA", "OUTPUT"]


def test_every_page_scrolls_instead_of_squeezing(dlg):
    """Regression: fixed-height tabs squeezed View Controls until its spin
    boxes and combos overlapped. Pages now live in scroll areas."""
    from PySide6.QtWidgets import QScrollArea

    for i, page in enumerate(dlg._pages):
        holder = dlg._stack.widget(i)
        assert isinstance(holder, QScrollArea)
        assert holder.widget() is page
        assert holder.widgetResizable()


def test_nav_click_switches_page_and_header(dlg):
    from pycsamt.app.desktop.dialogs.settings_dialog import _PAGE_ROLE

    item = next(
        dlg._nav.item(r) for r in range(dlg._nav.count())
        if dlg._nav.item(r).data(_PAGE_ROLE) == 8
    )
    dlg._nav.setCurrentItem(item)
    assert dlg._stack.currentIndex() == 8
    assert dlg.current_key() == "pipeline"
    assert dlg._page_title.text() == "Pipeline"


def test_search_filters_sidebar_by_setting_labels(dlg):
    from pycsamt.app.desktop.dialogs.settings_dialog import _PAGE_ROLE

    dlg._search.setText("dpi")  # a field label, not a page title
    shown = [
        dlg._nav.item(r).text().strip()
        for r in range(dlg._nav.count())
        if not dlg._nav.item(r).isHidden()
        and dlg._nav.item(r).data(_PAGE_ROLE) is not None
    ]
    assert "Pipeline" in shown and "Plot & Export" in shown
    assert "Topography" not in shown
    dlg._search.clear()
    assert all(not dlg._nav.item(r).isHidden() for r in range(dlg._nav.count()))


def test_snapshot_taken_at_construction(dlg, ctrl):
    assert dlg._snapshot == ctrl.snapshot()


# ── open_tab ──────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "key,idx",
    [
        ("pseudosections", 0),
        ("view_controls", 1),
        ("display", 2),
        ("rendering", 3),
        ("topography", 4),
        ("interpretation", 5),
        ("ordering", 6),
        ("output", 7),
        ("pipeline", 8),
    ],
)
def test_open_tab_valid_keys(dlg, key, idx):
    dlg.open_tab(key)
    assert dlg._stack.currentIndex() == idx
    assert dlg.current_key() == key


def test_open_tab_invalid_key_noop(dlg):
    dlg.open_tab("view_controls")
    dlg.open_tab("not_a_real_tab")
    assert dlg._stack.currentIndex() == 1


def test_constructor_open_tab_kwarg(qapp, ctrl):
    from pycsamt.app.desktop.dialogs.settings_dialog import APIConfigDialog

    d = APIConfigDialog(ctrl, open_tab="display")
    assert d._stack.currentIndex() == 2
    d.close()


def test_constructor_no_open_tab_defaults_to_first(dlg):
    assert dlg._stack.currentIndex() == 0


# ── _on_apply ─────────────────────────────────────────────────────────────


def test_on_apply_no_fields_emits_nothing(dlg):
    dlg._pages = [_StubPage(collect_return={}) for _ in range(5)]
    spy = mock.Mock()
    dlg.settings_changed.connect(spy)
    dlg._on_apply()
    spy.assert_not_called()


def test_on_apply_touches_matching_method(dlg):
    dlg._pages = [_StubPage() for _ in range(5)]
    dlg._pages[1] = _StubPage(collect_return={"view_controls": {"some_flag": True}})
    dlg._ctrl.apply_view_controls = mock.Mock()
    spy = mock.Mock()
    dlg.settings_changed.connect(spy)
    dlg._on_apply()
    dlg._ctrl.apply_view_controls.assert_called_once_with(some_flag=True)
    spy.assert_called_once_with(["view_controls"])


def test_on_apply_skips_empty_fields_dict(dlg):
    dlg._pages = [_StubPage() for _ in range(5)]
    dlg._pages[0] = _StubPage(collect_return={"pseudosections": {}})
    dlg._ctrl.apply_pseudosections = mock.Mock()
    dlg._on_apply()
    dlg._ctrl.apply_pseudosections.assert_not_called()


def test_on_apply_unknown_apply_key_skipped(dlg):
    dlg._pages = [_StubPage() for _ in range(5)]
    dlg._pages[0] = _StubPage(collect_return={"not_a_real_key": {"x": 1}})
    spy = mock.Mock()
    dlg.settings_changed.connect(spy)
    dlg._on_apply()
    spy.assert_not_called()


def test_on_apply_method_raising_is_swallowed(dlg):
    dlg._pages = [_StubPage() for _ in range(5)]
    dlg._pages[2] = _StubPage(collect_return={"display": {"x": 1}})
    dlg._ctrl.apply_display = mock.Mock(side_effect=RuntimeError("boom"))
    spy = mock.Mock()
    dlg.settings_changed.connect(spy)
    dlg._on_apply()  # must not raise
    spy.assert_not_called()


def test_on_apply_collect_raising_is_swallowed(dlg):
    dlg._pages = [_StubPage() for _ in range(5)]
    dlg._pages[0] = _StubPage(collect_raises=RuntimeError("bad collect"))
    dlg._on_apply()  # must not raise


def test_on_apply_multiple_pages_touched(dlg):
    dlg._pages = [_StubPage() for _ in range(5)]
    dlg._pages[0] = _StubPage(collect_return={"pseudosections": {"a": 1}})
    dlg._pages[3] = _StubPage(collect_return={"topography": {"b": 2}})
    dlg._ctrl.apply_pseudosections = mock.Mock()
    dlg._ctrl.apply_topography = mock.Mock()
    spy = mock.Mock()
    dlg.settings_changed.connect(spy)
    dlg._on_apply()
    assert spy.call_args[0][0] == ["pseudosections", "topography"]


# ── _on_ok / _on_cancel ───────────────────────────────────────────────────


def test_on_ok_applies_then_accepts(dlg):
    from PySide6.QtWidgets import QDialog

    dlg._pages = [_StubPage() for _ in range(5)]
    dlg._on_ok()
    assert dlg.result() == QDialog.DialogCode.Accepted


def test_on_cancel_restores_snapshot_and_rejects(dlg):
    from PySide6.QtWidgets import QDialog

    dlg._ctrl.restore = mock.Mock()
    dlg._on_cancel()
    dlg._ctrl.restore.assert_called_once_with(dlg._snapshot)
    assert dlg.result() == QDialog.DialogCode.Rejected


def test_on_cancel_restore_raising_is_swallowed_and_still_rejects(dlg):
    from PySide6.QtWidgets import QDialog

    dlg._ctrl.restore = mock.Mock(side_effect=RuntimeError("boom"))
    dlg._on_cancel()  # must not raise
    assert dlg.result() == QDialog.DialogCode.Rejected


# ── _on_reset_tab ─────────────────────────────────────────────────────────


def test_on_reset_tab_calls_page_reset_and_emits(dlg):
    stub = _StubPage()
    dlg._pages = [stub] + [_StubPage() for _ in range(4)]
    dlg._stack.setCurrentIndex(0)
    spy = mock.Mock()
    dlg.settings_changed.connect(spy)
    dlg._on_reset_tab()
    assert stub.reset_called is True
    spy.assert_called_once_with(["station", "section"])


def test_on_reset_tab_reset_raising_is_swallowed(dlg):
    stub = _StubPage(reset_raises=RuntimeError("boom"))
    dlg._pages = [_StubPage(), stub] + [_StubPage() for _ in range(3)]
    dlg._stack.setCurrentIndex(1)
    spy = mock.Mock()
    dlg.settings_changed.connect(spy)
    dlg._on_reset_tab()  # must not raise
    spy.assert_called_once_with(["view_controls"])


# ── Button wiring ─────────────────────────────────────────────────────────


def test_reset_button_triggers_on_reset_tab(dlg):
    dlg._pages = [_StubPage() for _ in range(5)]
    with mock.patch.object(dlg, "_on_reset_tab") as m:
        dlg._reset_btn.click()
        m.assert_called_once()
