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


def test_five_tabs_created(dlg):
    assert len(dlg._pages) == 5
    assert dlg._tabs.count() == 5


def test_tab_labels(dlg):
    labels = [dlg._tabs.tabText(i) for i in range(dlg._tabs.count())]
    assert labels == [
        "Pseudosections",
        "View Controls",
        "Display",
        "Topography",
        "Interpretation",
    ]


def test_snapshot_taken_at_construction(dlg, ctrl):
    assert dlg._snapshot == ctrl.snapshot()


# ── open_tab ──────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "key,idx",
    [
        ("pseudosections", 0),
        ("view_controls", 1),
        ("display", 2),
        ("topography", 3),
        ("interpretation", 4),
    ],
)
def test_open_tab_valid_keys(dlg, key, idx):
    dlg.open_tab(key)
    assert dlg._tabs.currentIndex() == idx


def test_open_tab_invalid_key_noop(dlg):
    dlg.open_tab("view_controls")
    dlg.open_tab("not_a_real_tab")
    assert dlg._tabs.currentIndex() == 1


def test_constructor_open_tab_kwarg(qapp, ctrl):
    from pycsamt.app.desktop.dialogs.settings_dialog import APIConfigDialog

    d = APIConfigDialog(ctrl, open_tab="display")
    assert d._tabs.currentIndex() == 2
    d.close()


def test_constructor_no_open_tab_defaults_to_first(dlg):
    assert dlg._tabs.currentIndex() == 0


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
    dlg._tabs.setCurrentIndex(0)
    spy = mock.Mock()
    dlg.settings_changed.connect(spy)
    dlg._on_reset_tab()
    assert stub.reset_called is True
    spy.assert_called_once_with(["pseudosections"])


def test_on_reset_tab_reset_raising_is_swallowed(dlg):
    stub = _StubPage(reset_raises=RuntimeError("boom"))
    dlg._pages = [_StubPage(), stub] + [_StubPage() for _ in range(3)]
    dlg._tabs.setCurrentIndex(1)
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
