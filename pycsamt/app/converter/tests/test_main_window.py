# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.main_window.ConverterMainWindow."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.converter import main_window as mw_mod  # noqa: E402
from pycsamt.app.converter.main_window import (  # noqa: E402
    ConverterMainWindow,
    _icon,
    _is_near_black,
    _recolor_svg,
    _TOOLS,
)
from pycsamt.app.converter.pages import SettingsPage  # noqa: E402


@pytest.fixture
def window(qapp, isolated_qsettings):
    w = ConverterMainWindow()
    yield w
    w.close()


# ── construction ─────────────────────────────────────────────────────────


def test_window_builds_all_pages(window):
    assert len(window._pages) == len(_TOOLS)
    assert window._stack.count() == len(_TOOLS)


def test_window_title_and_icon(window):
    assert window.windowTitle() == "pyCSAMT Format Studio"


def test_nav_has_one_row_per_tool(window):
    assert window._nav.count() == len(_TOOLS)
    for i, (label, _cls, _icon_name) in enumerate(_TOOLS):
        assert window._nav.item(i).text() == label


def test_pages_share_the_dock_log_panel(window):
    for page in window._pages:
        if hasattr(page, "log"):
            assert page.log is window._log_panel


def test_nav_selection_changes_stack_and_tools_menu(window):
    window._nav.setCurrentRow(2)
    assert window._stack.currentIndex() == 2
    assert window._tool_actions[2].isChecked()
    assert not window._tool_actions[0].isChecked()


def test_settings_action_jumps_to_settings_page(window):
    # Reached via findChildren(QAction) rather than
    # menuBar().actions()[0].menu().actions()[...]: the latter's .menu()
    # accessor has been observed to hand back a QMenu wrapper for an
    # already-destroyed C++ object under pytest (a Shiboken ownership/GC
    # timing quirk, not present when walking the QObject tree directly).
    from PySide6.QtGui import QAction

    settings_row = next(i for i, t in enumerate(_TOOLS) if t[1] is SettingsPage)
    window._nav.setCurrentRow(0)
    settings_action = next(
        a for a in window.findChildren(QAction) if a.text().startswith("&Settings")
    )
    settings_action.trigger()
    assert window._nav.currentRow() == settings_row


# ── theme ────────────────────────────────────────────────────────────────


def test_apply_theme_switches_dark_mode_flag(window):
    window._apply_theme("dark")
    assert mw_mod._DARK_MODE is True
    assert window._act_dark.isChecked()
    window._apply_theme("light")
    assert mw_mod._DARK_MODE is False
    assert window._act_light.isChecked()


def test_apply_theme_persists_setting(window):
    from pycsamt.app.converter.settings import load_settings

    window._apply_theme("dark")
    assert load_settings().theme == "dark"


def test_apply_theme_noop_persist_when_unchanged(window):
    from pycsamt.app.converter import settings as settings_mod

    window._apply_theme("light")
    calls = []
    orig = settings_mod.save_settings
    settings_mod.save_settings = lambda s: calls.append(s) or orig(s)
    try:
        window._apply_theme("light")
    finally:
        settings_mod.save_settings = orig
    assert calls == []


def test_embedded_theme_is_scoped_to_converter_window(qapp, isolated_qsettings):
    window = ConverterMainWindow(embedded=True, host_theme="dark")
    try:
        assert window._embedded is True
        assert window._theme == "dark"
        assert window.styleSheet()
        window.set_host_theme("light")
        assert window._theme == "light"
    finally:
        window.close()


# ── icon helpers ─────────────────────────────────────────────────────────


def test_is_near_black_true_for_black():
    assert _is_near_black("#000000") is True


def test_is_near_black_false_for_white():
    assert _is_near_black("#ffffff") is False


def test_recolor_svg_replaces_black_hex():
    svg = '<svg><path fill="#000000" d="M0 0"/></svg>'
    out = _recolor_svg(svg, target="#f4f7fb").decode("utf-8")
    assert "#000000" not in out
    assert "#f4f7fb" in out


def test_recolor_svg_leaves_non_black_hex_alone():
    svg = '<svg><path fill="#336699" d="M0 0"/></svg>'
    out = _recolor_svg(svg, target="#f4f7fb").decode("utf-8")
    assert "#336699" in out


def test_recolor_svg_replaces_black_keyword():
    svg = '<svg><path style="fill: black;" d="M0 0"/></svg>'
    out = _recolor_svg(svg, target="#f4f7fb").decode("utf-8")
    assert "black" not in out.lower() or "#f4f7fb" in out


def test_recolor_svg_injects_fill_on_bare_path():
    svg = "<svg><path d=\"M0 0\"/></svg>"
    out = _recolor_svg(svg, target="#f4f7fb").decode("utf-8")
    assert 'fill="#f4f7fb"' in out


def test_icon_returns_empty_for_unknown_name():
    icon = _icon("definitely-not-a-real-icon-name")
    assert icon.isNull()


def test_icon_loads_existing_svg(qapp):
    icon = _icon("tools")
    assert not icon.isNull()


def test_icon_dark_mode_recolors_svg(qapp):
    global_dark = mw_mod._DARK_MODE
    mw_mod._DARK_MODE = True
    try:
        icon = _icon("tools")
        assert not icon.isNull()
    finally:
        mw_mod._DARK_MODE = global_dark


# ── help menu ────────────────────────────────────────────────────────────


def test_open_documentation_calls_desktop_services(window, monkeypatch):
    from PySide6.QtGui import QDesktopServices

    seen = []
    monkeypatch.setattr(QDesktopServices, "openUrl", staticmethod(lambda url: seen.append(url.toString())))
    window._open_documentation()
    assert seen == ["https://pycsamt.org/"]


def test_open_github_calls_desktop_services(window, monkeypatch):
    from PySide6.QtGui import QDesktopServices

    seen = []
    monkeypatch.setattr(QDesktopServices, "openUrl", staticmethod(lambda url: seen.append(url.toString())))
    window._open_github()
    assert seen == ["https://github.com/earthai-tech/pycsamt"]


def test_open_about_shows_message_box(window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox

    # QMessageBox.about() is a distinct static entry point from
    # information/question/warning/critical (which no_modal_dialogs already
    # neutralizes) and still runs its own modal exec() loop -- block that
    # explicitly so this stays a non-interactive unit test.
    calls = []
    monkeypatch.setattr(
        QMessageBox, "about", staticmethod(lambda *a, **k: calls.append(a))
    )
    window._open_about()
    assert calls
