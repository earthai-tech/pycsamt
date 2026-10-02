# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Menu icons must stay visible on the dark theme.

Edit ▸ Stations / Frequencies / Tensor, View ▸ PCSF 3-D / Airborne / Log
and Tools ▸ Solver Builder used to stay black on the dark theme: they were
missing from the re-icon list, and ``log.svg`` / ``frequency-editor.svg``
write ``rgb(0, 0, 0)``, which the recolourer did not handle.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtCore import QSize  # noqa: E402


def _luminance(icon) -> float:
    img = icon.pixmap(QSize(32, 32)).toImage()
    tot = lum = 0.0
    for y in range(img.height()):
        for x in range(img.width()):
            c = img.pixelColor(x, y)
            if c.alpha() > 100:
                tot += 1
                lum += 0.299 * c.red() + 0.587 * c.green() + 0.114 * c.blue()
    return lum / tot if tot else 255.0


def _menu_actions(menu, path=""):
    for act in menu.actions():
        label = f"{path} ▸ {act.text().replace('&', '')}"
        if act.menu() is not None:
            yield label, act
            yield from _menu_actions(act.menu(), label)
        elif not act.isSeparator():
            yield label, act


@pytest.fixture
def mw(qapp):
    from pycsamt.app.desktop.main_window import MainWindow

    w = MainWindow()
    yield w
    w._apply_theme("light")
    w.close()


def _dark_icons(mw):
    bar = list(mw.menuBar().actions())  # keep the list: see QMenu gotcha
    out = []
    for top in bar:
        if top.menu() is None:
            continue
        for label, act in _menu_actions(top.menu(), top.text()):
            if not act.icon().isNull():
                out.append((label, _luminance(act.icon())))
    return out


def test_recolor_handles_rgb_and_short_hex():
    from pycsamt.app.desktop.main_window import _recolor_svg

    out = _recolor_svg('<path style="fill: rgb(0, 0, 0); stroke: #111"/>'
                       '<path fill="rgba(10,10,10,0.5)"/>'
                       '<path fill="rgb(200, 30, 30)"/>').decode()
    assert "rgb(0" not in out and "#111" not in out and "rgba(10" not in out
    assert "rgb(200, 30, 30)" in out  # coloured parts are kept


def test_every_menu_icon_is_bright_on_dark(mw):
    mw._apply_theme("dark")
    icons = _dark_icons(mw)
    labels = " ".join(label for label, _l in icons)
    for needle in ("Stations", "Frequencies", "Tensor", "PCSF", "Airborne",
                   "Log", "Solver Builder", "Theme"):
        assert needle in labels, needle
    dim = [(label, round(lum)) for label, lum in icons if lum < 110]
    assert not dim, dim


def test_back_to_light_restores_dark_icons(mw):
    mw._apply_theme("dark")
    mw._apply_theme("light")
    icons = dict(_dark_icons(mw))
    stations = next(v for k, v in icons.items() if "Stations" in k)
    assert stations < 110  # the original (dark) artwork again


def test_stroke_only_svg_is_not_filled():
    """depth.svg inherits fill="none" from <svg>; injecting a fill turned
    it into a solid square on the dark theme."""
    from pycsamt.app.desktop.main_window import _recolor_svg

    out = _recolor_svg('<svg fill="none"><path d="M0 0L1 1" '
                       'stroke="#000000"/></svg>').decode()
    assert 'fill="#cdd6f4"' not in out and "#000000" not in out
    bare = _recolor_svg('<svg><path d="M0 0L1 1Z"/></svg>').decode()
    assert 'fill="#cdd6f4"' in bare  # default-black shapes still recoloured


def test_edit_menu_uses_its_own_icons(mw):
    from pycsamt.app.desktop.main_window import _ICON_NAMES

    names = {_ICON_NAMES.get(a.icon().cacheKey())
             for a in (mw._act_history, mw._act_find)}
    assert names == {"history", "find-station"}
    points = next(a for a in mw._edit_freq_menu.actions()
                  if "Point Editor" in a.text())
    assert _ICON_NAMES.get(points.icon().cacheKey()) == "point-editor"


def test_panel_window_icons_follow_the_theme(qapp):
    """Advanced ▸ Analyses / Utilities icons stayed black on dark."""
    from pycsamt.app.desktop.windows.advanced_window import (
        AdvancedToolsWindow,
    )

    w = AdvancedToolsWindow()
    try:
        items = [w._nav.item(i) for i in range(w._nav.count())
                 if not w._nav.item(i).icon().isNull()]
        assert items
        w.set_dark_mode(True)
        assert all(_luminance(it.icon()) > 110 for it in items)
        w.set_dark_mode(False)
        assert all(_luminance(it.icon()) < 110 for it in items)
    finally:
        w.close()
