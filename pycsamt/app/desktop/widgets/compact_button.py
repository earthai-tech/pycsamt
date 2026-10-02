# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
compact_button — make small glyph buttons (+, −, ↑, ↓, …, 📂, ▶) legible.

Both theme stylesheets give every ``QPushButton`` ``padding: 5px 14px``,
i.e. 28 px of horizontal padding. A button fixed at 26-36 px wide is then
left with zero (or negative) room for its label, so Qt clips the glyph
away and the button renders as an empty rounded square.

:func:`compact_button` drops the horizontal padding for one button only,
keeping every other theme rule (colours, border, radius, hover), and fixes
its size. Use it for every icon-sized push button instead of a bare
``setFixedWidth``.
"""

from __future__ import annotations

from PySide6.QtWidgets import QPushButton

# Local rule, not a theme selector: it must win over the app stylesheet's
# generic ``QPushButton { padding: 5px 14px; }`` regardless of which theme
# (or none) is loaded. The size is pinned *in the rule* too: QSS
# min-/max-width are applied at polish time, after (and over) any
# setFixedWidth() call, so a bare "min-width: 0" let layouts squeeze the
# button to ~10 px.
_COMPACT_QSS = (
    "QPushButton {{ padding: 0px; min-width: {w}px; max-width: {w}px;{h} }}"
)


def compact_button(
    btn: QPushButton, width: int = 28, height: int | None = None
) -> QPushButton:
    """Size *btn* to ``width`` (× ``height``) with its glyph fully visible."""
    btn.setProperty("compact", True)
    h = (
        ""
        if height is None
        else f" min-height: {height}px; max-height: {height}px;"
    )
    btn.setStyleSheet(_COMPACT_QSS.format(w=width, h=h))
    if height is None:
        btn.setFixedWidth(width)
    else:
        btn.setFixedSize(width, height)
    return btn


__all__ = ["compact_button"]
