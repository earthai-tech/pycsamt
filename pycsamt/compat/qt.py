# File: pycsamt/compat/qt.py
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

"""
Compatibility helpers for Qt/PySide6 fragility that varies by
Python version.

Constructing many top-level Qt widgets (each often holding a
matplotlib ``Figure``/``FigureCanvasQTAgg``) inside one long-lived
offscreen ``QApplication`` leaves cyclic reference garbage that
CPython's generational collector reclaims at an unpredictable point
-- *any* later allocation, anywhere, can trigger a collection pass.
Reclaiming a widget then runs its C++/Shiboken destructor, and on
Python 3.9 that destructor call has been observed to corrupt the
interpreter and segfault, with the crash surfacing far from the
widget that actually went uncollected (see
``pycsamt/app/tests/conftest.py``'s ``_stub_mpl_qt_toolbar`` for the
sibling issue this complements).

:func:`settle_qt_gc` forces that collection to happen at a
known-safe point -- right after a widget is closed and discarded --
instead of leaving it to fire mid-construction of some unrelated
widget later on.
"""

from __future__ import annotations

import gc

__all__ = ["settle_qt_gc"]


def settle_qt_gc() -> None:
    """Deterministically reclaim cyclic Qt/Shiboken garbage.

    Call this right after a top-level Qt widget is closed and its
    last Python reference dropped, so any pending destructor calls
    happen now, at this known location, rather than during an
    unrelated widget's construction later in the process.
    """
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        gc.collect()
        return

    app = QApplication.instance()
    if app is not None:
        app.processEvents()
    gc.collect()
    if app is not None:
        app.processEvents()
