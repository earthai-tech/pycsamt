# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SettingsPage — abstract base widget for one tab in APIConfigDialog.

Subclasses must implement:
    populate()  — read from PYCSAMT_* singletons → populate widgets
    collect()   — read widgets → return {"tab_key": {field: value, …}, …}
    reset()     — reset singletons then re-populate
"""

from __future__ import annotations

from PySide6.QtWidgets import QWidget


class SettingsPage(QWidget):
    """
    Abstract base for a settings tab page.

    Each concrete subclass owns one or more PYCSAMT_* singletons and
    exposes their fields as Qt widgets.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

    # ── Interface ─────────────────────────────────────────────────────────────

    def populate(self) -> None:
        """Read from the relevant PYCSAMT_* singletons and fill widgets."""

    def collect(self) -> dict:
        """
        Return a dict mapping controller apply-method keys to field dicts.

        Example::

            {
                "station": {"side": "top", "marker_symbol": "v"},
                "section": {"y_direction": "down"},
            }

        Return an empty dict (or sub-dicts) to signal no changes.
        """
        return {}

    def reset(self) -> None:
        """Reset the page's singletons to defaults, then re-populate."""
