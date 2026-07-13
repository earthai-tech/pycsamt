# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
APIConfigDialog — tabbed Settings dialog for all PYCSAMT_* API singletons.

Layout
------
┌──────────────────────────────────────────────┐
│  [Pseudosections] [View Controls] [Display]  │
│  [Topography]     [Interpretation]           │
├──────────────────────────────────────────────┤
│                                              │
│              <active page>                   │
│                                              │
├──────────────────────────────────────────────┤
│  [Reset Tab]  ·  [Apply]  [OK]  [Cancel]     │
└──────────────────────────────────────────────┘

Signals
-------
settings_changed(list[str])
    Emitted after Apply or OK with a list of tab keys that were written
    (e.g. ``["station", "section", "view_controls"]``).  MainWindow
    connects to this to refresh only the affected panels.

Deep-linking
------------
Call ``dialog.open_tab("pseudosections")`` before ``exec()`` to open the
dialog with a specific tab pre-selected.  Any panel in the app can use this
to let users jump straight to the relevant settings:

    from pycsamt.app.desktop.dialogs.settings_dialog import APIConfigDialog
    dlg = APIConfigDialog(ctrl, parent=self)
    dlg.open_tab("pseudosections")
    dlg.settings_changed.connect(handler)
    dlg.exec()
"""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.settings_controller import (
    SettingsController,
)
from pycsamt.app.desktop.widgets.settings_pages import (
    DisplayPage,
    InterpretationPage,
    PseudosectionsPage,
    TopographyPage,
    ViewControlsPage,
)

# Ordered tab registry: (tab_key, display_label, PageClass)
_TABS = [
    ("pseudosections", "Pseudosections", PseudosectionsPage),
    ("view_controls", "View Controls", ViewControlsPage),
    ("display", "Display", DisplayPage),
    ("topography", "Topography", TopographyPage),
    ("interpretation", "Interpretation", InterpretationPage),
]

# Map tab_key → index (built at module level for open_tab look-ups)
_TAB_INDEX: dict[str, int] = {key: i for i, (key, _, _) in enumerate(_TABS)}


class APIConfigDialog(QDialog):
    """
    Tabbed dialog for configuring all PYCSAMT_* API singletons.

    Parameters
    ----------
    ctrl : SettingsController
        Shared controller instance (owned by MainWindow).
    parent : QWidget, optional
    open_tab : str, optional
        If given, the dialog opens with this tab selected.
        Must be one of ``"pseudosections"``, ``"view_controls"``,
        ``"display"``, ``"topography"``, ``"interpretation"``.
    """

    settings_changed = Signal(list)  # list[str] of apply-method keys touched

    def __init__(
        self,
        ctrl: SettingsController,
        parent: QWidget | None = None,
        open_tab: str | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("API Configuration")
        self.setMinimumSize(560, 480)
        self._ctrl = ctrl
        self._snapshot = ctrl.snapshot()  # for Cancel
        self._pages: list = []
        self._build_ui()
        if open_tab:
            self.open_tab(open_tab)

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(8)

        # Tab widget
        self._tabs = QTabWidget()
        for _key, label, PageClass in _TABS:
            page = PageClass(self)
            self._pages.append(page)
            self._tabs.addTab(page, label)
        root.addWidget(self._tabs)

        # Button row
        btn_row = QHBoxLayout()

        self._reset_btn = QPushButton("Reset Tab")
        self._reset_btn.setToolTip(
            "Reset the current tab's settings to package defaults"
        )
        self._reset_btn.clicked.connect(self._on_reset_tab)
        btn_row.addWidget(self._reset_btn)
        btn_row.addStretch()

        box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Apply
            | QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
        )
        box.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(
            self._on_apply
        )
        box.accepted.connect(self._on_ok)
        box.rejected.connect(self._on_cancel)
        btn_row.addWidget(box)

        root.addLayout(btn_row)

    # ── Public API ────────────────────────────────────────────────────────────

    def open_tab(self, key: str) -> None:
        """Switch to the tab identified by *key* (e.g. ``"pseudosections"``)."""
        idx = _TAB_INDEX.get(key)
        if idx is not None:
            self._tabs.setCurrentIndex(idx)

    # ── Slots ─────────────────────────────────────────────────────────────────

    def _on_apply(self) -> None:
        """Apply all pages and emit settings_changed with the touched keys."""
        touched: list[str] = []
        for i, (_tab_key, _, _) in enumerate(_TABS):
            page = self._pages[i]
            try:
                kw_map = page.collect()
            except Exception:
                kw_map = {}
            for apply_key, fields in kw_map.items():
                if not fields:
                    continue
                method = getattr(self._ctrl, f"apply_{apply_key}", None)
                if method:
                    try:
                        method(**fields)
                        if apply_key not in touched:
                            touched.append(apply_key)
                    except Exception:
                        pass
        if touched:
            self.settings_changed.emit(touched)

    def _on_ok(self) -> None:
        self._on_apply()
        self.accept()

    def _on_cancel(self) -> None:
        """Restore all singletons from the pre-dialog snapshot."""
        try:
            self._ctrl.restore(self._snapshot)
        except Exception:
            pass
        self.reject()

    def _on_reset_tab(self) -> None:
        """Reset the active tab's singletons and re-populate its page."""
        idx = self._tabs.currentIndex()
        key = _TABS[idx][0]
        page = self._pages[idx]
        try:
            page.reset()
        except Exception:
            pass
        # Emit so live panels refresh
        self.settings_changed.emit([key])
