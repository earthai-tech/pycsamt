# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
LicenseExpiredDialog -- the hard-block trial/license gate.

Product decision (PYCSAMT-DESKTOP-V2.6-MODERNIZATION-PLAN.md, Phase 11 /
Section 4): when no active trial or license is found, the app does not open
in a degraded/read-only mode -- it shows this modal dialog instead and, if
the user doesn't activate a key, quits entirely. Wiring this into the
startup sequence (calling it from the splash screen before ``MainWindow``
is shown) is Phase 12's job; this module only builds the dialog itself, so
Phase 12 can consume it as ``dlg.exec() == QDialog.DialogCode.Accepted``
== "proceed to MainWindow", anything else == "quit".
"""

from __future__ import annotations

from PySide6.QtCore import QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .. import branding
from ..licensing.interfaces import LicenseManager, LicenseStatus
from ..widgets.license_page import LicensePage


class LicenseExpiredDialog(QDialog):
    """Blocking gate shown when :meth:`LicenseManager.status` is not
    ``TRIAL_ACTIVE``/``LICENSED``.

    Embeds :class:`~pycsamt.app.desktop.widgets.license_page.LicensePage`
    for the actual key-entry/activation UI (no duplicated logic), adds a
    "Buy a license…" link and an explicit "Exit" button, and auto-accepts
    itself the moment activation succeeds
    (:attr:`LicensePage.activated <pycsamt.app.desktop.widgets.license_page.LicensePage.activated>`).
    """

    def __init__(
        self,
        license_manager: LicenseManager,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(f"{branding.APP_DISPLAY_NAME} — License required")
        self.setMinimumWidth(440)
        self.setModal(True)
        self._manager = license_manager
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(12)

        status = self._manager.status()
        if status == LicenseStatus.TRIAL_EXPIRED:
            headline = "Your free trial has ended."
        elif status == LicenseStatus.INVALID:
            headline = "Your license could not be verified."
        else:
            headline = "A license or active trial is required to continue."
        self._title_lbl = QLabel(f"<b>{headline}</b>")
        self._title_lbl.setWordWrap(True)
        root.addWidget(self._title_lbl)

        note = QLabel(
            "Enter a license key below to continue, or get one from the "
            f"link below. {branding.APP_DISPLAY_NAME} will close if you "
            "exit without activating."
        )
        note.setWordWrap(True)
        note.setObjectName("InfoLabel")
        root.addWidget(note)

        self._license_page = LicensePage(license_manager=self._manager)
        self._license_page.activated.connect(self.accept)
        root.addWidget(self._license_page)

        btn_row = QHBoxLayout()
        self._buy_btn = QPushButton("Buy a license…")
        self._buy_btn.clicked.connect(self._open_buy_page)
        btn_row.addWidget(self._buy_btn)
        btn_row.addStretch()
        self._exit_btn = QPushButton("Exit")
        self._exit_btn.clicked.connect(self.reject)
        btn_row.addWidget(self._exit_btn)
        root.addLayout(btn_row)

    def _open_buy_page(self) -> None:
        QDesktopServices.openUrl(QUrl(branding.URL_DOCS))
