# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
LicensePage — license activation / trial status widget.

Embedded as a tab in ``PreferencesDialog`` (the app's actual "Settings"
dialog — not ``APIConfigDialog``'s ``settings_pages/``, which is scoped
to plot-rendering/display preferences, an unrelated concern) and reused
by ``AboutDialog``'s compact status panel.

Depends on the ``LicenseManager`` Protocol
(``pycsamt.app.desktop.licensing.interfaces``), not a concrete class.
``main_window.py`` passes the real
``pycsamt.app.desktop.licensing.manager.get_default_manager()`` in; callers
that omit ``license_manager`` (tests, ``PreferencesDialog``'s own default)
fall back to ``NullLicenseManager`` -- an unconditional, non-persisting
trial, never the real license/trial state.
"""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.licensing.interfaces import LicenseManager, LicenseStatus
from pycsamt.app.desktop.licensing.null_manager import NullLicenseManager

_STATUS_TEXT = {
    LicenseStatus.TRIAL_ACTIVE: "Trial active",
    LicenseStatus.TRIAL_EXPIRED: "Trial expired",
    LicenseStatus.LICENSED: "Licensed",
    LicenseStatus.INVALID: "Invalid license",
    LicenseStatus.UNKNOWN: "Unknown",
}

_STATUS_COLOR = {
    LicenseStatus.TRIAL_ACTIVE: "#3b82f6",
    LicenseStatus.TRIAL_EXPIRED: "#ef4444",
    LicenseStatus.LICENSED: "#22c55e",
    LicenseStatus.INVALID: "#ef4444",
    LicenseStatus.UNKNOWN: "#6b7280",
}


def status_badge_html(status: LicenseStatus) -> str:
    """A small colour-coded HTML badge for *status* -- shared by
    ``LicensePage`` and ``AboutDialog``'s compact panel so the two never
    drift into inconsistent wording/colour for the same status."""
    text = _STATUS_TEXT.get(status, "Unknown")
    color = _STATUS_COLOR.get(status, "#6b7280")
    return (
        f"<span style='background:{color}; color:white; padding:2px 8px; "
        f"border-radius:8px; font-size:11px; font-weight:bold;'>{text}</span>"
    )


class LicensePage(QWidget):
    """License activation / trial-status tab.

    Parameters
    ----------
    license_manager : LicenseManager, optional
        Defaults to :class:`NullLicenseManager` (always an unlimited
        trial, rejects every key) when the caller doesn't pass one.
        ``main_window.py`` passes
        ``licensing.manager.get_default_manager()`` for the real,
        persisted license/trial state.
    """

    #: Emitted after :meth:`_on_activate` results in ``LicenseStatus.LICENSED``.
    #: The trial-expiry gate dialog connects to this to auto-close itself
    #: once a real key is accepted, instead of polling ``refresh()`` state.
    activated = Signal()

    def __init__(
        self,
        license_manager: LicenseManager | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._manager: LicenseManager = license_manager or NullLicenseManager()
        self._build_ui()
        self.refresh()

    # ── UI ────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(12)

        status_row = QHBoxLayout()
        status_row.addWidget(QLabel("Status:"))
        self._status_lbl = QLabel()
        self._status_lbl.setTextFormat(self._status_lbl.textFormat().RichText)
        status_row.addWidget(self._status_lbl)
        status_row.addStretch()
        root.addLayout(status_row)

        self._trial_lbl = QLabel()
        self._trial_lbl.setWordWrap(True)
        self._trial_lbl.setObjectName("InfoLabel")
        root.addWidget(self._trial_lbl)

        key_row = QHBoxLayout()
        key_row.addWidget(QLabel("License key:"))
        self._key_edit = QLineEdit()
        self._key_edit.setPlaceholderText("Enter your license key…")
        key_row.addWidget(self._key_edit)
        self._btn_activate = QPushButton("Activate")
        self._btn_activate.clicked.connect(self._on_activate)
        key_row.addWidget(self._btn_activate)
        root.addLayout(key_row)

        self._activate_msg = QLabel("")
        self._activate_msg.setWordWrap(True)
        root.addWidget(self._activate_msg)

        root.addStretch()

    # ── Behaviour ─────────────────────────────────────────────────────

    def refresh(self) -> None:
        """Re-read the license manager and update every widget."""
        status = self._manager.status()
        self._status_lbl.setText(status_badge_html(status))

        trial = self._manager.trial_state()
        if trial is not None:
            if trial.is_expired:
                self._trial_lbl.setText("Your trial has ended.")
            else:
                self._trial_lbl.setText(
                    f"{trial.days_remaining} day(s) remaining in your trial."
                )
        elif status == LicenseStatus.LICENSED:
            self._trial_lbl.setText("Thank you for licensing pyCSAMT.")
        else:
            self._trial_lbl.setText("")

        self._btn_activate.setEnabled(status != LicenseStatus.LICENSED)
        self._key_edit.setEnabled(status != LicenseStatus.LICENSED)

    def _on_activate(self) -> None:
        key = self._key_edit.text().strip()
        if not key:
            self._activate_msg.setText("Enter a license key first.")
            return
        result = self._manager.activate(key)
        if result == LicenseStatus.LICENSED:
            self._activate_msg.setText("License activated — thank you!")
            self._key_edit.clear()
            self.activated.emit()
        else:
            self._activate_msg.setText(
                "That key could not be activated. Check it and try again."
            )
        self.refresh()
