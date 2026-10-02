"""Reusable native empty-state view for unavailable analytical results."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)


class UnavailableResultView(QWidget):
    """Centered explanation shown when a result cannot be produced."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("UnavailableResultView")
        self.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )

        outer = QVBoxLayout(self)
        outer.setContentsMargins(28, 28, 28, 28)
        outer.addStretch(1)

        row = QHBoxLayout()
        row.addStretch(1)
        card = QFrame()
        card.setObjectName("UnavailableResultCard")
        card.setMaximumWidth(660)
        card_layout = QVBoxLayout(card)
        card_layout.setContentsMargins(34, 30, 34, 30)
        card_layout.setSpacing(12)

        self._badge = QLabel("RESULT UNAVAILABLE")
        self._badge.setObjectName("UnavailableResultBadge")
        self._badge.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._badge.setFixedHeight(28)
        self._badge.setMaximumWidth(190)
        card_layout.addWidget(
            self._badge, alignment=Qt.AlignmentFlag.AlignHCenter
        )

        self._title = QLabel()
        self._title.setObjectName("UnavailableResultTitle")
        self._title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._title.setWordWrap(True)
        card_layout.addWidget(self._title)

        self._reason = QLabel()
        self._reason.setObjectName("UnavailableResultReason")
        self._reason.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._reason.setWordWrap(True)
        card_layout.addWidget(self._reason)

        self._guidance = QLabel()
        self._guidance.setObjectName("UnavailableResultGuidance")
        self._guidance.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._guidance.setWordWrap(True)
        card_layout.addWidget(self._guidance)

        row.addWidget(card, 1)
        row.addStretch(1)
        outer.addLayout(row)
        outer.addStretch(1)

        self.set_content(
            "Load survey data to begin",
            "No stations are currently available to the QC dashboard.",
            "Load EDI or EMTF-XML files, then select a diagnostic.",
        )

    def set_content(self, title: str, reason: str, guidance: str = "") -> None:
        """Replace the explanatory copy without rebuilding the template."""
        self._title.setText(title)
        self._reason.setText(reason)
        self._guidance.setText(guidance)
        self._guidance.setVisible(bool(guidance.strip()))


__all__ = ["UnavailableResultView"]
