# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
NoDataDialog — shown whenever a tool or panel needs survey data but none is loaded.

Usage (from MainWindow or any panel window)
-------------------------------------------
    from pycsamt.app.desktop.dialogs.no_data_dialog import NoDataDialog

    # Returns True if the user clicked "Load EDI files…" (caller should then
    # trigger the file-open flow).  Returns False on Cancel.
    if NoDataDialog.require(parent=self, tool_name="Phase Tensor Map"):
        self._act_open.trigger()

Or instantiate directly if you need finer control:

    dlg = NoDataDialog("Strike Analyzer", parent=self)
    if dlg.exec() == QDialog.DialogCode.Accepted:
        ...
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QStyle,
    QVBoxLayout,
    QWidget,
)


class NoDataDialog(QDialog):
    """
    Modal dialog explaining that survey / EDI data must be loaded first.

    Parameters
    ----------
    tool_name : str
        Name of the tool or panel that was requested (shown in the message).
    parent : QWidget, optional
    """

    def __init__(
        self, tool_name: str = "", parent: QWidget | None = None
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("No Survey Data Loaded")
        self.setFixedWidth(420)
        self.setSizePolicy(
            self.sizePolicy().horizontalPolicy(),
            self.sizePolicy().verticalPolicy(),
        )
        self._tool_name = tool_name
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(14)
        root.setContentsMargins(20, 20, 20, 16)

        # ── Icon + text row ───────────────────────────────────────────────
        row = QHBoxLayout()
        row.setSpacing(16)

        icon_lbl = QLabel()
        icon_lbl.setFixedSize(48, 48)
        pix = (
            self.style()
            .standardIcon(QStyle.StandardPixmap.SP_MessageBoxWarning)
            .pixmap(48, 48)
        )
        icon_lbl.setPixmap(pix)
        icon_lbl.setAlignment(Qt.AlignmentFlag.AlignTop)
        row.addWidget(icon_lbl)

        text_col = QVBoxLayout()
        text_col.setSpacing(6)

        title = QLabel("<b>No survey data loaded</b>")
        title.setTextFormat(Qt.TextFormat.RichText)
        text_col.addWidget(title)

        if self._tool_name:
            body_text = (
                f"<b>{self._tool_name}</b> requires survey data "
                f"(EDI files or a compatible format) to be loaded first."
            )
        else:
            body_text = (
                "This action requires survey data "
                "(EDI files or a compatible format) to be loaded first."
            )
        body = QLabel(body_text + "<br><br>Would you like to load files now?")
        body.setTextFormat(Qt.TextFormat.RichText)
        body.setWordWrap(True)
        text_col.addWidget(body)

        row.addLayout(text_col, stretch=1)
        root.addLayout(row)

        # ── Buttons ───────────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        btn_row.addStretch()

        self._load_btn = QPushButton("Load EDI files…")
        self._load_btn.setDefault(True)
        self._load_btn.setMinimumWidth(130)
        self._load_btn.clicked.connect(self.accept)

        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)

        btn_row.addWidget(self._load_btn)
        btn_row.addWidget(cancel_btn)
        root.addLayout(btn_row)

    # ── Convenience class-method ──────────────────────────────────────────────

    @classmethod
    def require(cls, parent: QWidget | None, tool_name: str = "") -> bool:
        """
        Show the dialog; return ``True`` if the user clicked *Load EDI files*.

        Typical caller pattern in MainWindow::

            if NoDataDialog.require(self, "Strike Analyzer"):
                self._act_open.trigger()   # open the file-load dialog
        """
        return cls(tool_name, parent).exec() == QDialog.DialogCode.Accepted
