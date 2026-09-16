# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.widgets.drop_zone
==========================================

A drag-and-drop-or-browse file/folder picker used at the top of every
conversion page. It only collects a path -- pages that care what the
path *is* (an Occam2D directory? a .pcsf file?) call
:func:`pycsamt.app.converter.jobs.detect_source_job` themselves and
show the result next to this widget, keeping ``DropZone`` reusable for
pages (PCGL/PCGS/PCPT builders) that don't need auto-detection.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


class DropZone(QWidget):
    """A labeled path field that accepts drag-and-drop or Browse.

    Parameters
    ----------
    label : str
        Field label, e.g. ``"Source"``.
    accept_dirs, accept_files : bool
        What kind of path this field accepts. When both are true, a
        drop of either kind is accepted, and Browse offers a file
        picker (folders can be typed/dropped, or picked with a
        secondary "Folder…" button).
    file_filter : str
        Qt file-dialog filter string, e.g. ``"EDI files (*.edi)"``.
    """

    pathChanged = Signal(str)

    def __init__(
        self,
        label: str = "Source",
        *,
        accept_dirs: bool = True,
        accept_files: bool = True,
        file_filter: str = "All files (*.*)",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._accept_dirs = accept_dirs
        self._accept_files = accept_files
        self._file_filter = file_filter
        self.setAcceptDrops(True)
        self.setObjectName("DropZone")
        # Plain QWidget subclasses don't paint a stylesheet background/
        # border by default (only widgets with their own paintEvent, like
        # QFrame, do) -- without this, DropZone's #DropZone QSS rule is
        # silently ignored and the "drop box" has no visible border at all.
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self._build_ui(label)

    def _build_ui(self, label: str) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(10, 8, 10, 10)
        root.setSpacing(6)

        caption = QLabel(label)
        caption.setObjectName("DropZoneLabel")
        root.addWidget(caption)

        row = QHBoxLayout()
        row.setSpacing(6)
        self._edit = QLineEdit()
        self._edit.setPlaceholderText("Drop a file or folder here, or Browse…")
        self._edit.textChanged.connect(self._on_text_changed)
        row.addWidget(self._edit)

        if self._accept_files:
            btn_file = QPushButton("Browse…")
            btn_file.clicked.connect(self._browse_file)
            row.addWidget(btn_file)
        if self._accept_dirs:
            btn_dir = QPushButton("Folder…")
            btn_dir.clicked.connect(self._browse_dir)
            row.addWidget(btn_dir)

        root.addLayout(row)

    # -- Qt drag/drop -----------------------------------------------------

    def dragEnterEvent(self, event) -> None:  # noqa: N802 - Qt override
        if event.mimeData().hasUrls():
            self._set_drag_active(True)
            event.acceptProposedAction()

    def dragLeaveEvent(self, event) -> None:  # noqa: N802 - Qt override
        self._set_drag_active(False)

    def dropEvent(self, event) -> None:  # noqa: N802 - Qt override
        self._set_drag_active(False)
        urls = event.mimeData().urls()
        if not urls:
            return
        path = Path(urls[0].toLocalFile())
        if path.is_dir() and not self._accept_dirs:
            return
        if path.is_file() and not self._accept_files:
            return
        self._edit.setText(str(path))
        event.acceptProposedAction()

    def _set_drag_active(self, active: bool) -> None:
        self.setProperty("dragActive", "true" if active else "false")
        self.style().unpolish(self)
        self.style().polish(self)

    # -- Browse -------------------------------------------------------------

    def _browse_file(self) -> None:
        # QFileDialog rejects a Path object outright (TypeError) -- always
        # pass str, or Browse crashes as soon as the field already holds a
        # value from a previous pick/drop.
        start = str(self.path() or Path.home())
        chosen, _ = QFileDialog.getOpenFileName(self, "Select file", start, self._file_filter)
        if chosen:
            self._edit.setText(chosen)

    def _browse_dir(self) -> None:
        start = str(self.path() or Path.home())
        chosen = QFileDialog.getExistingDirectory(self, "Select folder", start)
        if chosen:
            self._edit.setText(chosen)

    # -- API ------------------------------------------------------------

    def _on_text_changed(self, text: str) -> None:
        self.pathChanged.emit(text)

    def path(self) -> Path | None:
        text = self._edit.text().strip()
        return Path(text) if text else None

    def set_path(self, path: str | Path) -> None:
        self._edit.setText(str(path))

    def clear(self) -> None:
        self._edit.clear()
