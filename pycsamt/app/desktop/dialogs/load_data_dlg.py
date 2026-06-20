# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
LoadDataDialog — file-open dialog for EDI / AVG / J survey data.

Supports drag-and-drop of files *and* folders (recurses, filters by
selected format), Browse Files / Browse Folder buttons, per-file
removal, and a Clear All action.  Returns accepted file paths via
``selected_paths``.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QDragEnterEvent, QDragLeaveEvent, QDropEvent
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

_FORMAT_MAP = {
    "EDI":          ("*.edi",),
    "AVG":          ("*.avg",),
    "J / ModEM":    ("*.j",),
    "All supported": ("*.edi", "*.avg", "*.j"),
}

_EXT_FROM_GLOB = {
    glob.lstrip("*."): glob.lstrip("*")
    for globs in _FORMAT_MAP.values()
    for glob in globs
}


class _DropZone(QLabel):
    """Styled drop area — emits raw dropped paths for the dialog to process."""

    raw_paths_dropped = Signal(list)

    _TEXT_IDLE  = "⬇   Drop EDI files or a folder here"
    _TEXT_HOVER = "  Release to add files"

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("DropZone")
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setAcceptDrops(True)
        self.setMinimumHeight(88)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setProperty("drag_over", "false")
        self.setText(self._TEXT_IDLE)

    # ── Qt drag events ─────────────────────────────────────────────

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
            self._set_drag_over(True)

    def dragLeaveEvent(self, event: QDragLeaveEvent) -> None:
        self._set_drag_over(False)

    def dropEvent(self, event: QDropEvent) -> None:
        urls  = event.mimeData().urls()
        paths = [u.toLocalFile() for u in urls if u.isLocalFile()]
        if paths:
            self.raw_paths_dropped.emit(paths)
        self._set_drag_over(False)
        event.acceptProposedAction()

    # ── Helpers ────────────────────────────────────────────────────

    def _set_drag_over(self, state: bool) -> None:
        self.setProperty("drag_over", "true" if state else "false")
        self.style().unpolish(self)
        self.style().polish(self)
        self.setText(self._TEXT_HOVER if state else self._TEXT_IDLE)


class LoadDataDialog(QDialog):
    """
    Modal dialog for selecting survey data files.

    After ``exec()`` returns ``QDialog.Accepted``, read ``selected_paths``
    for the confirmed file paths.

    Parameters
    ----------
    recomputed_dir : Path or str, optional
        If provided and the directory exists, a *Load Recomputed EDIs* button
        is shown so the user can instantly load the output of the last
        EDIRecomputer run without navigating manually.
    """

    def __init__(
        self,
        parent: QWidget | None = None,
        last_dir: str = "",
        recomputed_dir=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Open Survey Data")
        self.setMinimumSize(580, 460)
        self._last_dir = last_dir or str(Path.home())
        self._recomputed_dir = Path(recomputed_dir) if recomputed_dir else None
        self.selected_paths: List[str] = []
        self._build_ui()

    # ── UI construction ────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(16, 14, 16, 12)
        root.setSpacing(10)

        # ── Format selector ───────────────────────────────────────
        fmt_row = QHBoxLayout()
        fmt_lbl = QLabel("Format:")
        fmt_lbl.setObjectName("DialogLabel")
        fmt_row.addWidget(fmt_lbl)
        self._fmt_combo = QComboBox()
        self._fmt_combo.addItems(list(_FORMAT_MAP.keys()))
        self._fmt_combo.setCurrentText("EDI")
        self._fmt_combo.setFixedWidth(170)
        fmt_row.addWidget(self._fmt_combo)
        fmt_row.addStretch()
        root.addLayout(fmt_row)

        # ── Drop zone ─────────────────────────────────────────────
        self._drop_zone = _DropZone(self)
        self._drop_zone.raw_paths_dropped.connect(self._on_dropped)
        root.addWidget(self._drop_zone)

        # ── Browse buttons ────────────────────────────────────────
        browse_row = QHBoxLayout()
        browse_row.setSpacing(8)
        btn_files = QPushButton("Browse Files…")
        btn_files.setObjectName("BrowseButton")
        btn_files.clicked.connect(self._browse_files)
        btn_folder = QPushButton("Browse Folder…")
        btn_folder.setObjectName("BrowseButton")
        btn_folder.clicked.connect(self._browse_folder)
        browse_row.addWidget(btn_files)
        browse_row.addWidget(btn_folder)

        # Show shortcut only when a previous recompute output folder exists.
        if self._recomputed_dir and self._recomputed_dir.is_dir():
            btn_recomp = QPushButton("◈  Load Recomputed EDIs")
            btn_recomp.setObjectName("BrowseButton")
            btn_recomp.setToolTip(
                f"Load all EDI files from the last recomputed output:\n{self._recomputed_dir}"
            )
            btn_recomp.clicked.connect(self._load_recomputed)
            browse_row.addWidget(btn_recomp)

        browse_row.addStretch()
        root.addLayout(browse_row)

        # ── Separator ─────────────────────────────────────────────
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setObjectName("Separator")
        root.addWidget(sep)

        # ── File list header: count + action buttons ──────────────
        list_hdr = QHBoxLayout()
        list_hdr.setSpacing(6)
        self._count_lbl = QLabel("Selected files: 0")
        self._count_lbl.setObjectName("FileCountLabel")
        list_hdr.addWidget(self._count_lbl)
        list_hdr.addStretch()

        self._btn_remove = QPushButton("Remove Selected")
        self._btn_remove.setObjectName("FileListBtn")
        self._btn_remove.setEnabled(False)
        self._btn_remove.clicked.connect(self._remove_selected)

        self._btn_clear = QPushButton("Clear All")
        self._btn_clear.setObjectName("FileListBtn")
        self._btn_clear.setEnabled(False)
        self._btn_clear.clicked.connect(self._clear_all)

        list_hdr.addWidget(self._btn_remove)
        list_hdr.addWidget(self._btn_clear)
        root.addLayout(list_hdr)

        # ── File list ─────────────────────────────────────────────
        self._file_list = QListWidget()
        self._file_list.setObjectName("FileList")
        self._file_list.setSelectionMode(
            QListWidget.SelectionMode.ExtendedSelection
        )
        self._file_list.setMinimumHeight(150)
        self._file_list.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        self._file_list.itemSelectionChanged.connect(self._on_sel_changed)
        root.addWidget(self._file_list)

        # ── OK / Cancel ───────────────────────────────────────────
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
        )
        self._ok_btn = buttons.button(QDialogButtonBox.StandardButton.Ok)
        self._ok_btn.setText("Load Data")
        self._ok_btn.setEnabled(False)
        buttons.accepted.connect(self._on_accepted)
        buttons.rejected.connect(self.reject)
        root.addWidget(buttons)

    # ── Drag-and-drop handler ──────────────────────────────────────

    def _on_dropped(self, raw_paths: List[str]) -> None:
        """Expand folders, filter by current format, append to list."""
        exts = _FORMAT_MAP[self._fmt_combo.currentText()]
        suffix_set = {e.lstrip("*.").lower() for e in exts}
        found: List[str] = []
        for p in raw_paths:
            path = Path(p)
            if path.is_dir():
                for ext in exts:
                    found.extend(str(x) for x in sorted(path.rglob(ext)))
            elif path.is_file():
                if path.suffix.lower().lstrip(".") in suffix_set:
                    found.append(str(path))
        if found:
            self._add_paths(found)
        elif raw_paths:
            # User dropped something but nothing matched — give visual cue
            self._drop_zone.setText(
                f"⚠  No {self._fmt_combo.currentText()} files found in drop"
            )

    # ── Browse slots ───────────────────────────────────────────────

    def _load_recomputed(self) -> None:
        """Load all EDI files from the last EDIRecomputer output folder."""
        if not (self._recomputed_dir and self._recomputed_dir.is_dir()):
            return
        found = sorted(str(p) for p in self._recomputed_dir.rglob("*.edi"))
        if found:
            self._set_paths(found)
        else:
            self._drop_zone.setText("⚠  No EDI files found in recomputed folder")

    def _browse_files(self) -> None:
        exts     = " ".join(_FORMAT_MAP[self._fmt_combo.currentText()])
        fmt_name = self._fmt_combo.currentText()
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            f"Select {fmt_name} files",
            self._last_dir,
            f"{fmt_name} files ({exts});;All files (*)",
        )
        if paths:
            self._last_dir = str(Path(paths[0]).parent)
            self._set_paths(paths)

    def _browse_folder(self) -> None:
        folder = QFileDialog.getExistingDirectory(
            self, "Select survey folder", self._last_dir
        )
        if not folder:
            return
        self._last_dir = folder
        exts  = _FORMAT_MAP[self._fmt_combo.currentText()]
        found: List[str] = []
        for ext in exts:
            found.extend(str(p) for p in Path(folder).rglob(ext))
        found.sort()
        self._set_paths(found)

    # ── File-list mutations ────────────────────────────────────────

    def _set_paths(self, paths: List[str]) -> None:
        """Replace the entire file list."""
        self._file_list.clear()
        for p in paths:
            self._file_list.addItem(p)
        self._refresh_ui()

    def _add_paths(self, paths: List[str]) -> None:
        """Append paths, skipping duplicates already in the list."""
        existing = {
            self._file_list.item(i).text()
            for i in range(self._file_list.count())
        }
        for p in paths:
            if p not in existing:
                self._file_list.addItem(p)
                existing.add(p)
        self._refresh_ui()

    def _remove_selected(self) -> None:
        for item in list(self._file_list.selectedItems()):
            self._file_list.takeItem(self._file_list.row(item))
        self._refresh_ui()

    def _clear_all(self) -> None:
        self._file_list.clear()
        self._refresh_ui()
        self._drop_zone.setText(_DropZone._TEXT_IDLE)

    # ── State refresh ──────────────────────────────────────────────

    def _refresh_ui(self) -> None:
        n = self._file_list.count()
        self._count_lbl.setText(f"Selected files: {n}")
        self._btn_clear.setEnabled(n > 0)
        self._ok_btn.setEnabled(n > 0)
        self._on_sel_changed()

    def _on_sel_changed(self) -> None:
        self._btn_remove.setEnabled(bool(self._file_list.selectedItems()))

    # ── Accept ────────────────────────────────────────────────────

    def _on_accepted(self) -> None:
        self.selected_paths = [
            self._file_list.item(i).text()
            for i in range(self._file_list.count())
        ]
        self.accept()
