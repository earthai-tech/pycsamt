# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
LoadDataDialog — file-open dialog for survey data (EDI, EMTF-XML, AVG, J).

Supports drag-and-drop of files *and* folders (recurses, filters by
selected format), Browse Files / Browse Folder buttons, per-file
removal, and a Clear All action.  Returns accepted file paths via
``selected_paths``.

Everything that names the data follows the selected format
(:data:`FORMATS`): the drop zone ("Drop AVG files (.avg) or a folder
here"), its hover text, the "nothing found" message and the file-browser
title.  Dropping files of *another* supported format switches the format
to them instead of ignoring the drop.
"""

from __future__ import annotations

import os
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import (
    QDragEnterEvent,
    QDragLeaveEvent,
    QDropEvent,
)
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

_FORMAT_MAP = {
    "EDI": ("*.edi",),
    "EMTF XML": ("*.xml",),
    "AVG": ("*.avg",),
    "J / ModEM": ("*.j",),
    "All supported": ("*.edi", "*.xml", "*.avg", "*.j"),
}

# format -> (what the files are called, extensions shown to the user)
FORMATS = {
    "EDI": ("EDI files", ".edi"),
    "EMTF XML": ("EMTF-XML transfer functions", ".xml"),
    "AVG": ("AVG files", ".avg"),
    "J / ModEM": ("J files", ".j"),
    "All supported": ("EDI, EMTF-XML, AVG or J files",
                      ".edi · .xml · .avg · .j"),
}


def detect_format(paths) -> str | None:
    """The single format *paths* (files or folders) contain, if any."""
    found = set()
    single = {k: v for k, v in _FORMAT_MAP.items() if k != "All supported"}
    for p in paths:
        path = Path(p)
        files = ([x for x in path.rglob("*") if x.is_file()]
                 if path.is_dir() else [path])
        for f in files:
            suffix = f.suffix.lower()
            for name, globs in single.items():
                if any(suffix == g.lstrip("*").lower() for g in globs):
                    found.add(name)
    if not found:
        return None
    return found.pop() if len(found) == 1 else "All supported"


_EXT_FROM_GLOB = {
    glob.lstrip("*."): glob.lstrip("*")
    for globs in _FORMAT_MAP.values()
    for glob in globs
}


class _DropZone(QLabel):
    """Styled drop area — emits raw dropped paths for the dialog to process."""

    raw_paths_dropped = Signal(list)

    # EDI defaults; set_format() words them for the selected format
    _TEXT_IDLE = "⬇   Drop EDI files (.edi) or a folder here"
    _TEXT_HOVER = "  Release to add EDI files"

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("DropZone")
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setAcceptDrops(True)
        self.setMinimumHeight(88)
        self.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.setProperty("drag_over", "false")
        self.idle_text = self._TEXT_IDLE
        self.hover_text = self._TEXT_HOVER
        self.setText(self.idle_text)

    def set_format(self, fmt: str) -> None:
        """Word the drop zone for *fmt* (a key of :data:`FORMATS`)."""
        noun, exts = FORMATS.get(fmt, FORMATS["EDI"])
        self.idle_text = f"⬇   Drop {noun} ({exts}) or a folder here"
        self.hover_text = f"  Release to add {noun}"
        if self.property("drag_over") != "true":
            self.setText(self.idle_text)

    # ── Qt drag events ─────────────────────────────────────────────

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
            self._set_drag_over(True)

    def dragLeaveEvent(self, event: QDragLeaveEvent) -> None:
        self._set_drag_over(False)

    def dropEvent(self, event: QDropEvent) -> None:
        urls = event.mimeData().urls()
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
        self.setText(self.hover_text if state else self.idle_text)


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

    Signals
    -------
    open_format_studio_requested()
        Emitted when the user clicks "Open in Format Studio…" — this
        dialog only loads EDI/EMTF-XML/AVG/J survey data; anything else
        (PCSF/PCSM/PCBH/PCGL/PCGS/PCPT, or converting between formats)
        needs the full converter app. The dialog closes itself
        (rejected) so the caller can open that app without two modal
        windows competing.
    """

    open_format_studio_requested = Signal()

    def __init__(
        self,
        parent: QWidget | None = None,
        last_dir: str = "",
        recomputed_dir=None,
        existing_count: int = 0,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Open Survey Data")
        self.setMinimumSize(580, 460)
        self._last_dir = last_dir or str(Path.home())
        self._recomputed_dir = Path(recomputed_dir) if recomputed_dir else None
        self.selected_paths: list[str] = []
        self._existing_count = existing_count
        self._build_ui()

    # ── UI construction ────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(16, 14, 16, 12)
        root.setSpacing(10)

        self._mode_combo = QComboBox()
        self._mode_combo.addItem("Add to current survey", "append")
        self._mode_combo.addItem("Replace current survey", "replace")
        self._mode_combo.setCurrentIndex(0 if self._existing_count else 1)
        self._mode_combo.setVisible(self._existing_count > 0)
        self._mode_hint = QLabel()
        self._mode_hint.setWordWrap(True)
        self._mode_hint.setVisible(self._existing_count > 0)
        root.addWidget(self._mode_combo)
        root.addWidget(self._mode_hint)
        self._mode_combo.currentIndexChanged.connect(self._update_load_mode)

        # ── Format selector ───────────────────────────────────────
        fmt_row = QHBoxLayout()
        fmt_lbl = QLabel("Format:")
        fmt_lbl.setObjectName("DialogLabel")
        fmt_row.addWidget(fmt_lbl)
        self._fmt_combo = QComboBox()
        self._fmt_combo.addItems(list(_FORMAT_MAP.keys()))
        self._fmt_combo.setCurrentText("EDI")
        self._fmt_combo.setFixedWidth(170)
        self._fmt_combo.setToolTip(
            "Which files to pick up from drops and folders")
        fmt_row.addWidget(self._fmt_combo)
        self._fmt_hint = QLabel("")
        self._fmt_hint.setObjectName("InfoLabel")
        fmt_row.addWidget(self._fmt_hint)
        fmt_row.addStretch()
        root.addLayout(fmt_row)

        # ── Drop zone ─────────────────────────────────────────────
        self._drop_zone = _DropZone(self)
        self._drop_zone.raw_paths_dropped.connect(self._on_dropped)
        root.addWidget(self._drop_zone)
        self._fmt_combo.currentTextChanged.connect(self._on_format)
        self._on_format(self._fmt_combo.currentText())

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
                "Load all EDI files from the last recomputed output:\n"
                f"{self._recomputed_dir}"
            )
            btn_recomp.clicked.connect(self._load_recomputed)
            browse_row.addWidget(btn_recomp)

        browse_row.addStretch()
        btn_format_studio = QPushButton("Open in Format Studio…")
        btn_format_studio.setObjectName("BrowseButton")
        btn_format_studio.setToolTip(
            "Need PCSF/PCSM/PCBH/PCGL/PCGS/PCPT, or a format not listed "
            "above? Open the full pyCSAMT Format Studio."
        )
        btn_format_studio.clicked.connect(self._on_open_format_studio)
        browse_row.addWidget(btn_format_studio)
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
        self._update_load_mode()

    @property
    def load_mode(self) -> str:
        return self._mode_combo.currentData()

    def _update_load_mode(self) -> None:
        append = self.load_mode == "append"
        self._ok_btn.setText("Add Data" if append else "Load Data")
        self._mode_hint.setText(
            f"Keep the {self._existing_count} loaded stations "
            "and add new ones. "
            "Matching station IDs are skipped; existing data is kept."
            if append else
            f"Replace the {self._existing_count} loaded stations "
            "after the new files load successfully."
        )

    # ── Format ─────────────────────────────────────────────────────

    def _on_format(self, fmt: str) -> None:
        self._drop_zone.set_format(fmt)
        noun, exts = FORMATS.get(fmt, FORMATS["EDI"])
        self._fmt_hint.setText(exts)

    # ── Drag-and-drop handler ──────────────────────────────────────

    def _on_dropped(self, raw_paths: list[str]) -> None:
        """Expand folders, filter by current format, append to list.

        When nothing matches the selected format but the drop holds
        another supported format, switch to it (the user dropped AVG
        files while EDI was selected) rather than ignoring the drop.
        """
        found = self._matching(raw_paths)
        if not found and raw_paths:
            other = detect_format(raw_paths)
            if other and other != self._fmt_combo.currentText():
                self._fmt_combo.setCurrentText(other)
                found = self._matching(raw_paths)
        if found:
            self._add_paths(found)
        elif raw_paths:
            noun, exts = FORMATS.get(self._fmt_combo.currentText(),
                                     FORMATS["EDI"])
            # User dropped something but nothing matched — give visual cue
            self._drop_zone.setText(f"⚠  No {noun} ({exts}) found in the "
                                    "drop")

    def _matching(self, raw_paths: list[str]) -> list[str]:
        exts = _FORMAT_MAP[self._fmt_combo.currentText()]
        suffix_set = {e.lstrip("*.").lower() for e in exts}
        found: list[str] = []
        for p in raw_paths:
            path = Path(p)
            if path.is_dir():
                for ext in exts:
                    found.extend(str(x) for x in sorted(path.rglob(ext)))
            elif path.is_file():
                if path.suffix.lower().lstrip(".") in suffix_set:
                    found.append(str(path))
        return found

    # ── Browse slots ───────────────────────────────────────────────

    def _on_open_format_studio(self) -> None:
        self.open_format_studio_requested.emit()
        self.reject()

    def _load_recomputed(self) -> None:
        """Load all EDI files from the last EDIRecomputer output folder."""
        if not (self._recomputed_dir and self._recomputed_dir.is_dir()):
            return
        found = sorted(str(p) for p in self._recomputed_dir.rglob("*.edi"))
        if found:
            self._mode_combo.setCurrentIndex(1)
            self._set_paths(found)
        else:
            self._drop_zone.setText(
                "⚠  No EDI files found in recomputed folder"
            )

    def _browse_files(self) -> None:
        exts = " ".join(_FORMAT_MAP[self._fmt_combo.currentText()])
        noun, _e = FORMATS.get(self._fmt_combo.currentText(), FORMATS["EDI"])
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            f"Select {noun}",
            self._last_dir,
            f"{noun} ({exts});;All files (*)",
        )
        if paths:
            self._last_dir = str(Path(paths[0]).parent)
            self._add_paths(paths)

    def _browse_folder(self) -> None:
        folder = QFileDialog.getExistingDirectory(
            self, "Select survey folder", self._last_dir
        )
        if not folder:
            return
        self._last_dir = folder
        exts = _FORMAT_MAP[self._fmt_combo.currentText()]
        found: list[str] = []
        for ext in exts:
            found.extend(str(p) for p in Path(folder).rglob(ext))
        found.sort()
        self._add_paths(found)

    # ── File-list mutations ────────────────────────────────────────

    def _set_paths(self, paths: list[str]) -> None:
        """Replace the entire file list."""
        self._file_list.clear()
        for p in paths:
            self._file_list.addItem(p)
        self._refresh_ui()

    def _add_paths(self, paths: list[str]) -> None:
        """Append paths, skipping duplicates already in the list."""
        existing = {
            self._path_key(self._file_list.item(i).text())
            for i in range(self._file_list.count())
        }
        for p in paths:
            key = self._path_key(p)
            if key not in existing:
                self._file_list.addItem(p)
                existing.add(key)
        self._refresh_ui()

    @staticmethod
    def _path_key(path: str) -> str:
        return os.path.normcase(str(Path(path).resolve()))

    def _remove_selected(self) -> None:
        for item in list(self._file_list.selectedItems()):
            self._file_list.takeItem(self._file_list.row(item))
        self._refresh_ui()

    def _clear_all(self) -> None:
        self._file_list.clear()
        self._refresh_ui()
        self._drop_zone.setText(self._drop_zone.idle_text)

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
