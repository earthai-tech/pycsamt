# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages.batch_page
==========================================

A queue of mixed conversion jobs -- inversion results / AI arrays /
.pcsf / .pcsm (-> PCSF or PCSM) and .edi / EMTF-XML files (converted
each to the other) -- auto-classified on add and run sequentially with
a per-row status, via :func:`pycsamt.app.converter.jobs.run_batch_job`.

PCBH/PCGL/PCGS/PCPT builders are not queued here: each needs its own
structured per-file options (collar coordinates, planar/linear/fault
groupings, ...) that don't fit one generic row, so they stay on their
own dedicated pages.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
)

from pycsamt.app.converter import jobs
from pycsamt.app.converter.pages._base import ConverterPage, mark_primary, page_root_layout
from pycsamt.app.converter.widgets import DropZone

_COLUMNS = ("Path", "Kind", "Status")


class BatchPage(ConverterPage):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._rows: list[Path] = []
        self._build_ui()

    def _build_ui(self) -> None:
        root = page_root_layout(self)

        add_row = QHBoxLayout()
        add_files_btn = QPushButton("Add files…")
        add_files_btn.clicked.connect(self._add_files)
        add_row.addWidget(add_files_btn)

        add_dir_btn = QPushButton("Add folder…")
        add_dir_btn.clicked.connect(self._add_folder)
        add_row.addWidget(add_dir_btn)

        remove_btn = QPushButton("Remove selected")
        remove_btn.clicked.connect(self._remove_selected)
        add_row.addWidget(remove_btn)

        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self._clear)
        add_row.addWidget(clear_btn)
        add_row.addStretch()
        root.addLayout(add_row)

        self._table = QTableWidget(0, len(_COLUMNS))
        self._table.setHorizontalHeaderLabels(_COLUMNS)
        self._table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setAlternatingRowColors(True)
        self._table.verticalHeader().setVisible(False)
        root.addWidget(self._table)

        self._output_dir = DropZone(
            "Output folder", accept_files=False, accept_dirs=True
        )
        root.addWidget(self._output_dir)

        opts = QHBoxLayout()
        opts.addWidget(QLabel("Inversion target format:"))
        self._to_format = QComboBox()
        self._to_format.addItems(["pcsf", "pcsm"])
        opts.addWidget(self._to_format)
        self._overwrite = QCheckBox("Overwrite existing output")
        opts.addWidget(self._overwrite)
        opts.addStretch()
        root.addLayout(opts)

        self._run_button = QPushButton("Run all")
        self._run_button.clicked.connect(self._on_run_all)
        mark_primary(self._run_button)
        root.addWidget(self._run_button)

    # -- queue management -------------------------------------------------

    def _add_paths(self, paths: list[Path]) -> None:
        for path in paths:
            if path in self._rows:
                continue
            self._rows.append(path)
            row = self._table.rowCount()
            self._table.insertRow(row)
            self._table.setItem(row, 0, QTableWidgetItem(str(path)))
            kind = jobs.classify_batch_item(path)
            self._table.setItem(row, 1, QTableWidgetItem(kind))
            self._table.setItem(row, 2, QTableWidgetItem("pending"))

    def _add_files(self) -> None:
        files, _ = QFileDialog.getOpenFileNames(self, "Add files")
        self._add_paths([Path(f) for f in files])

    def _add_folder(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Add folder")
        if not folder:
            return
        root = Path(folder)
        candidates = sorted(root.glob("*.edi")) + sorted(root.glob("*.xml"))
        if not candidates:
            candidates = [root]
        self._add_paths(candidates)

    def _remove_selected(self) -> None:
        rows = sorted({idx.row() for idx in self._table.selectedIndexes()}, reverse=True)
        for row in rows:
            del self._rows[row]
            self._table.removeRow(row)

    def _clear(self) -> None:
        self._rows.clear()
        self._table.setRowCount(0)

    # -- run ----------------------------------------------------------

    def _on_run_all(self) -> None:
        if not self._rows:
            self.log.clear()
            self.log.append("Add at least one file first.")
            return
        output_dir = self._output_dir.path()
        if output_dir is None:
            self.log.clear()
            self.log.append("Pick an output folder first.")
            return

        for row in range(self._table.rowCount()):
            self._table.item(row, 2).setText("pending")

        items = [
            {"path": path, "kind": self._table.item(row, 1).text()}
            for row, path in enumerate(self._rows)
        ]

        def _on_success(results: list[dict]) -> None:
            for row, result in enumerate(results):
                self._table.item(row, 2).setText(result["status"])
            n_done = sum(1 for r in results if r["status"] == "done")
            n_error = sum(1 for r in results if r["status"] == "error")
            self.log.append(f"{n_done} done, {n_error} failed, {len(results) - n_done - n_error} skipped.")

        self.run_job(
            jobs.run_batch_job,
            items,
            output_dir,
            to_format=self._to_format.currentText(),
            overwrite=self._overwrite.isChecked(),
            supports_progress=True,
            on_success=_on_success,
            busy_message=f"Running {len(items)} job(s)…",
            success_message="Batch complete.",
        )
