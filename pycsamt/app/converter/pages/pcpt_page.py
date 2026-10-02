# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages.pcpt_page
=========================================

Build a PCPT (targets / points of interest, ``.pcpt.json``) document
from a CSV or XLSX table -- the GUI face of ``pycsamt format
build-pcpt`` (:mod:`pycsamt.format.pointset`).
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QPushButton,
)

from pycsamt.app.converter import jobs
from pycsamt.app.converter.pages._base import ConverterPage, mark_primary, page_root_layout
from pycsamt.app.converter.widgets import DropZone


class PcptPage(ConverterPage):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        root = page_root_layout(self)

        hint = QLabel("CSV/XLSX columns: name, x, y[, z] at minimum.")
        hint.setObjectName("PageHint")
        hint.setWordWrap(True)
        root.addWidget(hint)

        self._source = DropZone(
            "Points CSV/XLSX", accept_dirs=False,
            file_filter="CSV/XLSX files (*.csv *.xlsx *.xlsm)",
        )
        root.addWidget(self._source)

        self._output = DropZone("Output .pcpt.json", accept_dirs=False)
        root.addWidget(self._output)

        form = QFormLayout()
        self._sheet = QLineEdit()
        self._sheet.setPlaceholderText("XLSX only: sheet name/index")
        form.addRow("Sheet:", self._sheet)
        self._crs = QLineEdit()
        self._crs.setPlaceholderText("e.g. EPSG:4326")
        form.addRow("CRS:", self._crs)
        self._document_id = QLineEdit()
        form.addRow("Document id (optional):", self._document_id)
        root.addLayout(form)

        self._overwrite = QCheckBox("Overwrite existing output")
        root.addWidget(self._overwrite)

        self._run_button = QPushButton("Build PCPT")
        self._run_button.clicked.connect(self._on_build)
        mark_primary(self._run_button)
        root.addWidget(self._run_button)
        root.addStretch()

    def _on_build(self) -> None:
        source = self._source.path()
        output = self._output.path()
        if source is None or not source.exists():
            self.log.clear()
            self.log.append("Pick a points CSV/XLSX first.")
            return
        if output is None:
            self.log.clear()
            self.log.append("Pick an output .pcpt.json path first.")
            return

        sheet_text = self._sheet.text().strip()
        sheet: str | int | None = sheet_text or None
        if sheet is not None and sheet.isdigit():
            sheet = int(sheet)

        self.run_job(
            jobs.build_pcpt_job,
            source,
            output,
            sheet=sheet,
            crs=self._crs.text().strip() or None,
            document_id=self._document_id.text().strip() or None,
            overwrite=self._overwrite.isChecked(),
            busy_message="Building PCPT…",
            success_message=lambda r: f"Wrote {r['file']}\n  points: {r['n_points']}",
        )
