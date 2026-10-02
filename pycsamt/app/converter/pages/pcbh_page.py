# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages.pcbh_page
=========================================

Build a PCBH (borehole, ``.pcbh.json``) document from a combined CSV, a
directory of relational tables, an XLSX workbook, or a single LAS 2.0
log -- the GUI face of ``pycsamt format build-pcbh``
(:mod:`pycsamt.format.borehole`).
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QLineEdit,
    QPushButton,
)

from pycsamt.app.converter import jobs
from pycsamt.app.converter.pages._base import (
    ConverterPage,
    mark_primary,
    page_root_layout,
    two_col_form,
)
from pycsamt.app.converter.widgets import DropZone


class PcbhPage(ConverterPage):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        root = page_root_layout(self)

        self._source = DropZone("Source (CSV / folder / XLSX / LAS)")
        root.addWidget(self._source)

        self._output = DropZone("Output .pcbh.json", accept_dirs=False)
        root.addWidget(self._output)

        form = QFormLayout()
        self._kind = QComboBox()
        self._kind.addItems(["auto", "csv", "csv-dir", "xlsx", "las"])
        self._kind.currentTextChanged.connect(self._on_kind_changed)
        form.addRow("Source kind:", self._kind)

        self._document_id = QLineEdit()
        form.addRow("Document id (optional):", self._document_id)

        self._created_by = QLineEdit("pyCSAMT Format Studio")
        form.addRow("Created by:", self._created_by)
        root.addLayout(form)

        self._las_box = QGroupBox("LAS collar (required when source kind is LAS)")
        self._collar_id = QLineEdit()
        self._x = QDoubleSpinBox()
        self._x.setRange(-1e9, 1e9)
        self._x.setDecimals(3)
        self._y = QDoubleSpinBox()
        self._y.setRange(-1e9, 1e9)
        self._y.setDecimals(3)
        self._z = QDoubleSpinBox()
        self._z.setRange(-1e6, 1e6)
        self._z.setDecimals(3)
        self._crs = QLineEdit()
        self._crs.setPlaceholderText("e.g. EPSG:32650")
        las_grid = two_col_form(
            [
                ("Collar id:", self._collar_id),
                ("X:", self._x),
                ("Y:", self._y),
                ("Z (elevation):", self._z),
                ("Horizontal CRS:", self._crs),
            ]
        )
        self._las_box.setLayout(las_grid)
        root.addWidget(self._las_box)
        self._las_box.setVisible(False)

        self._overwrite = QCheckBox("Overwrite existing output")
        root.addWidget(self._overwrite)

        self._run_button = QPushButton("Build PCBH")
        self._run_button.clicked.connect(self._on_build)
        mark_primary(self._run_button)
        root.addWidget(self._run_button)
        root.addStretch()

    def _on_kind_changed(self, kind: str) -> None:
        self._las_box.setVisible(kind == "las")

    def _on_build(self) -> None:
        source = self._source.path()
        output = self._output.path()
        if source is None or not source.exists():
            self.log.clear()
            self.log.append("Pick a source file or folder first.")
            return
        if output is None:
            self.log.clear()
            self.log.append("Pick an output .pcbh.json path first.")
            return

        kind = self._kind.currentText()
        if kind == "las" or (kind == "auto" and source.suffix.lower() == ".las"):
            if not self._collar_id.text().strip() or not self._crs.text().strip():
                self.log.clear()
                self.log.append("LAS import needs a collar id and CRS at minimum.")
                return

        self.run_job(
            jobs.build_pcbh_job,
            source,
            output,
            source_kind=kind,
            collar_id=self._collar_id.text().strip() or None,
            x=self._x.value(),
            y=self._y.value(),
            z=self._z.value(),
            crs_horizontal=self._crs.text().strip() or None,
            document_id=self._document_id.text().strip() or None,
            created_by=self._created_by.text() or "pyCSAMT Format Studio",
            overwrite=self._overwrite.isChecked(),
            busy_message="Building PCBH…",
            success_message=lambda r: (
                f"Wrote {r['file']}\n"
                f"  boreholes: {r['n_boreholes']}"
                + (f"   issues: {len(r['issues'])}" if "issues" in r else "")
            ),
        )
