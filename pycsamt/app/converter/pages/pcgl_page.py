# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages.pcgl_page
=========================================

Build a PCGL (geology legend, ``.pcgl.json``) document from a CSV of
named resistivity units -- the GUI face of ``pycsamt format
build-pcgl`` (:mod:`pycsamt.format.geology`).
"""

from __future__ import annotations

from PySide6.QtWidgets import QCheckBox, QFormLayout, QLineEdit, QPushButton

from pycsamt.app.converter import jobs
from pycsamt.app.converter.pages._base import ConverterPage, mark_primary, page_root_layout
from pycsamt.app.converter.widgets import DropZone


class PcglPage(ConverterPage):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        root = page_root_layout(self)

        root.addWidget(
            _hint(
                "CSV columns: name, rho_min, rho_max required; color, "
                "description, code, pattern_id optional."
            )
        )

        self._source = DropZone(
            "Units CSV", accept_dirs=False, file_filter="CSV files (*.csv)"
        )
        root.addWidget(self._source)

        self._output = DropZone("Output .pcgl.json", accept_dirs=False)
        root.addWidget(self._output)

        form = QFormLayout()
        self._title = QLineEdit()
        form.addRow("Title:", self._title)
        self._document_id = QLineEdit()
        form.addRow("Document id (optional):", self._document_id)
        self._created_by = QLineEdit("pyCSAMT Format Studio")
        form.addRow("Created by:", self._created_by)
        root.addLayout(form)

        self._overwrite = QCheckBox("Overwrite existing output")
        root.addWidget(self._overwrite)

        self._run_button = QPushButton("Build PCGL")
        self._run_button.clicked.connect(self._on_build)
        mark_primary(self._run_button)
        root.addWidget(self._run_button)
        root.addStretch()

    def _on_build(self) -> None:
        source = self._source.path()
        output = self._output.path()
        if source is None or not source.exists():
            self.log.clear()
            self.log.append("Pick a units CSV first.")
            return
        if output is None:
            self.log.clear()
            self.log.append("Pick an output .pcgl.json path first.")
            return

        self.run_job(
            jobs.build_pcgl_job,
            source,
            output,
            title=self._title.text(),
            document_id=self._document_id.text().strip() or None,
            created_by=self._created_by.text() or "pyCSAMT Format Studio",
            overwrite=self._overwrite.isChecked(),
            busy_message="Building PCGL…",
            success_message=lambda r: f"Wrote {r['file']}\n  entries: {r['n_entries']}",
        )


def _hint(text: str):
    from PySide6.QtWidgets import QLabel

    label = QLabel(text)
    label.setObjectName("PageHint")
    label.setWordWrap(True)
    return label
