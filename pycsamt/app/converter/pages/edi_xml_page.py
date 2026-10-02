# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages.edi_xml_page
============================================

Lossless EDI <-> EMTF-XML conversion, one file or a whole directory at
a time -- the GUI face of ``pycsamt format edi-to-xml``/``xml-to-edi``
(:mod:`pycsamt.emtf`).
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QPushButton,
)

from pycsamt.app.converter import jobs
from pycsamt.app.converter.pages._base import ConverterPage, mark_primary, page_root_layout
from pycsamt.app.converter.widgets import DropZone


class EdiXmlPage(ConverterPage):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        root = page_root_layout(self)

        form = QFormLayout()
        self._direction = QComboBox()
        self._direction.addItems(["EDI → EMTF-XML", "EMTF-XML → EDI"])
        self._direction.currentIndexChanged.connect(self._on_direction_changed)
        form.addRow("Direction:", self._direction)
        root.addLayout(form)

        self._source = DropZone("Source (.edi file or folder)")
        root.addWidget(self._source)

        self._output_dir = DropZone(
            "Output folder", accept_files=False, accept_dirs=True
        )
        root.addWidget(self._output_dir)

        opts = QHBoxLayout()
        self._prefer_spectra = QCheckBox("Prefer EDI SPECTRA blocks")
        self._prefer_spectra.setChecked(True)
        opts.addWidget(self._prefer_spectra)

        self._on_loss = QComboBox()
        self._on_loss.addItems(["warn", "raise", "ignore"])
        opts.addWidget(self._on_loss)

        self._strict = QCheckBox("Strict XML parsing")
        self._strict.setChecked(True)
        opts.addWidget(self._strict)

        self._overwrite = QCheckBox("Overwrite existing output")
        opts.addWidget(self._overwrite)
        root.addLayout(opts)

        self._run_button = QPushButton("Convert")
        self._run_button.clicked.connect(self._on_convert)
        mark_primary(self._run_button)
        root.addWidget(self._run_button)
        root.addStretch()
        self._on_direction_changed(0)

    def _on_direction_changed(self, index: int) -> None:
        is_edi_to_xml = index == 0
        self._prefer_spectra.setVisible(is_edi_to_xml)
        self._on_loss.setVisible(not is_edi_to_xml)
        self._strict.setVisible(not is_edi_to_xml)
        self._source.setToolTip(
            "A .edi file or a folder of them" if is_edi_to_xml
            else "A .xml file or a folder of them"
        )

    def _on_convert(self) -> None:
        source = self._source.path()
        output_dir = self._output_dir.path()
        if source is None or not source.exists():
            self.log.clear()
            self.log.append("Pick a source file or folder first.")
            return
        if output_dir is None:
            self.log.clear()
            self.log.append("Pick an output folder first.")
            return

        if self._direction.currentIndex() == 0:
            self.run_job(
                jobs.edi_to_xml_job,
                source,
                output_dir,
                prefer_spectra=self._prefer_spectra.isChecked(),
                overwrite=self._overwrite.isChecked(),
                supports_progress=True,
                busy_message="Converting EDI → EMTF-XML…",
                success_message=lambda r: f"Wrote {len(r)} EMTF-XML file(s) to {output_dir}",
            )
        else:
            self.run_job(
                jobs.xml_to_edi_job,
                source,
                output_dir,
                on_loss=self._on_loss.currentText(),
                strict=self._strict.isChecked(),
                overwrite=self._overwrite.isChecked(),
                supports_progress=True,
                busy_message="Converting EMTF-XML → EDI…",
                success_message=lambda r: f"Wrote {len(r)} EDI file(s) to {output_dir}",
            )
