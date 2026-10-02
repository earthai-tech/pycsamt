# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages.pcgs_page
=========================================

Build a PCGS (structural evidence, ``.pcgs.json``) document from up to
three CSVs (planar / linear / faults) -- the GUI face of ``pycsamt
format build-pcgs`` (:mod:`pycsamt.format.structure`).
"""

from __future__ import annotations

from PySide6.QtWidgets import QCheckBox, QFormLayout, QLabel, QLineEdit, QPushButton

from pycsamt.app.converter import jobs
from pycsamt.app.converter.pages._base import ConverterPage, mark_primary, page_root_layout
from pycsamt.app.converter.widgets import DropZone


class PcgsPage(ConverterPage):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        root = page_root_layout(self)

        hint = QLabel(
            "At least one CSV is required. Planar: x, kind, strike_deg, "
            "dip_deg, dip_direction_deg. Linear: x, kind, trend_deg, "
            "plunge_deg. Faults: x, dip_deg, downthrown_side."
        )
        hint.setObjectName("PageHint")
        hint.setWordWrap(True)
        root.addWidget(hint)

        self._planar = DropZone(
            "Planar measurements CSV (optional)", accept_dirs=False,
            file_filter="CSV files (*.csv)",
        )
        root.addWidget(self._planar)
        self._linear = DropZone(
            "Linear measurements CSV (optional)", accept_dirs=False,
            file_filter="CSV files (*.csv)",
        )
        root.addWidget(self._linear)
        self._faults = DropZone(
            "Fault traces CSV (optional)", accept_dirs=False,
            file_filter="CSV files (*.csv)",
        )
        root.addWidget(self._faults)

        self._output = DropZone("Output .pcgs.json", accept_dirs=False)
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

        self._run_button = QPushButton("Build PCGS")
        self._run_button.clicked.connect(self._on_build)
        mark_primary(self._run_button)
        root.addWidget(self._run_button)
        root.addStretch()

    def _on_build(self) -> None:
        planar = self._planar.path()
        linear = self._linear.path()
        faults = self._faults.path()
        output = self._output.path()
        if not any(p and p.exists() for p in (planar, linear, faults)):
            self.log.clear()
            self.log.append("Pick at least one of the planar/linear/faults CSVs.")
            return
        if output is None:
            self.log.clear()
            self.log.append("Pick an output .pcgs.json path first.")
            return

        self.run_job(
            jobs.build_pcgs_job,
            output,
            planar_path=planar if planar and planar.exists() else None,
            linear_path=linear if linear and linear.exists() else None,
            faults_path=faults if faults and faults.exists() else None,
            title=self._title.text(),
            document_id=self._document_id.text().strip() or None,
            created_by=self._created_by.text() or "pyCSAMT Format Studio",
            overwrite=self._overwrite.isChecked(),
            busy_message="Building PCGS…",
            success_message=lambda r: (
                f"Wrote {r['file']}\n"
                f"  planar: {r['n_planar']}   linear: {r['n_linear']}   "
                f"faults: {r['n_faults']}"
            ),
        )
