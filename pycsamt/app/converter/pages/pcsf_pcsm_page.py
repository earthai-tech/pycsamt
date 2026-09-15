# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages.pcsf_pcsm_page
==============================================

Direct PCSF <-> PCSM transcoding, plus Validate and Info actions on an
existing .pcsf/.pcsm file -- the GUI face of ``pycsamt format
validate``/``info`` and the PCSF<->PCSM half of ``pycsamt format convert``.
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


class PcsfPcsmPage(ConverterPage):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        root = page_root_layout(self)

        self._source = DropZone(
            "PCSF or PCSM file", accept_dirs=False, file_filter="PCSF/PCSM (*.pcsf *.pcsm *.pcsm.gz)"
        )
        root.addWidget(self._source)

        form = QFormLayout()
        self._target_format = QComboBox()
        self._target_format.addItems(["pcsm", "pcsf"])
        form.addRow("Transcode to:", self._target_format)
        root.addLayout(form)

        self._output = DropZone("Output file", accept_dirs=False)
        root.addWidget(self._output)

        opts = QHBoxLayout()
        self._log10_view = QCheckBox("PCSM in log10(rho)")
        opts.addWidget(self._log10_view)
        self._overwrite = QCheckBox("Overwrite existing output")
        opts.addWidget(self._overwrite)
        root.addLayout(opts)

        buttons = QHBoxLayout()
        self._transcode_btn = QPushButton("Transcode")
        self._transcode_btn.clicked.connect(self._on_transcode)
        mark_primary(self._transcode_btn)
        buttons.addWidget(self._transcode_btn)

        self._validate_btn = QPushButton("Validate")
        self._validate_btn.clicked.connect(self._on_validate)
        buttons.addWidget(self._validate_btn)

        self._info_btn = QPushButton("Info")
        self._info_btn.clicked.connect(self._on_info)
        buttons.addWidget(self._info_btn)
        root.addLayout(buttons)
        root.addStretch()

    def _on_transcode(self) -> None:
        source = self._source.path()
        output = self._output.path()
        if source is None or not source.exists():
            self.log.clear()
            self.log.append("Pick a .pcsf/.pcsm source file first.")
            return
        if output is None:
            self.log.clear()
            self.log.append("Pick an output file path first.")
            return
        self.run_job(
            jobs.transcode_pcsf_pcsm_job,
            source,
            output,
            self._target_format.currentText(),
            log10_view=self._log10_view.isChecked(),
            overwrite=self._overwrite.isChecked(),
            run_button=self._transcode_btn,
            busy_message="Transcoding…",
            success_message=lambda r: f"Wrote {r['file']}",
        )

    def _on_validate(self) -> None:
        source = self._source.path()
        if source is None or not source.exists():
            self.log.clear()
            self.log.append("Pick a .pcsf/.pcsm file first.")
            return
        self.run_job(
            jobs.validate_pcsf_job,
            source,
            run_button=self._validate_btn,
            busy_message="Validating…",
            success_message=self._format_validation,
        )

    def _on_info(self) -> None:
        source = self._source.path()
        if source is None or not source.exists():
            self.log.clear()
            self.log.append("Pick a .pcsf/.pcsm file first.")
            return
        self.run_job(
            jobs.info_pcsf_job,
            source,
            run_button=self._info_btn,
            busy_message="Reading…",
            success_message=self._format_info,
        )

    @staticmethod
    def _format_validation(result: dict) -> str:
        lines = [f"{'VALID' if result['valid'] else 'INVALID'} — {result['file']}"]
        for check in result["checks"]:
            mark = "✓" if check["ok"] else "✗"
            lines.append(f"  {mark} {check['step']}: {check['detail'] or ('pass' if check['ok'] else 'FAIL')}")
        return "\n".join(lines)

    @staticmethod
    def _format_info(result: dict) -> str:
        lines = [
            f"{result['file']}  ({result['size_bytes']:,} bytes)",
            f"  container: {result['container']}   geometry: {result['geometry_kind']}",
            f"  backend: {result['source_backend']}   crs: {result['crs'] or '-'}",
        ]
        rho = result.get("resistivity")
        if rho:
            rng = f"{rho['min']:.4g} … {rho['max']:.4g} (median {rho['median']:.4g})" if "min" in rho else "all NaN"
            lines.append(f"  resistivity {rho['shape']}  {rng}")
        if "stations" in result:
            lines.append(f"  stations: {result['stations']['n']}")
        return "\n".join(lines)
