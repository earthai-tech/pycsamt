# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages.inversion_page
==============================================

ModEM 3-D / Occam2D / MARE2DEM inversion result, or any AI/DL array
bundle (.npz/.npy) -> PCSF/PCSM. Auto-detects the source the same way
``pycsamt format convert`` does (:func:`pycsamt.format.detect.detect_source`).
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
)

from pycsamt.app.converter import jobs
from pycsamt.app.converter.pages._base import (
    ConverterPage,
    mark_primary,
    page_root_layout,
    two_col_form,
)
from pycsamt.app.converter.settings import load_settings
from pycsamt.app.converter.widgets import DropZone


class InversionPage(ConverterPage):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        root = page_root_layout(self)

        self._source = DropZone("Inversion result / AI array bundle (source)")
        self._source.pathChanged.connect(self._on_source_changed)
        root.addWidget(self._source)

        self._detected = QLabel("Drop a file or folder to auto-detect it.")
        self._detected.setObjectName("PageHint")
        self._detected.setWordWrap(True)
        root.addWidget(self._detected)

        self._output_dir = DropZone(
            "Output folder", accept_files=False, accept_dirs=True
        )
        root.addWidget(self._output_dir)

        basic = QFormLayout()
        self._to_format = QComboBox()
        self._to_format.addItems(["pcsf", "pcsm", "pcsm.gz"])
        basic.addRow("Output format:", self._to_format)
        root.addLayout(basic)

        adv_box = QGroupBox("Advanced options")
        adv = QVBoxLayout(adv_box)
        adv.setSpacing(10)

        self._solver = QComboBox()
        self._solver.addItems(["auto", "occam2d", "modem", "mare2dem"])

        self._iteration = QSpinBox()
        self._iteration.setRange(0, 9999)
        self._iteration.setSpecialValueText("final (default)")

        self._epsg = QSpinBox()
        self._epsg.setRange(0, 99999)
        self._epsg.setSpecialValueText("(none)")

        self._utm_zone = QLineEdit()
        self._utm_zone.setPlaceholderText("e.g. 48N")

        self._encoding = QComboBox()
        self._encoding.addItems(["auto", "linear", "log10", "ln"])

        self._station_z = QComboBox()
        self._station_z.addItems(["auto", "elevation", "depth_down"])

        self._air_threshold = QDoubleSpinBox()
        self._air_threshold.setRange(0.0, 1e12)
        self._air_threshold.setDecimals(0)
        self._air_threshold.setValue(1e8)

        self._created_by = QLineEdit("pyCSAMT Format Studio")

        self._description = QLineEdit()

        # Short scalar fields two-up, so 9 fields take 5 rows instead of 9.
        adv.addLayout(
            two_col_form(
                [
                    ("Force solver:", self._solver),
                    ("Occam2D iteration:", self._iteration),
                    ("EPSG:", self._epsg),
                    ("UTM zone:", self._utm_zone),
                    ("AI array encoding:", self._encoding),
                    ("ModEM station-Z:", self._station_z),
                    ("ModEM air threshold (Ω·m):", self._air_threshold),
                    ("Created by:", self._created_by),
                    ("Description:", self._description),
                ]
            )
        )

        self._topo = DropZone("Topography source (optional)")
        adv.addWidget(self._topo)

        self._poly = DropZone(
            "MARE2DEM .poly (optional)", accept_dirs=False, accept_files=True
        )
        adv.addWidget(self._poly)

        self._log10_view = QCheckBox("Write PCSM block in log10(rho)")
        adv.addWidget(self._log10_view)

        self._overwrite = QCheckBox("Overwrite existing output")
        adv.addWidget(self._overwrite)

        root.addWidget(adv_box)

        self._run_button = QPushButton("Convert")
        self._run_button.clicked.connect(self._on_convert)
        mark_primary(self._run_button)
        root.addWidget(self._run_button)
        root.addStretch()

    # -- behaviour --------------------------------------------------------

    def _on_source_changed(self, _text: str) -> None:
        source = self._source.path()
        if source is None or not source.exists():
            self._detected.setText("Drop a file or folder to auto-detect it.")
            return
        try:
            info = jobs.detect_source_job(source)
        except Exception as exc:  # noqa: BLE001
            self._detected.setText(f"Could not classify this source: {exc}")
            return
        self._detected.setText(
            f"Detected: {info['category']}"
            + (f" ({info['backend']})" if info.get("backend") else "")
            + f"  —  {info['detail']}"
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

        settings = load_settings()
        solver = None if self._solver.currentText() == "auto" else self._solver.currentText()
        iteration = None if self._iteration.value() == 0 else self._iteration.value()
        topo = self._topo.path()
        epsg = self._epsg.value() or settings.default_epsg or None
        utm_zone = self._utm_zone.text().strip() or settings.default_utm_zone or None
        encoding = None if self._encoding.currentText() == "auto" else self._encoding.currentText()
        poly = self._poly.path()

        self.run_job(
            jobs.convert_to_pcsf_job,
            source,
            None,
            self._to_format.currentText(),
            output_dir,
            overwrite=self._overwrite.isChecked() or settings.overwrite_without_asking,
            solver=solver,
            iteration=iteration,
            topo=topo,
            epsg=epsg,
            utm_zone=utm_zone,
            encoding=encoding,
            station_z_convention=self._station_z.currentText(),
            air_threshold_ohm_m=self._air_threshold.value() or None,
            poly=poly,
            created_by=self._created_by.text() or "pyCSAMT Format Studio",
            description=self._description.text(),
            log10_view=self._log10_view.isChecked() or settings.log10_view,
            busy_message="Converting…",
            success_message=lambda r: (
                f"Wrote {r['file']}\n"
                f"  resistivity shape: {r.get('resistivity_shape')}\n"
                f"  stations: {r.get('n_stations', 0)}"
            ),
        )
