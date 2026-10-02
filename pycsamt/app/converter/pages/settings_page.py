# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages.settings_page
=============================================

Persisted defaults applied to new conversion jobs across every page --
see :mod:`pycsamt.app.converter.settings`. Other pages call
:func:`~pycsamt.app.converter.settings.load_settings` fresh each time
Convert is clicked, so a change here takes effect immediately without
needing to reopen the app.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.converter.pages._base import mark_primary, page_root_layout, two_col_form
from pycsamt.app.converter.settings import ConverterSettings, load_settings, save_settings


class SettingsPage(QWidget):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build_ui()
        self._load()

    def _build_ui(self) -> None:
        root = page_root_layout(self)

        self._epsg = QSpinBox()
        self._epsg.setRange(0, 99999)
        self._epsg.setSpecialValueText("(none)")

        self._utm_zone = QLineEdit()
        self._utm_zone.setPlaceholderText("e.g. 48N")

        self._station_z = QComboBox()
        self._station_z.addItems(["auto", "elevation", "depth_down"])

        self._air_threshold = QDoubleSpinBox()
        self._air_threshold.setRange(0.0, 1e12)
        self._air_threshold.setDecimals(0)

        self._on_loss = QComboBox()
        self._on_loss.addItems(["warn", "raise", "ignore"])

        root.addLayout(
            two_col_form(
                [
                    ("Default EPSG:", self._epsg),
                    ("Default UTM zone:", self._utm_zone),
                    ("ModEM station-Z convention:", self._station_z),
                    ("ModEM air threshold (Ω·m):", self._air_threshold),
                    ("EDI/XML data-loss policy:", self._on_loss),
                ]
            )
        )

        checks = QVBoxLayout()
        checks.setSpacing(6)

        self._prefer_spectra = QCheckBox("Prefer EDI SPECTRA blocks")
        checks.addWidget(self._prefer_spectra)

        self._xml_strict = QCheckBox("Strict EMTF-XML parsing")
        checks.addWidget(self._xml_strict)

        self._log10_view = QCheckBox("Write PCSM blocks in log10(rho) by default")
        checks.addWidget(self._log10_view)

        self._open_after = QCheckBox("Open output folder after each conversion")
        checks.addWidget(self._open_after)

        self._overwrite = QCheckBox("Overwrite existing output without asking")
        checks.addWidget(self._overwrite)

        root.addLayout(checks)

        self._status = QLabel("")
        self._status.setObjectName("PageHint")
        root.addWidget(self._status)

        self._save_button = QPushButton("Save settings")
        self._save_button.clicked.connect(self._on_save)
        mark_primary(self._save_button)
        root.addWidget(self._save_button)
        root.addStretch()

    def _load(self) -> None:
        s = load_settings()
        self._epsg.setValue(s.default_epsg or 0)
        self._utm_zone.setText(s.default_utm_zone)
        self._station_z.setCurrentText(s.station_z_convention)
        self._air_threshold.setValue(s.air_threshold_ohm_m)
        self._on_loss.setCurrentText(s.on_loss)
        self._prefer_spectra.setChecked(s.prefer_spectra)
        self._xml_strict.setChecked(s.xml_strict)
        self._log10_view.setChecked(s.log10_view)
        self._open_after.setChecked(s.open_output_folder_after_convert)
        self._overwrite.setChecked(s.overwrite_without_asking)

    def _on_save(self) -> None:
        settings = ConverterSettings(
            default_epsg=self._epsg.value() or None,
            default_utm_zone=self._utm_zone.text().strip(),
            station_z_convention=self._station_z.currentText(),
            air_threshold_ohm_m=self._air_threshold.value(),
            on_loss=self._on_loss.currentText(),
            prefer_spectra=self._prefer_spectra.isChecked(),
            xml_strict=self._xml_strict.isChecked(),
            log10_view=self._log10_view.isChecked(),
            open_output_folder_after_convert=self._open_after.isChecked(),
            overwrite_without_asking=self._overwrite.isChecked(),
        )
        save_settings(settings)
        self._status.setText("Settings saved.")
