# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
CoordTransformDialog — UTM ↔ Geographic (Lat / Lon) coordinate converter.

Back-end priority
-----------------
1. pyproj  (installed — v3.6.1)
2. pycsamt.gis.utils.ll_to_utm / utm_to_ll  (pure-Python fallback)

Features
--------
• Bidirectional: Lat/Lon → UTM   or   UTM → Lat/Lon
• Zone auto-detection from longitude (or override)
• Hemisphere N / S toggle
• Datum: WGS84 (default) or NAD83
• Paste tabular input  (one point per line: lat lon  or  easting northing)
• Copy button on output
• Survey-mode: load station list from loaded survey and show all coordinates
"""

from __future__ import annotations

import re

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QRadioButton,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

_DATUMS = ["WGS84", "NAD83", "GRS80"]


# ── Coordinate math ───────────────────────────────────────────────────────────


def _ll_to_utm(lat: float, lon: float, zone: int | None, hem: str, datum: str):
    """Return (easting, northing, zone) via pyproj or fallback."""
    try:
        from pyproj import Proj

        if zone is None:
            zone = int((lon + 180) / 6) + 1
        proj = Proj(
            proj="utm",
            zone=zone,
            datum=datum,
            south=(hem == "S"),
            ellps=datum,
        )
        e, n = proj(lon, lat)
        return e, n, zone
    except Exception:
        from pycsamt.gis.utils import ll_to_utm

        res = ll_to_utm(lat, lon)
        return (
            res["easting"],
            res["northing"],
            res.get("zone_number", zone or 0),
        )


def _utm_to_ll(
    easting: float, northing: float, zone: int, hem: str, datum: str
):
    """Return (lat, lon) via pyproj or fallback."""
    try:
        from pyproj import Proj

        proj = Proj(
            proj="utm",
            zone=zone,
            datum=datum,
            south=(hem == "S"),
            ellps=datum,
        )
        lon, lat = proj(easting, northing, inverse=True)
        return lat, lon
    except Exception:
        from pycsamt.gis.utils import utm_to_ll

        res = utm_to_ll(easting, northing, zone, hem)
        return res["lat"], res["lon"]


def _parse_pairs(text: str) -> list[tuple[float, float]]:
    """Parse whitespace / comma / tab separated pairs from multi-line text."""
    pairs = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"[,\s\t]+", line)
        if len(parts) >= 2:
            try:
                pairs.append((float(parts[0]), float(parts[1])))
            except ValueError:
                pass
    return pairs


# ── Dialog ────────────────────────────────────────────────────────────────────


class CoordTransformDialog(QDialog):
    """
    UTM ↔ Geographic coordinate converter.

    Parameters
    ----------
    sites : any, optional
        If supplied, the Survey tab is pre-populated with station coordinates.
    parent : QWidget, optional
    """

    def __init__(self, sites=None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Coordinate Transformer")
        self.setMinimumSize(700, 560)
        self._sites = sites
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(10)

        # ── Direction + options row ───────────────────────────────────────────
        opts_grp = QGroupBox("Options")
        opts_form = QFormLayout(opts_grp)
        opts_form.setSpacing(8)

        dir_row = QHBoxLayout()
        self._ll_to_utm_rb = QRadioButton("Lat / Lon  →  UTM")
        self._utm_to_ll_rb = QRadioButton("UTM  →  Lat / Lon")
        self._ll_to_utm_rb.setChecked(True)
        dir_row.addWidget(self._ll_to_utm_rb)
        dir_row.addWidget(self._utm_to_ll_rb)
        dir_row.addStretch()
        opts_form.addRow("Direction:", dir_row)

        zone_row = QHBoxLayout()
        self._zone_spin = QSpinBox()
        self._zone_spin.setRange(1, 60)
        self._zone_spin.setValue(30)
        self._auto_zone_lbl = QLabel("(auto-detected from longitude)")
        zone_row.addWidget(self._zone_spin)
        zone_row.addWidget(self._auto_zone_lbl)
        opts_form.addRow("UTM zone:", zone_row)

        hem_row = QHBoxLayout()
        self._hem_n = QRadioButton("N  (northern)")
        self._hem_s = QRadioButton("S  (southern)")
        self._hem_n.setChecked(True)
        hem_row.addWidget(self._hem_n)
        hem_row.addWidget(self._hem_s)
        opts_form.addRow("Hemisphere:", hem_row)

        self._datum_combo = QComboBox()
        self._datum_combo.addItems(_DATUMS)
        opts_form.addRow("Datum:", self._datum_combo)

        root.addWidget(opts_grp)

        # ── Manual / batch input ──────────────────────────────────────────────
        split = QSplitter(Qt.Orientation.Horizontal)

        in_widget = QWidget()
        in_layout = QVBoxLayout(in_widget)
        in_layout.addWidget(QLabel("<b>Input</b>  (one point per line)"))
        self._in_edit = QTextEdit()
        self._in_edit.setPlaceholderText(
            "Lat/Lon mode:\n  51.5074   -0.1278\n  48.8566    2.3522\n\n"
            "UTM mode:\n  699369   5710671\n  453000   5411000"
        )
        mono = QFont("Monospace")
        mono.setStyleHint(QFont.StyleHint.Monospace)
        self._in_edit.setFont(mono)
        in_layout.addWidget(self._in_edit)
        split.addWidget(in_widget)

        out_widget = QWidget()
        out_layout = QVBoxLayout(out_widget)
        hdr = QHBoxLayout()
        hdr.addWidget(QLabel("<b>Output</b>"))
        copy_btn = QPushButton("Copy")
        copy_btn.setFixedWidth(60)
        copy_btn.clicked.connect(self._copy_output)
        hdr.addWidget(copy_btn)
        out_layout.addLayout(hdr)
        self._out_edit = QTextEdit()
        self._out_edit.setReadOnly(True)
        self._out_edit.setFont(mono)
        out_layout.addWidget(self._out_edit)
        split.addWidget(out_widget)

        root.addWidget(split, stretch=1)

        # ── Survey table ──────────────────────────────────────────────────────
        if self._sites is not None:
            self._build_survey_table(root)

        # ── Buttons ───────────────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        run_btn = QPushButton("Transform")
        run_btn.clicked.connect(self._on_transform)
        btn_row.addWidget(run_btn)
        btn_row.addStretch()
        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        btn_row.addWidget(box)
        root.addLayout(btn_row)

    def _build_survey_table(self, root: QVBoxLayout) -> None:
        """Build the station-list section when a survey is loaded."""
        grp = QGroupBox("Survey Stations (from loaded data)")
        vlay = QVBoxLayout(grp)
        self._surv_table = QTableWidget(0, 5)
        self._surv_table.setHorizontalHeaderLabels(
            ["Station", "Lat (°)", "Lon (°)", "UTM E (m)", "UTM N (m)"]
        )
        self._surv_table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        self._surv_table.setEditTriggers(
            QTableWidget.EditTrigger.NoEditTriggers
        )
        self._surv_table.setMaximumHeight(140)
        vlay.addWidget(self._surv_table)
        root.addWidget(grp)
        self._populate_survey_table()

    def _populate_survey_table(self) -> None:
        try:
            from pycsamt.emtools._core import (
                _iter_items,
                _unwrap,
            )

            items = list(_iter_items(self._sites))
        except Exception:
            return
        self._zone_spin.value()
        hem = "N" if self._hem_n.isChecked() else "S"
        datum = self._datum_combo.currentText()
        self._surv_table.setRowCount(0)
        for ed in items:
            try:
                ed = _unwrap(ed)
            except Exception:
                pass
            name = getattr(ed, "station", None) or getattr(ed, "id", "?")
            lat = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
            lon = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
            r = self._surv_table.rowCount()
            self._surv_table.insertRow(r)
            self._surv_table.setItem(r, 0, QTableWidgetItem(str(name)))
            if lat is not None and lon is not None:
                try:
                    la, lo = float(lat), float(lon)
                    e, n, z = _ll_to_utm(la, lo, None, hem, datum)
                    self._surv_table.setItem(
                        r, 1, QTableWidgetItem(f"{la:.6f}")
                    )
                    self._surv_table.setItem(
                        r, 2, QTableWidgetItem(f"{lo:.6f}")
                    )
                    self._surv_table.setItem(
                        r, 3, QTableWidgetItem(f"{e:.1f}")
                    )
                    self._surv_table.setItem(
                        r, 4, QTableWidgetItem(f"{n:.1f}")
                    )
                except Exception:
                    pass
            else:
                for c in range(1, 5):
                    self._surv_table.setItem(r, c, QTableWidgetItem("—"))

    # ── Transform ─────────────────────────────────────────────────────────────

    def _on_transform(self) -> None:
        text = self._in_edit.toPlainText()
        pairs = _parse_pairs(text)
        if not pairs:
            self._out_edit.setPlainText("No valid input pairs found.")
            return

        zone = self._zone_spin.value()
        hem = "N" if self._hem_n.isChecked() else "S"
        datum = self._datum_combo.currentText()
        ll2utm = self._ll_to_utm_rb.isChecked()

        lines = []
        for a, b in pairs:
            try:
                if ll2utm:
                    e, n, z = _ll_to_utm(a, b, zone, hem, datum)
                    lines.append(f"{e:>14.3f}   {n:>14.3f}   zone {z}{hem}")
                else:
                    la, lo = _utm_to_ll(a, b, zone, hem, datum)
                    lines.append(f"{la:>12.6f}   {lo:>13.6f}")
            except Exception as exc:
                lines.append(f"ERROR: {exc}")

        self._out_edit.setPlainText("\n".join(lines))

    def _copy_output(self) -> None:
        from PySide6.QtWidgets import QApplication

        QApplication.clipboard().setText(self._out_edit.toPlainText())
