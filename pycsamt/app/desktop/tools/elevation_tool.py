# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ElevationEnrichDialog — fetch elevation for all loaded stations via open APIs.

Reads lat/lon from every station in the survey, queries an elevation API
(open-meteo or open-topo-data), and shows the result in a table.  Requires an
active internet connection and the ``requests`` library.

Calls ``pycsamt.gis.utils.get_elevation_from_api``.

Usage
-----
    from pycsamt.app.desktop.tools.elevation_tool import ElevationEnrichDialog
    dlg = ElevationEnrichDialog(sites, parent=self)
    dlg.exec()
"""
from __future__ import annotations

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QProgressBar,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

_APIS = ["open_meteo", "open_topo_data"]
_C_OK   = QColor("#c8e6c9")
_C_WARN = QColor("#fff9c4")
_C_ERR  = QColor("#ffcdd2")


# ── Worker ────────────────────────────────────────────────────────────────────

class _ElevWorker(QThread):
    progress = Signal(int, int, str, float)  # cur, total, station, elevation
    done     = Signal(list)                  # list of (name, lat, lon, elev)
    error    = Signal(str)

    def __init__(self, stations: list, api: str):
        super().__init__()
        self._stations = stations   # list of (name, lat, lon)
        self._api      = api

    def run(self):
        try:
            from pycsamt.gis.utils import (
                get_elevation_from_api,
            )
        except ImportError as exc:
            self.error.emit(f"Cannot import gis.utils: {exc}")
            return

        results = []
        total = len(self._stations)
        for idx, (name, lat, lon) in enumerate(self._stations):
            try:
                if lat is None or lon is None:
                    raise ValueError("no coordinates")
                elev = float(get_elevation_from_api(float(lat), float(lon),
                                                    api_name=self._api))
                self.progress.emit(idx + 1, total, str(name), elev)
                results.append((str(name), lat, lon, elev))
            except Exception:
                self.progress.emit(idx + 1, total, str(name), float("nan"))
                results.append((str(name), lat, lon, float("nan")))

        self.done.emit(results)


# ── Dialog ────────────────────────────────────────────────────────────────────

class ElevationEnrichDialog(QDialog):
    """
    Fetch and display elevation for every station in the loaded survey.

    Parameters
    ----------
    sites : any
        Loaded survey.
    parent : QWidget, optional
    """

    def __init__(self, sites=None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Elevation Enrichment")
        self.setMinimumSize(700, 520)
        self._sites   = sites
        self._worker  = None
        self._station_list: list[tuple] = []   # (name, lat, lon)
        self._results: list[tuple]      = []   # (name, lat, lon, elev)
        self._build_ui()
        self._populate_stations()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(8)

        # ── Controls ──────────────────────────────────────────────────────
        ctrl = QHBoxLayout()
        grp = QGroupBox("API")
        form = QFormLayout(grp)
        self._api_combo = QComboBox()
        self._api_combo.addItems(_APIS)
        form.addRow("Service:", self._api_combo)
        ctrl.addWidget(grp)
        ctrl.addStretch()

        btn_col = QVBoxLayout()
        self._run_btn = QPushButton("Fetch elevations")
        self._run_btn.clicked.connect(self._on_run)
        btn_col.addWidget(self._run_btn)
        self._export_btn = QPushButton("Export CSV…")
        self._export_btn.setEnabled(False)
        self._export_btn.clicked.connect(self._on_export)
        btn_col.addWidget(self._export_btn)
        ctrl.addLayout(btn_col)
        root.addLayout(ctrl)

        # ── Progress ──────────────────────────────────────────────────────
        self._progress = QProgressBar()
        self._progress.setValue(0)
        root.addWidget(self._progress)

        self._status_lbl = QLabel("")
        root.addWidget(self._status_lbl)

        # ── Table ─────────────────────────────────────────────────────────
        self._table = QTableWidget(0, 4)
        self._table.setHorizontalHeaderLabels(
            ["Station", "Lat (°)", "Lon (°)", "Elevation (m)"]
        )
        self._table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        root.addWidget(self._table, stretch=1)

        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        root.addWidget(box)

    # ── Populate station list ──────────────────────────────────────────────────

    def _populate_stations(self) -> None:
        if self._sites is None:
            self._status_lbl.setText("No survey loaded.")
            self._run_btn.setEnabled(False)
            return
        try:
            from pycsamt.emtools._core import _iter_items
            for ed in _iter_items(self._sites):
                # unwrap TFBundle / container objects if needed
                for attr in ("edi", "site", "data"):
                    inner = getattr(ed, attr, None)
                    if inner is not None and hasattr(inner, "Z"):
                        ed = inner
                        break
                name = (
                    getattr(ed, "station", None)
                    or getattr(ed, "id", "?")
                )
                lat = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
                lon = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
                try:
                    lat = float(lat) if lat is not None else None
                    lon = float(lon) if lon is not None else None
                except Exception:
                    lat = lon = None
                self._station_list.append((str(name), lat, lon))
        except Exception as exc:
            self._status_lbl.setText(f"Cannot read stations: {exc}")
            self._run_btn.setEnabled(False)
            return

        n_with_coords = sum(
            1 for _, la, lo in self._station_list
            if la is not None and lo is not None
        )
        self._status_lbl.setText(
            f"{len(self._station_list)} stations found, "
            f"{n_with_coords} with coordinates."
        )
        if n_with_coords == 0:
            self._run_btn.setEnabled(False)

    # ── Run ───────────────────────────────────────────────────────────────────

    def _on_run(self) -> None:
        self._run_btn.setEnabled(False)
        self._export_btn.setEnabled(False)
        self._table.setRowCount(0)
        self._results = []
        self._progress.setValue(0)
        self._status_lbl.setText("Fetching elevations…  (requires internet)")

        self._worker = _ElevWorker(
            self._station_list,
            api=self._api_combo.currentText(),
        )
        self._worker.progress.connect(self._on_progress)
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_progress(self, cur: int, total: int, name: str, elev: float) -> None:
        self._progress.setMaximum(total)
        self._progress.setValue(cur)
        import math
        ok = not math.isnan(elev)
        r = self._table.rowCount()
        self._table.insertRow(r)
        bg = _C_OK if ok else _C_ERR

        for col, val in enumerate((name, "—", "—", f"{elev:.1f}" if ok else "ERROR")):
            it = QTableWidgetItem(val)
            it.setBackground(bg)
            it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self._table.setItem(r, col, it)

    def _on_done(self, results: list) -> None:
        self._run_btn.setEnabled(True)
        self._results = results
        import math
        n_ok  = sum(1 for _, _, _, e in results if not math.isnan(e))
        n_err = len(results) - n_ok
        self._status_lbl.setText(
            f"Done — {n_ok} elevations fetched, {n_err} errors."
        )
        self._export_btn.setEnabled(bool(self._results))
        # Fill lat/lon columns from station list
        for r, (_name, lat, lon, _) in enumerate(results):
            self._table.setItem(
                r, 1,
                QTableWidgetItem(f"{lat:.6f}" if lat is not None else "—")
            )
            self._table.setItem(
                r, 2,
                QTableWidgetItem(f"{lon:.6f}" if lon is not None else "—")
            )

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._status_lbl.setText(f"Error: {msg}")

    # ── Export ────────────────────────────────────────────────────────────────

    def _on_export(self) -> None:
        if not self._results:
            return
        from pathlib import Path
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Elevations",
            str(Path.home() / "station_elevations.csv"),
            "CSV (*.csv)",
        )
        if not path:
            return
        try:
            import csv
            import math
            with open(path, "w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["station", "lat", "lon", "elevation_m"])
                for name, lat, lon, elev in self._results:
                    writer.writerow([
                        name,
                        f"{lat:.6f}" if lat is not None else "",
                        f"{lon:.6f}" if lon is not None else "",
                        f"{elev:.2f}" if not math.isnan(elev) else "",
                    ])
            self._status_lbl.setText(f"Saved → {path}")
        except Exception as exc:
            self._status_lbl.setText(f"Export error: {exc}")
