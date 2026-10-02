# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
StationEditorDialog — Edit ▸ Stations (coordinates, names, lines, header).

One dialog, four tabs over the same station list:

* **Coordinates** — latitude / longitude / elevation, or easting /
  northing in any EPSG (UTM zone suggested).  Paste a block copied from a
  spreadsheet at the current cell (Ctrl+V), import a CSV, or type.
* **Names** — prefix / suffix, find & replace (optionally a regular
  expression), upper / lower case, zero-padding of the station number;
  the new names can also be typed.  Empty or duplicate names block Apply.
* **Lines** — auto-detect from the station names, type a line, or rename a
  whole line at once.
* **Header** — EDI HEAD fields (acquired by, dates, project, survey,
  prospect, location, country, datum, declination); a value can be filled
  into every station, or only the selected rows.

Changed cells are highlighted and counted.  **Apply** reduces everything
to the cells that changed (:func:`~pycsamt.app.desktop.controllers.
station_edits.table_changes`) and hands them to the main window, which
applies them as one validated, undoable edit.
"""

from __future__ import annotations

import math

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers import station_edits as se

TABS = ("coordinates", "names", "lines", "header")
_CHANGED = QColor(255, 214, 102, 110)  # soft amber, readable in both themes


class _Table(QTableWidget):
    """Station table whose edits are compared with the original."""

    def __init__(self, columns: list[tuple[str, str]], *, editable: set,
                 parent=None) -> None:
        super().__init__(parent)
        self.keys = [k for k, _l in columns]
        self.editable = set(editable)
        self.setColumnCount(len(columns))
        self.setHorizontalHeaderLabels([lbl for _k, lbl in columns])
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.setAlternatingRowColors(True)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self._original: pd.DataFrame | None = None
        self.itemChanged.connect(self._mark)
        QShortcut(QKeySequence.StandardKey.Paste, self,
                  activated=self.paste_block)

    def fill(self, frame: pd.DataFrame) -> None:
        self._original = frame.reset_index(drop=True).copy()
        self._shown: dict[tuple[int, int], str] = {}  # text as displayed
        self.blockSignals(True)
        self.setRowCount(len(frame))
        for r, (_i, row) in enumerate(self._original.iterrows()):
            for c, key in enumerate(self.keys):
                v = row[key]
                text = "" if (isinstance(v, float) and math.isnan(v)) else (
                    f"{v:.7g}" if isinstance(v, float) else str(v))
                it = QTableWidgetItem(text)
                self._shown[(r, c)] = text
                if key not in self.editable:
                    it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.setItem(r, c, it)
        self.blockSignals(False)
        self.resizeColumnsToContents()

    def frame(self) -> pd.DataFrame:
        """The table as edited (numbers parsed; blanks stay NaN)."""
        out = self._original.copy()
        for r in range(self.rowCount()):
            for c, key in enumerate(self.keys):
                if key not in self.editable:
                    continue
                text = self.item(r, c).text().strip()
                if text == self._shown.get((r, c)):
                    continue  # untouched: keep the exact original value
                if key in se.COORD_COLUMNS or key in ("easting", "northing"):
                    try:
                        out.at[r, key] = float(text) if text else math.nan
                    except ValueError:
                        out.at[r, key] = text  # reported by apply
                else:
                    out.at[r, key] = text
        return out

    def _mark(self, item: QTableWidgetItem) -> None:
        if self._original is None:
            return
        key = self.keys[item.column()]
        old = self._original.at[item.row(), key]
        new = item.text().strip()
        if new == self._shown.get((item.row(), item.column())):
            same = True
        elif isinstance(old, float):
            try:
                same = (math.isnan(old) and not new) or (
                    new and abs(float(new) - old) <= 1e-9 * max(1, abs(old)))
            except ValueError:
                same = False
        else:
            same = str(old) == new
        item.setBackground(QBrush() if same else QBrush(_CHANGED))
        parent = self.parent()
        while parent is not None and not isinstance(parent,
                                                    StationEditorDialog):
            parent = parent.parent()
        if parent is not None:
            parent.update_summary()

    def paste_block(self) -> None:
        """Paste a spreadsheet block at the current cell."""
        cells = se.parse_pasted(QApplication.clipboard().text())
        if not cells:
            return
        r0 = max(self.currentRow(), 0)
        c0 = max(self.currentColumn(), 0)
        for dr, row in enumerate(cells):
            for dc, text in enumerate(row):
                r, c = r0 + dr, c0 + dc
                if r >= self.rowCount() or c >= self.columnCount():
                    continue
                if self.keys[c] in self.editable:
                    self.item(r, c).setText(text)

    def selected_rows(self) -> list[int]:
        return sorted({i.row() for i in self.selectedIndexes()})


class StationEditorDialog(QDialog):
    """Edit station coordinates, names, lines and header metadata."""

    def __init__(self, sites, lines: dict | None = None, *,
                 tab: str = "coordinates", parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Station Editor")
        self.resize(900, 560)
        self._sites = sites
        self._lines = dict(lines or {})
        self._base = se.station_table(sites, self._lines)
        self.changes: dict = {}  # filled on Apply
        v = QVBoxLayout(self)
        self._tabs = QTabWidget()
        v.addWidget(self._tabs, 1)
        self._tabs.addTab(self._build_coords(), "Coordinates")
        self._tabs.addTab(self._build_names(), "Names")
        self._tabs.addTab(self._build_lines(), "Lines")
        self._tabs.addTab(self._build_header(), "Header")
        self._tabs.setCurrentIndex(TABS.index(tab) if tab in TABS else 0)
        bottom = QHBoxLayout()
        self._summary = QLabel("")
        self._summary.setObjectName("InfoLabel")
        bottom.addWidget(self._summary, 1)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel)
        self._btn_apply = bb.addButton("Apply",
                                       QDialogButtonBox.ButtonRole.AcceptRole)
        self._btn_apply.clicked.connect(self._on_apply)
        bb.rejected.connect(self.reject)
        bottom.addWidget(bb)
        v.addLayout(bottom)
        self.update_summary()

    # ── coordinates ───────────────────────────────────────────────────
    def _build_coords(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        row = QHBoxLayout()
        row.addWidget(QLabel("Coordinates:"))
        self._crs = QComboBox()
        self._crs.addItem("Geographic (lat / lon, WGS84)", 0)
        lat = self._base["lat"].dropna()
        lon = self._base["lon"].dropna()
        if len(lat) and len(lon):
            utm = se.utm_epsg(float(lat.median()), float(lon.median()))
            self._crs.addItem(f"UTM (EPSG:{utm}, suggested)", utm)
        self._crs.addItem("Other EPSG…", -1)
        self._crs.currentIndexChanged.connect(self._on_crs)
        row.addWidget(self._crs)
        self._epsg = QSpinBox()
        self._epsg.setRange(1024, 999999)
        self._epsg.setPrefix("EPSG:")
        self._epsg.setValue(32650)
        self._epsg.setVisible(False)
        self._epsg.editingFinished.connect(lambda: self._on_crs(0))
        row.addWidget(self._epsg)
        row.addStretch(1)
        b = QPushButton("Import CSV…")
        b.setToolTip("A table with a station column and lat/lon/elev or "
                     "easting/northing columns")
        b.clicked.connect(self._import_csv)
        row.addWidget(b)
        v.addLayout(row)
        self._coords = _Table([("station", "Station"), ("line", "Line"),
                               ("lat", "Latitude °"), ("lon", "Longitude °"),
                               ("elev", "Elevation (m)")],
                              editable={"lat", "lon", "elev"})
        self._coords.fill(self._base)
        v.addWidget(self._coords, 1)
        tip = QLabel("Tip: copy a block of cells in Excel and press Ctrl+V "
                     "on the first cell to fill them all.")
        tip.setObjectName("InfoLabel")
        v.addWidget(tip)
        self._projected_epsg = 0
        return w

    def _current_epsg(self) -> int:
        data = self._crs.currentData()
        return int(self._epsg.value()) if data == -1 else int(data or 0)

    def _on_crs(self, _i: int) -> None:
        self._epsg.setVisible(self._crs.currentData() == -1)
        geo = self._coords_geographic()  # keep edits across a switch
        epsg = self._current_epsg()
        try:
            if epsg:
                shown = se.project_table(geo, epsg)
                cols = [("station", "Station"), ("line", "Line"),
                        ("easting", "Easting (m)"),
                        ("northing", "Northing (m)"),
                        ("elev", "Elevation (m)")]
                editable = {"easting", "northing", "elev"}
            else:
                shown, cols = geo, [("station", "Station"), ("line", "Line"),
                                    ("lat", "Latitude °"),
                                    ("lon", "Longitude °"),
                                    ("elev", "Elevation (m)")]
                editable = {"lat", "lon", "elev"}
        except Exception as exc:
            QMessageBox.warning(self, "Coordinates", str(exc))
            return
        self._projected_epsg = epsg
        base = (se.project_table(self._base, epsg) if epsg else self._base)
        self._coords.keys = [k for k, _l in cols]
        self._coords.editable = editable
        self._coords.setHorizontalHeaderLabels([lbl for _k, lbl in cols])
        self._coords.fill(base)
        # re-apply edits made before the switch
        self._coords.blockSignals(True)
        for r in range(len(shown)):
            for c, key in enumerate(self._coords.keys):
                if key in editable:
                    v = shown.at[r, key]
                    self._coords.item(r, c).setText(
                        "" if v != v else f"{float(v):.7g}")
        self._coords.blockSignals(False)
        for r in range(self._coords.rowCount()):
            for c in range(self._coords.columnCount()):
                self._coords._mark(self._coords.item(r, c))

    def _coords_geographic(self) -> pd.DataFrame:
        f = self._coords.frame()
        if self._projected_epsg:
            base = se.project_table(self._base, self._projected_epsg)
            moved = ~(f["easting"].astype(float).eq(base["easting"])
                      & f["northing"].astype(float).eq(base["northing"]))
            geo = se.unproject_table(f, self._projected_epsg)
            # untouched rows keep their exact lat / lon (no round trip)
            f["lat"] = geo["lat"].where(moved, self._base["lat"])
            f["lon"] = geo["lon"].where(moved, self._base["lon"])
        return f[["station", "line", "lat", "lon", "elev"]]

    def _import_csv(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Import coordinates", "", "Tables (*.csv *.txt);;All (*)")
        if path:
            self.import_coordinates(path)

    def import_coordinates(self, source) -> int:
        """Fill the table from a CSV / DataFrame (matched by station);
        returns the number of stations matched."""
        t = source if isinstance(source, pd.DataFrame) else pd.read_csv(
            source, sep=None, engine="python")
        cols = {c.lower().strip(): c for c in t.columns}
        key = next((cols[c] for c in ("station", "name", "site", "id")
                    if c in cols), None)
        if key is None:
            raise ValueError("no station / name / site / id column")
        pick = {k: next((cols[a] for a in al if a in cols), None) for k, al in
                (("lat", ("lat", "latitude")),
                 ("lon", ("lon", "long", "longitude")),
                 ("elev", ("elev", "elevation", "z", "altitude")),
                 ("easting", ("easting", "east", "x")),
                 ("northing", ("northing", "north", "y")))}
        if pick["lat"] is None and pick["easting"] is not None:
            epsg = self._current_epsg()
            if not epsg:
                raise ValueError("the table has easting / northing: choose "
                                 "its projection (UTM / EPSG) first")
            t = se.unproject_table(t.rename(columns={
                pick["easting"]: "easting", pick["northing"]: "northing"}),
                epsg)
            pick["lat"], pick["lon"] = "lat", "lon"
        by = {str(s).strip().casefold(): i for i, s in
              enumerate(self._base["station"])}
        geo = self._coords_geographic()
        matched = 0
        for _i, row in t.iterrows():
            r = by.get(str(row[key]).strip().casefold())
            if r is None:
                continue
            matched += 1
            for col in ("lat", "lon", "elev"):
                if pick[col] is not None and pd.notna(row[pick[col]]):
                    geo.at[r, col] = float(row[pick[col]])
        self._set_geographic(geo)
        return matched

    def _set_geographic(self, geo: pd.DataFrame) -> None:
        epsg = self._projected_epsg
        shown = se.project_table(geo, epsg) if epsg else geo
        for r in range(len(shown)):
            for c, key in enumerate(self._coords.keys):
                if key in self._coords.editable:
                    v = shown.at[r, key]
                    self._coords.item(r, c).setText(
                        "" if v != v else f"{float(v):.7g}")

    # ── names ─────────────────────────────────────────────────────────
    def _build_names(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        form = QFormLayout()
        row = QHBoxLayout()
        self._prefix = QLineEdit()
        self._prefix.setPlaceholderText("prefix")
        self._suffix = QLineEdit()
        self._suffix.setPlaceholderText("suffix")
        row.addWidget(self._prefix)
        row.addWidget(self._suffix)
        form.addRow("Add:", row)
        row = QHBoxLayout()
        self._find = QLineEdit()
        self._find.setPlaceholderText("find")
        self._replace = QLineEdit()
        self._replace.setPlaceholderText("replace with")
        self._regex = QCheckBox("regular expression")
        row.addWidget(self._find)
        row.addWidget(self._replace)
        row.addWidget(self._regex)
        form.addRow("Replace:", row)
        row = QHBoxLayout()
        self._case = QComboBox()
        for data, label in (("keep", "Keep case"), ("upper", "UPPER"),
                            ("lower", "lower")):
            self._case.addItem(label, data)
        self._pad = QSpinBox()
        self._pad.setRange(0, 6)
        self._pad.setSpecialValueText("no padding")
        self._pad.setSuffix(" digits")
        self._pad.setToolTip("Zero-pad the station number (L1-3 → L1-003)")
        row.addWidget(self._case)
        row.addWidget(self._pad)
        btn = QPushButton("Preview rules")
        btn.clicked.connect(self._preview_names)
        row.addWidget(btn)
        row.addStretch(1)
        form.addRow("Format:", row)
        v.addLayout(form)
        self._names = _Table([("station", "Current name"),
                              ("new", "New name")], editable={"new"})
        frame = self._base[["station"]].copy()
        frame["new"] = frame["station"]
        self._names.fill(frame)
        v.addWidget(self._names, 1)
        self._names_problems = QLabel("")
        self._names_problems.setStyleSheet("color: #c0392b;")
        v.addWidget(self._names_problems)
        return w

    def _preview_names(self) -> None:
        new = se.rename_preview(
            list(self._base["station"]), prefix=self._prefix.text(),
            suffix=self._suffix.text(), find=self._find.text(),
            replace=self._replace.text(), regex=self._regex.isChecked(),
            case=self._case.currentData(), pad=self._pad.value())
        for r, name in enumerate(new):
            self._names.item(r, 1).setText(name)

    # ── lines ─────────────────────────────────────────────────────────
    def _build_lines(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        row = QHBoxLayout()
        b = QPushButton("Auto-detect from names")
        b.setToolTip("L18-001, L18-002 → line L18 (numeric prefixes get 'L')")
        b.clicked.connect(self._auto_lines)
        row.addWidget(b)
        row.addSpacing(16)
        row.addWidget(QLabel("Rename line"))
        self._line_from = QComboBox()
        self._line_to = QLineEdit()
        self._line_to.setPlaceholderText("new line name")
        btn = QPushButton("Rename")
        btn.clicked.connect(self._rename_line)
        row.addWidget(self._line_from)
        row.addWidget(QLabel("to"))
        row.addWidget(self._line_to)
        row.addWidget(btn)
        row.addStretch(1)
        v.addLayout(row)
        self._lines_tbl = _Table([("station", "Station"), ("line", "Line")],
                                 editable={"line"})
        self._lines_tbl.fill(self._base[["station", "line"]])
        self._refresh_line_names()
        v.addWidget(self._lines_tbl, 1)
        return w

    def _refresh_line_names(self) -> None:
        names = sorted({self._lines_tbl.item(r, 1).text()
                        for r in range(self._lines_tbl.rowCount())} - {""})
        self._line_from.clear()
        self._line_from.addItems(names)

    def _auto_lines(self) -> None:
        found = se.detect_lines(self._base["station"])
        for r, st in enumerate(self._base["station"]):
            if st in found:
                self._lines_tbl.item(r, 1).setText(found[st])
        self._refresh_line_names()

    def _rename_line(self) -> None:
        old, new = self._line_from.currentText(), self._line_to.text().strip()
        if not old or not new:
            return
        for r in range(self._lines_tbl.rowCount()):
            it = self._lines_tbl.item(r, 1)
            if it.text() == old:
                it.setText(new)
        self._refresh_line_names()

    # ── header ────────────────────────────────────────────────────────
    def _build_header(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        row = QHBoxLayout()
        row.addWidget(QLabel("Set"))
        self._fill_field = QComboBox()
        for label, key in se.HEAD_FIELDS:
            self._fill_field.addItem(label, key)
        self._fill_value = QLineEdit()
        self._fill_value.setPlaceholderText("value")
        self._fill_scope = QComboBox()
        self._fill_scope.addItem("for all stations", "all")
        self._fill_scope.addItem("for selected rows", "selected")
        btn = QPushButton("Fill")
        btn.clicked.connect(self._fill_header)
        for x in (self._fill_field, QLabel("to"), self._fill_value,
                  self._fill_scope, btn):
            row.addWidget(x)
        row.addStretch(1)
        v.addLayout(row)
        cols = [("station", "Station")] + [(k, lbl) for lbl, k in
                                           se.HEAD_FIELDS]
        self._header = _Table(cols, editable={k for _l, k in
                                              se.HEAD_FIELDS})
        self._header.fill(self._base[[k for k, _l in cols]])
        v.addWidget(self._header, 1)
        return w

    def _fill_header(self) -> None:
        key = self._fill_field.currentData()
        c = self._header.keys.index(key)
        rows = (self._header.selected_rows()
                if self._fill_scope.currentData() == "selected"
                else range(self._header.rowCount()))
        for r in rows:
            self._header.item(r, c).setText(self._fill_value.text())

    # ── result ────────────────────────────────────────────────────────
    def edited_table(self) -> pd.DataFrame:
        """The base table with every tab's edits folded in."""
        out = self._base.copy()
        geo = self._coords_geographic()
        for col in ("lat", "lon", "elev"):
            out[col] = geo[col].to_numpy()
        out["station"] = self._names.frame()["new"].astype(str).str.strip()
        out["line"] = self._lines_tbl.frame()["line"].astype(str).str.strip()
        head = self._header.frame()
        for _label, key in se.HEAD_FIELDS:
            out[key] = head[key].to_numpy()
        return out

    def compute_changes(self) -> dict:
        return se.table_changes(self._base, self.edited_table())

    def update_summary(self) -> None:
        try:
            changes = self.compute_changes()
            problems = se.check_names(self._names.frame()["new"])
        except Exception as exc:  # a half-typed number, a bad EPSG
            self._summary.setText(f"⚠ {exc}")
            return
        self._names_problems.setText("⚠ " + "; ".join(problems)
                                     if problems else "")
        n_cells = sum(len(d) for d in changes.values())
        self._summary.setText(
            f"{n_cells} change(s) in {len(changes)} station(s)" if changes
            else "No changes yet")
        self._btn_apply.setEnabled(bool(changes) and not problems)

    def _on_apply(self) -> None:
        try:
            self.changes = self.compute_changes()
        except Exception as exc:
            QMessageBox.warning(self, "Station Editor", str(exc))
            return
        problems = se.check_names(self.edited_table()["station"])
        if problems:
            QMessageBox.warning(self, "Station Editor", "; ".join(problems))
            return
        self.accept()

    def change_label(self) -> str:
        """A short undo label for the changes ("Edit coordinates, names")."""
        kinds = []
        cols = {k for d in self.changes.values() for k in d}
        if cols & {"lat", "lon", "elev"}:
            kinds.append("coordinates")
        if "name" in cols:
            kinds.append("names")
        if "line" in cols:
            kinds.append("lines")
        if cols & {k for _l, k in se.HEAD_FIELDS}:
            kinds.append("header")
        return "Edit stations: " + (", ".join(kinds) or "no change")


__all__ = ["StationEditorDialog", "TABS"]
