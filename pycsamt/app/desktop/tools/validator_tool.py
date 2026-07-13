# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
EDIValidatorDialog — per-station data-quality checklist + Station Manager.

Quality checks
--------------
For every station in the loaded survey:
  ✓ Station name present
  ✓ Coordinates (lat / lon) present
  ✓ Z-tensor has at least 2 non-NaN components (xy, yx)
  ✓ Impedance error estimates (z_err) present
  ✓ Frequency count ≥ 1
  ✓ Frequency range (min T … max T)
  ⚠ Any frequency with all-NaN Z row

Colour coding
-------------
  Green  — all checks passed
  Yellow — minor warnings (no error estimates, few frequencies)
  Red    — critical failure (no Z, no coordinates)
  Gray   — excluded by user (will be removed on Apply)

Station management toolbar
--------------------------
  ✏ Rename        — double-click station name or use button / context menu
  🗑 Delete        — permanently remove from working list (with confirmation)
  ⊘ Exclude       — soft-hide; station is kept in display but filtered on Apply
  ↑ / ↓ Move      — reorder stations along profile
  Export CSV      — save current table as CSV
  Apply to Survey — push changes (renames + filtered list) to MainWindow

Signal
------
``modified_sites`` attribute is populated when the user clicks *Apply to Survey*
and ``exec()`` returns ``Accepted``.  MainWindow reads it to update the survey.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QIcon
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QMenu,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

# ── Colours ───────────────────────────────────────────────────────────────────
_GREEN = QColor("#c8e6c9")
_YELLOW = QColor("#fff9c4")
_RED = QColor("#ffcdd2")
_EXCLUDED = QColor("#e0e0e0")

# Foreground colours — always dark so text is visible over light pastels
# regardless of whether the OS is in dark or light mode.
_FG_NORMAL = QColor("#1a1a1a")
_FG_EXCLUDED = QColor("#666666")
_ICONS = Path(__file__).resolve().parents[1] / "resources" / "icons"


def _icon(name: str) -> QIcon:
    for candidate in (name, f"{name}.svg", f"{name}.png"):
        path = _ICONS / candidate
        if path.exists():
            return QIcon(str(path))
    return QIcon()


# ── Column layout ─────────────────────────────────────────────────────────────
_COLS = [
    "Station",
    "Lat",
    "Lon",
    "Z ok",
    "Errors",
    "N freq",
    "T min (s)",
    "T max (s)",
    "NaN rows",
    "Status",
]
_COL_STATION = 0
_COL_STATUS = len(_COLS) - 1


# ── QC helper ─────────────────────────────────────────────────────────────────


def _check_site(ed) -> dict:
    """Run all QC checks on one EDI/site object and return a result dict."""
    result = {c: "—" for c in _COLS}
    result["Status"] = "PASS"

    # Station name
    name = (
        getattr(ed, "station", None)
        or getattr(ed, "Site", {}).get("Name", "")
        or getattr(ed, "id", "?")
    )
    result["Station"] = str(name)

    # Coordinates
    lat = getattr(ed, "lat", None) or getattr(ed, "latitude", None)
    lon = getattr(ed, "lon", None) or getattr(ed, "longitude", None)
    if lat is None:
        try:
            lat = ed.Header.lat
        except Exception:
            pass
    if lon is None:
        try:
            lon = ed.Header.lon
        except Exception:
            pass

    result["Lat"] = f"{float(lat):.4f}" if lat is not None else "MISSING"
    result["Lon"] = f"{float(lon):.4f}" if lon is not None else "MISSING"
    if lat is None or lon is None:
        result["Status"] = "WARN"

    # Z tensor
    Z_obj = getattr(ed, "Z", None) or getattr(
        getattr(ed, "edi", None), "Z", None
    )
    z = z_err = freqs = None
    if Z_obj is not None:
        z = getattr(Z_obj, "z", None)
        z_err = getattr(Z_obj, "z_err", None)
        freqs = getattr(Z_obj, "freq", None)

    if z is None or (hasattr(z, "size") and z.size == 0):
        result["Z ok"] = "NO"
        result["Status"] = "FAIL"
    else:
        z_arr = np.asarray(z)
        xy_ok = (
            np.any(np.isfinite(z_arr[:, 0, 1])) if z_arr.ndim == 3 else False
        )
        yx_ok = (
            np.any(np.isfinite(z_arr[:, 1, 0])) if z_arr.ndim == 3 else False
        )
        result["Z ok"] = "YES" if (xy_ok and yx_ok) else "PARTIAL"
        if not (xy_ok and yx_ok) and result["Status"] == "PASS":
            result["Status"] = "WARN"
        if z_arr.ndim == 3:
            nan_rows = int(
                np.sum(
                    np.all(
                        ~np.isfinite(z_arr.reshape(z_arr.shape[0], -1)),
                        axis=1,
                    )
                )
            )
            result["NaN rows"] = str(nan_rows)
            if nan_rows > 0 and result["Status"] == "PASS":
                result["Status"] = "WARN"

    # Error estimates
    if z_err is None or (
        hasattr(z_err, "size") and not np.any(np.isfinite(z_err))
    ):
        result["Errors"] = "NO"
        if result["Status"] == "PASS":
            result["Status"] = "WARN"
    else:
        result["Errors"] = "YES"

    # Frequency
    if freqs is not None and hasattr(freqs, "__len__") and len(freqs) > 0:
        fa = np.asarray(freqs, dtype=float)
        fa = fa[fa > 0]
        result["N freq"] = str(len(fa))
        if fa.size > 0:
            p = 1.0 / fa
            result["T min (s)"] = f"{p.min():.4g}"
            result["T max (s)"] = f"{p.max():.4g}"
        if len(fa) < 3 and result["Status"] == "PASS":
            result["Status"] = "WARN"
    else:
        result["N freq"] = "0"
        if result["Status"] == "PASS":
            result["Status"] = "FAIL"

    return result


# ── Dialog ────────────────────────────────────────────────────────────────────


class EDIValidatorDialog(QDialog):
    """
    Per-station QC checklist + interactive Station Manager.

    Attributes
    ----------
    modified_sites : list or None
        Populated when the user clicks *Apply to Survey* and the dialog
        accepts.  Contains the filtered / renamed list of site objects.

    Signals
    -------
    open_recompute_requested
        Emitted when the user clicks *Recompute EDIs…* inside this dialog.
        MainWindow connects this to ``_open_recompute()``.

    Parameters
    ----------
    sites : any
        Loaded survey (from MainWindow._controller.sites).
    parent : QWidget, optional
    """

    open_recompute_requested = Signal()

    def __init__(self, sites, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("EDI Validator & Station Manager")
        self.setMinimumSize(1020, 580)

        self._sites = sites
        self._items: list = []  # raw EDI objects, parallel to rows
        self._results: list = []  # QC result dicts, parallel to rows
        self._excluded: set = set()  # row indices soft-excluded
        self._populating: bool = False
        self.modified_sites = None

        self._build_ui()
        self._run_checks()

    # ── Build UI ──────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(6)
        root.setContentsMargins(8, 8, 8, 8)

        # ── Summary + legend ──────────────────────────────────────────────
        self._summary_lbl = QLabel("Running checks…")
        self._summary_lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
        root.addWidget(self._summary_lbl)

        leg = QHBoxLayout()
        for color, fg, text in (
            (_GREEN, _FG_NORMAL, "Pass"),
            (_YELLOW, _FG_NORMAL, "Warning"),
            (_RED, _FG_NORMAL, "Fail"),
            (_EXCLUDED, _FG_EXCLUDED, "Excluded"),
        ):
            lbl = QLabel(f"  {text}  ")
            lbl.setAutoFillBackground(True)
            p = lbl.palette()
            p.setColor(lbl.backgroundRole(), color)
            p.setColor(lbl.foregroundRole(), fg)
            lbl.setPalette(p)
            leg.addWidget(lbl)
        leg.addStretch()
        root.addLayout(leg)

        # ── Management toolbar ────────────────────────────────────────────
        tb = QHBoxLayout()
        tb.setSpacing(4)

        def _tbtn(label: str, tip: str, slot) -> QPushButton:
            b = QPushButton(label)
            b.setToolTip(tip)
            b.setFixedHeight(26)
            b.clicked.connect(slot)
            return b

        tb.addWidget(
            _tbtn(
                "✏  Rename",
                "Rename selected station  (or double-click name)",
                self._on_rename,
            )
        )
        tb.addWidget(
            _tbtn(
                "🗑  Delete",
                "Remove selected station(s) from the survey",
                self._on_delete,
            )
        )
        tb.addSpacing(6)
        self._excl_btn = _tbtn(
            "⊘  Exclude",
            "Toggle exclude / include for selected stations",
            self._on_toggle_exclude,
        )
        tb.addWidget(self._excl_btn)
        tb.addSpacing(6)
        tb.addWidget(
            _tbtn(
                "↑  Move Up",
                "Move selected station up in profile order",
                self._on_move_up,
            )
        )
        tb.addWidget(
            _tbtn(
                "↓  Move Down",
                "Move selected station down in profile order",
                self._on_move_down,
            )
        )
        tb.addStretch()
        tb.addWidget(
            _tbtn("Export CSV…", "Save current table as CSV", self._on_export)
        )
        self._apply_btn = _tbtn(
            "Apply to Survey",
            "Push renames + filtered list back to the loaded survey",
            self._on_apply,
        )
        self._apply_btn.setEnabled(False)
        self._apply_btn.setStyleSheet("font-weight: bold;")
        tb.addWidget(self._apply_btn)

        root.addLayout(tb)

        # ── Table ─────────────────────────────────────────────────────────
        self._table = QTableWidget(0, len(_COLS))
        self._table.setHorizontalHeaderLabels(_COLS)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(_COL_STATION, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(
            _COL_STATUS, QHeaderView.ResizeMode.ResizeToContents
        )
        self._table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self._table.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection
        )
        self._table.setEditTriggers(QTableWidget.EditTrigger.DoubleClicked)
        self._table.setAlternatingRowColors(False)
        self._table.itemChanged.connect(self._on_item_changed)
        self._table.setContextMenuPolicy(
            Qt.ContextMenuPolicy.CustomContextMenu
        )
        self._table.customContextMenuRequested.connect(
            self._show_context_menu
        )
        root.addWidget(self._table, stretch=1)

        # ── Bottom buttons ────────────────────────────────────────────────
        bottom = QHBoxLayout()
        bottom.setSpacing(6)
        btn_recompute = QPushButton("Recompute EDIs…")
        btn_recompute.setIcon(_icon("recompute"))
        btn_recompute.setToolTip(
            "Open the Recompute EDIs panel to rotate, filter and rewrite the survey EDIs.\n"
            "Recomputed stations are marked with ◈ in the main station list."
        )
        btn_recompute.setObjectName("RecomputeShortcutBtn")
        btn_recompute.clicked.connect(self.open_recompute_requested)
        bottom.addWidget(btn_recompute)
        bottom.addStretch()
        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        bottom.addWidget(box)
        root.addLayout(bottom)

    # ── QC run ────────────────────────────────────────────────────────────────

    def _run_checks(self) -> None:
        if self._sites is None:
            self._summary_lbl.setText("No survey loaded.")
            return
        try:
            from pycsamt.emtools._core import (
                _iter_items,
                _unwrap,
            )

            items = list(_iter_items(self._sites))
        except Exception:
            try:
                items = (
                    list(self._sites)
                    if hasattr(self._sites, "__iter__")
                    else [self._sites]
                )
            except Exception:
                self._summary_lbl.setText("Cannot iterate over loaded data.")
                return

        self._items = []
        self._results = []
        for ed in items:
            try:
                ed = (
                    _unwrap(ed)
                    if "pycsamt" in str(type(ed).__module__)
                    else ed
                )
            except Exception:
                pass
            self._items.append(ed)
            try:
                self._results.append(_check_site(ed))
            except Exception as exc:
                self._results.append(
                    {
                        "Station": "?",
                        "Status": "FAIL",
                        **{
                            c: "—"
                            for c in _COLS
                            if c not in ("Station", "Status")
                        },
                        "NaN rows": str(exc),
                    }
                )

        self._excluded = set()
        self._fill_table()
        self._apply_btn.setEnabled(bool(self._items))

    # ── Table fill ────────────────────────────────────────────────────────────

    def _fill_table(self) -> None:
        self._populating = True
        self._table.setRowCount(0)
        n_pass = n_warn = n_fail = n_excl = 0

        for r, res in enumerate(self._results):
            self._table.insertRow(r)
            excluded = r in self._excluded
            status = res.get("Status", "PASS")

            if excluded:
                bg = _EXCLUDED
                n_excl += 1
            elif status == "PASS":
                bg, n_pass = _GREEN, n_pass + 1
            elif status == "WARN":
                bg, n_warn = _YELLOW, n_warn + 1
            else:
                bg, n_fail = _RED, n_fail + 1

            fg = _FG_EXCLUDED if excluded else _FG_NORMAL
            for c_idx, col in enumerate(_COLS):
                it = QTableWidgetItem(res.get(col, "—"))
                it.setBackground(bg)
                it.setForeground(fg)
                it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                # Only Station column is editable
                if c_idx != _COL_STATION:
                    it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if excluded:
                    f = it.font()
                    f.setStrikeOut(True)
                    it.setFont(f)
                self._table.setItem(r, c_idx, it)

        total = len(self._results)
        excl_note = f"  ·  {n_excl} excluded" if n_excl else ""
        self._summary_lbl.setText(
            f"{total} stations — "
            f"{n_pass} pass  |  {n_warn} warn  |  {n_fail} fail{excl_note}"
        )
        self._populating = False

    # ── Item change (inline rename) ───────────────────────────────────────────

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        if self._populating or item.column() != _COL_STATION:
            return
        row = item.row()
        if 0 <= row < len(self._results):
            self._results[row]["Station"] = item.text().strip()

    # ── Selection helpers ─────────────────────────────────────────────────────

    def _selected_rows(self) -> list[int]:
        return sorted({idx.row() for idx in self._table.selectedIndexes()})

    # ── Toolbar actions ───────────────────────────────────────────────────────

    def _on_rename(self) -> None:
        rows = self._selected_rows()
        if not rows:
            return
        row = rows[0]
        old_name = self._results[row]["Station"]
        new_name, ok = QInputDialog.getText(
            self,
            "Rename Station",
            f"New name for  '{old_name}':",
            text=old_name,
        )
        if ok and new_name.strip() and new_name.strip() != old_name:
            self._results[row]["Station"] = new_name.strip()
            self._fill_table()
            self._table.selectRow(row)

    def _on_delete(self) -> None:
        rows = self._selected_rows()
        if not rows:
            return
        names = ", ".join(self._results[r]["Station"] for r in rows[:4])
        if len(rows) > 4:
            names += f" … (+{len(rows) - 4} more)"
        reply = QMessageBox.question(
            self,
            "Delete stations?",
            f"Permanently remove {len(rows)} station(s) from the working list?\n\n"
            f"  {names}\n\n"
            "This does not affect files on disk.  "
            "Click <b>Apply to Survey</b> to push the change to the loaded survey.",
            QMessageBox.StandardButton.Yes
            | QMessageBox.StandardButton.Cancel,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return

        for r in sorted(rows, reverse=True):
            del self._items[r]
            del self._results[r]

        # Remap excluded indices after deletion
        removed = set(rows)
        new_excluded: set[int] = set()
        for ex in self._excluded:
            if ex in removed:
                continue
            shift = sum(1 for r in rows if r < ex)
            new_excluded.add(ex - shift)
        self._excluded = new_excluded

        self._fill_table()

    def _on_toggle_exclude(self) -> None:
        rows = self._selected_rows()
        if not rows:
            return
        if any(r not in self._excluded for r in rows):
            self._excluded.update(rows)  # exclude all selected
        else:
            self._excluded -= set(rows)  # include all selected
        self._fill_table()
        for r in rows:
            self._table.selectRow(r)

    def _on_move_up(self) -> None:
        rows = self._selected_rows()
        if not rows or rows[0] == 0:
            return
        for r in rows:
            if r == 0:
                continue
            self._items[r], self._items[r - 1] = (
                self._items[r - 1],
                self._items[r],
            )
            self._results[r], self._results[r - 1] = (
                self._results[r - 1],
                self._results[r],
            )
            self._excluded = _swap_excluded(self._excluded, r, r - 1)
        self._fill_table()
        for r in rows:
            if r > 0:
                self._table.selectRow(r - 1)

    def _on_move_down(self) -> None:
        rows = self._selected_rows()
        n = len(self._items)
        if not rows or rows[-1] >= n - 1:
            return
        for r in reversed(rows):
            if r >= n - 1:
                continue
            self._items[r], self._items[r + 1] = (
                self._items[r + 1],
                self._items[r],
            )
            self._results[r], self._results[r + 1] = (
                self._results[r + 1],
                self._results[r],
            )
            self._excluded = _swap_excluded(self._excluded, r, r + 1)
        self._fill_table()
        for r in rows:
            if r < n - 1:
                self._table.selectRow(r + 1)

    # ── Context menu ──────────────────────────────────────────────────────────

    def _show_context_menu(self, pos) -> None:
        rows = self._selected_rows()
        if not rows:
            return
        menu = QMenu(self)
        menu.addAction("✏  Rename…", self._on_rename)
        menu.addAction("🗑  Delete", self._on_delete)
        menu.addSeparator()
        any_included = any(r not in self._excluded for r in rows)
        excl_label = "⊘  Exclude" if any_included else "✓  Include"
        menu.addAction(excl_label, self._on_toggle_exclude)
        menu.addSeparator()
        menu.addAction("↑  Move Up", self._on_move_up)
        menu.addAction("↓  Move Down", self._on_move_down)
        menu.exec(self._table.viewport().mapToGlobal(pos))

    # ── Apply to Survey ───────────────────────────────────────────────────────

    def _on_apply(self) -> None:
        # Push renames to the actual EDI objects
        for row, ed in enumerate(self._items):
            new_name = self._results[row].get("Station", "")
            if not new_name:
                continue
            for attr in ("station", "id", "name"):
                try:
                    if getattr(ed, attr, None) is not None:
                        setattr(ed, attr, new_name)
                except Exception:
                    pass

        # Build filtered list (excluded rows removed)
        self.modified_sites = [
            ed for i, ed in enumerate(self._items) if i not in self._excluded
        ]

        n_kept = len(self.modified_sites)
        n_excl = len(self._excluded)
        self._summary_lbl.setText(
            f"Applied — {n_kept} station(s) kept"
            + (f", {n_excl} excluded / removed." if n_excl else ".")
        )
        self._apply_btn.setEnabled(False)
        self.accept()

    # ── Export CSV ────────────────────────────────────────────────────────────

    def _on_export(self) -> None:
        if not self._results:
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Validation Report",
            str(Path.home() / "survey_validation.csv"),
            "CSV (*.csv)",
        )
        if not path:
            return
        try:
            import csv

            with open(path, "w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=_COLS + ["Excluded"])
                writer.writeheader()
                for r, res in enumerate(self._results):
                    writer.writerow(
                        {
                            **res,
                            "Excluded": "yes"
                            if r in self._excluded
                            else "no",
                        }
                    )
            self._summary_lbl.setText(
                self._summary_lbl.text() + f"  →  saved {path}"
            )
        except Exception as exc:
            self._summary_lbl.setText(f"Export error: {exc}")


# ── Internal helper ───────────────────────────────────────────────────────────


def _swap_excluded(excluded: set, a: int, b: int) -> set:
    """Swap indices a and b inside the excluded set."""
    new = set(excluded)
    a_in = a in excluded
    b_in = b in excluded
    if a_in and not b_in:
        new.discard(a)
        new.add(b)
    elif b_in and not a_in:
        new.discard(b)
        new.add(a)
    return new
