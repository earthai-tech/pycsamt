# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SurveyOpsDialog — Edit ▸ Frequencies / Tensor whole-survey operations.

Left: the operation (trim band, regrid, decimate, align to a common grid,
remove duplicates, fill gaps, rotate, rotate to strike).  Right: its
settings and a live preview -- per station, the number of frequencies,
the band and the masked points before -> after (and the strike angle when
rotating to strike).  **Apply** hands the result back as one undoable
edit.  The operations themselves are in
:mod:`pycsamt.app.desktop.controllers.survey_ops`.
"""

from __future__ import annotations

from dataclasses import replace

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QSplitter,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers import survey_ops as so
from pycsamt.app.desktop.windows.inversion.forms import SettingsForm

_KEY = Qt.ItemDataRole.UserRole
_LESS = QColor(255, 214, 102, 110)   # fewer frequencies / more masked
_MORE = QColor(116, 192, 252, 90)    # more frequencies


class SurveyOpsDialog(QDialog):
    """Preview and apply a whole-survey frequency or tensor operation."""

    def __init__(self, sites, *, op: str = "trim", parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Frequency & Tensor Tools")
        self.resize(980, 600)
        self._sites = sites
        self.result = None  # the operated survey, set on Apply
        self._forms: dict[str, SettingsForm] = {}
        names = [str(s.name) for s in sites]
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(300)
        self._timer.timeout.connect(self.preview)

        v = QVBoxLayout(self)
        split = QSplitter(Qt.Orientation.Horizontal)
        self._list = QListWidget()
        self._list.setMaximumWidth(250)
        self._stack = QStackedWidget()
        for group, title in (("frequencies", "FREQUENCIES"),
                             ("tensor", "TENSOR")):
            head = QListWidgetItem(title)
            head.setFlags(Qt.ItemFlag.NoItemFlags)
            self._list.addItem(head)
            for o in so.OPS:
                if o.menu != group:
                    continue
                it = QListWidgetItem("   " + o.label)
                it.setData(_KEY, o.key)
                it.setToolTip(o.help)
                self._list.addItem(it)
                fields = [replace(f, choices=(("", "—"),) + tuple(
                    (n, n) for n in names)) if f.key == "ref" else f
                    for f in o.fields]
                page = QWidget()
                pv = QVBoxLayout(page)
                pv.setContentsMargins(0, 0, 0, 0)
                lbl = QLabel(f"<b>{o.label}</b><br><span style='color:"
                             f"#6b7280'>{o.help}</span>")
                lbl.setWordWrap(True)
                pv.addWidget(lbl)
                form = SettingsForm(fields)
                form.changed.connect(self._timer.start)
                pv.addWidget(form)
                self._forms[o.key] = form
                self._stack.addWidget(page)
        self._list.currentItemChanged.connect(self._on_op)
        split.addWidget(self._list)

        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(6, 0, 0, 0)
        rv.addWidget(self._stack)
        self._table = QTableWidget()
        self._table.setAlternatingRowColors(True)
        self._table.verticalHeader().setVisible(False)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        rv.addWidget(self._table, 1)
        self._status = QLabel("")
        self._status.setObjectName("InfoLabel")
        self._status.setWordWrap(True)
        rv.addWidget(self._status)
        split.addWidget(right)
        split.setStretchFactor(1, 1)
        v.addWidget(split, 1)

        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel)
        self._btn_apply = bb.addButton("Apply to survey",
                                       QDialogButtonBox.ButtonRole.AcceptRole)
        self._btn_apply.clicked.connect(self._on_apply)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)
        self.select_op(op)

    # ── operation ─────────────────────────────────────────────────────
    @property
    def op_key(self) -> str:
        it = self._list.currentItem()
        return it.data(_KEY) if it is not None else ""

    def select_op(self, key: str) -> None:
        for r in range(self._list.count()):
            if self._list.item(r).data(_KEY) == key:
                self._list.setCurrentRow(r)
                return

    def _on_op(self, cur, _prev) -> None:
        if cur is None or cur.data(_KEY) is None:
            return
        keys = [o.key for o in so.OPS]
        self._stack.setCurrentIndex(keys.index(cur.data(_KEY)))
        self.preview()

    def values(self) -> dict:
        return self._forms[self.op_key].values()

    def set_values(self, values: dict) -> None:
        self._forms[self.op_key].set_values(values)
        self._timer.stop()
        self.preview()

    # ── preview / apply ───────────────────────────────────────────────
    def preview(self) -> bool:
        """Run the operation on a copy and show the per-station effect."""
        self._pending = None
        self._btn_apply.setEnabled(False)
        try:
            out = so.run(self.op_key, self._sites, self.values())
            table = so.preview_table(self._sites, out,
                                     with_strike=self.op_key == "strike")
        except Exception as exc:
            self._table.setRowCount(0)
            self._status.setText(f"⚠ {exc}")
            return False
        self._pending = out
        self._fill(table)
        n0, n1 = int(table["n_before"].sum()), int(table["n_after"].sum())
        m0 = int(table["masked_before"].sum())
        m1 = int(table["masked_after"].sum())
        changed = (table["n_before"] != table["n_after"]) | (
            table["band_before"] != table["band_after"]) | (
            table["masked_before"] != table["masked_after"])
        rotated = self.op_key in ("rotate", "strike")
        self._status.setText(
            f"{len(table)} stations · frequencies {n0} → {n1} · masked "
            f"points {m0} → {m1}"
            + (" · every tensor rotated" if rotated else "")
            + ("" if (changed.any() or rotated)
               else " — nothing would change"))
        empty = table.index[table["n_after"] == 0].tolist()
        if empty:
            self._status.setText(self._status.text() + f" · ⚠ "
                                 f"{len(empty)} station(s) left with no "
                                 "frequency")
        self._btn_apply.setEnabled(bool(changed.any() or rotated)
                                   and not empty)
        return True

    def _fill(self, table) -> None:
        cols = [("station", "Station"), ("n_before", "Freqs before"),
                ("n_after", "Freqs after"), ("band_before", "Band before"),
                ("band_after", "Band after"),
                ("masked_before", "Masked before"),
                ("masked_after", "Masked after")]
        if "strike" in table:
            cols.append(("strike", "Strike °"))
        self._table.setColumnCount(len(cols))
        self._table.setHorizontalHeaderLabels([lbl for _k, lbl in cols])
        self._table.setRowCount(len(table))
        for r, (_i, row) in enumerate(table.iterrows()):
            for c, (key, _l) in enumerate(cols):
                val = row[key]
                text = f"{val:.1f}" if key == "strike" and val == val else (
                    "—" if key == "strike" else str(val))
                it = QTableWidgetItem(text)
                if key == "n_after" and row["n_after"] != row["n_before"]:
                    it.setBackground(QBrush(_LESS if row["n_after"] <
                                            row["n_before"] else _MORE))
                if key == "masked_after" and row["masked_after"] > \
                        row["masked_before"]:
                    it.setBackground(QBrush(_LESS))
                self._table.setItem(r, c, it)
        self._table.resizeColumnsToContents()

    def _on_apply(self) -> None:
        if getattr(self, "_pending", None) is None and not self.preview():
            return
        self.result = self._pending
        self.accept()

    def label(self) -> str:
        return so.op(self.op_key).label


__all__ = ["SurveyOpsDialog"]
