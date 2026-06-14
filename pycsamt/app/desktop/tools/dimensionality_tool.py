# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
DimensionalityDialog — classify each station×frequency as 1D / 2D / 3D.

Upper panel: colour-coded QTableWidget (station × period rows).
Lower panel: bar chart summary of 1D/2D/3D counts per station.

Calls:
  ``pycsamt.emtools.dimensionality.classify_dimensionality`` → DataFrame with
  column ``dim`` (0=1D, 1=2D, 2=3D).

  ``pycsamt.emtools.tensor.plot_dim_confidence_grid`` (optional) → map figure.

Usage
-----
    from pycsamt.app.desktop.tools.dimensionality_tool import DimensionalityDialog
    dlg = DimensionalityDialog(sites, parent=self)
    dlg.exec()
"""
from __future__ import annotations

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

_C_1D = QColor("#c8e6c9")   # green
_C_2D = QColor("#fff9c4")   # yellow
_C_3D = QColor("#ffcdd2")   # red
_DIM_LABEL = {0: "1D", 1: "2D", 2: "3D"}
_DIM_COLOR = {0: _C_1D, 1: _C_2D, 2: _C_3D}


# ── Worker ────────────────────────────────────────────────────────────────────

class _DimWorker(QThread):
    done  = Signal(object, object)   # (DataFrame, Figure|None)
    error = Signal(str)

    def __init__(self, sites, skew_th: float, ellipt_th: float, show_map: bool):
        super().__init__()
        self._sites    = sites
        self._skew_th  = skew_th
        self._ellipt_th = ellipt_th
        self._show_map = show_map

    def run(self):
        try:
            from pycsamt.emtools.dimensionality import classify_dimensionality
            df = classify_dimensionality(
                self._sites,
                skew_th=self._skew_th,
                ellipt_th=self._ellipt_th,
                verbose=0,
                api=False,
            )
            # Unwrap TFBundle / APIFrame if needed
            if hasattr(df, "frame"):
                df = df.frame
            if hasattr(df, "_frame"):
                df = df._frame

            map_fig = None
            if self._show_map:
                try:
                    from pycsamt.emtools.tensor import plot_dim_confidence_grid
                    before = set(plt.get_fignums())
                    ax = plot_dim_confidence_grid(self._sites, verbose=0)
                    after = set(plt.get_fignums())
                    new = after - before
                    map_fig = ax.figure if hasattr(ax, "figure") else (
                        plt.figure(max(new)) if new else None
                    )
                except Exception:
                    pass

            self.done.emit(df, map_fig)
        except Exception as exc:
            self.error.emit(str(exc))


# ── Dialog ────────────────────────────────────────────────────────────────────

class DimensionalityDialog(QDialog):
    """
    1D / 2D / 3D classification table + bar-chart summary.

    Parameters
    ----------
    sites : any
        Loaded survey.
    parent : QWidget, optional
    """

    def __init__(self, sites=None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Dimensionality Classifier")
        self.setMinimumSize(980, 640)
        self._sites  = sites
        self._worker = None
        self._df     = None
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(8)

        # ── Top: controls bar ─────────────────────────────────────────────
        ctrl_bar = QHBoxLayout()

        grp_thresh = QGroupBox("Thresholds")
        form_t = QFormLayout(grp_thresh)
        self._skew_spin = QDoubleSpinBox()
        self._skew_spin.setRange(0.1, 45.0)
        self._skew_spin.setDecimals(1)
        self._skew_spin.setValue(3.0)
        self._skew_spin.setSuffix("°")
        form_t.addRow("Skew β max:", self._skew_spin)
        self._ellipt_spin = QDoubleSpinBox()
        self._ellipt_spin.setRange(0.01, 1.0)
        self._ellipt_spin.setDecimals(2)
        self._ellipt_spin.setValue(0.20)
        form_t.addRow("Ellipticity max:", self._ellipt_spin)
        ctrl_bar.addWidget(grp_thresh)

        grp_opt = QGroupBox("Options")
        form_o = QFormLayout(grp_opt)
        self._map_cb = QCheckBox("Show dimensionality map")
        self._map_cb.setChecked(True)
        form_o.addRow(self._map_cb)
        ctrl_bar.addWidget(grp_opt)

        ctrl_bar.addStretch()

        btn_col = QVBoxLayout()
        self._run_btn = QPushButton("Classify")
        self._run_btn.setFixedWidth(100)
        self._run_btn.clicked.connect(self._on_run)
        btn_col.addWidget(self._run_btn)
        self._status_lbl = QLabel("")
        btn_col.addWidget(self._status_lbl)
        ctrl_bar.addLayout(btn_col)

        root.addLayout(ctrl_bar)

        if self._sites is None:
            self._run_btn.setEnabled(False)
            self._status_lbl.setText("No survey loaded.")

        # ── Legend ────────────────────────────────────────────────────────
        leg = QHBoxLayout()
        for color, text in ((_C_1D, "1D"), (_C_2D, "2D"), (_C_3D, "3D")):
            lbl = QLabel(f"  {text}  ")
            lbl.setAutoFillBackground(True)
            p = lbl.palette()
            p.setColor(lbl.backgroundRole(), color)
            lbl.setPalette(p)
            leg.addWidget(lbl)
        leg.addStretch()
        root.addLayout(leg)

        # ── Tabs: table + map ─────────────────────────────────────────────
        self._tabs = QTabWidget()

        # Table tab
        self._table = QTableWidget(0, 4)
        self._table.setHorizontalHeaderLabels(
            ["Station", "Period (s)", "Skew β (°)", "Dim"]
        )
        self._table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setAlternatingRowColors(False)
        self._tabs.addTab(self._table, "Classification Table")

        # Bar chart tab
        self._bar_canvas = MplCanvas(parent=self)
        self._tabs.addTab(self._bar_canvas, "Summary Chart")

        # Map tab
        self._map_canvas = MplCanvas(parent=self)
        self._tabs.addTab(self._map_canvas, "Dimensionality Map")

        root.addWidget(self._tabs, stretch=1)

        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        root.addWidget(box)

    # ── Run ───────────────────────────────────────────────────────────────────

    def _on_run(self) -> None:
        self._run_btn.setEnabled(False)
        self._status_lbl.setText("Classifying…")
        self._table.setRowCount(0)

        self._worker = _DimWorker(
            self._sites,
            self._skew_spin.value(),
            self._ellipt_spin.value(),
            self._map_cb.isChecked(),
        )
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, df, map_fig) -> None:
        self._run_btn.setEnabled(True)
        self._df = df
        if df is None or df.empty:
            self._status_lbl.setText("No results.")
            return
        self._fill_table(df)
        self._draw_bar(df)
        if map_fig is not None:
            self._map_canvas.show_figure(map_fig)
        total = len(df)
        counts = df["dim"].value_counts().to_dict() if "dim" in df.columns else {}
        self._status_lbl.setText(
            f"{total} rows — "
            f"1D: {counts.get(0,0)}  2D: {counts.get(1,0)}  3D: {counts.get(2,0)}"
        )

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._status_lbl.setText(f"Error: {msg}")

    # ── Fill table ────────────────────────────────────────────────────────────

    def _fill_table(self, df) -> None:
        self._table.setRowCount(0)
        station_col = next(
            (c for c in ("station", "Station", "id") if c in df.columns), None
        )
        period_col = next(
            (c for c in ("period", "Period", "T") if c in df.columns), None
        )
        skew_col = next(
            (c for c in ("beta_abs", "skew", "beta") if c in df.columns), None
        )
        dim_col = "dim" if "dim" in df.columns else None

        for _, row in df.iterrows():
            r = self._table.rowCount()
            self._table.insertRow(r)
            dim = int(row[dim_col]) if dim_col else 2
            bg = _DIM_COLOR.get(dim, _C_3D)

            def _item(val):
                it = QTableWidgetItem(str(val) if val is not None else "—")
                it.setBackground(bg)
                it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                return it

            self._table.setItem(r, 0, _item(row[station_col] if station_col else "—"))
            self._table.setItem(r, 1, _item(
                f"{row[period_col]:.4g}" if period_col else "—"
            ))
            self._table.setItem(r, 2, _item(
                f"{row[skew_col]:.2f}" if skew_col else "—"
            ))
            self._table.setItem(r, 3, _item(_DIM_LABEL.get(dim, "3D")))

    # ── Bar chart ─────────────────────────────────────────────────────────────

    def _draw_bar(self, df) -> None:
        if "dim" not in df.columns:
            return
        station_col = next(
            (c for c in ("station", "Station", "id") if c in df.columns), None
        )
        if station_col is None:
            return

        stations = df[station_col].unique()
        counts_1d, counts_2d, counts_3d = [], [], []
        for st in stations:
            sub = df[df[station_col] == st]["dim"]
            counts_1d.append((sub == 0).sum())
            counts_2d.append((sub == 1).sum())
            counts_3d.append((sub == 2).sum())

        fig, ax = plt.subplots(figsize=(max(6, len(stations) * 0.5 + 2), 3.5))
        x = np.arange(len(stations))
        w = 0.26
        ax.bar(x - w, counts_1d, width=w, label="1D", color="#66bb6a")
        ax.bar(x,     counts_2d, width=w, label="2D", color="#ffd54f")
        ax.bar(x + w, counts_3d, width=w, label="3D", color="#ef9a9a")
        ax.set_xticks(x)
        ax.set_xticklabels([str(s) for s in stations],
                           rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("# frequency rows")
        ax.set_title("1D / 2D / 3D counts per station")
        ax.legend(fontsize=8)
        fig.tight_layout()
        self._bar_canvas.show_figure(fig)
