# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
StrikeAnalyzerDialog — regional-strike and dimensionality quick analysis.

Backed by
---------
:func:`pycsamt.emtools.strike.estimate_strike_consensus`

Layout
------
  Controls (left)          │  Rose diagram (right)
  ─────────────────────────┼───────────────────────
  Period band  min / max   │  polar histogram of
  Weights  sweep / PT      │  per-station strikes
  [Run]                    │
  ─────────────────────────┴───────────────────────
  Per-station table (station | strike° | IQR | n)
"""
from __future__ import annotations

import numpy as np

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas


# ── Background worker ─────────────────────────────────────────────────────────

class _StrikeWorker(QThread):
    done   = Signal(object)   # DataFrame
    error  = Signal(str)

    def __init__(self, sites, band, w_sweep, w_pt):
        super().__init__()
        self._sites   = sites
        self._band    = band
        self._w_sweep = w_sweep
        self._w_pt    = w_pt

    def run(self):
        try:
            from pycsamt.emtools.strike import estimate_strike_consensus
            df = estimate_strike_consensus(
                self._sites,
                band=self._band,
                w_sweep=self._w_sweep,
                w_pt=self._w_pt,
                verbose=0,
            )
            self.done.emit(df)
        except Exception as exc:
            self.error.emit(str(exc))


# ── Dialog ────────────────────────────────────────────────────────────────────

class StrikeAnalyzerDialog(QDialog):
    """
    Quick regional-strike analysis dialog.

    Parameters
    ----------
    sites : any
        Loaded survey (passed from MainWindow._controller.sites).
    parent : QWidget, optional
    """

    def __init__(self, sites, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Strike Analyzer")
        self.setMinimumSize(860, 560)
        self._sites  = sites
        self._worker = None
        self._df     = None
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        # ── Top split: controls + canvas ──────────────────────────────────────
        top = QSplitter(Qt.Orientation.Horizontal)

        # Controls panel
        ctrl_widget = QWidget()
        ctrl_layout = QVBoxLayout(ctrl_widget)
        ctrl_layout.setSpacing(10)

        grp = QGroupBox("Parameters")
        form = QFormLayout(grp)
        form.setSpacing(8)

        self._band_min = QDoubleSpinBox()
        self._band_min.setRange(0.0001, 100_000.0)
        self._band_min.setDecimals(4)
        self._band_min.setValue(0.001)
        self._band_min.setSuffix(" s")
        form.addRow("Period band  min:", self._band_min)

        self._band_max = QDoubleSpinBox()
        self._band_max.setRange(0.0001, 100_000.0)
        self._band_max.setDecimals(1)
        self._band_max.setValue(1000.0)
        self._band_max.setSuffix(" s")
        form.addRow("Period band  max:", self._band_max)

        self._w_sweep = QDoubleSpinBox()
        self._w_sweep.setRange(0.0, 1.0)
        self._w_sweep.setSingleStep(0.1)
        self._w_sweep.setDecimals(2)
        self._w_sweep.setValue(0.5)
        form.addRow("Weight  sweep:", self._w_sweep)

        self._w_pt = QDoubleSpinBox()
        self._w_pt.setRange(0.0, 1.0)
        self._w_pt.setSingleStep(0.1)
        self._w_pt.setDecimals(2)
        self._w_pt.setValue(0.5)
        form.addRow("Weight  phase-tensor:", self._w_pt)

        ctrl_layout.addWidget(grp)

        self._run_btn = QPushButton("Run Analysis")
        self._run_btn.clicked.connect(self._on_run)
        ctrl_layout.addWidget(self._run_btn)

        self._status_lbl = QLabel("")
        self._status_lbl.setWordWrap(True)
        ctrl_layout.addWidget(self._status_lbl)
        ctrl_layout.addStretch()

        top.addWidget(ctrl_widget)
        top.setStretchFactor(0, 0)

        # Rose diagram canvas
        self._canvas = MplCanvas(self, toolbar=False)
        top.addWidget(self._canvas)
        top.setStretchFactor(1, 1)

        root.addWidget(top, stretch=2)

        # ── Per-station table ──────────────────────────────────────────────────
        self._table = QTableWidget(0, 4)
        self._table.setHorizontalHeaderLabels(
            ["Station", "Strike (°)", "IQR (°)", "n obs"]
        )
        self._table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setAlternatingRowColors(True)
        root.addWidget(self._table, stretch=1)

        # ── Close button ──────────────────────────────────────────────────────
        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        root.addWidget(box)

    # ── Slots ─────────────────────────────────────────────────────────────────

    def _on_run(self) -> None:
        if self._sites is None:
            QMessageBox.warning(self, "No data", "Load survey data first.")
            return
        self._run_btn.setEnabled(False)
        self._status_lbl.setText("Running…")

        band = (self._band_min.value(), self._band_max.value())
        self._worker = _StrikeWorker(
            self._sites, band,
            self._w_sweep.value(), self._w_pt.value(),
        )
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, df) -> None:
        self._df = df
        self._run_btn.setEnabled(True)
        if df.empty:
            self._status_lbl.setText("No strike estimates — check data coverage.")
            return
        self._status_lbl.setText(f"Done — {len(df)} stations.")
        self._fill_table(df)
        self._draw_rose(df)

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._status_lbl.setText(f"Error: {msg}")

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _fill_table(self, df) -> None:
        self._table.setRowCount(0)
        for _, row in df.iterrows():
            r = self._table.rowCount()
            self._table.insertRow(r)
            ang = row.get("ang", float("nan"))
            iqr = row.get("iqr", float("nan"))
            n   = int(row.get("n",  0))
            self._table.setItem(r, 0, QTableWidgetItem(str(row.get("station", ""))))
            self._table.setItem(r, 1, QTableWidgetItem(
                f"{ang:.1f}" if np.isfinite(ang) else "—"
            ))
            self._table.setItem(r, 2, QTableWidgetItem(
                f"{iqr:.1f}" if np.isfinite(iqr) else "—"
            ))
            self._table.setItem(r, 3, QTableWidgetItem(str(n)))

    def _draw_rose(self, df) -> None:
        angles = df["ang"].dropna().values
        if angles.size == 0:
            return
        fig = self._canvas.figure
        fig.clear()
        ax = fig.add_subplot(111, projection="polar")

        # Rose histogram (0-180° because strike is undirected)
        bins = np.linspace(0, np.pi, 19)
        rad  = np.deg2rad(angles % 180)
        counts, edges = np.histogram(rad, bins=bins)
        # Mirror for full rose (0-360)
        width = edges[1] - edges[0]
        ax.bar(edges[:-1], counts, width=width, color="#1e66f5",
               alpha=0.75, align="edge")
        ax.bar(edges[:-1] + np.pi, counts, width=width, color="#1e66f5",
               alpha=0.75, align="edge")

        consensus = float(np.nanmedian(angles % 180))
        ax.set_title(
            f"Regional strike  —  consensus {consensus:.1f}°",
            fontsize=9, pad=12,
        )
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        self._canvas.draw()
