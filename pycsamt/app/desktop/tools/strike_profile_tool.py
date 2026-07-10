# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
StrikeProfileDialog — strike angle vs. station-position line plot.

Calls ``pycsamt.emtools.strike.plot_strike_profile`` which renders a line plot
of structural strike angle (°) versus station position with an IQR confidence
ribbon.  Three estimation methods are supported: consensus (default), impedance
sweep, and phase-tensor azimuth.

Usage
-----
    from pycsamt.app.desktop.tools.strike_profile_tool import StrikeProfileDialog
    dlg = StrikeProfileDialog(sites, parent=self)
    dlg.exec()
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

_METHODS = ["consensus", "sweep", "pt"]
_SORT_BY = ["auto", "lon", "lat", "name"]


# ── Worker ────────────────────────────────────────────────────────────────────

class _ProfileWorker(QThread):
    done  = Signal(object)   # Figure
    error = Signal(str)

    def __init__(self, sites, method: str, band, sort_by: str):
        super().__init__()
        self._sites   = sites
        self._method  = method
        self._band    = band
        self._sort_by = sort_by

    def run(self):
        try:
            from pycsamt.emtools.strike import (
                plot_strike_profile,
            )
            set(plt.get_fignums())
            ax = plot_strike_profile(
                self._sites,
                method=self._method,
                band=self._band,
                sort_by=self._sort_by,
                verbose=0,
            )
            fig = ax.figure
            self.done.emit(fig)
        except Exception as exc:
            self.error.emit(str(exc))


# ── Dialog ────────────────────────────────────────────────────────────────────

class StrikeProfileDialog(QDialog):
    """
    Strike angle vs. station-position line plot.

    Parameters
    ----------
    sites : any
        Loaded survey.
    parent : QWidget, optional
    """

    def __init__(self, sites=None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Strike Profile Viewer")
        self.setMinimumSize(860, 480)
        self._sites  = sites
        self._worker = None
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(8)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── Left: controls ────────────────────────────────────────────────
        ctrl = QWidget()
        ctrl.setFixedWidth(220)
        ctrl_lay = QVBoxLayout(ctrl)
        ctrl_lay.setContentsMargins(6, 6, 6, 6)
        ctrl_lay.setSpacing(10)

        grp_method = QGroupBox("Estimation")
        form_m = QFormLayout(grp_method)
        self._method_combo = QComboBox()
        self._method_combo.addItems(_METHODS)
        form_m.addRow("Method:", self._method_combo)
        self._sort_combo = QComboBox()
        self._sort_combo.addItems(_SORT_BY)
        form_m.addRow("Sort by:", self._sort_combo)
        ctrl_lay.addWidget(grp_method)

        grp_band = QGroupBox("Period band (s)  — leave 0 for full range")
        form_b = QFormLayout(grp_band)
        self._tmin_spin = QDoubleSpinBox()
        self._tmin_spin.setRange(0, 1e6)
        self._tmin_spin.setDecimals(4)
        self._tmin_spin.setValue(0.0)
        self._tmin_spin.setSuffix(" s")
        self._tmax_spin = QDoubleSpinBox()
        self._tmax_spin.setRange(0, 1e6)
        self._tmax_spin.setDecimals(2)
        self._tmax_spin.setValue(0.0)
        self._tmax_spin.setSuffix(" s")
        form_b.addRow("T min:", self._tmin_spin)
        form_b.addRow("T max:", self._tmax_spin)
        ctrl_lay.addWidget(grp_band)

        self._run_btn = QPushButton("Plot Profile")
        self._run_btn.clicked.connect(self._on_plot)
        ctrl_lay.addWidget(self._run_btn)

        self._status_lbl = QLabel("")
        self._status_lbl.setWordWrap(True)
        ctrl_lay.addWidget(self._status_lbl)
        ctrl_lay.addStretch()

        if self._sites is None:
            self._run_btn.setEnabled(False)
            self._status_lbl.setText("No survey loaded.")

        splitter.addWidget(ctrl)

        # ── Right: canvas ─────────────────────────────────────────────────
        self._canvas = MplCanvas(parent=self)
        splitter.addWidget(self._canvas)
        splitter.setStretchFactor(1, 1)

        root.addWidget(splitter, stretch=1)

        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        root.addWidget(box)

    # ── Plot ──────────────────────────────────────────────────────────────────

    def _on_plot(self) -> None:
        tmin = self._tmin_spin.value()
        tmax = self._tmax_spin.value()
        band = (tmin, tmax) if (tmin > 0 and tmax > tmin) else None

        self._run_btn.setEnabled(False)
        self._status_lbl.setText("Computing strike profile…")

        self._worker = _ProfileWorker(
            self._sites,
            self._method_combo.currentText(),
            band,
            self._sort_combo.currentText(),
        )
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, fig) -> None:
        self._run_btn.setEnabled(True)
        self._status_lbl.setText("Done.")
        self._canvas.show_figure(fig)

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._status_lbl.setText(f"Error: {msg}")
