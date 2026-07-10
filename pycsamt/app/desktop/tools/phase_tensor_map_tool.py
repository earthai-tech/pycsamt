# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PhaseTensorMapDialog — geographic map of phase-tensor ellipses at one period.

Calls ``pycsamt.emtools.tensor.plot_phase_tensor_map`` which positions each
station as a tensor ellipse (shape = φ_max/φ_min, orientation = θ, colour =
chosen scalar).  Induction arrows are overlaid when tipper data are available.

Usage
-----
    from pycsamt.app.desktop.tools.phase_tensor_map_tool import PhaseTensorMapDialog
    dlg = PhaseTensorMapDialog(sites, parent=self)
    dlg.exec()
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QCheckBox,
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

_COLOR_BY = ["skew", "ellipt", "theta", "alpha", "s1", "s2"]
_TIPPER_CONV = ["parkinson", "wiese"]
_TIPPER_COMP = ["real", "imag", "amplitude"]


# ── Worker ────────────────────────────────────────────────────────────────────

class _MapWorker(QThread):
    done  = Signal(object)   # Figure
    error = Signal(str)

    def __init__(self, sites, period: float, c_by: str,
                 show_tipper: bool, tipper_conv: str, tipper_comp: str,
                 station_labels: bool):
        super().__init__()
        self._sites         = sites
        self._period        = period
        self._c_by          = c_by
        self._show_tipper   = show_tipper
        self._tipper_conv   = tipper_conv
        self._tipper_comp   = tipper_comp
        self._station_labels = station_labels

    def run(self):
        try:
            from pycsamt.emtools.tensor import (
                plot_phase_tensor_map,
            )
            set(plt.get_fignums())
            ax = plot_phase_tensor_map(
                self._sites,
                period=self._period,
                c_by=self._c_by,
                show_tipper=self._show_tipper,
                tipper_convention=self._tipper_conv,
                tipper_component=self._tipper_comp,
                station_labels=self._station_labels,
                verbose=0,
            )
            fig = ax.figure
            self.done.emit(fig)
        except Exception as exc:
            self.error.emit(str(exc))


# ── Dialog ────────────────────────────────────────────────────────────────────

class PhaseTensorMapDialog(QDialog):
    """
    Geographic phase-tensor ellipse map at a single period.

    Parameters
    ----------
    sites : any
        Loaded survey.
    parent : QWidget, optional
    """

    def __init__(self, sites=None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Phase Tensor Map")
        self.setMinimumSize(920, 640)
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
        ctrl.setFixedWidth(230)
        ctrl_lay = QVBoxLayout(ctrl)
        ctrl_lay.setContentsMargins(6, 6, 6, 6)
        ctrl_lay.setSpacing(10)

        grp_period = QGroupBox("Period")
        form_p = QFormLayout(grp_period)
        self._period_spin = QDoubleSpinBox()
        self._period_spin.setRange(1e-5, 1e5)
        self._period_spin.setDecimals(4)
        self._period_spin.setValue(10.0)
        self._period_spin.setSuffix(" s")
        self._period_spin.setSingleStep(1.0)
        form_p.addRow("T:", self._period_spin)
        ctrl_lay.addWidget(grp_period)

        grp_style = QGroupBox("Ellipse style")
        form_s = QFormLayout(grp_style)
        self._cby_combo = QComboBox()
        self._cby_combo.addItems(_COLOR_BY)
        form_s.addRow("Colour by:", self._cby_combo)
        ctrl_lay.addWidget(grp_style)

        grp_tip = QGroupBox("Tipper arrows")
        form_t = QFormLayout(grp_tip)
        self._show_tipper_cb = QCheckBox("Show tipper")
        self._show_tipper_cb.setChecked(True)
        form_t.addRow(self._show_tipper_cb)
        self._tipper_conv_combo = QComboBox()
        self._tipper_conv_combo.addItems(_TIPPER_CONV)
        form_t.addRow("Convention:", self._tipper_conv_combo)
        self._tipper_comp_combo = QComboBox()
        self._tipper_comp_combo.addItems(_TIPPER_COMP)
        form_t.addRow("Component:", self._tipper_comp_combo)
        ctrl_lay.addWidget(grp_tip)

        grp_map = QGroupBox("Map")
        form_map = QFormLayout(grp_map)
        self._labels_cb = QCheckBox("Station labels")
        self._labels_cb.setChecked(True)
        form_map.addRow(self._labels_cb)
        ctrl_lay.addWidget(grp_map)

        self._run_btn = QPushButton("Draw Map")
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
        self._run_btn.setEnabled(False)
        T = self._period_spin.value()
        self._status_lbl.setText(f"Drawing map at T = {T:.4g} s…")

        self._worker = _MapWorker(
            self._sites,
            period=T,
            c_by=self._cby_combo.currentText(),
            show_tipper=self._show_tipper_cb.isChecked(),
            tipper_conv=self._tipper_conv_combo.currentText(),
            tipper_comp=self._tipper_comp_combo.currentText(),
            station_labels=self._labels_cb.isChecked(),
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
