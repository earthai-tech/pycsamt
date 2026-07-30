# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PhaseTensorStripGridDialog — multi-profile phase-tensor ellipse-strip grid.

Calls ``pycsamt.emtools.tensor.plot_phase_tensor_strip_grid``, which tiles one
single-station ellipse strip (:func:`~pycsamt.emtools.tensor.plot_phase_tensor_strip`)
per selected station (rows) across survey lines/profiles (columns), sharing a
single colour scale and one colorbar for the whole figure.

Survey lines are auto-detected from station IDs using the same detector the
web app's Lines panel uses (:func:`pycsamt.site.lines.detect_lines_from_station_ids`),
so this works out of the box for any loaded survey without extra metadata.

Usage
-----
    from pycsamt.app.desktop.tools.phase_tensor_strip_grid_tool import (
        PhaseTensorStripGridDialog,
    )
    dlg = PhaseTensorStripGridDialog(sites, parent=self)
    dlg.exec()
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QPushButton,
    QSpinBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.site.lines import pick_representative_stations

_COLOR_BY = ["skew", "ellipt", "theta", "alpha", "phi_min", "phi_max"]


# ── Worker ────────────────────────────────────────────────────────────────────


class _StripGridWorker(QThread):
    done = Signal(object)  # Figure
    error = Signal(str)

    def __init__(self, sites, profiles: dict, c_by: str):
        super().__init__()
        self._sites = sites
        self._profiles = profiles
        self._c_by = c_by

    def run(self):
        try:
            from pycsamt.emtools.tensor import (
                plot_phase_tensor_strip_grid,
            )

            fig = plot_phase_tensor_strip_grid(
                self._sites,
                profiles=self._profiles,
                c_by=self._c_by,
                verbose=0,
            )
            self.done.emit(fig)
        except Exception as exc:
            self.error.emit(str(exc))


# ── Dialog ────────────────────────────────────────────────────────────────────


class PhaseTensorStripGridDialog(QDialog):
    """
    Multi-profile phase-tensor ellipse-strip grid.

    Parameters
    ----------
    sites : any
        Loaded survey.
    parent : QWidget, optional
    """

    def __init__(self, sites=None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Phase Tensor Strip Grid")
        self.setMinimumSize(1000, 640)
        self._sites = sites
        self._worker = None
        self._lines: dict[str, list] = {}
        self._build_ui()
        self._detect_lines()

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

        grp_lines = QGroupBox("Survey lines (auto-detected)")
        lay_lines = QVBoxLayout(grp_lines)
        self._lines_lbl = QLabel("—")
        self._lines_lbl.setWordWrap(True)
        lay_lines.addWidget(self._lines_lbl)
        ctrl_lay.addWidget(grp_lines)

        grp_sel = QGroupBox("Selection")
        form_sel = QFormLayout(grp_sel)
        self._per_line_spin = QSpinBox()
        self._per_line_spin.setRange(1, 20)
        self._per_line_spin.setValue(4)
        form_sel.addRow("Stations / line:", self._per_line_spin)
        ctrl_lay.addWidget(grp_sel)

        grp_style = QGroupBox("Ellipse style")
        form_s = QFormLayout(grp_style)
        self._cby_combo = QComboBox()
        self._cby_combo.addItems(_COLOR_BY)
        form_s.addRow("Colour by:", self._cby_combo)
        ctrl_lay.addWidget(grp_style)

        self._run_btn = QPushButton("Draw Grid")
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

    # ── Line detection ────────────────────────────────────────────────────────

    def _detect_lines(self) -> None:
        if self._sites is None:
            return
        try:
            from pycsamt.emtools._core import (
                _iter_items,
                _unwrap,
            )
            from pycsamt.site.lines import (
                detect_lines_from_station_ids,
            )

            names = []
            for ed in _iter_items(self._sites):
                try:
                    ed = _unwrap(ed)
                except Exception:
                    pass
                nm = (
                    getattr(ed, "station", None)
                    or getattr(ed, "id", None)
                    or "?"
                )
                names.append(str(nm))
            self._lines = detect_lines_from_station_ids(names)
            summary = ", ".join(
                f"{k} ({len(v)})" for k, v in self._lines.items()
            )
            self._lines_lbl.setText(summary or "No stations found.")
        except Exception as exc:
            self._lines_lbl.setText(f"Line detection failed: {exc}")
            self._lines = {}

    # ── Plot ──────────────────────────────────────────────────────────────────

    def _on_plot(self) -> None:
        if not self._lines:
            self._status_lbl.setText("No survey lines detected.")
            return
        k = self._per_line_spin.value()
        profiles = {
            ln: pick_representative_stations(names, k)
            for ln, names in self._lines.items()
        }

        self._run_btn.setEnabled(False)
        self._status_lbl.setText("Drawing ellipse-strip grid…")

        self._worker = _StripGridWorker(
            self._sites,
            profiles=profiles,
            c_by=self._cby_combo.currentText(),
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
