# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PhaseTensorMapDialog — geographic map of phase-tensor ellipses.

Two modes:

  Single period  Calls ``pycsamt.emtools.tensor.plot_phase_tensor_map`` --
                 one panel, each station drawn as a tensor ellipse (shape =
                 φ_max/φ_min, orientation = θ, colour = chosen scalar).
                 Induction arrows are overlaid when tipper data are
                 available.
  Multi-period   Calls ``pycsamt.emtools.tensor.plot_phase_tensor_map_grid``
  grid           -- the same ellipse map tiled across several periods, one
                 panel each, sharing a single colour scale/colorbar so the
                 panels are directly comparable ("skew at 30/3/0.3/0.03 Hz").

Both auto-select the ``pt_skew`` / ``pt_skew_abs`` diverging/sequential
colormaps from :mod:`pycsamt.api.style` whenever ``c_by`` is skew-like and
no explicit ``cmap`` is passed -- neither mode needs to choose a colormap.

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
    QLineEdit,
    QPushButton,
    QSpinBox,
    QSplitter,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView

_COLOR_BY = ["skew", "ellipt", "theta", "alpha", "s1", "s2"]
_TIPPER_CONV = ["parkinson", "wiese"]
_TIPPER_COMP = ["real", "imag", "amplitude"]
_MODES = ["Single period", "Multi-period grid"]


def _parse_periods(text: str) -> list[float]:
    """Parse a comma/whitespace-separated period list, e.g. ``"1, 10 100"``."""
    raw = [tok for tok in text.replace(",", " ").split() if tok]
    values = [float(tok) for tok in raw]
    if not values:
        raise ValueError("Enter at least one period (seconds).")
    return values


# ── Worker ────────────────────────────────────────────────────────────────────


class _MapWorker(QThread):
    done = Signal(object)  # Figure
    error = Signal(str)

    def __init__(
        self,
        sites,
        period: float,
        c_by: str,
        show_tipper: bool,
        tipper_conv: str,
        tipper_comp: str,
        station_labels: bool,
        *,
        grid: bool = False,
        periods: list[float] | None = None,
        n_cols: int | None = None,
        panel_labels: bool = True,
        abs_skew: bool = False,
    ):
        super().__init__()
        self._sites = sites
        self._period = period
        self._c_by = c_by
        self._show_tipper = show_tipper
        self._tipper_conv = tipper_conv
        self._tipper_comp = tipper_comp
        self._station_labels = station_labels
        self._grid = grid
        self._periods = periods
        self._n_cols = n_cols
        self._panel_labels = panel_labels
        self._abs_skew = abs_skew

    def run(self):
        try:
            if self._grid:
                from pycsamt.emtools.tensor import (
                    plot_phase_tensor_map_grid,
                )

                fig = plot_phase_tensor_map_grid(
                    self._sites,
                    periods=self._periods,
                    n_cols=self._n_cols,
                    panel_labels=self._panel_labels,
                    c_by=self._c_by,
                    abs_skew=self._abs_skew,
                    station_labels=self._station_labels,
                    show_tipper=self._show_tipper,
                    tipper_convention=self._tipper_conv,
                    tipper_component=self._tipper_comp,
                    verbose=0,
                )
            else:
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
    Geographic phase-tensor ellipse map — one period, or a multi-period grid.

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
        self._sites = sites
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

        grp_mode = QGroupBox("Mode")
        form_mode = QFormLayout(grp_mode)
        self._mode_combo = QComboBox()
        self._mode_combo.addItems(_MODES)
        self._mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        form_mode.addRow(self._mode_combo)
        ctrl_lay.addWidget(grp_mode)

        # Single-period vs. multi-period grid inputs share one stack so
        # switching modes never leaves a stale/hidden control's value applied.
        self._period_stack = QStackedWidget()

        grp_period = QGroupBox("Period")
        form_p = QFormLayout(grp_period)
        self._period_spin = QDoubleSpinBox()
        self._period_spin.setRange(1e-5, 1e5)
        self._period_spin.setDecimals(4)
        self._period_spin.setValue(10.0)
        self._period_spin.setSuffix(" s")
        self._period_spin.setSingleStep(1.0)
        form_p.addRow("T:", self._period_spin)
        self._period_stack.addWidget(grp_period)

        grp_grid = QGroupBox("Periods")
        form_g = QFormLayout(grp_grid)
        self._periods_edit = QLineEdit("1, 10, 100")
        self._periods_edit.setToolTip(
            "Comma/space-separated periods in seconds, one panel each"
            " -- e.g. the classic 30/3/0.3/0.03 Hz set is periods"
            " 0.033, 0.33, 3.3, 33."
        )
        form_g.addRow("T list (s):", self._periods_edit)
        self._ncols_spin = QSpinBox()
        self._ncols_spin.setRange(0, 6)
        self._ncols_spin.setValue(2)
        self._ncols_spin.setSpecialValueText("Auto")
        form_g.addRow("Columns:", self._ncols_spin)
        self._panel_labels_cb = QCheckBox("Panel labels (a, b, …)")
        self._panel_labels_cb.setChecked(True)
        form_g.addRow(self._panel_labels_cb)
        self._abs_skew_cb = QCheckBox("Use |skew| (sequential colour)")
        self._abs_skew_cb.setToolTip(
            "Colour by |β| with the pt_skew_abs colormap instead of signed"
            " β with the diverging pt_skew colormap. Only applies when"
            " 'Colour by' is skew."
        )
        form_g.addRow(self._abs_skew_cb)
        self._period_stack.addWidget(grp_grid)

        ctrl_lay.addWidget(self._period_stack)

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
        self._canvas_view = CanvasResultView(
            self,
            toolbar=True,
            empty_title="No map yet",
            empty_reason="Click Draw Map to render the phase-tensor map.",
        )
        self._canvas = self._canvas_view.canvas
        self._canvas.set_refresh_callback(
            self._on_plot, tooltip="Redraw the phase-tensor map"
        )
        splitter.addWidget(self._canvas_view)
        splitter.setStretchFactor(1, 1)

        root.addWidget(splitter, stretch=1)

        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        root.addWidget(box)

    # ── Slots ─────────────────────────────────────────────────────────────────

    def _on_mode_changed(self, index: int) -> None:
        self._period_stack.setCurrentIndex(index)

    # ── Plot ──────────────────────────────────────────────────────────────────

    def _on_plot(self) -> None:
        is_grid = self._mode_combo.currentIndex() == 1

        if is_grid:
            try:
                periods = _parse_periods(self._periods_edit.text())
            except ValueError as exc:
                self._status_lbl.setText(f"Error: {exc}")
                return
            self._status_lbl.setText(
                f"Drawing grid over {len(periods)} period(s)…"
            )
            n_cols = self._ncols_spin.value()
            worker_kwargs = dict(
                grid=True,
                periods=periods,
                n_cols=(n_cols if n_cols > 0 else None),
                panel_labels=self._panel_labels_cb.isChecked(),
                abs_skew=self._abs_skew_cb.isChecked(),
            )
            T = periods[0]  # unused by the grid path; kept for a valid call
        else:
            T = self._period_spin.value()
            self._status_lbl.setText(f"Drawing map at T = {T:.4g} s…")
            worker_kwargs = {}

        self._run_btn.setEnabled(False)
        self._worker = _MapWorker(
            self._sites,
            period=T,
            c_by=self._cby_combo.currentText(),
            show_tipper=self._show_tipper_cb.isChecked(),
            tipper_conv=self._tipper_conv_combo.currentText(),
            tipper_comp=self._tipper_comp_combo.currentText(),
            station_labels=self._labels_cb.isChecked(),
            **worker_kwargs,
        )
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, fig) -> None:
        self._run_btn.setEnabled(True)
        self._status_lbl.setText("Done.")
        self._canvas.show_figure(fig)
        self._canvas_view.show_canvas()

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._status_lbl.setText(f"Error: {msg}")
        self._canvas_view.show_unavailable("Phase tensor map unavailable", msg)
