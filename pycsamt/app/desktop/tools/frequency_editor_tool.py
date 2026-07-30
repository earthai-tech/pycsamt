# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
FrequencyEditorDialog — confidence-based frequency QC workflow.

Three tabs:
  1. Settings  — mode (recover / drop / mask), confidence threshold,
                 CI bounds, method, interpolation strategy.
  2. Decision table — colour-coded station × period × action table
                      (green=keep, yellow=recover, red=drop/mask).
  3. Summary chart  — before / after frequency count per station.

After the user is satisfied they can click "Apply to survey" to push
the edited sites back to MainWindow.

Calls:
  ``pycsamt.emtools.frequency.edit_frequencies_by_confidence``
  ``pycsamt.emtools.frequency.plot_frequency_edit_summary``

Usage
-----
    from pycsamt.app.desktop.tools.frequency_editor_tool import FrequencyEditorDialog
    dlg = FrequencyEditorDialog(sites, parent=self)
    if dlg.exec() and dlg.edited_sites is not None:
        # push dlg.edited_sites back to MainWindow
        ...
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

_MODES = ["recover", "drop", "mask"]
_METHODS = ["composite", "snr", "phase_slope", "coherence"]
_INTERPS = ["linear", "nearest"]
_REJECTS = ["drop", "mask", "keep"]
_ALSOS = ["both", "z", "tipper"]

_C_KEEP = QColor("#c8e6c9")
_C_RECOVER = QColor("#fff9c4")
_C_DROP = QColor("#ffcdd2")
_FG_DARK = QColor(
    "#1a1a1a"
)  # forced dark text so pastels stay readable in dark mode
_ACTION_COLOR = {
    "kept": _C_KEEP,
    "recovered": _C_RECOVER,
    "dropped": _C_DROP,
    "masked": _C_DROP,
}


# ── Worker ────────────────────────────────────────────────────────────────────


class _EditWorker(QThread):
    done = Signal(object, object, object)  # (result, decisions_df, fig)
    error = Signal(str)

    def __init__(
        self,
        sites,
        mode: str,
        method: str,
        threshold: float,
        ci_hi: float,
        ci_lo: float,
        interp: str,
        reject: str,
        also: str,
    ):
        super().__init__()
        self._sites = sites
        self._mode = mode
        self._method = method
        self._threshold = threshold
        self._ci_hi = ci_hi
        self._ci_lo = ci_lo
        self._interp = interp
        self._reject = reject
        self._also = also

    def run(self):
        try:
            from pycsamt.emtools.frequency import (
                edit_frequencies_by_confidence,
                plot_frequency_edit_summary,
            )

            result = edit_frequencies_by_confidence(
                self._sites,
                mode=self._mode,
                method=self._method,
                threshold=self._threshold,
                ci_hi=self._ci_hi,
                ci_lo=self._ci_lo,
                interpolation=self._interp,
                reject=self._reject,
                also=self._also,
                inplace=False,
                verbose=0,
            )

            decisions = getattr(result, "decisions", None)
            if hasattr(decisions, "frame"):
                decisions = decisions.frame
            if hasattr(decisions, "_frame"):
                decisions = decisions._frame

            # Summary figure
            before = set(plt.get_fignums())
            try:
                ax = plot_frequency_edit_summary(
                    self._sites,
                    result.sites if hasattr(result, "sites") else self._sites,
                    method=self._method,
                    ci_hi=self._ci_hi,
                    ci_lo=self._ci_lo,
                )
                fig = ax.figure
            except Exception:
                after_nums = set(plt.get_fignums()) - before
                fig = plt.figure(max(after_nums)) if after_nums else None

            self.done.emit(result, decisions, fig)
        except Exception as exc:
            self.error.emit(str(exc))


# ── Dialog ────────────────────────────────────────────────────────────────────


class FrequencyEditorDialog(QDialog):
    """
    Confidence-based frequency QC workflow dialog.

    Attributes
    ----------
    edited_sites : any or None
        Set after a successful run; contains the edited survey.

    Parameters
    ----------
    sites : any
        Loaded survey.
    parent : QWidget, optional
    """

    def __init__(self, sites=None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Frequency Editor")
        self.setMinimumSize(960, 620)
        self._sites = sites
        self._worker = None
        self._result = None
        self.edited_sites = None
        self._build_ui()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(8)

        # ── Settings bar ─────────────────────────────────────────────────
        bar = QHBoxLayout()

        grp_mode = QGroupBox("Edit strategy")
        form_m = QFormLayout(grp_mode)
        self._mode_combo = QComboBox()
        self._mode_combo.addItems(_MODES)
        form_m.addRow("Mode:", self._mode_combo)
        self._method_combo = QComboBox()
        self._method_combo.addItems(_METHODS)
        form_m.addRow("CI method:", self._method_combo)
        bar.addWidget(grp_mode)

        grp_thresh = QGroupBox("Thresholds")
        form_t = QFormLayout(grp_thresh)
        self._thresh_spin = QDoubleSpinBox()
        self._thresh_spin.setRange(0.0, 1.0)
        self._thresh_spin.setDecimals(2)
        self._thresh_spin.setValue(0.50)
        self._thresh_spin.setSingleStep(0.05)
        form_t.addRow("Threshold:", self._thresh_spin)
        self._cihi_spin = QDoubleSpinBox()
        self._cihi_spin.setRange(0.0, 1.0)
        self._cihi_spin.setDecimals(2)
        self._cihi_spin.setValue(0.90)
        self._cihi_spin.setSingleStep(0.05)
        form_t.addRow("CI high:", self._cihi_spin)
        self._cilo_spin = QDoubleSpinBox()
        self._cilo_spin.setRange(0.0, 1.0)
        self._cilo_spin.setDecimals(2)
        self._cilo_spin.setValue(0.50)
        self._cilo_spin.setSingleStep(0.05)
        form_t.addRow("CI low:", self._cilo_spin)
        bar.addWidget(grp_thresh)

        grp_adv = QGroupBox("Advanced")
        form_a = QFormLayout(grp_adv)
        self._interp_combo = QComboBox()
        self._interp_combo.addItems(_INTERPS)
        form_a.addRow("Interpolation:", self._interp_combo)
        self._reject_combo = QComboBox()
        self._reject_combo.addItems(_REJECTS)
        form_a.addRow("Reject rows:", self._reject_combo)
        self._also_combo = QComboBox()
        self._also_combo.addItems(_ALSOS)
        form_a.addRow("Apply to:", self._also_combo)
        bar.addWidget(grp_adv)

        bar.addStretch()

        btn_col = QVBoxLayout()
        self._run_btn = QPushButton("Run")
        self._run_btn.setFixedWidth(90)
        self._run_btn.clicked.connect(self._on_run)
        btn_col.addWidget(self._run_btn)
        self._apply_btn = QPushButton("Apply to survey")
        self._apply_btn.setFixedWidth(120)
        self._apply_btn.setEnabled(False)
        self._apply_btn.clicked.connect(self._on_apply)
        btn_col.addWidget(self._apply_btn)
        bar.addLayout(btn_col)

        root.addLayout(bar)

        self._status_lbl = QLabel("")
        root.addWidget(self._status_lbl)

        if self._sites is None:
            self._run_btn.setEnabled(False)
            self._status_lbl.setText("No survey loaded.")

        # ── Tabs ──────────────────────────────────────────────────────────
        self._tabs = QTabWidget()

        # Decision table
        self._table = QTableWidget(0, 4)
        self._table.setHorizontalHeaderLabels(
            ["Station", "Period (s)", "Confidence", "Action"]
        )
        self._table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._tabs.addTab(self._table, "Decision Table")

        # Summary chart
        self._chart_canvas = MplCanvas(parent=self)
        self._tabs.addTab(self._chart_canvas, "Before / After Chart")

        root.addWidget(self._tabs, stretch=1)

        # Legend
        leg = QHBoxLayout()
        for color, text in (
            (_C_KEEP, "Kept"),
            (_C_RECOVER, "Recovered"),
            (_C_DROP, "Dropped / Masked"),
        ):
            lbl = QLabel(f"  {text}  ")
            lbl.setAutoFillBackground(True)
            p = lbl.palette()
            p.setColor(lbl.backgroundRole(), color)
            p.setColor(lbl.foregroundRole(), _FG_DARK)
            lbl.setPalette(p)
            leg.addWidget(lbl)
        leg.addStretch()
        root.addLayout(leg)

        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        root.addWidget(box)

    # ── Run ───────────────────────────────────────────────────────────────────

    def _on_run(self) -> None:
        self._run_btn.setEnabled(False)
        self._apply_btn.setEnabled(False)
        self._status_lbl.setText("Running frequency editor…")
        self._table.setRowCount(0)

        self._worker = _EditWorker(
            self._sites,
            mode=self._mode_combo.currentText(),
            method=self._method_combo.currentText(),
            threshold=self._thresh_spin.value(),
            ci_hi=self._cihi_spin.value(),
            ci_lo=self._cilo_spin.value(),
            interp=self._interp_combo.currentText(),
            reject=self._reject_combo.currentText(),
            also=self._also_combo.currentText(),
        )
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, result, decisions, fig) -> None:
        self._run_btn.setEnabled(True)
        self._result = result
        self.edited_sites = result.sites if hasattr(result, "sites") else None
        self._apply_btn.setEnabled(self.edited_sites is not None)

        n_drop = getattr(result, "n_dropped", 0)
        n_mask = getattr(result, "n_masked", 0)
        n_recv = getattr(result, "n_recovered", 0)
        self._status_lbl.setText(
            f"Done — dropped: {n_drop}  masked: {n_mask}  recovered: {n_recv}"
        )

        if decisions is not None and not decisions.empty:
            self._fill_table(decisions)
        if fig is not None:
            self._chart_canvas.show_figure(fig)
            self._tabs.setCurrentIndex(1)

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
        conf_col = next(
            (c for c in ("confidence", "ci", "score") if c in df.columns),
            None,
        )
        action_col = "action" if "action" in df.columns else None

        for _, row in df.iterrows():
            r = self._table.rowCount()
            self._table.insertRow(r)
            action = str(row[action_col]) if action_col else "kept"
            bg = _ACTION_COLOR.get(action, _C_KEEP)

            def _item(val):
                it = QTableWidgetItem(str(val) if val is not None else "—")
                it.setBackground(bg)
                it.setForeground(_FG_DARK)
                it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                return it

            self._table.setItem(
                r, 0, _item(row[station_col] if station_col else "—")
            )
            self._table.setItem(
                r, 1, _item(f"{row[period_col]:.4g}" if period_col else "—")
            )
            self._table.setItem(
                r, 2, _item(f"{row[conf_col]:.3f}" if conf_col else "—")
            )
            self._table.setItem(r, 3, _item(action))

    # ── Apply ─────────────────────────────────────────────────────────────────

    def _on_apply(self) -> None:
        self._status_lbl.setText(
            "Edited sites ready — close this dialog to continue."
        )
        self.accept()
