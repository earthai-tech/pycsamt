# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PipelineWindow — 3-panel processing pipeline window.

┌──────────────┬──────────────────────────┬───────────────────┐
│  STEPPER     │  STEP CONFIGURATION      │  LOG + PREVIEW    │
│  (left)      │  (centre)                │  (right)          │
│              │                          │                   │
│ ① Load   ✅  │  Method  [AMA ▼]         │  timestamped log  │
│ ② QC     ✅  │  Window  [3  ]           │                   │
│ ③ SS     🔄  │  Weights [tri▼]          │  small canvas     │
│ ④ …      ⬜  │                          │  (step preview)   │
│              │  [▶ Run Step] [⏭ Skip]  │                   │
│ [▶ Quick]   │  ─────────────────────── │  [📊 Full view]  │
│ [▶ Run All] │  Result:                 │                   │
│ [⏹ Stop]   │  28/28 corrected          │                   │
│ [↺ Reset]  │  0.3 s                   │                   │
│ [💾 Save]  │                          │                   │
│             │                          │                   │
│ [████] 37%  │                          │                   │
└──────────────┴──────────────────────────┴───────────────────┘
"""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QByteArray, Qt, Signal, Slot
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTabWidget,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.pipeline_controller import (
    STATUS_ICON,
    ParamSpec,
    PipelineController,
    PipelineStep,
    StepStatus,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

_ICONS = Path(__file__).parent.parent / "resources" / "icons"


def _icon(name: str) -> QIcon:
    for c in (name, f"{name}.svg", f"{name}.png"):
        p = _ICONS / c
        if p.exists():
            return QIcon(str(p))
    return QIcon()


def _btn(
    text: str, icon_name: str = "", tip: str = "", fixed_w: int = 0
) -> QPushButton:
    b = QPushButton(text)
    if icon_name:
        ic = _icon(icon_name)
        if not ic.isNull():
            b.setIcon(ic)
    if tip:
        b.setToolTip(tip)
    if fixed_w:
        b.setFixedWidth(fixed_w)
    return b


# ── PipelineWindow ────────────────────────────────────────────────────────────


class PipelineWindow(QWidget):
    """
    Standalone 3-panel pipeline window.

    Signals
    -------
    pipeline_finished(object)
        Emitted when all steps complete, carrying the output Sites.
    """

    pipeline_finished = Signal(object)

    def __init__(self, parent: QWidget | None = None) -> None:
        flags = (
            Qt.WindowType.Window
            | Qt.WindowType.WindowCloseButtonHint
            | Qt.WindowType.WindowMinimizeButtonHint
            | Qt.WindowType.WindowMaximizeButtonHint
        )
        super().__init__(parent, flags)
        self.setWindowTitle("pycsamt — Processing Pipeline")
        ic = _icon("pipeline")
        if not ic.isNull():
            self.setWindowIcon(ic)
        self.resize(1300, 800)

        self._ctrl = PipelineController()
        self._worker = None
        self._selected_step: int = 0
        self._param_widgets: dict = {}

        self._build_ui()
        self._refresh_stepper()
        self._select_step(0)

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(4)

        splitter.addWidget(self._build_stepper_panel())
        splitter.addWidget(self._build_config_panel())
        splitter.addWidget(self._build_output_panel())

        splitter.setSizes([200, 420, 300])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 1)

        root.addWidget(splitter)

    # ── LEFT PANEL: Stepper ───────────────────────────────────────────────────

    def _build_stepper_panel(self) -> QWidget:
        w = QWidget()
        w.setObjectName("StepperPanel")
        w.setMinimumWidth(170)
        w.setMaximumWidth(240)
        v = QVBoxLayout(w)
        v.setContentsMargins(6, 8, 6, 8)
        v.setSpacing(6)

        title = QLabel("PIPELINE")
        title.setObjectName("PanelTitle")
        title.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        v.addWidget(title)

        # Step list
        self._step_list = QListWidget()
        self._step_list.setObjectName("StepperList")
        self._step_list.currentRowChanged.connect(self._select_step)
        v.addWidget(self._step_list, stretch=1)

        # Overall progress
        sep = QLabel("─" * 22)
        sep.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        sep.setObjectName("InfoLabel")
        v.addWidget(sep)

        self._overall_progress = QProgressBar()
        self._overall_progress.setRange(0, 100)
        self._overall_progress.setValue(0)
        self._overall_progress.setMaximumHeight(14)
        v.addWidget(self._overall_progress)

        self._progress_lbl = QLabel("0 / 8 steps done")
        self._progress_lbl.setObjectName("InfoLabel")
        self._progress_lbl.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        v.addWidget(self._progress_lbl)

        sep2 = QLabel("─" * 22)
        sep2.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        sep2.setObjectName("InfoLabel")
        v.addWidget(sep2)

        # Control buttons
        self._btn_quick = _btn(
            "⚡ Quick Pipeline",
            "pipeline",
            "Run all 8 steps with default parameters",
        )
        self._btn_run = _btn(
            "▶  Run All", "processing", "Run all pending steps sequentially"
        )
        self._btn_stop = _btn(
            "⏹  Stop", "", "Interrupt after the current step"
        )
        self._btn_reset = _btn("↺  Reset", "", "Reset all steps to Pending")
        self._btn_save = _btn(
            "💾  Save config", "save", "Save pipeline configuration to JSON"
        )
        self._btn_load = _btn(
            "📂  Load config", "open", "Load pipeline configuration from JSON"
        )

        self._btn_stop.setEnabled(False)
        for b in (
            self._btn_quick,
            self._btn_run,
            self._btn_stop,
            self._btn_reset,
            self._btn_save,
            self._btn_load,
        ):
            v.addWidget(b)

        self._btn_quick.clicked.connect(self._on_quick_pipeline)
        self._btn_run.clicked.connect(self._on_run_all)
        self._btn_stop.clicked.connect(self._on_stop)
        self._btn_reset.clicked.connect(self._on_reset)
        self._btn_save.clicked.connect(self._on_save_config)
        self._btn_load.clicked.connect(self._on_load_config)

        return w

    # ── CENTRE PANEL: Step configuration ─────────────────────────────────────

    def _build_config_panel(self) -> QWidget:
        w = QWidget()
        w.setObjectName("ConfigPanel")
        v = QVBoxLayout(w)
        v.setContentsMargins(10, 10, 10, 10)
        v.setSpacing(8)

        # Step title
        self._step_title = QLabel("")
        self._step_title.setObjectName("StepTitle")
        self._step_title.setWordWrap(True)
        v.addWidget(self._step_title)

        # Step description
        self._step_desc = QLabel("")
        self._step_desc.setObjectName("InfoLabel")
        self._step_desc.setWordWrap(True)
        v.addWidget(self._step_desc)

        # Method selector
        grp_method = QGroupBox("Method")
        grp_method.setObjectName("ParamsGroup")
        m_lay = QVBoxLayout(grp_method)
        self._combo_method = QComboBox()
        self._combo_method.currentTextChanged.connect(self._on_method_changed)
        m_lay.addWidget(self._combo_method)
        v.addWidget(grp_method)

        # Dynamic parameters
        self._grp_params = QGroupBox("Parameters")
        self._grp_params.setObjectName("ParamsGroup")
        self._form_params = QFormLayout(self._grp_params)
        self._form_params.setSpacing(5)
        v.addWidget(self._grp_params)

        # Step result info (after run)
        self._grp_result = QGroupBox("Last result")
        self._grp_result.setObjectName("ParamsGroup")
        r_lay = QVBoxLayout(self._grp_result)
        self._result_lbl = QLabel("—")
        self._result_lbl.setWordWrap(True)
        self._result_lbl.setObjectName("InfoLabel")
        self._elapsed_lbl = QLabel("")
        self._elapsed_lbl.setObjectName("InfoLabel")
        r_lay.addWidget(self._result_lbl)
        r_lay.addWidget(self._elapsed_lbl)
        self._grp_result.setVisible(False)
        v.addWidget(self._grp_result)

        v.addStretch(1)

        # Run / Skip buttons
        btn_row = QHBoxLayout()
        self._btn_run_step = _btn(
            "▶  Run step", "processing", "Execute only this step"
        )
        self._btn_skip = _btn("⏭  Skip", "", "Mark this step as skipped")
        btn_row.addWidget(self._btn_run_step)
        btn_row.addWidget(self._btn_skip)
        v.addLayout(btn_row)

        self._btn_run_step.clicked.connect(self._on_run_step)
        self._btn_skip.clicked.connect(self._on_skip_step)

        # Step progress bar (per-step)
        self._step_progress = QProgressBar()
        self._step_progress.setRange(0, 0)  # indeterminate
        self._step_progress.setMaximumHeight(10)
        self._step_progress.setTextVisible(False)
        self._step_progress.setVisible(False)
        v.addWidget(self._step_progress)

        return w

    # ── RIGHT PANEL: Log + preview ────────────────────────────────────────────

    def _build_output_panel(self) -> QWidget:
        w = QWidget()
        w.setObjectName("OutputPanel")
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)

        tabs = QTabWidget()
        tabs.setDocumentMode(True)

        # Log tab
        self._log_text = QPlainTextEdit()
        self._log_text.setReadOnly(True)
        self._log_text.setMaximumBlockCount(5000)
        self._log_text.setPlaceholderText("Pipeline log will appear here…")
        tabs.addTab(self._log_text, "Log")

        # Preview tab (small canvas + full-view button)
        preview_w = QWidget()
        pv = QVBoxLayout(preview_w)
        pv.setContentsMargins(4, 4, 4, 4)
        self._preview_canvas = MplCanvas(preview_w, toolbar=False)
        self._preview_canvas.setMinimumHeight(200)
        pv.addWidget(self._preview_canvas)
        self._btn_full_view = _btn(
            "📊  View full plot",
            "profile-view",
            "Open this plot in a larger window",
        )
        self._btn_full_view.clicked.connect(self._on_full_view)
        pv.addWidget(self._btn_full_view)
        tabs.addTab(preview_w, "Preview")

        # Summary tab (HTML)
        self._summary_browser = QTextBrowser()
        self._summary_browser.setPlaceholderText(
            "Pipeline summary will appear here…"
        )
        tabs.addTab(self._summary_browser, "Summary")

        v.addWidget(tabs)
        return w

    # ── Public API ────────────────────────────────────────────────────────────

    def set_input_sites(self, sites) -> None:
        self._ctrl.set_input_sites(sites)
        step0 = self._ctrl.steps[0]
        if step0.active_method == "current" and sites is not None:
            from pycsamt.app.desktop.controllers.pipeline_controller import (
                StepStatus,
            )

            step0.status = StepStatus.DONE
            step0.output_sites = sites
            step0.result_info = f"{self._ctrl._n(sites)} stations ready"
            self._ctrl._sites_chain[0] = sites
        self._refresh_stepper()
        if self._selected_step == 0:
            self._select_step(0)

    def set_dark_mode(self, dark: bool) -> None:
        pass  # theme applied via QSS

    # ── Stepper refresh ───────────────────────────────────────────────────────

    def _refresh_stepper(self) -> None:
        self._step_list.blockSignals(True)
        cur = self._step_list.currentRow()
        self._step_list.clear()
        done = 0
        for step in self._ctrl.steps:
            icon = STATUS_ICON[step.status]
            item = QListWidgetItem(f"{icon}  {step.id + 1}. {step.name}")
            item.setTextAlignment(Qt.AlignmentFlag.AlignVCenter)
            # Colour by status
            if step.status == StepStatus.DONE:
                item.setForeground(Qt.GlobalColor.green)
                done += 1
            elif step.status == StepStatus.ERROR:
                item.setForeground(Qt.GlobalColor.red)
            elif step.status == StepStatus.RUNNING:
                item.setForeground(Qt.GlobalColor.yellow)
            elif step.status == StepStatus.SKIPPED:
                item.setForeground(Qt.GlobalColor.gray)
            self._step_list.addItem(item)

        self._step_list.blockSignals(False)
        if cur >= 0:
            self._step_list.setCurrentRow(cur)

        pct = int(done / len(self._ctrl.steps) * 100)
        self._overall_progress.setValue(pct)
        self._progress_lbl.setText(
            f"{done} / {len(self._ctrl.steps)} steps done"
        )

    # ── Step selection ────────────────────────────────────────────────────────

    def _select_step(self, row: int) -> None:
        if row < 0 or row >= len(self._ctrl.steps):
            return
        self._selected_step = row
        step = self._ctrl.steps[row]

        self._step_title.setText(
            f"{STATUS_ICON[step.status]}  Step {row + 1}:  {step.name}"
        )
        self._step_desc.setText(step.description)

        # Populate method combo
        self._combo_method.blockSignals(True)
        self._combo_method.clear()
        for m in step.methods:
            self._combo_method.addItem(m.label, userData=m.name)
        # Select the active method
        for i in range(self._combo_method.count()):
            if self._combo_method.itemData(i) == step.active_method:
                self._combo_method.setCurrentIndex(i)
                break
        self._combo_method.blockSignals(False)
        self._rebuild_params_form(step)

        # Result info
        if step.status in (StepStatus.DONE, StepStatus.ERROR):
            self._grp_result.setVisible(True)
            self._result_lbl.setText(step.result_info or step.error_msg or "—")
            if step.elapsed_s > 0:
                self._elapsed_lbl.setText(f"Time: {step.elapsed_s:.2f}s")
            # Show preview if available
            self._draw_step_preview(step)
        else:
            self._grp_result.setVisible(False)

        # Enable/disable Run step based on state
        can_run = (
            self._worker is None or not self._worker.isRunning()
        ) and step.status != StepStatus.RUNNING
        self._btn_run_step.setEnabled(can_run)
        self._btn_skip.setEnabled(can_run and step.can_skip)

    # ── Params form ───────────────────────────────────────────────────────────

    def _on_method_changed(self, label: str) -> None:
        if self._selected_step < 0:
            return
        step = self._ctrl.steps[self._selected_step]
        # Find method by label
        for m in step.methods:
            if m.label == label:
                step.active_method = m.name
                step.reset_params()
                break
        self._rebuild_params_form(step)

    def _rebuild_params_form(self, step: PipelineStep) -> None:
        while self._form_params.rowCount():
            self._form_params.removeRow(0)
        self._param_widgets.clear()

        spec = step.active_method_spec
        if spec is None or not spec.params:
            self._grp_params.setVisible(False)
            return

        self._grp_params.setVisible(True)
        for param_name, ps in spec.params.items():
            w = self._make_param_widget(
                ps, step.params.get(param_name, ps.default)
            )
            # Wire change → save to step.params
            self._connect_param_widget(w, ps, param_name, step)
            self._form_params.addRow(f"{ps.label}:", w)
            self._param_widgets[param_name] = w

    def _make_param_widget(self, ps: ParamSpec, value) -> QWidget:
        if ps.kind == "combo":
            w = QComboBox()
            for opt in ps.options:
                w.addItem(str(opt))
            idx = w.findText(str(value))
            if idx >= 0:
                w.setCurrentIndex(idx)
            return w
        if ps.kind == "float":
            w = QDoubleSpinBox()
            w.setRange(ps.min_val, ps.max_val)
            w.setDecimals(ps.decimals)
            w.setValue(float(value) if value is not None else ps.default)
            return w
        if ps.kind == "int":
            w = QSpinBox()
            w.setRange(int(ps.min_val), int(ps.max_val))
            w.setValue(int(value) if value is not None else int(ps.default))
            return w
        if ps.kind == "bool":
            w = QCheckBox()
            w.setChecked(bool(value))
            return w
        if ps.kind == "path":
            container = QWidget()
            row = QHBoxLayout(container)
            row.setContentsMargins(0, 0, 0, 0)
            edit = QLineEdit(str(value or ""))
            edit.setObjectName("path_edit")
            browse = QPushButton("…")
            browse.setFixedWidth(26)
            browse.clicked.connect(lambda: self._browse_path(edit))
            row.addWidget(edit)
            row.addWidget(browse)
            return container
        w = QLineEdit(str(value or ""))
        return w

    def _connect_param_widget(
        self, w, ps: ParamSpec, param_name: str, step: PipelineStep
    ) -> None:
        def _save(val):
            step.params[param_name] = val

        if ps.kind == "combo":
            w.currentTextChanged.connect(_save)
        elif ps.kind in ("float",):
            w.valueChanged.connect(_save)
        elif ps.kind in ("int",):
            w.valueChanged.connect(_save)
        elif ps.kind == "bool":
            w.toggled.connect(_save)
        elif ps.kind == "path":
            edit = w.findChild(QLineEdit)
            if edit:
                edit.textChanged.connect(_save)
        else:
            if hasattr(w, "textChanged"):
                w.textChanged.connect(_save)

    @staticmethod
    def _browse_path(edit: QLineEdit) -> None:
        path = QFileDialog.getExistingDirectory(None, "Select folder")
        if path:
            edit.setText(path)

    # ── Pipeline execution ────────────────────────────────────────────────────

    def _on_quick_pipeline(self) -> None:
        """Reset all steps to defaults then run all."""
        self._ctrl.reset()
        # Keep step 0 pre-loaded if we have sites
        if self._ctrl._sites_input is not None:
            step0 = self._ctrl.steps[0]
            step0.active_method = "current"
            step0.status = StepStatus.DONE
            step0.output_sites = self._ctrl._sites_input
            step0.result_info = (
                f"{self._ctrl._n(self._ctrl._sites_input)} stations"
            )
            self._ctrl._sites_chain[0] = self._ctrl._sites_input
        self._refresh_stepper()
        self._log(
            "⚡ Quick Pipeline — running all steps with default parameters."
        )
        self._start_worker(list(range(len(self._ctrl.steps))))

    def _on_run_all(self) -> None:
        """Run all steps that are not yet DONE or SKIPPED."""
        pending = [
            s.id
            for s in self._ctrl.steps
            if s.status not in (StepStatus.DONE, StepStatus.SKIPPED)
        ]
        if not pending:
            self._log("All steps already complete.")
            return
        self._start_worker(pending)

    def _on_run_step(self) -> None:
        """Run only the currently selected step."""
        self._start_worker([self._selected_step])

    def _on_skip_step(self) -> None:
        step = self._ctrl.steps[self._selected_step]
        step.status = StepStatus.SKIPPED
        step.result_info = "Skipped by user."
        self._refresh_stepper()
        self._select_step(self._selected_step)
        self._log(
            f"Step {self._selected_step + 1} ({step.name}) marked as skipped."
        )

    def _on_stop(self) -> None:
        if self._worker and self._worker.isRunning():
            self._worker.requestInterruption()
            self._log("Stop requested — waiting for current step to finish…")

    def _on_reset(self) -> None:
        if self._worker and self._worker.isRunning():
            return
        self._ctrl.reset()
        # If sites are available, keep step 0 pre-done
        if self._ctrl._sites_input is not None:
            self.set_input_sites(self._ctrl._sites_input)
        else:
            self._refresh_stepper()
        self._log("Pipeline reset.")
        self._select_step(0)

    def _start_worker(self, step_ids: list[int]) -> None:
        if self._worker and self._worker.isRunning():
            return
        if not self._ctrl.is_ready():
            self._log("No survey data loaded — run step 0 (Load) first.")
            return

        from pycsamt.app.desktop.workers.pipeline_worker import (
            PipelineWorker,
        )

        self._worker = PipelineWorker(self._ctrl, step_ids, parent=self)
        self._worker.step_started.connect(self._on_step_started)
        self._worker.step_done.connect(self._on_step_done)
        self._worker.step_error.connect(self._on_step_error)
        self._worker.step_skipped.connect(self._on_step_skipped)
        self._worker.log_line.connect(self._on_log_line)
        self._worker.progress.connect(self._overall_progress.setValue)
        self._worker.all_done.connect(self._on_all_done)

        self._btn_run.setEnabled(False)
        self._btn_quick.setEnabled(False)
        self._btn_run_step.setEnabled(False)
        self._btn_stop.setEnabled(True)
        self._worker.start()

    # ── Worker slots ──────────────────────────────────────────────────────────

    @Slot(int)
    def _on_step_started(self, idx: int) -> None:
        self._ctrl.steps[idx].status = StepStatus.RUNNING
        self._refresh_stepper()
        self._step_progress.setVisible(True)
        if idx == self._selected_step:
            self._select_step(idx)

    @Slot(int, object)
    def _on_step_done(self, idx: int, result) -> None:
        self._refresh_stepper()
        self._step_progress.setVisible(False)
        if idx == self._selected_step:
            self._select_step(idx)
        # Jump stepper to next step
        next_idx = idx + 1
        if next_idx < self._step_list.count():
            self._step_list.setCurrentRow(next_idx)

    @Slot(int, str)
    def _on_step_error(self, idx: int, msg: str) -> None:
        self._refresh_stepper()
        self._step_progress.setVisible(False)
        if idx == self._selected_step:
            self._select_step(idx)

    @Slot(int)
    def _on_step_skipped(self, idx: int) -> None:
        self._refresh_stepper()

    @Slot(str)
    def _on_log_line(self, line: str) -> None:
        self._log_text.appendPlainText(line)

    @Slot()
    def _on_all_done(self) -> None:
        self._btn_run.setEnabled(True)
        self._btn_quick.setEnabled(True)
        self._btn_run_step.setEnabled(True)
        self._btn_stop.setEnabled(False)
        self._step_progress.setVisible(False)
        self._refresh_stepper()
        self._update_summary()

        result = self._ctrl.get_output_sites()
        if result is not None:
            self.pipeline_finished.emit(result)

    # ── Preview ───────────────────────────────────────────────────────────────

    def _draw_step_preview(self, step: PipelineStep) -> None:
        if step.diag_fn is None or step.output_sites is None:
            return
        import inspect

        import pycsamt.emtools as et

        fn = getattr(et, step.diag_fn, None)
        if fn is None:
            return
        fig = self._preview_canvas.figure
        fig.clear()
        try:
            sig = inspect.signature(fn)
            ax = fig.add_subplot(111)
            if "ax" in sig.parameters:
                fn(step.output_sites, ax=ax, verbose=0)
            else:
                fn(step.output_sites, verbose=0)
        except Exception:
            pass
        self._preview_canvas.draw()

    def _on_full_view(self) -> None:
        """Open the preview plot in a larger standalone window."""
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))
        step = self._ctrl.steps[self._selected_step]
        if step.diag_fn and step.output_sites:
            import inspect

            import pycsamt.emtools as et

            fn = getattr(et, step.diag_fn, None)
            if fn:
                try:
                    sig = inspect.signature(fn)
                    if "ax" in sig.parameters:
                        fn(step.output_sites, ax=ax, verbose=0)
                    else:
                        fn(step.output_sites, verbose=0)
                    plt.show()
                    return
                except Exception:
                    pass
        plt.close(fig)

    # ── Config save / load ────────────────────────────────────────────────────

    def _on_save_config(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save pipeline config", "", "JSON files (*.json)"
        )
        if not path:
            return
        cfg = self._ctrl.to_config()
        with open(path, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=2)
        self._log(f"Pipeline config saved to {path}")

    def _on_load_config(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load pipeline config", "", "JSON files (*.json)"
        )
        if not path:
            return
        with open(path, encoding="utf-8") as f:
            cfg = json.load(f)
        self._ctrl.from_config(cfg)
        self._refresh_stepper()
        self._select_step(self._selected_step)
        self._log(f"Pipeline config loaded from {path}")

    # ── Summary ───────────────────────────────────────────────────────────────

    def _update_summary(self) -> None:
        rows = []
        for step in self._ctrl.steps:
            icon = STATUS_ICON[step.status]
            info = step.result_info or step.error_msg or "—"
            t = f"{step.elapsed_s:.2f}s" if step.elapsed_s > 0 else "—"
            rows.append(
                f"<tr><td>{icon}</td><td><b>{step.name}</b></td>"
                f"<td>{info}</td><td>{t}</td></tr>"
            )
        html = (
            "<html><body style='font-family:monospace;font-size:11px'>"
            "<table border='0' cellspacing='4'>"
            "<tr><th></th><th>Step</th><th>Result</th><th>Time</th></tr>"
            + "".join(rows)
            + "</table></body></html>"
        )
        self._summary_browser.setHtml(html)

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _log(self, msg: str) -> None:
        import datetime

        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self._log_text.appendPlainText(f"[{ts}]  {msg}")

    # ── Geometry persistence ──────────────────────────────────────────────────

    def save_geometry_to(self, store: dict) -> None:
        store["pipeline_window"] = {
            "geometry": self.saveGeometry().toBase64().data().decode(),
            "visible": self.isVisible(),
        }

    def restore_geometry_from(self, store: dict) -> None:
        entry = store.get("pipeline_window")
        if not entry:
            return
        geo = entry.get("geometry")
        if geo:
            try:
                self.restoreGeometry(QByteArray.fromBase64(geo.encode()))
            except Exception:
                pass

    def closeEvent(self, event) -> None:
        self.hide()
        event.ignore()
