# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PipelineWindow — Pipeline Studio: build, run and review processing workflows.

A desktop front-end to the library pipeline engine (:mod:`pycsamt.pipeline`,
55+ registry steps, presets, YAML/JSON workflows, reports, run history).

┌──────────────────────────────────────────────────────────────────────────┐
│ Input: 28 stations · main window  [Load EDI/XML…]    │ Preset [▼][Apply] │
│                                      [Open…] [Save…]          [Library ▸]│
├───────────────────┬───────────────────────────────────┬──────────────────┤
│ WORKFLOW          │ [Step] [Dashboard] [Log] [History]│ STEP LIBRARY     │
│ 1 ☑ notch   Done  │  Power-line Harmonic Notch         │ 🔍 search…        │
│ 2 ☑ band  Running │  NR001 · noise_removal             │ ▸ noise_removal  │
│ 3 ☐ skew     Off  │  Parameters / More parameters     │ ▸ frequency      │
│ [+ Add][↑][↓][✕]  │  Result card · QC preview         │ …                │
│ [▶ Run all] …     │                                   │ [Add to workflow]│
│ [✓ Apply to main] │                                   │ ☐ AI steps       │
└───────────────────┴───────────────────────────────────┴──────────────────┘

Status badges are filled pills with ≥ 4.5 : 1 contrast on white text, so they
read on both themes (the old stepper coloured *text* green/yellow; yellow was
near-invisible in light mode).  The library drawer collapses (toggle or
Ctrl+L) so the step inspector can use the full width, like Forward's.

The run does not call ``Pipeline.run``: the worker executes the library's
own ``Step`` objects one by one, which gives a live "Running" state, Stop
between steps, and a per-step snapshot for QC previews.  Results reach the
main window only through the explicit **Apply to main data** button.
"""

from __future__ import annotations

from pycsamt.app.desktop.windows._base import refresh_icons

from pathlib import Path

from PySide6.QtCore import QByteArray, QSize, Qt, Signal, Slot
from PySide6.QtGui import QIcon, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.workflow_controller import (
    STATUS_STYLE,
    ParamField,
    RunStatus,
    WorkflowController,
    count_sites,
    parse_literal,
)
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.widgets.compact_button import compact_button

_ICONS = Path(__file__).parent.parent / "resources" / "icons"
_DEFAULT_PRESET = "full_processing"
_CODE_ROLE = Qt.ItemDataRole.UserRole
_CATEGORY_ROLE = Qt.ItemDataRole.UserRole + 1
_ACRONYMS = {"qc": "QC", "ai": "AI"}


def _category_title(cat: str) -> str:
    """'noise_removal' -> 'Noise Removal', 'qc' -> 'QC'."""
    return " ".join(_ACRONYMS.get(w, w.capitalize()) for w in cat.split("_"))


def _icon(name: str) -> QIcon:
    for c in (name, f"{name}.svg", f"{name}.png"):
        p = _ICONS / c
        if p.exists():
            return QIcon(str(p))
    return QIcon()


def _pill_qss(status: RunStatus) -> str:
    _text, colour = STATUS_STYLE[status]
    return (
        f"QLabel {{ background: {colour}; color: white; border-radius: 8px;"
        " padding: 1px 8px; font-weight: 600; font-size: 11px; }"
    )


class StatusPill(QLabel):
    """Filled, high-contrast status badge."""

    def __init__(self, status: RunStatus = RunStatus.PENDING) -> None:
        super().__init__()
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.set_status(status)

    def set_status(self, status: RunStatus) -> None:
        self.status = status
        self.setText(STATUS_STYLE[status][0])
        self.setStyleSheet(_pill_qss(status))


class _StepRow(QWidget):
    """One workflow row: index · enabled box · label/code · status pill."""

    def __init__(self, index: int, ws, on_toggle) -> None:
        super().__init__()
        h = QHBoxLayout(self)
        h.setContentsMargins(4, 3, 6, 3)
        h.setSpacing(6)
        self.num = QLabel(f"{index + 1}")
        self.num.setFixedWidth(18)
        self.num.setAlignment(Qt.AlignmentFlag.AlignRight
                              | Qt.AlignmentFlag.AlignVCenter)
        h.addWidget(self.num)
        self.check = QCheckBox()
        self.check.setChecked(ws.enabled)
        self.check.setToolTip("Include this step in runs")
        self.check.toggled.connect(lambda on: on_toggle(index, on))
        h.addWidget(self.check)
        text = QVBoxLayout()
        text.setSpacing(0)
        self.title = QLabel(ws.label)
        self.title.setStyleSheet("font-weight: 600;")
        self.sub = QLabel(f"{ws.code} · {ws.spec.label}")
        self.sub.setObjectName("InfoLabel")
        self.sub.setToolTip(f"{ws.code} · {ws.spec.label}")
        # Text yields width to the status pill (never the other way round):
        # long labels are clipped instead of pushing the pill off-screen.
        for lbl in (self.title, self.sub):
            lbl.setSizePolicy(QSizePolicy.Policy.Ignored,
                              QSizePolicy.Policy.Preferred)
            lbl.setMinimumWidth(0)
        text.addWidget(self.title)
        text.addWidget(self.sub)
        h.addLayout(text, 1)
        self.pill = StatusPill(ws.status)
        self.pill.setSizePolicy(QSizePolicy.Policy.Fixed,
                                QSizePolicy.Policy.Fixed)
        h.addWidget(self.pill)


# ── PipelineWindow ────────────────────────────────────────────────────────────


class PipelineWindow(QWidget):
    """
    Pipeline Studio window.

    Signals
    -------
    pipeline_finished(object)
        Emitted by **Apply to main data** with the workflow's output Sites.
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
        self.setWindowTitle("pycsamt — Pipeline Studio")
        ic = _icon("pipeline")
        if not ic.isNull():
            self.setWindowIcon(ic)
        self.resize(1280, 800)

        self._ctrl = WorkflowController()
        self._worker = None
        self._selected = -1
        self._rows: list[_StepRow] = []
        self._planned = 0
        self._finished = 0
        self._dark = False

        self._build_ui()
        try:
            self._ctrl.load_preset(_DEFAULT_PRESET)
            self._combo_preset.setCurrentIndex(
                max(self._combo_preset.findData(_DEFAULT_PRESET), 0)
            )
        except Exception:
            pass
        self._ctrl.reset_run()
        self._rebuild_rows(select=0)
        self._update_input_label()
        self._update_actions()

    # ══ UI construction ══════════════════════════════════════════════════════

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)
        root.addWidget(self._build_header())

        split = QSplitter(Qt.Orientation.Horizontal)
        split.setHandleWidth(4)
        split.setChildrenCollapsible(False)
        split.addWidget(self._build_workflow_panel())
        split.addWidget(self._build_center())
        self._lib_panel = self._build_library_panel()
        split.addWidget(self._lib_panel)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setStretchFactor(2, 0)
        split.setSizes([290, 700, 250])
        root.addWidget(split, 1)

        QShortcut(QKeySequence("Ctrl+L"), self, activated=self._toggle_library)
        QShortcut(QKeySequence("Ctrl+Up"), self,
                  activated=lambda: self._on_move(-1))
        QShortcut(QKeySequence("Ctrl+Down"), self,
                  activated=lambda: self._on_move(+1))

    def _build_header(self) -> QWidget:
        bar = QWidget()
        h = QHBoxLayout(bar)
        h.setContentsMargins(2, 0, 2, 0)
        h.setSpacing(6)

        self._input_lbl = QLabel("")
        self._input_lbl.setObjectName("InfoLabel")
        h.addWidget(self._input_lbl)
        btn_load = QPushButton("Load EDI/XML folder…")
        btn_load.setToolTip("Use a folder of EDI or EMTF-XML transfer "
                            "functions as the workflow input instead of the "
                            "main-window data (processed stations are "
                            "written as EDI)")
        btn_load.clicked.connect(self._on_load_folder)
        h.addWidget(btn_load)

        h.addWidget(self._vline())
        h.addWidget(QLabel("Preset:"))
        self._combo_preset = QComboBox()
        self._combo_preset.setMinimumContentsLength(16)
        self._combo_preset.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        for i, p in enumerate(self._ctrl.presets()):
            self._combo_preset.addItem(p.name.replace("_", " "), userData=p.name)
            self._combo_preset.setItemData(
                i, p.description, Qt.ItemDataRole.ToolTipRole
            )
        h.addWidget(self._combo_preset)
        btn_preset = QPushButton("Apply")
        btn_preset.setToolTip("Replace the workflow with this preset")
        btn_preset.clicked.connect(self._on_apply_preset)
        h.addWidget(btn_preset)

        h.addWidget(self._vline())
        btn_open = QPushButton("Open…")
        btn_open.setToolTip("Open a workflow (YAML/JSON — same format as the "
                            "CLI and Pipeline.from_yaml)")
        btn_open.clicked.connect(self._on_open)
        btn_save = QPushButton("Save…")
        btn_save.setToolTip("Save the workflow as YAML/JSON")
        btn_save.clicked.connect(self._on_save)
        h.addWidget(btn_open)
        h.addWidget(btn_save)
        h.addStretch(1)

        self._btn_library = QToolButton()
        self._btn_library.setCheckable(True)
        self._btn_library.setChecked(True)
        self._btn_library.setText("Library ▸")
        self._btn_library.setToolTip("Show / hide the step library  (Ctrl+L)")
        self._btn_library.toggled.connect(self._set_library_visible)
        h.addWidget(self._btn_library)
        self._header_buttons = [btn_load, btn_preset, btn_open, btn_save]
        return bar

    @staticmethod
    def _vline() -> QFrame:
        f = QFrame()
        f.setFrameShape(QFrame.Shape.VLine)
        f.setObjectName("Separator")
        return f

    # ── Left: workflow ────────────────────────────────────────────────────

    def _build_workflow_panel(self) -> QWidget:
        w = QWidget()
        w.setMinimumWidth(260)
        w.setMaximumWidth(360)
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(6)

        title = QLabel("WORKFLOW")
        title.setObjectName("PanelTitle")
        v.addWidget(title)

        self._step_list = QListWidget()
        self._step_list.setObjectName("StepperList")
        self._step_list.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection
        )
        self._step_list.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._step_list.currentRowChanged.connect(self._on_row_selected)
        v.addWidget(self._step_list, 1)

        edit_row = QHBoxLayout()
        self._btn_add = QPushButton("+ Add step")
        self._btn_add.setToolTip("Add the step selected in the Library below "
                                 "the current step")
        self._btn_add.clicked.connect(self._on_add_from_library)
        edit_row.addWidget(self._btn_add, 1)
        self._btn_up = compact_button(QPushButton("↑"))
        self._btn_up.setToolTip("Move step up  (Ctrl+Up)")
        self._btn_up.clicked.connect(lambda: self._on_move(-1))
        self._btn_dn = compact_button(QPushButton("↓"))
        self._btn_dn.setToolTip("Move step down  (Ctrl+Down)")
        self._btn_dn.clicked.connect(lambda: self._on_move(+1))
        self._btn_rm = compact_button(QPushButton("×"))
        self._btn_rm.setToolTip("Remove step")
        self._btn_rm.clicked.connect(self._on_remove)
        for b in (self._btn_up, self._btn_dn, self._btn_rm):
            edit_row.addWidget(b)
        v.addLayout(edit_row)

        run = QGroupBox("Run")
        rv = QVBoxLayout(run)
        self._btn_run_all = QPushButton("▶  Run all")
        self._btn_run_all.setObjectName("ComputeButton")
        self._btn_run_all.clicked.connect(lambda: self._start_run("all"))
        rv.addWidget(self._btn_run_all)
        r2 = QHBoxLayout()
        self._btn_run_from = QPushButton("Run from here")
        self._btn_run_from.setToolTip("Run the selected step and every "
                                      "step after it")
        self._btn_run_from.clicked.connect(lambda: self._start_run("from"))
        self._btn_run_one = QPushButton("Run step")
        self._btn_run_one.setToolTip("Run only the selected step")
        self._btn_run_one.clicked.connect(lambda: self._start_run("one"))
        r2.addWidget(self._btn_run_from)
        r2.addWidget(self._btn_run_one)
        rv.addLayout(r2)
        r3 = QHBoxLayout()
        self._btn_stop = QPushButton("■  Stop")
        self._btn_stop.setToolTip("Stop after the current step")
        self._btn_stop.clicked.connect(self._on_stop)
        self._btn_reset = QPushButton("Reset")
        self._btn_reset.setToolTip("Clear all run results")
        self._btn_reset.clicked.connect(self._on_reset)
        r3.addWidget(self._btn_stop)
        r3.addWidget(self._btn_reset)
        rv.addLayout(r3)
        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setMaximumHeight(12)
        self._progress.setTextVisible(False)
        rv.addWidget(self._progress)
        self._progress_lbl = QLabel("Ready")
        self._progress_lbl.setObjectName("InfoLabel")
        rv.addWidget(self._progress_lbl)
        v.addWidget(run)

        out = QGroupBox("Result")
        ov = QVBoxLayout(out)
        self._btn_apply_main = QPushButton("✓  Apply to main data")
        self._btn_apply_main.setObjectName("CommitButton")
        self._btn_apply_main.setToolTip(
            "Replace the main-window dataset with this workflow's output"
        )
        self._btn_apply_main.clicked.connect(self._on_apply_main)
        ov.addWidget(self._btn_apply_main)
        self._btn_export = QPushButton("Export results…")
        self._btn_export.setToolTip(
            "Write processed EDIs, reports and QC figures "
            "(formats: Settings ▸ API Configuration ▸ Pipeline)"
        )
        self._btn_export.clicked.connect(self._on_export)
        ov.addWidget(self._btn_export)
        self._chk_history = QCheckBox("Record runs in history")
        self._chk_history.setToolTip(
            "Append a summary of each finished run to the run-history file"
        )
        ov.addWidget(self._chk_history)
        v.addWidget(out)
        return w

    # ── Centre: inspector / dashboard / log / history ─────────────────────

    def _build_center(self) -> QWidget:
        self._tabs = QTabWidget()
        self._tabs.setDocumentMode(True)
        self._tabs.addTab(self._build_inspector(), "Step")

        self._dash_view = CanvasResultView(
            None, toolbar=True, empty_title="No run yet",
            empty_reason="Run the workflow to see step status, timing and "
                         "station flow.",
        )
        self._tabs.addTab(self._dash_view, "Dashboard")

        self._log_text = QPlainTextEdit()
        self._log_text.setReadOnly(True)
        self._log_text.setMaximumBlockCount(5000)
        self._log_text.setPlaceholderText("Run log will appear here…")
        self._tabs.addTab(self._log_text, "Log")

        hist = QWidget()
        hv = QVBoxLayout(hist)
        hv.setContentsMargins(4, 4, 4, 4)
        self._hist_table = QTableWidget(0, 6)
        self._hist_table.setHorizontalHeaderLabels(
            ["When (UTC)", "Workflow", "Steps", "Errors", "Stations",
             "Duration"]
        )
        self._hist_table.verticalHeader().setVisible(False)
        self._hist_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self._hist_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents
        )
        self._hist_table.horizontalHeader().setStretchLastSection(True)
        hv.addWidget(self._hist_table, 1)
        self._hist_empty = QLabel(
            "No recorded runs. Tick “Record runs in history” to log runs."
        )
        self._hist_empty.setObjectName("InfoLabel")
        hv.addWidget(self._hist_empty)
        btn_ref = QPushButton("Refresh")
        btn_ref.clicked.connect(self._refresh_history)
        hv.addWidget(btn_ref, 0, Qt.AlignmentFlag.AlignRight)
        self._tabs.addTab(hist, "History")
        self._tabs.currentChanged.connect(self._on_tab_changed)
        return self._tabs

    def _build_inspector(self) -> QWidget:
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        inner = QWidget()
        v = QVBoxLayout(inner)
        v.setContentsMargins(10, 8, 10, 8)
        v.setSpacing(8)

        self._insp_title = QLabel("")
        self._insp_title.setObjectName("StepTitle")
        self._insp_title.setStyleSheet(
            "QLabel#StepTitle { font-size: 15px; font-weight: 600; }"
        )
        self._insp_title.setWordWrap(True)
        v.addWidget(self._insp_title)
        self._insp_meta = QLabel("")
        self._insp_meta.setObjectName("InfoLabel")
        self._insp_meta.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
        )
        v.addWidget(self._insp_meta)
        self._insp_desc = QLabel("")
        self._insp_desc.setWordWrap(True)
        v.addWidget(self._insp_desc)

        self._grp_params = QGroupBox("Parameters")
        pv = QVBoxLayout(self._grp_params)
        self._form_basic = QFormLayout()
        pv.addLayout(self._form_basic)
        self._no_params = QLabel("This step has no default parameters.")
        self._no_params.setObjectName("InfoLabel")
        pv.addWidget(self._no_params)
        self._btn_adv = QToolButton()
        self._btn_adv.setCheckable(True)
        self._btn_adv.setText("More parameters")
        self._btn_adv.setToolButtonStyle(
            Qt.ToolButtonStyle.ToolButtonTextBesideIcon
        )
        self._btn_adv.setArrowType(Qt.ArrowType.RightArrow)
        self._btn_adv.toggled.connect(self._on_adv_toggled)
        pv.addWidget(self._btn_adv)
        self._adv_box = QWidget()
        self._form_adv = QFormLayout(self._adv_box)
        self._form_adv.setContentsMargins(0, 0, 0, 0)
        self._adv_box.setVisible(False)
        pv.addWidget(self._adv_box)
        self._btn_reset_params = QPushButton("Reset parameters")
        self._btn_reset_params.clicked.connect(self._on_reset_params)
        pv.addWidget(self._btn_reset_params, 0, Qt.AlignmentFlag.AlignLeft)
        v.addWidget(self._grp_params)

        self._grp_result = QGroupBox("Last run")
        rv = QVBoxLayout(self._grp_result)
        rr = QHBoxLayout()
        self._res_pill = StatusPill()
        rr.addWidget(self._res_pill)
        self._res_lbl = QLabel("")
        rr.addWidget(self._res_lbl, 1)
        rv.addLayout(rr)
        self._res_err = QLabel("")
        self._res_err.setWordWrap(True)
        self._res_err.setStyleSheet("color: #c92a2a;")
        self._res_err.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
        )
        rv.addWidget(self._res_err)
        v.addWidget(self._grp_result)

        self._grp_qc = QGroupBox("QC figure")
        qv = QVBoxLayout(self._grp_qc)
        qr = QHBoxLayout()
        self._combo_qc = QComboBox()
        qr.addWidget(self._combo_qc, 1)
        self._btn_qc = QPushButton("Show")
        self._btn_qc.setToolTip("Draw this QC figure for the step's output")
        self._btn_qc.clicked.connect(self._on_show_qc)
        qr.addWidget(self._btn_qc)
        qv.addLayout(qr)
        self._qc_view = CanvasResultView(
            None, toolbar=True, empty_title="No QC figure yet",
            empty_reason="Run the step, then choose a QC figure and click Show.",
        )
        self._qc_view.setMinimumHeight(320)
        qv.addWidget(self._qc_view)
        v.addWidget(self._grp_qc)

        self._insp_empty = QLabel(
            "The workflow is empty — choose a preset above, or add steps "
            "from the Library."
        )
        self._insp_empty.setWordWrap(True)
        self._insp_empty.setObjectName("InfoLabel")
        v.addWidget(self._insp_empty)
        v.addStretch(1)
        scroll.setWidget(inner)
        return scroll

    # ── Right: step library drawer ────────────────────────────────────────

    def _build_library_panel(self) -> QWidget:
        w = QWidget()
        w.setObjectName("LibraryPanel")
        w.setMinimumWidth(200)
        w.setMaximumWidth(320)
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(6)
        title = QLabel("STEP LIBRARY")
        title.setObjectName("PanelTitle")
        v.addWidget(title)
        self._lib_search = QLineEdit()
        self._lib_search.setPlaceholderText("Search steps…")
        self._lib_search.setClearButtonEnabled(True)
        self._lib_search.textChanged.connect(self._filter_library)
        v.addWidget(self._lib_search)
        self._lib_tree = QTreeWidget()
        self._lib_tree.setHeaderHidden(True)
        self._lib_tree.itemDoubleClicked.connect(
            lambda item, _c: self._on_add_from_library()
        )
        self._lib_tree.currentItemChanged.connect(self._on_lib_selected)
        v.addWidget(self._lib_tree, 1)
        self._lib_desc = QLabel("")
        self._lib_desc.setWordWrap(True)
        self._lib_desc.setObjectName("InfoLabel")
        self._lib_desc.setMinimumHeight(48)
        v.addWidget(self._lib_desc)
        btn = QPushButton("Add to workflow")
        btn.clicked.connect(self._on_add_from_library)
        v.addWidget(btn)
        self._chk_ai = QCheckBox("Show AI steps (experimental)")
        self._chk_ai.toggled.connect(self._on_ai_toggled)
        v.addWidget(self._chk_ai)
        self._populate_library()
        return w

    def _populate_library(self) -> None:
        self._lib_tree.clear()
        for cat in self._ctrl.categories():
            parent = QTreeWidgetItem([_category_title(cat)])
            parent.setData(0, _CATEGORY_ROLE, cat)
            parent.setFlags(Qt.ItemFlag.ItemIsEnabled)
            for spec in self._ctrl.catalogue(cat):
                item = QTreeWidgetItem([f"{spec.code}  {spec.label}"])
                item.setData(0, _CODE_ROLE, spec.code)
                item.setToolTip(0, f"{spec.name}()")
                parent.addChild(item)
            self._lib_tree.addTopLevelItem(parent)
        self._filter_library(self._lib_search.text())

    # ══ Public API ════════════════════════════════════════════════════════════

    def set_input_sites(self, sites) -> None:
        """Main-window data becomes the workflow input."""
        self._ctrl.set_input(sites, "main window")
        self._refresh_statuses()
        self._update_input_label()
        self._update_actions()

    def set_dark_mode(self, dark: bool) -> None:
        # Status pills are theme-independent (white on saturated colour).
        self._dark = dark
        refresh_icons(self, dark)

    # ══ Workflow list ═════════════════════════════════════════════════════════

    def _rebuild_rows(self, select: int | None = None) -> None:
        self._step_list.blockSignals(True)
        self._step_list.clear()
        self._rows = []
        for i, ws in enumerate(self._ctrl.steps):
            row = _StepRow(i, ws, self._on_toggle_step)
            item = QListWidgetItem()
            # height only: the list sizes row width to its viewport
            item.setSizeHint(QSize(1, row.sizeHint().height()))
            self._step_list.addItem(item)
            self._step_list.setItemWidget(item, row)
            self._rows.append(row)
        self._step_list.blockSignals(False)
        n = len(self._ctrl.steps)
        if n:
            target = min(max(select if select is not None else self._selected,
                             0), n - 1)
            self._step_list.setCurrentRow(target)
            self._on_row_selected(target)
        else:
            self._selected = -1
            self._show_inspector(None)
        self._update_actions()

    def _refresh_statuses(self) -> None:
        for row, ws in zip(self._rows, self._ctrl.steps):
            row.pill.set_status(ws.status)
        if 0 <= self._selected < len(self._ctrl.steps):
            self._show_result(self._ctrl.steps[self._selected])

    def _on_row_selected(self, row: int) -> None:
        if not 0 <= row < len(self._ctrl.steps):
            return
        self._selected = row
        self._show_inspector(row)
        self._update_actions()

    def _on_toggle_step(self, index: int, on: bool) -> None:
        self._ctrl.set_enabled(index, on)
        self._refresh_statuses()
        self._update_actions()

    def _on_move(self, delta: int) -> None:
        if self._running() or self._selected < 0:
            return
        new = self._ctrl.move_step(self._selected, delta)
        self._rebuild_rows(select=new)

    def _on_remove(self) -> None:
        if self._running() or self._selected < 0:
            return
        self._ctrl.remove_step(self._selected)
        self._rebuild_rows(select=self._selected)

    def _on_add_from_library(self) -> None:
        if self._running():
            return
        item = self._lib_tree.currentItem()
        code = item.data(0, _CODE_ROLE) if item is not None else None
        if not code:
            if not self._lib_panel.isVisible():
                self._set_library_visible(True)
            self._log("Select a step in the Library, then click Add.")
            return
        at = self._selected + 1 if self._selected >= 0 else None
        idx = self._ctrl.add_step(code, at)
        self._rebuild_rows(select=idx)
        self._log(f"Added [{code}] {self._ctrl.steps[idx].spec.label}.")

    # ══ Inspector ═════════════════════════════════════════════════════════════

    def _show_inspector(self, index: int | None) -> None:
        has = index is not None and 0 <= index < len(self._ctrl.steps)
        for w in (self._insp_title, self._insp_meta, self._insp_desc,
                  self._grp_params, self._grp_result, self._grp_qc):
            w.setVisible(has)
        self._insp_empty.setVisible(not has)
        if not has:
            return
        ws = self._ctrl.steps[index]
        spec = ws.spec
        self._insp_title.setText(f"{index + 1}.  {spec.label}")
        self._insp_meta.setText(
            f"{spec.code}  ·  {spec.category.replace('_', ' ')}  ·  "
            f"{spec.name}()   —   label “{ws.label}”"
        )
        self._insp_desc.setText(self._ctrl.describe(spec))
        self._build_param_forms(index)
        self._combo_qc.clear()
        for mod, fn in getattr(spec, "qc_defs", None) or []:
            self._combo_qc.addItem(fn, userData=(mod, fn))
        self._show_result(ws)
        self._qc_view.show_unavailable(
            "No QC figure yet",
            "Choose a QC figure and click Show." if self._combo_qc.count()
            else "This step defines no QC figure.",
        )

    def _show_result(self, ws) -> None:
        self._res_pill.set_status(ws.status)
        if ws.status in (RunStatus.DONE, RunStatus.ERROR):
            self._res_lbl.setText(
                f"{ws.n_in} → {ws.n_out} stations  ·  {ws.elapsed:.2f} s"
            )
        elif ws.status is RunStatus.OUTDATED:
            self._res_lbl.setText("Settings changed since the last run — "
                                  "run again to update.")
        elif ws.status is RunStatus.DISABLED:
            self._res_lbl.setText("Excluded from runs (unticked).")
        else:
            self._res_lbl.setText("Not run yet.")
        self._res_err.setText(ws.error)
        self._res_err.setVisible(bool(ws.error))
        can_qc = ws.status is RunStatus.DONE and self._combo_qc.count() > 0
        self._btn_qc.setEnabled(can_qc and not self._running())

    def _build_param_forms(self, index: int) -> None:
        for form in (self._form_basic, self._form_adv):
            while form.rowCount():
                form.removeRow(0)
        fields = self._ctrl.param_fields(index)
        basic = [f for f in fields if not f.advanced]
        adv = [f for f in fields if f.advanced]
        for f in basic:
            self._form_basic.addRow(f"{f.name}:", self._param_widget(index, f))
        for f in adv:
            self._form_adv.addRow(f"{f.name}:", self._param_widget(index, f))
        self._no_params.setVisible(not basic)
        self._btn_adv.setVisible(bool(adv))
        self._btn_adv.setText(f"More parameters ({len(adv)})")
        self._adv_box.setVisible(bool(adv) and self._btn_adv.isChecked())

    def _param_widget(self, index: int, f: ParamField) -> QWidget:
        def commit(value):
            self._ctrl.set_param(index, f.name, value)
            self._refresh_statuses()

        if f.kind == "bool":
            w = QCheckBox()
            w.setChecked(bool(f.value))
            w.toggled.connect(commit)
        elif f.kind == "int":
            w = QSpinBox()
            w.setRange(-1_000_000_000, 1_000_000_000)
            w.setValue(int(f.value))
            w.valueChanged.connect(commit)
        else:
            # float/str/literal: free text keeps 1e-3 and tuples exact
            w = QLineEdit("" if f.value is None else repr(f.value)
                          if f.kind != "str" else str(f.value))
            w.setPlaceholderText("None" if f.value is None else "")
            w.setToolTip(f"default: {f.default!r}")

            def on_edit(edit=w, kind=f.kind):
                text = edit.text()
                if kind == "float":
                    try:
                        value = float(text)
                    except ValueError:
                        edit.setStyleSheet("border: 1px solid #c92a2a;")
                        return
                elif kind == "str":
                    value = text
                else:
                    value = parse_literal(text)
                edit.setStyleSheet("")
                commit(value)

            w.editingFinished.connect(on_edit)
        w.setEnabled(not self._running())
        return w

    def _on_adv_toggled(self, on: bool) -> None:
        self._btn_adv.setArrowType(
            Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow
        )
        self._adv_box.setVisible(on)

    def _on_reset_params(self) -> None:
        if self._selected < 0 or self._running():
            return
        self._ctrl.reset_params(self._selected)
        self._build_param_forms(self._selected)
        self._refresh_statuses()

    def _on_show_qc(self) -> None:
        if self._selected < 0:
            return
        ws = self._ctrl.steps[self._selected]
        data = self._combo_qc.currentData()
        if ws.sites_after is None or not data:
            return
        mod, fn_name = data
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            import importlib

            from pycsamt.pipeline._steps import _to_figure, call_qc_fn

            fn = getattr(importlib.import_module(mod), fn_name)
            # before = this step's input: comparison plots need both states
            fig = _to_figure(call_qc_fn(
                fn, ws.sites_after, before=self._ctrl.input_for(self._selected)
            ))
            if fig is None:
                raise ValueError("the function returned no figure")
            self._qc_view.canvas.show_figure(fig)
            self._qc_view.show_canvas()
        except Exception as exc:
            self._qc_view.show_unavailable(
                "QC figure unavailable", f"{fn_name}: {exc}"
            )
        finally:
            QApplication.restoreOverrideCursor()

    # ══ Library drawer ════════════════════════════════════════════════════════

    @property
    def library_visible(self) -> bool:
        return self._btn_library.isChecked()

    def _toggle_library(self) -> None:
        self._btn_library.toggle()

    def _set_library_visible(self, visible: bool) -> None:
        if self._btn_library.isChecked() != visible:
            self._btn_library.setChecked(visible)
            return
        self._lib_panel.setVisible(visible)
        self._btn_library.setText("Library ▸" if visible else "◂ Library")

    def _filter_library(self, text: str) -> None:
        q = text.strip().lower()
        for i in range(self._lib_tree.topLevelItemCount()):
            parent = self._lib_tree.topLevelItem(i)
            any_shown = False
            for j in range(parent.childCount()):
                child = parent.child(j)
                hit = not q or q in child.text(0).lower() or q in (
                    child.toolTip(0).lower()
                )
                child.setHidden(not hit)
                any_shown |= hit
            parent.setHidden(not any_shown)
            parent.setExpanded(bool(q) and any_shown)

    def _on_lib_selected(self, item, _prev) -> None:
        code = item.data(0, _CODE_ROLE) if item is not None else None
        if not code:
            self._lib_desc.setText("")
            return
        from pycsamt.pipeline import lookup_step

        spec = lookup_step(code)
        desc = self._ctrl.describe(spec)
        self._lib_desc.setText(
            f"<b>{spec.label}</b><br>{desc[:220]}" if desc else spec.label
        )

    def _on_ai_toggled(self, on: bool) -> None:
        if on:
            codes = self._ctrl.enable_ai_steps()
            if codes:
                self._log("AI steps enabled: " + ", ".join(codes))
        self._populate_library()
        # hide/show the AI category (registration itself is permanent)
        for i in range(self._lib_tree.topLevelItemCount()):
            parent = self._lib_tree.topLevelItem(i)
            if parent.data(0, _CATEGORY_ROLE) == "ai":
                parent.setHidden(not on)

    # ══ Header actions ════════════════════════════════════════════════════════

    def _update_input_label(self) -> None:
        n = count_sites(self._ctrl.input_sites)
        if self._ctrl.input_sites is None:
            self._input_lbl.setText("Input: none — load data in the main "
                                    "window or a folder")
        else:
            src = self._ctrl.input_label or "main window"
            self._input_lbl.setText(f"Input: {n} stations · {src}")

    def _on_load_folder(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select a folder of EDI or EMTF-XML files")
        if not path:
            return
        try:
            from pycsamt.emtools import ensure_sites

            sites = ensure_sites(path)
        except Exception as exc:
            QMessageBox.warning(self, "Load failed", str(exc))
            return
        self._ctrl.set_input(sites, Path(path).name)
        self._refresh_statuses()
        self._update_input_label()
        self._update_actions()
        self._log(f"Input: {count_sites(sites)} stations from {path}")

    def _on_apply_preset(self) -> None:
        if self._running():
            return
        name = self._combo_preset.currentData()
        if not name:
            return
        if self._ctrl.steps and QMessageBox.question(
            self, "Replace workflow",
            f"Replace the current workflow with the “{name}” preset?",
        ) != QMessageBox.StandardButton.Yes:
            return
        self._ctrl.load_preset(name)
        self._ctrl.reset_run()
        self._rebuild_rows(select=0)
        self._log(f"Preset “{name}” loaded ({len(self._ctrl.steps)} steps).")

    def _on_open(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open workflow", "",
            "Workflows (*.yaml *.yml *.json);;All files (*)",
        )
        if not path:
            return
        try:
            self._ctrl.load(path)
        except Exception as exc:
            QMessageBox.warning(self, "Open failed", str(exc))
            return
        self._ctrl.reset_run()
        self._rebuild_rows(select=0)
        self._log(f"Workflow opened: {path}")

    def _on_save(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save workflow", f"{self._ctrl.name}.yaml",
            "YAML (*.yaml *.yml);;JSON (*.json)",
        )
        if not path:
            return
        try:
            self._ctrl.save(path)
        except Exception as exc:
            QMessageBox.warning(self, "Save failed", str(exc))
            return
        self._log(f"Workflow saved: {path}")

    # ══ Running ═══════════════════════════════════════════════════════════════

    def _running(self) -> bool:
        return self._worker is not None and self._worker.isRunning()

    def _start_run(self, mode: str) -> None:
        if self._running():
            return
        if self._ctrl.input_sites is None:
            self._log("No input data — load data in the main window or "
                      "click “Load EDI/XML folder…”.")
            self._tabs.setCurrentIndex(2)
            return
        indices = self._ctrl.plan(mode, max(self._selected, 0))
        if not indices:
            self._log("Nothing to run (no enabled steps selected).")
            return
        if mode != "all" and self._ctrl.input_for(indices[0]) is None:
            self._log("The earlier steps have not been run — use "
                      "“Run all” first.")
            self._tabs.setCurrentIndex(2)
            return
        from pycsamt.app.desktop.workers.workflow_worker import WorkflowWorker

        for i in indices:
            self._ctrl.steps[i].status = RunStatus.QUEUED
        # steps after a partial run no longer follow from fresh results
        if mode == "one":
            self._ctrl._invalidate_from(indices[0] + 1)
        self._planned, self._finished = len(indices), 0
        self._progress.setValue(0)
        self._progress_lbl.setText(f"0 / {self._planned} steps")
        self._refresh_statuses()
        self._log(f"Run started ({mode}): {self._planned} step(s).")

        self._worker = WorkflowWorker(self._ctrl, indices, parent=self)
        self._worker.step_started.connect(self._on_step_started)
        self._worker.step_finished.connect(self._on_step_finished)
        self._worker.log_line.connect(self._log)
        self._worker.run_finished.connect(self._on_run_finished)
        self._update_actions(running=True)
        self._worker.start()

    @Slot(int)
    def _on_step_started(self, idx: int) -> None:
        # Cross-thread signals are queued: a fast step may already be DONE
        # (set by the worker) when this arrives -- never overwrite that.
        ws = self._ctrl.steps[idx]
        if ws.status is RunStatus.QUEUED:
            ws.status = RunStatus.RUNNING
        self._refresh_statuses()

    @Slot(int)
    def _on_step_finished(self, idx: int) -> None:
        # Authoritative final state, on the GUI thread: "finished" is always
        # delivered after "started" for the same step, so this also repairs
        # a RUNNING written in the (tiny) window after the worker's DONE.
        ws = self._ctrl.steps[idx]
        if ws.result is not None:
            ws.status = (RunStatus.DONE if ws.result.error is None
                         else RunStatus.ERROR)
        self._finished += 1
        pct = int(100 * self._finished / max(self._planned, 1))
        self._progress.setValue(pct)
        self._progress_lbl.setText(f"{self._finished} / {self._planned} steps")
        self._refresh_statuses()

    @Slot(str)
    def _on_run_finished(self, outcome: str) -> None:
        # anything still queued did not run (the worker stops only between
        # steps, and sets DONE/ERROR itself, so RUNNING cannot be left over
        # except from a queued "started" signal -- see _on_step_started)
        for ws in self._ctrl.steps:
            if ws.status is RunStatus.QUEUED:
                ws.status = RunStatus.PENDING
        self._refresh_statuses()
        n_err = sum(ws.status is RunStatus.ERROR for ws in self._ctrl.steps)
        msg = {
            "completed": "Run finished",
            "stopped": "Run stopped",
            "aborted": "Run stopped at an error (policy: stop at first error)",
            "blocked": "Run could not continue",
        }.get(outcome, outcome)
        summary = (f"{msg} — {self._finished}/{self._planned} steps, "
                   f"{n_err} error(s), {self._ctrl.run_elapsed:.1f} s.")
        self._progress_lbl.setText(summary)
        self._log(summary)
        self._draw_dashboard()
        if self._chk_history.isChecked() and self._finished:
            try:
                p = self._ctrl.record_history()
                self._log(f"Run recorded in history: {p}")
            except Exception as exc:
                self._log(f"Could not record history: {exc}")
        self._worker = None
        self._update_actions()
        if 0 <= self._selected < len(self._ctrl.steps):
            self._show_result(self._ctrl.steps[self._selected])

    def _on_stop(self) -> None:
        if self._running():
            self._worker.requestInterruption()
            self._log("Stop requested — finishing the current step…")

    def _on_reset(self) -> None:
        if self._running():
            return
        self._ctrl.reset_run()
        self._refresh_statuses()
        self._progress.setValue(0)
        self._progress_lbl.setText("Ready")
        self._dash_view.show_unavailable(
            "No run yet", "Run the workflow to see the dashboard."
        )
        self._update_actions()
        self._log("Run results cleared.")

    # ══ Results ═══════════════════════════════════════════════════════════════

    def _draw_dashboard(self) -> None:
        result = self._ctrl.build_result()
        if not result.step_results:
            self._dash_view.show_unavailable(
                "No run yet", "Run the workflow to see the dashboard."
            )
            return
        from pycsamt.pipeline import (
            plot_pipeline_status,
            plot_pipeline_timing,
            plot_site_count_flow,
        )

        fig = self._dash_view.canvas.figure
        fig.clear()
        fig.set_layout_engine("constrained")
        try:
            # status | timing on top, station flow across the bottom; the
            # library titles embed the workflow name and collided when the
            # three panels sat side by side, so use short titles here.
            gs = fig.add_gridspec(2, 2)
            a1 = fig.add_subplot(gs[0, 0])
            a2 = fig.add_subplot(gs[0, 1])
            a3 = fig.add_subplot(gs[1, :])
            plot_pipeline_status(result, ax=a1)
            plot_pipeline_timing(result, ax=a2)
            plot_site_count_flow(result, ax=a3)
            for ax, title in ((a1, "Step status"), (a2, "Time per step"),
                              (a3, "Stations in → out")):
                ax.set_title(title, fontsize=10)
            fig.suptitle(f"Run of “{result.pipeline_name}”  ·  "
                         f"{result.elapsed_sec:.1f} s", fontsize=11)
            self._dash_view.canvas.draw()
            self._dash_view.show_canvas()
        except Exception as exc:
            fig.clear()
            self._dash_view.show_unavailable("Dashboard unavailable", str(exc))

    def _on_apply_main(self) -> None:
        out = self._ctrl.output_sites
        if out is None:
            return
        self.pipeline_finished.emit(out)
        self._log(f"Applied to main data: {count_sites(out)} stations.")

    def _on_export(self) -> None:
        if self._ctrl.output_sites is None:
            return
        path = QFileDialog.getExistingDirectory(self, "Export results to…")
        if not path:
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            written = self._ctrl.export(path)
        except Exception as exc:
            QApplication.restoreOverrideCursor()
            QMessageBox.warning(self, "Export failed", str(exc))
            return
        QApplication.restoreOverrideCursor()
        self._log(
            f"Exported to {written['root']}: {len(written['edis'])} EDI "
            f"file(s), {len(written['reports'])} report(s), "
            f"{len(written['figures'])} figure(s)."
        )

    def _refresh_history(self) -> None:
        rows = list(reversed(self._ctrl.load_history(last=100)))
        self._hist_table.setRowCount(len(rows))
        for r, rec in enumerate(rows):
            steps = rec.get("steps") or []
            cells = [
                str(rec.get("timestamp", "")).replace("T", " ").rstrip("Z"),
                str(rec.get("pipeline_name", "")),
                str(len(steps)),
                str(rec.get("n_errors", 0)),
                f"{rec.get('n_sites_in', 0)} → {rec.get('n_sites_out', 0)}",
                f"{float(rec.get('elapsed_sec', 0.0)):.1f} s",
            ]
            for c, text in enumerate(cells):
                self._hist_table.setItem(r, c, QTableWidgetItem(text))
        self._hist_empty.setVisible(not rows)

    def _on_tab_changed(self, index: int) -> None:
        if self._tabs.tabText(index) == "History":
            self._refresh_history()

    # ══ State ═════════════════════════════════════════════════════════════════

    def _update_actions(self, running: bool | None = None) -> None:
        running = self._running() if running is None else running
        has_steps = bool(self._ctrl.steps)
        sel = 0 <= self._selected < len(self._ctrl.steps)
        has_input = self._ctrl.input_sites is not None
        for b in self._header_buttons + [self._btn_add, self._btn_reset,
                                         self._btn_reset_params]:
            b.setEnabled(not running)
        self._step_list.setEnabled(not running)
        self._btn_up.setEnabled(not running and sel and self._selected > 0)
        self._btn_dn.setEnabled(
            not running and sel and self._selected < len(self._ctrl.steps) - 1
        )
        self._btn_rm.setEnabled(not running and sel)
        can_run = not running and has_steps and has_input
        self._btn_run_all.setEnabled(can_run)
        self._btn_run_from.setEnabled(can_run and sel)
        self._btn_run_one.setEnabled(can_run and sel)
        self._btn_stop.setEnabled(running)
        has_out = self._ctrl.output_sites is not None
        self._btn_apply_main.setEnabled(not running and has_out)
        self._btn_export.setEnabled(not running and has_out)
        for form in (self._form_basic, self._form_adv):
            for r in range(form.rowCount()):
                item = form.itemAt(r, QFormLayout.ItemRole.FieldRole)
                if item is not None and item.widget() is not None:
                    item.widget().setEnabled(not running)

    def _log(self, msg: str) -> None:
        import datetime

        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self._log_text.appendPlainText(f"[{ts}]  {msg}")

    # ── Geometry persistence ──────────────────────────────────────────────

    def save_geometry_to(self, store: dict) -> None:
        store["pipeline_window"] = {
            "geometry": self.saveGeometry().toBase64().data().decode(),
            "visible": self.isVisible(),
            "library_visible": self.library_visible,
            "record_history": self._chk_history.isChecked(),
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
        if "library_visible" in entry:
            self._set_library_visible(bool(entry["library_visible"]))
        if "record_history" in entry:
            self._chk_history.setChecked(bool(entry["record_history"]))

    def closeEvent(self, event) -> None:
        self.hide()
        event.ignore()


__all__ = ["PipelineWindow", "StatusPill"]
