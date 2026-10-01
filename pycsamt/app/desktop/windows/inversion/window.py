# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
InversionWindow — the Inversion Studio (desktop v2.6).

┌ Inversion Studio ─────────────────────────────────────────────────────────┐
│ [Classical|AI]  [Build & Run|View Results]      28 stations   [Library ▸] │
├──────────────┬──────────────────────────────────────────────┬─────────────┤
│ ENGINES      │ ① Data  ② Mesh  ③ Settings  ④ Run            │ PRESETS     │
│ Occam1D  ●   │  (step page)                                 │ RECENT RUNS │
│ Occam2D  ●   │                                              │             │
│ ModEM 2D ●   ├──────────────────────────────────────────────┤             │
│ ModEM 3D ●   │ [Running] Station 3/10  Iter 7/30  RMS 1.84  │             │
│ MARE2DEM ●   │ ▁▂▃▅                    02:13 · ~05:40 [■]   │             │
├──────────────┴──────────────────────────────────────────────┴─────────────┤
│ Console (terminal theme, searchable, ⇱ pops out into its own window)      │
└───────────────────────────────────────────────────────────────────────────┘

Every engine goes through the same Qt-free interface
(:mod:`pycsamt.app.desktop.controllers.inversion_engines`): its settings
become forms, **Build** writes the solver inputs into the run folder and
previews the mesh actually written, **Run** streams the solver console and
tracks iterations (Occam1D reports them itself; external solvers are
followed through their log files), and **View Results** reopens any run
folder -- from this window, a script, or the solver itself.

Closing the window only hides it: a run keeps going, and the main window's
status bar shows its progress.
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

from PySide6.QtCore import QByteArray, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QStackedWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers import inversion_runs as runs_store
from pycsamt.app.desktop.controllers.inversion_engines import (
    ENGINES,
    BuildInfo,
    LoadedRun,
    detect_engine,
    select_sites,
    site_names,
)
from pycsamt.app.desktop.panels.section_panel import SectionPanel
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.widgets.compact_button import compact_button
from pycsamt.app.desktop.windows._base import _icon, make_group
from pycsamt.app.desktop.windows._base import refresh_icons
from pycsamt.app.desktop.windows.inversion.ai_page import (
    AIInversionPage,
    draw_into,
)
from pycsamt.app.desktop.windows.inversion.console import (
    ConsolePanel,
    ConsoleWindow,
)
from pycsamt.app.desktop.windows.inversion.forms import SettingsForm
from pycsamt.app.desktop.windows.inversion.monitor import RunMonitor

STEPS = ("Data", "Mesh", "Settings", "Run")
ENGINE_ORDER = ("occam1d", "occam2d", "modem2d", "modem3d", "mare2dem")
_DIM_ENGINE = {"1D": "occam1d", "2D": "occam2d", "3D": "modem3d"}
_KEY_ROLE = Qt.ItemDataRole.UserRole

_SEG_QSS = (
    "QPushButton { border: 1px solid #8a94a3; padding: 4px 14px; "
    "background: transparent; font-weight: 600; }"
    "QPushButton:checked { background: #1864ab; color: white; "
    "border-color: #1864ab; }"
    "QPushButton#SegFirst { border-top-left-radius: 6px; "
    "border-bottom-left-radius: 6px; }"
    "QPushButton#SegLast { border-top-right-radius: 6px; "
    "border-bottom-right-radius: 6px; border-left: none; }"
)
_PILL = {
    "builtin": ("Built-in", "#2a7f3f"),
    "ready": ("Ready", "#2a7f3f"),
    "wsl": ("WSL", "#1864ab"),
    "missing": ("No binary", "#9c4f00"),
    "checking": ("Checking…", "#5f6b7a"),
}


def _segmented(labels: tuple[str, ...], parent: QWidget):
    box = QWidget(parent)
    h = QHBoxLayout(box)
    h.setContentsMargins(0, 0, 0, 0)
    h.setSpacing(0)
    group = QButtonGroup(box)
    group.setExclusive(True)
    buttons = []
    for i, label in enumerate(labels):
        b = QPushButton(label)
        b.setCheckable(True)
        b.setObjectName("SegFirst" if i == 0 else "SegLast"
                        if i == len(labels) - 1 else "SegMid")
        b.setStyleSheet(_SEG_QSS)
        group.addButton(b, i)
        h.addWidget(b)
        buttons.append(b)
    buttons[0].setChecked(True)
    return box, group, buttons


def _pill_qss(colour: str) -> str:
    return (f"QLabel {{ background: {colour}; color: white; border-radius: "
            "7px; padding: 0px 7px; font-weight: 600; font-size: 10px; }")


def _vline() -> QFrame:
    f = QFrame()
    f.setFrameShape(QFrame.Shape.VLine)
    f.setObjectName("Separator")
    return f


class _EngineRow(QWidget):
    def __init__(self, eng) -> None:
        super().__init__()
        v = QVBoxLayout(self)
        v.setContentsMargins(6, 5, 6, 5)
        v.setSpacing(1)
        top = QHBoxLayout()
        top.setSpacing(6)
        self.title = QLabel(eng.label)
        self.title.setStyleSheet("font-weight: 600;")
        top.addWidget(self.title)
        dim = QLabel(eng.dim)
        dim.setObjectName("InfoLabel")
        top.addWidget(dim)
        top.addStretch(1)
        self.pill = QLabel()
        top.addWidget(self.pill)
        v.addLayout(top)
        self.set_status("builtin" if eng.binary_key is None else "checking")

    def set_status(self, status: str, tip: str = "") -> None:
        self.status = status
        text, colour = _PILL[status]
        self.pill.setText(text)
        self.pill.setStyleSheet(_pill_qss(colour))
        self.pill.setToolTip(tip)


class InversionWindow(QWidget):
    """Inversion Studio: classical and AI inversion, build, run, review.

    Signals
    -------
    result_ready(dict)
        ``{"engine": key, "result": <library result>}`` when a run finishes
        or a result is sent to Interpretation.
    build_solver_requested(str)
        Solver Builder key ("occam2d", "modem2d", "modem3d", "mare2dem").
    run_state(str, int, bool)
        Status text, percent (-1 = busy) and whether a run is active; the
        main window mirrors it in its status bar.
    panel_closed()
        The window was closed (hidden; runs keep going).
    """

    result_ready = Signal(dict)
    build_solver_requested = Signal(str)
    run_state = Signal(str, int, bool)
    panel_closed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        flags = (Qt.WindowType.Window | Qt.WindowType.WindowCloseButtonHint
                 | Qt.WindowType.WindowMinimizeButtonHint
                 | Qt.WindowType.WindowMaximizeButtonHint)
        super().__init__(parent, flags)
        self.setWindowTitle("pycsamt — Inversion Studio")
        ic = _icon("inversion")
        if not ic.isNull():
            self.setWindowIcon(ic)
        self.resize(1240, 800)
        self._session_key = "inversion_window"
        self._sites = None
        self._dark = False
        self._starting_model: dict | None = None
        self._engine_key = "occam2d"
        self._forms: dict[str, dict[str, SettingsForm]] = {}
        self._build: BuildInfo | None = None
        self._build_stale = True
        self._worker = None  # build or run worker
        self._worker_kind = ""
        self._run_after_build = False
        self._run_engine = ""
        self._run_workdir: Path | None = None
        self._loaded: LoadedRun | None = None
        self._binaries: dict[str, list] = {}
        self._user_binaries: dict[str, str] = {}  # per solver key
        self._detector = None

        self._build_ui()
        self._select_engine("occam2d")
        self._poll = QTimer(self)
        self._poll.setInterval(2000)
        self._poll.timeout.connect(self._poll_history)
        self.refresh_solver_binaries()
        # The WSL probe (MARE2DEM) runs on first show, not at app start.
        self._probed = False

    # ══ UI ═══════════════════════════════════════════════════════════════
    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)
        root.addWidget(self._build_header())

        self._vsplit = QSplitter(Qt.Orientation.Vertical)
        self._vsplit.setChildrenCollapsible(False)
        body = QSplitter(Qt.Orientation.Horizontal)
        body.setChildrenCollapsible(False)
        self._mode_stack = QStackedWidget()
        self._classical_stack = QStackedWidget()
        self._classical_stack.addWidget(self._build_run_page())
        self._classical_stack.addWidget(self._build_results_page())
        self._mode_stack.addWidget(self._classical_stack)
        self._ai_page = AIInversionPage()
        self._ai_page.log.connect(self._log)
        self._ai_page.run_state.connect(self.run_state)
        self._ai_page.result_ready.connect(self.result_ready)
        self._mode_stack.addWidget(self._ai_page)
        body.addWidget(self._mode_stack)
        self._lib_panel = self._build_library()
        body.addWidget(self._lib_panel)
        body.setStretchFactor(0, 1)
        body.setStretchFactor(1, 0)
        body.setSizes([1000, 230])
        self._vsplit.addWidget(body)

        self._console = ConsolePanel()
        self._console.pop_out_requested.connect(self._pop_out_console)
        self._console_win = ConsoleWindow(self)
        self._console_win.dock_requested.connect(self._dock_console)
        self._vsplit.addWidget(self._console)
        self._vsplit.setStretchFactor(0, 3)
        self._vsplit.setStretchFactor(1, 1)
        self._vsplit.setSizes([560, 180])
        root.addWidget(self._vsplit, 1)
        # Back-compat alias: the pre-2.6 window's log widget
        self._log_edit = self._console.view

        QShortcut(QKeySequence("Ctrl+L"), self,
                  activated=lambda: self._btn_library.toggle())
        QShortcut(QKeySequence("Ctrl+R"), self, activated=self._on_run)

    def _build_header(self) -> QWidget:
        bar = QWidget()
        h = QHBoxLayout(bar)
        h.setContentsMargins(2, 0, 2, 0)
        h.setSpacing(8)
        box, self._mode_group, _ = _segmented(("Classical", "AI"), bar)
        self._mode_group.idClicked.connect(self._on_mode)
        h.addWidget(box)
        h.addWidget(_vline())
        self._page_box, self._page_group, _ = _segmented(
            ("Build && Run", "View Results"), bar)
        self._page_group.idClicked.connect(self._show_classical_page)
        h.addWidget(self._page_box)
        h.addStretch(1)
        self._data_lbl = QLabel("No data loaded")
        self._data_lbl.setObjectName("InfoLabel")
        h.addWidget(self._data_lbl)
        self._btn_library = QToolButton()
        self._btn_library.setCheckable(True)
        self._btn_library.setChecked(True)
        self._btn_library.setText("Library ▸")
        self._btn_library.setToolTip("Show / hide presets and recent runs "
                                     "(Ctrl+L)")
        self._btn_library.toggled.connect(self._set_library_visible)
        h.addWidget(self._btn_library)
        return bar

    # ── Build & Run page ──────────────────────────────────────────────
    def _build_run_page(self) -> QWidget:
        page = QSplitter(Qt.Orientation.Horizontal)
        page.setChildrenCollapsible(False)

        left = QWidget()
        left.setMinimumWidth(190)
        left.setMaximumWidth(260)
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.setSpacing(4)
        cap = QLabel("ENGINES")
        cap.setObjectName("InfoLabel")
        lv.addWidget(cap)
        self._engine_list = QListWidget()
        self._engine_list.setObjectName("EngineList")
        self._engine_rows: dict[str, _EngineRow] = {}
        for key in ENGINE_ORDER:
            eng = ENGINES[key]
            item = QListWidgetItem()
            item.setData(_KEY_ROLE, key)
            item.setToolTip(eng.description)
            row = _EngineRow(eng)
            item.setSizeHint(row.sizeHint())
            self._engine_list.addItem(item)
            self._engine_list.setItemWidget(item, row)
            self._engine_rows[key] = row
        self._engine_list.currentItemChanged.connect(
            lambda cur, _prev: cur and self._select_engine(
                cur.data(_KEY_ROLE)))
        lv.addWidget(self._engine_list, 1)
        self._engine_desc = QLabel("")
        self._engine_desc.setWordWrap(True)
        self._engine_desc.setObjectName("InfoLabel")
        lv.addWidget(self._engine_desc)
        page.addWidget(left)

        center = QWidget()
        cv = QVBoxLayout(center)
        cv.setContentsMargins(0, 0, 0, 0)
        cv.setSpacing(4)
        steps = QHBoxLayout()
        steps.setSpacing(4)
        self._step_group = QButtonGroup(self)
        self._step_buttons: list[QPushButton] = []
        for i, name in enumerate(STEPS):
            b = QPushButton(f"{i + 1}  {name}")
            b.setCheckable(True)
            b.setStyleSheet(
                "QPushButton { padding: 5px 12px; border-radius: 12px; "
                "border: 1px solid #8a94a3; background: transparent; }"
                "QPushButton:checked { background: #1864ab; color: white; "
                "border-color: #1864ab; font-weight: 600; }")
            self._step_group.addButton(b, i)
            steps.addWidget(b)
            self._step_buttons.append(b)
        self._step_buttons[0].setChecked(True)
        self._step_group.idClicked.connect(self._go_step)
        steps.addStretch(1)
        self._build_state_lbl = QLabel("")
        self._build_state_lbl.setObjectName("InfoLabel")
        steps.addWidget(self._build_state_lbl)
        cv.addLayout(steps)

        self._step_stack = QStackedWidget()
        self._step_stack.addWidget(self._build_data_step())
        self._step_stack.addWidget(self._build_mesh_step())
        self._step_stack.addWidget(self._build_settings_step())
        self._step_stack.addWidget(self._build_launch_step())
        cv.addWidget(self._step_stack, 1)

        nav = QHBoxLayout()
        self._btn_back = QPushButton("← Back")
        self._btn_back.clicked.connect(lambda: self._go_step(
            self._step_stack.currentIndex() - 1))
        self._btn_next = QPushButton("Next →")
        self._btn_next.clicked.connect(lambda: self._go_step(
            self._step_stack.currentIndex() + 1))
        nav.addWidget(self._btn_back)
        nav.addStretch(1)
        nav.addWidget(self._btn_next)
        cv.addLayout(nav)

        self._monitor = RunMonitor()
        self._monitor.stop_requested.connect(self._on_stop)
        cv.addWidget(self._monitor)
        self._btn_stop = self._monitor.btn_stop
        page.addWidget(center)
        page.setStretchFactor(1, 1)
        page.setSizes([220, 900])
        return page

    @staticmethod
    def _scroll(widget: QWidget) -> QScrollArea:
        s = QScrollArea()
        s.setWidgetResizable(True)
        s.setFrameShape(QScrollArea.Shape.NoFrame)
        s.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        s.setWidget(widget)
        return s

    def _section_host(self, section: str) -> QStackedWidget:
        host = QStackedWidget()
        host.setProperty("section", section)
        host.setSizePolicy(QSizePolicy.Policy.Preferred,
                           QSizePolicy.Policy.Maximum)
        return host

    def _build_data_step(self) -> QWidget:
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        grp, gl = make_group("Stations")
        self._station_list = QListWidget()
        self._station_list.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection)
        self._station_list.itemChanged.connect(self._on_stations_changed)
        gl.addWidget(self._station_list, 1)
        row = QHBoxLayout()
        b_all = QPushButton("All")
        b_none = QPushButton("None")
        b_sel = QPushButton("Only selected")
        b_sel.setToolTip("Tick only the highlighted rows")
        for b in (b_all, b_none, b_sel):
            row.addWidget(b)
        b_all.clicked.connect(lambda: self._tick_all(True))
        b_none.clicked.connect(lambda: self._tick_all(False))
        b_sel.clicked.connect(self._tick_highlighted)
        gl.addLayout(row)
        self._station_count = QLabel("")
        self._station_count.setObjectName("InfoLabel")
        gl.addWidget(self._station_count)
        h.addWidget(grp, 1)

        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        grp, gl = make_group("Frequency band")
        self._chk_band = QCheckBox("Limit the frequency band")
        self._chk_band.toggled.connect(self._on_band_toggled)
        gl.addWidget(self._chk_band)
        f = QFormLayout()
        self._f_min = QDoubleSpinBox()
        self._f_min.setDecimals(5)
        self._f_min.setRange(1e-6, 1e6)
        self._f_min.setValue(1e-3)
        self._f_min.setSuffix(" Hz")
        self._f_max = QDoubleSpinBox()
        self._f_max.setDecimals(2)
        self._f_max.setRange(1e-5, 1e7)
        self._f_max.setValue(1e4)
        self._f_max.setSuffix(" Hz")
        for sb in (self._f_min, self._f_max):
            sb.setEnabled(False)
            sb.valueChanged.connect(self._mark_stale)
        f.addRow("f min:", self._f_min)
        f.addRow("f max:", self._f_max)
        gl.addLayout(f)
        rv.addWidget(grp)
        grp, gl = make_group("Data weighting")
        self._data_host = self._section_host("Data")
        gl.addWidget(self._data_host)
        rv.addWidget(grp)
        grp, gl = make_group("Starting model")
        self._fwd_model_label = QLabel("(no model from Forward)")
        self._fwd_model_label.setWordWrap(True)
        self._fwd_model_label.setObjectName("FwdModelLabel")
        gl.addWidget(self._fwd_model_label)
        btn_clear = QPushButton("Clear Forward model")
        btn_clear.clicked.connect(self._clear_starting_model)
        gl.addWidget(btn_clear)
        rv.addWidget(grp)
        self._data_check = QLabel("")
        self._data_check.setWordWrap(True)
        rv.addWidget(self._data_check)
        rv.addStretch(1)
        h.addWidget(self._scroll(right), 1)
        return w

    def _build_mesh_step(self) -> QWidget:
        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 4, 0)
        grp, gl = make_group("Mesh")
        self._mesh_host = self._section_host("Mesh")
        gl.addWidget(self._mesh_host)
        lv.addWidget(grp)
        self._btn_build = QPushButton("Build inputs && preview mesh")
        self._btn_build.setToolTip("Write the solver input files into the "
                                   "run folder and draw the mesh they use")
        self._btn_build.clicked.connect(self._on_build)
        lv.addWidget(self._btn_build)
        grp, gl = make_group("Build summary")
        self._summary_form = QFormLayout()
        self._summary_form.setSpacing(3)
        gl.addLayout(self._summary_form)
        self._build_warn = QLabel("")
        self._build_warn.setWordWrap(True)
        self._build_warn.setStyleSheet("color: #9c4f00;")
        gl.addWidget(self._build_warn)
        lv.addWidget(grp)
        lv.addStretch(1)
        scroll = self._scroll(left)
        scroll.setMinimumWidth(260)
        scroll.setMaximumWidth(360)
        split.addWidget(scroll)
        self._mesh_view = CanvasResultView(
            empty_title="No mesh yet",
            empty_reason="Click “Build inputs & preview mesh” to write the "
                         "solver files and see the grid they define.")
        split.addWidget(self._mesh_view)
        split.setStretchFactor(1, 1)
        return split

    def _build_settings_step(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        row.addWidget(QLabel("Preset:"))
        self._preset_combo = QComboBox()
        self._preset_combo.setMinimumContentsLength(14)
        row.addWidget(self._preset_combo, 1)
        b = QPushButton("Apply")
        b.clicked.connect(lambda: self._apply_preset(
            self._preset_combo.currentText()))
        row.addWidget(b)
        b = QPushButton("Save as…")
        b.clicked.connect(self._save_preset)
        row.addWidget(b)
        b = QPushButton("Defaults")
        b.setToolTip("Reset every setting of this engine")
        b.clicked.connect(self._reset_engine_settings)
        row.addWidget(b)
        v.addLayout(row)
        grp, gl = make_group("Inversion settings")
        self._settings_host = self._section_host("Settings")
        gl.addWidget(self._settings_host)
        v.addWidget(grp)
        v.addStretch(1)
        return self._scroll(w)

    def _build_launch_step(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        grp, gl = make_group("Run folder")
        row = QHBoxLayout()
        self._workdir_edit = QLineEdit()
        self._workdir_edit.setPlaceholderText(
            "Empty = a new folder under ~/pycsamt-inversions")
        self._workdir_edit.editingFinished.connect(self._mark_stale)
        row.addWidget(self._workdir_edit, 1)
        b = QPushButton("Browse…")
        b.clicked.connect(self._browse_workdir)
        row.addWidget(b)
        b = QPushButton("Open")
        b.setToolTip("Open the run folder in the file browser")
        b.clicked.connect(lambda: self._reveal(self._current_workdir(False)))
        row.addWidget(b)
        gl.addLayout(row)
        v.addWidget(grp)

        self._bin_group, gl = make_group("Solver binary")
        row = QHBoxLayout()
        self._binary_combo = QComboBox()
        self._binary_combo.setEditable(True)
        self._binary_combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self._binary_combo.lineEdit().setPlaceholderText(
            "Path to the solver executable")
        self._binary_combo.currentIndexChanged.connect(self._on_binary_pick)
        self._binary_combo.editTextChanged.connect(self._on_binary_edited)
        row.addWidget(self._binary_combo, 1)
        b = QPushButton("…")
        compact_button(b, 28)
        b.setToolTip("Browse for the executable")
        b.clicked.connect(self._browse_binary)
        row.addWidget(b)
        self._btn_detect = QPushButton("Detect")
        self._btn_detect.setToolTip("Search the Solver Builder registry, "
                                    "the source tree, PATH and WSL")
        self._btn_detect.clicked.connect(self._detect_binaries_background)
        row.addWidget(self._btn_detect)
        self._btn_build_solver = QPushButton("Build…")
        self._btn_build_solver.setToolTip("Compile this solver with the "
                                          "Solver Builder")
        self._btn_build_solver.clicked.connect(
            lambda: self.build_solver_requested.emit(
                ENGINES[self._engine_key].binary_key or ""))
        row.addWidget(self._btn_build_solver)
        gl.addLayout(row)
        self._binary_origin = QLabel("")
        self._binary_origin.setObjectName("InfoLabel")
        gl.addWidget(self._binary_origin)
        row = QHBoxLayout()
        self._chk_stage = QCheckBox("Copy the binary into the run folder")
        self._chk_stage.setToolTip("Makes the run folder self-contained: "
                                   "it can be re-run or archived without "
                                   "this machine's toolchain")
        row.addWidget(self._chk_stage)
        self._rb_copy = QRadioButton("copy")
        self._rb_move = QRadioButton("move")
        self._rb_copy.setChecked(True)
        self._rb_move.setToolTip("Moves the file; the Solver Builder "
                                 "registry is updated to the new place")
        for rb in (self._rb_copy, self._rb_move):
            rb.setEnabled(False)
            row.addWidget(rb)
        self._chk_stage.toggled.connect(self._rb_copy.setEnabled)
        self._chk_stage.toggled.connect(self._rb_move.setEnabled)
        row.addStretch(1)
        gl.addLayout(row)
        v.addWidget(self._bin_group)

        grp, gl = make_group("Parallel run")
        self._run_host = self._section_host("Run")
        gl.addWidget(self._run_host)
        self._run_group = grp
        v.addWidget(grp)

        grp, gl = make_group("Summary")
        self._run_summary = QLabel("")
        self._run_summary.setWordWrap(True)
        self._run_summary.setTextFormat(Qt.TextFormat.RichText)
        gl.addWidget(self._run_summary)
        v.addWidget(grp)
        row = QHBoxLayout()
        row.addStretch(1)
        self._btn_run = QPushButton("▶  Run inversion")
        self._btn_run.setObjectName("RunButton")
        self._btn_run.setMinimumWidth(180)
        self._btn_run.setToolTip("Build the inputs if needed, then run "
                                 "(Ctrl+R)")
        self._btn_run.clicked.connect(self._on_run)
        row.addWidget(self._btn_run)
        v.addLayout(row)
        self._run_status = QLabel("")
        self._run_status.setObjectName("RunStatusLabel")
        v.addWidget(self._run_status)
        v.addStretch(1)
        return self._scroll(w)

    # ── View Results page ─────────────────────────────────────────────
    def _build_results_page(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        row = QHBoxLayout()
        self._btn_open_run = QPushButton("Open run folder…")
        self._btn_open_run.clicked.connect(self._browse_run_folder)
        row.addWidget(self._btn_open_run)
        self._res_path = QLabel("No run loaded")
        self._res_path.setObjectName("InfoLabel")
        self._res_path.setSizePolicy(QSizePolicy.Policy.Ignored,
                                     QSizePolicy.Policy.Preferred)
        row.addWidget(self._res_path, 1)
        self._res_station = QComboBox()
        self._res_station.setToolTip("Station (Occam1D inverts each station "
                                     "separately)")
        self._res_station.currentIndexChanged.connect(self._reload_result)
        self._res_iter = QComboBox()
        self._res_iter.setToolTip("Iteration to show")
        self._res_iter.currentIndexChanged.connect(self._reload_result)
        row.addWidget(QLabel("Station:"))
        row.addWidget(self._res_station)
        row.addWidget(QLabel("Iteration:"))
        row.addWidget(self._res_iter)
        b = QPushButton("Reload")
        b.setToolTip("Read the folder again (a run may still be writing)")
        b.clicked.connect(lambda: self._loaded and self.open_run(
            self._loaded.path, self._loaded.engine))
        row.addWidget(b)
        self._btn_export = QPushButton("Export PCSF…")
        self._btn_export.setToolTip("Save the model in the pycsamt PCSF "
                                    "format (Map View, conversions)")
        self._btn_export.clicked.connect(self._export_pcsf)
        row.addWidget(self._btn_export)
        self._btn_send = QPushButton("Send to Interpretation")
        self._btn_send.clicked.connect(self._send_to_interpretation)
        row.addWidget(self._btn_send)
        v.addLayout(row)

        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)
        left = QWidget()
        left.setMinimumWidth(190)
        left.setMaximumWidth(260)
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        grp, gl = make_group("Run")
        self._res_summary = QFormLayout()
        self._res_summary.setSpacing(3)
        gl.addLayout(self._res_summary)
        lv.addWidget(grp)
        cap = QLabel("VIEWS")
        cap.setObjectName("InfoLabel")
        lv.addWidget(cap)
        self._view_list = QListWidget()
        self._view_list.currentItemChanged.connect(
            lambda cur, _p: cur and self._render_view(cur.data(_KEY_ROLE)))
        lv.addWidget(self._view_list, 1)
        split.addWidget(left)
        self._res_stack = QStackedWidget()
        self._result_view = CanvasResultView(
            empty_title="No run loaded",
            empty_reason="Open a run folder (Occam1D, Occam2D, ModEM, "
                         "MARE2DEM) or finish a run in Build & Run.")
        self._tab_section = SectionPanel()
        self._res_stack.addWidget(self._result_view)
        self._res_stack.addWidget(self._tab_section)
        split.addWidget(self._res_stack)
        split.setStretchFactor(1, 1)
        split.setSizes([220, 900])
        v.addWidget(split, 1)
        self._btn_export.setEnabled(False)
        self._btn_send.setEnabled(False)
        return w

    # ── Library drawer ────────────────────────────────────────────────
    def _build_library(self) -> QWidget:
        w = QWidget()
        w.setMinimumWidth(200)
        w.setMaximumWidth(300)
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        cap = QLabel("PRESETS")
        cap.setObjectName("InfoLabel")
        v.addWidget(cap)
        self._preset_list = QListWidget()
        self._preset_list.itemDoubleClicked.connect(
            lambda it: self._apply_preset(it.data(_KEY_ROLE)))
        v.addWidget(self._preset_list, 1)
        row = QHBoxLayout()
        b = QPushButton("Apply")
        b.clicked.connect(lambda: self._preset_list.currentItem() and
                          self._apply_preset(self._preset_list.currentItem()
                                             .data(_KEY_ROLE)))
        row.addWidget(b)
        self._btn_del_preset = QPushButton("Delete")
        self._btn_del_preset.clicked.connect(self._delete_preset)
        row.addWidget(self._btn_del_preset)
        v.addLayout(row)
        cap = QLabel("RECENT RUNS")
        cap.setObjectName("InfoLabel")
        v.addWidget(cap)
        self._recent_list = QListWidget()
        self._recent_list.setToolTip("Double-click to open in View Results")
        self._recent_list.itemDoubleClicked.connect(
            lambda it: self.open_run(it.data(_KEY_ROLE)["workdir"],
                                     it.data(_KEY_ROLE)["engine"]))
        v.addWidget(self._recent_list, 2)
        row = QHBoxLayout()
        b = QPushButton("Open")
        b.clicked.connect(lambda: self._recent_list.currentItem() and
                          self.open_run(self._recent_list.currentItem()
                                        .data(_KEY_ROLE)["workdir"]))
        row.addWidget(b)
        b = QPushButton("Folder")
        b.clicked.connect(lambda: self._recent_list.currentItem() and
                          self._reveal(self._recent_list.currentItem()
                                       .data(_KEY_ROLE)["workdir"]))
        row.addWidget(b)
        v.addLayout(row)
        self._refresh_recent()
        return w

    # ══ engine & forms ═══════════════════════════════════════════════════
    def _engine(self):
        return ENGINES[self._engine_key]

    def _engine_forms(self, key: str) -> dict[str, SettingsForm]:
        if key not in self._forms:
            eng = ENGINES[key]
            forms = {}
            for host in (self._data_host, self._mesh_host,
                         self._settings_host, self._run_host):
                section = host.property("section")
                fields = [f for f in eng.fields() if f.section == section]
                form = SettingsForm(fields)
                form.changed.connect(self._mark_stale)
                host.addWidget(form)
                forms[section] = form
            self._forms[key] = forms
        return self._forms[key]

    def _select_engine(self, key: str) -> None:
        if key not in ENGINES:
            return
        self._engine_key = key
        eng = ENGINES[key]
        forms = self._engine_forms(key)
        for host in (self._data_host, self._mesh_host, self._settings_host,
                     self._run_host):
            form = forms[host.property("section")]
            host.setCurrentWidget(form)
            host.setVisible(bool(form.fields))
        self._run_group.setVisible(bool(forms["Run"].fields))
        for i in range(self._engine_list.count()):
            item = self._engine_list.item(i)
            if item.data(_KEY_ROLE) == key:
                self._engine_list.blockSignals(True)
                self._engine_list.setCurrentItem(item)
                self._engine_list.blockSignals(False)
        self._engine_desc.setText(eng.description)
        self._bin_group.setVisible(eng.binary_key is not None)
        self._build = None
        self._build_stale = True
        self._mesh_view.show_unavailable(
            "No mesh yet", f"Build the {eng.label} inputs to preview the "
                           "mesh.")
        self._populate_presets()
        self._fill_binary_combo()
        self._update_data_check()
        self._update_run_summary()

    def values(self) -> dict:
        """All settings of the current engine (+ data window)."""
        out = {}
        for form in self._engine_forms(self._engine_key).values():
            out.update(form.values())
        if self._chk_band.isChecked():
            out["freq_min"] = float(self._f_min.value())
            out["freq_max"] = float(self._f_max.value())
        return out

    def _mark_stale(self, *_) -> None:
        if self._build is not None:
            self._build_stale = True
            self._build_state_lbl.setText("Settings changed — inputs will "
                                          "be rebuilt on Run")
        self._update_run_summary()

    def _on_band_toggled(self, on: bool) -> None:
        self._f_min.setEnabled(on)
        self._f_max.setEnabled(on)
        self._mark_stale()

    # ── navigation ────────────────────────────────────────────────────
    def _go_step(self, i: int) -> None:
        i = max(0, min(i, len(STEPS) - 1))
        self._step_stack.setCurrentIndex(i)
        self._step_buttons[i].setChecked(True)
        self._btn_back.setEnabled(i > 0)
        self._btn_next.setEnabled(i < len(STEPS) - 1)
        if i == 3:
            self._update_run_summary()

    def _on_mode(self, i: int) -> None:
        self._mode_stack.setCurrentIndex(i)
        self._page_box.setVisible(i == 0)
        if i == 1:
            self._ai_page.set_sites(self._sites, self._ticked_names())

    def _show_classical_page(self, i: int) -> None:
        self._classical_stack.setCurrentIndex(i)
        self._page_group.button(i).setChecked(True)

    def _set_library_visible(self, on: bool) -> None:
        self._lib_panel.setVisible(on)
        self._btn_library.setText("Library ▸" if on else "◂ Library")

    @property
    def library_visible(self) -> bool:
        return self._lib_panel.isVisible() or self._btn_library.isChecked()

    # ══ data ═════════════════════════════════════════════════════════════
    def set_sites(self, sites) -> None:
        """Called by MainWindow when data is loaded or changed."""
        self._sites = sites
        self._station_list.blockSignals(True)
        self._station_list.clear()
        for name in site_names(sites):
            item = QListWidgetItem(name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked)
            self._station_list.addItem(item)
        self._station_list.blockSignals(False)
        self._on_stations_changed()
        self._ai_page.set_sites(sites, self._ticked_names())

    def _ticked_names(self) -> list[str]:
        return [self._station_list.item(i).text()
                for i in range(self._station_list.count())
                if self._station_list.item(i).checkState()
                == Qt.CheckState.Checked]

    def _selected_sites(self):
        names = self._ticked_names()
        if self._sites is None:
            return None
        if len(names) == self._station_list.count():
            return self._sites
        return select_sites(self._sites, names) if names else None

    def _tick_all(self, on: bool) -> None:
        state = Qt.CheckState.Checked if on else Qt.CheckState.Unchecked
        self._station_list.blockSignals(True)
        for i in range(self._station_list.count()):
            self._station_list.item(i).setCheckState(state)
        self._station_list.blockSignals(False)
        self._on_stations_changed()

    def _tick_highlighted(self) -> None:
        chosen = {it.text() for it in self._station_list.selectedItems()}
        if not chosen:
            return
        self._station_list.blockSignals(True)
        for i in range(self._station_list.count()):
            it = self._station_list.item(i)
            it.setCheckState(Qt.CheckState.Checked if it.text() in chosen
                             else Qt.CheckState.Unchecked)
        self._station_list.blockSignals(False)
        self._on_stations_changed()

    def _on_stations_changed(self, *_) -> None:
        n_all = self._station_list.count()
        n = len(self._ticked_names())
        self._station_count.setText(f"{n} of {n_all} stations used")
        self._data_lbl.setText(f"{n_all} stations loaded" if n_all
                               else "No data loaded")
        self._mark_stale()
        self._update_data_check()
        self._ai_page.set_sites(self._sites, self._ticked_names())

    def _update_data_check(self) -> None:
        msg = self._engine().check_sites(self._selected_sites())
        if msg:
            self._data_check.setText(f"⚠ {msg}")
            self._data_check.setStyleSheet("color: #9c4f00;")
        else:
            self._data_check.setText(f"✓ Ready for {self._engine().label}.")
            self._data_check.setStyleSheet("color: #2a7f3f;")

    # ── starting model from Forward ───────────────────────────────────
    def load_starting_model(self, payload: dict) -> None:
        """Called by MainWindow when ForwardModelWindow sends a model."""
        self._starting_model = payload
        dim = payload.get("dim", "?")
        rho_list = payload.get("resistivity", []) or []
        label = f"Forward model: dim={dim}"
        if rho_list:
            label += f", {len(rho_list)} layers"
        self._fwd_model_label.setText(label)
        key = _DIM_ENGINE.get(str(dim).upper())
        if key:
            self._select_engine(key)
        if rho_list:
            rho0 = float(rho_list[0])
            for forms in self._forms.values():
                for form in forms.values():
                    for k in ("starting_resistivity", "initial_rho"):
                        if k in form.widgets:
                            form.set_values({k: rho0})
        self._mode_group.button(0).setChecked(True)
        self._on_mode(0)
        self._show_classical_page(0)
        self._go_step(0)
        self._log(f"Starting model received from Forward ({label}).")

    def _clear_starting_model(self) -> None:
        self._starting_model = None
        self._fwd_model_label.setText("(no model from Forward)")

    # ══ binaries ═════════════════════════════════════════════════════════
    def refresh_solver_binaries(self) -> None:
        """Pick up binaries (Solver Builder registry, source tree, PATH).

        Fast (no WSL call); :meth:`_detect_binaries_background` adds the
        WSL build.  A path the user typed or browsed to is never replaced.
        """
        from pycsamt.models.solver_build import discover_binaries

        for key in ENGINE_ORDER:
            bkey = ENGINES[key].binary_key
            if bkey is None:
                continue
            found = discover_binaries(bkey, probe_wsl=False)
            known = self._binaries.get(bkey, [])
            merged = list(found) + [c for c in known
                                    if all(c.path != f.path for f in found)]
            self._binaries[bkey] = merged
        self._update_engine_pills()
        self._fill_binary_combo()

    def _detect_binaries_background(self) -> None:
        if self._detector is not None and self._detector.isRunning():
            return
        from pycsamt.app.desktop.workers.inversion_worker import (
            EngineRunWorker,
        )
        from pycsamt.models.solver_build import discover_binaries

        keys = sorted({e.binary_key for e in ENGINES.values()
                       if e.binary_key})

        def task(_rep):
            return {k: discover_binaries(k, probe_wsl=True) for k in keys}

        self._btn_detect.setEnabled(False)
        self._binary_origin.setText("Searching for solver binaries…")
        self._detector = EngineRunWorker(task, parent=self)
        self._detector.finished.connect(self._on_binaries_detected)
        self._detector.error.connect(lambda _m: self._btn_detect.setEnabled(
            True))
        self._detector.start()

    def _on_binaries_detected(self, found: dict) -> None:
        self._btn_detect.setEnabled(True)
        for k, cands in found.items():
            self._binaries[k] = list(cands)
        self._update_engine_pills()
        self._fill_binary_combo()

    def _update_engine_pills(self) -> None:
        for key, row in self._engine_rows.items():
            bkey = ENGINES[key].binary_key
            if bkey is None:
                row.set_status("builtin", "Pure Python — no binary needed")
                continue
            cands = self._binaries.get(bkey, [])
            if not cands:
                row.set_status("missing", "No executable found — use "
                                          "Build… or browse to one")
            elif cands[0].path.startswith("wsl:"):
                row.set_status("wsl", f"{cands[0].path} (runs in WSL)")
            else:
                row.set_status("ready", f"{cands[0].path} ({cands[0].origin})")

    def _fill_binary_combo(self) -> None:
        """Candidates for the current solver; a path the user chose for
        *this* solver is kept (each solver remembers its own)."""
        bkey = self._engine().binary_key
        combo = self._binary_combo
        if bkey is None:
            return
        user = self._user_binaries.get(bkey, "")
        combo.blockSignals(True)
        combo.clear()
        for c in self._binaries.get(bkey, []):
            combo.addItem(c.path, c.origin)
        if user:
            i = combo.findText(user)
            if i < 0:
                combo.insertItem(0, user, "your path")
                i = 0
            combo.setCurrentIndex(i)
        elif combo.count():
            combo.setCurrentIndex(0)
        else:
            combo.setEditText("")
        combo.blockSignals(False)
        self._on_binary_pick()

    def _on_binary_edited(self, text: str) -> None:
        """Typed or browsed path: remember it for this solver only."""
        bkey = self._engine().binary_key
        if bkey is None:
            return
        text = text.strip()
        cands = self._binaries.get(bkey, [])
        best = cands[0].path if cands else ""
        if text and text != best:
            self._user_binaries[bkey] = text  # typed, browsed or 2nd pick
        else:
            self._user_binaries.pop(bkey, None)
        self._update_run_summary()

    def _on_binary_pick(self, *_) -> None:
        origin = self._binary_combo.currentData()
        path = self._binary_combo.currentText().strip()
        wsl = path.startswith("wsl:")
        if not path:
            self._binary_origin.setText("No executable found — click "
                                        "“Build…” or browse to one.")
        else:
            self._binary_origin.setText(
                f"Source: {origin or 'your path'}"
                + (" · runs inside WSL (Linux)" if wsl else ""))
        self._chk_stage.setEnabled(bool(path) and not wsl)
        if wsl:
            self._chk_stage.setChecked(False)
            self._chk_stage.setToolTip("A WSL binary runs inside Linux and "
                                       "cannot be copied into a Windows "
                                       "folder.")
        self._update_run_summary()

    def _current_binary(self) -> str | None:
        text = self._binary_combo.currentText().strip()
        return text or None

    def _browse_binary(self) -> None:
        p, _ = QFileDialog.getOpenFileName(self, "Select solver executable",
                                           str(Path.home()))
        if p:
            bkey = self._engine().binary_key
            if bkey:
                self._user_binaries[bkey] = p
            self._fill_binary_combo()

    # ══ build ════════════════════════════════════════════════════════════
    def _current_workdir(self, create: bool = True) -> Path:
        text = self._workdir_edit.text().strip()
        if text:
            p = Path(text).expanduser()
        else:
            p = (Path.home() / "pycsamt-inversions" /
                 f"{self._engine_key}-{time.strftime('%Y%m%d-%H%M%S')}")
            if create:
                self._workdir_edit.setText(str(p))
        if create:
            p.mkdir(parents=True, exist_ok=True)
        return p

    def _busy(self) -> bool:
        return self._worker is not None and self._worker.isRunning()

    def _on_build(self) -> None:
        self._start_build(then_run=False)

    def _start_build(self, *, then_run: bool) -> bool:
        if self._busy():
            return False
        eng = self._engine()
        sites = self._selected_sites()
        reason = eng.check_sites(sites)
        if reason:
            self._go_step(0)
            self._monitor.finish("error", reason)
            return False
        from pycsamt.app.desktop.workers.inversion_worker import (
            EngineRunWorker,
        )

        workdir = self._current_workdir()
        values = self.values()
        self._run_after_build = then_run
        self._console.append(f"── Build {eng.label} inputs → {workdir} ──")
        self._monitor.set_busy(f"Building {eng.label} inputs…")
        self._worker_kind = "build"
        self._worker = EngineRunWorker(
            lambda _rep: eng.build(sites, workdir, values), parent=self)
        self._worker.finished.connect(self._on_build_done)
        self._worker.error.connect(self._on_build_error)
        self._set_running_ui(True)
        self._worker.start()
        return True

    def _on_build_done(self, info: BuildInfo) -> None:
        self._worker = None
        self._set_running_ui(False)
        self._build = info
        self._build_stale = False
        self._build_state_lbl.setText(f"Inputs built · {info.workdir.name}")
        while self._summary_form.rowCount():
            self._summary_form.removeRow(0)
        for k, v in info.summary:
            lbl = QLabel(v)
            lbl.setWordWrap(True)
            lbl.setTextInteractionFlags(
                Qt.TextInteractionFlag.TextSelectableByMouse)
            self._summary_form.addRow(f"{k}:", lbl)
        self._build_warn.setText("\n".join(f"⚠ {w}" for w in info.warnings))
        for k, v in info.summary:
            self._console.append(f"  {k}: {v}")
        for w in info.warnings:
            self._console.append(f"  warning: {w}")
        draw_into(self._mesh_view,
                  lambda fig: self._engine().plot_mesh(fig, info))
        self._monitor.finish("ready", "Inputs written — ready to run")
        self._update_run_summary()
        if self._run_after_build:
            self._run_after_build = False
            self._launch()
        elif self._step_stack.currentIndex() == 0:
            self._go_step(1)

    def _on_build_error(self, msg: str) -> None:
        self._worker = None
        self._run_after_build = False
        self._set_running_ui(False)
        self._console.append(f"ERROR building inputs: {msg}")
        self._monitor.finish("error", "Build failed — " + msg.splitlines()[0]
                             [:140])
        self._mesh_view.show_unavailable("Build failed", msg)

    # ══ run ══════════════════════════════════════════════════════════════
    def _on_run(self) -> None:
        if self._mode_stack.currentIndex() == 1:
            self._ai_page.run()
            return
        if self._busy():
            return
        eng = self._engine()
        if eng.binary_key and not self._current_binary():
            self._go_step(3)
            self._monitor.finish("error", f"No {eng.label} executable — "
                                 "click “Build…” or browse to one.")
            return
        if self._build is None or self._build_stale:
            self._start_build(then_run=True)
        else:
            self._launch()

    def _launch(self) -> None:
        from pycsamt.app.desktop.workers.inversion_worker import (
            EngineRunWorker,
        )
        from pycsamt.models.solver_build import stage_binary

        eng = self._engine()
        build = self._build
        values = dict(build.values)
        binary = self._current_binary() if eng.binary_key else None
        if binary and self._chk_stage.isChecked() and \
                not binary.startswith("wsl:"):
            try:
                binary = stage_binary(binary, build.workdir,
                                      move=self._rb_move.isChecked())
                self._console.append(f"Binary {'moved' if self._rb_move.isChecked() else 'copied'} "
                                     f"to {binary}")
                if self._rb_move.isChecked():
                    self.refresh_solver_binaries()
            except Exception as exc:
                self._console.append(f"warning: could not stage the binary "
                                     f"({exc}); using it in place")
        max_it = int(values.get("max_iterations", 0) or 0)
        target = values.get("target_misfit", values.get("target_rms"))
        n_units = len(build.stations) if eng.key == "occam1d" else 1
        self._console.append(f"── Run {eng.label} · {build.workdir} ──")
        if binary:
            self._console.append(f"$ binary: {binary}")
        self._run_engine = eng.key
        self._run_workdir = build.workdir
        self._worker_kind = "run"
        self._worker = EngineRunWorker(eng.make_task(build, values, binary),
                                       parent=self)
        w = self._worker
        w.line.connect(self._console.append)
        w.stage.connect(self._monitor.set_stage)
        w.iteration.connect(self._on_iteration)
        w.finished.connect(self._on_run_done)
        w.error.connect(self._on_run_error)
        w.cancelled.connect(self._on_run_cancelled)
        self._monitor.start(f"Running {eng.label}", max_iter=max_it,
                            target=float(target) if target else None,
                            n_units=n_units)
        self._set_running_ui(True)
        if eng.binary_key:
            self._poll.start()
        self.run_state.emit(f"Inversion: {eng.label}", 0, True)
        w.start()

    def _on_iteration(self, label: str, n: int, rms: float) -> None:
        self._monitor.add_iteration(label, n, rms)
        self.run_state.emit(f"Inversion: {self._engine().label}",
                            self._monitor.percent(), True)

    def _poll_history(self) -> None:
        if self._run_workdir is None or not self._run_engine:
            return
        try:
            pairs = ENGINES[self._run_engine].history(self._run_workdir)
        except Exception:
            return
        if pairs:
            self._monitor.set_history(pairs)
            self.run_state.emit(f"Inversion: "
                                f"{ENGINES[self._run_engine].label}",
                                self._monitor.percent(), True)

    def _on_stop(self) -> None:
        if self._mode_stack.currentIndex() == 1 and self._ai_page.running:
            self._ai_page.stop()
            return
        if self._busy():
            self._console.append("Stopping…")
            self._worker.cancel()

    def _end_run(self) -> None:
        self._poll.stop()
        self._poll_history()
        self._worker = None
        self._set_running_ui(False)

    def _on_run_done(self, result) -> None:
        key, wd = self._run_engine, self._run_workdir
        self._end_run()
        rms = self._monitor.history[-1][2] if self._monitor.history else None
        self._monitor.finish("done", f"{ENGINES[key].label} finished"
                             + (f" — RMS {rms:.3f}" if rms else ""))
        self._run_status.setText("Done.")
        self._console.append("── Run finished ──")
        self.run_state.emit("", 100, False)
        try:
            runs_store.record_run(key, wd, stations=len(
                self._build.stations if self._build else []),
                final_rms=rms, status="done")
        except Exception:
            pass
        self._refresh_recent()
        self.open_run(wd, key)
        self._show_classical_page(1)
        self.result_ready.emit({"engine": key, "result": result,
                               "path": str(wd)})

    def _on_run_error(self, msg: str) -> None:
        key, wd = self._run_engine, self._run_workdir
        self._end_run()
        self._console.append(f"ERROR: {msg}")
        self._monitor.finish("error", msg.splitlines()[0][:160])
        self._run_status.setText("Error.")
        self.run_state.emit("", 0, False)
        try:
            runs_store.record_run(key, wd, status="error")
            self._refresh_recent()
        except Exception:
            pass

    def _on_run_cancelled(self) -> None:
        self._end_run()
        self._console.append("── Run stopped by the user ──")
        self._monitor.finish("stopped", "Stopped — files written so far are "
                             "kept in the run folder")
        self._run_status.setText("Stopped.")
        self.run_state.emit("", 0, False)

    def _set_running_ui(self, running: bool) -> None:
        self._btn_run.setEnabled(not running)
        self._btn_build.setEnabled(not running)
        self._engine_list.setEnabled(not running)
        for forms in self._forms.values():
            for form in forms.values():
                form.set_enabled(not running)
        self._btn_stop.setEnabled(running and self._worker_kind == "run")

    # ── run summary ───────────────────────────────────────────────────
    def _update_run_summary(self, *_) -> None:
        if not hasattr(self, "_run_summary"):
            return
        eng = self._engine()
        n = len(self._ticked_names())
        v = self.values()
        rows = [("Engine", f"{eng.label} ({eng.dim})"),
                ("Stations", str(n)),
                ("Band", f"{v['freq_min']:g} – {v['freq_max']:g} Hz"
                 if "freq_min" in v else "all frequencies"),
                ("Iterations", str(v.get("max_iterations", "–"))),
                ("Target RMS", str(v.get("target_misfit",
                                         v.get("target_rms", "–")))),
                ("Inputs", "built" if self._build and not self._build_stale
                 else "will be built on Run")]
        if eng.binary_key:
            b = self._current_binary()
            rows.append(("Binary", b or "<span style='color:#c92a2a'>"
                                       "missing</span>"))
        self._run_summary.setText("<table>" + "".join(
            f"<tr><td style='padding-right:12px'><b>{k}</b></td>"
            f"<td>{val}</td></tr>" for k, val in rows) + "</table>")

    # ══ results ══════════════════════════════════════════════════════════
    def _browse_run_folder(self) -> None:
        start = str(self._loaded.path) if self._loaded else str(
            Path.home() / "pycsamt-inversions")
        d = QFileDialog.getExistingDirectory(self, "Open inversion run "
                                             "folder", start)
        if d:
            self.open_run(d)

    def open_run(self, path, engine_key: str | None = None, *,
                 iteration=None, station=None) -> bool:
        """Open a run folder in View Results; returns success."""
        path = Path(path)
        key = engine_key if engine_key in ENGINES else detect_engine(path)
        self._show_classical_page(1)
        self._mode_group.button(0).setChecked(True)
        self._on_mode(0)
        if key is None:
            self._result_view.show_unavailable(
                "Not an inversion run folder",
                f"No Occam1D, Occam2D, ModEM or MARE2DEM run was recognised "
                f"in {path}.")
            self._res_path.setText(str(path))
            return False
        try:
            run = ENGINES[key].load(path, iteration=iteration,
                                    station=station)
        except Exception as exc:
            self._result_view.show_unavailable(
                "Could not load the run", f"{type(exc).__name__}: {exc}")
            self._res_stack.setCurrentWidget(self._result_view)
            return False
        self._set_loaded(run)
        return True

    def _set_loaded(self, run: LoadedRun) -> None:
        self._loaded = run
        eng = ENGINES[run.engine]
        self._res_path.setText(f"{eng.label} · {run.path}")
        self._res_path.setToolTip(str(run.path))
        while self._res_summary.rowCount():
            self._res_summary.removeRow(0)
        self._res_summary.addRow("Engine:", QLabel(eng.label))
        for k, v in run.summary:
            lbl = QLabel(str(v))
            lbl.setWordWrap(True)
            self._res_summary.addRow(f"{k}:", lbl)
        for combo in (self._res_station, self._res_iter):
            combo.blockSignals(True)
            combo.clear()
        for s in run.stations:
            self._res_station.addItem(s)
        self._res_station.setCurrentText(run.station or "")
        self._res_station.setEnabled(len(run.stations) > 1)
        selectable = run.engine in ("occam1d", "occam2d")
        for it in run.iterations:
            self._res_iter.addItem(str(it), it)
        if run.iteration is not None:
            self._res_iter.setCurrentIndex(max(
                self._res_iter.findData(run.iteration), 0))
        self._res_iter.setEnabled(selectable and len(run.iterations) > 1)
        self._res_iter.setToolTip(
            "Iteration to show" if selectable else
            f"{eng.label} results show the final iteration")
        for combo in (self._res_station, self._res_iter):
            combo.blockSignals(False)
        self._btn_export.setEnabled(run.engine != "occam1d")
        self._btn_export.setToolTip(
            "Save the model as PCSF" if run.engine != "occam1d" else
            "PCSF export covers 2-D/3-D models; Occam1D results are saved "
            "as text in each station's model-text folder.")
        self._btn_send.setEnabled(hasattr(run.result, "rho_2d"))
        current = (self._view_list.currentItem().data(_KEY_ROLE)
                   if self._view_list.currentItem() else None)
        self._view_list.blockSignals(True)
        self._view_list.clear()
        views = eng.views(run)
        for key, label in views:
            item = QListWidgetItem(label)
            item.setData(_KEY_ROLE, key)
            self._view_list.addItem(item)
        self._view_list.blockSignals(False)
        keys = [k for k, _ in views]
        pick = keys.index(current) if current in keys else 0
        if keys:
            self._view_list.setCurrentRow(pick)
            self._render_view(keys[pick])
        if run.engine == "occam1d":
            self._ai_page.set_classical_run(run)

    def _reload_result(self, *_) -> None:
        run = self._loaded
        if run is None:
            return
        station = self._res_station.currentText() or None
        iteration = self._res_iter.currentData()
        if station != run.station:
            iteration = None  # each station has its own history
        try:
            new = ENGINES[run.engine].load(run.path, iteration=iteration,
                                           station=station)
        except Exception as exc:
            self._result_view.show_unavailable(
                "Could not load", f"{type(exc).__name__}: {exc}")
            return
        self._set_loaded(new)

    def _render_view(self, key: str) -> None:
        run = self._loaded
        if run is None or not key:
            return
        if key == "section" and run.engine == "occam2d":
            self._res_stack.setCurrentWidget(self._tab_section)
            try:
                self._tab_section.set_result(run.result)
            except Exception as exc:
                self._console.append(f"Section plot error: {exc}")
            return
        self._res_stack.setCurrentWidget(self._result_view)
        eng = ENGINES[run.engine]
        draw_into(self._result_view, lambda fig: eng.render(key, fig, run))

    def _export_pcsf(self) -> None:
        run = self._loaded
        if run is None:
            return
        dst, _ = QFileDialog.getSaveFileName(
            self, "Export PCSF", str(run.path / f"{run.path.name}.pcsf"),
            "PCSF (*.pcsf)")
        if dst:
            self.export_pcsf(dst)

    def export_pcsf(self, dst) -> Path | None:
        from pycsamt.format import convert_engine as ce

        run = self._loaded
        try:
            sk = ce.detect(run.path, None)
            model = ce.build_model(sk, ce.ConvertOptions(
                created_by="pycsamt desktop — Inversion Studio"))
            ce.write_model(model, Path(dst), "pcsf", False)
        except Exception as exc:
            self._console.append(f"PCSF export failed: {exc}")
            QMessageBox.warning(self, "PCSF export failed", str(exc))
            return None
        self._console.append(f"Exported {dst}")
        return Path(dst)

    def _send_to_interpretation(self) -> None:
        if self._loaded is not None:
            self.result_ready.emit({"engine": self._loaded.engine,
                                    "result": self._loaded.result,
                                    "path": str(self._loaded.path)})
            self._console.append("Model sent to the Interpretation Studio.")

    # ══ library ══════════════════════════════════════════════════════════
    def _populate_presets(self) -> None:
        items = runs_store.presets(self._engine_key)
        self._preset_list.clear()
        self._preset_combo.clear()
        for p in items:
            label = p.name + ("" if p.builtin else "  (yours)")
            it = QListWidgetItem(label)
            it.setData(_KEY_ROLE, p.name)
            it.setToolTip(p.description or ", ".join(
                f"{k}={v}" for k, v in p.values.items()))
            self._preset_list.addItem(it)
            self._preset_combo.addItem(p.name)

    def _apply_preset(self, name: str) -> None:
        for p in runs_store.presets(self._engine_key):
            if p.name == name:
                forms = self._engine_forms(self._engine_key)
                for form in forms.values():
                    form.set_values({k: v for k, v in p.values.items()
                                     if k in form.widgets})
                self._console.append(f"Preset “{name}” applied to "
                                     f"{self._engine().label}.")
                return

    def _save_preset(self) -> None:
        name, ok = QInputDialog.getText(self, "Save preset",
                                        "Preset name:")
        if ok and name.strip():
            vals = {k: v for k, v in self.values().items()
                    if not k.startswith("freq_")}
            runs_store.save_preset(self._engine_key, name.strip(), vals)
            self._populate_presets()

    def _delete_preset(self) -> None:
        it = self._preset_list.currentItem()
        if it is None:
            return
        name = it.data(_KEY_ROLE)
        if any(p.name == name and p.builtin
               for p in runs_store.presets(self._engine_key)):
            return
        runs_store.delete_preset(self._engine_key, name)
        self._populate_presets()

    def _reset_engine_settings(self) -> None:
        for form in self._engine_forms(self._engine_key).values():
            form.reset()

    def _refresh_recent(self) -> None:
        self._recent_list.clear()
        for r in runs_store.recent_runs():
            eng = ENGINES.get(r.get("engine"))
            rms = r.get("final_rms")
            text = (f"{r.get('label')}\n{eng.label if eng else r['engine']}"
                    f" · {r.get('when', '')}"
                    + (f" · RMS {rms:.2f}" if rms else "")
                    + ("" if r.get("status") == "done"
                       else f" · {r.get('status')}"))
            it = QListWidgetItem(text)
            it.setData(_KEY_ROLE, r)
            it.setToolTip(r["workdir"])
            self._recent_list.addItem(it)

    # ══ console / misc ═══════════════════════════════════════════════════
    def _log(self, msg: str) -> None:
        self._console.append(msg)

    def _pop_out_console(self) -> None:
        if self._console.parent() is self._console_win:
            self._console_win.raise_()
            return
        self._console_win.hold(self._console)
        self._console.btn_pop.setText("Dock")
        self._console.btn_pop.setToolTip("Dock the console back")
        self._console.pop_out_requested.disconnect()
        self._console.pop_out_requested.connect(self._dock_console)
        self._console_win.show()

    def _dock_console(self) -> None:
        if self._console.parent() is not self._console_win:
            return
        self._vsplit.addWidget(self._console)
        self._console.btn_pop.setText("Pop out")
        self._console.btn_pop.setToolTip("Open the console in its own window")
        self._console.pop_out_requested.disconnect()
        self._console.pop_out_requested.connect(self._pop_out_console)
        self._console_win.hide()
        self._vsplit.setSizes([560, 180])

    @property
    def console_popped_out(self) -> bool:
        return self._console.parent() is self._console_win

    def _browse_workdir(self) -> None:
        d = QFileDialog.getExistingDirectory(
            self, "Select run folder",
            self._workdir_edit.text() or str(Path.home()))
        if d:
            self._workdir_edit.setText(d)
            self._mark_stale()

    @staticmethod
    def _reveal(path) -> None:
        if path and Path(path).exists():
            if sys.platform == "win32":
                os.startfile(str(path))  # noqa: S606
            else:
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))

    def set_dark_mode(self, dark: bool) -> None:
        self._dark = dark
        refresh_icons(self, dark)
        from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

        for canvas in self.findChildren(MplCanvas):
            canvas.apply_theme(dark)
        try:
            self._tab_section.set_dark_mode(dark)
        except Exception:
            pass

    # ── session ───────────────────────────────────────────────────────
    def save_geometry_to(self, store: dict) -> None:
        store[self._session_key] = {
            "geometry": self.saveGeometry().toBase64().data().decode(),
            "visible": self.isVisible(),
            "library_visible": self._btn_library.isChecked(),
            "console_theme": self._console.view.theme,
            "engine": self._engine_key,
            "mode": self._mode_group.checkedId(),
        }

    def restore_geometry_from(self, store: dict) -> None:
        entry = store.get(self._session_key)
        if not entry:
            return
        geo = entry.get("geometry")
        if geo:
            try:
                self.restoreGeometry(QByteArray.fromBase64(geo.encode()))
            except Exception:
                pass
        if "library_visible" in entry:
            self._btn_library.setChecked(bool(entry["library_visible"]))
        if entry.get("console_theme"):
            self._console.set_theme(entry["console_theme"])
        if entry.get("engine") in ENGINES:
            self._select_engine(entry["engine"])

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        if not self._probed:
            self._probed = True
            QTimer.singleShot(0, self._detect_binaries_background)

    def closeEvent(self, event) -> None:  # noqa: N802
        """Hide instead of destroying: a run keeps going."""
        self.hide()
        event.ignore()
        self.panel_closed.emit()


__all__ = ["InversionWindow"]
