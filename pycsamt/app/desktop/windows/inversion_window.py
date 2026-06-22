# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
InversionWindow — full 1-D / 2-D / 3-D MT inversion panel.

Three-panel layout (left params | right content)
────────────────────────────────────────────────
Left (290 px)
  • Dimension      : 1D / 2D / 3D radio buttons
  • Block A        : Traditional solver (Occam2D · ModEM · MARE2DEM)
                     config QStackedWidget per solver
  • Block B        : AI solver (EMInverter1D · EMInverter2D · GCNInverter3D)
                     config QStackedWidget per solver
  • Starting model : initial ρ + "← From Forward" load button
  • Data           : station list, TE/TM modes, freq range
  • Output         : working directory + Browse
  • [▶ Run]  [■ Stop]

Right (tabs)
  • Log         : live QPlainTextEdit + QProgressBar
  • Model       : resistivity result canvas
  • Fit         : observed vs. predicted response canvas
  • Convergence : RMS / loss vs. iteration canvas
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
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
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QStackedWidget,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base import PanelWindow, make_group, icon_button

logger = logging.getLogger(__name__)

# ── Dimension availability maps ───────────────────────────────────────────────
# Which solvers are valid for each dimension
_TRAD_FOR_DIM = {
    "1D": [],                           # no traditional 1-D solver
    "2D": ["occam2d", "modem", "mare2dem"],
    "3D": ["modem"],
}
_AI_FOR_DIM = {
    "1D": ["inv1d"],
    "2D": ["inv2d"],
    "3D": ["inv3d"],
}


def _spin(value, lo, hi, decimals=0, step=1.0) -> QDoubleSpinBox:
    sb = QDoubleSpinBox()
    sb.setRange(lo, hi)
    sb.setDecimals(decimals)
    sb.setSingleStep(step)
    sb.setValue(value)
    sb.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    return sb


def _ispin(value, lo, hi) -> QSpinBox:
    sb = QSpinBox()
    sb.setRange(lo, hi)
    sb.setValue(value)
    sb.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    return sb


# ═════════════════════════════════════════════════════════════════════════════
# InversionWindow
# ═════════════════════════════════════════════════════════════════════════════

class InversionWindow(PanelWindow):
    """
    Full 1-D / 2-D / 3-D MT inversion panel window.

    Replaces the old ``InversionWizardDialog``.  Stays open as a floating
    window so the user can monitor long runs while using other panels.
    """

    # Emitted when a result is ready (payload for external consumers)
    result_ready = Signal(dict)

    def __init__(self, parent: QWidget | None = None) -> None:
        self._sites          = None
        self._worker         = None
        self._starting_model: dict | None = None

        super().__init__(
            title="Inversion",
            session_key="inversion_window",
            params_width=295,
            icon_name="inversion",
            parent=parent,
        )
        self.resize(1260, 820)

        # Wire dim + engine selectors (widgets already built by super)
        self._dim_group.buttonClicked.connect(self._on_dim_changed)
        self._engine_group.buttonClicked.connect(self._on_engine_changed)
        self._on_dim_changed()

    # =========================================================================
    # PanelWindow overrides
    # =========================================================================

    def _build_params(self, layout: QVBoxLayout) -> None:
        self._build_dim_selector(layout)
        self._build_engine_selector(layout)
        self._build_starting_model_group(layout)
        self._build_data_group(layout)
        self._build_output_group(layout)
        self._build_run_buttons(layout)

    def _build_content(self, layout: QVBoxLayout) -> None:
        self._tabs = QTabWidget()
        self._tabs.setObjectName("InversionTabWidget")

        self._tab_log  = self._build_log_tab()
        self._tab_model       = MplCanvas(toolbar=True)
        self._tab_fit         = MplCanvas(toolbar=True)
        self._tab_convergence = MplCanvas(toolbar=True)

        self._tabs.addTab(self._tab_log,         "Log")
        self._tabs.addTab(self._tab_model,       "Model")
        self._tabs.addTab(self._tab_fit,         "Fit")
        self._tabs.addTab(self._tab_convergence, "Convergence")

        layout.addWidget(self._tabs)

    # =========================================================================
    # Left-panel builders
    # =========================================================================

    # ── Dimension ──────────────────────────────────────────────────────────

    def _build_dim_selector(self, layout: QVBoxLayout) -> None:
        grp, lay = make_group("Dimension")
        row = QHBoxLayout()
        self._dim_group = QButtonGroup(self)
        for label in ("1D", "2D", "3D"):
            rb = QRadioButton(label)
            rb.setProperty("dim", label)
            if label == "2D":
                rb.setChecked(True)
            self._dim_group.addButton(rb)
            row.addWidget(rb)
        lay.addLayout(row)
        layout.addWidget(grp)

    # ── Engine selectors (Traditional + AI) ────────────────────────────────

    def _build_engine_selector(self, layout: QVBoxLayout) -> None:
        self._engine_group = QButtonGroup(self)    # global: only one engine active

        # ── Block A: Traditional ──────────────────────────────────────
        grp_a, lay_a = make_group("Traditional Solvers")

        self._rb_occam2d  = QRadioButton("Occam2D  (2D)")
        self._rb_modem    = QRadioButton("ModEM    (2D / 3D)")
        self._rb_mare2dem = QRadioButton("MARE2DEM (2D)")
        self._rb_occam2d.setProperty("engine", "occam2d")
        self._rb_modem.setProperty("engine",   "modem")
        self._rb_mare2dem.setProperty("engine","mare2dem")
        self._rb_modem.setChecked(True)

        for rb in (self._rb_occam2d, self._rb_modem, self._rb_mare2dem):
            self._engine_group.addButton(rb)
            lay_a.addWidget(rb)

        # Stacked config per traditional solver
        self._trad_stack = QStackedWidget()
        self._trad_stack.addWidget(self._build_cfg_occam2d())    # 0
        self._trad_stack.addWidget(self._build_cfg_modem())      # 1
        self._trad_stack.addWidget(self._build_cfg_mare2dem())   # 2
        lay_a.addWidget(self._trad_stack)
        layout.addWidget(grp_a)

        # ── Block B: AI ────────────────────────────────────────────────
        grp_b, lay_b = make_group("AI Inversion")

        self._rb_inv1d = QRadioButton("EMInverter1D  ResNet/CNN/FCN")
        self._rb_inv2d = QRadioButton("EMInverter2D  U-Net")
        self._rb_inv3d = QRadioButton("GCNInverter3D  Graph NN")
        self._rb_inv1d.setProperty("engine", "inv1d")
        self._rb_inv2d.setProperty("engine", "inv2d")
        self._rb_inv3d.setProperty("engine", "inv3d")

        for rb in (self._rb_inv1d, self._rb_inv2d, self._rb_inv3d):
            self._engine_group.addButton(rb)
            lay_b.addWidget(rb)

        self._ai_stack = QStackedWidget()
        self._ai_stack.addWidget(self._build_cfg_inv1d())   # 0
        self._ai_stack.addWidget(self._build_cfg_inv2d())   # 1
        self._ai_stack.addWidget(self._build_cfg_inv3d())   # 2
        lay_b.addWidget(self._ai_stack)
        layout.addWidget(grp_b)

        # Wire engine → stacked pages
        self._rb_occam2d.toggled.connect( lambda c: self._trad_stack.setCurrentIndex(0) if c else None)
        self._rb_modem.toggled.connect(   lambda c: self._trad_stack.setCurrentIndex(1) if c else None)
        self._rb_mare2dem.toggled.connect(lambda c: self._trad_stack.setCurrentIndex(2) if c else None)
        self._rb_inv1d.toggled.connect(   lambda c: self._ai_stack.setCurrentIndex(0)   if c else None)
        self._rb_inv2d.toggled.connect(   lambda c: self._ai_stack.setCurrentIndex(1)   if c else None)
        self._rb_inv3d.toggled.connect(   lambda c: self._ai_stack.setCurrentIndex(2)   if c else None)

    # ── Traditional config panels ───────────────────────────────────────────

    def _build_cfg_occam2d(self) -> QWidget:
        w = QWidget(); f = QFormLayout(w); f.setSpacing(3)
        self._occ_n_layers    = _ispin(30, 5, 100)
        self._occ_n_air       = _ispin(5,  1, 20)
        self._occ_cell_h      = _spin(100.0, 10.0, 10000.0, 0, 10.0)
        self._occ_depth_scale = _spin(1.2,   1.01, 2.0,     2, 0.01)
        self._occ_max_iter    = _ispin(100, 1, 500)
        self._occ_target_rms  = _spin(1.0,   0.1, 50.0, 2, 0.1)
        self._occ_init_rho    = _spin(100.0, 1.0, 1e5,  1, 10.0)
        self._occ_err_rho     = _spin(0.05,  0.001, 1.0, 3, 0.01)
        self._occ_err_phase   = _spin(0.5,   0.001, 10.0, 2, 0.1)
        self._occ_binary      = QLineEdit()
        self._occ_binary.setPlaceholderText("occam2d binary path (optional)")
        btn_bin = QPushButton("…"); btn_bin.setFixedWidth(26)
        btn_bin.clicked.connect(lambda: self._browse_binary(self._occ_binary))
        bin_row = QHBoxLayout(); bin_row.addWidget(self._occ_binary); bin_row.addWidget(btn_bin)
        bin_w = QWidget(); bin_w.setLayout(bin_row)
        f.addRow("Layers:",         self._occ_n_layers)
        f.addRow("Air layers:",     self._occ_n_air)
        f.addRow("Cell size h (m):",self._occ_cell_h)
        f.addRow("Depth scale:",    self._occ_depth_scale)
        f.addRow("Max iter:",       self._occ_max_iter)
        f.addRow("Target RMS:",     self._occ_target_rms)
        f.addRow("Init ρ (Ω·m):",  self._occ_init_rho)
        f.addRow("Err floor ρ:",    self._occ_err_rho)
        f.addRow("Err floor φ (°):",self._occ_err_phase)
        f.addRow("Binary:",         bin_w)
        return w

    def _build_cfg_modem(self) -> QWidget:
        w = QWidget(); f = QFormLayout(w); f.setSpacing(3)
        self._mod_mode      = QComboBox()
        self._mod_mode.addItems(["2D", "3D"])
        self._mod_nx        = _ispin(30, 5, 200)
        self._mod_nz        = _ispin(30, 5, 200)
        self._mod_n_air     = _ispin(5,  1, 20)
        self._mod_cell_h    = _spin(500.0, 10.0, 50000.0, 0, 50.0)
        self._mod_max_iter  = _ispin(100, 1, 500)
        self._mod_target    = _spin(1.0, 0.1, 50.0, 2, 0.1)
        self._mod_init_rho  = _spin(100.0, 1.0, 1e5, 1, 10.0)
        self._mod_smooth_x  = _spin(0.3, 0.01, 1.0, 2, 0.05)
        self._mod_smooth_z  = _spin(0.3, 0.01, 1.0, 2, 0.05)
        self._mod_binary    = QLineEdit()
        self._mod_binary.setPlaceholderText("ModEM binary path (optional)")
        btn_bin = QPushButton("…"); btn_bin.setFixedWidth(26)
        btn_bin.clicked.connect(lambda: self._browse_binary(self._mod_binary))
        bin_row = QHBoxLayout(); bin_row.addWidget(self._mod_binary); bin_row.addWidget(btn_bin)
        bin_w = QWidget(); bin_w.setLayout(bin_row)
        f.addRow("Mode:",           self._mod_mode)
        f.addRow("nx (x cells):",   self._mod_nx)
        f.addRow("nz (z cells):",   self._mod_nz)
        f.addRow("Air layers:",     self._mod_n_air)
        f.addRow("Cell size h (m):",self._mod_cell_h)
        f.addRow("Max iter:",       self._mod_max_iter)
        f.addRow("Target RMS:",     self._mod_target)
        f.addRow("Init ρ (Ω·m):",  self._mod_init_rho)
        f.addRow("Smooth X:",       self._mod_smooth_x)
        f.addRow("Smooth Z:",       self._mod_smooth_z)
        f.addRow("Binary:",         bin_w)
        return w

    def _build_cfg_mare2dem(self) -> QWidget:
        w = QWidget(); f = QFormLayout(w); f.setSpacing(3)
        self._m2d_max_iter  = _ispin(30, 1, 200)
        self._m2d_target    = _spin(1.0, 0.1, 50.0, 2, 0.1)
        self._m2d_init_rho  = _spin(100.0, 1.0, 1e5, 1, 10.0)
        self._m2d_use_mpi   = QCheckBox("Use MPI")
        self._m2d_n_procs   = _ispin(4, 1, 128)
        self._m2d_binary    = QLineEdit()
        self._m2d_binary.setPlaceholderText("MARE2DEM binary path (optional)")
        btn_bin = QPushButton("…"); btn_bin.setFixedWidth(26)
        btn_bin.clicked.connect(lambda: self._browse_binary(self._m2d_binary))
        bin_row = QHBoxLayout(); bin_row.addWidget(self._m2d_binary); bin_row.addWidget(btn_bin)
        bin_w = QWidget(); bin_w.setLayout(bin_row)
        f.addRow("Max iter:",       self._m2d_max_iter)
        f.addRow("Target RMS:",     self._m2d_target)
        f.addRow("Init ρ (Ω·m):",  self._m2d_init_rho)
        f.addRow("",                self._m2d_use_mpi)
        f.addRow("MPI procs:",      self._m2d_n_procs)
        f.addRow("Binary:",         bin_w)
        return w

    # ── AI config panels ────────────────────────────────────────────────────

    def _build_cfg_inv1d(self) -> QWidget:
        w = QWidget(); f = QFormLayout(w); f.setSpacing(3)
        self._ai1_arch      = QComboBox()
        self._ai1_arch.addItems(["resnet", "cnn1d", "fcn"])
        self._ai1_solver    = QComboBox()
        self._ai1_solver.addItems(["mt1d", "csamt1d"])
        self._ai1_n_layers  = _ispin(5, 2, 20)
        self._ai1_n_samples = _ispin(2000, 100, 50000)
        self._ai1_epochs    = _ispin(60, 5, 500)
        self._ai1_batch     = _ispin(256, 8, 2048)
        self._ai1_lr        = _spin(1e-3, 1e-5, 1e-1, 5, 1e-4)
        self._ai1_noise     = _spin(0.05, 0.0, 0.5, 3, 0.01)
        self._ai1_geology   = QComboBox()
        self._ai1_geology.addItems(
            ["(none)", "sedimentary", "crystalline", "geothermal",
             "marine", "permafrost", "basement", "coastal",
             "evaporite", "hydrothermal", "laterite",
             "mineralized", "porphyry", "volcanic"]
        )
        f.addRow("Architecture:",   self._ai1_arch)
        f.addRow("Solver:",         self._ai1_solver)
        f.addRow("n_layers:",       self._ai1_n_layers)
        f.addRow("Train samples:",  self._ai1_n_samples)
        f.addRow("Epochs:",         self._ai1_epochs)
        f.addRow("Batch size:",     self._ai1_batch)
        f.addRow("Learning rate:",  self._ai1_lr)
        f.addRow("Noise level:",    self._ai1_noise)
        f.addRow("Geology prior:",  self._ai1_geology)
        return w

    def _build_cfg_inv2d(self) -> QWidget:
        w = QWidget(); f = QFormLayout(w); f.setSpacing(3)
        self._ai2_n_comp    = _ispin(4, 1, 8)
        self._ai2_n_depth   = _ispin(40, 10, 200)
        self._ai2_n_sta     = _ispin(20, 2, 100)
        self._ai2_n_freq    = _ispin(32, 8, 128)
        self._ai2_n_samples = _ispin(500, 50, 10000)
        self._ai2_epochs    = _ispin(40, 5, 300)
        self._ai2_batch     = _ispin(16, 2, 128)
        self._ai2_lr        = _spin(1e-3, 1e-5, 1e-1, 5, 1e-4)
        f.addRow("n_components:",   self._ai2_n_comp)
        f.addRow("n_depth:",        self._ai2_n_depth)
        f.addRow("n_stations:",     self._ai2_n_sta)
        f.addRow("n_freqs:",        self._ai2_n_freq)
        f.addRow("Train samples:",  self._ai2_n_samples)
        f.addRow("Epochs:",         self._ai2_epochs)
        f.addRow("Batch size:",     self._ai2_batch)
        f.addRow("Learning rate:",  self._ai2_lr)
        return w

    def _build_cfg_inv3d(self) -> QWidget:
        w = QWidget(); f = QFormLayout(w); f.setSpacing(3)
        self._ai3_n_feat    = _ispin(40, 10, 200)
        self._ai3_n_layers  = _ispin(5,  2, 20)
        self._ai3_hidden    = QLineEdit("256,128,64")
        self._ai3_hidden.setToolTip("Hidden layer sizes, comma-separated")
        self._ai3_dropout   = _spin(0.1, 0.0, 0.9, 2, 0.05)
        self._ai3_n_sta     = _ispin(16, 4, 100)
        self._ai3_n_samples = _ispin(300, 50, 5000)
        self._ai3_epochs    = _ispin(40, 5, 300)
        self._ai3_batch     = _ispin(16, 2, 128)
        self._ai3_lr        = _spin(1e-3, 1e-5, 1e-1, 5, 1e-4)
        self._ai3_radius    = _spin(5000.0, 100.0, 1e6, 0, 500.0)
        f.addRow("n_features:",     self._ai3_n_feat)
        f.addRow("n_layers:",       self._ai3_n_layers)
        f.addRow("Hidden layers:",  self._ai3_hidden)
        f.addRow("Dropout:",        self._ai3_dropout)
        f.addRow("Stations:",       self._ai3_n_sta)
        f.addRow("Train samples:",  self._ai3_n_samples)
        f.addRow("Epochs:",         self._ai3_epochs)
        f.addRow("Batch size:",     self._ai3_batch)
        f.addRow("Learning rate:",  self._ai3_lr)
        f.addRow("Graph radius (m):",self._ai3_radius)
        return w

    # ── Starting model ──────────────────────────────────────────────────────

    def _build_starting_model_group(self, layout: QVBoxLayout) -> None:
        grp, lay = make_group("Starting Model")
        f = QFormLayout()
        f.setSpacing(3)
        self._init_rho = _spin(100.0, 0.1, 1e5, 1, 10.0)
        f.addRow("Init ρ (Ω·m):", self._init_rho)
        lay.addLayout(f)

        self._fwd_model_label = QLabel("(no model from Forward)")
        self._fwd_model_label.setWordWrap(True)
        self._fwd_model_label.setObjectName("FwdModelLabel")
        lay.addWidget(self._fwd_model_label)

        btn_clear = QPushButton("✕ Clear Forward model")
        btn_clear.clicked.connect(self._clear_starting_model)
        lay.addWidget(btn_clear)
        layout.addWidget(grp)

    # ── Data ────────────────────────────────────────────────────────────────

    def _build_data_group(self, layout: QVBoxLayout) -> None:
        grp, lay = make_group("Data")

        # Station list
        self._station_list = QListWidget()
        self._station_list.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection
        )
        self._station_list.setMaximumHeight(110)
        lay.addWidget(QLabel("Stations (Ctrl-click to select):"))
        lay.addWidget(self._station_list)

        # Modes
        mode_row = QHBoxLayout()
        self._mode_te = QCheckBox("TE (xy)")
        self._mode_te.setChecked(True)
        self._mode_tm = QCheckBox("TM (yx)")
        self._mode_tm.setChecked(True)
        mode_row.addWidget(self._mode_te)
        mode_row.addWidget(self._mode_tm)
        lay.addLayout(mode_row)

        # Frequency range
        f = QFormLayout(); f.setSpacing(3)
        self._f_min = _spin(1e-3, 1e-4, 1e3, 4, 1e-4)
        self._f_max = _spin(1e3,  1e-3, 1e6, 1, 10.0)
        f.addRow("f_min (Hz):", self._f_min)
        f.addRow("f_max (Hz):", self._f_max)
        lay.addLayout(f)

        layout.addWidget(grp)

    # ── Output dir ──────────────────────────────────────────────────────────

    def _build_output_group(self, layout: QVBoxLayout) -> None:
        grp, lay = make_group("Output Directory")
        row = QHBoxLayout()
        self._workdir_edit = QLineEdit()
        self._workdir_edit.setPlaceholderText("Select working directory…")
        btn = QPushButton("Browse…")
        btn.clicked.connect(self._browse_workdir)
        row.addWidget(self._workdir_edit)
        row.addWidget(btn)
        lay.addLayout(row)
        layout.addWidget(grp)

    # ── Run / Stop buttons ──────────────────────────────────────────────────

    def _build_run_buttons(self, layout: QVBoxLayout) -> None:
        row = QHBoxLayout()
        self._btn_run  = icon_button("▶  Run",  "results", "Start inversion")
        self._btn_stop = icon_button("■  Stop", "export",  "Stop inversion")
        self._btn_run.setObjectName("RunButton")
        self._btn_stop.setEnabled(False)
        self._btn_run.clicked.connect(self._on_run)
        self._btn_stop.clicked.connect(self._on_stop)
        row.addWidget(self._btn_run)
        row.addWidget(self._btn_stop)
        layout.addLayout(row)

        self._run_status = QLabel("")
        self._run_status.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._run_status.setObjectName("RunStatusLabel")
        layout.addWidget(self._run_status)

    # =========================================================================
    # Right panel: log tab
    # =========================================================================

    def _build_log_tab(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(4, 4, 4, 4)
        v.setSpacing(4)

        self._log_edit = QPlainTextEdit()
        self._log_edit.setReadOnly(True)
        self._log_edit.setMaximumBlockCount(4000)
        self._log_edit.setObjectName("InversionLog")
        v.addWidget(self._log_edit)

        self._progress_bar = QProgressBar()
        self._progress_bar.setRange(0, 100)
        self._progress_bar.setVisible(False)
        v.addWidget(self._progress_bar)
        return w

    # =========================================================================
    # Dimension switching
    # =========================================================================

    def _current_dim(self) -> str:
        btn = self._dim_group.checkedButton()
        return btn.property("dim") if btn else "2D"

    def _on_dim_changed(self, *_) -> None:
        dim = self._current_dim()
        trad_ok = _TRAD_FOR_DIM[dim]
        ai_ok   = _AI_FOR_DIM[dim]

        for rb, key in [
            (self._rb_occam2d,  "occam2d"),
            (self._rb_modem,    "modem"),
            (self._rb_mare2dem, "mare2dem"),
        ]:
            rb.setEnabled(key in trad_ok)

        for rb, key in [
            (self._rb_inv1d, "inv1d"),
            (self._rb_inv2d, "inv2d"),
            (self._rb_inv3d, "inv3d"),
        ]:
            rb.setEnabled(key in ai_ok)

        # Auto-select a valid engine if current is now disabled
        checked = self._engine_group.checkedButton()
        if checked and not checked.isEnabled():
            # pick first available
            for rb in self._engine_group.buttons():
                if rb.isEnabled():
                    rb.setChecked(True)
                    break
        self._on_engine_changed()

    def _on_engine_changed(self, *_) -> None:
        btn = self._engine_group.checkedButton()
        if btn is None:
            return
        engine = btn.property("engine")
        # Show/hide stacked config widgets
        trad_engines = {"occam2d", "modem", "mare2dem"}
        if engine in trad_engines:
            self._trad_stack.setVisible(True)
            self._ai_stack.setVisible(False)
        else:
            self._trad_stack.setVisible(False)
            self._ai_stack.setVisible(True)

    def _current_engine(self) -> str:
        btn = self._engine_group.checkedButton()
        return btn.property("engine") if btn else "modem"

    # =========================================================================
    # Starting model (from ForwardModelWindow)
    # =========================================================================

    def load_starting_model(self, payload: dict) -> None:
        """Called by MainWindow when ForwardModelWindow.send_to_inversion fires."""
        self._starting_model = payload
        dim = payload.get("dim", "?")
        rho_list = payload.get("resistivity", [])
        label = f"Forward model: dim={dim}"
        if rho_list:
            label += f", {len(rho_list)} layers"
        self._fwd_model_label.setText(label)
        # Propagate init_rho from first layer if 1D
        if rho_list:
            self._init_rho.setValue(float(rho_list[0]))
        # Switch inversion dim to match
        for btn in self._dim_group.buttons():
            if btn.property("dim") == dim:
                btn.setChecked(True)
        self._on_dim_changed()
        self._log(f"Starting model received from Forward ({label}).")
        self._tabs.setCurrentIndex(0)  # show log tab

    def _clear_starting_model(self) -> None:
        self._starting_model = None
        self._fwd_model_label.setText("(no model from Forward)")

    # =========================================================================
    # Run / Stop
    # =========================================================================

    def _on_run(self) -> None:
        if self._worker and self._worker.isRunning():
            return
        engine = self._current_engine()
        if engine in {"occam2d", "modem", "mare2dem"}:
            self._run_traditional(engine)
        else:
            self._run_ai(engine)

    def _on_stop(self) -> None:
        if self._worker:
            try:
                self._worker.cancel()
            except AttributeError:
                self._worker.terminate()
            self._worker.quit()
        self._btn_run.setEnabled(True)
        self._btn_stop.setEnabled(False)
        self._progress_bar.setVisible(False)
        self._log("Stopped.")
        self._run_status.setText("Stopped.")

    # ── Traditional run ─────────────────────────────────────────────────────

    def _run_traditional(self, engine: str) -> None:
        from pycsamt.app.desktop.workers.inversion_worker import InversionWorker

        workdir = self._workdir_edit.text().strip()
        if not workdir:
            QMessageBox.warning(self, "No workdir", "Please set an output directory.")
            return
        Path(workdir).mkdir(parents=True, exist_ok=True)

        if engine == "occam2d":
            ok = self._launch_occam2d(workdir)
        elif engine == "modem":
            ok = self._launch_modem(workdir)
        else:
            ok = self._launch_mare2dem(workdir)

        if not ok:
            return

    def _launch_occam2d(self, workdir: str) -> bool:
        from pycsamt.models.occam2d import InputBuilder, OccamConfig
        from pycsamt.app.desktop.workers.inversion_worker import InversionWorker

        modes = []
        if self._mode_te.isChecked(): modes.append("te")
        if self._mode_tm.isChecked(): modes.append("tm")
        cfg = OccamConfig(
            modes=modes,
            n_layers=self._occ_n_layers.value(),
            n_airlayers=self._occ_n_air.value(),
            cell_size_horizontal=self._occ_cell_h.value(),
            depth_scale=self._occ_depth_scale.value(),
            max_iterations=self._occ_max_iter.value(),
            target_misfit=self._occ_target_rms.value(),
            initial_rho=self._init_rho.value(),
            error_floor_rho=self._occ_err_rho.value(),
            error_floor_phase=self._occ_err_phase.value(),
        )
        sites = self._selected_sites()
        try:
            self._log("Building Occam2D input files…")
            InputBuilder(source=sites, workdir=workdir, config=cfg).build()
            self._log("Input files written.")
        except Exception as exc:
            self._log(f"InputBuilder error: {exc}")
            return False

        binary = self._occ_binary.text().strip() or None
        worker = InversionWorker(
            workdir=workdir, binary_path=binary,
            max_iter=self._occ_max_iter.value(),
            target_misfit=self._occ_target_rms.value(),
            parent=self,
        )
        self._attach_trad_worker(worker)
        return True

    def _launch_modem(self, workdir: str) -> bool:
        from pycsamt.models.modem import ModEmConfig, InputBuilder, ModEmRunner
        from pycsamt.app.desktop.workers.inversion_worker import InversionWorker

        mode_str = self._mod_mode.currentText().lower()  # "2d" or "3d"
        cfg = ModEmConfig(
            mode=mode_str,
            nx_2d=self._mod_nx.value(),
            nz_2d=self._mod_nz.value(),
            n_airlayers_2d=self._mod_n_air.value(),
            cell_size_h_2d=self._mod_cell_h.value(),
            max_iterations=self._mod_max_iter.value(),
            target_rms=self._mod_target.value(),
            initial_rho=self._init_rho.value(),
            smooth_x=self._mod_smooth_x.value(),
            smooth_z=self._mod_smooth_z.value(),
            binary_2d=self._mod_binary.text().strip() or None,
        )
        sites = self._selected_sites()
        try:
            self._log(f"Building ModEM {mode_str.upper()} input files…")
            InputBuilder(source=sites, workdir=workdir, config=cfg).build()
            self._log("Input files written.")
        except Exception as exc:
            self._log(f"InputBuilder error: {exc}")
            return False

        binary = (
            self._mod_binary.text().strip() or
            (cfg.binary_2d if mode_str == "2d" else getattr(cfg, "binary_3d", None))
        )
        worker = InversionWorker(
            workdir=workdir, binary_path=binary,
            max_iter=self._mod_max_iter.value(),
            target_misfit=self._mod_target.value(),
            parent=self,
        )
        self._attach_trad_worker(worker)
        return True

    def _launch_mare2dem(self, workdir: str) -> bool:
        from pycsamt.models.mare2dem import Mare2DEMConfig, InputBuilder, Mare2DEMRunner
        from pycsamt.app.desktop.workers.inversion_worker import InversionWorker

        cfg = Mare2DEMConfig(
            max_iterations=self._m2d_max_iter.value(),
            target_rms=self._m2d_target.value(),
            initial_rho=self._init_rho.value(),
            use_mpi=self._m2d_use_mpi.isChecked(),
            n_procs=self._m2d_n_procs.value(),
            binary=self._m2d_binary.text().strip() or None,
        )
        sites = self._selected_sites()
        try:
            self._log("Building MARE2DEM input files…")
            InputBuilder(source=sites, workdir=workdir, config=cfg).build()
            self._log("Input files written.")
        except Exception as exc:
            self._log(f"InputBuilder error: {exc}")
            return False

        binary = self._m2d_binary.text().strip() or None
        worker = InversionWorker(
            workdir=workdir, binary_path=binary,
            max_iter=self._m2d_max_iter.value(),
            target_misfit=self._m2d_target.value(),
            parent=self,
        )
        self._attach_trad_worker(worker)
        return True

    def _attach_trad_worker(self, worker) -> None:
        self._worker = worker
        worker.stdout_line.connect(self._log)
        worker.progress.connect(self._progress_bar.setValue)
        worker.finished.connect(self._on_trad_finished)
        worker.error.connect(self._on_error)
        self._start_worker()

    # ── AI run ───────────────────────────────────────────────────────────────

    def _run_ai(self, engine: str) -> None:
        from pycsamt.app.desktop.workers.ai_inversion_worker import AIInversionWorker

        dim    = self._current_dim()
        params = self._build_ai_params(engine, dim)
        self._worker = AIInversionWorker(params, parent=self)
        self._worker.log_line.connect(self._log)
        self._worker.progress.connect(self._progress_bar.setValue)
        self._worker.finished.connect(self._on_ai_finished)
        self._worker.error.connect(self._on_error)
        self._start_worker()

    def _build_ai_params(self, engine: str, dim: str) -> dict:
        p: dict = {
            "dim":   dim,
            "f_min": self._f_min.value(),
            "f_max": self._f_max.value(),
        }
        if engine == "inv1d":
            geo = self._ai1_geology.currentText()
            p.update({
                "arch":        self._ai1_arch.currentText(),
                "solver":      self._ai1_solver.currentText(),
                "n_layers":    self._ai1_n_layers.value(),
                "n_samples":   self._ai1_n_samples.value(),
                "epochs":      self._ai1_epochs.value(),
                "batch_size":  self._ai1_batch.value(),
                "lr":          self._ai1_lr.value(),
                "noise_level": self._ai1_noise.value(),
                "geology":     None if geo == "(none)" else geo,
                "n_freq":      30,
            })
        elif engine == "inv2d":
            p.update({
                "n_components": self._ai2_n_comp.value(),
                "n_depth":      self._ai2_n_depth.value(),
                "n_stations":   self._ai2_n_sta.value(),
                "n_freq":       self._ai2_n_freq.value(),
                "n_samples":    self._ai2_n_samples.value(),
                "epochs":       self._ai2_epochs.value(),
                "batch_size":   self._ai2_batch.value(),
                "lr":           self._ai2_lr.value(),
            })
        else:
            hidden_txt = self._ai3_hidden.text().strip()
            try:
                hidden = [int(x.strip()) for x in hidden_txt.split(",") if x.strip()]
            except ValueError:
                hidden = [256, 128, 64]
            p.update({
                "n_features":  self._ai3_n_feat.value(),
                "n_layers":    self._ai3_n_layers.value(),
                "hidden":      hidden,
                "dropout":     self._ai3_dropout.value(),
                "n_sta":       self._ai3_n_sta.value(),
                "n_samples":   self._ai3_n_samples.value(),
                "epochs":      self._ai3_epochs.value(),
                "batch_size":  self._ai3_batch.value(),
                "lr":          self._ai3_lr.value(),
                "radius":      self._ai3_radius.value(),
                "n_freq":      20,
            })
        # Attach observed data if sites are loaded
        p["X_obs"] = self._build_X_obs(engine, dim, p)
        return p

    def _build_X_obs(self, engine: str, dim: str, params: dict):
        """Convert loaded Sites to X_obs array expected by the AI inverter."""
        if self._sites is None:
            return None
        try:
            import numpy as np
            sites_list = self._selected_sites()
            if sites_list is None:
                return None
            f_min  = max(float(params.get("f_min", 1e-3)), 1e-6)
            f_max  = max(float(params.get("f_max", 1e3)),  f_min * 10)
            n_freq = int(params.get("n_freq", 30))
            freqs  = np.logspace(np.log10(f_min), np.log10(f_max), n_freq)
            rows   = []
            for site in sites_list.as_list():
                rho_a = site.interpolate_rho_a(freqs)
                phase = site.interpolate_phase(freqs)
                rows.append(np.concatenate([
                    np.log10(np.maximum(rho_a, 1e-12)),
                    phase,
                ]))
            return np.array(rows, dtype=float) if rows else None
        except Exception:
            return None

    # ── Shared worker start ──────────────────────────────────────────────────

    def _start_worker(self) -> None:
        self._btn_run.setEnabled(False)
        self._btn_stop.setEnabled(True)
        self._progress_bar.setValue(0)
        self._progress_bar.setVisible(True)
        self._run_status.setText("Running…")
        self._tabs.setCurrentIndex(0)   # switch to Log tab
        self._worker.start()

    # =========================================================================
    # Result handlers
    # =========================================================================

    def _on_trad_finished(self, result) -> None:
        self._finish_run()
        self._log("Inversion finished.")
        self._run_status.setText("Done.")
        try:
            self._plot_trad_result(result)
        except Exception as exc:
            self._log(f"Plot error: {exc}")
        self.result_ready.emit({"engine": self._current_engine(), "result": result})

    def _on_ai_finished(self, result: dict) -> None:
        self._finish_run()
        self._log("AI inversion finished.")
        self._run_status.setText("Done.")
        try:
            self._plot_ai_result(result)
        except Exception as exc:
            self._log(f"Plot error: {exc}")
        self.result_ready.emit({"engine": self._current_engine(), "result": result})

    def _on_error(self, msg: str) -> None:
        self._finish_run()
        self._log(f"ERROR: {msg}")
        self._run_status.setText("Error.")
        QMessageBox.critical(self, "Inversion Error", msg)

    def _finish_run(self) -> None:
        self._btn_run.setEnabled(True)
        self._btn_stop.setEnabled(False)
        self._progress_bar.setVisible(False)

    # =========================================================================
    # Plotting helpers
    # =========================================================================

    def _plot_trad_result(self, result) -> None:
        """Plot Occam2D / ModEM / MARE2DEM result."""
        # Model tab
        self._tab_model.figure.clear()
        ax = self._tab_model.figure.add_subplot(111)
        try:
            result.plot_model(ax=ax)
        except AttributeError:
            try:
                from pycsamt.models.occam2d import PlotModel
                PlotModel(result).plot(ax=ax)
            except Exception as exc:
                ax.set_title(f"Model plot unavailable: {exc}")
        self._tab_model.figure.tight_layout()
        self._tab_model.draw()
        self._tabs.setCurrentIndex(1)

        # Fit tab — predicted vs observed response
        self._tab_fit.figure.clear()
        ax = self._tab_fit.figure.add_subplot(111)
        try:
            result.plot_response(ax=ax)
        except Exception as exc:
            ax.set_title(f"Response plot unavailable: {exc}")
        self._tab_fit.figure.tight_layout()
        self._tab_fit.draw()

        # Convergence tab — RMS vs iteration
        self._tab_convergence.figure.clear()
        ax = self._tab_convergence.figure.add_subplot(111)
        try:
            result.plot_misfit(ax=ax)
        except Exception as exc:
            ax.set_title(f"Convergence plot unavailable: {exc}")
        self._tab_convergence.figure.tight_layout()
        self._tab_convergence.draw()

    def _plot_ai_result(self, result: dict) -> None:
        dim    = result.get("dim", "1D")
        y_pred = result.get("y_pred")
        freqs  = result.get("freqs")

        # ── Model tab: predicted resistivity profile(s) ─────────────
        self._tab_model.figure.clear()
        ax = self._tab_model.figure.add_subplot(111)
        if y_pred is not None:
            try:
                self._plot_ai_model(ax, y_pred, dim, result)
            except Exception as exc:
                ax.set_title(f"AI model plot error: {exc}")
        self._tab_model.figure.tight_layout()
        self._tab_model.draw()
        self._tabs.setCurrentIndex(1)

        # ── Fit tab: forward re-predict vs observed ──────────────────
        self._tab_fit.figure.clear()
        ax = self._tab_fit.figure.add_subplot(111)
        try:
            self._plot_ai_fit(ax, result, freqs)
        except Exception as exc:
            ax.set_title(f"Fit plot unavailable: {exc}")
        self._tab_fit.figure.tight_layout()
        self._tab_fit.draw()

        # ── Convergence tab: training loss ───────────────────────────
        self._tab_convergence.figure.clear()
        ax = self._tab_convergence.figure.add_subplot(111)
        try:
            inv = result.get("inverter")
            if inv is not None and hasattr(inv, "loss_history_"):
                history = inv.loss_history_
                ax.plot(history, "b-", lw=1.5)
                ax.set_xlabel("Epoch")
                ax.set_ylabel("Loss")
                ax.set_title("Training loss")
                ax.set_yscale("log")
            else:
                ax.set_title("Loss history not available")
        except Exception as exc:
            ax.set_title(f"Convergence: {exc}")
        self._tab_convergence.figure.tight_layout()
        self._tab_convergence.draw()

    def _plot_ai_model(self, ax, y_pred, dim: str, result: dict) -> None:
        import numpy as np
        if dim == "1D":
            n_layers = result.get("n_layers", 5)
            # y_pred shape: (n_stations, 2*n_layers-1)
            # first n_layers cols = log10(rho), rest = thickness
            for i, row in enumerate(y_pred[:min(5, len(y_pred))]):
                rho   = 10.0 ** row[:n_layers]
                thick = row[n_layers:]
                depths = np.concatenate([[0], np.cumsum(thick), [np.sum(thick) * 1.5]])
                ax.step(np.append(rho, rho[-1]), depths, where="post",
                        label=f"sta {i+1}", alpha=0.8)
            ax.set_xscale("log")
            ax.invert_yaxis()
            ax.set_xlabel("Resistivity (Ω·m)")
            ax.set_ylabel("Depth (m)")
            ax.set_title("AI 1-D predicted models")
            if len(y_pred) <= 5:
                ax.legend(fontsize=8)
        elif dim == "2D":
            # y_pred shape: (n_profiles, n_stations, n_params)
            n_stations = result.get("n_stations", 20)
            n_depth    = result.get("n_depth",    40)
            rho_2d     = y_pred[0, :, :n_depth] if y_pred.ndim == 3 else y_pred[:n_stations, :n_depth]
            im = ax.imshow(
                rho_2d.T,
                aspect="auto", cmap="jet_r",
                extent=[0, n_stations, n_depth, 0],
            )
            self._tab_model.figure.colorbar(im, ax=ax, label="log₁₀ ρ (Ω·m)")
            ax.set_xlabel("Station index")
            ax.set_ylabel("Depth index")
            ax.set_title("AI 2-D predicted section (log₁₀ ρ)")
        else:
            coords = result.get("coords")
            n_params = y_pred.shape[-1] if y_pred.ndim > 1 else 1
            mid_col = n_params // 2
            vals = y_pred[:, mid_col] if y_pred.ndim > 1 else y_pred
            if coords is not None and len(coords) == len(vals):
                sc = ax.scatter(coords[:, 0], coords[:, 1], c=vals, cmap="jet_r", s=60)
                self._tab_model.figure.colorbar(sc, ax=ax, label="log₁₀ ρ mid-depth")
                ax.set_xlabel("X (m)")
                ax.set_ylabel("Y (m)")
            else:
                ax.plot(vals, "o-")
                ax.set_xlabel("Station")
                ax.set_ylabel("log₁₀ ρ mid-depth")
            ax.set_title("AI 3-D predicted resistivity")

    def _plot_ai_fit(self, ax, result: dict, freqs) -> None:
        import numpy as np
        X_obs = result.get("X_obs")
        if X_obs is None or freqs is None:
            ax.set_title("No observed data to compare")
            return
        n_freq = len(freqs)
        period = 1.0 / freqs
        # X_obs rows = [log_rho_a (n_freq), phase (n_freq)]
        for i, row in enumerate(X_obs[:min(3, len(X_obs))]):
            rho_a  = 10.0 ** row[:n_freq]
            phase  = row[n_freq:2*n_freq]
            ax.loglog(period, rho_a, "o-", ms=3, label=f"obs {i+1}")
        ax.set_xlabel("Period (s)")
        ax.set_ylabel("ρₐ (Ω·m)")
        ax.set_title("Observed apparent resistivity")
        ax.legend(fontsize=8)

    # =========================================================================
    # Helpers
    # =========================================================================

    def _selected_sites(self):
        if self._sites is None:
            return None
        selected_names = {
            item.text() for item in self._station_list.selectedItems()
        }
        if not selected_names:
            return self._sites
        try:
            return self._sites.select(list(selected_names))
        except Exception:
            return self._sites

    def set_sites(self, sites) -> None:
        self._sites = sites
        self._station_list.clear()
        if sites is None:
            return
        try:
            for site in sites.as_list():
                self._station_list.addItem(QListWidgetItem(site.name))
        except Exception:
            pass

    def _browse_workdir(self) -> None:
        d = QFileDialog.getExistingDirectory(
            self, "Select output directory",
            self._workdir_edit.text() or str(Path.home()),
        )
        if d:
            self._workdir_edit.setText(d)

    def _browse_binary(self, line_edit: QLineEdit) -> None:
        p, _ = QFileDialog.getOpenFileName(self, "Select binary", str(Path.home()))
        if p:
            line_edit.setText(p)

    def _log(self, msg: str) -> None:
        self._log_edit.appendPlainText(msg)
