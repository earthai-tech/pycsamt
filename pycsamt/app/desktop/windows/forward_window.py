# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ForwardModelWindow — interactive 1-D / 2-D / 3-D MT forward modelling.

Three-panel layout
──────────────────
  ┌──────────────┬──────────────────────────┬──────────────┐
  │ MODEL BUILDER│ RESULTS  (tab canvas)    │ LIBRARY      │
  │  ~280 px     │  stretches               │  ~210 px     │
  └──────────────┴──────────────────────────┴──────────────┘

Model Builder (left)
    • Dimension selector  1D / 2D / 3D  → switches QStackedWidget
    • 1D: layer table (ρ, h) + CRUD row buttons
    • 2D: grid params (nx, nz, dx, n_pad), background ρ, optional anomaly
    • 3D: grid params (nx, ny, nz, dx, dy, n_pad), background ρ, optional anomaly
    • Shared: frequencies (n_freq, f_min, f_max, x-axis), mode, noise
    • [▶ Compute]   [⬆ Export]

Results (centre tabs)
    1D: ρₐ/φ Curves | Model Profile | Sensitivity | vs Observed
    2D: Pseudosection | 2D Model | Profile Responses
    3D: 3D Model | Response Map | Section | Tensors

Library (right)
    • Saved-model list  → [Save] [Rename] [Delete] [Load]
    • Presets           → one button per geology prior
    • [→ Send to Inversion]
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.forward_controller import (
    GEOLOGY_PRESET_NAMES,
    ForwardController,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    icon_button,
    make_group,
)

logger = logging.getLogger(__name__)


# ── tiny helpers ──────────────────────────────────────────────────────────────


def _spin(
    value: float, lo: float, hi: float, decimals: int = 0, step: float = 1.0
) -> QDoubleSpinBox:
    sb = QDoubleSpinBox()
    sb.setRange(lo, hi)
    sb.setDecimals(decimals)
    sb.setSingleStep(step)
    sb.setValue(value)
    sb.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    return sb


def _ispin(value: int, lo: int, hi: int) -> QSpinBox:
    sb = QSpinBox()
    sb.setRange(lo, hi)
    sb.setValue(value)
    sb.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    return sb


def _label(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
    return lbl


# ── ForwardModelWindow ────────────────────────────────────────────────────────


class ForwardModelWindow(PanelWindow):
    """
    3-panel floating window for interactive MT forward modelling.

    The window overrides PanelWindow's 2-panel shell: the right
    "content" area is itself split (centre canvas tabs + right library).
    """

    # Emitted when "Send to Inversion" is clicked; carries a params dict
    send_to_inversion = Signal(dict)

    def __init__(self, parent: QWidget | None = None) -> None:
        # Must be set before super().__init__() because _build_content is
        # called inside PanelWindow._build_shell() during super().__init__().
        self._ctrl = ForwardController()
        self._worker = None
        self._last_result: Any = None

        super().__init__(
            title="Forward Modelling",
            session_key="forward_model",
            params_width=285,
            icon_name="forward",
            parent=parent,
        )
        self.resize(1300, 820)

        # Wire dim selector after all widgets exist
        self._dim_group.buttonClicked.connect(self._on_dim_changed)
        self._on_dim_changed()

    # =========================================================================
    # PanelWindow overrides
    # =========================================================================

    def _build_params(self, layout: QVBoxLayout) -> None:
        """Left panel — model builder."""
        self._build_dim_selector(layout)
        self._build_stacked_model(layout)
        self._build_freq_group(layout)
        self._build_run_buttons(layout)

    def _build_content(self, layout: QVBoxLayout) -> None:
        """Right area — inner splitter: centre tabs | library."""
        inner = QSplitter(Qt.Orientation.Horizontal)
        inner.setHandleWidth(4)

        # ── Centre: tab widget with result canvases ────────────────────
        self._tab_widget = QTabWidget()
        self._tab_widget.setObjectName("ForwardTabWidget")
        self._build_result_tabs()
        inner.addWidget(self._tab_widget)

        # ── Right: library panel ───────────────────────────────────────
        lib_widget = self._build_library_widget()
        lib_widget.setMinimumWidth(190)
        lib_widget.setMaximumWidth(240)
        inner.addWidget(lib_widget)

        inner.setStretchFactor(0, 1)
        inner.setStretchFactor(1, 0)
        inner.setSizes([900, 210])

        layout.addWidget(inner)

    # =========================================================================
    # Left panel builders
    # =========================================================================

    def _build_dim_selector(self, layout: QVBoxLayout) -> None:
        grp, lay = make_group("Dimension")
        row = QHBoxLayout()
        self._dim_group = QButtonGroup(self)
        for label, val in [("1D", "1D"), ("2D", "2D"), ("3D", "3D")]:
            rb = QRadioButton(label)
            rb.setProperty("dim", val)
            if val == "1D":
                rb.setChecked(True)
            self._dim_group.addButton(rb)
            row.addWidget(rb)
        lay.addLayout(row)
        layout.addWidget(grp)

    def _build_stacked_model(self, layout: QVBoxLayout) -> None:
        self._stacked = QStackedWidget()
        self._stacked.addWidget(self._build_1d_page())  # index 0
        self._stacked.addWidget(self._build_2d_page())  # index 1
        self._stacked.addWidget(self._build_3d_page())  # index 2
        layout.addWidget(self._stacked)

    # ── 1D page ───────────────────────────────────────────────────────────────

    def _build_1d_page(self) -> QWidget:
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)

        grp, lay = make_group("Layer Model")

        # Table: # | ρ (Ω·m) | h (m)
        self._layer_table = QTableWidget(0, 3)
        self._layer_table.setHorizontalHeaderLabels(["#", "ρ (Ω·m)", "h (m)"])
        self._layer_table.horizontalHeader().setStretchLastSection(True)
        self._layer_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self._layer_table.setMinimumHeight(140)
        self._layer_table.setMaximumHeight(220)
        self._layer_table.verticalHeader().setVisible(False)
        self._layer_table.setColumnWidth(0, 28)
        self._layer_table.setColumnWidth(1, 80)
        lay.addWidget(self._layer_table)

        # CRUD row
        btn_row = QHBoxLayout()
        self._btn_add_layer = QPushButton("+")
        self._btn_add_layer.setFixedWidth(28)
        self._btn_add_layer.setToolTip("Add layer below")
        self._btn_rem_layer = QPushButton("−")
        self._btn_rem_layer.setFixedWidth(28)
        self._btn_rem_layer.setToolTip("Remove selected layer")
        self._btn_up_layer = QPushButton("↑")
        self._btn_up_layer.setFixedWidth(28)
        self._btn_up_layer.setToolTip("Move layer up")
        self._btn_dn_layer = QPushButton("↓")
        self._btn_dn_layer.setFixedWidth(28)
        self._btn_dn_layer.setToolTip("Move layer down")
        for b in (
            self._btn_add_layer,
            self._btn_rem_layer,
            self._btn_up_layer,
            self._btn_dn_layer,
        ):
            btn_row.addWidget(b)
        btn_row.addStretch()
        lay.addLayout(btn_row)

        v.addWidget(grp)

        # Wire layer CRUD
        self._btn_add_layer.clicked.connect(self._add_layer)
        self._btn_rem_layer.clicked.connect(self._remove_layer)
        self._btn_up_layer.clicked.connect(self._move_layer_up)
        self._btn_dn_layer.clicked.connect(self._move_layer_down)

        # Default 2-layer model: 100 Ω·m / 300 m  +  1000 Ω·m halfspace
        self._populate_default_1d_model()

        return page

    def _populate_default_1d_model(self) -> None:
        rows = [
            (100.0, 300.0),
            (1000.0, None),  # halfspace — no thickness
        ]
        self._layer_table.setRowCount(0)
        for rho, h in rows:
            self._insert_layer_row(rho, h)

    def _insert_layer_row(self, rho: float, h: float | None) -> None:
        r = self._layer_table.rowCount()
        self._layer_table.insertRow(r)
        idx_item = QTableWidgetItem(str(r + 1))
        idx_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
        self._layer_table.setItem(r, 0, idx_item)
        self._layer_table.setItem(r, 1, QTableWidgetItem(f"{rho:.1f}"))
        h_str = f"{h:.1f}" if h is not None else "∞"
        self._layer_table.setItem(r, 2, QTableWidgetItem(h_str))
        self._renumber_layers()

    def _renumber_layers(self) -> None:
        n = self._layer_table.rowCount()
        for i in range(n):
            item = self._layer_table.item(i, 0)
            if item:
                item.setText(str(i + 1))
            # Halfspace always has ∞ thickness
            if i == n - 1:
                h_item = self._layer_table.item(i, 2)
                if h_item and h_item.text() not in ("∞", "inf"):
                    h_item.setText("∞")

    def _add_layer(self) -> None:
        sel = self._layer_table.currentRow()
        n = self._layer_table.rowCount()
        # Insert before halfspace (last row)
        insert_at = max(0, min(n - 1, sel + 1) if sel >= 0 else n - 1)
        self._layer_table.insertRow(insert_at)
        idx = QTableWidgetItem(str(insert_at + 1))
        idx.setFlags(Qt.ItemFlag.ItemIsEnabled)
        self._layer_table.setItem(insert_at, 0, idx)
        self._layer_table.setItem(insert_at, 1, QTableWidgetItem("100.0"))
        self._layer_table.setItem(insert_at, 2, QTableWidgetItem("200.0"))
        self._renumber_layers()

    def _remove_layer(self) -> None:
        if self._layer_table.rowCount() <= 2:
            return
        sel = self._layer_table.currentRow()
        if sel < 0:
            return
        # Can't remove halfspace
        if sel == self._layer_table.rowCount() - 1:
            return
        self._layer_table.removeRow(sel)
        self._renumber_layers()

    def _move_layer_up(self) -> None:
        sel = self._layer_table.currentRow()
        if sel <= 0:
            return
        self._swap_rows(sel - 1, sel)
        self._layer_table.setCurrentCell(sel - 1, 1)

    def _move_layer_down(self) -> None:
        n = self._layer_table.rowCount()
        sel = self._layer_table.currentRow()
        if sel < 0 or sel >= n - 2:  # can't move halfspace
            return
        self._swap_rows(sel, sel + 1)
        self._layer_table.setCurrentCell(sel + 1, 1)

    def _swap_rows(self, a: int, b: int) -> None:
        t = self._layer_table
        for col in (1, 2):
            ia = t.item(a, col)
            ib = t.item(b, col)
            ta = ia.text() if ia else ""
            tb = ib.text() if ib else ""
            if ia:
                ia.setText(tb)
            if ib:
                ib.setText(ta)

    def _read_1d_model(self) -> tuple[np.ndarray, np.ndarray]:
        n = self._layer_table.rowCount()
        rho = []
        h = []
        for i in range(n):
            r_item = self._layer_table.item(i, 1)
            h_item = self._layer_table.item(i, 2)
            r_val = float(r_item.text() if r_item else "100")
            rho.append(max(r_val, 1e-3))
            if i < n - 1:
                h_txt = h_item.text() if h_item else "100"
                h_val = float(h_txt) if h_txt not in ("∞", "inf") else 100.0
                h.append(max(h_val, 1.0))
        return np.array(rho), np.array(h)

    def _set_1d_model(self, rho: list[float], thick: list[float]) -> None:
        self._layer_table.setRowCount(0)
        for i, r in enumerate(rho):
            h = thick[i] if i < len(thick) else None
            self._insert_layer_row(r, h)

    # ── 2D page ───────────────────────────────────────────────────────────────

    def _build_2d_page(self) -> QWidget:
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)

        # Grid group
        grp_g, lay_g = make_group("Grid")
        form_g = QFormLayout()
        form_g.setSpacing(4)
        self._2d_nx = _ispin(30, 5, 200)
        self._2d_nz = _ispin(20, 5, 100)
        self._2d_dx = _spin(500.0, 10.0, 50000.0, 0, 50.0)
        self._2d_dzmin = _spin(50.0, 1.0, 1000.0, 0, 10.0)
        self._2d_dzmax = _spin(1000.0, 10.0, 50000.0, 0, 100.0)
        self._2d_npad = _ispin(5, 2, 20)
        form_g.addRow("nx cells:", self._2d_nx)
        form_g.addRow("nz cells:", self._2d_nz)
        form_g.addRow("dx (m):", self._2d_dx)
        form_g.addRow("dz min (m):", self._2d_dzmin)
        form_g.addRow("dz max (m):", self._2d_dzmax)
        form_g.addRow("n_pad:", self._2d_npad)
        lay_g.addLayout(form_g)
        v.addWidget(grp_g)

        # Background
        grp_b, lay_b = make_group("Background")
        form_b = QFormLayout()
        form_b.setSpacing(4)
        self._2d_bgrho = _spin(100.0, 0.1, 100000.0, 1, 10.0)
        form_b.addRow("ρ (Ω·m):", self._2d_bgrho)
        self._2d_nsta = _ispin(10, 2, 100)
        form_b.addRow("Stations:", self._2d_nsta)
        lay_b.addLayout(form_b)
        v.addWidget(grp_b)

        # Anomaly
        grp_a, lay_a = make_group("Anomaly (optional)")
        self._2d_anom_cb = QCheckBox("Include anomaly")
        lay_a.addWidget(self._2d_anom_cb)
        self._2d_anom_widget = QWidget()
        form_a = QFormLayout(self._2d_anom_widget)
        form_a.setSpacing(4)
        self._2d_ax = _spin(7500.0, 0.0, 1e6, 0, 500.0)
        self._2d_az = _spin(500.0, 0.0, 1e5, 0, 50.0)
        self._2d_aw = _spin(2000.0, 10.0, 1e5, 0, 100.0)
        self._2d_ah = _spin(1000.0, 10.0, 1e5, 0, 100.0)
        self._2d_arho = _spin(10.0, 0.1, 1e5, 1, 5.0)
        form_a.addRow("Centre x (m):", self._2d_ax)
        form_a.addRow("Top z (m):", self._2d_az)
        form_a.addRow("Width (m):", self._2d_aw)
        form_a.addRow("Height (m):", self._2d_ah)
        form_a.addRow("ρ (Ω·m):", self._2d_arho)
        self._2d_anom_widget.setVisible(False)
        lay_a.addWidget(self._2d_anom_widget)
        self._2d_anom_cb.toggled.connect(self._2d_anom_widget.setVisible)
        v.addWidget(grp_a)

        v.addStretch(1)
        return page

    # ── 3D page ───────────────────────────────────────────────────────────────

    def _build_3d_page(self) -> QWidget:
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)

        # Grid group
        grp_g, lay_g = make_group("Grid")
        form_g = QFormLayout()
        form_g.setSpacing(4)
        self._3d_nx = _ispin(12, 4, 80)
        self._3d_ny = _ispin(12, 4, 80)
        self._3d_nz = _ispin(10, 4, 50)
        self._3d_dx = _spin(1000.0, 50.0, 50000.0, 0, 100.0)
        self._3d_dy = _spin(1000.0, 50.0, 50000.0, 0, 100.0)
        self._3d_npad = _ispin(4, 2, 12)
        form_g.addRow("nx:", self._3d_nx)
        form_g.addRow("ny:", self._3d_ny)
        form_g.addRow("nz:", self._3d_nz)
        form_g.addRow("dx (m):", self._3d_dx)
        form_g.addRow("dy (m):", self._3d_dy)
        form_g.addRow("n_pad:", self._3d_npad)
        lay_g.addLayout(form_g)
        v.addWidget(grp_g)

        # Background + stations
        grp_b, lay_b = make_group("Background & Stations")
        form_b = QFormLayout()
        form_b.setSpacing(4)
        self._3d_bgrho = _spin(100.0, 0.1, 1e5, 1, 10.0)
        self._3d_nsx = _ispin(4, 1, 20)
        self._3d_nsy = _ispin(4, 1, 20)
        self._3d_stasp = _spin(1000.0, 50.0, 1e5, 0, 100.0)
        form_b.addRow("Background ρ:", self._3d_bgrho)
        form_b.addRow("Stations X:", self._3d_nsx)
        form_b.addRow("Stations Y:", self._3d_nsy)
        form_b.addRow("Spacing (m):", self._3d_stasp)
        lay_b.addLayout(form_b)
        v.addWidget(grp_b)

        # Anomaly
        grp_a, lay_a = make_group("Anomaly (optional)")
        self._3d_anom_cb = QCheckBox("Include anomaly")
        lay_a.addWidget(self._3d_anom_cb)
        self._3d_anom_widget = QWidget()
        form_a = QFormLayout(self._3d_anom_widget)
        form_a.setSpacing(4)
        self._3d_ax = _spin(6000.0, 0.0, 1e6, 0, 500.0)
        self._3d_ay = _spin(6000.0, 0.0, 1e6, 0, 500.0)
        self._3d_az = _spin(500.0, 0.0, 1e5, 0, 50.0)
        self._3d_arho = _spin(10.0, 0.1, 1e5, 1, 5.0)
        form_a.addRow("Centre x (m):", self._3d_ax)
        form_a.addRow("Centre y (m):", self._3d_ay)
        form_a.addRow("Top z (m):", self._3d_az)
        form_a.addRow("ρ (Ω·m):", self._3d_arho)
        self._3d_anom_widget.setVisible(False)
        lay_a.addWidget(self._3d_anom_widget)
        self._3d_anom_cb.toggled.connect(self._3d_anom_widget.setVisible)
        v.addWidget(grp_a)

        v.addStretch(1)
        return page

    # ── Shared frequency / mode / noise group ─────────────────────────────────

    def _build_freq_group(self, layout: QVBoxLayout) -> None:
        grp, lay = make_group("Frequencies")
        form = QFormLayout()
        form.setSpacing(4)

        self._nfreq = _ispin(30, 5, 200)
        self._fmin = _spin(1e-3, 1e-4, 1e3, 4, 1e-4)  # 4 decimals, min 0.0001
        self._fmax = _spin(1e3, 1e-3, 1e6, 1, 10.0)

        self._axis_combo = QComboBox()
        self._axis_combo.addItems(["Period (s)", "Frequency (Hz)"])

        self._mode_combo = QComboBox()
        self._mode_combo.addItems(["MT", "CSAMT", "TEM"])

        self._noise_combo = QComboBox()
        self._noise_combo.addItems(["None", "Gaussian 5%", "Multiplicative 5%"])

        form.addRow("n_freq:", self._nfreq)
        form.addRow("f_min (Hz):", self._fmin)
        form.addRow("f_max (Hz):", self._fmax)
        form.addRow("x-axis:", self._axis_combo)
        form.addRow("Mode:", self._mode_combo)
        form.addRow("Noise:", self._noise_combo)
        lay.addLayout(form)
        layout.addWidget(grp)

    def _build_run_buttons(self, layout: QVBoxLayout) -> None:
        self._btn_compute = icon_button(
            "▶  Compute", "results", "Run the forward solver"
        )
        self._btn_compute.setObjectName("ComputeButton")
        self._btn_export = icon_button("⬆  Export", "export", "Export current figure")
        layout.addWidget(self._btn_compute)
        layout.addWidget(self._btn_export)
        self._btn_compute.clicked.connect(self._on_compute)
        self._btn_export.clicked.connect(self._on_export)

        self._compute_label = QLabel("")
        self._compute_label.setObjectName("ComputeLabel")
        self._compute_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self._compute_label)

    # =========================================================================
    # Centre: result tabs
    # =========================================================================

    def _build_result_tabs(self) -> None:
        """Create all canvas tabs for 1D, 2D, 3D results."""
        # ── 1D canvases ───────────────────────────────────────────────
        self._c1d_curves = MplCanvas(toolbar=True)
        self._c1d_model = MplCanvas(toolbar=True)
        self._c1d_sensitivity = MplCanvas(toolbar=True)
        self._c1d_observed = MplCanvas(toolbar=True)

        # ── 2D canvases ───────────────────────────────────────────────
        self._c2d_pseudo = MplCanvas(toolbar=True)
        self._c2d_model = MplCanvas(toolbar=True)
        self._c2d_profiles = MplCanvas(toolbar=True)

        # ── 3D canvases ───────────────────────────────────────────────
        self._c3d_model = MplCanvas(toolbar=True)
        self._c3d_map = MplCanvas(toolbar=True)
        self._c3d_section = MplCanvas(toolbar=True)
        self._c3d_tensors = MplCanvas(toolbar=True)

        # Initially show 1D tabs
        self._set_tabs_for_dim("1D")

    def _set_tabs_for_dim(self, dim: str) -> None:
        self._tab_widget.clear()
        if dim == "1D":
            self._tab_widget.addTab(self._c1d_curves, "ρₐ / φ Curves")
            self._tab_widget.addTab(self._c1d_model, "Model Profile")
            self._tab_widget.addTab(self._c1d_sensitivity, "Sensitivity")
            self._tab_widget.addTab(self._c1d_observed, "vs Observed")
        elif dim == "2D":
            self._tab_widget.addTab(self._c2d_pseudo, "Pseudosection")
            self._tab_widget.addTab(self._c2d_model, "2D Model")
            self._tab_widget.addTab(self._c2d_profiles, "Profile Responses")
        else:
            self._tab_widget.addTab(self._c3d_model, "3D Model")
            self._tab_widget.addTab(self._c3d_map, "Response Map")
            self._tab_widget.addTab(self._c3d_section, "Section")
            self._tab_widget.addTab(self._c3d_tensors, "Tensors")

    # =========================================================================
    # Right: library panel
    # =========================================================================

    def _build_library_widget(self) -> QWidget:
        widget = QWidget()
        widget.setObjectName("LibraryPanel")
        v = QVBoxLayout(widget)
        v.setContentsMargins(4, 4, 4, 4)
        v.setSpacing(4)

        # ── Saved models ──────────────────────────────────────────────
        saved_grp, saved_lay = make_group("Saved Models")

        self._lib_list = QListWidget()
        self._lib_list.setObjectName("LibraryList")
        self._lib_list.setMinimumHeight(100)
        self._lib_list.setMaximumHeight(160)
        self._lib_list.itemDoubleClicked.connect(self._load_model_from_lib)
        saved_lay.addWidget(self._lib_list)
        self._refresh_lib_list()

        btn_save = QPushButton("Save current")
        btn_rename = QPushButton("Rename")
        btn_delete = QPushButton("Delete")
        btn_load = QPushButton("Load →")
        for b in (btn_save, btn_rename, btn_delete, btn_load):
            b.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            saved_lay.addWidget(b)
        btn_save.clicked.connect(self._save_current_model)
        btn_rename.clicked.connect(self._rename_selected_model)
        btn_delete.clicked.connect(self._delete_selected_model)
        btn_load.clicked.connect(self._load_selected_model)
        v.addWidget(saved_grp)

        # ── Presets ───────────────────────────────────────────────────
        preset_grp, preset_lay = make_group("Presets (1D)")

        scroll_inner = QWidget()
        scroll_v = QVBoxLayout(scroll_inner)
        scroll_v.setContentsMargins(0, 0, 0, 0)
        scroll_v.setSpacing(2)

        for name in GEOLOGY_PRESET_NAMES:
            btn = QPushButton(name.capitalize())
            btn.setToolTip(f"Load {name} geological prior")
            btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            btn.clicked.connect(lambda _, n=name: self._load_preset(n))
            scroll_v.addWidget(btn)
        scroll_v.addStretch(1)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setWidget(scroll_inner)
        scroll.setMaximumHeight(240)
        preset_lay.addWidget(scroll)
        v.addWidget(preset_grp)

        # ── Send to Inversion ─────────────────────────────────────────
        self._btn_send_inv = QPushButton("→ Send to Inversion")
        self._btn_send_inv.setObjectName("SendToInversionButton")
        self._btn_send_inv.setToolTip(
            "Send current model as starting model for the Inversion window"
        )
        self._btn_send_inv.clicked.connect(self._on_send_to_inversion)
        v.addWidget(self._btn_send_inv)

        v.addStretch(1)
        return widget

    # =========================================================================
    # Dimension switching
    # =========================================================================

    def _current_dim(self) -> str:
        btn = self._dim_group.checkedButton()
        return btn.property("dim") if btn else "1D"

    def _on_dim_changed(self, *_) -> None:
        dim = self._current_dim()
        idx = {"1D": 0, "2D": 1, "3D": 2}[dim]
        self._stacked.setCurrentIndex(idx)
        self._set_tabs_for_dim(dim)
        # 3D is slow — warn the user
        if dim == "3D":
            self._compute_label.setText("3D can take several minutes.")
        else:
            self._compute_label.setText("")

    # =========================================================================
    # Compute
    # =========================================================================

    def _build_params_dict(self) -> dict:
        dim = self._current_dim()
        noise_map = {
            "None": "none",
            "Gaussian 5%": "gaussian",
            "Multiplicative 5%": "multiplicative",
        }
        p: dict = {
            "dim": dim,
            "n_freq": self._nfreq.value(),
            "f_min": self._fmin.value(),
            "f_max": self._fmax.value(),
            "noise": noise_map.get(self._noise_combo.currentText(), "none"),
        }
        if dim == "1D":
            rho, h = self._read_1d_model()
            p["resistivity"] = rho.tolist()
            p["thickness"] = h.tolist()
        elif dim == "2D":
            p.update(
                {
                    "nx": self._2d_nx.value(),
                    "nz": self._2d_nz.value(),
                    "dx": self._2d_dx.value(),
                    "dz_min": self._2d_dzmin.value(),
                    "dz_max": self._2d_dzmax.value(),
                    "n_pad": self._2d_npad.value(),
                    "bg_rho": self._2d_bgrho.value(),
                    "n_stations": int(self._2d_nsta.value()),
                    "anomaly": self._2d_anom_cb.isChecked(),
                    "anom_x": self._2d_ax.value(),
                    "anom_z": self._2d_az.value(),
                    "anom_w": self._2d_aw.value(),
                    "anom_h": self._2d_ah.value(),
                    "anom_rho": self._2d_arho.value(),
                }
            )
        else:
            p.update(
                {
                    "nx": self._3d_nx.value(),
                    "ny": self._3d_ny.value(),
                    "nz": self._3d_nz.value(),
                    "dx": self._3d_dx.value(),
                    "dy": self._3d_dy.value(),
                    "n_pad": self._3d_npad.value(),
                    "bg_rho": self._3d_bgrho.value(),
                    "n_sx": self._3d_nsx.value(),
                    "n_sy": self._3d_nsy.value(),
                    "sta_spacing": self._3d_stasp.value(),
                    "anomaly": self._3d_anom_cb.isChecked(),
                    "anom_x": self._3d_ax.value(),
                    "anom_y": self._3d_ay.value(),
                    "anom_z": self._3d_az.value(),
                    "anom_rho": self._3d_arho.value(),
                }
            )
        return p

    def _on_compute(self) -> None:
        if self._worker and self._worker.isRunning():
            return
        from pycsamt.app.desktop.workers.forward_worker import (
            ForwardWorker,
        )

        params = self._build_params_dict()
        self._btn_compute.setEnabled(False)
        self._compute_label.setText("Computing…")
        self._worker = ForwardWorker(params, parent=self)
        self._worker.finished.connect(self._on_result)
        self._worker.error.connect(self._on_error)
        self._worker.progress.connect(
            lambda v: self._compute_label.setText(f"Computing… {v}%")
        )
        self._worker.start()

    def _on_result(self, result) -> None:
        self._btn_compute.setEnabled(True)
        self._compute_label.setText("Done.")
        self._last_result = result
        dim = self._current_dim()
        try:
            if dim == "1D":
                self._plot_1d(result)
            elif dim == "2D":
                self._plot_2d(result)
            else:
                self._plot_3d(result)
        except Exception as exc:
            logger.exception("Plot error")
            self._compute_label.setText(f"Plot error: {exc}")

    def _on_error(self, msg: str) -> None:
        self._btn_compute.setEnabled(True)
        self._compute_label.setText(f"Error: {msg}")
        QMessageBox.critical(self, "Forward Error", msg)

    # =========================================================================
    # Plotting helpers
    # =========================================================================

    def _plot_1d(self, resp) -> None:
        use_period = self._axis_combo.currentIndex() == 0
        x = 1.0 / resp.freqs if use_period else resp.freqs
        xlabel = "Period (s)" if use_period else "Frequency (Hz)"

        # Tab 0: ρₐ/φ curves
        ax1, ax2 = self._reset_canvas_2ax(self._c1d_curves)
        ax1.loglog(x, resp.rho_a, "b-o", markersize=4, label="ρₐ")
        ax1.set_ylabel("ρₐ (Ω·m)")
        ax1.set_xlabel(xlabel)
        ax1.legend(fontsize=8)
        ax2.semilogx(x, resp.phase, "r-s", markersize=4, label="φ")
        ax2.set_ylabel("Phase (°)")
        ax2.set_xlabel(xlabel)
        ax2.legend(fontsize=8)
        ax2.set_ylim(0, 90)
        self._c1d_curves.figure.tight_layout()
        self._c1d_curves.draw()

        # Tab 1: model profile
        rho, h = self._read_1d_model()
        self._c1d_model.figure.clear()
        ax = self._c1d_model.figure.add_subplot(111)
        from pycsamt.forward.plot import plot_model_1d

        try:
            plot_model_1d(ax=ax, resistivity=rho, thickness=h, label="model")
        except Exception:
            depths = np.concatenate([[0], np.cumsum(h), [np.sum(h) * 1.5]])
            for i, _r in enumerate(rho):
                ax.axhline(y=depths[i], color="gray", lw=0.5, ls="--")
            np.repeat(rho, 2)[1:-1]
            np.repeat(depths[:-1], 2)[1:]
            ax.plot(rho, depths[:-1], drawstyle="steps-post")
            ax.set_xscale("log")
            ax.invert_yaxis()
            ax.set_xlabel("ρ (Ω·m)")
            ax.set_ylabel("Depth (m)")
        self._c1d_model.figure.tight_layout()
        self._c1d_model.draw()

        # Tab 2: sensitivity (approximate — dρₐ/dρ per layer, numeric)
        try:
            self._plot_1d_sensitivity(resp)
        except Exception:
            pass

    def _reset_canvas_2ax(self, canvas: MplCanvas):
        canvas.figure.clear()
        ax1 = canvas.figure.add_subplot(211)
        ax2 = canvas.figure.add_subplot(212)
        return ax1, ax2

    def _plot_1d_sensitivity(self, resp) -> None:
        from pycsamt.forward import LayeredModel, MT1DForward

        rho, h = self._read_1d_model()
        n_layers = len(rho)
        n_freq = len(resp.freqs)
        sens = np.zeros((n_layers, n_freq))
        eps = 0.05  # 5% perturbation

        for i in range(n_layers):
            rho_p = rho.copy()
            rho_p[i] *= 1 + eps
            resp_p = MT1DForward(resp.freqs).run(
                LayeredModel(resistivity=rho_p, thickness=h)
            )
            sens[i] = (np.log(resp_p.rho_a) - np.log(resp.rho_a)) / eps

        self._c1d_sensitivity.figure.clear()
        ax = self._c1d_sensitivity.figure.add_subplot(111)
        im = ax.imshow(
            sens,
            aspect="auto",
            cmap="RdBu_r",
            interpolation="nearest",
            extent=[0, n_freq, n_layers, 0],
        )
        self._c1d_sensitivity.figure.colorbar(im, ax=ax, label="d ln ρₐ / d ln ρ")
        ax.set_xlabel("Frequency index (low → high)")
        ax.set_ylabel("Layer index (top → bottom)")
        ax.set_title("Jacobian sensitivity")
        self._c1d_sensitivity.figure.tight_layout()
        self._c1d_sensitivity.draw()

    def _plot_2d(self, resp) -> None:
        from pycsamt.forward.plot import (
            plot_model_2d,
            plot_pseudosection_2d,
            plot_response_profiles,
        )

        # Tab 0: pseudosection
        self._c2d_pseudo.figure.clear()
        ax = self._c2d_pseudo.figure.add_subplot(111)
        try:
            plot_pseudosection_2d(resp, ax=ax, mode="TE")
        except Exception as exc:
            ax.set_title(f"Pseudosection unavailable: {exc}")
        self._c2d_pseudo.figure.tight_layout()
        self._c2d_pseudo.draw()

        # Tab 1: 2D model
        self._c2d_model.figure.clear()
        ax = self._c2d_model.figure.add_subplot(111)
        try:
            plot_model_2d(resp.grid, ax=ax)
        except Exception as exc:
            ax.set_title(f"Model plot unavailable: {exc}")
        self._c2d_model.figure.tight_layout()
        self._c2d_model.draw()

        # Tab 2: station profiles
        self._c2d_profiles.figure.clear()
        ax = self._c2d_profiles.figure.add_subplot(111)
        try:
            plot_response_profiles(resp, ax=ax, mode="TE")
        except Exception as exc:
            ax.set_title(f"Profile plot unavailable: {exc}")
        self._c2d_profiles.figure.tight_layout()
        self._c2d_profiles.draw()

    def _plot_3d(self, resp) -> None:
        import io

        import matplotlib.pyplot as plt
        from matplotlib.image import imread as _imread

        from pycsamt.forward.plot import (
            plot_model_3d,
            plot_response_map_3d,
            plot_response_section_3d,
            plot_tensor_components_3d,
        )

        def _render_own_figure(canvas, func, *args, label="", **kwargs):
            """Call a plot function that creates its own figure, then rasterise
            it into the canvas so it lives inside the existing MplCanvas."""
            canvas.figure.clear()
            ax_img = canvas.figure.add_subplot(111)
            try:
                result = func(*args, **kwargs)
                # result is an ndarray of Axes — grab their parent figure
                if hasattr(result, "flat"):
                    src_fig = next(iter(result.flat)).figure
                else:
                    src_fig = result.figure
                buf = io.BytesIO()
                src_fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
                plt.close(src_fig)
                buf.seek(0)
                img = _imread(buf)
                ax_img.imshow(img)
                ax_img.axis("off")
            except Exception as exc:
                ax_img.set_title(f"{label}: {exc}")
            canvas.figure.tight_layout(pad=0)
            canvas.draw()

        # ── Tab 0: 3D Model (no ax, takes grid3d) ──────────────────────
        _render_own_figure(
            self._c3d_model,
            plot_model_3d,
            resp.grid,
            label="3D Model",
        )

        # ── Tab 1: Response Map (accepts ax) ───────────────────────────
        self._c3d_map.figure.clear()
        ax = self._c3d_map.figure.add_subplot(111)
        try:
            plot_response_map_3d(resp, ax=ax)
        except Exception as exc:
            ax.set_title(f"Response Map: {exc}")
        self._c3d_map.figure.tight_layout()
        self._c3d_map.draw()

        # ── Tab 2: Section (accepts ax) ────────────────────────────────
        self._c3d_section.figure.clear()
        ax = self._c3d_section.figure.add_subplot(111)
        try:
            plot_response_section_3d(resp, ax=ax)
        except Exception as exc:
            ax.set_title(f"Section: {exc}")
        self._c3d_section.figure.tight_layout()
        self._c3d_section.draw()

        # ── Tab 3: Tensor components (no ax, multi-panel) ──────────────
        _render_own_figure(
            self._c3d_tensors,
            plot_tensor_components_3d,
            resp,
            label="Tensors",
        )

    # =========================================================================
    # Library actions
    # =========================================================================

    def _refresh_lib_list(self) -> None:
        self._lib_list.clear()
        for name in self._ctrl.model_names:
            self._lib_list.addItem(QListWidgetItem(name))

    def _save_current_model(self) -> None:
        dim = self._current_dim()
        if dim != "1D":
            QMessageBox.information(
                self,
                "Save Model",
                "Library save is supported for 1D models.\n"
                "2D/3D parameters are reproduced from the builder controls.",
            )
            return
        rho, h = self._read_1d_model()
        name, ok = QInputDialog.getText(self, "Save Model", "Model name:")
        if not ok or not name.strip():
            return
        record = ForwardController.model_to_record(name.strip(), "1D", rho, h)
        self._ctrl.save_model(name.strip(), record)
        self._refresh_lib_list()

    def _selected_lib_name(self) -> str | None:
        item = self._lib_list.currentItem()
        return item.text() if item else None

    def _rename_selected_model(self) -> None:
        name = self._selected_lib_name()
        if not name:
            return
        new_name, ok = QInputDialog.getText(
            self, "Rename Model", "New name:", text=name
        )
        if ok and new_name.strip():
            self._ctrl.rename_model(name, new_name.strip())
            self._refresh_lib_list()

    def _delete_selected_model(self) -> None:
        name = self._selected_lib_name()
        if not name:
            return
        r = QMessageBox.question(
            self,
            "Delete Model",
            f"Delete '{name}'?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if r == QMessageBox.StandardButton.Yes:
            self._ctrl.delete_model(name)
            self._refresh_lib_list()

    def _load_selected_model(self) -> None:
        name = self._selected_lib_name()
        if name:
            self._load_model_by_name(name)

    def _load_model_from_lib(self, item: QListWidgetItem) -> None:
        self._load_model_by_name(item.text())

    def _load_model_by_name(self, name: str) -> None:
        record = self._ctrl.get_model_record(name)
        if not record:
            return
        dim = record.get("dim", "1D")
        if dim == "1D":
            rho, h = ForwardController.record_to_arrays(record)
            # Switch to 1D view
            for btn in self._dim_group.buttons():
                if btn.property("dim") == "1D":
                    btn.setChecked(True)
            self._on_dim_changed()
            self._set_1d_model(rho.tolist(), h.tolist())

    # ── Geology presets ───────────────────────────────────────────────────────

    def _load_preset(self, name: str) -> None:
        try:
            data = ForwardController.build_preset_1d(name)
            # Switch to 1D if not already
            for btn in self._dim_group.buttons():
                if btn.property("dim") == "1D":
                    btn.setChecked(True)
            self._on_dim_changed()
            self._set_1d_model(data["resistivity"], data["thickness"])
            self._compute_label.setText(f"Preset loaded: {name}")
        except Exception as exc:
            QMessageBox.warning(self, "Preset Error", str(exc))

    # =========================================================================
    # Send to Inversion
    # =========================================================================

    def _on_send_to_inversion(self) -> None:
        """Serialise current model and emit signal for MainWindow to open Inversion."""
        dim = self._current_dim()
        payload: dict = {"dim": dim}
        if dim == "1D":
            rho, h = self._read_1d_model()
            payload["resistivity"] = rho.tolist()
            payload["thickness"] = h.tolist()
        elif dim == "2D":
            payload["bg_rho"] = self._2d_bgrho.value()
        else:
            payload["bg_rho"] = self._3d_bgrho.value()
        if self._last_result is not None:
            payload["has_result"] = True
        self.send_to_inversion.emit(payload)

    # =========================================================================
    # Export
    # =========================================================================

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )

        current = self._tab_widget.currentWidget()
        if not isinstance(current, MplCanvas):
            return
        ExportDialog(figure=current.figure, parent=self).exec()

    # =========================================================================
    # Public API (called by MainWindow)
    # =========================================================================

    def set_observed_sites(self, sites) -> None:
        """Feed real EDI data for the 'vs Observed' overlay tab."""
        self._sites = sites

    def load_starting_model(self, payload: dict) -> None:
        """Called when InversionWindow feeds a starting model back here."""
        dim = payload.get("dim", "1D")
        if dim == "1D" and "resistivity" in payload:
            for btn in self._dim_group.buttons():
                if btn.property("dim") == "1D":
                    btn.setChecked(True)
            self._on_dim_changed()
            self._set_1d_model(payload["resistivity"], payload.get("thickness", []))
