# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
InterpretationWindow — floating Interpretation Studio panel.

Layout
──────
    ┌─────────────────────┬───────────────────────────────────────────────┐
    │  WORKFLOW  (300px)  │  PLOT STUDIO  (expands)                       │
    │  ┌─────────────────┐│  ┌──────────────────────────────────────────┐ │
    │  │ Model Status    ││  │  Tab 1  │  Tab 2  │  …  │  +             │ │
    │  └─────────────────┘│  │  MplCanvas + NavigationToolbar            │ │
    │  Category list       │  └──────────────────────────────────────────┘ │
    │  Plot selector combo │  ┌── Gallery strip ────────────────────────┐ │
    │  Param form          │  │  [thumb] [thumb] [thumb] …              │ │
    │  ─────────────────── │  └────────────────────────────────────────── │
    │  Run actions         │                                               │
    └─────────────────────┴───────────────────────────────────────────────┘

Every "Generate" click adds a new tab in the plot studio and a thumbnail in
the gallery.  Tabs can be closed, renamed (double-click), or exported.
"""
from __future__ import annotations

import io
from pathlib import Path
from typing import List

from PySide6.QtCore  import Qt, QSize, QThread, Signal
from PySide6.QtGui   import QIcon, QPixmap
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QTabBar,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.interp_controller import (
    InterpController,
    WORKFLOW_CATALOGUE,
    CATEGORIES,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base      import PanelWindow, make_group, icon_button


# ── Background workers ────────────────────────────────────────────────────────

class _RunWorker(QThread):
    """Runs a heavy InterpController computation off the GUI thread."""
    finished = Signal(str)

    def __init__(self, fn, parent=None):
        super().__init__(parent)
        self._fn = fn

    def run(self):
        try:
            msg = self._fn()
            self.finished.emit(msg or "Done.")
        except Exception as exc:
            self.finished.emit(f"Error: {exc}")


class _GeneratePlotWorker(QThread):
    """Generates a matplotlib Figure off the GUI thread, then signals it."""
    finished = Signal(object, str)   # (Figure, label)
    failed   = Signal(str)

    def __init__(self, ctrl, method_name: str, label: str, kwargs: dict,
                 parent=None):
        super().__init__(parent)
        self._ctrl        = ctrl
        self._method_name = method_name
        self._label       = label
        self._kwargs      = kwargs

    def run(self):
        import matplotlib
        matplotlib.use("Agg")   # ensure non-interactive backend in thread
        try:
            fig = self._ctrl.generate(self._method_name, **self._kwargs)
            self.finished.emit(fig, self._label)
        except Exception as exc:
            self.failed.emit(str(exc))


# ── Interpretation Window ─────────────────────────────────────────────────────

class InterpretationWindow(PanelWindow):
    """
    Floating Interpretation Studio.

    Left:  workflow navigator · parameters · run actions
    Right: tabbed plot studio + gallery strip
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(
            title="Interpretation",
            session_key="interpretation",
            params_width=300,
            icon_name="interpret",
            parent=parent,
        )
        self.resize(1400, 860)
        self._ctrl = InterpController()
        self._worker: QThread | None = None
        self._plot_worker: _GeneratePlotWorker | None = None
        self._populate_category_combo()
        self._on_category_changed(0)

    # ── Left: params panel ────────────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:

        # ── Model status card ─────────────────────────────────────────────
        grp_status, lay_status = make_group("Model Status")
        self._status_card = QLabel("No model loaded\nOpen EDI data first")
        self._status_card.setObjectName("ModelStatusCard")
        self._status_card.setWordWrap(True)
        self._status_card.setAlignment(Qt.AlignmentFlag.AlignTop)
        lay_status.addWidget(self._status_card)

        load_row = QHBoxLayout()
        btn_from_inv = QPushButton("← From Inversion")
        btn_from_inv.setObjectName("FileListBtn")
        btn_from_inv.setToolTip("Use the model from the Inversion window")
        btn_from_inv.clicked.connect(self._load_from_inversion)
        btn_load_occ = QPushButton("Load Occam2D…")
        btn_load_occ.setObjectName("FileListBtn")
        btn_load_occ.clicked.connect(self._load_occam2d)
        load_row.addWidget(btn_from_inv)
        load_row.addWidget(btn_load_occ)
        lay_status.addLayout(load_row)
        layout.addWidget(grp_status)

        # ── Category navigator ────────────────────────────────────────────
        grp_nav, lay_nav = make_group("Workflow")
        self._combo_category = QComboBox()
        self._combo_category.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._combo_category.currentIndexChanged.connect(self._on_category_changed)
        lay_nav.addWidget(self._combo_category)

        self._combo_plot = QComboBox()
        self._combo_plot.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        lay_nav.addWidget(self._combo_plot)

        self._plot_desc = QLabel("")
        self._plot_desc.setObjectName("InfoLabel")
        self._plot_desc.setWordWrap(True)
        lay_nav.addWidget(self._plot_desc)
        layout.addWidget(grp_nav)

        # ── Context-sensitive parameters ──────────────────────────────────
        self._grp_params, self._lay_params_inner = make_group("Parameters")
        self._param_form = QFormLayout()
        self._param_form.setSpacing(4)
        self._param_form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        self._lay_params_inner.addLayout(self._param_form)
        layout.addWidget(self._grp_params)

        # ── Borehole management ───────────────────────────────────────────
        grp_bh, lay_bh = make_group("Boreholes")
        self._bh_list = QListWidget()
        self._bh_list.setObjectName("StackList")
        self._bh_list.setMaximumHeight(70)
        bh_row = QHBoxLayout()
        btn_add_csv = QPushButton("+ CSV")
        btn_add_csv.setObjectName("FileListBtn")
        btn_add_csv.clicked.connect(self._add_borehole_csv)
        btn_add_las = QPushButton("+ LAS")
        btn_add_las.setObjectName("FileListBtn")
        btn_add_las.clicked.connect(self._add_borehole_las)
        btn_rem_bh  = QPushButton("Remove")
        btn_rem_bh.setObjectName("FileListBtn")
        btn_rem_bh.clicked.connect(self._remove_borehole)
        bh_row.addWidget(btn_add_csv)
        bh_row.addWidget(btn_add_las)
        bh_row.addWidget(btn_rem_bh)
        lay_bh.addWidget(self._bh_list)
        lay_bh.addLayout(bh_row)
        layout.addWidget(grp_bh)

        # ── Run actions ───────────────────────────────────────────────────
        grp_run, lay_run = make_group("Actions")

        self._btn_generate = icon_button("▶  Generate Plot", tooltip="Render the selected plot")
        self._btn_generate.clicked.connect(self._on_generate)
        lay_run.addWidget(self._btn_generate)

        self._btn_run_geo = icon_button("⚙  Run Geological",
                                         tooltip="Classify resistivity → strat. logs")
        self._btn_run_geo.clicked.connect(self._on_run_geological)
        lay_run.addWidget(self._btn_run_geo)

        self._btn_run_hydro = icon_button("⚙  Run Hydrology",
                                           tooltip="Compute K, Sw, water table")
        self._btn_run_hydro.clicked.connect(self._on_run_hydro)
        lay_run.addWidget(self._btn_run_hydro)

        self._btn_run_mc = icon_button("🎲  Run Monte Carlo",
                                        tooltip="Uncertainty propagation")
        self._btn_run_mc.clicked.connect(self._on_run_mc)
        lay_run.addWidget(self._btn_run_mc)

        self._run_status = QLabel("")
        self._run_status.setObjectName("InfoLabel")
        self._run_status.setWordWrap(True)
        lay_run.addWidget(self._run_status)
        layout.addWidget(grp_run)

        # ── Export ────────────────────────────────────────────────────────
        grp_exp, lay_exp = make_group("Export")
        exp_row1 = QHBoxLayout()
        for label, slot in [("XYZ", self._export_xyz), ("LAS", self._export_las)]:
            b = QPushButton(label)
            b.setObjectName("FileListBtn")
            b.clicked.connect(slot)
            exp_row1.addWidget(b)
        exp_row2 = QHBoxLayout()
        for label, slot in [("CSV", self._export_csv), ("VTK", self._export_vtk)]:
            b = QPushButton(label)
            b.setObjectName("FileListBtn")
            b.clicked.connect(slot)
            exp_row2.addWidget(b)
        lay_exp.addLayout(exp_row1)
        lay_exp.addLayout(exp_row2)

        btn_export_fig = icon_button("⬆  Export Current Plot", "export",
                                      "Save the active plot to file")
        btn_export_fig.clicked.connect(self._export_current_figure)
        lay_exp.addWidget(btn_export_fig)
        layout.addWidget(grp_exp)

    # ── Right: plot studio ────────────────────────────────────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        # Studio splitter: tabs on top, gallery strip below
        studio_splitter = QSplitter(Qt.Orientation.Vertical)
        studio_splitter.setHandleWidth(4)

        # ── Tab widget (plot tabs) ─────────────────────────────────────────
        self._tab_widget = QTabWidget()
        self._tab_widget.setObjectName("PlotStudio")
        self._tab_widget.setTabsClosable(True)
        self._tab_widget.setMovable(True)
        self._tab_widget.setDocumentMode(True)
        self._tab_widget.tabCloseRequested.connect(self._close_tab)
        self._tab_widget.tabBarDoubleClicked.connect(self._rename_tab)
        self._tab_widget.setTabBarAutoHide(False)

        # "+" button to add a new blank tab
        add_btn = QToolButton()
        add_btn.setText("+")
        add_btn.setToolTip("Generate another plot")
        add_btn.clicked.connect(self._on_generate)
        self._tab_widget.setCornerWidget(add_btn, Qt.Corner.TopRightCorner)

        # Placeholder when no plots exist
        self._placeholder = QWidget()
        ph_lay = QVBoxLayout(self._placeholder)
        ph_lbl = QLabel(
            "▶  Select a workflow category\n"
            "    choose a plot and click  Generate\n\n"
            "    Each result opens as a new tab."
        )
        ph_lbl.setObjectName("SurveyPlaceholder")
        ph_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        ph_lay.addWidget(ph_lbl)
        self._tab_widget.addTab(self._placeholder, "  Welcome  ")
        self._tab_widget.tabBar().setTabButton(0, QTabBar.ButtonPosition.RightSide, None)

        studio_splitter.addWidget(self._tab_widget)

        # ── Gallery strip ──────────────────────────────────────────────────
        gallery_frame = QWidget()
        gallery_frame.setObjectName("GalleryFrame")
        gallery_frame.setFixedHeight(84)
        gallery_lay = QVBoxLayout(gallery_frame)
        gallery_lay.setContentsMargins(0, 0, 0, 0)
        gallery_lay.setSpacing(0)

        gallery_header = QLabel("  Plot Gallery")
        gallery_header.setObjectName("CanvasLabel")
        gallery_header.setFixedHeight(18)
        gallery_lay.addWidget(gallery_header)

        self._gallery = QListWidget()
        self._gallery.setObjectName("GalleryList")
        self._gallery.setFlow(QListWidget.Flow.LeftToRight)
        self._gallery.setViewMode(QListWidget.ViewMode.IconMode)
        self._gallery.setIconSize(QSize(80, 56))
        self._gallery.setSpacing(3)
        self._gallery.setResizeMode(QListWidget.ResizeMode.Adjust)
        self._gallery.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self._gallery.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._gallery.itemClicked.connect(self._on_gallery_click)
        gallery_lay.addWidget(self._gallery)
        studio_splitter.addWidget(gallery_frame)

        studio_splitter.setSizes([720, 84])
        studio_splitter.setStretchFactor(0, 1)
        studio_splitter.setStretchFactor(1, 0)

        layout.addWidget(studio_splitter)

    # ── Category / plot navigation ────────────────────────────────────────

    def _populate_category_combo(self) -> None:
        icons = {
            "Setup & Model":       "📐",
            "Geological":          "🪨",
            "Hydrology":           "💧",
            "Field Constraints":   "📍",
            "EM Diagnostics":      "📡",
            "Uncertainty (MC)":    "📊",
            "Advanced Plots":      "🔬",
            "Fusion & Time-Lapse": "🔀",
            "Export":              "💾",
        }
        self._combo_category.blockSignals(True)
        for cat in CATEGORIES:
            icon = icons.get(cat, "•")
            self._combo_category.addItem(f"{icon}  {cat}")
        self._combo_category.blockSignals(False)

    def _on_category_changed(self, row: int) -> None:
        if row < 0 or row >= len(CATEGORIES):
            return
        cat = CATEGORIES[row]
        plots = WORKFLOW_CATALOGUE[cat]
        self._combo_plot.blockSignals(True)
        self._combo_plot.clear()
        for label, _, _desc in plots:
            self._combo_plot.addItem(label)
        self._combo_plot.blockSignals(False)
        self._combo_plot.setCurrentIndex(0)
        self._update_plot_desc(row, 0)
        self._rebuild_param_form(cat, 0)
        self._combo_plot.currentIndexChanged.connect(
            lambda idx, r=row: self._on_plot_changed(r, idx)
        )

    def _on_plot_changed(self, cat_row: int, plot_idx: int) -> None:
        self._update_plot_desc(cat_row, plot_idx)
        cat = CATEGORIES[cat_row]
        self._rebuild_param_form(cat, plot_idx)

    def _update_plot_desc(self, cat_row: int, plot_idx: int) -> None:
        try:
            cat = CATEGORIES[cat_row]
            label, fn, desc = WORKFLOW_CATALOGUE[cat][plot_idx]
            self._plot_desc.setText(
                f"<small style='color:#888'>{desc}</small>"
            )
        except (IndexError, KeyError):
            self._plot_desc.setText("")

    def _rebuild_param_form(self, cat: str, plot_idx: int) -> None:
        """Build dynamic parameter widgets for the selected plot."""
        while self._param_form.rowCount():
            self._param_form.removeRow(0)
        self._param_widgets = {}

        # ── Category-specific global params ───────────────────────────────
        if cat == "Geological":
            w = QComboBox()
            stations = self._station_names()
            w.addItems(["(all)"] + stations)
            self._param_widgets["station"] = w
            self._param_form.addRow("Station:", w)

        elif cat == "Hydrology":
            for name, label, default, lo, hi, step in [
                ("m",     "Archie m",    2.0, 1.0, 3.5, 0.1),
                ("n",     "Archie n",    2.0, 1.0, 3.5, 0.1),
                ("rho_w", "ρ_w (Ω·m)",  10.0, 0.1, 1000.0, 0.1),
                ("phi",   "Porosity",   0.35, 0.01, 0.80, 0.01),
            ]:
                w = QDoubleSpinBox()
                w.setRange(lo, hi)
                w.setSingleStep(step)
                w.setDecimals(2)
                w.setValue(default)
                self._param_widgets[name] = w
                self._param_form.addRow(label + ":", w)

        elif cat == "Uncertainty (MC)":
            n_w = QSpinBox()
            n_w.setRange(50, 2000)
            n_w.setSingleStep(50)
            n_w.setValue(300)
            self._param_widgets["n_samples"] = n_w
            self._param_form.addRow("Samples:", n_w)

        elif cat == "Export":
            pass   # export uses file dialogs, no inline params

    # ── Generate plot ─────────────────────────────────────────────────────

    def _on_generate(self) -> None:
        cat_row  = self._combo_category.currentIndex()
        plot_idx = self._combo_plot.currentIndex()
        if cat_row < 0 or plot_idx < 0:
            return
        cat = CATEGORIES[cat_row]
        label, method_name, _desc = WORKFLOW_CATALOGUE[cat][plot_idx]

        # Collect params
        kwargs = {}
        for key, widget in getattr(self, "_param_widgets", {}).items():
            if isinstance(widget, (QSpinBox, QDoubleSpinBox)):
                kwargs[key] = widget.value()
            elif isinstance(widget, QComboBox):
                v = widget.currentText()
                kwargs[key] = "" if v == "(all)" else v
            else:
                kwargs[key] = widget.text() if hasattr(widget, "text") else None

        self._btn_generate.setEnabled(False)
        self._run_status.setText(f"Generating {label}…")

        self._plot_worker = _GeneratePlotWorker(
            self._ctrl, method_name, label, kwargs, parent=self
        )
        self._plot_worker.finished.connect(self._on_plot_ready)
        self._plot_worker.failed.connect(self._on_plot_failed)
        self._plot_worker.start()

    def _on_plot_ready(self, fig, label: str) -> None:
        self._add_plot_tab(fig, label)
        self._run_status.setText(f"✓  {label}")
        self._btn_generate.setEnabled(True)

    def _on_plot_failed(self, msg: str) -> None:
        self._run_status.setText(f"✕  {msg}")
        self._btn_generate.setEnabled(True)

    def _add_plot_tab(self, fig, title: str) -> None:
        """Add fig as a new tab and a thumbnail to the gallery."""
        # Remove placeholder tab if still there
        if self._tab_widget.count() == 1 and self._tab_widget.widget(0) is self._placeholder:
            self._tab_widget.clear()

        canvas = MplCanvas(self, toolbar=True)
        canvas.figure = fig
        canvas.draw()

        idx = self._tab_widget.addTab(canvas, f"  {title}  ")
        self._tab_widget.setCurrentIndex(idx)

        # Thumbnail
        thumb = self._make_thumbnail(fig, title)
        self._gallery.addItem(thumb)

        # Store (tab_index → gallery_index mapping via item data)
        thumb.setData(Qt.ItemDataRole.UserRole, idx)

    def _make_thumbnail(self, fig, title: str) -> QListWidgetItem:
        """Render a 80×56 px thumbnail from the figure."""
        try:
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=30,
                        bbox_inches="tight", facecolor=fig.get_facecolor())
            buf.seek(0)
            pix = QPixmap()
            pix.loadFromData(buf.read())
            icon = QIcon(pix.scaled(80, 56, Qt.AspectRatioMode.KeepAspectRatio,
                                    Qt.TransformationMode.SmoothTransformation))
        except Exception:
            icon = QIcon()
        item = QListWidgetItem(icon, "")
        item.setToolTip(title)
        return item

    # ── Tab management ────────────────────────────────────────────────────

    def _close_tab(self, idx: int) -> None:
        widget = self._tab_widget.widget(idx)
        if widget is self._placeholder:
            return
        self._tab_widget.removeTab(idx)
        # Remove gallery item at same index (best-effort)
        if idx < self._gallery.count():
            self._gallery.takeItem(idx)
        if self._tab_widget.count() == 0:
            self._tab_widget.addTab(self._placeholder, "  Welcome  ")
            self._tab_widget.tabBar().setTabButton(
                0, QTabBar.ButtonPosition.RightSide, None
            )

    def _rename_tab(self, idx: int) -> None:
        if idx < 0 or self._tab_widget.widget(idx) is self._placeholder:
            return
        from PySide6.QtWidgets import QInputDialog
        current = self._tab_widget.tabText(idx).strip()
        name, ok = QInputDialog.getText(self, "Rename Tab", "New name:", text=current)
        if ok and name:
            self._tab_widget.setTabText(idx, f"  {name}  ")

    def _on_gallery_click(self, item: QListWidgetItem) -> None:
        idx = item.data(Qt.ItemDataRole.UserRole)
        if idx is not None and 0 <= idx < self._tab_widget.count():
            self._tab_widget.setCurrentIndex(idx)

    # ── Action runners ────────────────────────────────────────────────────

    def _on_run_geological(self) -> None:
        self._run_status.setText("Running geological classification…")
        self._btn_run_geo.setEnabled(False)
        def _fn():
            return self._ctrl.run_geological()
        self._worker = _RunWorker(_fn, parent=self)
        self._worker.finished.connect(self._on_run_finished)
        self._worker.finished.connect(lambda _: self._btn_run_geo.setEnabled(True))
        self._worker.start()

    def _on_run_hydro(self) -> None:
        # Push petro config from UI
        kw = {}
        for k, w in getattr(self, "_param_widgets", {}).items():
            if isinstance(w, QDoubleSpinBox):
                kw[k] = w.value()
        if kw:
            self._ctrl.set_petro_config(**kw)
        self._run_status.setText("Running hydrology estimation…")
        self._btn_run_hydro.setEnabled(False)
        self._worker = _RunWorker(self._ctrl.run_hydro, parent=self)
        self._worker.finished.connect(self._on_run_finished)
        self._worker.finished.connect(lambda _: self._btn_run_hydro.setEnabled(True))
        self._worker.start()

    def _on_run_mc(self) -> None:
        n = self._param_widgets.get("n_samples")
        n_samples = n.value() if n else 200
        self._run_status.setText(f"Running MC ({n_samples} samples)…")
        self._btn_run_mc.setEnabled(False)
        def _fn():
            return self._ctrl.run_monte_carlo(n_samples=n_samples)
        self._worker = _RunWorker(_fn, parent=self)
        self._worker.finished.connect(self._on_run_finished)
        self._worker.finished.connect(lambda _: self._btn_run_mc.setEnabled(True))
        self._worker.start()

    def _on_run_finished(self, msg: str) -> None:
        self._run_status.setText(msg)
        self._update_status_card()

    # ── Model loading ─────────────────────────────────────────────────────

    def _load_from_inversion(self) -> None:
        parent = self.parent()
        if parent is not None and hasattr(parent, "_inversion_win"):
            inv_win = parent._inversion_win
            model = getattr(inv_win, "_result_model", None)
            if model is not None:
                self._ctrl.set_model(model)
                self._update_status_card()
                self._run_status.setText("Model loaded from Inversion window.")
                return
        self._run_status.setText("No model available from Inversion window.")

    def _load_occam2d(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select Occam2D result directory"
        )
        if not path:
            return
        try:
            self._ctrl.set_model_from_occam2d(path)
            self._update_status_card()
            self._run_status.setText("Occam2D model loaded.")
        except Exception as exc:
            self._run_status.setText(f"Load failed: {exc}")

    # ── Borehole management ───────────────────────────────────────────────

    def _add_borehole_csv(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load borehole CSV", "", "CSV files (*.csv);;All files (*)"
        )
        if path:
            try:
                name = self._ctrl.add_borehole_csv(path)
                self._bh_list.addItem(name)
            except Exception as exc:
                self._run_status.setText(f"Borehole load failed: {exc}")

    def _add_borehole_las(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load borehole LAS", "", "LAS files (*.las *.LAS);;All files (*)"
        )
        if path:
            try:
                name = self._ctrl.add_borehole_las(path)
                self._bh_list.addItem(name)
            except Exception as exc:
                self._run_status.setText(f"LAS load failed: {exc}")

    def _remove_borehole(self) -> None:
        item = self._bh_list.currentItem()
        if item:
            self._ctrl.remove_borehole(item.text())
            self._bh_list.takeItem(self._bh_list.row(item))

    # ── Export ────────────────────────────────────────────────────────────

    def _export_xyz(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Oasis Montaj XYZ", "", "XYZ (*.xyz);;All files (*)"
        )
        if path:
            self._run_status.setText(self._ctrl.export_xyz(path))

    def _export_las(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Export LAS 2.0", "", "LAS (*.las);;All files (*)"
        )
        if path:
            station = ""
            w = getattr(self, "_param_widgets", {}).get("station")
            if w and isinstance(w, QComboBox):
                station = w.currentText()
            self._run_status.setText(self._ctrl.export_las(path, station=station))

    def _export_csv(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Export CSV table", "", "CSV (*.csv);;All files (*)"
        )
        if path:
            self._run_status.setText(self._ctrl.export_csv(path))

    def _export_vtk(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Export VTK", "", "VTK (*.vtk);;All files (*)"
        )
        if path:
            self._run_status.setText(self._ctrl.export_vtk(path))

    def _export_current_figure(self) -> None:
        canvas = self._tab_widget.currentWidget()
        if not isinstance(canvas, MplCanvas):
            self._run_status.setText("No plot to export.")
            return
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog
        ExportDialog(figure=canvas.figure, parent=self).exec()

    # ── Public API ────────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        super().set_sites(sites)
        self._ctrl.set_sites(sites)
        self._update_status_card()

    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._ctrl.dark = dark

    # ── Helpers ───────────────────────────────────────────────────────────

    def _update_status_card(self) -> None:
        st = self._ctrl.model_status
        if not st["loaded"]:
            sites = self._ctrl.state.sites
            n_st  = ""
            if sites is not None:
                try:
                    from pycsamt.emtools._core import _iter_items
                    n_st = f"  {len(list(_iter_items(sites)))} stations loaded\n"
                except Exception:
                    pass
            self._status_card.setText(f"{n_st}No resistivity model\nLoad from inversion or file")
            return
        depth = st["depth_max"]
        prof  = st["profile_m"]
        rms   = st["rms"]
        txt   = (
            f"✓ {st['method']}  ·  {st['n_x']} × {st['n_z']}\n"
            f"Depth: {depth:.0f} m   Profile: {(prof or 0)/1e3:.1f} km\n"
            f"RMS: {rms or '—'}"
        )
        self._status_card.setText(txt)

    def _station_names(self) -> List[str]:
        model = self._ctrl.state.model
        if model and hasattr(model, "station_names") and model.station_names:
            return list(model.station_names)
        sites = self._ctrl.state.sites
        if sites is not None:
            try:
                from pycsamt.emtools._core import _iter_items, _name
                return [_name(ed, i) for i, ed in enumerate(_iter_items(sites))]
            except Exception:
                pass
        return []
