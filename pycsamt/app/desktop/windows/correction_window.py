# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
CorrectionWindow — interactive data-conditioning panel for AMT/CSAMT surveys.

Left panel  — category selector, correction chooser, dynamic parameter form,
              Preview/Apply actions, correction stack with undo/remove,
              Commit-to-Main and Revert-to-Raw controls.

Right panel — dual MplCanvas (Before / After) with view-mode switcher:
              Before/After  → two stacked canvases
              Overlay       → both datasets on one canvas (dashed=before, solid=after)
              Diff          → relative change section (ΔΩ/Ω %)

The correction stack is non-destructive: raw Sites are never modified.
``corrections_committed`` signal carries the final corrected Sites so
MainWindow can replace the global dataset.

``corrections_reverted`` notifies MainWindow that no corrections are active.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QMenu,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QRadioButton,
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

# Small symbol prefix per category for the dropdown — purely cosmetic
_CAT_ICON = {
    "Static Shift": "⇅",
    "Noise Removal": "∿",
    "Source Effects": "⊕",
    "Tensor Rotation": "↻",
    "Coordinates": "⊙",
    "Stratagem": "✦",
}

from pycsamt.app.desktop.controllers.correction_controller import (
    CATALOGUE,
    CATEGORIES,
    COORD_CATEGORIES,
    ROTATION_CATEGORIES,
    STATIC_SHIFT_CATEGORIES,
    STRATAGEM_CATEGORIES,
    CorrectionController,
    ParamSpec,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    icon_button,
    make_group,
)


class CorrectionWindow(PanelWindow):
    """
    Floating data-correction panel.

    Signals
    -------
    corrections_committed(object)
        Emitted when the user clicks "Commit to Main".
        Payload is the final corrected Sites object.
    corrections_reverted()
        Emitted when the user clicks "Revert to Raw".
    """

    corrections_committed = Signal(object)
    corrections_reverted = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(
            title="Data Corrections",
            session_key="correction_window",
            params_width=290,
            icon_name="sites-correction",
            parent=parent,
        )
        self.resize(1300, 780)
        self._ctrl = CorrectionController()
        self._preview_sites = None  # result of last Preview (not in stack)

        self._populate_category_combo()
        self._on_category_changed(
            0
        )  # populate _combo_correction on first open
        self._refresh_all()

    # ── Left params panel ─────────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:

        # ── 1. Correction Type ─────────────────────────────────────────
        grp_cat, lay_cat = make_group("Correction Type")

        # Category — single dropdown (one click, select, done)
        cat_lbl = QLabel("Category")
        cat_lbl.setObjectName("FieldLabel")
        lay_cat.addWidget(cat_lbl)
        self._combo_category = QComboBox()
        self._combo_category.setObjectName("CategoryCombo")
        self._combo_category.setToolTip("Select the correction family")
        self._combo_category.currentIndexChanged.connect(
            self._on_category_changed
        )
        lay_cat.addWidget(self._combo_category)

        # Sub-correction within the chosen category
        corr_lbl = QLabel("Correction")
        corr_lbl.setObjectName("FieldLabel")
        lay_cat.addWidget(corr_lbl)
        self._combo_correction = QComboBox()
        self._combo_correction.setObjectName("CorrectionCombo")
        self._combo_correction.currentIndexChanged.connect(
            self._on_correction_changed
        )
        lay_cat.addWidget(self._combo_correction)

        # Description chip — subtle left-accent border via stylesheet
        self._desc_lbl = QLabel("")
        self._desc_lbl.setWordWrap(True)
        self._desc_lbl.setObjectName("DescChip")
        self._desc_lbl.setAlignment(Qt.AlignmentFlag.AlignTop)
        self._desc_lbl.setVisible(False)
        lay_cat.addWidget(self._desc_lbl)

        layout.addWidget(grp_cat)

        # ── 2. Stratagem Source (hidden until Stratagem is selected) ───
        self._grp_strat, lay_strat = make_group("Stratagem Source")
        self._grp_strat.setVisible(False)

        self._radio_use_current = QRadioButton("Use currently loaded data")
        self._radio_load_dir = QRadioButton("Load from EDI directory")
        self._radio_use_current.setChecked(True)
        self._strat_radio_grp = QButtonGroup(self)
        self._strat_radio_grp.addButton(self._radio_use_current, 0)
        self._strat_radio_grp.addButton(self._radio_load_dir, 1)
        lay_strat.addWidget(self._radio_use_current)
        lay_strat.addWidget(self._radio_load_dir)

        # Path row (enabled only when "Load from dir" radio is active)
        self._strat_dir_row = QWidget()
        dir_h = QHBoxLayout(self._strat_dir_row)
        dir_h.setContentsMargins(0, 2, 0, 0)
        dir_h.setSpacing(4)
        self._edi_dir_label = QLabel("(no directory selected)")
        self._edi_dir_label.setObjectName("InfoLabel")
        self._edi_dir_label.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        btn_browse = QPushButton("📂")
        btn_browse.setFixedSize(28, 26)
        btn_browse.setToolTip("Browse for EDI directory")
        btn_browse.clicked.connect(self._on_browse_edi_dir)
        dir_h.addWidget(self._edi_dir_label)
        dir_h.addWidget(btn_browse)
        lay_strat.addWidget(self._strat_dir_row)

        self._btn_load_edi = QPushButton("⬆  Load EDI Dir")
        self._btn_load_edi.setObjectName("CommitButton")
        self._btn_load_edi.setToolTip(
            "Load all EDI files from the selected directory via EDIBatch.\n"
            "Switches the correction pipeline to Stratagem mode."
        )
        self._btn_load_edi.clicked.connect(self._on_load_edi_dir)
        lay_strat.addWidget(self._btn_load_edi)

        self._edi_load_status = QLabel("")
        self._edi_load_status.setObjectName("InfoLabel")
        self._edi_load_status.setWordWrap(True)
        lay_strat.addWidget(self._edi_load_status)

        # Wire radio → enable/disable the dir row
        self._strat_radio_grp.idToggled.connect(self._on_strat_source_toggled)
        self._strat_dir_row.setEnabled(False)
        self._btn_load_edi.setEnabled(False)

        layout.addWidget(self._grp_strat)

        # ── 3. Parameters ──────────────────────────────────────────────
        self._grp_params, self._lay_params_gb = make_group("Parameters")
        self._param_form = QFormLayout()
        self._param_form.setSpacing(5)
        self._param_form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        self._lay_params_gb.addLayout(self._param_form)
        self._param_widgets: dict = {}

        self._no_params_lbl = QLabel("Nothing to configure")
        self._no_params_lbl.setObjectName("InfoLabel")
        self._no_params_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._no_params_lbl.setVisible(False)
        self._lay_params_gb.addWidget(self._no_params_lbl)

        layout.addWidget(self._grp_params)

        # ── 3b. Affected Stations (Static Shift only) ──────────────────
        self._grp_ss_affected, lay_ss = make_group("Affected Stations")
        self._grp_ss_affected.setVisible(False)

        ss_hint = QLabel(
            "Station names that carry static shift\n"
            "(comma or newline separated):"
        )
        ss_hint.setObjectName("InfoLabel")
        ss_hint.setWordWrap(True)
        lay_ss.addWidget(ss_hint)

        self._txt_ss_stations = QPlainTextEdit()
        self._txt_ss_stations.setPlaceholderText("e.g. S001, S002, S003")
        self._txt_ss_stations.setMaximumHeight(72)
        self._txt_ss_stations.setObjectName("SSStationsEdit")
        lay_ss.addWidget(self._txt_ss_stations)

        layout.addWidget(self._grp_ss_affected)

        # ── 4. Actions — flat row, no group box border ─────────────────
        act_row = QHBoxLayout()
        act_row.setSpacing(6)
        self._btn_preview = QPushButton("▶  Preview")
        self._btn_preview.setToolTip("Compute result without saving to stack")
        self._btn_apply = QPushButton("✓  Apply")
        self._btn_apply.setObjectName("CommitButton")
        self._btn_apply.setToolTip("Apply correction and push to stack")
        self._btn_preview.clicked.connect(self._on_preview)
        self._btn_apply.clicked.connect(self._on_apply)
        act_row.addWidget(self._btn_preview)
        act_row.addWidget(self._btn_apply)
        layout.addLayout(act_row)

        self._action_status = QLabel("")
        self._action_status.setObjectName("InfoLabel")
        self._action_status.setWordWrap(True)
        self._action_status.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self._action_status)

        # ── 5. Applied corrections stack ───────────────────────────────
        self._grp_stk, lay_stk = make_group("Applied")
        self._stack_list = QListWidget()
        self._stack_list.setObjectName("StackList")
        self._stack_list.setMaximumHeight(120)
        self._stack_list.setContextMenuPolicy(
            Qt.ContextMenuPolicy.CustomContextMenu
        )
        self._stack_list.customContextMenuRequested.connect(
            self._on_stack_ctx_menu
        )
        lay_stk.addWidget(self._stack_list)

        undo_row = QHBoxLayout()
        undo_row.setSpacing(6)
        self._btn_undo = QPushButton("↩  Undo")
        self._btn_undo.setObjectName("FileListBtn")
        self._btn_undo.clicked.connect(self._on_undo)
        self._btn_clr = QPushButton("✗  Clear")
        self._btn_clr.setObjectName("FileListBtn")
        self._btn_clr.clicked.connect(self._on_clear_stack)
        undo_row.addWidget(self._btn_undo)
        undo_row.addWidget(self._btn_clr)
        lay_stk.addLayout(undo_row)

        self._stack_status = QLabel("No corrections applied")
        self._stack_status.setObjectName("InfoLabel")
        self._stack_status.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay_stk.addWidget(self._stack_status)
        layout.addWidget(self._grp_stk)

        # ── 6. Output ──────────────────────────────────────────────────
        grp_commit, lay_commit = make_group("Output")

        self._btn_commit = QPushButton("⬆  Commit to Main")
        self._btn_commit.setObjectName("CommitButton")
        self._btn_commit.setToolTip(
            "Replace the loaded dataset with the current corrected Sites.\n"
            "All panels (Profile, Map, QC, …) will update."
        )
        self._btn_commit.clicked.connect(self._on_commit)

        self._btn_revert = QPushButton("↺  Revert to Raw")
        self._btn_revert.setObjectName("RevertButton")
        self._btn_revert.setToolTip(
            "Clear all corrections in this panel (does not affect main data)."
        )
        self._btn_revert.clicked.connect(self._on_revert)

        self._btn_export_strat = QPushButton("💾  Export Stratagem EDIs…")
        self._btn_export_strat.setObjectName("FileListBtn")
        self._btn_export_strat.setToolTip(
            "Write Stratagem-corrected EDI files to a directory."
        )
        self._btn_export_strat.clicked.connect(self._on_strat_export)
        self._btn_export_strat.setVisible(False)

        lay_commit.addWidget(self._btn_commit)
        lay_commit.addWidget(self._btn_revert)
        lay_commit.addWidget(self._btn_export_strat)
        layout.addWidget(grp_commit)

    # ── Right content panel ───────────────────────────────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        # ── Toolbar ───────────────────────────────────────────────────
        bar = QHBoxLayout()
        bar.setContentsMargins(8, 4, 8, 0)
        bar.setSpacing(8)

        bar.addWidget(QLabel("View:"))
        self._combo_mode = QComboBox()
        self._combo_mode.addItems(
            ["Before / After", "Overlay", "Diff", "2D Section"]
        )
        self._combo_mode.setFixedWidth(150)
        self._combo_mode.currentIndexChanged.connect(self._on_mode_changed)
        bar.addWidget(self._combo_mode)

        # Sub-view for Coordinates category
        self._combo_coord_view = QComboBox()
        self._combo_coord_view.addItems(["Position map", "Elevation profile"])
        self._combo_coord_view.setFixedWidth(140)
        self._combo_coord_view.currentIndexChanged.connect(self._refresh_plots)
        self._combo_coord_view.setVisible(False)
        bar.addWidget(self._combo_coord_view)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.VLine)
        sep.setObjectName("Separator")
        bar.addWidget(sep)

        self._btn_export = icon_button(
            "⬆  Export…", "export", "Export current figure"
        )
        self._btn_export.setFixedWidth(110)
        self._btn_export.clicked.connect(self._on_export)
        bar.addWidget(self._btn_export)

        bar.addStretch()
        self._view_status = QLabel("")
        self._view_status.setObjectName("InfoLabel")
        bar.addWidget(self._view_status)

        bar_w = QWidget()
        bar_w.setLayout(bar)
        layout.addWidget(bar_w)

        # ── Stacked view: page 0 = split, page 1 = single (overlay/diff) ──
        self._view_stack = QStackedWidget()
        layout.addWidget(self._view_stack)

        # ── Page 0: vertical split before / after ─────────────────────
        split_page = QWidget()
        split_v = QVBoxLayout(split_page)
        split_v.setContentsMargins(0, 0, 0, 0)
        split_v.setSpacing(0)

        self._canvas_splitter = QSplitter(Qt.Orientation.Vertical)
        self._canvas_splitter.setHandleWidth(4)

        def _canvas_pane(title: str) -> tuple[QWidget, MplCanvas]:
            pane = QWidget()
            pane.setObjectName("CanvasPane")
            v = QVBoxLayout(pane)
            v.setContentsMargins(0, 0, 0, 0)
            v.setSpacing(0)
            lbl = QLabel(f"  {title}")
            lbl.setObjectName("CanvasLabel")
            lbl.setFixedHeight(22)
            v.addWidget(lbl)
            canvas = MplCanvas(pane, toolbar=False)
            v.addWidget(canvas)
            return pane, canvas

        before_pane, self._canvas_before = _canvas_pane("Before")
        after_pane, self._canvas_after = _canvas_pane("After / Preview")

        self._canvas_splitter.addWidget(before_pane)
        self._canvas_splitter.addWidget(after_pane)
        self._canvas_splitter.setSizes([360, 360])
        split_v.addWidget(self._canvas_splitter)
        self._view_stack.addWidget(split_page)

        # ── Page 1: single canvas for overlay / diff ──────────────────
        single_page = QWidget()
        single_v = QVBoxLayout(single_page)
        single_v.setContentsMargins(0, 0, 0, 0)
        single_v.setSpacing(0)
        self._canvas_single = MplCanvas(single_page, toolbar=True)
        single_v.addWidget(self._canvas_single)
        self._view_stack.addWidget(single_page)

        # ── Page 2: Stratagem Studio ───────────────────────────────────
        strat_page = QWidget()
        strat_v = QVBoxLayout(strat_page)
        strat_v.setContentsMargins(0, 0, 0, 0)
        strat_v.setSpacing(0)

        # Progress bar (hidden by default; visible while a worker runs)
        self._strat_progress = QProgressBar()
        self._strat_progress.setRange(0, 0)  # indeterminate pulse
        self._strat_progress.setFixedHeight(6)
        self._strat_progress.setVisible(False)
        strat_v.addWidget(self._strat_progress)

        self._strat_tabs = QTabWidget()
        self._strat_tabs.setDocumentMode(True)
        strat_v.addWidget(self._strat_tabs)

        # Tab 0 — QC Report (horizontal bar charts)
        qc_page = QWidget()
        qc_v = QVBoxLayout(qc_page)
        qc_v.setContentsMargins(0, 0, 0, 0)
        self._canvas_strat_qc = MplCanvas(qc_page, toolbar=True)
        qc_v.addWidget(self._canvas_strat_qc)
        self._strat_tabs.addTab(qc_page, "QC Report")

        # Tab 1 — Before / After impedance
        ba_page = QWidget()
        ba_v = QVBoxLayout(ba_page)
        ba_v.setContentsMargins(0, 0, 0, 0)
        ba_splitter = QSplitter(Qt.Orientation.Vertical)
        ba_splitter.setHandleWidth(4)

        def _strat_pane(title):
            p = QWidget()
            p.setObjectName("CanvasPane")
            v = QVBoxLayout(p)
            v.setContentsMargins(0, 0, 0, 0)
            v.setSpacing(0)
            lbl = QLabel(f"  {title}")
            lbl.setObjectName("CanvasLabel")
            lbl.setFixedHeight(22)
            v.addWidget(lbl)
            c = MplCanvas(p, toolbar=False)
            v.addWidget(c)
            return p, c

        ba_before_pane, self._canvas_strat_before = _strat_pane("Before (raw)")
        ba_after_pane, self._canvas_strat_after = _strat_pane(
            "After (corrected)"
        )
        ba_splitter.addWidget(ba_before_pane)
        ba_splitter.addWidget(ba_after_pane)
        ba_splitter.setSizes([350, 350])
        ba_v.addWidget(ba_splitter)
        self._strat_tabs.addTab(ba_page, "Before / After")

        # Tab 2 — Static-shift factors bar chart
        ss_page = QWidget()
        ss_v = QVBoxLayout(ss_page)
        ss_v.setContentsMargins(0, 0, 0, 0)
        self._canvas_strat_ss = MplCanvas(ss_page, toolbar=True)
        ss_v.addWidget(self._canvas_strat_ss)
        self._strat_tabs.addTab(ss_page, "SS Factors")

        # Tab 3 — Station QC table (raw numbers)
        tbl_page = QWidget()
        tbl_v = QVBoxLayout(tbl_page)
        self._strat_qc_table = QTableWidget()
        self._strat_qc_table.setObjectName("StationTable")
        self._strat_qc_table.setAlternatingRowColors(True)
        self._strat_qc_table.horizontalHeader().setStretchLastSection(True)
        tbl_v.addWidget(self._strat_qc_table)
        self._strat_tabs.addTab(tbl_page, "QC Table")

        self._view_stack.addWidget(strat_page)

        # ── Page 3: 2-D ρ_a pseudosection (Static Shift) ─────────────
        pseudo_page = QWidget()
        pseudo_v = QVBoxLayout(pseudo_page)
        pseudo_v.setContentsMargins(0, 0, 0, 0)
        pseudo_v.setSpacing(0)
        self._canvas_pseudo = MplCanvas(pseudo_page, toolbar=True)
        pseudo_v.addWidget(self._canvas_pseudo)
        self._view_stack.addWidget(pseudo_page)

        self._view_stack.setCurrentIndex(0)

    # ── Populate category combo ────────────────────────────────────────

    def _populate_category_combo(self) -> None:
        self._combo_category.blockSignals(True)
        for cat in CATEGORIES:
            icon = _CAT_ICON.get(cat, "·")
            self._combo_category.addItem(f"{icon}  {cat}")
        self._combo_category.blockSignals(False)

    # ── Public API ────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        super().set_sites(sites)
        self._ctrl.set_raw_sites(sites)
        self._preview_sites = None
        self._refresh_all()

    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._ctrl.dark = dark
        self._refresh_plots()

    # ── Slots: left panel ─────────────────────────────────────────────

    def _on_category_changed(self, row: int) -> None:
        if row < 0 or row >= len(CATEGORIES):
            return
        cat = CATEGORIES[row]
        corrections = list(CATALOGUE[cat].keys())
        self._combo_correction.blockSignals(True)
        self._combo_correction.clear()
        self._combo_correction.addItems(corrections)
        self._combo_correction.blockSignals(False)
        self._combo_correction.setCurrentIndex(0)

        is_strat = cat in STRATAGEM_CATEGORIES
        is_ss = cat in STATIC_SHIFT_CATEGORIES

        # Show / hide category-specific side panels
        self._grp_strat.setVisible(is_strat)
        self._grp_ss_affected.setVisible(is_ss)
        self._combo_mode.setVisible(not is_strat)
        self._combo_coord_view.setVisible(
            self._is_coord_category() and not is_strat
        )
        self._btn_export_strat.setVisible(is_strat)

        if is_strat:
            self._view_stack.setCurrentIndex(2)
            self._refresh_strat_plots()
        else:
            # Auto-switch view mode
            if is_ss:
                self._combo_mode.blockSignals(True)
                self._combo_mode.setCurrentText("2D Section")
                self._combo_mode.blockSignals(False)
                self._view_stack.setCurrentIndex(3)
            else:
                mode = self._combo_mode.currentText()
                if mode == "2D Section":
                    # 2D Section requires impedance data — not for coord corrections
                    # or when switching away from Static Shift to another category
                    self._combo_mode.blockSignals(True)
                    self._combo_mode.setCurrentText("Before / After")
                    self._combo_mode.blockSignals(False)
                    self._view_stack.setCurrentIndex(0)
                elif mode == "Before / After":
                    self._view_stack.setCurrentIndex(0)
                else:
                    self._view_stack.setCurrentIndex(1)

        self._on_correction_changed(0)
        self._preview_sites = None
        if not is_strat:
            self._refresh_plots()

    def _on_correction_changed(self, idx: int) -> None:
        cat_row = self._combo_category.currentIndex()
        if cat_row < 0:
            return
        cat = CATEGORIES[cat_row]
        corrections = list(CATALOGUE[cat].keys())
        if idx < 0 or idx >= len(corrections):
            return
        label = corrections[idx]
        info = CATALOGUE[cat][label]
        desc = info.get("desc", "")
        self._desc_lbl.setVisible(bool(desc))
        self._desc_lbl.setText(
            f"<small style='color:#888'>{desc}</small>" if desc else ""
        )
        self._rebuild_param_form(info.get("params", []))
        # Auto-refresh view when the user picks a different correction sub-type
        if cat not in STRATAGEM_CATEGORIES:
            self._refresh_plots()

    def _rebuild_param_form(self, params: list) -> None:
        """Clear and recreate the dynamic parameter form."""
        while self._param_form.rowCount():
            self._param_form.removeRow(0)
        self._param_widgets.clear()

        for spec in params:
            widget = self._make_widget(spec)
            self._param_widgets[spec.name] = widget
            if spec.tip:
                widget.setToolTip(spec.tip)
            self._param_form.addRow(spec.label + ":", widget)

        self._no_params_lbl.setVisible(len(params) == 0)

    def _make_widget(self, spec: ParamSpec) -> QWidget:
        if spec.kind == "spin":
            w = QSpinBox()
            lo, hi, step = spec.opts
            w.setRange(lo, hi)
            w.setSingleStep(step)
            w.setValue(int(spec.default))
            return w
        if spec.kind == "dspin":
            w = QDoubleSpinBox()
            lo, hi, step = spec.opts
            decimals = max(
                0, -int(np.floor(np.log10(step))) if step < 1 else 1
            )
            w.setRange(lo, hi)
            w.setSingleStep(step)
            w.setDecimals(decimals)
            w.setValue(float(spec.default))
            return w
        if spec.kind == "combo":
            w = QComboBox()
            w.addItems(spec.opts)
            idx = (
                spec.opts.index(spec.default)
                if spec.default in spec.opts
                else 0
            )
            w.setCurrentIndex(idx)
            return w
        if spec.kind == "check":
            from PySide6.QtWidgets import QCheckBox

            w = QCheckBox()
            w.setChecked(bool(spec.default))
            return w
        # fallback: line edit
        from PySide6.QtWidgets import QLineEdit

        w = QLineEdit(str(spec.default))
        return w

    def _get_param_values(self) -> dict:
        from PySide6.QtWidgets import QLineEdit

        vals = {}
        for name, widget in self._param_widgets.items():
            if isinstance(widget, QSpinBox):
                vals[name] = widget.value()
            elif isinstance(widget, QDoubleSpinBox):
                vals[name] = widget.value()
            elif isinstance(widget, QComboBox):
                vals[name] = widget.currentText()
            elif isinstance(widget, QCheckBox):
                vals[name] = widget.isChecked()
            elif isinstance(widget, QLineEdit):
                try:
                    vals[name] = float(widget.text())
                except Exception:
                    vals[name] = widget.text()
        # Include affected-station names for Static Shift
        cat_row = self._combo_category.currentIndex()
        if 0 <= cat_row < len(CATEGORIES):
            if CATEGORIES[cat_row] in STATIC_SHIFT_CATEGORIES:
                vals["affected_stations"] = self._get_affected_stations()
        return vals

    def _get_affected_stations(self) -> list[str]:
        """Parse station names from the Affected Stations text field."""
        import re

        text = self._txt_ss_stations.toPlainText().strip()
        if not text:
            return []
        names = [n.strip() for n in re.split(r"[,;\n]+", text)]
        return [n for n in names if n]

    def _current_fn_label(self) -> tuple[str, str]:
        """Return (fn_name, human_label) for the currently selected correction."""
        cat_row = self._combo_category.currentIndex()
        if cat_row < 0:
            return "", ""
        cat = CATEGORIES[cat_row]
        corrections = list(CATALOGUE[cat].keys())
        cidx = self._combo_correction.currentIndex()
        if cidx < 0 or cidx >= len(corrections):
            return "", ""
        label = corrections[cidx]
        fn = CATALOGUE[cat][label]["fn"]
        return fn, label

    # ── Action slots ──────────────────────────────────────────────────

    def _on_preview(self) -> None:
        if not self._ctrl.has_data:
            self._action_status.setText("Load survey data first.")
            return
        fn_name, label = self._current_fn_label()
        if not fn_name:
            return
        kwargs = self._get_param_values()
        self._action_status.setText("Computing preview…")
        self._btn_preview.setEnabled(False)
        try:
            result = self._ctrl.preview(fn_name, kwargs)
            if result is not None:
                self._preview_sites = (
                    result  # may be DataFrame for coord corrections
                )
                self._action_status.setText(f"Preview: {label}")
                self._refresh_plots()
            else:
                self._action_status.setText(
                    "Preview failed — check parameters."
                )
        except Exception as exc:
            self._action_status.setText(f"Error: {exc}")
        finally:
            self._btn_preview.setEnabled(True)

    def _on_apply(self) -> None:
        if not self._ctrl.has_data:
            self._action_status.setText("Load survey data first.")
            return
        fn_name, label = self._current_fn_label()
        if not fn_name:
            return
        kwargs = self._get_param_values()
        params_str = ", ".join(f"{k}={v}" for k, v in kwargs.items())
        full_label = f"{label}  ({params_str})" if params_str else label
        self._action_status.setText(f"Applying {label}…")
        self._btn_apply.setEnabled(False)
        try:
            step = self._ctrl.apply(fn_name, kwargs, full_label)
            if step is not None:
                self._preview_sites = None
                self._action_status.setText(f"Applied: {label}")
                self._refresh_stack_list()
                # Stratagem page needs its own refresh; normal page otherwise
                if self._view_stack.currentIndex() == 2:
                    self._refresh_strat_plots()
                else:
                    self._refresh_plots()
            else:
                self._action_status.setText("Apply failed — check parameters.")
        except Exception as exc:
            self._action_status.setText(f"Error: {exc}")
        finally:
            self._btn_apply.setEnabled(True)

    def _on_undo(self) -> None:
        self._ctrl.undo_last()
        self._preview_sites = None
        self._refresh_stack_list()
        self._refresh_plots()

    def _on_clear_stack(self) -> None:
        self._ctrl.revert_all()
        self._preview_sites = None
        self._refresh_stack_list()
        self._refresh_plots()

    def _on_stack_ctx_menu(self, pos) -> None:
        item = self._stack_list.itemAt(pos)
        if item is None:
            return
        row = self._stack_list.row(item)
        menu = QMenu(self)
        act_remove = menu.addAction(f"Remove step {row + 1}")
        act_view = menu.addAction("Show 'After' for this step")
        chosen = menu.exec(self._stack_list.mapToGlobal(pos))
        if chosen == act_remove:
            self._ctrl.remove_step(row)
            self._preview_sites = None
            self._refresh_stack_list()
            self._refresh_plots()
        elif chosen == act_view:
            stack = self._ctrl.stack
            if row < len(stack):
                self._preview_sites = stack[row].sites_after
                self._refresh_plots()

    # ── Stratagem slots ───────────────────────────────────────────────

    def _on_strat_source_toggled(self, btn_id: int, checked: bool) -> None:
        """Enable/disable EDI dir row based on radio selection."""
        is_load = self._strat_radio_grp.checkedId() == 1
        self._strat_dir_row.setEnabled(is_load)
        self._btn_load_edi.setEnabled(is_load)

    def _on_browse_edi_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select EDI Directory", "", QFileDialog.Option.ShowDirsOnly
        )
        if path:
            self._edi_dir_label.setText(path)
            self._edi_dir_label.setToolTip(path)

    def _on_load_edi_dir(self) -> None:
        path = self._edi_dir_label.text()
        if not path or path == "(no directory selected)":
            self._edi_load_status.setText("⚠  Select a directory first.")
            return
        self._strat_progress.setVisible(True)
        self._btn_load_edi.setEnabled(False)
        self._edi_load_status.setText("Loading…")
        try:
            n = self._ctrl.load_edi_dir(path)
            self._edi_load_status.setText(
                f"✓  {n} station{'s' if n != 1 else ''} loaded."
            )
            self._refresh_all()
            # Switch to Stratagem category automatically
            strat_row = CATEGORIES.index("Stratagem")
            self._combo_category.setCurrentIndex(strat_row)
        except Exception as exc:
            self._edi_load_status.setText(f"✕  {exc}")
        finally:
            self._strat_progress.setVisible(False)
            self._btn_load_edi.setEnabled(True)

    def _refresh_strat_plots(self) -> None:
        """Refresh all three Stratagem canvases from controller state."""
        # QC bar chart
        fig_qc = self._canvas_strat_qc.figure
        try:
            self._ctrl.plot_strat_qc(fig_qc)
            self._canvas_strat_qc.draw()
        except Exception:
            pass
        # SS factors
        ax_ss = self._canvas_strat_ss.axes
        try:
            self._ctrl.plot_strat_ss_factors(ax_ss)
            self._canvas_strat_ss.draw()
        except Exception:
            pass
        # Before / After impedance curves
        try:
            raw_sites = self._ctrl.raw_sites
            curr_sites = self._ctrl.current_sites
            if raw_sites is not None:
                ax_b = self._canvas_strat_before.axes
                ax_b.cla()
                self._ctrl.plot_rho_curves(raw_sites, ax_b, "Before (raw)")
                self._canvas_strat_before.draw()
            if curr_sites is not None and curr_sites is not raw_sites:
                ax_a = self._canvas_strat_after.axes
                ax_a.cla()
                self._ctrl.plot_rho_curves(
                    curr_sites, ax_a, "After (corrected)"
                )
                self._canvas_strat_after.draw()
        except Exception:
            pass
        # QC data table
        self._populate_strat_qc_table()

    def _populate_strat_qc_table(self) -> None:
        df = self._ctrl._strat_qc_report
        tbl = self._strat_qc_table
        tbl.clearContents()
        if df is None or df.empty:
            tbl.setRowCount(0)
            tbl.setColumnCount(0)
            return
        cols = list(df.columns)
        tbl.setColumnCount(len(cols))
        tbl.setRowCount(len(df))
        tbl.setHorizontalHeaderLabels(cols)
        for r, row in df.iterrows():
            for c, col in enumerate(cols):
                val = row[col]
                item = QTableWidgetItem(
                    f"{val:.4f}" if isinstance(val, float) else str(val)
                )
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                # Flag rows in red
                if "flagged" in df.columns and bool(row.get("flagged", False)):
                    from PySide6.QtGui import QColor

                    item.setForeground(QColor("#f38ba8"))
                tbl.setItem(r, c, item)
        tbl.resizeColumnsToContents()

    def _on_strat_export(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self,
            "Export Stratagem EDIs to…",
            "",
            QFileDialog.Option.ShowDirsOnly,
        )
        if not path:
            return
        try:
            n = self._ctrl.export_stratagem(path)
            self._action_status.setText(
                f"✓  Exported {n} EDI file(s) to {path}"
            )
        except Exception as exc:
            self._action_status.setText(f"✕  Export failed: {exc}")

    # ── View-mode slot ────────────────────────────────────────────────

    def _on_mode_changed(self, _idx: int) -> None:
        mode = self._combo_mode.currentText()
        if mode == "Before / After":
            self._view_stack.setCurrentIndex(0)
        elif mode == "2D Section":
            self._view_stack.setCurrentIndex(3)
        else:
            self._view_stack.setCurrentIndex(1)
        self._refresh_plots()

    # ── Commit / Revert ───────────────────────────────────────────────

    def _on_commit(self) -> None:
        if not self._ctrl.has_corrections:
            self._view_status.setText("No corrections to commit.")
            return
        sites = self._ctrl.current_sites
        # Apply any coordinate corrections to the EDI headers (only at commit time)
        if self._ctrl.has_coord_corrections:
            from pycsamt.gis.coord_correction import (
                apply_coords_df_to_sites,
            )

            final_coords = self._ctrl.current_coords_df()
            if final_coords is not None:
                apply_coords_df_to_sites(sites, final_coords)
        self.corrections_committed.emit(sites)
        n = self._ctrl.n_steps
        self._view_status.setText(
            f"✓ Committed {n} correction{'s' if n != 1 else ''} to main data."
        )

    def _on_revert(self) -> None:
        self._ctrl.revert_all()
        self._preview_sites = None
        self._refresh_stack_list()
        self._refresh_plots()
        self.corrections_reverted.emit()
        self._view_status.setText("Reverted to raw data in this panel.")

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )

        mode = self._combo_mode.currentText()
        fig = (
            self._canvas_single.figure
            if mode != "Before / After"
            else self._canvas_after.figure
        )
        ExportDialog(figure=fig, parent=self).exec()

    # ── View-mode helpers ─────────────────────────────────────────────

    def _is_coord_category(self) -> bool:
        row = self._combo_category.currentIndex()
        if row < 0 or row >= len(CATEGORIES):
            return False
        return CATEGORIES[row] in COORD_CATEGORIES

    # ── Refresh helpers ───────────────────────────────────────────────

    def _refresh_all(self) -> None:
        self._refresh_stack_list()
        self._refresh_plots()
        enabled = self._ctrl.has_data
        self._btn_preview.setEnabled(enabled)
        self._btn_apply.setEnabled(enabled)
        self._btn_commit.setEnabled(enabled and self._ctrl.has_corrections)
        self._btn_revert.setEnabled(enabled and self._ctrl.has_corrections)

    def _refresh_stack_list(self) -> None:
        self._stack_list.clear()
        stack = self._ctrl.stack
        for i, step in enumerate(stack):
            self._stack_list.addItem(f"  {i + 1}.  {step.label}")
        n = len(stack)
        self._grp_stk.setTitle(f"Applied  ({n})" if n > 0 else "Applied")
        if n == 0:
            self._stack_status.setText("No corrections applied")
        else:
            self._stack_status.setText(
                f"{n} correction{'s' if n != 1 else ''} in stack"
            )
        has_corr = n > 0
        self._btn_undo.setEnabled(has_corr)
        self._btn_clr.setEnabled(has_corr)
        self._btn_commit.setEnabled(has_corr)
        self._btn_revert.setEnabled(has_corr)

    def _is_ss_category(self) -> bool:
        row = self._combo_category.currentIndex()
        if row < 0 or row >= len(CATEGORIES):
            return False
        return CATEGORIES[row] in STATIC_SHIFT_CATEGORIES

    def _refresh_plots(self) -> None:
        mode = self._combo_mode.currentText()
        is_coord = self._is_coord_category()
        _crow = self._combo_category.currentIndex()
        is_rotation = (
            CATEGORIES[_crow] in ROTATION_CATEGORIES
            if 0 <= _crow < len(CATEGORIES)
            else False
        )
        elev_view = (
            self._combo_coord_view.currentText() == "Elevation profile"
            if is_coord
            else False
        )

        # ── 2D pseudosection view (Static Shift only) ────────────────────────
        if mode == "2D Section":
            self._view_stack.setCurrentIndex(3)
            fig_ps = self._canvas_pseudo.figure
            if is_coord:
                # Coordinates category has no impedance data — show message
                fig_ps.clear()
                ax_msg = fig_ps.add_subplot(111)
                ax_msg.set_facecolor(
                    "#181825" if self._ctrl.dark else "#eff1f5"
                )
                fig_ps.patch.set_facecolor(
                    "#1e1e2e" if self._ctrl.dark else "#e6e9ef"
                )
                tc = "#a6adc8" if self._ctrl.dark else "#6c6f85"
                ax_msg.text(
                    0.5,
                    0.5,
                    "2D pseudosection requires impedance (Z) data.\n"
                    "Not applicable for coordinate corrections.\n\n"
                    "Switch to  Before / After  or  Overlay  to compare positions.",
                    transform=ax_msg.transAxes,
                    ha="center",
                    va="center",
                    fontsize=10,
                    color=tc,
                    multialignment="center",
                )
                ax_msg.axis("off")
                self._canvas_pseudo.draw()
                return
            current_data = (
                self._preview_sites
                if self._preview_sites is not None
                else self._ctrl.current_sites
            )
            affected = self._get_affected_stations()
            try:
                self._ctrl.plot_rho_pseudosection(
                    current_data,
                    fig_ps,
                    affected_stations=affected or None,
                    title="ρ_a Pseudosection  —  Static Shift",
                )
                try:
                    fig_ps.tight_layout(pad=1.0)
                except Exception:
                    pass
            except Exception:
                pass
            self._canvas_pseudo.draw()
            return

        # ── Determine before / after data sources ────────────────────────────
        if is_coord:
            # Coords: always use DataFrame snapshots (never EDI head objects)
            before_data = self._ctrl.raw_coords_df
            after_data = (
                self._preview_sites  # may be a DataFrame from preview()
                if self._preview_sites is not None
                else self._ctrl.current_coords_df()
            )
            after_label = (
                "Preview"
                if self._preview_sites is not None
                else f"After  ({self._ctrl.n_steps} step{'s' if self._ctrl.n_steps != 1 else ''})"
            )
        else:
            before_data = self._ctrl.raw_sites
            after_data = (
                self._preview_sites
                if self._preview_sites is not None
                else self._ctrl.current_sites
            )
            after_label = (
                "Preview"
                if self._preview_sites is not None
                else f"After  ({self._ctrl.n_steps} step{'s' if self._ctrl.n_steps != 1 else ''})"
            )

        def _draw(fn, *args):
            try:
                fn(*args)
            except Exception:
                pass

        # ── Before / After split view ─────────────────────────────────────────
        if mode == "Before / After":
            self._view_stack.setCurrentIndex(0)

            self._canvas_before.figure.clear()
            ax_b = self._canvas_before.figure.add_subplot(111)
            if is_coord and elev_view:
                _draw(
                    self._ctrl.plot_station_elevation,
                    before_data,
                    ax_b,
                    "Before — elevation",
                )
            elif is_coord:
                _draw(
                    self._ctrl.plot_station_map,
                    before_data,
                    ax_b,
                    "Before — positions",
                )
            else:
                _draw(
                    self._ctrl.plot_rho_curves,
                    before_data,
                    ax_b,
                    "Before (raw)",
                )
            try:
                self._canvas_before.figure.tight_layout(pad=1.0)
            except Exception:
                pass
            self._canvas_before.draw()

            self._canvas_after.figure.clear()
            ax_a = self._canvas_after.figure.add_subplot(111)
            if is_coord and elev_view:
                _draw(
                    self._ctrl.plot_station_elevation,
                    after_data,
                    ax_a,
                    after_label,
                )
            elif is_coord:
                _draw(
                    self._ctrl.plot_station_map, after_data, ax_a, after_label
                )
            else:
                _draw(
                    self._ctrl.plot_rho_curves, after_data, ax_a, after_label
                )
            try:
                self._canvas_after.figure.tight_layout(pad=1.0)
            except Exception:
                pass
            self._canvas_after.draw()

        # ── Overlay / Rose / Diff single-canvas view ──────────────────────────
        else:
            self._view_stack.setCurrentIndex(1)
            fig = self._canvas_single.figure
            fig.clear()

            if mode == "Overlay":
                if is_coord and elev_view:
                    # True overlay: both profiles on one axis, distinct colours
                    ax = fig.add_subplot(111)
                    _draw(
                        self._ctrl.plot_station_elevation_overlay,
                        before_data,
                        after_data,
                        ax,
                    )
                elif is_coord:
                    ax = fig.add_subplot(111)
                    _draw(
                        self._ctrl.plot_station_map_overlay,
                        before_data,
                        after_data,
                        ax,
                    )
                elif is_rotation:
                    # Rose diagram for tensor rotation
                    _draw(
                        self._ctrl.plot_rotation_rose,
                        before_data,
                        after_data,
                        fig,
                    )
                else:
                    ax = fig.add_subplot(111)
                    _draw(self._ctrl.plot_overlay, before_data, after_data, ax)

            else:  # Diff
                if is_coord:
                    ax = fig.add_subplot(111)
                    _draw(
                        self._ctrl.plot_displacement_diff,
                        before_data,
                        after_data,
                        ax,
                    )
                elif is_rotation:
                    _draw(
                        self._ctrl.plot_rotation_rose,
                        before_data,
                        after_data,
                        fig,
                    )
                else:
                    ax = fig.add_subplot(111)
                    _draw(self._ctrl.plot_diff, before_data, after_data, ax)

            try:
                fig.tight_layout(pad=1.2)
            except Exception:
                pass
            self._canvas_single.draw()


# ── Needed for type hint in _make_widget ──────────────────────────────────────
import numpy as np
