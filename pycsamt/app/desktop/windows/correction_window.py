# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
CorrectionWindow — interactive data-conditioning panel for AMT/CSAMT surveys.

Left panel  — category selector, correction chooser, dynamic parameter form,
              Preview/Apply actions, correction stack with undo/remove,
              Commit-to-Main and Revert-to-Raw controls.

Right panel — one comparison canvas driven by two orthogonal choices:
              Compare  → Before / After · Overlay · Diff
              Display  → Curves (1-D) · Pseudosection (2-D) · Strike rose
                         (rotation) · Position map / Elevation profile
                         (coordinates)
              Every Compare mode works with every Display, so e.g. a
              static-shift correction can be judged as a 2-D before/after
              section, a 1-D overlay for one station, or a ρ_a-ratio / Δφ
              diff.  Rendering lives in ``controllers/correction_views.py``.
              When a view cannot be drawn the canvas is replaced by a card
              explaining why — never left as empty axes.

The correction stack is non-destructive: raw Sites are never modified.
``corrections_committed`` signal carries the final corrected Sites so
MainWindow can replace the global dataset.

``corrections_reverted`` notifies MainWindow that no corrections are active.
"""

from __future__ import annotations

import warnings

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
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
    "Distortion": "◈",
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
from pycsamt.app.desktop.controllers.correction_views import (
    ALL_STATIONS,
    COMPARE_MODES,
    COMPONENT_CHOICES,
    QUANTITY_CHOICES,
    PlotUnavailable,
    extract_responses,
    figure_blank_reason,
    render_curves,
    render_section,
)
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    icon_button,
    make_group,
)
from pycsamt.app.desktop.widgets.compact_button import compact_button


# Stacked-widget pages of the right-hand panel
_PAGE_COMPARE = 0
_PAGE_STRAT = 1

# Display choices per category family
DISPLAY_CURVES = "Curves (1-D)"
DISPLAY_SECTION = "Pseudosection (2-D)"
DISPLAY_ROSE = "Strike rose"
DISPLAY_MAP = "Position map"
DISPLAY_ELEV = "Elevation profile"
# Static shift is a lateral (station-to-station) effect: open it as a section.
_DEFAULT_DISPLAY = {"Static Shift": DISPLAY_SECTION}


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
            params_width=270,
            icon_name="sites-correction",
            parent=parent,
        )
        # MainWindow clamps panels to 75 % of the screen (960 px on a
        # 1280-px-wide laptop display); keep the minimum well below that.
        self.resize(1180, 760)
        self._ctrl = CorrectionController()
        self._preview_sites = None  # result of last Preview (not in stack)
        # Last Display chosen per category, restored when switching back
        self._display_memory: dict[str, str] = dict(_DEFAULT_DISPLAY)

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
        compact_button(btn_browse, 28, 26)
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
        # Two orthogonal choices drive every comparison view:
        #   Compare  — HOW before and after are contrasted
        #              (Before / After, Overlay, Diff)
        #   Display  — WHAT is drawn (1-D curves, 2-D pseudosection,
        #              strike rose, station map, ...), per category
        # plus quantity / component / station refinements for Z views.
        #
        # Laid out on two short rows with fixed-length combos: one long row
        # of combos sized to their longest entry (e.g. station names) forced
        # a ~1700 px minimum window width that did not fit on laptop screens.
        bar = QHBoxLayout()
        bar.setContentsMargins(8, 4, 8, 2)
        bar.setSpacing(6)

        self._view_controls = QWidget()
        vc_rows = QVBoxLayout(self._view_controls)
        vc_rows.setContentsMargins(0, 0, 0, 0)
        vc_rows.setSpacing(3)
        vc = QHBoxLayout()
        vc.setSpacing(6)
        vc_rows.addLayout(vc)

        def _labelled(text: str, widget: QWidget, tip: str) -> None:
            lbl = QLabel(text)
            lbl.setObjectName("FieldLabel")
            widget.setToolTip(tip)
            vc.addWidget(lbl)
            vc.addWidget(widget)

        self._combo_mode = QComboBox()
        self._combo_mode.addItems(list(COMPARE_MODES))
        self._combo_mode.currentIndexChanged.connect(self._on_mode_changed)
        _labelled(
            "Compare:",
            self._combo_mode,
            "Before / After — raw and corrected side by side, same scale\n"
            "Overlay — both states on one plot (before dashed/grey or as "
            "contours)\n"
            "Diff — what the correction changed: ρ_a ratio (×) and Δφ (°)",
        )

        self._combo_display = QComboBox()
        self._combo_display.currentIndexChanged.connect(
            self._on_display_changed
        )
        _labelled("Display:", self._combo_display, "What to draw")

        self._z_controls = QWidget()
        zc = QHBoxLayout(self._z_controls)
        zc.setContentsMargins(0, 0, 0, 0)
        zc.setSpacing(6)
        self._combo_quantity = QComboBox()
        self._combo_quantity.addItems(list(QUANTITY_CHOICES))
        self._combo_quantity.setToolTip(
            "ρ_a + φ is recommended: a correction that alters phase (e.g. "
            "rotation, distortion removal) is only visible in φ, while a "
            "pure static shift must leave φ unchanged."
        )
        self._combo_component = QComboBox()
        self._combo_component.addItems(list(COMPONENT_CHOICES))
        self._combo_component.setToolTip("Impedance component(s) to show")
        self._combo_station = QComboBox()
        self._combo_station.setMaxVisibleItems(20)
        self._combo_station.addItem(ALL_STATIONS)
        self._combo_station.setToolTip(
            "Curves: show only this station's before/after.\n"
            "Pseudosection: outline this station's column."
        )
        # One "Show:" label for the row: the entries ("ρ_a + φ", "XY + YX",
        # "All stations" / a station name) are self-describing, and three
        # separate labels cost ~150 px of window width.
        lbl = QLabel("Show:")
        lbl.setObjectName("FieldLabel")
        zc.addWidget(lbl)
        for w in (self._combo_quantity, self._combo_component,
                  self._combo_station):
            zc.addWidget(w)
            w.currentIndexChanged.connect(self._refresh_plots)
        zc.addStretch(1)
        vc.addStretch(1)
        vc_rows.addWidget(self._z_controls)
        for combo, chars in (
            (self._combo_mode, 12),
            (self._combo_display, 15),
            (self._combo_quantity, 7),
            (self._combo_component, 7),
            (self._combo_station, 9),
        ):
            _compact_combo(combo, chars)
        bar.addWidget(self._view_controls, 1)

        side = QVBoxLayout()
        side.setSpacing(3)
        self._btn_export = icon_button(
            "⬆  Export…", "export", "Export current figure"
        )
        self._btn_export.setFixedWidth(110)
        self._btn_export.clicked.connect(self._on_export)
        side.addWidget(self._btn_export, 0, Qt.AlignmentFlag.AlignRight)
        self._view_status = QLabel("")
        self._view_status.setObjectName("InfoLabel")
        self._view_status.setAlignment(Qt.AlignmentFlag.AlignRight)
        # Long status messages must not widen the window
        self._view_status.setSizePolicy(
            QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred
        )
        side.addWidget(self._view_status)
        side_w = QWidget()
        side_w.setLayout(side)
        side_w.setFixedWidth(120)
        bar.addWidget(side_w)

        bar_w = QWidget()
        bar_w.setLayout(bar)
        layout.addWidget(bar_w)

        # ── Stacked view: page 0 = comparison canvas, page 1 = Stratagem ──
        self._view_stack = QStackedWidget()
        layout.addWidget(self._view_stack)

        # ── Page 0: one canvas for every comparison view ──────────────
        # A single figure (not separate Before/After canvases) keeps both
        # states on a shared axis scale and exports as one image.
        cmp_page = QWidget()
        cmp_v = QVBoxLayout(cmp_page)
        cmp_v.setContentsMargins(0, 0, 0, 0)
        cmp_v.setSpacing(0)
        self._plot_view = CanvasResultView(
            cmp_page,
            toolbar=True,
            empty_title="No data loaded",
            empty_reason="Load survey data, then Preview or Apply a correction.",
        )
        self._canvas = self._plot_view.canvas
        self._canvas.set_refresh_callback(
            self._refresh_plots, tooltip="Redraw current view"
        )
        cmp_v.addWidget(self._plot_view)
        self._view_stack.addWidget(cmp_page)

        # ── Page 1: Stratagem Studio ───────────────────────────────────
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
        self._canvas_strat_qc_view = CanvasResultView(
            qc_page,
            toolbar=True,
            empty_title="No QC report yet",
            empty_reason="Load Stratagem EDI data to generate a QC report.",
        )
        self._canvas_strat_qc = self._canvas_strat_qc_view.canvas
        self._canvas_strat_qc.set_refresh_callback(
            self._refresh_strat_plots, tooltip="Redraw Stratagem plots"
        )
        qc_v.addWidget(self._canvas_strat_qc_view)
        self._strat_tabs.addTab(qc_page, "QC Report")

        # Tab 1 — Before / After impedance
        ba_page = QWidget()
        ba_v = QVBoxLayout(ba_page)
        ba_v.setContentsMargins(0, 0, 0, 0)
        ba_splitter = QSplitter(Qt.Orientation.Vertical)
        ba_splitter.setHandleWidth(4)

        def _strat_pane(title, empty_title, empty_reason):
            p = QWidget()
            p.setObjectName("CanvasPane")
            v = QVBoxLayout(p)
            v.setContentsMargins(0, 0, 0, 0)
            v.setSpacing(0)
            lbl = QLabel(f"  {title}")
            lbl.setObjectName("CanvasLabel")
            lbl.setFixedHeight(22)
            v.addWidget(lbl)
            view = CanvasResultView(
                p, toolbar=False,
                empty_title=empty_title, empty_reason=empty_reason,
            )
            v.addWidget(view)
            return p, view

        ba_before_pane, self._canvas_strat_before_view = _strat_pane(
            "Before (raw)",
            "No data loaded",
            "Load Stratagem EDI data to see the raw curves.",
        )
        ba_after_pane, self._canvas_strat_after_view = _strat_pane(
            "After (corrected)",
            "No corrections applied yet",
            "Apply a Stratagem correction to see the corrected curves.",
        )
        self._canvas_strat_before = self._canvas_strat_before_view.canvas
        self._canvas_strat_after = self._canvas_strat_after_view.canvas
        for view in (self._canvas_strat_before_view, self._canvas_strat_after_view):
            view.canvas.set_refresh_callback(
                self._refresh_strat_plots, tooltip="Redraw Stratagem plots"
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
        self._canvas_strat_ss_view = CanvasResultView(
            ss_page,
            toolbar=True,
            empty_title="No static-shift factors yet",
            empty_reason="Load Stratagem EDI data to compute static-shift factors.",
        )
        self._canvas_strat_ss = self._canvas_strat_ss_view.canvas
        self._canvas_strat_ss.set_refresh_callback(
            self._refresh_strat_plots, tooltip="Redraw Stratagem plots"
        )
        ss_v.addWidget(self._canvas_strat_ss_view)
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

        self._view_stack.setCurrentIndex(_PAGE_COMPARE)

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
        self._populate_station_combo()
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
        self._view_controls.setVisible(not is_strat)
        self._btn_export_strat.setVisible(is_strat)

        if is_strat:
            self._view_stack.setCurrentIndex(_PAGE_STRAT)
            self._refresh_strat_plots()
        else:
            self._view_stack.setCurrentIndex(_PAGE_COMPARE)
            self._populate_display_combo(cat)

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
                if self._view_stack.currentIndex() == _PAGE_STRAT:
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
            self._canvas_strat_qc_view.show_canvas()
        except Exception as exc:
            self._canvas_strat_qc_view.show_unavailable(
                "No QC report yet",
                str(exc) or "Load Stratagem EDI data to generate a QC report.",
            )
        # SS factors
        ax_ss = self._canvas_strat_ss.axes
        try:
            self._ctrl.plot_strat_ss_factors(ax_ss)
            self._canvas_strat_ss.draw()
            self._canvas_strat_ss_view.show_canvas()
        except Exception as exc:
            self._canvas_strat_ss_view.show_unavailable(
                "No static-shift factors yet",
                str(exc) or "Load Stratagem EDI data to compute static-shift factors.",
            )
        # Before / After impedance curves
        raw_sites = self._ctrl.raw_sites
        curr_sites = self._ctrl.current_sites
        if raw_sites is not None:
            try:
                ax_b = self._canvas_strat_before.axes
                ax_b.cla()
                self._ctrl.plot_rho_curves(raw_sites, ax_b, "Before (raw)")
                self._canvas_strat_before.draw()
                self._canvas_strat_before_view.show_canvas()
            except Exception as exc:
                self._canvas_strat_before_view.show_unavailable(
                    "Plot unavailable", str(exc)
                )
        else:
            self._canvas_strat_before_view.show_unavailable(
                "No data loaded", "Load Stratagem EDI data to see the raw curves."
            )
        if curr_sites is not None and curr_sites is not raw_sites:
            try:
                ax_a = self._canvas_strat_after.axes
                ax_a.cla()
                self._ctrl.plot_rho_curves(
                    curr_sites, ax_a, "After (corrected)"
                )
                self._canvas_strat_after.draw()
                self._canvas_strat_after_view.show_canvas()
            except Exception as exc:
                self._canvas_strat_after_view.show_unavailable(
                    "Plot unavailable", str(exc)
                )
        else:
            self._canvas_strat_after_view.show_unavailable(
                "No corrections applied yet",
                "Apply a Stratagem correction to see the corrected curves.",
            )
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
        self._refresh_plots()

    def _on_display_changed(self, _idx: int) -> None:
        cat = self._current_category()
        display = self._combo_display.currentText()
        if cat and display:
            self._display_memory[cat] = display
        # Quantity / component / station only refine impedance views
        self._z_controls.setVisible(display in (DISPLAY_CURVES, DISPLAY_SECTION))
        self._refresh_plots()

    def _populate_display_combo(self, cat: str) -> None:
        """Offer the Display choices that make sense for *cat*."""
        if cat in COORD_CATEGORIES:
            choices = [DISPLAY_MAP, DISPLAY_ELEV]
        else:
            choices = [DISPLAY_CURVES, DISPLAY_SECTION]
            if cat in ROTATION_CATEGORIES:
                choices.append(DISPLAY_ROSE)
        wanted = self._display_memory.get(cat, choices[0])
        self._combo_display.blockSignals(True)
        self._combo_display.clear()
        self._combo_display.addItems(choices)
        self._combo_display.setCurrentText(
            wanted if wanted in choices else choices[0]
        )
        self._combo_display.blockSignals(False)
        _fit_popup(self._combo_display)
        self._z_controls.setVisible(
            self._combo_display.currentText()
            in (DISPLAY_CURVES, DISPLAY_SECTION)
        )

    def _populate_station_combo(self) -> None:
        """Refill the Station picker from the raw dataset, keeping the
        current choice when that station still exists."""
        current = self._combo_station.currentText()
        names = []
        if self._ctrl.raw_sites is not None:
            try:
                names = list(extract_responses(self._ctrl.raw_sites))
            except Exception:
                names = []
        self._combo_station.blockSignals(True)
        self._combo_station.clear()
        self._combo_station.addItem(ALL_STATIONS)
        self._combo_station.addItems(names)
        self._combo_station.setCurrentText(
            current if current in names else ALL_STATIONS
        )
        self._combo_station.blockSignals(False)
        _fit_popup(self._combo_station)

    def _current_category(self) -> str:
        row = self._combo_category.currentIndex()
        return CATEGORIES[row] if 0 <= row < len(CATEGORIES) else ""

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

        ExportDialog(figure=self._canvas.figure, parent=self).exec()

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
        """Redraw the comparison canvas for the current Compare × Display.

        Never leaves empty axes on screen: anything that prevents a real
        plot (no data, nothing to compare yet, missing component, a
        plotting error) swaps the canvas for a card stating the reason.
        """
        if self._current_category() in STRATAGEM_CATEGORIES:
            return
        view = self._plot_view
        fig = self._canvas.figure
        if not self._ctrl.has_data:
            view.show_unavailable(
                "No data loaded",
                "There is no survey in this panel yet.",
                "Load EDI / EMTF-XML data in the main window; it is sent "
                "here automatically.",
            )
            return

        with warnings.catch_warnings():
            # Clearing shared log axes briefly resets limits to (0, 1)
            warnings.simplefilter("ignore", UserWarning)
            fig.clear()
        # Constrained layout keeps shared colour bars and side-by-side
        # panels aligned at any canvas size.
        fig.set_layout_engine("constrained")
        try:
            self._render_comparison(fig)
            reason = figure_blank_reason(fig)
            if reason:
                raise PlotUnavailable(
                    "Nothing to display for this view",
                    reason,
                    "Try another Display or Compare mode, or check the "
                    "correction parameters.",
                )
        except PlotUnavailable as exc:
            fig.clear()
            view.show_unavailable(exc.title, exc.reason, exc.guidance)
            return
        except Exception as exc:  # a plotting bug must not blank the UI
            fig.clear()
            view.show_unavailable(
                "This view could not be drawn",
                f"{type(exc).__name__}: {exc}",
                "Try another Display or Compare mode. If it persists, "
                "please report it with the data that triggers it.",
            )
            return
        self._canvas.draw()
        view.show_canvas()

    def _after_state(self):
        """Return (after_data, after_title) for the current comparison."""
        n = self._ctrl.n_steps
        if self._preview_sites is not None:
            fn, label = self._current_fn_label()
            return self._preview_sites, f"Preview: {label}" if label else "Preview"
        if n == 0:
            title = "After (no corrections yet)"
        else:
            title = f"After ({n} step{'s' if n != 1 else ''})"
        if self._is_coord_category():
            return self._ctrl.current_coords_df(), title
        return self._ctrl.current_sites, title

    def _render_comparison(self, fig) -> None:
        mode = self._combo_mode.currentText()
        display = self._combo_display.currentText()
        after_data, after_title = self._after_state()

        if self._is_coord_category():
            self._render_coords(fig, mode, display, after_data, after_title)
            return

        before = self._ctrl.raw_sites
        if display == DISPLAY_ROSE:
            if mode == "Diff":
                raise PlotUnavailable(
                    "No Diff view for a strike rose",
                    "A rose diagram compares two direction distributions; "
                    "subtracting them bin by bin has no physical meaning.",
                    "Use Before / After or Overlay for the rose, or switch "
                    "Display to Curves / Pseudosection to see how the "
                    "rotation changed ρ_a and φ.",
                )
            if mode == "Overlay":
                self._ctrl.plot_rotation_rose_overlay(before, after_data, fig)
            else:
                self._ctrl.plot_rotation_rose(before, after_data, fig)
            return

        theme = _theme(self._ctrl.dark)
        station = self._combo_station.currentText()
        quantities = QUANTITY_CHOICES[self._combo_quantity.currentText()]
        components = COMPONENT_CHOICES[self._combo_component.currentText()]
        if display == DISPLAY_SECTION:
            render_section(
                fig,
                before,
                after_data,
                mode=mode,
                station=station,
                quantities=quantities,
                components=components,
                affected_stations=(
                    self._get_affected_stations()
                    if self._is_ss_category()
                    else None
                ),
                theme=theme,
                after_title=after_title,
            )
        else:
            render_curves(
                fig,
                before,
                after_data,
                mode=mode,
                station=station,
                quantities=quantities,
                components=components,
                theme=theme,
                after_title=after_title,
            )

    def _render_coords(self, fig, mode, display, after, after_title) -> None:
        before = self._ctrl.raw_coords_df
        elev = display == DISPLAY_ELEV
        if mode == "Diff":
            if _coords_equal(before, after):
                raise PlotUnavailable(
                    "Nothing to compare yet",
                    "No coordinate correction has been previewed or "
                    "applied, so every station is still at its raw "
                    "position.",
                    "Choose a correction and click Preview or Apply.",
                )
            self._ctrl.plot_displacement_diff(before, after, fig.add_subplot(111))
        elif mode == "Overlay":
            ax = fig.add_subplot(111)
            if elev:
                self._ctrl.plot_station_elevation_overlay(before, after, ax)
            else:
                self._ctrl.plot_station_map_overlay(before, after, ax)
        else:
            # Not sharex/sharey: the station map uses an equal aspect with
            # adjustable="datalim", which matplotlib forbids on axes shared
            # in both directions. Identical limits are imposed afterwards.
            ax_b, ax_a = fig.subplots(1, 2)
            plot = (
                self._ctrl.plot_station_elevation
                if elev
                else self._ctrl.plot_station_map
            )
            what = "elevation" if elev else "positions"
            plot(before, ax_b, f"Before (raw) — {what}")
            plot(after, ax_a, f"{after_title} — {what}")
            for get, set_ in (("get_xlim", "set_xlim"), ("get_ylim", "set_ylim")):
                lims = [getattr(ax, get)() for ax in (ax_b, ax_a)]
                lo = min(min(l) for l in lims)
                hi = max(max(l) for l in lims)
                for ax in (ax_b, ax_a):
                    getattr(ax, set_)(lo, hi)


def _compact_combo(combo: QComboBox, chars: int) -> None:
    """Size *combo* for ~*chars* characters instead of its longest item.

    The drop-down list still shows full item text; only the closed box is
    capped, so long station names cannot inflate the window's minimum width.
    """
    combo.setSizeAdjustPolicy(
        QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
    )
    combo.setMinimumContentsLength(chars)
    _fit_popup(combo)


def _fit_popup(combo: QComboBox) -> None:
    """Let the drop-down list be as wide as its longest item."""
    combo.view().setMinimumWidth(combo.view().sizeHintForColumn(0) + 24)


def _theme(dark: bool) -> dict:
    from pycsamt.app.desktop.controllers.correction_controller import (
        _DARK,
        _LIGHT,
    )

    return _DARK if dark else _LIGHT


def _coords_equal(a, b) -> bool:
    """True when two coordinate tables hold the same station positions."""
    try:
        if a is None or b is None:
            return a is b
        cols = [c for c in ("station", "lat", "lon", "elev") if c in a.columns]
        return a[cols].reset_index(drop=True).equals(
            b[cols].reset_index(drop=True)
        )
    except Exception:
        return False


# ── Needed for type hint in _make_widget ──────────────────────────────────────
import numpy as np
