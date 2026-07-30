# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AdvancedToolsWindow — independent floating window for advanced emtools analyses,
topography configuration, and format conversion.

Left params panel
─────────────────
  Category     Sidebar list (Strike · Phase Tensor · Induction · Impedance ·
                              Depth Imaging · Survey Tools · Topography · Conversion)
  [Stacked pages based on category selection]

Right content
─────────────
  Stacked pages:
    Page 0 — single large MplCanvas (for emtools plots)
    Page 1 — topography preview canvas + view selector
    Page 2 — conversion results tabs (table + curves + map)
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.advanced_controller import (
    ADVANCED_GROUP_ICONS,
    ADVANCED_GROUPS,
    CONV_INDEX,
    TOPO_INDEX,
    AdvancedController,
    ConversionController,
    ConversionWorker,
    DimModelWorker,
    TopoPreviewController,
    describe_advanced_plot,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    _icon,
    icon_button,
    make_group,
)


class AdvancedToolsWindow(PanelWindow):
    """
    Floating Advanced Tools window.

    Pages
    -----
    0-5  emtools plots (Strike / Phase Tensor / Induction / Impedance /
                        Depth / Survey)
    6    Topography configuration + preview
    7    Format conversion (AVG / J / Spectra -> EDI)
    """

    # Emitted when the user commits a conversion result as the main dataset
    conversion_committed = Signal(object)  # EDICollection

    def __init__(self, parent: QWidget | None = None) -> None:
        # Controllers created before super().__init__ so they are available
        # inside _build_params / _build_content which are called from __init__.
        self._ctrl = AdvancedController()
        self._topo_ctrl = TopoPreviewController()
        self._conv_ctrl = ConversionController()

        # Topo color state
        self._topo_fill_color = "#a89070"
        self._topo_line_color = "#6b4e2a"

        # Conversion worker handle + result rows
        self._conv_worker: ConversionWorker | None = None
        self._conv_running = False
        self._conv_result_stats: list = []

        # Dictionary-model training worker
        self._dim_worker: DimModelWorker | None = None

        super().__init__(
            title="Advanced Tools",
            session_key="advanced_tools",
            params_width=310,
            icon_name="advanced-tools",
            parent=parent,
        )
        self.resize(1240, 820)
        self._auto_rendered = False
        self._populate_category_combo()
        self._on_category_changed(0)

    # ── Params panel (left) ───────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:
        # ── Category selector ─────────────────────────────────────────
        grp_cat, lay_cat = make_group("Category")
        self._combo_category = QComboBox()
        self._combo_category.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._combo_category.currentIndexChanged.connect(
            self._on_category_changed
        )
        lay_cat.addWidget(self._combo_category)
        layout.addWidget(grp_cat)

        # ── Stacked pages for per-category params ─────────────────────
        self._params_stack = QStackedWidget()
        self._params_stack.addWidget(self._build_plot_params_page())  # 0
        self._params_stack.addWidget(self._build_topo_params_page())  # 1
        self._params_stack.addWidget(self._build_conv_params_page())  # 2
        layout.addWidget(self._params_stack)

    # ── Plot params page (page 0) ─────────────────────────────────────

    def _build_plot_params_page(self) -> QWidget:
        page = QWidget()
        vlay = QVBoxLayout(page)
        vlay.setContentsMargins(0, 0, 0, 0)
        vlay.setSpacing(5)

        # Plot combo
        grp_plot, lay_plot = make_group("Plot")
        self._combo_plot = QComboBox()
        self._combo_plot.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._combo_plot.currentIndexChanged.connect(self._on_plot_changed)
        lay_plot.addWidget(self._combo_plot)
        vlay.addWidget(grp_plot)

        # Description
        self._desc_lbl = QLabel("")
        self._desc_lbl.setWordWrap(True)
        self._desc_lbl.setObjectName("InfoLabel")
        self._desc_lbl.setAlignment(Qt.AlignmentFlag.AlignTop)
        vlay.addWidget(self._desc_lbl)

        # ── Dictionary Model group (shown only for ATOM psection) ─────
        grp_model, lay_model = make_group("Dictionary Model")

        row_natoms = QWidget()
        h_na = QHBoxLayout(row_natoms)
        h_na.setContentsMargins(0, 0, 0, 0)
        h_na.addWidget(QLabel("n_atoms"))
        self._spin_n_atoms = QDoubleSpinBox()
        self._spin_n_atoms.setDecimals(0)
        self._spin_n_atoms.setRange(2, 20)
        self._spin_n_atoms.setSingleStep(1)
        self._spin_n_atoms.setValue(6)
        h_na.addWidget(self._spin_n_atoms)
        lay_model.addWidget(row_natoms)

        row_niter = QWidget()
        h_ni = QHBoxLayout(row_niter)
        h_ni.setContentsMargins(0, 0, 0, 0)
        h_ni.addWidget(QLabel("n_iter"))
        self._spin_n_iter = QDoubleSpinBox()
        self._spin_n_iter.setDecimals(0)
        self._spin_n_iter.setRange(10, 200)
        self._spin_n_iter.setSingleStep(10)
        self._spin_n_iter.setValue(40)
        h_ni.addWidget(self._spin_n_iter)
        lay_model.addWidget(row_niter)

        self._btn_train_model = QPushButton("⚙  Train Model from Survey")
        self._btn_train_model.setToolTip(
            "Learn a sparse-coding dictionary from the phase-tensor features\n"
            "of the loaded sites.  Required before plotting ATOM psection."
        )
        self._btn_train_model.clicked.connect(self._on_train_model)
        lay_model.addWidget(self._btn_train_model)

        self._model_status_lbl = QLabel("Not trained")
        self._model_status_lbl.setObjectName("InfoLabel")
        self._model_status_lbl.setWordWrap(True)
        lay_model.addWidget(self._model_status_lbl)

        self._grp_model = grp_model
        self._grp_model.setVisible(False)
        vlay.addWidget(grp_model)

        # Actions
        grp_act, lay_act = make_group("Actions")
        self._btn_run = icon_button(
            "↻  Run / Refresh", "advanced-tools", "Render selected plot"
        )
        self._btn_export = icon_button(
            "⬆  Export…", "export", "Save figure to file"
        )
        self._btn_run.clicked.connect(self._on_run)
        self._btn_export.clicked.connect(self._on_export)
        lay_act.addWidget(self._btn_run)
        lay_act.addWidget(self._btn_export)
        vlay.addWidget(grp_act)

        # Status
        self._status_lbl = QLabel("")
        self._status_lbl.setObjectName("InfoLabel")
        self._status_lbl.setWordWrap(True)
        vlay.addWidget(self._status_lbl)

        vlay.addStretch(1)
        return page

    # ── Topo params page (page 1) ─────────────────────────────────────

    def _build_topo_params_page(self) -> QWidget:
        page = QWidget()
        vlay = QVBoxLayout(page)
        vlay.setContentsMargins(0, 0, 0, 0)
        vlay.setSpacing(4)

        # ── Configuration group ───────────────────────────────────────
        grp_cfg, lay_cfg = make_group("Configuration")

        self._chk_topo_enabled = QCheckBox("Enable topography")
        lay_cfg.addWidget(self._chk_topo_enabled)

        row_src = QWidget()
        h_src = QHBoxLayout(row_src)
        h_src.setContentsMargins(0, 0, 0, 0)
        h_src.addWidget(QLabel("Source"))
        self._combo_topo_source = QComboBox()
        self._combo_topo_source.addItems(["sites", "file", "array"])
        self._combo_topo_source.currentTextChanged.connect(
            self._on_topo_source_changed
        )
        h_src.addWidget(self._combo_topo_source)
        lay_cfg.addWidget(row_src)

        # File row (shown only when source == "file")
        self._topo_file_row = QWidget()
        h_file = QHBoxLayout(self._topo_file_row)
        h_file.setContentsMargins(0, 0, 0, 0)
        self._edit_topo_file = QLineEdit()
        self._edit_topo_file.setPlaceholderText("Elevation file path…")
        btn_browse_topo = QPushButton("…")
        btn_browse_topo.setFixedWidth(28)
        btn_browse_topo.clicked.connect(self._browse_topo_file)
        h_file.addWidget(self._edit_topo_file)
        h_file.addWidget(btn_browse_topo)
        self._topo_file_row.setVisible(False)
        lay_cfg.addWidget(self._topo_file_row)

        row_interp = QWidget()
        h_interp = QHBoxLayout(row_interp)
        h_interp.setContentsMargins(0, 0, 0, 0)
        h_interp.addWidget(QLabel("Interpolation"))
        self._combo_topo_interp = QComboBox()
        self._combo_topo_interp.addItems(["linear", "cubic", "nearest"])
        h_interp.addWidget(self._combo_topo_interp)
        lay_cfg.addWidget(row_interp)

        row_exag = QWidget()
        h_exag = QHBoxLayout(row_exag)
        h_exag.setContentsMargins(0, 0, 0, 0)
        h_exag.addWidget(QLabel("Exaggeration"))
        self._spin_topo_exag = QDoubleSpinBox()
        self._spin_topo_exag.setRange(0.1, 20.0)
        self._spin_topo_exag.setSingleStep(0.1)
        self._spin_topo_exag.setValue(1.0)
        h_exag.addWidget(self._spin_topo_exag)
        lay_cfg.addWidget(row_exag)

        vlay.addWidget(grp_cfg)

        # ── Style group ───────────────────────────────────────────────
        grp_style, lay_style = make_group("Style")

        row_fc = QWidget()
        h_fc = QHBoxLayout(row_fc)
        h_fc.setContentsMargins(0, 0, 0, 0)
        h_fc.addWidget(QLabel("Fill color"))
        self._btn_topo_fill_color = QPushButton()
        self._btn_topo_fill_color.setFixedWidth(36)
        self._apply_topo_btn_color(
            self._btn_topo_fill_color, self._topo_fill_color
        )
        self._btn_topo_fill_color.clicked.connect(
            lambda: self._pick_color("fill")
        )
        h_fc.addWidget(self._btn_topo_fill_color)
        lay_style.addWidget(row_fc)

        row_fa = QWidget()
        h_fa = QHBoxLayout(row_fa)
        h_fa.setContentsMargins(0, 0, 0, 0)
        h_fa.addWidget(QLabel("Fill α"))
        self._spin_topo_fill_alpha = QDoubleSpinBox()
        self._spin_topo_fill_alpha.setRange(0.0, 1.0)
        self._spin_topo_fill_alpha.setSingleStep(0.05)
        self._spin_topo_fill_alpha.setValue(0.4)
        h_fa.addWidget(self._spin_topo_fill_alpha)
        lay_style.addWidget(row_fa)

        row_lc = QWidget()
        h_lc = QHBoxLayout(row_lc)
        h_lc.setContentsMargins(0, 0, 0, 0)
        h_lc.addWidget(QLabel("Line color"))
        self._btn_topo_line_color = QPushButton()
        self._btn_topo_line_color.setFixedWidth(36)
        self._apply_topo_btn_color(
            self._btn_topo_line_color, self._topo_line_color
        )
        self._btn_topo_line_color.clicked.connect(
            lambda: self._pick_color("line")
        )
        h_lc.addWidget(self._btn_topo_line_color)
        lay_style.addWidget(row_lc)

        row_lw = QWidget()
        h_lw = QHBoxLayout(row_lw)
        h_lw.setContentsMargins(0, 0, 0, 0)
        h_lw.addWidget(QLabel("Line width"))
        self._spin_topo_line_w = QDoubleSpinBox()
        self._spin_topo_line_w.setRange(0.1, 5.0)
        self._spin_topo_line_w.setSingleStep(0.1)
        self._spin_topo_line_w.setValue(1.2)
        h_lw.addWidget(self._spin_topo_line_w)
        lay_style.addWidget(row_lw)

        self._chk_surface_line = QCheckBox("Surface line")
        self._chk_surface_line.setChecked(True)
        lay_style.addWidget(self._chk_surface_line)

        self._chk_clip_below = QCheckBox("Clip below surface")
        self._chk_clip_below.setChecked(True)
        lay_style.addWidget(self._chk_clip_below)

        self._chk_pins_at_surface = QCheckBox("Station pins at surface")
        self._chk_pins_at_surface.setChecked(True)
        lay_style.addWidget(self._chk_pins_at_surface)

        vlay.addWidget(grp_style)

        # ── Pseudosection strip group ─────────────────────────────────
        grp_strip, lay_strip = make_group("Pseudosection Strip")

        self._chk_topo_strip = QCheckBox("Show topo strip")
        self._chk_topo_strip.setChecked(True)
        lay_strip.addWidget(self._chk_topo_strip)

        row_sh = QWidget()
        h_sh = QHBoxLayout(row_sh)
        h_sh.setContentsMargins(0, 0, 0, 0)
        h_sh.addWidget(QLabel("Strip height"))
        self._spin_strip_h = QDoubleSpinBox()
        self._spin_strip_h.setRange(0.05, 0.5)
        self._spin_strip_h.setSingleStep(0.02)
        self._spin_strip_h.setValue(0.18)
        h_sh.addWidget(self._spin_strip_h)
        lay_strip.addWidget(row_sh)

        vlay.addWidget(grp_strip)

        # ── Actions group ─────────────────────────────────────────────
        grp_ta, lay_ta = make_group("Actions")
        self._btn_topo_apply = icon_button(
            "↻  Apply to Global Config", "", "Apply settings to PYCSAMT_TOPO"
        )
        self._btn_topo_reset = icon_button(
            "↩  Reset to Defaults", "", "Reset PYCSAMT_TOPO to defaults"
        )
        self._btn_topo_preview = icon_button(
            "▶  Preview / Refresh", "", "Refresh the preview canvas"
        )
        self._btn_topo_apply.clicked.connect(self._on_topo_apply)
        self._btn_topo_reset.clicked.connect(self._on_topo_reset)
        self._btn_topo_preview.clicked.connect(self._refresh_topo_preview)
        lay_ta.addWidget(self._btn_topo_apply)
        lay_ta.addWidget(self._btn_topo_reset)
        lay_ta.addWidget(self._btn_topo_preview)
        vlay.addWidget(grp_ta)

        self._topo_status = QLabel("")
        self._topo_status.setObjectName("InfoLabel")
        self._topo_status.setWordWrap(True)
        vlay.addWidget(self._topo_status)

        vlay.addStretch(1)
        return page

    # ── Conversion params page (page 2) ──────────────────────────────

    def _build_conv_params_page(self) -> QWidget:
        page = QWidget()
        vlay = QVBoxLayout(page)
        vlay.setContentsMargins(0, 0, 0, 0)
        vlay.setSpacing(4)

        # ── Input group ───────────────────────────────────────────────
        grp_inp, lay_inp = make_group("Input")

        row_type = QWidget()
        h_type = QHBoxLayout(row_type)
        h_type.setContentsMargins(0, 0, 0, 0)
        h_type.addWidget(QLabel("Type"))
        self._combo_conv_type = QComboBox()
        self._combo_conv_type.addItems(
            ["AVG → EDI", "J → EDI", "Spectra → EDI"]
        )
        self._combo_conv_type.currentIndexChanged.connect(
            self._on_conv_type_changed
        )
        h_type.addWidget(self._combo_conv_type)
        lay_inp.addWidget(row_type)

        self._conv_file_row = QWidget()
        h_cfile = QHBoxLayout(self._conv_file_row)
        h_cfile.setContentsMargins(0, 0, 0, 0)
        self._edit_conv_path = QLineEdit()
        self._edit_conv_path.setPlaceholderText("File or directory path…")
        self._edit_conv_path.textChanged.connect(self._update_conv_run_state)
        btn_browse_conv = QPushButton("📂")
        btn_browse_conv.setFixedWidth(28)
        btn_browse_conv.clicked.connect(self._browse_conv_path)
        h_cfile.addWidget(self._edit_conv_path)
        h_cfile.addWidget(btn_browse_conv)
        lay_inp.addWidget(self._conv_file_row)

        self._conv_dir_hint = QLabel("Select file or directory")
        self._conv_dir_hint.setObjectName("InfoLabel")
        lay_inp.addWidget(self._conv_dir_hint)

        vlay.addWidget(grp_inp)

        # ── Options group (stacked per type) ──────────────────────────
        grp_opts, lay_opts = make_group("Options")
        self._conv_opts_stack = QStackedWidget()

        # Page 0 — AVG opts
        avg_page = QWidget()
        avg_form = QFormLayout(avg_page)
        avg_form.setSpacing(5)
        self._avg_freq_order = QComboBox()
        self._avg_freq_order.addItems(["ascending", "descending"])
        avg_form.addRow("Freq order:", self._avg_freq_order)
        self._avg_freq_tol = QDoubleSpinBox()
        self._avg_freq_tol.setDecimals(9)
        self._avg_freq_tol.setRange(1e-9, 1e-2)
        self._avg_freq_tol.setSingleStep(1e-6)
        self._avg_freq_tol.setValue(1e-5)
        avg_form.addRow("Freq tol:", self._avg_freq_tol)
        self._avg_compute_z = QCheckBox()
        avg_form.addRow("Compute Z from ρ/φ:", self._avg_compute_z)
        self._avg_compute_rho = QCheckBox()
        avg_form.addRow("Compute ρ/φ from Z:", self._avg_compute_rho)
        self._conv_opts_stack.addWidget(avg_page)

        # Page 1 — J opts
        j_page = QWidget()
        j_form = QFormLayout(j_page)
        j_form.setSpacing(5)
        self._j_freq_order = QComboBox()
        self._j_freq_order.addItems(["ascending", "descending"])
        j_form.addRow("Freq order:", self._j_freq_order)
        self._j_station_name = QLineEdit()
        self._j_station_name.setPlaceholderText("Optional single-station name")
        j_form.addRow("Station name:", self._j_station_name)
        self._conv_opts_stack.addWidget(j_page)

        # Page 2 — Spectra opts
        sp_page = QWidget()
        sp_form = QFormLayout(sp_page)
        sp_form.setSpacing(5)
        self._sp_e_labels = QLineEdit("EX,EY")
        sp_form.addRow("E labels:", self._sp_e_labels)
        self._sp_h_labels = QLineEdit("HX,HY")
        sp_form.addRow("H labels:", self._sp_h_labels)
        self._sp_estimate_errors = QCheckBox()
        sp_form.addRow("Estimate errors:", self._sp_estimate_errors)
        self._sp_remote_ref = QCheckBox()
        sp_form.addRow("Use remote ref:", self._sp_remote_ref)
        self._sp_station_suffix = QLineEdit()
        sp_form.addRow("Station suffix:", self._sp_station_suffix)
        self._sp_skip_errors = QCheckBox()
        self._sp_skip_errors.setChecked(True)
        sp_form.addRow("Skip errors:", self._sp_skip_errors)
        self._conv_opts_stack.addWidget(sp_page)

        lay_opts.addWidget(self._conv_opts_stack)
        vlay.addWidget(grp_opts)

        # ── AVG topography / station profile group ───────────────────
        self._avg_topo_group, lay_avg_topo = make_group(
            "Add topography / station profile"
        )
        topo_hint = QLabel(
            "Optional: attach a Zonge .stn station profile so converted EDI "
            "files receive station elevation, longitude, and latitude."
        )
        topo_hint.setObjectName("InfoLabel")
        topo_hint.setWordWrap(True)
        lay_avg_topo.addWidget(topo_hint)

        row_stn = QWidget()
        h_stn = QHBoxLayout(row_stn)
        h_stn.setContentsMargins(0, 0, 0, 0)
        h_stn.addWidget(QLabel("File"))
        self._avg_stn_path = QLineEdit()
        self._avg_stn_path.setPlaceholderText(
            "K1.stn or station profile file..."
        )
        btn_stn = QPushButton("Browse...")
        btn_stn.setFixedWidth(76)
        btn_stn.clicked.connect(self._browse_avg_stn_path)
        h_stn.addWidget(self._avg_stn_path)
        h_stn.addWidget(btn_stn)
        lay_avg_topo.addWidget(row_stn)

        avg_topo_form = QFormLayout()
        avg_topo_form.setSpacing(5)
        self._avg_convert_coords = QCheckBox()
        self._avg_convert_coords.setChecked(True)
        avg_topo_form.addRow(
            "Convert UTM → lat/lon:", self._avg_convert_coords
        )
        self._avg_epsg = QLineEdit()
        self._avg_epsg.setPlaceholderText("EPSG code, e.g. 32650")
        avg_topo_form.addRow("EPSG:", self._avg_epsg)
        self._avg_utm_zone = QLineEdit()
        self._avg_utm_zone.setPlaceholderText("UTM zone, e.g. 50N")
        avg_topo_form.addRow("UTM zone:", self._avg_utm_zone)
        lay_avg_topo.addLayout(avg_topo_form)
        coord_hint = QLabel(
            "If conversion is enabled, provide either EPSG or UTM zone. "
            "Disable it only when the profile already contains latitude and "
            "longitude columns."
        )
        coord_hint.setObjectName("InfoLabel")
        coord_hint.setWordWrap(True)
        lay_avg_topo.addWidget(coord_hint)
        vlay.addWidget(self._avg_topo_group)

        # ── Output group (optional) ───────────────────────────────────
        grp_out, lay_out = make_group("Output")
        self._chk_write_edis = QCheckBox("Write EDI files to disk")
        lay_out.addWidget(self._chk_write_edis)
        row_outdir = QWidget()
        h_outdir = QHBoxLayout(row_outdir)
        h_outdir.setContentsMargins(0, 0, 0, 0)
        h_outdir.addWidget(QLabel("Output dir"))
        self._edit_out_dir = QLineEdit()
        self._edit_out_dir.setPlaceholderText("Output directory…")
        btn_browse_out = QPushButton("…")
        btn_browse_out.setFixedWidth(28)
        btn_browse_out.clicked.connect(self._browse_out_dir)
        h_outdir.addWidget(self._edit_out_dir)
        h_outdir.addWidget(btn_browse_out)
        lay_out.addWidget(row_outdir)
        vlay.addWidget(grp_out)

        # ── Actions group ─────────────────────────────────────────────
        grp_ca, lay_ca = make_group("Actions")

        self._btn_conv_run = icon_button(
            "▶▶  Run Transform", "", "Run the conversion"
        )
        self._btn_conv_run.setObjectName("CommitButton")
        self._btn_conv_run.clicked.connect(self._on_conv_run)

        self._btn_conv_commit = icon_button(
            "⬆  Commit as Main Dataset",
            "",
            "Use converted data as main dataset",
        )
        self._btn_conv_commit.setObjectName("CommitButton")
        self._btn_conv_commit.setEnabled(False)
        self._btn_conv_commit.clicked.connect(self._on_conv_commit)

        self._btn_conv_export = icon_button(
            "💾  Export EDIs…", "", "Write EDI files to disk"
        )
        self._btn_conv_export.setEnabled(False)
        self._btn_conv_export.clicked.connect(self._on_conv_export)

        self._btn_conv_clear = icon_button(
            "✗  Clear", "", "Clear conversion result"
        )
        self._btn_conv_clear.clicked.connect(self._on_conv_clear)

        lay_ca.addWidget(self._btn_conv_run)
        lay_ca.addWidget(self._btn_conv_commit)
        lay_ca.addWidget(self._btn_conv_export)
        lay_ca.addWidget(self._btn_conv_clear)
        vlay.addWidget(grp_ca)
        self._update_conv_run_state()

        self._conv_progress = QProgressBar()
        self._conv_progress.setRange(0, 0)  # indeterminate
        self._conv_progress.setVisible(False)
        vlay.addWidget(self._conv_progress)

        self._conv_status = QLabel("")
        self._conv_status.setObjectName("InfoLabel")
        self._conv_status.setWordWrap(True)
        vlay.addWidget(self._conv_status)

        vlay.addStretch(1)
        return page

    # ── Content panel (right) ─────────────────────────────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        self._content_stack = QStackedWidget()

        # ── Page 0: plot canvas ───────────────────────────────────────
        page0 = QWidget()
        v0 = QVBoxLayout(page0)
        v0.setContentsMargins(0, 0, 0, 0)
        self._canvas = MplCanvas(page0, toolbar=True)
        v0.addWidget(self._canvas)
        self._content_stack.addWidget(page0)

        # ── Page 1: topo content ──────────────────────────────────────
        page1 = QWidget()
        v1 = QVBoxLayout(page1)
        v1.setContentsMargins(4, 4, 4, 4)

        bar1 = QHBoxLayout()
        bar1.addWidget(QLabel("Preview:"))
        self._combo_topo_view = QComboBox()
        self._combo_topo_view.addItems(
            [
                "Elevation Profile",
                "Terrain Fill Preview",
                "Elevation Histogram",
            ]
        )
        self._combo_topo_view.currentIndexChanged.connect(
            self._on_topo_view_changed
        )
        bar1.addWidget(self._combo_topo_view)
        bar1.addStretch()
        self._topo_stats_lbl = QLabel("")
        self._topo_stats_lbl.setObjectName("InfoLabel")
        bar1.addWidget(self._topo_stats_lbl)
        bar_w1 = QWidget()
        bar_w1.setLayout(bar1)
        v1.addWidget(bar_w1)

        self._canvas_topo = MplCanvas(page1, toolbar=True)
        v1.addWidget(self._canvas_topo)
        self._content_stack.addWidget(page1)

        # ── Page 2: conversion content ────────────────────────────────
        page2 = QWidget()
        v2 = QVBoxLayout(page2)
        v2.setContentsMargins(0, 0, 0, 0)

        self._conv_tabs = QTabWidget()
        self._conv_tabs.setDocumentMode(True)

        # Tab 0 — Results table
        results_page = QWidget()
        rp_v = QVBoxLayout(results_page)
        self._conv_table = QTableWidget()
        self._conv_table.setAlternatingRowColors(True)
        self._conv_table.horizontalHeader().setStretchLastSection(True)
        self._conv_table.setEditTriggers(
            QTableWidget.EditTrigger.NoEditTriggers
        )
        rp_v.addWidget(self._conv_table)
        self._conv_tabs.addTab(results_page, "Results")

        # Tab 1 — Impedance Curves
        curves_page = QWidget()
        cp_v = QVBoxLayout(curves_page)
        self._canvas_conv_curves = MplCanvas(curves_page, toolbar=False)
        cp_v.addWidget(self._canvas_conv_curves)
        self._conv_tabs.addTab(curves_page, "Impedance Curves")

        # Tab 2 — Station Map
        map_page = QWidget()
        mp_v = QVBoxLayout(map_page)
        self._canvas_conv_map = MplCanvas(map_page, toolbar=False)
        mp_v.addWidget(self._canvas_conv_map)
        self._conv_tabs.addTab(map_page, "Station Map")

        v2.addWidget(self._conv_tabs)
        self._content_stack.addWidget(page2)

        layout.addWidget(self._content_stack)

    # ── Public API ────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        super().set_sites(sites)
        self._ctrl.set_sites(sites)
        self._topo_ctrl.set_sites(sites)
        self._auto_rendered = False
        # Refresh topo preview if that page is active
        if self._combo_category.currentIndex() == TOPO_INDEX:
            self._refresh_topo_preview()
        else:
            self._auto_render_if_ready()

    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._ctrl.dark = dark
        self._topo_ctrl.dark = dark
        self._conv_ctrl.dark = dark

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._auto_render_if_ready()

    # ── Category slot ─────────────────────────────────────────────────

    def _on_category_changed(self, row: int) -> None:
        if row < 0 or row >= len(ADVANCED_GROUPS):
            return

        if row == TOPO_INDEX:
            self._params_stack.setCurrentIndex(1)
            self._content_stack.setCurrentIndex(1)
            self._refresh_topo_preview()
            return

        if row == CONV_INDEX:
            self._params_stack.setCurrentIndex(2)
            self._content_stack.setCurrentIndex(2)
            return

        # Regular emtools plot category
        self._params_stack.setCurrentIndex(0)
        self._content_stack.setCurrentIndex(0)

        _label, plots = ADVANCED_GROUPS[row]
        self._combo_plot.blockSignals(True)
        self._combo_plot.clear()
        for label, _fn, _has_ax in plots:
            self._combo_plot.addItem(label)
        self._combo_plot.blockSignals(False)
        self._combo_plot.setCurrentIndex(0)
        self._update_desc(row, 0)
        self._update_model_group_visibility()

    # ── Run / Export slots ────────────────────────────────────────────

    def _on_run(self) -> None:
        cat_row = self._combo_category.currentIndex()
        plot_row = self._combo_plot.currentIndex()
        if cat_row < 0 or plot_row < 0:
            return
        _label, plots = ADVANCED_GROUPS[cat_row]
        if plot_row >= len(plots):
            return
        _plot_label, fn_name, has_ax = plots[plot_row]

        if self._ctrl._sites is None:
            self._status_lbl.setText("Load survey data first.")
            return

        self._status_lbl.setText(f"Running {fn_name}…")
        self._btn_run.setEnabled(False)
        try:
            new_fig = self._ctrl.draw(fn_name, has_ax, self._canvas.figure)
            if new_fig is not None:
                self._canvas.show_figure(new_fig)
            else:
                self._canvas.draw()
            self._status_lbl.setText("Done.")
        except Exception as exc:
            self._status_lbl.setText(f"Error: {exc}")
        finally:
            self._btn_run.setEnabled(True)

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )

        ExportDialog(figure=self._canvas.figure, parent=self).exec()

    def _auto_render_if_ready(self) -> None:
        if self._auto_rendered or self._ctrl._sites is None:
            return
        if not self.isVisible():
            return
        cat_row = self._combo_category.currentIndex()
        if cat_row < 0 or cat_row >= len(ADVANCED_GROUPS):
            return
        _label, plots = ADVANCED_GROUPS[cat_row]
        if not plots:
            return
        self._auto_rendered = True
        QTimer.singleShot(0, self._on_run)

    # ── Plot-selection slot ───────────────────────────────────────────

    def _on_plot_changed(self, _: int) -> None:
        cat_row = self._combo_category.currentIndex()
        plot_row = self._combo_plot.currentIndex()
        self._update_desc(cat_row, plot_row)
        self._update_model_group_visibility()

    def _update_model_group_visibility(self) -> None:
        self._grp_model.setVisible(self._is_atom_psection_selected())

    def _is_atom_psection_selected(self) -> bool:
        cat = self._combo_category.currentIndex()
        if cat < 0 or cat >= len(ADVANCED_GROUPS):
            return False
        _, plots = ADVANCED_GROUPS[cat]
        pi = self._combo_plot.currentIndex()
        if pi < 0 or pi >= len(plots):
            return False
        return plots[pi][1] == "plot_atom_psection"

    # ── Train-model slots ─────────────────────────────────────────────

    def _on_train_model(self) -> None:
        if self._ctrl._sites is None:
            self._model_status_lbl.setText("Load survey data first.")
            return
        n_atoms = int(self._spin_n_atoms.value())
        n_iter = int(self._spin_n_iter.value())
        self._btn_train_model.setEnabled(False)
        self._model_status_lbl.setText("Training…")
        self._dim_worker = DimModelWorker(self._ctrl, n_atoms, n_iter)
        self._dim_worker.finished.connect(self._on_model_trained)
        self._dim_worker.error.connect(self._on_model_train_error)
        self._dim_worker.start()

    def _on_model_trained(self, model: dict) -> None:
        self._btn_train_model.setEnabled(True)
        meta = model.get("meta", {})
        n_atoms = model["D"].shape[1] if model.get("D") is not None else "?"
        n_samples = meta.get("samples", "?")
        self._model_status_lbl.setText(
            f"Ready: {n_atoms} atoms · {n_samples} samples"
        )
        self._status_lbl.setText("Model trained. Click Run to plot.")

    def _on_model_train_error(self, msg: str) -> None:
        self._btn_train_model.setEnabled(True)
        self._model_status_lbl.setText(f"Error: {msg}")

    # ── Topo slots ────────────────────────────────────────────────────

    def _on_topo_source_changed(self, text: str) -> None:
        self._topo_file_row.setVisible(text == "file")

    def _browse_topo_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select elevation file",
            "",
            "CSV / Parquet (*.csv *.parquet *.txt);;All Files (*)",
        )
        if path:
            self._edit_topo_file.setText(path)

    def _pick_color(self, which: str) -> None:
        current = (
            self._topo_fill_color if which == "fill" else self._topo_line_color
        )
        color = QColorDialog.getColor(current, self, f"Choose {which} color")
        if color.isValid():
            hex_color = color.name()
            if which == "fill":
                self._topo_fill_color = hex_color
                self._apply_topo_btn_color(
                    self._btn_topo_fill_color, hex_color
                )
            else:
                self._topo_line_color = hex_color
                self._apply_topo_btn_color(
                    self._btn_topo_line_color, hex_color
                )

    def _apply_topo_btn_color(self, btn: QPushButton, hex_color: str) -> None:
        btn.setStyleSheet(f"background:{hex_color}; border:1px solid #888")

    def _on_topo_apply(self) -> None:
        try:
            from pycsamt.topo.config import configure_topo

            configure_topo(
                enabled=self._chk_topo_enabled.isChecked(),
                source=self._combo_topo_source.currentText(),
                elev_file=(self._edit_topo_file.text().strip() or None),
                interp_method=self._combo_topo_interp.currentText(),
                exaggeration=self._spin_topo_exag.value(),
                fill_color=self._topo_fill_color,
                fill_alpha=self._spin_topo_fill_alpha.value(),
                line_color=self._topo_line_color,
                line_width=self._spin_topo_line_w.value(),
                show_surface_line=self._chk_surface_line.isChecked(),
                clip_below_surface=self._chk_clip_below.isChecked(),
                station_pins_at_surface=self._chk_pins_at_surface.isChecked(),
                show_topo_strip=self._chk_topo_strip.isChecked(),
                strip_height_ratio=self._spin_strip_h.value(),
            )
            from pycsamt.topo.config import PYCSAMT_TOPO

            self._topo_status.setText(PYCSAMT_TOPO.summary())
        except Exception as exc:
            self._topo_status.setText(f"Error: {exc}")
        self._refresh_topo_preview()

    def _on_topo_reset(self) -> None:
        try:
            from pycsamt.topo.config import reset_topo

            reset_topo()
            self._sync_topo_widgets_from_config()
            from pycsamt.topo.config import PYCSAMT_TOPO

            self._topo_status.setText("Reset. " + PYCSAMT_TOPO.summary())
        except Exception as exc:
            self._topo_status.setText(f"Error: {exc}")
        self._refresh_topo_preview()

    def _sync_topo_widgets_from_config(self) -> None:
        try:
            from pycsamt.topo.config import PYCSAMT_TOPO as T

            self._chk_topo_enabled.setChecked(T.enabled)
            idx = self._combo_topo_source.findText(T.source)
            if idx >= 0:
                self._combo_topo_source.setCurrentIndex(idx)
            self._topo_file_row.setVisible(T.source == "file")
            if T.elev_file:
                self._edit_topo_file.setText(T.elev_file)
            idx_interp = self._combo_topo_interp.findText(T.interp_method)
            if idx_interp >= 0:
                self._combo_topo_interp.setCurrentIndex(idx_interp)
            self._spin_topo_exag.setValue(T.exaggeration)
            self._topo_fill_color = T.fill_color
            self._apply_topo_btn_color(self._btn_topo_fill_color, T.fill_color)
            self._spin_topo_fill_alpha.setValue(T.fill_alpha)
            self._topo_line_color = T.line_color
            self._apply_topo_btn_color(self._btn_topo_line_color, T.line_color)
            self._spin_topo_line_w.setValue(T.line_width)
            self._chk_surface_line.setChecked(T.show_surface_line)
            self._chk_clip_below.setChecked(T.clip_below_surface)
            self._chk_pins_at_surface.setChecked(T.station_pins_at_surface)
            self._chk_topo_strip.setChecked(T.show_topo_strip)
            self._spin_strip_h.setValue(T.strip_height_ratio)
        except Exception:
            pass

    def _refresh_topo_preview(self) -> None:
        if self._topo_ctrl._sites is None:
            self._topo_stats_lbl.setText("No data loaded")
        else:
            try:
                stats = self._topo_ctrl.get_stats()
                self._topo_stats_lbl.setText(
                    f"{stats['n_stations']} stations  |  "
                    f"elev {stats['elev_min']:.0f}–{stats['elev_max']:.0f} m"
                    if stats["n_stations"] > 0
                    else "No elevation data"
                )
            except Exception:
                self._topo_stats_lbl.setText("")

        view_idx = self._combo_topo_view.currentIndex()
        fig = self._canvas_topo.figure
        try:
            if view_idx == 0:
                self._topo_ctrl.plot_elevation_profile(fig)
            elif view_idx == 1:
                self._topo_ctrl.plot_fill_preview(fig)
            else:
                self._topo_ctrl.plot_elevation_histogram(fig)
        except Exception as exc:
            fig.clear()
            ax = fig.add_subplot(111)
            ax.text(
                0.5,
                0.5,
                f"Preview error: {exc}",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=9,
            )
        self._canvas_topo.draw()

    def _on_topo_view_changed(self, _: int) -> None:
        self._refresh_topo_preview()

    # ── Conversion slots ──────────────────────────────────────────────

    def _on_conv_type_changed(self, idx: int) -> None:
        self._conv_opts_stack.setCurrentIndex(idx)
        if hasattr(self, "_avg_topo_group"):
            self._avg_topo_group.setVisible(idx == 0)
        self._update_conv_run_state()

    def _update_conv_run_state(self) -> None:
        if not hasattr(self, "_btn_conv_run"):
            return
        has_source = bool(
            getattr(self, "_edit_conv_path", None)
            and self._edit_conv_path.text().strip()
        )
        self._btn_conv_run.setEnabled(has_source and not self._conv_running)

    def _browse_conv_path(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select input file", "", "All Files (*)"
        )
        if not path:
            # fallback to directory
            path = QFileDialog.getExistingDirectory(
                self, "Select input directory", ""
            )
        if path:
            self._edit_conv_path.setText(path)

    def _browse_out_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select output directory", ""
        )
        if path:
            self._edit_out_dir.setText(path)

    def _browse_avg_stn_path(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select AVG topography / station profile",
            "",
            "Station files (*.stn *.txt *.csv);;All Files (*)",
        )
        if path:
            self._avg_stn_path.setText(path)

    def _on_conv_run(self) -> None:
        path = self._edit_conv_path.text().strip()
        if not path:
            self._conv_status.setText(
                "Please provide an input file/directory."
            )
            return

        type_str = self._combo_conv_type.currentText()
        self._conv_ctrl.set_source(type_str, path)

        # Collect options from the active options page
        options: dict = {}
        idx = self._combo_conv_type.currentIndex()
        try:
            if idx == 0:  # AVG
                options["freq_order"] = self._avg_freq_order.currentText()
                options["freq_tol"] = self._avg_freq_tol.value()
                if self._avg_compute_z.isChecked():
                    options["compute_z"] = True
                if self._avg_compute_rho.isChecked():
                    options["compute_rho_phi"] = True
                stn_path = self._avg_stn_path.text().strip()
                if stn_path:
                    options["stn_path"] = stn_path
                    options["convert_stn_coords"] = (
                        self._avg_convert_coords.isChecked()
                    )
                epsg = self._avg_epsg.text().strip()
                if epsg:
                    options["epsg"] = epsg
                utm_zone = self._avg_utm_zone.text().strip()
                if utm_zone:
                    options["utm_zone"] = utm_zone
            elif idx == 1:  # J
                options["freq_order"] = self._j_freq_order.currentText()
                station_name = self._j_station_name.text().strip()
                if station_name:
                    options["name"] = station_name
            elif idx == 2:  # Spectra
                options["e_labels"] = self._sp_e_labels.text()
                options["h_labels"] = self._sp_h_labels.text()
                options["estimate_error"] = (
                    self._sp_estimate_errors.isChecked()
                )
                options["use_remote"] = self._sp_remote_ref.isChecked()
                suf = self._sp_station_suffix.text().strip()
                if suf:
                    options["station_suffix"] = suf
                options["skip_errors"] = self._sp_skip_errors.isChecked()
        except Exception:
            pass

        if self._chk_write_edis.isChecked():
            out_dir = self._edit_out_dir.text().strip()
            if out_dir:
                options["output_dir"] = out_dir

        self._conv_progress.setVisible(True)
        self._conv_running = True
        self._update_conv_run_state()
        self._conv_status.setText(f"Running {type_str}…")

        self._conv_worker = ConversionWorker(self._conv_ctrl, options)
        self._conv_worker.finished.connect(self._on_conv_finished)
        self._conv_worker.error.connect(self._on_conv_error)
        self._conv_worker.start()

    def _on_conv_finished(self, collection, failures: list) -> None:
        self._conv_progress.setVisible(False)
        self._conv_running = False
        self._update_conv_run_state()

        self._conv_ctrl._result = collection
        stats = self._conv_ctrl.build_stats(collection, failures)
        self._conv_ctrl._stats = stats
        self._conv_result_stats = stats.get("rows", [])

        # Populate table
        rows_data = stats.get("rows", [])
        cols = [
            "Station",
            "N Freqs",
            "F min (Hz)",
            "F max (Hz)",
            "Latitude",
            "Longitude",
            "Elevation",
            "Has Z",
            "Has Tipper",
        ]
        self._conv_table.setColumnCount(len(cols))
        self._conv_table.setHorizontalHeaderLabels(cols)
        self._conv_table.setRowCount(len(rows_data))

        def _fmt_coord(value, digits=6):
            try:
                value = float(value)
                if value == value:
                    return f"{value:.{digits}f}"
            except Exception:
                pass
            return "—"

        for r, row in enumerate(rows_data):
            self._conv_table.setItem(r, 0, QTableWidgetItem(row["station"]))
            self._conv_table.setItem(
                r, 1, QTableWidgetItem(str(row["n_freqs"]))
            )
            self._conv_table.setItem(
                r,
                2,
                QTableWidgetItem(
                    f"{row['f_min']:.4g}"
                    if row["f_min"] == row["f_min"]
                    else "—"
                ),
            )
            self._conv_table.setItem(
                r,
                3,
                QTableWidgetItem(
                    f"{row['f_max']:.4g}"
                    if row["f_max"] == row["f_max"]
                    else "—"
                ),
            )
            self._conv_table.setItem(
                r, 4, QTableWidgetItem(_fmt_coord(row.get("lat")))
            )
            self._conv_table.setItem(
                r, 5, QTableWidgetItem(_fmt_coord(row.get("lon")))
            )
            self._conv_table.setItem(
                r, 6, QTableWidgetItem(_fmt_coord(row.get("elev"), digits=2))
            )
            self._conv_table.setItem(
                r, 7, QTableWidgetItem("Yes" if row["has_Z"] else "No")
            )
            self._conv_table.setItem(
                r, 8, QTableWidgetItem("Yes" if row["has_tipper"] else "No")
            )

        # Refresh plots
        self._conv_ctrl.plot_impedance_curves(self._canvas_conv_curves.figure)
        self._canvas_conv_curves.draw()
        self._conv_ctrl.plot_station_map(self._canvas_conv_map.figure)
        self._canvas_conv_map.draw()

        n_ok = stats.get("n_total", 0)
        n_fail = stats.get("n_failures", 0)
        self._conv_status.setText(
            f"Done: {n_ok} stations"
            + (f", {n_fail} failures" if n_fail else "")
        )
        self._btn_conv_commit.setEnabled(True)
        self._btn_conv_export.setEnabled(True)

    def _on_conv_error(self, msg: str) -> None:
        self._conv_progress.setVisible(False)
        self._conv_running = False
        self._update_conv_run_state()
        self._conv_status.setText(f"Error: {msg}")

    def _on_conv_commit(self) -> None:
        if self._conv_ctrl.result is not None:
            self.conversion_committed.emit(self._conv_ctrl.result)

    def _on_conv_export(self) -> None:
        if not self._conv_ctrl.has_result:
            return
        out_dir = QFileDialog.getExistingDirectory(
            self, "Select output directory for EDI files", ""
        )
        if not out_dir:
            return
        try:
            from pycsamt.emtools._core import _iter_items

            written = 0
            for ed in _iter_items(self._conv_ctrl.result):
                try:
                    write_fn = getattr(ed, "write_edifile", None) or getattr(
                        ed, "write", None
                    )
                    if write_fn is not None:
                        write_fn(save_dir=out_dir)
                        written += 1
                except Exception:
                    pass
            self._conv_status.setText(
                f"Exported {written} EDI files to {out_dir}"
            )
        except Exception as exc:
            self._conv_status.setText(f"Export error: {exc}")

    def _on_conv_clear(self) -> None:
        self._conv_ctrl._result = None
        self._conv_ctrl._stats = {}
        self._conv_ctrl._failures = []
        self._conv_result_stats = []
        self._conv_table.clearContents()
        self._conv_table.setRowCount(0)
        self._conv_table.setColumnCount(0)
        # Clear canvases
        for canvas in (self._canvas_conv_curves, self._canvas_conv_map):
            canvas.figure.clear()
            canvas.draw()
        self._btn_conv_commit.setEnabled(False)
        self._btn_conv_export.setEnabled(False)
        self._conv_status.setText("Cleared.")
        self._conv_progress.setVisible(False)

    # ── Helpers ───────────────────────────────────────────────────────

    def _populate_category_combo(self) -> None:
        self._combo_category.blockSignals(True)
        for group_label, _plots in ADVANCED_GROUPS:
            icon_name = ADVANCED_GROUP_ICONS.get(group_label, "advanced-tools")
            icon = _icon(icon_name)
            if icon.isNull():
                self._combo_category.addItem(group_label)
            else:
                self._combo_category.addItem(icon, group_label)
        self._combo_category.blockSignals(False)

    def _update_desc(self, cat_row: int, plot_row: int) -> None:
        try:
            _label, plots = ADVANCED_GROUPS[cat_row]
            if not plots:
                self._desc_lbl.setText("")
                return
            plot_label, fn_name, has_ax = plots[plot_row]
            proj = (
                " · polar"
                if fn_name
                in __import__(
                    "pycsamt.app.desktop.controllers.advanced_controller",
                    fromlist=["_POLAR_FNS"],
                )._POLAR_FNS
                else ""
            )
            desc = describe_advanced_plot(fn_name)
            self._desc_lbl.setText(
                f"<b>{plot_label}</b>{proj}<br/>"
                f"<small style='color:#888'>{desc}</small>"
            )
        except Exception:
            self._desc_lbl.setText("")
