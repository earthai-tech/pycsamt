# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
MainWindow — redesigned host window for the pycsamt desktop application.

Philosophy
──────────
The main window is intentionally lean: it owns only the station list and the
station detail card.  Every scientific task (profile curves, map, QC,
inversion, agents) opens as an independent floating window so the user can
arrange them freely — on multiple monitors, side by side, or overlapping.

Layout
──────
    ┌──────────────────────────────────────────────────────────┐
    │  Toolbar:  Open  Load EDI  │  Profile  Map  QC           │
    │            Inversion       │  Agents   │  Export  Theme  │
    ├──────────────────┬───────────────────────────────────────┤
    │  Station List    │  Station Detail Card                  │
    │  (search + table)│  (name, coords, quality, actions)    │
    ├──────────────────┴───────────────────────────────────────┤
    │  Log  (collapsible, 80 px)                               │
    └──────────────────────────────────────────────────────────┘

Menu bar
────────
    File     — Open, Load EDI, Recent Files, Save Session, Quit
    View     — toggle each panel window + theme
    Settings — API Configuration, Reset to Defaults, Save/Load Profile
    Help     — Documentation, About
"""

from __future__ import annotations

import re
from pathlib import Path

from PySide6.QtCore import QByteArray, QSize, Qt
from PySide6.QtGui  import QAction, QActionGroup, QIcon, QKeySequence, QPainter, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QDockWidget,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMenu,
    QProgressBar,
    QSizePolicy,
    QSplitter,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.app_controller import AppController
from pycsamt.app.desktop.agent_master_bridge     import launch_agent_master
from pycsamt.app.desktop.models.session            import SessionState
from pycsamt.app.desktop.panels.log_panel          import LogPanel
from pycsamt.app.desktop.panels.station_panel      import StationPanel
from pycsamt.app.desktop.widgets.mpl_canvas        import (
    apply_mpl_dark_theme,
    apply_mpl_light_theme,
)
from pycsamt.app.desktop.widgets.station_detail    import StationDetailCard
from pycsamt.app.desktop.widgets.survey_overview   import SurveyOverviewWidget
from pycsamt.app.desktop.windows                   import (
    AdvancedToolsWindow,
    CorrectionWindow,
    ForwardModelWindow,
    InterpretationWindow,
    InversionWindow,
    MapViewerWindow,
    ProfileViewerWindow,
    QCDashboardWindow,
    TDEMWindow,
)
from pycsamt.app.desktop.windows.pipeline_window   import PipelineWindow

_RESOURCES = Path(__file__).parent / "resources"
_ICONS     = _RESOURCES / "icons"

# Tracks current theme so _icon() can recolor SVGs without being passed a flag.
_DARK_MODE: bool = False

# Matches any 6-digit hex colour in an SVG source string.
_HEX6_RE = re.compile(r'#[0-9a-fA-F]{6}', re.IGNORECASE)


def _is_near_black(hex6: str) -> bool:
    r, g, b = int(hex6[1:3], 16), int(hex6[3:5], 16), int(hex6[5:7], 16)
    return (0.299 * r + 0.587 * g + 0.114 * b) < 80


def _recolor_svg(text: str, target: str = "#cdd6f4") -> bytes:
    # Pass 1: replace near-black hex colours
    result = _HEX6_RE.sub(
        lambda m: target if _is_near_black(m.group(0)) else m.group(0),
        text,
    )
    # Pass 2: replace the named keyword "black" in attribute values / inline styles
    result = re.sub(r'(?<=[";\s:])black(?=[";\s])', target, result, flags=re.IGNORECASE)
    # Pass 3: inject fill on individual shape elements that carry no explicit fill,
    # so icons that use the SVG default (black) fill — e.g. interpret.svg — become
    # visible on a dark background.  Elements with an existing fill= or fill: are
    # left untouched, and CSS class rules (fill:none) still override the injected
    # presentation attribute via normal specificity.
    def _maybe_fill(m: re.Match) -> str:
        tag = m.group(0)
        if "fill=" not in tag and "fill:" not in tag:
            return (tag[:-2] + f' fill="{target}"/>'
                    if tag.endswith("/>")
                    else tag[:-1] + f' fill="{target}">')
        return tag
    result = re.sub(
        r'<(?:path|rect|circle|ellipse|polygon|polyline|line)\b[^>]*?(?:/>|>)',
        _maybe_fill,
        result,
    )
    return result.encode("utf-8")


def _icon(name: str) -> QIcon:
    """Load an icon by name, recolouring dark SVG strokes for dark mode."""
    for c in (name, f"{name}.svg", f"{name}.png"):
        p = _ICONS / c
        if p.exists():
            if _DARK_MODE and str(c).endswith(".svg"):
                try:
                    from PySide6.QtSvg import QSvgRenderer
                    svg_bytes = _recolor_svg(p.read_text("utf-8"))
                    renderer  = QSvgRenderer(QByteArray(svg_bytes))
                    icon = QIcon()
                    for s in (16, 22, 32, 48):
                        pm = QPixmap(s, s)
                        pm.fill(Qt.GlobalColor.transparent)
                        painter = QPainter(pm)
                        renderer.render(painter)
                        painter.end()
                        icon.addPixmap(pm)
                    return icon
                except Exception:
                    pass  # fall through to plain load
            return QIcon(str(p))
    return QIcon()


# ── Main window ───────────────────────────────────────────────────────────────

class MainWindow(QMainWindow):
    """Lean host window — station list + detail card + panel launchers."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._session    = SessionState.load()
        self._controller = AppController(self._session)
        self._loader     = None      # LoaderWorker kept alive while running
        self._recomputed_ids: set[str] = set()   # stations marked recomputed
        self._last_recompute_output = None        # Path to last recomputed EDI folder

        # Settings controller — load saved profile (silently no-ops if absent)
        from pycsamt.app.desktop.controllers.settings_controller import SettingsController
        self._settings_ctrl = SettingsController()
        self._settings_ctrl.load()

        self._setup_window()
        self._create_panel_windows()          # windows must exist before theme is applied
        self._apply_theme(self._session.theme)
        self._create_central_widget()
        self._create_log_dock()
        self._create_menu_bar()
        self._create_tool_bar()
        self._create_status_bar()
        self._restore_layout()
        self._wire_signals()
        self._log("pycsamt ready — load EDI files to begin.")

    # ── Window setup ──────────────────────────────────────────────────

    def _setup_window(self) -> None:
        self.setWindowTitle("pycsamt")
        self.setMinimumSize(1000, 600)
        self.resize(1200, 720)
        ico = _ICONS / "pycsamt.ico"
        if ico.exists():
            self.setWindowIcon(QIcon(str(ico)))
        # Tracks every icon-bearing QAction so _apply_theme can re-ice them all.
        # Populated by _create_menu_bar and _create_tool_bar (both run after _setup_window).
        self._all_icon_actions: list = []

    # ── Theme ──────────────────────────────────────────────────────────

    def _apply_theme(self, theme: str) -> None:
        global _DARK_MODE
        _DARK_MODE = (theme == "dark")

        qss = _RESOURCES / f"{theme}_theme.qss"
        if qss.exists():
            QApplication.instance().setStyleSheet(
                qss.read_text(encoding="utf-8")
            )
        (apply_mpl_dark_theme if theme == "dark" else apply_mpl_light_theme)()
        self._session.theme = theme

        # Sync the View > Theme checkboxes (created after first call)
        if hasattr(self, "_act_dark"):
            self._act_dark.setChecked(theme == "dark")
            self._act_light.setChecked(theme == "light")

        # Re-apply ALL icon-bearing actions (menu + toolbar + More menu)
        for action, icon_name in getattr(self, "_all_icon_actions", []):
            action.setIcon(_icon(icon_name))
        if hasattr(self, "_more_btn"):
            self._more_btn.setIcon(_icon("more"))

        for win in self._panel_windows():
            try:
                win.set_dark_mode(theme == "dark")
            except Exception:
                pass


    def _toggle_theme(self) -> None:
        new = "light" if self._session.theme == "dark" else "dark"
        self._apply_theme(new)
        # Show the next target: now in *new*, so label points to the other theme.
        self._act_theme.setText("☾  Dark" if new == "light" else "☀  Light")

    # ── Independent panel windows (created once, shown/hidden on demand) ──

    def _create_panel_windows(self) -> None:
        self._profile_win    = ProfileViewerWindow(parent=self)
        self._map_win        = MapViewerWindow(parent=self)
        self._qc_win         = QCDashboardWindow(parent=self)
        self._correction_win = CorrectionWindow(parent=self)
        self._advanced_win   = AdvancedToolsWindow(parent=self)
        self._tdem_win       = TDEMWindow(parent=self)
        self._pipeline_win   = PipelineWindow(parent=self)
        self._forward_win    = ForwardModelWindow(parent=self)
        self._inversion_win  = InversionWindow(parent=self)
        self._interp_win     = InterpretationWindow(parent=self)

        # Wire forward → inversion bridge
        self._forward_win.send_to_inversion.connect(self._on_forward_send_to_inversion)
        # Wire inversion result → interpretation model
        self._inversion_win.result_ready.connect(self._on_inversion_result_ready)

        # Wire correction window → main data
        self._correction_win.corrections_committed.connect(self._on_corrections_committed)
        # Wire advanced tools conversion → main data
        self._advanced_win.conversion_committed.connect(self._on_conversion_committed)

        # Restore positions from session
        geo = self._session.window_geometries
        for win in self._panel_windows():
            win.restore_geometry_from(geo)

    # ── Central widget — station list + detail card ───────────────────

    def _create_central_widget(self) -> None:
        container = QWidget(self)
        container.setObjectName("CentralContainer")
        h_layout = QVBoxLayout(container)
        h_layout.setContentsMargins(0, 0, 0, 0)
        h_layout.setSpacing(0)

        splitter = QSplitter(Qt.Orientation.Horizontal, container)
        splitter.setHandleWidth(3)

        # ── Left: station list ────────────────────────────────────────
        left = QWidget()
        left.setObjectName("StationListPane")
        left.setMinimumWidth(220)
        left.setMaximumWidth(400)
        left_v = QVBoxLayout(left)
        left_v.setContentsMargins(4, 4, 4, 4)
        left_v.setSpacing(4)

        # Search / filter bar
        self._search_bar = QLineEdit()
        self._search_bar.setPlaceholderText("🔍  Filter stations…")
        self._search_bar.setObjectName("SearchBar")
        self._search_bar.textChanged.connect(self._on_filter_changed)
        left_v.addWidget(self._search_bar)

        # Station table
        self._station_panel = StationPanel(left)
        left_v.addWidget(self._station_panel)

        # ── Right: vertical splitter — survey overview + detail card ──
        right_splitter = QSplitter(Qt.Orientation.Vertical, container)
        right_splitter.setHandleWidth(3)

        self._survey_overview = SurveyOverviewWidget(container)
        self._survey_overview.setObjectName("SurveyOverview")

        self._detail_card = StationDetailCard(container)
        self._detail_card.setObjectName("StationDetailCard")

        right_splitter.addWidget(self._survey_overview)
        right_splitter.addWidget(self._detail_card)
        right_splitter.setStretchFactor(0, 0)
        right_splitter.setStretchFactor(1, 1)
        right_splitter.setSizes([220, 360])

        splitter.addWidget(left)
        splitter.addWidget(right_splitter)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([280, 720])

        h_layout.addWidget(splitter)
        self.setCentralWidget(container)

    # ── Log dock (bottom, thin) ───────────────────────────────────────

    def _create_log_dock(self) -> None:
        self._log_panel = LogPanel(self)
        self._log_dock  = QDockWidget("Log", self)
        self._log_dock.setObjectName("LogDock")
        self._log_dock.setAllowedAreas(Qt.DockWidgetArea.BottomDockWidgetArea)
        self._log_dock.setWidget(self._log_panel)
        self._log_dock.setMaximumHeight(120)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self._log_dock)

    # ── Menu bar (File | View | Help) ─────────────────────────────────

    def _create_menu_bar(self) -> None:
        mb = self.menuBar()

        # ── File ──────────────────────────────────────────────────────
        file_menu = mb.addMenu("&File")

        self._act_open = QAction(_icon("open"), "&Open / Load EDI…", self)
        self._act_open.setShortcut(QKeySequence.StandardKey.Open)
        self._act_open.setStatusTip("Open EDI / AVG / J survey files")
        self._act_open.triggered.connect(self._open_files)

        self._act_save = QAction(_icon("save-session"), "&Save Session", self)
        self._act_save.setShortcut(QKeySequence.StandardKey.Save)
        self._act_save.triggered.connect(self._on_save_session)

        self._act_recent = file_menu.addMenu("Recent Files")
        self._rebuild_recent_menu()

        act_quit = QAction(_icon("quit"), "&Quit", self)
        act_quit.setShortcut(QKeySequence.StandardKey.Quit)
        act_quit.triggered.connect(self.close)

        file_menu.addAction(self._act_open)
        file_menu.addAction(self._act_save)
        file_menu.addSeparator()
        file_menu.addMenu(self._act_recent)
        file_menu.addSeparator()
        act_prefs = QAction(_icon("tools"), "&Preferences…", self)
        act_prefs.setShortcut(QKeySequence.StandardKey.Preferences)
        act_prefs.triggered.connect(self._open_preferences)
        file_menu.addAction(act_prefs)
        file_menu.addSeparator()
        file_menu.addAction(act_quit)

        # ── View ──────────────────────────────────────────────────────
        view_menu = mb.addMenu("&View")

        act_profile = QAction(_icon("profile-view"), "&Profile Viewer", self)
        act_profile.setShortcut("Ctrl+P")
        act_profile.triggered.connect(lambda: self._show_window(self._profile_win))
        view_menu.addAction(act_profile)

        act_map = QAction(_icon("map-view"), "&Map Viewer", self)
        act_map.setShortcut("Ctrl+M")
        act_map.triggered.connect(lambda: self._show_window(self._map_win))
        view_menu.addAction(act_map)

        act_qc = QAction(_icon("qc"), "&QC Dashboard", self)
        act_qc.setShortcut("Ctrl+Q")
        act_qc.triggered.connect(lambda: self._show_window(self._qc_win))
        view_menu.addAction(act_qc)

        act_corr = QAction(_icon("sites-correction"), "&Data Corrections", self)
        act_corr.setShortcut("Ctrl+R")
        act_corr.triggered.connect(lambda: self._show_window(self._correction_win))
        view_menu.addAction(act_corr)

        act_fwd = QAction(_icon("forward"), "&Forward Modelling", self)
        act_fwd.setShortcut("Ctrl+F")
        act_fwd.triggered.connect(lambda: self._show_window(self._forward_win))
        view_menu.addAction(act_fwd)

        act_inv = QAction(_icon("inversion"), "&Inversion Wizard…", self)
        act_inv.setShortcut("Ctrl+I")
        act_inv.triggered.connect(self._open_inversion_wizard)
        view_menu.addAction(act_inv)

        act_interp = QAction(_icon("interpret"), "&Interpretation Studio", self)
        act_interp.setShortcut("Ctrl+Shift+I")
        act_interp.triggered.connect(lambda: self._show_window(self._interp_win))
        view_menu.addAction(act_interp)

        act_pipe = QAction(_icon("pipeline"), "&Processing Pipeline", self)
        act_pipe.setShortcut("Ctrl+Shift+P")
        act_pipe.triggered.connect(lambda: self._show_window(self._pipeline_win))
        view_menu.addAction(act_pipe)

        act_tdem = QAction(_icon("tdem"), "&TDEM Analysis", self)
        act_tdem.setShortcut("Ctrl+T")
        act_tdem.triggered.connect(lambda: self._show_window(self._tdem_win))
        view_menu.addAction(act_tdem)

        act_agents = QAction(_icon("agents"), "&Agent Master", self)
        act_agents.setShortcut("Ctrl+Shift+A")
        act_agents.triggered.connect(self._open_agent_master)
        view_menu.addAction(act_agents)

        act_adv = QAction(_icon("advanced-tools"), "&Advanced Tools", self)
        act_adv.setShortcut("Ctrl+Shift+T")
        act_adv.triggered.connect(lambda: self._show_window(self._advanced_win))
        view_menu.addAction(act_adv)

        view_menu.addSeparator()
        _log_tva = self._log_dock.toggleViewAction()
        _log_tva.setIcon(_icon("log"))
        view_menu.addAction(_log_tva)
        view_menu.addSeparator()

        theme_menu = view_menu.addMenu("Theme")
        self._act_dark  = QAction("☾  Dark",  self, checkable=True)
        self._act_light = QAction("☀  Light", self, checkable=True)
        # QActionGroup enforces mutual exclusivity: checking one unchecks the other.
        _theme_grp = QActionGroup(self)
        _theme_grp.setExclusive(True)
        _theme_grp.addAction(self._act_dark)
        _theme_grp.addAction(self._act_light)
        self._act_dark.setChecked(self._session.theme  == "dark")
        self._act_light.setChecked(self._session.theme == "light")
        self._act_dark.triggered.connect(lambda: self._apply_theme("dark"))
        self._act_light.triggered.connect(lambda: self._apply_theme("light"))
        theme_menu.addAction(self._act_dark)
        theme_menu.addAction(self._act_light)

        # ── Tools ─────────────────────────────────────────────────────
        tools_menu = mb.addMenu("&Tools")

        act_strike = QAction(_icon("strike-analyzer"), "&Strike Analyzer…", self)
        act_strike.setShortcut("Ctrl+Shift+S")
        act_strike.setStatusTip("Regional-strike and dimensionality analysis")
        act_strike.triggered.connect(self._open_strike_analyzer)
        tools_menu.addAction(act_strike)

        act_valid = QAction(_icon("EDI-validator"), "ED&I Validator…", self)
        act_valid.setShortcut("Ctrl+Shift+V")
        act_valid.setStatusTip("Per-station data-quality checklist")
        act_valid.triggered.connect(self._open_edi_validator)
        tools_menu.addAction(act_valid)

        tools_menu.addSeparator()

        act_recompute = QAction(_icon("recompute"), "&Recompute EDIs…", self)
        act_recompute.setShortcut("Ctrl+Shift+X")
        act_recompute.setStatusTip(
            "Recompute EDI files: rotate, filter frequencies, fill missing, rewrite"
        )
        act_recompute.triggered.connect(self._open_recompute)
        tools_menu.addAction(act_recompute)

        tools_menu.addSeparator()

        act_conv = QAction(_icon("format-converter"), "&Format Converter…", self)
        act_conv.setShortcut("Ctrl+Shift+C")
        act_conv.setStatusTip("Export survey to EDI / CSV / JSON")
        act_conv.triggered.connect(self._open_format_converter)
        tools_menu.addAction(act_conv)

        act_batch = QAction(_icon("batch-export"), "&Batch Export Plots…", self)
        act_batch.setShortcut("Ctrl+Shift+E")
        act_batch.setStatusTip("Save every open canvas figure to a folder")
        act_batch.triggered.connect(self._open_batch_export)
        tools_menu.addAction(act_batch)

        tools_menu.addSeparator()

        act_coord = QAction(_icon("coordinate-transformer"), "&Coordinate Transformer…", self)
        act_coord.setShortcut("Ctrl+Shift+G")
        act_coord.setStatusTip("Convert UTM ↔ Lat/Lon (pyproj-backed)")
        act_coord.triggered.connect(self._open_coord_transformer)
        tools_menu.addAction(act_coord)

        tools_menu.addSeparator()

        act_station_resp = QAction(_icon("station-response"), "Station &Response Inspector…", self)
        act_station_resp.setShortcut("Ctrl+Shift+R")
        act_station_resp.setStatusTip("Plot full impedance tensor response for a selected station")
        act_station_resp.triggered.connect(self._open_station_response)
        tools_menu.addAction(act_station_resp)

        act_strike_profile = QAction(_icon("strike-profile"), "Strike &Profile Viewer…", self)
        act_strike_profile.setShortcut("Ctrl+Shift+P")
        act_strike_profile.setStatusTip("Strike angle vs. station-position line plot with IQR ribbon")
        act_strike_profile.triggered.connect(self._open_strike_profile)
        tools_menu.addAction(act_strike_profile)

        act_pt_map = QAction(_icon("phase-tensor"), "Phase &Tensor Map…", self)
        act_pt_map.setShortcut("Ctrl+Shift+T")
        act_pt_map.setStatusTip("Geographic map of phase-tensor ellipses at a chosen period")
        act_pt_map.triggered.connect(self._open_phase_tensor_map)
        tools_menu.addAction(act_pt_map)

        tools_menu.addSeparator()

        act_dim = QAction(_icon("dimensionnality"), "&Dimensionality Classifier…", self)
        act_dim.setShortcut("Ctrl+Shift+D")
        act_dim.setStatusTip("Classify each station × frequency as 1D / 2D / 3D")
        act_dim.triggered.connect(self._open_dimensionality)
        tools_menu.addAction(act_dim)

        act_freq_ed = QAction(_icon("frequency-editor"), "&Frequency Editor…", self)
        act_freq_ed.setShortcut("Ctrl+Shift+F")
        act_freq_ed.setStatusTip("Confidence-based frequency QC: drop / mask / recover frequency bands")
        act_freq_ed.triggered.connect(self._open_frequency_editor)
        tools_menu.addAction(act_freq_ed)

        tools_menu.addSeparator()

        act_lm = QAction(_icon("layered-model"), "&Layered Model Builder…", self)
        act_lm.setShortcut("Ctrl+Shift+L")
        act_lm.setStatusTip("Build and preview a 1-D layered earth model")
        act_lm.triggered.connect(self._open_layered_model)
        tools_menu.addAction(act_lm)

        act_elev = QAction(_icon("elevation"), "&Elevation Enrichment…", self)
        act_elev.setShortcut("Ctrl+Shift+H")
        act_elev.setStatusTip("Fetch elevation for all loaded stations via open elevation API")
        act_elev.triggered.connect(self._open_elevation_enrichment)
        tools_menu.addAction(act_elev)

        self._all_icon_actions.extend([
            (act_strike,         "strike-analyzer"),
            (act_valid,          "EDI-validator"),
            (act_recompute,      "recompute"),
            (act_conv,           "format-converter"),
            (act_batch,          "batch-export"),
            (act_coord,          "coordinate-transformer"),
            (act_station_resp,   "station-response"),
            (act_strike_profile, "strike-profile"),
            (act_pt_map,         "phase-tensor"),
            (act_dim,            "dimensionnality"),
            (act_freq_ed,        "frequency-editor"),
            (act_lm,             "layered-model"),
            (act_elev,           "elevation"),
        ])

        # ── Settings ──────────────────────────────────────────────────
        settings_menu = mb.addMenu("&Settings")

        act_api = QAction(_icon("tools"), "&API Configuration…", self)
        act_api.setShortcut("Ctrl+,")
        act_api.setStatusTip(
            "Configure PYCSAMT_* singletons: pseudosections, view controls, display, topography"
        )
        act_api.triggered.connect(self._open_api_config)
        settings_menu.addAction(act_api)

        settings_menu.addSeparator()

        act_reset_all = QAction("Reset All to Defaults", self)
        act_reset_all.setStatusTip("Reset all API singletons to package defaults")
        act_reset_all.triggered.connect(self._reset_all_settings)
        settings_menu.addAction(act_reset_all)

        settings_menu.addSeparator()

        act_save_profile = QAction("Save Profile…", self)
        act_save_profile.setStatusTip("Save current API configuration to a JSON profile")
        act_save_profile.triggered.connect(self._save_settings_profile)
        settings_menu.addAction(act_save_profile)

        act_load_profile = QAction("Load Profile…", self)
        act_load_profile.setStatusTip("Load an API configuration profile from JSON")
        act_load_profile.triggered.connect(self._load_settings_profile)
        settings_menu.addAction(act_load_profile)

        self._all_icon_actions.append((act_api, "tools"))

        # ── Help ──────────────────────────────────────────────────────
        help_menu = mb.addMenu("&Help")
        act_docs  = QAction(_icon("docs"), "&Documentation", self)
        act_docs.setStatusTip("Open pycsamt documentation in your browser")
        act_docs.triggered.connect(self._open_documentation)
        help_menu.addAction(act_docs)

        act_gh = QAction(_icon("github"), "pycsamt on &GitHub", self)
        act_gh.setStatusTip("Open the pycsamt GitHub repository in your browser")
        act_gh.triggered.connect(self._open_github)
        help_menu.addAction(act_gh)

        help_menu.addSeparator()

        act_about = QAction(_icon("help"), "&About pycsamt", self)
        act_about.setStatusTip("About pycsamt v2 — version, author and links")
        act_about.triggered.connect(self._open_about)
        help_menu.addAction(act_about)

        # Register every icon-bearing action created in this method so that
        # _apply_theme can re-ice all of them when the user switches theme.
        self._all_icon_actions.extend([
            (self._act_open,  "open"),
            (self._act_save,  "save-session"),
            (act_prefs,       "tools"),
            (act_profile,     "profile-view"),
            (act_map,         "map-view"),
            (act_qc,          "qc"),
            (act_corr,        "sites-correction"),
            (act_fwd,         "forward"),
            (act_inv,         "inversion"),
            (act_interp,      "interpret"),
            (act_pipe,        "pipeline"),
            (act_tdem,        "tdem"),
            (act_agents,      "agents"),
            (act_adv,         "advanced-tools"),
            (_log_tva,        "log"),
            (act_docs,        "docs"),
            (act_gh,          "github"),
            (act_about,       "help"),
        ])

    # ── Toolbar ───────────────────────────────────────────────────────

    def _create_tool_bar(self) -> None:
        tb = self.addToolBar("Main")
        tb.setObjectName("MainToolBar")
        tb.setIconSize(QSize(22, 22))
        tb.setMovable(False)
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)

        tb.addAction(self._act_open)
        tb.addAction(self._act_save)
        tb.addSeparator()

        def _tb(icon_name, label, tip, slot):
            a = QAction(_icon(icon_name), label, self)
            a.setStatusTip(tip)
            a.triggered.connect(slot)
            tb.addAction(a)
            self._all_icon_actions.append((a, icon_name))
            return a

        # ── Primary workflow — max 10 visible actions (+ More) ───────────
        # Slots 1-2: Open / Save already added above.
        # Slots 3-10: the eight core scientific tools.
        _tb("profile-view",     "Profile",     "Open Profile Viewer",
            lambda: self._show_window(self._profile_win))
        _tb("map-view",         "Map",         "Open Map Viewer",
            lambda: self._show_window(self._map_win))
        tb.addSeparator()
        _tb("qc",               "QC",          "Open QC Dashboard",
            lambda: self._show_window(self._qc_win))
        _tb("sites-correction", "Corrections", "Open Data Corrections",
            lambda: self._show_window(self._correction_win))
        _tb("forward",          "Forward",     "Open Forward Modelling",
            lambda: self._show_window(self._forward_win))
        _tb("inversion",        "Inversion",   "Open Inversion Wizard",
            self._open_inversion_wizard)
        _tb("interpret",        "Interpret",   "Open Interpretation Studio",
            lambda: self._show_window(self._interp_win))
        _tb("pipeline",         "Pipeline",    "Open Processing Pipeline",
            lambda: self._show_window(self._pipeline_win))
        tb.addSeparator()
        _tb("agents",           "Agents",      "Open Agent Master",
            self._open_agent_master)

        # ── Secondary-tools dropdown — sits directly next to Agents ──────────
        # One click opens the menu; items are grouped into three sections so
        # the user can find them at a glance without hunting through a flat list.
        # Stored as instance variables to prevent Python GC from collecting them.
        self._more_menu = QMenu(self)
        sec_menu = self._more_menu

        # ── Section 1: time-domain & frequency-domain extensions ─────────────
        sec_menu.addSection("EM Extensions")
        _act_tdem = QAction(_icon("tdem"), "TDEM Analysis", self)
        _act_tdem.setStatusTip("Open time-domain EM (TDEM) analysis panel")
        _act_tdem.triggered.connect(lambda: self._show_window(self._tdem_win))
        sec_menu.addAction(_act_tdem)

        # ── Section 2: advanced processing / diagnostics ─────────────────────
        sec_menu.addSection("Processing")
        _act_adv = QAction(_icon("advanced-tools"), "Advanced Tools", self)
        _act_adv.setStatusTip("Open advanced EM processing and diagnostics")
        _act_adv.triggered.connect(lambda: self._show_window(self._advanced_win))
        sec_menu.addAction(_act_adv)

        # ── Section 3: output ─────────────────────────────────────────────────
        sec_menu.addSection("Output")
        _act_exp = QAction(_icon("export"), "Export Figure…", self)
        _act_exp.setStatusTip("Save the active figure to PNG / PDF / SVG")
        _act_exp.triggered.connect(self._on_export_figure)
        sec_menu.addAction(_act_exp)

        # Track More-menu actions so _apply_theme can re-ice them on theme switch
        self._all_icon_actions.extend([
            (_act_tdem,  "tdem"),
            (_act_adv,   "advanced-tools"),
            (_act_exp,   "export"),
        ])

        self._more_btn = QToolButton(self)
        self._more_btn.setObjectName("MoreToolButton")
        self._more_btn.setText("More ▾")
        self._more_btn.setIcon(_icon("more"))
        self._more_btn.setIconSize(QSize(22, 22))
        self._more_btn.setToolTip("TDEM · Advanced Tools · Export")
        self._more_btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        self._more_btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self._more_btn.setMenu(self._more_menu)
        self._more_btn.setMinimumWidth(58)
        tb.addWidget(self._more_btn)
        tb.addSeparator()

        # Label shows the TARGET theme (what you get when you click),
        # not the current theme — so in dark mode it reads "☀ Light", etc.
        self._act_theme = QAction(
            "☀  Light" if self._session.theme == "dark" else "☾  Dark", self
        )
        self._act_theme.setStatusTip("Toggle dark / light theme")
        self._act_theme.triggered.connect(self._toggle_theme)
        tb.addAction(self._act_theme)

    # ── Status bar ────────────────────────────────────────────────────

    def _create_status_bar(self) -> None:
        sb = self.statusBar()
        self._status_file_lbl = QLabel("No data loaded")
        self._status_freq_lbl = QLabel("")
        self._status_freq_lbl.setObjectName("StatusFreqLabel")
        self._status_ready_lbl = QLabel("Ready  ●")
        self._status_ready_lbl.setObjectName("StatusReadyLabel")
        self._progress_bar = QProgressBar()
        self._progress_bar.setFixedWidth(160)
        self._progress_bar.setMaximumHeight(14)
        self._progress_bar.setRange(0, 100)
        self._progress_bar.setVisible(False)
        sb.addWidget(self._status_file_lbl)
        sb.addWidget(self._status_freq_lbl)
        sb.addPermanentWidget(self._progress_bar)
        sb.addPermanentWidget(self._status_ready_lbl)

    # ── Signal wiring ─────────────────────────────────────────────────

    def _wire_signals(self) -> None:
        ctrl = self._controller

        # Data loaded → update everything
        ctrl.on_data_loaded(self._on_data_loaded)

        # Status messages
        ctrl.on_status_message(
            lambda msg: self.statusBar().showMessage(msg, 4000)
        )

        # Station selected → detail card + profile window (if open)
        ctrl.on_station_selected(self._on_station_selected)

        # Station panel clicks → controller
        self._station_panel.station_selected.connect(ctrl.select_station)
        # Double-click on station → open Profile Viewer for that station
        self._station_panel.station_selected.connect(self._on_station_double_clicked)

        # Detail card action buttons
        self._detail_card.open_profile_requested.connect(
            self._on_open_profile_for_station
        )
        self._detail_card.show_on_map_requested.connect(
            self._on_show_on_map
        )

        # Map window station pick → controller
        self._map_win.station_selected.connect(ctrl.select_station)

    # ── Data loading ──────────────────────────────────────────────────

    def _open_files(self) -> None:
        from pycsamt.app.desktop.dialogs.load_data_dlg import LoadDataDialog
        dlg = LoadDataDialog(
            self,
            last_dir=self._session.last_data_dir,
            recomputed_dir=self._last_recompute_output,
        )
        if dlg.exec() != LoadDataDialog.DialogCode.Accepted:
            return
        paths = dlg.selected_paths
        if not paths:
            return
        self._session.last_data_dir = str(Path(paths[0]).parent)
        for p in paths:
            self._controller.add_recent_file(p)
        self._rebuild_recent_menu()
        self._start_loading(paths)

    def _start_loading(self, paths: list) -> None:
        from pycsamt.app.desktop.workers.loader_worker import LoaderWorker
        self._loaded_paths = list(paths)
        self._loader = LoaderWorker(paths, parent=self)
        self._loader.progress.connect(self._progress_bar.setValue)
        self._loader.finished.connect(self._controller.set_sites)
        self._loader.error.connect(self._on_load_error)
        self._progress_bar.setValue(0)
        self._progress_bar.setVisible(True)
        self._status_ready_lbl.setText("Loading…")
        self._log(f"Loading {len(paths)} file(s)…")
        self._loader.start()

    def _on_data_loaded(self, sites) -> None:
        # Fresh load → clear the local set; StationModel.set_dataframe() will
        # clear its own _recomputed_ids automatically.
        self._recomputed_ids.clear()

        ctrl = self._controller
        df   = None
        if self._loader is not None:
            df = self._loader.data_controller.dataframe

        try:
            if df is not None:
                self._station_panel.set_dataframe(df)
                self._survey_overview.update_survey(
                    df, getattr(self, "_loaded_paths", None)
                )
                # Feed all panel windows
                for _call in [
                    lambda: self._correction_win.set_sites(sites),
                    lambda: self._profile_win.set_sites(sites),
                    lambda: self._map_win.set_dataframe(df),
                    lambda: self._map_win.set_sites(sites),
                    lambda: self._qc_win.set_sites(sites),
                    lambda: self._advanced_win.set_sites(sites),
                    lambda: self._pipeline_win.set_input_sites(sites),
                    lambda: self._forward_win.set_observed_sites(sites),
                    lambda: self._inversion_win.set_sites(sites),
                    lambda: self._interp_win.set_sites(sites),
                ]:
                    try:
                        _call()
                    except Exception:
                        pass   # never let a panel window block the status update
        finally:
            # Always update the status bar — even if a panel window throws.
            n = ctrl.n_stations
            self._status_file_lbl.setText(f"{n} stations loaded")
            self._progress_bar.setVisible(False)
            self._status_ready_lbl.setText("Ready  ●")
            self._log(f"Loaded {n} stations.")

    # ── Station selection ──────────────────────────────────────────────

    def _on_station_selected(self, station_id: str) -> None:
        """Update detail card; propagate to open windows."""
        sites = self._controller.sites
        if sites and station_id:
            self._detail_card.update_station(station_id, sites)
        if self._profile_win.isVisible():
            self._profile_win.set_station(station_id)
        if self._map_win.isVisible():
            self._map_win.highlight_station(station_id)

    def _on_station_double_clicked(self, station_id: str) -> None:
        """Open Profile Viewer and jump to that station."""
        self._on_open_profile_for_station(station_id)

    def _on_open_profile_for_station(self, station_id: str) -> None:
        self._show_window(self._profile_win)
        self._profile_win.set_station(station_id)

    def _on_show_on_map(self, station_id: str) -> None:
        self._show_window(self._map_win)
        self._map_win.highlight_station(station_id)

    # ── Panel window helpers ───────────────────────────────────────────

    def _show_window(self, win: QWidget) -> None:
        """Show, raise and activate an independent panel window."""
        if not win.isVisible():
            # Default position: offset from main window
            geo = self.geometry()
            win.move(geo.left() + 60, geo.top() + 60)
        win.show()
        win.raise_()
        win.activateWindow()

    def _open_agent_master(self) -> None:
        """Launch Agent Master in the user's default browser."""
        try:
            result = launch_agent_master(open_browser=True)
        except Exception as exc:  # noqa: BLE001 - surface GUI errors
            msg = f"Could not launch Agent Master: {exc}"
            self._log(msg)
            self.statusBar().showMessage(msg, 7000)
            return

        state = "Starting" if result.started else "Opening"
        msg = f"{state} Agent Master at {result.url}"
        self._log(msg)
        self.statusBar().showMessage(msg, 6000)

    def _panel_windows(self) -> list:
        wins = []
        for attr in ("_profile_win", "_map_win", "_qc_win",
                     "_correction_win",
                     "_advanced_win", "_tdem_win",
                     "_pipeline_win", "_forward_win", "_inversion_win",
                     "_interp_win"):
            w = getattr(self, attr, None)
            if w is not None:
                wins.append(w)
        return wins

    # ── Correction → Main data bridge ────────────────────────────────

    def _on_corrections_committed(self, corrected_sites) -> None:
        """Replace global sites with the corrected dataset from CorrectionWindow."""
        self._controller.set_sites(corrected_sites)
        # Rebuild the dataframe view from the corrected sites
        from pycsamt.app.desktop.controllers.data_controller import DataController
        dc = DataController()
        dc._sites = corrected_sites
        dc._df    = dc._build_dataframe()
        df = dc.dataframe
        if df is not None and not df.empty:
            self._station_panel.set_dataframe(df)
            self._survey_overview.update_survey(df, getattr(self, "_loaded_paths", None))
            self._profile_win.set_sites(corrected_sites)
            self._map_win.set_dataframe(df)
            self._map_win.set_sites(corrected_sites)
            self._qc_win.set_sites(corrected_sites)
            self._advanced_win.set_sites(corrected_sites)
            self._pipeline_win.set_input_sites(corrected_sites)
            self._forward_win.set_observed_sites(corrected_sites)
            self._inversion_win.set_sites(corrected_sites)
            self._interp_win.set_sites(corrected_sites)
        n = self._controller.n_stations
        self._status_file_lbl.setText(f"{n} stations (corrected)")
        self._log(f"Corrected dataset committed — {n} stations.")

    # ── Conversion → Main data bridge ────────────────────────────────

    def _on_conversion_committed(self, converted_sites) -> None:
        """Replace global sites with an EDICollection converted by AdvancedToolsWindow."""
        try:
            from pycsamt.site.base import to_sites
            converted_sites = to_sites(converted_sites)
        except Exception as exc:
            self._log(f"Could not commit converted dataset: {exc}")
            return

        self._controller.set_sites(converted_sites)
        from pycsamt.app.desktop.controllers.data_controller import DataController
        dc = DataController()
        dc._sites = converted_sites
        dc._df    = dc._build_dataframe()
        df = dc.dataframe
        if df is not None and not df.empty:
            self._station_panel.set_dataframe(df)
            self._survey_overview.update_survey(df, getattr(self, "_loaded_paths", None))
            self._correction_win.set_sites(converted_sites)
            self._profile_win.set_sites(converted_sites)
            self._map_win.set_dataframe(df)
            self._map_win.set_sites(converted_sites)
            self._qc_win.set_sites(converted_sites)
            self._pipeline_win.set_input_sites(converted_sites)
            self._forward_win.set_observed_sites(converted_sites)
            self._inversion_win.set_sites(converted_sites)
            self._interp_win.set_sites(converted_sites)
        n = self._controller.n_stations
        self._status_file_lbl.setText(f"{n} stations (converted)")
        self._log(f"Converted dataset committed — {n} stations.")

    # ── Forward → Inversion bridge ────────────────────────────────────

    def _on_forward_send_to_inversion(self, payload: dict) -> None:
        """Show the InversionWindow pre-loaded with the forward model."""
        self._log(
            f"Forward model sent to Inversion (dim={payload.get('dim','1D')})."
        )
        self._inversion_win.load_starting_model(payload)
        self._show_window(self._inversion_win)

    # ── Inversion ─────────────────────────────────────────────────────

    def _open_inversion_wizard(self) -> None:
        self._show_window(self._inversion_win)

    def _on_inversion_result_ready(self, payload: dict) -> None:
        """Forward a completed inversion model to the Interpretation Studio."""
        result = payload.get("result")
        if result is None:
            return
        try:
            from pycsamt.interp._base import ResistivityModel
            if isinstance(result, ResistivityModel) or hasattr(result, "rho_2d"):
                self._interp_win._ctrl.set_model(result)
            else:
                result_dir = getattr(result, "result_dir", None)
                if result_dir:
                    self._interp_win._ctrl.set_model_from_occam2d(result_dir)
        except Exception:
            pass
        self._interp_win._update_status_card()
        self._log("Inversion model forwarded to Interpretation Studio.")

    # ── Export ────────────────────────────────────────────────────────

    def _on_export_figure(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog
        # Try to grab figure from the most recently active panel window
        fig = None
        for win in reversed(self._panel_windows()):
            if not win.isVisible():
                continue
            for attr in ("_canvas", "_result_canvas"):
                canvas = getattr(win, attr, None)
                if canvas is not None:
                    fig = canvas.figure
                    break
            if fig is not None:
                break
        if fig is None:
            self.statusBar().showMessage("No figure to export.", 3000)
            return
        ExportDialog(figure=fig, parent=self).exec()

    # ── Preferences ───────────────────────────────────────────────────

    def _open_preferences(self) -> None:
        from pycsamt.app.desktop.dialogs.preferences_dlg import PreferencesDialog
        dlg = PreferencesDialog(session=self._session, parent=self)
        if dlg.exec() == PreferencesDialog.DialogCode.Accepted:
            self._apply_theme(self._session.theme)
            self._act_theme.setText(
                "☀  Light" if self._session.theme == "dark" else "☾  Dark"
            )
            self._session.save()
            self._log("Preferences saved.")

    # ── API Configuration / Settings ──────────────────────────────────

    def _open_api_config(self, tab: str | None = None) -> None:
        """Open the API Configuration dialog, optionally pre-selecting *tab*."""
        from pycsamt.app.desktop.dialogs.settings_dialog import APIConfigDialog
        dlg = APIConfigDialog(self._settings_ctrl, parent=self, open_tab=tab)
        dlg.settings_changed.connect(self._on_settings_changed)
        dlg.exec()
        # Persist automatically after every dialog session
        try:
            self._settings_ctrl.save()
        except Exception:
            pass

    def _on_settings_changed(self, touched: list) -> None:
        """Refresh open panels that are affected by the changed setting keys."""
        # Pseudosection rendering or view controls → redraw Profile Viewer
        if any(k in touched for k in ("station", "section", "view_controls")):
            try:
                if self._profile_win.isVisible():
                    self._profile_win._on_refresh()
            except Exception:
                pass

        # Pseudosection rendering → redraw QC Dashboard current plot
        if any(k in touched for k in ("station", "section")):
            try:
                if self._qc_win.isVisible():
                    self._qc_win._on_run()
            except Exception:
                pass

        keys_str = ", ".join(touched)
        self._log(f"Settings applied: {keys_str}")

    def _reset_all_settings(self) -> None:
        """Reset all PYCSAMT_* singletons to defaults and refresh panels."""
        from PySide6.QtWidgets import QMessageBox
        reply = QMessageBox.question(
            self,
            "Reset All Settings",
            "Reset ALL API configuration to package defaults?\n\n"
            "This will revert all pseudosection, view, display, and topography "
            "settings to their original values.\n\nThis action cannot be undone.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        self._settings_ctrl.reset_all()
        self._on_settings_changed(["station", "section", "view_controls", "topography"])
        self._log("All API settings reset to package defaults.")

    def _save_settings_profile(self) -> None:
        from PySide6.QtWidgets import QFileDialog
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Settings Profile",
            str(self._settings_ctrl.SETTINGS_PATH.parent),
            "JSON files (*.json)",
        )
        if path:
            try:
                self._settings_ctrl.save(path)
                self._log(f"Settings profile saved → {path}")
            except Exception as exc:
                self._log(f"Save profile error: {exc}")

    def _load_settings_profile(self) -> None:
        from PySide6.QtWidgets import QFileDialog
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Settings Profile",
            str(self._settings_ctrl.SETTINGS_PATH.parent),
            "JSON files (*.json)",
        )
        if path:
            ok = self._settings_ctrl.load(path)
            if ok:
                self._on_settings_changed(
                    ["station", "section", "view_controls", "topography"]
                )
                self._log(f"Settings profile loaded ← {path}")
            else:
                self._log(f"Could not load settings profile: {path}")

    # ── Guard: require loaded survey data ────────────────────────────────

    def _require_sites(self, tool_name: str = "") -> bool:
        """Return True if survey data is loaded; else show NoDataDialog.

        If the user clicks *Load EDI files* in the dialog, the open-file
        action is triggered automatically and this method returns False
        (the caller should abort — the user will retry after loading).
        """
        sites = getattr(self._controller, "sites", None)
        if sites is not None:
            return True
        from pycsamt.app.desktop.dialogs.no_data_dialog import NoDataDialog
        if NoDataDialog.require(self, tool_name):
            self._act_open.trigger()
        return False

    # ── Tools menu handlers ───────────────────────────────────────────

    def _open_strike_analyzer(self) -> None:
        if not self._require_sites("Strike Analyzer"):
            return
        from pycsamt.app.desktop.tools.strike_tool import StrikeAnalyzerDialog
        StrikeAnalyzerDialog(getattr(self._controller, "sites", None), parent=self).exec()

    def _open_edi_validator(self) -> None:
        if not self._require_sites("EDI Validator & Station Manager"):
            return
        from pycsamt.app.desktop.tools.validator_tool import EDIValidatorDialog
        dlg = EDIValidatorDialog(getattr(self._controller, "sites", None), parent=self)
        dlg.open_recompute_requested.connect(self._open_recompute)
        if dlg.exec() and dlg.modified_sites is not None:
            self._apply_modified_sites(dlg.modified_sites, source="EDI Validator")

    def _open_format_converter(self) -> None:
        if not self._require_sites("Format Converter"):
            return
        from pycsamt.app.desktop.tools.converter_tool import FormatConverterDialog
        FormatConverterDialog(getattr(self._controller, "sites", None), parent=self).exec()

    def _open_batch_export(self) -> None:
        from pycsamt.app.desktop.tools.batch_export_tool import BatchExportDialog
        figures = self._collect_figures()
        BatchExportDialog(figures, parent=self).exec()

    def _open_coord_transformer(self) -> None:
        from pycsamt.app.desktop.tools.coord_tool import CoordTransformDialog
        CoordTransformDialog(getattr(self._controller, "sites", None), parent=self).exec()

    def _collect_figures(self) -> list:
        """Return [(label, Figure)] for every visible canvas in all panel windows."""
        figures = []
        _LABEL = {
            "_profile_win":    "Profile",
            "_map_win":        "Map",
            "_qc_win":         "QC",
            "_correction_win": "Correction",
            "_advanced_win":   "Advanced",
            "_tdem_win":       "TDEM",
            "_pipeline_win":   "Pipeline",
            "_forward_win":    "Forward",
            "_inversion_win":  "Inversion",
            "_interp_win":     "Interpretation",
        }
        for attr, label in _LABEL.items():
            win = getattr(self, attr, None)
            if win is None or not win.isVisible():
                continue
            for canvas_attr in ("_canvas", "_result_canvas"):
                canvas = getattr(win, canvas_attr, None)
                if canvas is not None and hasattr(canvas, "figure"):
                    figures.append((label, canvas.figure))
        return figures

    def _open_station_response(self) -> None:
        if not self._require_sites("Station Response Inspector"):
            return
        from pycsamt.app.desktop.tools.station_response_tool import StationResponseDialog
        StationResponseDialog(getattr(self._controller, "sites", None), parent=self).exec()

    def _open_strike_profile(self) -> None:
        if not self._require_sites("Strike Profile Viewer"):
            return
        from pycsamt.app.desktop.tools.strike_profile_tool import StrikeProfileDialog
        StrikeProfileDialog(getattr(self._controller, "sites", None), parent=self).exec()

    def _open_phase_tensor_map(self) -> None:
        if not self._require_sites("Phase Tensor Map"):
            return
        from pycsamt.app.desktop.tools.phase_tensor_map_tool import PhaseTensorMapDialog
        PhaseTensorMapDialog(getattr(self._controller, "sites", None), parent=self).exec()

    def _open_dimensionality(self) -> None:
        if not self._require_sites("Dimensionality Classifier"):
            return
        from pycsamt.app.desktop.tools.dimensionality_tool import DimensionalityDialog
        DimensionalityDialog(getattr(self._controller, "sites", None), parent=self).exec()

    def _open_frequency_editor(self) -> None:
        if not self._require_sites("Frequency Editor"):
            return
        from pycsamt.app.desktop.tools.frequency_editor_tool import FrequencyEditorDialog
        dlg = FrequencyEditorDialog(getattr(self._controller, "sites", None), parent=self)
        if dlg.exec() and dlg.edited_sites is not None:
            self._apply_modified_sites(dlg.edited_sites, source="Frequency Editor")

    def _open_layered_model(self) -> None:
        from pycsamt.app.desktop.tools.layered_model_tool import LayeredModelDialog
        LayeredModelDialog(parent=self).exec()

    def _open_elevation_enrichment(self) -> None:
        if not self._require_sites("Elevation Enrichment"):
            return
        from pycsamt.app.desktop.tools.elevation_tool import ElevationEnrichDialog
        ElevationEnrichDialog(getattr(self._controller, "sites", None), parent=self).exec()

    def _open_recompute(self) -> None:
        from PySide6.QtWidgets import QMessageBox
        from pycsamt.app.desktop.dialogs.recompute_dlg import RecomputeDialog

        sites = getattr(self._controller, "sites", None)

        # If already recomputed — warn and ask for confirmation before re-opening
        if self._recomputed_ids:
            n   = len(self._recomputed_ids)
            ans = QMessageBox.question(
                self,
                "Already Recomputed",
                f"{n} station(s) in this survey have already been recomputed.\n\n"
                "Do you want to recompute again with new settings?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if ans != QMessageBox.StandardButton.Yes:
                return

        dlg = RecomputeDialog(sites=sites, parent=self)
        dlg.recompute_committed.connect(self._on_recompute_committed)
        dlg.exec()

    def _on_recompute_committed(self, result) -> None:
        """Apply a completed recompute result back into the app."""
        new_ids = {rec.station for rec in result.records if rec.status == "ok"}
        self._recomputed_ids.update(new_ids)
        self._last_recompute_output = result.output_root

        # Load first: set_dataframe() clears model's _recomputed_ids internally.
        # Badges must be applied AFTER the new DataFrame is in place so the
        # station IDs from the fresh data are the ones being matched.
        recomputed_sites = result.sites
        self._apply_modified_sites(recomputed_sites, source="Recompute")

        # Now stamp the badges onto the freshly-loaded model rows.
        self._station_panel._table._model.mark_recomputed(self._recomputed_ids)

        n = len(new_ids)
        self._log(f"Recompute committed — {n} station(s) marked with ◈ badge.")
        self.statusBar().showMessage(
            f"Recomputed {n} station(s). Marked with ◈ in the station list.", 6000
        )

    def _apply_modified_sites(self, new_sites, *, source: str = "") -> None:
        """Push a modified site list back into the controller and refresh all panels."""
        self._controller.set_sites(new_sites)
        try:
            from pycsamt.app.desktop.controllers.data_controller import DataController
            dc = DataController()
            dc._sites = new_sites
            dc._df    = dc._build_dataframe()
            df = dc.dataframe
            if df is not None and not df.empty:
                self._station_panel.set_dataframe(df)
                self._survey_overview.update_survey(
                    df, getattr(self, "_loaded_paths", None)
                )
                try:
                    self._map_win.set_dataframe(df)
                    self._map_win.set_sites(new_sites)
                except Exception:
                    pass
                try:
                    self._profile_win.set_sites(new_sites)
                except Exception:
                    pass
        except Exception:
            pass
        n = len(new_sites) if hasattr(new_sites, "__len__") else "?"
        tag = f"{source}: " if source else ""
        self._log(f"{tag}{n} station(s) applied to survey.")

    # ── Help menu handlers ────────────────────────────────────────────

    def _open_documentation(self) -> None:
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices
        QDesktopServices.openUrl(QUrl("http://pycsamt.readthedocs.io/"))

    def _open_github(self) -> None:
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices
        QDesktopServices.openUrl(QUrl("https://github.com/earthai-tech/pycsamt"))

    def _open_about(self) -> None:
        from pycsamt.app.desktop.dialogs.about_dialog import AboutDialog
        AboutDialog(parent=self).exec()

    # ── Recent files ──────────────────────────────────────────────────

    def _rebuild_recent_menu(self) -> None:
        self._act_recent.clear()
        if not self._session.recent_files:
            a = QAction("(none)", self)
            a.setEnabled(False)
            self._act_recent.addAction(a)
            return
        for path in self._session.recent_files[:20]:
            a = QAction(path, self)
            a.triggered.connect(
                lambda _=False, p=path: self._start_loading([p])
            )
            self._act_recent.addAction(a)

    # ── Search / filter ───────────────────────────────────────────────

    def _on_filter_changed(self, text: str) -> None:
        try:
            self._station_panel.filter(text)
        except Exception:
            pass

    # ── Session ───────────────────────────────────────────────────────

    def _on_save_session(self) -> None:
        self._save_layout()
        self._session.save()
        self._log("Session saved.")
        self.statusBar().showMessage("Session saved.", 3000)

    def _save_layout(self) -> None:
        self._session.dock_geometry = (
            self.saveGeometry().toBase64().data().decode()
        )
        self._session.dock_state = (
            self.saveState().toBase64().data().decode()
        )
        for win in self._panel_windows():
            win.save_geometry_to(self._session.window_geometries)

    def _restore_layout(self) -> None:
        if self._session.dock_geometry:
            try:
                self.restoreGeometry(
                    QByteArray.fromBase64(
                        self._session.dock_geometry.encode()
                    )
                )
            except Exception:
                pass
        if self._session.dock_state:
            try:
                self.restoreState(
                    QByteArray.fromBase64(
                        self._session.dock_state.encode()
                    )
                )
            except Exception:
                pass

    # ── Error handling ────────────────────────────────────────────────

    def _on_load_error(self, message: str) -> None:
        self._progress_bar.setVisible(False)
        self._status_ready_lbl.setText("Error  ✕")
        self._log(f"ERROR: {message}")
        self.statusBar().showMessage(f"Load failed: {message}", 6000)

    # ── Logging helper ────────────────────────────────────────────────

    def _log(self, text: str) -> None:
        self._log_panel.append_line(text)

    def log(self, text: str) -> None:
        """Public logging entry point — delegates to the log panel."""
        self._log_panel.append_line(text)

    # ── Public helpers ────────────────────────────────────────────────

    def set_status(self, message: str, timeout_ms: int = 4000) -> None:
        self.statusBar().showMessage(message, timeout_ms)

    def set_file_label(self, text: str) -> None:
        """Update the status-bar file/dataset label."""
        self._status_file_lbl.setText(text)

    def set_freq_label(self, text: str) -> None:
        """Update the status-bar frequency info label."""
        self._status_freq_lbl.setText(text)

    # ── Qt overrides ──────────────────────────────────────────────────

    def closeEvent(self, event) -> None:
        from PySide6.QtWidgets import QMessageBox
        reply = QMessageBox.question(
            self,
            "Quit pycsamt",
            "Are you sure you want to quit?\n\n"
            "Any unsaved session data will be lost.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel,
        )
        if reply != QMessageBox.StandardButton.Yes:
            event.ignore()
            return
        self._save_layout()
        self._session.save()
        for win in self._panel_windows():
            win.close()
        self._log("Goodbye.")
        super().closeEvent(event)
