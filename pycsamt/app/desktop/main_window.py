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
    │  Toolbar:  Open  Save  │  Profile  Map  QC ...           │
    │            Inversion       │  Agents   │  Export  Theme  │
    ├──────────────────┬───────────────────────────────────────┤
    │  Station List    │  Station Detail Card                  │
    │  (search + table)│  (name, coords, quality, actions)    │
    ├──────────────────┴───────────────────────────────────────┤
    │  Log  (collapsible, 80 px)                               │
    └──────────────────────────────────────────────────────────┘

Menu bar
────────
    File     — Open / Load Data, Save Session, Save Edited Survey,
               Recent Files, Quit
    View     — toggle each panel window + theme
    Settings — API Configuration, Reset to Defaults, Save/Load Profile
    Help     — Documentation, About
"""

from __future__ import annotations

import re
from pathlib import Path

from PySide6.QtCore import QByteArray, QSize, Qt
from PySide6.QtGui import (
    QAction,
    QActionGroup,
    QIcon,
    QKeySequence,
    QPainter,
    QPixmap,
)
from PySide6.QtWidgets import (
    QApplication,
    QDockWidget,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMenu,
    QMenuBar,
    QProgressBar,
    QPushButton,
    QSplitter,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop import branding
from pycsamt.app.desktop.agent_master_bridge import (
    launch_agent_master,
)
from pycsamt.app.desktop.controllers.app_controller import (
    AppController,
)
from pycsamt.app.desktop.models.session import SessionState
from pycsamt.app.desktop.panels.log_panel import LogPanel
from pycsamt.app.desktop.panels.station_panel import (
    StationPanel,
)
from pycsamt.app.desktop.widgets.mpl_canvas import (
    apply_mpl_dark_theme,
    apply_mpl_light_theme,
)
from pycsamt.app.desktop.widgets.station_detail import (
    StationDetailCard,
)
from pycsamt.app.desktop.widgets.survey_overview import (
    SurveyOverviewWidget,
)
from pycsamt.app.desktop.windows import (
    AdvancedToolsWindow,
    AirborneWindow,
    CorrectionWindow,
    ForwardModelWindow,
    InterpretationWindow,
    InversionWindow,
    MapViewerWindow,
    Pcsf3DWindow,
    ProfileViewerWindow,
    QCDashboardWindow,
    TDEMWindow,
)
from pycsamt.app.desktop.windows.pipeline_window import (
    PipelineWindow,
)
from pycsamt.app.desktop.windows.solver_builder_window import (
    SolverBuilderWindow,
)

_RESOURCES = branding.RESOURCES_DIR
_ICONS = branding.ICONS_DIR

# Tracks current theme so _icon() can recolor SVGs without being passed a flag.
_DARK_MODE: bool = False

# Matches any 6-digit hex colour in an SVG source string.
_HEX6_RE = re.compile(r"#[0-9a-fA-F]{6}", re.IGNORECASE)


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
    result = re.sub(
        r'(?<=[";\s:])black(?=[";\s])', target, result, flags=re.IGNORECASE
    )

    # Pass 3: inject fill on individual shape elements that carry no explicit fill,
    # so icons that use the SVG default (black) fill — e.g. interpret.svg — become
    # visible on a dark background.  Elements with an existing fill= or fill: are
    # left untouched, and CSS class rules (fill:none) still override the injected
    # presentation attribute via normal specificity.
    def _maybe_fill(m: re.Match) -> str:
        tag = m.group(0)
        if "fill=" not in tag and "fill:" not in tag:
            return (
                tag[:-2] + f' fill="{target}"/>'
                if tag.endswith("/>")
                else tag[:-1] + f' fill="{target}">'
            )
        return tag

    result = re.sub(
        r"<(?:path|rect|circle|ellipse|polygon|polyline|line)\b[^>]*?(?:/>|>)",
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
                    renderer = QSvgRenderer(QByteArray(svg_bytes))
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
        self._session = SessionState.load()
        self._controller = AppController(self._session)
        self._loader = None  # LoaderWorker kept alive while running
        self._load_in_progress = False
        self._append_on_commit = False
        self._pending_load_paths = []
        self._recomputed_ids: set[str] = set()  # stations marked recomputed
        self._last_recompute_output = (
            None  # Path to last recomputed EDI folder
        )
        # The unfiltered survey remains available while a smaller, global
        # active-line scope is distributed to every tool window.
        self._all_sites = None
        from pycsamt.app.desktop.controllers.edit_history import EditHistory

        self._history = EditHistory()  # Edit ▸ Undo / Redo / History
        self._all_dataframe = None
        self._line_actions: dict[str, QAction] = {}
        self._converter_app_window = None

        # Settings controller — load saved profile (silently no-ops if absent)
        from pycsamt.app.desktop.controllers.settings_controller import (
            SettingsController,
        )

        self._settings_ctrl = SettingsController()
        self._settings_ctrl.load()

        self._setup_window()
        self._create_panel_windows()  # windows must exist before theme is applied
        self._apply_theme(self._session.theme)
        self._create_central_widget()
        self._create_log_dock()
        self._create_menu_bar()
        self._create_tool_bar()
        self._create_status_bar()
        self._restore_layout()
        self._wire_signals()
        self._log("pycsamt ready — load survey data (EDI, EMTF-XML, AVG, "
                  "J) to begin.")

    # ── Window setup ──────────────────────────────────────────────────

    def _setup_window(self) -> None:
        self.setWindowTitle("pycsamt")
        self.setMinimumSize(760, 480)
        self._fit_to_screen(force_default=True)
        ico = _ICONS / "pycsamt.logo.ico"
        if ico.exists():
            self.setWindowIcon(QIcon(str(ico)))
        # Tracks every icon-bearing QAction so _apply_theme can re-ice them all.
        # Populated by _create_menu_bar and _create_tool_bar (both run after _setup_window).
        self._all_icon_actions: list = []

    def _fit_to_screen(self, *, force_default: bool = False) -> None:
        """Keep the main window comfortably inside the active desktop."""
        screen = self.screen() or QApplication.primaryScreen()
        if screen is None:
            if force_default:
                self.resize(1180, 700)
            return

        available = screen.availableGeometry()
        max_width = max(
            self.minimumWidth(),
            min(1280, int(available.width() * 0.75), available.width() - 40),
        )
        max_height = max(
            self.minimumHeight(),
            min(760, int(available.height() * 0.75), available.height() - 40),
        )

        if force_default:
            width, height = min(1180, max_width), min(700, max_height)
        else:
            width = min(max(self.width(), self.minimumWidth()), max_width)
            height = min(max(self.height(), self.minimumHeight()), max_height)

        self.setWindowState(Qt.WindowState.WindowNoState)
        x = available.x() + (available.width() - width) // 2
        y = available.y() + (available.height() - height) // 2
        self.setGeometry(x, y, width, height)

    # ── Theme ──────────────────────────────────────────────────────────

    def _apply_theme(self, theme: str) -> None:
        global _DARK_MODE
        _DARK_MODE = theme == "dark"

        qss = _RESOURCES / f"{theme}_theme.qss"
        if qss.exists():
            QApplication.instance().setStyleSheet(
                qss.read_text(encoding="utf-8")
            )
        (apply_mpl_dark_theme if theme == "dark" else apply_mpl_light_theme)()
        self._session.theme = theme
        if hasattr(self, "_log_panel"):
            self._log_panel.set_dark(theme == "dark")
        self._style_toolbar_overflow()

        # Sync the View > Theme checkboxes (created after first call)
        if hasattr(self, "_act_dark"):
            self._act_dark.setChecked(theme == "dark")
            self._act_light.setChecked(theme == "light")

        # Re-apply ALL icon-bearing actions (menu + toolbar)
        for action, icon_name in getattr(self, "_all_icon_actions", []):
            action.setIcon(_icon(icon_name))

        for win in self._panel_windows():
            try:
                win.set_dark_mode(theme == "dark")
            except Exception:
                pass
        converter = getattr(self, "_converter_app_window", None)
        if converter is not None:
            converter.set_host_theme(theme)

    def _toggle_theme(self) -> None:
        new = "light" if self._session.theme == "dark" else "dark"
        self._apply_theme(new)
        # Show the next target: now in *new*, so label points to the other theme.
        self._act_theme.setText("☾  Dark" if new == "light" else "☀  Light")

    # ── Independent panel windows (created once, shown/hidden on demand) ──

    def _create_panel_windows(self) -> None:
        self._profile_win = ProfileViewerWindow(parent=self)
        self._map_win = MapViewerWindow(parent=self)
        self._pcsf3d_win = Pcsf3DWindow(parent=self)
        self._qc_win = QCDashboardWindow(parent=self)
        self._correction_win = CorrectionWindow(parent=self)
        self._advanced_win = AdvancedToolsWindow(parent=self)
        self._tdem_win = TDEMWindow(parent=self)
        self._pipeline_win = PipelineWindow(parent=self)
        # "Apply to main data" in Pipeline Studio (was never connected, so a
        # pipeline result used to stay trapped in the pipeline window)
        self._pipeline_win.pipeline_finished.connect(self._on_pipeline_applied)
        self._forward_win = ForwardModelWindow(parent=self)
        self._inversion_win = InversionWindow(parent=self)
        self._interp_win = InterpretationWindow(parent=self)
        self._airborne_win = AirborneWindow(parent=self)
        # Tools ▸ Solver Builder: compiled binaries flow into Inversion
        self._solver_builder_win = SolverBuilderWindow(parent=self)
        self._solver_builder_win.binary_built.connect(
            lambda _key, _binary: self._inversion_win.refresh_solver_binaries()
        )
        self._inversion_win.build_solver_requested.connect(
            self._open_solver_builder
        )

        # Wire forward → inversion bridge
        self._forward_win.send_to_inversion.connect(
            self._on_forward_send_to_inversion
        )
        # Wire inversion result → interpretation model
        self._inversion_win.result_ready.connect(
            self._on_inversion_result_ready
        )
        self._inversion_win.run_state.connect(self._on_inversion_run_state)

        # Wire correction window → main data
        self._correction_win.corrections_committed.connect(
            self._on_corrections_committed
        )
        # Wire advanced tools conversion → main data
        self._advanced_win.conversion_committed.connect(
            self._on_conversion_committed
        )

        # Restore positions from session
        geo = self._session.window_geometries
        for win in self._panel_windows():
            restore_geometry = getattr(win, "restore_geometry_from", None)
            if callable(restore_geometry):
                restore_geometry(geo)

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

        line_row = QHBoxLayout()
        line_row.setSpacing(4)
        self._line_scope_btn = QToolButton(left)
        self._line_scope_btn.setText("Lines: all")
        self._line_scope_btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self._line_scope_btn.setToolTip(
            "Choose survey lines included throughout the application"
        )
        self._line_scope_menu = QMenu(self._line_scope_btn)
        self._line_scope_btn.setMenu(self._line_scope_menu)
        self._primary_line_combo = QComboBox(left)
        self._primary_line_combo.setToolTip(
            "Primary line used by operations that require one profile"
        )
        self._primary_line_combo.currentTextChanged.connect(
            self._on_primary_line_changed
        )
        line_row.addWidget(self._line_scope_btn, 1)
        primary_label = QLabel("Primary:", left)
        primary_label.setObjectName("LineScopeLabel")
        line_row.addWidget(primary_label)
        line_row.addWidget(self._primary_line_combo, 1)
        left_v.addLayout(line_row)

        # Station table
        self._station_panel = StationPanel(left)
        left_v.addWidget(self._station_panel)

        # ── Right: vertical splitter — survey overview + detail card ──
        self._stats_tabs = QTabWidget(container)
        self._stats_tabs.setDocumentMode(True)

        self._survey_overview = SurveyOverviewWidget(container)
        self._survey_overview.setObjectName("SurveyOverview")

        self._detail_card = StationDetailCard(container)
        self._detail_card.setObjectName("StationDetailCard")
        self._detail_card.set_animations(
            getattr(self._session, "ui_animations", True))

        self._stats_tabs.addTab(self._survey_overview, "Survey Overview")
        self._stats_tabs.addTab(self._detail_card, "Station Statistics")

        splitter.addWidget(left)
        splitter.addWidget(self._stats_tabs)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([280, 720])

        h_layout.addWidget(splitter)
        self.setCentralWidget(container)

    # ── Log dock (bottom, thin) ───────────────────────────────────────

    def _create_log_dock(self) -> None:
        """The Log dock: hidden by default (status-bar "Log" chip, View ▸
        Log); keeps an unread count while hidden."""
        self._log_panel = LogPanel(self)
        self._log_panel.set_dark(self._session.theme == "dark")
        self._log_unread = {"n": 0, "level": "info"}
        self._log_dock = QDockWidget("Log", self)
        self._log_dock.setObjectName("LogDock")
        self._log_dock.setAllowedAreas(Qt.DockWidgetArea.BottomDockWidgetArea)
        self._log_dock.setWidget(self._log_panel)
        self._log_dock.setMinimumHeight(90)
        self.addDockWidget(
            Qt.DockWidgetArea.BottomDockWidgetArea, self._log_dock
        )
        self.resizeDocks([self._log_dock], [160], Qt.Orientation.Vertical)
        self._log_dock.setVisible(bool(self._session.log_visible))
        self._log_panel.message_logged.connect(self._on_log_message)
        self._log_dock.visibilityChanged.connect(self._on_log_visibility)

    def _on_log_message(self, level: str, text: str) -> None:
        if level == "error":
            self.statusBar().showMessage(text, 6000)
        if not self._log_dock.isHidden():  # open (even if window minimised)
            return
        u = self._log_unread
        u["n"] += 1
        rank = {"info": 0, "ok": 0, "warn": 1, "error": 2}
        if rank.get(level, 0) > rank.get(u["level"], 0):
            u["level"] = level
        self._update_log_chip()

    def set_log_visible(self, on: bool) -> None:
        """Show / hide the log and remember it (user choice)."""
        self._session.log_visible = bool(on)
        self._log_dock.setVisible(bool(on))
        self._on_log_visibility(bool(on))

    def _on_log_visibility(self, visible: bool) -> None:
        # also fires when the whole window hides/closes: the remembered
        # choice changes only through set_log_visible / View ▸ Log
        if visible:
            self._log_unread.update(n=0, level="info")
        self._update_log_chip()

    def _update_log_chip(self) -> None:
        chip = getattr(self, "_log_chip", None)
        if chip is None:
            return
        u = self._log_unread
        shown = not self._log_dock.isHidden()
        chip.setText("Log" + (f"  {u['n']}" if u["n"] and not shown else ""))
        chip.setChecked(shown)
        level = u["level"] if u["n"] else "info"
        chip.setProperty("level", level)
        chip.setToolTip(
            ("Hide" if shown else "Show") + " the log  (Ctrl+Alt+L)"
            + (f" — {u['n']} new message(s)" if u["n"] and not shown else ""))
        chip.style().unpolish(chip)
        chip.style().polish(chip)

    # ── Menu bar (File | View | Help) ─────────────────────────────────

    def _create_menu_bar(self) -> None:
        mb = self.menuBar()

        # ── File ──────────────────────────────────────────────────────
        file_menu = mb.addMenu("&File")

        self._act_open = QAction(_icon("open"), "&Open / Load Data…", self)
        self._act_open.setShortcut(QKeySequence.StandardKey.Open)
        self._act_open.setStatusTip(
            "Open survey data: EDI, EMTF-XML transfer functions, AVG or J "
            "files")
        self._act_open.triggered.connect(self._open_files)

        self._act_save = QAction(_icon("save-session"), "&Save Session", self)
        self._act_save.setShortcut(QKeySequence.StandardKey.Save)
        self._act_save.triggered.connect(self._on_save_session)

        self._act_recent = file_menu.addMenu("Recent Files")
        self._rebuild_recent_menu()

        act_quit = QAction(_icon("quit"), "&Quit", self)
        act_quit.setShortcut(QKeySequence.StandardKey.Quit)
        act_quit.triggered.connect(self.close)

        self._act_save_survey = QAction(_icon("save"), "Save &Edited Survey…",
                                        self)
        self._act_save_survey.setShortcut("Ctrl+Alt+S")
        self._act_save_survey.setStatusTip(
            "Write the edited survey as EDI or EMTF-XML files (the loaded "
            "files are never overwritten by an edit)")
        self._act_save_survey.triggered.connect(self._save_edited_survey)

        file_menu.addAction(self._act_open)
        file_menu.addAction(self._act_save)
        file_menu.addAction(self._act_save_survey)
        file_menu.addSeparator()
        file_menu.addMenu(self._act_recent)
        file_menu.addSeparator()
        file_menu.addAction(act_quit)

        act_prefs = QAction(_icon("tools"), "&Preferences…", self)
        act_prefs.setShortcut(QKeySequence.StandardKey.Preferences)
        act_prefs.triggered.connect(self._open_preferences)
        self._build_edit_menu(mb, act_prefs)

        # ── View ──────────────────────────────────────────────────────
        view_menu = mb.addMenu("&View")

        act_profile = QAction(_icon("profile-view"), "&Profile Viewer", self)
        act_profile.setShortcut("Ctrl+P")
        act_profile.triggered.connect(
            lambda: self._show_window(self._profile_win)
        )
        view_menu.addAction(act_profile)

        self._act_profile_on_select = QAction(
            "Open Profile when a station is selected", self, checkable=True
        )
        self._act_profile_on_select.setStatusTip(
            "Automatically open Profile Viewer when selecting a station"
        )
        self._act_profile_on_select.setChecked(
            self._session.open_profile_on_station_select
        )
        self._act_profile_on_select.toggled.connect(
            self._on_profile_on_select_toggled
        )
        view_menu.addAction(self._act_profile_on_select)
        view_menu.addSeparator()

        act_map = QAction(_icon("map-view"), "&Map Viewer", self)
        act_map.setShortcut("Ctrl+M")
        act_map.triggered.connect(lambda: self._show_window(self._map_win))
        view_menu.addAction(act_map)

        act_pcsf3d = QAction(_icon("3d"), "&PCSF 3D Viewer", self)
        act_pcsf3d.setShortcut("Ctrl+Shift+M")
        act_pcsf3d.triggered.connect(
            lambda: self._show_window(self._pcsf3d_win)
        )
        view_menu.addAction(act_pcsf3d)

        act_qc = QAction(_icon("qc"), "&QC Dashboard", self)
        act_qc.setShortcut("Ctrl+Shift+Q")  # Ctrl+Q is Quit
        act_qc.triggered.connect(lambda: self._show_window(self._qc_win))
        view_menu.addAction(act_qc)

        act_corr = QAction(
            _icon("sites-correction"), "&Data Corrections", self
        )
        act_corr.setShortcut("Ctrl+R")
        act_corr.triggered.connect(
            lambda: self._show_window(self._correction_win)
        )
        view_menu.addAction(act_corr)

        act_fwd = QAction(_icon("forward"), "&Forward Modelling", self)
        act_fwd.setShortcut("Ctrl+Alt+F")  # Ctrl+F is Find station
        act_fwd.triggered.connect(lambda: self._show_window(self._forward_win))
        view_menu.addAction(act_fwd)

        act_inv = QAction(_icon("inversion"), "&Inversion Wizard…", self)
        act_inv.setShortcut("Ctrl+I")
        act_inv.triggered.connect(self._open_inversion_wizard)
        view_menu.addAction(act_inv)

        act_interp = QAction(
            _icon("interpret"), "&Interpretation Studio", self
        )
        act_interp.setShortcut("Ctrl+Shift+I")
        act_interp.triggered.connect(
            lambda: self._show_window(self._interp_win)
        )
        view_menu.addAction(act_interp)

        act_airborne = QAction(
            _icon("induction"), "&Airborne EM…", self
        )
        act_airborne.setShortcut("Ctrl+Shift+B")
        act_airborne.triggered.connect(
            lambda: self._show_window(self._airborne_win)
        )
        view_menu.addAction(act_airborne)

        act_pipe = QAction(_icon("pipeline"), "&Processing Pipeline", self)
        act_pipe.setShortcut("Ctrl+Shift+P")
        act_pipe.triggered.connect(
            lambda: self._show_window(self._pipeline_win)
        )
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
        act_adv.triggered.connect(
            lambda: self._show_window(self._advanced_win)
        )
        view_menu.addAction(act_adv)

        view_menu.addSeparator()
        _log_tva = self._log_dock.toggleViewAction()
        _log_tva.setIcon(_icon("log"))
        _log_tva.setShortcut("Ctrl+Alt+L")
        _log_tva.triggered.connect(
            lambda on: setattr(self._session, "log_visible", bool(on)))
        view_menu.addAction(_log_tva)
        view_menu.addSeparator()

        theme_menu = view_menu.addMenu("Theme")
        self._act_dark = QAction("☾  Dark", self, checkable=True)
        self._act_light = QAction("☀  Light", self, checkable=True)
        # QActionGroup enforces mutual exclusivity: checking one unchecks the other.
        _theme_grp = QActionGroup(self)
        _theme_grp.setExclusive(True)
        _theme_grp.addAction(self._act_dark)
        _theme_grp.addAction(self._act_light)
        self._act_dark.setChecked(self._session.theme == "dark")
        self._act_light.setChecked(self._session.theme == "light")
        self._act_dark.triggered.connect(lambda: self._apply_theme("dark"))
        self._act_light.triggered.connect(lambda: self._apply_theme("light"))
        theme_menu.addAction(self._act_dark)
        theme_menu.addAction(self._act_light)

        # ── Tools ─────────────────────────────────────────────────────
        tools_menu = mb.addMenu("&Tools")

        act_strike = QAction(
            _icon("strike-analyzer"), "&Strike Analyzer…", self
        )
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


        act_conv = QAction(
            _icon("format-converter"), "&Format Converter…", self
        )
        act_conv.setShortcut("Ctrl+Shift+C")
        act_conv.setStatusTip("Export survey to EDI / CSV / JSON")
        act_conv.triggered.connect(self._open_format_converter)
        tools_menu.addAction(act_conv)

        act_batch = QAction(
            _icon("batch-export"), "&Batch Export Plots…", self
        )
        act_batch.setShortcut("Ctrl+Shift+E")
        act_batch.setStatusTip("Save every open canvas figure to a folder")
        act_batch.triggered.connect(self._open_batch_export)
        tools_menu.addAction(act_batch)

        tools_menu.addSeparator()

        act_coord = QAction(
            _icon("coordinate-transformer"), "&Coordinate Transformer…", self
        )
        act_coord.setShortcut("Ctrl+Shift+G")
        act_coord.setStatusTip("Convert UTM ↔ Lat/Lon (pyproj-backed)")
        act_coord.triggered.connect(self._open_coord_transformer)
        tools_menu.addAction(act_coord)

        tools_menu.addSeparator()

        act_station_resp = QAction(
            _icon("station-response"), "Station &Response Inspector…", self
        )
        act_station_resp.setShortcut("Ctrl+Shift+R")
        act_station_resp.setStatusTip(
            "Plot full impedance tensor response for a selected station"
        )
        act_station_resp.triggered.connect(self._open_station_response)
        tools_menu.addAction(act_station_resp)

        act_strike_profile = QAction(
            _icon("strike-profile"), "Strike &Profile Viewer…", self
        )
        act_strike_profile.setShortcut("Ctrl+Alt+P")  # ≠ Pipeline
        act_strike_profile.setStatusTip(
            "Strike angle vs. station-position line plot with IQR ribbon"
        )
        act_strike_profile.triggered.connect(self._open_strike_profile)
        tools_menu.addAction(act_strike_profile)

        act_pt_map = QAction(_icon("phase-tensor"), "Phase &Tensor Map…", self)
        act_pt_map.setShortcut("Ctrl+Alt+T")  # ≠ Advanced Tools
        act_pt_map.setStatusTip(
            "Geographic map of phase-tensor ellipses at a chosen period"
        )
        act_pt_map.triggered.connect(self._open_phase_tensor_map)
        tools_menu.addAction(act_pt_map)

        act_pt_strip_grid = QAction(
            _icon("phase-tensor"), "Phase Tensor Strip &Grid…", self
        )
        act_pt_strip_grid.setShortcut("Ctrl+Shift+N")
        act_pt_strip_grid.setStatusTip(
            "Ellipse-strip vs. period, tiled by survey line (multi-profile view)"
        )
        act_pt_strip_grid.triggered.connect(self._open_phase_tensor_strip_grid)
        tools_menu.addAction(act_pt_strip_grid)

        tools_menu.addSeparator()

        act_dim = QAction(
            _icon("dimensionnality"), "&Dimensionality Classifier…", self
        )
        act_dim.setShortcut("Ctrl+Shift+D")
        act_dim.setStatusTip(
            "Classify each station × frequency as 1D / 2D / 3D"
        )
        act_dim.triggered.connect(self._open_dimensionality)
        tools_menu.addAction(act_dim)

        act_freq_ed = QAction(
            _icon("frequency-editor"), "&Frequency Editor…", self
        )
        act_freq_ed.setShortcut("Ctrl+Shift+F")
        act_freq_ed.setStatusTip(
            "Confidence-based frequency QC: drop / mask / recover frequency bands"
        )
        act_freq_ed.triggered.connect(self._open_frequency_editor)

        tools_menu.addSeparator()

        act_lm = QAction(
            _icon("layered-model"), "&Layered Model Builder…", self
        )
        act_lm.setShortcut("Ctrl+Shift+L")
        act_lm.setStatusTip("Build and preview a 1-D layered earth model")
        act_lm.triggered.connect(self._open_layered_model)
        tools_menu.addAction(act_lm)

        act_elev = QAction(_icon("elevation"), "&Elevation Enrichment…", self)
        act_elev.setShortcut("Ctrl+Shift+H")
        act_elev.setStatusTip(
            "Fetch elevation for all loaded stations via open elevation API"
        )
        act_elev.triggered.connect(self._open_elevation_enrichment)
        # Frequency Editor, Recompute EDIs and Elevation Enrichment live
        # in the Edit menu (they change the survey)
        self._edit_freq_menu.addAction(act_freq_ed)
        self._edit_tensor_menu.addAction(act_recompute)
        self._edit_station_menu.addAction(act_elev)

        # Last: solver compilation belongs to inversion, which follows the
        # processing tools above.
        tools_menu.addSeparator()
        act_solvers = QAction(_icon("tools"), "&Solver Builder…", self)
        act_solvers.setStatusTip(
            "Compile Occam2D, ModEM 2-D/3-D or MARE2DEM (installs the "
            "compilers if needed)"
        )
        act_solvers.triggered.connect(lambda: self._open_solver_builder(None))
        tools_menu.addAction(act_solvers)

        self._all_icon_actions.extend(
            [
                (act_strike, "strike-analyzer"),
                (act_valid, "EDI-validator"),
                (act_recompute, "recompute"),
                (act_conv, "format-converter"),
                (act_batch, "batch-export"),
                (act_coord, "coordinate-transformer"),
                (act_station_resp, "station-response"),
                (act_strike_profile, "strike-profile"),
                (act_pt_map, "phase-tensor"),
                (act_pt_strip_grid, "phase-tensor"),
                (act_dim, "dimensionnality"),
                (act_freq_ed, "frequency-editor"),
                (act_lm, "layered-model"),
                (act_elev, "elevation"),
            ]
        )

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
        act_reset_all.setStatusTip(
            "Reset all API singletons to package defaults"
        )
        act_reset_all.triggered.connect(self._reset_all_settings)
        settings_menu.addAction(act_reset_all)

        settings_menu.addSeparator()

        act_save_profile = QAction("Save Profile…", self)
        act_save_profile.setStatusTip(
            "Save current API configuration to a JSON profile"
        )
        act_save_profile.triggered.connect(self._save_settings_profile)
        settings_menu.addAction(act_save_profile)

        act_load_profile = QAction("Load Profile…", self)
        act_load_profile.setStatusTip(
            "Load an API configuration profile from JSON"
        )
        act_load_profile.triggered.connect(self._load_settings_profile)
        settings_menu.addAction(act_load_profile)

        self._all_icon_actions.append((act_api, "tools"))

        # ── Help ──────────────────────────────────────────────────────
        help_menu = QMenu("&Help", self)
        self._help_menu = help_menu
        help_bar = QMenuBar(mb)
        help_bar.addMenu(help_menu)
        mb.setCornerWidget(help_bar, Qt.Corner.TopRightCorner)
        act_docs = QAction(_icon("docs"), "&Documentation", self)
        act_docs.setStatusTip("Open pycsamt documentation in your browser")
        act_docs.triggered.connect(self._open_documentation)
        help_menu.addAction(act_docs)

        act_gh = QAction(_icon("github"), "pycsamt on &GitHub", self)
        act_gh.setStatusTip(
            "Open the pycsamt GitHub repository in your browser"
        )
        act_gh.triggered.connect(self._open_github)
        help_menu.addAction(act_gh)

        help_menu.addSeparator()

        act_about = QAction(_icon("help"), "&About pycsamt", self)
        act_about.setStatusTip("About pycsamt v2 — version, author and links")
        act_about.triggered.connect(self._open_about)
        help_menu.addAction(act_about)

        # Register every icon-bearing action created in this method so that
        # _apply_theme can re-ice all of them when the user switches theme.
        self._all_icon_actions.extend(
            [
                (self._act_open, "open"),
                (self._act_save, "save-session"),
                (act_prefs, "tools"),
                (act_profile, "profile-view"),
                (act_map, "map-view"),
                (act_qc, "qc"),
                (act_corr, "sites-correction"),
                (act_fwd, "forward"),
                (act_inv, "inversion"),
                (act_interp, "interpret"),
                (act_pipe, "pipeline"),
                (act_tdem, "tdem"),
                (act_agents, "agents"),
                (act_adv, "advanced-tools"),
                (_log_tva, "log"),
                (act_docs, "docs"),
                (act_gh, "github"),
                (act_about, "help"),
            ]
        )

    # ── Toolbar ───────────────────────────────────────────────────────

    def _style_toolbar_overflow(self) -> None:
        """Give the toolbar's overflow button ("…", shown when the tools do
        not fit) an icon in the theme's text colour: Qt's own glyph is
        drawn dark and vanished on the dark theme."""
        tb = getattr(self, "_main_toolbar", None)
        if tb is None:
            return
        from PySide6.QtCore import QPointF
        from PySide6.QtGui import QColor, QIcon, QPainter, QPixmap

        btn = tb.findChild(QToolButton, "qt_toolbar_ext_button")
        if btn is None:
            return
        dark = self._session.theme == "dark"
        colour = QColor("#cdd6f4" if dark else "#4c4f69")
        # The button is narrow (about 12 px) and as tall as the toolbar:
        # a vertical ellipsis fits it; a horizontal row shrank to nothing.
        ratio = max(1.0, float(self.devicePixelRatioF()))
        w_px, h_px = 6, 20
        pix = QPixmap(int(w_px * ratio), int(h_px * ratio))
        pix.setDevicePixelRatio(ratio)
        pix.fill(Qt.GlobalColor.transparent)
        p = QPainter(pix)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(colour)
        for cy in (4.0, 10.0, 16.0):
            p.drawEllipse(QPointF(3.0, cy), 2.0, 2.0)
        p.end()
        btn.setIconSize(QSize(w_px, h_px))
        btn.setIcon(QIcon(pix))
        btn.setToolTip("More tools")

    def _create_tool_bar(self) -> None:
        tb = self.addToolBar("Main")
        self._main_toolbar = tb
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
        _tb(
            "profile-view",
            "Profile",
            "Open Profile Viewer",
            lambda: self._show_window(self._profile_win),
        )
        _tb(
            "map-view",
            "Map",
            "Open Map Viewer",
            lambda: self._show_window(self._map_win),
        )
        tb.addSeparator()
        _tb(
            "qc",
            "QC",
            "Open QC Dashboard",
            lambda: self._show_window(self._qc_win),
        )
        _tb(
            "sites-correction",
            "Corrections",
            "Open Data Corrections",
            lambda: self._show_window(self._correction_win),
        )
        _tb(
            "forward",
            "Forward",
            "Open Forward Modelling",
            lambda: self._show_window(self._forward_win),
        )
        _tb(
            "inversion",
            "Inversion",
            "Open Inversion Wizard",
            self._open_inversion_wizard,
        )
        _tb(
            "interpret",
            "Interpret",
            "Open Interpretation Studio",
            lambda: self._show_window(self._interp_win),
        )
        _tb(
            "induction",
            "Airborne",
            "Open Airborne EM",
            lambda: self._show_window(self._airborne_win),
        )
        _tb(
            "pipeline",
            "Pipeline",
            "Open Processing Pipeline",
            lambda: self._show_window(self._pipeline_win),
        )
        tb.addSeparator()
        _tb("agents", "Agents", "Open Agent Master", self._open_agent_master)
        tb.addSeparator()

        # ── Secondary tools — listed directly rather than tucked behind a
        # "More" dropdown: the toolbar still has room for them, and an
        # overflow menu only earns its keep once space actually runs out.
        _tb(
            "tdem",
            "TDEM",
            "Open time-domain EM (TDEM) analysis panel",
            lambda: self._show_window(self._tdem_win),
        )
        _tb(
            "advanced-tools",
            "Advanced",
            "Open advanced EM processing and diagnostics",
            lambda: self._show_window(self._advanced_win),
        )
        # Full converter application. Unlike Tools > Format Converter, this
        # does not require a survey to be loaded in the desktop application.
        self._act_converter_app = _tb(
            "format-converter",
            "Converter",
            "Open the complete pycsamt format-conversion application",
            self._open_converter_app,
        )
        _tb(
            "export",
            "Export",
            "Save the active figure to PNG / PDF / SVG",
            self._on_export_figure,
        )
        tb.addSeparator()

        # Label shows the TARGET theme (what you get when you click),
        # not the current theme — so in dark mode it reads "☀ Light", etc.
        self._act_theme = QAction(
            "☀  Light" if self._session.theme == "dark" else "☾  Dark", self
        )
        self._act_theme.setStatusTip("Toggle dark / light theme")
        self._act_theme.triggered.connect(self._toggle_theme)
        tb.addAction(self._act_theme)
        self._style_toolbar_overflow()

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
        # Inversion runs keep going when their window is hidden; this chip
        # shows their progress and reopens the window on click.
        self._inv_chip = QPushButton("")
        self._inv_chip.setObjectName("InversionRunChip")
        self._inv_chip.setFlat(True)
        self._inv_chip.setCursor(Qt.CursorShape.PointingHandCursor)
        self._inv_chip.setToolTip("Inversion running — click to open the "
                                  "Inversion Studio")
        self._inv_chip.setStyleSheet(
            "QPushButton#InversionRunChip { background: #1864ab; color: "
            "white; border-radius: 8px; padding: 1px 10px; "
            "font-weight: 600; }")
        self._inv_chip.setVisible(False)
        self._inv_chip.clicked.connect(
            lambda: self._show_window(self._inversion_win))
        sb.addPermanentWidget(self._inv_chip)
        # Log toggle with an unread count (amber: warnings, red: errors)
        self._log_chip = QPushButton("Log")
        self._log_chip.setObjectName("LogChip")
        self._log_chip.setCheckable(True)
        self._log_chip.setFlat(True)
        self._log_chip.setCursor(Qt.CursorShape.PointingHandCursor)
        self._log_chip.setStyleSheet(
            "QPushButton#LogChip { border: 1px solid #8a94a3; "
            "border-radius: 8px; padding: 0px 10px; background: transparent;"
            " }"
            "QPushButton#LogChip:checked { background: #1864ab; color: white;"
            " border-color: #1864ab; }"
            "QPushButton#LogChip[level=\"warn\"] { background: #b7791f; "
            "color: white; border-color: #b7791f; }"
            "QPushButton#LogChip[level=\"error\"] { background: #c0392b; "
            "color: white; border-color: #c0392b; font-weight: 600; }")
        self._log_chip.clicked.connect(self.set_log_visible)
        sb.addPermanentWidget(self._log_chip)
        self._update_log_chip()
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

        # Detail card action buttons
        self._detail_card.open_profile_requested.connect(
            self._on_open_profile_for_station
        )
        self._detail_card.show_on_map_requested.connect(self._on_show_on_map)

        # Map window station pick → controller
        self._map_win.station_selected.connect(ctrl.select_station)

    # ── Data loading ──────────────────────────────────────────────────

    def _open_files(self) -> None:
        from pycsamt.app.desktop.dialogs.load_data_dlg import (
            LoadDataDialog,
        )

        if self._load_in_progress:
            self.statusBar().showMessage("A survey is already loading", 3500)
            return
        dlg = LoadDataDialog(
            self,
            last_dir=self._session.last_data_dir,
            recomputed_dir=self._last_recompute_output,
            existing_count=(
                len(self._all_dataframe)
                if self._all_dataframe is not None else 0
            ),
        )
        dlg.open_format_studio_requested.connect(self._open_converter_app)
        if dlg.exec() != LoadDataDialog.DialogCode.Accepted:
            return
        paths = dlg.selected_paths
        if not paths:
            return
        self._start_loading(paths, mode=getattr(dlg, "load_mode", "replace"))

    def _start_loading(self, paths: list, *, mode: str | None = None) -> None:
        from pycsamt.app.desktop.workers.loader_worker import (
            LoaderWorker,
        )

        if self._load_in_progress or (
            self._loader is not None
            and getattr(self._loader, "isRunning", lambda: False)()
        ):
            self.statusBar().showMessage("A survey is already loading", 3500)
            return
        if not paths:
            return
        if mode is None:
            mode = "append" if self._all_sites is not None else "replace"
        if mode not in {"append", "replace"}:
            raise ValueError("Load mode must be append or replace")
        self._load_mode = mode
        self._pending_load_paths = list(paths)
        self._load_in_progress = True
        self._loader = LoaderWorker(paths, parent=self)
        self._loader.progress.connect(self._progress_bar.setValue)
        self._loader.finished.connect(self._on_loader_finished)
        self._loader.error.connect(self._on_load_error)
        self._progress_bar.setValue(0)
        self._progress_bar.setVisible(True)
        self._status_ready_lbl.setText("Loading…")
        self._log(f"Loading {len(paths)} file(s)…")
        self._loader.start()

    def _on_loader_finished(self, sites) -> None:
        """Commit a successful batch without replacing state on read failure."""
        incoming = self._loader.data_controller
        append = self._load_mode == "append" and self._all_sites is not None
        try:
            if append:
                report = incoming.prepend_existing(
                    self._all_sites, self._all_dataframe
                )
                if report["added"]:
                    sites = incoming.sites
                    # Preserve old line choices and activate newly added lines.
                    old_lines = set(self._all_dataframe["Line"].astype(str))
                    new_lines = set(incoming.dataframe["Line"].astype(str))
                    self._session.active_lines = sorted(
                        self._active_line_names() | (new_lines - old_lines)
                    )
                    self._loaded_paths = list(dict.fromkeys(
                        getattr(self, "_loaded_paths", [])
                        + incoming.source_paths
                    ))
                added = report["added"]
                message = (
                    f"Added {added} station{'s' if added != 1 else ''}; "
                    f"skipped {len(report['skipped'])} matching station IDs."
                )
                if report["skipped"]:
                    self._log(
                        "Kept existing stations: "
                        + ", ".join(report["skipped"])
                    )
            else:
                self._loaded_paths = list(self._pending_load_paths)
                message = f"Loaded {len(sites)} stations."
            self._append_on_commit = append
            if not append or report["added"]:
                self._controller.set_sites(sites)
            self._session.last_data_dir = str(
                Path(self._pending_load_paths[0]).parent
            )
            for path in self._pending_load_paths:
                self._controller.add_recent_file(path)
            self._rebuild_recent_menu()
            self._progress_bar.setVisible(False)
            self._status_ready_lbl.setText("Ready  ●")
            self._log(message)
            self.statusBar().showMessage(message, 8000)
        except Exception as exc:
            self._on_load_error(str(exc))
        finally:
            self._append_on_commit = False
            self._load_in_progress = False
            self._pending_load_paths = []

    def _on_data_loaded(self, sites) -> None:
        # Fresh load → clear the local set; StationModel.set_dataframe() will
        # clear its own _recomputed_ids automatically.
        if not self._append_on_commit:
            self._recomputed_ids.clear()

        ctrl = self._controller
        df = None
        if self._loader is not None:
            df = self._loader.data_controller.dataframe
        before = self._survey_state()

        try:
            if df is not None:
                self._all_sites = sites
                if self._append_on_commit and before[0] is not None:
                    self._history.record("Add stations", before,
                                         (sites, self._lines_of(df)))
                else:
                    self._history.reset((sites, self._lines_of(df)))
                self._refresh_edit_state()
                self._all_dataframe = df.copy()
                self._configure_line_scope()
                self._apply_active_line_scope(initial=True)
                if self._append_on_commit:
                    self._station_panel._table._model.mark_recomputed(
                        self._recomputed_ids
                    )
                    if ctrl.selected_station:
                        ctrl.select_station(ctrl.selected_station)
        finally:
            # Always update the status bar — even if a panel window throws.
            n = ctrl.n_stations
            total = len(self._all_dataframe) if self._all_dataframe is not None else n
            self._status_file_lbl.setText(f"{n}/{total} stations active")
            self._progress_bar.setVisible(False)
            self._status_ready_lbl.setText("Ready  ●")
            self._log(f"Loaded {n} stations.")

    def _configure_line_scope(self) -> None:
        """Build the global active-line menu for the newly loaded survey."""
        df = self._all_dataframe
        if df is None or "Line" not in df.columns:
            lines = ["Survey"]
        else:
            line_names = set()
            for value in df["Line"].dropna():
                name = str(value).strip()
                if (
                    name
                    and name.lower() not in {"nan", "none"}
                    and any(character.isalnum() for character in name)
                ):
                    line_names.add(name)
            lines = sorted(line_names) or ["Survey"]

        restored = set(self._session.active_lines) & set(lines)
        active = restored or set(lines)
        self._line_scope_menu.clear()
        self._line_actions.clear()

        all_action = QAction("All lines", self, checkable=True)
        all_action.setChecked(active == set(lines))
        all_action.triggered.connect(self._on_all_lines_toggled)
        self._line_scope_menu.addAction(all_action)
        self._all_lines_action = all_action
        self._line_scope_menu.addSeparator()

        counts = (
            df["Line"].astype(str).value_counts().to_dict()
            if df is not None and "Line" in df.columns
            else {"Survey": len(df) if df is not None else 0}
        )
        for line in lines:
            action = QAction(
                f"{line}  ({counts.get(line, 0)} stations)",
                self,
                checkable=True,
            )
            action.setData(line)
            action.setChecked(line in active)
            action.toggled.connect(self._on_line_scope_toggled)
            self._line_scope_menu.addAction(action)
            self._line_actions[line] = action

        self._primary_line_combo.blockSignals(True)
        self._primary_line_combo.clear()
        self._primary_line_combo.addItems(sorted(active))
        primary = self._session.primary_line
        if primary not in active:
            primary = sorted(active)[0]
        self._primary_line_combo.setCurrentText(primary)
        self._primary_line_combo.blockSignals(False)

    def _active_line_names(self) -> set[str]:
        return {
            line for line, action in self._line_actions.items() if action.isChecked()
        }

    def _on_all_lines_toggled(self, checked: bool) -> None:
        if not checked:
            # "All lines" is a select-all command, not a way to disable the
            # complete survey. Individual entries control exclusions.
            self._all_lines_action.blockSignals(True)
            self._all_lines_action.setChecked(
                len(self._active_line_names()) == len(self._line_actions)
            )
            self._all_lines_action.blockSignals(False)
            return
        for action in self._line_actions.values():
            action.blockSignals(True)
            action.setChecked(True)
            action.blockSignals(False)
        self._apply_active_line_scope()

    def _on_line_scope_toggled(self, checked: bool) -> None:
        active = self._active_line_names()
        if not active:
            action = self.sender()
            if isinstance(action, QAction):
                action.blockSignals(True)
                action.setChecked(True)
                action.blockSignals(False)
            self.statusBar().showMessage("At least one survey line must remain active", 3500)
            return
        self._all_lines_action.blockSignals(True)
        self._all_lines_action.setChecked(len(active) == len(self._line_actions))
        self._all_lines_action.blockSignals(False)
        self._apply_active_line_scope()

    def _on_primary_line_changed(self, line: str) -> None:
        if not line:
            return
        self._session.primary_line = line
        self._controller.primary_line = line
        for window in self._panel_windows():
            window.setProperty("primaryLine", line)
        self.statusBar().showMessage(f"Primary survey line: {line}", 3000)

    def _apply_active_line_scope(self, *, initial: bool = False) -> None:
        """Filter once, then distribute the same active dataset everywhere."""
        if self._all_dataframe is None or self._all_sites is None:
            return
        active = self._active_line_names()
        df = self._all_dataframe
        if "Line" in df.columns and "Survey" not in active:
            active_df = df[df["Line"].astype(str).isin(active)].copy()
        else:
            active_df = df.copy()
        wanted = set(active_df["ID"].astype(str))
        from pycsamt.site.base import Sites

        active_sites = Sites(
            [
                site
                for site in self._all_sites
                if str(getattr(site, "name", getattr(site, "station", "")))
                in wanted
            ]
        )
        primary = self._primary_line_combo.currentText()
        if primary not in active:
            primary = sorted(active)[0]
            self._primary_line_combo.blockSignals(True)
            self._primary_line_combo.clear()
            self._primary_line_combo.addItems(sorted(active))
            self._primary_line_combo.setCurrentText(primary)
            self._primary_line_combo.blockSignals(False)
        else:
            current = self._primary_line_combo.currentText()
            self._primary_line_combo.blockSignals(True)
            self._primary_line_combo.clear()
            self._primary_line_combo.addItems(sorted(active))
            self._primary_line_combo.setCurrentText(current)
            self._primary_line_combo.blockSignals(False)

        self._controller.set_active_scope(active_sites, active, primary)
        if self._controller.selected_station not in wanted:
            self._controller.selected_station = None
            self._session.selected_station = None
            self._detail_card.clear()
            self._stats_tabs.setCurrentWidget(self._survey_overview)
        self._station_panel.set_dataframe(active_df)
        self._survey_overview.update_survey(
            active_df, getattr(self, "_loaded_paths", None)
        )
        if initial:
            self._stats_tabs.setCurrentWidget(self._survey_overview)
        for call in (
            lambda: self._correction_win.set_sites(active_sites),
            lambda: self._profile_win.set_sites(active_sites, active_df),
            lambda: self._map_win.set_dataframe(active_df),
            lambda: self._map_win.set_sites(active_sites),
            lambda: self._pcsf3d_win.set_sites(active_sites),
            lambda: self._qc_win.set_sites(active_sites),
            lambda: self._advanced_win.set_sites(active_sites),
            lambda: self._pipeline_win.set_input_sites(active_sites),
            lambda: self._forward_win.set_observed_sites(active_sites),
            lambda: self._inversion_win.set_sites(active_sites),
            lambda: self._interp_win.set_sites(active_sites),
        ):
            try:
                call()
            except Exception:
                pass
        for window in self._panel_windows():
            window.setProperty("primaryLine", primary)
        total = len(self._all_dataframe)
        count = len(active_df)
        self._line_scope_btn.setText(
            f"Lines: {len(active)}/{len(self._line_actions)}"
        )
        self._status_file_lbl.setText(f"{count}/{total} stations active")
        self._log(
            f"Active lines: {', '.join(sorted(active))} â€” {count}/{total} stations."
        )

    # ── Station selection ──────────────────────────────────────────────

    def _on_station_selected(self, station_id: str) -> None:
        """Update detail card; propagate to open windows."""
        sites = self._controller.sites
        if sites and station_id:
            self._detail_card.update_station(station_id, sites)
            self._stats_tabs.setCurrentWidget(self._detail_card)
        if self._profile_win.isVisible():
            self._profile_win.set_station(station_id)
        if self._map_win.isVisible():
            self._map_win.highlight_station(station_id)
        if (
            station_id
            and self._act_profile_on_select.isChecked()
            and not self._profile_win.isVisible()
        ):
            self._on_open_profile_for_station(station_id)

    def _on_profile_on_select_toggled(self, checked: bool) -> None:
        """Persist the optional station-click Profile Viewer behaviour."""
        self._session.open_profile_on_station_select = checked

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
            self._fit_panel_to_screen(win)
        win.show()
        win.raise_()
        win.activateWindow()

    def _fit_panel_to_screen(self, win: QWidget) -> None:
        """Clamp and center a panel using the main-window sizing policy."""
        screen = win.screen() or self.screen() or QApplication.primaryScreen()
        if screen is None:
            win.resize(min(win.width(), 1180), min(win.height(), 700))
            return

        available = screen.availableGeometry()
        max_width = max(
            win.minimumWidth(),
            min(1280, int(available.width() * 0.75), available.width() - 40),
        )
        max_height = max(
            win.minimumHeight(),
            min(760, int(available.height() * 0.75), available.height() - 40),
        )
        preferred_min_width = min(760, max_width)
        preferred_min_height = min(480, max_height)
        width = min(max(win.width(), preferred_min_width), max_width)
        height = min(max(win.height(), preferred_min_height), max_height)

        win.setWindowState(Qt.WindowState.WindowNoState)
        x = available.x() + (available.width() - width) // 2
        y = available.y() + (available.height() - height) // 2
        win.setGeometry(x, y, width, height)

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
        for attr in (
            "_profile_win",
            "_map_win",
            "_pcsf3d_win",
            "_qc_win",
            "_correction_win",
            "_advanced_win",
            "_tdem_win",
            "_pipeline_win",
            "_forward_win",
            "_inversion_win",
            "_interp_win",
            "_airborne_win",
            "_solver_builder_win",
        ):
            w = getattr(self, attr, None)
            if w is not None:
                wins.append(w)
        return wins

    # ── Correction → Main data bridge ────────────────────────────────

    @staticmethod
    def _site_name(site) -> str:
        return str(getattr(site, "name", getattr(site, "station", "")))

    def _merge_active_result(self, changed_sites):
        """Replace active stations while retaining globally disabled lines."""
        if self._all_sites is None or not self._controller.active_lines:
            return changed_sites
        from pycsamt.site.base import Sites

        changed = {self._site_name(site): site for site in changed_sites}
        merged = []
        for site in self._all_sites:
            name = self._site_name(site)
            merged.append(changed.pop(name, site))
        merged.extend(changed.values())
        return Sites(merged)

    @staticmethod
    def _lines_of(frame) -> dict:
        from pycsamt.app.desktop.controllers.station_edits import (
            lines_from_frame,
        )

        return lines_from_frame(frame)

    def _survey_state(self) -> tuple:
        """``(sites, {station: line})`` -- what an undo step restores
        (lines live in the station table, not in the sites)."""
        return self._all_sites, self._lines_of(self._all_dataframe)

    def _adopt_full_dataset(self, sites, *, label: str = "Edit",
                            record: bool = True, lines: dict | None = None,
                            merge: bool = True) -> bool:
        """Install a changed dataset and preserve its existing line labels.

        Every change goes through here, so it is recorded as one Edit ▸
        Undo step (*label*); undo / redo themselves pass ``record=False``.
        *lines* replaces the station -> line map (renames, line edits);
        ``merge=False`` when *sites* is already the whole survey.
        """
        before = self._survey_state()
        from pycsamt.app.desktop.controllers.data_controller import DataController

        old_lines = {}
        if self._all_dataframe is not None and "Line" in self._all_dataframe:
            old_lines = dict(
                zip(
                    self._all_dataframe["ID"].astype(str),
                    self._all_dataframe["Line"].astype(str),
                )
            )
        if merge:
            sites = self._merge_active_result(sites)
        dc = DataController()
        dc._sites = sites
        dc._df = dc._build_dataframe()
        df = dc.dataframe
        if df is None or df.empty:
            return False
        if lines is not None:
            old_lines = dict(lines)
        if old_lines:
            restored = df["ID"].astype(str).map(old_lines)
            df.loc[restored.notna(), "Line"] = restored[restored.notna()]
        self._controller.all_sites = sites
        self._all_sites = sites
        self._all_dataframe = df.copy()
        self._configure_line_scope()
        self._apply_active_line_scope()
        if record:
            self._history.record(label, before, self._survey_state())
        self._refresh_edit_state()
        return True

    def _on_corrections_committed(self, corrected_sites) -> None:
        """Replace global sites with the corrected dataset from CorrectionWindow."""
        self._adopt_full_dataset(corrected_sites, label="Data corrections")
        n = self._controller.n_stations
        self._status_file_lbl.setText(f"{n} stations (corrected)")
        self._log(f"Corrected dataset committed — {n} stations.")

    def _open_solver_builder(self, key: str | None = None) -> None:
        """Show the Solver Builder, optionally on a given solver."""
        self._show_window(self._solver_builder_win)
        if key:
            self._solver_builder_win.open_solver(key)

    def _on_pipeline_applied(self, processed_sites) -> None:
        """Adopt the Pipeline Studio output as the global dataset."""
        self._adopt_full_dataset(processed_sites,
                                 label="Processing pipeline")
        n = self._controller.n_stations
        self._status_file_lbl.setText(f"{n} stations (pipeline)")
        self._log(f"Pipeline output applied to main data — {n} stations.")

    # ── Conversion → Main data bridge ────────────────────────────────

    def _on_conversion_committed(self, converted_sites) -> None:
        """Replace global sites with an EDICollection converted by AdvancedToolsWindow."""
        try:
            from pycsamt.site.base import to_sites

            converted_sites = to_sites(converted_sites)
        except Exception as exc:
            self._log(f"Could not commit converted dataset: {exc}")
            return

        self._adopt_full_dataset(converted_sites, label="Format conversion")
        n = self._controller.n_stations
        self._status_file_lbl.setText(f"{n} stations (converted)")
        self._log(f"Converted dataset committed — {n} stations.")

    # ── Forward → Inversion bridge ────────────────────────────────────

    def _on_forward_send_to_inversion(self, payload: dict) -> None:
        """Show the InversionWindow pre-loaded with the forward model."""
        self._log(
            f"Forward model sent to Inversion (dim={payload.get('dim', '1D')})."
        )
        self._inversion_win.load_starting_model(payload)
        self._show_window(self._inversion_win)

    # ── Inversion ─────────────────────────────────────────────────────

    def _open_inversion_wizard(self) -> None:
        self._show_window(self._inversion_win)

    def _on_inversion_run_state(self, text: str, percent: int,
                                running: bool) -> None:
        """Mirror an Inversion Studio run in the status bar chip."""
        chip = getattr(self, "_inv_chip", None)
        if chip is None:
            return
        if not running:
            chip.setVisible(False)
            return
        pct = "" if percent < 0 else f"  {percent}%"
        chip.setText(f"⟳ {text or 'Inversion'}{pct}")
        chip.setVisible(True)

    def _on_inversion_result_ready(self, payload: dict) -> None:
        """Forward a completed inversion model to the Interpretation Studio."""
        if payload.get("result") is None and not payload.get("path"):
            return
        self._interp_win.receive_inversion(payload)
        self._log("Inversion model forwarded to Interpretation Studio.")

    # ── Export ────────────────────────────────────────────────────────

    def _on_export_figure(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )

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
        from pycsamt.app.desktop.dialogs.preferences_dlg import (
            PreferencesDialog,
        )
        from pycsamt.app.desktop.licensing import get_default_manager

        dlg = PreferencesDialog(
            session=self._session,
            parent=self,
            license_manager=get_default_manager(),
        )
        if dlg.exec() == PreferencesDialog.DialogCode.Accepted:
            self._apply_theme(self._session.theme)
            self._detail_card.set_animations(self._session.ui_animations)
            self._act_theme.setText(
                "☀  Light" if self._session.theme == "dark" else "☾  Dark"
            )
            self._session.save()
            self._log("Preferences saved.")

    # ── API Configuration / Settings ──────────────────────────────────

    def _open_api_config(self, tab: str | None = None) -> None:
        """Open the API Configuration dialog, optionally pre-selecting *tab*."""
        from pycsamt.app.desktop.dialogs.settings_dialog import (
            APIConfigDialog,
        )

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
        if any(
            k in touched
            for k in (
                "station",
                "section",
                "view_controls",
                "contour",
                "mesh",
                "style",
                "interpretation",
            )
        ):
            try:
                if self._profile_win.isVisible():
                    self._profile_win._on_refresh()
            except Exception:
                pass

        # Pseudosection rendering → redraw QC Dashboard current plot
        if any(
            k in touched
            for k in ("station", "section", "contour", "style", "interpretation")
        ):
            try:
                if self._qc_win.isVisible():
                    self._qc_win._on_run()
            except Exception:
                pass

        if "view_controls" in touched:
            from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

            for canvas in self.findChildren(MplCanvas):
                canvas.draw()

        keys_str = ", ".join(touched)
        self._log(f"Settings applied: {keys_str}")

    def _reset_all_settings(self) -> None:
        """Reset all PYCSAMT_* singletons to defaults and refresh panels."""
        from PySide6.QtWidgets import QMessageBox

        reply = QMessageBox.question(
            self,
            "Reset All Settings",
            "Reset ALL API configuration to package defaults?\n\n"
            "This will revert all display, export, rendering, ordering, and "
            "topography "
            "settings to their original values.\n\nThis action cannot be undone.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        self._settings_ctrl.reset_all()
        self._on_settings_changed(
            [
                "station", "section", "view_controls", "topography", "style",
                "interpretation", "plot", "contour", "mesh", "ordering",
            ]
        )
        self._log("All API settings reset to package defaults.")

    def _save_settings_profile(self) -> None:
        from PySide6.QtWidgets import QFileDialog

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Settings Profile",
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
            self,
            "Load Settings Profile",
            str(self._settings_ctrl.SETTINGS_PATH.parent),
            "JSON files (*.json)",
        )
        if path:
            ok = self._settings_ctrl.load(path)
            if ok:
                self._on_settings_changed(
                    [
                        "station", "section", "view_controls", "topography",
                        "style", "interpretation", "plot", "contour", "mesh",
                        "ordering",
                    ]
                )
                self._log(f"Settings profile loaded ← {path}")
            else:
                self._log(f"Could not load settings profile: {path}")

    # ── Guard: require loaded survey data ────────────────────────────────

    def _require_sites(self, tool_name: str = "") -> bool:
        """Return True if survey data is loaded; else show NoDataDialog.

        If the user clicks *Load Data* in the dialog, the open-file
        action is triggered automatically and this method returns False
        (the caller should abort — the user will retry after loading).
        """
        sites = getattr(self._controller, "sites", None)
        if sites is not None:
            return True
        from pycsamt.app.desktop.dialogs.no_data_dialog import (
            NoDataDialog,
        )

        if NoDataDialog.require(self, tool_name):
            self._act_open.trigger()
        return False

    # ── Tools menu handlers ───────────────────────────────────────────

    def _open_strike_analyzer(self) -> None:
        if not self._require_sites("Strike Analyzer"):
            return
        from pycsamt.app.desktop.tools.strike_tool import (
            StrikeAnalyzerDialog,
        )

        StrikeAnalyzerDialog(
            getattr(self._controller, "sites", None), parent=self
        ).exec()

    def _open_edi_validator(self) -> None:
        if not self._require_sites("EDI Validator & Station Manager"):
            return
        from pycsamt.app.desktop.tools.validator_tool import (
            EDIValidatorDialog,
        )

        dlg = EDIValidatorDialog(
            getattr(self._controller, "sites", None), parent=self
        )
        dlg.open_recompute_requested.connect(self._open_recompute)
        if dlg.exec() and dlg.modified_sites is not None:
            self._apply_modified_sites(
                dlg.modified_sites, source="EDI Validator"
            )

    def _open_format_converter(self) -> None:
        if not self._require_sites("Format Converter"):
            return
        from pycsamt.app.desktop.tools.converter_tool import (
            FormatConverterDialog,
        )

        dlg = FormatConverterDialog(
            getattr(self._controller, "sites", None), parent=self
        )
        dlg.open_format_studio_requested.connect(self._open_converter_app)
        dlg.exec()

    def _open_converter_app(self) -> None:
        """Open the complete :mod:`pycsamt.app.converter` application."""
        try:
            if self._converter_app_window is None:
                from pycsamt.app.converter.main_window import (
                    ConverterMainWindow,
                )

                self._converter_app_window = ConverterMainWindow(
                    embedded=True,
                    host_theme=self._session.theme,
                )
            self._show_window(self._converter_app_window)
            self._log("Converter App opened.")
        except Exception as exc:
            message = f"Could not open Converter App: {exc}"
            self._log(message)
            self.statusBar().showMessage(message, 7000)

    def _open_batch_export(self) -> None:
        from pycsamt.app.desktop.tools.batch_export_tool import (
            BatchExportDialog,
        )

        figures = self._collect_figures()
        BatchExportDialog(figures, parent=self).exec()

    def _open_coord_transformer(self) -> None:
        from pycsamt.app.desktop.tools.coord_tool import (
            CoordTransformDialog,
        )

        CoordTransformDialog(
            getattr(self._controller, "sites", None), parent=self
        ).exec()

    def _collect_figures(self) -> list:
        """Return [(label, Figure)] for every visible canvas in all panel windows."""
        figures = []
        _LABEL = {
            "_profile_win": "Profile",
            "_map_win": "Map",
            "_pcsf3d_win": "PCSF 3D",
            "_qc_win": "QC",
            "_correction_win": "Correction",
            "_advanced_win": "Advanced",
            "_tdem_win": "TDEM",
            "_pipeline_win": "Pipeline",
            "_forward_win": "Forward",
            "_inversion_win": "Inversion",
            "_interp_win": "Interpretation",
            "_airborne_win": "Airborne",
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
        from pycsamt.app.desktop.tools.station_response_tool import (
            StationResponseDialog,
        )

        StationResponseDialog(
            getattr(self._controller, "sites", None), parent=self
        ).exec()

    def _open_strike_profile(self) -> None:
        if not self._require_sites("Strike Profile Viewer"):
            return
        from pycsamt.app.desktop.tools.strike_profile_tool import (
            StrikeProfileDialog,
        )

        StrikeProfileDialog(
            getattr(self._controller, "sites", None), parent=self
        ).exec()

    def _open_phase_tensor_map(self) -> None:
        if not self._require_sites("Phase Tensor Map"):
            return
        from pycsamt.app.desktop.tools.phase_tensor_map_tool import (
            PhaseTensorMapDialog,
        )

        PhaseTensorMapDialog(
            getattr(self._controller, "sites", None), parent=self
        ).exec()

    def _open_phase_tensor_strip_grid(self) -> None:
        if not self._require_sites("Phase Tensor Strip Grid"):
            return
        from pycsamt.app.desktop.tools.phase_tensor_strip_grid_tool import (
            PhaseTensorStripGridDialog,
        )

        PhaseTensorStripGridDialog(
            getattr(self._controller, "sites", None), parent=self
        ).exec()

    def _open_dimensionality(self) -> None:
        if not self._require_sites("Dimensionality Classifier"):
            return
        from pycsamt.app.desktop.tools.dimensionality_tool import (
            DimensionalityDialog,
        )

        DimensionalityDialog(
            getattr(self._controller, "sites", None), parent=self
        ).exec()

    def _open_frequency_editor(self) -> None:
        if not self._require_sites("Frequency Editor"):
            return
        from pycsamt.app.desktop.tools.frequency_editor_tool import (
            FrequencyEditorDialog,
        )

        dlg = FrequencyEditorDialog(
            getattr(self._controller, "sites", None), parent=self
        )
        if dlg.exec() and dlg.edited_sites is not None:
            self._apply_modified_sites(
                dlg.edited_sites, source="Frequency Editor"
            )

    def _open_layered_model(self) -> None:
        from pycsamt.app.desktop.tools.layered_model_tool import (
            LayeredModelDialog,
        )

        LayeredModelDialog(parent=self).exec()

    def _open_elevation_enrichment(self) -> None:
        if not self._require_sites("Elevation Enrichment"):
            return
        from pycsamt.app.desktop.tools.elevation_tool import (
            ElevationEnrichDialog,
        )

        ElevationEnrichDialog(
            getattr(self._controller, "sites", None), parent=self
        ).exec()

    def _open_recompute(self) -> None:
        from PySide6.QtWidgets import QMessageBox

        from pycsamt.app.desktop.dialogs.recompute_dlg import (
            RecomputeDialog,
        )

        sites = getattr(self._controller, "sites", None)

        # If already recomputed — warn and ask for confirmation before re-opening
        if self._recomputed_ids:
            n = len(self._recomputed_ids)
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
            f"Recomputed {n} station(s). Marked with ◈ in the station list.",
            6000,
        )

    def _apply_modified_sites(self, new_sites, *, source: str = "") -> None:
        """Push a modified site list back into the controller and refresh all panels."""
        tag = f"{source}: " if source else ""
        if new_sites is None or (hasattr(new_sites, "__len__")
                                 and len(new_sites) == 0):
            # (an empty result merged with the untouched stations and
            # was logged as "0 station(s) applied")
            self._log(f"Warning: {tag}the edited survey is empty; "
                      "nothing was applied.")
            return
        try:
            ok = self._adopt_full_dataset(new_sites, label=source or "Edit")
        except Exception as exc:  # was swallowed: say it did not apply
            self._log(f"ERROR: {tag}could not apply the edit — {exc}")
            return
        if not ok:
            self._log(f"Warning: {tag}the edited survey is empty; "
                      "nothing was applied.")
            return
        n = len(new_sites) if hasattr(new_sites, "__len__") else "?"
        self._log(f"{tag}{n} station(s) applied to survey.")

    # ── Edit menu ─────────────────────────────────────────────────────

    def _build_edit_menu(self, mb, act_prefs) -> None:
        """Edit: undo / redo / history, find, and the survey editors."""
        edit = mb.addMenu("&Edit")
        self._edit_menu = edit
        self._act_undo = QAction(_icon("reset"), "&Undo", self)
        self._act_undo.setShortcuts([QKeySequence.StandardKey.Undo])
        self._act_undo.triggered.connect(self.undo_edit)
        self._act_redo = QAction("&Redo", self)
        self._act_redo.setShortcuts(
            [QKeySequence("Ctrl+Y"), QKeySequence("Ctrl+Shift+Z")])
        self._act_redo.triggered.connect(self.redo_edit)
        self._act_history = QAction(_icon("log"), "&History…", self)
        self._act_history.setShortcut("Ctrl+H")
        self._act_history.setStatusTip("Every change to the survey; jump "
                                       "back to any of them")
        self._act_history.triggered.connect(self._open_history)
        self._act_revert = QAction("Re&vert to As Loaded", self)
        self._act_revert.setStatusTip(
            "Back to the survey as it was loaded (can itself be undone)")
        self._act_revert.triggered.connect(self.revert_to_loaded)
        for act in (self._act_undo, self._act_redo, self._act_history):
            edit.addAction(act)
        edit.addAction(self._act_revert)
        edit.addSeparator()
        self._act_find = QAction(_icon("tools"), "&Find Station…", self)
        self._act_find.setShortcut(QKeySequence.StandardKey.Find)
        self._act_find.setStatusTip("Jump to a station by name")
        self._act_find.triggered.connect(self._find_station)
        edit.addAction(self._act_find)
        edit.addSeparator()
        self._edit_station_menu = edit.addMenu(_icon("station-response"),
                                               "&Stations")
        for text, tab, key, tip in (
                ("Edit &Coordinates…", "coordinates", "Ctrl+Alt+C",
                 "Latitude / longitude / elevation or easting / northing; "
                 "paste from a spreadsheet or import a CSV"),
                ("&Rename Stations…", "names", "Ctrl+Alt+R",
                 "Prefix, find & replace, case, zero-padding"),
                ("Assign &Lines…", "lines", "Ctrl+Alt+N",
                 "Detect lines from station names, rename lines"),
                ("&Header && Metadata…", "header", "Ctrl+Alt+M",
                 "EDI header: acquired by, dates, project, location…")):
            act = QAction(text, self)
            act.setShortcut(key)
            act.setStatusTip(tip)
            act.triggered.connect(lambda _c=False, t=tab:
                                  self.open_station_editor(t))
            self._edit_station_menu.addAction(act)
        self._edit_station_menu.addSeparator()
        self._edit_freq_menu = edit.addMenu(_icon("frequency-editor"),
                                            "F&requencies")
        act_points = QAction(_icon("station-response"), "&Point Editor…",
                             self)
        act_points.setShortcut("Ctrl+Alt+E")
        act_points.setStatusTip(
            "Mask, delete, interpolate or restore data points and apply "
            "static shifts, station by station")
        act_points.triggered.connect(self.open_point_editor)
        self._edit_freq_menu.addAction(act_points)
        self._edit_freq_menu.addSeparator()
        self._edit_tensor_menu = edit.addMenu(_icon("impendance-tensor"),
                                              "&Tensor")
        from pycsamt.app.desktop.controllers.survey_ops import OPS

        for o in OPS:
            act = QAction(o.label + "…", self)
            act.setStatusTip(o.help)
            act.triggered.connect(lambda _c=False, k=o.key:
                                  self.open_survey_ops(k))
            (self._edit_freq_menu if o.menu == "frequencies"
             else self._edit_tensor_menu).addAction(act)
        self._edit_freq_menu.addSeparator()
        self._edit_tensor_menu.addSeparator()
        edit.addSeparator()
        edit.addAction(act_prefs)
        self._all_icon_actions.extend([(self._act_undo, "reset"),
                                       (self._act_history, "log"),
                                       (self._act_find, "tools")])
        self._refresh_edit_state()

    def open_station_editor(self, tab: str = "coordinates") -> None:
        if not self._require_sites("Station Editor"):
            return
        from pycsamt.app.desktop.dialogs.station_editor import (
            StationEditorDialog,
        )

        sites, lines = self._survey_state()
        dlg = StationEditorDialog(sites, lines, tab=tab, parent=self)
        if dlg.exec() and dlg.changes:
            self.apply_station_changes(dlg.changes, dlg.change_label())

    def open_survey_ops(self, op: str = "trim") -> None:
        if not self._require_sites("Frequency & Tensor Tools"):
            return
        from pycsamt.app.desktop.dialogs.survey_ops import SurveyOpsDialog

        sites, _lines = self._survey_state()
        dlg = SurveyOpsDialog(sites, op=op, parent=self)
        if dlg.exec() and dlg.result is not None:
            self.apply_survey_op(dlg.result, dlg.label())

    def apply_survey_op(self, sites, label: str) -> bool:
        """Adopt a whole-survey operation's result as one undoable edit."""
        _old, lines = self._survey_state()
        if not self._adopt_full_dataset(sites, label=label, lines=lines,
                                        merge=False):
            self._log(f"Warning: {label} — nothing was applied.")
            return False
        self._log(f"{label} — applied to {self._controller.n_stations} "
                  "station(s).")
        return True

    def open_point_editor(self, station: str | None = None) -> None:
        if not self._require_sites("Point Editor"):
            return
        from pycsamt.app.desktop.dialogs.point_editor import (
            PointEditorDialog,
        )

        sites, _lines = self._survey_state()
        dlg = PointEditorDialog(
            sites, station=station or self._controller.selected_station,
            parent=self)
        if dlg.exec():
            self.apply_point_edits(dlg.edited_sites(),
                                   dlg.edited_stations())

    def apply_point_edits(self, sites, stations: list[str]) -> bool:
        """Adopt the Point Editor result as one undoable edit."""
        if not stations:
            return False
        _old, lines = self._survey_state()
        label = (f"Point Editor: {stations[0]}" if len(stations) == 1
                 else f"Point Editor: {len(stations)} stations")
        if not self._adopt_full_dataset(sites, label=label, lines=lines,
                                        merge=False):
            self._log(f"Warning: {label} — nothing was applied.")
            return False
        self._log(f"{label} — edits written to each station's INFO "
                  "(PYCSAMT_POINT_EDITS).")
        return True

    def apply_station_changes(self, changes: dict,
                              label: str = "Edit stations") -> bool:
        """Apply Station Editor changes as one undoable edit."""
        from pycsamt.app.desktop.controllers.station_edits import (
            apply_changes,
            carry_lines,
        )

        sites, lines = self._survey_state()
        try:
            new_sites, n, _line_changes = apply_changes(sites, changes)
        except Exception as exc:
            self._log(f"ERROR: {label} — {exc}")
            return False
        new_lines = carry_lines(lines, changes)
        if not self._adopt_full_dataset(new_sites, label=label,
                                        lines=new_lines, merge=False):
            self._log(f"Warning: {label} — nothing was applied.")
            return False
        self._log(f"{label} — {n} station(s) changed.")
        return True

    def _refresh_edit_state(self) -> None:
        """Undo / Redo labels and enabled state, and the title's ●."""
        h = self._history
        if hasattr(self, "_act_undo"):
            self._act_undo.setEnabled(h.can_undo)
            self._act_undo.setText(
                f"&Undo {h.undo_label}" if h.can_undo else "&Undo")
            self._act_redo.setEnabled(h.can_redo)
            self._act_redo.setText(
                f"&Redo {h.redo_label}" if h.can_redo else "&Redo")
            self._act_history.setEnabled(len(h) > 0)
            self._act_revert.setEnabled(len(h) > 0)
            has = self._all_sites is not None
            self._act_save_survey.setEnabled(has)
            self._act_find.setEnabled(has)
        self.setWindowTitle("pycsamt" + ("  — ● edited" if h.dirty
                                          and self._all_sites is not None
                                          else ""))

    def _restore_survey(self, state, message: str) -> None:
        sites, lines = state if isinstance(state, tuple) else (state, None)
        if sites is None or not self._adopt_full_dataset(
                sites, record=False, lines=lines, merge=False):
            self._log(f"Warning: {message} — nothing to restore.")
            return
        n = self._controller.n_stations
        self._status_file_lbl.setText(f"{n} stations")
        self._log(message)
        self.statusBar().showMessage(message, 3000)

    def undo_edit(self) -> None:
        if not self._history.can_undo:
            return
        label = self._history.undo_label
        self._restore_survey(self._history.undo(), f"Undo: {label}")

    def redo_edit(self) -> None:
        if not self._history.can_redo:
            return
        label = self._history.redo_label
        self._restore_survey(self._history.redo(), f"Redo: {label}")

    def jump_to_history(self, n_applied: int) -> None:
        """Go to the survey with the first *n_applied* steps applied."""
        self._restore_survey(self._history.jump(n_applied),
                             f"History: {n_applied} step(s) applied")

    def revert_to_loaded(self) -> None:
        original = self._history.original()
        if original is None:
            return
        sites, lines = (original if isinstance(original, tuple)
                        else (original, None))
        # recorded, so the revert itself can be undone
        if self._adopt_full_dataset(sites, label="Revert to as loaded",
                                    lines=lines, merge=False):
            self._log("Reverted to the survey as loaded.")

    def _open_history(self) -> None:
        from PySide6.QtWidgets import (
            QDialog,
            QDialogButtonBox,
            QListWidget,
            QListWidgetItem,
        )

        dlg = QDialog(self)
        dlg.setWindowTitle("Edit history")
        dlg.resize(460, 380)
        v = QVBoxLayout(dlg)
        info = QLabel("Every change to the survey, oldest first. Pick a step "
                      "to go back (or forward) to it; later steps stay "
                      "available until you make a new change.")
        info.setWordWrap(True)
        info.setObjectName("InfoLabel")
        v.addWidget(info)
        lst = QListWidget()
        start = QListWidgetItem("◦  As loaded")
        lst.addItem(start)
        for time, label, applied in self._history.entries():
            it = QListWidgetItem(f"{'●' if applied else '○'}  {time}  {label}")
            if not applied:
                it.setForeground(Qt.GlobalColor.gray)
            lst.addItem(it)
        lst.setCurrentRow(self._history.position)
        v.addWidget(lst, 1)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        go = bb.addButton("Go to step", QDialogButtonBox.ButtonRole.AcceptRole)
        bb.rejected.connect(dlg.reject)

        def apply():
            row = lst.currentRow()
            if row >= 0 and row != self._history.position:
                self.jump_to_history(row)
            dlg.accept()

        go.clicked.connect(apply)
        lst.itemDoubleClicked.connect(lambda _it: apply())
        v.addWidget(bb)
        dlg.exec()

    def _station_ids(self) -> list[str]:
        df = self._all_dataframe
        if df is None or "ID" not in df:
            return []
        return [str(x) for x in df["ID"].tolist()]

    def find_station(self, name: str) -> bool:
        """Select *name* (case-insensitive) everywhere; ``True`` if found."""
        ids = self._station_ids()
        match = next((i for i in ids if i == name), None) or next(
            (i for i in ids if i.lower() == name.strip().lower()), None)
        if match is None:
            return False
        self._controller.select_station(match)
        try:
            self._station_panel.highlight_station(match)
        except Exception:
            pass
        return True

    def _find_station(self) -> None:
        from PySide6.QtWidgets import QInputDialog

        ids = self._station_ids()
        if not ids:
            return
        name, ok = QInputDialog.getItem(self, "Find station", "Station:",
                                        ids, 0, True)
        if ok and name and not self.find_station(name):
            self.statusBar().showMessage(f"No station named {name!r}", 4000)

    def save_edited_survey(self, folder, fmt: str = "edi") -> list:
        """Write the current survey to *folder* (``edi`` or ``xml``)."""
        from pathlib import Path

        from pycsamt.site.base import to_sites

        sites = to_sites(self._all_sites)
        folder = Path(folder)
        paths = (sites.write_xml(folder) if fmt == "xml"
                 else sites.write(folder, exist_ok=True))
        self._history.mark_saved()
        self._refresh_edit_state()
        self._log(f"Edited survey saved — {len(paths)} {fmt.upper()} "
                  f"file(s) in {folder}")
        return paths

    def _save_edited_survey(self) -> None:
        from PySide6.QtWidgets import QFileDialog, QMessageBox

        if self._all_sites is None:
            return
        folder = QFileDialog.getExistingDirectory(
            self, "Save edited survey to a folder")
        if not folder:
            return
        box = QMessageBox(self)
        box.setWindowTitle("Save edited survey")
        box.setText("Write the survey as:")
        edi = box.addButton("EDI files", QMessageBox.ButtonRole.AcceptRole)
        xml = box.addButton("EMTF-XML files", QMessageBox.ButtonRole.AcceptRole)
        box.addButton(QMessageBox.StandardButton.Cancel)
        box.exec()
        if box.clickedButton() not in (edi, xml):
            return
        try:
            self.save_edited_survey(folder,
                                    "xml" if box.clickedButton() is xml
                                    else "edi")
        except Exception as exc:
            self._log(f"ERROR: could not save the edited survey — {exc}")

    # ── Help menu handlers ────────────────────────────────────────────

    def _open_documentation(self) -> None:
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        QDesktopServices.openUrl(QUrl("https://pycsamt.org/"))

    def _open_github(self) -> None:
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        QDesktopServices.openUrl(
            QUrl("https://github.com/earthai-tech/pycsamt")
        )

    def _open_about(self) -> None:
        from pycsamt.app.desktop.dialogs.about_dialog import (
            AboutDialog,
        )
        from pycsamt.app.desktop.licensing import get_default_manager

        AboutDialog(parent=self, license_manager=get_default_manager()).exec()

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
        self._session.dock_state = self.saveState().toBase64().data().decode()
        for win in self._panel_windows():
            save_geometry = getattr(win, "save_geometry_to", None)
            if callable(save_geometry):
                save_geometry(self._session.window_geometries)

    def _restore_layout(self) -> None:
        if self._session.dock_geometry:
            try:
                self.restoreGeometry(
                    QByteArray.fromBase64(self._session.dock_geometry.encode())
                )
            except Exception:
                pass
        # Old sessions may contain maximized or nearly full-screen geometry.
        # Clamp after restoring so upgrades immediately get a modern default.
        self._fit_to_screen()
        if self._session.dock_state:
            try:
                self.restoreState(
                    QByteArray.fromBase64(self._session.dock_state.encode())
                )
            except Exception:
                pass
        # The saved dock state of older sessions always had the log open;
        # the log's own flag (default: hidden) decides.
        self._log_dock.setVisible(bool(self._session.log_visible))

    # ── Error handling ────────────────────────────────────────────────

    def _on_load_error(self, message: str) -> None:
        self._load_in_progress = False
        self._pending_load_paths = []
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
            + ("The survey has edits that are not saved to disk "
               "(File ▸ Save Edited Survey).\n\n"
               if self._history.dirty and self._all_sites is not None
               else "")
            + "Any unsaved session data will be lost.",
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
        if self._converter_app_window is not None:
            self._converter_app_window.close()
        self._log("Goodbye.")
        super().closeEvent(event)
