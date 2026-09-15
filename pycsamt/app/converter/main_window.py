# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.main_window
====================================

The converter app's top-level window: a branded sidebar of conversion
tools next to a stacked page area -- same nav idea as before -- now
wrapped in real application chrome mirroring the full desktop app's
``MainWindow`` (menu bar, toolbar, status bar) plus one shared,
dockable log panel that every page reports to instead of each
embedding its own.
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
    QFrame,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.converter.pages import (
    BatchPage,
    EdiXmlPage,
    InversionPage,
    PcbhPage,
    PcglPage,
    PcgsPage,
    PcptPage,
    PcsfPcsmPage,
    SettingsPage,
)
from pycsamt.app.converter.settings import load_settings, save_settings
from pycsamt.app.converter.widgets import LogPanel

_RESOURCES = Path(__file__).parent / "resources"
_ICONS = _RESOURCES / "icons"

_TOOLS = (
    ("Inversion → PCSF/PCSM", InversionPage, "inversion.svg"),
    ("PCSF ⇄ PCSM / Validate / Info", PcsfPcsmPage, "format-converter.svg"),
    ("EDI ⇄ EMTF-XML", EdiXmlPage, "EDI-validator.svg"),
    ("Build PCBH (borehole)", PcbhPage, "layered-model.svg"),
    ("Build PCGL (geology legend)", PcglPage, "topography.svg"),
    ("Build PCGS (structure)", PcgsPage, "strike-analyzer.svg"),
    ("Build PCPT (points)", PcptPage, "coverage.svg"),
    ("Batch queue", BatchPage, "batch-export.svg"),
    ("Settings", SettingsPage, "tools.svg"),
)

_SETTINGS_ROW = next(i for i, t in enumerate(_TOOLS) if t[1] is SettingsPage)

# Tracks current theme so _icon() can recolor SVGs without being passed a flag.
_DARK_MODE: bool = False

# Matches any 6-digit hex colour in an SVG source string.
_HEX6_RE = re.compile(r"#[0-9a-fA-F]{6}", re.IGNORECASE)


def _is_near_black(hex6: str) -> bool:
    r, g, b = int(hex6[1:3], 16), int(hex6[3:5], 16), int(hex6[5:7], 16)
    return (0.299 * r + 0.587 * g + 0.114 * b) < 80


def _recolor_svg(text: str, target: str = "#f4f7fb") -> bytes:
    # Pass 1: replace near-black hex colours.
    result = _HEX6_RE.sub(
        lambda m: target if _is_near_black(m.group(0)) else m.group(0),
        text,
    )
    # Pass 2: replace the named keyword "black" in attribute values / inline styles.
    result = re.sub(
        r'(?<=[";\s:])black(?=[";\s])', target, result, flags=re.IGNORECASE
    )

    # Pass 3: inject fill on shape elements that carry no explicit fill, so
    # icons relying on the SVG default (black) become visible on a dark
    # background. Elements with an existing fill= or fill: are left alone.
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


class ConverterMainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self._settings = load_settings()
        self._theme = "light"
        self._tool_actions: list = []
        self._all_icon_actions: list = []

        self.setWindowTitle("pyCSAMT Format Studio")
        self.resize(920, 620)
        self.setMinimumSize(760, 520)
        ico = _ICONS / "pycsamt.ico"
        if ico.exists():
            self.setWindowIcon(QIcon(str(ico)))

        self._build_ui()

        theme = self._settings.theme if self._settings.theme in ("light", "dark") else "light"
        self._apply_theme(theme)

    def _build_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        sidebar = self._build_sidebar()
        root.addWidget(sidebar)

        self._stack = QStackedWidget()
        self._stack.setObjectName("PageSurface")
        self._pages = [page_cls() for _label, page_cls, _icon_name in _TOOLS]
        for page in self._pages:
            self._stack.addWidget(page)
        root.addWidget(self._stack, stretch=1)

        self._create_log_dock()
        self._create_menu_bar()
        self._create_status_bar()

        self._nav.currentRowChanged.connect(self._stack.setCurrentIndex)
        self._nav.currentRowChanged.connect(self._sync_tools_menu)
        self._nav.setCurrentRow(0)

        # Every page reports to one shared dock instead of an embedded
        # LogPanel -- avoids duplicating the same ~120px log block on every
        # single page.
        for page in self._pages:
            if hasattr(page, "log"):
                page.log = self._log_panel

    def _build_sidebar(self) -> QWidget:
        # No brand banner here -- the window title bar already reads
        # "pyCSAMT Format Studio", so repeating it atop the nav list would
        # just be redundant chrome.
        sidebar = QFrame()
        sidebar.setObjectName("Sidebar")
        sidebar.setFixedWidth(250)
        layout = QVBoxLayout(sidebar)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self._nav = QListWidget()
        self._nav.setObjectName("ToolNav")
        self._nav.setFrameShape(QFrame.Shape.NoFrame)
        self._nav.setIconSize(QSize(18, 18))
        self._nav.setSpacing(2)
        for label, _page_cls, icon_name in _TOOLS:
            item = QListWidgetItem(label)
            item.setIcon(_icon(icon_name))
            self._nav.addItem(item)
        layout.addWidget(self._nav, stretch=1)

        return sidebar

    # ── Log dock (bottom, thin, shared by every page) ──────────────────

    def _create_log_dock(self) -> None:
        self._log_panel = LogPanel()
        self._log_dock = QDockWidget("Log", self)
        self._log_dock.setObjectName("LogDock")
        self._log_dock.setAllowedAreas(Qt.DockWidgetArea.BottomDockWidgetArea)
        self._log_dock.setWidget(self._log_panel)
        self._log_dock.setMaximumHeight(140)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self._log_dock)

    # ── Menu bar (File | Tools | View | Help) ───────────────────────────

    def _create_menu_bar(self) -> None:
        mb = self.menuBar()

        # ── File ─────────────────────────────────────────────────────
        file_menu = mb.addMenu("&File")

        act_settings = QAction(_icon("tools"), "&Settings…", self)
        act_settings.setShortcut("Ctrl+,")
        act_settings.setStatusTip("Persisted defaults for new conversion jobs")
        act_settings.triggered.connect(lambda: self._nav.setCurrentRow(_SETTINGS_ROW))
        file_menu.addAction(act_settings)

        file_menu.addSeparator()

        act_quit = QAction("&Quit", self)
        act_quit.setShortcut(QKeySequence.StandardKey.Quit)
        act_quit.triggered.connect(self.close)
        file_menu.addAction(act_quit)

        self._all_icon_actions.append((act_settings, "tools"))

        # ── Tools -- one entry per page, mirrors the sidebar nav ────────
        tools_menu = mb.addMenu("&Tools")
        for i, (label, _page_cls, icon_name) in enumerate(_TOOLS):
            act = QAction(_icon(icon_name), label, self)
            act.setCheckable(True)
            act.setShortcut(f"Ctrl+{i + 1}")
            act.triggered.connect(lambda _=False, row=i: self._nav.setCurrentRow(row))
            tools_menu.addAction(act)
            self._tool_actions.append(act)
            self._all_icon_actions.append((act, icon_name))
        self._tool_actions[0].setChecked(True)

        # ── View ─────────────────────────────────────────────────────
        view_menu = mb.addMenu("&View")

        _log_tva = self._log_dock.toggleViewAction()
        view_menu.addAction(_log_tva)
        view_menu.addSeparator()

        theme_menu = view_menu.addMenu("Theme")
        self._act_light = QAction("☀  Light", self, checkable=True)
        self._act_dark = QAction("☾  Dark", self, checkable=True)
        theme_grp = QActionGroup(self)
        theme_grp.setExclusive(True)
        theme_grp.addAction(self._act_light)
        theme_grp.addAction(self._act_dark)
        self._act_light.triggered.connect(lambda: self._apply_theme("light"))
        self._act_dark.triggered.connect(lambda: self._apply_theme("dark"))
        theme_menu.addAction(self._act_light)
        theme_menu.addAction(self._act_dark)

        # ── Help ─────────────────────────────────────────────────────
        help_menu = mb.addMenu("&Help")

        act_docs = QAction("&Documentation", self)
        act_docs.setStatusTip("Open pycsamt.org in your browser")
        act_docs.triggered.connect(self._open_documentation)
        help_menu.addAction(act_docs)

        act_gh = QAction("pyCSAMT on &GitHub", self)
        act_gh.setStatusTip("Open the pycsamt GitHub repository in your browser")
        act_gh.triggered.connect(self._open_github)
        help_menu.addAction(act_gh)

        help_menu.addSeparator()

        act_about = QAction("&About pyCSAMT Format Studio", self)
        act_about.triggered.connect(self._open_about)
        help_menu.addAction(act_about)

    # ── Status bar ────────────────────────────────────────────────────

    def _create_status_bar(self) -> None:
        sb = self.statusBar()
        sb.addPermanentWidget(QLabel("Ready  ●"))

    # ── Nav <-> Tools menu sync ──────────────────────────────────────

    def _sync_tools_menu(self, row: int) -> None:
        for i, act in enumerate(self._tool_actions):
            act.setChecked(i == row)

    # ── Theme ────────────────────────────────────────────────────────

    def _apply_theme(self, theme: str) -> None:
        global _DARK_MODE
        _DARK_MODE = theme == "dark"
        self._theme = theme

        qss = _RESOURCES / f"{theme}_theme.qss"
        if qss.exists():
            QApplication.instance().setStyleSheet(qss.read_text(encoding="utf-8"))

        if hasattr(self, "_act_dark"):
            self._act_dark.setChecked(theme == "dark")
            self._act_light.setChecked(theme == "light")

        # Re-ice every icon-bearing action (File/Tools menus) + sidebar nav.
        for action, icon_name in self._all_icon_actions:
            action.setIcon(_icon(icon_name))
        if hasattr(self, "_nav"):
            for row, (_label, _page_cls, icon_name) in enumerate(_TOOLS):
                item = self._nav.item(row)
                if item is not None:
                    item.setIcon(_icon(icon_name))

        if self._settings.theme != theme:
            self._settings.theme = theme
            save_settings(self._settings)

    # ── Help menu handlers ───────────────────────────────────────────

    def _open_documentation(self) -> None:
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        QDesktopServices.openUrl(QUrl("https://pycsamt.org/"))

    def _open_github(self) -> None:
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        QDesktopServices.openUrl(QUrl("https://github.com/earthai-tech/pycsamt"))

    def _open_about(self) -> None:
        QMessageBox.about(
            self,
            "About pyCSAMT Format Studio",
            "<b>pyCSAMT Format Studio</b><br>"
            "Version 2.0<br><br>"
            "Convert inversion results, AI/DL array bundles, EDI/EMTF-XML "
            "files, and build PCBH/PCGL/PCGS/PCPT documents -- the "
            "standalone GUI face of <code>pycsamt format</code>.<br><br>"
            "© earthai-tech — "
            '<a href="https://pycsamt.org/">pycsamt.org</a>',
        )
