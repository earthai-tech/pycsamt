# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
APIConfigDialog — Settings dialog for all PYCSAMT_* API singletons.

Layout
------
┌──────────────────┬──────────────────────────────────────────┐
│ 🔍 Search…        │  Pseudosections                          │
│                  │  Station markers and section axis style. │
│ PLOTS            │ ──────────────────────────────────────── │
│   Pseudosections │                                          │
│   View Controls  │        <active page, scrollable>         │
│   Display        │                                          │
│   Contours & Mesh│                                          │
│ DATA             │                                          │
│   Topography     │                                          │
│   …              │                                          │
├──────────────────┴──────────────────────────────────────────┤
│  [Reset page]                          [OK] [Cancel] [Apply]│
└─────────────────────────────────────────────────────────────┘

Sidebar navigation replaced a tab strip: with nine sections the tabs no
longer fit (the last ones hid behind scroll arrows with clipped titles),
and fixed-height tabs squeezed tall pages until spin boxes and combos
overlapped. Every page now sits in its own scroll area, so no page is
ever compressed, whatever the dialog or screen size.

Signals
-------
settings_changed(list[str])
    Emitted after Apply or OK with a list of apply-keys that were written
    (e.g. ``["station", "section", "view_controls"]``).  MainWindow
    connects to this to refresh only the affected panels.

Deep-linking
------------
Call ``dialog.open_tab("pseudosections")`` (or pass ``open_tab=``) to open
the dialog on a specific page:

    from pycsamt.app.desktop.dialogs.settings_dialog import APIConfigDialog
    dlg = APIConfigDialog(ctrl, parent=self, open_tab="pipeline")
    dlg.settings_changed.connect(handler)
    dlg.exec()
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont, QGuiApplication
from PySide6.QtWidgets import (
    QAbstractButton,
    QDialog,
    QDialogButtonBox,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QScrollArea,
    QSplitter,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.settings_controller import (
    SettingsController,
)
from pycsamt.app.desktop.widgets.settings_pages import (
    DisplayPage,
    InterpretationPage,
    OrderingPage,
    OutputPage,
    PipelinePage,
    PseudosectionsPage,
    RenderingPage,
    TopographyPage,
    ViewControlsPage,
)

# Ordered page registry: (key, label, PageClass, group, one-line description)
_TABS = [
    ("pseudosections", "Pseudosections", PseudosectionsPage, "Plots",
     "Station markers and pseudosection axis orientation."),
    ("view_controls", "View Controls", ViewControlsPage, "Plots",
     "Resistivity/phase display, profile x-axis and plot-panel behaviour."),
    ("display", "Display", DisplayPage, "Plots",
     "Colours and line widths of MT components and corrections."),
    ("rendering", "Contours & Mesh", RenderingPage, "Plots",
     "Contour overlays and mesh-review styling."),
    ("topography", "Topography", TopographyPage, "Data",
     "Terrain-following geometry for sections and plots."),
    ("interpretation", "Interpretation", InterpretationPage, "Data",
     "Colour maps and markers of interpreted sections."),
    ("ordering", "Data Ordering", OrderingPage, "Data",
     "How loaders order survey sites when no order is requested."),
    ("output", "Plot & Export", OutputPage, "Output",
     "Default figure format, resolution and export folder."),
    ("pipeline", "Pipeline", PipelinePage, "Output",
     "Processing-pipeline results, error handling, reports, cache and "
     "run history."),
]

# Map key → page index (for open_tab look-ups)
_TAB_INDEX: dict[str, int] = {key: i for i, (key, *_rest) in enumerate(_TABS)}

_RESET_KEYS: dict[str, list[str]] = {
    "pseudosections": ["station", "section"],
    "view_controls": ["view_controls"],
    "display": ["style"],
    "topography": ["topography"],
    "interpretation": ["interpretation"],
    "output": ["plot"],
    "rendering": ["contour", "mesh"],
    "ordering": ["ordering"],
    "pipeline": ["pipe"],
}

_PAGE_ROLE = Qt.ItemDataRole.UserRole  # nav item -> page index (None: header)
_GROUP_ROLE = Qt.ItemDataRole.UserRole + 1  # header item -> group name


class APIConfigDialog(QDialog):
    """
    Sidebar-navigated dialog for configuring all PYCSAMT_* API singletons.

    Parameters
    ----------
    ctrl : SettingsController
        Shared controller instance (owned by MainWindow).
    parent : QWidget, optional
    open_tab : str, optional
        If given, the dialog opens on this page (a key of ``_TABS``, e.g.
        ``"pseudosections"`` or ``"pipeline"``).
    """

    settings_changed = Signal(list)  # list[str] of apply-method keys touched

    def __init__(
        self,
        ctrl: SettingsController,
        parent: QWidget | None = None,
        open_tab: str | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("API Configuration")
        self.setMinimumSize(640, 440)
        self._ctrl = ctrl
        self._snapshot = ctrl.snapshot()  # for Cancel
        self._pages: list = []
        self._build_ui()
        self._fit_to_screen(860, 600)
        self.open_tab(open_tab or _TABS[0][0])

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(8)

        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)
        split.setHandleWidth(4)

        # ── Sidebar: search + grouped navigation ──────────────────────
        side = QWidget()
        sv = QVBoxLayout(side)
        sv.setContentsMargins(0, 0, 0, 0)
        sv.setSpacing(6)
        self._search = QLineEdit()
        self._search.setPlaceholderText("Search settings…")
        self._search.setClearButtonEnabled(True)
        self._search.textChanged.connect(self._filter_nav)
        sv.addWidget(self._search)

        self._nav = QListWidget()
        self._nav.setObjectName("SettingsNav")
        self._nav.setMinimumWidth(170)
        self._nav.setMaximumWidth(230)
        # Roomier rows than the app-wide list style: this is a navigation
        # list, read at a glance (local rule, wins over the theme QSS).
        self._nav.setStyleSheet(
            "QListWidget#SettingsNav::item { padding: 5px 4px; }"
        )
        sv.addWidget(self._nav, 1)
        split.addWidget(side)

        # ── Content: header + stacked scrollable pages ────────────────
        content = QWidget()
        cv = QVBoxLayout(content)
        cv.setContentsMargins(6, 0, 0, 0)
        cv.setSpacing(4)
        self._page_title = QLabel("")
        self._page_title.setObjectName("SettingsPageTitle")
        # stylesheet, not setFont(): the theme's QLabel rule overrides fonts
        self._page_title.setStyleSheet(
            "QLabel#SettingsPageTitle { font-size: 15px; font-weight: 600; }"
        )
        cv.addWidget(self._page_title)
        self._page_desc = QLabel("")
        self._page_desc.setObjectName("InfoLabel")
        self._page_desc.setWordWrap(True)
        cv.addWidget(self._page_desc)
        rule = QFrame()
        rule.setFrameShape(QFrame.Shape.HLine)
        rule.setObjectName("Separator")
        cv.addWidget(rule)

        self._stack = QStackedWidget()
        cv.addWidget(self._stack, 1)
        split.addWidget(content)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([190, 640])
        root.addWidget(split, 1)

        group = None
        for i, (_key, label, PageClass, grp, _desc) in enumerate(_TABS):
            if grp != group:
                group = grp
                self._nav.addItem(self._header_item(grp))
            page = PageClass(self)
            self._pages.append(page)
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QFrame.Shape.NoFrame)
            scroll.setHorizontalScrollBarPolicy(
                Qt.ScrollBarPolicy.ScrollBarAlwaysOff
            )
            scroll.setWidget(page)
            self._stack.addWidget(scroll)
            item = QListWidgetItem(f"   {label}")
            item.setData(_PAGE_ROLE, i)
            self._nav.addItem(item)

        self._nav.currentItemChanged.connect(self._on_nav_changed)
        self._stack.currentChanged.connect(self._sync_header)

        # ── Button row ────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        self._reset_btn = QPushButton("Reset page")
        self._reset_btn.setToolTip(
            "Reset this page's settings to package defaults"
        )
        self._reset_btn.clicked.connect(self._on_reset_tab)
        btn_row.addWidget(self._reset_btn)
        btn_row.addStretch()

        box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Apply
        )
        box.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(
            self._on_apply
        )
        box.accepted.connect(self._on_ok)
        box.rejected.connect(self._on_cancel)
        btn_row.addWidget(box)
        root.addLayout(btn_row)

    @staticmethod
    def _header_item(text: str) -> QListWidgetItem:
        item = QListWidgetItem(text.upper())
        item.setFlags(Qt.ItemFlag.NoItemFlags)  # not selectable
        f = QFont()
        f.setBold(True)
        f.setPointSizeF(max(f.pointSizeF() * 0.85, 7.0))
        item.setFont(f)
        item.setData(_PAGE_ROLE, None)
        item.setData(_GROUP_ROLE, text)
        return item

    def _fit_to_screen(self, width: int, height: int) -> None:
        screen = self.screen() or QGuiApplication.primaryScreen()
        if screen is not None:
            avail = screen.availableGeometry()
            width = min(width, int(avail.width() * 0.9))
            height = min(height, int(avail.height() * 0.9))
        self.resize(max(width, self.minimumWidth()),
                    max(height, self.minimumHeight()))

    # ── Navigation ────────────────────────────────────────────────────────────

    def _nav_item_for(self, index: int) -> QListWidgetItem | None:
        for row in range(self._nav.count()):
            item = self._nav.item(row)
            if item.data(_PAGE_ROLE) == index:
                return item
        return None

    def _on_nav_changed(self, current, _previous) -> None:
        if current is None:
            return
        index = current.data(_PAGE_ROLE)
        if index is not None:
            self._stack.setCurrentIndex(index)

    def _sync_header(self, index: int) -> None:
        if not 0 <= index < len(_TABS):
            return
        _key, label, _cls, _grp, desc = _TABS[index]
        self._page_title.setText(label)
        self._page_desc.setText(desc)
        item = self._nav_item_for(index)
        if item is not None and self._nav.currentItem() is not item:
            self._nav.blockSignals(True)
            self._nav.setCurrentItem(item)
            self._nav.blockSignals(False)

    def _filter_nav(self, text: str) -> None:
        """Show only pages whose title or any visible setting label matches."""
        query = text.strip().lower()
        visible_groups: set[str] = set()
        for i, (_key, label, _cls, grp, desc) in enumerate(_TABS):
            match = not query or query in label.lower() or query in desc.lower()
            if not match:
                match = any(query in t.lower() for t in self._page_texts(i))
            self._nav_item_for(i).setHidden(not match)
            if match:
                visible_groups.add(grp)
        for row in range(self._nav.count()):
            item = self._nav.item(row)
            if item.data(_PAGE_ROLE) is None:
                item.setHidden(item.data(_GROUP_ROLE) not in visible_groups)
        # Jump to the first match so the results are immediately visible
        current = self._nav.currentItem()
        if current is None or current.isHidden():
            for row in range(self._nav.count()):
                item = self._nav.item(row)
                if item.data(_PAGE_ROLE) is not None and not item.isHidden():
                    self._nav.setCurrentItem(item)
                    break

    def _page_texts(self, index: int) -> list[str]:
        page = self._pages[index]
        texts = [w.text() for w in page.findChildren(QLabel)]
        texts += [w.text() for w in page.findChildren(QAbstractButton)]
        texts += [w.title() for w in page.findChildren(QGroupBox)]
        return [t for t in texts if t]

    # ── Public API ────────────────────────────────────────────────────────────

    def open_tab(self, key: str) -> None:
        """Switch to the page identified by *key* (e.g. ``"pipeline"``)."""
        idx = _TAB_INDEX.get(key)
        if idx is not None:
            self._stack.setCurrentIndex(idx)
            self._sync_header(idx)

    def current_key(self) -> str:
        """Key of the page currently shown."""
        return _TABS[self._stack.currentIndex()][0]

    def page_labels(self) -> list[str]:
        """Page titles in navigation order."""
        return [label for _key, label, *_rest in _TABS]

    # ── Slots ─────────────────────────────────────────────────────────────────

    def _on_apply(self) -> None:
        """Apply all pages and emit settings_changed with the touched keys."""
        touched: list[str] = []
        for page in self._pages:
            try:
                kw_map = page.collect()
            except Exception:
                kw_map = {}
            for apply_key, fields in kw_map.items():
                if not fields:
                    continue
                method = getattr(self._ctrl, f"apply_{apply_key}", None)
                if method:
                    try:
                        method(**fields)
                        if apply_key not in touched:
                            touched.append(apply_key)
                    except Exception:
                        pass
        if touched:
            self.settings_changed.emit(touched)

    def _on_ok(self) -> None:
        self._on_apply()
        self.accept()

    def _on_cancel(self) -> None:
        """Restore all singletons from the pre-dialog snapshot."""
        try:
            self._ctrl.restore(self._snapshot)
        except Exception:
            pass
        self.reject()

    def _on_reset_tab(self) -> None:
        """Reset the active page's singletons and re-populate it."""
        idx = self._stack.currentIndex()
        key = _TABS[idx][0]
        page = self._pages[idx]
        try:
            page.reset()
        except Exception:
            pass
        # Emit so live panels refresh
        self.settings_changed.emit(_RESET_KEYS.get(key, [key]))
