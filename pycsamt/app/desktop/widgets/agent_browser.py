# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AgentBrowserWidget — category-tabbed card browser for pycsamt agents.

Layout
------
  ┌─────────────────────────────────────────┐
  │  🔍  Search agents…                     │  ← QLineEdit
  │  [All][QC][Pre-proc][Inv.][AI Inv.][…]  │  ← QTabBar (scrollable chips)
  ├─────────────────────────────────────────┤
  │  ▌ QC Analysis              [LLM]       │  ←  _AgentCard (54 px)
  │  ▌ Quality-control scoring via…         │
  │  ▌ Anomaly Detection        [LLM]       │
  │    Unsupervised CAE anomaly…            │
  │  ▌ QC Quicklook             [FAST]      │
  │    Fast QC overview plot…               │
  │  …  (scrollable)                        │
  └─────────────────────────────────────────┘

Signals
-------
agent_selected(str)   — emitted on single-click / keyboard selection
agent_activated(str)  — emitted on double-click (Run immediately)
"""

from __future__ import annotations

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QSizePolicy,
    QTabBar,
    QVBoxLayout,
    QWidget,
)

# ── colour constants ──────────────────────────────────────────────────────────

_LLM_COLOR  = "#f0a500"   # amber  — LLM agents
_PROC_COLOR = "#00c896"   # teal   — fast processing agents

# Short label displayed inside each category tab
_CAT_ABBREV: dict[str, str] = {
    "Data Loading":               "Load",
    "Quality Control":            "QC",
    "Pre-processing":             "Pre-proc",
    "Phase & Dimensionality":     "Dimension",
    "Forward Modelling":          "Forward",
    "Inversion":                  "Inversion",
    "AI Inversion":               "AI Inv",
    "Post-processing":            "Post",
    "Interpretation & Reporting": "Interp",
    "Workflow":                   "Workflow",
    "⚡ Processing":               "⚡ Fast",
}

# ── Card stylesheet (light + dark safe, uses relative colours) ────────────────

_BROWSER_QSS = """
/* ── search bar ──────────────────────────────────────────────── */
QLineEdit#AgentSearch {
    border-radius: 4px;
    padding: 4px 8px;
    font-size: 12px;
}

/* ── category tab bar ────────────────────────────────────────── */
QTabBar#CategoryTabBar {
    font-size: 11px;
}
QTabBar#CategoryTabBar::tab {
    padding: 3px 10px;
    margin-right: 2px;
    border-radius: 10px;
    min-width: 40px;
}
QTabBar#CategoryTabBar::tab:selected {
    background: #4c7af0;
    color: #ffffff;
    font-weight: bold;
}

/* ── agent card list ─────────────────────────────────────────── */
QListWidget#AgentCardList {
    border: none;
    background: transparent;
    outline: none;
}
QListWidget#AgentCardList::item {
    border-radius: 5px;
    padding: 0px;
}
QListWidget#AgentCardList::item:selected {
    background: rgba(100, 160, 255, 0.14);
    border: 1px solid rgba(100, 160, 255, 0.35);
}
QListWidget#AgentCardList::item:hover:!selected {
    background: rgba(255, 255, 255, 0.05);
}

/* ── card labels ─────────────────────────────────────────────── */
QLabel#CardName {
    font-weight: bold;
    font-size: 12px;
}
QLabel#CardDesc {
    font-size: 10px;
    color: #9a9aaa;
}
QLabel#LLMBadge {
    background: #f0a500;
    color: #1e1e2e;
    font-size: 9px;
    font-weight: bold;
    padding: 1px 6px;
    border-radius: 4px;
}
QLabel#ProcBadge {
    background: #00c896;
    color: #1e1e2e;
    font-size: 9px;
    font-weight: bold;
    padding: 1px 6px;
    border-radius: 4px;
}
"""


# ═════════════════════════════════════════════════════════════════════════════
# _AgentCard
# ═════════════════════════════════════════════════════════════════════════════

class _AgentCard(QFrame):
    """
    Single-row agent card.

    ┌────┬─────────────────────────────────────────────┐
    │ ▌  │  Agent Name                       [BADGE]   │  ← 54 px tall
    │    │  Short description, elided to fit           │
    └────┴─────────────────────────────────────────────┘
    """

    def __init__(
        self,
        name: str,
        description: str,
        is_llm: bool,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._name   = name
        self._is_llm = is_llm
        self.setObjectName("AgentCard")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFixedHeight(54)
        self.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._build(name, description, is_llm)

    # ── build ─────────────────────────────────────────────────────────────────

    def _build(self, name: str, desc: str, is_llm: bool) -> None:
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(0)

        # Left accent bar (4 px)
        bar = QFrame()
        bar.setFixedWidth(4)
        color = _LLM_COLOR if is_llm else _PROC_COLOR
        bar.setStyleSheet(f"background: {color}; border: none; border-radius: 0px;")
        row.addWidget(bar)

        # Content area
        content = QVBoxLayout()
        content.setContentsMargins(8, 6, 8, 6)
        content.setSpacing(2)

        # Name + badge row
        name_row = QHBoxLayout()
        name_row.setSpacing(6)
        name_lbl = QLabel(name)
        name_lbl.setObjectName("CardName")
        name_row.addWidget(name_lbl, 1)

        badge_text = "LLM" if is_llm else "FAST"
        badge = QLabel(badge_text)
        badge.setObjectName("LLMBadge" if is_llm else "ProcBadge")
        name_row.addWidget(badge)
        content.addLayout(name_row)

        # Description (single line, elided)
        desc_lbl = QLabel()
        desc_lbl.setObjectName("CardDesc")
        # Elide eagerly — 220 px is safe for a 300-px panel
        fm     = desc_lbl.fontMetrics()
        elided = fm.elidedText(desc, Qt.TextElideMode.ElideRight, 220)
        desc_lbl.setText(elided)
        desc_lbl.setToolTip(desc)
        content.addWidget(desc_lbl)

        row.addLayout(content)

    # ── public API ────────────────────────────────────────────────────────────

    @property
    def agent_name(self) -> str:
        return self._name


# ═════════════════════════════════════════════════════════════════════════════
# AgentBrowserWidget
# ═════════════════════════════════════════════════════════════════════════════

class AgentBrowserWidget(QWidget):
    """
    Search bar + category tab bar + scrollable agent card list.

    Signals
    -------
    agent_selected(str)   — fired when a card row is highlighted
    agent_activated(str)  — fired on double-click (caller may auto-run)
    """

    agent_selected  = Signal(str)
    agent_activated = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._cat_names:  list[str]              = []
        self._all_agents: dict[str, list[tuple]] = {}
        self._current_agent: str                 = ""
        self._build_ui()
        self._populate()
        self.setStyleSheet(_BROWSER_QSS)

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        # Search bar
        self._search = QLineEdit()
        self._search.setObjectName("AgentSearch")
        self._search.setPlaceholderText("🔍  Search agents…")
        self._search.setClearButtonEnabled(True)
        self._search.textChanged.connect(self._on_search)
        layout.addWidget(self._search)

        # Category tab bar (horizontal, scrollable chips)
        self._cat_bar = QTabBar()
        self._cat_bar.setObjectName("CategoryTabBar")
        self._cat_bar.setExpanding(False)
        self._cat_bar.setDocumentMode(False)
        self._cat_bar.setUsesScrollButtons(True)
        self._cat_bar.currentChanged.connect(self._on_cat_changed)
        layout.addWidget(self._cat_bar)

        # Agent card list
        self._card_list = QListWidget()
        self._card_list.setObjectName("AgentCardList")
        self._card_list.setSpacing(2)
        self._card_list.setFrameShape(QListWidget.Shape.NoFrame)
        self._card_list.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._card_list.setVerticalScrollMode(
            QListWidget.ScrollMode.ScrollPerPixel
        )
        self._card_list.currentRowChanged.connect(self._on_row_changed)
        self._card_list.itemDoubleClicked.connect(self._on_double_click)
        layout.addWidget(self._card_list, 1)

    # ── data population ───────────────────────────────────────────────────────

    def _populate(self) -> None:
        from pycsamt.app.desktop.agent_registry import (
            AGENT_REGISTRY,
            agents_by_category,
        )

        # Build category → [(name, entry)] ordered map
        cat_map: dict[str, list[tuple[str, dict]]] = {"All": []}

        for cat, names in agents_by_category().items():
            entries = [(n, AGENT_REGISTRY[n]) for n in names]
            cat_map[cat]    = entries
            cat_map["All"].extend(entries)

        proc = [
            (n, e) for n, e in AGENT_REGISTRY.items()
            if e.get("type") == "processing"
        ]
        if proc:
            cat_map["⚡ Processing"] = proc
            cat_map["All"].extend(proc)

        self._all_agents = cat_map
        self._cat_names  = list(cat_map.keys())

        # Add tabs
        for cat in self._cat_names:
            short = _CAT_ABBREV.get(cat, cat)
            count = len(cat_map[cat])
            self._cat_bar.addTab(short)
            idx = self._cat_bar.count() - 1
            self._cat_bar.setTabToolTip(idx, f"{cat}  ({count} agents)")

        # Trigger initial card refresh for "All"
        if self._cat_names:
            self._cat_bar.setCurrentIndex(0)
            self._refresh_cards("All", "")

    # ── card rendering ────────────────────────────────────────────────────────

    def _refresh_cards(self, cat: str, search: str) -> None:
        current = self._current_agent
        self._card_list.clear()

        entries = self._all_agents.get(cat, [])
        search  = search.strip().lower()

        for name, entry in entries:
            if search and (
                search not in name.lower()
                and search not in entry.get("description", "").lower()
            ):
                continue

            is_llm = entry.get("type") == "llm"
            desc   = entry.get("description", "")
            card   = _AgentCard(name, desc, is_llm)

            item = QListWidgetItem()
            item.setSizeHint(QSize(0, 56))   # height enforced; width flexible
            item.setData(Qt.ItemDataRole.UserRole, name)
            item.setToolTip(desc)

            self._card_list.addItem(item)
            self._card_list.setItemWidget(item, card)

            # Re-select the previously active agent
            if name == current:
                self._card_list.setCurrentItem(item)

    # ── slots ─────────────────────────────────────────────────────────────────

    def _on_cat_changed(self, idx: int) -> None:
        if 0 <= idx < len(self._cat_names):
            cat = self._cat_names[idx]
            self._refresh_cards(cat, self._search.text())

    def _on_search(self, text: str) -> None:
        # Jump to "All" tab while searching so nothing is hidden
        if text and self._cat_bar.currentIndex() != 0:
            self._cat_bar.blockSignals(True)
            self._cat_bar.setCurrentIndex(0)
            self._cat_bar.blockSignals(False)
        cat = self._cat_names[self._cat_bar.currentIndex()] \
              if self._cat_names else "All"
        self._refresh_cards(cat, text)

    def _on_row_changed(self, row: int) -> None:
        item = self._card_list.item(row)
        if item is None:
            return
        name = item.data(Qt.ItemDataRole.UserRole)
        if name:
            self._current_agent = name
            self.agent_selected.emit(name)

    def _on_double_click(self, item: QListWidgetItem) -> None:
        name = item.data(Qt.ItemDataRole.UserRole)
        if name:
            self.agent_activated.emit(name)

    # ── public API ────────────────────────────────────────────────────────────

    @property
    def current_agent(self) -> str:
        return self._current_agent

    def select_agent(self, name: str) -> None:
        """Programmatically select an agent and scroll it into view."""
        for i in range(self._card_list.count()):
            item = self._card_list.item(i)
            if item and item.data(Qt.ItemDataRole.UserRole) == name:
                self._card_list.setCurrentItem(item)
                self._card_list.scrollToItem(item)
                return

        # Agent not visible in current tab — switch to All and retry
        if self._cat_bar.currentIndex() != 0:
            self._cat_bar.setCurrentIndex(0)   # triggers _refresh_cards
            self.select_agent(name)

    def apply_theme(self, dark: bool) -> None:
        """Re-apply the stylesheet so badge/search colours adapt to theme."""
        # Currently the QSS works for both modes; call update to repaint
        self.update()
