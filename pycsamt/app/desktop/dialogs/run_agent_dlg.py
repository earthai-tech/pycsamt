# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
RunAgentDialog — agent picker with a dynamically generated parameter form.

Phase 4: Two-column layout:
  Left   — QListWidget listing all registered agents, grouped by type
  Right  — QStackedWidget with one QFormLayout page per agent
  Bottom — API-key field (for LLM agents), OK / Cancel
"""

from __future__ import annotations

from typing import Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPlainTextEdit,
    QSpinBox,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.agent_registry import (
    agent_names,
    agents_by_category,
    get_entry,
)

# ──────────────────────────────────────────────────────────────────────────────
# Per-agent parameter page
# ──────────────────────────────────────────────────────────────────────────────

class _ParamPage(QWidget):
    """A QFormLayout page for one agent's parameters."""

    def __init__(self, entry: dict, parent=None) -> None:
        super().__init__(parent)
        self._entry    = entry
        self._widgets: dict[str, QWidget] = {}
        self._build()

    def _build(self) -> None:
        form = QFormLayout(self)
        form.setContentsMargins(8, 8, 8, 8)
        form.setSpacing(8)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        # Description
        desc = self._entry.get("description", "")
        if desc:
            lbl = QLabel(desc)
            lbl.setWordWrap(True)
            lbl.setObjectName("AgentDescLabel")
            form.addRow(lbl)

        for param_name, spec in self._entry.get("params", {}).items():
            label = spec.get("label", param_name)
            tip   = spec.get("tip", "")
            w     = self._make_widget(spec)
            if tip:
                w.setToolTip(tip)
            self._widgets[param_name] = w
            form.addRow(f"{label}:", w)

        if not self._entry.get("params"):
            form.addRow(QLabel("<i>No parameters for this operation.</i>"))

    @staticmethod
    def _make_widget(spec: dict) -> QWidget:
        kind = spec.get("type", "str")
        if kind == "combo":
            w = QComboBox()
            w.addItems(spec.get("options", []))
            default = spec.get("default", "")
            idx = w.findText(str(default))
            if idx >= 0:
                w.setCurrentIndex(idx)
            return w
        if kind == "float":
            lo, hi = spec.get("range", (0.0, 1e9))
            w = QDoubleSpinBox()
            w.setRange(lo, hi)
            w.setSingleStep(spec.get("step", 0.1))
            w.setDecimals(3)
            w.setValue(float(spec.get("default", lo)))
            return w
        if kind == "int":
            lo, hi = spec.get("range", (0, 10_000))
            w = QSpinBox()
            w.setRange(lo, hi)
            w.setSingleStep(spec.get("step", 1))
            w.setValue(int(spec.get("default", lo)))
            return w
        if kind == "bool":
            w = QCheckBox()
            w.setChecked(bool(spec.get("default", False)))
            return w
        # default: string
        w = QLineEdit()
        w.setText(str(spec.get("default", "")))
        return w

    def current_params(self) -> dict[str, Any]:
        """Return {param_name: current_value} dict."""
        result: dict[str, Any] = {}
        for name, widget in self._widgets.items():
            if isinstance(widget, QComboBox):
                result[name] = widget.currentText()
            elif isinstance(widget, (QDoubleSpinBox, QSpinBox)):
                result[name] = widget.value()
            elif isinstance(widget, QCheckBox):
                result[name] = widget.isChecked()
            elif isinstance(widget, QLineEdit):
                result[name] = widget.text()
            else:
                result[name] = None
        return result


# ──────────────────────────────────────────────────────────────────────────────
# Main dialog
# ──────────────────────────────────────────────────────────────────────────────

class RunAgentDialog(QDialog):
    """
    Modal dialog for selecting and configuring a pycsamt agent.

    After ``exec()`` returns ``QDialog.Accepted``, read:
      ``selected_agent`` — display name string
      ``params``         — {param_name: value} dict
      ``api_key``        — API key string (may be empty)
    """

    def __init__(
        self,
        api_key: str | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Run Agent / Processing")
        self.setMinimumSize(680, 440)

        self.selected_agent: str = ""
        self.params:         dict[str, Any] = {}
        self.api_key:        str = api_key or ""

        self._pages: dict[str, _ParamPage] = {}
        self._build_ui()
        # Select first item by default
        if self._list.count() > 0:
            self._list.setCurrentRow(0)

    # ── UI construction ───────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        # ── Top: list + stacked page ───────────────────────────────
        top_row = QHBoxLayout()

        # Left: agent list
        self._list = QListWidget()
        self._list.setFixedWidth(200)
        self._list.setObjectName("AgentList")
        self._list.currentRowChanged.connect(self._on_agent_selected)
        top_row.addWidget(self._list)

        # Right: parameter pages
        self._stack = QStackedWidget()
        top_row.addWidget(self._stack)
        root.addLayout(top_row)

        # Populate list + stack
        self._populate_agents()

        # ── Geological context (LLM agents only) ──────────────────
        self._grp_context = QGroupBox("Geological Context  (optional — sent to LLM)")
        self._grp_context.setVisible(False)
        ctx_lay = QVBoxLayout(self._grp_context)
        ctx_lay.setContentsMargins(8, 6, 8, 6)
        ctx_hint = QLabel(
            "Describe the survey setting or pose a specific question. "
            "The LLM will use it alongside the analysis results."
        )
        ctx_hint.setWordWrap(True)
        ctx_lay.addWidget(ctx_hint)
        self._ctx_edit = QPlainTextEdit()
        self._ctx_edit.setPlaceholderText(
            "e.g. 'Semi-arid Precambrian terrain, looking for aquifer zones. "
            "What do low-resistivity anomalies below 200 m indicate?'"
        )
        self._ctx_edit.setMaximumHeight(72)
        ctx_lay.addWidget(self._ctx_edit)
        root.addWidget(self._grp_context)

        # ── API key row ────────────────────────────────────────────
        key_row = QHBoxLayout()
        key_row.addWidget(QLabel("API key (LLM agents):"))
        self._api_key_edit = QLineEdit(self.api_key)
        self._api_key_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self._api_key_edit.setPlaceholderText(
            "Leave blank to use ANTHROPIC_API_KEY / OPENAI_API_KEY env var"
        )
        key_row.addWidget(self._api_key_edit)
        root.addLayout(key_row)

        # ── OK / Cancel ────────────────────────────────────────────
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._on_accepted)
        buttons.rejected.connect(self.reject)
        root.addWidget(buttons)

    def _populate_agents(self) -> None:
        def _header(text: str) -> QListWidgetItem:
            item = QListWidgetItem(text)
            item.setFlags(Qt.ItemFlag.NoItemFlags)
            item.setForeground(Qt.GlobalColor.gray)
            return item

        # ── LLM agents grouped by category ────────────────────────────────
        self._list.addItem(_header("─── LLM Agents (API key required) ───"))
        for category, names in agents_by_category().items():
            self._list.addItem(_header(f"  {category}"))
            for name in names:
                entry = get_entry(name)
                item = QListWidgetItem(f"    {name}")
                item.setData(Qt.ItemDataRole.UserRole, name)
                item.setToolTip(entry.get("description", ""))
                self._list.addItem(item)
                page = _ParamPage(entry)
                self._pages[name] = page
                self._stack.addWidget(page)

        # ── Processing operations ──────────────────────────────────────────
        self._list.addItem(_header("─── Processing (no API key) ───"))
        for name in agent_names():
            entry = get_entry(name)
            if entry["type"] != "processing":
                continue
            item = QListWidgetItem(f"  {name}")
            item.setData(Qt.ItemDataRole.UserRole, name)
            item.setToolTip(entry.get("description", ""))
            self._list.addItem(item)
            page = _ParamPage(entry)
            self._pages[name] = page
            self._stack.addWidget(page)

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_agent_selected(self, row: int) -> None:
        item = self._list.item(row)
        if item is None:
            return
        name = item.data(Qt.ItemDataRole.UserRole)
        if not name:
            # header item — skip to next
            next_row = row + 1
            if next_row < self._list.count():
                self._list.setCurrentRow(next_row)
            return
        page = self._pages.get(name)
        if page:
            self._stack.setCurrentWidget(page)
        # Show context box only for LLM agents
        entry  = get_entry(name) or {}
        is_llm = entry.get("type") == "llm"
        self._grp_context.setVisible(is_llm)

    def _on_accepted(self) -> None:
        item = self._list.currentItem()
        if item is None:
            return
        name = item.data(Qt.ItemDataRole.UserRole)
        if not name:
            return
        self.selected_agent = name
        self.params         = self._pages[name].current_params()
        self.api_key        = self._api_key_edit.text().strip()
        ctx = self._ctx_edit.toPlainText().strip()
        if ctx:
            self.params["user_prompt"] = ctx
            self.params["context"]     = ctx   # compat alias
        self.accept()
