# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PreferencesDialog — full application preferences (Phase 5).

Four tabs:
  General   — theme, default data directory, max recent files
  Solvers   — paths to Occam2D / ModEM / MARE2DEM executables + workdir
  AI / LLM  — API key (Anthropic / OpenAI), provider, model
  Advanced  — log level, map tile provider, cache directory

Changes are written back to the SessionState and applied immediately
(theme, log level).  The caller is responsible for calling session.save().
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.models.session import SessionState

_THEMES = ["dark", "light"]
_LOG_LEVELS = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
_LLM_PROVIDERS = ["claude", "openai"]
_CLAUDE_MODELS = [
    "claude-opus-4-8",
    "claude-sonnet-4-6",
    "claude-haiku-4-5-20251001",
]
_OPENAI_MODELS = ["gpt-4o", "gpt-4o-mini", "gpt-3.5-turbo"]
_TILE_PROVIDERS = [
    "OpenStreetMap.Mapnik",
    "Esri.WorldTopoMap",
    "Esri.WorldImagery",
    "CartoDB.Positron",
    "CartoDB.DarkMatter",
]


# ──────────────────────────────────────────────────────────────────────────────
# Helper: labelled path-browse row
# ──────────────────────────────────────────────────────────────────────────────


def _path_row(
    label: str, default: str, parent
) -> tuple[QHBoxLayout, QLineEdit]:
    """Return (layout, line_edit) for a path + Browse button."""
    row = QHBoxLayout()
    edit = QLineEdit(default)
    edit.setPlaceholderText(label)
    btn = QPushButton("Browse…")
    btn.setFixedWidth(70)

    def _browse() -> None:
        if "directory" in label.lower() or "workdir" in label.lower():
            p = QFileDialog.getExistingDirectory(
                parent, f"Select {label}", default
            )
        else:
            p, _ = QFileDialog.getOpenFileName(
                parent, f"Select {label}", default
            )
        if p:
            edit.setText(p)

    btn.clicked.connect(_browse)
    row.addWidget(edit)
    row.addWidget(btn)
    return row, edit


# ──────────────────────────────────────────────────────────────────────────────
# PreferencesDialog
# ──────────────────────────────────────────────────────────────────────────────


class PreferencesDialog(QDialog):
    """
    Four-tab preferences dialog.

    Pass the current ``SessionState``; after ``exec()`` returns ``Accepted``
    the session is already updated — the caller just needs to ``save()`` it.
    """

    def __init__(
        self,
        session: SessionState,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Preferences")
        self.setMinimumSize(540, 380)
        self._session = session
        self._build_ui()

    # ── Build ──────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        self._tabs = QTabWidget()
        root.addWidget(self._tabs)

        self._build_general_tab()
        self._build_solvers_tab()
        self._build_llm_tab()
        self._build_advanced_tab()

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._on_accepted)
        buttons.rejected.connect(self.reject)
        root.addWidget(buttons)

    # ── Tab 0: General ────────────────────────────────────────────────

    def _build_general_tab(self) -> None:
        w = QWidget()
        form = QFormLayout(w)
        form.setSpacing(10)

        self._theme_combo = QComboBox()
        self._theme_combo.addItems(_THEMES)
        self._theme_combo.setCurrentText(self._session.theme)
        form.addRow("Theme:", self._theme_combo)

        data_row, self._data_dir_edit = _path_row(
            "Default data directory", self._session.last_data_dir, self
        )
        form.addRow("Data directory:", data_row)

        self._max_recent_spin = QSpinBox()
        self._max_recent_spin.setRange(1, 100)
        self._max_recent_spin.setValue(self._session.max_recent_files)
        form.addRow("Max recent files:", self._max_recent_spin)

        self._tabs.addTab(w, "General")

    # ── Tab 1: Solvers ────────────────────────────────────────────────

    def _build_solvers_tab(self) -> None:
        w = QWidget()
        form = QFormLayout(w)
        form.setSpacing(10)

        oc_row, self._occam2d_edit = _path_row(
            "Occam2D binary", self._session.occam2d_binary, self
        )
        form.addRow("Occam2D binary:", oc_row)

        mo_row, self._modem_edit = _path_row(
            "ModEM binary", self._session.modem_binary, self
        )
        form.addRow("ModEM binary:", mo_row)

        mr_row, self._mare2dem_edit = _path_row(
            "MARE2DEM binary", self._session.mare2dem_binary, self
        )
        form.addRow("MARE2DEM binary:", mr_row)

        wd_row, self._workdir_edit = _path_row(
            "Inversion working directory",
            self._session.inversion_workdir,
            self,
        )
        form.addRow("Default workdir:", wd_row)

        self._tabs.addTab(w, "Solvers")

    # ── Tab 2: AI / LLM ──────────────────────────────────────────────

    def _build_llm_tab(self) -> None:
        w = QWidget()
        form = QFormLayout(w)
        form.setSpacing(10)

        self._api_key_edit = QLineEdit(self._session.api_key)
        self._api_key_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self._api_key_edit.setPlaceholderText(
            "sk-ant-… or sk-… (leave blank to use env var)"
        )
        form.addRow("API key:", self._api_key_edit)

        self._provider_combo = QComboBox()
        self._provider_combo.addItems(_LLM_PROVIDERS)
        self._provider_combo.currentTextChanged.connect(
            self._on_provider_changed
        )
        form.addRow("Provider:", self._provider_combo)

        self._model_combo = QComboBox()
        self._model_combo.addItems(_CLAUDE_MODELS)
        form.addRow("Model:", self._model_combo)

        note = QLabel(
            "<i>API key is stored in ~/.pycsamt/session.json.<br>"
            "Alternatively set ANTHROPIC_API_KEY or OPENAI_API_KEY env var.</i>"
        )
        note.setWordWrap(True)
        form.addRow(note)

        self._tabs.addTab(w, "AI / LLM")

    def _on_provider_changed(self, provider: str) -> None:
        self._model_combo.clear()
        if provider == "claude":
            self._model_combo.addItems(_CLAUDE_MODELS)
        else:
            self._model_combo.addItems(_OPENAI_MODELS)

    # ── Tab 3: Advanced ───────────────────────────────────────────────

    def _build_advanced_tab(self) -> None:
        w = QWidget()
        form = QFormLayout(w)
        form.setSpacing(10)

        self._log_level_combo = QComboBox()
        self._log_level_combo.addItems(_LOG_LEVELS)
        self._log_level_combo.setCurrentText(self._session.log_level)
        form.addRow("Log level:", self._log_level_combo)

        self._tile_combo = QComboBox()
        self._tile_combo.addItems(_TILE_PROVIDERS)
        self._tile_combo.setCurrentText(self._session.tile_provider)
        form.addRow("Map tile provider:", self._tile_combo)

        cache_row, self._cache_edit = _path_row(
            "Tile cache directory",
            str(Path.home() / ".pycsamt" / "tile_cache"),
            self,
        )
        form.addRow("Tile cache dir:", cache_row)

        self._tabs.addTab(w, "Advanced")

    # ── Accept ────────────────────────────────────────────────────────

    def _on_accepted(self) -> None:
        s = self._session

        # General
        s.theme = self._theme_combo.currentText()
        s.last_data_dir = self._data_dir_edit.text().strip()
        s.max_recent_files = self._max_recent_spin.value()

        # Solvers
        s.occam2d_binary = self._occam2d_edit.text().strip()
        s.modem_binary = self._modem_edit.text().strip()
        s.mare2dem_binary = self._mare2dem_edit.text().strip()
        s.inversion_workdir = self._workdir_edit.text().strip()

        # AI / LLM
        s.api_key = self._api_key_edit.text().strip()

        # Advanced
        s.log_level = self._log_level_combo.currentText()
        s.tile_provider = self._tile_combo.currentText()

        # Apply log level immediately
        import logging

        logging.getLogger().setLevel(
            getattr(logging, s.log_level, logging.WARNING)
        )

        self.accept()
