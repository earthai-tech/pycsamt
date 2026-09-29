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
    QCheckBox,
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
_LLM_PROVIDERS = ["claude", "openai", "gemini"]
_CLAUDE_MODELS = [
    "claude-opus-4-8",
    "claude-sonnet-4-6",
    "claude-haiku-4-5-20251001",
]
_OPENAI_MODELS = ["gpt-4o", "gpt-4o-mini", "gpt-3.5-turbo"]
_GEMINI_MODELS = ["gemini-2.0-flash", "gemini-1.5-pro", "gemini-1.5-flash"]
_MODELS_BY_PROVIDER = {
    "claude": _CLAUDE_MODELS,
    "openai": _OPENAI_MODELS,
    "gemini": _GEMINI_MODELS,
}
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
    Five-tab preferences dialog (General / Solvers / AI-LLM / License /
    Advanced).

    Pass the current ``SessionState``; after ``exec()`` returns ``Accepted``
    the session is already updated — the caller just needs to ``save()`` it.

    ``license_manager`` is forwarded to the License tab's
    ``LicensePage`` unchanged; omit it to fall back to an unconditional
    ``NullLicenseManager`` trial (tests, or any caller that deliberately
    wants no persisted license state).
    """

    def __init__(
        self,
        session: SessionState,
        parent: QWidget | None = None,
        license_manager=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Preferences")
        self.setMinimumSize(540, 380)
        self._session = session
        self._license_manager = license_manager
        self._build_ui()

    # ── Build ──────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        self._tabs = QTabWidget()
        root.addWidget(self._tabs)

        self._build_general_tab()
        self._build_solvers_tab()
        self._build_llm_tab()
        self._build_license_tab()
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

        self._anim_chk = QCheckBox("Animate station statistics")
        self._anim_chk.setToolTip(
            "Draw the response preview and frequency coverage with a short "
            "left-to-right animation when a station is selected")
        self._anim_chk.setChecked(bool(getattr(self._session,
                                               "ui_animations", True)))
        form.addRow("Motion:", self._anim_chk)

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

        # Restore the saved provider/model now that both combos exist.
        # setCurrentText on the provider combo fires _on_provider_changed
        # (rebuilding the model list for that provider) before the saved
        # model is selected below.
        saved_provider = self._session.llm_provider
        if saved_provider in _LLM_PROVIDERS:
            self._provider_combo.setCurrentText(saved_provider)
        if self._session.llm_model:
            idx = self._model_combo.findText(self._session.llm_model)
            if idx >= 0:
                self._model_combo.setCurrentIndex(idx)

        note = QLabel(
            "<i>API key is stored in ~/.pycsamt/session.json.<br>"
            "Alternatively set ANTHROPIC_API_KEY, OPENAI_API_KEY, or "
            "GEMINI_API_KEY / GOOGLE_API_KEY env var.</i>"
        )
        note.setWordWrap(True)
        form.addRow(note)

        self._tabs.addTab(w, "AI / LLM")

    def _on_provider_changed(self, provider: str) -> None:
        self._model_combo.clear()
        self._model_combo.addItems(
            _MODELS_BY_PROVIDER.get(provider, _CLAUDE_MODELS)
        )

    # ── Tab 3: License ────────────────────────────────────────────────

    def _build_license_tab(self) -> None:
        from pycsamt.app.desktop.widgets.license_page import LicensePage

        self._license_page = LicensePage(license_manager=self._license_manager)
        self._tabs.addTab(self._license_page, "License")

    # ── Tab 4: Advanced ───────────────────────────────────────────────

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
        s.ui_animations = self._anim_chk.isChecked()

        # Solvers
        s.occam2d_binary = self._occam2d_edit.text().strip()
        s.modem_binary = self._modem_edit.text().strip()
        s.mare2dem_binary = self._mare2dem_edit.text().strip()
        s.inversion_workdir = self._workdir_edit.text().strip()

        # AI / LLM
        s.api_key = self._api_key_edit.text().strip()
        s.llm_provider = self._provider_combo.currentText()
        s.llm_model = self._model_combo.currentText()

        # Advanced
        s.log_level = self._log_level_combo.currentText()
        s.tile_provider = self._tile_combo.currentText()

        # Apply log level immediately
        import logging

        logging.getLogger().setLevel(
            getattr(logging, s.log_level, logging.WARNING)
        )

        self.accept()
