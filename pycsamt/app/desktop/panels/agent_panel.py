# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AgentPanel — right dock widget for AI workflow and agent output.

Phase 4: Contains:
  * Run Agent… button → opens RunAgentDialog
  * Export Figure… button → opens ExportDialog on the current result canvas
  * Progress bar + Stop button
  * Three-tab output viewer:
      Log      — streaming lines from AgentWorker.log_line
      Result   — MplCanvas showing the agent result figure
      Summary  — formatted AgentResult summary and LLM interpretation
"""

from __future__ import annotations

import datetime
from typing import Optional

from PySide6.QtCore import Qt, Signal, Slot
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QTabWidget,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas


class AgentPanel(QWidget):
    """
    Right-dock panel for running pycsamt agents and viewing results.

    Parameters
    ----------
    app_controller : AppController
        Shared application state — used to access the loaded Sites.
    parent : QWidget, optional
    """

    # Emitted just before the worker thread starts so the dock can show itself
    agent_started = Signal(str)   # agent_name

    def __init__(self, app_controller=None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._ctrl   = app_controller
        self._worker: Optional[object] = None   # AgentWorker (set on run)
        self._result_figure = None              # last matplotlib Figure
        self._build_ui()

    # ── Construction ──────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(4, 4, 4, 4)
        root.setSpacing(4)

        # ── Button bar ─────────────────────────────────────────────
        btn_row = QHBoxLayout()

        self._btn_run = QPushButton("▶  Run Agent…")
        self._btn_run.setObjectName("RunAgentButton")
        self._btn_run.setStatusTip("Open the agent picker dialog")
        self._btn_run.clicked.connect(self._on_run)
        btn_row.addWidget(self._btn_run)

        self._btn_export = QPushButton("⬆  Export…")
        self._btn_export.setObjectName("ExportButton")
        self._btn_export.setStatusTip("Export the current result figure")
        self._btn_export.setEnabled(False)
        self._btn_export.clicked.connect(self._on_export)
        btn_row.addWidget(self._btn_export)

        self._btn_stop = QPushButton("■  Stop")
        self._btn_stop.setObjectName("StopButton")
        self._btn_stop.setEnabled(False)
        self._btn_stop.clicked.connect(self._on_stop)
        btn_row.addWidget(self._btn_stop)

        root.addLayout(btn_row)

        # ── Progress bar ───────────────────────────────────────────
        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setMaximumHeight(12)
        self._progress.setVisible(False)
        self._progress.setTextVisible(False)
        root.addWidget(self._progress)

        # ── Status label ───────────────────────────────────────────
        self._status_lbl = QLabel("No agent running.")
        self._status_lbl.setObjectName("AgentStatusLabel")
        root.addWidget(self._status_lbl)

        # ── Output tabs ────────────────────────────────────────────
        self._tabs = QTabWidget()
        self._tabs.setDocumentMode(True)
        root.addWidget(self._tabs)

        # Tab 0: Log
        self._log_text = QPlainTextEdit()
        self._log_text.setReadOnly(True)
        self._log_text.setMaximumBlockCount(3000)
        self._log_text.setPlaceholderText("Agent log output will appear here…")
        self._tabs.addTab(self._log_text, "Log")

        # Tab 1: Result figure
        self._result_canvas = MplCanvas(self, toolbar=True)
        self._tabs.addTab(self._result_canvas, "Result")

        # Tab 2: Summary / interpretation
        self._summary_browser = QTextBrowser()
        self._summary_browser.setPlaceholderText(
            "Agent summary and LLM interpretation will appear here…"
        )
        self._tabs.addTab(self._summary_browser, "Summary")

    # ── Public API (called from MainWindow / AppController) ───────────

    def set_app_controller(self, ctrl) -> None:
        self._ctrl = ctrl

    def append_log(self, text: str) -> None:
        """Append a timestamped line to the Log tab."""
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self._log_text.appendPlainText(f"[{ts}]  {text}")

    # ── Run / Stop ────────────────────────────────────────────────────

    def _on_run(self) -> None:
        from pycsamt.app.desktop.dialogs.run_agent_dlg import RunAgentDialog

        if self._ctrl is None or self._ctrl.sites is None:
            self._status_lbl.setText("Load survey data before running an agent.")
            return

        # Get API key from session if available
        api_key: Optional[str] = None
        if self._ctrl is not None and hasattr(self._ctrl.session, "api_key"):
            api_key = self._ctrl.session.api_key or None

        dlg = RunAgentDialog(api_key=api_key, parent=self)
        if dlg.exec() != RunAgentDialog.DialogCode.Accepted:
            return

        self._start_worker(
            agent_name=dlg.selected_agent,
            params=dlg.params,
            api_key=dlg.api_key,
        )

    def _start_worker(
        self,
        agent_name: str,
        params: dict,
        api_key: Optional[str] = None,
    ) -> None:
        from pycsamt.app.desktop.workers.agent_worker import AgentWorker

        self.agent_started.emit(agent_name)   # let the dock show itself first

        self._log_text.clear()
        self._summary_browser.clear()
        self._progress.setValue(0)
        self._progress.setVisible(True)
        self._btn_stop.setEnabled(True)
        self._btn_run.setEnabled(False)
        self._status_lbl.setText(f"Running: {agent_name}…")
        self._tabs.setCurrentIndex(0)   # show Log tab

        self._worker = AgentWorker(
            agent_name=agent_name,
            sites=self._ctrl.sites,
            params=params,
            api_key=api_key,
            parent=self,
        )
        self._worker.log_line.connect(self._on_log_line)
        self._worker.progress.connect(self._progress.setValue)
        self._worker.result.connect(self._on_result)
        self._worker.error.connect(self._on_error)
        self._worker.finished.connect(self._on_worker_finished)
        self._worker.start()

    def _on_stop(self) -> None:
        if self._worker is not None:
            self._worker.cancel()
            self._worker.quit()
            self._status_lbl.setText("Stopped.")
            self._on_worker_finished()

    # ── Worker signals ────────────────────────────────────────────────

    @Slot(str)
    def _on_log_line(self, text: str) -> None:
        self.append_log(text)

    @Slot(object)
    def _on_result(self, result) -> None:
        """Handle AgentResult (LLM) or raw function return value."""
        self._render_result(result)
        self._btn_export.setEnabled(self._result_figure is not None)

    @Slot(str)
    def _on_error(self, message: str) -> None:
        self.append_log(f"ERROR: {message}")
        self._status_lbl.setText(f"Failed: {message[:60]}")

    def _on_worker_finished(self) -> None:
        self._progress.setVisible(False)
        self._btn_stop.setEnabled(False)
        self._btn_run.setEnabled(True)
        if "Running" in self._status_lbl.text():
            self._status_lbl.setText("Completed.")

    # ── Result rendering ──────────────────────────────────────────────

    def _render_result(self, result) -> None:
        """Display result in the Result and Summary tabs."""
        renderable = _extract_renderable(result)
        if renderable is not None:
            self._result_figure = (
                renderable.get_figure()
                if hasattr(renderable, "get_figure")
                else renderable
            )
            try:
                self._result_canvas.show_figure(renderable)
            except Exception:
                pass
            self._tabs.setCurrentIndex(1)

        self._render_summary(result)

    def _render_summary(self, result) -> None:
        """Format and display summary / LLM interpretation."""
        lines = []

        if hasattr(result, "status"):
            lines.append(f"<b>Status:</b> {result.status}")
        if hasattr(result, "summary") and result.summary:
            lines.append(f"<b>Summary:</b> {result.summary}")
        if hasattr(result, "elapsed_seconds"):
            lines.append(f"<b>Elapsed:</b> {result.elapsed_seconds:.1f} s")
        if hasattr(result, "cost_estimate") and result.cost_estimate:
            lines.append(f"<b>Cost estimate:</b> {result.cost_estimate}")
        if hasattr(result, "llm_interpretation") and result.llm_interpretation:
            lines.append("<hr>")
            lines.append("<b>LLM Interpretation:</b>")
            lines.append(result.llm_interpretation.replace("\n", "<br>"))
        if hasattr(result, "warnings") and result.warnings:
            lines.append("<hr>")
            lines.append("<b>Warnings:</b>")
            for w in result.warnings:
                lines.append(f"• {w}")

        if not lines:
            lines.append("<i>Processing result — see Log tab for details.</i>")

        self._summary_browser.setHtml("<br>".join(lines))
        if lines:
            self._tabs.setCurrentIndex(2)

    # ── Export ────────────────────────────────────────────────────────

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog
        fig = self._result_figure or self._result_canvas.figure
        if fig is None:
            return
        dlg = ExportDialog(figure=fig, parent=self)
        dlg.exec()


# ── Module-level helper (shared with agent_window) ────────────────────────────

def _is_renderable(obj) -> bool:
    """Return True if *obj* is a matplotlib Figure or Axes."""
    return hasattr(obj, "savefig") or hasattr(obj, "get_figure")


def _extract_renderable(result):
    """Return the first Figure or Axes found in *result*, or None."""
    figs = _extract_all_renderables(result)
    return figs[0][1] if figs else None


def _extract_all_renderables(result) -> list[tuple[str, object]]:
    """Return all (label, Figure-or-Axes) pairs from *result*, in order.

    Handles:
      • Direct Figure / Axes
      • AgentResult.data[key]  → Figure/Axes
      • AgentResult.data["figures"][key] → Figure/Axes  (nested dict)
      • AgentResult.figure attribute
    """
    found: list[tuple[str, object]] = []

    if _is_renderable(result):
        return [("Figure", result)]

    if hasattr(result, "data") and isinstance(result.data, dict):
        for key, val in result.data.items():
            if _is_renderable(val):
                found.append((_label(key), val))
            elif isinstance(val, dict):
                # one level of nesting — e.g. data["figures"] = {"ss_summary": fig, …}
                for inner_key, inner_val in val.items():
                    if _is_renderable(inner_val):
                        found.append((_label(inner_key), inner_val))

    fig = getattr(result, "figure", None)
    if fig is not None and _is_renderable(fig):
        if not any(v is fig for _, v in found):
            found.append(("Figure", fig))

    return found


def _label(key: str) -> str:
    """Convert a snake_case key to a Title Case display label."""
    return key.replace("_", " ").title()
