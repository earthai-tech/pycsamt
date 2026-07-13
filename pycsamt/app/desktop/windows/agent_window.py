# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AgentRunnerWindow — floating window for AI / processing agents.

Layout
------
  ┌──────────────────────────┬────────────────────────────────────────────┐
  │  BROWSER (top-left)      │  CONTENT (right)                           │
  │  🔍 Search agents…       │  ┌─ 📋 Log ─ 📊 Result ─ 📝 Summary ───┐  │
  │  [All][QC][Pre-proc][…]  │  │                                       │  │
  │  ┌────────────────────┐  │  │  (agent output / figures)             │  │
  │  │ ▌ QC Analysis [LLM]│  │  └───────────────────────────────────────┘  │
  │  │ ▌ Anomaly Det [LLM]│  │═══════════ splitter ════════════════════════│
  │  │ ▌ QC Quicklook     │  │  💬 Chat with Agent          [🗑  Clear]    │
  │  └────────────────────┘  │  ┌───────────────────────────────────────┐  │
  ├──────────────────────────┤  │  welcome / chat history (QTextBrowser) │  │
  │  DETAIL (bottom-left)    │  └───────────────────────────────────────┘  │
  │  ── QC Analysis ──────   │  [Ask a follow-up question…         ] [▶]   │
  │  Geological Context…     │   Enter to send · Shift+Enter for new line   │
  │  [sk-… 🔒]               │                                             │
  │  [▶ Run Agent] [■ Stop]  │                                             │
  └──────────────────────────┴─────────────────────────────────────────────┘

Left column : QSplitter(Vertical) — browser / detail proportions resizable.
Right column: QSplitter(Vertical) — tabs / chat proportions resizable.
"""

from __future__ import annotations

import datetime

from PySide6.QtCore import QEvent, Qt, Slot
from PySide6.QtGui import QKeyEvent, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QTabWidget,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.agent_browser import (
    AgentBrowserWidget,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    _icon,
)

# ─────────────────────────────────────────────────────────────────────────────
# Figure pop-out: hover-reveal overlay button + publication-ready viewer
# ─────────────────────────────────────────────────────────────────────────────


class _HoverRevealButton(QPushButton):
    """
    Transparent expand button overlaid on a result container widget.

    Invisible by default; fades in when the mouse enters the container
    (including any of its child widgets).  Click opens FigureViewerWindow.
    """

    _SS_HIDDEN = """QPushButton {
        background: transparent;
        border: none;
        color:       transparent;
    }"""
    _SS_SHOWN = """QPushButton {
        background:   rgba(18, 20, 32, 0.72);
        border:       1px solid rgba(255, 255, 255, 0.22);
        border-radius: 7px;
        color:        rgba(255, 255, 255, 0.90);
        font-size:    14px;
    }
    QPushButton:hover {
        background:   rgba(55, 60, 90, 0.90);
        border-color: rgba(255, 255, 255, 0.45);
        color:        white;
    }
    QPushButton:pressed {
        background:   rgba(80, 86, 130, 0.95);
    }"""

    def __init__(
        self,
        container: QWidget,
        get_canvas_fn,  # () → MplCanvas | None
        parent=None,
    ) -> None:
        super().__init__("⤢", parent or container)
        self._container = container
        self._get_canvas = get_canvas_fn
        self._shown = False

        self.setFixedSize(30, 30)
        self.setToolTip("Open in full-size viewer (publication-ready)")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setStyleSheet(self._SS_HIDDEN)

        self.clicked.connect(self._on_click)
        container.installEventFilter(self)
        self._reposition()
        self.raise_()

    # ── Event filter on parent container ──────────────────────────────────────

    def eventFilter(self, obj, event) -> bool:
        if obj is self._container:
            t = event.type()
            if t == QEvent.Type.Resize:
                self._reposition()
                self.raise_()
            elif t == QEvent.Type.Enter:
                self._set_shown(True)
            elif t == QEvent.Type.Leave:
                # Leave only fires when mouse exits the whole container tree
                self._set_shown(False)
        return super().eventFilter(obj, event)

    def _set_shown(self, shown: bool) -> None:
        if shown == self._shown:
            return
        self._shown = shown
        self.setStyleSheet(self._SS_SHOWN if shown else self._SS_HIDDEN)
        if shown:
            self.raise_()

    def _reposition(self) -> None:
        m = 8
        self.move(self._container.width() - self.width() - m, m)

    # ── Click → open viewer ───────────────────────────────────────────────────

    def _on_click(self) -> None:
        canvas = self._get_canvas()
        if canvas is None:
            return
        fig = getattr(canvas, "figure", None)
        if fig is None:
            return
        # Label = current inner-tab text (e.g. "Ss Summary")
        inner_tab = self._container.findChild(QTabWidget)
        label = (
            inner_tab.tabText(inner_tab.currentIndex())
            if inner_tab
            else "Figure Viewer"
        )
        viewer = FigureViewerWindow(
            figure=fig,
            title=label,
            parent=self._container.window(),
        )
        viewer.show()


class FigureViewerWindow(QDialog):
    """
    Large, scrollable, publication-ready figure viewer.

    Renders the matplotlib Figure to an in-memory PNG pixmap via
    ``fig.savefig()`` — which does NOT mutate the figure's size.
    This means:
      • Re-opening the viewer always starts from the original figure size.
      • Closing the viewer leaves the main Result canvas unchanged.
      • A QScrollArea lets the user pan around the full-resolution image.

    Controls
    --------
    DPI combo   — choose render resolution (72 / 100 / 150 / 200 / 300)
    Ctrl+scroll — zoom in / out (25 DPI steps, range 50–400)
    Export …    — save PNG / PDF / SVG at the currently selected DPI
    """

    _DPIS = (72, 100, 150, 200, 300)

    def __init__(
        self,
        figure,
        title: str = "Figure Viewer",
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._figure = figure
        self._dpi = 150

        self.setWindowTitle(title)
        self.setWindowFlags(
            Qt.WindowType.Window
            | Qt.WindowType.WindowMinMaxButtonsHint
            | Qt.WindowType.WindowCloseButtonHint
        )
        self.resize(1280, 920)
        self.setMinimumSize(800, 600)

        vlay = QVBoxLayout(self)
        vlay.setContentsMargins(0, 0, 0, 0)
        vlay.setSpacing(0)

        # ── Scrollable image area ─────────────────────────────────────────
        self._scroll = QScrollArea()
        self._scroll.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._scroll.setWidgetResizable(False)  # widget drives size, not area
        self._scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        self._scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        self._img_label = QLabel()
        self._img_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._scroll.setWidget(self._img_label)
        vlay.addWidget(self._scroll)

        # ── Bottom action bar ─────────────────────────────────────────────
        bar = QWidget()
        bar.setObjectName("SummaryToolbar")
        bar_lay = QHBoxLayout(bar)
        bar_lay.setContentsMargins(8, 4, 8, 4)
        bar_lay.setSpacing(8)

        bar_lay.addWidget(QLabel("DPI:"))
        self._dpi_combo = QComboBox()
        for d in self._DPIS:
            self._dpi_combo.addItem(str(d))
        self._dpi_combo.setCurrentText("150")
        self._dpi_combo.currentTextChanged.connect(self._on_dpi_changed)
        bar_lay.addWidget(self._dpi_combo)

        _hint = QLabel("Ctrl+scroll to zoom  ·  scrollbars to pan")
        _hint.setObjectName("InfoLabel")
        bar_lay.addWidget(_hint)
        bar_lay.addStretch()

        for fmt, lbl in (
            ("png", "Export PNG…"),
            ("pdf", "Export PDF…"),
            ("svg", "Export SVG…"),
        ):
            btn = QPushButton(lbl)
            btn.clicked.connect(lambda _=False, f=fmt: self._export(f))
            bar_lay.addWidget(btn)

        vlay.addWidget(bar)

        # Initial render
        self._render()

    # ── Rendering ─────────────────────────────────────────────────────────

    def _render(self) -> None:
        """Rasterise figure at self._dpi → pixmap in QLabel.  Figure unchanged."""
        import io

        buf = io.BytesIO()
        try:
            self._figure.savefig(
                buf, format="png", dpi=self._dpi, bbox_inches="tight"
            )
        except Exception as exc:
            self._img_label.setText(f"Render error:\n{exc}")
            return
        buf.seek(0)
        px = QPixmap()
        px.loadFromData(buf.getvalue())
        self._img_label.setPixmap(px)
        self._img_label.resize(px.size())  # let scrollbars appear as needed

    def _on_dpi_changed(self, text: str) -> None:
        try:
            self._dpi = int(text)
        except ValueError:
            return
        self._render()

    # ── Ctrl+scroll zoom ──────────────────────────────────────────────────

    def wheelEvent(self, event) -> None:
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            delta = event.angleDelta().y()
            step = 25 if delta > 0 else -25
            new_dpi = max(50, min(400, self._dpi + step))
            if new_dpi != self._dpi:
                self._dpi = new_dpi
                # Sync combo without triggering a second _render call
                self._dpi_combo.blockSignals(True)
                idx = self._dpi_combo.findText(str(new_dpi))
                if idx >= 0:
                    self._dpi_combo.setCurrentIndex(idx)
                self._dpi_combo.blockSignals(False)
                self._render()
        else:
            super().wheelEvent(event)

    # ── Export ────────────────────────────────────────────────────────────

    def _export(self, fmt: str) -> None:
        _filters = {
            "png": "PNG Image (*.png)",
            "pdf": "PDF Document (*.pdf)",
            "svg": "SVG Vector (*.svg)",
        }
        path, _ = QFileDialog.getSaveFileName(
            self,
            f"Export figure as {fmt.upper()}",
            f"figure.{fmt}",
            _filters[fmt],
        )
        if not path:
            return
        try:
            self._figure.savefig(path, dpi=self._dpi, bbox_inches="tight")
        except Exception as exc:
            QMessageBox.warning(self, "Export Failed", str(exc))


# ─────────────────────────────────────────────────────────────────────────────
# Module-level constants / helpers
# ─────────────────────────────────────────────────────────────────────────────

_CHAT_WELCOME_HTML = (
    "<html><body style='font-family:system-ui,sans-serif;"
    "font-size:13px;padding:24px;'>"
    "<p style='color:#555;text-align:center;font-size:22px;margin:0'>💬</p>"
    "<p style='color:#888;text-align:center;font-size:12px;margin:8px 0 0'>"
    "<b>Chat with the agent</b></p>"
    "<p style='color:#666;text-align:center;font-size:11px;margin:4px 0 0'>"
    "Run an LLM agent, then ask follow-up questions here."
    "</p></body></html>"
)


def _format_result_html(result) -> str:
    """Build structured HTML for the Summary tab from an AgentResult."""

    def _esc(s: str) -> str:
        return (
            str(s)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
        )

    # Figure-only result (processing agents returning Axes / Figure)
    if hasattr(result, "savefig") or hasattr(result, "get_figure"):
        return (
            "<html><body style='font-family:system-ui,sans-serif;"
            "font-size:13px;padding:24px;'>"
            "<p style='color:#888;text-align:center;font-size:20px'>📊</p>"
            "<p style='color:#888;text-align:center;font-size:12px;'>"
            "Figure result — see the <b>Result</b> tab for the plot.</p>"
            "</body></html>"
        )

    sections: list[str] = []

    # ── Status badge + metadata row ────────────────────────────────────────
    status = getattr(result, "status", None)
    elapsed = getattr(result, "elapsed_seconds", None)
    cost = getattr(result, "cost_estimate_usd", None)
    if status or elapsed is not None or cost:
        if status == "success":
            badge_bg = "#a6e3a1"
        elif status in ("error", "failed"):
            badge_bg = "#f38ba8"
        else:
            badge_bg = "#6c7086"
        meta_parts: list[str] = []
        if elapsed is not None:
            meta_parts.append(f"{elapsed:.1f} s")
        if cost:
            meta_parts.append(f"${cost:.4f}")
        meta_str = "  ·  ".join(meta_parts)
        sections.append(
            "<table width='100%' cellpadding='0' cellspacing='0'><tr>"
            f"<td><span style='background:{badge_bg};color:#1e1e2e;"
            f"padding:2px 10px;border-radius:4px;font-weight:bold;"
            f"font-size:11px'>&nbsp;{_esc(str(status or '?').upper())}"
            f"&nbsp;</span></td>"
            f"<td align='right' style='color:#888;font-size:11px'>"
            f"{meta_str}</td></tr></table>"
        )

    # ── Overview / summary text ────────────────────────────────────────────
    summary = getattr(result, "summary", None)
    if summary:
        sections.append(
            "<p style='margin:12px 0 4px;font-weight:bold;"
            "color:#cdd6f4;font-size:12px;'>Overview</p>"
            "<p style='margin:0 0 10px;font-size:12px;line-height:1.5;"
            f"color:#cdd6f4;'>{_esc(summary)}</p>"
        )

    # ── LLM Interpretation (primary content, left-border accent) ──────────
    interp = getattr(result, "llm_interpretation", None)
    if interp:
        safe_interp = _esc(interp).replace("\n", "<br/>")
        sections.append(
            "<div style='border-left:3px solid #89b4fa;padding:8px 14px;"
            "margin:10px 0;background:#1a2535;border-radius:0 6px 6px 0;'>"
            "<div style='font-weight:bold;color:#89b4fa;margin-bottom:8px;"
            "font-size:12px;'>🤖&nbsp; LLM Interpretation</div>"
            "<div style='font-size:12px;line-height:1.65;color:#cdd6f4'>"
            f"{safe_interp}</div></div>"
        )

    # ── Warnings ──────────────────────────────────────────────────────────
    warnings = getattr(result, "warnings", None)
    if warnings:
        rows = "".join(
            f"<div style='color:#e0af68;font-size:11px;margin:3px 0'>"
            f"• {_esc(w)}</div>"
            for w in warnings
        )
        sections.append(
            "<div style='background:#2a1500;border:1px solid #fab387;"
            "border-radius:6px;padding:10px;margin:10px 0;'>"
            "<div style='color:#fab387;font-weight:bold;margin-bottom:6px;"
            f"font-size:11px;'>⚠&nbsp; Warnings ({len(warnings)})</div>"
            + rows
            + "</div>"
        )

    if not sections:
        sections.append(
            f"<pre style='font-size:11px;color:#888;"
            f"white-space:pre-wrap;'>{_esc(str(result))}</pre>"
        )

    body = "\n".join(sections)
    return (
        "<html><body style='font-family:system-ui,sans-serif;"
        "font-size:13px;padding:12px;'>" + body + "</body></html>"
    )


# ═════════════════════════════════════════════════════════════════════════════
# AgentRunnerWindow
# ═════════════════════════════════════════════════════════════════════════════


class AgentRunnerWindow(PanelWindow):
    """Floating Agent Runner with category-card browser on the left."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(
            title="Agent Runner",
            session_key="agent_runner",
            params_width=320,
            icon_name="agents",
            parent=parent,
        )
        self.resize(1400, 840)
        self._worker = None
        self._ctrl = None
        self._current_agent: str = ""
        self._params_widgets: dict[str, QWidget] = {}
        self._last_agent = None  # last-run LLM agent instance
        self._last_result_context = ""  # text context for follow-up chat
        self._last_interpretation = ""  # raw LLM interpretation (for copy)
        self._chat_worker = None  # ChatWorker for follow-up queries

    # ── Shell override ────────────────────────────────────────────────────────

    def _build_shell(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self._splitter = QSplitter(Qt.Orientation.Horizontal, self)
        self._splitter.setHandleWidth(4)

        # ── Left: vertical splitter (browser | detail) ───────────────────────
        self._left_vsplit = QSplitter(Qt.Orientation.Vertical)
        self._left_vsplit.setHandleWidth(3)
        self._left_vsplit.setObjectName("ParamsPanel")
        self._left_vsplit.setMinimumWidth(220)
        self._left_vsplit.setMaximumWidth(480)

        self._browser = AgentBrowserWidget()
        self._browser.agent_selected.connect(self._on_agent_selected)
        self._browser.agent_activated.connect(self._on_agent_activated)
        self._left_vsplit.addWidget(self._browser)

        detail_outer = QWidget()
        detail_outer.setObjectName("ParamsPanel")
        detail_vlay = QVBoxLayout(detail_outer)
        detail_vlay.setContentsMargins(0, 0, 0, 0)
        detail_vlay.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)

        inner = QWidget()
        inner.setObjectName("ParamsInner")
        self._params_layout = QVBoxLayout(inner)
        self._params_layout.setContentsMargins(8, 6, 8, 8)
        self._params_layout.setSpacing(6)
        self._build_params(self._params_layout)
        self._params_layout.addStretch(1)

        scroll.setWidget(inner)
        detail_vlay.addWidget(scroll)
        self._left_vsplit.addWidget(detail_outer)
        self._left_vsplit.setSizes([440, 300])

        # ── Right: content ───────────────────────────────────────────────────
        self._content_widget = QWidget()
        self._content_widget.setObjectName("ContentPanel")
        self._content_layout = QVBoxLayout(self._content_widget)
        self._content_layout.setContentsMargins(0, 0, 0, 0)
        self._content_layout.setSpacing(0)
        self._build_content(self._content_layout)

        self._splitter.addWidget(self._left_vsplit)
        self._splitter.addWidget(self._content_widget)
        self._splitter.setStretchFactor(0, 0)
        self._splitter.setStretchFactor(1, 1)
        self._splitter.setSizes([260, 900])

        outer.addWidget(self._splitter)

    # ── Detail panel (bottom-left) ────────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:
        self._detail_header = QLabel("No agent selected")
        self._detail_header.setObjectName("DetailHeader")
        self._detail_header.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._detail_header.setStyleSheet(
            "font-weight: bold; font-size: 12px; padding: 4px 0;"
        )
        layout.addWidget(self._detail_header)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setObjectName("Separator")
        layout.addWidget(sep)

        # Dynamic params form
        self._grp_params = QGroupBox("Parameters")
        self._grp_params.setObjectName("ParamsGroup")
        self._form_params = QFormLayout(self._grp_params)
        self._form_params.setSpacing(5)
        self._form_params.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        self._grp_params.setVisible(False)
        layout.addWidget(self._grp_params)

        # Geological context (LLM only)
        self._grp_context = QGroupBox("Geological Context")
        self._grp_context.setObjectName("ParamsGroup")
        ctx_lay = QVBoxLayout(self._grp_context)
        ctx_lay.setContentsMargins(6, 8, 6, 6)
        ctx_lay.setSpacing(4)
        ctx_hint = QLabel(
            "Optional: describe the survey setting or ask a specific question. "
            "Sent to the LLM alongside the analysis results."
        )
        ctx_hint.setWordWrap(True)
        ctx_hint.setObjectName("InfoLabel")
        ctx_lay.addWidget(ctx_hint)
        self._ctx_edit = QPlainTextEdit()
        self._ctx_edit.setPlaceholderText(
            "e.g. 'Semi-arid Precambrian terrain, looking for aquifer zones. "
            "What do low-resistivity anomalies below 200 m indicate?'"
        )
        self._ctx_edit.setMaximumHeight(80)
        ctx_lay.addWidget(self._ctx_edit)
        self._grp_context.setVisible(False)
        layout.addWidget(self._grp_context)

        # API Key (LLM only)
        self._grp_apikey = QGroupBox("API Key")
        self._grp_apikey.setObjectName("ParamsGroup")
        api_lay = QVBoxLayout(self._grp_apikey)
        api_lay.setContentsMargins(6, 8, 6, 6)
        self._edit_apikey = QLineEdit()
        self._edit_apikey.setPlaceholderText(
            "sk-…  (leave blank for env var)"
        )
        self._edit_apikey.setEchoMode(QLineEdit.EchoMode.Password)
        api_lay.addWidget(self._edit_apikey)
        self._grp_apikey.setVisible(False)
        layout.addWidget(self._grp_apikey)

        # Run / Stop + progress
        grp_run = QGroupBox("Run")
        grp_run.setObjectName("ParamsGroup")
        run_lay = QVBoxLayout(grp_run)
        run_lay.setContentsMargins(6, 8, 6, 6)
        run_lay.setSpacing(5)

        btn_row = QHBoxLayout()
        self._btn_run = QPushButton("▶  Run Agent")
        self._btn_run.setIcon(_icon("agents"))
        self._btn_run.setEnabled(False)
        self._btn_stop = QPushButton("■  Stop")
        self._btn_stop.setEnabled(False)
        self._btn_run.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._btn_stop.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._btn_run.clicked.connect(self._on_run)
        self._btn_stop.clicked.connect(self._on_stop)
        btn_row.addWidget(self._btn_run)
        btn_row.addWidget(self._btn_stop)
        run_lay.addLayout(btn_row)

        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setMaximumHeight(10)
        self._progress.setTextVisible(False)
        self._progress.setVisible(False)
        run_lay.addWidget(self._progress)

        self._status_lbl = QLabel("Select an agent above, then click Run.")
        self._status_lbl.setWordWrap(True)
        self._status_lbl.setObjectName("InfoLabel")
        run_lay.addWidget(self._status_lbl)

        layout.addWidget(grp_run)

    # ── Right panel: tabs (top) + persistent chat (bottom) ───────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        # Vertical splitter lets the user resize the result vs. chat areas
        self._content_vsplit = QSplitter(Qt.Orientation.Vertical)
        self._content_vsplit.setHandleWidth(4)

        # ── Top: tab widget (Log | Result | Summary) ────────────────────────
        self._tabs = QTabWidget()
        self._tabs.setDocumentMode(True)

        self._log_text = QPlainTextEdit()
        self._log_text.setReadOnly(True)
        self._log_text.setMaximumBlockCount(5000)
        self._log_text.setPlaceholderText("Agent output will appear here…")
        self._tabs.addTab(self._log_text, _icon("log"), "Log")

        # Result tab: inner QTabWidget so multiple figures each get their own tab
        _result_container = QWidget()
        _result_vlay = QVBoxLayout(_result_container)
        _result_vlay.setContentsMargins(0, 0, 0, 0)
        _result_vlay.setSpacing(0)
        self._result_inner_tabs = QTabWidget()
        self._result_inner_tabs.setDocumentMode(True)
        self._result_inner_tabs.setTabBarAutoHide(
            True
        )  # hidden when only 1 figure
        # Seed with a placeholder canvas so the widget is never empty
        self._result_canvas = MplCanvas(self, toolbar=True)
        self._result_inner_tabs.addTab(self._result_canvas, "Result")
        _result_vlay.addWidget(self._result_inner_tabs)
        self._tabs.addTab(_result_container, _icon("results"), "Result")

        # Summary tab = thin copy-button toolbar + browser
        _sum_widget = QWidget()
        _sum_vlay = QVBoxLayout(_sum_widget)
        _sum_vlay.setContentsMargins(0, 0, 0, 0)
        _sum_vlay.setSpacing(0)

        _sum_toolbar = QWidget()
        _sum_toolbar.setObjectName("SummaryToolbar")
        _sum_tlay = QHBoxLayout(_sum_toolbar)
        _sum_tlay.setContentsMargins(6, 3, 6, 3)
        _sum_tlay.addStretch()
        self._btn_copy_interp = QPushButton("📋  Copy Interpretation")
        self._btn_copy_interp.setFlat(True)
        self._btn_copy_interp.setEnabled(False)
        self._btn_copy_interp.setToolTip(
            "Copy LLM interpretation text to clipboard"
        )
        self._btn_copy_interp.clicked.connect(self._on_copy_interpretation)
        _sum_tlay.addWidget(self._btn_copy_interp)
        _sum_vlay.addWidget(_sum_toolbar)

        self._summary_browser = QTextBrowser()
        self._summary_browser.setPlaceholderText("Summary will appear here…")
        self._summary_browser.setOpenExternalLinks(False)
        _sum_vlay.addWidget(self._summary_browser)

        self._tabs.addTab(_sum_widget, _icon("summary"), "Summary")

        # Hover-reveal pop-out button — floats over the result container
        _HoverRevealButton(
            container=_result_container,
            get_canvas_fn=lambda: (
                w
                if isinstance(
                    w := self._result_inner_tabs.currentWidget(), MplCanvas
                )
                else None
            ),
        )

        self._content_vsplit.addWidget(self._tabs)

        # ── Bottom: chat panel (always visible, resizable) ───────────────────
        _chat_outer = QFrame()
        _chat_outer.setObjectName("ChatPanel")
        _chat_outer.setFrameShape(QFrame.Shape.StyledPanel)

        chat_vlay = QVBoxLayout(_chat_outer)
        chat_vlay.setContentsMargins(8, 6, 8, 8)
        chat_vlay.setSpacing(4)

        # Header row
        chat_hdr = QHBoxLayout()
        _chat_ic = _icon("chat")
        _ic_lbl = QLabel()
        if not _chat_ic.isNull():
            _ic_lbl.setPixmap(_chat_ic.pixmap(16, 16))
        else:
            _ic_lbl.setText("💬")
        chat_hdr.addWidget(_ic_lbl)

        _chat_title = QLabel("Chat with Agent")
        _f = _chat_title.font()
        _f.setBold(True)
        _chat_title.setFont(_f)
        chat_hdr.addWidget(_chat_title)
        chat_hdr.addStretch()

        self._btn_chat_clear = QPushButton("🗑  Clear")
        self._btn_chat_clear.setFlat(True)
        self._btn_chat_clear.setFixedHeight(22)
        self._btn_chat_clear.setToolTip("Clear chat history")
        self._btn_chat_clear.clicked.connect(self._on_chat_clear)
        chat_hdr.addWidget(self._btn_chat_clear)
        chat_vlay.addLayout(chat_hdr)

        _sep = QFrame()
        _sep.setFrameShape(QFrame.Shape.HLine)
        _sep.setObjectName("Separator")
        chat_vlay.addWidget(_sep)

        # Chat history browser
        self._chat_history = QTextBrowser()
        self._chat_history.setOpenExternalLinks(False)
        self._chat_history.setObjectName("ChatHistory")
        self._chat_history.setHtml(_CHAT_WELCOME_HTML)
        chat_vlay.addWidget(self._chat_history)

        # Multi-line input + Send button
        chat_inp_row = QHBoxLayout()
        self._chat_input = QPlainTextEdit()
        self._chat_input.setPlaceholderText("Ask a follow-up question…")
        self._chat_input.setMaximumHeight(72)
        self._chat_input.setMinimumHeight(36)
        self._chat_input.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        self._chat_input.installEventFilter(
            self
        )  # Enter=send, Shift+Enter=newline
        chat_inp_row.addWidget(self._chat_input)

        self._btn_chat_send = QPushButton("▶")
        self._btn_chat_send.setFixedSize(36, 36)
        self._btn_chat_send.setToolTip("Send  (Enter)")
        self._btn_chat_send.clicked.connect(self._on_chat_send)
        chat_inp_row.addWidget(
            self._btn_chat_send, 0, Qt.AlignmentFlag.AlignBottom
        )
        chat_vlay.addLayout(chat_inp_row)

        _hint = QLabel("Enter to send  ·  Shift+Enter for new line")
        _hint.setObjectName("InfoLabel")
        chat_vlay.addWidget(_hint)

        self._chat_status = QLabel("")
        self._chat_status.setObjectName("InfoLabel")
        self._chat_status.setWordWrap(True)
        chat_vlay.addWidget(self._chat_status)

        self._content_vsplit.addWidget(_chat_outer)

        # Default split: ~65 % tabs / ~35 % chat
        self._content_vsplit.setSizes([560, 280])

        layout.addWidget(self._content_vsplit)

    # ── Public API ─────────────────────────────────────────────────────────────

    def set_app_controller(self, ctrl) -> None:
        self._ctrl = ctrl
        api_key = getattr(getattr(ctrl, "session", None), "api_key", "") or ""
        self._edit_apikey.setText(api_key)

    def append_log(self, text: str) -> None:
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self._log_text.appendPlainText(f"[{ts}]  {text}")

    @property
    def result_figure(self):
        return getattr(self, "_last_figure", None)

    # ── Agent selection ────────────────────────────────────────────────────────

    @Slot(str)
    def _on_agent_selected(self, name: str) -> None:
        self._current_agent = name
        self._detail_header.setText(name or "No agent selected")
        self._rebuild_params_form(name)
        self._btn_run.setEnabled(bool(name))

    @Slot(str)
    def _on_agent_activated(self, name: str) -> None:
        self._on_agent_selected(name)
        self._on_run()

    # ── Run / Stop ─────────────────────────────────────────────────────────────

    def _on_run(self) -> None:
        from pycsamt.app.desktop.workers.agent_worker import (
            AgentWorker,
        )

        if self._ctrl is None or getattr(self._ctrl, "sites", None) is None:
            self._status_lbl.setText("Load survey data first.")
            return

        agent_name = self._current_agent
        if not agent_name:
            self._status_lbl.setText("Select an agent first.")
            return

        params = self._collect_params()
        api_key = self._edit_apikey.text().strip() or None

        self._log_text.clear()
        self._summary_browser.clear()
        self._progress.setValue(0)
        self._progress.setVisible(True)
        self._btn_run.setEnabled(False)
        self._btn_stop.setEnabled(True)
        self._status_lbl.setText(f"Running: {agent_name}…")
        self._tabs.setCurrentIndex(0)

        # Reset chat to welcome state for each new run
        self._chat_history.setHtml(_CHAT_WELCOME_HTML)
        self._chat_status.setText("")

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
        self._worker.finished.connect(self._on_finished)
        self._worker.agent_ready.connect(self._on_agent_ready)
        self._worker.start()

    def _on_stop(self) -> None:
        if self._worker and self._worker.isRunning():
            self._worker.cancel()
            self._worker.quit()
            self._worker.wait(2000)
        self._reset_run_ui()

    # ── Worker signal handlers ─────────────────────────────────────────────────

    @Slot(str)
    def _on_log_line(self, line: str) -> None:
        self.append_log(line)

    @Slot(object)
    def _on_result(self, result) -> None:
        from pycsamt.app.desktop.panels.agent_panel import (
            _extract_all_renderables,
        )

        figs = _extract_all_renderables(result)
        if figs:
            # Clear previous figure tabs (keep the placeholder canvas aside)
            while self._result_inner_tabs.count():
                self._result_inner_tabs.removeTab(0)

            for i, (label, renderable) in enumerate(figs):
                canvas = (
                    self._result_canvas
                    if i == 0
                    else MplCanvas(self, toolbar=True)
                )
                try:
                    canvas.show_figure(renderable)
                except Exception:
                    pass
                self._result_inner_tabs.addTab(canvas, label)

            self._last_figure = (
                figs[0][1].get_figure()
                if hasattr(figs[0][1], "get_figure")
                else figs[0][1]
            )
            self._tabs.setCurrentIndex(1)  # Result tab
        else:
            self._tabs.setCurrentIndex(2)  # Summary tab

        self._summary_browser.setHtml(_format_result_html(result))

        # Store LLM interpretation for copy button
        self._last_interpretation = (
            getattr(result, "llm_interpretation", None) or ""
        )
        self._btn_copy_interp.setEnabled(bool(self._last_interpretation))

        # Build chat context
        ctx_parts: list[str] = []
        summary = getattr(result, "summary", None)
        if summary:
            ctx_parts.append(f"Summary: {summary}")
        if self._last_interpretation:
            ctx_parts.append(f"LLM analysis:\n{self._last_interpretation}")
        if ctx_parts:
            self._last_result_context = "\n\n".join(ctx_parts)

        # Drop a system message into chat when an LLM agent completed
        if self._last_agent is not None:
            self._append_chat_message(
                "system",
                "Analysis complete — ask me a follow-up question about the results.",
            )

    @Slot(str)
    def _on_error(self, msg: str) -> None:
        self._status_lbl.setText(f"ERROR: {msg}")
        self.append_log(f"ERROR: {msg}")

    @Slot()
    def _on_finished(self) -> None:
        self._reset_run_ui()
        self._status_lbl.setText("Finished.")

    # ── Form helpers ──────────────────────────────────────────────────────────

    def _rebuild_params_form(self, agent_name: str) -> None:
        from pycsamt.app.desktop.agent_registry import (
            AGENT_REGISTRY,
        )

        while self._form_params.rowCount():
            self._form_params.removeRow(0)
        self._params_widgets.clear()

        meta = AGENT_REGISTRY.get(agent_name, {})
        params_spec = meta.get("params", {})
        is_llm = meta.get("type") == "llm"

        for param_name, spec in params_spec.items():
            w = self._make_param_widget(spec)
            if spec.get("tip"):
                w.setToolTip(spec["tip"])
            label = spec.get("label", param_name)
            self._form_params.addRow(f"{label}:", w)
            self._params_widgets[param_name] = w

        self._grp_params.setVisible(bool(params_spec))
        self._grp_context.setVisible(is_llm)
        self._grp_apikey.setVisible(is_llm)

    @staticmethod
    def _make_param_widget(spec: dict) -> QWidget:
        kind = spec.get("type", "str")
        default = spec.get("default")

        if kind == "combo":
            w = QComboBox()
            for opt in spec.get("options", []):
                w.addItem(str(opt))
            if default is not None:
                idx = w.findText(str(default))
                if idx >= 0:
                    w.setCurrentIndex(idx)
            return w

        if kind == "float":
            lo, hi = spec.get("range", (0.0, 1e9))
            w = QDoubleSpinBox()
            w.setRange(lo, hi)
            w.setSingleStep(spec.get("step", 0.1))
            w.setDecimals(spec.get("decimals", 3))
            if default is not None:
                w.setValue(float(default))
            return w

        if kind == "int":
            lo, hi = spec.get("range", (0, 10_000))
            w = QSpinBox()
            w.setRange(int(lo), int(hi))
            w.setSingleStep(spec.get("step", 1))
            if default is not None:
                w.setValue(int(default))
            return w

        if kind == "bool":
            w = QCheckBox()
            w.setChecked(bool(default))
            return w

        w = QLineEdit()
        if default is not None:
            w.setText(str(default))
        return w

    def _collect_params(self) -> dict:
        from pycsamt.app.desktop.agent_registry import (
            AGENT_REGISTRY,
        )

        meta = AGENT_REGISTRY.get(self._current_agent, {})
        params_spec = meta.get("params", {})
        result: dict = {}

        for param_name, _spec in params_spec.items():
            w = self._params_widgets.get(param_name)
            if w is None:
                continue
            if isinstance(w, QComboBox):
                result[param_name] = w.currentText()
            elif isinstance(w, (QDoubleSpinBox, QSpinBox)):
                result[param_name] = w.value()
            elif isinstance(w, QCheckBox):
                result[param_name] = w.isChecked()
            elif isinstance(w, QLineEdit):
                result[param_name] = w.text()

        ctx = self._ctx_edit.toPlainText().strip()
        if ctx:
            result["user_prompt"] = ctx
            result["context"] = ctx
        return result

    def _reset_run_ui(self) -> None:
        self._progress.setVisible(False)
        self._btn_run.setEnabled(bool(self._current_agent))
        self._btn_stop.setEnabled(False)

    # ── Agent-ready ───────────────────────────────────────────────────────────

    @Slot(object)
    def _on_agent_ready(self, agent) -> None:
        self._last_agent = agent

    # ── Summary: copy interpretation ─────────────────────────────────────────

    def _on_copy_interpretation(self) -> None:
        if self._last_interpretation:
            QApplication.clipboard().setText(self._last_interpretation)

    # ── Chat ──────────────────────────────────────────────────────────────────

    def eventFilter(self, obj, event) -> bool:
        """Intercept Enter in the chat input to send instead of inserting a newline."""
        if obj is self._chat_input and event.type() == QEvent.Type.KeyPress:
            ke = QKeyEvent(event)
            if ke.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter) and not (
                ke.modifiers() & Qt.KeyboardModifier.ShiftModifier
            ):
                self._on_chat_send()
                return True
        return super().eventFilter(obj, event)

    def _on_chat_clear(self) -> None:
        self._chat_history.setHtml(_CHAT_WELCOME_HTML)
        self._chat_status.setText("")

    def _on_chat_send(self) -> None:
        from pycsamt.app.desktop.workers.agent_worker import (
            ChatWorker,
        )

        question = self._chat_input.toPlainText().strip()
        if not question:
            return
        if self._last_agent is None:
            self._chat_status.setText("Run an LLM agent first.")
            return
        if not getattr(self._last_agent, "api_key", None):
            self._chat_status.setText(
                "No API key — configure one in the API Key field and re-run the agent."
            )
            return

        self._append_chat_message("user", question)
        self._chat_input.clear()
        self._btn_chat_send.setEnabled(False)
        self._chat_input.setEnabled(False)
        self._chat_status.setText("Asking LLM…")

        self._chat_worker = ChatWorker(
            self._last_agent,
            question,
            self._last_result_context,
            parent=self,
        )
        self._chat_worker.reply_done.connect(self._on_chat_reply)
        self._chat_worker.error.connect(self._on_chat_error)
        self._chat_worker.start()

    @Slot(str)
    def _on_chat_reply(self, text: str) -> None:
        self._btn_chat_send.setEnabled(True)
        self._chat_input.setEnabled(True)
        self._chat_status.setText("")
        self._append_chat_message("agent", text)
        self._chat_input.setFocus()

    @Slot(str)
    def _on_chat_error(self, msg: str) -> None:
        self._btn_chat_send.setEnabled(True)
        self._chat_input.setEnabled(True)
        self._chat_status.setText(f"Error: {msg}")

    def _append_chat_message(self, role: str, text: str) -> None:
        safe = (
            text.replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace("\n", "<br/>")
        )
        if role == "system":
            html = (
                f"<p style='color:#666;text-align:center;font-size:11px;"
                f"font-style:italic;margin:6px 0;'>{safe}</p>"
            )
        elif role == "user":
            html = (
                "<div style='background:#2a3a52;color:#cdd6f4;"
                "text-align:right;margin:4px 2px;padding:8px 12px;"
                "border-radius:12px 12px 4px 12px;font-size:12px;'>"
                f"<small><b>You</b></small><br/>{safe}</div>"
            )
        else:  # agent
            html = (
                "<div style='background:#1e3428;color:#cdd6f4;"
                "text-align:left;margin:4px 2px;padding:8px 12px;"
                "border-radius:12px 12px 12px 4px;font-size:12px;'>"
                f"<small><b>Agent</b></small><br/>{safe}</div>"
            )
        self._chat_history.append(html)
        self._chat_history.ensureCursorVisible()

    # ── Theme ──────────────────────────────────────────────────────────────────

    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._browser.apply_theme(dark)
