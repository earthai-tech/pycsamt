# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Solver console: a styled, searchable, detachable scrolling log.

``ConsolePanel`` shows solver output as it streams (ModEM, Occam2D,
MARE2DEM print hundreds of lines per iteration).  Lines are queued and
flushed every 100 ms, so a chatty solver never floods the GUI thread.

Themes: *Terminal* (green on near-black, the default), *Amber*, and
*Follow app* (the application palette).  Errors, warnings and
iteration/RMS lines are coloured.  The pop-out button moves the panel into
its own window (``ConsoleWindow``) that keeps streaming while the user
works elsewhere; closing that window docks the panel back.
"""

from __future__ import annotations

import re
from pathlib import Path

from PySide6.QtCore import QRegularExpression, Qt, QTimer, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QFontDatabase,
    QKeySequence,
    QShortcut,
    QSyntaxHighlighter,
    QTextCharFormat,
    QTextDocument,
)
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

# name -> (background, text, error, warning, highlight, dim)
THEMES: dict[str, tuple[str, ...]] = {
    "terminal": ("#0b0f0c", "#3ddc72", "#ff6b6b", "#ffd43b", "#74c0fc",
                 "#6c8a74"),
    "amber": ("#140e02", "#ffb000", "#ff6b6b", "#ffe066", "#fff3bf",
              "#8c6d1f"),
}
THEME_LABELS = {"terminal": "Terminal", "amber": "Amber", "app": "Follow app"}

_ERR = re.compile(r"\b(error|fatal|failed|segmentation|traceback|abort)\b",
                  re.I)
_WARN = re.compile(r"\b(warn(ing)?|cannot|could not|not found)\b", re.I)
# The application stylesheet sets a proportional font on every widget;
# the console pins a monospace family in its own stylesheet.
_MONO = ("font-family: Consolas, 'Cascadia Mono', 'DejaVu Sans Mono', "
         "'Courier New', monospace; font-size: 9pt;")
_ITER = re.compile(r"\b(iter(ation)?\b|rms|misfit)", re.I)


class _Highlighter(QSyntaxHighlighter):
    def __init__(self, doc: QTextDocument) -> None:
        super().__init__(doc)
        self._fmt: dict[str, QTextCharFormat] = {}
        self.set_colors("#ff6b6b", "#ffd43b", "#74c0fc", "#6c8a74")

    def set_colors(self, err, warn, hi, dim) -> None:
        def f(color, bold=False):
            fmt = QTextCharFormat()
            fmt.setForeground(QColor(color))
            if bold:
                fmt.setFontWeight(QFont.Weight.Bold)
            return fmt

        self._fmt = {"err": f(err, True), "warn": f(warn), "hi": f(hi, True),
                     "dim": f(dim)}
        self.rehighlight()

    def highlightBlock(self, text: str) -> None:  # noqa: N802
        if not text:
            return
        if text.startswith(("$ ", "──", "==")):
            kind = "dim"
        elif _ERR.search(text):
            kind = "err"
        elif _WARN.search(text):
            kind = "warn"
        elif _ITER.search(text):
            kind = "hi"
        else:
            return
        self.setFormat(0, len(text), self._fmt[kind])


class ConsoleView(QPlainTextEdit):
    """Read-only monospace log with batched appends."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("SolverConsole")
        self.setReadOnly(True)
        self.setMaximumBlockCount(20000)
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        font.setPointSize(9)
        self.setFont(font)
        self.follow = True
        self._pending: list[str] = []
        self._timer = QTimer(self)
        self._timer.setInterval(100)
        self._timer.timeout.connect(self.flush)
        self._timer.start()
        self._hl = _Highlighter(self.document())
        self.theme = ""
        self.set_theme("terminal")

    def append_line(self, line: str) -> None:
        self._pending.append(str(line))

    def flush(self) -> None:
        if not self._pending:
            return
        text = "\n".join(self._pending)
        self._pending.clear()
        self.appendPlainText(text)
        if self.follow:
            bar = self.verticalScrollBar()
            bar.setValue(bar.maximum())

    def set_theme(self, name: str) -> None:
        self.theme = name if name in THEME_LABELS else "terminal"
        if self.theme == "app":
            self.setStyleSheet(f"QPlainTextEdit#SolverConsole {{ {_MONO} }}")
            self._hl.set_colors("#c92a2a", "#9c4f00", "#1864ab", "#5f6b7a")
            return
        bg, fg, err, warn, hi, dim = THEMES[self.theme]
        self.setStyleSheet(
            f"QPlainTextEdit#SolverConsole {{ background: {bg}; color: {fg};"
            f" border: none; selection-background-color: {dim}; {_MONO} }}")
        self._hl.set_colors(err, warn, hi, dim)

    def lines(self) -> list[str]:
        self.flush()
        return self.toPlainText().splitlines()


class ConsolePanel(QWidget):
    """Toolbar + :class:`ConsoleView`; can be popped out of its host."""

    pop_out_requested = Signal()
    theme_changed = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(2)
        bar = QHBoxLayout()
        bar.setContentsMargins(4, 2, 4, 0)
        bar.setSpacing(4)
        self.title = QLabel("Console")
        self.title.setStyleSheet("font-weight: 600;")
        bar.addWidget(self.title)
        self.status = QLabel("")
        self.status.setObjectName("InfoLabel")
        bar.addWidget(self.status, 1)

        self.search = QLineEdit()
        self.search.setPlaceholderText("Find…  (Enter = next)")
        self.search.setClearButtonEnabled(True)
        self.search.setMaximumWidth(190)
        self.search.returnPressed.connect(self.find_next)
        bar.addWidget(self.search)

        def tool(text, tip, checkable=False, checked=False):
            b = QToolButton()
            b.setText(text)
            b.setToolTip(tip)
            b.setCheckable(checkable)
            b.setChecked(checked)
            b.setAutoRaise(True)
            bar.addWidget(b)
            return b

        self.btn_follow = tool("Follow", "Keep scrolling to the newest line",
                               True, True)
        self.btn_follow.toggled.connect(self._set_follow)
        self.btn_wrap = tool("Wrap", "Wrap long lines", True, False)
        self.btn_wrap.toggled.connect(self._set_wrap)
        self.theme_combo = QComboBox()
        for key, label in THEME_LABELS.items():
            self.theme_combo.addItem(label, key)
        self.theme_combo.setToolTip("Console colours")
        self.theme_combo.currentIndexChanged.connect(self._on_theme)
        bar.addWidget(self.theme_combo)
        self.btn_copy = tool("Copy", "Copy the whole log")
        self.btn_copy.clicked.connect(self.copy_all)
        self.btn_save = tool("Save…", "Save the log to a text file")
        self.btn_save.clicked.connect(self.save_dialog)
        self.btn_clear = tool("Clear", "Clear the console")
        self.btn_clear.clicked.connect(self.clear)
        self.btn_pop = tool("Pop out", "Open the console in its own window "
                                       "— it keeps streaming while you work "
                                       "elsewhere")
        from pycsamt.app.desktop.windows._base import _icon

        ic = _icon("pop-out")
        if not ic.isNull():
            self.btn_pop.setIcon(ic)
            self.btn_pop.setToolButtonStyle(
                Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.btn_pop.clicked.connect(self.pop_out_requested)
        v.addLayout(bar)

        self.view = ConsoleView(self)
        v.addWidget(self.view, 1)
        QShortcut(QKeySequence.StandardKey.Find, self,
                  activated=lambda: self.search.setFocus())

    # ── API ───────────────────────────────────────────────────────────
    def append(self, line: str) -> None:
        self.view.append_line(line)

    def clear(self) -> None:
        self.view._pending.clear()
        self.view.clear()

    def text(self) -> str:
        self.view.flush()
        return self.view.toPlainText()

    def set_theme(self, name: str) -> None:
        i = self.theme_combo.findData(name)
        self.theme_combo.setCurrentIndex(max(i, 0))

    def copy_all(self) -> None:
        QApplication.clipboard().setText(self.text())

    def save_to(self, path) -> Path:
        p = Path(path)
        p.write_text(self.text() + "\n", encoding="utf-8")
        return p

    def save_dialog(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save console log", "run.log", "Log files (*.log *.txt)")
        if path:
            self.save_to(path)

    def find_next(self) -> bool:
        text = self.search.text()
        if not text:
            return False
        self.view.flush()  # include lines still queued for display
        rx = QRegularExpression(QRegularExpression.escape(text),
                                QRegularExpression.PatternOption
                                .CaseInsensitiveOption)
        found = self.view.find(rx)
        if not found:  # wrap around
            cur = self.view.textCursor()
            cur.movePosition(cur.MoveOperation.Start)
            self.view.setTextCursor(cur)
            found = self.view.find(rx)
        if found:
            self.btn_follow.setChecked(False)
        return bool(found)

    # ── internals ─────────────────────────────────────────────────────
    def _set_follow(self, on: bool) -> None:
        self.view.follow = on
        if on:
            bar = self.view.verticalScrollBar()
            bar.setValue(bar.maximum())

    def _set_wrap(self, on: bool) -> None:
        self.view.setLineWrapMode(
            QPlainTextEdit.LineWrapMode.WidgetWidth if on
            else QPlainTextEdit.LineWrapMode.NoWrap)

    def _on_theme(self, _i: int) -> None:
        key = self.theme_combo.currentData()
        self.view.set_theme(key)
        self.theme_changed.emit(key)


class ConsoleWindow(QWidget):
    """Top-level home of a popped-out :class:`ConsolePanel`."""

    dock_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        flags = (Qt.WindowType.Window | Qt.WindowType.WindowCloseButtonHint
                 | Qt.WindowType.WindowMinimizeButtonHint
                 | Qt.WindowType.WindowMaximizeButtonHint)
        super().__init__(parent, flags)
        self.setWindowTitle("pycsamt — Inversion console")
        self.resize(900, 420)
        self._lay = QVBoxLayout(self)
        self._lay.setContentsMargins(4, 4, 4, 4)

    def hold(self, panel: ConsolePanel) -> None:
        self._lay.addWidget(panel)
        panel.show()

    def closeEvent(self, event) -> None:  # noqa: N802
        event.ignore()
        self.hide()
        self.dock_requested.emit()


__all__ = ["THEMES", "ConsolePanel", "ConsoleView", "ConsoleWindow"]
