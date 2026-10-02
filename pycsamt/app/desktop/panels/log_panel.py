# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
LogPanel — the main window's execution log (bottom dock, hidden by default).

Every line gets a timestamp and a level -- ``info``, ``ok``, ``warn`` or
``error`` -- either passed by the caller or read from the text ("ERROR:",
"failed", "warning", "✓", "✕", ...).  Levels are coloured with a light or
a dark palette (:meth:`LogPanel.set_dark`), can be filtered, searched,
copied and saved.  :attr:`LogPanel.message_logged` lets the main window
keep an unread count on its status-bar chip while the dock is hidden.
"""

from __future__ import annotations

import datetime
import re
from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtGui import QColor, QFont, QTextCharFormat, QTextCursor
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

LEVELS = ("info", "ok", "warn", "error")

# level -> (light, dark) text colours; the background follows the theme
_COLOURS = {
    "time": ("#6b7280", "#7f8796"),
    "info": ("#1f2937", "#d8dee9"),
    "ok": ("#1b7f3b", "#6fd08c"),
    "warn": ("#9a5b00", "#f2c14e"),
    "error": ("#b42318", "#ff7b72"),
}
_BACKGROUND = ("#ffffff", "#161b22")
_BORDER = ("#d0d7de", "#30363d")

_ERROR_RE = re.compile(r"\b(error|failed|fatal|exception|traceback)\b|✕",
                       re.I)
_WARN_RE = re.compile(r"\b(warn(ing)?|could not|cannot|not found|skipp)",
                      re.I)
_OK_RE = re.compile(r"✓|\b(done|saved|loaded|finished|complete[d]?|"
                    r"exported)\b", re.I)
_MONO = ("Consolas, 'Cascadia Mono', 'DejaVu Sans Mono', 'Courier New', "
         "monospace")


def classify(text: str) -> str:
    """The level a log line reads as (``info`` when nothing matches)."""
    if _ERROR_RE.search(text):
        return "error"
    if _WARN_RE.search(text):
        return "warn"
    if _OK_RE.search(text):
        return "ok"
    return "info"


class LogPanel(QWidget):
    """Levelled, themed, filterable execution log."""

    message_logged = Signal(str, str)  # level, text

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._dark = False
        self._records: list[tuple[str, str, str]] = []  # (time, level, text)
        self._max = 2000
        self._build_ui()
        self._apply_style()

    # ── UI ────────────────────────────────────────────────────────────
    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(4, 2, 4, 4)
        root.setSpacing(2)

        bar = QHBoxLayout()
        bar.setContentsMargins(0, 0, 0, 0)
        bar.setSpacing(4)
        self._level = QComboBox()
        for key, label in (("all", "All messages"), ("warn", "Warnings + "
                           "errors"), ("error", "Errors only")):
            self._level.addItem(label, key)
        self._level.setToolTip("Show messages of this level and above")
        self._level.currentIndexChanged.connect(lambda _i: self._rebuild())
        bar.addWidget(self._level)
        self._find = QLineEdit()
        self._find.setPlaceholderText("Find…")
        self._find.setClearButtonEnabled(True)
        self._find.setMaximumWidth(220)
        self._find.textChanged.connect(lambda _t: self._rebuild())
        bar.addWidget(self._find)
        bar.addStretch(1)
        self._counts = QLabel("")
        self._counts.setObjectName("InfoLabel")
        bar.addWidget(self._counts)

        def tool(text, tip, slot, checkable=False):
            b = QToolButton()
            b.setText(text)
            b.setToolTip(tip)
            b.setAutoRaise(True)
            b.setCheckable(checkable)
            if checkable:
                b.toggled.connect(slot)
            else:
                b.clicked.connect(slot)
            bar.addWidget(b)
            return b

        self._btn_wrap = tool("Wrap", "Wrap long lines", self._set_wrap,
                              checkable=True)
        tool("Copy", "Copy the shown lines", self.copy_all)
        tool("Save…", "Save the whole log to a text file", self.save_dialog)
        tool("Clear", "Clear the log", self.clear)
        root.addLayout(bar)

        self._text = QPlainTextEdit(self)
        self._text.setReadOnly(True)
        self._text.setMaximumBlockCount(self._max)
        self._text.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        font = QFont()
        font.setStyleHint(QFont.StyleHint.TypeWriter)
        font.setFamilies([f.strip(" '") for f in _MONO.split(",")])
        font.setPointSize(9)
        self._text.setFont(font)
        root.addWidget(self._text, 1)

    def _apply_style(self) -> None:
        i = 1 if self._dark else 0
        self._text.setStyleSheet(
            f"QPlainTextEdit {{ background: {_BACKGROUND[i]}; color: "
            f"{_COLOURS['info'][i]}; border: 1px solid {_BORDER[i]}; "
            f"border-radius: 4px; font-family: {_MONO}; font-size: 9pt; }}")

    # ── public API ────────────────────────────────────────────────────
    def append_line(self, text: str, level: str | None = None) -> None:
        """Log *text* (level read from the text unless given)."""
        text = str(text)
        level = level if level in LEVELS else classify(text)
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self._records.append((ts, level, text))
        if len(self._records) > self._max:
            del self._records[: len(self._records) - self._max]
        if self._shown(level, text):
            self._write(ts, level, text)
        self._update_counts()
        self.message_logged.emit(level, text)

    def set_dark(self, dark: bool) -> None:
        """Recolour for the dark or the light theme."""
        self._dark = bool(dark)
        self._apply_style()
        self._rebuild()

    def clear(self) -> None:
        self._records.clear()
        self._text.clear()
        self._update_counts()

    def text(self) -> str:
        return self._text.toPlainText()

    def records(self) -> list[tuple[str, str, str]]:
        return list(self._records)

    def copy_all(self) -> None:
        QApplication.clipboard().setText(self.text())

    def save_to(self, path) -> Path:
        path = Path(path)
        path.write_text("\n".join(f"[{t}] {lv.upper():5s} {x}"
                                  for t, lv, x in self._records) + "\n",
                        encoding="utf-8")
        return path

    def save_dialog(self) -> None:
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save log", f"pycsamt_log_{stamp}.txt",
            "Text files (*.txt);;All files (*)")
        if path:
            self.save_to(path)

    # ── internals ─────────────────────────────────────────────────────
    def _shown(self, level: str, text: str) -> bool:
        floor = self._level.currentData() or "all"
        if floor == "error" and level != "error":
            return False
        if floor == "warn" and level not in ("warn", "error"):
            return False
        needle = self._find.text().strip().lower()
        return not needle or needle in text.lower()

    def _write(self, ts: str, level: str, text: str) -> None:
        i = 1 if self._dark else 0
        cur = self._text.textCursor()
        cur.movePosition(QTextCursor.MoveOperation.End)
        if not self._text.document().isEmpty():
            cur.insertBlock()
        fmt = QTextCharFormat()
        fmt.setForeground(QColor(_COLOURS["time"][i]))
        cur.insertText(f"[{ts}]  ", fmt)
        fmt = QTextCharFormat()
        fmt.setForeground(QColor(_COLOURS[level][i]))
        if level == "error":
            fmt.setFontWeight(QFont.Weight.Bold)
        cur.insertText(text, fmt)
        self._text.setTextCursor(cur)
        self._text.ensureCursorVisible()

    def _rebuild(self) -> None:
        self._text.clear()
        for ts, level, text in self._records:
            if self._shown(level, text):
                self._write(ts, level, text)

    def _update_counts(self) -> None:
        n_err = sum(1 for _t, lv, _x in self._records if lv == "error")
        n_warn = sum(1 for _t, lv, _x in self._records if lv == "warn")
        bits = [f"{len(self._records)} lines"]
        if n_warn:
            bits.append(f"{n_warn} warning{'s' if n_warn > 1 else ''}")
        if n_err:
            bits.append(f"{n_err} error{'s' if n_err > 1 else ''}")
        self._counts.setText(" · ".join(bits))

    def _set_wrap(self, on: bool) -> None:
        self._text.setLineWrapMode(
            QPlainTextEdit.LineWrapMode.WidgetWidth if on
            else QPlainTextEdit.LineWrapMode.NoWrap)


__all__ = ["LEVELS", "LogPanel", "classify"]
