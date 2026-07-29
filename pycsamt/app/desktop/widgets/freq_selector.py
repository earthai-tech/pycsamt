# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
FreqSelector — dual-handle log-scale range slider for frequency / period.

Phase 3: QWidget showing a horizontal track with two draggable handles
spanning the full period/frequency range present in the loaded data.
Emits ``range_changed(f_min_hz, f_max_hz)`` whenever either handle moves.
A toggle button switches the axis label and tick display between Hz and
seconds (period).
"""

from __future__ import annotations

import math
from pathlib import Path

from PySide6.QtCore import QPoint, Qt, Signal
from PySide6.QtGui import QColor, QFont, QIcon, QPainter, QPen
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

_ICONS = Path(__file__).parent.parent / "resources" / "icons"


def _icon(name: str) -> QIcon:
    for c in (name, f"{name}.svg", f"{name}.png"):
        p = _ICONS / c
        if p.exists():
            return QIcon(str(p))
    return QIcon()


# ──────────────────────────────────────────────────────────────────────────────
# Colours (match dark_theme.qss palette)
# ──────────────────────────────────────────────────────────────────────────────
_C_TRACK = QColor("#313244")
_C_RANGE = QColor("#89b4fa")
_C_HANDLE = QColor("#cdd6f4")
_C_BORDER = QColor("#89b4fa")
_C_TICK = QColor("#585b70")
_C_TEXT = QColor("#a6adc8")
_C_RANGE_A = QColor(137, 180, 250, 60)  # translucent fill

_LOG_DECADES = [0.001, 0.01, 0.1, 1, 10, 100, 1_000, 10_000]
_HANDLE_R = 7
_TRACK_H = 5
_TICK_H = 6


# ──────────────────────────────────────────────────────────────────────────────
# Internal range-slider widget
# ──────────────────────────────────────────────────────────────────────────────


class _RangeSlider(QWidget):
    """Low-level dual-handle log-frequency slider (no labels)."""

    range_changed = Signal(float, float)  # (lo_hz, hi_hz)

    def __init__(
        self,
        f_min: float = 0.001,
        f_max: float = 10_000.0,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._f_min = max(f_min, 1e-9)
        self._f_max = max(f_max, self._f_min * 10)
        self._lo_hz = self._f_min
        self._hi_hz = self._f_max
        self._drag = None  # 'lo' | 'hi'
        self.setMouseTracking(True)
        self.setMinimumHeight(36)
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    # ── Coordinate helpers ─────────────────────────────────────────────

    def _pad(self) -> int:
        return _HANDLE_R + 4

    def _track_width(self) -> int:
        return max(1, self.width() - 2 * self._pad())

    def _hz_to_frac(self, hz: float) -> float:
        lo = math.log10(self._f_min)
        hi = math.log10(self._f_max)
        return (math.log10(max(hz, self._f_min)) - lo) / (hi - lo)

    def _frac_to_hz(self, frac: float) -> float:
        lo = math.log10(self._f_min)
        hi = math.log10(self._f_max)
        return 10 ** (lo + max(0.0, min(1.0, frac)) * (hi - lo))

    def _px(self, hz: float) -> int:
        return self._pad() + int(self._hz_to_frac(hz) * self._track_width())

    def _px_to_hz(self, x: float) -> float:
        frac = (x - self._pad()) / max(1, self._track_width())
        return self._frac_to_hz(frac)

    # ── Public API ─────────────────────────────────────────────────────

    def set_range(self, f_min: float, f_max: float) -> None:
        self._f_min = max(f_min, 1e-9)
        self._f_max = max(f_max, self._f_min * 10)
        self._lo_hz = self._f_min
        self._hi_hz = self._f_max
        self.update()

    def set_selection(self, lo_hz: float, hi_hz: float) -> None:
        self._lo_hz = max(self._f_min, min(lo_hz, self._f_max))
        self._hi_hz = max(self._f_min, min(hi_hz, self._f_max))
        self.update()

    @property
    def lo_hz(self) -> float:
        return self._lo_hz

    @property
    def hi_hz(self) -> float:
        return self._hi_hz

    # ── Paint ──────────────────────────────────────────────────────────

    def paintEvent(self, event) -> None:
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)

        mid_y = self.height() // 2
        lo_px = self._px(self._lo_hz)
        hi_px = self._px(self._hi_hz)
        pad = self._pad()
        tw = self._track_width()

        # Background track
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(_C_TRACK)
        p.drawRoundedRect(pad, mid_y - _TRACK_H // 2, tw, _TRACK_H, 3, 3)

        # Selected range (fill)
        p.setBrush(_C_RANGE_A)
        p.drawRect(lo_px, mid_y - _TRACK_H // 2, hi_px - lo_px, _TRACK_H)

        # Selected range (top line)
        p.setPen(QPen(_C_RANGE, 2))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawLine(lo_px, mid_y - _TRACK_H // 2, hi_px, mid_y - _TRACK_H // 2)

        # Decade tick marks
        p.setPen(QPen(_C_TICK, 1))
        font = QFont()
        font.setPointSize(7)
        p.setFont(font)
        for decade in _LOG_DECADES:
            if self._f_min <= decade <= self._f_max:
                tx = self._px(decade)
                p.drawLine(tx, mid_y + _TRACK_H, tx, mid_y + _TRACK_H + _TICK_H)

        # Handles
        for px in (lo_px, hi_px):
            p.setPen(QPen(_C_BORDER, 2))
            p.setBrush(_C_HANDLE)
            p.drawEllipse(QPoint(px, mid_y), _HANDLE_R, _HANDLE_R)

    # ── Mouse interaction ──────────────────────────────────────────────

    def mousePressEvent(self, event) -> None:
        x = event.position().x()
        lo_px = self._px(self._lo_hz)
        hi_px = self._px(self._hi_hz)
        self._drag = "lo" if abs(x - lo_px) <= abs(x - hi_px) else "hi"
        self._move(x)

    def mouseMoveEvent(self, event) -> None:
        if self._drag:
            self._move(event.position().x())

    def mouseReleaseEvent(self, event) -> None:
        self._drag = None

    def _move(self, x: float) -> None:
        hz = self._px_to_hz(x)
        if self._drag == "lo":
            self._lo_hz = max(self._f_min, min(hz, self._hi_hz))
        else:
            self._hi_hz = min(self._f_max, max(hz, self._lo_hz))
        self.update()
        self.range_changed.emit(self._lo_hz, self._hi_hz)


# ──────────────────────────────────────────────────────────────────────────────
# Public FreqSelector widget
# ──────────────────────────────────────────────────────────────────────────────


class FreqSelector(QWidget):
    """
    Dual-handle log-scale range selector for frequency / period.

    Parameters
    ----------
    f_min, f_max : float
        Full available frequency range in Hz.
    parent : QWidget, optional

    Signals
    -------
    range_changed(float, float)
        Emitted on handle move: ``(lo_hz, hi_hz)`` always in Hz.
    """

    range_changed = Signal(float, float)  # (lo_hz, hi_hz)

    def __init__(
        self,
        f_min: float = 0.001,
        f_max: float = 10_000.0,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._show_period = False
        self._build_ui(f_min, f_max)

    # ── Construction ──────────────────────────────────────────────────

    def _build_ui(self, f_min: float, f_max: float) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(4, 2, 4, 2)
        root.setSpacing(0)

        # Top row: label + log₁₀ value readout + unit toggle + reset
        top_row = QHBoxLayout()
        top_row.setContentsMargins(0, 0, 0, 0)

        # "Frequency range (log₁₀ Hz):" — updated by _toggle_unit
        self._unit_lbl = QLabel("Frequency range (log₁₀ Hz):")
        self._unit_lbl.setObjectName("FreqSelectorLabel")
        top_row.addWidget(self._unit_lbl)

        top_row.addStretch()

        self._lo_lbl = QLabel()
        self._hi_lbl = QLabel()
        self._lo_lbl.setObjectName("FreqValueLabel")
        self._hi_lbl.setObjectName("FreqValueLabel")
        top_row.addWidget(self._lo_lbl)
        top_row.addWidget(QLabel(" – "))
        top_row.addWidget(self._hi_lbl)

        # Toggle: "→ log₁₀(T)" switches to period view
        self._toggle_btn = QPushButton("→ log₁₀(T)")
        self._toggle_btn.setFixedWidth(90)
        self._toggle_btn.setObjectName("FreqToggleButton")
        self._toggle_btn.clicked.connect(self._toggle_unit)
        top_row.addWidget(self._toggle_btn)

        # Reset button — use reset.svg icon if available, else Unicode fallback
        self._reset_btn = QPushButton()
        ico = _icon("reset")
        if ico.isNull():
            self._reset_btn.setText("↺")
        else:
            self._reset_btn.setIcon(ico)
        self._reset_btn.setFixedWidth(28)
        self._reset_btn.setObjectName("FreqResetButton")
        self._reset_btn.setToolTip("Reset to full frequency / period range")
        self._reset_btn.clicked.connect(self._on_reset)
        top_row.addWidget(self._reset_btn)

        root.addLayout(top_row)

        # Slider
        self._slider = _RangeSlider(f_min, f_max, self)
        self._slider.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._slider.range_changed.connect(self._on_slider_changed)
        root.addWidget(self._slider)

        self._refresh_labels(f_min, f_max)

    # ── Public API ─────────────────────────────────────────────────────

    def set_freq_range(self, f_min: float, f_max: float) -> None:
        """Update the available frequency range (Hz)."""
        self._slider.set_range(f_min, f_max)
        self._refresh_labels(f_min, f_max)

    def set_selection(self, lo_hz: float, hi_hz: float) -> None:
        self._slider.set_selection(lo_hz, hi_hz)
        self._refresh_labels(lo_hz, hi_hz)

    @property
    def lo_hz(self) -> float:
        return self._slider.lo_hz

    @property
    def hi_hz(self) -> float:
        return self._slider.hi_hz

    # ── Internal ──────────────────────────────────────────────────────

    def _toggle_unit(self) -> None:
        self._show_period = not self._show_period
        if self._show_period:
            self._toggle_btn.setText("→ log₁₀(Hz)")
            self._unit_lbl.setText("Period range (log₁₀ s):")
        else:
            self._toggle_btn.setText("→ log₁₀(T)")
            self._unit_lbl.setText("Frequency range (log₁₀ Hz):")
        self._refresh_labels(self._slider.lo_hz, self._slider.hi_hz)

    def _on_reset(self) -> None:
        """Reset both handles to the full available range and notify listeners."""
        lo = self._slider._f_min
        hi = self._slider._f_max
        self._slider.set_selection(lo, hi)
        self._slider.update()
        self._refresh_labels(lo, hi)
        self.range_changed.emit(lo, hi)

    def _on_slider_changed(self, lo_hz: float, hi_hz: float) -> None:
        self._refresh_labels(lo_hz, hi_hz)
        self.range_changed.emit(lo_hz, hi_hz)

    def _refresh_labels(self, lo_hz: float, hi_hz: float) -> None:
        if self._show_period:
            # display log₁₀(T) — smaller T (higher f) is the left handle
            lo_T = 1.0 / hi_hz if hi_hz > 0 else float("inf")
            hi_T = 1.0 / lo_hz if lo_hz > 0 else float("inf")
            self._lo_lbl.setText(self._fmt_log10(lo_T))
            self._hi_lbl.setText(self._fmt_log10(hi_T))
        else:
            # display log₁₀(Hz)
            self._lo_lbl.setText(self._fmt_log10(lo_hz))
            self._hi_lbl.setText(self._fmt_log10(hi_hz))

    @staticmethod
    def _fmt_log10(val: float) -> str:
        """Format a positive real value as its log₁₀, e.g. 1000 → '3.00'."""
        if val <= 0 or not math.isfinite(val):
            return "—"
        v = math.log10(val)
        # integer decade: no decimal; otherwise 2 d.p.
        return f"{int(round(v))}" if abs(v - round(v)) < 0.05 else f"{v:.2f}"
