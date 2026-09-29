# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Run monitor: what an inversion is doing right now.

One compact strip under the stepper:

``[Running]  Station 3/10 · S07   Iter 7 / 30   RMS 1.84 → 1.00
 ▁▂▃▅ (live RMS curve)   02:13 elapsed · ~05:40 left   [■ Stop]``

Iterations come either from the engine itself (Occam1D callback) or from
the solver's log, polled by the window while an external binary runs.
Progress is iteration / max-iterations; for Occam1D it is spread over the
stations (station k of n contributes 1/n).
"""

from __future__ import annotations

import math
import time

from PySide6.QtCore import QPointF, QRectF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

# Same filled, ≥ 4.5:1 pills as Pipeline Studio
STATES = {
    "idle": ("Idle", "#5f6b7a"),
    "building": ("Building", "#3b5bdb"),
    "ready": ("Ready", "#3b5bdb"),
    "running": ("Running", "#1864ab"),
    "done": ("Done", "#2a7f3f"),
    "error": ("Error", "#c92a2a"),
    "stopped": ("Stopped", "#9c4f00"),
}


def fmt_seconds(s: float | None) -> str:
    if s is None or not math.isfinite(s):
        return "–"
    s = int(round(s))
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    return f"{h}:{m:02d}:{sec:02d}" if h else f"{m:02d}:{sec:02d}"


class StatePill(QLabel):
    def __init__(self) -> None:
        super().__init__()
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setMinimumWidth(72)
        self.set_state("idle")

    def set_state(self, state: str) -> None:
        self.state = state
        text, colour = STATES.get(state, STATES["idle"])
        self.setText(text)
        self.setStyleSheet(
            f"QLabel {{ background: {colour}; color: white; "
            "border-radius: 8px; padding: 1px 10px; font-weight: 600; "
            "font-size: 11px; }")


class Sparkline(QWidget):
    """Tiny RMS-vs-iteration curve (log y) with the target as a dash."""

    def __init__(self) -> None:
        super().__init__()
        self.values: list[float] = []
        self.target: float | None = None
        self.setMinimumSize(140, 34)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Fixed)
        self.setToolTip("RMS misfit per iteration (log scale); dashed line "
                        "= target")

    def set_data(self, values, target=None) -> None:
        self.values = [float(v) for v in values if v and math.isfinite(v)
                       and v > 0]
        self.target = target
        self.update()

    def paintEvent(self, _event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        r = QRectF(self.rect()).adjusted(3, 4, -3, -4)
        p.setPen(QPen(QColor(128, 128, 128, 90), 1))
        p.drawRoundedRect(QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5),
                          4, 4)
        vals = self.values
        if not vals:
            p.setPen(QColor(128, 128, 128))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, "no iterations yet")
            return
        logs = [math.log10(v) for v in vals]
        all_logs = logs + ([math.log10(self.target)]
                           if self.target and self.target > 0 else [])
        lo, hi = min(all_logs), max(all_logs)
        if hi - lo < 1e-6:
            lo, hi = lo - 0.5, hi + 0.5

        def pt(i, lv):
            x = r.left() + (r.width() * i / max(len(logs) - 1, 1))
            y = r.bottom() - (lv - lo) / (hi - lo) * r.height()
            return QPointF(x, y)

        if self.target and self.target > 0:
            y = pt(0, math.log10(self.target)).y()
            p.setPen(QPen(QColor("#2a7f3f"), 1, Qt.PenStyle.DashLine))
            p.drawLine(QPointF(r.left(), y), QPointF(r.right(), y))
        path = QPainterPath(pt(0, logs[0]))
        for i, lv in enumerate(logs[1:], 1):
            path.lineTo(pt(i, lv))
        p.setPen(QPen(QColor("#1c7ed6"), 1.8))
        p.drawPath(path)
        p.setBrush(QColor("#1c7ed6"))
        p.drawEllipse(pt(len(logs) - 1, logs[-1]), 2.5, 2.5)


class RunMonitor(QWidget):
    """Live status strip for one run."""

    stop_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("RunMonitor")
        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 4, 6, 4)
        outer.setSpacing(3)
        top = QHBoxLayout()
        top.setSpacing(10)
        self.pill = StatePill()
        top.addWidget(self.pill)
        self.stage = QLabel("No run yet")
        self.stage.setStyleSheet("font-weight: 600;")
        self.stage.setSizePolicy(QSizePolicy.Policy.Ignored,
                                 QSizePolicy.Policy.Preferred)
        top.addWidget(self.stage, 2)
        self.iter_lbl = QLabel("Iter –")
        self.rms_lbl = QLabel("RMS –")
        self.time_lbl = QLabel("")
        self.time_lbl.setObjectName("InfoLabel")
        for w in (self.iter_lbl, self.rms_lbl, self.time_lbl):
            top.addWidget(w)
        self.btn_stop = QPushButton("■ Stop")
        self.btn_stop.setToolTip("Stop the solver (the process tree is "
                                 "terminated; files written so far stay)")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.stop_requested)
        top.addWidget(self.btn_stop)
        outer.addLayout(top)
        low = QHBoxLayout()
        low.setSpacing(8)
        self.bar = QProgressBar()
        self.bar.setRange(0, 1000)
        self.bar.setTextVisible(False)
        self.bar.setFixedHeight(8)
        low.addWidget(self.bar, 3)
        self.spark = Sparkline()
        low.addWidget(self.spark, 2)
        outer.addLayout(low)

        self._t0: float | None = None
        self._t_end: float | None = None
        self.max_iter = 0
        self.target: float | None = None
        self.n_units = 1  # stations for Occam1D
        self.unit = 0
        self.history: list[tuple[str, int, float]] = []
        self._clock = QTimer(self)
        self._clock.setInterval(1000)
        self._clock.timeout.connect(self._tick)

    # ── lifecycle ─────────────────────────────────────────────────────
    def start(self, stage: str, *, max_iter: int, target: float | None,
              n_units: int = 1) -> None:
        self.history = []
        self.max_iter = max(int(max_iter), 0)
        self.target = target
        self.n_units = max(int(n_units), 1)
        self.unit = 0
        self._t0 = time.monotonic()
        self._t_end = None
        self.pill.set_state("running")
        self.stage.setText(stage)
        self.iter_lbl.setText(f"Iter 0 / {self.max_iter}"
                              if self.max_iter else "Iter –")
        self.rms_lbl.setText("RMS –" + (f" → {target:g}" if target else ""))
        self.bar.setRange(0, 1000)
        self.bar.setValue(0)
        self.spark.set_data([], target)
        self.btn_stop.setEnabled(True)
        self._clock.start()
        self._tick()

    def set_busy(self, stage: str, state: str = "building") -> None:
        """Indeterminate activity (building inputs, loading results)."""
        self.pill.set_state(state)
        self.stage.setText(stage)
        self.bar.setRange(0, 0)

    def set_stage(self, text: str) -> None:
        self.stage.setText(text)
        m = _station_index(text)
        if m is not None:
            self.unit, self.n_units = m

    def add_iteration(self, label: str, n: int, rms: float) -> None:
        self.history.append((label, int(n), float(rms)))
        cur = [h for h in self.history if h[0] == label]
        self.spark.set_data([h[2] for h in cur], self.target)
        its = max(n, 0)
        self.iter_lbl.setText(f"Iter {its} / {self.max_iter}"
                              if self.max_iter else f"Iter {its}")
        self.rms_lbl.setText(f"RMS {rms:.3f}"
                             + (f" → {self.target:g}" if self.target else ""))
        if self.max_iter:
            frac = (min(its, self.max_iter) / self.max_iter)
            unit = max(self.unit - 1, 0) if self.n_units > 1 else 0
            total = (unit + frac) / self.n_units
            self.bar.setRange(0, 1000)
            self.bar.setValue(int(1000 * min(total, 1.0)))
        self._tick()

    def set_history(self, pairs: list[tuple[int, float]]) -> None:
        """Iterations parsed from a solver log (replaces the history)."""
        known = {(h[1], h[2]) for h in self.history}
        for n, rms in pairs:
            if (n, rms) not in known:
                self.add_iteration("", n, rms)

    def finish(self, state: str, text: str) -> None:
        self._t_end = time.monotonic()
        self._clock.stop()
        self._tick()
        self.pill.set_state(state)
        self.stage.setText(text)
        self.bar.setRange(0, 1000)
        if state == "done":
            self.bar.setValue(1000)
        self.btn_stop.setEnabled(False)

    def reset(self, text: str = "No run yet") -> None:
        self._clock.stop()
        self.pill.set_state("idle")
        self.stage.setText(text)
        self.iter_lbl.setText("Iter –")
        self.rms_lbl.setText("RMS –")
        self.time_lbl.setText("")
        self.bar.setRange(0, 1000)
        self.bar.setValue(0)
        self.spark.set_data([])
        self.btn_stop.setEnabled(False)

    # ── derived ───────────────────────────────────────────────────────
    @property
    def running(self) -> bool:
        return self.pill.state == "running"

    def percent(self) -> int:
        if self.bar.maximum() == 0:
            return -1
        return int(round(100 * self.bar.value() / 1000))

    def eta(self) -> float | None:
        if self._t0 is None or not self.max_iter:
            return None
        frac = self.bar.value() / 1000
        if frac <= 0.02 or frac >= 1:
            return None
        elapsed = (self._t_end or time.monotonic()) - self._t0
        return elapsed * (1 - frac) / frac

    def _tick(self) -> None:
        if self._t0 is None:
            return
        elapsed = (self._t_end or time.monotonic()) - self._t0
        txt = f"{fmt_seconds(elapsed)} elapsed"
        eta = self.eta() if self._t_end is None else None
        if eta is not None:
            txt += f" · ~{fmt_seconds(eta)} left"
        self.time_lbl.setText(txt)


def _station_index(text: str):
    """(k, n) from "Station k/n · NAME", else None."""
    import re

    m = re.match(r"\s*Station\s+(\d+)\s*/\s*(\d+)", text)
    return (int(m.group(1)), int(m.group(2))) if m else None


__all__ = ["STATES", "RunMonitor", "Sparkline", "StatePill", "fmt_seconds"]
