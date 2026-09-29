"""Compact, animated survey overview for the desktop workspace."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from PySide6.QtCore import (
    QEasingCurve,
    QPointF,
    QRectF,
    Qt,
    QVariantAnimation,
)
from PySide6.QtGui import QColor, QFont, QLinearGradient, QPainter, QPen
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

_COLORS = ["#3b82f6", "#14b8a6", "#a78bfa", "#f59e0b", "#ec4899", "#06b6d4"]


def _numeric(df, name):
    values = pd.to_numeric(
        df.get(name, pd.Series(np.nan, index=df.index)), errors="coerce"
    )
    return values.where(np.isfinite(values))


@dataclass
class _SurveySummary:
    """Row-aligned metadata that preserves coordinate pairs."""

    frame: pd.DataFrame
    positioned: pd.Series
    frequency: pd.Series
    elevation: pd.Series
    tipper: pd.Series
    lines: pd.Series

    @classmethod
    def from_frame(cls, frame):
        df = frame.copy().reset_index(drop=True)
        lat, lon = _numeric(df, "Latitude"), _numeric(df, "Longitude")
        df["Latitude"], df["Longitude"] = lat, lon
        positioned = lat.between(-90, 90) & lon.between(-180, 180)
        nf = _numeric(df, "N_freq")
        nf = nf.where((nf >= 0) & (nf % 1 == 0))
        tipper = df.get(
            "Tipper", pd.Series(None, index=df.index, dtype=object)
        )
        tipper = tipper.map(
            lambda v: {
                "true": True,
                "1": True,
                "1.0": True,
                "yes": True,
                "false": False,
                "0": False,
                "0.0": False,
                "no": False,
            }.get(str(v).strip().lower(), None)
        )
        lines = (
            df.get("Line", pd.Series("Unassigned", index=df.index))
            .fillna("Unassigned")
            .astype(str)
            .str.strip()
        )
        lines = lines.replace({"": "Unassigned", "—": "Unassigned"})
        return cls(
            df, positioned, nf, _numeric(df, "Elevation"), tipper, lines
        )


class _AnimatedNumber(QLabel):
    """Short eased count-up, settling on the exact value without looping."""

    def __init__(self, parent=None):
        super().__init__("—", parent)
        self.setObjectName("SurveyKpiValue")
        self._value = 0.0
        self._decimals = 0
        self._animation = QVariantAnimation(self)
        self._animation.setDuration(700)
        self._animation.setEasingCurve(QEasingCurve.Type.OutCubic)
        self._animation.valueChanged.connect(self._advance)

    def set_number(self, value):
        self._animation.stop()
        if value is None:
            self._value = 0.0
            self.setText("—")
            return
        self._target = float(value)
        self._decimals = int(not float(value).is_integer())
        self._animation.setStartValue(self._value)
        self._animation.setEndValue(self._target)
        self._animation.start()

    def _advance(self, value):
        self._value = float(value)
        self.setText(f"{self._value:,.{self._decimals}f}")

    def hideEvent(self, event):  # noqa: N802
        if self._animation.state() == QVariantAnimation.State.Running:
            self._animation.stop()
            self._advance(self._target)
        super().hideEvent(event)


class _SurveyCharts(QWidget):
    """A compact row of painter charts, with details available on hover."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(285)
        self.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        self.setMouseTracking(True)
        self.setAccessibleName("Survey overview charts")
        self._progress = 1.0
        self._hovered = -1
        self._cards = []
        self._points = []
        self._animation = QVariantAnimation(self)
        self._animation.setDuration(850)
        self._animation.setEasingCurve(QEasingCurve.Type.OutCubic)
        self._animation.valueChanged.connect(self._advance)
        self.set_dataframe(None)

    def _advance(self, value):
        self._progress = float(value)
        self.update()

    def set_dataframe(self, df):
        self._animation.stop()
        self._summary = _SurveySummary.from_frame(
            df if df is not None else pd.DataFrame()
        )
        s = self._summary
        self._total = len(s.frame)
        self._tipper = int(s.tipper.eq(True).sum())
        self._kinds = (["tipper"] if self._tipper else []) + [
            "frequency",
            "footprint",
        ]
        self._cards = []
        self._points = []
        self._hovered = -1
        nf = s.frequency.dropna()
        self._hist = np.array([])
        if len(nf):
            bins = min(12, max(3, int(np.sqrt(len(nf)))))
            if nf.nunique() == 1:
                bins = [nf.iloc[0] - 0.5, nf.iloc[0] + 0.5]
            self._hist, self._edges = np.histogram(nf, bins=bins)
        self._animation.setStartValue(0.0)
        self._animation.setEndValue(1.0)
        self._progress = 0.0 if self._total else 1.0
        self.setAccessibleDescription(
            f"{self._total} stations; {self._tipper} with tipper; "
            f"{int(s.positioned.sum())} with valid coordinate pairs."
        )
        if self.isVisible() and self._total:
            self._animation.start()
        self.update()

    def showEvent(self, event):  # noqa: N802
        super().showEvent(event)
        if self._total and self._progress < 1:
            self._animation.start()

    def hideEvent(self, event):  # noqa: N802
        self._animation.stop()
        self._progress = 1.0
        super().hideEvent(event)

    def _text(self, p, rect, text, size=9, bold=False, color=None, align=None):
        font = QFont(self.font())
        font.setPointSizeF(size)
        font.setBold(bold)
        p.setFont(font)
        p.setPen(color or self.palette().color(self.palette().ColorRole.Text))
        text = p.fontMetrics().elidedText(
            text, Qt.TextElideMode.ElideRight, max(1, int(rect.width()))
        )
        p.drawText(rect, align or Qt.AlignmentFlag.AlignCenter, text)

    def paintEvent(self, event):  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        pal = self.palette()
        dark = pal.color(pal.ColorRole.Window).lightness() < 128
        surface = QColor("#24273a" if dark else "#ffffff")
        border = QColor("#3c435b" if dark else "#dce5f1")
        muted = QColor("#a6adc8" if dark else "#73839a")
        count = len(self._kinds)
        gap, margin = 10.0, 1.0
        width = max(
            1.0, (self.width() - margin * 2 - gap * (count - 1)) / count
        )
        self._cards = []
        self._points = []
        titles = {
            "tipper": "Tipper coverage",
            "frequency": "Frequency sampling",
            "footprint": "Survey footprint",
        }
        for i, kind in enumerate(self._kinds):
            rect = QRectF(
                margin + i * (width + gap), 4, width, self.height() - 8
            )
            self._cards.append((kind, rect))
            p.setBrush(surface)
            p.setPen(
                QPen(QColor("#7aa7ed") if i == self._hovered else border, 1)
            )
            p.drawRoundedRect(rect, 12, 12)
            self._text(
                p,
                rect.adjusted(12, 12, -12, -rect.height() + 38),
                titles[kind],
                9,
                True,
                align=Qt.AlignmentFlag.AlignLeft
                | Qt.AlignmentFlag.AlignVCenter,
            )
            plot = rect.adjusted(18, 53, -18, -40)
            if plot.width() <= 2 or plot.height() <= 2:
                continue
            p.save()
            p.setClipRect(plot.adjusted(-8, -8, 8, 8))
            if kind == "tipper":
                footer = self._donut(p, plot, dark, muted)
            elif kind == "frequency":
                footer = self._frequency(p, plot, muted)
            else:
                footer = self._footprint(p, plot, muted)
            p.restore()
            self._text(
                p,
                QRectF(
                    rect.x() + 10, rect.bottom() - 32, rect.width() - 20, 22
                ),
                footer,
                8,
                color=muted,
            )
        p.end()

    def _donut(self, p, plot, dark, muted):
        size = max(8, min(plot.width(), plot.height()) * 0.77)
        circle = QRectF(
            plot.center().x() - size / 2,
            plot.center().y() - size / 2,
            size,
            size,
        )
        thickness = max(7, min(16, size * 0.12))
        p.setBrush(Qt.BrushStyle.NoBrush)
        pen = QPen(QColor("#363c52" if dark else "#edf2f8"), thickness)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        p.setPen(pen)
        p.drawArc(circle, 0, 360 * 16)
        percent = self._tipper / self._total
        gradient = QLinearGradient(circle.topLeft(), circle.bottomRight())
        gradient.setColorAt(0, QColor("#60a5fa"))
        gradient.setColorAt(1, QColor("#14b8a6"))
        pen.setBrush(gradient)
        p.setPen(pen)
        p.drawArc(circle, 90 * 16, -int(360 * 16 * percent * self._progress))
        self._text(
            p,
            circle,
            f"{percent * self._progress:.0%}",
            min(24, size * 0.2),
            True,
        )
        return f"{self._tipper} of {self._total} stations"

    def _frequency(self, p, plot, muted):
        if not len(self._hist):
            self._text(p, plot, "No frequency data", color=muted)
            return "—"
        bw = plot.width() / len(self._hist)
        peak = max(1, self._hist.max())
        for i, count in enumerate(self._hist):
            # Slight stagger gives a soft wave without a continuous animation.
            progress = min(
                1, max(0, self._progress * 1.15 - i / len(self._hist) * 0.15)
            )
            height = plot.height() * count / peak * progress
            bar = QRectF(
                plot.left() + i * bw + 2,
                plot.bottom() - height,
                max(1, bw - 4),
                height,
            )
            gradient = QLinearGradient(bar.topLeft(), bar.bottomLeft())
            gradient.setColorAt(0, QColor("#60a5fa"))
            gradient.setColorAt(1, QColor("#2563eb"))
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(gradient)
            p.drawRoundedRect(bar, min(4, bw / 4), min(4, height / 2))
        nf = self._summary.frequency.dropna()
        return f"{nf.min():g}–{nf.max():g} frequencies / station"

    def _footprint(self, p, plot, muted):
        s = self._summary
        valid = s.frame.loc[s.positioned]
        if valid.empty:
            self._text(p, plot, "No coordinates", color=muted)
            return "—"
        # Bound paint cost for large surveys; stats still include every row.
        points = valid.iloc[
            np.linspace(0, len(valid) - 1, min(len(valid), 2000), dtype=int)
        ]
        lon, lat = points.Longitude.to_numpy(), points.Latitude.to_numpy()
        if np.ptp(lon) > 180:
            lon = lon % 360
        x = (lon - lon.mean()) * max(0.15, np.cos(np.deg2rad(lat.mean())))
        y = lat - lat.mean()
        scale = (
            min(
                plot.width() / max(np.ptp(x), 1e-9),
                plot.height() / max(np.ptp(y), 1e-9),
            )
            * 0.86
        )
        xs = plot.center().x() + (x - (x.max() + x.min()) / 2) * scale
        ys = plot.center().y() - (y - (y.max() + y.min()) / 2) * scale
        lines = list(s.lines.unique())
        palette = {
            line: _COLORS[i % len(_COLORS)] for i, line in enumerate(lines)
        }
        for row, px, py in zip(points.index, xs, ys):
            color = QColor(palette[s.lines.loc[row]])
            color.setAlphaF(0.85 * self._progress)
            p.setBrush(color)
            p.setPen(Qt.PenStyle.NoPen)
            radius = 2.0 + 1.5 * self._progress
            point = QPointF(float(px), float(py))
            p.drawEllipse(point, radius, radius)
            self._points.append((point, row))
        return f"{len(valid):,} positioned stations"

    def mouseMoveEvent(self, event):  # noqa: N802
        position = event.position()
        hovered = next(
            (
                i
                for i, (_, rect) in enumerate(self._cards)
                if rect.contains(position)
            ),
            -1,
        )
        if hovered != self._hovered:
            self._hovered = hovered
            self.update()
        if hovered < 0:
            self.setToolTip("")
            return
        kind = self._cards[hovered][0]
        s = self._summary
        if kind == "tipper":
            message = (
                f"{self._tipper} of {self._total} stations with tipper\n"
                f"{s.tipper.isna().sum()} unknown"
            )
        elif kind == "frequency":
            nf = s.frequency.dropna()
            message = (
                (
                    f"Median: {nf.median():g} frequency samples per station\n"
                    f"{s.frequency.isna().sum()} unknown counts"
                )
                if len(nf)
                else "Frequency counts unavailable"
            )
        else:
            message = (
                f"{int(s.positioned.sum())} valid coordinate pairs\n"
                f"{int((~s.positioned).sum())} missing or invalid\n"
                "Colours group survey lines; "
                "approximate local geographic aspect."
            )
            nearest = min(
                self._points,
                key=lambda item: (item[0] - position).manhattanLength(),
                default=None,
            )
            if nearest and (nearest[0] - position).manhattanLength() < 12:
                row = nearest[1]
                record = s.frame.loc[row]
                message = (
                    f"{record.get('ID', 'Station')} · {s.lines.loc[row]}\n"
                    f"{record.Latitude:.5f}°, {record.Longitude:.5f}°"
                )
        self.setToolTip(message)
        super().mouseMoveEvent(event)

    def leaveEvent(self, event):  # noqa: N802
        self._hovered = -1
        self.setToolTip("")
        self.update()
        super().leaveEvent(event)


class SurveyOverviewWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("SurveyOverview")
        self._build_ui()
        self.clear()

    def _kpi(self, title):
        card = QFrame()
        card.setObjectName("SurveyKpiCard")
        layout = QVBoxLayout(card)
        layout.setContentsMargins(12, 8, 12, 8)
        label = QLabel(title)
        label.setObjectName("SurveyKpiLabel")
        label.setWordWrap(True)
        value = _AnimatedNumber(card)
        sub = QLabel("—")
        sub.setObjectName("SurveyKpiSubtext")
        sub.setWordWrap(True)
        for widget in (label, value, sub):
            layout.addWidget(widget)
        return card, value, sub

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)
        self._placeholder = QLabel(
            "Open EDI or EMTF XML files to view survey analytics"
        )
        self._placeholder.setObjectName("SurveyPlaceholder")
        self._placeholder.setWordWrap(True)
        self._placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        root.addWidget(self._placeholder, 1)
        self._scroll = QScrollArea(self)
        self._scroll.setObjectName("SurveyDashboardScroll")
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._dashboard = QWidget()
        self._dashboard.setObjectName("SurveyDashboardBody")
        self._dashboard.setMinimumHeight(420)
        dash = QVBoxLayout(self._dashboard)
        dash.setContentsMargins(0, 0, 0, 0)
        dash.setSpacing(10)
        kpis = QHBoxLayout()
        kpis.setSpacing(10)
        self._kpi_cards = {}
        for name, title in [
            ("stations", "STATIONS"),
            ("nfreq", "FREQUENCIES"),
            ("elev", "ELEVATION"),
            ("tipper", "TIPPER COVERAGE"),
        ]:
            card, value, sub = self._kpi(title)
            setattr(self, "_v_" + name, value)
            setattr(self, "_s_" + name, sub)
            self._kpi_cards[name] = card
            kpis.addWidget(card, 1)
        dash.addLayout(kpis)
        self._charts = _SurveyCharts(self._dashboard)
        dash.addWidget(self._charts, 1)
        self._meta = QLabel()
        self._meta.setTextFormat(Qt.TextFormat.PlainText)
        self._meta.setObjectName("SurveyMeta")
        self._meta.setWordWrap(True)
        dash.addWidget(self._meta)
        self._scroll.setWidget(self._dashboard)
        root.addWidget(self._scroll, 1)

    def update_survey(self, df, paths=None):
        if df is None or df.empty:
            self.clear()
            return
        s = _SurveySummary.from_frame(df)
        n = len(s.frame)
        nf, elev = s.frequency.dropna(), s.elevation.dropna()
        tipper = int(s.tipper.eq(True).sum())
        self._placeholder.hide()
        self._scroll.show()
        self._dashboard.show()
        self._v_stations.set_number(n)
        self._s_stations.setText("survey sites loaded")
        self._v_nfreq.set_number(float(nf.median()) if len(nf) else None)
        self._s_nfreq.setText(
            f"median · {nf.min():g}–{nf.max():g}"
            if len(nf)
            else "no frequency data"
        )
        self._v_elev.set_number(
            round(float(elev.max() - elev.min())) if len(elev) else None
        )
        self._s_elev.setText(
            f"metres relief · {elev.min():.0f}–{elev.max():.0f} m"
            if len(elev)
            else "no elevation data"
        )
        self._kpi_cards["tipper"].setVisible(tipper > 0)
        self._v_tipper.set_number(round(100 * tipper / n))
        self._s_tipper.setText(f"percent · {tipper} of {n} stations")
        self._charts.set_dataframe(df)
        formats = sorted(
            {
                Path(p).suffix.lstrip(".").upper()
                for p in paths or []
                if Path(p).suffix
            }
        )
        directories = sorted({str(Path(p).parent) for p in paths or []})
        source = (
            directories[0]
            if len(directories) == 1
            else f"{len(directories)} source directories"
            if directories
            else "—"
        )
        if len(directories) == 1 and len(Path(source).parts) > 4:
            source = str(Path("…", *Path(source).parts[-3:]))
        self._meta.setText(
            f"FORMAT  {' / '.join(formats) or '—'}     ·     SOURCE  {source}"
        )

    def clear(self):
        self._scroll.hide()
        self._dashboard.hide()
        self._placeholder.show()
        self._kpi_cards["tipper"].hide()
        self._charts.set_dataframe(None)
