# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
StationDetailCard — the main window's "Station Statistics" tab.

    ┌──────────────────────────────────────────────┐
    │  18-023A                          [ B ]      │  grade badge (A–E)
    │  Line L18 · 53 frequencies · tipper          │
    │  ── Quality ──────────────────────────────── │
    │  Completeness    100 %  ██████████▏          │  bar + survey median ▏
    │  Signal / noise   39    ███████▌  ▏          │
    │  2-D consistency  55 %  █████▌   ▏           │
    │  ── Response ─────────────────────────────── │
    │  ρa  ╲___╱‾‾   (XY ●  YX ●)                  │  animated preview
    │  φ   ‾‾╲__                                   │
    │  ── Frequency coverage ───────────────────── │
    │  ▕▕▕▕ ▕▕▕▕▕▕▕  ▕▕   ticks coloured by error,  │
    │  0.001 Hz ─────── 1000 Hz  survey band, gaps │
    │  ── Location ─────────────────────────────── │
    │  Latitude / Longitude / Elevation / Period…  │
    ├──────────────────────────────────────────────┤
    │  [Open Profile]        [Show on Map]         │
    └──────────────────────────────────────────────┘

The grade and bars are measured (see
:mod:`pycsamt.app.desktop.controllers.station_stats`); they replace the
old "quality dots", which were ``n_frequencies // 10``.  On a new station
the curves and the coverage ticks draw left to right (~0.35 s) unless
animations are off (Preferences ▸ General).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from PySide6.QtCore import (
    Property,
    QEasingCurve,
    QPointF,
    QPropertyAnimation,
    QRectF,
    Qt,
    Signal,
)
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.station_stats import (
    GRADE_COLOURS,
    StationStats,
    station_stats,
    survey_medians,
)

_ICONS = Path(__file__).parent.parent / "resources" / "icons"
_XY = QColor("#d6336c")
_YX = QColor("#1c7ed6")


def _lbl(text: str, obj_name: str = "") -> QLabel:
    lbl = QLabel(text)
    if obj_name:
        lbl.setObjectName(obj_name)
    return lbl


def _section(title: str) -> QLabel:
    lbl = QLabel(title.upper())
    lbl.setObjectName("InfoLabel")
    lbl.setStyleSheet("font-weight: 600; letter-spacing: 1px; "
                      "padding-top: 6px;")
    return lbl


_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _pow10(e: int) -> str:
    return "10" + str(int(e)).translate(_SUP)


def _err_colour(rel: float) -> QColor:
    """Tick colour for a relative impedance error."""
    if not math.isfinite(rel):
        return QColor("#3b82f6")  # no error data: neutral blue
    if rel <= 0.05:
        return QColor("#2a7f3f")
    if rel <= 0.15:
        return QColor("#d08700")
    return QColor("#c92a2a")


class _Animated(QWidget):
    """Base for widgets revealed left to right (``progress`` 0 → 1)."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._progress = 1.0
        self._anim = QPropertyAnimation(self, b"progress", self)
        self._anim.setDuration(350)
        self._anim.setEasingCurve(QEasingCurve.Type.OutCubic)

    def _get_progress(self) -> float:
        return self._progress

    def _set_progress(self, v: float) -> None:
        self._progress = float(v)
        self.update()

    progress = Property(float, _get_progress, _set_progress)

    def reveal(self, animate: bool) -> None:
        self._anim.stop()
        if animate:
            self._anim.setStartValue(0.0)
            self._anim.setEndValue(1.0)
            self._anim.start()
        else:
            self._set_progress(1.0)


class FrequencyCoverageGraphic(_Animated):
    """Log-frequency ticks coloured by relative error, the survey's band
    behind them, and highlighted gaps (spacing > 2.5 x the median)."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._values: list[float] = []
        self._errors: list[float] = []
        self._band: tuple[float, float] | None = None
        self.setMinimumHeight(86)
        self.setMaximumHeight(100)
        self.setToolTip("Frequency sampling (log scale). Tick colour = "
                        "relative impedance error: green ≤ 5 %, amber ≤ "
                        "15 %, red > 15 %, blue = no errors. Grey band = "
                        "the survey's range; shaded = gaps.")

    def set_frequencies(self, frequencies, rel_err=None, band=None) -> None:
        values = np.asarray(frequencies if frequencies is not None else [],
                            dtype=float)
        good = np.isfinite(values) & (values > 0)
        errs = (np.asarray(rel_err, float)[good] if rel_err is not None
                and len(rel_err) == len(values) else
                np.full(good.sum(), np.nan))
        order = np.argsort(values[good])
        self._values = np.log10(values[good][order]).tolist()
        self._errors = errs[order].tolist()
        self._band = band
        self.update()

    def gaps(self) -> list[tuple[float, float]]:
        v = np.asarray(self._values)
        if v.size < 3:
            return []
        d = np.diff(v)
        med = float(np.median(d))
        return [(float(v[i]), float(v[i + 1])) for i in np.flatnonzero(
            d > 2.5 * max(med, 1e-9))]

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        pal = self.palette()
        fg = pal.color(pal.ColorRole.Text)
        muted = pal.color(pal.ColorRole.Mid)
        area = QRectF(8, 20, max(1, self.width() - 16), self.height() - 40)
        y = area.center().y()
        if not self._values:
            p.setPen(fg)
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                       "No frequency data")
            return
        lo, hi = self._values[0], self._values[-1]
        if self._band is not None and all(
                math.isfinite(b) and b > 0 for b in self._band):
            lo = min(lo, math.log10(self._band[0]))
            hi = max(hi, math.log10(self._band[1]))
        span = max(hi - lo, 1e-6)

        def X(v):
            return area.left() + (v - lo) / span * area.width()

        if self._band is not None and all(
                math.isfinite(b) and b > 0 for b in self._band):
            x0, x1 = X(math.log10(self._band[0])), X(math.log10(self._band[1]))
            p.fillRect(QRectF(x0, y - 16, x1 - x0, 32),
                       QColor(128, 128, 128, 38))
        for a, b in self.gaps():
            p.fillRect(QRectF(X(a), y - 16, X(b) - X(a), 32),
                       QColor(201, 42, 42, 40))
        p.setPen(QPen(muted, 1))
        p.drawLine(QPointF(area.left(), y), QPointF(area.right(), y))
        cut = area.left() + self._progress * area.width()
        for v, e in zip(self._values, self._errors):
            x = X(v)
            if x > cut + 0.5:
                break
            p.setPen(QPen(_err_colour(e), 2.2))
            p.drawLine(QPointF(x, y - 13), QPointF(x, y + 13))
        p.setPen(fg)
        p.drawText(QRectF(area.left(), 0, area.width(), 18),
                   Qt.AlignmentFlag.AlignCenter,
                   f"{len(self._values)} frequencies  ·  "
                   f"{10 ** self._values[0]:.3g}–{10 ** self._values[-1]:.3g}"
                   " Hz" + (f"  ·  {len(self.gaps())} gap(s)"
                            if self.gaps() else ""))
        # axis ends (the survey band may extend beyond the station)
        p.drawText(QRectF(area.left(), area.bottom() + 4, area.width(), 16),
                   Qt.AlignmentFlag.AlignLeft, f"{10 ** lo:.3g} Hz")
        p.drawText(QRectF(area.left(), area.bottom() + 4, area.width(), 16),
                   Qt.AlignmentFlag.AlignRight, f"{10 ** hi:.3g} Hz")


class ResponsePreview(_Animated):
    """Apparent resistivity (log-log) and phase against period, XY and
    YX, drawn with QPainter so it is instant and theme-aware."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._st: StationStats | None = None
        self.setMinimumHeight(170)
        self.setMaximumHeight(210)
        self.setToolTip("Apparent resistivity (top, log scale) and phase "
                        "(bottom) against period. Magenta = XY, blue = YX.")

    def set_stats(self, st: StationStats | None) -> None:
        self._st = st
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        pal = self.palette()
        fg = pal.color(pal.ColorRole.Text)
        grid = pal.color(pal.ColorRole.Mid)
        st = self._st
        if st is None or st.period.size < 2:
            p.setPen(fg)
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                       "No response to preview")
            return
        left, right = 44.0, 8.0
        w = max(10.0, self.width() - left - right)
        h = self.height()
        top_r = QRectF(left, 8, w, h * 0.55 - 12)
        top_p = QRectF(left, h * 0.55 + 4, w, h * 0.45 - 22)
        lt = np.log10(st.period)
        t0, t1 = float(np.nanmin(lt)), float(np.nanmax(lt))
        tspan = max(t1 - t0, 1e-6)
        rhos = np.concatenate([st.rho_xy, st.rho_yx])
        rhos = rhos[np.isfinite(rhos) & (rhos > 0)]
        if not rhos.size:
            p.setPen(fg)
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                       "No finite impedances")
            return
        r0 = math.floor(np.log10(rhos.min()))
        r1 = math.ceil(np.log10(rhos.max()))
        r1 = r1 if r1 > r0 else r0 + 1
        cut = self._progress

        def xpos(v):
            return left + (v - t0) / tspan * w

        for rect in (top_r, top_p):
            p.setPen(QPen(grid, 1))
            p.drawRect(rect)
        p.setPen(fg)
        f = p.font()
        f.setPointSizeF(max(6.5, f.pointSizeF() - 1.5))
        p.setFont(f)
        p.drawText(QRectF(0, top_r.top(), left - 4, 14),
                   Qt.AlignmentFlag.AlignRight, _pow10(r1))
        p.drawText(QRectF(0, top_r.bottom() - 14, left - 4, 14),
                   Qt.AlignmentFlag.AlignRight, _pow10(r0))
        # colour key
        kx = top_r.right() - 70
        for i, (lab, col) in enumerate((("XY", _XY), ("YX", _YX))):
            if lab == "YX" and st.scalar:
                continue
            p.setPen(QPen(col, 2.2))
            p.drawLine(QPointF(kx + i * 36, top_r.top() + 9),
                       QPointF(kx + i * 36 + 12, top_r.top() + 9))
            p.setPen(fg)
            p.drawText(QRectF(kx + i * 36 + 15, top_r.top() + 2, 22, 14),
                       Qt.AlignmentFlag.AlignLeft, lab)
        p.drawText(QRectF(0, top_r.center().y() - 7, left - 4, 14),
                   Qt.AlignmentFlag.AlignRight, "ρa")
        p.drawText(QRectF(0, top_p.top(), left - 4, 14),
                   Qt.AlignmentFlag.AlignRight, "90°")
        p.drawText(QRectF(0, top_p.bottom() - 14, left - 4, 14),
                   Qt.AlignmentFlag.AlignRight, "0°")
        p.drawText(QRectF(left, h - 16, w, 14), Qt.AlignmentFlag.AlignLeft,
                   f"{10 ** t0:.3g} s")
        p.drawText(QRectF(left, h - 16, w, 14), Qt.AlignmentFlag.AlignRight,
                   f"{10 ** t1:.3g} s")

        def curve(xv, yv, rect, lo, hi, colour, logy):
            path = QPainterPath()
            started = False
            # Reveal from short to long period; periods come high -> low
            # (frequencies are ascending), so sort, and allow for rounding
            # (the last point used to be dropped -> nothing drawn).
            limit = t0 + cut * tspan + 1e-9 * max(tspan, 1.0)
            order = np.argsort(xv)
            xv, yv = np.asarray(xv)[order], np.asarray(yv)[order]
            for x, y in zip(xv, yv):
                if not (np.isfinite(y) and (y > 0 or not logy)):
                    started = False
                    continue
                if x > limit:
                    break
                yy = math.log10(y) if logy else y
                py = rect.bottom() - (yy - lo) / (hi - lo) * rect.height()
                pt = QPointF(xpos(x), min(max(py, rect.top()),
                                          rect.bottom()))
                if started:
                    path.lineTo(pt)
                else:
                    path.moveTo(pt)
                    started = True
            p.setPen(QPen(colour, 1.8))
            p.drawPath(path)

        curve(lt, st.rho_xy, top_r, r0, r1, _XY, True)
        curve(lt, st.phi_xy, top_p, 0.0, 90.0, _XY, False)
        if not st.scalar:
            curve(lt, st.rho_yx, top_r, r0, r1, _YX, True)
            curve(lt, st.phi_yx, top_p, 0.0, 90.0, _YX, False)


class _MetricBar(QWidget):
    """Score bar (0–1) with an optional survey-median marker."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.score = 0.0
        self.survey: float | None = None
        self.setMinimumSize(80, 12)
        self.setMaximumHeight(14)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Fixed)

    def set_values(self, score: float, survey: float | None) -> None:
        self.score = float(np.clip(score, 0.0, 1.0)) \
            if math.isfinite(score) else 0.0
        self.survey = survey if survey is not None and math.isfinite(
            survey) else None
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        r = QRectF(0, 2, self.width() - 1, self.height() - 4)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(128, 128, 128, 55))
        p.drawRoundedRect(r, 4, 4)
        colour = QColor("#2a7f3f" if self.score >= 0.7 else
                        "#d08700" if self.score >= 0.45 else "#c92a2a")
        p.setBrush(colour)
        p.drawRoundedRect(QRectF(r.left(), r.top(), r.width() * self.score,
                                 r.height()), 4, 4)
        if self.survey is not None:
            x = r.left() + r.width() * float(np.clip(self.survey, 0, 1))
            p.setPen(QPen(self.palette().color(
                self.palette().ColorRole.Text), 2))
            p.drawLine(QPointF(x, 0), QPointF(x, self.height()))


class StationDetailCard(QWidget):
    """
    Measured statistics for the station selected in the main table.

    Signals
    -------
    open_profile_requested(station_id)
    show_on_map_requested(station_id)
    """

    open_profile_requested = Signal(str)
    show_on_map_requested = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._station_id: str = ""
        self._animate = True
        self._survey_key = None
        self._survey: dict = {}
        self.stats: StationStats | None = None
        self._build_ui()
        self.clear()

    # ── construction ──────────────────────────────────────────────
    def _build_ui(self) -> None:
        self.setObjectName("StationDetailCard")
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = QScrollArea(self)
        scroll.setObjectName("DetailScrollArea")
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        body = QWidget()
        body.setObjectName("DetailBody")
        v = QVBoxLayout(body)
        v.setContentsMargins(12, 12, 12, 8)
        v.setSpacing(5)

        head = QHBoxLayout()
        name_box = QVBoxLayout()
        name_box.setSpacing(0)
        self._lbl_name = QLabel("No station selected")
        self._lbl_name.setObjectName("StationName")
        self._lbl_name.setWordWrap(True)
        name_box.addWidget(self._lbl_name)
        self._lbl_sub = _lbl("", "InfoLabel")
        name_box.addWidget(self._lbl_sub)
        head.addLayout(name_box, 1)
        self._badge = QLabel("–")
        self._badge.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._badge.setFixedSize(44, 44)
        head.addWidget(self._badge, 0, Qt.AlignmentFlag.AlignTop)
        v.addLayout(head)

        v.addWidget(_section("Quality"))
        qgrid = QGridLayout()
        qgrid.setHorizontalSpacing(8)
        qgrid.setVerticalSpacing(3)
        self._metric_rows = []
        for i, name in enumerate(("Completeness", "Signal / noise",
                                  "2-D consistency")):
            k = _lbl(name, "DetailKey")
            val = _lbl("—", "DetailValue")
            bar = _MetricBar()
            qgrid.addWidget(k, i, 0)
            qgrid.addWidget(val, i, 1, Qt.AlignmentFlag.AlignRight)
            qgrid.addWidget(bar, i, 2)
            self._metric_rows.append((val, bar))
        qgrid.setColumnStretch(2, 1)
        v.addLayout(qgrid)
        self._lbl_vs = _lbl("", "InfoLabel")
        self._lbl_vs.setWordWrap(True)
        v.addWidget(self._lbl_vs)

        v.addWidget(_section("Response"))
        self._response = ResponsePreview(body)
        v.addWidget(self._response)
        v.addWidget(_section("Frequency coverage"))
        self._frequency_graph = FrequencyCoverageGraphic(body)
        v.addWidget(self._frequency_graph)

        v.addWidget(_section("Location & data"))
        grid = QGridLayout()
        grid.setSpacing(4)
        grid.setColumnMinimumWidth(0, 100)

        def row(r: int, key: str) -> QLabel:
            val = _lbl("—", "DetailValue")
            val.setTextInteractionFlags(
                Qt.TextInteractionFlag.TextSelectableByMouse)
            val.setWordWrap(True)
            grid.addWidget(_lbl(key, "DetailKey"), r, 0,
                           Qt.AlignmentFlag.AlignLeft)
            grid.addWidget(val, r, 1, Qt.AlignmentFlag.AlignLeft)
            return val

        self._lbl_lat = row(0, "Latitude")
        self._lbl_lon = row(1, "Longitude")
        self._lbl_elev = row(2, "Elevation")
        self._lbl_nfreq = row(3, "Frequencies")
        self._lbl_frange = row(4, "Period range")
        self._lbl_tipper = row(5, "Tipper")
        self._lbl_comps = row(6, "Components")
        v.addLayout(grid)
        # Back-compat: the grade replaces the old dots label.
        self._lbl_quality = self._badge
        v.addStretch(1)
        scroll.setWidget(body)
        outer.addWidget(scroll, stretch=1)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setObjectName("Separator")
        outer.addWidget(sep)
        btn_row = QHBoxLayout()
        btn_row.setContentsMargins(12, 6, 12, 8)
        btn_row.setSpacing(6)
        self._btn_profile = QPushButton("Open Profile")
        self._btn_map = QPushButton("Show on Map")
        for b in (self._btn_profile, self._btn_map):
            b.setObjectName("CardButton")
            b.setSizePolicy(QSizePolicy.Policy.Expanding,
                            QSizePolicy.Policy.Fixed)
            btn_row.addWidget(b)
        self._btn_profile.clicked.connect(
            lambda: self.open_profile_requested.emit(self._station_id))
        self._btn_map.clicked.connect(
            lambda: self.show_on_map_requested.emit(self._station_id))
        outer.addLayout(btn_row)

    # ── public API ────────────────────────────────────────────────
    def set_animations(self, on: bool) -> None:
        """Preferences ▸ General ▸ "Animate station statistics"."""
        self._animate = bool(on)

    def _set_badge(self, grade: str, tip: str = "") -> None:
        colour = GRADE_COLOURS.get(grade, GRADE_COLOURS["–"])
        self._badge.setText(grade)
        self._badge.setToolTip(tip)
        self._badge.setStyleSheet(
            f"QLabel {{ background: {colour}; color: white; "
            "border-radius: 22px; font-size: 20px; font-weight: 700; }")

    def _survey_stats(self, sites) -> dict:
        key = (id(sites), len(sites) if hasattr(sites, "__len__") else 0)
        if key != self._survey_key:
            try:
                self._survey = survey_medians(sites)
            except Exception:
                self._survey = {}
            self._survey_key = key
        return self._survey

    def update_station(self, station_id: str, sites) -> None:
        """Populate the card from a Sites collection for *station_id*."""
        self._station_id = station_id
        try:
            site = sites.get(station_id)
            s = site.summary()
        except Exception:
            self._lbl_name.setText(station_id)
            return
        lat, lon, elev = s.get("lat"), s.get("lon"), s.get("elev")
        nfreq = s.get("nfreq", 0)
        has_tip = bool(s.get("tipper", False))
        comps = s.get("components", [])
        self._lbl_name.setText(s.get("name", station_id))
        self._lbl_lat.setText(f"{lat:+.4f} °" if lat is not None else "—")
        self._lbl_lon.setText(f"{lon:+.4f} °" if lon is not None else "—")
        self._lbl_elev.setText(f"{elev:.0f} m" if elev is not None else "—")
        self._lbl_nfreq.setText(str(nfreq))
        self._lbl_tipper.setText("Yes" if has_tip else "No")
        self._lbl_comps.setText("  ".join(comps) if comps else "—")

        st = station_stats(site)
        self.stats = st
        survey = self._survey_stats(sites)
        if st.period.size:
            self._lbl_frange.setText(
                f"{st.period.min():.2e} – {st.period.max():.2e} s")
        else:
            self._lbl_frange.setText("—")
        self._lbl_sub.setText(
            f"{nfreq} frequencies · {'tipper' if has_tip else 'no tipper'}"
            + (" · scalar (Zxy only)" if st.scalar else ""))

        # quality rows + survey medians
        from pycsamt.app.desktop.controllers.station_stats import _snr_score

        medians = [survey.get("completeness"),
                   _snr_score(survey.get("snr", float("nan"))),
                   survey.get("skew_ok")]
        tip_lines = []
        for (val, bar), (label, text, score), med in zip(
                self._metric_rows, st.breakdown(), medians):
            val.setText(text)
            bar.set_values(score, med)
            tip_lines.append(f"{label}: {text}")
        self._set_badge(st.grade, "Grade " + st.grade + "\n"
                        + "\n".join(tip_lines)
                        + "\n(completeness, SNR and phase-tensor "
                          "|β| ≤ 3° averaged)")
        if survey:
            better = (st.score >= survey.get("score", float("nan")))
            self._lbl_vs.setText(
                f"▏ = survey median ({survey['n_stations']} stations). "
                + ("Above" if better else "Below") + " the survey's "
                "median quality.")
        else:
            self._lbl_vs.setText("")

        band = (survey.get("fmin"), survey.get("fmax")) if survey else None
        self._frequency_graph.set_frequencies(st.freq, st.rel_err, band)
        self._response.set_stats(st)
        self._frequency_graph.reveal(self._animate)
        self._response.reveal(self._animate)
        self._btn_profile.setEnabled(True)
        self._btn_map.setEnabled(True)

    def clear(self) -> None:
        """Reset to empty / no-selection state."""
        self._station_id = ""
        self.stats = None
        self._lbl_name.setText("No station selected")
        self._lbl_sub.setText("Select a station in the table")
        for lbl in (self._lbl_lat, self._lbl_lon, self._lbl_elev,
                    self._lbl_nfreq, self._lbl_frange, self._lbl_tipper,
                    self._lbl_comps):
            lbl.setText("—")
        for val, bar in self._metric_rows:
            val.setText("—")
            bar.set_values(0.0, None)
        self._lbl_vs.setText("")
        self._set_badge("–")
        self._btn_profile.setEnabled(False)
        self._btn_map.setEnabled(False)
        self._frequency_graph.set_frequencies([])
        self._response.set_stats(None)
