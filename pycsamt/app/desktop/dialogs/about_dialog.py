# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AboutDialog — branded "About pycsamt v2" window.

Features
--------
* Hero banner painted with QPainter:
    - Dark-blue linear gradient background
    - pycsamt SVG logo centred and scaled
    - Version badge (amber)
    - Tagline text (steel blue)
    - Decorative dual EM-wave sinusoids at the bottom
    - Soft radial corner glow for depth
* Content panel:
    - One-line description + 6 feature bullets
    - Metadata table (version / Python / license / platform / author)
    - Clickable Documentation and GitHub links (open browser)
    - Built-with credits
* Clean Close button
"""

from __future__ import annotations

import math
import platform
import sys

from PySide6.QtCore import QRectF, Qt, QUrl
from PySide6.QtGui import (
    QColor,
    QDesktopServices,
    QFont,
    QLinearGradient,
    QPainter,
    QPainterPath,
    QPen,
    QRadialGradient,
)
from PySide6.QtWidgets import (
    QDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop import branding

# ── Paths ─────────────────────────────────────────────────────────────────────
_LOGO_SVG = branding.LOGO_SVG

# ── External links ────────────────────────────────────────────────────────────
_URL_DOCS = branding.URL_DOCS
_URL_GH = branding.URL_GITHUB

# ── Brand palette ─────────────────────────────────────────────────────────────
_C_BG_DARK = QColor("#0c1f4a")
_C_BG_MID = QColor("#1a3a7c")
_C_BG_EDGE = QColor("#0a1830")
_C_AMBER = QColor("#fbb040")
_C_STEEL = QColor("#8cb4de")
_C_WHITE = QColor("#ffffff")
_C_LINK = "#1a5cb8"


# ── Hero banner ───────────────────────────────────────────────────────────────


class _HeroBanner(QWidget):
    """
    Custom-painted hero widget:
      • Gradient background  (dark-navy → mid-blue → dark-navy)
      • Radial corner glow for depth
      • pycsamt SVG logo (centred, scaled)
      • Version badge (amber)
      • Tagline (steel-blue)
      • Two overlapping EM sinusoidal waves (decorative, bottom strip)
    """

    _LOGO_ASPECT = 400.29 / 66.91  # from SVG viewBox

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(175)
        self.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )

        self._renderer = None
        if _LOGO_SVG.exists():
            try:
                from PySide6.QtSvg import QSvgRenderer

                r = QSvgRenderer(str(_LOGO_SVG))
                if r.isValid():
                    self._renderer = r
            except Exception:
                pass

    def paintEvent(self, event) -> None:
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        rect = self.rect()

        # ── Background gradient ────────────────────────────────────────
        bg = QLinearGradient(0, 0, w, h)
        bg.setColorAt(0.00, _C_BG_EDGE)
        bg.setColorAt(0.40, _C_BG_MID)
        bg.setColorAt(1.00, _C_BG_DARK)
        p.fillRect(rect, bg)

        # ── Radial glow (top-left corner, brand blue-cyan) ─────────────
        glow = QRadialGradient(w * 0.12, h * 0.2, w * 0.45)
        glow.setColorAt(0.0, QColor(45, 120, 210, 55))
        glow.setColorAt(1.0, QColor(45, 120, 210, 0))
        p.fillRect(rect, glow)

        # ── EM-wave decoration (bottom strip) ─────────────────────────
        wave_y = h - 32
        amp1 = 9.0
        amp2 = 5.5
        freq = 5.5 * math.pi / max(w, 1)

        path1 = QPainterPath()
        path2 = QPainterPath()
        path1.moveTo(0, wave_y)
        path2.moveTo(0, wave_y + 4)
        for x in range(1, w + 1):
            path1.lineTo(x, wave_y - amp1 * math.sin(freq * x))
            path2.lineTo(
                x, wave_y + 4 - amp2 * math.sin(freq * x + math.pi * 0.65)
            )

        pen1 = QPen(QColor(255, 255, 255, 38), 1.6)
        pen2 = QPen(QColor(251, 176, 64, 48), 1.2)
        p.setPen(pen1)
        p.drawPath(path1)
        p.setPen(pen2)
        p.drawPath(path2)

        # ── Logo SVG ───────────────────────────────────────────────────
        if self._renderer:
            lh = 54
            lw = int(lh * self._LOGO_ASPECT)
            lx = (w - lw) // 2
            ly = 22
            self._renderer.render(p, QRectF(lx, ly, lw, lh))

        # ── Version badge ──────────────────────────────────────────────
        ver = branding.get_version()

        vfont = QFont()
        vfont.setFamily("Arial")
        vfont.setPixelSize(15)
        vfont.setBold(True)
        p.setFont(vfont)
        p.setPen(_C_AMBER)
        p.drawText(
            QRectF(0, 88, w, 26),
            Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter,
            f"v {ver}",
        )

        # ── Tagline ────────────────────────────────────────────────────
        tfont = QFont()
        tfont.setFamily("Arial")
        tfont.setPixelSize(10)
        p.setFont(tfont)
        p.setPen(_C_STEEL)
        p.drawText(
            QRectF(0, 114, w, 18),
            Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter,
            branding.TAGLINE,
        )

        p.end()


# ── Link button (styled QPushButton that opens a URL) ─────────────────────────


class _LinkButton(QPushButton):
    def __init__(
        self,
        icon_char: str,
        label: str,
        url: str,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(f"  {icon_char}  {label}", parent)
        self._url = url
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFlat(True)
        self.setStyleSheet(
            f"QPushButton {{"
            f"  color: {_C_LINK};"
            f"  text-align: left;"
            f"  border: 1px solid #b0c8e8;"
            f"  border-radius: 4px;"
            f"  padding: 5px 14px;"
            f"  font-size: 12px;"
            f"}}"
            f"QPushButton:hover {{"
            f"  background: #e8f0fb;"
            f"  border-color: {_C_LINK};"
            f"}}"
        )
        self.clicked.connect(self._open)

    def _open(self) -> None:
        QDesktopServices.openUrl(QUrl(self._url))


# ── About dialog ──────────────────────────────────────────────────────────────


class AboutDialog(QDialog):
    """
    Branded "About pycsamt v2" window.

    Parameters
    ----------
    parent : QWidget, optional
    license_manager : LicenseManager, optional
        Defaults to :class:`~pycsamt.app.desktop.licensing.null_manager
        .NullLicenseManager` (unlimited trial) when the caller doesn't
        pass one. ``main_window.py`` passes
        ``licensing.manager.get_default_manager()`` for the real,
        persisted license/trial state.
    """

    def __init__(
        self,
        parent: QWidget | None = None,
        license_manager=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("About pycsamt")
        self.setFixedWidth(530)
        self.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum
        )
        self._build_ui(license_manager)

    def _build_ui(self, license_manager=None) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(0)
        root.setContentsMargins(0, 0, 0, 0)

        # ── Hero ──────────────────────────────────────────────────────
        root.addWidget(_HeroBanner())

        # ── Content ───────────────────────────────────────────────────
        body = QWidget()
        lay = QVBoxLayout(body)
        lay.setSpacing(10)
        lay.setContentsMargins(24, 16, 24, 8)

        # Description
        desc = QLabel(
            "A comprehensive open-source Python library for processing, analysing, "
            "and inverting audio-frequency Magnetotelluric, CSAMT/CSEM, and "
            "TEM/TDEM electromagnetic data — with integrated AI/ML capabilities."
        )
        desc.setWordWrap(True)
        desc.setAlignment(Qt.AlignmentFlag.AlignJustify)
        lay.addWidget(desc)

        # Feature bullets
        feat = QLabel(
            "<ul style='margin:4px 0 4px 0; padding-left:20px; line-height:185%;'>"
            "<li>Phase-tensor, strike analysis and dimensionality classification</li>"
            "<li>Impedance tensor correction &amp; static-shift removal</li>"
            "<li>2-D / 3-D inversion pipeline (Occam-2D, GCN-3D inverter)</li>"
            "<li>AI / ML resistivity prediction and deep-learning inversion</li>"
            "<li>TDEM / TEM waveform transformation and sounding analysis</li>"
            "<li>Geological cross-section interpretation &amp; spatial interpolation</li>"
            "</ul>"
        )
        feat.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(feat)

        # Divider
        div = QFrame()
        div.setFrameShape(QFrame.Shape.HLine)
        div.setFrameShadow(QFrame.Shadow.Sunken)
        lay.addWidget(div)

        # Metadata
        ver = branding.get_version()

        py_ver = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
        plat = platform.system()

        meta = QLabel(
            "<table cellspacing='5' style='font-size:12px;'>"
            f"<tr>"
            f"  <td><b>Version</b></td><td>&nbsp;{ver}</td>"
            f"  <td>&nbsp;&nbsp;&nbsp;<b>Python</b></td><td>&nbsp;{py_ver}</td>"
            f"</tr>"
            f"<tr>"
            f"  <td><b>License</b></td><td>&nbsp;{branding.LICENSE_SPDX}</td>"
            f"  <td>&nbsp;&nbsp;&nbsp;<b>Platform</b></td><td>&nbsp;{plat}</td>"
            f"</tr>"
            f"<tr>"
            f"  <td><b>Author&nbsp;</b></td>"
            f"  <td colspan='3'>&nbsp;{branding.AUTHOR_NAME} &nbsp;"
            f"<a href='mailto:{branding.AUTHOR_EMAIL}' style='color:{_C_LINK};'>"
            f"{branding.AUTHOR_EMAIL}</a></td>"
            f"</tr>"
            f"</table>"
        )
        meta.setTextFormat(Qt.TextFormat.RichText)
        meta.setOpenExternalLinks(True)
        lay.addWidget(meta)

        # License / trial status — see widgets/license_page.py for the
        # full activation UI (Preferences ▸ License tab); this is a
        # compact read-only summary sharing the same status wording.
        from pycsamt.app.desktop.licensing.null_manager import (
            NullLicenseManager,
        )
        from pycsamt.app.desktop.widgets.license_page import (
            status_badge_html,
        )

        manager = license_manager or NullLicenseManager()
        status = manager.status()
        trial = manager.trial_state()
        trial_note = ""
        if trial is not None and not trial.is_expired:
            trial_note = f" &nbsp; {trial.days_remaining} day(s) remaining"
        license_row = QLabel(
            f"<b>License status</b>&nbsp; {status_badge_html(status)}{trial_note}"
        )
        license_row.setObjectName("LicenseStatusRow")
        license_row.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(license_row)

        # Link buttons row
        link_row = QHBoxLayout()
        link_row.setSpacing(10)
        link_row.addWidget(_LinkButton("📄", "Documentation", _URL_DOCS))
        link_row.addWidget(_LinkButton("⑂", "GitHub Repository", _URL_GH))
        link_row.addStretch()
        lay.addLayout(link_row)

        # Built-with credits
        built = QLabel(
            "<div style='color:#888; font-size:10px; text-align:center;'>"
            "Built with &nbsp; NumPy · SciPy · Matplotlib · PySide6 · "
            "scikit-learn · PyTorch · pyproj"
            "</div>"
        )
        built.setTextFormat(Qt.TextFormat.RichText)
        built.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay.addWidget(built)

        root.addWidget(body)

        # ── Close button ──────────────────────────────────────────────
        btn_area = QWidget()
        btn_lay = QHBoxLayout(btn_area)
        btn_lay.setContentsMargins(24, 4, 24, 14)
        btn_lay.addStretch()
        close_btn = QPushButton("Close")
        close_btn.setDefault(True)
        close_btn.setFixedWidth(90)
        close_btn.clicked.connect(self.accept)
        btn_lay.addWidget(close_btn)
        root.addWidget(btn_area)
