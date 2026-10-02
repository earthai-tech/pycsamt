# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AboutDialog — the pyCSAMT Help & About panel (Help, F1).

Layout
------
* **Hero banner** (QPainter): navy gradient, logo, version badge, tagline,
  EM-wave decoration.
* Tabs:

  - **Overview** -- what pyCSAMT does and what is new in 2.6;
  - **Resources** -- documentation, desktop guide, tutorials, release
    notes, GitHub, issue tracker, and how to cite;
  - **System** -- version, Python, Qt, platform, license status and a
    *Copy system info* button for bug reports;
  - **Author** -- Laurent Kouadio, website and e-mail.

* Footer: copyright, license, *Close*.
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
    QGuiApplication,
    QLinearGradient,
    QPainter,
    QPainterPath,
    QPen,
    QRadialGradient,
)
from PySide6.QtWidgets import (
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QTabWidget,
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

_WIDTH = 600

#: what 2.6 brought to the desktop application (Overview tab)
_NEW_IN_26 = (
    ("Inversion Studio", "Occam-2D, ModEM 2-D/3-D and MARE2DEM runs with "
                         "live progress and a Solver Builder."),
    ("Interpretation Studio", "Lithology, boreholes, hydrogeology and "
                              "time-lapse views from one place."),
    ("QC Studio", "Per-line diagnostics and a pass / warn / fail station "
                  "summary."),
    ("Airborne & TDEM studios", "ZTEM, AFMAG, MobileMT; TEM soundings "
                                "converted to EDI."),
    ("PCSF 3-D viewer", "Fence, block and depth-slice scenes, handed to "
                        "Map View with their settings."),
    ("Edit menu", "Undo history, station and frequency editors, survey "
                  "operations."),
)


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
        *,
        detail: str = "",
        icon_name: str = "",
        color: str = _C_LINK,
        hover: str = "#e8f0fb",
        dark: bool = False,
    ) -> None:
        # an app SVG icon when given (emoji glyphs are missing from some
        # system fonts and drew as boxes)
        svg = branding.ICONS_DIR / f"{icon_name}.svg" if icon_name else None
        text = (f"  {label}" if svg is not None and svg.exists()
                else f"  {icon_char}  {label}")
        super().__init__(text, parent)
        if svg is not None and svg.exists():
            from PySide6.QtCore import QSize
            from PySide6.QtGui import QIcon

            icon = QIcon(str(svg))
            if dark:  # black artwork: recolour it like the menu icons
                try:
                    from PySide6.QtCore import QByteArray
                    from PySide6.QtGui import QPixmap
                    from PySide6.QtSvg import QSvgRenderer

                    from pycsamt.app.desktop.main_window import _recolor_svg

                    r = QSvgRenderer(QByteArray(_recolor_svg(
                        svg.read_text("utf-8"), target=color)))
                    pm = QPixmap(36, 36)
                    pm.fill(Qt.GlobalColor.transparent)
                    painter = QPainter(pm)
                    r.render(painter)
                    painter.end()
                    icon = QIcon(pm)
                except Exception:
                    pass
            self.setIcon(icon)
            self.setIconSize(QSize(18, 18))
        self._url = url
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFlat(True)
        self.setToolTip(detail or url)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Fixed)
        self.setStyleSheet(
            f"QPushButton {{"
            f"  color: {color};"
            f"  text-align: left;"
            f"  border: 1px solid #b0c8e8;"
            f"  border-radius: 6px;"
            f"  padding: 7px 14px;"
            f"  font-size: 12px;"
            f"}}"
            f"QPushButton:hover {{"
            f"  background: {hover};"
            f"  border-color: {color};"
            f"}}"
        )
        self.clicked.connect(self._open)

    def _open(self) -> None:
        QDesktopServices.openUrl(QUrl(self._url))


def _rich(html: str, *, links: bool = True) -> QLabel:
    lab = QLabel(html)
    lab.setTextFormat(Qt.TextFormat.RichText)
    lab.setWordWrap(True)
    lab.setOpenExternalLinks(links)
    lab.setTextInteractionFlags(
        Qt.TextInteractionFlag.TextBrowserInteraction)
    return lab


def _page() -> tuple[QWidget, QVBoxLayout]:
    w = QWidget()
    lay = QVBoxLayout(w)
    lay.setContentsMargins(22, 14, 22, 10)
    lay.setSpacing(8)
    return w, lay


def system_info(license_manager=None) -> dict[str, str]:
    """Version / environment facts shown on the System tab (and copied
    for bug reports)."""
    try:
        import PySide6
        from PySide6.QtCore import qVersion

        qt = f"{qVersion()} (PySide6 {PySide6.__version__})"
    except Exception:
        qt = "unknown"
    try:
        import pycsamt

        where = str(__import__("pathlib").Path(pycsamt.__file__).parent)
    except Exception:
        where = "unknown"
    status = ""
    if license_manager is not None:
        try:
            status = str(license_manager.status().name).replace("_", " ")
            status = status.capitalize()
        except Exception:
            status = ""
    info = {
        "pyCSAMT": branding.get_version(),
        "Python": platform.python_version(),
        "Qt": qt,
        "Platform": f"{platform.system()} {platform.release()} "
                    f"({platform.machine()})",
        "License": branding.LICENSE_SPDX,
        "Installed at": where,
    }
    if status:
        info["License status"] = status
    return info


class _Monogram(QWidget):
    """Initials in a brand-gradient disc (Author tab)."""

    def __init__(self, name: str, size: int = 58,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._text = "".join(part[0] for part in name.split()[:2]).upper()
        self.setFixedSize(size, size)

    def paintEvent(self, event) -> None:
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        r = QRectF(self.rect()).adjusted(1, 1, -1, -1)
        g = QLinearGradient(0, 0, r.width(), r.height())
        g.setColorAt(0.0, _C_BG_MID)
        g.setColorAt(1.0, _C_BG_DARK)
        p.setBrush(g)
        p.setPen(QPen(_C_AMBER, 2))
        p.drawEllipse(r)
        f = QFont("Arial")
        f.setPixelSize(int(r.height() * 0.38))
        f.setBold(True)
        p.setFont(f)
        p.setPen(_C_WHITE)
        p.drawText(r, Qt.AlignmentFlag.AlignCenter, self._text)
        p.end()


# ── About dialog ──────────────────────────────────────────────────────────────


class AboutDialog(QDialog):
    """
    The pyCSAMT Help & About panel.

    Parameters
    ----------
    parent : QWidget, optional
    license_manager : LicenseManager, optional
        Defaults to :class:`~pycsamt.app.desktop.licensing.null_manager
        .NullLicenseManager` (unlimited trial) when the caller doesn't
        pass one. ``main_window.py`` passes
        ``licensing.manager.get_default_manager()`` for the real,
        persisted license/trial state.
    tab : str, optional
        Tab to open on: ``"overview"`` (default), ``"resources"``,
        ``"system"`` or ``"author"``.
    dark : bool, optional
        Dark theme colours (detected from the application stylesheet when
        omitted).
    """

    TABS = ("overview", "resources", "system", "author")

    def __init__(
        self,
        parent: QWidget | None = None,
        license_manager=None,
        tab: str = "overview",
        dark: bool | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("About pycsamt")
        self.setFixedWidth(_WIDTH)
        self.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum
        )
        from pycsamt.app.desktop.licensing.null_manager import (
            NullLicenseManager,
        )

        self._manager = license_manager or NullLicenseManager()
        if dark is None:  # the app theme is a stylesheet, not a palette
            from PySide6.QtWidgets import QApplication

            qss = QApplication.instance().styleSheet() if                 QApplication.instance() else ""
            dark = "#1e1e2e" in qss or                 self.palette().window().color().lightness() < 128
        self._dark = bool(dark)
        self._c_head = "#8cb4de" if dark else "#1a3a7c"
        self._c_muted = "#bac2de" if dark else "#555"
        self._c_faint = "#9399b2" if dark else "#888"
        self._c_link = "#89b4fa" if dark else _C_LINK
        self._build_ui()
        if tab in self.TABS:
            self.tabs.setCurrentIndex(self.TABS.index(tab))

    # ── layout ────────────────────────────────────────────────────────
    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(0)
        root.setContentsMargins(0, 0, 0, 0)
        root.addWidget(_HeroBanner())

        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.addTab(self._overview(), "Overview")
        self.tabs.addTab(self._resources(), "Resources")
        self.tabs.addTab(self._system(), "System")
        self.tabs.addTab(self._author(), "Author")
        root.addWidget(self.tabs)

        foot = QWidget()
        fl = QHBoxLayout(foot)
        fl.setContentsMargins(22, 6, 22, 14)
        year = __import__("datetime").date.today().year
        copy = QLabel(f"<span style='color:{self._c_faint}; font-size:10px;'>© "
                      f"{year} {branding.COPYRIGHT_HOLDER} · "
                      f"{branding.LICENSE_SPDX} · open source</span>")
        copy.setTextFormat(Qt.TextFormat.RichText)
        fl.addWidget(copy, 1)
        close_btn = QPushButton("Close")
        close_btn.setDefault(True)
        close_btn.setFixedWidth(90)
        close_btn.clicked.connect(self.accept)
        fl.addWidget(close_btn)
        root.addWidget(foot)
        close_btn.setFocus()  # not a selectable label (drew highlighted)

    def _overview(self) -> QWidget:
        w, lay = _page()
        lay.addWidget(_rich(
            "An open-source Python suite for processing, analysing, "
            "modelling and interpreting magnetotelluric, AMT, CSAMT/CSEM, "
            "airborne and TEM electromagnetic data &mdash; from raw EDI to "
            "inversion and geological interpretation, with integrated "
            "AI/ML tools."))
        rows = "".join(
            f"<tr><td style='padding:3px 10px 3px 0; white-space:nowrap;'>"
            f"<b>{name}</b></td><td style='padding:3px 0; color:{self._c_muted};'>"
            f"{text}</td></tr>" for name, text in _NEW_IN_26)
        lay.addWidget(_rich(
            f"<div style='margin-top:4px;'><b style='color:{self._c_head};'>"
            "New in 2.6</b></div>"
            f"<table cellspacing='0' style='font-size:12px;'>{rows}</table>"))
        lay.addWidget(_rich(
            f"<span style='color:{self._c_faint}; font-size:11px;'>Press "
            "<b>F1</b> anywhere to reopen this panel · full notes under "
            "<i>Resources ▸ Release notes</i>.</span>"))
        lay.addStretch(1)
        return w

    def _resources(self) -> QWidget:
        w, lay = _page()
        grid = QGridLayout()
        grid.setSpacing(8)
        base = _URL_DOCS.rstrip("/")
        links = (
            ("docs", "Documentation", _URL_DOCS, "Guides, theory and API"),
            ("overview", "Desktop guide",
             f"{base}/applications/desktop/index.html",
             "Every window of this application"),
            ("summary", "Tutorials", f"{base}/tutorials/index.html",
             "Worked, end-to-end examples"),
            ("log", "Release notes", f"{base}/release_notes/index.html",
             "What changed in each version"),
            ("github", "GitHub repository", _URL_GH, "Source code"),
            ("chat", "Report an issue", branding.URL_ISSUES,
             "Bugs and feature requests"),
        )
        for i, (icon, label, url, detail) in enumerate(links):
            grid.addWidget(_LinkButton(
                "•", label, url, detail=detail, icon_name=icon,
                color=self._c_link, dark=self._dark,
                hover="#313244" if self._dark else "#e8f0fb"),
                i // 2, i % 2)
        lay.addLayout(grid)
        lay.addWidget(_rich(
            f"<div style='margin-top:8px;'><b style='color:{self._c_head};'>"
            "Cite pyCSAMT</b></div><div style='font-size:11px;'>"
            f"{branding.CITATION} &nbsp;<a href='{branding.URL_CITATION}' "
            f"style='color:{self._c_link};'>doi:{branding.CITATION_DOI}</a>"
            "</div>"))
        lay.addStretch(1)
        return w

    def _system(self) -> QWidget:
        w, lay = _page()
        info = system_info()
        rows = "".join(
            f"<tr><td style='padding:2px 14px 2px 0;'><b>{k}</b></td>"
            f"<td>{v}</td></tr>" for k, v in info.items())
        lay.addWidget(_rich(f"<table style='font-size:12px;'>{rows}</table>",
                            links=False))

        # License / trial status -- the full activation UI is
        # Preferences ▸ License (widgets/license_page.py); this is a
        # compact read-only summary sharing the same wording.
        from pycsamt.app.desktop.widgets.license_page import (
            status_badge_html,
        )

        status = self._manager.status()
        trial = self._manager.trial_state()
        trial_note = ""
        if trial is not None and not trial.is_expired:
            trial_note = f" &nbsp; {trial.days_remaining} day(s) remaining"
        license_row = QLabel(
            f"<b>License status</b>&nbsp; {status_badge_html(status)}"
            f"{trial_note}")
        license_row.setObjectName("LicenseStatusRow")
        license_row.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(license_row)

        row = QHBoxLayout()
        self.btn_copy = QPushButton("Copy system info")
        self.btn_copy.setToolTip("Copy these details to paste into a bug "
                                 "report")
        self.btn_copy.clicked.connect(self.copy_system_info)
        row.addWidget(self.btn_copy)
        self._copied = QLabel("")
        self._copied.setObjectName("InfoLabel")
        row.addWidget(self._copied, 1)
        lay.addLayout(row)
        lay.addWidget(_rich(
            f"<span style='color:{self._c_faint}; font-size:10px;'>Built with NumPy · "
            "SciPy · pandas · Matplotlib · Plotly · PySide6 · scikit-learn "
            "· PyTorch · pyproj</span>", links=False))
        lay.addStretch(1)
        return w

    def copy_system_info(self) -> str:
        info = system_info(self._manager)
        text = "\n".join(f"{k}: {v}" for k, v in info.items())
        QGuiApplication.clipboard().setText(text)
        self._copied.setText("Copied to the clipboard.")
        return text

    def _author(self) -> QWidget:
        w, lay = _page()
        site = branding.AUTHOR_WEBSITE
        head = QHBoxLayout()
        head.setSpacing(14)
        head.addWidget(_Monogram(branding.AUTHOR_NAME))
        head.addWidget(_rich(
            f"<div style='font-size:16px; color:{self._c_head};'><b>"
            f"{branding.AUTHOR_NAME}</b></div>"
            f"<div style='color:{self._c_muted};'>{branding.AUTHOR_ROLE}"
            "</div>"), 1)
        lay.addLayout(head)
        rows = (
            ("Website", f"<a href='{site}' style='color:{self._c_link};'>"
                        f"{site.split('//')[-1]}</a>"),
            ("E-mail", f"<a href='mailto:{branding.AUTHOR_EMAIL}' "
                       f"style='color:{self._c_link};'>{branding.AUTHOR_EMAIL}</a>"),
            ("GitHub", f"<a href='{_URL_GH}' style='color:{self._c_link};'>"
                       f"{_URL_GH.split('github.com/')[-1]}</a>"),
        )
        lay.addWidget(_rich(
            "<table style='font-size:12px; margin-top:6px;'>" + "".join(
                f"<tr><td style='padding:3px 16px 3px 0;'><b>{k}</b></td>"
                f"<td>{v}</td></tr>" for k, v in rows) + "</table>"))
        lay.addWidget(_rich(
            f"<div style='color:{self._c_faint}; font-size:11px; margin-top:6px;'>"
            "Questions, collaborations and feedback are welcome — by e-mail "
            "or on the GitHub issue tracker.</div>"))
        lay.addStretch(1)
        return w


__all__ = ["AboutDialog", "system_info"]
