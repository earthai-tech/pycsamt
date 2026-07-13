# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PanelWindow — base class for all independent pycsamt panel windows.

Layout (always):
    ┌──────────────┬──────────────────────────────────────┐
    │  PARAMETERS  │  VISUALISATION / RESULTS             │
    │  (fixed      │  (stretches to fill remaining space) │
    │   ~260 px)   │                                      │
    └──────────────┴──────────────────────────────────────┘

Subclasses override:
    _build_params(layout)  → populate the left QVBoxLayout
    _build_content(layout) → populate the right QVBoxLayout
    _apply_icon()          → optional window icon
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QByteArray, Qt, Signal
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import (
    QGroupBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
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


# ── Shared helpers ────────────────────────────────────────────────────────────


def make_group(title: str) -> tuple[QGroupBox, QVBoxLayout]:
    """Return (QGroupBox, inner_layout) ready to populate."""
    gb = QGroupBox(title)
    gb.setObjectName("ParamsGroup")
    lay = QVBoxLayout(gb)
    lay.setContentsMargins(6, 10, 6, 6)
    lay.setSpacing(5)
    return gb, lay


def icon_button(
    text: str, icon_name: str = "", tooltip: str = ""
) -> QPushButton:
    btn = QPushButton(text)
    if icon_name:
        ic = _icon(icon_name)
        if not ic.isNull():
            btn.setIcon(ic)
    if tooltip:
        btn.setToolTip(tooltip)
    btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    return btn


# ── Base window ───────────────────────────────────────────────────────────────


class PanelWindow(QWidget):
    """
    Base class for every independent pycsamt panel window.

    Parameters
    ----------
    title : str
        Window title-bar text.
    session_key : str
        Key used to persist geometry inside SessionState.window_geometries.
    params_width : int
        Initial fixed width of the left params panel (default 270).
    icon_name : str, optional
        Icon filename (without extension) from resources/icons/.
    parent : QWidget, optional
    """

    # Emitted when the user closes the window (window is hidden, not destroyed)
    panel_closed = Signal()

    def __init__(
        self,
        title: str,
        session_key: str,
        params_width: int = 270,
        icon_name: str = "",
        parent: QWidget | None = None,
    ) -> None:
        flags = (
            Qt.WindowType.Window
            | Qt.WindowType.WindowCloseButtonHint
            | Qt.WindowType.WindowMinimizeButtonHint
            | Qt.WindowType.WindowMaximizeButtonHint
        )
        super().__init__(parent, flags)
        self.setWindowTitle(f"pycsamt — {title}")
        if icon_name:
            ic = _icon(icon_name)
            if not ic.isNull():
                self.setWindowIcon(ic)

        self._session_key = session_key
        self._params_width = params_width
        self._sites = None
        self._dark: bool = True

        self._build_shell()

    # ── Shell layout ──────────────────────────────────────────────────────────

    def _build_shell(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self._splitter = QSplitter(Qt.Orientation.Horizontal, self)
        self._splitter.setHandleWidth(4)

        # ── Left: scrollable params panel ─────────────────────────────
        params_outer = QWidget()
        params_outer.setObjectName("ParamsPanel")
        params_outer.setMinimumWidth(200)
        params_outer.setMaximumWidth(360)
        params_outer.setFixedWidth(self._params_width)

        params_v = QVBoxLayout(params_outer)
        params_v.setContentsMargins(0, 0, 0, 0)
        params_v.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)

        inner = QWidget()
        inner.setObjectName("ParamsInner")
        self._params_layout = QVBoxLayout(inner)
        self._params_layout.setContentsMargins(8, 8, 8, 8)
        self._params_layout.setSpacing(6)

        self._build_params(self._params_layout)
        self._params_layout.addStretch(1)

        scroll.setWidget(inner)
        params_v.addWidget(scroll)

        # ── Right: content panel ───────────────────────────────────────
        self._content_widget = QWidget()
        self._content_widget.setObjectName("ContentPanel")
        self._content_layout = QVBoxLayout(self._content_widget)
        self._content_layout.setContentsMargins(0, 0, 0, 0)
        self._content_layout.setSpacing(0)
        self._build_content(self._content_layout)

        self._splitter.addWidget(params_outer)
        self._splitter.addWidget(self._content_widget)
        self._splitter.setStretchFactor(0, 0)
        self._splitter.setStretchFactor(1, 1)

        outer.addWidget(self._splitter)

    # ── Subclass interface ────────────────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:
        """Populate the left params panel.  Override in subclass."""

    def _build_content(self, layout: QVBoxLayout) -> None:
        """Populate the right content area.  Override in subclass."""

    # ── Public API ────────────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        """Called by MainWindow after EDI data is loaded."""
        self._sites = sites

    def set_dark_mode(self, dark: bool) -> None:
        self._dark = dark
        # Re-style any MplCanvas instances that already hold a rendered figure
        # so the theme applies immediately without requiring a manual redraw.
        from pycsamt.app.desktop.widgets.mpl_canvas import (
            MplCanvas,
        )

        for canvas in self.findChildren(MplCanvas):
            canvas.apply_theme(dark)

    # ── Session geometry ──────────────────────────────────────────────────────

    def save_geometry_to(self, store: dict) -> None:
        store[self._session_key] = {
            "geometry": self.saveGeometry().toBase64().data().decode(),
            "visible": self.isVisible(),
        }

    def restore_geometry_from(self, store: dict) -> None:
        entry = store.get(self._session_key)
        if not entry:
            return
        geo = entry.get("geometry")
        if geo:
            try:
                self.restoreGeometry(QByteArray.fromBase64(geo.encode()))
            except Exception:
                pass

    # ── Qt overrides ──────────────────────────────────────────────────────────

    def closeEvent(self, event) -> None:
        """Hide instead of destroying so state is preserved."""
        self.hide()
        event.ignore()
        self.panel_closed.emit()
