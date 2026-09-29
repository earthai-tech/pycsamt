# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ColorbarWidget — matplotlib colorbar embedded in a narrow QWidget.

Phase 3: Wraps a slim matplotlib Figure that renders only a colourbar.
Call ``update_colorbar(cmap, vmin, vmax, label)`` to refresh without
recreating the widget.
"""

from __future__ import annotations

import matplotlib.cm as mcm
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
from PySide6.QtWidgets import (
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from pycsamt.compat.matplotlib import get_cmap


def _qt_canvas_class():
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

    return FigureCanvasQTAgg


class ColorbarWidget(QWidget):
    """
    A thin QWidget that displays a matplotlib colorbar.

    Parameters
    ----------
    orientation : 'vertical' | 'horizontal'
    parent : QWidget, optional
    """

    _THEME = {
        True: dict(fig_bg="#1e1e2e", tick="#a6adc8", label="#cdd6f4"),
        False: dict(fig_bg="#e6e9ef", tick="#6c6f85", label="#4c4f69"),
    }

    def __init__(
        self,
        orientation: str = "vertical",
        parent: QWidget | None = None,
        dark: bool = True,
    ) -> None:
        super().__init__(parent)
        self._orientation = orientation
        self._dark = dark
        self._cb = None
        # Cached params so set_dark_mode() can redraw without the caller
        # needing to re-supply the last cmap/vmin/vmax/label.
        self._last_cmap: str | mcolors.Colormap = "plasma"
        self._last_vmin = 0.0
        self._last_vmax = 1.0
        self._last_label = ""
        self._build_ui()

    # ── Construction ──────────────────────────────────────────────────

    def _build_ui(self) -> None:
        FigureCanvasQTAgg = _qt_canvas_class()

        if self._orientation == "vertical":
            figsize = (0.7, 3.5)
        else:
            figsize = (4.0, 0.5)

        self._fig = Figure(
            figsize=figsize, facecolor=self._THEME[self._dark]["fig_bg"]
        )
        self._ax = self._fig.add_axes([0.1, 0.05, 0.4, 0.9])
        self._canvas = FigureCanvasQTAgg(self._fig)
        self._canvas.setSizePolicy(
            QSizePolicy.Policy.Preferred,
            QSizePolicy.Policy.Expanding,
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._canvas)

        # Draw a default colorbar
        self.update_colorbar("plasma", 0, 1, "")

    # ── Public API ─────────────────────────────────────────────────────

    def set_dark_mode(self, dark: bool) -> None:
        """Re-theme the colorbar figure to match the app's dark/light mode.

        Previously this widget's figure background was hardcoded to the
        dark palette regardless of the active theme, so in light mode it
        stood out as a mismatched black box next to an otherwise
        light-themed panel.
        """
        if dark == self._dark:
            return
        self._dark = dark
        self.update_colorbar(
            self._last_cmap, self._last_vmin, self._last_vmax,
            self._last_label,
        )

    def clear(self) -> None:
        """Blank the colorbar (no gradient/ticks/label) for an empty state.

        Used instead of leaving a stale gradient visible when the parent
        panel has no data loaded — a colorbar with no data behind it is
        just as unprofessional as an empty axes with default 0..1 ticks.
        """
        self._fig.clear()
        self._fig.set_facecolor(self._THEME[self._dark]["fig_bg"])
        self._cb = None
        self._canvas.draw_idle()

    def update_colorbar(
        self,
        cmap: str | mcolors.Colormap = "plasma",
        vmin: float = 0.0,
        vmax: float = 1.0,
        label: str = "",
    ) -> None:
        """Redraw the colorbar with new range and label."""
        self._last_cmap, self._last_vmin = cmap, vmin
        self._last_vmax, self._last_label = vmax, label

        # Recreate axes each call to avoid matplotlib figure-ownership issues
        # when a previous colorbar is removed and the cax becomes orphaned.
        self._fig.clear()
        self._fig.set_facecolor(self._THEME[self._dark]["fig_bg"])
        if self._orientation == "vertical":
            self._ax = self._fig.add_axes([0.15, 0.05, 0.35, 0.90])
        else:
            self._ax = self._fig.add_axes([0.05, 0.20, 0.90, 0.35])

        if isinstance(cmap, str):
            cmap = get_cmap(cmap)
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        sm = mcm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        theme = self._THEME[self._dark]
        self._cb = self._fig.colorbar(
            sm,
            cax=self._ax,
            orientation=self._orientation,
        )
        self._cb.ax.tick_params(colors=theme["tick"], labelsize=8)
        if label:
            self._cb.set_label(label, color=theme["label"], fontsize=9)

        self._canvas.draw_idle()
