# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ColorbarWidget — matplotlib colorbar embedded in a narrow QWidget.

Phase 3: Wraps a slim matplotlib Figure that renders only a colourbar.
Call ``update_colorbar(cmap, vmin, vmax, label)`` to refresh without
recreating the widget.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("QtAgg")

import matplotlib.cm as mcm
import matplotlib.colors as mcolors
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg,
)
from matplotlib.figure import Figure
from PySide6.QtWidgets import (
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)


class ColorbarWidget(QWidget):
    """
    A thin QWidget that displays a matplotlib colorbar.

    Parameters
    ----------
    orientation : 'vertical' | 'horizontal'
    parent : QWidget, optional
    """

    def __init__(
        self,
        orientation: str = "vertical",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._orientation = orientation
        self._cb = None
        self._build_ui()

    # ── Construction ──────────────────────────────────────────────────

    def _build_ui(self) -> None:
        if self._orientation == "vertical":
            figsize = (0.7, 3.5)
        else:
            figsize = (4.0, 0.5)

        self._fig = Figure(figsize=figsize, facecolor="#1e1e2e")
        self._ax  = self._fig.add_axes([0.1, 0.05, 0.4, 0.9])
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

    def update_colorbar(
        self,
        cmap: str | mcolors.Colormap = "plasma",
        vmin: float = 0.0,
        vmax: float = 1.0,
        label: str = "",
    ) -> None:
        """Redraw the colorbar with new range and label."""
        # Recreate axes each call to avoid matplotlib figure-ownership issues
        # when a previous colorbar is removed and the cax becomes orphaned.
        self._fig.clear()
        if self._orientation == "vertical":
            self._ax = self._fig.add_axes([0.15, 0.05, 0.35, 0.90])
        else:
            self._ax = self._fig.add_axes([0.05, 0.20, 0.90, 0.35])

        if isinstance(cmap, str):
            cmap = mcm.get_cmap(cmap)
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        sm   = mcm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        self._cb = self._fig.colorbar(
            sm,
            cax=self._ax,
            orientation=self._orientation,
        )
        self._cb.ax.tick_params(colors="#a6adc8", labelsize=8)
        if label:
            self._cb.set_label(label, color="#cdd6f4", fontsize=9)

        self._canvas.draw_idle()
