# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
MplCanvas — FigureCanvasQTAgg wrapper with an optional NavigationToolbar.

Phase 1/2: Drop-in widget that hosts a matplotlib Figure.  All pycsamt
plot functions accept a Figure or Axes; this widget exposes both.

Usage::

    canvas = MplCanvas(parent=self)
    canvas.axes.plot(x, y, "o")
    canvas.draw()
"""

from __future__ import annotations

import matplotlib

matplotlib.use("QtAgg")  # must be set before importing Figure / canvas

from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT,
)
from matplotlib.figure import Figure
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)


class MplCanvas(QWidget):
    """
    A QWidget containing a matplotlib Figure (+ optional toolbar).

    Parameters
    ----------
    parent : QWidget, optional
    toolbar : bool
        Show the matplotlib navigation toolbar below the canvas.
    """

    def __init__(
        self,
        parent: QWidget | None = None,
        toolbar: bool = True,
    ) -> None:
        super().__init__(parent)

        self.figure = Figure(tight_layout=True)
        self.axes = self.figure.add_subplot(111)

        self._canvas = FigureCanvasQTAgg(self.figure)
        self._canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._canvas)

        if toolbar:
            self._toolbar = NavigationToolbar2QT(self._canvas, self)
            layout.addWidget(self._toolbar)
        else:
            self._toolbar = None

    # ── Convenience pass-throughs ─────────────────────────────────────

    def draw(self) -> None:
        self._canvas.draw_idle()

    def clear_axes(self) -> None:
        self.axes.cla()
        self.draw()

    def show_figure(self, obj) -> None:
        """Display a Figure or Axes returned by an agent/emtools function.

        Replaces the FigureCanvasQTAgg entirely so there are no stale pixel
        buffers or canvas/figure geometry mismatches (which caused hatched
        artifacts and green-fill artefacts around axes).
        """
        # Resolve Axes → Figure
        if hasattr(obj, "get_figure"):
            fig = obj.get_figure()
        elif hasattr(obj, "savefig"):
            fig = obj
        else:
            return

        layout = self.layout()

        # ── Remove and discard the old canvas ─────────────────────────
        layout.removeWidget(self._canvas)
        self._canvas.setParent(None)
        self._canvas.deleteLater()

        # ── Resize figure to fill the available widget area ────────────
        # Use the widget's current pixel size; fall back to window size.
        avail_w = self.width() or 900
        avail_h = self.height() or 600
        toolbar_h = self._toolbar.sizeHint().height() if self._toolbar else 0
        canvas_h = max(avail_h - toolbar_h - 4, 100)
        dpi = fig.get_dpi() or 100
        fig.set_size_inches(avail_w / dpi, canvas_h / dpi)

        # ── Create a fresh canvas for this figure ──────────────────────
        self._canvas = FigureCanvasQTAgg(fig)
        self._canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        # Insert above the toolbar (index 0)
        layout.insertWidget(0, self._canvas)

        # ── Update toolbar ─────────────────────────────────────────────
        if self._toolbar is not None:
            self._toolbar.canvas = self._canvas
            self._toolbar.update()

        # ── Store references ───────────────────────────────────────────
        self.figure = fig
        self.axes = fig.axes[0] if fig.axes else fig.add_subplot(111)

        # Draw after the layout has settled
        QTimer.singleShot(0, self._canvas.draw_idle)

    def apply_theme(self, dark: bool) -> None:
        """Re-style an already-rendered figure when the app theme changes.

        `Figure()` bakes rcParams at construction time, so switching themes
        via rcParams alone leaves existing figures untouched.  This method
        explicitly updates every axes in the figure and triggers a redraw.
        """
        c = _DARK_STYLE if dark else _LIGHT_STYLE
        self.figure.patch.set_facecolor(c["fig_bg"])
        for ax in self.figure.get_axes():
            try:
                ax.set_facecolor(c["axes_bg"])
                ax.tick_params(colors=c["tick"])
                ax.xaxis.label.set_color(c["fg"])
                ax.yaxis.label.set_color(c["fg"])
                ax.title.set_color(c["fg"])
                for spine in ax.spines.values():
                    spine.set_edgecolor(c["spine"])
            except Exception:
                pass
        self._canvas.draw_idle()


# ──────────────────────────────────────────────────────────────────────────────
# Global matplotlib theme helpers — call from MainWindow._apply_theme()
# ──────────────────────────────────────────────────────────────────────────────

import matplotlib as _mpl

# Theme colour tables — kept in sync with the QSS palettes.
_DARK_STYLE = dict(
    fig_bg="#1e1e2e",
    axes_bg="#181825",
    fg="#cdd6f4",
    tick="#a6adc8",
    spine="#45475a",
)
_LIGHT_STYLE = dict(
    fig_bg="#e6e9ef",
    axes_bg="#eff1f5",
    fg="#4c4f69",
    tick="#6c6f85",
    spine="#bcc0cc",
)


def apply_mpl_dark_theme() -> None:
    """Set matplotlib rcParams to match the dark QSS palette."""
    _mpl.rcParams.update(
        {
            "axes.facecolor": "#181825",
            "figure.facecolor": "#1e1e2e",
            "savefig.facecolor": "#1e1e2e",
            "axes.edgecolor": "#45475a",
            "axes.labelcolor": "#cdd6f4",
            "text.color": "#cdd6f4",
            "xtick.color": "#a6adc8",
            "ytick.color": "#a6adc8",
            "grid.color": "#313244",
            "axes.grid": True,
            "grid.linestyle": "--",
            "grid.alpha": 0.35,
            "lines.color": "#89b4fa",
            "patch.edgecolor": "#45475a",
            "legend.facecolor": "#1e1e2e",
            "legend.edgecolor": "#45475a",
            "legend.labelcolor": "#cdd6f4",
        }
    )


def apply_mpl_light_theme() -> None:
    """Set matplotlib rcParams to match the light QSS palette."""
    _mpl.rcParams.update(
        {
            "axes.facecolor": "#eff1f5",
            "figure.facecolor": "#e6e9ef",
            "savefig.facecolor": "#e6e9ef",
            "axes.edgecolor": "#bcc0cc",
            "axes.labelcolor": "#4c4f69",
            "text.color": "#4c4f69",
            "xtick.color": "#6c6f85",
            "ytick.color": "#6c6f85",
            "grid.color": "#ccd0da",
            "axes.grid": True,
            "grid.linestyle": "--",
            "grid.alpha": 0.5,
            "lines.color": "#1e66f5",
            "patch.edgecolor": "#bcc0cc",
            "legend.facecolor": "#eff1f5",
            "legend.edgecolor": "#bcc0cc",
            "legend.labelcolor": "#4c4f69",
        }
    )
