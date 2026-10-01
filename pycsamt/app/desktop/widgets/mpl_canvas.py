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

import warnings

import matplotlib as _mpl
from matplotlib.figure import Figure
from matplotlib.layout_engine import PlaceHolderLayoutEngine
from PySide6.QtCore import QSize, QTimer
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QSizePolicy,
    QStyle,
    QToolButton,
    QVBoxLayout,
    QWidget,
)


def _qt_mpl_classes():
    from matplotlib.backends.backend_qtagg import (
        FigureCanvasQTAgg,
        NavigationToolbar2QT,
    )

    return FigureCanvasQTAgg, NavigationToolbar2QT


class MplCanvas(QWidget):
    """
    A QWidget containing a matplotlib Figure (+ optional toolbar).

    Parameters
    ----------
    parent : QWidget, optional
    toolbar : bool
        Show compact navigation and export controls above the canvas.
    """

    def __init__(
        self,
        parent: QWidget | None = None,
        toolbar: bool = True,
    ) -> None:
        super().__init__(parent)
        FigureCanvasQTAgg, _ = _qt_mpl_classes()

        # Embedded plots manage their own margins. Automatic tight layout runs
        # at every Qt paint and is costly for colorbars and station labels.
        # Only one detached snapshot window is kept open per canvas at a
        # time; the toolbar button is disabled while it is open so clicking
        # it again re-raises that window instead of stacking up duplicates.
        self._detached_dialog = None
        self._detach_action = None
        self._refresh_btn = None
        self._refresh_connected = False
        self._layout_positions = None
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)

        self._canvas = FigureCanvasQTAgg(self.figure)
        self._canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._canvas, 1)

        self._toolbar = None
        self._tools = QWidget(self)
        tools_layout = QHBoxLayout(self._tools)
        tools_layout.setContentsMargins(4, 2, 4, 2)
        tools_layout.setSpacing(2)
        tools_layout.addStretch()
        layout.insertWidget(0, self._tools)
        self._has_toolbar = toolbar
        if toolbar:
            self._build_toolbar()
        self._tools.setVisible(toolbar)

        self._fit_timer = QTimer(self)
        self._fit_timer.setSingleShot(True)
        self._fit_timer.timeout.connect(self.fit_to_view)

    def _build_toolbar(self) -> None:
        """Recreate navigation bindings whenever the underlying canvas changes."""
        if self._toolbar is not None:
            self._tools.layout().removeWidget(self._toolbar)
            self._toolbar.deleteLater()
        _, NavigationToolbar2QT = _qt_mpl_classes()
        self._toolbar = NavigationToolbar2QT(self._canvas, self._tools)
        if hasattr(self._toolbar, "setIconSize"):
            self._toolbar.setIconSize(QSize(18, 18))
            self._toolbar.setMovable(False)
            # Keep navigation and axis editing, replace Matplotlib's save UI.
            for action in list(self._toolbar.actions()):
                if action.text() == "Save":
                    self._toolbar.removeAction(action)
            if hasattr(self._toolbar, "locLabel"):
                self._toolbar.locLabel.hide()
            export = self._toolbar.addAction(
                self.style().standardIcon(QStyle.StandardPixmap.SP_DialogSaveButton),
                "Export figure (PNG, SVG, PDF, EPS, TIFF)", self.export_figure,
            )
            export.setToolTip("Export publication figure…")
            detached = self._toolbar.addAction(
                self.style().standardIcon(QStyle.StandardPixmap.SP_TitleBarMaxButton),
                "Open in separate plot window", self.open_detached,
            )
            detached.setToolTip("Open an independent plot window")
            # _build_toolbar() reruns on every show_figure() (new toolbar
            # action objects each time), so re-apply the disabled state
            # here rather than only at the moment the dialog is opened.
            detached.setEnabled(self._detached_dialog is None)
            self._detach_action = detached
        self._toolbar.setStyleSheet("""
            QToolBar { background: transparent; border: none; spacing: 2px; }
            QToolButton { background: transparent; border: none;
                          border-radius: 5px; padding: 5px; color: #334155; }
            QToolButton:hover { background: #e2e8f0; }
            QToolButton:pressed, QToolButton:checked { background: #bfdbfe; }
            QToolButton:focus { border: 1px solid #3b82f6; }
        """)
        self._tools.layout().addWidget(self._toolbar)

    def _apply_panel_style(self) -> None:
        from matplotlib.text import Text
        from pycsamt.api.control import PYCSAMT_CONTROL

        settings = PYCSAMT_CONTROL.panel
        transparent = settings.background == "transparent"
        bg = "none" if transparent else "white"
        self.figure.set_facecolor(bg)
        self._canvas.setStyleSheet(
            "background: transparent;" if transparent else "background: white;"
        )
        self._tools.setStyleSheet("background: white;")
        self._tools.setVisible(self._has_toolbar and settings.toolbar)
        self._position_refresh_btn()
        for ax in self.figure.axes:
            ax.set_facecolor(bg)
            ax.tick_params(colors="#475569")
            for spine in ax.spines.values():
                spine.set_color("#94a3b8")
            # Theme-generated labels must remain legible on publication white.
            for text in ax.findobj(Text):
                text.set_color("#334155")
                # ... and so must their boxes: a dark-theme restyle had
                # painted them near-black, leaving the (now slate) text
                # unreadable -- the phase-tensor |β| legend, size reference
                # and PT-strip station label on dark.
                box = text.get_bbox_patch()
                if box is not None:
                    box.set_facecolor("white")
                    box.set_edgecolor("#cbd5e1")
                    box.set_alpha(0.92)
            legend = ax.get_legend()
            if legend is not None:
                legend.get_frame().set_facecolor(bg)
        for text in self.figure.texts:
            text.set_color("#334155")

    def export_figure(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog

        ExportDialog(figure=self.figure, parent=self).exec()

    def open_detached(self) -> None:
        """Open an independent snapshot with its own navigation and export.

        Only one detached window is kept per canvas: while one is open the
        toolbar button is disabled, so a second click can't happen and
        stack up duplicate windows — clicking is simply not possible until
        the existing one is closed.
        """
        if self._detached_dialog is not None:
            self._detached_dialog.raise_()
            self._detached_dialog.activateWindow()
            return

        import copy
        from PySide6.QtWidgets import QMessageBox

        try:
            figure = copy.deepcopy(self.figure)
        except Exception as exc:
            QMessageBox.warning(self, "Open plot", str(exc))
            return
        dialog = QDialog(self.window())
        dialog.setWindowTitle("Plot view")
        dialog.resize(1000, 720)
        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(0, 0, 0, 0)
        view = MplCanvas(dialog)
        layout.addWidget(view)
        view.show_figure(figure)
        self._detached_dialog = dialog
        if self._detach_action is not None:
            self._detach_action.setEnabled(False)

        def release_dialog():
            view._fit_timer.stop()
            self._detached_dialog = None
            if self._detach_action is not None:
                self._detach_action.setEnabled(True)
            dialog.deleteLater()

        dialog.finished.connect(release_dialog)
        dialog.show()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._fit_timer.start(0)

    # ── Convenience pass-throughs ─────────────────────────────────────

    def draw(self) -> None:
        self._layout_positions = None
        self._apply_panel_style()
        self._fit_timer.start(0)

    def fit_to_view(self) -> None:
        """Fit the embedded figure to the current viewport without resizing it.

        Matplotlib helpers often create figures at a fixed inch size.  An
        embedded canvas must instead follow the Qt viewport, otherwise titles
        and axis labels can sit outside the visible area after the figure is
        transferred or the panel is resized.
        """
        # Too small to hold a plot (hidden panel, collapsed splitter, first
        # paint): laying out would only collapse the axes; resizeEvent
        # fits it again once the canvas has a real size.
        if self._canvas.width() < 40 or self._canvas.height() < 40:
            return
        # Qt geometry is logical pixels; Matplotlib's DPI includes the screen
        # device ratio. Mixing the two shrinks plots on Windows at 150/200%.
        ratio = self._canvas.device_pixel_ratio
        width_px = max(self._canvas.width(), 1) * ratio
        height_px = max(self._canvas.height(), 1) * ratio
        dpi = self.figure.get_dpi() or 100
        self.figure.set_size_inches(
            width_px / dpi,
            height_px / dpi,
            forward=False,
        )

        # A single axes is safe to lay out responsively. Multi-axes figures
        # may contain manually positioned colorbars/insets and retain the
        # layout defined by their plotting function.
        visible_axes = [ax for ax in self.figure.axes if ax.get_visible()]
        layout_engine = self.figure.get_layout_engine()
        from pycsamt.api.control import PYCSAMT_CONTROL

        # Subplot-backed colorbars participate in tight_layout. Manually
        # positioned insets retain their geometry and scientific aspect ratio.
        compatible = all(
            getattr(ax, "get_subplotspec", lambda: None)() is not None
            for ax in visible_axes
        )
        if not compatible and self._layout_positions:
            for ax, position in self._layout_positions.items():
                if ax in visible_axes and ax.get_axes_locator() is None:
                    in_layout = ax.get_in_layout()
                    ax.set_position(position)
                    ax.set_in_layout(in_layout)
        if (PYCSAMT_CONTROL.panel.fit_layout and compatible
                and visible_axes
                and (layout_engine is None
                     or isinstance(layout_engine, PlaceHolderLayoutEngine))):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                try:
                    self.figure.tight_layout(pad=1.4)
                except Exception:
                    pass
        self._apply_panel_style()
        if PYCSAMT_CONTROL.panel.fit_layout:
            self._reserve_label_space(visible_axes)
        self._canvas.draw_idle()

    def _reserve_label_space(self, axes) -> None:
        """Measure actual decorations, including axes with manual colorbars.

        A fixed bottom percentage cannot accommodate multiline or rotated
        labels. Fit the axes group into the remaining space after measuring
        overflow in display pixels, preserving relative panel positions.
        """
        from matplotlib.transforms import Bbox
        import numpy as np

        if not axes:
            return
        if self._layout_positions is None:
            self._layout_positions = {
                ax: ax.get_position(original=True).frozen() for ax in axes
            }
        for _ in range(3):
            self._canvas.draw()
            renderer = self._canvas.get_renderer()
            bounds = [ax.get_tightbbox(renderer) for ax in axes]
            bounds = [b for b in bounds if b is not None and np.isfinite(b.extents).all()]
            if not bounds:
                return
            tight = Bbox.union(bounds)
            width, height = self.figure.bbox.size
            padding = 10 * self._canvas.device_pixel_ratio
            left = max(0, padding - tight.x0) / width
            bottom = max(0, padding - tight.y0) / height
            right = max(0, tight.x1 + padding - width) / width
            top = max(0, tight.y1 + padding - height) / height
            if max(left, bottom, right, top) < 0.001:
                return
            group = Bbox.union([ax.get_position(original=True) for ax in axes])
            sx = max(0.2, (group.width - left - right) / group.width)
            sy = max(0.2, (group.height - bottom - top) / group.height)
            for ax in axes:
                # Divider/inset axes follow their parent during the next draw.
                if ax.get_axes_locator() is not None:
                    continue
                pos = ax.get_position(original=True)
                in_layout = ax.get_in_layout()
                ax.set_position([
                    group.x0 + left + (pos.x0 - group.x0) * sx,
                    group.y0 + bottom + (pos.y0 - group.y0) * sy,
                    pos.width * sx, pos.height * sy,
                ])
                ax.set_in_layout(in_layout)

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
        FigureCanvasQTAgg, _NavigationToolbar2QT = _qt_mpl_classes()

        # ── Remove and discard the old canvas ─────────────────────────
        layout.removeWidget(self._canvas)
        self._canvas.setParent(None)
        self._canvas.deleteLater()

        # ── Resize figure to fill the available widget area ────────────
        # Use the widget's current pixel size; fall back to window size.
        avail_w = self.width() or 900
        avail_h = self.height() or 600
        toolbar_h = self._tools.sizeHint().height() if self._has_toolbar else 0
        canvas_h = max(avail_h - toolbar_h - 4, 100)
        dpi = fig.get_dpi() or 100
        fig.set_size_inches(avail_w / dpi, canvas_h / dpi)

        # ── Create a fresh canvas for this figure ──────────────────────
        self._canvas = FigureCanvasQTAgg(fig)
        self._canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        # The canvas stretches below the compact control strip.
        layout.addWidget(self._canvas, 1)

        # ── Update toolbar ─────────────────────────────────────────────
        if self._has_toolbar:
            self._build_toolbar()

        # ── Store references ───────────────────────────────────────────
        self.figure = fig
        self._layout_positions = None
        self.axes = fig.axes[0] if fig.axes else fig.add_subplot(111)

        # Fit and draw after the replacement canvas has joined the layout.
        self._fit_timer.start(0)
        self._position_refresh_btn()

    def closeEvent(self, event):  # noqa: N802
        self._fit_timer.stop()
        super().closeEvent(event)

    def resizeEvent(self, event) -> None:  # noqa: N802
        """Reflow labels after the canvas viewport changes size."""
        super().resizeEvent(event)
        if hasattr(self, "_fit_timer"):
            self._fit_timer.start(60)
        self._position_refresh_btn()

    # ── Floating hover-reveal refresh button ────────────────────────────

    def set_refresh_callback(
        self,
        callback,
        tooltip: str = "Hard refresh — recompute this plot immediately",
    ) -> None:
        """Attach (or reassign) a floating refresh button to this canvas.

        The button floats over the top-left corner of the *drawing area*
        — below the navigation toolbar strip when one is shown — so it
        never overlaps the toolbar's own icons, in particular "Open in
        separate plot window" which sits at the toolbar's right edge.
        It stays nearly transparent until hovered (see the
        ``PlotRefreshOverlay`` QSS rules), matching the existing pattern
        rather than a hard show/hide on canvas hover.
        """
        if self._refresh_btn is None:
            self._refresh_btn = QToolButton(self)
            self._refresh_btn.setObjectName("PlotRefreshOverlay")
            self._refresh_btn.setText("↻")
            self._refresh_btn.setFixedSize(32, 32)
        elif self._refresh_connected:
            self._refresh_btn.clicked.disconnect()
        self._refresh_btn.setToolTip(tooltip)
        self._refresh_btn.setAccessibleName(tooltip)
        self._refresh_btn.clicked.connect(callback)
        self._refresh_connected = True
        self._refresh_btn.show()
        self._refresh_btn.raise_()
        self._position_refresh_btn()

    @property
    def refresh_button(self) -> QToolButton | None:
        """The floating refresh button, once set_refresh_callback() has run."""
        return self._refresh_btn

    def _position_refresh_btn(self) -> None:
        if self._refresh_btn is None:
            return
        margin = 6
        top = margin
        if self._has_toolbar and self._tools.isVisible():
            top += self._tools.height()
        self._refresh_btn.move(margin, top)
        self._refresh_btn.raise_()

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
        self._apply_panel_style()
        self._canvas.draw_idle()


# ──────────────────────────────────────────────────────────────────────────────
# Global matplotlib theme helpers — call from MainWindow._apply_theme()
# ──────────────────────────────────────────────────────────────────────────────

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


# Global rcParams stay *publication* neutral whatever the app theme.
# They used to carry the theme colours (#e6e9ef / #1e1e2e), and because
# rcParams are process-wide every figure made in the desktop -- including
# the ones the library saves to disk (pipeline QC plots, agent reports,
# batch exports) and every ``savefig`` via ``savefig.facecolor`` -- came
# out grey (light theme) or navy (dark theme).  On-screen canvases style
# themselves (``MplCanvas._apply_panel_style``), so nothing needs them.
_PUBLICATION_RC = {
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "text.color": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.grid": False,
    "legend.facecolor": "inherit",
    "legend.edgecolor": "0.8",
    "legend.labelcolor": "None",
    "patch.edgecolor": "black",
}


def apply_publication_rc() -> None:
    """Reset the colour-related rcParams to publication (white) defaults."""
    _mpl.rcParams.update(_PUBLICATION_RC)


def apply_mpl_dark_theme() -> None:
    """Kept for compatibility: the app theme no longer leaks into
    matplotlib's global rcParams (see :data:`_PUBLICATION_RC`)."""
    apply_publication_rc()


def apply_mpl_light_theme() -> None:
    """Kept for compatibility: see :func:`apply_mpl_dark_theme`."""
    apply_publication_rc()
