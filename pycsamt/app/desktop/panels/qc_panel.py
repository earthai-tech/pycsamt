# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
QCPanel — quality-control diagnostic panel.

Six QC categories, each with multiple selectable plots:

  Overview       — QC quicklook, coverage heatmap, confidence
  Coverage       — Frequency presence, polar coverage
  Noise / SNR    — SNR histogram, harmonic waterfall, off-diagonal psections
  Skew / Dim     — Skew traffic, dimensionality psection, phase-tensor skew
  Static Shift   — SS QC psection / profile / station curves, SS radar
  Distortion     — Near-surface, overprint, field zones, strike

Interaction model (mirrors correction window):
    [Category ▾] [Plot ▾] [↻ Refresh] [⬆ Export]
    ─────────────────────────────────────────────
    [                   MplCanvas                ]

Changing category populates the plot combo; changing the plot or pressing
↻ re-renders the canvas.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.qc_controller import (
    ALL_GROUPS,
    QCController,
)
from pycsamt.app.desktop.widgets.category_bar import (
    CategoryComboBar,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas


class QCPanel(QWidget):
    """
    QC diagnostic panel — category + plot dropdown, single canvas.

    Call :meth:`set_sites` after loading data to enable all plots.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._ctrl = QCController()
        self._build_ui()

    # ── Construction ──────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── Category + plot selector bar ──────────────────────────────
        self._bar = CategoryComboBar(show_item=True, parent=self)
        cats = [group_label for group_label, _ in ALL_GROUPS]
        self._bar.set_categories(cats)
        self._bar.category_changed.connect(self._on_category_changed)
        self._bar.item_changed.connect(self._on_item_changed)
        self._bar.refresh_clicked.connect(self.refresh)
        self._bar.export_clicked.connect(self._on_export)
        root.addWidget(self._bar)

        # ── Single canvas ─────────────────────────────────────────────
        self._canvas = MplCanvas(self, toolbar=False)
        root.addWidget(self._canvas)

        # Populate the item combo for the first category without triggering
        # a double-render (category_changed fires before items are set).
        self._populate_items(0)
        self._draw_empty()

    # ── Public API ────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        """Load a Sites collection and refresh the current plot."""
        self._ctrl.set_sites(sites)
        self.refresh()

    def set_dark_mode(self, dark: bool) -> None:
        self._ctrl.dark = dark
        self.refresh()

    def refresh(self) -> None:
        """Re-render the currently selected plot."""
        cat_idx  = self._bar.current_category()
        item_idx = self._bar.current_item()
        if cat_idx < 0 or cat_idx >= len(ALL_GROUPS):
            return
        _, plot_list = ALL_GROUPS[cat_idx]
        if item_idx < 0 or item_idx >= len(plot_list):
            return
        _label, fn_name, has_ax = plot_list[item_idx]
        new_fig = self._ctrl.draw(fn_name, has_ax, self._canvas.figure)
        if new_fig is not None:
            self._canvas.show_figure(new_fig)
        else:
            self._canvas.draw()

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_category_changed(self, idx: int) -> None:
        self._populate_items(idx)
        # Render first item of the new category
        self._on_item_changed(0)

    def _on_item_changed(self, _idx: int) -> None:
        self.refresh()

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )
        ExportDialog(figure=self._canvas.figure, parent=self).exec()

    # ── Helpers ───────────────────────────────────────────────────────

    def _populate_items(self, cat_idx: int) -> None:
        if cat_idx < 0 or cat_idx >= len(ALL_GROUPS):
            return
        _, plot_list = ALL_GROUPS[cat_idx]
        self._bar.set_items([label for label, _fn, _ax in plot_list])

    def _draw_empty(self) -> None:
        from pycsamt.app.desktop.controllers.plot_controller import (
            _annotate_empty,
            style_axes,
        )
        fig = self._canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        _annotate_empty(ax, "Load survey data to enable QC plots")
        style_axes(ax, self._ctrl.dark)
        self._canvas.draw()
