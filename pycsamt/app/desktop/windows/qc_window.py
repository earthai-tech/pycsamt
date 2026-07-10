# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
QCDashboardWindow — independent floating QC diagnostics window.

Left params panel
─────────────────
  Category     clickable list (Overview · Coverage · Noise/SNR ·
                               Skew/Dim · Static Shift · Distortion)
  Plot         ComboBox — changes with selected category
  Description  small info label
  [↻ Run]      [⬆ Export]

Right content
─────────────
  Single large MplCanvas + NavigationToolbar
  (one focused plot at a time — much cleaner than tabbed combos)
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QComboBox,
    QLabel,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.qc_controller import (
    ALL_GROUPS,
    GROUP_ICONS,
    QCController,
    describe_plot,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    _icon,
    icon_button,
    make_group,
)


class QCDashboardWindow(PanelWindow):
    """Floating QC Dashboard — category nav + plot selector left, canvas right."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(
            title="QC Dashboard",
            session_key="qc_dashboard",
            params_width=260,
            icon_name="qc",
            parent=parent,
        )
        self.resize(1200, 800)
        self._ctrl = QCController()
        self._auto_rendered = False
        # Populate after UI is built
        self._populate_category_combo()
        self._on_category_changed(0)

    # ── Params panel ──────────────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:
        # ── Category selector ─────────────────────────────────────────
        grp_cat, lay_cat = make_group("Category")
        self._combo_category = QComboBox()
        self._combo_category.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._combo_category.currentIndexChanged.connect(self._on_category_changed)
        lay_cat.addWidget(self._combo_category)
        layout.addWidget(grp_cat)

        # ── Plot selector ─────────────────────────────────────────────
        grp_plot, lay_plot = make_group("Plot")
        self._combo_plot = QComboBox()
        self._combo_plot.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._combo_plot.currentIndexChanged.connect(self._on_plot_changed)
        lay_plot.addWidget(self._combo_plot)
        layout.addWidget(grp_plot)

        # ── Info / description ────────────────────────────────────────
        self._desc_lbl = QLabel("")
        self._desc_lbl.setWordWrap(True)
        self._desc_lbl.setObjectName("InfoLabel")
        self._desc_lbl.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._desc_lbl)

        # ── Actions ───────────────────────────────────────────────────
        grp_act, lay_act = make_group("Actions")
        self._btn_run    = icon_button("↻  Run / Refresh", "qc",
                                       "Render the selected QC plot")
        self._btn_export = icon_button("⬆  Export…", "export",
                                       "Save figure to file")
        self._btn_run.clicked.connect(self._on_run)
        self._btn_export.clicked.connect(self._on_export)
        lay_act.addWidget(self._btn_run)
        lay_act.addWidget(self._btn_export)
        layout.addWidget(grp_act)

        # Status
        self._status_lbl = QLabel("")
        self._status_lbl.setObjectName("InfoLabel")
        layout.addWidget(self._status_lbl)

    # ── Content panel ─────────────────────────────────────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        self._canvas = MplCanvas(self, toolbar=True)
        layout.addWidget(self._canvas)

    # ── Public API ────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        super().set_sites(sites)
        self._ctrl.set_sites(sites)
        self._auto_rendered = False
        self._auto_render_if_ready()

    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._ctrl.dark = dark

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._auto_render_if_ready()

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_category_changed(self, row: int) -> None:
        if row < 0 or row >= len(ALL_GROUPS):
            return
        _label, plots = ALL_GROUPS[row]
        self._combo_plot.blockSignals(True)
        self._combo_plot.clear()
        for label, _fn, _has_ax in plots:
            self._combo_plot.addItem(label)
        self._combo_plot.blockSignals(False)
        self._combo_plot.setCurrentIndex(0)
        self._update_desc(row, 0)

    def _on_plot_changed(self, row: int) -> None:
        self._update_desc(self._combo_category.currentIndex(), row)

    def _on_run(self) -> None:
        cat_row  = self._combo_category.currentIndex()
        plot_row = self._combo_plot.currentIndex()
        if cat_row < 0 or plot_row < 0:
            return
        _label, plots = ALL_GROUPS[cat_row]
        if plot_row >= len(plots):
            return
        _plot_label, fn_name, has_ax = plots[plot_row]
        if self._ctrl._sites is None:
            self._status_lbl.setText("Load survey data first.")
            return
        self._status_lbl.setText(f"Running {fn_name}…")
        self._btn_run.setEnabled(False)
        try:
            new_fig = self._ctrl.draw(fn_name, has_ax, self._canvas.figure)
            if new_fig is not None:
                self._canvas.show_figure(new_fig)
            else:
                self._canvas.draw()
            self._status_lbl.setText("Done.")
        except Exception as exc:
            self._status_lbl.setText(f"Error: {exc}")
        finally:
            self._btn_run.setEnabled(True)

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )
        ExportDialog(figure=self._canvas.figure, parent=self).exec()

    # ── Helpers ───────────────────────────────────────────────────────

    def _auto_render_if_ready(self) -> None:
        if self._auto_rendered or self._ctrl._sites is None:
            return
        if not self.isVisible():
            return
        self._auto_rendered = True
        QTimer.singleShot(0, self._on_run)

    def _populate_category_combo(self) -> None:
        self._combo_category.blockSignals(True)
        for group_label, _plots in ALL_GROUPS:
            icon_name = GROUP_ICONS.get(group_label, "qc")
            icon = _icon(icon_name)
            if icon.isNull():
                self._combo_category.addItem(group_label)
            else:
                self._combo_category.addItem(icon, group_label)
        self._combo_category.blockSignals(False)

    def _update_desc(self, cat_row: int, plot_row: int) -> None:
        try:
            _label, plots = ALL_GROUPS[cat_row]
            plot_label, fn_name, _has_ax = plots[plot_row]
            desc = describe_plot(fn_name)
            self._desc_lbl.setText(
                f"<b>{plot_label}</b><br/>"
                f"<small style='color:#888'>{desc}</small>"
            )
        except (IndexError, Exception):
            self._desc_lbl.setText("")
