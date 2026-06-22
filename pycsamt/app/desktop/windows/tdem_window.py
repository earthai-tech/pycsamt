# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TDEMWindow — independent floating window for Time Domain EM analysis.

Left params panel
─────────────────
  Load data    Folder browse button → reads all .AVG/.Z/.LOG files
  Data info    Shows loaded survey summary (n stations, n files, etc.)
  Category     Clickable list: Decay/Rho · Section · Map & Overview · Dashboard
  Plot         ComboBox — changes with category
  [↻ Run]      [⬆ Export]

Right content
─────────────
  Single large MplCanvas + NavigationToolbar.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QLabel,
    QProgressBar,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.tdem_controller import (
    TDEMController,
    TDEM_GROUPS,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base import PanelWindow, make_group, icon_button


class TDEMWindow(PanelWindow):
    """
    Floating TDEM Analysis window — data loader + plot controls left,
    large MplCanvas right.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(
            title="TDEM Analysis",
            session_key="tdem_window",
            params_width=265,
            icon_name="tdem",
            parent=parent,
        )
        self.resize(1240, 820)
        self._ctrl = TDEMController()
        self._populate_category_combo()
        self._on_category_changed(0)

    # ── Params panel (left) ───────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:
        # ── Data loading ──────────────────────────────────────────────
        grp_load, lay_load = make_group("Data")
        self._btn_browse = icon_button("📂  Browse folder…", "open",
                                       "Select folder with .AVG / .Z / .LOG files")
        self._btn_browse.clicked.connect(self._on_browse)
        lay_load.addWidget(self._btn_browse)

        self._load_progress = QProgressBar()
        self._load_progress.setRange(0, 100)
        self._load_progress.setMaximumHeight(10)
        self._load_progress.setTextVisible(False)
        self._load_progress.setVisible(False)
        lay_load.addWidget(self._load_progress)

        self._info_lbl = QLabel("No TDEM data loaded.\nBrowse to a folder\ncontaining .AVG files.")
        self._info_lbl.setObjectName("InfoLabel")
        self._info_lbl.setWordWrap(True)
        self._info_lbl.setAlignment(Qt.AlignmentFlag.AlignTop)
        lay_load.addWidget(self._info_lbl)
        layout.addWidget(grp_load)

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
        lay_plot.addWidget(self._combo_plot)
        layout.addWidget(grp_plot)

        # ── Description ───────────────────────────────────────────────
        self._desc_lbl = QLabel("")
        self._desc_lbl.setWordWrap(True)
        self._desc_lbl.setObjectName("InfoLabel")
        self._desc_lbl.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._desc_lbl)

        # ── Actions ───────────────────────────────────────────────────
        grp_act, lay_act = make_group("Actions")
        self._btn_run    = icon_button("↻  Run / Refresh", "tdem",
                                       "Render the selected TDEM plot")
        self._btn_export = icon_button("⬆  Export…", "export",
                                       "Save figure to file")
        self._btn_run.clicked.connect(self._on_run)
        self._btn_export.clicked.connect(self._on_export)
        lay_act.addWidget(self._btn_run)
        lay_act.addWidget(self._btn_export)
        layout.addWidget(grp_act)

        self._status_lbl = QLabel("")
        self._status_lbl.setObjectName("InfoLabel")
        self._status_lbl.setWordWrap(True)
        layout.addWidget(self._status_lbl)

    # ── Content panel (right) ─────────────────────────────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        self._canvas = MplCanvas(self, toolbar=True)
        layout.addWidget(self._canvas)

    # ── Public API ────────────────────────────────────────────────────

    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._ctrl.dark = dark

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_browse(self) -> None:
        folder = QFileDialog.getExistingDirectory(
            self,
            "Select TDEM survey folder",
            "",
            QFileDialog.Option.ShowDirsOnly,
        )
        if not folder:
            return

        self._status_lbl.setText("Loading…")
        self._load_progress.setVisible(True)
        self._load_progress.setValue(0)
        self._btn_run.setEnabled(False)

        ok = self._ctrl.load_folder(folder, progress_cb=self._load_progress.setValue)

        self._load_progress.setVisible(False)
        self._btn_run.setEnabled(True)

        if ok:
            self._info_lbl.setText(self._ctrl.summary)
            self._status_lbl.setText("Data loaded — click Run to render.")
        else:
            self._info_lbl.setText(self._ctrl.summary or "No TDEM files found.")
            self._status_lbl.setText("Load failed.")

    def _on_category_changed(self, row: int) -> None:
        if row < 0 or row >= len(TDEM_GROUPS):
            return
        _label, plots = TDEM_GROUPS[row]
        self._combo_plot.blockSignals(True)
        self._combo_plot.clear()
        for label, _cls, _has_ax, _key in plots:
            self._combo_plot.addItem(label)
        self._combo_plot.blockSignals(False)
        self._combo_plot.setCurrentIndex(0)
        self._update_desc(row, 0)

    def _on_run(self) -> None:
        cat_row  = self._combo_category.currentIndex()
        plot_row = self._combo_plot.currentIndex()
        if cat_row < 0 or plot_row < 0:
            return
        _label, plots = TDEM_GROUPS[cat_row]
        if plot_row >= len(plots):
            return
        _plot_label, class_name, has_ax, data_key = plots[plot_row]

        self._status_lbl.setText(f"Running {class_name}…")
        self._btn_run.setEnabled(False)
        try:
            new_fig = self._ctrl.draw(class_name, has_ax, data_key, self._canvas.figure)
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
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog
        ExportDialog(figure=self._canvas.figure, parent=self).exec()

    # ── Helpers ───────────────────────────────────────────────────────

    def _populate_category_combo(self) -> None:
        self._combo_category.blockSignals(True)
        for group_label, _plots in TDEM_GROUPS:
            self._combo_category.addItem(group_label)
        self._combo_category.blockSignals(False)

    def _update_desc(self, cat_row: int, plot_row: int) -> None:
        try:
            _label, plots = TDEM_GROUPS[cat_row]
            plot_label, class_name, has_ax, data_key = plots[plot_row]
            self._desc_lbl.setText(
                f"<b>{plot_label}</b><br/>"
                f"<small style='color:#888'>tdem.{class_name}({data_key}).plot()</small>"
            )
        except Exception:
            self._desc_lbl.setText("")
