# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AdvancedToolsWindow — the Advanced Tools studio (desktop v2.6).

┌ Advanced Tools ───────────────────────────────────────────────────────────┐
│ Phase Tensor — ellipses, roses, skew…          26 stations  Export▾ Lib ▸ │
├────────────┬──────────────┬─────────────────────────────────┬─────────────┤
│ ANALYSES   │ PLOTS        │                                 │ PINNED      │
│ ▣ Strike   │ ● PT psection│   figure (or "why not" card)    │ [thumb]     │
│ ▣ Phase T. │ ● PT map     │                                 │             │
│ ▣ Induction│ ○ PT strip   │                                 │             │
│ ▣ Impedance│ OPTIONS      │                                 │             │
│ ▣ Depth    │ Colour map ▾ │                                 │             │
│ ▣ Survey   │ Station  ▾   │                                 │             │
│ UTILITIES  │ [↻ Draw]     │           [📌 Pin] [Export…]    │             │
│ ▣ Topo     │              │                                 │             │
│ ▣ Convert  │              │                                 │             │
└────────────┴──────────────┴─────────────────────────────────┴─────────────┘

The navigation rail picks a section: six emtools analysis sections (one
plot page) and two utilities with their own pages (Topography,
Conversion).  Each plot's options come from its signature
(:mod:`pycsamt.app.desktop.controllers.advanced_studio`), so only options
the plot takes are shown; plots that need a station or line grouping get
it from the Options box.  Figures are always publication-white; pinned
ones collect in the Library for side-by-side review and batch export.
"""

from __future__ import annotations

import io
from pathlib import Path

from PySide6.QtCore import QByteArray, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QIcon, QKeySequence, QPixmap, QShortcut
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSplitter,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
    QLineEdit,
)

from pycsamt.app.desktop.controllers.advanced_controller import (
    AdvancedController,
    ConversionController,
    ConversionWorker,
    DimModelWorker,
    TopoPreviewController,
    describe_advanced_plot,
)
from pycsamt.app.desktop.controllers.advanced_studio import (
    LINE_MODES,
    SECTIONS,
    is_polar,
    line_groups,
    plot_inputs,
    plot_kwargs,
    plot_options,
    site_names,
)
from pycsamt.app.desktop.controllers.correction_views import (
    figure_blank_reason,
)
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.widgets.compact_button import compact_button
from pycsamt.app.desktop.windows._base import _icon, icon_button, make_group
from pycsamt.app.desktop.windows._base import refresh_icons
from pycsamt.app.desktop.windows.inversion.forms import SettingsForm

_KEY_ROLE = Qt.ItemDataRole.UserRole
_MISSING = QColor("#8a94a3")


def _thumbnail(fig) -> QIcon:
    try:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=28, facecolor="white")
        pix = QPixmap()
        pix.loadFromData(buf.getvalue())
        return QIcon(pix.scaled(120, 80, Qt.AspectRatioMode.KeepAspectRatio,
                                Qt.TransformationMode.SmoothTransformation))
    except Exception:
        return QIcon()


def _caption(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setObjectName("InfoLabel")
    return lbl


class AdvancedToolsWindow(QWidget):
    """Advanced Tools studio: emtools analyses, topography, conversion."""

    # the user commits a conversion result as the main dataset
    conversion_committed = Signal(object)  # EDICollection
    panel_closed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        flags = (Qt.WindowType.Window | Qt.WindowType.WindowCloseButtonHint
                 | Qt.WindowType.WindowMinimizeButtonHint
                 | Qt.WindowType.WindowMaximizeButtonHint)
        super().__init__(parent, flags)
        self.setWindowTitle("pycsamt — Advanced Tools")
        ic = _icon("advanced-tools")
        if not ic.isNull():
            self.setWindowIcon(ic)
        self.resize(1320, 840)
        self._session_key = "advanced_tools"
        self._sites = None
        self._dark = False
        self._ctrl = AdvancedController()
        self._ctrl.dark = False  # figures are for publication: white
        self._topo_ctrl = TopoPreviewController()
        self._conv_ctrl = ConversionController()
        self._topo_fill_color = "#a89070"
        self._topo_line_color = "#6b4e2a"
        self._conv_worker: ConversionWorker | None = None
        self._conv_running = False
        self._conv_result_stats: list = []
        self._dim_worker: DimModelWorker | None = None
        self._section = 0
        self._plot_row: dict[int, int] = {}  # last plot per section
        self._forms: dict[str, SettingsForm] = {}
        self._figure = None
        self._figure_label = ""
        self._pinned: list = []  # (label, Figure)
        self._auto_rendered = False
        self._build_ui()
        self.select_section(0)

    # ══ UI ═══════════════════════════════════════════════════════════════
    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)
        root.addWidget(self._build_header())
        body = QSplitter(Qt.Orientation.Horizontal)
        body.setChildrenCollapsible(False)
        body.addWidget(self._build_nav())
        self._page_stack = QStackedWidget()
        self._page_stack.addWidget(self._build_plot_page())
        self._page_stack.addWidget(self._build_utility_page(
            self._build_topo_params_page(), self._build_topo_content()))
        self._page_stack.addWidget(self._build_utility_page(
            self._build_conv_params_page(), self._build_conv_content()))
        body.addWidget(self._page_stack)
        self._lib_panel = self._build_library()
        body.addWidget(self._lib_panel)
        body.setStretchFactor(1, 1)
        body.setSizes([190, 940, 190])
        root.addWidget(body, 1)
        QShortcut(QKeySequence("Ctrl+L"), self,
                  activated=lambda: self._btn_library.toggle())
        QShortcut(QKeySequence("Ctrl+R"), self, activated=self._on_run)
        QShortcut(QKeySequence("F5"), self, activated=self._on_run)

    def _build_header(self) -> QWidget:
        bar = QWidget()
        h = QHBoxLayout(bar)
        h.setContentsMargins(2, 0, 2, 0)
        h.setSpacing(8)
        self._title_lbl = QLabel("")
        self._title_lbl.setTextFormat(Qt.TextFormat.RichText)
        h.addWidget(self._title_lbl, 1)
        self._data_lbl = _caption("No EDI/XML data")
        h.addWidget(self._data_lbl)
        self._btn_export_menu = QToolButton()
        self._btn_export_menu.setText("Export  ▾")
        self._btn_export_menu.setPopupMode(
            QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QMenu(self._btn_export_menu)
        menu.addAction("Current figure…", self._on_export)
        menu.addAction("Pinned figures to a folder…", self._export_pinned)
        self._btn_export_menu.setMenu(menu)
        h.addWidget(self._btn_export_menu)
        self._btn_library = QToolButton()
        self._btn_library.setCheckable(True)
        self._btn_library.setChecked(True)
        self._btn_library.setText("Library ▸")
        self._btn_library.setToolTip("Show / hide pinned figures (Ctrl+L)")
        self._btn_library.toggled.connect(
            lambda on: self._lib_panel.setVisible(on))
        h.addWidget(self._btn_library)
        return bar

    def _build_nav(self) -> QWidget:
        w = QWidget()
        w.setMinimumWidth(170)
        w.setMaximumWidth(230)
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        v.addWidget(_caption("ANALYSES  ·  UTILITIES"))
        self._nav = QListWidget()
        self._nav.setObjectName("EngineList")
        self._nav.setIconSize(QSize(20, 20))
        for i, sec in enumerate(SECTIONS):
            if i and sec.key != "plots" and SECTIONS[i - 1].key == "plots":
                sep = QListWidgetItem("")
                sep.setFlags(Qt.ItemFlag.NoItemFlags)
                sep.setSizeHint(QSize(10, 10))
                self._nav.addItem(sep)
            item = QListWidgetItem(_icon(sec.icon), sec.label)
            item.setData(_KEY_ROLE, i)
            item.setToolTip(sec.help)
            item.setSizeHint(QSize(160, 30))
            self._nav.addItem(item)
        self._nav.currentItemChanged.connect(
            lambda cur, _p: cur is not None and cur.data(_KEY_ROLE)
            is not None and self.select_section(cur.data(_KEY_ROLE)))
        v.addWidget(self._nav, 1)
        return w

    def _build_plot_page(self) -> QWidget:
        page = QSplitter(Qt.Orientation.Horizontal)
        page.setChildrenCollapsible(False)
        left = QWidget()
        left.setMinimumWidth(230)
        left.setMaximumWidth(320)
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.setSpacing(4)
        lv.addWidget(_caption("PLOTS"))
        self._plot_list = QListWidget()
        self._plot_list.setObjectName("EngineList")
        self._plot_list.setStyleSheet("QListWidget::item { padding: 4px 3px; }")
        self._plot_list.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._plot_list.currentRowChanged.connect(self._on_plot_changed)
        self._plot_list.itemDoubleClicked.connect(lambda *_: self._on_run())
        lv.addWidget(self._plot_list, 3)
        self._desc_lbl = QLabel("")
        self._desc_lbl.setWordWrap(True)
        self._desc_lbl.setObjectName("InfoLabel")
        self._desc_lbl.setTextFormat(Qt.TextFormat.RichText)
        lv.addWidget(self._desc_lbl)

        # options: survey inputs + the plot's own options
        self._grp_opts, gl = make_group("Options")
        inputs = QFormLayout()
        inputs.setSpacing(5)
        self._station_combo = QComboBox()
        self._station_combo.setToolTip("Station for single-station plots")
        self._station_combo.activated.connect(lambda *_: self._on_run())
        self._station_row_lbl = QLabel("Station:")
        inputs.addRow(self._station_row_lbl, self._station_combo)
        self._lines_combo = QComboBox()
        for value, label in LINE_MODES:
            self._lines_combo.addItem(label, value)
        self._lines_combo.setToolTip("How stations are grouped into lines")
        self._lines_combo.activated.connect(self._on_lines_changed)
        self._lines_row_lbl = QLabel("Lines:")
        inputs.addRow(self._lines_row_lbl, self._lines_combo)
        self._lines_info = _caption("")
        self._lines_info.setWordWrap(True)
        inputs.addRow(self._lines_info)
        gl.addLayout(inputs)
        self._form_stack = QStackedWidget()
        self._no_opts = _caption("This plot has no options.")
        self._form_stack.addWidget(self._no_opts)
        gl.addWidget(self._form_stack)
        lv.addWidget(self._grp_opts)

        # dictionary model (ATOM pseudosection)
        grp_model, lay_model = make_group("Dictionary model")
        row = QFormLayout()
        self._spin_n_atoms = QDoubleSpinBox()
        self._spin_n_atoms.setDecimals(0)
        self._spin_n_atoms.setRange(2, 20)
        self._spin_n_atoms.setValue(6)
        row.addRow("Atoms:", self._spin_n_atoms)
        self._spin_n_iter = QDoubleSpinBox()
        self._spin_n_iter.setDecimals(0)
        self._spin_n_iter.setRange(10, 200)
        self._spin_n_iter.setSingleStep(10)
        self._spin_n_iter.setValue(40)
        row.addRow("Iterations:", self._spin_n_iter)
        lay_model.addLayout(row)
        self._btn_train_model = QPushButton("⚙  Train from survey")
        self._btn_train_model.setToolTip(
            "Learn a sparse-coding dictionary from the phase-tensor features\n"
            "of the loaded sites.  Required before plotting ATOM psection.")
        self._btn_train_model.clicked.connect(self._on_train_model)
        lay_model.addWidget(self._btn_train_model)
        self._model_status_lbl = _caption("Not trained")
        self._model_status_lbl.setWordWrap(True)
        lay_model.addWidget(self._model_status_lbl)
        self._grp_model = grp_model
        self._grp_model.setVisible(False)
        lv.addWidget(grp_model)

        run = QHBoxLayout()
        self._btn_run = QPushButton("↻  Draw")
        self._btn_run.setToolTip("Render the selected plot (Ctrl+R / F5)")
        self._btn_run.clicked.connect(self._on_run)
        run.addWidget(self._btn_run)
        self._chk_auto = QCheckBox("Auto")
        self._chk_auto.setChecked(True)
        self._chk_auto.setToolTip("Redraw when the plot or an option "
                                  "changes")
        run.addWidget(self._chk_auto)
        lv.addLayout(run)
        self._status_lbl = _caption("")
        self._status_lbl.setWordWrap(True)
        lv.addWidget(self._status_lbl)
        page.addWidget(left)

        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.setSpacing(2)
        self._canvas_view = CanvasResultView(
            right, toolbar=True, empty_title="No plot yet",
            empty_reason="Load EDI/XML data in the main window, then pick a "
                         "plot.")
        self._canvas = self._canvas_view.canvas
        self._canvas.set_refresh_callback(self._on_run,
                                          tooltip="Render selected plot")
        rv.addWidget(self._canvas_view, 1)
        tools = QHBoxLayout()
        tools.addStretch(1)
        self._btn_pin = QPushButton("📌  Pin to library")
        self._btn_pin.setObjectName("FileListBtn")
        self._btn_pin.clicked.connect(self._pin_current)
        self._btn_export = QPushButton("⬆  Export…")
        self._btn_export.setObjectName("FileListBtn")
        self._btn_export.clicked.connect(self._on_export)
        tools.addWidget(self._btn_pin)
        tools.addWidget(self._btn_export)
        rv.addLayout(tools)
        page.addWidget(right)
        page.setStretchFactor(1, 1)
        page.setSizes([260, 700])
        return page

    @staticmethod
    def _build_utility_page(params: QWidget, content: QWidget) -> QWidget:
        page = QSplitter(Qt.Orientation.Horizontal)
        page.setChildrenCollapsible(False)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(params)
        scroll.setMinimumWidth(280)
        scroll.setMaximumWidth(360)
        page.addWidget(scroll)
        page.addWidget(content)
        page.setStretchFactor(1, 1)
        page.setSizes([310, 700])
        return page

    def _build_topo_content(self) -> QWidget:
        page1 = QWidget()
        v1 = QVBoxLayout(page1)
        v1.setContentsMargins(0, 0, 0, 0)
        bar1 = QHBoxLayout()
        bar1.addWidget(QLabel("Preview:"))
        self._combo_topo_view = QComboBox()
        self._combo_topo_view.addItems(
            ["Elevation Profile", "Terrain Fill Preview",
             "Elevation Histogram"])
        self._combo_topo_view.currentIndexChanged.connect(
            self._on_topo_view_changed)
        bar1.addWidget(self._combo_topo_view)
        bar1.addStretch()
        self._topo_stats_lbl = _caption("")
        bar1.addWidget(self._topo_stats_lbl)
        v1.addLayout(bar1)
        self._canvas_topo_view = CanvasResultView(
            page1, toolbar=True, empty_title="No topography preview yet",
            empty_reason="Load survey data, then click Preview / Refresh.")
        self._canvas_topo = self._canvas_topo_view.canvas
        self._canvas_topo.set_refresh_callback(
            self._refresh_topo_preview, tooltip="Refresh the preview canvas")
        v1.addWidget(self._canvas_topo_view)
        return page1

    def _build_conv_content(self) -> QWidget:
        self._conv_tabs = QTabWidget()
        self._conv_tabs.setDocumentMode(True)
        results_page = QWidget()
        rp_v = QVBoxLayout(results_page)
        self._conv_table = QTableWidget()
        self._conv_table.setAlternatingRowColors(True)
        self._conv_table.horizontalHeader().setStretchLastSection(True)
        self._conv_table.setEditTriggers(
            QTableWidget.EditTrigger.NoEditTriggers)
        rp_v.addWidget(self._conv_table)
        self._conv_tabs.addTab(results_page, "Results")
        for attr, title, empty in (
                ("_canvas_conv_curves", "Impedance Curves",
                 "No impedance curves yet"),
                ("_canvas_conv_map", "Station Map", "No station map yet")):
            tab = QWidget()
            tv = QVBoxLayout(tab)
            view = CanvasResultView(
                tab, toolbar=True, empty_title=empty,
                empty_reason="Run a conversion to see it.")
            view.canvas.set_refresh_callback(self._on_conv_run,
                                             tooltip="Run the conversion")
            setattr(self, f"{attr}_view", view)
            setattr(self, attr, view.canvas)
            tv.addWidget(view)
            self._conv_tabs.addTab(tab, title)
        return self._conv_tabs

    def _build_library(self) -> QWidget:
        w = QWidget()
        w.setMinimumWidth(160)
        w.setMaximumWidth(250)
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        v.addWidget(_caption("PINNED FIGURES"))
        self._gallery = QListWidget()
        self._gallery.setObjectName("GalleryList")
        self._gallery.setViewMode(QListWidget.ViewMode.IconMode)
        self._gallery.setIconSize(QSize(120, 80))
        self._gallery.setResizeMode(QListWidget.ResizeMode.Adjust)
        self._gallery.setMovement(QListWidget.Movement.Static)
        self._gallery.setWordWrap(True)
        self._gallery.setSpacing(4)
        self._gallery.itemClicked.connect(self._show_pinned)
        v.addWidget(self._gallery, 1)
        row = QHBoxLayout()
        b = QPushButton("Remove")
        b.setObjectName("FileListBtn")
        b.clicked.connect(self._unpin)
        row.addWidget(b)
        b = QPushButton("Export all…")
        b.setObjectName("FileListBtn")
        b.clicked.connect(self._export_pinned)
        row.addWidget(b)
        v.addLayout(row)
        return w

    # ══ sections & plots ═════════════════════════════════════════════════
    def select_section(self, i: int) -> None:
        if not 0 <= i < len(SECTIONS):
            return
        if i != self._section and self._plot_list.currentRow() >= 0 \
                and SECTIONS[self._section].key == "plots":
            self._plot_row[self._section] = self._plot_list.currentRow()
        self._section = i
        sec = SECTIONS[i]
        for r in range(self._nav.count()):
            if self._nav.item(r).data(_KEY_ROLE) == i:
                self._nav.blockSignals(True)
                self._nav.setCurrentRow(r)
                self._nav.blockSignals(False)
        self._title_lbl.setText(
            f"<b>{sec.label}</b> <span style='color:#6b7280'>— "
            f"{sec.help}</span>")
        if sec.key == "topo":
            self._page_stack.setCurrentIndex(1)
            self._refresh_topo_preview()
            return
        if sec.key == "conv":
            self._page_stack.setCurrentIndex(2)
            return
        self._page_stack.setCurrentIndex(0)
        self._plot_list.blockSignals(True)
        self._plot_list.clear()
        for label, fn, _has_ax in sec.plots:
            tags = " · polar" if is_polar(fn) else ""
            item = QListWidgetItem(label)
            item.setData(_KEY_ROLE, fn)
            item.setToolTip(describe_advanced_plot(fn) + tags)
            self._plot_list.addItem(item)
        self._plot_list.blockSignals(False)
        self._plot_list.setCurrentRow(self._plot_row.get(i, 0))
        self._on_plot_changed(self._plot_list.currentRow())

    @property
    def current_plot(self) -> tuple:
        sec = SECTIONS[self._section]
        row = self._plot_list.currentRow()
        if sec.key != "plots" or not 0 <= row < len(sec.plots):
            return ()
        return sec.plots[row]

    def select_plot(self, fn_name: str) -> None:
        """Switch to the section holding *fn_name* and select it."""
        for i, sec in enumerate(SECTIONS):
            fns = [p[1] for p in sec.plots]
            if fn_name in fns:
                self._plot_row[i] = fns.index(fn_name)
                if i == self._section:
                    self._plot_list.setCurrentRow(fns.index(fn_name))
                else:
                    self.select_section(i)
                return

    def _form_for(self, fn: str) -> SettingsForm | None:
        if fn not in self._forms:
            fields = plot_options(fn)
            if not fields:
                return None
            form = SettingsForm(fields)
            form.changed.connect(self._on_option_changed)
            self._forms[fn] = form
            self._form_stack.addWidget(form)
        return self._forms[fn]

    def _on_plot_changed(self, _row: int) -> None:
        plot = self.current_plot
        if not plot:
            return
        label, fn, _has_ax = plot
        self._update_desc(label, fn)
        form = self._form_for(fn)
        self._form_stack.setCurrentWidget(form or self._no_opts)
        need = plot_inputs(fn)
        for w in (self._station_row_lbl, self._station_combo):
            w.setVisible("station" in need)
        for w in (self._lines_row_lbl, self._lines_combo, self._lines_info):
            w.setVisible("lines" in need)
        self._update_lines_info()
        self._grp_model.setVisible("model" in need)
        if self._chk_auto.isChecked():
            self._on_run()

    def _update_desc(self, label: str, fn: str) -> None:
        polar = " · polar" if is_polar(fn) else ""
        self._desc_lbl.setText(
            f"<b>{label}</b><span style='color:#888'>{polar}</span><br/>"
            f"<small style='color:#6b7280'>{describe_advanced_plot(fn)}"
            f"</small>")

    def _on_option_changed(self) -> None:
        if self._chk_auto.isChecked():
            QTimer.singleShot(0, self._on_run)

    def _lines(self) -> dict:
        return line_groups(site_names(self._ctrl._sites),
                           self._lines_combo.currentData() or "prefix")

    def _update_lines_info(self) -> None:
        g = self._lines()
        self._lines_info.setText(
            ", ".join(f"{k} ({len(v)})" for k, v in list(g.items())[:6])
            + (" …" if len(g) > 6 else "") if g else "")

    def _on_lines_changed(self, *_) -> None:
        self._update_lines_info()
        self._on_run()

    # ══ drawing ══════════════════════════════════════════════════════════
    def _on_run(self) -> None:
        plot = self.current_plot
        if not plot or not self._page_stack.currentIndex() == 0:
            return
        label, fn_name, has_ax = plot
        if self._ctrl._sites is None:
            self._status_lbl.setText("Load survey data first.")
            self._show_card(label, "Load EDI/XML data in the main window.")
            return
        need = plot_inputs(fn_name)
        form = self._forms.get(fn_name)
        kw = plot_kwargs(
            fn_name, form.values() if form else {},
            station=self._station_combo.currentText()
            if "station" in need else "",
            lines=self._lines() if "lines" in need else None)
        self._status_lbl.setText(f"Drawing {label}…")
        self._btn_run.setEnabled(False)
        self.setCursor(Qt.CursorShape.WaitCursor)
        import matplotlib.pyplot as plt

        target = plt.figure()
        try:
            new_fig = self._ctrl.draw(fn_name, has_ax, target, **kw)
        except Exception as exc:
            new_fig, target = None, None
            self._show_card(label, f"{fn_name} failed: {exc}")
            self._status_lbl.setText(f"Error: {exc}")
            return
        finally:
            self._btn_run.setEnabled(True)
            self.unsetCursor()
        fig = new_fig if new_fig is not None else target
        if new_fig is not None and target is not new_fig:
            plt.close(target)
        why = figure_blank_reason(fig)
        if why is not None:
            plt.close(fig)
            self._show_card(label, why.replace("error:", "error —"))
            self._status_lbl.setText("Nothing drawn — see the card.")
            return
        self._close_figure(self._figure)
        self._figure, self._figure_label = fig, label
        self._canvas.show_figure(fig)
        self._canvas_view.show_canvas()
        self._btn_pin.setEnabled(True)
        self._status_lbl.setText("Done.")

    def _show_card(self, title: str, reason: str) -> None:
        self._close_figure(self._figure)
        self._figure = None
        self._btn_pin.setEnabled(False)
        self._canvas_view.show_unavailable(title, reason)

    def _close_figure(self, fig) -> None:
        if fig is None or any(f is fig for _l, f in self._pinned):
            return
        import matplotlib.pyplot as plt

        plt.close(fig)

    def _on_export(self) -> None:
        if self._page_stack.currentIndex() == 1:
            fig = self._canvas_topo.figure
        elif self._figure is not None:
            fig = self._figure
        else:
            self._status_lbl.setText("No figure to export.")
            return
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog

        ExportDialog(figure=fig, parent=self).exec()

    def _auto_render_if_ready(self) -> None:
        if self._auto_rendered or self._ctrl._sites is None:
            return
        if not self.isVisible() or not self.current_plot:
            return
        self._auto_rendered = True
        QTimer.singleShot(0, self._on_run)

    # ══ library ══════════════════════════════════════════════════════════
    def _pin_current(self) -> None:
        if self._figure is None:
            return
        label = f"{SECTIONS[self._section].label} · {self._figure_label}"
        self._pinned.append((label, self._figure))
        item = QListWidgetItem(_thumbnail(self._figure), label)
        item.setToolTip(label)
        self._gallery.addItem(item)
        self._btn_library.setChecked(True)

    def _show_pinned(self, item: QListWidgetItem) -> None:
        i = self._gallery.row(item)
        if 0 <= i < len(self._pinned):
            label, fig = self._pinned[i]
            if self._page_stack.currentIndex() != 0:
                self._page_stack.setCurrentIndex(0)
            self._close_figure(self._figure)
            self._figure, self._figure_label = fig, label
            self._canvas.show_figure(fig)
            self._canvas_view.show_canvas()

    def _unpin(self) -> None:
        i = self._gallery.currentRow()
        if 0 <= i < len(self._pinned):
            _label, fig = self._pinned.pop(i)
            self._gallery.takeItem(i)
            if fig is not self._figure:
                self._close_figure(fig)

    def pinned_labels(self) -> list[str]:
        return [label for label, _f in self._pinned]

    def export_pinned(self, folder, fmt: str = "png", dpi: int = 300) -> list:
        """Save every pinned figure into *folder*; returns the paths."""
        out = []
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        for i, (label, fig) in enumerate(self._pinned, 1):
            stem = "".join(ch if ch.isalnum() else "_" for ch in label)
            path = folder / f"{i:02d}_{stem.strip('_')}.{fmt}"
            fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
            out.append(path)
        return out

    def _export_pinned(self) -> None:
        if not self._pinned:
            self._status_lbl.setText("Pin figures first (📌 under the plot).")
            return
        d = QFileDialog.getExistingDirectory(self, "Export pinned figures")
        if d:
            paths = self.export_pinned(d)
            self._status_lbl.setText(f"Exported {len(paths)} figure(s).")

    # ══ public API (main window) ═════════════════════════════════════════
    def set_sites(self, sites) -> None:
        self._sites = sites
        self._ctrl.set_sites(sites)
        self._topo_ctrl.set_sites(sites)
        self._auto_rendered = False
        names = site_names(sites)
        self._data_lbl.setText(f"{len(names)} EDI/XML stations" if names
                               else "No EDI/XML data")
        cur = self._station_combo.currentText()
        self._station_combo.clear()
        self._station_combo.addItems(names)
        if cur in names:
            self._station_combo.setCurrentText(cur)
        self._update_lines_info()
        if SECTIONS[self._section].key == "topo":
            self._refresh_topo_preview()
        elif self.isVisible() and self.current_plot:
            self._on_run()

    def set_dark_mode(self, dark: bool) -> None:
        # the UI follows the app theme; figures stay publication-white
        self._dark = dark
        refresh_icons(self, dark)
        self._ctrl.dark = False
        self._topo_ctrl.dark = False
        self._conv_ctrl.dark = False

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._auto_render_if_ready()

    def save_geometry_to(self, store: dict) -> None:
        store[self._session_key] = {
            "geometry": self.saveGeometry().toBase64().data().decode(),
            "visible": self.isVisible(),
            "library_visible": self._btn_library.isChecked(),
            "section": self._section,
            "auto": self._chk_auto.isChecked(),
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
        if "library_visible" in entry:
            self._btn_library.setChecked(bool(entry["library_visible"]))
        if "auto" in entry:
            self._chk_auto.setChecked(bool(entry["auto"]))
        sec = entry.get("section")
        if isinstance(sec, int) and 0 <= sec < len(SECTIONS):
            self.select_section(sec)

    def closeEvent(self, event) -> None:  # noqa: N802
        """Hide instead of destroying so state is preserved."""
        self.hide()
        event.ignore()
        self.panel_closed.emit()

    # ══ topography & conversion (unchanged from v2.5) ═══════════════════

    def _build_topo_params_page(self) -> QWidget:
        page = QWidget()
        vlay = QVBoxLayout(page)
        vlay.setContentsMargins(0, 0, 0, 0)
        vlay.setSpacing(4)

        # ── Configuration group ───────────────────────────────────────
        grp_cfg, lay_cfg = make_group("Configuration")

        self._chk_topo_enabled = QCheckBox("Enable topography")
        lay_cfg.addWidget(self._chk_topo_enabled)

        row_src = QWidget()
        h_src = QHBoxLayout(row_src)
        h_src.setContentsMargins(0, 0, 0, 0)
        h_src.addWidget(QLabel("Source"))
        self._combo_topo_source = QComboBox()
        self._combo_topo_source.addItems(["sites", "file", "array"])
        self._combo_topo_source.currentTextChanged.connect(
            self._on_topo_source_changed
        )
        h_src.addWidget(self._combo_topo_source)
        lay_cfg.addWidget(row_src)

        # File row (shown only when source == "file")
        self._topo_file_row = QWidget()
        h_file = QHBoxLayout(self._topo_file_row)
        h_file.setContentsMargins(0, 0, 0, 0)
        self._edit_topo_file = QLineEdit()
        self._edit_topo_file.setPlaceholderText("Elevation file path…")
        btn_browse_topo = QPushButton("…")
        compact_button(btn_browse_topo)
        btn_browse_topo.clicked.connect(self._browse_topo_file)
        h_file.addWidget(self._edit_topo_file)
        h_file.addWidget(btn_browse_topo)
        self._topo_file_row.setVisible(False)
        lay_cfg.addWidget(self._topo_file_row)

        row_interp = QWidget()
        h_interp = QHBoxLayout(row_interp)
        h_interp.setContentsMargins(0, 0, 0, 0)
        h_interp.addWidget(QLabel("Interpolation"))
        self._combo_topo_interp = QComboBox()
        self._combo_topo_interp.addItems(["linear", "cubic", "nearest"])
        h_interp.addWidget(self._combo_topo_interp)
        lay_cfg.addWidget(row_interp)

        row_exag = QWidget()
        h_exag = QHBoxLayout(row_exag)
        h_exag.setContentsMargins(0, 0, 0, 0)
        h_exag.addWidget(QLabel("Exaggeration"))
        self._spin_topo_exag = QDoubleSpinBox()
        self._spin_topo_exag.setRange(0.1, 20.0)
        self._spin_topo_exag.setSingleStep(0.1)
        self._spin_topo_exag.setValue(1.0)
        h_exag.addWidget(self._spin_topo_exag)
        lay_cfg.addWidget(row_exag)

        vlay.addWidget(grp_cfg)

        # ── Style group ───────────────────────────────────────────────
        grp_style, lay_style = make_group("Style")

        row_fc = QWidget()
        h_fc = QHBoxLayout(row_fc)
        h_fc.setContentsMargins(0, 0, 0, 0)
        h_fc.addWidget(QLabel("Fill color"))
        self._btn_topo_fill_color = QPushButton()
        self._btn_topo_fill_color.setFixedWidth(36)
        self._apply_topo_btn_color(
            self._btn_topo_fill_color, self._topo_fill_color
        )
        self._btn_topo_fill_color.clicked.connect(
            lambda: self._pick_color("fill")
        )
        h_fc.addWidget(self._btn_topo_fill_color)
        lay_style.addWidget(row_fc)

        row_fa = QWidget()
        h_fa = QHBoxLayout(row_fa)
        h_fa.setContentsMargins(0, 0, 0, 0)
        h_fa.addWidget(QLabel("Fill α"))
        self._spin_topo_fill_alpha = QDoubleSpinBox()
        self._spin_topo_fill_alpha.setRange(0.0, 1.0)
        self._spin_topo_fill_alpha.setSingleStep(0.05)
        self._spin_topo_fill_alpha.setValue(0.4)
        h_fa.addWidget(self._spin_topo_fill_alpha)
        lay_style.addWidget(row_fa)

        row_lc = QWidget()
        h_lc = QHBoxLayout(row_lc)
        h_lc.setContentsMargins(0, 0, 0, 0)
        h_lc.addWidget(QLabel("Line color"))
        self._btn_topo_line_color = QPushButton()
        self._btn_topo_line_color.setFixedWidth(36)
        self._apply_topo_btn_color(
            self._btn_topo_line_color, self._topo_line_color
        )
        self._btn_topo_line_color.clicked.connect(
            lambda: self._pick_color("line")
        )
        h_lc.addWidget(self._btn_topo_line_color)
        lay_style.addWidget(row_lc)

        row_lw = QWidget()
        h_lw = QHBoxLayout(row_lw)
        h_lw.setContentsMargins(0, 0, 0, 0)
        h_lw.addWidget(QLabel("Line width"))
        self._spin_topo_line_w = QDoubleSpinBox()
        self._spin_topo_line_w.setRange(0.1, 5.0)
        self._spin_topo_line_w.setSingleStep(0.1)
        self._spin_topo_line_w.setValue(1.2)
        h_lw.addWidget(self._spin_topo_line_w)
        lay_style.addWidget(row_lw)

        self._chk_surface_line = QCheckBox("Surface line")
        self._chk_surface_line.setChecked(True)
        lay_style.addWidget(self._chk_surface_line)

        self._chk_clip_below = QCheckBox("Clip below surface")
        self._chk_clip_below.setChecked(True)
        lay_style.addWidget(self._chk_clip_below)

        self._chk_pins_at_surface = QCheckBox("Station pins at surface")
        self._chk_pins_at_surface.setChecked(True)
        lay_style.addWidget(self._chk_pins_at_surface)

        vlay.addWidget(grp_style)

        # ── Pseudosection strip group ─────────────────────────────────
        grp_strip, lay_strip = make_group("Pseudosection Strip")

        self._chk_topo_strip = QCheckBox("Show topo strip")
        self._chk_topo_strip.setChecked(True)
        lay_strip.addWidget(self._chk_topo_strip)

        row_sh = QWidget()
        h_sh = QHBoxLayout(row_sh)
        h_sh.setContentsMargins(0, 0, 0, 0)
        h_sh.addWidget(QLabel("Strip height"))
        self._spin_strip_h = QDoubleSpinBox()
        self._spin_strip_h.setRange(0.05, 0.5)
        self._spin_strip_h.setSingleStep(0.02)
        self._spin_strip_h.setValue(0.18)
        h_sh.addWidget(self._spin_strip_h)
        lay_strip.addWidget(row_sh)

        vlay.addWidget(grp_strip)

        # ── Actions group ─────────────────────────────────────────────
        grp_ta, lay_ta = make_group("Actions")
        self._btn_topo_apply = icon_button(
            "↻  Apply to Global Config", "", "Apply settings to PYCSAMT_TOPO"
        )
        self._btn_topo_reset = icon_button(
            "↩  Reset to Defaults", "", "Reset PYCSAMT_TOPO to defaults"
        )
        self._btn_topo_preview = icon_button(
            "▶  Preview / Refresh", "", "Refresh the preview canvas"
        )
        self._btn_topo_apply.clicked.connect(self._on_topo_apply)
        self._btn_topo_reset.clicked.connect(self._on_topo_reset)
        self._btn_topo_preview.clicked.connect(self._refresh_topo_preview)
        lay_ta.addWidget(self._btn_topo_apply)
        lay_ta.addWidget(self._btn_topo_reset)
        lay_ta.addWidget(self._btn_topo_preview)
        vlay.addWidget(grp_ta)

        self._topo_status = QLabel("")
        self._topo_status.setObjectName("InfoLabel")
        self._topo_status.setWordWrap(True)
        vlay.addWidget(self._topo_status)

        vlay.addStretch(1)
        return page

    # ── Conversion params page (page 2) ──────────────────────────────

    def _build_conv_params_page(self) -> QWidget:
        page = QWidget()
        vlay = QVBoxLayout(page)
        vlay.setContentsMargins(0, 0, 0, 0)
        vlay.setSpacing(4)

        # ── Input group ───────────────────────────────────────────────
        grp_inp, lay_inp = make_group("Input")

        row_type = QWidget()
        h_type = QHBoxLayout(row_type)
        h_type.setContentsMargins(0, 0, 0, 0)
        h_type.addWidget(QLabel("Type"))
        self._combo_conv_type = QComboBox()
        self._combo_conv_type.addItems(
            ["AVG → EDI", "J → EDI", "Spectra → EDI"]
        )
        self._combo_conv_type.currentIndexChanged.connect(
            self._on_conv_type_changed
        )
        h_type.addWidget(self._combo_conv_type)
        lay_inp.addWidget(row_type)

        self._conv_file_row = QWidget()
        h_cfile = QHBoxLayout(self._conv_file_row)
        h_cfile.setContentsMargins(0, 0, 0, 0)
        self._edit_conv_path = QLineEdit()
        self._edit_conv_path.setPlaceholderText("File or directory path…")
        self._edit_conv_path.textChanged.connect(self._update_conv_run_state)
        btn_browse_conv = QPushButton("📂")
        compact_button(btn_browse_conv)
        btn_browse_conv.clicked.connect(self._browse_conv_path)
        h_cfile.addWidget(self._edit_conv_path)
        h_cfile.addWidget(btn_browse_conv)
        lay_inp.addWidget(self._conv_file_row)

        self._conv_dir_hint = QLabel("Select file or directory")
        self._conv_dir_hint.setObjectName("InfoLabel")
        lay_inp.addWidget(self._conv_dir_hint)

        vlay.addWidget(grp_inp)

        # ── Options group (stacked per type) ──────────────────────────
        grp_opts, lay_opts = make_group("Options")
        self._conv_opts_stack = QStackedWidget()

        # Page 0 — AVG opts
        avg_page = QWidget()
        avg_form = QFormLayout(avg_page)
        avg_form.setSpacing(5)
        self._avg_freq_order = QComboBox()
        self._avg_freq_order.addItems(["ascending", "descending"])
        avg_form.addRow("Freq order:", self._avg_freq_order)
        self._avg_freq_tol = QDoubleSpinBox()
        self._avg_freq_tol.setDecimals(9)
        self._avg_freq_tol.setRange(1e-9, 1e-2)
        self._avg_freq_tol.setSingleStep(1e-6)
        self._avg_freq_tol.setValue(1e-5)
        avg_form.addRow("Freq tol:", self._avg_freq_tol)
        self._avg_compute_z = QCheckBox()
        avg_form.addRow("Compute Z from ρ/φ:", self._avg_compute_z)
        self._avg_compute_rho = QCheckBox()
        avg_form.addRow("Compute ρ/φ from Z:", self._avg_compute_rho)
        self._conv_opts_stack.addWidget(avg_page)

        # Page 1 — J opts
        j_page = QWidget()
        j_form = QFormLayout(j_page)
        j_form.setSpacing(5)
        self._j_freq_order = QComboBox()
        self._j_freq_order.addItems(["ascending", "descending"])
        j_form.addRow("Freq order:", self._j_freq_order)
        self._j_station_name = QLineEdit()
        self._j_station_name.setPlaceholderText("Optional single-station name")
        j_form.addRow("Station name:", self._j_station_name)
        self._conv_opts_stack.addWidget(j_page)

        # Page 2 — Spectra opts
        sp_page = QWidget()
        sp_form = QFormLayout(sp_page)
        sp_form.setSpacing(5)
        self._sp_e_labels = QLineEdit("EX,EY")
        sp_form.addRow("E labels:", self._sp_e_labels)
        self._sp_h_labels = QLineEdit("HX,HY")
        sp_form.addRow("H labels:", self._sp_h_labels)
        self._sp_estimate_errors = QCheckBox()
        sp_form.addRow("Estimate errors:", self._sp_estimate_errors)
        self._sp_remote_ref = QCheckBox()
        sp_form.addRow("Use remote ref:", self._sp_remote_ref)
        self._sp_station_suffix = QLineEdit()
        sp_form.addRow("Station suffix:", self._sp_station_suffix)
        self._sp_skip_errors = QCheckBox()
        self._sp_skip_errors.setChecked(True)
        sp_form.addRow("Skip errors:", self._sp_skip_errors)
        self._conv_opts_stack.addWidget(sp_page)

        lay_opts.addWidget(self._conv_opts_stack)
        vlay.addWidget(grp_opts)

        # ── AVG topography / station profile group ───────────────────
        self._avg_topo_group, lay_avg_topo = make_group(
            "Add topography / station profile"
        )
        topo_hint = QLabel(
            "Optional: attach a Zonge .stn station profile so converted EDI "
            "files receive station elevation, longitude, and latitude."
        )
        topo_hint.setObjectName("InfoLabel")
        topo_hint.setWordWrap(True)
        lay_avg_topo.addWidget(topo_hint)

        row_stn = QWidget()
        h_stn = QHBoxLayout(row_stn)
        h_stn.setContentsMargins(0, 0, 0, 0)
        h_stn.addWidget(QLabel("File"))
        self._avg_stn_path = QLineEdit()
        self._avg_stn_path.setPlaceholderText(
            "K1.stn or station profile file..."
        )
        btn_stn = QPushButton("Browse...")
        btn_stn.setFixedWidth(76)
        btn_stn.clicked.connect(self._browse_avg_stn_path)
        h_stn.addWidget(self._avg_stn_path)
        h_stn.addWidget(btn_stn)
        lay_avg_topo.addWidget(row_stn)

        avg_topo_form = QFormLayout()
        avg_topo_form.setSpacing(5)
        self._avg_convert_coords = QCheckBox()
        self._avg_convert_coords.setChecked(True)
        avg_topo_form.addRow(
            "Convert UTM → lat/lon:", self._avg_convert_coords
        )
        self._avg_epsg = QLineEdit()
        self._avg_epsg.setPlaceholderText("EPSG code, e.g. 32650")
        avg_topo_form.addRow("EPSG:", self._avg_epsg)
        self._avg_utm_zone = QLineEdit()
        self._avg_utm_zone.setPlaceholderText("UTM zone, e.g. 50N")
        avg_topo_form.addRow("UTM zone:", self._avg_utm_zone)
        lay_avg_topo.addLayout(avg_topo_form)
        coord_hint = QLabel(
            "If conversion is enabled, provide either EPSG or UTM zone. "
            "Disable it only when the profile already contains latitude and "
            "longitude columns."
        )
        coord_hint.setObjectName("InfoLabel")
        coord_hint.setWordWrap(True)
        lay_avg_topo.addWidget(coord_hint)
        vlay.addWidget(self._avg_topo_group)

        # ── Output group (optional) ───────────────────────────────────
        grp_out, lay_out = make_group("Output")
        self._chk_write_edis = QCheckBox("Write EDI files to disk")
        lay_out.addWidget(self._chk_write_edis)
        row_outdir = QWidget()
        h_outdir = QHBoxLayout(row_outdir)
        h_outdir.setContentsMargins(0, 0, 0, 0)
        h_outdir.addWidget(QLabel("Output dir"))
        self._edit_out_dir = QLineEdit()
        self._edit_out_dir.setPlaceholderText("Output directory…")
        btn_browse_out = QPushButton("…")
        compact_button(btn_browse_out)
        btn_browse_out.clicked.connect(self._browse_out_dir)
        h_outdir.addWidget(self._edit_out_dir)
        h_outdir.addWidget(btn_browse_out)
        lay_out.addWidget(row_outdir)
        vlay.addWidget(grp_out)

        # ── Actions group ─────────────────────────────────────────────
        grp_ca, lay_ca = make_group("Actions")

        self._btn_conv_run = icon_button(
            "▶▶  Run Transform", "", "Run the conversion"
        )
        self._btn_conv_run.setObjectName("CommitButton")
        self._btn_conv_run.clicked.connect(self._on_conv_run)

        self._btn_conv_commit = icon_button(
            "⬆  Commit as Main Dataset",
            "",
            "Use converted data as main dataset",
        )
        self._btn_conv_commit.setObjectName("CommitButton")
        self._btn_conv_commit.setEnabled(False)
        self._btn_conv_commit.clicked.connect(self._on_conv_commit)

        self._btn_conv_export = icon_button(
            "💾  Export EDIs…", "", "Write EDI files to disk"
        )
        self._btn_conv_export.setEnabled(False)
        self._btn_conv_export.clicked.connect(self._on_conv_export)

        self._btn_conv_clear = icon_button(
            "✗  Clear", "", "Clear conversion result"
        )
        self._btn_conv_clear.clicked.connect(self._on_conv_clear)

        lay_ca.addWidget(self._btn_conv_run)
        lay_ca.addWidget(self._btn_conv_commit)
        lay_ca.addWidget(self._btn_conv_export)
        lay_ca.addWidget(self._btn_conv_clear)
        vlay.addWidget(grp_ca)
        self._update_conv_run_state()

        self._conv_progress = QProgressBar()
        self._conv_progress.setRange(0, 0)  # indeterminate
        self._conv_progress.setVisible(False)
        vlay.addWidget(self._conv_progress)

        self._conv_status = QLabel("")
        self._conv_status.setObjectName("InfoLabel")
        self._conv_status.setWordWrap(True)
        vlay.addWidget(self._conv_status)

        vlay.addStretch(1)
        return page

    # ── Train-model slots ─────────────────────────────────────────────

    def _on_train_model(self) -> None:
        if self._ctrl._sites is None:
            self._model_status_lbl.setText("Load survey data first.")
            return
        n_atoms = int(self._spin_n_atoms.value())
        n_iter = int(self._spin_n_iter.value())
        self._btn_train_model.setEnabled(False)
        self._model_status_lbl.setText("Training…")
        self._dim_worker = DimModelWorker(self._ctrl, n_atoms, n_iter)
        self._dim_worker.finished.connect(self._on_model_trained)
        self._dim_worker.error.connect(self._on_model_train_error)
        self._dim_worker.start()

    def _on_model_trained(self, model: dict) -> None:
        self._btn_train_model.setEnabled(True)
        meta = model.get("meta", {})
        n_atoms = model["D"].shape[1] if model.get("D") is not None else "?"
        n_samples = meta.get("samples", "?")
        self._model_status_lbl.setText(
            f"Ready: {n_atoms} atoms · {n_samples} samples"
        )
        self._status_lbl.setText("Model trained.")
        self._on_run()

    def _on_model_train_error(self, msg: str) -> None:
        self._btn_train_model.setEnabled(True)
        self._model_status_lbl.setText(f"Error: {msg}")

    # ── Topo slots ────────────────────────────────────────────────────

    def _on_topo_source_changed(self, text: str) -> None:
        self._topo_file_row.setVisible(text == "file")

    def _browse_topo_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select elevation file",
            "",
            "CSV / Parquet (*.csv *.parquet *.txt);;All Files (*)",
        )
        if path:
            self._edit_topo_file.setText(path)

    def _pick_color(self, which: str) -> None:
        current = (
            self._topo_fill_color if which == "fill" else self._topo_line_color
        )
        color = QColorDialog.getColor(current, self, f"Choose {which} color")
        if color.isValid():
            hex_color = color.name()
            if which == "fill":
                self._topo_fill_color = hex_color
                self._apply_topo_btn_color(
                    self._btn_topo_fill_color, hex_color
                )
            else:
                self._topo_line_color = hex_color
                self._apply_topo_btn_color(
                    self._btn_topo_line_color, hex_color
                )

    def _apply_topo_btn_color(self, btn: QPushButton, hex_color: str) -> None:
        btn.setStyleSheet(f"background:{hex_color}; border:1px solid #888")

    def _on_topo_apply(self) -> None:
        try:
            from pycsamt.topo.config import configure_topo

            configure_topo(
                enabled=self._chk_topo_enabled.isChecked(),
                source=self._combo_topo_source.currentText(),
                elev_file=(self._edit_topo_file.text().strip() or None),
                interp_method=self._combo_topo_interp.currentText(),
                exaggeration=self._spin_topo_exag.value(),
                fill_color=self._topo_fill_color,
                fill_alpha=self._spin_topo_fill_alpha.value(),
                line_color=self._topo_line_color,
                line_width=self._spin_topo_line_w.value(),
                show_surface_line=self._chk_surface_line.isChecked(),
                clip_below_surface=self._chk_clip_below.isChecked(),
                station_pins_at_surface=self._chk_pins_at_surface.isChecked(),
                show_topo_strip=self._chk_topo_strip.isChecked(),
                strip_height_ratio=self._spin_strip_h.value(),
            )
            from pycsamt.topo.config import PYCSAMT_TOPO

            self._topo_status.setText(PYCSAMT_TOPO.summary())
        except Exception as exc:
            self._topo_status.setText(f"Error: {exc}")
        self._refresh_topo_preview()

    def _on_topo_reset(self) -> None:
        try:
            from pycsamt.topo.config import reset_topo

            reset_topo()
            self._sync_topo_widgets_from_config()
            from pycsamt.topo.config import PYCSAMT_TOPO

            self._topo_status.setText("Reset. " + PYCSAMT_TOPO.summary())
        except Exception as exc:
            self._topo_status.setText(f"Error: {exc}")
        self._refresh_topo_preview()

    def _sync_topo_widgets_from_config(self) -> None:
        try:
            from pycsamt.topo.config import PYCSAMT_TOPO as T

            self._chk_topo_enabled.setChecked(T.enabled)
            idx = self._combo_topo_source.findText(T.source)
            if idx >= 0:
                self._combo_topo_source.setCurrentIndex(idx)
            self._topo_file_row.setVisible(T.source == "file")
            if T.elev_file:
                self._edit_topo_file.setText(T.elev_file)
            idx_interp = self._combo_topo_interp.findText(T.interp_method)
            if idx_interp >= 0:
                self._combo_topo_interp.setCurrentIndex(idx_interp)
            self._spin_topo_exag.setValue(T.exaggeration)
            self._topo_fill_color = T.fill_color
            self._apply_topo_btn_color(self._btn_topo_fill_color, T.fill_color)
            self._spin_topo_fill_alpha.setValue(T.fill_alpha)
            self._topo_line_color = T.line_color
            self._apply_topo_btn_color(self._btn_topo_line_color, T.line_color)
            self._spin_topo_line_w.setValue(T.line_width)
            self._chk_surface_line.setChecked(T.show_surface_line)
            self._chk_clip_below.setChecked(T.clip_below_surface)
            self._chk_pins_at_surface.setChecked(T.station_pins_at_surface)
            self._chk_topo_strip.setChecked(T.show_topo_strip)
            self._spin_strip_h.setValue(T.strip_height_ratio)
        except Exception:
            pass

    def _refresh_topo_preview(self) -> None:
        if self._topo_ctrl._sites is None:
            self._topo_stats_lbl.setText("No data loaded")
        else:
            try:
                stats = self._topo_ctrl.get_stats()
                self._topo_stats_lbl.setText(
                    f"{stats['n_stations']} stations  |  "
                    f"elev {stats['elev_min']:.0f}–{stats['elev_max']:.0f} m"
                    if stats["n_stations"] > 0
                    else "No elevation data"
                )
            except Exception:
                self._topo_stats_lbl.setText("")

        view_idx = self._combo_topo_view.currentIndex()
        fig = self._canvas_topo.figure
        try:
            if view_idx == 0:
                self._topo_ctrl.plot_elevation_profile(fig)
            elif view_idx == 1:
                self._topo_ctrl.plot_fill_preview(fig)
            else:
                self._topo_ctrl.plot_elevation_histogram(fig)
            self._canvas_topo.draw()
            self._canvas_topo_view.show_canvas()
        except Exception as exc:
            self._canvas_topo_view.show_unavailable(
                "Preview unavailable", f"Preview error: {exc}"
            )

    def _on_topo_view_changed(self, _: int) -> None:
        self._refresh_topo_preview()

    # ── Conversion slots ──────────────────────────────────────────────

    def _on_conv_type_changed(self, idx: int) -> None:
        self._conv_opts_stack.setCurrentIndex(idx)
        if hasattr(self, "_avg_topo_group"):
            self._avg_topo_group.setVisible(idx == 0)
        self._update_conv_run_state()

    def _update_conv_run_state(self) -> None:
        if not hasattr(self, "_btn_conv_run"):
            return
        has_source = bool(
            getattr(self, "_edit_conv_path", None)
            and self._edit_conv_path.text().strip()
        )
        self._btn_conv_run.setEnabled(has_source and not self._conv_running)

    def _browse_conv_path(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select input file", "", "All Files (*)"
        )
        if not path:
            # fallback to directory
            path = QFileDialog.getExistingDirectory(
                self, "Select input directory", ""
            )
        if path:
            self._edit_conv_path.setText(path)

    def _browse_out_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select output directory", ""
        )
        if path:
            self._edit_out_dir.setText(path)

    def _browse_avg_stn_path(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select AVG topography / station profile",
            "",
            "Station files (*.stn *.txt *.csv);;All Files (*)",
        )
        if path:
            self._avg_stn_path.setText(path)

    def _on_conv_run(self) -> None:
        path = self._edit_conv_path.text().strip()
        if not path:
            self._conv_status.setText(
                "Please provide an input file/directory."
            )
            return

        type_str = self._combo_conv_type.currentText()
        self._conv_ctrl.set_source(type_str, path)

        # Collect options from the active options page
        options: dict = {}
        idx = self._combo_conv_type.currentIndex()
        try:
            if idx == 0:  # AVG
                options["freq_order"] = self._avg_freq_order.currentText()
                options["freq_tol"] = self._avg_freq_tol.value()
                if self._avg_compute_z.isChecked():
                    options["compute_z"] = True
                if self._avg_compute_rho.isChecked():
                    options["compute_rho_phi"] = True
                stn_path = self._avg_stn_path.text().strip()
                if stn_path:
                    options["stn_path"] = stn_path
                    options["convert_stn_coords"] = (
                        self._avg_convert_coords.isChecked()
                    )
                epsg = self._avg_epsg.text().strip()
                if epsg:
                    options["epsg"] = epsg
                utm_zone = self._avg_utm_zone.text().strip()
                if utm_zone:
                    options["utm_zone"] = utm_zone
            elif idx == 1:  # J
                options["freq_order"] = self._j_freq_order.currentText()
                station_name = self._j_station_name.text().strip()
                if station_name:
                    options["name"] = station_name
            elif idx == 2:  # Spectra
                options["e_labels"] = self._sp_e_labels.text()
                options["h_labels"] = self._sp_h_labels.text()
                options["estimate_error"] = (
                    self._sp_estimate_errors.isChecked()
                )
                options["use_remote"] = self._sp_remote_ref.isChecked()
                suf = self._sp_station_suffix.text().strip()
                if suf:
                    options["station_suffix"] = suf
                options["skip_errors"] = self._sp_skip_errors.isChecked()
        except Exception:
            pass

        if self._chk_write_edis.isChecked():
            out_dir = self._edit_out_dir.text().strip()
            if out_dir:
                options["output_dir"] = out_dir

        self._conv_progress.setVisible(True)
        self._conv_running = True
        self._update_conv_run_state()
        self._conv_status.setText(f"Running {type_str}…")

        self._conv_worker = ConversionWorker(self._conv_ctrl, options)
        self._conv_worker.finished.connect(self._on_conv_finished)
        self._conv_worker.error.connect(self._on_conv_error)
        self._conv_worker.start()

    def _on_conv_finished(self, collection, failures: list) -> None:
        self._conv_progress.setVisible(False)
        self._conv_running = False
        self._update_conv_run_state()

        self._conv_ctrl._result = collection
        stats = self._conv_ctrl.build_stats(collection, failures)
        self._conv_ctrl._stats = stats
        self._conv_result_stats = stats.get("rows", [])

        # Populate table
        rows_data = stats.get("rows", [])
        cols = [
            "Station",
            "N Freqs",
            "F min (Hz)",
            "F max (Hz)",
            "Latitude",
            "Longitude",
            "Elevation",
            "Has Z",
            "Has Tipper",
        ]
        self._conv_table.setColumnCount(len(cols))
        self._conv_table.setHorizontalHeaderLabels(cols)
        self._conv_table.setRowCount(len(rows_data))

        def _fmt_coord(value, digits=6):
            try:
                value = float(value)
                if value == value:
                    return f"{value:.{digits}f}"
            except Exception:
                pass
            return "—"

        for r, row in enumerate(rows_data):
            self._conv_table.setItem(r, 0, QTableWidgetItem(row["station"]))
            self._conv_table.setItem(
                r, 1, QTableWidgetItem(str(row["n_freqs"]))
            )
            self._conv_table.setItem(
                r,
                2,
                QTableWidgetItem(
                    f"{row['f_min']:.4g}"
                    if row["f_min"] == row["f_min"]
                    else "—"
                ),
            )
            self._conv_table.setItem(
                r,
                3,
                QTableWidgetItem(
                    f"{row['f_max']:.4g}"
                    if row["f_max"] == row["f_max"]
                    else "—"
                ),
            )
            self._conv_table.setItem(
                r, 4, QTableWidgetItem(_fmt_coord(row.get("lat")))
            )
            self._conv_table.setItem(
                r, 5, QTableWidgetItem(_fmt_coord(row.get("lon")))
            )
            self._conv_table.setItem(
                r, 6, QTableWidgetItem(_fmt_coord(row.get("elev"), digits=2))
            )
            self._conv_table.setItem(
                r, 7, QTableWidgetItem("Yes" if row["has_Z"] else "No")
            )
            self._conv_table.setItem(
                r, 8, QTableWidgetItem("Yes" if row["has_tipper"] else "No")
            )

        # Refresh plots
        try:
            self._conv_ctrl.plot_impedance_curves(
                self._canvas_conv_curves.figure
            )
            self._canvas_conv_curves.draw()
            self._canvas_conv_curves_view.show_canvas()
        except Exception as exc:
            self._canvas_conv_curves_view.show_unavailable(
                "Impedance curves unavailable", str(exc)
            )
        try:
            self._conv_ctrl.plot_station_map(self._canvas_conv_map.figure)
            self._canvas_conv_map.draw()
            self._canvas_conv_map_view.show_canvas()
        except Exception as exc:
            self._canvas_conv_map_view.show_unavailable(
                "Station map unavailable", str(exc)
            )

        n_ok = stats.get("n_total", 0)
        n_fail = stats.get("n_failures", 0)
        self._conv_status.setText(
            f"Done: {n_ok} stations"
            + (f", {n_fail} failures" if n_fail else "")
        )
        self._btn_conv_commit.setEnabled(True)
        self._btn_conv_export.setEnabled(True)

    def _on_conv_error(self, msg: str) -> None:
        self._conv_progress.setVisible(False)
        self._conv_running = False
        self._update_conv_run_state()
        self._conv_status.setText(f"Error: {msg}")

    def _on_conv_commit(self) -> None:
        if self._conv_ctrl.result is not None:
            self.conversion_committed.emit(self._conv_ctrl.result)

    def _on_conv_export(self) -> None:
        if not self._conv_ctrl.has_result:
            return
        out_dir = QFileDialog.getExistingDirectory(
            self, "Select output directory for EDI files", ""
        )
        if not out_dir:
            return
        try:
            from pycsamt.emtools._core import _iter_items

            written = 0
            for ed in _iter_items(self._conv_ctrl.result):
                try:
                    write_fn = getattr(ed, "write_edifile", None) or getattr(
                        ed, "write", None
                    )
                    if write_fn is not None:
                        write_fn(save_dir=out_dir)
                        written += 1
                except Exception:
                    pass
            self._conv_status.setText(
                f"Exported {written} EDI files to {out_dir}"
            )
        except Exception as exc:
            self._conv_status.setText(f"Export error: {exc}")

    def _on_conv_clear(self) -> None:
        self._conv_ctrl._result = None
        self._conv_ctrl._stats = {}
        self._conv_ctrl._failures = []
        self._conv_result_stats = []
        self._conv_table.clearContents()
        self._conv_table.setRowCount(0)
        self._conv_table.setColumnCount(0)
        # Clear canvases
        for canvas in (self._canvas_conv_curves, self._canvas_conv_map):
            canvas.figure.clear()
            canvas.draw()
        self._canvas_conv_curves_view.show_unavailable(
            "No impedance curves yet", "Run a conversion to see impedance curves."
        )
        self._canvas_conv_map_view.show_unavailable(
            "No station map yet", "Run a conversion to see the station map."
        )
        self._btn_conv_commit.setEnabled(False)
        self._btn_conv_export.setEnabled(False)
        self._conv_status.setText("Cleared.")
        self._conv_progress.setVisible(False)



__all__ = ["AdvancedToolsWindow"]
