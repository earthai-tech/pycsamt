# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AirborneWindow — ZTEM / AFMAG / MobileMT diagnostics and map plotting.

Net-new floating panel (Phase 8 of the desktop modernization plan).
Loads an :class:`~pycsamt.airborne.site.AirborneSites` from an EMTF-XML
file or directory via :func:`~pycsamt.airborne.site.ensure_asites`, then
dispatches Category -> Plot combos into
:class:`~pycsamt.app.desktop.controllers.airborne_controller
.AirborneController`.

MobileMT is generic-adapter-only: only already-decoded EMTF-XML is ever
read (raw vendor MobileMT files are a permanent restriction, not a gap
this window works around — see the controller module docstring and
``MEMORY.md``'s ``project_mobilemt_vendor_data_blocked.md``). The "Load"
button and the MobileMT category description say so explicitly rather
than implying raw-format import is supported or forthcoming.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.airborne_controller import (
    CATALOGUE,
    CATEGORIES,
    AirborneController,
    ParamSpec,
)
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    icon_button,
    make_group,
)

_CATEGORY_NOTE = {
    "ZTEM": "Along-profile divergence, phase rotation, and usable-band "
    "diagnostics for Geotech ZTEM-style tipper data.",
    "AFMAG": "Tilt-angle profiles/pseudosections and aircraft-motion "
    "susceptibility diagnostics.",
    "MobileMT": "Generic-adapter-only: reads already-decoded EMTF-XML. "
    "pyCSAMT never reads raw vendor MobileMT files — this is a "
    "permanent restriction, not a temporary gap.",
}


class AirborneWindow(PanelWindow):
    """Floating ZTEM / AFMAG / MobileMT diagnostics panel."""

    def __init__(self, parent: QWidget | None = None) -> None:
        self._ctrl = AirborneController()
        super().__init__(
            title="Airborne EM",
            session_key="airborne_window",
            params_width=300,
            icon_name="induction",
            parent=parent,
        )
        self.resize(1100, 720)
        self._populate_category_combo()
        self._on_category_changed(0)

    # =========================================================================
    # Left panel
    # =========================================================================

    def _build_params(self, layout: QVBoxLayout) -> None:
        # ── Data ─────────────────────────────────────────────────────────
        grp_data, lay_data = make_group("Data")
        self._data_status = QLabel("No airborne data loaded")
        self._data_status.setObjectName("InfoLabel")
        self._data_status.setWordWrap(True)
        lay_data.addWidget(self._data_status)

        row = QHBoxLayout()
        btn_load = QPushButton("Load EMTF-XML…")
        btn_load.setToolTip(
            "Select a directory of EMTF-XML files (one per station/"
            "line) or a single EMTF-XML file. Raw vendor MobileMT "
            "files are never read — see the MobileMT category note."
        )
        btn_load.clicked.connect(self._on_load)
        btn_clear = QPushButton("Clear")
        btn_clear.clicked.connect(self._on_clear)
        row.addWidget(btn_load)
        row.addWidget(btn_clear)
        lay_data.addLayout(row)
        layout.addWidget(grp_data)

        # ── Category / plot navigation ──────────────────────────────────
        grp_nav, lay_nav = make_group("Diagnostic")
        self._combo_category = QComboBox()
        self._combo_category.currentIndexChanged.connect(
            self._on_category_changed
        )
        lay_nav.addWidget(QLabel("Category:"))
        lay_nav.addWidget(self._combo_category)

        self._category_note = QLabel("")
        self._category_note.setWordWrap(True)
        self._category_note.setObjectName("InfoLabel")
        lay_nav.addWidget(self._category_note)

        self._combo_plot = QComboBox()
        self._combo_plot.currentIndexChanged.connect(self._on_plot_changed)
        lay_nav.addWidget(QLabel("Plot:"))
        lay_nav.addWidget(self._combo_plot)

        self._plot_desc = QLabel("")
        self._plot_desc.setWordWrap(True)
        self._plot_desc.setObjectName("InfoLabel")
        lay_nav.addWidget(self._plot_desc)
        layout.addWidget(grp_nav)

        # ── Dynamic parameters ───────────────────────────────────────────
        self._grp_params, lay_params = make_group("Parameters")
        self._param_form = QFormLayout()
        self._param_form.setSpacing(4)
        lay_params.addLayout(self._param_form)
        self._no_params_lbl = QLabel("(no parameters)")
        self._no_params_lbl.setObjectName("InfoLabel")
        lay_params.addWidget(self._no_params_lbl)
        layout.addWidget(self._grp_params)

        # ── Draw ──────────────────────────────────────────────────────────
        self._btn_draw = icon_button(
            "▶  Draw", "results", "Render the selected diagnostic"
        )
        self._btn_draw.clicked.connect(self._on_draw)
        layout.addWidget(self._btn_draw)

    # =========================================================================
    # Right panel
    # =========================================================================

    def _build_content(self, layout: QVBoxLayout) -> None:
        self._canvas_view = CanvasResultView(
            toolbar=True,
            empty_title="No plot yet",
            empty_reason="Pick a category and plot, then click Draw.",
        )
        self._canvas = self._canvas_view.canvas
        self._canvas.set_refresh_callback(
            self._on_draw, tooltip="Render the selected diagnostic"
        )
        layout.addWidget(self._canvas_view)

    # =========================================================================
    # Data loading
    # =========================================================================

    def _on_load(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select EMTF-XML directory"
        )
        if not path:
            path, _ = QFileDialog.getOpenFileName(
                self, "Select EMTF-XML file", "",
                "EMTF-XML (*.xml);;All files (*)",
            )
        if not path:
            return
        try:
            n = self._ctrl.load(path)
            self._data_status.setText(f"{n} station(s) loaded.")
        except Exception as exc:
            self._data_status.setText(f"Load failed: {exc}")

    def _on_clear(self) -> None:
        self._ctrl.clear()
        self._data_status.setText("No airborne data loaded")

    # =========================================================================
    # Category / plot navigation
    # =========================================================================

    def _populate_category_combo(self) -> None:
        self._combo_category.blockSignals(True)
        self._combo_category.addItems(CATEGORIES)
        self._combo_category.blockSignals(False)

    def _on_category_changed(self, row: int) -> None:
        if row < 0 or row >= len(CATEGORIES):
            return
        cat = CATEGORIES[row]
        self._category_note.setText(
            f"<small style='color:#888'>{_CATEGORY_NOTE.get(cat, '')}</small>"
        )
        entries = CATALOGUE[cat]
        self._combo_plot.blockSignals(True)
        self._combo_plot.clear()
        for label, fn_name, desc, params, multi in entries:
            self._combo_plot.addItem(label)
        self._combo_plot.blockSignals(False)
        self._combo_plot.setCurrentIndex(0)
        self._on_plot_changed(0)

    def _on_plot_changed(self, idx: int) -> None:
        cat = CATEGORIES[self._combo_category.currentIndex()]
        entries = CATALOGUE[cat]
        if idx < 0 or idx >= len(entries):
            return
        label, fn_name, desc, params, multi = entries[idx]
        self._plot_desc.setText(f"<small style='color:#888'>{desc}</small>")
        self._rebuild_param_form(params)

    def _current_entry(self):
        cat = CATEGORIES[self._combo_category.currentIndex()]
        idx = self._combo_plot.currentIndex()
        entries = CATALOGUE[cat]
        if idx < 0 or idx >= len(entries):
            return cat, None
        return cat, entries[idx]

    # =========================================================================
    # Dynamic parameter form (same widget kinds as CorrectionWindow's)
    # =========================================================================

    def _rebuild_param_form(self, params: list) -> None:
        while self._param_form.rowCount():
            self._param_form.removeRow(0)
        self._param_widgets: dict[str, QWidget] = {}

        for spec in params:
            widget = self._make_widget(spec)
            self._param_widgets[spec.name] = widget
            if spec.tip:
                widget.setToolTip(spec.tip)
            self._param_form.addRow(spec.label + ":", widget)

        self._no_params_lbl.setVisible(len(params) == 0)

    def _make_widget(self, spec: ParamSpec) -> QWidget:
        if spec.kind == "spin":
            w = QSpinBox()
            lo, hi, step = spec.opts
            w.setRange(lo, hi)
            w.setSingleStep(step)
            w.setValue(int(spec.default))
            return w
        if spec.kind == "dspin":
            w = QDoubleSpinBox()
            lo, hi, step = spec.opts
            decimals = max(
                0, -int(np.floor(np.log10(step))) if step < 1 else 1
            )
            w.setRange(lo, hi)
            w.setSingleStep(step)
            w.setDecimals(decimals)
            w.setValue(float(spec.default))
            return w
        if spec.kind == "combo":
            w = QComboBox()
            w.addItems(spec.opts)
            idx = (
                spec.opts.index(spec.default)
                if spec.default in spec.opts
                else 0
            )
            w.setCurrentIndex(idx)
            return w
        if spec.kind == "check":
            w = QCheckBox()
            w.setChecked(bool(spec.default))
            return w
        w = QLineEdit(str(spec.default))
        return w

    def _get_param_values(self) -> dict:
        vals: dict = {}
        for name, widget in self._param_widgets.items():
            if isinstance(widget, QSpinBox):
                vals[name] = widget.value()
            elif isinstance(widget, QDoubleSpinBox):
                vals[name] = widget.value()
            elif isinstance(widget, QComboBox):
                vals[name] = widget.currentText()
            elif isinstance(widget, QCheckBox):
                vals[name] = widget.isChecked()
            elif isinstance(widget, QLineEdit):
                vals[name] = widget.text()
        return vals

    # =========================================================================
    # Draw
    # =========================================================================

    def _on_draw(self) -> None:
        cat, entry = self._current_entry()
        if entry is None:
            return
        label, fn_name, desc, params, multi = entry
        kwargs = self._get_param_values()
        try:
            fig = self._ctrl.generate(cat, fn_name, **kwargs)
            self._canvas.show_figure(fig)
            self._canvas_view.show_canvas()
        except Exception as exc:
            self._canvas_view.show_unavailable(
                "Plot unavailable", f"{label} failed: {exc}"
            )

    # =========================================================================
    # Theme
    # =========================================================================

    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._ctrl.dark = dark
