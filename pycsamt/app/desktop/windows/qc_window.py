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
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGridLayout,
    QLabel,
    QLineEdit,
    QSizePolicy,
    QSpinBox,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.qc_controller import (
    ALL_GROUPS,
    GROUP_ICONS,
    QCController,
    QC_STATIC_SHIFT_METHODS,
    QC_STATIC_SHIFT_PLOTS,
    describe_plot,
    qc_parameter_specs,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.widgets.unavailable_view import (
    UnavailableResultView,
)
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
            params_width=300,
            icon_name="qc",
            parent=parent,
        )
        self.resize(1200, 800)
        self._ctrl = QCController()
        self._render_timer = QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.timeout.connect(self._on_run)
        self._parameter_widgets: dict[str, tuple[QWidget, str]] = {}
        self._parameter_cache: dict[str, dict] = {}
        self._method_parameter_cache: dict[tuple[str, str], dict] = {}
        self._parameter_method: str | None = None
        self._parameter_plot: str | None = None
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
        self._combo_category.currentIndexChanged.connect(
            self._on_category_changed
        )
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

        self._controls_group, self._controls_layout = make_group("Controls")
        layout.addWidget(self._controls_group)

        self._view_group, self._view_layout = make_group("Plot view")
        layout.addWidget(self._view_group)

        self._status_lbl = QLabel("")
        self._status_lbl.setObjectName("InfoLabel")
        layout.addWidget(self._status_lbl)

        # ── Actions ───────────────────────────────────────────────────
        self._actions_group, lay_act = make_group("Actions")
        self._btn_run = icon_button(
            "↻  Hard refresh",
            "qc",
            "Recompute the selected plot immediately",
        )
        self._btn_export = icon_button(
            "⬆  Export…", "export", "Save figure to file"
        )
        self._btn_run.clicked.connect(self._on_hard_refresh)
        self._btn_export.clicked.connect(self._on_export)
        lay_act.addWidget(self._btn_run)
        lay_act.addWidget(self._btn_export)
        layout.addWidget(self._actions_group)

    # ── Content panel ─────────────────────────────────────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        self._result_stack = QStackedWidget()
        self._plot_page = QWidget(self)
        plot_layout = QGridLayout(self._plot_page)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.setSpacing(0)

        self._canvas = MplCanvas(self._plot_page, toolbar=True)
        plot_layout.addWidget(self._canvas, 0, 0)

        # The overlay floats over the canvas's own drawing area (below its
        # navigation toolbar strip), so it never overlaps toolbar icons
        # like "Open in separate plot window" — see MplCanvas.
        self._canvas.set_refresh_callback(
            self._on_hard_refresh,
            tooltip="Hard refresh — recompute this plot immediately",
        )
        self._btn_overlay_refresh = self._canvas.refresh_button

        self._unavailable_view = UnavailableResultView(self)
        self._result_stack.addWidget(self._plot_page)
        self._result_stack.addWidget(self._unavailable_view)
        self._result_stack.setCurrentWidget(self._unavailable_view)
        self._btn_export.setEnabled(False)
        layout.addWidget(self._result_stack)

    # ── Public API ────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        super().set_sites(sites)
        self._ctrl.set_sites(sites)
        try:
            category = self._combo_category.currentIndex()
            plot = self._combo_plot.currentIndex()
            fn_name = ALL_GROUPS[category][1][plot][1]
            self._rebuild_parameter_editors(fn_name)
        except (IndexError, TypeError):
            pass
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
        self._on_plot_changed(0)

    def _on_plot_changed(self, row: int) -> None:
        category = self._combo_category.currentIndex()
        self._update_desc(category, row)
        try:
            _label, plots = ALL_GROUPS[category]
            _plot_label, fn_name, _has_ax = plots[row]
        except (IndexError, TypeError):
            return
        self._rebuild_parameter_editors(fn_name)
        self._request_render()

    def _clear_parameter_layout(self, layout: QVBoxLayout) -> None:
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.hide()
                widget.deleteLater()

    def _rebuild_parameter_editors(
        self, fn_name: str, method: str | None = None,
    ) -> None:
        """Build scientific and view controls for the selected function."""
        self._store_parameter_values()
        self._clear_parameter_layout(self._controls_layout)
        self._clear_parameter_layout(self._view_layout)
        self._parameter_widgets = {}
        self._parameter_plot = fn_name
        saved = self._parameter_cache.get(fn_name, {})
        method = method or saved.get("method", "ama")
        specs = qc_parameter_specs(fn_name, method)
        if fn_name in QC_STATIC_SHIFT_PLOTS:
            method_saved = self._method_parameter_cache.get(
                (fn_name, method), {}
            )
            saved = {
                spec.name: (method_saved if spec.method_specific else saved)
                .get(spec.name, spec.default)
                for spec in specs
            }
            saved["method"] = method
            self._parameter_method = method
        else:
            self._parameter_method = None
        control_specs = [spec for spec in specs if not spec.view]
        view_specs = [spec for spec in specs if spec.view]

        self._add_parameter_form(
            self._controls_layout,
            control_specs,
            saved,
            empty_text="This plot uses its standard analytical defaults.",
        )
        self._add_parameter_form(self._view_layout, view_specs, saved)
        self._view_group.setVisible(bool(view_specs))

    def _add_parameter_form(
        self,
        host: QVBoxLayout,
        specs,
        saved: dict,
        empty_text: str = "",
    ) -> None:
        if not specs:
            if empty_text:
                label = QLabel(empty_text)
                label.setObjectName("InfoLabel")
                label.setWordWrap(True)
                host.addWidget(label)
            return
        container = QWidget()
        form = QFormLayout(container)
        form.setContentsMargins(0, 0, 0, 0)
        form.setSpacing(6)
        form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        for spec in specs:
            widget = self._make_parameter_widget(spec, saved.get(spec.name))
            self._parameter_widgets[spec.name] = (widget, spec.kind)
            form.addRow(f"{spec.label}:", widget)
        host.addWidget(container)

    def _make_parameter_widget(self, spec, saved_value):
        value = spec.default if saved_value is None else saved_value
        if spec.kind == "bool":
            widget = QCheckBox()
            widget.setChecked(bool(value))
            widget.toggled.connect(self._on_parameter_changed)
        elif spec.kind == "int":
            widget = QSpinBox()
            widget.setRange(0, 1_000_000)
            if spec.name == "half_window":
                widget.setMinimum(1)
            elif spec.name == "poly":
                widget.setRange(0, 1)
            widget.setValue(int(value))
            widget.valueChanged.connect(self._on_parameter_changed)
        elif spec.kind == "float":
            widget = QDoubleSpinBox()
            lower = -1_000_000_000.0 if spec.name == "rotate_deg" else 0.0
            widget.setRange(lower, 1_000_000_000.0)
            widget.setDecimals(5)
            widget.setSingleStep(max(abs(float(value)) / 10.0, 0.01))
            widget.setValue(float(value))
            widget.valueChanged.connect(self._on_parameter_changed)
        elif spec.choices:
            widget = QComboBox()
            widget.setMinimumContentsLength(14)
            widget.setSizeAdjustPolicy(
                QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
            )
            if spec.default is None:
                widget.addItem("Use default", None)
            for choice in spec.choices:
                if choice == "auto" and spec.default is None:
                    continue
                label = choice.replace("_", " ").title()
                if (self._parameter_plot in QC_STATIC_SHIFT_PLOTS
                        and spec.name == "method"):
                    label = QC_STATIC_SHIFT_METHODS[choice]
                widget.addItem(label, choice)
            index = widget.findData(value)
            widget.setCurrentIndex(max(index, 0))
            widget.currentIndexChanged.connect(self._on_parameter_changed)
        elif spec.name == "station":
            widget = QComboBox()
            widget.addItem("Automatic", None)
            for station in self._station_names():
                widget.addItem(station, station)
            index = widget.findData(value)
            widget.setCurrentIndex(max(index, 0))
            widget.currentIndexChanged.connect(self._on_parameter_changed)
        else:
            widget = QLineEdit()
            if spec.kind in {"sequence", "optional_sequence"}:
                if value is not None:
                    widget.setText(
                        value
                        if isinstance(value, str)
                        else ", ".join(str(part) for part in value)
                    )
                widget.setPlaceholderText("Comma-separated values")
            elif value is not None:
                widget.setText(str(value))
            if spec.name == "source_offset":
                widget.setPlaceholderText("Use metadata or enter metres")
            elif spec.name in {"sig_dist", "sig_val"}:
                widget.setPlaceholderText("Automatic (leave blank)")
            widget.textChanged.connect(self._on_parameter_changed)
        widget.setToolTip(
            f"Function parameter: {spec.name} (default: {spec.default!r})"
        )
        return widget

    def _station_names(self) -> list[str]:
        if self._ctrl._sites is None:
            return []
        try:
            from pycsamt.emtools._core import _iter_items, _name

            return [
                _name(item, index)
                for index, item in enumerate(_iter_items(self._ctrl._sites))
            ]
        except Exception:
            return []

    @staticmethod
    def _raw_parameter_value(widget: QWidget):
        if isinstance(widget, QLineEdit):
            return widget.text()
        if isinstance(widget, QComboBox):
            return widget.currentData()
        if isinstance(widget, QCheckBox):
            return widget.isChecked()
        if isinstance(widget, (QSpinBox, QDoubleSpinBox)):
            return widget.value()
        return None

    def _store_parameter_values(self) -> None:
        if not getattr(self, "_parameter_plot", None):
            return
        self._parameter_cache[self._parameter_plot] = {
            name: self._raw_parameter_value(widget)
            for name, (widget, _kind) in self._parameter_widgets.items()
        }
        if self._parameter_method is not None:
            self._method_parameter_cache[
                (self._parameter_plot, self._parameter_method)
            ] = self._parameter_cache[self._parameter_plot].copy()

    def _collect_plot_kwargs(self) -> tuple[dict, str | None]:
        kwargs: dict = {}
        try:
            for name, (widget, kind) in self._parameter_widgets.items():
                raw = self._raw_parameter_value(widget)
                if kind in {"optional_float", "optional_int"}:
                    text = str(raw).strip()
                    if not text:
                        continue
                    value = float(text)
                    if kind == "optional_int":
                        value = int(value)
                elif kind in {"sequence", "optional_sequence"}:
                    text = str(raw).strip().strip("()[]")
                    if not text:
                        if kind == "optional_sequence":
                            continue
                        return {}, f"{name.replace('_', ' ')} cannot be empty."
                    value = tuple(
                        float(part.strip())
                        for part in text.split(",")
                        if part.strip()
                    )
                elif kind == "optional_str":
                    if raw is None or not str(raw).strip():
                        continue
                    value = str(raw).strip()
                else:
                    value = raw
                if name == "source_offset" and value is not None and value <= 0:
                    return {}, "Source offset must be greater than zero metres."
                if name in {"sig_dist", "sig_val"} and value <= 0:
                    return {}, "Filter bandwidths must be greater than zero."
                kwargs[name] = value
        except ValueError:
            return {}, "Enter valid numeric values in the highlighted controls."

        if {"near_threshold", "far_threshold"} <= kwargs.keys():
            if kwargs["near_threshold"] >= kwargs["far_threshold"]:
                return {}, "Near threshold must be lower than far threshold."
        if {"ci_lo", "ci_hi"} <= kwargs.keys():
            if kwargs["ci_lo"] >= kwargs["ci_hi"]:
                return {}, "CI low must be lower than CI high."
        self._store_parameter_values()
        return kwargs, None

    def _on_parameter_changed(self, *_args) -> None:
        if self._parameter_plot in QC_STATIC_SHIFT_PLOTS:
            widget, _ = self._parameter_widgets["method"]
            method = widget.currentData()
            if method != self._parameter_method:
                self._rebuild_parameter_editors(self._parameter_plot, method)
        self._store_parameter_values()
        self._request_render(350)

    def _request_render(self, delay_ms: int = 0) -> None:
        """Debounce automatic redraws caused by selection/control changes."""
        if not hasattr(self, "_render_timer") or self._ctrl._sites is None:
            return
        if not self.isVisible():
            self._auto_rendered = False
            return
        self._render_timer.start(max(0, delay_ms))

    def _on_hard_refresh(self) -> None:
        self._render_timer.stop()
        self._on_run()

    def _on_run(self) -> None:
        cat_row = self._combo_category.currentIndex()
        plot_row = self._combo_plot.currentIndex()
        if cat_row < 0 or plot_row < 0:
            return
        _label, plots = ALL_GROUPS[cat_row]
        if plot_row >= len(plots):
            return
        _plot_label, fn_name, has_ax = plots[plot_row]
        if self._ctrl._sites is None:
            self._status_lbl.setText("Load survey data first.")
            self._show_unavailable(
                "Load survey data to begin",
                "No stations are currently available to the QC dashboard.",
                "Load EDI or EMTF-XML files, then select a diagnostic.",
            )
            return
        kwargs, validation_error = self._collect_plot_kwargs()
        if validation_error:
            self._status_lbl.setText("Waiting for valid parameters.")
            self._show_unavailable(
                "Adjust the selected plot controls",
                validation_error,
                (
                    "Correct the value in Controls or Plot view; the plot "
                    "will update automatically."
                ),
            )
            return
        self._status_lbl.setText(f"Running {fn_name}…")
        self._btn_run.setEnabled(False)
        self._btn_overlay_refresh.setEnabled(False)
        try:
            new_fig = self._ctrl.draw(
                fn_name,
                has_ax,
                self._canvas.figure,
                **kwargs,
            )
            unavailable = self._ctrl.last_unavailable
            if unavailable is not None:
                self._show_unavailable(
                    unavailable.title,
                    unavailable.reason,
                    unavailable.guidance,
                )
                self._status_lbl.setText("Result unavailable.")
                return
            if new_fig is not None:
                self._canvas.show_figure(new_fig)
            else:
                self._canvas.fit_to_view()
            self._result_stack.setCurrentWidget(self._plot_page)
            self._btn_export.setEnabled(True)
            self._status_lbl.setText("Done.")
        except Exception as exc:
            self._status_lbl.setText(f"Error: {exc}")
            self._show_unavailable(
                "The diagnostic could not be displayed",
                str(exc) or "An unexpected rendering error occurred.",
                "Review the selected data and try the diagnostic again.",
            )
        finally:
            self._btn_run.setEnabled(True)
            self._btn_overlay_refresh.setEnabled(True)

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )

        ExportDialog(figure=self._canvas.figure, parent=self).exec()

    # ── Helpers ───────────────────────────────────────────────────────

    def _show_unavailable(
        self, title: str, reason: str, guidance: str = ""
    ) -> None:
        """Show a useful explanation instead of a blank plot and toolbar."""
        self._unavailable_view.set_content(title, reason, guidance)
        self._result_stack.setCurrentWidget(self._unavailable_view)
        self._btn_export.setEnabled(False)

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
