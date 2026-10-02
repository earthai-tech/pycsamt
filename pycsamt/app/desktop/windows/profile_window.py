# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ProfileViewerWindow — independent floating window for MT response visualisation.

Left params panel
─────────────────
  Station        ComboBox — select station (or all)
  Period range   two QDoubleSpinBox  (min / max in seconds)
  Components     checkboxes XY · YX · XX · YY
  Error bars     checkbox
  [↻ Refresh]   [⬆ Export]

Right content
─────────────
  ProfilePanel tabs:
    ρₐ / φ  |  Pseudosection ρₐ  |  Pseudosection φ  |
    Tipper (when available)  |  Phase Tensor  |  PT Strip
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.panels.profile_panel import (
    ProfilePanel,
)
from pycsamt.app.desktop.widgets.searchable_combo import (
    SearchableComboBox,
)
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    icon_button,
    make_group,
)


class ProfileViewerWindow(PanelWindow):
    """Floating Profile Viewer with parameters and scientific plots."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(
            title="Profile Viewer",
            session_key="profile_viewer",
            params_width=260,
            icon_name="profile-view",
            parent=parent,
        )
        self.resize(1280, 780)

    # ── Params panel (left) ───────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:
        grp_line, lay_line = make_group("Survey lines")
        self._line_list = QListWidget(self)
        self._line_list.setSelectionMode(
            QListWidget.SelectionMode.MultiSelection
        )
        self._line_list.setMaximumHeight(112)
        self._line_list.setToolTip(
            "Select one or more survey lines. All lines is the default."
        )
        self._line_list.itemSelectionChanged.connect(
            self._on_line_selection_changed
        )
        self._line_list.itemClicked.connect(self._on_line_item_clicked)
        lay_line.addWidget(self._line_list)
        layout.addWidget(grp_line)

        # ── Station ───────────────────────────────────────────────────
        grp_sta, lay_sta = make_group("Station")

        # SearchableComboBox: closed = shows selected name / placeholder;
        # open = popup with a live-filter search row + scrollable list (12 max).
        self._combo_station = SearchableComboBox(max_visible=12, parent=self)
        self._combo_station.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._combo_station.station_selected.connect(self._on_station_picked)

        lay_sta.addWidget(self._combo_station)
        layout.addWidget(grp_sta)

        # ── Period range ──────────────────────────────────────────────
        grp_per, lay_per = make_group("Period range (s)")
        form = QFormLayout()
        form.setSpacing(4)
        self._spin_tmin = QDoubleSpinBox()
        self._spin_tmin.setRange(1e-6, 1e6)
        self._spin_tmin.setDecimals(4)
        self._spin_tmin.setValue(1e-4)
        self._spin_tmin.setSingleStep(0.001)
        self._spin_tmax = QDoubleSpinBox()
        self._spin_tmax.setRange(1e-6, 1e6)
        self._spin_tmax.setDecimals(2)
        self._spin_tmax.setValue(1000.0)
        self._spin_tmax.setSingleStep(10.0)
        form.addRow("Min T:", self._spin_tmin)
        form.addRow("Max T:", self._spin_tmax)
        lay_per.addLayout(form)
        layout.addWidget(grp_per)

        # ── Components ────────────────────────────────────────────────
        grp_cmp, lay_cmp = make_group("Components")
        row1 = QHBoxLayout()
        row2 = QHBoxLayout()
        self._chk_xy = QCheckBox("XY")
        self._chk_xy.setChecked(True)
        self._chk_yx = QCheckBox("YX")
        self._chk_yx.setChecked(True)
        self._chk_xx = QCheckBox("XX")
        self._chk_xx.setChecked(False)
        self._chk_yy = QCheckBox("YY")
        self._chk_yy.setChecked(False)
        row1.addWidget(self._chk_xy)
        row1.addWidget(self._chk_yx)
        row2.addWidget(self._chk_xx)
        row2.addWidget(self._chk_yy)
        lay_cmp.addLayout(row1)
        lay_cmp.addLayout(row2)
        # Auto-refresh when any component checkbox changes
        for chk in (self._chk_xy, self._chk_yx, self._chk_xx, self._chk_yy):
            chk.toggled.connect(self._on_component_changed)
        layout.addWidget(grp_cmp)

        # ── Phase range ───────────────────────────────────────────────
        grp_phase, lay_phase = make_group("Phase range")
        self._combo_phase = QComboBox()
        # (label, ymin, ymax) — None means auto
        self._PHASE_RANGES = [
            ("Auto (data range)", None, None),
            ("−90° to 90°    (−π/2 to π/2)", -90.0, 90.0),
            ("−180° to 180°  (−π to π)", -180.0, 180.0),
            ("−360° to 360°  (−2π to 2π)", -360.0, 360.0),
            ("0° to 90°", 0.0, 90.0),
        ]
        for label, _, _ in self._PHASE_RANGES:
            self._combo_phase.addItem(label)
        self._combo_phase.currentIndexChanged.connect(
            self._on_phase_range_changed
        )
        lay_phase.addWidget(self._combo_phase)
        layout.addWidget(grp_phase)

        # ── Display options ───────────────────────────────────────────
        grp_opt, lay_opt = make_group("Display")
        self._chk_errbar = QCheckBox("Error bars")
        self._chk_errbar.setChecked(True)
        self._chk_errbar.toggled.connect(self._on_errbar_toggled)
        self._chk_legend = QCheckBox("Legend")
        self._chk_legend.setChecked(True)
        self._chk_bw = QCheckBox("B/W (black lines)")
        self._chk_bw.setChecked(False)
        self._chk_bw.setToolTip(
            "Paint all response curves black — useful for raw-data\n"
            "publication figures where colour would not reproduce."
        )
        self._chk_bw.toggled.connect(self._on_bw_toggled)
        lay_opt.addWidget(self._chk_errbar)
        lay_opt.addWidget(self._chk_legend)
        lay_opt.addWidget(self._chk_bw)
        skew_row = QFormLayout()
        self._combo_pt_skew = QComboBox()
        self._combo_pt_skew.addItems(["Signed skew β", "Absolute skew |β|"])
        self._combo_pt_skew.setToolTip(
            "Choose signed or absolute skew colouring — shared by the\n"
            "Phase Tensor and PT Strip tabs. Absolute (|β|) puts the\n"
            "1-D/2-D vs 3-D structure threshold (|β| = 3°) on a plain\n"
            "0..max colour scale instead of splitting it across a\n"
            "signed -10..10° range."
        )
        self._combo_pt_skew.currentIndexChanged.connect(
            self._on_pt_skew_mode_changed
        )
        skew_row.addRow("Skew colour:", self._combo_pt_skew)
        lay_opt.addLayout(skew_row)
        self._chk_pt_labels = QCheckBox("Phase-tensor in-plot labels")
        self._chk_pt_labels.setChecked(True)
        self._chk_pt_labels.setToolTip(
            "Show the |β| 1-D/2-D vs 3-D legend and the size-reference\n"
            "ellipse (Phase Tensor tab) and the station name (PT Strip)."
        )
        self._chk_pt_labels.toggled.connect(self._on_pt_labels_toggled)
        lay_opt.addWidget(self._chk_pt_labels)
        layout.addWidget(grp_opt)

        # ── Topography ────────────────────────────────────────────────
        grp_topo, lay_topo = make_group("Topography")
        self._chk_topo = QCheckBox("Include in 2D sections")
        self._chk_topo.setChecked(False)
        self._chk_topo.setToolTip(
            "Enable terrain-following topography for pseudosections\n"
            "and the 2-D inversion section.\n"
            "Updates the global PYCSAMT_TOPO setting."
        )
        self._chk_topo.toggled.connect(self._on_topo_toggled)
        lay_topo.addWidget(self._chk_topo)

        form_topo = QFormLayout()
        form_topo.setSpacing(4)
        self._spin_exag = QDoubleSpinBox()
        self._spin_exag.setRange(0.1, 20.0)
        self._spin_exag.setSingleStep(0.5)
        self._spin_exag.setValue(1.0)
        self._spin_exag.setDecimals(1)
        self._spin_exag.setSuffix(" ×")
        self._spin_exag.setEnabled(False)
        self._spin_exag.setToolTip(
            "Vertical exaggeration of terrain surface (1 = true scale)"
        )
        self._spin_exag.valueChanged.connect(self._on_exag_changed)
        form_topo.addRow("Exaggeration:", self._spin_exag)
        lay_topo.addLayout(form_topo)
        layout.addWidget(grp_topo)

        # ── Action buttons ────────────────────────────────────────────
        grp_act, lay_act = make_group("Actions")
        self._btn_refresh = icon_button(
            "↻  Refresh", "profile-view", "Redraw current tab"
        )
        self._btn_refresh.clicked.connect(self._on_refresh)
        self._btn_export = icon_button(
            "⬆  Export…", "export", "Export figure to file"
        )
        self._btn_export.clicked.connect(self._on_export)
        self._btn_pub = icon_button(
            "📐  Publication…",
            "profile-view",
            "Open a standalone publication-quality multi-panel view\n"
            "for the selected station (does not overwrite the main panel).",
        )
        self._btn_pub.clicked.connect(self._on_pub_view)
        lay_act.addWidget(self._btn_refresh)
        lay_act.addWidget(self._btn_export)
        lay_act.addWidget(self._btn_pub)
        layout.addWidget(grp_act)

        # Info label at bottom
        self._info_lbl = QLabel("")
        self._info_lbl.setWordWrap(True)
        self._info_lbl.setObjectName("InfoLabel")
        layout.addWidget(self._info_lbl)

    # ── Content panel (right) ─────────────────────────────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        self._profile_panel = ProfilePanel(self)
        layout.addWidget(self._profile_panel)

        # Floating hover-reveal "hard refresh" button on every tab's own
        # canvas, alongside the sidebar Refresh button — see
        # MplCanvas.set_refresh_callback() for why it never overlaps the
        # canvas's own "open in separate window" toolbar icon.
        for canvas in (
            self._profile_panel._canvas_rho_phi,
            self._profile_panel._canvas_rho_ps,
            self._profile_panel._canvas_ph_ps,
            self._profile_panel._canvas_tipper,
            self._profile_panel._canvas_pt,
            self._profile_panel._canvas_pt_strip,
        ):
            canvas.set_refresh_callback(
                self._on_refresh,
                tooltip="Hard refresh — recompute the current tab",
            )

    # ── Public API ────────────────────────────────────────────────────

    def set_sites(self, sites, dataframe=None) -> None:
        super().set_sites(sites)
        self._all_sites = sites
        try:
            self._station_to_line = self._build_line_map(sites, dataframe)
            self._populate_line_list()
        except Exception:
            self._station_to_line = {}
            self._line_list.clear()
            self._line_list.setEnabled(False)
        try:
            self._profile_panel.set_sites(sites)
        except Exception:
            pass  # panel redraw errors must not block combo population
        self._populate_station_combo(sites)
        self._update_period_range(sites)

    def _build_line_map(self, sites, dataframe=None) -> dict[str, str]:
        """Map station names to survey lines from metadata or source paths."""
        mapping: dict[str, str] = {}
        if dataframe is not None and {"ID", "Line"}.issubset(dataframe.columns):
            for row in dataframe[["ID", "Line"]].itertuples(index=False):
                line = str(row.Line)
                if line and line != "—":
                    mapping[str(row.ID)] = line
        for site in sites:
            name = str(getattr(site, "name", ""))
            if name in mapping:
                continue
            try:
                path = getattr(site.edi, "path", None)
                if path is not None:
                    mapping[name] = str(path.parent.name)
            except Exception:
                pass
        return mapping

    def _populate_line_list(self) -> None:
        counts: dict[str, int] = {}
        for line in self._station_to_line.values():
            counts[line] = counts.get(line, 0) + 1

        self._line_list.blockSignals(True)
        self._line_list.clear()
        all_item = QListWidgetItem(
            f"All lines ({len(self._all_sites)} stations)"
        )
        all_item.setData(Qt.ItemDataRole.UserRole, None)
        self._line_list.addItem(all_item)
        all_item.setSelected(True)
        for line in sorted(counts):
            item = QListWidgetItem(f"{line} ({counts[line]} stations)")
            item.setData(Qt.ItemDataRole.UserRole, line)
            self._line_list.addItem(item)
        self._line_list.setEnabled(bool(counts))
        self._line_list.blockSignals(False)

    def _on_line_selection_changed(self) -> None:
        """Apply selected survey lines only to the Profile Viewer."""
        if not hasattr(self, "_all_sites") or self._line_list.count() == 0:
            return
        selected_items = self._line_list.selectedItems()
        specific = [
            item.data(Qt.ItemDataRole.UserRole)
            for item in selected_items
            if item.data(Qt.ItemDataRole.UserRole) is not None
        ]
        all_item = self._line_list.item(0)
        self._line_list.blockSignals(True)
        if specific:
            all_item.setSelected(False)
        elif all_item not in selected_items:
            all_item.setSelected(True)
        self._line_list.blockSignals(False)

        if specific:
            from pycsamt.site.base import Sites

            wanted = set(specific)
            sites = Sites(
                [
                    site
                    for site in self._all_sites
                    if self._station_to_line.get(site.name) in wanted
                ]
            )
        else:
            sites = self._all_sites

        self._profile_panel.set_sites(sites)
        names = [site.name for site in sites]
        current = self._profile_panel._ctrl._station_id
        if current not in names:
            self._profile_panel._ctrl.set_station(names[0] if names else None)
        self._populate_station_combo(sites)
        self._update_period_range(sites)
        label = ", ".join(specific) if specific else "All lines"
        self._info_lbl.setText(f"Lines: {label} · {len(sites)} stations")

    def _on_line_item_clicked(self, item: QListWidgetItem) -> None:
        """Give the synthetic All-lines row radio-button-like behavior."""
        is_all = item.data(Qt.ItemDataRole.UserRole) is None
        self._line_list.blockSignals(True)
        if is_all:
            self._line_list.clearSelection()
            item.setSelected(True)
        elif self._line_list.count():
            self._line_list.item(0).setSelected(False)
        self._line_list.blockSignals(False)
        self._on_line_selection_changed()

    def set_station(self, station_id: str) -> None:
        """Select a station by ID — called from main window on double-click."""
        self._apply_station(station_id)

    # ── Station selection (internal) ──────────────────────────────────

    def _on_station_picked(self, name: str) -> None:
        """Fired by SearchableComboBox when the user confirms a choice."""
        self._apply_station(name)

    def _apply_station(self, name: str) -> None:
        """Push *name* to the panel and sync the combo without re-emitting."""
        if not name:
            return
        self._combo_station.select_station(name)
        self._profile_panel.set_selected_station(name)
        self._info_lbl.setText(f"Station: {name}")

    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._profile_panel.set_dark_mode(dark)

    # ── Topo slots ────────────────────────────────────────────────────

    def _on_topo_toggled(self, checked: bool) -> None:
        from pycsamt.topo import configure_topo

        configure_topo(enabled=checked)
        self._profile_panel._mark_all_dirty()
        self._spin_exag.setEnabled(checked)
        # Redraw pseudosection tabs
        tab = self._profile_panel._tabs.currentIndex()
        if tab in (1, 2):
            self._profile_panel._redraw_current_tab()

    def _on_exag_changed(self, value: float) -> None:
        from pycsamt.topo import configure_topo

        configure_topo(exaggeration=value)
        self._profile_panel._mark_all_dirty()
        if self._chk_topo.isChecked():
            tab = self._profile_panel._tabs.currentIndex()
            if tab in (1, 2):
                self._profile_panel._redraw_current_tab()

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_errbar_toggled(self, checked: bool) -> None:
        self._profile_panel._ctrl.set_show_errbar(checked)
        if self._profile_panel._tabs.currentIndex() == 0:
            self._profile_panel._redraw_rho_phi()

    def _on_bw_toggled(self, checked: bool) -> None:
        self._profile_panel._ctrl.set_bw_mode(checked)
        if self._profile_panel._tabs.currentIndex() == 0:
            self._profile_panel._redraw_rho_phi()

    def _on_pt_labels_toggled(self, checked: bool) -> None:
        panel = self._profile_panel
        panel._ctrl.set_pt_annotations(checked)
        panel._dirty_canvases.add(panel._canvas_pt)
        panel._dirty_canvases.add(panel._canvas_pt_strip)
        current = panel._tabs.currentWidget()
        if current in (panel._canvas_pt, panel._canvas_pt_strip):
            panel._redraw_current_tab(force=True)

    def _on_pt_skew_mode_changed(self, index: int) -> None:
        """Switch Phase Tensor / PT Strip colouring between signed and
        absolute skew — one shared setting drives both tabs."""
        panel = self._profile_panel
        panel._ctrl.set_pt_absolute_skew(index == 1)
        panel._dirty_canvases.add(panel._canvas_pt)
        panel._dirty_canvases.add(panel._canvas_pt_strip)
        current = panel._tabs.currentWidget()
        if current in (panel._canvas_pt, panel._canvas_pt_strip):
            panel._redraw_current_tab(force=True)

    def _on_component_changed(self) -> None:
        """Push new component selection to controller and redraw ρₐ/φ tab."""
        self._apply_components()
        self._profile_panel._mark_all_dirty()
        # Only re-render if ρₐ/φ tab (index 0) is active or
        # the user is on that tab — always safe to refresh it
        tab = self._profile_panel._tabs.currentIndex()
        if tab == 0:
            self._profile_panel._redraw_rho_phi()
        # Also refresh if on any pseudosection tab (components affect all)
        elif tab in (1, 2):
            self._profile_panel._redraw_current_tab()

    def _on_phase_range_changed(self, idx: int) -> None:
        """Apply the selected phase y-limit to the controller and redraw."""
        _label, ymin, ymax = self._PHASE_RANGES[idx]
        self._profile_panel._ctrl.set_phase_ylim(ymin, ymax)
        # Refresh ρₐ/φ tab if active; otherwise mark stale for next switch
        if self._profile_panel._tabs.currentIndex() == 0:
            self._profile_panel._redraw_rho_phi()

    def _apply_components(self) -> None:
        """Collect checked components and push to PlotController."""
        comps = []
        if self._chk_xy.isChecked():
            comps.append("xy")
        if self._chk_yx.isChecked():
            comps.append("yx")
        if self._chk_xx.isChecked():
            comps.append("xx")
        if self._chk_yy.isChecked():
            comps.append("yy")
        if not comps:
            comps = ["xy"]  # always show at least one component
            self._chk_xy.setChecked(True)
        self._profile_panel._ctrl.set_components(tuple(comps))

    def _on_refresh(self) -> None:
        T_min = self._spin_tmin.value()
        T_max = self._spin_tmax.value()
        if T_min >= T_max:
            T_min, T_max = None, None
        # Push all current param state to controller before redrawing
        self._apply_components()
        idx = self._combo_phase.currentIndex()
        _label, ymin, ymax = self._PHASE_RANGES[idx]
        self._profile_panel._ctrl.set_phase_ylim(ymin, ymax)
        self._profile_panel._ctrl.set_period_range(T_min, T_max)
        # Explicit refresh → force the Phase Tensor tab to fully recompute,
        # even if the panel-level key would otherwise say "nothing changed".
        self._profile_panel.invalidate_phase_tensor()
        self._profile_panel._mark_all_dirty()
        self._profile_panel._redraw_current_tab()

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )

        canvas = self._profile_panel.current_canvas()
        fig = canvas.figure if canvas is not None else None
        if fig:
            ExportDialog(
                figure=fig,
                figure_factory=self._build_publication_export_figure,
                parent=self,
            ).exec()

    def _build_publication_export_figure(self):
        """Redraw the active plot as a white publication-quality figure."""
        from matplotlib.figure import Figure

        panel = self._profile_panel
        widget = panel._tabs.currentWidget()
        try:
            n_stations = len(panel._ctrl._sites or ())
        except Exception:
            n_stations = 1

        multi_station = widget in (
            panel._canvas_rho_ps,
            panel._canvas_ph_ps,
            panel._canvas_pt,
        )
        if multi_station:
            width = max(10.0, min(24.0, 4.0 + 0.16 * n_stations))
            figsize = (width, 7.5)
        elif widget is panel._canvas_pt_strip:
            figsize = (11.0, 4.2)
        else:
            figsize = (9.0, 7.0)

        fig = Figure(figsize=figsize, facecolor="white")
        old_dark = panel._ctrl.dark
        panel._ctrl.dark = False
        try:
            if widget is panel._canvas_rho_phi:
                panel._ctrl.draw_rho_phi(fig)
            else:
                ax = fig.add_subplot(111)
                draw = {
                    panel._canvas_rho_ps: panel._ctrl.draw_rho_pseudosection,
                    panel._canvas_ph_ps: panel._ctrl.draw_phase_pseudosection,
                    panel._canvas_tipper: panel._ctrl.draw_tipper,
                    panel._canvas_pt: panel._ctrl.draw_phase_tensor,
                    panel._canvas_pt_strip: panel._ctrl.draw_phase_tensor_strip,
                }.get(widget)
                if draw is not None:
                    draw(ax)
        finally:
            panel._ctrl.dark = old_dark

        fig.set_facecolor("white")
        for ax in fig.axes:
            ax.set_facecolor("white")
            ax.tick_params(colors="#222222", labelcolor="#222222")
            ax.xaxis.label.set_color("#222222")
            ax.yaxis.label.set_color("#222222")
            ax.title.set_color("#111111")
            for spine in ax.spines.values():
                spine.set_color("#444444")
            for text in ax.texts:
                text.set_color("#222222")
        if multi_station:
            fig.subplots_adjust(left=0.07, right=0.94, bottom=0.10, top=0.82)
        return fig

    def _on_pub_view(self) -> None:
        """Open a standalone publication-quality figure in a new dialog."""
        from pycsamt.app.desktop.windows.publication_view_dialog import (
            PublicationViewDialog,
        )

        ctrl = self._profile_panel._ctrl
        name = self._combo_station.current_station()
        if not name:
            self._info_lbl.setText(
                "Select a station before opening Publication View."
            )
            return
        # Push current component / errbar state
        self._apply_components()
        dlg = PublicationViewDialog(
            controller=ctrl,
            station_name=name,
            dark=ctrl.dark,
            parent=self,
        )
        dlg.show()  # non-modal: user can keep interacting with the main window

    # ── Helpers ───────────────────────────────────────────────────────

    def _populate_station_combo(self, sites) -> None:
        """Populate the SearchableComboBox with station names from *sites*."""
        names: list[str] = []
        try:
            for site in sites:
                names.append(site.name)
        except Exception:
            pass
        self._combo_station.set_names(names)

    def _update_period_range(self, sites) -> None:
        try:
            all_f = []
            for site in sites:
                f = site.freq
                if f is not None:
                    all_f.extend(f.ravel().tolist())
            if all_f:
                f_min = max(float(min(all_f)), 1e-6)
                f_max = float(max(all_f))
                self._spin_tmin.setValue(round(1.0 / f_max, 5))
                self._spin_tmax.setValue(round(1.0 / f_min, 1))
        except Exception:
            pass
