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
  7-tab ProfilePanel:
    ρₐ / φ  |  Pseudosection ρₐ  |  Pseudosection φ  |
    Tipper  |  Phase Tensor       |  PT Strip  |  2D Section
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
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
    """Floating Profile Viewer — params left, 6-tab plot panel right."""

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
        self._chk_xy = QCheckBox("XY"); self._chk_xy.setChecked(True)
        self._chk_yx = QCheckBox("YX"); self._chk_yx.setChecked(True)
        self._chk_xx = QCheckBox("XX"); self._chk_xx.setChecked(False)
        self._chk_yy = QCheckBox("YY"); self._chk_yy.setChecked(False)
        row1.addWidget(self._chk_xy); row1.addWidget(self._chk_yx)
        row2.addWidget(self._chk_xx); row2.addWidget(self._chk_yy)
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
            ("Auto (data range)",          None,   None),
            ("−90° to 90°    (−π/2 to π/2)",  -90.0,  90.0),
            ("−180° to 180°  (−π to π)",   -180.0, 180.0),
            ("−360° to 360°  (−2π to 2π)", -360.0, 360.0),
            ("0° to 90°",                    0.0,   90.0),
        ]
        for label, _, _ in self._PHASE_RANGES:
            self._combo_phase.addItem(label)
        self._combo_phase.currentIndexChanged.connect(self._on_phase_range_changed)
        lay_phase.addWidget(self._combo_phase)
        layout.addWidget(grp_phase)

        # ── Display options ───────────────────────────────────────────
        grp_opt, lay_opt = make_group("Display")
        self._chk_errbar = QCheckBox("Error bars")
        self._chk_errbar.setChecked(True)
        self._chk_errbar.toggled.connect(self._on_errbar_toggled)
        self._chk_legend  = QCheckBox("Legend")
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
        self._spin_exag.setToolTip("Vertical exaggeration of terrain surface (1 = true scale)")
        self._spin_exag.valueChanged.connect(self._on_exag_changed)
        form_topo.addRow("Exaggeration:", self._spin_exag)
        lay_topo.addLayout(form_topo)
        layout.addWidget(grp_topo)

        # ── Action buttons ────────────────────────────────────────────
        grp_act, lay_act = make_group("Actions")
        self._btn_refresh = icon_button("↻  Refresh", "profile-view", "Redraw current tab")
        self._btn_refresh.clicked.connect(self._on_refresh)
        self._btn_export  = icon_button("⬆  Export…", "export", "Export figure to file")
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

    # ── Public API ────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        super().set_sites(sites)
        try:
            self._profile_panel.set_sites(sites)
        except Exception:
            pass   # panel redraw errors must not block combo population
        self._populate_station_combo(sites)
        self._update_period_range(sites)

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
        self._spin_exag.setEnabled(checked)
        # Sync the section_panel checkbox if it exists
        try:
            sp = self._profile_panel._section_panel
            sp._chk_topo.blockSignals(True)
            sp._chk_topo.setChecked(checked)
            sp._chk_topo.blockSignals(False)
        except Exception:
            pass
        # Redraw pseudosection tabs
        tab = self._profile_panel._tabs.currentIndex()
        if tab in (1, 2):
            self._profile_panel._redraw_current_tab()
        elif tab == 6:
            self._profile_panel._section_panel._redraw()

    def _on_exag_changed(self, value: float) -> None:
        from pycsamt.topo import configure_topo
        configure_topo(exaggeration=value)
        if self._chk_topo.isChecked():
            tab = self._profile_panel._tabs.currentIndex()
            if tab in (1, 2):
                self._profile_panel._redraw_current_tab()
            elif tab == 5:
                self._profile_panel._section_panel._redraw()

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_errbar_toggled(self, checked: bool) -> None:
        self._profile_panel._ctrl.set_show_errbar(checked)
        if self._profile_panel._tabs.currentIndex() == 0:
            self._profile_panel._redraw_rho_phi()

    def _on_bw_toggled(self, checked: bool) -> None:
        self._profile_panel._ctrl.set_bw_mode(checked)
        if self._profile_panel._tabs.currentIndex() == 0:
            self._profile_panel._redraw_rho_phi()

    def _on_component_changed(self) -> None:
        """Push new component selection to controller and redraw ρₐ/φ tab."""
        self._apply_components()
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
        if self._chk_xy.isChecked(): comps.append("xy")
        if self._chk_yx.isChecked(): comps.append("yx")
        if self._chk_xx.isChecked(): comps.append("xx")
        if self._chk_yy.isChecked(): comps.append("yy")
        if not comps:
            comps = ["xy"]   # always show at least one component
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
        self._profile_panel._redraw_current_tab()

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )
        tab  = self._profile_panel._tabs.currentIndex()
        canvases = [
            self._profile_panel._canvas_rho_phi,
            self._profile_panel._canvas_rho_ps,
            self._profile_panel._canvas_ph_ps,
            self._profile_panel._canvas_tipper,
            self._profile_panel._canvas_pt,
            self._profile_panel._canvas_pt_strip,
        ]
        fig = canvases[tab].figure if tab < len(canvases) else None
        if fig:
            ExportDialog(figure=fig, parent=self).exec()

    def _on_pub_view(self) -> None:
        """Open a standalone publication-quality figure in a new dialog."""
        from pycsamt.app.desktop.windows.publication_view_dialog import (
            PublicationViewDialog,
        )
        ctrl  = self._profile_panel._ctrl
        name  = self._combo_station.current_station()
        if not name:
            self._info_lbl.setText("Select a station before opening Publication View.")
            return
        # Push current component / errbar state
        self._apply_components()
        dlg = PublicationViewDialog(
            controller=ctrl,
            station_name=name,
            dark=ctrl.dark,
            parent=self,
        )
        dlg.show()     # non-modal: user can keep interacting with the main window

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
