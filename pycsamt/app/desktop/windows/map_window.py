# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
MapViewerWindow — enhanced floating map viewer.

Left param panel (scrollable)
─────────────────────────────
  Map Type         station | elevation | depth | resistivity
  ──────────── (conditional) ─────────────────────────────────
  Freq / Period    combo populated from loaded EDI data; Hz ↔ s toggle
  Component        XY | YX | Det
  ──────────── (depth only) ────────────────────────────────────
  Depth Settings   target depth (m), info label showing ≈ period
  ──────────── ──────────────────────────────────────────────────
  Input CRS        Geographic (EPSG:4326) | UTM Zone (1-60, N/S) | Custom EPSG
  Display          colormap, marker size, opacity, colorbar orientation
  Contour          none | lines | filled | filled+labels; N levels
  Overlay          profile lines, station labels, grid
  Basemap          provider combo, opacity, install hint
  Actions          Refresh | Export

Right content
─────────────
  MapPanel — interactive matplotlib map
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore    import Qt, Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QSlider,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.panels.map_panel   import MapPanel
from pycsamt.app.desktop.windows._base      import (
    PanelWindow,
    make_group,
    icon_button,
    _icon,
)


# ── available colormaps ───────────────────────────────────────────────────────

_CMAPS = [
    "plasma", "viridis", "magma", "inferno", "jet",
    "RdBu_r", "seismic", "coolwarm", "terrain",
    "hot", "YlOrRd", "copper",
]

# ── basemap providers (must match MapPanel._BASEMAP_PROVIDERS keys) ───────────

_BASEMAP_ITEMS = [
    "None",
    "OpenStreetMap",
    "Stamen Terrain",
    "Satellite (ESRI)",
    "CartoDB Positron",
    "CartoDB DarkMatter",
]

_MAP_TYPE_ICONS = {
    "Station": "station-response",
    "Elevation": "elevation",
    "Depth": "depth",
    "Resistivity": "resistivity",
}


def _add_icon_item(combo: QComboBox, label: str, icon_name: str) -> None:
    icon = _icon(icon_name)
    if icon.isNull():
        combo.addItem(label)
    else:
        combo.addItem(icon, label)


class _VSep(QFrame):
    """Thin vertical separator for use in horizontal toolbar rows."""

    def __init__(self) -> None:
        super().__init__()
        self.setFrameShape(QFrame.Shape.VLine)
        self.setFrameShadow(QFrame.Shadow.Sunken)
        self.setFixedWidth(2)


class MapViewerWindow(PanelWindow):
    """Floating Map Viewer with elevation, depth, and resistivity map types."""

    station_selected = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(
            title="Map Viewer",
            session_key="map_viewer",
            params_width=270,
            icon_name="map-view",
            parent=parent,
        )
        self.resize(1100, 740)
        self._freq_list: list[float] = []

    # ── params panel ──────────────────────────────────────────────────────────

    def _build_params(self, layout: QVBoxLayout) -> None:

        # ── 1. Map type ───────────────────────────────────────────────────────
        grp_type, lay_type = make_group("Map Type")
        self._combo_type = QComboBox()
        for item in ("Station", "Elevation", "Depth", "Resistivity"):
            _add_icon_item(
                self._combo_type,
                item,
                _MAP_TYPE_ICONS.get(item, "map-view"),
            )
        self._combo_type.currentTextChanged.connect(self._on_type_changed)
        lay_type.addWidget(self._combo_type)
        layout.addWidget(grp_type)

        # ── 2. Frequency / Period (depth + resistivity) ───────────────────────
        self._grp_freq, lay_freq = make_group("Frequency / Period")
        form_f = QFormLayout()
        form_f.setSpacing(4)

        # Hz ↔ s radio buttons
        radio_row = QHBoxLayout()
        self._radio_hz = QRadioButton("Hz")
        self._radio_hz.setChecked(True)
        self._radio_s  = QRadioButton("s (period)")
        self._radio_grp = QButtonGroup(self)
        self._radio_grp.addButton(self._radio_hz)
        self._radio_grp.addButton(self._radio_s)
        self._radio_hz.toggled.connect(self._on_freq_unit_toggled)
        radio_row.addWidget(self._radio_hz)
        radio_row.addWidget(self._radio_s)
        lay_freq.addLayout(radio_row)

        # Frequency combo (populated when sites load)
        self._combo_freq = QComboBox()
        self._combo_freq.setEditable(True)
        self._combo_freq.setMaxVisibleItems(20)
        form_f.addRow("Value:", self._combo_freq)

        # Component
        self._combo_comp = QComboBox()
        for c in ("XY", "YX", "Det"):
            self._combo_comp.addItem(c)
        form_f.addRow("Component:", self._combo_comp)
        lay_freq.addLayout(form_f)
        layout.addWidget(self._grp_freq)

        # ── 3. Depth settings (depth type only) ──────────────────────────────
        self._grp_depth, lay_depth = make_group("Depth Settings")
        form_d = QFormLayout()
        form_d.setSpacing(4)
        self._spin_depth = QDoubleSpinBox()
        self._spin_depth.setRange(10, 100_000)
        self._spin_depth.setValue(500)
        self._spin_depth.setSuffix(" m")
        self._spin_depth.setSingleStep(100)
        form_d.addRow("Target depth:", self._spin_depth)
        self._lbl_depth_period = QLabel("≈ — (load data first)")
        self._lbl_depth_period.setWordWrap(True)
        self._lbl_depth_period.setObjectName("InfoLabel")
        form_d.addRow("", self._lbl_depth_period)
        lay_depth.addLayout(form_d)
        self._spin_depth.valueChanged.connect(self._update_depth_label)
        layout.addWidget(self._grp_depth)

        # ── 4. CRS / Coordinate System ───────────────────────────────────────
        grp_crs, lay_crs = make_group("Input CRS")

        self._combo_crs_mode = QComboBox()
        for m in ("Geographic (lat/lon)", "UTM Zone", "Custom EPSG"):
            self._combo_crs_mode.addItem(m)
        self._combo_crs_mode.setToolTip(
            "Coordinate system of the loaded station data.\n"
            "Standard EDI files use Geographic (EPSG:4326).\n"
            "Some surveys deliver UTM or national-grid coordinates."
        )
        lay_crs.addWidget(self._combo_crs_mode)

        # UTM sub-row (zone + hemisphere)
        self._wgt_utm = QWidget()
        utm_row = QHBoxLayout(self._wgt_utm)
        utm_row.setContentsMargins(0, 0, 0, 0)
        utm_row.setSpacing(4)
        utm_row.addWidget(QLabel("Zone:"))
        self._spin_utm_zone = QSpinBox()
        self._spin_utm_zone.setRange(1, 60)
        self._spin_utm_zone.setValue(50)
        self._spin_utm_zone.setToolTip("UTM zone number (1–60)")
        utm_row.addWidget(self._spin_utm_zone)
        self._radio_utm_n = QRadioButton("N")
        self._radio_utm_n.setChecked(True)
        self._radio_utm_s = QRadioButton("S")
        _utm_hem_grp = QButtonGroup(self)
        _utm_hem_grp.addButton(self._radio_utm_n)
        _utm_hem_grp.addButton(self._radio_utm_s)
        self._utm_hem_grp = _utm_hem_grp   # keep alive
        utm_row.addWidget(self._radio_utm_n)
        utm_row.addWidget(self._radio_utm_s)
        lay_crs.addWidget(self._wgt_utm)

        # Custom EPSG sub-row
        self._wgt_epsg = QWidget()
        epsg_row = QHBoxLayout(self._wgt_epsg)
        epsg_row.setContentsMargins(0, 0, 0, 0)
        epsg_row.setSpacing(4)
        epsg_row.addWidget(QLabel("EPSG:"))
        self._edit_epsg = QLineEdit("4326")
        self._edit_epsg.setPlaceholderText("e.g. 32650")
        self._edit_epsg.setMaximumWidth(80)
        epsg_row.addWidget(self._edit_epsg)
        self._btn_epsg_validate = QPushButton("✓")
        self._btn_epsg_validate.setMaximumWidth(28)
        self._btn_epsg_validate.setToolTip("Validate EPSG code")
        self._btn_epsg_validate.clicked.connect(self._validate_epsg)
        epsg_row.addWidget(self._btn_epsg_validate)
        lay_crs.addWidget(self._wgt_epsg)

        # Resolved CRS info label
        self._lbl_crs_info = QLabel("")
        self._lbl_crs_info.setObjectName("InfoLabel")
        self._lbl_crs_info.setWordWrap(True)
        lay_crs.addWidget(self._lbl_crs_info)

        layout.addWidget(grp_crs)

        # Wire CRS controls
        self._combo_crs_mode.currentTextChanged.connect(self._on_crs_mode_changed)
        self._spin_utm_zone.valueChanged.connect(self._update_crs_info)
        self._radio_utm_n.toggled.connect(self._update_crs_info)
        self._edit_epsg.textChanged.connect(self._update_crs_info)
        self._on_crs_mode_changed(self._combo_crs_mode.currentText())  # init visibility

        # ── 5. Display ────────────────────────────────────────────────────────
        grp_disp, lay_disp = make_group("Display")
        form_disp = QFormLayout()
        form_disp.setSpacing(4)

        self._combo_cmap = QComboBox()
        self._combo_cmap.addItems(_CMAPS)
        form_disp.addRow("Colormap:", self._combo_cmap)

        self._spin_ms = QSpinBox()
        self._spin_ms.setRange(2, 30)
        self._spin_ms.setValue(8)
        form_disp.addRow("Marker size:", self._spin_ms)

        self._slider_alpha = QSlider(Qt.Orientation.Horizontal)
        self._slider_alpha.setRange(20, 100)
        self._slider_alpha.setValue(90)
        self._slider_alpha.setToolTip("Marker opacity (%)")
        form_disp.addRow("Opacity:", self._slider_alpha)

        lay_disp.addLayout(form_disp)

        # Colorbar sub-section
        self._chk_cbar = QCheckBox("Show colorbar")
        self._chk_cbar.setChecked(True)
        self._chk_cbar.toggled.connect(self._on_cbar_toggle)
        lay_disp.addWidget(self._chk_cbar)

        cbar_row = QHBoxLayout()
        self._radio_cbar_v = QRadioButton("Vertical")
        self._radio_cbar_v.setChecked(True)
        self._radio_cbar_h = QRadioButton("Horizontal")
        cbar_orient_grp = QButtonGroup(self)
        cbar_orient_grp.addButton(self._radio_cbar_v)
        cbar_orient_grp.addButton(self._radio_cbar_h)
        cbar_row.addWidget(QLabel("  Bar:"))
        cbar_row.addWidget(self._radio_cbar_v)
        cbar_row.addWidget(self._radio_cbar_h)
        lay_disp.addLayout(cbar_row)
        self._cbar_orient_grp = cbar_orient_grp   # keep alive

        layout.addWidget(grp_disp)

        # ── 5. Contour ────────────────────────────────────────────────────────
        grp_cont, lay_cont = make_group("Contour")
        form_c = QFormLayout()
        form_c.setSpacing(4)

        self._combo_contour = QComboBox()
        for m in ("None", "Lines", "Filled", "Filled + labels"):
            self._combo_contour.addItem(m)
        form_c.addRow("Mode:", self._combo_contour)

        self._spin_levels = QSpinBox()
        self._spin_levels.setRange(3, 30)
        self._spin_levels.setValue(10)
        form_c.addRow("Levels:", self._spin_levels)

        lay_cont.addLayout(form_c)
        layout.addWidget(grp_cont)

        # ── 6. Overlay ────────────────────────────────────────────────────────
        grp_over, lay_over = make_group("Overlay")
        self._chk_profile = QCheckBox("Profile lines")
        self._chk_profile.setChecked(True)
        self._chk_labels  = QCheckBox("Station labels")
        self._chk_labels.setChecked(True)
        self._chk_grid    = QCheckBox("Grid")
        self._chk_grid.setChecked(True)
        lay_over.addWidget(self._chk_profile)
        lay_over.addWidget(self._chk_labels)
        lay_over.addWidget(self._chk_grid)
        layout.addWidget(grp_over)

        # ── 7. Basemap ────────────────────────────────────────────────────────
        grp_base, lay_base = make_group("Basemap")
        self._combo_basemap = QComboBox()
        self._combo_basemap.addItems(_BASEMAP_ITEMS)
        lay_base.addWidget(self._combo_basemap)

        form_base = QFormLayout()
        form_base.setSpacing(4)
        self._slider_base_alpha = QSlider(Qt.Orientation.Horizontal)
        self._slider_base_alpha.setRange(10, 100)
        self._slider_base_alpha.setValue(70)
        self._slider_base_alpha.setToolTip("Basemap opacity (%)")
        form_base.addRow("Opacity:", self._slider_base_alpha)
        lay_base.addLayout(form_base)

        self._lbl_ctx_hint = QLabel(
            "contextily not installed — basemap tiles will trigger an install "
            "prompt when a provider is selected."
        )
        self._lbl_ctx_hint.setWordWrap(True)
        self._lbl_ctx_hint.setObjectName("InfoLabel")
        self._lbl_ctx_hint.setVisible(not self._contextily_available())
        lay_base.addWidget(self._lbl_ctx_hint)

        layout.addWidget(grp_base)

        # ── 8. Actions ────────────────────────────────────────────────────────
        grp_act, lay_act = make_group("Actions")
        self._btn_refresh = icon_button("↻  Refresh", "map-view")
        self._btn_refresh.clicked.connect(self._on_refresh)
        self._btn_export  = icon_button("⬆  Export…", "export")
        self._btn_export.clicked.connect(self._on_export)
        lay_act.addWidget(self._btn_refresh)
        lay_act.addWidget(self._btn_export)
        layout.addWidget(grp_act)

        # Initial visibility
        self._on_type_changed(self._combo_type.currentText())

    # ── content panel ─────────────────────────────────────────────────────────

    def _build_content(self, layout: QVBoxLayout) -> None:
        self._map_panel = MapPanel(self)
        self._map_panel.station_selected.connect(self.station_selected)
        self._map_panel.freq_list_ready.connect(self._populate_freq_combo)
        layout.addWidget(self._map_panel)

    # ── public API ────────────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        super().set_sites(sites)
        self._map_panel.set_sites(sites)

    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._map_panel.set_dark_mode(dark)

    def set_dataframe(self, df) -> None:
        self._map_panel.set_dataframe(df)

    def highlight_station(self, station_id: str) -> None:
        self._map_panel.highlight_station(station_id)

    # ── internal helpers ──────────────────────────────────────────────────────

    @staticmethod
    def _contextily_available() -> bool:
        try:
            import contextily  # noqa: F401
            return True
        except ImportError:
            return False

    def _populate_freq_combo(self, freq_list: list[float]) -> None:
        """Fill the frequency combo from EDI data (called via signal)."""
        self._freq_list = freq_list
        self._combo_freq.clear()
        use_hz = self._radio_hz.isChecked()
        for f in freq_list:
            if use_hz:
                self._combo_freq.addItem(f"{f:.5g} Hz")
            else:
                T = 1.0 / f if f > 0 else 0.0
                self._combo_freq.addItem(f"{T:.5g} s")
        self._update_depth_label()

    def _on_freq_unit_toggled(self, checked: bool) -> None:
        if self._freq_list:
            self._populate_freq_combo(self._freq_list)

    def _current_freq_hz(self) -> float:
        """Parse current combo text → Hz."""
        txt = self._combo_freq.currentText().strip()
        try:
            num = float(txt.split()[0])
        except (ValueError, IndexError):
            return 1.0
        if "s" in txt.lower() and "hz" not in txt.lower():
            return 1.0 / num if num > 0 else 1.0
        return num

    def _update_depth_label(self) -> None:
        depth = self._spin_depth.value()
        # Estimate median ρ̄ from freq list length as a rough proxy
        rho_bar = self._estimate_median_rho()
        if rho_bar > 0:
            # T ≈ depth² / (503² × ρ̄)
            T = depth**2 / (503.0**2 * rho_bar)
            f = 1.0 / T if T > 0 else 0.0
            self._lbl_depth_period.setText(
                f"≈ T = {T:.4g} s  (f = {f:.4g} Hz)\n"
                f"using ρ̄ ≈ {rho_bar:.3g} Ω·m"
            )
        else:
            self._lbl_depth_period.setText("Load data for period estimate")

    def _estimate_median_rho(self) -> float:
        """Rough geometric mean ρ̄ from the loaded sites (for depth label)."""
        sites = self._sites
        if sites is None:
            return 0.0
        vals: list[float] = []
        for site in sites:
            try:
                from pycsamt.site.base import _extract_z_arrays
                arrs = _extract_z_arrays(site.edi)
                z = arrs.get("z")
                f = arrs.get("freq")
                if z is None or f is None:
                    continue
                z = np.asarray(z)
                f = np.asarray(f, float)
                if z.ndim == 2 and z.shape[1] >= 4:
                    z = z[:, :4].reshape(-1, 2, 2)
                if z.ndim == 3 and z.shape[1] == 2:
                    zxy = z[:, 0, 1]
                    rho = 0.2 * np.abs(zxy)**2 / (f + 1e-24)
                    rho = rho[np.isfinite(rho) & (rho > 0)]
                    if rho.size:
                        vals.append(float(np.exp(np.log(rho).mean())))
            except Exception:
                pass
        if not vals:
            return 0.0
        return float(np.exp(np.log(np.asarray(vals)).mean()))

    def _on_type_changed(self, text: str) -> None:
        t = text.lower()
        self._grp_freq.setVisible(t in ("depth", "resistivity"))
        self._grp_depth.setVisible(t == "depth")

    def _on_cbar_toggle(self, checked: bool) -> None:
        self._radio_cbar_v.setEnabled(checked)
        self._radio_cbar_h.setEnabled(checked)

    # ── CRS helpers ───────────────────────────────────────────────────────────

    def _on_crs_mode_changed(self, text: str) -> None:
        is_utm    = "UTM" in text
        is_custom = "Custom" in text
        self._wgt_utm.setVisible(is_utm)
        self._wgt_epsg.setVisible(is_custom)
        self._update_crs_info()

    def _resolve_source_crs(self) -> str:
        """Return the EPSG string for the currently selected input CRS."""
        mode = self._combo_crs_mode.currentText()
        if "UTM" in mode:
            zone = self._spin_utm_zone.value()
            hemi = "N" if self._radio_utm_n.isChecked() else "S"
            epsg = 32600 + zone if hemi == "N" else 32700 + zone
            return f"EPSG:{epsg}"
        if "Custom" in mode:
            txt = self._edit_epsg.text().strip().upper()
            if txt.startswith("EPSG:"):
                return txt
            try:
                int(txt)
                return f"EPSG:{txt}"
            except ValueError:
                return "EPSG:4326"
        return "EPSG:4326"

    def _update_crs_info(self) -> None:
        """Update the resolved CRS info label and validate silently."""
        crs = self._resolve_source_crs()
        try:
            from pyproj import CRS
            c = CRS.from_user_input(crs)
            geo = "geographic" if c.is_geographic else "projected"
            self._lbl_crs_info.setText(f"→ {crs}  ({geo})")
            self._lbl_crs_info.setStyleSheet("")
        except Exception as exc:
            self._lbl_crs_info.setText(f"Invalid: {exc}")
            self._lbl_crs_info.setStyleSheet("color: #f38ba8;")  # red hint

    def _validate_epsg(self) -> None:
        """Explicit validate button: show a richer message in the info label."""
        crs = self._resolve_source_crs()
        try:
            from pyproj import CRS, Transformer
            c   = CRS.from_user_input(crs)
            geo = "geographic" if c.is_geographic else "projected"
            # Quick sanity-check: can we project to EPSG:3857?
            Transformer.from_crs(crs, "EPSG:3857", always_xy=True)
            self._lbl_crs_info.setText(
                f"✓ {crs}  ({geo})\n{c.name}"
            )
            self._lbl_crs_info.setStyleSheet("color: #a6e3a1;")  # green
        except Exception as exc:
            self._lbl_crs_info.setText(f"✗ {exc}")
            self._lbl_crs_info.setStyleSheet("color: #f38ba8;")

    # ── contour mode string ───────────────────────────────────────────────────

    def _contour_mode_str(self) -> str:
        lbl = self._combo_contour.currentText()
        return {
            "None":           "none",
            "Lines":          "lines",
            "Filled":         "filled",
            "Filled + labels": "filled_labels",
        }.get(lbl, "none")

    # ── pop-out sync ──────────────────────────────────────────────────────────

    def sync_from_panel(self, panel) -> None:
        """
        Mirror all settings from an existing MapPanel into this window's
        controls, then refresh.  Called by _MapPopOutButton._on_click().
        """
        # Map type
        _type_labels = {
            "station": "Station", "elevation": "Elevation",
            "depth": "Depth",     "resistivity": "Resistivity",
        }
        idx = self._combo_type.findText(_type_labels.get(panel._map_type, "Station"))
        if idx >= 0:
            self._combo_type.setCurrentIndex(idx)

        # CRS
        crs = panel._source_crs.upper()
        if crs in ("EPSG:4326", ""):
            self._combo_crs_mode.setCurrentText("Geographic (lat/lon)")
        else:
            code_str = crs.replace("EPSG:", "")
            try:
                code = int(code_str)
                if 32601 <= code <= 32660:
                    self._combo_crs_mode.setCurrentText("UTM Zone")
                    self._spin_utm_zone.setValue(code - 32600)
                    self._radio_utm_n.setChecked(True)
                elif 32701 <= code <= 32760:
                    self._combo_crs_mode.setCurrentText("UTM Zone")
                    self._spin_utm_zone.setValue(code - 32700)
                    self._radio_utm_s.setChecked(True)
                else:
                    self._combo_crs_mode.setCurrentText("Custom EPSG")
                    self._edit_epsg.setText(code_str)
            except ValueError:
                self._combo_crs_mode.setCurrentText("Custom EPSG")
                self._edit_epsg.setText(code_str)

        # Display
        idx = self._combo_cmap.findText(panel._cmap_name)
        if idx >= 0:
            self._combo_cmap.setCurrentIndex(idx)
        self._spin_ms.setValue(panel._marker_size)
        self._slider_alpha.setValue(int(panel._marker_alpha * 100))

        # Colorbar
        self._chk_cbar.setChecked(panel._show_cbar)
        if panel._cbar_orient == "horizontal":
            self._radio_cbar_h.setChecked(True)
        else:
            self._radio_cbar_v.setChecked(True)

        # Contour
        _contour_labels = {
            "none": "None", "lines": "Lines",
            "filled": "Filled", "filled_labels": "Filled + labels",
        }
        idx = self._combo_contour.findText(
            _contour_labels.get(panel._contour_mode, "None")
        )
        if idx >= 0:
            self._combo_contour.setCurrentIndex(idx)
        self._spin_levels.setValue(panel._contour_levels)

        # Overlay
        self._chk_profile.setChecked(panel._show_profile)
        self._chk_labels.setChecked(panel._show_labels)
        self._chk_grid.setChecked(panel._show_grid)

        # Basemap
        idx = self._combo_basemap.findText(panel._provider)
        if idx >= 0:
            self._combo_basemap.setCurrentIndex(idx)
        self._slider_base_alpha.setValue(int(panel._basemap_alpha * 100))

        # Frequency / component (set text directly; combo may be empty if no
        # sites are loaded yet — redraw() on map_panel handles it gracefully)
        comp_labels = {"xy": "XY", "yx": "YX", "det": "Det"}
        idx = self._combo_comp.findText(comp_labels.get(panel._component, "XY"))
        if idx >= 0:
            self._combo_comp.setCurrentIndex(idx)
        self._spin_depth.setValue(panel._target_depth_m)
        if panel._target_freq_hz > 0 and self._combo_freq.count() == 0:
            # No freq combo items yet — add a placeholder so the value isn't lost
            self._combo_freq.addItem(f"{panel._target_freq_hz:.5g} Hz")
            self._combo_freq.setCurrentIndex(0)

        # Trigger a full redraw with all synced settings
        self._on_refresh()

    # ── slots ─────────────────────────────────────────────────────────────────

    def _on_refresh(self) -> None:
        map_type   = self._combo_type.currentText().lower()
        freq_hz    = self._current_freq_hz()
        depth_m    = self._spin_depth.value()
        component  = self._combo_comp.currentText().lower()
        cmap       = self._combo_cmap.currentText()
        ms         = self._spin_ms.value()
        alpha      = self._slider_alpha.value() / 100.0
        c_mode     = self._contour_mode_str()
        c_levels   = self._spin_levels.value()
        show_prof  = self._chk_profile.isChecked()
        show_lbl   = self._chk_labels.isChecked()
        show_grid  = self._chk_grid.isChecked()
        provider   = self._combo_basemap.currentText()
        base_a     = self._slider_base_alpha.value() / 100.0
        show_cbar   = self._chk_cbar.isChecked()
        cbar_orient = "vertical" if self._radio_cbar_v.isChecked() else "horizontal"
        source_crs  = self._resolve_source_crs()

        try:
            self._map_panel.redraw(
                map_type       = map_type,
                target_freq_hz = freq_hz,
                target_depth_m = depth_m,
                component      = component,
                cmap_name      = cmap,
                marker_size    = ms,
                marker_alpha   = alpha,
                contour_mode   = c_mode,
                contour_levels = c_levels,
                show_profile   = show_prof,
                show_labels    = show_lbl,
                show_grid      = show_grid,
                provider       = provider,
                basemap_alpha  = base_a,
                show_cbar      = show_cbar,
                cbar_orient    = cbar_orient,
                source_crs     = source_crs,
            )
        except Exception:
            try:
                self._map_panel._draw_map()
            except Exception:
                pass

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog
        ExportDialog(figure=self._map_panel._canvas.figure, parent=self).exec()


# ── MapDetailWindow ───────────────────────────────────────────────────────────

class MapDetailWindow(QDialog):
    """
    Full-detail pop-out map viewer opened by the ⤢ button on MapPanel.

    Layout: two compact horizontal control strips at the top, then the map
    panel filling all remaining space.  This is intentionally different from
    MapViewerWindow (which has a left sidebar) so the pop-out feels like a
    genuine zoom/expand rather than a duplicate panel.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("Map Viewer — Full Detail")
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self._freq_list: list[float] = []
        self._build_ui()
        self.resize(1440, 920)

    # ── build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 4)
        root.setSpacing(4)

        self._build_row1(root)
        self._build_row2(root)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setObjectName("Separator")
        root.addWidget(sep)

        self._map_panel = MapPanel(self)
        self._map_panel.freq_list_ready.connect(self._on_freq_list)
        root.addWidget(self._map_panel, stretch=1)

        # Initialise conditional visibility
        self._on_type_changed(self._combo_type.currentText())
        self._on_crs_mode(self._combo_crs.currentText())

    def _build_row1(self, parent: QVBoxLayout) -> None:
        """Row 1 — map type / freq / depth / CRS / colormap."""
        row = QHBoxLayout()
        row.setSpacing(6)

        # Map type
        row.addWidget(QLabel("Type:"))
        self._combo_type = QComboBox()
        for t in ("Station", "Elevation", "Depth", "Resistivity"):
            _add_icon_item(
                self._combo_type,
                t,
                _MAP_TYPE_ICONS.get(t, "map-view"),
            )
        self._combo_type.currentTextChanged.connect(self._on_type_changed)
        row.addWidget(self._combo_type)

        row.addWidget(_VSep())

        # Frequency / Period (hidden for Station & Elevation)
        self._wgt_freq = QWidget()
        _fl = QHBoxLayout(self._wgt_freq)
        _fl.setContentsMargins(0, 0, 0, 0)
        _fl.setSpacing(4)
        _fl.addWidget(QLabel("Freq:"))
        self._combo_freq = QComboBox()
        self._combo_freq.setEditable(True)
        self._combo_freq.setMinimumWidth(110)
        _fl.addWidget(self._combo_freq)
        self._btn_hz = QRadioButton("Hz")
        self._btn_hz.setChecked(True)
        self._btn_s  = QRadioButton("s")
        _fg = QButtonGroup(self)
        _fg.addButton(self._btn_hz)
        _fg.addButton(self._btn_s)
        self._btn_hz.toggled.connect(self._on_freq_unit)
        self._freq_btn_grp = _fg
        _fl.addWidget(self._btn_hz)
        _fl.addWidget(self._btn_s)
        row.addWidget(self._wgt_freq)

        # Component (hidden for Station & Elevation)
        self._wgt_comp = QWidget()
        _cl = QHBoxLayout(self._wgt_comp)
        _cl.setContentsMargins(0, 0, 0, 0)
        _cl.setSpacing(4)
        _cl.addWidget(QLabel("Comp:"))
        self._combo_comp = QComboBox()
        for c in ("XY", "YX", "Det"):
            self._combo_comp.addItem(c)
        _cl.addWidget(self._combo_comp)
        row.addWidget(self._wgt_comp)

        # Depth target (visible only for Depth type)
        self._wgt_depth = QWidget()
        _dl = QHBoxLayout(self._wgt_depth)
        _dl.setContentsMargins(0, 0, 0, 0)
        _dl.setSpacing(4)
        _dl.addWidget(QLabel("Depth:"))
        self._spin_depth = QDoubleSpinBox()
        self._spin_depth.setRange(10, 100_000)
        self._spin_depth.setValue(500)
        self._spin_depth.setSuffix(" m")
        self._spin_depth.setSingleStep(100)
        self._spin_depth.setFixedWidth(100)
        _dl.addWidget(self._spin_depth)
        row.addWidget(self._wgt_depth)

        row.addWidget(_VSep())

        # CRS
        row.addWidget(QLabel("CRS:"))
        self._combo_crs = QComboBox()
        for m in ("Geographic (lat/lon)", "UTM Zone", "Custom EPSG"):
            self._combo_crs.addItem(m)
        self._combo_crs.currentTextChanged.connect(self._on_crs_mode)
        row.addWidget(self._combo_crs)

        # UTM zone sub-widget (conditional)
        self._wgt_utm = QWidget()
        _ul = QHBoxLayout(self._wgt_utm)
        _ul.setContentsMargins(0, 0, 0, 0)
        _ul.setSpacing(2)
        _ul.addWidget(QLabel("Zone:"))
        self._spin_utm = QSpinBox()
        self._spin_utm.setRange(1, 60)
        self._spin_utm.setValue(50)
        self._spin_utm.setFixedWidth(48)
        _ul.addWidget(self._spin_utm)
        self._radio_utm_n = QRadioButton("N")
        self._radio_utm_n.setChecked(True)
        self._radio_utm_s = QRadioButton("S")
        _uh = QButtonGroup(self)
        _uh.addButton(self._radio_utm_n)
        _uh.addButton(self._radio_utm_s)
        self._utm_hem_grp = _uh
        _ul.addWidget(self._radio_utm_n)
        _ul.addWidget(self._radio_utm_s)
        row.addWidget(self._wgt_utm)

        # Custom EPSG sub-widget (conditional)
        self._wgt_epsg = QWidget()
        _el = QHBoxLayout(self._wgt_epsg)
        _el.setContentsMargins(0, 0, 0, 0)
        _el.setSpacing(2)
        _el.addWidget(QLabel("EPSG:"))
        self._edit_epsg = QLineEdit("4326")
        self._edit_epsg.setFixedWidth(60)
        _el.addWidget(self._edit_epsg)
        row.addWidget(self._wgt_epsg)

        row.addStretch()

        # Colormap
        row.addWidget(QLabel("Cmap:"))
        self._combo_cmap = QComboBox()
        self._combo_cmap.addItems(_CMAPS)
        self._combo_cmap.setFixedWidth(100)
        row.addWidget(self._combo_cmap)

        row.addWidget(_VSep())

        # Log-scale toggle (for resistivity / depth maps)
        self._chk_log = QCheckBox("Log ρ")
        self._chk_log.setChecked(True)
        row.addWidget(self._chk_log)

        parent.addLayout(row)

    def _build_row2(self, parent: QVBoxLayout) -> None:
        """Row 2 — marker display, contour, basemap, overlays, colorbar, actions."""
        row = QHBoxLayout()
        row.setSpacing(6)

        # Marker size & opacity
        row.addWidget(QLabel("Size:"))
        self._spin_ms = QSpinBox()
        self._spin_ms.setRange(2, 30)
        self._spin_ms.setValue(8)
        self._spin_ms.setFixedWidth(50)
        row.addWidget(self._spin_ms)
        row.addWidget(QLabel("Opacity:"))
        self._slider_alpha = QSlider(Qt.Orientation.Horizontal)
        self._slider_alpha.setRange(20, 100)
        self._slider_alpha.setValue(90)
        self._slider_alpha.setFixedWidth(80)
        self._slider_alpha.setToolTip("Marker opacity (%)")
        row.addWidget(self._slider_alpha)

        row.addWidget(_VSep())

        # Contour
        row.addWidget(QLabel("Contour:"))
        self._combo_contour = QComboBox()
        for m in ("None", "Lines", "Filled", "Filled + labels"):
            self._combo_contour.addItem(m)
        row.addWidget(self._combo_contour)
        row.addWidget(QLabel("Lvl:"))
        self._spin_levels = QSpinBox()
        self._spin_levels.setRange(3, 30)
        self._spin_levels.setValue(10)
        self._spin_levels.setFixedWidth(45)
        row.addWidget(self._spin_levels)

        row.addWidget(_VSep())

        # Basemap
        row.addWidget(QLabel("Basemap:"))
        self._combo_basemap = QComboBox()
        self._combo_basemap.addItems(_BASEMAP_ITEMS)
        row.addWidget(self._combo_basemap)
        row.addWidget(QLabel("α:"))
        self._slider_base_alpha = QSlider(Qt.Orientation.Horizontal)
        self._slider_base_alpha.setRange(10, 100)
        self._slider_base_alpha.setValue(70)
        self._slider_base_alpha.setFixedWidth(70)
        self._slider_base_alpha.setToolTip("Basemap opacity (%)")
        row.addWidget(self._slider_base_alpha)

        row.addWidget(_VSep())

        # Overlay toggles
        self._chk_profile = QCheckBox("Profile")
        self._chk_profile.setChecked(True)
        self._chk_labels  = QCheckBox("Labels")
        self._chk_labels.setChecked(True)
        self._chk_grid    = QCheckBox("Grid")
        self._chk_grid.setChecked(True)
        row.addWidget(self._chk_profile)
        row.addWidget(self._chk_labels)
        row.addWidget(self._chk_grid)

        row.addWidget(_VSep())

        # Colorbar
        self._chk_cbar = QCheckBox("Colorbar")
        self._chk_cbar.setChecked(True)
        row.addWidget(self._chk_cbar)
        self._radio_cbar_v = QRadioButton("V")
        self._radio_cbar_v.setChecked(True)
        self._radio_cbar_h = QRadioButton("H")
        _cg = QButtonGroup(self)
        _cg.addButton(self._radio_cbar_v)
        _cg.addButton(self._radio_cbar_h)
        self._cbar_grp = _cg
        row.addWidget(self._radio_cbar_v)
        row.addWidget(self._radio_cbar_h)

        row.addStretch()

        # Action buttons
        self._btn_refresh = QPushButton("↻  Refresh")
        self._btn_refresh.clicked.connect(self._on_refresh)
        row.addWidget(self._btn_refresh)
        btn_export = QPushButton("⬆  Export…")
        btn_export.clicked.connect(self._on_export)
        row.addWidget(btn_export)

        parent.addLayout(row)

    # ── public API ────────────────────────────────────────────────────────────

    def set_dataframe(self, df) -> None:
        self._map_panel.set_dataframe(df)

    def set_sites(self, sites) -> None:
        self._map_panel.set_sites(sites)

    def set_dark_mode(self, dark: bool) -> None:
        self._map_panel.set_dark_mode(dark)

    def sync_from_panel(self, panel) -> None:
        """
        Mirror all settings from an existing MapPanel into this window's
        controls, then refresh.  Called by _MapPopOutButton._on_click().
        """
        # Map type
        _type_labels = {
            "station": "Station", "elevation": "Elevation",
            "depth": "Depth",     "resistivity": "Resistivity",
        }
        idx = self._combo_type.findText(_type_labels.get(panel._map_type, "Station"))
        if idx >= 0:
            self._combo_type.setCurrentIndex(idx)

        # CRS
        crs = panel._source_crs.upper()
        if crs in ("EPSG:4326", ""):
            self._combo_crs.setCurrentText("Geographic (lat/lon)")
        else:
            code_str = crs.replace("EPSG:", "")
            try:
                code = int(code_str)
                if 32601 <= code <= 32660:
                    self._combo_crs.setCurrentText("UTM Zone")
                    self._spin_utm.setValue(code - 32600)
                    self._radio_utm_n.setChecked(True)
                elif 32701 <= code <= 32760:
                    self._combo_crs.setCurrentText("UTM Zone")
                    self._spin_utm.setValue(code - 32700)
                    self._radio_utm_s.setChecked(True)
                else:
                    self._combo_crs.setCurrentText("Custom EPSG")
                    self._edit_epsg.setText(code_str)
            except ValueError:
                self._combo_crs.setCurrentText("Custom EPSG")
                self._edit_epsg.setText(crs.replace("EPSG:", ""))

        # Display
        idx = self._combo_cmap.findText(panel._cmap_name)
        if idx >= 0:
            self._combo_cmap.setCurrentIndex(idx)
        self._spin_ms.setValue(panel._marker_size)
        self._slider_alpha.setValue(int(panel._marker_alpha * 100))
        self._chk_log.setChecked(panel._log_scale)

        # Colorbar
        self._chk_cbar.setChecked(panel._show_cbar)
        if panel._cbar_orient == "horizontal":
            self._radio_cbar_h.setChecked(True)
        else:
            self._radio_cbar_v.setChecked(True)

        # Contour
        _contour_labels = {
            "none": "None", "lines": "Lines",
            "filled": "Filled", "filled_labels": "Filled + labels",
        }
        idx = self._combo_contour.findText(
            _contour_labels.get(panel._contour_mode, "None")
        )
        if idx >= 0:
            self._combo_contour.setCurrentIndex(idx)
        self._spin_levels.setValue(panel._contour_levels)

        # Overlay
        self._chk_profile.setChecked(panel._show_profile)
        self._chk_labels.setChecked(panel._show_labels)
        self._chk_grid.setChecked(panel._show_grid)

        # Basemap
        idx = self._combo_basemap.findText(panel._provider)
        if idx >= 0:
            self._combo_basemap.setCurrentIndex(idx)
        self._slider_base_alpha.setValue(int(panel._basemap_alpha * 100))

        # Component / depth / freq
        _comp_labels = {"xy": "XY", "yx": "YX", "det": "Det"}
        idx = self._combo_comp.findText(_comp_labels.get(panel._component, "XY"))
        if idx >= 0:
            self._combo_comp.setCurrentIndex(idx)
        self._spin_depth.setValue(panel._target_depth_m)
        if panel._target_freq_hz > 0 and self._combo_freq.count() == 0:
            self._combo_freq.addItem(f"{panel._target_freq_hz:.5g} Hz")
            self._combo_freq.setCurrentIndex(0)

        self._on_refresh()

    # ── internal helpers ──────────────────────────────────────────────────────

    def _on_freq_list(self, freq_list: list[float]) -> None:
        self._freq_list = freq_list
        self._populate_freq_combo()

    def _populate_freq_combo(self) -> None:
        self._combo_freq.clear()
        use_hz = self._btn_hz.isChecked()
        for f in self._freq_list:
            if use_hz:
                self._combo_freq.addItem(f"{f:.5g} Hz")
            else:
                T = 1.0 / f if f > 0 else 0.0
                self._combo_freq.addItem(f"{T:.5g} s")

    def _on_freq_unit(self, _checked: bool) -> None:
        if self._freq_list:
            self._populate_freq_combo()

    def _on_type_changed(self, text: str) -> None:
        t = text.lower()
        self._wgt_freq.setVisible(t in ("depth", "resistivity"))
        self._wgt_comp.setVisible(t in ("depth", "resistivity"))
        self._wgt_depth.setVisible(t == "depth")

    def _on_crs_mode(self, text: str) -> None:
        self._wgt_utm.setVisible("UTM" in text)
        self._wgt_epsg.setVisible("Custom" in text)

    def _resolve_crs(self) -> str:
        mode = self._combo_crs.currentText()
        if "UTM" in mode:
            zone = self._spin_utm.value()
            hemi = "N" if self._radio_utm_n.isChecked() else "S"
            epsg = 32600 + zone if hemi == "N" else 32700 + zone
            return f"EPSG:{epsg}"
        if "Custom" in mode:
            txt = self._edit_epsg.text().strip().upper()
            if txt.startswith("EPSG:"):
                return txt
            try:
                int(txt)
                return f"EPSG:{txt}"
            except ValueError:
                return "EPSG:4326"
        return "EPSG:4326"

    def _current_freq_hz(self) -> float:
        txt = self._combo_freq.currentText().strip()
        try:
            num = float(txt.split()[0])
        except (ValueError, IndexError):
            return 1.0
        if "s" in txt.lower() and "hz" not in txt.lower():
            return 1.0 / num if num > 0 else 1.0
        return num

    def _contour_mode_str(self) -> str:
        return {
            "None": "none", "Lines": "lines",
            "Filled": "filled", "Filled + labels": "filled_labels",
        }.get(self._combo_contour.currentText(), "none")

    def _on_refresh(self) -> None:
        self._map_panel.redraw(
            map_type       = self._combo_type.currentText().lower(),
            target_freq_hz = self._current_freq_hz(),
            target_depth_m = self._spin_depth.value(),
            component      = self._combo_comp.currentText().lower(),
            cmap_name      = self._combo_cmap.currentText(),
            marker_size    = self._spin_ms.value(),
            marker_alpha   = self._slider_alpha.value() / 100.0,
            contour_mode   = self._contour_mode_str(),
            contour_levels = self._spin_levels.value(),
            show_profile   = self._chk_profile.isChecked(),
            show_labels    = self._chk_labels.isChecked(),
            show_grid      = self._chk_grid.isChecked(),
            provider       = self._combo_basemap.currentText(),
            basemap_alpha  = self._slider_base_alpha.value() / 100.0,
            show_cbar      = self._chk_cbar.isChecked(),
            cbar_orient    = "vertical" if self._radio_cbar_v.isChecked() else "horizontal",
            source_crs     = self._resolve_crs(),
            log_scale      = self._chk_log.isChecked(),
        )

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog
        ExportDialog(figure=self._map_panel._canvas.figure, parent=self).exec()
