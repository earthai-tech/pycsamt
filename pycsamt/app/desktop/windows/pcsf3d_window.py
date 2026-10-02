# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Pcsf3DWindow — the desktop app's 3-D/volume panel.

Renders any ``.pcsf``/``.pcsm`` inversion-result file (Occam2D/DUHI
``grid2d``, ModEM ``grid3d``, MARE2DEM ``mesh_unstructured``, or a
``multiline`` stack) as a real fence, block, depth-slice, or surface scene
— the same views and geometry Map View builds, via
:meth:`pycsamt.map.MapView.from_pcsf` and :meth:`~pycsamt.map.MapView.map3d`.
Neither is reimplemented here: both are backend-neutral (plain Plotly
figures, no Dash dependency), so this window's own job is only to load a
file, expose the render options, and display the returned figure through
:class:`~pycsamt.app.desktop.widgets.plotly_view.PlotlyView`.

Layout (v2.6, modelled on Map View's 3-D view)
──────────────────────────────────────────────
::

    ┌ File / Topography / Geology / Boreholes / Structure ┬ toolbar ──────┐
    │ [Load .pcsf…]                                       │ Fence Block … │
    │ Topography: Auto ▾  [Fill gaps] [Online]            │ Depth ▾ Colour│
    │ Geology legend, boreholes, structure overlays       │ ⛰ Topo ▦ Bar …│
    │                                                     ├───────────────┤
    │                                                     │  3-D scene    │
    └─────────────────────────────────────────────────────┴───────────────┘

The toolbar above the scene holds what Map View's 3-D toolbar holds —
mode, depth window, colour map (incl. ``jet_r``), vertical exaggeration,
opacity, topography / terrain / colour bar / legend / stations / labels
toggles, background, camera presets, spin — plus PNG / HTML export and
**Open in Map View**, which hands Map View the *same scene* (mode,
colours, depth, topography, overlays, spin, camera) so nothing has to be
rebuilt there (:func:`~pycsamt.app.desktop.controllers.pcsf_scene.
mapview_state`).  Option changes re-render in place (``Plotly.react``), so
the camera is kept.

Topography: the file's elevations come from its station table, its
per-station block *and* a topography raster sampled at the stations
(:func:`~pycsamt.app.desktop.controllers.pcsf_scene.file_elevations`);
**Auto** takes the file's when it has any, else the survey loaded in the
main window; stations still missing an elevation can be filled along
their line.

Optional overlays (geology legend ``.pcgl.json``/``.csv``, boreholes
``.pcbh.json``, structure ``.pcgs.json``) are drawn into the same scene;
see the overlay methods below for how each is placed.
"""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSize, Qt, QTimer
from PySide6.QtGui import QAction, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QMenu,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.pcsf_scene import (
    CMAPS,
    DEPTH_PRESETS,
    MODES,
    file_elevations,
    fill_line_gaps,
    mapview_state,
)
from pycsamt.app.desktop.widgets.plotly_view import PlotlyView
from pycsamt.app.desktop.windows._base import (
    PanelWindow,
    _icon,
    icon_button,
    make_group,
)

_CMAPS = list(CMAPS)
_MODES = [m for m, _l in MODES]
_BOREHOLE_FAMILIES = ["lithology", "structure", "alteration"]
# station markers (triangles are pycsamt geometry: Plotly 3-D has none)
_STATION_SYMBOLS = [
    ("triangle-down", "▼ Triangle (filled)"),
    ("triangle-down-open", "▽ Triangle (open)"),
    ("diamond", "◆ Diamond"), ("diamond-open", "◇ Diamond (open)"),
    ("circle", "● Circle"), ("circle-open", "○ Circle (open)"),
    ("square", "■ Square"), ("cross", "✚ Cross"), ("x", "✕ X"),
]
# Scene backgrounds: (page, scene) colours and text colour.
_BACKGROUNDS = {
    "White": ("#ffffff", "#ffffff", "#1f2937"),
    "Light grey": ("#f3f4f6", "#f8f9fb", "#1f2937"),
    "Dark": ("#111827", "#0b1220", "#e5e7eb"),
}
_TOPO_SOURCES = [
    ("auto", "Auto (file, else loaded stations)"),
    ("file", "From the file"),
    ("stations", "Loaded stations (main window)"),
    ("topo_file", "Topography file…"),
    ("online", "Fetch online"),
    ("flat", "Flat (no topography)"),
]
# camera eyes (Plotly scene units)
_CAMERAS = {
    "reset": {"eye": {"x": 1.25, "y": 1.25, "z": 1.25}},
    "top": {"eye": {"x": 0.0, "y": 0.0, "z": 2.5},
            "up": {"x": 0, "y": 1, "z": 0}},
    "front": {"eye": {"x": 0.0, "y": -2.4, "z": 0.25}},
    "side": {"eye": {"x": 2.4, "y": 0.0, "z": 0.25}},
    "iso": {"eye": {"x": 1.6, "y": -1.6, "z": 1.0}},
}
# Rotates the camera about the vertical axis while window.__pycsamtSpin,
# and remembers the camera the user leaves (for the Map View hand-off).
_SPIN_JS = """
(function () {
  var gd = document.getElementById('{plot_id}');
  if (!gd) { return; }
  window.__pycsamtSpin = %s;
  gd.on('plotly_relayout', function () {
    try { window.__pycsamtCamera = gd._fullLayout.scene.camera; } catch (e) {}
  });
  var angle = null, radius = null, height = null;
  setInterval(function () {
    if (!window.__pycsamtSpin || !gd._fullLayout || !gd._fullLayout.scene) {
      angle = null; return;
    }
    var eye = gd._fullLayout.scene.camera.eye;
    if (angle === null) {
      radius = Math.sqrt(eye.x * eye.x + eye.y * eye.y) || 1.8;
      angle = Math.atan2(eye.y, eye.x); height = eye.z;
    }
    angle += 0.012;
    Plotly.relayout(gd, {'scene.camera.eye': {
      x: radius * Math.cos(angle), y: radius * Math.sin(angle), z: height}});
  }, 40);
})();
"""

_TB_QSS = (
    "QToolButton { padding: 2px 6px; border: 1px solid #c3c9d4; "
    "border-radius: 5px; background: transparent; }"
    "QToolButton:checked { background: #1864ab; color: white; "
    "border-color: #1864ab; }"
    "QToolButton:hover { border-color: #1864ab; }"
    "QToolButton::menu-indicator { image: none; width: 0px; }"
)


def _tb(text: str, tip: str, *, checkable: bool = False,
        checked: bool = False) -> QToolButton:
    b = QToolButton()
    b.setText(text)
    b.setToolTip(tip)
    b.setCheckable(checkable)
    if checkable:
        b.setChecked(checked)
    b.setStyleSheet(_TB_QSS)
    b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
    return b


def _sep() -> QFrame:
    f = QFrame()
    f.setFrameShape(QFrame.Shape.VLine)
    f.setObjectName("Separator")
    return f


def _cap(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setObjectName("InfoLabel")
    return lbl


class Pcsf3DWindow(PanelWindow):
    """Floating window that loads and renders a ``.pcsf`` file in 3-D."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(
            title="PCSF 3D Viewer",
            session_key="pcsf3d",
            params_width=250,
            icon_name="3d",
            parent=parent,
        )
        self.resize(1180, 780)
        self._view = None  # MapView, built by load_pcsf()
        self._path: Path | None = None
        self._legend = None  # GeologyLegend, built by _on_load_legend()
        self._pcbh_document = None  # PCBHDocument, built by _on_load_boreholes()
        self._pcgs_document = None  # StructModel, built by _on_load_structure()
        self._file_elev: dict = {}  # station elevations stored in the file
        self._file_elev_how = ""
        self._topo_elev: dict = {}  # from a topography file
        self._online_elev: dict = {}
        self._elev_used: dict = {}  # elevations of the last render
        self._camera: dict | None = None  # last camera seen on the page
        self._cam_timer = QTimer(self)
        self._cam_timer.setInterval(1500)
        self._cam_timer.timeout.connect(self._poll_camera)
        self._install_shortcuts()
        self._set_panel_visible(False)  # render controls: toolbar "Controls"

    # ── Controls panel (left, hidden by default: toolbar "Controls") ──

    def _section(self, layout: QVBoxLayout, title: str, *,
                 open_: bool = False) -> QVBoxLayout:
        """A collapsible section; returns the layout to fill."""
        head = QToolButton()
        head.setText(title.replace("&", "&&"))  # "&" alone is a mnemonic
        head.setCheckable(True)
        head.setChecked(open_)
        head.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        head.setArrowType(Qt.ArrowType.DownArrow if open_
                          else Qt.ArrowType.RightArrow)
        head.setAutoRaise(True)
        head.setStyleSheet("QToolButton { font-weight: 600; border: none; "
                           "padding: 3px 0px; }")
        head.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Fixed)
        body = QWidget()
        lay = QVBoxLayout(body)
        lay.setContentsMargins(6, 0, 0, 6)
        lay.setSpacing(4)
        body.setVisible(open_)

        def flip(on):
            body.setVisible(on)
            head.setArrowType(Qt.ArrowType.DownArrow if on
                              else Qt.ArrowType.RightArrow)

        head.toggled.connect(flip)
        layout.addWidget(head)
        layout.addWidget(body)
        self._sections[title] = head
        return lay

    @staticmethod
    def _form(lay: QVBoxLayout) -> QFormLayout:
        f = QFormLayout()
        f.setSpacing(4)
        f.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        lay.addLayout(f)
        return f

    def _dspin(self, lo, hi, value, step, decimals=0, suffix="",
               special="", tip="") -> QDoubleSpinBox:
        w = QDoubleSpinBox()
        w.setRange(lo, hi)
        w.setDecimals(decimals)
        w.setSingleStep(step)
        w.setValue(value)
        if suffix:
            w.setSuffix(suffix)
        if special:
            w.setSpecialValueText(special)
        if tip:
            w.setToolTip(tip)
        w.valueChanged.connect(self._schedule_render)
        return w

    def _combo(self, items, value=None, tip="") -> QComboBox:
        w = QComboBox()
        for data, label in items:
            w.addItem(label, data)
        if value is not None:
            w.setCurrentIndex(max(w.findData(value), 0))
        if tip:
            w.setToolTip(tip)
        w.currentIndexChanged.connect(self._schedule_render)
        return w

    def _check(self, text, checked, tip="") -> QCheckBox:
        w = QCheckBox(text)
        w.setChecked(checked)
        if tip:
            w.setToolTip(tip)
        w.toggled.connect(self._schedule_render)
        return w

    def _build_params(self, layout: QVBoxLayout) -> None:
        self._sections: dict[str, QToolButton] = {}
        self._render_timer = QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.setInterval(250)
        self._render_timer.timeout.connect(self._rerender)

        grp_file, lay_file = make_group("File")
        self._btn_load = icon_button(
            "📂  Load .pcsf / .pcsm…", "3d",
            "Open a .pcsf/.pcsm inversion-result file (Ctrl+O)")
        self._btn_load.clicked.connect(self._on_load)
        lay_file.addWidget(self._btn_load)
        self._file_lbl = QLabel("No file loaded.")
        self._file_lbl.setWordWrap(True)
        self._file_lbl.setObjectName("InfoLabel")
        lay_file.addWidget(self._file_lbl)
        layout.addWidget(grp_file)

        # ── scene (what is drawn) ──────────────────────────────────────
        lay = self._section(layout, "Scene", open_=True)
        f = self._form(lay)
        self._combo_mode = QComboBox()
        for key, label in MODES:
            self._combo_mode.addItem(label, key)
        self._combo_mode.setToolTip("Scene type  (keys 1–4)")
        self._combo_mode.currentIndexChanged.connect(self._on_mode_changed)
        f.addRow("Mode:", self._combo_mode)
        self._combo_depth = QComboBox()
        for value, label in DEPTH_PRESETS:
            self._combo_depth.addItem(label, value)
        self._combo_depth.addItem("Custom…", -1.0)
        self._combo_depth.setToolTip("Deepest level shown")
        self._combo_depth.currentIndexChanged.connect(self._on_depth_preset)
        f.addRow("Depth:", self._combo_depth)
        self._combo_cmap = QComboBox()
        self._combo_cmap.addItems(_CMAPS)
        self._combo_cmap.setToolTip("Colour map (jet_r: conductive red)")
        self._combo_cmap.currentIndexChanged.connect(
            lambda _i: self._rerender())
        f.addRow("Colours:", self._combo_cmap)
        self._spin_ve = QDoubleSpinBox()
        self._spin_ve.setRange(0.0, 50.0)
        self._spin_ve.setDecimals(1)
        self._spin_ve.setSingleStep(0.5)
        self._spin_ve.setValue(0.0)
        self._spin_ve.setSpecialValueText("Auto")
        self._spin_ve.setSuffix(" ×")
        self._spin_ve.setToolTip(
            "Vertical exaggeration (keys +/−). A survey tens of km wide "
            "over a few km of model looks like a thin sheet at true scale; "
            "Auto makes the depth about half the horizontal extent.")
        self._spin_ve.editingFinished.connect(self._rerender)
        f.addRow("Vertical:", self._spin_ve)

        # ── stations & labels ──────────────────────────────────────────
        lay = self._section(layout, "Stations & labels")
        f = self._form(lay)
        self._combo_sta_symbol = self._combo(
            _STATION_SYMBOLS, "triangle-down",
            "Marker drawn at each station (triangles are real 3-D "
            "geometry, like a section's station marks)")
        f.addRow("Marker:", self._combo_sta_symbol)
        self._spin_sta_size = self._dspin(2, 30, 6, 1,
                                          tip="Marker size (px)")
        f.addRow("Size:", self._spin_sta_size)
        self._sta_color = "#111111"
        self._btn_sta_color = QPushButton()
        self._btn_sta_color.setToolTip("Marker colour")
        self._btn_sta_color.clicked.connect(self._pick_station_color)
        self._paint_color_button()
        f.addRow("Colour:", self._btn_sta_color)
        self._spin_sta_max = self._dspin(
            0, 500, 0, 1, special="All",
            tip="Markers per line, evenly spaced (first and last kept)")
        f.addRow("Per line:", self._spin_sta_max)
        self._spin_label_angle = self._dspin(
            -90, 90, 0, 15, suffix="°",
            tip="Rotate station labels (0 = horizontal; 45/90 for a "
                "crowded line)")
        f.addRow("Label angle:", self._spin_label_angle)
        self._combo_label_density = self._combo(
            [(v, f"{int(v * 100)} %") for v in (1.0, 0.75, 0.5, 0.25, 0.1,
                                               0.05)], 1.0,
            "Share of stations labelled per line (markers are unaffected)")
        f.addRow("Labelled:", self._combo_label_density)
        self._edit_label_names = QLineEdit()
        self._edit_label_names.setPlaceholderText("only these, e.g. S01, S12")
        self._edit_label_names.setToolTip(
            "Label only the named stations (comma separated)")
        self._edit_label_names.editingFinished.connect(self._schedule_render)
        f.addRow("Named:", self._edit_label_names)
        self._chk_line_labels = self._check("Line names", True,
                                            "Name each profile line")
        lay.addWidget(self._chk_line_labels)

        # ── depth & units ──────────────────────────────────────────────
        lay = self._section(layout, "Depth & units")
        f = self._form(lay)
        self._spin_depth_lo = self._dspin(0, 1e7, 0, 100, suffix=" m",
                                          tip="Shallowest level shown")
        f.addRow("From:", self._spin_depth_lo)
        self._spin_depth = self._dspin(
            0, 1e7, 0, 500, suffix=" m", special="Full model",
            tip="Deepest level shown (the toolbar presets set it)")
        self._spin_depth.valueChanged.connect(self._sync_depth_preset)
        f.addRow("To:", self._spin_depth)
        self._combo_depth_unit = self._combo(
            [("m", "metres"), ("km", "kilometres")], "m",
            "Unit of the depth axis")
        f.addRow("Depth axis:", self._combo_depth_unit)
        self._combo_x_unit = self._combo(
            [("m", "metres"), ("km", "kilometres")], "m",
            "Unit of the horizontal axes")
        f.addRow("Distance:", self._combo_x_unit)

        # ── resistivity & colours ──────────────────────────────────────
        lay = self._section(layout, "Resistivity & colours")
        f = self._form(lay)
        self._combo_scale = self._combo(
            [("log", "Logarithmic"), ("linear", "Linear")], "log",
            "Colour scale of resistivity")
        f.addRow("Scale:", self._combo_scale)
        self._spin_vmin = self._dspin(0, 1e9, 0, 1, 2, " Ω·m", "Auto",
                                      "Colour range minimum")
        self._spin_vmax = self._dspin(0, 1e9, 0, 100, 2, " Ω·m", "Auto",
                                      "Colour range maximum")
        f.addRow("Colour min:", self._spin_vmin)
        f.addRow("Colour max:", self._spin_vmax)
        self._spin_plo = self._dspin(0, 49, 2, 1, 0, " %",
                                     tip="Auto range: lower percentile")
        self._spin_phi = self._dspin(51, 100, 98, 1, 0, " %",
                                     tip="Auto range: upper percentile")
        f.addRow("Clip low:", self._spin_plo)
        f.addRow("Clip high:", self._spin_phi)
        self._spin_rho_lo = self._dspin(0, 1e9, 0, 10, 2, " Ω·m", "Off",
                                        "Show only cells at least this "
                                        "resistive")
        self._spin_rho_hi = self._dspin(0, 1e9, 0, 100, 2, " Ω·m", "Off",
                                        "Show only cells at most this "
                                        "resistive")
        f.addRow("Keep from:", self._spin_rho_lo)
        f.addRow("Keep to:", self._spin_rho_hi)
        self._spin_rho_cutoff = self._dspin(
            0, 1e12, 0, 1000, 0, " Ω·m", "Off",
            "Hide values above this (air / padding fill)")
        f.addRow("Cut above:", self._spin_rho_cutoff)
        self._chk_contours = self._check("Contour lines", False,
                                         "Iso-resistivity contours")
        lay.addWidget(self._chk_contours)

        # ── view & geometry ────────────────────────────────────────────
        lay = self._section(layout, "View & geometry")
        f = self._form(lay)
        self._combo_aspect = self._combo(
            [("ve", "Vertical exaggeration (toolbar)"),
             ("data", "True proportions"), ("cube", "Equal cube"),
             ("manual", "Stretch to box")], "ve",
            "How the scene box is proportioned")
        f.addRow("Proportions:", self._combo_aspect)
        self._spin_opacity = self._dspin(0.1, 1.0, 0.85, 0.05, 2,
                                         tip="Model opacity")
        f.addRow("Opacity:", self._spin_opacity)
        self._spin_nslices = self._dspin(2, 40, 8, 1,
                                         tip="Depth slices (depth mode)")
        f.addRow("Depth slices:", self._spin_nslices)
        self._spin_surfaces = self._dspin(
            2, 40, 12, 1, tip="Iso-surfaces (iso-surface mode)")
        f.addRow("Iso-surfaces:", self._spin_surfaces)
        self._spin_spacing = self._dspin(0.1, 10, 1.0, 0.1, 1, " ×",
                                         tip="Spread between lines")
        f.addRow("Line spacing:", self._spin_spacing)
        self._spin_azimuth = self._dspin(-180, 180, 0, 5, 0, "°",
                                         tip="Rotate the survey frame")
        f.addRow("Azimuth:", self._spin_azimuth)
        self._chk_smooth = self._check("Smooth sections", True,
                                       "Spline-resample fence panels")
        lay.addWidget(self._chk_smooth)
        self._combo_section_res = self._combo(
            [(60, "Coarse"), (100, "Normal"), (160, "Fine"),
             (240, "Ultra")], 100, "Resampled section resolution")
        f2 = self._form(lay)
        f2.addRow("Resolution:", self._combo_section_res)
        self._combo_vol_smooth = self._combo(
            [(0.0, "Off (raw cells)"), (0.8, "Light"), (1.5, "Medium"),
             (2.5, "Strong"), (4.0, "Very strong")], 0.0,
            "Smooth the block / iso volume")
        f2.addRow("Volume smooth:", self._combo_vol_smooth)
        self._combo_bg = QComboBox()
        self._combo_bg.addItems(list(_BACKGROUNDS))
        self._combo_bg.setToolTip("Scene background (white for figures)")
        self._combo_bg.currentIndexChanged.connect(lambda _i: self._rerender())
        f2.addRow("Background:", self._combo_bg)

        # ── topography ─────────────────────────────────────────────────
        lay = self._section(layout, "Topography", open_=True)
        self._combo_topo = QComboBox()
        for key, label in _TOPO_SOURCES:
            self._combo_topo.addItem(label, key)
        self._combo_topo.setToolTip(
            "Surface elevation of the stations, draped on the model:\n"
            "the file's own (station table, topography block or raster), "
            "the survey loaded in the main window (matched by station "
            "name), a .csv/.bln/.stn file, or an online lookup.")
        self._combo_topo.currentIndexChanged.connect(self._on_topo_source)
        lay.addWidget(self._combo_topo)
        self._btn_topo_file = icon_button(
            "📂  Load topography file…", "topography",
            "A .csv/.bln/.stn file of station elevations")
        self._btn_topo_file.clicked.connect(self._on_load_topo_file)
        self._btn_topo_file.setVisible(False)
        lay.addWidget(self._btn_topo_file)
        self._chk_fill_gaps = QCheckBox("Fill gaps along lines")
        self._chk_fill_gaps.setChecked(True)
        self._chk_fill_gaps.setToolTip(
            "Stations without an elevation take one interpolated between\n"
            "their neighbours on the same line (no holes in the surface).")
        self._chk_fill_gaps.toggled.connect(lambda _on: self._rerender())
        lay.addWidget(self._chk_fill_gaps)
        self._chk_fetch_elev = QCheckBox("Fetch online on load")
        self._chk_fetch_elev.setChecked(False)
        self._chk_fetch_elev.setToolTip(
            "A ModEM/MARE2DEM result carries no real elevation of its own.\n"
            "When checked, a best-effort online lookup runs on load — off by\n"
            "default so this window works without a network connection.")
        lay.addWidget(self._chk_fetch_elev)
        self._topo_lbl = QLabel("")
        self._topo_lbl.setWordWrap(True)
        self._topo_lbl.setObjectName("InfoLabel")
        lay.addWidget(self._topo_lbl)

        # ── overlays ───────────────────────────────────────────────────
        lay = self._section(layout, "Geology legend")
        self._btn_load_legend = icon_button(
            "📂  Load legend…", "interpret",
            "Open a .pcgl.json legend or a plain CSV"
            " (name, rho_min, rho_max[, color])")
        self._btn_load_legend.clicked.connect(self._on_load_legend)
        lay.addWidget(self._btn_load_legend)
        self._legend_lbl = QLabel("No legend loaded.")
        self._legend_lbl.setWordWrap(True)
        self._legend_lbl.setObjectName("InfoLabel")
        lay.addWidget(self._legend_lbl)
        self._chk_show_geology = QCheckBox("Show geology overlay")
        self._chk_show_geology.setChecked(True)
        self._chk_show_geology.setEnabled(False)
        self._chk_show_geology.toggled.connect(lambda _on: self._rerender())
        lay.addWidget(self._chk_show_geology)
        self._chk_pattern_fill = QCheckBox("Pattern fill (fence/depth only)")
        self._chk_pattern_fill.setChecked(False)
        self._chk_pattern_fill.setEnabled(False)
        self._chk_pattern_fill.setToolTip(
            "Texture each band with its assigned pattern instead of a flat\n"
            "colour. Block/surface modes always render solid — Plotly's\n"
            "Volume/Isosurface traces have no per-cell texture support.")
        self._chk_pattern_fill.toggled.connect(lambda _on: self._rerender())
        lay.addWidget(self._chk_pattern_fill)

        lay = self._section(layout, "Boreholes")
        self._btn_load_boreholes = icon_button(
            "📂  Load PCBH…", "3d", "Open a .pcbh.json borehole document")
        self._btn_load_boreholes.clicked.connect(self._on_load_boreholes)
        lay.addWidget(self._btn_load_boreholes)
        self._boreholes_lbl = QLabel("No boreholes loaded.")
        self._boreholes_lbl.setWordWrap(True)
        self._boreholes_lbl.setObjectName("InfoLabel")
        lay.addWidget(self._boreholes_lbl)
        self._combo_borehole_family = QComboBox()
        self._combo_borehole_family.addItems(_BOREHOLE_FAMILIES)
        self._combo_borehole_family.setToolTip("Interval family drawn")
        self._combo_borehole_family.currentIndexChanged.connect(
            lambda _i: self._rerender())
        lay.addWidget(self._combo_borehole_family)
        self._chk_show_boreholes = QCheckBox("Show boreholes")
        self._chk_show_boreholes.setChecked(True)
        self._chk_show_boreholes.setEnabled(False)
        self._chk_show_boreholes.toggled.connect(lambda _on: self._rerender())
        lay.addWidget(self._chk_show_boreholes)

        lay = self._section(layout, "Structure")
        self._btn_load_structure = icon_button(
            "📂  Load PCGS…", "3d",
            "Open a .pcgs.json structural-geology document")
        self._btn_load_structure.clicked.connect(self._on_load_structure)
        lay.addWidget(self._btn_load_structure)
        self._structure_lbl = QLabel("No structure loaded.")
        self._structure_lbl.setWordWrap(True)
        self._structure_lbl.setObjectName("InfoLabel")
        lay.addWidget(self._structure_lbl)
        self._chk_show_structure = QCheckBox("Show structure overlay")
        self._chk_show_structure.setChecked(True)
        self._chk_show_structure.setEnabled(False)
        self._chk_show_structure.toggled.connect(lambda _on: self._rerender())
        lay.addWidget(self._chk_show_structure)

        btn_reset = QPushButton("Reset render settings")
        btn_reset.setObjectName("FileListBtn")
        btn_reset.clicked.connect(self.reset_render_settings)
        layout.addWidget(btn_reset)
        keys = QLabel(
            "<small style='color:#6b7280'><b>Keys</b> 1–4 mode · Space spin"
            " · R reset view · T topography · C colour bar · L legend · "
            "S stations · +/− vertical · F5 render · F9 controls · "
            "Ctrl+P PNG · Ctrl+M Map View</small>")
        keys.setWordWrap(True)
        layout.addWidget(keys)

    def _schedule_render(self, *_a) -> None:
        if self._view is not None:
            self._render_timer.start()

    def _pick_station_color(self) -> None:
        from PySide6.QtGui import QColor
        from PySide6.QtWidgets import QColorDialog

        c = QColorDialog.getColor(QColor(self._sta_color), self,
                                  "Station marker colour")
        if c.isValid():
            self.set_station_color(c.name())

    def set_station_color(self, colour: str) -> None:
        self._sta_color = colour
        self._paint_color_button()
        self._schedule_render()

    def _paint_color_button(self) -> None:
        self._btn_sta_color.setText(self._sta_color)
        self._btn_sta_color.setStyleSheet(
            f"QPushButton {{ background: {self._sta_color}; color: "
            f"{'#ffffff' if self._is_dark(self._sta_color) else '#111111'};"
            " border: 1px solid #9ca3af; border-radius: 4px; padding: 2px; }")

    @staticmethod
    def _is_dark(colour: str) -> bool:
        try:
            r, g, b = (int(colour[i:i + 2], 16) for i in (1, 3, 5))
        except (ValueError, IndexError):
            return False
        return (0.299 * r + 0.587 * g + 0.114 * b) < 128

    def _sync_depth_preset(self, value: float) -> None:
        i = self._combo_depth.findData(float(value))
        self._combo_depth.blockSignals(True)
        self._combo_depth.setCurrentIndex(
            i if i >= 0 else self._combo_depth.findData(-1.0))
        self._combo_depth.blockSignals(False)

    def reset_render_settings(self) -> None:
        """Back to the defaults (the Map View look)."""
        defaults = (
            (self._spin_sta_size, 6), (self._spin_sta_max, 0),
            (self._spin_label_angle, 0), (self._spin_depth_lo, 0),
            (self._spin_depth, 0), (self._spin_vmin, 0), (self._spin_vmax, 0),
            (self._spin_plo, 2), (self._spin_phi, 98),
            (self._spin_rho_lo, 0), (self._spin_rho_hi, 0),
            (self._spin_rho_cutoff, 0), (self._spin_opacity, 0.85),
            (self._spin_nslices, 8), (self._spin_surfaces, 12),
            (self._spin_spacing, 1.0), (self._spin_azimuth, 0),
            (self._spin_ve, 0))
        for w, v in defaults:
            w.blockSignals(True)
            w.setValue(v)
            w.blockSignals(False)
        for combo, data in ((self._combo_sta_symbol, "triangle-down"),
                            (self._combo_label_density, 1.0),
                            (self._combo_depth_unit, "m"),
                            (self._combo_x_unit, "m"),
                            (self._combo_scale, "log"),
                            (self._combo_aspect, "ve"),
                            (self._combo_section_res, 100),
                            (self._combo_vol_smooth, 0.0),
                            (self._combo_depth, 0.0)):
            combo.blockSignals(True)
            combo.setCurrentIndex(max(combo.findData(data), 0))
            combo.blockSignals(False)
        for chk, on in ((self._chk_line_labels, True),
                        (self._chk_contours, False),
                        (self._chk_smooth, True)):
            chk.blockSignals(True)
            chk.setChecked(on)
            chk.blockSignals(False)
        self._edit_label_names.clear()
        self._sta_color = "#111111"
        self._paint_color_button()
        self._rerender()

    def render_options(self) -> dict:
        """The scene's ``map3d`` keyword arguments (also what Map View is
        handed, see :func:`~pycsamt.app.desktop.controllers.pcsf_scene.
        mapview_controls`)."""
        def pos(v):
            return float(v) if v and v > 0 else None

        lo, hi = float(self._spin_depth_lo.value()), self._depth_max()
        vmin, vmax = pos(self._spin_vmin.value()), pos(self._spin_vmax.value())
        rlo, rhi = pos(self._spin_rho_lo.value()), pos(self._spin_rho_hi.value())
        names = tuple(n.strip() for n in
                      self._edit_label_names.text().split(",") if n.strip())
        aspect = self._combo_aspect.currentData()
        opts = dict(
            mode=self._mode(),
            cmap=self._combo_cmap.currentText(),
            opacity=float(self._spin_opacity.value()),
            topography=self._topo_source() != "flat"
            and self._chk_topo.isChecked(),
            show_terrain=self._chk_terrain.isChecked(),
            show_stations=self._chk_stations.isChecked(),
            station_labels=self._chk_station_labels.isChecked(),
            station_symbol=self._combo_sta_symbol.currentData(),
            station_size=int(self._spin_sta_size.value()),
            station_color=self._sta_color,
            max_stations=int(self._spin_sta_max.value()) or None,
            station_label_angle=float(self._spin_label_angle.value()),
            station_label_fraction=float(
                self._combo_label_density.currentData()),
            station_label_names=names or None,
            show_labels=self._chk_line_labels.isChecked(),
            depth_range=(lo, hi) if hi > 0 else ((lo, 1e9) if lo > 0
                                                 else None),
            depth_unit=self._combo_depth_unit.currentData(),
            x_unit=self._combo_x_unit.currentData(),
            log_color=self._combo_scale.currentData() == "log",
            value_range=(vmin, vmax) if vmin and vmax else None,
            crange_percentile=(float(self._spin_plo.value()),
                               float(self._spin_phi.value())),
            rho_range=(rlo or 0.0, rhi or 1e12) if (rlo or rhi) else None,
            rho_display_max=pos(self._spin_rho_cutoff.value()),
            show_contours=self._chk_contours.isChecked(),
            aspectmode="data" if aspect == "ve" else aspect,
            vertical_exaggeration=float(self._spin_ve.value())
            if aspect == "ve" else None,
            n_slices=int(self._spin_nslices.value()),
            surface_count=int(self._spin_surfaces.value()),
            line_spacing=float(self._spin_spacing.value()),
            azimuth=float(self._spin_azimuth.value()),
            smooth_sections=self._chk_smooth.isChecked(),
            section_res=int(self._combo_section_res.currentData()),
            volume_smoothing=float(self._combo_vol_smooth.currentData()),
        )
        geo = self._geology_kwargs()
        if geo:
            opts.update(geo)
        opts["geology_legend"] = self._chk_legend.isChecked()
        return opts

    # ── Content (right): one compact toolbar + scene ──────────────────

    def _menu_button(self, icon: str, tip: str) -> tuple:
        b = self._icon_button(icon, tip)
        b.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QMenu(b)
        b.setMenu(menu)
        return b, menu

    @staticmethod
    def _icon_button(icon: str, tip: str, *, checkable: bool = False,
                     checked: bool = False, text: str = "") -> QToolButton:
        """A square icon shortcut (text only if the icon is missing)."""
        b = _tb(text, tip, checkable=checkable, checked=checked)
        ic = _icon(icon)
        if not ic.isNull():
            b.setIcon(ic)
            b.setIconSize(QSize(18, 18))
            b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        elif not text:
            b.setText(tip.split("(")[0].strip()[:10])
        return b

    def _build_content(self, layout: QVBoxLayout) -> None:
        bar = QWidget()
        bar.setObjectName("SceneToolbar")
        row = QHBoxLayout(bar)
        row.setContentsMargins(6, 3, 6, 3)
        row.setSpacing(3)

        self._btn_panel = self._icon_button(
            "tools", "Render controls: scene, stations, labels, depth, "
            "resistivity, view, topography, overlays  (F9)",
            checkable=True, checked=False)
        self._btn_panel.toggled.connect(self._set_panel_visible)
        row.addWidget(self._btn_panel)
        self._btn_open = self._icon_button(
            "open", "Open a .pcsf/.pcsm model  (Ctrl+O)")
        self._btn_open.clicked.connect(self._on_load)
        row.addWidget(self._btn_open)
        row.addWidget(_sep())

        # layers: every show/hide toggle in one menu
        self._btn_layers, menu = self._menu_button(
            "layered-model",
            "Layers: topography, terrain, colour bar, legend, stations")

        def toggle(text, tip, checked, shortcut=""):
            act = QAction(text, self)
            act.setCheckable(True)
            act.setChecked(checked)
            act.setToolTip(tip)
            if shortcut:
                act.setText(f"{text}\t{shortcut}")
            act.toggled.connect(lambda _on: self._rerender())
            menu.addAction(act)
            return act

        self._chk_topo = toggle("Topography", "Drape on the stations' "
                                "elevations", True, "T")
        self._chk_terrain = toggle("Terrain line", "Ground line on each "
                                   "section", True)
        menu.addSeparator()
        self._chk_colorbar = toggle("Colour bar", "Show the colour bar",
                                    True, "C")
        self._chk_legend = toggle("Legend", "Show the legend", True, "L")
        menu.addSeparator()
        self._chk_stations = toggle("Station markers", "Show stations",
                                    True, "S")
        self._chk_station_labels = toggle("Station labels", "Station names",
                                          False)
        row.addWidget(self._btn_layers)

        self._btn_camera, cmenu = self._menu_button(
            "3d", "Camera: reset, top, front, side, oblique")
        for key, label in (("reset", "Reset view\tR"), ("top", "Top"),
                           ("front", "Front (look north)"),
                           ("side", "Side (look west)"),
                           ("iso", "Oblique 3-D")):
            cmenu.addAction(label, lambda k=key: self.set_camera(k))
        row.addWidget(self._btn_camera)
        self._btn_spin = self._icon_button(
            "recompute", "Spin the scene  (Space)", checkable=True)
        self._btn_spin.toggled.connect(self._on_spin)
        row.addWidget(self._btn_spin)
        row.addWidget(_sep())

        self._btn_out, omenu = self._menu_button(
            "export", "Export: PNG (as shown) or interactive HTML")
        self._btn_png = omenu.addAction("PNG image (as shown)…\tCtrl+P",
                                        self._on_png)
        self._btn_export = omenu.addAction("Interactive HTML…",
                                           self._on_export)
        self._btn_png.setEnabled(False)
        self._btn_export.setEnabled(False)
        row.addWidget(self._btn_out)
        self._btn_mapview = self._icon_button(
            "map-view",
            "Continue in Map View (browser) with this exact scene: mode, "
            "colours, depth, topography, overlays, spin and camera  "
            "(Ctrl+M)")
        self._btn_mapview.clicked.connect(self._on_open_mapview)
        self._btn_mapview.setEnabled(False)
        row.addWidget(self._btn_mapview)
        row.addStretch(1)
        self._scene_lbl = _cap("")
        row.addWidget(self._scene_lbl)
        layout.addWidget(bar)

        self._plotly_view = PlotlyView(self)
        layout.addWidget(self._plotly_view, 1)
        # on-scene shortcut, like the plot canvases' refresh button
        self._btn_render = self._plotly_view.add_overlay_button(
            _icon("reset"), "Hard render: re-read the file and rebuild the "
            "scene  (F5)", self.hard_render, text="↻")
        self._btn_render.setEnabled(False)
        self._status_lbl = QLabel("")
        self._status_lbl.setObjectName("InfoLabel")
        self._status_lbl.setWordWrap(True)
        self._status_lbl.setContentsMargins(6, 2, 6, 2)
        layout.addWidget(self._status_lbl)

    def hard_render(self) -> None:
        """Re-read the file and rebuild the page from scratch (the normal
        render swaps the figure in place and keeps the camera)."""
        if self._path is None:
            return
        try:
            from pycsamt.map import MapView

            self._view = MapView.from_pcsf(
                str(self._path),
                fetch_elevation=self._chk_fetch_elev.isChecked())
        except Exception as exc:
            self._status_lbl.setText(f"Could not re-read the file: {exc}")
            return
        self._plotly_view.reload_next()
        self._on_render()

    def _set_panel_visible(self, on: bool) -> None:
        panel = self._splitter.widget(0)
        if panel is not None:
            panel.setVisible(on)

    def _install_shortcuts(self) -> None:
        def key(seq, fn):
            QShortcut(QKeySequence(seq), self, activated=fn)

        for i in range(len(MODES)):
            key(str(i + 1), lambda i=i: self._combo_mode.setCurrentIndex(i))
        key("Space", lambda: self._btn_spin.toggle())
        key("R", lambda: self.set_camera("reset"))
        key("T", lambda: self._chk_topo.toggle())
        key("C", lambda: self._chk_colorbar.toggle())
        key("L", lambda: self._chk_legend.toggle())
        key("S", lambda: self._chk_stations.toggle())
        key("+", lambda: self._step_ve(+1))
        key("=", lambda: self._step_ve(+1))
        key("-", lambda: self._step_ve(-1))
        key("F5", self.hard_render)
        key("F9", lambda: self._btn_panel.toggle())
        key("Ctrl+O", self._on_load)
        key("Ctrl+P", self._on_png)
        key("Ctrl+M", self._on_open_mapview)

    # ── Public API ────────────────────────────────────────────────────

    def load_pcsf(self, path: str) -> None:
        """Load and render *path* — the programmatic entry point behind
        the "Load .pcsf file…" button, also usable directly (e.g. by
        tests or a future "open with" integration)."""
        from pycsamt.map import MapView

        self._path = Path(path)
        self._view = MapView.from_pcsf(
            path, fetch_elevation=self._chk_fetch_elev.isChecked())
        try:
            file_elev, how = file_elevations(path)
        except Exception:
            file_elev, how = {}, ""
        # what MapView already read, plus a raster / block it leaves out
        own = {s.id: s.elevation for s in self._view.data.stations
               if s.elevation is not None}
        if any(abs(float(e)) > 0 for e in own.values()):
            file_elev = {**file_elev, **own}
        self._file_elev = {k: float(e) for k, e in file_elev.items()}
        self._file_elev_how = how
        self._topo_elev = {}
        self._online_elev = {}
        self._camera = None
        has_topo = any(abs(float(e)) > 0 for e in self._file_elev.values())
        self._file_lbl.setText(
            f"{self._path.name}\n"
            f"lines={len(self._view.lines)}  stations={self._view.n_stations}"
            f"  geo={'yes' if self._view.has_geo else 'no'}\n"
            f"topography: {how if has_topo else 'none in the file'}")
        self._combo_topo.blockSignals(True)
        self._combo_topo.setCurrentIndex(self._combo_topo.findData("auto"))
        self._combo_topo.blockSignals(False)
        self._btn_topo_file.setVisible(False)
        for b in (self._btn_render, self._btn_export, self._btn_png,
                  self._btn_mapview):
            b.setEnabled(True)
        self._on_render()

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_load(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load PCSF file", "",
            "PCSF / PCSM files (*.pcsf *.pcsm *.pcsm.gz)")
        if not path:
            return
        try:
            self.load_pcsf(path)
        except Exception as exc:
            self._status_lbl.setText(f"Load error: {exc}")
            self._btn_render.setEnabled(False)
            self._btn_export.setEnabled(False)

    def _on_load_legend(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load geology legend",
            "",
            "Geology legend (*.pcgl.json *.csv)",
        )
        if not path:
            return
        try:
            self._legend = self._read_legend(path)
            n = len(self._legend.entries)
            self._legend_lbl.setText(f"{Path(path).name}\n{n} band(s)")
            self._chk_show_geology.setEnabled(True)
            self._chk_pattern_fill.setEnabled(True)
            if self._view is not None:
                self._on_render()
        except Exception as exc:
            self._legend = None
            self._legend_lbl.setText(f"Legend load error: {exc}")
            self._chk_show_geology.setEnabled(False)
            self._chk_pattern_fill.setEnabled(False)

    @staticmethod
    def _read_legend(path: str):
        from pycsamt.format.geology import legend_from_csv, read_legend

        if str(path).lower().endswith(".csv"):
            return legend_from_csv(path)
        return read_legend(path)

    def _geology_kwargs(self) -> dict:
        """``VolumeMapOptions`` overrides for the loaded legend, or ``{}``."""
        if self._legend is None or not self._chk_show_geology.isChecked():
            return {}
        from pycsamt.app._geology import (
            geology_bands_from_store,
            geology_pattern_stencils_from_store,
        )

        bands = geology_bands_from_store(self._legend)
        if bands is None:
            return {}
        kwargs = {"geology": bands}
        if self._chk_pattern_fill.isChecked():
            kwargs["geology_fill"] = "pattern"
            kwargs["geology_patterns"] = geology_pattern_stencils_from_store(
                self._legend
            )
        return kwargs

    def _on_load_boreholes(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load PCBH document", "", "PCBH files (*.pcbh.json)"
        )
        if not path:
            return
        try:
            from pycsamt.format.borehole import read_pcbh

            self._pcbh_document = read_pcbh(path)
            n = len(self._pcbh_document.boreholes)
            self._boreholes_lbl.setText(f"{Path(path).name}\n{n} hole(s)")
            self._chk_show_boreholes.setEnabled(True)
            if self._view is not None:
                self._on_render()
        except Exception as exc:
            self._pcbh_document = None
            self._boreholes_lbl.setText(f"Load error: {exc}")
            self._chk_show_boreholes.setEnabled(False)

    def _borehole_alignment_kwargs(self, ids, lats, lons, lines, uv):
        """``(datum, surface, offset_shift)`` matching the fence renderer's
        own cross-strike normalisation -- ported from
        ``pycsamt.app.mapview._render._add_borehole_scene``, the only place
        this recipe existed before this window. A borehole collar's raw
        projected cross-strike position must go through the *same* shift
        the volume builder applied to its line panels
        (``pycsamt.map.geometry.normalize_offsets``) or the hole floats in
        front of / behind its own line.
        """
        import numpy as np

        from pycsamt.map.borehole_align import surface_from_sections

        elevs = [
            float(s.elevation) if s.elevation is not None else float("nan")
            for s in self._view.data.stations
            if s.id in ids
        ]
        paired = [
            (uv[i][0], e)
            for i, e in zip(ids, elevs)
            if i in uv and np.isfinite(e)
        ]
        surface = (
            surface_from_sections([([p[0] for p in paired], [p[1] for p in paired])])
            if len(paired) >= 2
            else None
        )
        datum = "surface" if surface is not None else "zero"

        line_v: dict = {}
        for sid, ln in zip(ids, lines):
            if sid in uv:
                line_v.setdefault(str(ln), []).append(uv[sid][1])
        medians = [float(np.nanmedian(v)) for v in line_v.values() if v]
        offset_shift = float(min(medians)) if medians else 0.0
        return datum, surface, offset_shift

    def _borehole_traces(self) -> list:
        """Plotly traces for the loaded PCBH document, or ``[]``.

        Only the "not enough geo-referenced stations to build a scene
        frame" case returns ``[]`` quietly (there is nothing wrong, just
        nothing placeable). Any other failure propagates to
        :meth:`_on_render`'s own try/except so a bad borehole document
        surfaces as a real status-bar error instead of silently rendering
        without boreholes.
        """
        if (
            self._pcbh_document is None
            or not self._chk_show_boreholes.isChecked()
            or self._view is None
        ):
            return []
        from pycsamt.app._borehole import scene_borehole_traces
        from pycsamt.map.borehole_align import align_boreholes_to_scene
        from pycsamt.map.geometry import survey_frame, survey_uv

        stations = [
            s
            for s in self._view.data.stations
            if s.latitude is not None and s.longitude is not None
        ]
        if len(stations) < 2:
            return []
        ids = [s.id for s in stations]
        lats = [float(s.latitude) for s in stations]
        lons = [float(s.longitude) for s in stations]
        lines = [s.line or "line" for s in stations]

        frame = survey_frame(lats, lons, lines)
        uv = survey_uv(ids, lats, lons, lines)
        datum, surface, offset_shift = self._borehole_alignment_kwargs(
            ids, lats, lons, lines, uv
        )

        alignment = align_boreholes_to_scene(
            self._pcbh_document,
            frame,
            family=self._combo_borehole_family.currentText(),
            datum=datum,
            surface=surface,
            offset_shift=offset_shift,
        )
        return scene_borehole_traces(
            alignment, as_tubes=True, opacity=0.9, show_labels=True
        )

    def _on_load_structure(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load PCGS document",
            "",
            "PCGS files (*.pcgs.json)",
        )
        if not path:
            return
        try:
            from pycsamt.format.structure import read_structure

            self._pcgs_document = read_structure(path)
            n_planar = len(self._pcgs_document.planar)
            n_linear = len(self._pcgs_document.linear)
            n_faults = len(self._pcgs_document.faults)
            self._structure_lbl.setText(
                f"{Path(path).name}\n"
                f"{n_faults} fault(s), {n_planar} planar, {n_linear} linear"
            )
            self._chk_show_structure.setEnabled(True)
            if self._view is not None:
                self._on_render()
        except Exception as exc:
            self._pcgs_document = None
            self._structure_lbl.setText(f"Load error: {exc}")
            self._chk_show_structure.setEnabled(False)

    def _structure_line_offsets(self, ids, lines, uv):
        """``(line_offsets, default_line)`` -- each line's own already-
        normalised cross-strike scene offset, the same per-line-median math
        the volume builder itself uses (``pycsamt.map.volume._line_offset``)
        and the borehole overlay's ``offset_shift`` derives from. Structural
        items need the *per-line* value directly (their ``x`` is already a
        profile position on their own line), not a single scalar shift.
        """
        import numpy as np

        line_v: dict = {}
        for sid, ln in zip(ids, lines):
            if sid in uv:
                line_v.setdefault(str(ln), []).append(uv[sid][1])
        medians = {
            ln: float(np.nanmedian(v)) for ln, v in line_v.items() if v
        }
        if not medians:
            return {}, None
        offset_shift = min(medians.values())
        line_offsets = {ln: (v - offset_shift) for ln, v in medians.items()}
        return line_offsets, next(iter(line_offsets))

    def _structure_traces(self) -> list:
        """Plotly traces for the loaded PCGS document, or ``[]``.

        Same "quiet vs. propagated" failure split as :meth:`_borehole_traces`.
        """
        if (
            self._pcgs_document is None
            or not self._chk_show_structure.isChecked()
            or self._view is None
        ):
            return []
        from pycsamt.app._structure import structure_scene_traces
        from pycsamt.map.geometry import survey_uv

        stations = [
            s
            for s in self._view.data.stations
            if s.latitude is not None and s.longitude is not None
        ]
        if len(stations) < 2:
            return []
        ids = [s.id for s in stations]
        lats = [float(s.latitude) for s in stations]
        lons = [float(s.longitude) for s in stations]
        lines = [s.line or "line" for s in stations]

        uv = survey_uv(ids, lats, lons, lines)
        _, surface, _ = self._borehole_alignment_kwargs(ids, lats, lons, lines, uv)
        line_offsets, default_line = self._structure_line_offsets(ids, lines, uv)
        if not line_offsets:
            return []

        return structure_scene_traces(
            self._pcgs_document.model,
            line_offsets=line_offsets,
            default_line=default_line,
            surface=surface,
        )

    def _on_render(self) -> None:
        if self._view is None:
            return
        try:
            view, topo_note = self._topography_view()
            dmax = self._depth_max()
            opts = self.render_options()
            # exaggerate once every overlay is in (they widen the scene)
            opts.pop("vertical_exaggeration", None)
            fig = view.map3d(**opts)
            for trace in self._borehole_traces():
                fig.add_trace(trace)
            for trace in self._structure_traces():
                fig.add_trace(trace)
            ve = self._style_scene(fig)
            self._plotly_view.show_figure(
                fig, post_script=_SPIN_JS % (
                    "true" if self._btn_spin.isChecked() else "false"))
            if self.isVisible():
                self._cam_timer.start()
            self._topo_lbl.setText(topo_note)
            self._scene_lbl.setText(
                f"{self._combo_mode.currentText()} · "
                f"{self._combo_cmap.currentText()} · ×{ve:.1f}")
            self._status_lbl.setText(
                f"{self._mode()} — "
                f"{self._view.n_stations} station(s) · vertical ×{ve:.1f}"
                f" · {self._combo_cmap.currentText()}"
                + (f" · to {dmax:,.0f} m" if dmax > 0 else ""))
        except Exception as exc:
            self._status_lbl.setText(f"Render error: {exc}")

    # ── toolbar helpers ───────────────────────────────────────────────
    def _rerender(self) -> None:
        if self._view is not None:
            self._on_render()

    def _on_mode_changed(self, _i: int) -> None:
        self._rerender()

    def _mode(self) -> str:
        return self._combo_mode.currentData() or _MODES[0]

    def _on_depth_preset(self, _i: int) -> None:
        value = self._combo_depth.currentData()
        if value is not None and value < 0:  # Custom: the panel's fields
            self._btn_panel.setChecked(True)
            self._sections["Depth & units"].setChecked(True)
            self._spin_depth.setFocus()
            return
        self._spin_depth.blockSignals(True)
        self._spin_depth.setValue(float(value or 0.0))
        self._spin_depth.blockSignals(False)
        self._rerender()

    def _depth_max(self) -> float:
        return float(self._spin_depth.value())

    def _step_ve(self, direction: int) -> None:
        ve = self._spin_ve.value() or 1.0
        steps = [1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0]
        if direction > 0:
            ve = next((s for s in steps if s > ve), steps[-1])
        else:
            ve = next((s for s in reversed(steps) if s < ve), 0.0)
        self._spin_ve.setValue(ve)
        self._rerender()

    def set_camera(self, key: str) -> None:
        """Move the camera to a preset: reset, top, front, side, iso."""
        cam = _CAMERAS.get(key)
        if cam is None:
            return
        cam = {"up": {"x": 0, "y": 0, "z": 1}, "center": {"x": 0, "y": 0,
                                                          "z": 0}, **cam}
        self._camera = cam
        if key == "reset" and self._btn_spin.isChecked():
            self._btn_spin.setChecked(False)
        self._plotly_view.run_js(
            "(function(){var gd=document.querySelector('.js-plotly-plot');"
            "if(gd){Plotly.relayout(gd,{'scene.camera':"
            + json.dumps(cam) + "});}})();")

    def _poll_camera(self) -> None:
        if not self.isVisible() or not self._plotly_view.plot_ready:
            return
        self._plotly_view.eval_js(
            "JSON.stringify(window.__pycsamtCamera || null)",
            self._got_camera)

    def _got_camera(self, value) -> None:
        try:
            cam = json.loads(value) if value else None
        except (TypeError, ValueError):
            cam = None
        if isinstance(cam, dict) and cam.get("eye"):
            self._camera = cam

    def _style_scene(self, fig) -> float:
        """Background, colour bar, legend and vertical exaggeration.

        The Map View themes paint the page ``#eff1f5`` (grey); figures
        for reports need white.  Returns the exaggeration used.
        """
        page, scene, text = _BACKGROUNDS[self._combo_bg.currentText()]
        axis = dict(backgroundcolor=scene, showbackground=True,
                    gridcolor="#d1d5db" if scene != "#0b1220" else "#374151",
                    color=text)
        fig.update_layout(paper_bgcolor=page, plot_bgcolor=page,
                          font_color=text,
                          showlegend=self._chk_legend.isChecked(),
                          scene=dict(bgcolor=scene, xaxis=axis, yaxis=axis,
                                     zaxis=axis))
        show_cb = self._chk_colorbar.isChecked()
        for tr in fig.data:
            if "showscale" in tr and tr.showscale is not None:
                tr.showscale = show_cb
            marker = getattr(tr, "marker", None)
            if marker is not None and getattr(marker, "showscale", None):
                marker.showscale = show_cb
        if self._camera:
            fig.update_layout(scene_camera=self._camera)
        if self._combo_aspect.currentData() != "ve":
            return 1.0  # the chosen Plotly aspect mode rules
        from pycsamt.map.volume import apply_vertical_exaggeration

        return apply_vertical_exaggeration(fig, self._spin_ve.value())

    def _on_spin(self, on: bool) -> None:
        self._plotly_view.run_js(
            f"window.__pycsamtSpin = {'true' if on else 'false'};")

    # ── topography ────────────────────────────────────────────────────
    def _topo_source(self) -> str:
        return self._combo_topo.currentData() or "auto"

    def _on_topo_source(self, _i: int) -> None:
        src = self._topo_source()
        self._btn_topo_file.setVisible(src == "topo_file")
        if src == "topo_file" and not self._topo_elev:
            self._on_load_topo_file()
            return
        self._rerender()

    def _on_load_topo_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load topography", "",
            "Topography (*.csv *.bln *.stn *.txt);;All files (*)")
        if path:
            self.load_topography(path)

    def load_topography(self, source) -> int:
        """Station elevations from a .csv/.bln/.stn file (or any
        :func:`pycsamt.format.topo_source.resolve_topo` source); returns
        the number of matched stations."""
        from pycsamt.format.topo_source import resolve_topo

        if self._view is None:
            return 0
        names = [s.id for s in self._view.data.stations]
        import warnings

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                att = resolve_topo(source, names, on_mismatch="warn")
            elev = dict(att.elevation)
        except Exception as exc:
            # A plain "station, elevation" table has no coordinates for
            # resolve_topo; match it by station name instead.
            elev = _name_elevation_table(source, names)
            if not elev:
                self._topo_lbl.setText(f"Topography error: {exc}")
                return 0
        self._topo_elev = elev
        self._combo_topo.blockSignals(True)
        self._combo_topo.setCurrentIndex(
            self._combo_topo.findData("topo_file"))
        self._combo_topo.blockSignals(False)
        self._btn_topo_file.setVisible(True)
        self._rerender()
        return len(self._topo_elev)

    def _source_elevations(self, src: str) -> tuple[dict, str]:
        """Raw ``({station: elevation}, note)`` of one source."""
        v = self._view
        if src == "file":
            good = {k: e for k, e in self._file_elev.items()
                    if e is not None}
            if not any(abs(float(e)) > 0 for e in good.values()):
                return {}, ("The file has no topography — choose loaded "
                            "stations, a topography file or online.")
            how = self._file_elev_how or "file"
            return good, f"Topography from the file ({how})"
        if src == "stations":
            if self._sites is None:
                return {}, ("No survey loaded in the main window — load "
                            "EDI/XML data there first.")
            from pycsamt.format.topo_source import resolve_topo

            import warnings

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    att = resolve_topo(self._sites,
                                       [s.id for s in v.data.stations],
                                       on_mismatch="warn")
            except Exception as exc:
                return {}, f"Topography from stations failed: {exc}"
            if not att.elevation:
                return {}, "No station name matches the loaded survey."
            return dict(att.elevation), "Topography from the loaded survey"
        if src == "topo_file":
            if not self._topo_elev:
                return {}, "Load a topography file."
            return dict(self._topo_elev), "Topography from file"
        if src == "online":
            if not self._online_elev:
                try:
                    self._online_elev = v.fetch_elevations()
                except Exception as exc:
                    return {}, f"Online elevation lookup failed: {exc}"
            return dict(self._online_elev), "Topography fetched online"
        return {}, ""

    def _topography_view(self):
        """(view with the chosen elevations, status note)."""
        v = self._view
        n = v.n_stations
        src = self._topo_source()
        if src == "flat":
            self._elev_used = {s.id: 0.0 for s in v.data.stations}
            return v.with_elevations(self._elev_used), \
                "Flat surface (elevation 0)."
        if src == "auto":
            elev, note = self._source_elevations("file")
            if not elev and self._sites is not None:
                elev, note = self._source_elevations("stations")
            if not elev:
                self._elev_used = {}
                return v, ("No topography: the file has none and no survey "
                           "is loaded — pick a topography file or online.")
        else:
            elev, note = self._source_elevations(src)
            if not elev:
                self._elev_used = {}
                return v, note
        matched = sum(1 for s in v.data.stations if s.id in elev)
        filled = 0
        if self._chk_fill_gaps.isChecked() and matched < n:
            elev, filled = fill_line_gaps(v.data.stations, elev)
        self._elev_used = elev
        extra = f", {filled} filled along lines" if filled else ""
        return v.with_elevations(elev), (
            f"{note} ({matched}/{n} stations{extra}).")

    def set_sites(self, sites) -> None:
        """The main window's survey (used for "Loaded stations" topo)."""
        super().set_sites(sites)
        if self._view is not None and self._topo_source() in ("stations",
                                                              "auto"):
            self._on_render()

    # ── output ────────────────────────────────────────────────────────
    def mapview_state(self) -> dict:
        """The scene as Map View should open it (see ``pcsf_scene``)."""
        render = self.render_options()
        render.pop("geology", None)  # bands travel as the legend store
        render.pop("geology_patterns", None)
        render.setdefault("geology_fill", "solid")
        return mapview_state(
            render,
            elevations=self._elev_used,
            spin=self._btn_spin.isChecked(),
            dark=self._combo_bg.currentText() == "Dark",
            camera=self._camera,
            geology=self._legend if self._chk_show_geology.isChecked()
            else None,
            boreholes=self._pcbh_document
            if self._chk_show_boreholes.isChecked() else None,
            structure=self._pcgs_document
            if self._chk_show_structure.isChecked() else None,
            source=str(self._path or ""),
        )

    def _on_open_mapview(self) -> None:
        if getattr(self, "_path", None) is None:
            return
        from pycsamt.app.desktop import mapview_bridge

        try:
            state = self.mapview_state()
        except Exception as exc:  # hand over the model at least
            state = None
            self._status_lbl.setText(f"Scene not handed over: {exc}")
        try:
            launch = mapview_bridge.launch_mapview(str(self._path),
                                                   state=state)
        except Exception as exc:
            self._status_lbl.setText(f"Map View could not start: {exc}")
            return
        self._status_lbl.setText(
            f"Opening Map View at {launch.url} with this scene (starts in "
            "a few seconds).")

    def _on_png(self) -> None:
        if self._view is None or not self._plotly_view.plot_ready:
            self._status_lbl.setText("Render a scene first.")
            return
        default = (self._path.with_suffix(".png").name if self._path
                   else "scene.png")
        path, _ = QFileDialog.getSaveFileName(self, "Save scene as PNG",
                                              default, "PNG (*.png)")
        if not path:
            return
        self._status_lbl.setText("Saving PNG…")
        self._plotly_view.snapshot_png(
            path, lambda p, err: self._status_lbl.setText(
                f"Saved {p}" if p else f"PNG failed: {err}"))

    def _on_export(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Export scene", "", "HTML files (*.html)")
        if not path:
            return
        try:
            self._plotly_view.export_html(path)
            self._status_lbl.setText(f"Exported to {path}")
        except Exception as exc:
            self._status_lbl.setText(f"Export error: {exc}")

    def hideEvent(self, event) -> None:  # noqa: N802
        self._cam_timer.stop()
        super().hideEvent(event)


def _name_elevation_table(source, names) -> dict:
    """{station: elevation} from a CSV/TXT with a name and an elevation
    column (matched case-insensitively to *names*); ``{}`` otherwise."""
    import pandas as pd

    try:
        df = pd.read_csv(source, sep=None, engine="python")
    except Exception:
        return {}
    cols = {c.lower().strip(): c for c in df.columns}
    name_col = next((cols[c] for c in ("station", "name", "site", "id",
                                       "station_id") if c in cols), None)
    elev_col = next((cols[c] for c in ("elevation", "elev", "z", "altitude",
                                       "height", "topo") if c in cols), None)
    if name_col is None or elev_col is None:
        return {}
    wanted = {str(n).lower(): n for n in names}
    out = {}
    for n, e in zip(df[name_col].astype(str), df[elev_col]):
        key = wanted.get(n.strip().lower())
        try:
            e = float(e)
        except (TypeError, ValueError):
            continue
        if key is not None and e == e:
            out[key] = e
    return out
