# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SectionPanel — 2-D resistivity section viewer (Phase 5).

Displays an Occam2D (or forward-model) 2-D section as a pcolormesh
with a log-resistivity colourbar.  Supports:

  * Loading an InversionResult from a working directory
  * Parameter controls (ρ range, depth clip, station markers, colormap)
  * Optional comparison overlay of a second result (dashed contours)
  * Export via ExportDialog

Wired as the content of Profile tab "2D Section" inside ProfilePanel.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.colors import LogNorm
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QPushButton,
    QSplitter,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.plot_controller import (
    style_axes,
)
from pycsamt.app.desktop.widgets.colorbar_widget import (
    ColorbarWidget,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

_COLORMAPS = ["jet_r", "RdBu_r", "viridis_r", "plasma_r", "rainbow_r", "bwr_r"]


class SectionPanel(QWidget):
    """
    2-D resistivity section viewer.

    Signals
    -------
    result_loaded : str
        Emitted with the workdir path after a result is loaded.
    """

    result_loaded = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._result        = None   # InversionResult (primary)
        self._result_ref    = None   # InversionResult (comparison overlay)
        self._dark: bool    = True
        self._build_ui()

    # ── Construction ──────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── Top button bar ─────────────────────────────────────────
        btn_bar = QHBoxLayout()
        btn_bar.setContentsMargins(4, 4, 4, 2)

        self._btn_load = QPushButton("Load Result…")
        self._btn_load.setStatusTip("Load an Occam2D inversion working directory")
        self._btn_load.clicked.connect(self._on_load)
        btn_bar.addWidget(self._btn_load)

        self._btn_compare = QPushButton("+ Compare…")
        self._btn_compare.setStatusTip("Overlay a second inversion result for comparison")
        self._btn_compare.setEnabled(False)
        self._btn_compare.clicked.connect(self._on_load_ref)
        btn_bar.addWidget(self._btn_compare)

        self._btn_clear_compare = QPushButton("Clear Compare")
        self._btn_clear_compare.setEnabled(False)
        self._btn_clear_compare.clicked.connect(self._on_clear_ref)
        btn_bar.addWidget(self._btn_clear_compare)

        btn_bar.addStretch()

        self._btn_export = QPushButton("⬆ Export…")
        self._btn_export.setEnabled(False)
        self._btn_export.clicked.connect(self._on_export)
        btn_bar.addWidget(self._btn_export)

        root.addLayout(btn_bar)

        # ── Main content: controls | canvas + colorbar ─────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(3)

        # Left: controls panel
        ctrl_widget = self._build_controls()
        ctrl_widget.setFixedWidth(200)
        splitter.addWidget(ctrl_widget)

        # Right: canvas + colorbar
        right = QWidget()
        right_layout = QHBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(0)

        self._canvas = MplCanvas(self, toolbar=True)
        right_layout.addWidget(self._canvas)

        self._cbar = ColorbarWidget(orientation="vertical", parent=self)
        self._cbar.setFixedWidth(60)
        right_layout.addWidget(self._cbar)

        splitter.addWidget(right)
        splitter.setStretchFactor(1, 1)
        root.addWidget(splitter)

        # ── Bottom: summary text ───────────────────────────────────
        self._summary = QTextBrowser()
        self._summary.setMaximumHeight(80)
        self._summary.setPlaceholderText("Load an inversion result to view the section.")
        root.addWidget(self._summary)

        # Draw initial empty state
        self._draw_empty()

    def _build_controls(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)

        # Resistivity range
        rho_box = QGroupBox("Resistivity range (Ω·m)")
        rho_form = QFormLayout(rho_box)

        self._rho_min = QDoubleSpinBox()
        self._rho_min.setRange(0.01, 1e6)
        self._rho_min.setDecimals(1)
        self._rho_min.setValue(1.0)
        self._rho_min.setSuffix(" Ω·m")
        rho_form.addRow("Min:", self._rho_min)

        self._rho_max = QDoubleSpinBox()
        self._rho_max.setRange(0.01, 1e8)
        self._rho_max.setDecimals(0)
        self._rho_max.setValue(1000.0)
        self._rho_max.setSuffix(" Ω·m")
        rho_form.addRow("Max:", self._rho_max)

        layout.addWidget(rho_box)

        # Display options
        disp_box = QGroupBox("Display")
        disp_form = QFormLayout(disp_box)

        self._depth_max = QDoubleSpinBox()
        self._depth_max.setRange(0.0, 100_000.0)
        self._depth_max.setValue(0.0)
        self._depth_max.setSpecialValueText("Auto")
        self._depth_max.setSuffix(" m")
        disp_form.addRow("Depth max:", self._depth_max)

        self._show_stations = QCheckBox()
        self._show_stations.setChecked(True)
        disp_form.addRow("Show stations:", self._show_stations)

        self._cmap_combo = QComboBox()
        self._cmap_combo.addItems(_COLORMAPS)
        disp_form.addRow("Colormap:", self._cmap_combo)

        layout.addWidget(disp_box)

        # Topography group
        topo_box = QGroupBox("Topography")
        topo_form = QFormLayout(topo_box)
        self._chk_topo = QCheckBox()
        self._chk_topo.setChecked(False)
        self._chk_topo.setToolTip(
            "Drape the section over real terrain elevation.\n"
            "Reads PYCSAMT_TOPO global setting; toggling this checkbox\n"
            "updates that setting for the whole session."
        )
        self._chk_topo.toggled.connect(self._on_topo_toggled)
        topo_form.addRow("Include topo:", self._chk_topo)

        self._spin_exag = QDoubleSpinBox()
        self._spin_exag.setRange(0.1, 20.0)
        self._spin_exag.setSingleStep(0.5)
        self._spin_exag.setValue(1.0)
        self._spin_exag.setDecimals(1)
        self._spin_exag.setSuffix(" ×")
        self._spin_exag.setEnabled(False)
        self._spin_exag.setToolTip("Vertical exaggeration of the terrain surface")
        self._spin_exag.valueChanged.connect(self._on_exag_changed)
        topo_form.addRow("Exaggeration:", self._spin_exag)

        layout.addWidget(topo_box)

        # Redraw button
        btn_redraw = QPushButton("↻ Redraw")
        btn_redraw.clicked.connect(self._redraw)
        layout.addWidget(btn_redraw)

        layout.addStretch()
        return w

    # ── Public API ─────────────────────────────────────────────────────

    def set_result(self, result) -> None:
        """Load an InversionResult programmatically (from wizard/worker)."""
        self._result = result
        self._btn_compare.setEnabled(True)
        self._btn_export.setEnabled(True)
        self._show_summary()
        self._redraw()
        try:
            self.result_loaded.emit(str(getattr(result, "workdir", "")))
        except Exception:
            pass

    def set_dark_mode(self, dark: bool) -> None:
        self._dark = dark
        if self._result is not None:
            self._redraw()

    # ── Topo slots ────────────────────────────────────────────────────

    def _on_topo_toggled(self, checked: bool) -> None:
        from pycsamt.topo import configure_topo
        configure_topo(enabled=checked)
        self._spin_exag.setEnabled(checked)
        if self._result is not None:
            self._redraw()

    def _on_exag_changed(self, value: float) -> None:
        from pycsamt.topo import configure_topo
        configure_topo(exaggeration=value)
        if self._result is not None and self._chk_topo.isChecked():
            self._redraw()

    def clear(self) -> None:
        self._result     = None
        self._result_ref = None
        self._btn_compare.setEnabled(False)
        self._btn_export.setEnabled(False)
        self._btn_clear_compare.setEnabled(False)
        self._summary.clear()
        self._draw_empty()

    # ── Load from disk ─────────────────────────────────────────────────

    def _on_load(self) -> None:
        workdir = QFileDialog.getExistingDirectory(
            self, "Select Occam2D working directory", str(Path.home())
        )
        if not workdir:
            return
        self._load_from_workdir(workdir, is_reference=False)

    def _on_load_ref(self) -> None:
        workdir = QFileDialog.getExistingDirectory(
            self, "Select reference Occam2D working directory", str(Path.home())
        )
        if not workdir:
            return
        self._load_from_workdir(workdir, is_reference=True)

    def _load_from_workdir(self, workdir: str, is_reference: bool) -> None:
        try:
            from pycsamt.models.occam2d import InversionResult
            result = InversionResult(workdir=workdir)
        except Exception as exc:
            self._summary.setPlainText(f"Failed to load '{workdir}':\n{exc}")
            return

        if is_reference:
            self._result_ref = result
            self._btn_clear_compare.setEnabled(True)
        else:
            self.set_result(result)

        self._redraw()

    def _on_clear_ref(self) -> None:
        self._result_ref = None
        self._btn_clear_compare.setEnabled(False)
        self._redraw()

    # ── Drawing ────────────────────────────────────────────────────────

    def _draw_empty(self) -> None:
        ax = self._canvas.axes
        ax.cla()
        ax.text(
            0.5, 0.5,
            "Load an inversion result\nor run an inversion from\n"
            "Processing → Inversion Wizard",
            transform=ax.transAxes,
            ha="center", va="center",
            fontsize=11, color="#585b70",
        )
        style_axes(ax, self._dark)
        self._canvas.draw()

    def _redraw(self) -> None:
        if self._result is None:
            self._draw_empty()
            return

        ax = self._canvas.axes
        ax.cla()

        try:
            self._draw_section(ax, self._result, alpha=1.0)
            if self._result_ref is not None:
                self._draw_section_contours(ax, self._result_ref)
        except Exception as exc:
            ax.cla()
            ax.text(
                0.5, 0.5, f"Render error:\n{exc}",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10, color="#f38ba8",
            )

        style_axes(ax, self._dark)
        self._canvas.draw()

        # Update colorbar
        self._cbar.update_colorbar(
            cmap=self._cmap_combo.currentText(),
            vmin=self._rho_min.value(),
            vmax=self._rho_max.value(),
            label="Resistivity (Ω·m)",
        )

    def _draw_section(self, ax, result, alpha: float = 1.0) -> None:
        """Draw pcolormesh 2-D resistivity section on *ax*."""
        rho_log10 = np.asarray(result.rho_2d)   # (nz, nx) in log10(Ω·m)
        mesh      = result.mesh
        xn = np.asarray(mesh.x_nodes)
        zn = np.asarray(mesh.z_nodes)

        rho_lin = np.where(np.isfinite(rho_log10), 10.0 ** rho_log10, np.nan)

        rho_min = self._rho_min.value()
        rho_max = self._rho_max.value()
        cmap    = self._cmap_combo.currentText()

        # Clip depth if requested
        depth_max = self._depth_max.value()
        if depth_max > 0:
            zn       = zn[: np.searchsorted(zn, depth_max) + 1]
            rho_lin  = rho_lin[: len(zn) - 1, :]

        norm = LogNorm(vmin=rho_min, vmax=rho_max, clip=True)

        # Centre x on profile midpoint for readability
        x_shift = float((xn.min() + xn.max()) / 2.0)
        xn_c    = (xn - x_shift) / 1000.0   # m → km
        zn_km   = zn / 1000.0

        # ── Topography draped rendering ────────────────────────────────
        from pycsamt.topo.config import PYCSAMT_TOPO
        use_topo = PYCSAMT_TOPO.enabled and self._chk_topo.isChecked()

        if use_topo:
            self._draw_section_topo(
                ax, xn_c, zn_km, rho_lin, norm, cmap, alpha, result, x_shift
            )
        else:
            # Standard flat-datum rendering
            ax.pcolormesh(
                xn_c, zn_km, rho_lin,
                norm=norm, cmap=cmap, alpha=alpha, shading="flat",
            )
            ax.set_ylim(zn_km.max(), 0)   # depth positive downward
            ax.set_ylabel("Depth (km)", fontsize=9)

            if self._show_stations.isChecked():
                try:
                    sx = (np.asarray(result.mesh.station_x) - x_shift) / 1000.0
                    ax.plot(sx, np.zeros_like(sx),
                            "v", color="#f5c542", ms=6, zorder=5)
                except Exception:
                    pass

        ax.set_xlabel("Distance (km)", fontsize=9)

        # RMS annotation
        try:
            rms_txt = f"RMS = {result.final_rms:.3f}  ({result.n_iterations} iter)"
            ax.set_title(rms_txt, fontsize=9, pad=3)
        except Exception:
            pass

    def _draw_section_topo(
        self, ax, xn_c, zn_km, rho_lin, norm, cmap, alpha, result, x_shift
    ) -> None:
        """Terrain-following variant of _draw_section."""
        from pycsamt.topo.config import PYCSAMT_TOPO
        from pycsamt.topo.drape import (
            drape_section,
            interp_elev,
        )
        from pycsamt.topo.overlay import draw_topo_section

        # Try to get station elevations from the result object
        try:
            sx_m      = np.asarray(result.mesh.station_x)
            sx_km     = (sx_m - x_shift) / 1000.0
            elev_m    = np.asarray(getattr(result.mesh, "station_elev",
                                           np.zeros(len(sx_m))))
            elev_km   = elev_m / 1000.0
        except Exception:
            # No elevation data → fall back to flat rendering
            ax.pcolormesh(xn_c, zn_km, rho_lin,
                          norm=norm, cmap=cmap, alpha=alpha, shading="flat")
            ax.set_ylim(zn_km.max(), 0)
            ax.set_ylabel("Depth (km)", fontsize=9)
            return

        # Cell-centre x positions for elevation interpolation
        x_centres = (xn_c[:-1] + xn_c[1:]) / 2.0
        elev_at_centres = interp_elev(sx_km, elev_km, x_centres,
                                      PYCSAMT_TOPO.interp_method)

        # Build terrain-following z grid
        _, z_draped, rho_draped = drape_section(
            xn_c, zn_km, rho_lin, elev_at_centres,
            exaggeration=PYCSAMT_TOPO.exaggeration,
            clip_above_surface=PYCSAMT_TOPO.clip_below_surface,
        )

        # pcolormesh with 2-D Y grid (requires Matplotlib >= 3.3)
        ax.pcolormesh(
            xn_c, z_draped, rho_draped,
            norm=norm, cmap=cmap, alpha=alpha, shading="flat",
        )

        # Y-axis: absolute elevation (surface at top)
        y_surface = np.max(elev_at_centres)
        y_deep    = y_surface - zn_km.max() * PYCSAMT_TOPO.exaggeration
        ax.set_ylim(y_deep, y_surface + 0.05 * abs(y_surface - y_deep))
        ax.set_ylabel("Elevation (km a.s.l.)", fontsize=9)

        # Terrain overlay
        if self._show_stations.isChecked():
            try:
                snames = [str(getattr(result.mesh, "station_names", [f"S{i}" for i in range(len(sx_km))])
                              [i] if hasattr(result.mesh, "station_names") else f"S{i}")
                          for i in range(len(sx_km))]
            except Exception:
                snames = None
            draw_topo_section(
                ax, sx_km, elev_m, snames,
                station_x_km=sx_km,
                dark=self._dark,
            )

    def _draw_section_contours(self, ax, result_ref) -> None:
        """Overlay dashed contour lines from the reference result."""
        try:
            rho_log10 = np.asarray(result_ref.rho_2d)
            mesh      = result_ref.mesh
            xn = np.asarray(mesh.x_nodes)
            zn = np.asarray(mesh.z_nodes)
            x_shift = float((xn.min() + xn.max()) / 2.0)
            xc = ((xn[:-1] + xn[1:]) / 2.0 - x_shift) / 1000.0
            zc = ((zn[:-1] + zn[1:]) / 2.0) / 1000.0
            ax.contour(
                xc, zc, rho_log10,
                levels=8,
                colors="white",
                linewidths=0.6,
                linestyles="--",
                alpha=0.55,
            )
        except Exception:
            pass

    # ── Summary ────────────────────────────────────────────────────────

    def _show_summary(self) -> None:
        if self._result is None:
            return
        try:
            self._summary.setPlainText(self._result.summary)
        except Exception:
            try:
                lines = [
                    f"Workdir:    {getattr(self._result, 'workdir', '?')}",
                    f"Iterations: {getattr(self._result, 'n_iterations', '?')}",
                    f"Final RMS:  {getattr(self._result, 'final_rms', '?')}",
                ]
                self._summary.setPlainText("\n".join(lines))
            except Exception:
                pass

    # ── Export ────────────────────────────────────────────────────────

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import (
            ExportDialog,
        )
        ExportDialog(figure=self._canvas.figure, parent=self).exec()
