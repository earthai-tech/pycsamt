# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ProfilePanel — tabbed scientific plot panel.

Six tabs driven by PlotController:

  [ρₐ / φ]           Single-station apparent-resistivity + phase curves
  [Pseudosection ρₐ]  ρₐ pseudosection (full profile)
  [Pseudosection φ]   φ pseudosection (full profile)
  [Tipper]            Tipper components
  [Phase Tensor]      Phase-tensor pseudosection (Caldwell 2004 style)
  [2D Section]        SectionPanel — inversion result viewer

Public API:
    set_sites(sites)
    set_selected_station(sid)
    set_period_range(lo_hz, hi_hz)
    set_dark_mode(bool)
    set_inversion_result(result)
"""

from __future__ import annotations

import numpy as np

from PySide6.QtWidgets import (
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.plot_controller import PlotController
from pycsamt.app.desktop.panels.section_panel import SectionPanel
from pycsamt.app.desktop.widgets.freq_selector import FreqSelector
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas


class ProfilePanel(QWidget):
    """Tabbed panel for MT response curves and pseudosections."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._ctrl = PlotController()
        self._build_ui()

    # ── Construction ──────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # Frequency selector bar
        self._freq_sel = FreqSelector(parent=self)
        self._freq_sel.setFixedHeight(56)
        self._freq_sel.range_changed.connect(self._on_freq_range_changed)
        root.addWidget(self._freq_sel)

        # Tab widget
        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        root.addWidget(self._tabs)

        # Tab 0: ρₐ / φ response curves
        self._canvas_rho_phi = MplCanvas(self, toolbar=True)
        self._tabs.addTab(self._canvas_rho_phi, "ρₐ / φ")

        # Tab 1: Apparent-resistivity pseudosection
        self._canvas_rho_ps = MplCanvas(self, toolbar=True)
        self._tabs.addTab(self._canvas_rho_ps, "Pseudosection ρₐ")

        # Tab 2: Phase pseudosection
        self._canvas_ph_ps = MplCanvas(self, toolbar=True)
        self._tabs.addTab(self._canvas_ph_ps, "Pseudosection φ")

        # Tab 3: Tipper
        self._canvas_tipper = MplCanvas(self, toolbar=True)
        self._tabs.addTab(self._canvas_tipper, "Tipper")

        # Tab 4: Phase tensor
        self._canvas_pt = MplCanvas(self, toolbar=True)
        self._tabs.addTab(self._canvas_pt, "Phase Tensor")

        # Tab 5: 2D section — SectionPanel
        self._section_panel = SectionPanel(self)
        self._tabs.addTab(self._section_panel, "2D Section")

        # Lazy redraw on tab switch
        self._tabs.currentChanged.connect(self._on_tab_changed)

        self._draw_empty_all()

    # ── Public API ─────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        """Load a Sites collection and redraw all tabs."""
        self._ctrl.set_sites(sites)
        try:
            freqs = []
            for site in sites:
                f = getattr(site, "freq", None)
                if f is not None:
                    freqs.extend(np.asarray(f).ravel().tolist())
            if freqs:
                f_min = max(float(min(freqs)), 1e-6)
                f_max = float(max(freqs))
                self._freq_sel.set_freq_range(f_min, f_max)
        except Exception:
            pass
        self._redraw_all()

    def set_selected_station(self, station_id: str) -> None:
        """Highlight a station; redraw ρₐ/φ tab and mark active pseudosection."""
        self._ctrl.set_station(station_id)
        self._redraw_rho_phi()
        tab = self._tabs.currentIndex()
        if tab == 1:
            self._redraw_rho_pseudosection()
        elif tab == 2:
            self._redraw_phase_pseudosection()
        elif tab == 3:
            self._redraw_tipper()

    def set_dark_mode(self, dark: bool) -> None:
        self._ctrl.dark = dark
        self._section_panel.set_dark_mode(dark)
        self._redraw_current_tab()

    def set_inversion_result(self, result) -> None:
        """Load an InversionResult into the 2D Section tab."""
        self._section_panel.set_result(result)
        self._tabs.setCurrentWidget(self._section_panel)

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_freq_range_changed(self, lo_hz: float, hi_hz: float) -> None:
        T_max = 1.0 / lo_hz if lo_hz > 0 else None
        T_min = 1.0 / hi_hz if hi_hz > 0 else None
        self._ctrl.set_period_range(T_min, T_max)
        self._redraw_current_tab()

    def _on_tab_changed(self, index: int) -> None:
        """Lazy-redraw: only draw the tab when it becomes visible."""
        _redraw = [
            self._redraw_rho_phi,
            self._redraw_rho_pseudosection,
            self._redraw_phase_pseudosection,
            self._redraw_tipper,
            self._redraw_phase_tensor,
        ]
        if index < len(_redraw):
            _redraw[index]()

    # ── Internal draw helpers ──────────────────────────────────────────

    def _draw_empty_all(self) -> None:
        from pycsamt.app.desktop.controllers.plot_controller import (
            style_axes, _annotate_empty,
        )
        for canvas, msg in (
            (self._canvas_rho_phi, "Select a station to view ρₐ / φ curves"),
            (self._canvas_rho_ps,  "Load survey data"),
            (self._canvas_ph_ps,   "Load survey data"),
            (self._canvas_tipper,  "Load survey data"),
            (self._canvas_pt,      "Load survey data"),
        ):
            canvas.figure.clear()
            ax = canvas.figure.add_subplot(111)
            _annotate_empty(ax, msg)
            style_axes(ax, self._ctrl.dark)
            canvas.draw()

    def _redraw_all(self) -> None:
        self._redraw_rho_phi()
        self._redraw_rho_pseudosection()
        self._redraw_phase_pseudosection()
        self._redraw_tipper()
        self._redraw_phase_tensor()

    def _redraw_current_tab(self) -> None:
        tab = self._tabs.currentIndex()
        _redraw = [
            self._redraw_rho_phi,
            self._redraw_rho_pseudosection,
            self._redraw_phase_pseudosection,
            self._redraw_tipper,
            self._redraw_phase_tensor,
        ]
        if tab < len(_redraw):
            _redraw[tab]()

    def _redraw_rho_phi(self) -> None:
        fig = self._canvas_rho_phi.figure
        self._ctrl.draw_rho_phi(fig)
        self._canvas_rho_phi.draw()

    def _redraw_rho_pseudosection(self) -> None:
        fig = self._canvas_rho_ps.figure
        fig.clear()
        ax = fig.add_subplot(111)
        self._canvas_rho_ps.axes = ax
        self._ctrl.draw_rho_pseudosection(ax)
        self._canvas_rho_ps.draw()

    def _redraw_phase_pseudosection(self) -> None:
        fig = self._canvas_ph_ps.figure
        fig.clear()
        ax = fig.add_subplot(111)
        self._canvas_ph_ps.axes = ax
        self._ctrl.draw_phase_pseudosection(ax)
        self._canvas_ph_ps.draw()

    def _redraw_tipper(self) -> None:
        fig = self._canvas_tipper.figure
        fig.clear()
        ax = fig.add_subplot(111)
        self._canvas_tipper.axes = ax
        self._ctrl.draw_tipper(ax)
        self._canvas_tipper.draw()

    def _redraw_phase_tensor(self) -> None:
        fig = self._canvas_pt.figure
        fig.clear()
        ax = fig.add_subplot(111)
        self._canvas_pt.axes = ax
        self._ctrl.draw_phase_tensor(ax)
        self._canvas_pt.draw()
