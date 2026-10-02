# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ProfilePanel — tabbed scientific plot panel.

Core tabs driven by PlotController, with Tipper added when available, plus
one self-contained tab that manages its own controls:

  [ρₐ / φ]           Single-station apparent-resistivity + phase curves
  [Pseudosection ρₐ]  ρₐ pseudosection (full profile)
  [Pseudosection φ]   φ pseudosection (full profile)
  [Res/Phase Section] MTPy-style stacked ρₐ/φ pseudo-section, one column
                       per component (PlotResPhasePseudoSection)
  [Tipper]            Tipper components
  [Phase Tensor]      Phase-tensor pseudosection (Caldwell 2004 style)
  [PT Strip]          Single-station phase-tensor ellipse strip vs period

Every tab here reads from raw survey soundings (``set_sites`` / EDI data).
A 2-D *inversion-result* section is a different kind of data (a loaded
``InversionResult``, not a Sites collection) and lives in
``InversionWindow`` instead, which is where a result is actually produced
and where every other result-viewing tab (Model/Fit/Convergence) already
is — see that window's own "2-D Section" tab (``SectionPanel``). It was
briefly wired in here during Phase 1 on the strength of a stale test/
docstring that anticipated it on this panel; moved once the mismatch was
caught (an inversion result has no natural link to ``set_sites``' survey
soundings at all).

Public API:
    set_sites(sites)
    set_selected_station(sid)
    set_period_range(lo_hz, hi_hz)
    set_dark_mode(bool)
"""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import (
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.plot_controller import (
    PlotController,
)
from pycsamt.app.desktop.panels.resphase_section_panel import (
    ResPhaseSectionPanel,
)
from pycsamt.app.desktop.widgets.freq_selector import (
    FreqSelector,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas


class ProfilePanel(QWidget):
    """Tabbed panel for MT response curves and pseudosections."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._ctrl = PlotController()
        self._dirty_canvases: set[MplCanvas] = set()
        # Phase-tensor tab cache key: tuple returned by PlotController.phase_tensor_key()
        # stored after the last successful draw.  When the current key matches,
        # skip the full matplotlib redraw and just repaint the existing canvas.
        self._pt_last_key: tuple | None = None
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

        # Optional Tipper tab (inserted at index 3 when data is available)
        self._canvas_tipper = MplCanvas(self, toolbar=True)
        self._canvas_tipper.hide()

        # Phase tensor
        self._canvas_pt = MplCanvas(self, toolbar=True)
        self._tabs.addTab(self._canvas_pt, "Phase Tensor")

        # Phase-tensor ellipse strip (single selected station)
        self._canvas_pt_strip = MplCanvas(self, toolbar=True)
        self._tabs.addTab(self._canvas_pt_strip, "PT Strip")

        # MTPy-style stacked res/phase pseudo-section (self-contained: owns
        # its own controls bar, so it is not part of the MplCanvas-only
        # dirty-tracking/lazy-redraw scheme below).
        self._resphase_section = ResPhaseSectionPanel(self)
        self._tabs.addTab(self._resphase_section, "Res/Phase Section")

        # Lazy redraw on tab switch
        self._tabs.currentChanged.connect(self._on_tab_changed)

        self._draw_empty_all()
        self._dirty_canvases.clear()

    # ── Public API ─────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        """Load a Sites collection and redraw all tabs."""
        self._pt_last_key = None  # new data always means a fresh PT draw
        self._ctrl.set_sites(sites)
        # A global/profile line change may remove the formerly selected
        # station. Never leave a stale ID that resolves to None in a
        # station-specific plot such as Tipper or PT Strip.
        try:
            names = [str(site.name) for site in sites]
        except Exception:
            names = []
        if self._ctrl._station_id not in names:
            self._ctrl.set_station(names[0] if names else None)
        self._set_tipper_tab_visible(self._sites_have_tipper(sites))
        try:
            freqs = []
            for site in sites:
                f = getattr(site, "freq", None)
                if f is not None:
                    freqs.extend(np.asarray(f).ravel().tolist())
            if freqs:
                f_min = max(float(min(freqs)), 1e-6)
                f_max = float(max(freqs))
                self._freq_sel.blockSignals(True)
                self._freq_sel.set_freq_range(f_min, f_max)
                self._freq_sel.blockSignals(False)
        except Exception:
            self._freq_sel.blockSignals(False)
            pass
        self._mark_all_dirty()
        try:
            self._redraw_current_tab(force=True)
        except Exception:
            pass  # a plot failure must not block line/station controls
        try:
            self._resphase_section.set_sites(sites)
        except Exception:
            pass  # a plot failure in one tab must not block the others

    def set_selected_station(self, station_id: str) -> None:
        """Highlight a station; redraw ρₐ/φ tab and mark active pseudosection."""
        self._ctrl.set_station(station_id)
        self._mark_all_dirty()
        self._redraw_current_tab(force=True)

    def set_dark_mode(self, dark: bool) -> None:
        self._ctrl.dark = dark
        # Dark mode changes the PT plot styling → force a full redraw next time
        self._pt_last_key = None
        self._mark_all_dirty()
        self._redraw_current_tab()
        try:
            self._resphase_section.set_dark_mode(dark)
        except Exception:
            pass

    def invalidate_phase_tensor(self) -> None:
        """Force the Phase Tensor tab to recompute and redraw on the next visit.

        Call this before an explicit user-triggered refresh so the DataFrame
        cache is also cleared (useful when the user changes settings and wants
        a fresh computation, not just a repaint of cached data).
        """
        self._pt_last_key = None
        self._ctrl.invalidate_phase_tensor()

    # ── Slots ─────────────────────────────────────────────────────────

    def _on_freq_range_changed(self, lo_hz: float, hi_hz: float) -> None:
        T_max = 1.0 / lo_hz if lo_hz > 0 else None
        T_min = 1.0 / hi_hz if hi_hz > 0 else None
        self._ctrl.set_period_range(T_min, T_max)
        self._mark_all_dirty()
        self._redraw_current_tab()

    def _on_tab_changed(self, index: int) -> None:
        """Lazy-redraw: only draw the tab when it becomes visible."""
        self._redraw_current_tab(force=False)
        current = self._tabs.widget(index)
        self._resphase_section.set_active(current is self._resphase_section)

    @staticmethod
    def _sites_have_tipper(sites) -> bool:
        try:
            for site in sites or ():
                tipper = getattr(site, "tipper", None)
                if tipper is not None and np.asarray(tipper).size:
                    return True
                if bool(site.summary().get("tipper", False)):
                    return True
        except Exception:
            return False
        return False

    def _set_tipper_tab_visible(self, visible: bool) -> None:
        index = self._tabs.indexOf(self._canvas_tipper)
        if visible and index < 0:
            self._canvas_tipper.show()
            self._tabs.insertTab(3, self._canvas_tipper, "Tipper")
        elif not visible and index >= 0:
            self._tabs.removeTab(index)
            self._canvas_tipper.hide()

    # ── Internal draw helpers ──────────────────────────────────────────

    def _draw_empty_all(self) -> None:
        from pycsamt.app.desktop.controllers.plot_controller import (
            _annotate_empty,
            style_axes,
        )

        for canvas, msg in (
            (self._canvas_rho_phi, "Select a station to view ρₐ / φ curves"),
            (self._canvas_rho_ps, "Load survey data"),
            (self._canvas_ph_ps, "Load survey data"),
            (self._canvas_tipper, "Load survey data"),
            (self._canvas_pt, "Load survey data"),
            (
                self._canvas_pt_strip,
                "Select a station to view its ellipse strip",
            ),
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
        self._redraw_phase_tensor_strip()
        self._dirty_canvases.clear()

    def _mark_all_dirty(self) -> None:
        self._dirty_canvases = {
            self._canvas_rho_phi,
            self._canvas_rho_ps,
            self._canvas_ph_ps,
            self._canvas_tipper,
            self._canvas_pt,
            self._canvas_pt_strip,
        }

    def _redraw_current_tab(self, force: bool = True) -> None:
        widget = self._tabs.currentWidget()
        if not force and widget not in self._dirty_canvases:
            return
        redraw = {
            self._canvas_rho_phi: self._redraw_rho_phi,
            self._canvas_rho_ps: self._redraw_rho_pseudosection,
            self._canvas_ph_ps: self._redraw_phase_pseudosection,
            self._canvas_tipper: self._redraw_tipper,
            self._canvas_pt: self._redraw_phase_tensor,
            self._canvas_pt_strip: self._redraw_phase_tensor_strip,
        }.get(widget)
        if redraw is not None:
            redraw()
            self._dirty_canvases.discard(widget)

    def current_canvas(self):
        """Return the canvas shown by the active tab."""
        widget = self._tabs.currentWidget()
        return widget if isinstance(widget, MplCanvas) else None

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
        """Redraw the Phase Tensor tab — skips the full matplotlib draw when
        the current plot key matches the last drawn key (tab revisit with no
        setting change).  Only the expensive initial draw and explicit refreshes
        trigger a full recompute + redraw cycle.
        """
        current_key = self._ctrl.phase_tensor_key()
        if (
            self._pt_last_key is not None
            and current_key == self._pt_last_key
            and len(self._canvas_pt.figure.axes) > 0
        ):
            # Nothing changed — just repaint the existing figure (fast path).
            self._canvas_pt.draw()
            return

        fig = self._canvas_pt.figure
        fig.clear()
        ax = fig.add_subplot(111)
        self._canvas_pt.axes = ax
        self._ctrl.draw_phase_tensor(ax)
        self._canvas_pt.draw()
        # Record the key so the next identical tab visit is free.
        self._pt_last_key = current_key

    def _redraw_phase_tensor_strip(self) -> None:
        fig = self._canvas_pt_strip.figure
        fig.clear()
        ax = fig.add_subplot(111)
        self._canvas_pt_strip.axes = ax
        self._ctrl.draw_phase_tensor_strip(ax)
        self._canvas_pt_strip.draw()
