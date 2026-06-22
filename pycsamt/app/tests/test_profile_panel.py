# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ProfilePanel (Phase 3)."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.panels.profile_panel import ProfilePanel
from pycsamt.app.desktop.widgets.freq_selector import FreqSelector


@pytest.fixture
def panel(qapp):
    p = ProfilePanel()
    yield p
    p.close()


# ── Construction ──────────────────────────────────────────────────────────

def test_profile_panel_creates(qapp):
    p = ProfilePanel()
    assert p is not None
    p.close()


def test_has_freq_selector(panel):
    assert isinstance(panel._freq_sel, FreqSelector)


def test_has_tab_widget(panel):
    assert panel._tabs is not None


def test_tab_count(panel):
    assert panel._tabs.count() == 6


def test_tab_labels(panel):
    labels = [panel._tabs.tabText(i) for i in range(panel._tabs.count())]
    assert "ρₐ / φ" in labels
    assert "Pseudosection ρₐ" in labels
    assert "Tipper" in labels
    assert "Phase Tensor" in labels
    assert "2D Section" in labels


# ── Canvases ──────────────────────────────────────────────────────────────

def test_rho_phi_canvas_exists(panel):
    assert panel._canvas_rho_phi is not None


def test_pseudosection_canvas_exists(panel):
    assert panel._canvas_rho_ps is not None


def test_phase_ps_canvas_exists(panel):
    assert panel._canvas_ph_ps is not None


def test_tipper_canvas_exists(panel):
    assert panel._canvas_tipper is not None


def test_phase_tensor_canvas_exists(panel):
    assert panel._canvas_pt is not None


# ── PlotController linkage ────────────────────────────────────────────────

def test_plot_controller_created(panel):
    from pycsamt.app.desktop.controllers.plot_controller import PlotController
    assert isinstance(panel._ctrl, PlotController)


def test_set_dark_mode_updates_controller(panel):
    panel.set_dark_mode(False)
    assert panel._ctrl.dark is False
    panel.set_dark_mode(True)
    assert panel._ctrl.dark is True


# ── No-data draw calls ────────────────────────────────────────────────────

def test_set_sites_none_does_not_raise(panel):
    panel._ctrl.set_sites(None)
    panel._redraw_all()   # must not raise


def test_set_selected_station_no_data(panel):
    panel.set_selected_station("S01")  # no sites loaded — must not raise


def test_redraw_rho_phi_does_not_raise(panel):
    panel._redraw_rho_phi()


def test_redraw_rho_pseudosection_does_not_raise(panel):
    panel._redraw_rho_pseudosection()


def test_redraw_phase_pseudosection_does_not_raise(panel):
    panel._redraw_phase_pseudosection()


def test_redraw_tipper_does_not_raise(panel):
    panel._redraw_tipper()


def test_redraw_phase_tensor_does_not_raise(panel):
    panel._redraw_phase_tensor()


# ── Freq range signal ─────────────────────────────────────────────────────

def test_freq_range_change_updates_controller(panel):
    panel._on_freq_range_changed(1.0, 1000.0)
    T_min = 1.0 / 1000.0
    T_max = 1.0 / 1.0
    assert panel._ctrl._period_range is not None
    assert panel._ctrl._period_range[0] == pytest.approx(T_min)
    assert panel._ctrl._period_range[1] == pytest.approx(T_max)


# ── Tab switching ─────────────────────────────────────────────────────────

def test_tab_switch_does_not_raise(panel):
    for i in range(5):
        panel._tabs.setCurrentIndex(i)
        panel._on_tab_changed(i)
