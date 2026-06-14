# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for PlotController (Phase 3) — no Qt required."""

from __future__ import annotations

import pytest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.app.desktop.controllers.plot_controller import (
    PlotController, style_axes, _annotate_empty,
)


# ── Helpers ───────────────────────────────────────────────────────────────

def make_fig_and_axes():
    fig, ax = plt.subplots()
    return fig, ax


# ── Construction ──────────────────────────────────────────────────────────

def test_plot_controller_creates():
    assert PlotController() is not None


def test_default_dark_mode():
    ctrl = PlotController()
    assert ctrl.dark is True


def test_default_sites_is_none():
    assert PlotController()._sites is None


def test_default_station_is_none():
    assert PlotController()._station_id is None


# ── State setters ─────────────────────────────────────────────────────────

def test_set_sites():
    ctrl = PlotController()
    ctrl.set_sites("mock_sites")
    assert ctrl._sites == "mock_sites"


def test_set_station():
    ctrl = PlotController()
    ctrl.set_station("S05")
    assert ctrl._station_id == "S05"


def test_set_period_range():
    ctrl = PlotController()
    ctrl.set_period_range(0.001, 10.0)
    assert ctrl._period_range == (0.001, 10.0)


def test_set_period_range_none():
    ctrl = PlotController()
    ctrl.set_period_range(None, None)
    assert ctrl._period_range is None


def test_clear():
    ctrl = PlotController()
    ctrl.set_sites("x")
    ctrl.set_station("S01")
    ctrl.clear()
    assert ctrl._sites is None
    assert ctrl._station_id is None


# ── draw_rho_phi — no data ────────────────────────────────────────────────

def test_draw_rho_phi_no_data_does_not_raise():
    ctrl = PlotController()
    fig, _ = make_fig_and_axes()
    ctrl.draw_rho_phi(fig)   # must not raise
    plt.close(fig)


def test_draw_rho_phi_creates_two_axes_when_no_data():
    ctrl = PlotController()
    fig = plt.figure()
    ctrl.draw_rho_phi(fig)
    assert len(fig.axes) == 2
    plt.close(fig)


def test_draw_rho_phi_no_station_shows_prompt():
    ctrl = PlotController()
    ctrl.set_sites(object())   # truthy mock
    fig = plt.figure()
    ctrl.draw_rho_phi(fig)
    texts = [t.get_text() for ax in fig.axes for t in ax.texts]
    assert any("Select a station" in t for t in texts)
    plt.close(fig)


# ── draw_rho_pseudosection — no data ─────────────────────────────────────

def test_draw_rho_pseudosection_no_data_does_not_raise():
    ctrl = PlotController()
    fig, ax = make_fig_and_axes()
    ctrl.draw_rho_pseudosection(ax)
    plt.close(fig)


def test_draw_rho_pseudosection_shows_load_message():
    ctrl = PlotController()
    fig, ax = make_fig_and_axes()
    ctrl.draw_rho_pseudosection(ax)
    texts = [t.get_text() for t in ax.texts]
    assert any("Load" in t for t in texts)
    plt.close(fig)


# ── draw_phase_pseudosection — no data ────────────────────────────────────

def test_draw_phase_pseudosection_no_data_does_not_raise():
    ctrl = PlotController()
    fig, ax = make_fig_and_axes()
    ctrl.draw_phase_pseudosection(ax)
    plt.close(fig)


# ── draw_tipper — no data ─────────────────────────────────────────────────

def test_draw_tipper_no_data_does_not_raise():
    ctrl = PlotController()
    fig, ax = make_fig_and_axes()
    ctrl.draw_tipper(ax)
    plt.close(fig)


# ── draw_phase_tensor — no data ───────────────────────────────────────────

def test_draw_phase_tensor_no_data_does_not_raise():
    ctrl = PlotController()
    fig, ax = make_fig_and_axes()
    ctrl.draw_phase_tensor(ax)
    plt.close(fig)


# ── style_axes ────────────────────────────────────────────────────────────

def test_style_axes_dark_sets_facecolor():
    fig, ax = make_fig_and_axes()
    style_axes(ax, dark=True)
    assert ax.get_facecolor() != (1, 1, 1, 1)   # not white
    plt.close(fig)


def test_style_axes_light_does_not_raise():
    fig, ax = make_fig_and_axes()
    style_axes(ax, dark=False)
    plt.close(fig)


# ── _annotate_empty ───────────────────────────────────────────────────────

def test_annotate_empty_adds_text():
    fig, ax = make_fig_and_axes()
    _annotate_empty(ax, "Test message")
    texts = [t.get_text() for t in ax.texts]
    assert "Test message" in texts
    plt.close(fig)
