# -*- coding: utf-8 -*-
"""Tests for pycsamt.emtools.tf (transfer-function / tipper)"""
from __future__ import annotations

import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.tf import (
    plot_tipper_hodograms,
    plot_induction_arrows,
    plot_induction_convention,
    plot_tipper_polar,
    plot_induction_rose,
)


# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

class _FakeZ:
    def __init__(self, z, freq):
        self.z    = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeTipper:
    def __init__(self, tipper, freq):
        # tipper shape: (n, 2) complex — t[:, 0] = Tx, t[:, 1] = Ty
        self.tipper = np.asarray(tipper, dtype=complex)
        self.freq   = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq, *, tipper=None,
                 east=None, north=None):
        self.station = station
        self.Z    = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if tipper is not None:
            self.Tipper = _FakeTipper(tipper, freq)
        if east is not None:
            self.east  = float(east)
            self.north = float(north)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 12, f_lo: float = 0.1, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] =  amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _tipper(freqs: np.ndarray, amp: float = 0.15) -> np.ndarray:
    """Synthetic (n, 2) complex tipper array."""
    t = np.zeros((freqs.size, 2), dtype=complex)
    t[:, 0] = amp * np.exp(1j * np.linspace(0, np.pi / 4, freqs.size))
    t[:, 1] = amp * 0.5 * np.exp(1j * np.linspace(0, np.pi / 6, freqs.size))
    return t


def _site(name: str, n: int = 12, *, with_tipper=True,
          east=None, north=None) -> _FakeSite:
    fr = _freqs(n)
    tip = _tipper(fr) if with_tipper else None
    return _FakeSite(name, _iso_z(fr), fr, tipper=tip,
                     east=east, north=north)


def _profile(n_sites: int = 4) -> list:
    return [_site(f"S{i:02d}", east=i * 200.0, north=0.0)
            for i in range(n_sites)]


# ─────────────────────────────────────────────────────────────────────────────
# plot_tipper_hodograms
# ─────────────────────────────────────────────────────────────────────────────

class TestPlotTipperHodograms:

    def test_returns_figure_with_tipper(self):
        sites = [_site("S00", with_tipper=True)]
        fig = plot_tipper_hodograms(sites)
        plt.close("all")
        assert fig is not None

    def test_no_tipper_graceful(self):
        sites = [_site("S00", with_tipper=False)]
        fig = plot_tipper_hodograms(sites)
        plt.close("all")
        assert fig is not None

    def test_external_axes(self):
        sites = [_site("S00")]
        fig_in, (ax1, ax2) = plt.subplots(1, 2)
        result = plot_tipper_hodograms(sites, axes=[ax1, ax2])
        plt.close("all")
        assert result is not None

    def test_multi_band(self):
        sites = [_site("S00")]
        fig = plot_tipper_hodograms(sites, n_bands=4)
        plt.close("all")
        assert fig is not None

    def test_station_selection(self):
        sites = [_site("S00"), _site("S01")]
        fig = plot_tipper_hodograms(sites, station="S00")
        plt.close("all")
        assert fig is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_induction_arrows
# ─────────────────────────────────────────────────────────────────────────────

class TestPlotInductionArrows:

    def test_returns_axes(self):
        sites = _profile(4)
        result = plot_induction_arrows(sites, periods=[1.0])
        plt.close("all")
        assert result is not None

    def test_external_ax(self):
        sites = _profile(4)
        fig, ax = plt.subplots()
        result = plot_induction_arrows(sites, periods=[1.0], ax=ax)
        plt.close("all")
        assert result is not None

    def test_no_tipper_graceful(self):
        sites = [_site(f"S{i}", with_tipper=False, east=i*200.0, north=0.0)
                 for i in range(3)]
        result = plot_induction_arrows(sites, periods=[1.0])
        plt.close("all")
        assert result is not None

    def test_convention_wiese(self):
        sites = _profile(3)
        result = plot_induction_arrows(sites, periods=[1.0], convention="wiese")
        plt.close("all")
        assert result is not None

    def test_multiple_periods(self):
        sites = _profile(3)
        result = plot_induction_arrows(sites, periods=[0.1, 1.0, 10.0])
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_induction_convention
# ─────────────────────────────────────────────────────────────────────────────

class TestPlotInductionConvention:

    def test_returns_figure(self):
        sites = _profile(3)
        result = plot_induction_convention(sites)
        plt.close("all")
        assert result is not None

    def test_external_axes(self):
        sites = _profile(3)
        fig, axes = plt.subplots(2, 2)
        result = plot_induction_convention(sites, axes=axes.ravel())
        plt.close("all")
        assert result is not None

    def test_custom_period(self):
        sites = _profile(3)
        result = plot_induction_convention(sites, period=0.1)
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_tipper_polar
# ─────────────────────────────────────────────────────────────────────────────

class TestPlotTipperPolar:

    def test_returns_figure(self):
        sites = [_site("S00")]
        result = plot_tipper_polar(sites)
        plt.close("all")
        assert result is not None

    def test_no_tipper_graceful(self):
        sites = [_site("S00", with_tipper=False)]
        result = plot_tipper_polar(sites)
        plt.close("all")
        assert result is not None

    def test_station_selection(self):
        sites = [_site("S00"), _site("S01")]
        result = plot_tipper_polar(sites, station="S01")
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_induction_rose
# ─────────────────────────────────────────────────────────────────────────────

class TestPlotInductionRose:

    def test_returns_figure(self):
        sites = _profile(4)
        result = plot_induction_rose(sites)
        plt.close("all")
        assert result is not None

    def test_no_tipper_graceful(self):
        sites = [_site(f"S{i}", with_tipper=False, east=i*200.0, north=0.0)
                 for i in range(3)]
        result = plot_induction_rose(sites)
        plt.close("all")
        assert result is not None

    def test_external_ax(self):
        sites = _profile(3)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        result = plot_induction_rose(sites, ax=ax)
        plt.close("all")
        assert result is not None
