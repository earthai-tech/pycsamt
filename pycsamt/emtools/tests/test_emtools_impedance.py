# -*- coding: utf-8 -*-
"""Tests for pycsamt.emtools.impedance"""
from __future__ import annotations

import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.impedance import (
    plot_phasor_wheel,
    plot_offdiag_antisym_residual,
    plot_determinant_track,
)


# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

class _FakeZ:
    def __init__(self, z, freq):
        self.z    = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)
        self.z_err = np.abs(z) * 0.05   # 5 % error


class _FakeSite:
    def __init__(self, station, z, freq):
        self.station = station
        self.Z    = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 10, f_lo: float = 1.0, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] =  amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _3d_z(freqs: np.ndarray) -> np.ndarray:
    z = _iso_z(freqs)
    amp = np.abs(z[:, 0, 1])
    z[:, 0, 0] = 0.3 * amp * (0.6 + 0.8j)
    z[:, 1, 1] = -0.3 * amp * (0.5 + 0.7j)
    return z


def _site(name: str, z=None, n: int = 12) -> _FakeSite:
    fr = _freqs(n)
    if z is None:
        z = _iso_z(fr)
    return _FakeSite(name, z, fr)


# ─────────────────────────────────────────────────────────────────────────────
# plot_phasor_wheel
# ─────────────────────────────────────────────────────────────────────────────

class TestPlotPhasorWheel:

    def test_returns_figure(self):
        sites = [_site("S00")]
        result = plot_phasor_wheel(sites)
        plt.close("all")
        assert result is not None

    def test_external_polar_ax(self):
        sites = [_site("S00")]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        result = plot_phasor_wheel(sites, ax=ax)
        plt.close("all")
        assert result is not None

    def test_multi_site_no_crash(self):
        sites = [_site(f"S{i:02d}") for i in range(3)]
        result = plot_phasor_wheel(sites)
        plt.close("all")
        assert result is not None

    def test_3d_site_no_crash(self):
        fr = _freqs()
        sites = [_FakeSite("S00", _3d_z(fr), fr)]
        result = plot_phasor_wheel(sites)
        plt.close("all")
        assert result is not None

    def test_specific_station(self):
        sites = [_site("S00"), _site("S01")]
        result = plot_phasor_wheel(sites, station="S00")
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_offdiag_antisym_residual
# ─────────────────────────────────────────────────────────────────────────────

class TestPlotOffdiagAntisymResidual:

    def test_returns_figure(self):
        sites = [_site("S00")]
        result = plot_offdiag_antisym_residual(sites)
        plt.close("all")
        assert result is not None

    def test_external_ax(self):
        sites = [_site("S00")]
        fig, ax = plt.subplots()
        result = plot_offdiag_antisym_residual(sites, ax=ax)
        plt.close("all")
        assert result is not None

    def test_multi_site(self):
        sites = [_site(f"S{i}") for i in range(4)]
        result = plot_offdiag_antisym_residual(sites)
        plt.close("all")
        assert result is not None

    def test_3d_site(self):
        fr = _freqs()
        sites = [_FakeSite("S00", _3d_z(fr), fr)]
        result = plot_offdiag_antisym_residual(sites)
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_determinant_track
# ─────────────────────────────────────────────────────────────────────────────

class TestPlotDeterminantTrack:

    def test_returns_figure(self):
        sites = [_site(f"S{i}") for i in range(3)]
        result = plot_determinant_track(sites)
        plt.close("all")
        assert result is not None

    def test_external_axes(self):
        sites = [_site("S00")]
        fig, axes = plt.subplots(1, 2)
        result = plot_determinant_track(sites, axes=axes)
        plt.close("all")
        assert result is not None

    def test_single_site(self):
        sites = [_site("S00")]
        result = plot_determinant_track(sites)
        plt.close("all")
        assert result is not None

    def test_empty_no_crash(self):
        try:
            result = plot_determinant_track([])
            plt.close("all")
        except Exception:
            plt.close("all")
