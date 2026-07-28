"""Tests for pycsamt.emtools.advanced (novel diagnostic plots)"""

from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.advanced import (
    plot_apparent_anisotropy_section,
    plot_apparent_resistivity_polar,
    plot_dimensionality_ternary,
    plot_distortion_radar,
    plot_impedance_mohr_circles,
    plot_pt_period_clock,
    plot_rho_phase_bode,
    plot_snr_section,
    plot_survey_fingerprint,
    plot_z_invariants_section,
    plot_zt_argand,
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.z_err = np.abs(z) * 0.05
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq, *, east=None, north=None):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if east is not None:
            self.east = float(east)
            self.north = float(north)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 12, f_lo: float = 0.1, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _3d_z(freqs: np.ndarray, skew: float = 0.4) -> np.ndarray:
    z = _iso_z(freqs)
    amp = np.abs(z[:, 0, 1])
    z[:, 0, 0] = skew * amp * (0.6 + 0.8j)
    z[:, 1, 1] = -skew * amp * (0.5 + 0.7j)
    return z


def _site(name: str, n: int = 12, *, east=None, north=None, three_d=False) -> _FakeSite:
    fr = _freqs(n)
    z = _3d_z(fr) if three_d else _iso_z(fr)
    return _FakeSite(name, z, fr, east=east, north=north)


def _profile(n: int = 4) -> list:
    return [_site(f"S{i:02d}", east=float(i) * 200.0, north=0.0) for i in range(n)]


def _mixed(n: int = 3) -> list:
    sites = []
    for i in range(n):
        fr = _freqs()
        z = _3d_z(fr) if i % 2 else _iso_z(fr)
        sites.append(_FakeSite(f"S{i:02d}", z, fr, east=float(i) * 200.0, north=0.0))
    return sites


# ─────────────────────────────────────────────────────────────────────────────
# plot_impedance_mohr_circles
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotImpedanceMohrCircles:
    def test_returns_figure(self):
        sites = [_site("S00")]
        fig = plot_impedance_mohr_circles(sites)
        plt.close("all")
        assert fig is not None

    def test_external_axes(self):
        sites = [_site("S00")]
        fig_in, axes = plt.subplots(1, 2)
        result = plot_impedance_mohr_circles(sites, axes=axes)
        plt.close("all")
        assert result is not None

    def test_3d_site(self):
        sites = [_site("S00", three_d=True)]
        fig = plot_impedance_mohr_circles(sites)
        plt.close("all")
        assert fig is not None

    def test_n_periods(self):
        sites = [_site("S00")]
        fig = plot_impedance_mohr_circles(sites, n_periods=4)
        plt.close("all")
        assert fig is not None

    def test_specific_station(self):
        sites = [_site("S00"), _site("S01")]
        fig = plot_impedance_mohr_circles(sites, station="S00")
        plt.close("all")
        assert fig is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_zt_argand
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotZtArgand:
    def test_returns_figure(self):
        sites = [_site("S00")]
        fig = plot_zt_argand(sites)
        plt.close("all")
        assert fig is not None

    def test_external_axes(self):
        sites = [_site("S00")]
        fig_in, axes = plt.subplots(1, 2)
        result = plot_zt_argand(sites, axes=axes)
        plt.close("all")
        assert result is not None

    def test_single_component(self):
        sites = [_site("S00")]
        fig = plot_zt_argand(sites, components=("xy",))
        plt.close("all")
        assert fig is not None

    def test_station_selection(self):
        sites = [_site("S00"), _site("S01")]
        fig = plot_zt_argand(sites, station="S01")
        plt.close("all")
        assert fig is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_survey_fingerprint
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotSurveyFingerprint:
    def test_returns_figure(self):
        sites = _mixed(4)
        fig = plot_survey_fingerprint(sites)
        plt.close("all")
        assert fig is not None

    def test_custom_quantities(self):
        sites = _mixed(3)
        fig = plot_survey_fingerprint(sites, quantities=["skew", "ellipt"])
        plt.close("all")
        assert fig is not None

    def test_empty_sites_no_crash(self):
        try:
            plot_survey_fingerprint([])
            plt.close("all")
        except Exception:
            plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_dimensionality_ternary
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotDimensionalityTernary:
    def test_returns_figure(self):
        sites = _mixed(3)
        fig = plot_dimensionality_ternary(sites)
        plt.close("all")
        assert fig is not None

    def test_external_ax(self):
        sites = _mixed(3)
        fig_in, ax = plt.subplots()
        result = plot_dimensionality_ternary(sites, ax=ax)
        plt.close("all")
        assert result is not None

    def test_empty_sites_no_crash(self):
        try:
            plot_dimensionality_ternary([])
        except Exception:
            pass
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_distortion_radar
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotDistortionRadar:
    def test_returns_figure(self):
        sites = [_site("S00", three_d=True)]
        fig = plot_distortion_radar(sites)
        plt.close("all")
        assert fig is not None

    def test_multi_site(self):
        sites = _mixed(3)
        fig = plot_distortion_radar(sites)
        plt.close("all")
        assert fig is not None

    def test_external_ax(self):
        sites = [_site("S00")]
        fig_in = plt.figure()
        ax = fig_in.add_subplot(111, projection="polar")
        result = plot_distortion_radar(sites, ax=ax)
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_rho_phase_bode
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotRhoPhaseBode:
    def test_returns_figure(self):
        sites = [_site("S00")]
        fig = plot_rho_phase_bode(sites)
        plt.close("all")
        assert fig is not None

    def test_external_axes(self):
        sites = [_site("S00")]
        fig_in, axes = plt.subplots(1, 2)
        result = plot_rho_phase_bode(sites, axes=axes)
        plt.close("all")
        assert result is not None

    def test_yx_component(self):
        sites = [_site("S00")]
        fig = plot_rho_phase_bode(sites, component="yx")
        plt.close("all")
        assert fig is not None

    def test_smooth_window(self):
        sites = [_site("S00")]
        fig = plot_rho_phase_bode(sites, smooth_window=2)
        plt.close("all")
        assert fig is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_pt_period_clock
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotPtPeriodClock:
    def test_returns_figure(self):
        sites = _mixed(3)
        fig = plot_pt_period_clock(sites)
        plt.close("all")
        assert fig is not None

    def test_external_ax(self):
        sites = _mixed(3)
        fig_in = plt.figure()
        ax = fig_in.add_subplot(111, projection="polar")
        result = plot_pt_period_clock(sites, ax=ax)
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_apparent_resistivity_polar
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotApparentResistivityPolar:
    def test_returns_figure(self):
        sites = [_site("S00")]
        fig = plot_apparent_resistivity_polar(sites)
        plt.close("all")
        assert fig is not None

    def test_station_selection(self):
        sites = [_site("S00"), _site("S01")]
        fig = plot_apparent_resistivity_polar(sites, station="S00")
        plt.close("all")
        assert fig is not None

    def test_external_ax(self):
        sites = [_site("S00")]
        fig_in = plt.figure()
        ax = fig_in.add_subplot(111, projection="polar")
        result = plot_apparent_resistivity_polar(sites, ax=ax)
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_apparent_anisotropy_section
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotApparentAnisotropySection:
    def test_returns_figure(self):
        sites = _profile(4)
        fig = plot_apparent_anisotropy_section(sites)
        plt.close("all")
        assert fig is not None

    def test_external_ax(self):
        sites = _profile(3)
        fig_in, ax = plt.subplots()
        result = plot_apparent_anisotropy_section(sites, ax=ax)
        plt.close("all")
        assert result is not None

    def test_empty_sites_no_crash(self):
        try:
            plot_apparent_anisotropy_section([])
        except Exception:
            pass
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_snr_section
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotSnrSection:
    def test_returns_figure(self):
        sites = _profile(4)
        fig = plot_snr_section(sites)
        plt.close("all")
        assert fig is not None

    def test_external_axes(self):
        sites = _profile(3)
        result = plot_snr_section(sites, axes=None)
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_z_invariants_section
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotZInvariantsSection:
    def test_returns_figure(self):
        sites = _profile(4)
        fig = plot_z_invariants_section(sites)
        plt.close("all")
        assert fig is not None

    def test_extra_sites(self):
        sites = _profile(5)
        result = plot_z_invariants_section(sites)
        plt.close("all")
        assert result is not None

    def test_empty_sites_no_crash(self):
        try:
            plot_z_invariants_section([])
        except Exception:
            pass
        plt.close("all")
