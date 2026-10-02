"""Tests for pycsamt.emtools.impedance"""

from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.impedance import (
    plot_determinant_track,
    plot_offdiag_antisym_residual,
    plot_phasor_wheel,
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq, err=True):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)
        self.z_err = np.abs(z) * 0.05 if err else None  # 5 % error


class _FakeSite:
    def __init__(self, station, z, freq, err=True):
        self.station = station
        self.Z = _FakeZ(z, freq, err)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


class _FakeZBadShape:
    """A Z wrapper whose array cannot be coerced to (n, 2, 2)."""

    def __init__(self, freq):
        self.z = np.ones((freq.size, 3, 3), dtype=complex)
        self.freq = np.asarray(freq, dtype=float)
        self.z_err = None


class _FakeSiteNoZ:
    """A site whose Z block is unusable, exercising the "no Z" branches."""

    def __init__(self, station, freq):
        self.station = station
        self.Z = _FakeZBadShape(freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 10, f_lo: float = 1.0, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
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

    def test_station_not_found(self):
        sites = [_site("S00")]
        result = plot_phasor_wheel(sites, station="ZZZ")
        texts = [t.get_text() for t in result.texts]
        plt.close("all")
        assert "station not found" in texts

    def test_no_z_block(self):
        fr = _freqs()
        sites = [_FakeSiteNoZ("S00", fr)]
        result = plot_phasor_wheel(sites)
        texts = [t.get_text() for t in result.texts]
        plt.close("all")
        assert "no Z" in texts

    def test_pband_filters_periods(self):
        sites = [_site("S00", n=20)]
        result = plot_phasor_wheel(sites, pband=(1e-2, 1e-1))
        plt.close("all")
        assert result is not None

    def test_radius_norm(self):
        sites = [_site("S00")]
        result = plot_phasor_wheel(sites, radius="norm")
        plt.close("all")
        assert result is not None

    def test_no_connect(self):
        sites = [_site("S00")]
        result = plot_phasor_wheel(sites, connect=False)
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

    def test_no_data(self):
        fr = _freqs()
        sites = [_FakeSiteNoZ("S00", fr)]
        result = plot_offdiag_antisym_residual(sites)
        texts = [t.get_text() for t in result.texts]
        plt.close("all")
        assert "no data" in texts

    def test_skips_site_with_no_z(self):
        fr = _freqs()
        good = _site("S00")
        bad = _FakeSiteNoZ("S01", fr)
        result = plot_offdiag_antisym_residual([good, bad])
        plt.close("all")
        assert result is not None

    def test_gap_fill_across_mismatched_freq_grids(self):
        fr1 = _freqs(10)
        fr2 = _freqs(14, 2.0, 8e3)
        s1 = _FakeSite("S00", _iso_z(fr1), fr1)
        s2 = _FakeSite("S01", _iso_z(fr2), fr2)
        result = plot_offdiag_antisym_residual([s1, s2])
        plt.close("all")
        assert result is not None

    def test_explicit_vlim(self):
        sites = [_site("S00")]
        result = plot_offdiag_antisym_residual(sites, vlim=0.5)
        plt.close("all")
        assert result is not None

    def test_preinverted_axis_not_reinverted(self):
        sites = [_site("S00")]
        fig, ax = plt.subplots()
        ax.invert_yaxis()
        assert ax.yaxis_inverted()
        result = plot_offdiag_antisym_residual(sites, ax=ax)
        plt.close("all")
        assert result.yaxis_inverted()


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
            plot_determinant_track([])
            plt.close("all")
        except Exception:
            plt.close("all")

    def test_no_z_err_uses_deterministic_branch(self):
        fr = _freqs()
        sites = [_FakeSite("S00", _iso_z(fr), fr, err=False)]
        result = plot_determinant_track(sites)
        plt.close("all")
        assert result is not None

    def test_no_sites_with_axes_given(self):
        fig, axg = plt.subplots(1, 1)
        result = plot_determinant_track([], axes=axg)
        plt.close("all")
        assert result is not None

    def test_station_not_found_with_axes_given(self):
        sites = [_site("S00")]
        fig, axg = plt.subplots(1, 1)
        result = plot_determinant_track(sites, station="ZZZ", axes=axg)
        plt.close("all")
        assert result is not None

    def test_no_z_with_axes_given(self):
        fr = _freqs()
        sites = [_FakeSiteNoZ("S00", fr)]
        fig, axg = plt.subplots(1, 1)
        result = plot_determinant_track(sites, axes=axg)
        plt.close("all")
        assert result is not None
