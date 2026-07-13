"""Tests for rho/phase trend smoothing."""

from __future__ import annotations

import numpy as np

from pycsamt.emtools import smooth_rho_phase
from pycsamt.emtools._core import _iter_items


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 16) -> np.ndarray:
    return np.logspace(-1.0, 3.0, n)


def _z_from_rho_phase(freq, rho, phase_deg, *, sign=1.0):
    amp = np.sqrt(5.0 * np.asarray(freq, float) * np.asarray(rho, float))
    return sign * amp * np.exp(1j * np.deg2rad(phase_deg))


def _rho_from_z(freq, z):
    return (np.abs(z) ** 2) / (5.0 * np.asarray(freq, float))


def _roughness(y):
    y = np.asarray(y, dtype=float)
    good = np.isfinite(y) & (y > 0.0)
    if np.count_nonzero(good) < 4:
        return np.nan
    return float(np.nanstd(np.diff(np.diff(np.log10(y[good])))))


def _noisy_site():
    fr = _freqs()
    x = np.log10(fr)
    trend = 80.0 * (1.0 + 0.18 * (x - x.mean()) ** 2)
    wiggle = np.array(
        [
            1.00,
            1.30,
            0.74,
            1.18,
            0.82,
            1.25,
            0.78,
            1.15,
            0.85,
            1.22,
            0.80,
            1.18,
            0.86,
            1.26,
            0.76,
            1.08,
        ]
    )
    rho = trend * wiggle
    phase = (
        43.0
        + 7.0 * (x - x.mean())
        + np.array(
            [
                0.0,
                9.0,
                -7.0,
                6.0,
                -8.0,
                7.0,
                -6.0,
                5.0,
                -4.0,
                6.0,
                -5.0,
                8.0,
                -6.0,
                7.0,
                -8.0,
                0.0,
            ]
        )
    )
    z = np.zeros((fr.size, 2, 2), dtype=complex)
    z[:, 0, 1] = _z_from_rho_phase(fr, rho, phase)
    z[:, 1, 0] = _z_from_rho_phase(fr, rho * 1.05, phase + 180.0, sign=-1.0)
    z[:, 0, 0] = _z_from_rho_phase(fr, 12.0 * np.ones_like(fr), 20.0)
    z[:, 1, 1] = _z_from_rho_phase(fr, 14.0 * np.ones_like(fr), -15.0)
    return _FakeSite("S00", z, fr)


def _first_z(sites):
    ed = next(_iter_items(sites))
    raw = getattr(ed, "edi", ed)
    return raw.Z.z, raw.Z.freq


def test_smooth_rho_phase_returns_sites_and_preserves_shape():
    site = _noisy_site()
    result = smooth_rho_phase([site], components="offdiag", degree=3)
    z, fr = _first_z(result)
    assert z.shape == site.Z.z.shape
    assert fr.shape == site.Z.freq.shape


def test_smooth_rho_phase_reduces_log_rho_roughness():
    site = _noisy_site()
    before = _rho_from_z(site.Z.freq, site.Z.z[:, 0, 1])
    result = smooth_rho_phase([site], components="xy", degree=3)
    z, fr = _first_z(result)
    after = _rho_from_z(fr, z[:, 0, 1])
    assert _roughness(after) < _roughness(before)


def test_smooth_rho_phase_respects_component_selection():
    site = _noisy_site()
    diag_before = site.Z.z[:, 0, 0].copy()
    result = smooth_rho_phase([site], components="xy", degree=2)
    z, _ = _first_z(result)
    np.testing.assert_allclose(z[:, 0, 0], diag_before)
    assert not np.allclose(z[:, 0, 1], site.Z.z[:, 0, 1])


def test_smooth_rho_phase_phase_only_keeps_rho_amplitude():
    site = _noisy_site()
    before_rho = _rho_from_z(site.Z.freq, site.Z.z[:, 0, 1])
    before_phase = np.rad2deg(np.angle(site.Z.z[:, 0, 1]))
    result = smooth_rho_phase(
        [site],
        components="xy",
        degree=2,
        smooth_rho=False,
        smooth_phase=True,
    )
    z, fr = _first_z(result)
    after_rho = _rho_from_z(fr, z[:, 0, 1])
    after_phase = np.rad2deg(np.angle(z[:, 0, 1]))
    np.testing.assert_allclose(after_rho, before_rho, rtol=1e-10, atol=1e-10)
    assert np.nanstd(after_phase - before_phase) > 0.0
