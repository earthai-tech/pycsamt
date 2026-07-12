# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.utils.mathext`.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

import pycsamt.utils.mathext as mathext
from pycsamt.exceptions import ResistivityError
from pycsamt.utils.mathext import (
    betaj,
    d_hanning_window,
    linkage_matrix,
    rhoa2z,
    rhophi2z,
    round_dipole_length,
    savitzky_golay1d,
    z2rhoa,
)


# ------------------------- round_dipole_length ------------------------

@pytest.mark.parametrize(
    "value, expected",
    [(12.0, 10.0), (13.0, 15.0), (15.0, 15.0), (0.0, 0.0), (-12.0, -10.0)],
)
def test_round_dipole_length(value, expected):
    assert round_dipole_length(value) == expected


def test_round_dipole_length_rejects_non_numeric():
    with pytest.raises(TypeError):
        round_dipole_length("50")


# --------------------------- compute_azimuth --------------------------

@pytest.fixture()
def stub_utm_to_ll(monkeypatch):
    """Deterministic UTM->lat/lon stub (1 km = 1 degree)."""

    def _stub(*, reference_ellipsoid, northing, easting, zone):
        lat = np.asarray(northing, dtype=float) / 1000.0 - 90.0
        lon = np.asarray(easting, dtype=float) / 1000.0 - 180.0
        return lat, lon

    monkeypatch.setattr(mathext, "utm_to_ll", _stub)


def test_compute_azimuth_due_north(stub_utm_to_ll):
    northing = np.array([80000.0, 81000.0, 82000.0])
    easting = np.full(3, 100000.0)
    az = mathext.compute_azimuth(easting, northing)
    assert az.shape == (2,)
    assert np.allclose(az, 0.0, atol=0.01)


def test_compute_azimuth_due_east(stub_utm_to_ll):
    northing = np.full(3, 90000.0)  # equator in stub coords
    easting = np.array([100000.0, 101000.0, 102000.0])
    az = mathext.compute_azimuth(easting, northing)
    assert np.allclose(az, 90.0, atol=0.5)


def test_compute_azimuth_extrapolate_prepends(stub_utm_to_ll):
    northing = np.array([80000.0, 81000.0, 82000.0])
    easting = np.full(3, 100000.0)
    az = mathext.compute_azimuth(easting, northing, extrapolate=True)
    assert az.shape == (3,)


def test_compute_azimuth_input_validation(stub_utm_to_ll):
    assert mathext.compute_azimuth([1.0], [2.0]).size == 0
    with pytest.raises(ValueError):
        mathext.compute_azimuth([1.0, 2.0], [1.0, 2.0, 3.0])


# ---------------------------- linkage_matrix --------------------------

def _cluster_frame():
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        rng.normal(size=(6, 3)), columns=list("abc")
    )


@pytest.mark.parametrize("kind", ["design", "squareform", "condense"])
def test_linkage_matrix_kinds_shape(kind):
    df = _cluster_frame()
    lm = linkage_matrix(df, kind=kind)
    assert lm.shape == (len(df) - 1, 4)


def test_linkage_matrix_as_frame_and_column_checks():
    df = _cluster_frame()
    lm = linkage_matrix(df.values, columns=list("abc"), as_frame=True)
    assert isinstance(lm, pd.DataFrame)
    assert lm.shape == (5, 4)

    with pytest.raises(TypeError):
        linkage_matrix(df.values, columns=["only-one"])
    with pytest.raises(ValueError):
        linkage_matrix(df, kind="unknown-kind")


# --------------------------- hanning / betaj --------------------------

def test_d_hanning_window_center_and_outside():
    W = 5
    assert d_hanning_window(0.0, 0.0, W) == pytest.approx(2.0 / W)
    assert d_hanning_window(3.0, 0.0, W) == 0.0


def test_betaj_reference_value():
    assert betaj(xj=2, L=1, W=5) == pytest.approx(
        0.35136534572813144, rel=1e-9
    )


def test_betaj_window_smaller_than_dipole_raises():
    with pytest.raises(ValueError):
        betaj(xj=0, L=10, W=5)


# ----------------------------- rhoa2z / z2rhoa ------------------------

def test_z2rhoa_documented_value():
    z = np.array([[2 + 3j]])
    f = np.array([1014.0])
    rhoa = z2rhoa(z, f)
    assert rhoa[0, 0] == pytest.approx(1623.73691735, rel=1e-6)


def test_rhoa2z_z2rhoa_roundtrip():
    z = np.array([[2 + 3j]])
    f = np.array([1014.0])
    rhoa = z2rhoa(z, f)
    phs = np.rad2deg(np.angle(z))
    z_back = rhoa2z(rhoa, phs, f)
    assert np.allclose(z_back, z)


def test_rhoa2z_and_z2rhoa_length_checks():
    with pytest.raises(ValueError):
        rhoa2z(np.ones((2, 1)), np.ones((3, 1)), np.ones(2))
    with pytest.raises(ValueError):
        z2rhoa(np.ones((2, 1), dtype=complex), np.ones(3))


# ------------------------------ rhophi2z ------------------------------

def test_rhophi2z_scalar_documented_example():
    z = rhophi2z(823, 25, 500)
    assert z.shape == (1,)
    assert z[0].real == pytest.approx(1300.00682824, rel=1e-6)
    assert z[0].imag == pytest.approx(606.20313966, rel=1e-6)


def test_rhophi2z_2x2_with_scalar_freq():
    rho = np.array([[823.0, 700.0], [723.0, 526.0]])
    phi = np.array([[45.0, 50.0], [90.0, 180.0]])
    z = rhophi2z(rho, phi, 500.0)
    assert z.shape == (2, 2)
    assert np.abs(z[0, 0]) == pytest.approx(
        math.sqrt(5 * 500 * 823), rel=1e-9
    )


def test_rhophi2z_grid_matches_freq_rows():
    rng = np.random.default_rng(1)
    rho = np.abs(rng.normal(size=(4, 3))) * 100 + 1
    phi = np.abs(rng.normal(size=(4, 3))) * 60
    freq = np.abs(rng.normal(size=4)) * 100 + 1
    z = rhophi2z(rho, phi, freq)
    assert z.shape == (4, 3)
    assert np.all(np.abs(z) > 0)


def test_rhophi2z_inconsistent_2x2_raises():
    with pytest.raises(ResistivityError):
        rhophi2z(np.ones((2, 2)), np.ones(3), 100.0)


# ---------------------------- savitzky_golay1d ------------------------

def test_savitzky_golay1d_smooths_noise():
    rng = np.random.default_rng(42)
    t = np.linspace(-4, 4, 200)
    clean = np.exp(-(t**2))
    noisy = clean + rng.normal(0, 0.05, t.shape)
    smoothed = savitzky_golay1d(noisy, window_size=31, order=4,
                                mode="valid")
    assert smoothed.shape == noisy.shape
    # smoothing should reduce the residual w.r.t. the clean signal
    assert (
        np.linalg.norm(smoothed - clean)
        < np.linalg.norm(noisy - clean)
    )


def test_savitzky_golay1d_derivative_of_line():
    y = 2.0 * np.arange(50.0)
    dy = savitzky_golay1d(y, window_size=7, order=2, deriv=1,
                          mode="valid")
    interior = dy[3:-3]
    assert np.allclose(interior, 2.0, atol=1e-8)


def test_savitzky_golay1d_input_validation():
    y = np.arange(10.0)
    with pytest.raises(TypeError):
        savitzky_golay1d(y, window_size=4, order=2)  # even window
    with pytest.raises(TypeError):
        savitzky_golay1d(y, window_size=3, order=5)  # order too high
    with pytest.raises(ValueError):
        savitzky_golay1d(y, window_size="w", order=2)
