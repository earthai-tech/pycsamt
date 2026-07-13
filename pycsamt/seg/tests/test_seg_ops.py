# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.seg.ops import (
    MU0,
    amp_or_psd,
    coherence_ms,
    pack_hermitian,
    rho_phi_to_z,
    rotate_impedance,
    rotate_spectra,
    rotate_tipper,
    time_vector,
    unpack_hermitian,
    z_to_rho_phi,
)


def _R(theta_deg: float) -> np.ndarray:
    th = np.deg2rad(theta_deg)
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, s], [-s, c]], float)


def test_rotate_impedance_basic() -> None:
    Z = np.array([[[1 + 2j, 3 + 4j], [5 + 6j, 7 + 8j]]], complex)
    th = 90.0
    R = _R(th)
    exp = (R @ Z[0] @ R.T)[None, ...]
    out = rotate_impedance(Z, th)
    assert np.allclose(out, exp)


def test_rotate_tipper_basic() -> None:
    t = np.array([1.0, 0.0])
    th = 90.0
    R = _R(th)
    exp = (R @ t).ravel()
    out = rotate_tipper(t, th)
    assert np.allclose(out, exp)


def test_z_to_rho_phi_roundtrip() -> None:
    f = np.array([1.0, 10.0, 100.0])
    Z = np.zeros((3, 2, 2), complex)
    Z[:, 0, 0] = 1.0 + 1.0j
    Z[:, 0, 1] = 0.5 - 0.2j
    Z[:, 1, 0] = -0.3 + 0.7j
    Z[:, 1, 1] = 2.0 - 1.0j

    rho, phi = z_to_rho_phi(Z, f, mu0=MU0)
    Z2 = rho_phi_to_z(rho, phi, f, mu0=MU0)

    # magnitude/phase must match
    m1 = np.abs(Z)
    m2 = np.abs(Z2)
    p1 = np.angle(Z)
    p2 = np.angle(Z2)
    assert np.allclose(m1, m2, rtol=1e-12, atol=1e-12)
    assert np.allclose(p1, p2, rtol=1e-12, atol=1e-12)


def test_time_vector() -> None:
    t = time_vector(5, 0.2)
    assert np.allclose(t, [0.0, 0.2, 0.4, 0.6, 0.8])


def test_coherence_ms_limits_and_value() -> None:
    # perfect coherence: |Sxy|^2 = Sxx*Syy
    Sxx = np.array([2.0, 3.0])
    Syy = np.array([4.0, 5.0])
    Sxy = np.sqrt(Sxx * Syy) * np.exp(1j * 0.3)
    coh = coherence_ms(Sxy, Sxx, Syy)
    assert np.allclose(coh, [1.0, 1.0])
    # zeros in denom → 0
    coh2 = coherence_ms(np.array([1.0]), np.array([0.0]), np.array([2.0]))
    assert float(coh2[0]) == 0.0


def test_pack_unpack_hermitian() -> None:
    # 3x3 Hermitian with real diagonal
    A = np.array(
        [
            [1.0, 1 + 2j, 0.5 - 1j],
            [1 - 2j, 3.0, -0.2 + 0.4j],
            [0.5 + 1j, -0.2 - 0.4j, 2.0],
        ],
        complex,
    )
    v = pack_hermitian(A)
    B = unpack_hermitian(v, 3)
    assert np.allclose(A, B)
    # diagonal must be real
    assert np.allclose(np.diag(B).imag, 0.0)


def test_rotate_spectra_against_manual() -> None:
    C = np.array([[2.0 + 0j, 1.0 + 1.0j], [1.0 - 1.0j, 1.5 + 0j]], complex)
    th = 45.0
    R = _R(th)
    exp = R @ C @ R.T
    out = rotate_spectra(C, th)
    assert np.allclose(out, exp)


def test_amp_or_psd_sine_peak() -> None:
    fs = 128.0
    n = 1024
    t = np.arange(n) / fs
    f0 = 8.0
    x = np.sin(2 * np.pi * f0 * t)

    f, A = amp_or_psd(x, fs, mode="amp")
    # peak bin near f0 has amplitude ~ 0.5
    k = int(np.argmin(np.abs(f - f0)))
    assert A[k] == pytest.approx(0.5, rel=0.15)

    f2, P = amp_or_psd(x, fs, mode="psd")
    assert f2.shape == P.shape
    assert np.all(P >= 0.0)
    assert P[int(np.argmin(np.abs(f2 - f0)))] > 0.0
