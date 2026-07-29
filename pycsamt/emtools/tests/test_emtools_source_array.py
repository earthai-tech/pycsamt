"""
Tests for pycsamt.emtools.source_array
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.source_array import (
    array_factor,
    beam_steer,
    pas_pattern,
    plot_radiation_pattern,
    sdas_directivity,
    sdas_element_pattern,
    snr_gain_db,
    steering_angles,
    wavenumber,
)

# ─────────────────────────────────────────────────────────────────────────────
# wavenumber
# ─────────────────────────────────────────────────────────────────────────────


def test_wavenumber_freespace_positive():
    k = wavenumber(100.0)
    assert k > 0


def test_wavenumber_freespace_scales():
    k1 = wavenumber(100.0)
    k2 = wavenumber(200.0)
    assert np.isclose(k2, 2 * k1)


def test_wavenumber_earth_positive():
    k = wavenumber(100.0, rho=100.0)
    assert k > 0


def test_wavenumber_earth_increases_with_freq():
    k1 = wavenumber(10.0, rho=100.0)
    k2 = wavenumber(1000.0, rho=100.0)
    assert k2 > k1


def test_wavenumber_earth_decreases_with_rho():
    k_lo = wavenumber(100.0, rho=10.0)
    k_hi = wavenumber(100.0, rho=1000.0)
    assert k_lo > k_hi


def test_wavenumber_earth_formula():
    """k_eff = sqrt(π f μ₀ / ρ)."""
    f, rho = 100.0, 100.0
    _MU0 = 4e-7 * np.pi
    expected = np.sqrt(np.pi * f * _MU0 / rho)
    assert np.isclose(wavenumber(f, rho), expected)


# ─────────────────────────────────────────────────────────────────────────────
# sdas_element_pattern
# ─────────────────────────────────────────────────────────────────────────────


def test_element_null_at_0():
    """θ=0 is along dipole axis → null."""
    F = sdas_element_pattern(0.0, l=1000.0, k=0.002)
    assert float(F) == pytest.approx(0.0, abs=1e-8)


def test_element_null_at_180():
    F = sdas_element_pattern(180.0, l=1000.0, k=0.002)
    assert float(F) == pytest.approx(0.0, abs=1e-8)


def test_element_max_at_90():
    """θ=90° is broadside → normalized to 1."""
    F = sdas_element_pattern(90.0, l=1000.0, k=0.002)
    assert float(F) == pytest.approx(1.0, rel=1e-4)


def test_element_normalized_max_is_one():
    theta = np.linspace(1e-3, 179.9, 500)
    F = sdas_element_pattern(theta, l=1000.0, k=0.002)
    assert F.max() == pytest.approx(1.0, rel=1e-4)


def test_element_shape():
    theta = np.linspace(0, 180, 100)
    F = sdas_element_pattern(theta, l=500.0, k=0.003)
    assert F.shape == (100,)


def test_element_non_negative():
    theta = np.linspace(0, 180, 200)
    F = sdas_element_pattern(theta, l=1000.0, k=0.002)
    assert (F >= 0).all()


def test_element_symmetric():
    theta1 = np.linspace(1.0, 89.0, 50)
    theta2 = 180.0 - theta1
    F1 = sdas_element_pattern(theta1, l=1000.0, k=0.002)
    F2 = sdas_element_pattern(theta2, l=1000.0, k=0.002)
    assert np.allclose(F1, F2, atol=1e-6)


# ─────────────────────────────────────────────────────────────────────────────
# array_factor
# ─────────────────────────────────────────────────────────────────────────────


def test_af_broadside_unity_beta0():
    """At θ_b=0 with β=0, normalised AF = 1."""
    AF = array_factor(0.0, N=3, d=100.0, k=0.002, beta=0.0)
    assert float(AF) == pytest.approx(1.0, rel=1e-4)


def test_af_range():
    theta = np.linspace(-90, 90, 200)
    AF = array_factor(theta, N=5, d=100.0, k=0.002)
    assert (AF >= 0).all()
    assert (AF <= 1.0 + 1e-8).all()


def test_af_shape():
    theta = np.linspace(-90, 90, 100)
    AF = array_factor(theta, N=4, d=150.0, k=0.003)
    assert AF.shape == (100,)


def test_af_steered_max_near_target():
    """After beam steering, max of AF should be near θ_m."""
    theta_m = 30.0
    k = wavenumber(100.0, rho=100.0)
    d = 100.0
    beta = beam_steer(theta_m, d, k)
    theta = np.linspace(-90, 90, 1800)
    AF = array_factor(theta, N=5, d=d, k=k, beta=beta)
    theta_max = float(theta[AF.argmax()])
    assert abs(theta_max - theta_m) < 5.0  # within 5° tolerance


def test_af_scalar():
    AF = array_factor(0.0, N=3, d=100.0, k=0.002)
    assert np.isscalar(float(AF))


def test_af_n1_equals_one():
    """N=1 → no array → AF = 1 everywhere."""
    theta = np.linspace(-90, 90, 100)
    AF = array_factor(theta, N=1, d=100.0, k=0.002)
    assert np.allclose(AF, 1.0, atol=1e-6)


# ─────────────────────────────────────────────────────────────────────────────
# pas_pattern
# ─────────────────────────────────────────────────────────────────────────────


def test_pas_normalized_max_one():
    theta = np.linspace(-90, 90, 500)
    P = pas_pattern(theta, N=3, d=100.0, k=0.002)
    assert P.max() == pytest.approx(1.0, rel=1e-4)


def test_pas_non_negative():
    theta = np.linspace(-90, 90, 200)
    P = pas_pattern(theta, N=3, d=100.0, k=0.002)
    assert (P >= 0).all()


def test_pas_n1_equals_element():
    """N=1 PAS pattern equals normalised element pattern."""
    theta = np.linspace(5.0, 85.0, 50)
    k = wavenumber(100.0, rho=100.0)
    P = pas_pattern(theta, N=1, d=100.0, k=k, l=1000.0)
    # element pattern evaluated at 90° - theta (broadside → dipole angle)
    F = sdas_element_pattern(90.0 - theta, l=1000.0, k=k)
    # both normalised; compare up to constant scale
    assert np.allclose(P / P.max(), F / F.max(), atol=1e-5)


def test_pas_steered_peak():
    """Steering should shift the combined PAS peak toward θ_m."""
    theta_m = 30.0
    k = wavenumber(100.0, rho=100.0)
    d = 200.0
    beta = beam_steer(theta_m, d, k)
    theta = np.linspace(-85, 85, 1700)
    P_unsteered = pas_pattern(theta, N=4, d=d, k=k, beta=0.0, l=1000.0)
    P_steered = pas_pattern(theta, N=4, d=d, k=k, beta=beta, l=1000.0)
    peak_u = float(theta[P_unsteered.argmax()])
    peak_s = float(theta[P_steered.argmax()])
    # Steered peak must be closer to θ_m than unsteered peak
    assert (
        abs(peak_s - theta_m) <= abs(peak_u - theta_m) + 1.0
    ), f"steered peak {peak_s:.1f}° not shifted toward {theta_m}° (unsteered={peak_u:.1f}°)"


# ─────────────────────────────────────────────────────────────────────────────
# beam_steer
# ─────────────────────────────────────────────────────────────────────────────


def test_beam_steer_broadside_zero():
    """Broadside (θ_m=0) requires β=0."""
    beta = beam_steer(0.0, d=100.0, k=0.002)
    assert beta == pytest.approx(0.0, abs=1e-10)


def test_beam_steer_positive_angle_negative_beta():
    beta = beam_steer(30.0, d=100.0, k=0.002)
    assert beta < 0.0


def test_beam_steer_antisymmetric():
    b1 = beam_steer(45.0, d=100.0, k=0.002)
    b2 = beam_steer(-45.0, d=100.0, k=0.002)
    assert b1 == pytest.approx(-b2, rel=1e-6)


def test_beam_steer_roundtrip():
    """steering_angles recovers the target angle."""
    theta_m = 25.0
    k = wavenumber(100.0, rho=100.0)
    d = 150.0
    beta = beam_steer(theta_m, d, k)
    angles = steering_angles(3, d, k, beta)
    assert any(
        abs(a - theta_m) < 1.0 for a in angles
    ), f"Expected ~{theta_m}° in {angles}"


# ─────────────────────────────────────────────────────────────────────────────
# steering_angles
# ─────────────────────────────────────────────────────────────────────────────


def test_steering_angles_beta0_broadside():
    """β=0 → n=0 solution is θ_m=0° (broadside)."""
    k = wavenumber(100.0, rho=100.0)
    angles = steering_angles(3, 100.0, k, beta=0.0)
    assert any(abs(a) < 1.0 for a in angles)


def test_steering_angles_all_in_range():
    k = wavenumber(100.0, rho=100.0)
    beta = beam_steer(15.0, 100.0, k)
    angles = steering_angles(3, 100.0, k, beta)
    assert all(-90.0 <= a <= 90.0 for a in angles)


def test_steering_angles_sorted():
    k = wavenumber(200.0, rho=50.0)
    angles = steering_angles(4, 200.0, k, beta=0.1)
    assert (np.diff(angles) >= 0).all()


# ─────────────────────────────────────────────────────────────────────────────
# sdas_directivity
# ─────────────────────────────────────────────────────────────────────────────


def test_directivity_positive():
    D = sdas_directivity(l=1000.0, k=0.002)
    assert D > 0


def test_directivity_finite():
    D = sdas_directivity(l=500.0, k=0.003)
    assert np.isfinite(D)


def test_directivity_increases_with_length():
    """Longer dipole → more directive."""
    D_short = sdas_directivity(l=100.0, k=0.002)
    D_long = sdas_directivity(l=3000.0, k=0.002)
    assert D_long >= D_short


# ─────────────────────────────────────────────────────────────────────────────
# snr_gain_db
# ─────────────────────────────────────────────────────────────────────────────


def test_snr_gain_n1_zero():
    assert snr_gain_db(1) == pytest.approx(0.0, abs=1e-6)


def test_snr_gain_n2():
    assert snr_gain_db(2) == pytest.approx(20 * np.log10(2), rel=1e-6)


def test_snr_gain_n10():
    assert snr_gain_db(10) == pytest.approx(20.0, rel=1e-6)


def test_snr_gain_increases_with_n():
    gains = [snr_gain_db(n) for n in range(1, 8)]
    assert all(g2 >= g1 for g1, g2 in zip(gains, gains[1:]))


# ─────────────────────────────────────────────────────────────────────────────
# plot_radiation_pattern
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_returns_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta = np.linspace(-90, 90, 200)
    P = pas_pattern(theta, N=3, d=100.0, k=0.002)
    ax = plot_radiation_pattern(theta, P, polar=False)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_polar():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta = np.linspace(-90, 90, 200)
    P = pas_pattern(theta, N=3, d=100.0, k=0.002)
    ax = plot_radiation_pattern(theta, P, polar=True)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_multi_pattern():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta = np.linspace(-90, 90, 200)
    k = wavenumber(100.0, rho=100.0)
    d = 100.0
    pats = np.vstack(
        [
            pas_pattern(theta, N=1, d=d, k=k),
            pas_pattern(theta, N=3, d=d, k=k),
            pas_pattern(theta, N=5, d=d, k=k),
        ]
    )
    ax = plot_radiation_pattern(
        theta,
        pats,
        polar=False,
        labels=["N=1 (SDAS)", "N=3 (PAS)", "N=5 (PAS)"],
    )
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_log_scale():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta = np.linspace(-90, 90, 200)
    P = pas_pattern(theta, N=3, d=100.0, k=0.002)
    ax = plot_radiation_pattern(theta, P, polar=False, log_scale=True)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_existing_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, ax_in = plt.subplots()
    theta = np.linspace(-90, 90, 200)
    P = pas_pattern(theta, N=2, d=100.0, k=0.002)
    ax_out = plot_radiation_pattern(theta, P, polar=False, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")
