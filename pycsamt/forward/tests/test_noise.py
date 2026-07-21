# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Physics/contract tests for :mod:`pycsamt.forward.noise`.

Exercises ``GaussianNoise`` on both the MT and TEM dispatch paths,
``MultiplicativeNoise`` (log-space perturbation), ``FieldRealisticNoise``
(frequency-dependent noise floor), the ``_BaseNoiseModel.__call__``
delegation, and the ``add_gaussian_noise`` / ``add_noise`` functional
helpers.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.forward.em1d import ForwardResponse, MT1DForward, TEM1DForward
from pycsamt.forward.noise import (
    FieldRealisticNoise,
    GaussianNoise,
    MultiplicativeNoise,
    add_gaussian_noise,
    add_noise,
)
from pycsamt.forward.synthetic import LayeredModel


# ─────────────────────────────────────────────────────────────────────────────
# Shared fixtures (plain helper functions, no fixture machinery needed)
# ─────────────────────────────────────────────────────────────────────────────


def _mt_response(n=2000):
    """A clean MT response with many points for stable noise statistics."""
    freqs = np.logspace(-3, 4, n)
    model = LayeredModel(resistivity=[100.0, 10.0, 500.0], thickness=[300.0, 800.0])
    return MT1DForward(freqs).run(model)


def _tem_response():
    """A clean TEM response (real solver, not a mock)."""
    times = np.logspace(-5, -2, 40)
    model = LayeredModel(resistivity=[100.0, 10.0], thickness=[300.0])
    return TEM1DForward(times).run(model)


# ─────────────────────────────────────────────────────────────────────────────
# GaussianNoise
# ─────────────────────────────────────────────────────────────────────────────


def test_gaussian_noise_invalid_level_raises():
    with pytest.raises(ValueError):
        GaussianNoise(level=0.0)
    with pytest.raises(ValueError):
        GaussianNoise(level=1.0)
    with pytest.raises(ValueError):
        GaussianNoise(level=-0.1)


def test_gaussian_noise_mt_dispatch_and_shapes():
    resp = _mt_response(n=50)
    noisy = GaussianNoise(level=0.05).apply(resp, seed=0)

    assert noisy.method == "MT1D"
    assert noisy.rho_a.shape == resp.rho_a.shape
    assert noisy.phase.shape == resp.phase.shape
    assert noisy.z.shape == resp.z.shape
    # original response must be untouched (apply returns a new object)
    assert not np.allclose(noisy.rho_a, resp.rho_a)


def test_gaussian_noise_mt_statistics_match_level():
    resp = _mt_response(n=4000)
    level = 0.05
    noisy = GaussianNoise(level=level, apply_to="rho").apply(resp, seed=1)

    diff = np.log10(noisy.rho_a) - np.log10(resp.rho_a)
    # sample std of a N(0, level) draw over 4000 points is tight
    assert diff.std() == pytest.approx(level, rel=0.15)
    assert diff.mean() == pytest.approx(0.0, abs=0.01)
    # phase must be untouched when apply_to="rho"
    assert np.array_equal(noisy.phase, resp.phase)


def test_gaussian_noise_mt_phase_only():
    resp = _mt_response(n=4000)
    phase_level = 3.0
    noisy = GaussianNoise(
        level=0.05, apply_to="phase", phase_level=phase_level
    ).apply(resp, seed=2)

    assert np.array_equal(noisy.rho_a, resp.rho_a)
    diff = noisy.phase - resp.phase
    assert diff.std() == pytest.approx(phase_level, rel=0.15)


def test_gaussian_noise_mt_z_recomputed_consistently():
    resp = _mt_response(n=200)
    noisy = GaussianNoise(level=0.05, apply_to="both").apply(resp, seed=3)

    omega = 2.0 * np.pi * resp.freqs
    mu0 = 4.0e-7 * np.pi
    expected_z = np.sqrt(noisy.rho_a * omega * mu0) * np.exp(
        1j * np.deg2rad(noisy.phase)
    )
    assert np.allclose(noisy.z, expected_z)


def test_gaussian_noise_apply_to_z_perturbs_both_rho_and_phase():
    """``apply_to='z'`` perturbs Z by jointly perturbing rho_a and phase
    (the two quantities that determine the complex impedance) -- it is
    equivalent to ``apply_to='both'`` given the same seed."""
    resp = _mt_response(n=100)
    noisy_z = GaussianNoise(level=0.2, apply_to="z").apply(resp, seed=4)
    noisy_both = GaussianNoise(level=0.2, apply_to="both").apply(resp, seed=4)

    assert not np.array_equal(noisy_z.rho_a, resp.rho_a)
    assert not np.array_equal(noisy_z.phase, resp.phase)
    assert not np.allclose(noisy_z.z, resp.z)
    assert np.array_equal(noisy_z.rho_a, noisy_both.rho_a)
    assert np.array_equal(noisy_z.phase, noisy_both.phase)
    assert np.array_equal(noisy_z.z, noisy_both.z)


def test_gaussian_noise_mt_reproducible_with_seed():
    resp = _mt_response(n=100)
    model = GaussianNoise(level=0.05)

    n1 = model.apply(resp, seed=42)
    n2 = model.apply(resp, seed=42)
    n3 = model.apply(resp, seed=43)

    assert np.array_equal(n1.rho_a, n2.rho_a)
    assert np.array_equal(n1.phase, n2.phase)
    assert not np.array_equal(n1.rho_a, n3.rho_a)


def test_gaussian_noise_tem_dispatch_path():
    resp = _tem_response()
    assert resp.method == "TEM1D"
    assert resp.freqs is None
    assert resp.dBz_dt is not None

    noisy = GaussianNoise(level=0.05).apply(resp, seed=0)

    assert noisy.method == "TEM1D"
    assert noisy.dBz_dt.shape == resp.dBz_dt.shape
    # sign must be preserved (only magnitude is perturbed in log space)
    assert np.array_equal(np.sign(noisy.dBz_dt), np.sign(resp.dBz_dt))
    assert not np.allclose(noisy.dBz_dt, resp.dBz_dt)


def test_gaussian_noise_tem_reproducible_with_seed():
    resp = _tem_response()
    model = GaussianNoise(level=0.08)

    n1 = model.apply(resp, seed=7)
    n2 = model.apply(resp, seed=7)
    n3 = model.apply(resp, seed=8)

    assert np.array_equal(n1.dBz_dt, n2.dBz_dt)
    assert not np.array_equal(n1.dBz_dt, n3.dBz_dt)


def test_gaussian_noise_call_delegates_to_apply():
    resp = _mt_response(n=30)
    model = GaussianNoise(level=0.05)

    via_call = model(resp, seed=11)
    via_apply = model.apply(resp, seed=11)

    assert np.array_equal(via_call.rho_a, via_apply.rho_a)
    assert np.array_equal(via_call.phase, via_apply.phase)


# ─────────────────────────────────────────────────────────────────────────────
# MultiplicativeNoise
# ─────────────────────────────────────────────────────────────────────────────


def test_multiplicative_noise_mt_log_space_statistics():
    resp = _mt_response(n=4000)
    sigma = 0.07
    noisy = MultiplicativeNoise(sigma_log10=sigma).apply(resp, seed=5)

    diff = np.log10(noisy.rho_a) - np.log10(resp.rho_a)
    assert diff.std() == pytest.approx(sigma, rel=0.15)
    assert diff.mean() == pytest.approx(0.0, abs=0.01)
    # MultiplicativeNoise never touches phase/z
    assert np.array_equal(noisy.phase, resp.phase)


def test_multiplicative_noise_manual_reproduction():
    """
    apply() is a deterministic function of (response, seed): reconstruct
    the exact rng draw sequence by hand and confirm bit-for-bit agreement,
    rather than merely re-running the same code path.
    """
    resp = _mt_response(n=64)
    sigma = 0.05
    noisy = MultiplicativeNoise(sigma_log10=sigma).apply(resp, seed=123)

    rng = np.random.default_rng(123)
    log_rho = np.log10(np.maximum(resp.rho_a, 1e-12))
    expected_rho = 10.0 ** (
        log_rho + rng.normal(0.0, sigma, log_rho.shape)
    )
    assert np.allclose(noisy.rho_a, expected_rho)


def test_multiplicative_noise_tem_path():
    resp = _tem_response()
    noisy = MultiplicativeNoise(sigma_log10=0.1).apply(resp, seed=6)

    assert np.array_equal(np.sign(noisy.dBz_dt), np.sign(resp.dBz_dt))
    assert not np.allclose(noisy.dBz_dt, resp.dBz_dt)


def test_multiplicative_noise_reproducible_with_seed():
    resp = _mt_response(n=50)
    model = MultiplicativeNoise(sigma_log10=0.05)

    n1 = model.apply(resp, seed=9)
    n2 = model.apply(resp, seed=9)
    n3 = model.apply(resp, seed=10)

    assert np.array_equal(n1.rho_a, n2.rho_a)
    assert not np.array_equal(n1.rho_a, n3.rho_a)


def test_multiplicative_noise_call_delegates_to_apply():
    resp = _mt_response(n=30)
    model = MultiplicativeNoise(sigma_log10=0.05)

    via_call = model(resp, seed=1)
    via_apply = model.apply(resp, seed=1)
    assert np.array_equal(via_call.rho_a, via_apply.rho_a)


# ─────────────────────────────────────────────────────────────────────────────
# FieldRealisticNoise
# ─────────────────────────────────────────────────────────────────────────────


def test_field_realistic_noise_profile_powerline_harmonics():
    model = FieldRealisticNoise(
        base_level=0.02, powerline_freq=50.0, powerline_level=0.30,
        dead_band=False,
    )
    freqs = np.array([50.0, 100.0, 150.0, 200.0, 10.0])
    sigma = model.noise_profile(freqs)

    assert sigma.shape == freqs.shape
    assert sigma[0] == pytest.approx(0.30)       # k=1
    assert sigma[1] == pytest.approx(0.15)       # k=2
    assert sigma[2] == pytest.approx(0.10)       # k=3
    assert sigma[3] == pytest.approx(0.075)      # k=4
    assert sigma[4] == pytest.approx(0.02)       # untouched baseline


def test_field_realistic_noise_profile_dead_band():
    model = FieldRealisticNoise(
        base_level=0.02, dead_band=True, dead_band_level=0.15,
        dead_band_freq_range=(3e-4, 1e-3),
    )
    freqs = np.array([5e-4, 10.0])
    sigma = model.noise_profile(freqs)

    assert sigma[0] == pytest.approx(0.15)   # inside dead band
    assert sigma[1] == pytest.approx(0.02)   # far from any harmonic

    # disabling the dead band drops back to the base level
    model_no_db = FieldRealisticNoise(base_level=0.02, dead_band=False)
    sigma_no_db = model_no_db.noise_profile(freqs)
    assert sigma_no_db[0] == pytest.approx(0.02)


def test_field_realistic_noise_apply_perturbs_and_recomputes_z():
    resp = _mt_response(n=200)
    model = FieldRealisticNoise()
    noisy = model.apply(resp, seed=0)

    assert not np.allclose(noisy.rho_a, resp.rho_a)
    assert not np.allclose(noisy.phase, resp.phase)

    omega = 2.0 * np.pi * resp.freqs
    mu0 = 4.0e-7 * np.pi
    expected_z = np.sqrt(noisy.rho_a * omega * mu0) * np.exp(
        1j * np.deg2rad(noisy.phase)
    )
    assert np.allclose(noisy.z, expected_z)


def test_field_realistic_noise_apply_manual_reproduction():
    """Reproduce the exact rng draw sequence (rho then phase) by hand."""
    resp = _mt_response(n=64)
    model = FieldRealisticNoise()
    noisy = model.apply(resp, seed=321)

    rng = np.random.default_rng(321)
    sigma = model.noise_profile(resp.freqs)
    log_rho = np.log10(np.maximum(resp.rho_a, 1e-12))
    log_rho = log_rho + rng.normal(0.0, 1.0, sigma.shape) * sigma
    expected_rho = 10.0**log_rho
    expected_phase = resp.phase + rng.normal(0.0, 1.0, sigma.shape) * sigma * 45.0

    assert np.allclose(noisy.rho_a, expected_rho)
    assert np.allclose(noisy.phase, expected_phase)


def test_field_realistic_noise_reproducible_with_seed():
    resp = _mt_response(n=80)
    model = FieldRealisticNoise()

    n1 = model.apply(resp, seed=17)
    n2 = model.apply(resp, seed=17)
    n3 = model.apply(resp, seed=18)

    assert np.array_equal(n1.rho_a, n2.rho_a)
    assert not np.array_equal(n1.rho_a, n3.rho_a)


def test_field_realistic_noise_requires_freqs():
    resp = _tem_response()
    assert resp.freqs is None
    with pytest.raises(ValueError, match="requires freqs"):
        FieldRealisticNoise().apply(resp, seed=0)


def test_field_realistic_noise_apply_skips_missing_rho_a():
    """
    Exercise the defensive `if response.rho_a is not None` branch
    (noise.py:340) with rho_a absent: no real MT/CSAMT solver produces
    this combination (rho_a and phase are always populated together),
    but ForwardResponse is a plain dataclass so this partial state is
    constructible directly without a mock/fake response class.
    """
    freqs = np.logspace(-3, 2, 20)
    phase = np.full(freqs.shape, 45.0)
    resp = ForwardResponse(method="MT1D", freqs=freqs, rho_a=None, phase=phase)

    noisy = FieldRealisticNoise().apply(resp, seed=0)

    assert noisy.rho_a is None
    assert not np.array_equal(noisy.phase, phase)
    # z cannot be recomputed without rho_a
    assert noisy.z is None


def test_field_realistic_noise_apply_skips_missing_phase():
    """Mirror case: phase absent, rho_a present (noise.py:345, 351)."""
    freqs = np.logspace(-3, 2, 20)
    rho_a = np.full(freqs.shape, 100.0)
    resp = ForwardResponse(method="MT1D", freqs=freqs, rho_a=rho_a, phase=None)

    noisy = FieldRealisticNoise().apply(resp, seed=0)

    assert noisy.phase is None
    assert not np.array_equal(noisy.rho_a, rho_a)
    assert noisy.z is None


def test_field_realistic_noise_call_delegates_to_apply():
    resp = _mt_response(n=30)
    model = FieldRealisticNoise()

    via_call = model(resp, seed=2)
    via_apply = model.apply(resp, seed=2)
    assert np.array_equal(via_call.rho_a, via_apply.rho_a)


def test_field_realistic_noise_plot_profile_returns_axes():
    freqs = np.logspace(-4, 3, 50)
    model = FieldRealisticNoise()

    ax = model.plot_profile(freqs)
    assert ax is not None
    assert ax.get_xlabel() == "Frequency (Hz)"
    assert ax.get_ylabel() == "Noise level (%)"
    line = ax.get_lines()[0]
    sigma = model.noise_profile(freqs)
    assert np.allclose(line.get_ydata(), sigma * 100.0)
    plt.close(ax.figure)


def test_field_realistic_noise_plot_profile_reuses_given_axes():
    freqs = np.logspace(-4, 3, 20)
    model = FieldRealisticNoise()
    fig, ax = plt.subplots()

    returned = model.plot_profile(freqs, ax=ax)
    assert returned is ax
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Functional helpers: add_gaussian_noise / add_noise
# ─────────────────────────────────────────────────────────────────────────────


def test_add_gaussian_noise_matches_class_directly():
    resp = _mt_response(n=100)
    via_fn = add_gaussian_noise(resp, level=0.07, seed=99)
    via_class = GaussianNoise(level=0.07).apply(resp, seed=99)

    assert np.array_equal(via_fn.rho_a, via_class.rho_a)
    assert np.array_equal(via_fn.phase, via_class.phase)


def test_add_noise_dispatch_gaussian_default():
    resp = _mt_response(n=100)
    via_fn = add_noise(resp, seed=99)  # default noise_model="gaussian"
    via_class = GaussianNoise(level=0.05).apply(resp, seed=99)

    assert np.array_equal(via_fn.rho_a, via_class.rho_a)


@pytest.mark.parametrize("alias", ["gaussian", "gauss", "Gaussian"])
def test_add_noise_gaussian_aliases(alias):
    resp = _mt_response(n=50)
    via_fn = add_noise(resp, alias, level=0.05, seed=1)
    via_class = GaussianNoise(level=0.05).apply(resp, seed=1)
    assert np.array_equal(via_fn.rho_a, via_class.rho_a)


@pytest.mark.parametrize("alias", ["multiplicative", "mult", "log"])
def test_add_noise_multiplicative_aliases(alias):
    resp = _mt_response(n=50)
    via_fn = add_noise(resp, alias, level=0.05, seed=1)
    via_class = MultiplicativeNoise(sigma_log10=0.05).apply(resp, seed=1)
    assert np.array_equal(via_fn.rho_a, via_class.rho_a)


@pytest.mark.parametrize("alias", ["field", "realistic", "field_realistic"])
def test_add_noise_field_aliases(alias):
    resp = _mt_response(n=50)
    via_fn = add_noise(resp, alias, level=0.02, seed=1)
    via_class = FieldRealisticNoise(base_level=0.02).apply(resp, seed=1)
    assert np.array_equal(via_fn.rho_a, via_class.rho_a)


def test_add_noise_accepts_prebuilt_instance():
    resp = _mt_response(n=50)
    prebuilt = GaussianNoise(level=0.1, apply_to="phase")
    via_fn = add_noise(resp, prebuilt, seed=3)
    via_class = prebuilt.apply(resp, seed=3)
    assert np.array_equal(via_fn.phase, via_class.phase)


def test_add_noise_forwards_kwargs_to_constructor():
    resp = _mt_response(n=50)
    via_fn = add_noise(resp, "gaussian", level=0.05, apply_to="phase", seed=4)
    via_class = GaussianNoise(level=0.05, apply_to="phase").apply(resp, seed=4)
    assert np.array_equal(via_fn.rho_a, via_class.rho_a)
    assert np.array_equal(via_fn.phase, via_class.phase)


def test_add_noise_unknown_model_raises():
    resp = _mt_response(n=10)
    with pytest.raises(ValueError, match="Unknown noise_model"):
        add_noise(resp, "not-a-real-model")
