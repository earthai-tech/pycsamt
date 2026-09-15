# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.zonge import ops


def test_calculate_rho():
    rho = ops.calculate_rho(mag_e=100.0, mag_h=50.0, asp=100.0, freq=1.0)
    assert rho > 0
    assert np.isfinite(rho)


def test_calculate_ip():
    assert ops.calculate_ip(phz_e=10.0, phz_h=3.0) == 7.0


def test_calculate_std_dev_normal_case():
    out = ops.calculate_std_dev([1.0, 2.0, 3.0, 4.0])
    assert np.isclose(out, np.std([1.0, 2.0, 3.0, 4.0], ddof=1))


def test_calculate_std_dev_fewer_than_two_values_returns_zero():
    assert ops.calculate_std_dev([1.0]) == 0.0
    assert ops.calculate_std_dev([]) == 0.0


def test_calculate_std_dev_negative_variance_clamped_to_zero(monkeypatch):
    # The shortcut formula is mathematically non-negative, so forcing
    # the negative-variance clamp requires simulating the floating-
    # point cancellation that can occur in practice.
    monkeypatch.setattr(ops.np, "sum", lambda *a, **k: 0.0)
    out = ops.calculate_std_dev([5.0, 5.0, 5.0])
    assert out == 0.0


def test_calculate_e_field_std_dev():
    out = ops.calculate_e_field_std_dev([10.0, 20.0, 30.0], asp=100.0, current=2.0)
    assert np.isfinite(out)


def test_calculate_h_field_std_dev():
    out = ops.calculate_h_field_std_dev([1.0, 2.0, 3.0], current=2.0)
    assert np.isfinite(out)


def test_calculate_c_var_normal():
    assert ops.calculate_c_var(sigma=2.0, average=4.0) == 50.0


def test_calculate_c_var_zero_average_returns_zero():
    assert ops.calculate_c_var(sigma=2.0, average=0.0) == 0.0


def test_calculate_std_dev_rho_p():
    out = ops.calculate_std_dev_rho_p([1.0, 2.0, 3.0])
    assert out == ops.calculate_std_dev([1.0, 2.0, 3.0])


def test_calculate_std_dev_rho_c_normal():
    out = ops.calculate_std_dev_rho_c(
        rho_c=10.0, e_avg=5.0, h_avg=2.0, sigma_e=0.5, sigma_h=0.2
    )
    assert out > 0


def test_calculate_std_dev_rho_c_zero_avg_returns_zero():
    assert ops.calculate_std_dev_rho_c(10.0, 0.0, 2.0, 0.5, 0.2) == 0.0
    assert ops.calculate_std_dev_rho_c(10.0, 5.0, 0.0, 0.5, 0.2) == 0.0


def test_calculate_avg_magnitude():
    assert ops.calculate_avg_magnitude([1.0, 2.0, 3.0]) == 2.0


def test_calculate_avg_phase():
    assert ops.calculate_avg_phase([10.0, 20.0]) == 15.0


def test_calculate_parameter_avg_rho():
    assert ops.calculate_parameter_avg_rho([1.0, 3.0]) == 2.0


def test_calculate_component_avg_rho():
    out = ops.calculate_component_avg_rho(e_mag_avg=10.0, h_mag_avg=5.0, freq=1.0)
    assert out > 0


def test_calculate_magnetic_induction():
    assert ops.calculate_magnetic_induction(h_mag=10.0, rho=2.0) == 5.0


def test_calculate_apparent_resistivity_default_geometric_factor():
    out = ops.calculate_apparent_resistivity(e_mag=2.0, h_mag=1.0)
    assert out == 10.0


def test_calculate_apparent_resistivity_custom_geometric_factor():
    out = ops.calculate_apparent_resistivity(e_mag=2.0, h_mag=1.0, geometric_factor=2.0)
    assert out == 20.0


def test_calculate_snr_normal():
    out = ops.calculate_snr([1.0, 2.0, 3.0], [1.0, 1.1, 0.9])
    assert np.isfinite(out)


def test_calculate_snr_zero_noise_std_returns_inf():
    out = ops.calculate_snr([1.0, 2.0, 3.0], [5.0])
    assert out == np.inf


def test_calculate_phase_error():
    assert ops.calculate_phase_error(3.0, 10.0) == 7.0


def test_propagate_resistivity_error():
    out = ops.propagate_resistivity_error(
        rho=10.0, e_avg=5.0, h_avg=2.0, sigma_e=0.5, sigma_h=0.1
    )
    assert out > 0


def test_calculate_avg_amplitude():
    assert ops.calculate_avg_amplitude([-1.0, 2.0, -3.0]) == 2.0


def test_calculate_relative_error_normal():
    assert ops.calculate_relative_error(rho=10.0, sigma_rho=1.0) == 10.0


def test_calculate_relative_error_zero_rho_returns_zero():
    assert ops.calculate_relative_error(rho=0.0, sigma_rho=1.0) == 0.0


def test_calculate_magnitude_ratio():
    assert ops.calculate_magnitude_ratio(e_mag=10.0, h_mag=5.0) == 2.0


def test_calculate_resistivity_phase():
    out = ops.calculate_resistivity_phase(rho=1.0, phase_e=1.0, phase_h=0.0)
    assert np.isclose(out, np.arctan2(1.0, 1.0))


def test_calculate_frequency_dependent_resistivity():
    out = ops.calculate_frequency_dependent_resistivity(e_mag=10.0, h_mag=5.0, freq=2.0)
    assert out == 1.0


def test_calculate_rho_correction():
    out = ops.calculate_rho_correction(
        rho=10.0, e_std=1.0, h_std=1.0, e_avg=5.0, h_avg=5.0
    )
    assert out == 10.0 * (1 + 0.2 + 0.2)


def test_calculate_averaged_magnitude():
    assert ops.calculate_averaged_magnitude([-2.0, 2.0]) == 2.0


def test_calculate_conductivity():
    assert ops.calculate_conductivity(rho=4.0) == 0.25


def test_calculate_error_propagation_amplitude():
    out = ops.calculate_error_propagation_amplitude(
        e_std=1.0, h_std=1.0, rho_std=1.0, e_avg=10.0, h_avg=10.0, rho=10.0
    )
    assert np.isclose(out, np.sqrt(3 * (0.1**2)))


def test_calculate_e_field_error():
    out = ops.calculate_e_field_error([10.0, 20.0, 30.0], asp=100.0, current=2.0)
    assert out >= 0
