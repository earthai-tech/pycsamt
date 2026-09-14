# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.exceptions import ZError
from pycsamt.z.utils import (
    align_frequency_stack,
    correct_for_sensor_orientation,
    ensure_z3,
    enforce_offdiag_antisymmetry,
    finite_mask,
    freq_from_periods,
    invert_z,
    periods_from_freq,
    rho_phi_from_z,
    rotate_z,
    sigma_clip_mask,
)


def _z2():
    return np.array([[0 + 0j, 1 + 1j], [-1 - 1j, 0 + 0j]])


def _z3(n=3):
    return np.stack([_z2() * (k + 1) for k in range(n)])


# ─────────────────────────────────────────────────────────────────────────
# correct_for_sensor_orientation
# ─────────────────────────────────────────────────────────────────────────


def test_correct_for_sensor_orientation_2d_input_identity_angles():
    Z, Zerr = correct_for_sensor_orientation(_z2())
    assert Z.shape == (2, 2)
    assert Zerr is None
    assert np.allclose(Z, _z2())


def test_correct_for_sensor_orientation_3d_input():
    Z, Zerr = correct_for_sensor_orientation(_z3())
    assert Z.shape == (3, 2, 2)


def test_correct_for_sensor_orientation_raises_bad_2d_shape():
    with pytest.raises(ZError, match="shape must be"):
        correct_for_sensor_orientation(np.zeros((3, 3)))


def test_correct_for_sensor_orientation_raises_bad_ndim():
    with pytest.raises(ZError, match="must have shape"):
        correct_for_sensor_orientation(np.zeros((2, 2, 2, 2)))


def test_correct_for_sensor_orientation_casts_real_input_to_complex():
    Z, _ = correct_for_sensor_orientation(np.eye(2))
    assert np.iscomplexobj(Z)


def test_correct_for_sensor_orientation_with_err_2d():
    Zp = _z2()
    err = np.ones((2, 2))
    Z, Zerr = correct_for_sensor_orientation(Zp, z_prime_err=err)
    assert Zerr is not None
    assert Zerr.shape == (2, 2)


def test_correct_for_sensor_orientation_with_err_3d():
    Zp = _z3()
    err = np.ones((3, 2, 2))
    Z, Zerr = correct_for_sensor_orientation(Zp, z_prime_err=err)
    assert Zerr.shape == (3, 2, 2)


def test_correct_for_sensor_orientation_err_raises_bad_2d_shape():
    with pytest.raises(ZError, match="z_prime_err.*must be"):
        correct_for_sensor_orientation(_z2(), z_prime_err=np.zeros((3, 3)))


def test_correct_for_sensor_orientation_err_raises_bad_ndim():
    with pytest.raises(ZError, match="z_prime_err.*must match"):
        correct_for_sensor_orientation(_z2(), z_prime_err=np.zeros((2, 2, 2, 2)))


def test_correct_for_sensor_orientation_err_raises_mismatched_nfreq():
    with pytest.raises(ZError, match="matching"):
        correct_for_sensor_orientation(_z3(3), z_prime_err=np.ones((2, 2, 2)))


def test_correct_for_sensor_orientation_err_casts_non_float():
    Zp = _z2()
    err = np.ones((2, 2), dtype=int)
    Z, Zerr = correct_for_sensor_orientation(Zp, z_prime_err=err)
    assert Zerr.dtype == float


def test_correct_for_sensor_orientation_raises_when_u_singular():
    # bx == by (mod 180) makes U singular
    with pytest.raises(ZError, match="singular"):
        correct_for_sensor_orientation(_z2(), bx=0.0, by=0.0)


# ─────────────────────────────────────────────────────────────────────────
# periods_from_freq / freq_from_periods
# ─────────────────────────────────────────────────────────────────────────


def test_periods_from_freq_basic():
    out = periods_from_freq([1.0, 2.0, 4.0])
    assert np.allclose(out, [1.0, 0.5, 0.25])


def test_periods_from_freq_log10():
    out = periods_from_freq([10.0], log10=True)
    assert np.isclose(out[0], -1.0)


def test_periods_from_freq_raises_on_nonpositive():
    with pytest.raises(ZError, match="strictly positive"):
        periods_from_freq([1.0, 0.0])


def test_freq_from_periods_basic():
    out = freq_from_periods([1.0, 0.5])
    assert np.allclose(out, [1.0, 2.0])


def test_freq_from_periods_log10():
    out = freq_from_periods([-1.0], log10=True)
    assert np.isclose(out[0], 10.0)


def test_freq_from_periods_raises_on_nonpositive():
    with pytest.raises(ZError, match="strictly positive"):
        freq_from_periods([0.0])


# ─────────────────────────────────────────────────────────────────────────
# ensure_z3
# ─────────────────────────────────────────────────────────────────────────


def test_ensure_z3_promotes_2d():
    out = ensure_z3(_z2())
    assert out.shape == (1, 2, 2)
    assert np.iscomplexobj(out)


def test_ensure_z3_keeps_3d():
    out = ensure_z3(_z3())
    assert out.shape == (3, 2, 2)


def test_ensure_z3_raises_for_bad_shape():
    with pytest.raises(ZError, match="Z must be shape"):
        ensure_z3(np.zeros((3, 3)))


# ─────────────────────────────────────────────────────────────────────────
# rho_phi_from_z
# ─────────────────────────────────────────────────────────────────────────


def test_rho_phi_from_z_without_err():
    rho, phi, rho_e, phi_e = rho_phi_from_z(_z3(2), [1.0, 2.0])
    assert rho.shape == (2, 2, 2)
    assert phi.shape == (2, 2, 2)
    assert rho_e is None and phi_e is None


def test_rho_phi_from_z_with_err():
    # Avoid zero-magnitude entries (the sample _z3 tensor has a zero
    # diagonal) so the relative-error computation stays well-defined.
    z = np.stack(
        [
            np.array([[1 + 1j, 2 + 2j], [-2 - 2j, 1 + 1j]]),
            np.array([[2 + 2j, 4 + 4j], [-4 - 4j, 2 + 2j]]),
        ]
    )
    err = np.ones_like(z, dtype=float)
    rho, phi, rho_e, phi_e = rho_phi_from_z(z, [1.0, 2.0], z_err=err)
    assert rho_e is not None and phi_e is not None
    assert rho_e.shape == rho.shape


def test_rho_phi_from_z_raises_on_freq_length_mismatch():
    with pytest.raises(ZError, match="Length of 'freq'"):
        rho_phi_from_z(_z3(2), [1.0, 2.0, 3.0])


# ─────────────────────────────────────────────────────────────────────────
# enforce_offdiag_antisymmetry
# ─────────────────────────────────────────────────────────────────────────


def test_enforce_offdiag_antisymmetry_2d_no_err():
    z = np.array([[1 + 0j, 2 + 0j], [3 + 0j, 1 + 0j]])
    out, err = enforce_offdiag_antisymmetry(z)
    assert out.shape == (2, 2)
    assert err is None
    assert np.isclose(out[0, 1], -out[1, 0])


def test_enforce_offdiag_antisymmetry_3d_no_err():
    z = _z3(2)
    out, err = enforce_offdiag_antisymmetry(z)
    assert out.shape == (2, 2, 2)
    assert err is None


def test_enforce_offdiag_antisymmetry_2d_with_err():
    z = np.array([[1 + 0j, 2 + 0j], [3 + 0j, 1 + 0j]])
    e = np.array([[0.1, 0.2], [0.3, 0.1]])
    out, err = enforce_offdiag_antisymmetry(z, z_err=e)
    assert out.shape == (2, 2)
    assert err.shape == (2, 2)
    assert np.isclose(err[0, 1], err[1, 0])


def test_enforce_offdiag_antisymmetry_3d_with_err():
    z = _z3(2)
    e = np.ones_like(z, dtype=float)
    out, err = enforce_offdiag_antisymmetry(z, z_err=e)
    assert out.shape == (2, 2, 2)
    assert err.shape == (2, 2, 2)


# ─────────────────────────────────────────────────────────────────────────
# rotate_z
# ─────────────────────────────────────────────────────────────────────────


def test_rotate_z_2d_scalar_angle_no_err():
    out, err = rotate_z(_z2(), 30.0)
    assert out.shape == (2, 2)
    assert err is None


def test_rotate_z_3d_scalar_angle():
    out, err = rotate_z(_z3(2), 15.0)
    assert out.shape == (2, 2, 2)


def test_rotate_z_list_of_one_angle():
    out, err = rotate_z(_z3(2), [10.0])
    assert out.shape == (2, 2, 2)


def test_rotate_z_per_frequency_angles():
    out, err = rotate_z(_z3(3), [0.0, 10.0, 20.0])
    assert out.shape == (3, 2, 2)


def test_rotate_z_raises_on_angle_length_mismatch():
    with pytest.raises(ZError, match="Expected"):
        rotate_z(_z3(3), [0.0, 10.0])


def test_rotate_z_with_err_2d():
    z = _z2()
    err = np.ones((2, 2))
    out, err_out = rotate_z(z, 30.0, z_err=err)
    assert err_out is not None
    assert err_out.shape == (2, 2)


def test_rotate_z_with_err_3d():
    z = _z3(2)
    err = np.ones((2, 2, 2))
    out, err_out = rotate_z(z, 10.0, z_err=err)
    assert err_out.shape == (2, 2, 2)


# ─────────────────────────────────────────────────────────────────────────
# invert_z
# ─────────────────────────────────────────────────────────────────────────


def test_invert_z_2d_no_err():
    z = np.array([[1 + 1j, 0.5j], [0.2, 1 - 1j]])
    out, err = invert_z(z)
    assert out.shape == (2, 2)
    assert err is None
    # round trip: Z @ Z^-1 ~ I
    ident = z @ out
    assert np.allclose(ident, np.eye(2), atol=1e-6)


def test_invert_z_3d_no_err():
    z = np.stack(
        [
            np.array([[1 + 1j, 0.5j], [0.2, 1 - 1j]]),
            np.array([[2 + 0j, 0.1j], [0.1, 2 - 0j]]),
        ]
    )
    out, err = invert_z(z)
    assert out.shape == (2, 2, 2)


def test_invert_z_with_err_2d():
    z = np.array([[1 + 1j, 0.5j], [0.2, 1 - 1j]])
    e = np.ones((2, 2)) * 0.01
    out, err = invert_z(z, z_err=e)
    assert err is not None


def test_invert_z_with_err_3d():
    z = np.stack(
        [
            np.array([[1 + 1j, 0.5j], [0.2, 1 - 1j]]),
            np.array([[2 + 0j, 0.1j], [0.1, 2 - 0j]]),
        ]
    )
    e = np.ones((2, 2, 2)) * 0.01
    out, err = invert_z(z, z_err=e)
    assert err.shape == (2, 2, 2)


def test_invert_z_raises_when_singular():
    z = np.array([[1 + 0j, 1 + 0j], [1 + 0j, 1 + 0j]])  # singular
    with pytest.raises(ZError, match="Singular"):
        invert_z(z)


# ─────────────────────────────────────────────────────────────────────────
# align_frequency_stack
# ─────────────────────────────────────────────────────────────────────────


def test_align_frequency_stack_basic():
    ref = [1.0, 2.0, 3.0]
    freq = [1.0, 3.0]
    z = np.array([10.0, 30.0])
    out = align_frequency_stack(ref, freq, z)
    assert out.shape == (3,)
    assert out[0] == 10.0
    assert np.isnan(out[1])
    assert out[2] == 30.0


def test_align_frequency_stack_with_tensor_data():
    ref = [1.0, 2.0]
    freq = [2.0]
    z = np.ones((1, 2, 2), dtype=complex)
    out = align_frequency_stack(ref, freq, z, fill_value=0.0)
    assert out.shape == (2, 2, 2)
    assert np.allclose(out[0], 0.0)
    assert np.allclose(out[1], 1.0)


def test_align_frequency_stack_raises_on_shape_mismatch():
    with pytest.raises(ZError, match="must match"):
        align_frequency_stack([1.0, 2.0], [1.0, 2.0], np.zeros(3))


def test_align_frequency_stack_raises_when_freq_not_in_ref():
    with pytest.raises(ZError, match="must be present"):
        align_frequency_stack([1.0, 2.0], [5.0], np.zeros(1))


# ─────────────────────────────────────────────────────────────────────────
# finite_mask / sigma_clip_mask
# ─────────────────────────────────────────────────────────────────────────


def test_finite_mask_real_array():
    out = finite_mask(np.array([1.0, np.nan, np.inf, 2.0]))
    assert list(out) == [True, False, False, True]


def test_finite_mask_complex_array():
    out = finite_mask(np.array([1 + 1j, np.nan + 1j, 1 + np.inf * 1j]))
    assert list(out) == [True, False, False]


def test_sigma_clip_mask_default_axis():
    a = np.array([1.0, 2.0, 3.0, 100.0])
    mask = sigma_clip_mask(a, nsigma=1.0)
    assert mask[-1] == np.False_ or mask[-1] is False or not mask[-1]


def test_sigma_clip_mask_custom_axis():
    a = np.array([[1.0, 100.0], [2.0, 3.0], [3.0, 4.0]])
    mask = sigma_clip_mask(a, axis=0, nsigma=1.0)
    assert mask.shape == a.shape
