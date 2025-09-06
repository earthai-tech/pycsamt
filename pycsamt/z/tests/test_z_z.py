# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

import numpy as np
import pytest

from pycsamt.z.z import Z
from pycsamt.exceptions import ZError


def _z2x2():
    return np.array(
        [[0.0 + 0.0j, 0.0 + 1.0j],
         [0.0 - 1.0j, 0.0 + 0.0j]]
    )


def test_init_from_2x2_and_freq_compute_rp():
    z = _z2x2()
    f = np.array([1.0])
    obj = Z(z_array=z, freq=f)
    assert obj.resistivity.shape == (1, 2, 2)
    assert obj.phase.shape == (1, 2, 2)
    #  |Z|=1 on off-diagonals -> rho = 0.2 / f = 0.2
    np.testing.assert_allclose(
        obj.resistivity[0, 0, 1], 0.2, rtol=1e-12
    )
    np.testing.assert_allclose(
        obj.resistivity[0, 1, 0], 0.2, rtol=1e-12
    )
    # phases: angle(1j)=+90, angle(-1j)=-90
    np.testing.assert_allclose(
        obj.phase[0, 0, 1], 90.0, rtol=1e-12
    )
    np.testing.assert_allclose(
        obj.phase[0, 1, 0], -90.0, rtol=1e-12
    )


def test_len_and_freq_mismatch_raises():
    z = np.zeros((2, 2, 2), complex)
    obj = Z(z_array=z)
    with pytest.raises(ZError):
        obj.freq = np.array([1.0])  # mismatch
    # len(obj) uses n_freq from BaseEM/ResPhase
    obj.freq = np.array([10.0, 1.0])
    assert len(obj) == 2


def test_z_setter_promotion_and_shape_guard():
    obj = Z()
    obj.z = _z2x2()  # promote to (1,2,2)
    assert obj.z.shape == (1, 2, 2)
    with pytest.raises(ZError):
        obj.z = np.zeros((3, 2, 3), complex)


# def test_z_err_setter_promotion_and_mismatch():
#     z = np.zeros((2, 2, 2), complex)
#     obj = Z(z_array=z)
#     obj.z_err = np.ones((2, 2)) * 0.1  # promote
#     assert obj.z_err.shape == (1, 2, 2) or obj.z_err.shape == (2, 2, 2)
#     # enforce mismatch error
#     with pytest.raises(ZError):
#         obj.z_err = np.ones((3, 2, 2))

def test_z_err_setter_promotion_for_single_freq():
    obj = Z(z_array=np.zeros((1, 2, 2), complex))
    obj.z_err = np.ones((2, 2)) * 0.1  # promote to (1, 2, 2)
    assert obj.z_err.shape == (1, 2, 2)

def test_z_err_setter_mismatch_raises_for_multi_freq():
    obj = Z(z_array=np.zeros((2, 2, 2), complex))
    with pytest.raises(ZError):
        obj.z_err = np.ones((2, 2)) * 0.1  # cannot broadcast to (2, 2, 2)

def test_inverse_no_error_and_with_error():
    z = np.array(
        [[[0.0 + 0.0j, 1.0 + 0.0j],
          [-1.0 + 0.0j, 0.0 + 0.0j]]]
    )
    obj = Z(z_array=z, freq=np.array([1.0]))
    inv = obj.inverse
    exp = np.linalg.inv(z[0])
    np.testing.assert_allclose(inv[0], exp)
    # with errors present (propagate internally)
    obj.z_err = np.full((1, 2, 2), 0.05)
    inv2 = obj.inverse
    np.testing.assert_allclose(inv2[0], exp, rtol=1e-6)


def test_inverse_singular_raises():
    z = np.array(
        [[[1.0 + 0.0j, 0.0 + 0.0j],
          [0.0 + 0.0j, 0.0 + 0.0j]]]
    )
    obj = Z(z_array=z, freq=np.array([1.0]))
    with pytest.raises(ZError):
        _ = obj.inverse


def test_rotate_updates_history_and_360_noop():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    f = np.array([10.0, 1.0])
    obj = Z(z_array=z, freq=f)
    z0 = obj.z.copy()
    obj.rotate(360.0)
    np.testing.assert_allclose(obj.z, z0, rtol=1e-12)
    # history accumulates modulo 360
    obj.rotate(90.0)
    obj.rotate(270.0)
    np.testing.assert_allclose(
        obj.rotation_angle, np.zeros(2), atol=1e-12
    )


def test_resphase_error_propagation_from_z_err():
    z = np.repeat(_z2x2()[None, ...], 1, axis=0)
    obj = Z(z_array=z, freq=np.array([1.0]))
    obj.z_err = np.full((1, 2, 2), 0.1)
    obj.compute_resistivity_phase()
    # For |Z|=1 and dz=0.1 -> rho_rel=0.2, rho=0.2 -> 0.04
    np.testing.assert_allclose(
        obj.resistivity_err[0, 0, 1], 0.04, rtol=1e-3
    )
    # phase_err = atan(0.1) deg ~ 5.7106
    np.testing.assert_allclose(
        obj.phase_err[0, 0, 1], 5.710593, rtol=1e-3
    )


def test_remove_static_shift_identity():
    z = np.repeat(_z2x2()[None, ...], 3, axis=0)
    f = np.array([10.0, 1.0, 0.1])
    obj = Z(z_array=z, freq=f)
    S, zc = obj.remove_static_shift(1.0, 1.0)
    np.testing.assert_allclose(S[:, 0, 0], 1.0)
    np.testing.assert_allclose(S[:, 1, 1], 1.0)
    np.testing.assert_allclose(zc, z)


def test_remove_distortion_identity():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    D = np.eye(2)
    D_used, Z0, Z0_err = obj.remove_distortion(D)
    np.testing.assert_allclose(D_used, D)
    np.testing.assert_allclose(Z0, z)
    assert Z0_err is None


def test_only_1d_and_only_2d_shapes_and_values():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    z1d = obj.only_1d
    z2d = obj.only_2d
    # diagonals zero in both
    np.testing.assert_allclose(z1d[:, 0, 0], 0.0)
    np.testing.assert_allclose(z1d[:, 1, 1], 0.0)
    np.testing.assert_allclose(z2d[:, 0, 0], 0.0)
    np.testing.assert_allclose(z2d[:, 1, 1], 0.0)
    # only_1d off-diagonals share mean magnitude of originals
    m = 0.5 * (np.abs(z[:, 0, 1]) + np.abs(z[:, 1, 0]))
    np.testing.assert_allclose(np.abs(z1d[:, 0, 1]), m)
    np.testing.assert_allclose(np.abs(z1d[:, 1, 0]), m)


def test_trace_skew_det_norm_and_errors():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    # basic values
    np.testing.assert_allclose(obj.trace, np.zeros(2))
    np.testing.assert_allclose(obj.skew, z[:, 0, 1] - z[:, 1, 0])
    detz = np.linalg.det(z)
    np.testing.assert_allclose(obj.det, detz)
    np.testing.assert_allclose(obj.norm,
                               np.linalg.norm(z, axis=(1, 2)))
    # with errors
    obj.z_err = np.full((2, 2, 2), 0.05)
    assert obj.trace_err.shape == (2,)
    assert obj.skew_err.shape == (2,)
    assert obj.det_err.shape == (2,)
    ne = obj.norm_err
    assert ne is not None and ne.shape == (2,)


def test_invariants_keys_and_shapes():
    z = np.repeat(_z2x2()[None, ...], 3, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0, 0.1]))
    inv = obj.invariants
    keys = {
        "z1", "det", "det_real", "det_imag",
        "trace", "skew", "norm",
        "lambda_plus", "lambda_minus",
        "sigma_plus", "sigma_minus",
    }
    assert keys.issubset(inv.keys())
    for k in inv:
        assert inv[k].shape == (3,)


def test_component_accessors_and_err_views():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    np.testing.assert_allclose(obj.z_xx, z[:, 0, 0])
    np.testing.assert_allclose(obj.z_xy, z[:, 0, 1])
    np.testing.assert_allclose(obj.z_yx, z[:, 1, 0])
    np.testing.assert_allclose(obj.z_yy, z[:, 1, 1])
    obj.z_err = np.ones_like(z, dtype=float) * 0.1
    np.testing.assert_allclose(obj.z_err_xx, obj.z_err[:, 0, 0])
    np.testing.assert_allclose(obj.z_err_xy, obj.z_err[:, 0, 1])
    np.testing.assert_allclose(obj.z_err_yx, obj.z_err[:, 1, 0])
    np.testing.assert_allclose(obj.z_err_yy, obj.z_err[:, 1, 1])


def test_from_res_phase_roundtrip_simple():
    f = np.array([5.0])
    rho = np.ones((1, 2, 2), float)  # 1 Ohm·m
    phi = np.zeros((1, 2, 2), float)  # 0 deg
    zobj = Z.from_res_phase(rho, phi, f)
    # |Z| = sqrt(5 f rho) = sqrt(25) = 5
    np.testing.assert_allclose(np.abs(zobj.z), 5.0, rtol=1e-12)
    np.testing.assert_allclose(np.angle(zobj.z, deg=True),
                               0.0, atol=1e-12)
