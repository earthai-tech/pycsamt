# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.exceptions import ZError
from pycsamt.z.z import Z


def _z2x2():
    return np.array([[0.0 + 0.0j, 0.0 + 1.0j], [0.0 - 1.0j, 0.0 + 0.0j]])


# ─────────────────────────────────────────────────────────────────────────
# __init__: z_err_array promotion at construction
# ─────────────────────────────────────────────────────────────────────────


def test_init_promotes_2d_z_err_array():
    z = np.repeat(_z2x2()[None, ...], 1, axis=0)
    obj = Z(z_array=z, z_err_array=np.ones((2, 2)) * 0.1, freq=np.array([1.0]))
    assert obj.z_err.shape == (1, 2, 2)


def test_init_accepts_3d_z_err_array():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(
        z_array=z,
        z_err_array=np.ones((2, 2, 2)) * 0.1,
        freq=np.array([10.0, 1.0]),
    )
    assert obj.z_err.shape == (2, 2, 2)


# ─────────────────────────────────────────────────────────────────────────
# freq setter branches
# ─────────────────────────────────────────────────────────────────────────


def test_freq_setter_none_clears_freq():
    obj = Z(z_array=np.repeat(_z2x2()[None, ...], 1, axis=0), freq=np.array([1.0]))
    obj.freq = None
    assert obj.freq is None


def test_freq_setter_raises_for_non_1d():
    obj = Z()
    with pytest.raises(ZError, match="1-D"):
        obj.freq = np.zeros((2, 2))


def test_freq_setter_raises_for_nonpositive():
    obj = Z()
    with pytest.raises(ZError, match="strictly positive"):
        obj.freq = np.array([1.0, -1.0])


def test_freq_setter_repeats_float_rotation_angle_when_z_already_set():
    obj = Z(z_array=np.repeat(_z2x2()[None, ...], 2, axis=0))
    # Force rotation_angle back to a bare float, as it would be before
    # any Z/freq interaction had a chance to promote it to an array.
    obj.rotation_angle = 0.0
    obj.freq = np.array([10.0, 1.0])
    assert isinstance(obj.rotation_angle, np.ndarray)
    assert obj.rotation_angle.shape == (2,)


# ─────────────────────────────────────────────────────────────────────────
# z_err setter: clearing recomputes rho/phi
# ─────────────────────────────────────────────────────────────────────────


def test_z_err_setter_none_recomputes_resphase_without_errors():
    z = np.repeat(_z2x2()[None, ...], 1, axis=0)
    obj = Z(z_array=z, freq=np.array([1.0]))
    obj.z_err = np.full((1, 2, 2), 0.1)
    assert obj.resistivity_err is not None
    obj.z_err = None
    assert obj.z_err is None
    assert obj.resistivity_err is None


# ─────────────────────────────────────────────────────────────────────────
# inverse
# ─────────────────────────────────────────────────────────────────────────


def test_inverse_raises_when_z_not_set():
    obj = Z()
    with pytest.raises(ZError, match="Z is not set"):
        _ = obj.inverse


def test_inverse_with_error_wraps_unexpected_exception(monkeypatch):
    import pycsamt.z.z as zmod

    z = np.array([[[0.0 + 0.0j, 1.0 + 0.0j], [-1.0 + 0.0j, 0.0 + 0.0j]]])
    obj = Z(z_array=z, freq=np.array([1.0]), z_err_array=np.full((1, 2, 2), 0.05))

    def _boom(*a, **k):
        raise RuntimeError("simulated inversion failure")

    monkeypatch.setattr(zmod, "invertmatrix_incl_errors", _boom)
    with pytest.raises(ZError, match="Failed to invert tensor"):
        _ = obj.inverse


# ─────────────────────────────────────────────────────────────────────────
# rotate
# ─────────────────────────────────────────────────────────────────────────


def test_rotate_raises_when_z_not_set():
    obj = Z()
    with pytest.raises(ZError, match="cannot rotate"):
        obj.rotate(10.0)


def test_rotate_raises_for_invalid_angle_count():
    z = np.repeat(_z2x2()[None, ...], 3, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0, 0.1]))
    with pytest.raises(ZError, match="Expected 1 angle"):
        obj.rotate([10.0, 20.0])  # neither 1 nor 3


def test_rotate_with_length_one_array_angle():
    z = np.repeat(_z2x2()[None, ...], 3, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0, 0.1]))
    obj.rotate(np.array([45.0]))  # not scalar, not list/tuple -> else branch
    assert obj.rotation_angle.shape == (3,)


def test_rotate_with_per_frequency_angles():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    obj.rotate([10.0, 20.0])
    assert np.allclose(obj.rotation_angle, [10.0, 20.0])


def test_rotate_resets_float_rotation_angle_when_manually_cleared():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    obj.rotation_angle = 0.0  # force back to float
    obj.rotate(15.0)
    assert isinstance(obj.rotation_angle, np.ndarray)


def test_rotate_propagates_z_err():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    obj.z_err = np.full((2, 2, 2), 0.05)
    obj.rotate(30.0)
    assert obj.z_err is not None
    assert obj.z_err.shape == (2, 2, 2)


# ─────────────────────────────────────────────────────────────────────────
# remove_static_shift
# ─────────────────────────────────────────────────────────────────────────


def test_remove_static_shift_raises_when_z_not_set():
    obj = Z()
    with pytest.raises(ZError, match="cannot remove static shift"):
        obj.remove_static_shift()


def test_remove_static_shift_raises_for_wrong_length_factor():
    z = np.repeat(_z2x2()[None, ...], 3, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0, 0.1]))
    with pytest.raises(ZError, match="scalar or length-3"):
        obj.remove_static_shift(reduce_res_factor_x=[1.0, 2.0])


def test_remove_static_shift_accepts_per_frequency_factor():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    S, zc = obj.remove_static_shift(
        reduce_res_factor_x=[1.0, 4.0], reduce_res_factor_y=[1.0, 1.0]
    )
    np.testing.assert_allclose(S[:, 0, 0], [1.0, 2.0])


def test_remove_static_shift_raises_for_nonpositive_factor():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    with pytest.raises(ZError, match="positive values"):
        obj.remove_static_shift(reduce_res_factor_x=-1.0)


# ─────────────────────────────────────────────────────────────────────────
# remove_distortion
# ─────────────────────────────────────────────────────────────────────────


def test_remove_distortion_raises_when_z_not_set():
    obj = Z()
    with pytest.raises(ZError, match="cannot remove distortion"):
        obj.remove_distortion(np.eye(2))


def test_remove_distortion_accepts_stacked_distortion_tensor(caplog):
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    D_stack = np.stack([np.eye(2), np.eye(2) * 2])
    D_used, Z0, Z0_err = obj.remove_distortion(D_stack)
    assert D_used.shape == (2, 2)
    np.testing.assert_allclose(D_used, np.eye(2))


def test_remove_distortion_raises_for_bad_distortion_shape():
    z = np.repeat(_z2x2()[None, ...], 1, axis=0)
    obj = Z(z_array=z, freq=np.array([1.0]))
    with pytest.raises(ZError, match="distortion_tensor must have shape"):
        obj.remove_distortion(np.zeros((3, 3)))


def test_remove_distortion_accepts_stacked_distortion_err_tensor():
    z = np.repeat(_z2x2()[None, ...], 1, axis=0)
    obj = Z(z_array=z, freq=np.array([1.0]))
    D = np.eye(2)
    D_err_stack = np.stack([np.ones((2, 2)) * 0.1, np.ones((2, 2)) * 0.2])
    D_used, Z0, Z0_err = obj.remove_distortion(D, D_err_stack)
    assert Z0_err is not None


def test_remove_distortion_raises_for_bad_distortion_err_shape():
    z = np.repeat(_z2x2()[None, ...], 1, axis=0)
    obj = Z(z_array=z, freq=np.array([1.0]))
    with pytest.raises(ZError, match="distortion_err_tensor must have shape"):
        obj.remove_distortion(np.eye(2), np.zeros((3, 3)))


def test_remove_distortion_raises_for_singular_distortion():
    z = np.repeat(_z2x2()[None, ...], 1, axis=0)
    obj = Z(z_array=z, freq=np.array([1.0]))
    singular_D = np.array([[1.0, 1.0], [1.0, 1.0]])
    with pytest.raises(ZError, match="singular"):
        obj.remove_distortion(singular_D)


def test_remove_distortion_with_z_err_propagates_errors():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    obj.z_err = np.full((2, 2, 2), 0.05)
    D = np.array([[1.2, 0.1], [0.05, 0.9]])
    D_used, Z0, Z0_err = obj.remove_distortion(D)
    assert Z0_err is not None
    assert Z0_err.shape == z.shape
    assert np.all(Z0_err >= 0.0)


def test_remove_distortion_with_distortion_err_only_propagates_errors():
    z = np.repeat(_z2x2()[None, ...], 1, axis=0)
    obj = Z(z_array=z, freq=np.array([1.0]))
    D = np.array([[1.2, 0.1], [0.05, 0.9]])
    D_err = np.full((2, 2), 0.02)
    D_used, Z0, Z0_err = obj.remove_distortion(D, D_err)
    assert Z0_err is not None


def test_remove_distortion_raises_when_zerr_shape_mismatch():
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    obj._z_err = np.ones((1, 2, 2))  # deliberately mismatched
    D = np.array([[1.2, 0.1], [0.05, 0.9]])
    with pytest.raises(ZError, match="must match 'z' shape"):
        obj.remove_distortion(D)


# ─────────────────────────────────────────────────────────────────────────
# Guard rails: every property that requires Z to be set
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "attr",
    [
        "only_1d",
        "only_2d",
        "trace",
        "skew",
        "det",
        "norm",
        "invariants",
        "z_xx",
        "z_xy",
        "z_yx",
        "z_yy",
    ],
)
def test_properties_raise_zerror_when_z_not_set(attr):
    obj = Z()
    with pytest.raises(ZError):
        getattr(obj, attr)


@pytest.mark.parametrize(
    "attr", ["trace_err", "skew_err", "det_err", "norm_err"]
)
def test_err_properties_return_none_when_zerr_not_set(attr):
    z = np.repeat(_z2x2()[None, ...], 2, axis=0)
    obj = Z(z_array=z, freq=np.array([10.0, 1.0]))
    assert getattr(obj, attr) is None


@pytest.mark.parametrize(
    "attr", ["z_err_xx", "z_err_xy", "z_err_yx", "z_err_yy"]
)
def test_z_err_component_views_none_when_unset(attr):
    obj = Z()
    assert getattr(obj, attr) is None
