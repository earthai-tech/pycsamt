import numpy as np
import pytest

from pycsamt.z.tipper import Tipper
from pycsamt.exceptions import ZError
from pycsamt.utils.zmath import rotatevector_incl_errors


def _cplx(a, b=0.0):
    return np.array(a, dtype=float) + 1j * np.array(b, dtype=float)


def test_init_shape_normalization():
    t1 = Tipper(tipper_array=np.array([1.0, 2.0]))
    assert t1.tipper.shape == (1, 1, 2)

    t2 = Tipper(tipper_array=np.ones((3, 2)))
    assert t2.tipper.shape == (3, 1, 2)

    t3 = Tipper(tipper_array=np.ones((1, 2)))
    assert t3.tipper.shape == (1, 1, 2)

    t4 = Tipper(tipper_array=np.ones((2, 1, 2)))
    assert t4.tipper.shape == (2, 1, 2)


def test_tipper_err_shapes_and_mismatch():
    t = Tipper(tipper_array=np.ones((2, 1, 2)))
    t.tipper_err = np.ones((2, 1, 2))
    assert t.tipper_err.shape == (2, 1, 2)

    with pytest.raises(ZError):
        t.tipper_err = np.ones((1, 1, 2))  # wrong n


def test_freq_validation_and_length():
    t = Tipper(tipper_array=np.ones((3, 1, 2)))
    t.freq = [10.0, 1.0, 0.1]

    with pytest.raises(ZError):
        t.freq = [10.0, -1.0, 0.1]

    with pytest.raises(ZError):
        t.freq = [1.0, 2.0]  # length mismatch


def test_compute_amp_phase_no_error():
    arr = np.zeros((2, 1, 2), dtype=complex)
    arr[0, 0, 0] = 1 + 1j
    arr[0, 0, 1] = 0 + 1j
    arr[1, 0, 0] = 1 + 0j
    arr[1, 0, 1] = 1 + 0j

    t = Tipper(tipper_array=arr)
    t.compute_amp_phase()

    amp = t.amplitude
    ph = t.phase

    exp_amp0 = np.array([np.sqrt(2.0), 1.0])[None, :]
    exp_amp1 = np.array([1.0, 1.0])[None, :]
    assert np.allclose(amp[0, 0], exp_amp0)
    assert np.allclose(amp[1, 0], exp_amp1)

    exp_ph0 = np.array([45.0, 90.0])[None, :]
    exp_ph1 = np.array([0.0, 0.0])[None, :]
    assert np.allclose(ph[0, 0], exp_ph0)
    assert np.allclose(ph[1, 0], exp_ph1)


def test_compute_amp_phase_with_error():
    arr = np.array([[[1 + 0j, 0 + 1j]]])
    err = np.array([[[0.1, 0.1]]], dtype=float)
    t = Tipper(tipper_array=arr, tipper_err_array=err)
    t.compute_amp_phase()

    assert t.amplitude_err is not None
    assert t.phase_err is not None
    assert t.amplitude_err.shape == arr.shape
    assert t.phase_err.shape == arr.shape
    assert np.all(t.amplitude_err >= 0.0)


def test_compute_mag_direction_known_angles():
    # k0: real vector = (-cos30, -sin30) -> θ=30°, |.|=1
    # k1: imag vector = (-cos60, -sin60) -> θ=60°, |.|=1
    th0 = np.deg2rad(30.0)
    th1 = np.deg2rad(60.0)

    arr = np.zeros((2, 1, 2), dtype=complex)
    arr[0, 0, 0] = -np.cos(th0) + 0j
    arr[0, 0, 1] = -np.sin(th0) + 0j
    arr[1, 0, 0] = 0.0 + 1j * (-np.cos(th1))
    arr[1, 0, 1] = 0.0 + 1j * (-np.sin(th1))

    t = Tipper(tipper_array=arr)
    t.compute_mag_direction()

    assert np.allclose(t.mag_real[0], 1.0, atol=1e-12)
    assert np.allclose(t.angle_real[0], 30.0, atol=1e-12)

    assert np.allclose(t.mag_imag[1], 1.0, atol=1e-12)
    assert np.allclose(t.angle_imag[1], 60.0, atol=1e-12)


def test_set_amp_phase_builds_complex():
    r = np.ones((2, 1, 2), dtype=float)
    phi = np.zeros_like(r)
    phi[:, 0, 1] = 90.0  # Ty at 90°
    t = Tipper()
    t.set_amp_phase(r, phi)

    T = t.tipper
    assert T.shape == (2, 1, 2)
    assert np.allclose(T[:, 0, 0].real, 1.0)
    assert np.allclose(T[:, 0, 0].imag, 0.0)
    assert np.allclose(T[:, 0, 1].real, 0.0, atol=1e-12)
    assert np.allclose(T[:, 0, 1].imag, 1.0, atol=1e-12)


def test_set_amp_phase_shape_mismatch_raises():
    r = np.ones((2, 1, 2))
    phi = np.ones((1, 1, 2))
    t = Tipper()
    with pytest.raises(ZError):
        t.set_amp_phase(r, phi)


def test_set_mag_direction_requires_initialized():
    t = Tipper()
    with pytest.raises(ZError):
        t.set_mag_direction(1.0, 0.0, 1.0, 0.0)


def test_set_mag_direction_reconstructs():
    t = Tipper(tipper_array=np.zeros((2, 1, 2), complex))
    mr = np.array([1.0, 2.0])
    ar = np.array([0.0, 90.0])
    mi = np.array([2.0, 1.0])
    ai = np.array([90.0, 0.0])

    t.set_mag_direction(mr, ar, mi, ai)

    T = t.tipper
    # k0: Tx=-1+0j, Ty=0-2j
    assert np.allclose(T[0, 0, 0], -1 + 0j)
    assert np.allclose(T[0, 0, 1], 0 - 2j)
    # k1: Tx=0-1j, Ty=-2+0j
    assert np.allclose(T[1, 0, 0], 0 - 1j)
    assert np.allclose(T[1, 0, 1], -2 + 0j)


def test_rotate_single_angle_and_history():
    arr = np.zeros((1, 1, 2), dtype=complex)
    arr[0, 0, 0] = 1.0 + 0.0j
    arr[0, 0, 1] = 0.0 + 0.0j

    err = np.full_like(arr, 0.1, dtype=float)

    t = Tipper(tipper_array=arr, tipper_err_array=err)
    t.rotate(90.0)

    exp, exp_e = rotatevector_incl_errors(arr[0], 90.0, err[0])
    assert np.allclose(t.tipper[0], exp)
    assert np.allclose(t.tipper_err[0], exp_e)
    assert np.allclose(t.rotation_angle[0], 90.0)


def test_rotate_wrong_len_raises():
    arr = np.ones((3, 1, 2), dtype=complex)
    t = Tipper(tipper_array=arr)
    with pytest.raises(ZError):
        t.rotate([0.0, 10.0])  # not 1 or n
