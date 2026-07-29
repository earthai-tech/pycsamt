# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.utils.zmath`.

Covers scalar helpers, period-list generation, matrix algebra with
error propagation, polar/rectangular error conversions, and the
rho-phi <-> Z conversion helpers.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pycsamt.utils.zmath import (
    EmptyArrayError,
    ImpedanceConversionError,
    InvalidSignificantFiguresError,
    LogArrayError,
    MatrixInversionError,
    MatrixMultiplyError,
    NearestIndexError,
    PeriodListError,
    ReorientDataError,
    RotationMatrixError,
    centre_point,
    compute_determinant_error,
    get_period_list,
    invertmatrix_incl_errors,
    make_log_increasing_array,
    mu0,
    multiplymatrices_incl_errors,
    nearest_index,
    old_z_error2r_phi_error,
    propagate_error_polar2rect,
    propagate_error_rect2polar,
    reorient_data2D,
    rhophi2z,
    rhophi_to_z,
    rotatematrix_incl_errors,
    rotatevector_incl_errors,
    roundsf,
    z_err_to_rphi_err,
    z_error2r_phi_error,
)

# --------------------------- centre_point -----------------------------


def test_centre_point_midrange():
    x = np.array([0.0, 2.0, 10.0])
    y = np.array([-4.0, 0.0, 4.0])
    assert centre_point(x, y) == (5.0, 0.0)


def test_centre_point_empty_raises():
    with pytest.raises(EmptyArrayError):
        centre_point(np.array([]), np.array([1.0]))


# ------------------------------ roundsf -------------------------------


@pytest.mark.parametrize(
    "number, sf, expected",
    [
        (123.456, 2, 120.0),
        (0.0012345, 2, 0.0012),
        (987.0, 1, 1000.0),
        (0.0, 5, 0.0),
        (-123.456, 2, -120.0),
    ],
)
def test_roundsf_values(number, sf, expected):
    assert roundsf(number, sf) == pytest.approx(expected)


def test_roundsf_invalid_sf_raises():
    with pytest.raises(InvalidSignificantFiguresError):
        roundsf(1.23, "two")


def test_roundsf_sf_below_one_clamps():
    # sf < 1 behaves like sf = 1
    assert roundsf(987.0, 0) == roundsf(987.0, 1)


# --------------------------- get_period_list --------------------------


def test_get_period_list_exact_decades():
    periods = get_period_list(1.0, 100.0, 4)
    assert periods.size == 9
    assert periods[0] == pytest.approx(1.0)
    assert periods[-1] == pytest.approx(100.0)
    # log-uniform spacing
    ratios = periods[1:] / periods[:-1]
    assert np.allclose(ratios, ratios[0])


def test_get_period_list_outside_range_covers_bounds():
    periods = get_period_list(2.0, 50.0, 4, include_outside_range=True)
    assert periods[0] <= 2.0
    assert periods[-1] >= 50.0


def test_get_period_list_inside_range_stays_within_bounds():
    periods = get_period_list(2.0, 50.0, 4, include_outside_range=False)
    assert periods[0] >= 2.0
    assert periods[-1] <= 50.0


@pytest.mark.parametrize("pmin, pmax", [(-1.0, 10.0), (0.0, 10.0), (10.0, 1.0)])
def test_get_period_list_invalid_bounds_raise(pmin, pmax):
    with pytest.raises(PeriodListError):
        get_period_list(pmin, pmax, 4)


# ---------------------------- nearest_index ---------------------------


def test_nearest_index_basic():
    arr = np.array([0.0, 1.0, 4.0, 9.0])
    assert nearest_index(3.9, arr) == 2
    assert nearest_index(-5.0, arr) == 0


def test_nearest_index_empty_raises():
    with pytest.raises(NearestIndexError):
        nearest_index(1.0, np.array([]))


# ---------------------- make_log_increasing_array ---------------------


def test_make_log_increasing_array_contract():
    z1, target, n = 1.0, 100.0, 10
    layers = make_log_increasing_array(z1, target, n)
    assert layers.size == n
    assert layers[0] == pytest.approx(z1)
    assert np.all(np.diff(layers) > 0)
    assert layers.sum() <= target


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(z1_layer=-1.0, target_depth=10.0, n_layers=5),
        dict(z1_layer=1.0, target_depth=0.0, n_layers=5),
        dict(z1_layer=1.0, target_depth=10.0, n_layers=0),
        dict(z1_layer=1.0, target_depth=10.0, n_layers=5, increment_factor=1.5),
    ],
)
def test_make_log_increasing_array_invalid_inputs(kwargs):
    with pytest.raises(LogArrayError):
        make_log_increasing_array(**kwargs)


# ---------------------- invertmatrix_incl_errors ----------------------


def test_invertmatrix_identity_roundtrip():
    m = np.array([[2.0, 1.0], [1.0, 2.0]])
    inv, err = invertmatrix_incl_errors(m)
    assert err is None
    assert np.allclose(m @ inv, np.eye(2))


def test_invertmatrix_error_propagation_shape_and_sign():
    m = np.array([[2.0, 0.0], [0.0, 4.0]])
    e = np.full((2, 2), 0.1)
    inv, ierr = invertmatrix_incl_errors(m, e)
    assert ierr.shape == (2, 2)
    assert np.all(ierr >= 0)


def test_invertmatrix_singular_raises():
    with pytest.raises(MatrixInversionError):
        invertmatrix_incl_errors(np.array([[1.0, 2.0], [2.0, 4.0]]))


def test_invertmatrix_bad_inputs_raise():
    with pytest.raises(MatrixInversionError):
        invertmatrix_incl_errors(None)
    with pytest.raises(MatrixInversionError):
        invertmatrix_incl_errors(np.ones((2, 3)))
    with pytest.raises(MatrixInversionError):
        invertmatrix_incl_errors(np.eye(2), np.ones((3, 3)))


# ------------------------------ rhophi2z ------------------------------


def test_rhophi2z_field_units_convention():
    # |Z| = sqrt(5 * f * rho) with phase preserved
    rho = np.full((2, 2), 823.0)
    phi = np.full((2, 2), 25.0)
    z = rhophi2z(rho, phi, 500.0)
    assert z.shape == (2, 2)
    assert np.allclose(np.abs(z), math.sqrt(5 * 500 * 823))
    assert np.allclose(np.rad2deg(np.angle(z)), 25.0)


def test_rhophi2z_requires_2x2():
    with pytest.raises(ImpedanceConversionError):
        rhophi2z(np.ones(3), np.ones(3), 100.0)


# ---------------------- compute_determinant_error ---------------------


def test_compute_determinant_error_theoretical():
    z = np.array([[1 + 1j, 2 + 0j], [0 + 1j, 3 + 3j]])
    err = np.array([[0.1, 0.2], [0.3, 0.4]])
    a, b, c, d = z[0, 0], z[0, 1], z[1, 0], z[1, 1]
    expected = abs(0.1 * abs(d) + 0.4 * abs(a) - 0.2 * abs(c) - 0.3 * abs(b))
    error, sqrt_error = compute_determinant_error(z, err)
    assert error == pytest.approx(expected)
    assert sqrt_error == pytest.approx(math.sqrt(expected))


def test_compute_determinant_error_stochastic_finite():
    rng_state = np.random.get_state()
    try:
        np.random.seed(0)
        z = np.array([[1 + 1j, 2 + 0j], [0 + 1j, 3 + 3j]])
        err = np.full((2, 2), 0.05)
        error, _ = compute_determinant_error(z, err, method="stochastic", repeats=200)
        assert np.isfinite(error)
        assert error > 0
    finally:
        np.random.set_state(rng_state)


def test_compute_determinant_error_shape_mismatch():
    from pycsamt.utils.zmath import DeterminantError

    with pytest.raises(DeterminantError):
        compute_determinant_error(np.eye(2), np.ones((3, 3)))


# ---------------------- polar/rect error transport --------------------


def test_propagate_error_polar2rect_pure_radial():
    x_err, y_err = propagate_error_polar2rect(1.0, 0.1, 0.0, 0.0)
    assert x_err == pytest.approx(0.1)
    assert y_err == pytest.approx(0.0)


def test_propagate_error_rect2polar_positive():
    rho_err, phi_err = propagate_error_rect2polar(1.0, 0.1, 0.0, 0.1)
    assert rho_err > 0
    assert phi_err > 0


def test_propagate_error_rect2polar_origin_inside():
    rho_err, phi_err = propagate_error_rect2polar(0.0, 0.1, 0.0, 0.1)
    assert phi_err == pytest.approx(180.0)
    assert rho_err == pytest.approx(math.hypot(0.1, 0.1))


def test_z_error2r_phi_error_scalar_and_cap():
    res_rel, phi = z_error2r_phi_error(3.0, 4.0, 0.5)
    assert res_rel == pytest.approx(2 * 0.5 / 5.0)
    assert phi == pytest.approx(math.degrees(math.atan(0.1)))
    # error larger than |Z| -> phase capped at 90 deg
    _, phi_cap = z_error2r_phi_error(0.3, 0.4, 10.0)
    assert phi_cap == 90.0


def test_z_error2r_phi_error_array_path():
    zr = np.array([3.0, 0.3])
    zi = np.array([4.0, 0.4])
    err = np.array([0.5, 10.0])
    res_rel, phi = z_error2r_phi_error(zr, zi, err)
    assert res_rel.shape == (2,)
    assert phi[1] == 90.0


def test_old_z_error2r_phi_error_basic():
    rho_err, phi_err = old_z_error2r_phi_error(3.0, 0.1, 4.0, 0.1)
    assert rho_err > 0
    assert 0 <= phi_err <= 90


# ----------------------- rotations and products -----------------------


def test_rotatematrix_identity_invariant():
    m = np.eye(2)
    mrot, merr = rotatematrix_incl_errors(m, 37.0)
    assert merr is None
    assert np.allclose(mrot, np.eye(2))


def test_rotatematrix_90deg_flips_diagonal():
    m = np.diag([1.0, -1.0])
    mrot, _ = rotatematrix_incl_errors(m, 90.0)
    assert np.allclose(mrot, np.diag([-1.0, 1.0]), atol=1e-12)


def test_rotatematrix_error_propagation_nonnegative():
    m = np.array([[1.0, 2.0], [3.0, 4.0]])
    e = np.full((2, 2), 0.1)
    _, erot = rotatematrix_incl_errors(m, 30.0, e)
    assert erot.shape == (2, 2)
    assert np.all(erot >= 0)


def test_rotatematrix_bad_inputs_raise():
    with pytest.raises(RotationMatrixError):
        rotatematrix_incl_errors(None, 10.0)
    with pytest.raises(RotationMatrixError):
        rotatematrix_incl_errors(np.ones(3), 10.0)
    with pytest.raises(RotationMatrixError):
        rotatematrix_incl_errors(np.eye(2), 10.0, np.ones((1, 2)))


def test_rotatevector_90deg():
    v = np.array([1.0, 0.0])
    vrot, verr = rotatevector_incl_errors(v, 90.0)
    assert verr is None
    assert np.allclose(vrot, [0.0, -1.0], atol=1e-12)


def test_rotatevector_error_stays_nonnegative():
    v = np.array([1.0, 2.0])
    e = np.array([0.1, 0.2])
    _, verr = rotatevector_incl_errors(v, 45.0, e)
    assert np.all(verr >= 0)


def test_multiplymatrices_product_and_errors():
    a = np.array([[1.0, 2.0], [3.0, 4.0]])
    b = np.array([[5.0, 6.0], [7.0, 8.0]])
    prod, perr = multiplymatrices_incl_errors(a, b)
    assert np.allclose(prod, a @ b)
    assert perr is None

    ea = np.full((2, 2), 0.1)
    eb = np.full((2, 2), 0.2)
    _, perr = multiplymatrices_incl_errors(a, b, ea, eb)
    assert perr.shape == (2, 2)
    assert np.all(perr > 0)


def test_multiplymatrices_bad_inputs_raise():
    with pytest.raises(MatrixMultiplyError):
        multiplymatrices_incl_errors(None, np.eye(2))
    with pytest.raises(MatrixMultiplyError):
        multiplymatrices_incl_errors(np.eye(2), np.ones((2, 3)))
    with pytest.raises(MatrixMultiplyError):
        multiplymatrices_incl_errors(np.eye(2), np.eye(2), np.ones((3, 3)), np.eye(2))


# ---------------------------- reorient_data2D -------------------------


def test_reorient_data2d_default_axes_identity():
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0, 6.0])
    xc, yc = reorient_data2D(x, y, 0.0, 90.0)
    assert np.allclose(xc, x)
    assert np.allclose(yc, y)


def test_reorient_data2d_swapped_sensors():
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0, 6.0])
    xc, yc = reorient_data2D(x, y, 90.0, 0.0)
    assert np.allclose(xc, y)
    assert np.allclose(yc, x)


def test_reorient_data2d_trims_to_common_length():
    x = np.arange(5.0)
    y = np.arange(3.0)
    xc, yc = reorient_data2D(x, y)
    assert xc.size == yc.size == 3


def test_reorient_data2d_parallel_sensors_raise():
    with pytest.raises(ReorientDataError):
        reorient_data2D(np.ones(3), np.ones(3), 45.0, 45.0)


# ------------------------------ rhophi_to_z ---------------------------


def test_rhophi_to_z_from_resistivity_si_units():
    freq = np.array([10.0])
    rho = np.array([100.0])
    phase = np.array([0.25])  # radians
    z_abs, z_re, z_im, z_c = rhophi_to_z(phase, freq, resistivity=rho)
    expected_abs = np.sqrt(mu0 * 2 * np.pi * freq * rho)
    assert np.allclose(z_abs, expected_abs)
    assert np.allclose(np.hypot(z_re, z_im), z_abs)
    assert np.allclose(np.angle(z_c), phase)


def test_rhophi_to_z_degree_flag_and_heuristic():
    freq = np.array([10.0, 10.0])
    rho = np.array([100.0, 100.0])
    deg_phase = np.array([45.0, 60.0])
    explicit = rhophi_to_z(deg_phase, freq, resistivity=rho, deg=True)
    auto = rhophi_to_z(deg_phase, freq, resistivity=rho)
    assert np.allclose(explicit[3], auto[3])
    assert np.allclose(np.angle(explicit[3]), np.deg2rad(deg_phase))


def test_rhophi_to_z_zabs_overrides_resistivity():
    z_abs, z_re, z_im, _ = rhophi_to_z(
        np.array([0.0]),
        np.array([1.0]),
        z_abs=np.array([7.0]),
        resistivity=np.array([1.0]),
    )
    assert z_abs[0] == pytest.approx(7.0)
    assert z_re[0] == pytest.approx(7.0)
    assert z_im[0] == pytest.approx(0.0)


def test_rhophi_to_z_missing_inputs_raise():
    with pytest.raises(ImpedanceConversionError):
        rhophi_to_z(np.array([0.1]), np.array([1.0]))
    with pytest.raises(EmptyArrayError):
        rhophi_to_z(np.array([]), np.array([1.0]), resistivity=np.array([1.0]))


# ---------------------------- z_err_to_rphi_err -----------------------


def test_z_err_to_rphi_err_relationship():
    rho_rel, phase_deg = z_err_to_rphi_err(
        np.array([3.0]), np.array([4.0]), np.array([0.5])
    )
    assert rho_rel[0] == pytest.approx(0.2)
    assert phase_deg[0] == pytest.approx(math.degrees(math.atan(0.1)))


def test_z_err_to_rphi_err_cap_and_units():
    rho_rel, phase_deg = z_err_to_rphi_err(
        np.array([0.3]), np.array([0.4]), np.array([10.0])
    )
    assert rho_rel[0] >= 1.0
    assert phase_deg[0] == pytest.approx(90.0)

    _, phase_mrad = z_err_to_rphi_err(
        np.array([3.0]),
        np.array([4.0]),
        np.array([0.5]),
        deg_out=False,
    )
    assert phase_mrad[0] == pytest.approx(math.atan(0.1) * 1e3)
