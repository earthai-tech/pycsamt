import numpy as np
import pytest

from pycsamt.exceptions import (
    PhaseError,
    ResistivityError,
    ZError,
)
from pycsamt.z.resphase import ResPhase


def _mk_stack(n=2, val=1.0 + 0.0j):
    z = np.full((n, 2, 2), val, dtype=complex)
    f = np.array([10.0, 1.0][:n], float)
    return z, f


def test_forward_basic_no_errors_shapes_and_values():
    z, f = _mk_stack()
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, freq=f)

    assert rp.resistivity.shape == (2, 2, 2)
    assert rp.phase.shape == (2, 2, 2)

    expect_rho = np.array([0.2 / 10.0, 0.2 / 1.0], float)
    expect_rho = expect_rho[:, None, None]
    assert np.allclose(rp.resistivity, expect_rho)

    assert np.allclose(rp.phase, 0.0)
    assert rp.resistivity_err is None
    assert rp.phase_err is None


def test_forward_with_errors_propagates_rho_and_phase():
    z, f = _mk_stack()
    e = 0.1
    ze = np.full_like(z, e, dtype=float)

    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, z_err_array=ze, freq=f)

    rho = rp.resistivity
    rho_e = rp.resistivity_err
    phi_e = rp.phase_err

    assert rho_e is not None and phi_e is not None
    # r_rel = 2 * (e / |Z|) = 0.2 → rho_err = 0.2 * rho
    assert np.allclose(rho_e, rho * 0.2, rtol=1e-6, atol=0.0)
    # phase_err ≈ atan(0.1) deg
    expect_phi = np.degrees(np.arctan(0.1))
    assert np.allclose(phi_e, expect_phi, rtol=1e-6, atol=0.0)


def test_forward_validates_shapes_and_frequency():
    z, f = _mk_stack()
    rp = ResPhase()

    with pytest.raises(ZError):
        rp.compute_resistivity_phase(z_array=z[:1], freq=f)

    with pytest.raises(ZError):
        rp.compute_resistivity_phase(z_array=z, freq=[10.0, -1.0])


def test_forward_preserves_nan_z_as_masked_values():
    z, f = _mk_stack()
    z[0, 0, 0] = np.nan
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, freq=f)

    assert np.isnan(rp.resistivity[0, 0, 0])
    assert np.isnan(rp.phase[0, 0, 0])
    assert np.isfinite(rp.resistivity[1, 0, 0])


def test_forward_allows_zerr_nan_where_z_is_masked():
    z, f = _mk_stack()
    ze = np.full_like(z, 0.1, dtype=float)
    z[0, 0, 0] = np.nan
    ze[0, 0, 0] = np.nan
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, z_err_array=ze, freq=f)

    assert np.isnan(rp.resistivity_err[0, 0, 0])
    assert np.isnan(rp.phase_err[0, 0, 0])
    assert np.isfinite(rp.resistivity_err[1, 0, 0])


def test_forward_allows_zerr_nan_where_z_is_valid():
    # A valid Z with an unavailable uncertainty (e.g. an EDI .VAR dummy
    # value at an otherwise-measured period) is legitimate real-world
    # data, not a caller error -- it must propagate as a missing error
    # estimate, not be rejected outright.
    z, f = _mk_stack()
    ze = np.full_like(z, 0.1, dtype=float)
    ze[0, 0, 0] = np.nan
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, z_err_array=ze, freq=f)

    assert np.isnan(rp.resistivity_err[0, 0, 0])
    assert np.isnan(rp.phase_err[0, 0, 0])
    assert np.isfinite(rp.resistivity_err[1, 0, 0])


def test_forward_rejects_negative_or_infinite_zerr_where_z_is_valid():
    z, f = _mk_stack()
    rp = ResPhase()

    ze_neg = np.full_like(z, 0.1, dtype=float)
    ze_neg[0, 0, 0] = -0.1
    with pytest.raises(ZError):
        rp.compute_resistivity_phase(z_array=z, z_err_array=ze_neg, freq=f)

    ze_inf = np.full_like(z, 0.1, dtype=float)
    ze_inf[0, 0, 0] = np.inf
    with pytest.raises(ZError):
        rp.compute_resistivity_phase(z_array=z, z_err_array=ze_inf, freq=f)


def test_inverse_builds_z_and_zerr_when_errors_given():
    _, f = _mk_stack()
    rho = np.array([0.2 / 10.0, 0.2 / 1.0])[:, None, None]
    rho = np.repeat(rho, 2, axis=1)
    rho = np.repeat(rho, 2, axis=2)
    phi = np.zeros_like(rho)

    rho_err = 0.2 * rho
    phi_err = np.full_like(rho, np.degrees(np.arctan(0.1)))

    rp = ResPhase()
    rp.set_res_phase(
        res_array=rho,
        phase_array=phi,
        freq=f,
        res_err_array=rho_err,
        phase_err_array=phi_err,
    )

    # |Z| = sqrt(5 f rho)
    z = rp._z
    z_abs = np.sqrt(5.0 * f[:, None, None] * rho)
    assert np.allclose(np.abs(z), z_abs, rtol=1e-12)

    # z_err computed and non-negative
    assert rp._z_err is not None
    assert rp._z_err.shape == z.shape
    assert np.all(rp._z_err >= 0.0)


def test_inverse_without_errors_leaves_z_err_none():
    _, f = _mk_stack()
    rho = np.array([0.2 / 10.0, 0.2 / 1.0])[:, None, None]
    rho = np.repeat(rho, 2, axis=1)
    rho = np.repeat(rho, 2, axis=2)
    phi = np.zeros_like(rho)

    rp = ResPhase()
    rp.set_res_phase(rho, phi, f)

    assert rp._z is not None
    assert rp._z_err is None


def test_component_views_and_guards():
    rp = ResPhase()
    with pytest.raises(ResistivityError):
        _ = rp.res_xx
    with pytest.raises(PhaseError):
        _ = rp.phase_xy

    z, f = _mk_stack()
    rp.compute_resistivity_phase(z_array=z, freq=f)

    assert rp.res_xx.shape == (2,)
    assert rp.res_xy.shape == (2,)
    assert rp.phase_yx.shape == (2,)
    assert rp.phase_yy.shape == (2,)


def test_determinant_metrics_shapes_and_finiteness():
    # Use identity so det = 1 → sqrt(det)=1 (phase 0)
    z = np.zeros((2, 2, 2), complex)
    z[:, 0, 0] = 1.0 + 0.0j
    z[:, 1, 1] = 1.0 + 0.0j
    f = np.array([10.0, 1.0], float)

    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, freq=f)

    ph_det = rp.phase_det
    rho_det = rp.res_det

    assert ph_det.shape == (2,)
    assert rho_det.shape == (2,)
    assert np.allclose(ph_det, 0.0)
    assert np.allclose(rho_det, np.array([0.02, 0.2]))


def test_freq_property_checks_via_forward_and_inverse():
    z, _ = _mk_stack()
    rp = ResPhase()
    with pytest.raises(ZError):
        rp.compute_resistivity_phase(z_array=z, freq=[1.0, -1.0])

    rho = np.ones((2, 2, 2), float)
    phi = np.zeros_like(rho)
    with pytest.raises(ZError):
        rp.set_res_phase(rho, phi, freq=[np.nan, 1.0])
