# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.exceptions import PhaseError, ResistivityError, ZError
from pycsamt.z.resphase import ResPhase


def _mk_stack(n=2, val=1.0 + 0.0j):
    z = np.full((n, 2, 2), val, dtype=complex)
    f = np.array([10.0, 1.0][:n], float)
    return z, f


# ─────────────────────────────────────────────────────────────────────────
# __init__ constructor path
# ─────────────────────────────────────────────────────────────────────────


def test_constructor_sets_z_zerr_and_freq_directly():
    z, f = _mk_stack()
    ze = np.full_like(z, 0.1, dtype=float)
    rp = ResPhase(z_array=z, z_err_array=ze, freq=f)
    assert rp._z is not None
    assert rp._z_err is not None
    assert np.allclose(rp.freq, f)


# ─────────────────────────────────────────────────────────────────────────
# resistivity / phase property getters/setters
# ─────────────────────────────────────────────────────────────────────────


def test_resistivity_getter_raises_when_unset():
    rp = ResPhase()
    with pytest.raises(ResistivityError):
        _ = rp.resistivity


def test_resistivity_setter_direct_assignment():
    rp = ResPhase()
    rp.resistivity = [[[1.0, 2.0], [3.0, 4.0]]]
    assert rp.resistivity.shape == (1, 2, 2)


def test_resistivity_err_setter_none_and_array():
    rp = ResPhase()
    rp.resistivity_err = None
    assert rp.resistivity_err is None
    rp.resistivity_err = [[[0.1, 0.1], [0.1, 0.1]]]
    assert rp.resistivity_err.shape == (1, 2, 2)


def test_phase_getter_raises_when_unset():
    rp = ResPhase()
    with pytest.raises(PhaseError):
        _ = rp.phase


def test_phase_setter_direct_assignment():
    rp = ResPhase()
    rp.phase = [[[0.0, 0.0], [0.0, 0.0]]]
    assert rp.phase.shape == (1, 2, 2)


def test_phase_err_setter_none_and_array():
    rp = ResPhase()
    rp.phase_err = None
    assert rp.phase_err is None
    rp.phase_err = [[[1.0, 1.0], [1.0, 1.0]]]
    assert rp.phase_err.shape == (1, 2, 2)


# ─────────────────────────────────────────────────────────────────────────
# compute_resistivity_phase: error branches
# ─────────────────────────────────────────────────────────────────────────


def test_compute_raises_when_z_and_freq_missing():
    rp = ResPhase()
    with pytest.raises(ZError, match="missing Z"):
        rp.compute_resistivity_phase()


def test_compute_raises_for_bad_z_shape():
    rp = ResPhase()
    with pytest.raises(ZError, match="must have shape"):
        rp.compute_resistivity_phase(z_array=np.zeros((2, 3)), freq=[1.0])


def test_compute_raises_for_freq_length_mismatch():
    z, _ = _mk_stack()
    rp = ResPhase()
    with pytest.raises(ZError, match="freq must be 1-D"):
        rp.compute_resistivity_phase(z_array=z, freq=[1.0])


def test_compute_raises_for_zerr_shape_mismatch():
    z, f = _mk_stack()
    rp = ResPhase()
    with pytest.raises(ZError, match="Z error must have same shape"):
        rp.compute_resistivity_phase(
            z_array=z, z_err_array=np.zeros((1, 2, 2)), freq=f
        )


# ─────────────────────────────────────────────────────────────────────────
# set_res_phase: error branches
# ─────────────────────────────────────────────────────────────────────────


def test_set_res_phase_raises_for_complex_rho():
    _, f = _mk_stack()
    rho = np.ones((2, 2, 2), dtype=complex)
    phi = np.zeros((2, 2, 2), dtype=float)
    rp = ResPhase()
    with pytest.raises(ResistivityError):
        rp.set_res_phase(rho, phi, f)


def test_set_res_phase_raises_for_complex_phase():
    _, f = _mk_stack()
    rho = np.ones((2, 2, 2), dtype=float)
    phi = np.zeros((2, 2, 2), dtype=complex)
    rp = ResPhase()
    with pytest.raises(PhaseError):
        rp.set_res_phase(rho, phi, f)


def test_set_res_phase_raises_for_freq_length_mismatch():
    rho = np.ones((2, 2, 2), float)
    phi = np.zeros_like(rho)
    rp = ResPhase()
    with pytest.raises(ZError, match="freq must be 1-D"):
        rp.set_res_phase(rho, phi, freq=[1.0])


def test_set_res_phase_raises_for_infinite_rho_or_phi():
    _, f = _mk_stack()
    rho = np.ones((2, 2, 2), float)
    rho[0, 0, 0] = np.inf
    phi = np.zeros_like(rho)
    rp = ResPhase()
    with pytest.raises(ZError, match="must be finite"):
        rp.set_res_phase(rho, phi, f)


def test_set_res_phase_raises_for_error_shape_mismatch():
    _, f = _mk_stack()
    rho = np.ones((2, 2, 2), float)
    phi = np.zeros_like(rho)
    rp = ResPhase()
    with pytest.raises(ZError, match="must match shapes"):
        rp.set_res_phase(
            rho, phi, f, res_err_array=np.zeros((1, 2, 2)), phase_err_array=phi
        )


def test_set_res_phase_raises_for_infinite_errors():
    _, f = _mk_stack()
    rho = np.ones((2, 2, 2), float)
    phi = np.zeros_like(rho)
    rho_err = np.zeros_like(rho)
    rho_err[0, 0, 0] = np.inf
    phi_err = np.zeros_like(rho)
    rp = ResPhase()
    with pytest.raises(ZError, match="must be finite"):
        rp.set_res_phase(rho, phi, f, res_err_array=rho_err, phase_err_array=phi_err)


def test_set_res_phase_raises_for_negative_errors():
    _, f = _mk_stack()
    rho = np.ones((2, 2, 2), float)
    phi = np.zeros_like(rho)
    rho_err = np.zeros_like(rho)
    rho_err[0, 0, 0] = -0.1
    phi_err = np.zeros_like(rho)
    rp = ResPhase()
    with pytest.raises(ZError, match="non-negative"):
        rp.set_res_phase(rho, phi, f, res_err_array=rho_err, phase_err_array=phi_err)


# ─────────────────────────────────────────────────────────────────────────
# Component views: full coverage of every slot
# ─────────────────────────────────────────────────────────────────────────


def test_all_resistivity_and_phase_component_views():
    z, f = _mk_stack()
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, freq=f)

    for attr in ("res_xx", "res_xy", "res_yx", "res_yy"):
        assert getattr(rp, attr).shape == (2,)
    for attr in ("phase_xx", "phase_xy", "phase_yx", "phase_yy"):
        assert getattr(rp, attr).shape == (2,)


def test_all_error_component_views_none_when_unset():
    z, f = _mk_stack()
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, freq=f)  # no z_err -> *_err None

    for attr in (
        "res_err_xx",
        "res_err_xy",
        "res_err_yx",
        "res_err_yy",
        "phase_err_xx",
        "phase_err_xy",
        "phase_err_yx",
        "phase_err_yy",
    ):
        assert getattr(rp, attr) is None


def test_all_error_component_views_present_when_set():
    z, f = _mk_stack()
    ze = np.full_like(z, 0.1, dtype=float)
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, z_err_array=ze, freq=f)

    for attr in (
        "res_err_xx",
        "res_err_xy",
        "res_err_yx",
        "res_err_yy",
        "phase_err_xx",
        "phase_err_xy",
        "phase_err_yx",
        "phase_err_yy",
    ):
        assert getattr(rp, attr).shape == (2,)


# ─────────────────────────────────────────────────────────────────────────
# Determinant-based metrics
# ─────────────────────────────────────────────────────────────────────────


def test_zdet_raises_when_z_not_set():
    rp = ResPhase()
    with pytest.raises(ZError, match="Z is not set"):
        _ = rp._zdet


def test_zdet_var_defaults_to_ones_when_no_zerr():
    z = np.zeros((2, 2, 2), complex)
    z[:, 0, 0] = 1.0
    z[:, 1, 1] = 1.0
    f = np.array([10.0, 1.0])
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, freq=f)
    assert np.allclose(rp._zdet_var, 1.0)


def test_zdet_var_uses_zerr_when_present():
    z = np.zeros((2, 2, 2), complex)
    z[:, 0, 0] = 1.0
    z[:, 1, 1] = 1.0
    ze = np.zeros((2, 2, 2), float)
    ze[:, 0, 0] = 0.1
    ze[:, 1, 1] = 0.1
    f = np.array([10.0, 1.0])
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, z_err_array=ze, freq=f)
    assert rp._zdet_var.shape == (2,)


def test_phase_det_err_and_res_det_err():
    z = np.zeros((2, 2, 2), complex)
    z[:, 0, 0] = 1.0 + 0.0j
    z[:, 1, 1] = 1.0 + 0.0j
    ze = np.zeros((2, 2, 2), float)
    ze[:, 0, 0] = 0.05
    ze[:, 1, 1] = 0.05
    f = np.array([10.0, 1.0])
    rp = ResPhase()
    rp.compute_resistivity_phase(z_array=z, z_err_array=ze, freq=f)

    pde = rp.phase_det_err
    rde = rp.res_det_err
    assert pde.shape == (2,)
    assert rde.shape == (2,)
    assert np.all(np.isfinite(pde))
    assert np.all(np.isfinite(rde))
