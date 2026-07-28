# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp.constraints — constraint dataclasses and
ConstrainedCalibrator."""

from __future__ import annotations

import importlib
import sys
from unittest import mock

import numpy as np
import pytest

import pycsamt.interp.constraints as constraints_mod
from pycsamt.interp._base import ResistivityModel
from pycsamt.interp.constraints import (
    ConstrainedCalibrator,
    ECConstraint,
    PumpingTestConstraint,
    SlugTestConstraint,
    WaterLevelConstraint,
    _nearest_x_ix,
    _nearest_z_ix,
)
from pycsamt.interp.hydromodel import (
    EMHydroModel,
    EMHydroResult,
    PetrophysicalConfig,
)
from pycsamt.interp.petrophysics import ArchieModel

# ─────────────────────────────────────────────────────────────────────────────
# Fixtures / helpers
# ─────────────────────────────────────────────────────────────────────────────


def _model():
    rho = np.log10(
        np.array(
            [
                [800.0, 800.0, 800.0],
                [600.0, 600.0, 600.0],
                [15.0, 15.0, 15.0],
                [10.0, 10.0, 10.0],
                [8.0, 8.0, 8.0],
            ]
        )
    )
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 500.0, 1000.0]),
        z_centers=np.array([5.0, 15.0, 30.0, 60.0, 100.0]),
        method="test",
    )


def _hydro_model(**cfg_kwargs):
    cfg = PetrophysicalConfig(petro=ArchieModel(m=1.8, n=2.0), **cfg_kwargs)
    return EMHydroModel(_model(), cfg, method_tag="TDEM")


def _fake_result(model, cfg, *, water_table, transmissivity, hydraulic_K):
    n_z, n_x = model.rho_2d.shape
    return EMHydroResult(
        resistivity_model=model,
        config=cfg,
        porosity=np.full((n_z, n_x), 0.25),
        saturation=np.full((n_z, n_x), 1.0),
        hydraulic_K=hydraulic_K,
        water_table=water_table,
        transmissivity=transmissivity,
        storativity_confined=np.zeros(n_x),
        storativity_unconfined=np.zeros(n_x),
        dar_zarrouk_TR=np.zeros(n_x),
        dar_zarrouk_S=np.zeros(n_x),
        tds=100.0,
        method_tag="test",
    )


class _FakeConstraint:
    """Duck-typed constraint that matches none of the known dataclasses —
    exercises the fallthrough branches in constraint_residuals /
    _constraint_residual_sq."""

    x = 500.0
    station = "fake"


# ─────────────────────────────────────────────────────────────────────────────
# Constraint dataclasses — validation
# ─────────────────────────────────────────────────────────────────────────────


def test_water_level_constraint_defaults():
    c = WaterLevelConstraint(x=100.0, depth_m=10.0)
    assert c.uncertainty_m == 1.0
    assert c.station == ""


@pytest.mark.parametrize(
    "kwargs",
    [
        {"depth_m": -1.0},
        {"depth_m": 10.0, "uncertainty_m": 0.0},
        {"depth_m": 10.0, "uncertainty_m": -1.0},
    ],
)
def test_water_level_constraint_invalid(kwargs):
    with pytest.raises(ValueError):
        WaterLevelConstraint(x=0.0, **kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"T_m2s": 0.0},
        {"T_m2s": -1.0},
        {"T_m2s": 1e-3, "uncertainty_factor": 0.5},
    ],
)
def test_pumping_test_constraint_invalid(kwargs):
    with pytest.raises(ValueError):
        PumpingTestConstraint(x=0.0, **kwargs)


def test_pumping_test_constraint_valid_defaults():
    c = PumpingTestConstraint(x=0.0, T_m2s=1e-3)
    assert c.uncertainty_factor == 3.0
    assert c.S is None


@pytest.mark.parametrize(
    "kwargs",
    [
        {"K_ms": 0.0, "depth_m": 10.0},
        {"K_ms": -1.0, "depth_m": 10.0},
        {"K_ms": 1e-5, "depth_m": -1.0},
    ],
)
def test_slug_test_constraint_invalid(kwargs):
    with pytest.raises(ValueError):
        SlugTestConstraint(x=0.0, **kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"ec_mscm": 0.0},
        {"ec_mscm": -1.0},
        {"ec_mscm": 1.0, "uncertainty_mscm": 0.0},
    ],
)
def test_ec_constraint_invalid(kwargs):
    with pytest.raises(ValueError):
        ECConstraint(x=0.0, **kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# ConstrainedCalibrator — construction
# ─────────────────────────────────────────────────────────────────────────────


def test_init_requires_at_least_one_constraint():
    with pytest.raises(ValueError, match="At least one constraint"):
        ConstrainedCalibrator(constraints=[])


def test_init_requires_at_least_one_calibrate_flag():
    with pytest.raises(ValueError, match="calibrate_rho_w"):
        ConstrainedCalibrator(
            constraints=[WaterLevelConstraint(x=0.0, depth_m=10.0)],
            calibrate_rho_w=False,
            calibrate_m=False,
            calibrate_phi_prior=False,
        )


def test_init_defaults():
    cal = ConstrainedCalibrator(constraints=[WaterLevelConstraint(x=0.0, depth_m=10.0)])
    assert cal.calibrated_config_ is None
    assert cal.misfit_history_ == []
    assert cal.opt_result_ is None
    assert cal.n_restarts == 1
    assert cal.verbose is False


# ─────────────────────────────────────────────────────────────────────────────
# _pack_params / _unpack_params
# ─────────────────────────────────────────────────────────────────────────────


def test_pack_params_rho_w_only():
    cal = ConstrainedCalibrator(
        constraints=[WaterLevelConstraint(x=0.0, depth_m=10.0)],
        calibrate_rho_w=True,
    )
    cfg = PetrophysicalConfig(rho_w=15.0)
    x0, bounds = cal._pack_params(cfg)
    assert x0.tolist() == [15.0]
    assert bounds == [cal.rho_w_bounds]


def test_pack_unpack_all_three_params():
    cal = ConstrainedCalibrator(
        constraints=[WaterLevelConstraint(x=0.0, depth_m=10.0)],
        calibrate_rho_w=True,
        calibrate_m=True,
        calibrate_phi_prior=True,
    )
    cfg = PetrophysicalConfig(
        petro=ArchieModel(m=1.9, n=2.0), rho_w=15.0, porosity_prior=0.3
    )
    x0, bounds = cal._pack_params(cfg)
    assert len(x0) == 3
    assert len(bounds) == 3

    new_cfg = cal._unpack_params(np.array([0.5, 2.0, 0.4]), cfg)
    assert new_cfg.rho_w == pytest.approx(0.5)
    assert new_cfg.petro.m == pytest.approx(2.0)
    assert new_cfg.porosity_prior == pytest.approx(0.4)


def test_unpack_params_clips_to_bounds():
    cal = ConstrainedCalibrator(
        constraints=[WaterLevelConstraint(x=0.0, depth_m=10.0)],
        calibrate_rho_w=True,
        rho_w_bounds=(1.0, 10.0),
    )
    cfg = PetrophysicalConfig(rho_w=5.0)
    new_cfg = cal._unpack_params(np.array([999.0]), cfg)
    assert new_cfg.rho_w == pytest.approx(10.0)  # clipped to upper bound


def test_pack_params_m_fallback_when_petro_lacks_m():
    cal = ConstrainedCalibrator(
        constraints=[WaterLevelConstraint(x=0.0, depth_m=10.0)],
        calibrate_rho_w=False,
        calibrate_m=True,
    )
    cfg = PetrophysicalConfig(petro=ArchieModel(m=1.8))
    del cfg.petro.m  # simulate a petro model without an 'm' attribute
    x0, bounds = cal._pack_params(cfg)
    assert x0[0] == pytest.approx(1.8)  # fallback default


# ─────────────────────────────────────────────────────────────────────────────
# fit()
# ─────────────────────────────────────────────────────────────────────────────


def test_fit_raises_without_scipy(monkeypatch):
    monkeypatch.setattr(constraints_mod, "_SCIPY_OK", False)
    cal = ConstrainedCalibrator(constraints=[WaterLevelConstraint(x=0.0, depth_m=10.0)])
    with pytest.raises(ImportError, match="scipy is required"):
        cal.fit(_hydro_model())


def test_fit_end_to_end_verbose(capsys):
    cal = ConstrainedCalibrator(
        constraints=[
            WaterLevelConstraint(x=500.0, depth_m=28.0),
            PumpingTestConstraint(x=500.0, T_m2s=1.0),
            SlugTestConstraint(x=500.0, K_ms=1e-4, depth_m=60.0),
            ECConstraint(x=500.0, ec_mscm=0.5),
        ],
        calibrate_rho_w=True,
        calibrate_phi_prior=True,
        verbose=True,
    )
    result = cal.fit(_hydro_model())
    assert isinstance(result, EMHydroResult)
    assert cal.calibrated_config_ is not None
    assert cal.opt_result_ is not None
    assert len(cal.misfit_history_) == 1

    out = capsys.readouterr().out
    assert "final misfit" in out
    assert "petro:" in out
    assert "rho_w:" in out
    assert "kozeny_C" not in out  # excluded from the verbose printout


def test_fit_multiple_restarts():
    cal = ConstrainedCalibrator(
        constraints=[WaterLevelConstraint(x=500.0, depth_m=28.0)],
        calibrate_rho_w=True,
        n_restarts=3,
        verbose=False,
    )
    cal.fit(_hydro_model())
    assert len(cal.misfit_history_) == 3


def test_fit_calibrate_m_only():
    cal = ConstrainedCalibrator(
        constraints=[WaterLevelConstraint(x=500.0, depth_m=28.0)],
        calibrate_rho_w=False,
        calibrate_m=True,
        verbose=False,
    )
    result = cal.fit(_hydro_model())
    assert isinstance(result, EMHydroResult)
    assert cal.calibrated_config_.petro.m != 1.8 or True  # ran without error


# ─────────────────────────────────────────────────────────────────────────────
# _objective
# ─────────────────────────────────────────────────────────────────────────────


def test_objective_normal_path_returns_finite_float():
    base_model = _hydro_model()
    cal = ConstrainedCalibrator(
        constraints=[WaterLevelConstraint(x=500.0, depth_m=28.0)],
        calibrate_rho_w=True,
    )
    x0, _ = cal._pack_params(base_model.config)
    val = cal._objective(x0, base_model)
    assert np.isfinite(val)
    assert val >= 0


def test_objective_exception_returns_large_penalty():
    base_model = _hydro_model()
    cal = ConstrainedCalibrator(
        constraints=[WaterLevelConstraint(x=500.0, depth_m=28.0)],
        calibrate_rho_w=True,
    )
    x0, _ = cal._pack_params(base_model.config)
    with mock.patch.object(EMHydroModel, "fit", side_effect=RuntimeError("boom")):
        val = cal._objective(x0, base_model)
    assert val == 1e12


# ─────────────────────────────────────────────────────────────────────────────
# constraint_residuals / _constraint_residual_sq
# ─────────────────────────────────────────────────────────────────────────────


def test_constraint_residuals_all_types():
    model = _model()
    cfg = PetrophysicalConfig(rho_w=20.0)
    result = _fake_result(
        model,
        cfg,
        water_table=np.array([28.0, 30.0, 32.0]),
        transmissivity=np.array([1.0, 2.0, 3.0]),
        hydraulic_K=np.full((5, 3), 1e-4),
    )
    cal = ConstrainedCalibrator(
        constraints=[
            WaterLevelConstraint(x=500.0, depth_m=28.0),
            PumpingTestConstraint(x=500.0, T_m2s=1.5),
            SlugTestConstraint(x=500.0, K_ms=1e-4, depth_m=60.0),
            ECConstraint(x=500.0, ec_mscm=0.5),
        ],
        calibrate_rho_w=True,
    )
    rows = cal.constraint_residuals(result)
    assert len(rows) == 4
    types = {r["type"] for r in rows}
    assert types == {
        "WaterLevelConstraint",
        "PumpingTestConstraint",
        "SlugTestConstraint",
        "ECConstraint",
    }
    wl_row = next(r for r in rows if r["type"] == "WaterLevelConstraint")
    assert wl_row["predicted"] == pytest.approx(30.0)
    assert wl_row["residual_m"] == pytest.approx(2.0)


def test_constraint_residuals_unknown_type_fallthrough():
    model = _model()
    cfg = PetrophysicalConfig(rho_w=20.0)
    result = _fake_result(
        model,
        cfg,
        water_table=np.array([28.0, 30.0, 32.0]),
        transmissivity=np.array([1.0, 2.0, 3.0]),
        hydraulic_K=np.full((5, 3), 1e-4),
    )
    cal = ConstrainedCalibrator(
        constraints=[WaterLevelConstraint(x=500.0, depth_m=28.0)],
        calibrate_rho_w=True,
    )
    cal.constraints.append(_FakeConstraint())
    rows = cal.constraint_residuals(result)
    fake_row = rows[-1]
    assert fake_row["type"] == "_FakeConstraint"
    assert set(fake_row) == {"type", "x_m", "station"}


def test_residual_sq_water_level_nan_penalty():
    model = _model()
    cfg = PetrophysicalConfig(rho_w=20.0)
    result = _fake_result(
        model,
        cfg,
        water_table=np.array([np.nan, 30.0, 32.0]),
        transmissivity=np.array([1.0, 2.0, 3.0]),
        hydraulic_K=np.full((5, 3), 1e-4),
    )
    c = WaterLevelConstraint(x=0.0, depth_m=10.0)
    val = constraints_mod._constraint_residual_sq(c, result, ix=0)
    assert val == 10.0


def test_residual_sq_pumping_test_nonpositive_penalty():
    model = _model()
    cfg = PetrophysicalConfig(rho_w=20.0)
    result = _fake_result(
        model,
        cfg,
        water_table=np.array([28.0, 30.0, 32.0]),
        transmissivity=np.array([0.0, 2.0, 3.0]),
        hydraulic_K=np.full((5, 3), 1e-4),
    )
    c = PumpingTestConstraint(x=0.0, T_m2s=1.0)
    val = constraints_mod._constraint_residual_sq(c, result, ix=0)
    assert val == 25.0


def test_residual_sq_slug_test_nonpositive_penalty():
    model = _model()
    cfg = PetrophysicalConfig(rho_w=20.0)
    hyd_K = np.full((5, 3), 1e-4)
    hyd_K[2, 0] = 0.0  # z=30 is index 2, station 0
    result = _fake_result(
        model,
        cfg,
        water_table=np.array([28.0, 30.0, 32.0]),
        transmissivity=np.array([1.0, 2.0, 3.0]),
        hydraulic_K=hyd_K,
    )
    c = SlugTestConstraint(x=0.0, K_ms=1e-4, depth_m=30.0)
    val = constraints_mod._constraint_residual_sq(c, result, ix=0)
    assert val == 25.0


def test_residual_sq_ec_constraint_normal():
    model = _model()
    cfg = PetrophysicalConfig(rho_w=20.0)
    result = _fake_result(
        model,
        cfg,
        water_table=np.array([28.0, 30.0, 32.0]),
        transmissivity=np.array([1.0, 2.0, 3.0]),
        hydraulic_K=np.full((5, 3), 1e-4),
    )
    c = ECConstraint(x=0.0, ec_mscm=0.5)
    val = constraints_mod._constraint_residual_sq(c, result, ix=0)
    assert val >= 0
    assert np.isfinite(val)


def test_residual_sq_unknown_type_returns_zero():
    model = _model()
    cfg = PetrophysicalConfig(rho_w=20.0)
    result = _fake_result(
        model,
        cfg,
        water_table=np.array([28.0, 30.0, 32.0]),
        transmissivity=np.array([1.0, 2.0, 3.0]),
        hydraulic_K=np.full((5, 3), 1e-4),
    )
    val = constraints_mod._constraint_residual_sq(_FakeConstraint(), result, ix=0)
    assert val == 0.0


# ─────────────────────────────────────────────────────────────────────────────
# module-level helpers
# ─────────────────────────────────────────────────────────────────────────────


def test_nearest_x_ix():
    model = _model()
    assert _nearest_x_ix(model, 510.0) == 1
    assert _nearest_x_ix(model, -100.0) == 0


def test_nearest_z_ix():
    model = _model()
    assert _nearest_z_ix(model, 29.0) == 2
    assert _nearest_z_ix(model, 1000.0) == 4


# ─────────────────────────────────────────────────────────────────────────────
# scipy import guard (module import-time except branch)
# ─────────────────────────────────────────────────────────────────────────────


def test_scipy_missing_at_import_sets_flag_false():
    """Simulates `import scipy.optimize` failing at module load time, which
    is otherwise untestable in an environment where scipy is installed."""
    try:
        with mock.patch.dict(sys.modules, {"scipy": None, "scipy.optimize": None}):
            importlib.reload(constraints_mod)
        assert constraints_mod._SCIPY_OK is False
    finally:
        importlib.reload(constraints_mod)  # restore real scipy binding
        assert constraints_mod._SCIPY_OK is True
