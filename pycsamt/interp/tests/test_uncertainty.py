# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp.uncertainty — UncertaintyBounds, UncertaintyResult,
MonteCarloHydro."""

from __future__ import annotations

import warnings
from unittest import mock

import numpy as np
import pytest

from pycsamt.interp._base import ResistivityModel
from pycsamt.interp.hydromodel import EMHydroModel, PetrophysicalConfig
from pycsamt.interp.petrophysics import ArchieModel, WaxmanSmitsModel
from pycsamt.interp.uncertainty import (
    MonteCarloHydro,
    UncertaintyBounds,
    UncertaintyResult,
    _draw,
    _extract_params,
    _perturb_petro,
)


def _model():
    rho = np.log10(
        np.array(
            [
                [800.0, 800.0],
                [600.0, 800.0],
                [15.0, 800.0],
                [10.0, 800.0],
                [8.0, 800.0],
            ]
        )
    )
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 500.0]),
        z_centers=np.array([5.0, 15.0, 30.0, 60.0, 100.0]),
        method="TDEM",
        rms=1.5,
    )


def _cfg():
    return PetrophysicalConfig(
        petro=ArchieModel(m=1.8, n=2.0), rho_w=20.0, porosity_prior=0.25
    )


# ─────────────────────────────────────────────────────────────────────────────
# UncertaintyBounds
# ─────────────────────────────────────────────────────────────────────────────


def test_bounds_invalid_dist_raises():
    with pytest.raises(ValueError, match="dist must be"):
        UncertaintyBounds(rho_w_range=(5.0, 80.0), dist="bogus")


def test_bounds_no_free_params_raises():
    with pytest.raises(ValueError, match="At least one"):
        UncertaintyBounds()


def test_bounds_n_free_and_free_names():
    b = UncertaintyBounds(rho_w_range=(5.0, 80.0), m_range=(1.4, 2.4))
    assert b.n_free == 2
    assert b.free_names == ["rho_w", "m"]


def test_bounds_free_names_without_rho_w():
    """rho_w_range=None exercises the False branch of that first check."""
    b = UncertaintyBounds(m_range=(1.4, 2.4))
    assert b.free_names == ["m"]
    assert b.n_free == 1


def test_bounds_n_free_all_four():
    b = UncertaintyBounds(
        rho_w_range=(5.0, 80.0),
        m_range=(1.4, 2.4),
        n_range=(1.8, 2.5),
        phi_prior_range=(0.1, 0.4),
    )
    assert b.n_free == 4
    assert b.free_names == ["rho_w", "m", "n", "phi_prior"]


def test_bounds_sample_uniform_shapes_and_types():
    b = UncertaintyBounds(rho_w_range=(5.0, 80.0), m_range=(1.4, 2.4))
    rng = np.random.default_rng(0)
    configs = b.sample(_cfg(), 10, rng)
    assert len(configs) == 10
    assert all(isinstance(c, PetrophysicalConfig) for c in configs)
    assert all(5.0 <= c.rho_w <= 80.0 for c in configs)
    assert all(1.4 <= c.petro.m <= 2.4 for c in configs)


def test_bounds_sample_normal_dist():
    b = UncertaintyBounds(rho_w_range=(20.0, 5.0), dist="normal")  # mean, std
    rng = np.random.default_rng(0)
    configs = b.sample(_cfg(), 20, rng)
    assert all(1e-3 <= c.rho_w <= 1e4 for c in configs)


def test_bounds_sample_fixed_params_stay_at_central_value():
    b = UncertaintyBounds(rho_w_range=(5.0, 80.0))  # only rho_w free
    cfg = _cfg()
    rng = np.random.default_rng(0)
    configs = b.sample(cfg, 5, rng)
    assert all(c.petro.m == pytest.approx(cfg.petro.m) for c in configs)
    assert all(c.porosity_prior == pytest.approx(cfg.porosity_prior) for c in configs)


# ─────────────────────────────────────────────────────────────────────────────
# UncertaintyResult.prob_wt_shallower_than
# ─────────────────────────────────────────────────────────────────────────────


def _bare_result(mean_wt, std_wt):
    return UncertaintyResult(
        resistivity_model=_model(),
        config=_cfg(),
        bounds=UncertaintyBounds(rho_w_range=(5.0, 80.0)),
        n_samples=10,
        mean_wt=mean_wt,
        std_wt=std_wt,
    )


def test_prob_wt_shallower_than_with_ensemble():
    result = _bare_result(np.array([30.0]), np.array([5.0]))
    ens = np.array([[10.0], [20.0], [40.0], [50.0]])
    p = result.prob_wt_shallower_than(25.0, wt_ensemble=ens)
    assert p[0] == pytest.approx(0.5)  # 2 of 4 samples < 25


def test_prob_wt_shallower_than_gaussian_approx_valid_std():
    result = _bare_result(np.array([30.0]), np.array([5.0]))
    p = result.prob_wt_shallower_than(30.0)
    assert p[0] == pytest.approx(0.5, abs=1e-6)  # at the mean


def test_prob_wt_shallower_than_zero_std_fallback():
    """std_wt == 0 -> deterministic 0/1 branch instead of the Gaussian CDF."""
    result = _bare_result(np.array([10.0, 50.0]), np.array([0.0, 0.0]))
    p = result.prob_wt_shallower_than(30.0)
    np.testing.assert_array_equal(p, [1.0, 0.0])


def test_prob_wt_shallower_than_nan_mean_stays_zero():
    result = _bare_result(np.array([np.nan]), np.array([0.0]))
    p = result.prob_wt_shallower_than(30.0)
    assert p[0] == 0.0


# ─────────────────────────────────────────────────────────────────────────────
# UncertaintyResult.station_report / to_csv
# ─────────────────────────────────────────────────────────────────────────────


def _full_result(n_x=2):
    model = _model()
    return UncertaintyResult(
        resistivity_model=model,
        config=_cfg(),
        bounds=UncertaintyBounds(rho_w_range=(5.0, 80.0)),
        n_samples=10,
        mean_wt=np.full(n_x, 30.0),
        std_wt=np.full(n_x, 5.0),
        p10_wt=np.full(n_x, 22.0),
        p50_wt=np.full(n_x, 30.0),
        p90_wt=np.full(n_x, 38.0),
        wt_detection_rate=np.full(n_x, 0.8),
        mean_T=np.full(n_x, 1e-3),
        std_T=np.full(n_x, 2e-4),
        p10_T=np.full(n_x, 5e-4),
        p50_T=np.full(n_x, 1e-3),
        p90_T=np.full(n_x, 1.5e-3),
    )


def test_station_report_contents():
    result = _full_result()
    rows = result.station_report()
    assert len(rows) == 2
    assert rows[0]["station"] == "S000"
    assert rows[0]["wt_range_m"] == pytest.approx(38.0 - 22.0)
    assert "log10_T_range" in rows[0]


def test_station_report_uses_station_names_when_present():
    model = ResistivityModel.from_array(
        np.log10(np.full((3, 2), 100.0)),
        x_centers=np.array([0.0, 100.0]),
        z_centers=np.array([5.0, 15.0, 30.0]),
        station_names=["Alpha", "Beta"],
    )
    result = UncertaintyResult(
        resistivity_model=model,
        config=_cfg(),
        bounds=UncertaintyBounds(rho_w_range=(5.0, 80.0)),
        n_samples=1,
        mean_wt=np.zeros(2),
        std_wt=np.zeros(2),
        p10_wt=np.zeros(2),
        p50_wt=np.zeros(2),
        p90_wt=np.zeros(2),
        wt_detection_rate=np.zeros(2),
        mean_T=np.zeros(2),
        std_T=np.zeros(2),
        p10_T=np.zeros(2),
        p90_T=np.zeros(2),
    )
    rows = result.station_report()
    assert rows[0]["station"] == "Alpha"
    assert rows[1]["station"] == "Beta"


def test_to_csv_writes_rows(tmp_path):
    result = _full_result()
    out = result.to_csv(tmp_path / "unc.csv")
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 1 + 2  # header + 2 stations


def test_to_csv_empty_rows_returns_without_writing(tmp_path):
    empty_model = ResistivityModel.from_array(
        np.zeros((3, 0)),
        x_centers=np.array([]),
        z_centers=np.array([5.0, 15.0, 30.0]),
        station_x=np.array([]),
        station_names=[],
    )
    result = UncertaintyResult(
        resistivity_model=empty_model,
        config=_cfg(),
        bounds=UncertaintyBounds(rho_w_range=(5.0, 80.0)),
        n_samples=1,
    )
    out = result.to_csv(tmp_path / "empty.csv")
    assert not out.exists()


# ─────────────────────────────────────────────────────────────────────────────
# MonteCarloHydro.run() / run_ensemble()
# ─────────────────────────────────────────────────────────────────────────────


def test_run_basic_shapes():
    bounds = UncertaintyBounds(rho_w_range=(15.0, 25.0), m_range=(1.6, 2.0))
    mc = MonteCarloHydro(_model(), _cfg(), bounds, n_samples=8, seed=0)
    unc = mc.run()
    assert isinstance(unc, UncertaintyResult)
    assert unc.mean_K.shape == (5, 2)
    assert unc.mean_wt.shape == (2,)
    assert unc.sampled_params.shape == (8, 2)
    assert unc.n_samples == 8


def test_run_verbose_prints_progress(capsys):
    bounds = UncertaintyBounds(rho_w_range=(15.0, 25.0))
    mc = MonteCarloHydro(_model(), _cfg(), bounds, n_samples=3, seed=0, verbose=True)
    mc.run()
    out = capsys.readouterr().out
    assert "MC sample 0/3" in out
    assert "MC complete: 3 samples" in out


def test_run_all_samples_fail_produces_all_nan_without_crash():
    bounds = UncertaintyBounds(rho_w_range=(15.0, 25.0))
    mc = MonteCarloHydro(_model(), _cfg(), bounds, n_samples=3, seed=0)
    with mock.patch.object(EMHydroModel, "fit", side_effect=RuntimeError("boom")):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            unc = mc.run()
    assert np.all(np.isnan(unc.mean_K))
    assert np.all(np.isnan(unc.mean_wt))


def test_run_ensemble_returns_raw_arrays():
    bounds = UncertaintyBounds(rho_w_range=(15.0, 25.0))
    mc = MonteCarloHydro(_model(), _cfg(), bounds, n_samples=5, seed=1)
    unc, wt_ens, T_ens = mc.run_ensemble()
    assert wt_ens.shape == (5, 2)
    assert T_ens.shape == (5, 2)
    assert isinstance(unc, UncertaintyResult)


def test_run_cv_K_is_nan_where_mean_K_nonpositive():
    bounds = UncertaintyBounds(rho_w_range=(15.0, 25.0))
    mc = MonteCarloHydro(_model(), _cfg(), bounds, n_samples=3, seed=0)
    with mock.patch.object(EMHydroModel, "fit", side_effect=RuntimeError("boom")):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            unc = mc.run()
    assert np.all(np.isnan(unc.cv_K))


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────


def test_draw_fixed_returns_central_repeated():
    rng = np.random.default_rng(0)
    vals = _draw(None, 42.0, 5, rng, "uniform", (0.0, 100.0))
    np.testing.assert_array_equal(vals, np.full(5, 42.0))


def test_draw_uniform_within_bounds():
    rng = np.random.default_rng(0)
    vals = _draw((10.0, 20.0), 15.0, 100, rng, "uniform", (0.0, 100.0))
    assert np.all(vals >= 10.0) and np.all(vals <= 20.0)


def test_draw_normal_clipped_to_physical_range():
    rng = np.random.default_rng(0)
    vals = _draw((1000.0, 500.0), 1000.0, 200, rng, "normal", (0.0, 1e4))
    assert vals.min() >= 0.0
    assert vals.max() <= 1e4


def test_perturb_petro_archie_preserves_other_fields():
    a = ArchieModel(m=1.8, n=2.0, a=0.9)
    new_a = _perturb_petro(a, m=2.1, n=1.9)
    assert new_a.m == pytest.approx(2.1)
    assert new_a.n == pytest.approx(1.9)
    assert new_a.a == pytest.approx(0.9)


def test_perturb_petro_waxman_smits_preserves_sigma_s():
    ws = WaxmanSmitsModel(m=1.8, n=2.0, sigma_s=0.02)
    new_ws = _perturb_petro(ws, m=2.0, n=2.2)
    assert new_ws.sigma_s == pytest.approx(0.02)


def test_extract_params_matches_free_names_order():
    bounds = UncertaintyBounds(
        rho_w_range=(5.0, 80.0),
        m_range=(1.4, 2.4),
        n_range=(1.8, 2.5),
        phi_prior_range=(0.1, 0.4),
    )
    cfg = PetrophysicalConfig(
        petro=ArchieModel(m=1.9, n=2.1), rho_w=33.0, porosity_prior=0.22
    )
    vals = _extract_params(cfg, bounds)
    np.testing.assert_allclose(vals, [33.0, 1.9, 2.1, 0.22])


def test_extract_params_only_rho_w_free():
    bounds = UncertaintyBounds(rho_w_range=(5.0, 80.0))
    cfg = _cfg()
    vals = _extract_params(cfg, bounds)
    assert vals.shape == (1,)
    assert vals[0] == pytest.approx(cfg.rho_w)


def test_extract_params_without_rho_w_free():
    """rho_w_range=None exercises the False branch of that first check."""
    bounds = UncertaintyBounds(m_range=(1.4, 2.4))
    cfg = _cfg()
    vals = _extract_params(cfg, bounds)
    assert vals.shape == (1,)
    assert vals[0] == pytest.approx(cfg.petro.m)


def test_run_ensemble_all_samples_fail_without_crash():
    bounds = UncertaintyBounds(rho_w_range=(15.0, 25.0))
    mc = MonteCarloHydro(_model(), _cfg(), bounds, n_samples=3, seed=0)
    with mock.patch.object(EMHydroModel, "fit", side_effect=RuntimeError("boom")):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            unc, wt_ens, T_ens = mc.run_ensemble()
    assert np.all(np.isnan(wt_ens))
    assert np.all(np.isnan(unc.mean_K))
