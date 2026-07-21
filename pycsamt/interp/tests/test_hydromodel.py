# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp.hydromodel — PetrophysicalConfig, EMHydroResult,
EMHydroModel."""

from __future__ import annotations

import sys
from unittest import mock

import numpy as np
import pytest

from pycsamt.interp._base import ResistivityModel
from pycsamt.interp.hydromodel import (
    EMHydroModel,
    EMHydroResult,
    PetrophysicalConfig,
    _archie_from_ws,
    _cell_thicknesses,
    _replace_config,
    _station_name,
)
from pycsamt.interp.petrophysics import ArchieModel, WaxmanSmitsModel


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────


def _model(station_names=None):
    # station 0: resistive -> conductive (water table detectable)
    # station 1: uniformly resistive (no water table)
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
        station_names=station_names,
        method="TDEM",
        rms=1.5,
    )


def _cfg(**kwargs):
    kwargs.setdefault("petro", ArchieModel(m=1.8, n=2.0))
    kwargs.setdefault("rho_w", 20.0)
    kwargs.setdefault("porosity_prior", 0.25)
    return PetrophysicalConfig(**kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# PetrophysicalConfig
# ─────────────────────────────────────────────────────────────────────────────


def test_config_defaults():
    cfg = PetrophysicalConfig()
    assert isinstance(cfg.petro, ArchieModel)
    assert cfg.rho_w == 20.0
    assert cfg.porosity_prior == 0.25


def test_config_rejects_bad_petro_type():
    with pytest.raises(TypeError, match="petro must be"):
        PetrophysicalConfig(petro="not-a-model")


@pytest.mark.parametrize("rho_w", [0.0, -1.0])
def test_config_rejects_nonpositive_rho_w(rho_w):
    with pytest.raises(ValueError, match="rho_w"):
        PetrophysicalConfig(rho_w=rho_w)


@pytest.mark.parametrize("phi", [0.0, 1.0, -0.1, 1.5])
def test_config_rejects_invalid_porosity_prior(phi):
    with pytest.raises(ValueError, match="porosity_prior"):
        PetrophysicalConfig(porosity_prior=phi)


def test_config_accepts_waxman_smits():
    cfg = PetrophysicalConfig(petro=WaxmanSmitsModel(sigma_s=0.01))
    assert isinstance(cfg.petro, WaxmanSmitsModel)


# ─────────────────────────────────────────────────────────────────────────────
# EMHydroModel — construction / overrides
# ─────────────────────────────────────────────────────────────────────────────


def test_model_init_default_config():
    m = EMHydroModel(_model())
    assert isinstance(m.config, PetrophysicalConfig)
    assert m.method_tag == ""


def test_model_init_overrides_petro_rho_w_porosity():
    base_cfg = _cfg()
    m = EMHydroModel(
        _model(),
        base_cfg,
        petro=ArchieModel(m=2.2, n=2.0),
        rho_w=5.0,
        porosity_prior=0.4,
    )
    assert m.config.petro.m == pytest.approx(2.2)
    assert m.config.rho_w == pytest.approx(5.0)
    assert m.config.porosity_prior == pytest.approx(0.4)


def test_model_init_no_overrides_keeps_config():
    base_cfg = _cfg(rho_w=15.0)
    m = EMHydroModel(_model(), base_cfg)
    assert m.config.rho_w == pytest.approx(15.0)


# ─────────────────────────────────────────────────────────────────────────────
# EMHydroModel.fit() — end to end
# ─────────────────────────────────────────────────────────────────────────────


def test_fit_returns_result_with_expected_shapes():
    model = _model()
    m = EMHydroModel(model, _cfg(), method_tag="TDEM")
    result = m.fit()
    assert isinstance(result, EMHydroResult)
    assert result.porosity.shape == (5, 2)
    assert result.saturation.shape == (5, 2)
    assert result.hydraulic_K.shape == (5, 2)
    assert result.water_table.shape == (2,)
    assert result.transmissivity.shape == (2,)
    assert result.method_tag == "TDEM"
    assert result.tds > 0


def test_fit_detects_water_table_on_station0_not_on_station1():
    model = _model()
    result = EMHydroModel(model, _cfg()).fit()
    assert np.isfinite(result.water_table[0])
    assert np.isnan(result.water_table[1])


def test_fit_method_tag_defaults_to_model_method_when_unset():
    model = _model()
    result = EMHydroModel(model, _cfg()).fit()  # method_tag="" default
    assert result.method_tag == "TDEM"  # falls back to model.method


def test_fit_with_waxman_smits_petro():
    model = _model()
    cfg = _cfg(petro=WaxmanSmitsModel(m=1.8, n=2.0, sigma_s=0.01))
    result = EMHydroModel(model, cfg).fit()
    assert result.porosity.shape == (5, 2)
    assert np.all(np.isfinite(result.porosity))


def test_fit_with_fracture_depth_blends_fracture_K():
    model = _model()
    cfg_no_frac = _cfg(fracture_depth_m=None)
    cfg_frac = _cfg(fracture_depth_m=50.0, fracture_rho_matrix=5000.0)
    r_no_frac = EMHydroModel(model, cfg_no_frac).fit()
    r_frac = EMHydroModel(model, cfg_frac).fit()
    # deep cells (below fracture_depth_m) should differ once fracture-K blends in
    assert not np.allclose(
        r_no_frac.hydraulic_K[3:, 0], r_frac.hydraulic_K[3:, 0]
    )


def test_fit_skips_nan_and_nonpositive_rho_cells():
    rho = np.log10(np.array([[800.0], [600.0], [15.0], [10.0], [8.0]]))
    rho[1, 0] = np.nan
    model = ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0]),
        z_centers=np.array([5.0, 15.0, 30.0, 60.0, 100.0]),
        method="TDEM",
    )
    result = EMHydroModel(model, _cfg()).fit()
    # porosity_prior is left in place for the skipped (nan) cell
    assert result.porosity[1, 0] == pytest.approx(0.25)


# ─────────────────────────────────────────────────────────────────────────────
# Internal helper methods — direct unit tests for edge branches
# ─────────────────────────────────────────────────────────────────────────────


def test_compute_transmissivity_storativity_no_saturated_cells():
    model = _model()
    m = EMHydroModel(model, _cfg())
    K_map = np.full((5, 2), 1e-4)
    phi_map = np.full((5, 2), 0.25)
    h = _cell_thicknesses(model.z_centers)
    wt = np.array([1e6, 1e6])  # deeper than every cell -> sat_mask all False
    T, Sc, Su = m._compute_transmissivity_storativity(K_map, phi_map, wt, h)
    assert T[0] == 0.0
    assert Sc[0] == 0.0
    assert Su[0] == 0.0


def test_compute_dar_zarrouk_direct():
    model = _model()
    m = EMHydroModel(model, _cfg())
    h = _cell_thicknesses(model.z_centers)
    TR, S = m._compute_dar_zarrouk(h)
    assert TR.shape == (2,)
    assert S.shape == (2,)
    assert np.all(TR > 0)
    assert np.all(S > 0)


# ─────────────────────────────────────────────────────────────────────────────
# EMHydroResult — station_report / to_csv / to_dataframe
# ─────────────────────────────────────────────────────────────────────────────


def test_station_report_normal():
    model = _model()
    result = EMHydroModel(model, _cfg()).fit()
    rows = result.station_report()
    assert len(rows) == 2
    assert rows[0]["station"] in ("S000", "0")
    assert set(rows[0]) >= {
        "station",
        "x_m",
        "water_table_m",
        "mean_porosity_sat",
        "mean_K_ms",
        "transmissivity_m2s",
        "tds_mg_per_L",
    }


def _fake_result(model, cfg, *, water_table):
    n_z, n_x = model.rho_2d.shape
    return EMHydroResult(
        resistivity_model=model,
        config=cfg,
        porosity=np.full((n_z, n_x), 0.25),
        saturation=np.ones((n_z, n_x)),
        hydraulic_K=np.full((n_z, n_x), 1e-4),
        water_table=water_table,
        transmissivity=np.ones(n_x),
        storativity_confined=np.full(n_x, 0.1),
        storativity_unconfined=np.full(n_x, 0.2),
        dar_zarrouk_TR=np.full(n_x, 100.0),
        dar_zarrouk_S=np.full(n_x, 0.5),
        tds=320.0,
        method_tag="test",
    )


def test_station_report_no_saturated_cells_gives_nan_means():
    model = _model()
    cfg = _cfg()
    # water table deeper than every cell -> sat_mask all False for station 1
    result = _fake_result(model, cfg, water_table=np.array([30.0, 1e6]))
    rows = result.station_report()
    assert np.isnan(rows[1]["mean_porosity_sat"])
    assert np.isnan(rows[1]["mean_K_ms"])
    assert np.isfinite(rows[0]["mean_porosity_sat"])


def test_to_csv_writes_all_cells(tmp_path):
    model = _model()
    result = EMHydroModel(model, _cfg()).fit()
    out = result.to_csv(tmp_path / "cells.csv")
    lines = out.read_text().strip().splitlines()
    assert lines[0].split(",")[0] == "station"
    assert len(lines) == 1 + 5 * 2  # header + n_z * n_x rows


def test_station_report_csv_normal(tmp_path):
    model = _model()
    result = EMHydroModel(model, _cfg()).fit()
    out = result.station_report_csv(tmp_path / "stations.csv")
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 1 + 2  # header + 2 stations


def test_station_report_csv_empty_rows_writes_nothing(tmp_path):
    empty_model = ResistivityModel.from_array(
        np.zeros((3, 0)),
        x_centers=np.array([]),
        z_centers=np.array([5.0, 15.0, 30.0]),
        station_x=np.array([]),
        station_names=[],
        method="empty",
    )
    result = _fake_result(empty_model, _cfg(), water_table=np.array([]))
    out = result.station_report_csv(tmp_path / "empty.csv")
    # no rows -> returns early without ever creating the file
    assert not out.exists()


def test_to_dataframe_success():
    model = _model()
    result = EMHydroModel(model, _cfg()).fit()
    df = result.to_dataframe()
    assert len(df) == 2


def test_to_dataframe_without_pandas_raises():
    model = _model()
    result = EMHydroModel(model, _cfg()).fit()
    with mock.patch.dict(sys.modules, {"pandas": None}):
        with pytest.raises(ImportError, match="pandas is required"):
            result.to_dataframe()


# ─────────────────────────────────────────────────────────────────────────────
# Module-level helpers
# ─────────────────────────────────────────────────────────────────────────────


def test_cell_thicknesses_multi_cell():
    z = np.array([10.0, 30.0, 60.0])
    h = _cell_thicknesses(z)
    assert h.shape == (3,)
    assert np.all(h > 0)


def test_cell_thicknesses_single_cell_large():
    h = _cell_thicknesses(np.array([50.0]))
    assert h[0] == pytest.approx(50.0)


def test_cell_thicknesses_single_cell_small_floors_at_one():
    h = _cell_thicknesses(np.array([0.2]))
    assert h[0] == pytest.approx(1.0)


def test_station_name_uses_model_names():
    model = _model(station_names=["Alpha", "Beta"])
    assert _station_name(model, 0) == "Alpha"
    assert _station_name(model, 1) == "Beta"


def test_station_name_falls_back_when_no_names():
    model = _model(station_names=[])
    assert _station_name(model, 0) == "S000"


def test_station_name_falls_back_when_index_out_of_range():
    model = _model(station_names=["OnlyOne"])
    assert _station_name(model, 5) == "S005"


def test_archie_from_ws_preserves_params():
    ws = WaxmanSmitsModel(m=1.9, n=2.1, a=1.05, sigma_s=0.02)
    a = _archie_from_ws(ws)
    assert isinstance(a, ArchieModel)
    assert a.m == pytest.approx(1.9)
    assert a.n == pytest.approx(2.1)
    assert a.a == pytest.approx(1.05)


def test_replace_config_overrides_selected_fields():
    cfg = _cfg(rho_w=20.0)
    new_cfg = _replace_config(cfg, rho_w=5.0)
    assert new_cfg.rho_w == pytest.approx(5.0)
    assert new_cfg.porosity_prior == cfg.porosity_prior
