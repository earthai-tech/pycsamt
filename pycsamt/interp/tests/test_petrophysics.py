# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp.petrophysics."""

from __future__ import annotations

from unittest import mock

import numpy as np
import pytest

from pycsamt.interp import petrophysics as pp
from pycsamt.interp.petrophysics import (
    ArchieModel,
    HashinShtrikmanBounds,
    WaxmanSmitsModel,
)

# ─────────────────────────────────────────────────────────────────────────────
# ArchieModel
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "kwargs",
    [
        {"m": 0.0},
        {"m": -1.0},
        {"n": 0.0},
        {"n": -1.0},
        {"a": 0.0},
        {"a": -1.0},
    ],
)
def test_archie_invalid_params(kwargs):
    with pytest.raises(ValueError):
        ArchieModel(**kwargs)


def test_archie_formation_factor_scalar_and_array():
    a = ArchieModel(m=1.8, n=2.0, a=1.0)
    assert a.formation_factor(0.30) == pytest.approx(1.0 * 0.30 ** (-1.8))
    F = a.formation_factor(np.array([0.1, 0.3, 0.5]))
    assert F.shape == (3,)


def test_archie_forward_scalar():
    a = ArchieModel(m=1.8, n=2.0, a=1.0)
    rho = a.forward(phi=0.30, Sw=1.0, rho_w=0.025)
    assert isinstance(rho, float)
    assert rho > 0


def test_archie_saturation_inverts_forward():
    a = ArchieModel(m=1.8, n=2.0, a=1.0)
    rho = a.forward(phi=0.30, Sw=0.6, rho_w=0.025)
    Sw = a.saturation(rho, phi=0.30, rho_w=0.025)
    assert Sw == pytest.approx(0.6, rel=1e-3)


def test_archie_porosity_inverts_forward():
    a = ArchieModel(m=1.8, n=2.0, a=1.0)
    rho = a.forward(phi=0.30, Sw=1.0, rho_w=0.025)
    phi = a.porosity(rho, Sw=1.0, rho_w=0.025)
    assert phi == pytest.approx(0.30, rel=1e-3)


def test_archie_fluid_resistivity_inverts_forward():
    a = ArchieModel(m=1.8, n=2.0, a=1.0)
    rho = a.forward(phi=0.30, Sw=1.0, rho_w=0.025)
    rho_w = a.fluid_resistivity(rho, phi=0.30, Sw=1.0)
    assert rho_w == pytest.approx(0.025, rel=1e-3)


def test_archie_water_content():
    a = ArchieModel()
    theta = a.water_content(phi=0.3, Sw=0.5)
    assert theta == pytest.approx(0.15)


def test_archie_array_inputs():
    a = ArchieModel()
    rho = a.forward(phi=np.array([0.2, 0.3]), Sw=np.array([1.0, 0.8]), rho_w=0.025)
    assert rho.shape == (2,)


def test_archie_repr():
    a = ArchieModel(m=1.8, n=2.0, a=1.0)
    assert repr(a) == "ArchieModel(m=1.8, n=2.0, a=1.0)"


# ─────────────────────────────────────────────────────────────────────────────
# WaxmanSmitsModel
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "kwargs", [{"m": 0.0}, {"n": 0.0}, {"a": 0.0}, {"sigma_s": -1.0}]
)
def test_waxman_smits_invalid_params(kwargs):
    with pytest.raises(ValueError):
        WaxmanSmitsModel(**kwargs)


def test_waxman_smits_forward_scalar():
    ws = WaxmanSmitsModel(m=1.8, n=2.0, sigma_s=0.01)
    rho = ws.forward(phi=0.30, Sw=0.70, sigma_w=40.0)
    assert rho > 0


def test_waxman_smits_forward_zero_sigma_s_matches_archie_like_shape():
    ws = WaxmanSmitsModel(m=1.8, n=2.0, sigma_s=0.0)
    rho = ws.forward(phi=0.30, Sw=1.0, sigma_w=40.0)
    assert rho > 0


def test_waxman_smits_saturation_inverts_forward_bracketed():
    ws = WaxmanSmitsModel(m=1.8, n=2.0, sigma_s=0.01)
    rho = ws.forward(phi=0.30, Sw=0.6, sigma_w=40.0)
    Sw = ws.saturation(rho, phi=0.30, sigma_w=40.0)
    assert Sw == pytest.approx(0.6, abs=1e-2)


def test_waxman_smits_saturation_array_input():
    ws = WaxmanSmitsModel(sigma_s=0.01)
    rho = ws.forward(phi=np.array([0.2, 0.3]), Sw=np.array([0.5, 0.9]), sigma_w=40.0)
    Sw = ws.saturation(rho, phi=np.array([0.2, 0.3]), sigma_w=40.0)
    assert Sw.shape == (2,)
    np.testing.assert_allclose(Sw, [0.5, 0.9], atol=1e-2)


def test_waxman_smits_saturation_no_bracket_uses_initial_guess():
    """When residual(1e-4) and residual(1.0) share sign, falls back to S0."""
    ws = WaxmanSmitsModel(sigma_s=0.01)
    # A very low observed resistivity (huge sigma_obs) relative to what any
    # Sw in [1e-4, 1] can produce pushes both endpoints to the same sign.
    Sw = ws.saturation(rho=1e-2, phi=0.30, sigma_w=40.0)
    assert 0.0 <= Sw <= 1.0


def test_waxman_smits_saturation_sigma_w_zero_uses_default_initial_guess():
    """sw_si <= 1e-12 takes the S0 = 0.5 fallback branch."""
    ws = WaxmanSmitsModel(sigma_s=0.01)
    Sw = ws.saturation(rho=100.0, phi=0.30, sigma_w=0.0)
    assert 0.0 <= Sw <= 1.0


def test_waxman_smits_saturation_brentq_exception_falls_back():
    # rho=100, phi=0.5 is a verified bracketed case (fa*fb < 0), so brentq
    # is actually invoked and its failure must be caught.
    ws = WaxmanSmitsModel(sigma_s=0.01)
    with mock.patch.object(pp, "brentq", side_effect=RuntimeError("boom")):
        Sw = ws.saturation(rho=100.0, phi=0.5, sigma_w=40.0)
    assert 0.0 <= Sw <= 1.0


def test_waxman_smits_porosity_inverts_forward():
    ws = WaxmanSmitsModel(m=1.8, n=2.0, sigma_s=0.01)
    rho = ws.forward(phi=0.30, Sw=1.0, sigma_w=40.0)
    phi = ws.porosity(rho, Sw=1.0, sigma_w=40.0)
    assert phi == pytest.approx(0.30, rel=1e-2)


def test_waxman_smits_repr():
    ws = WaxmanSmitsModel(m=1.8, n=2.0, a=1.0, sigma_s=0.01)
    assert repr(ws) == "WaxmanSmitsModel(m=1.8, n=2.0, a=1.0, sigma_s=0.01)"


# ─────────────────────────────────────────────────────────────────────────────
# HashinShtrikmanBounds
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "kwargs", [{"rho_matrix": 0.0}, {"rho_matrix": -1.0}, {"rho_fluid": 0.0}]
)
def test_hs_bounds_invalid_params(kwargs):
    defaults = {"rho_matrix": 1000.0, "rho_fluid": 25.0}
    defaults.update(kwargs)
    with pytest.raises(ValueError):
        HashinShtrikmanBounds(**defaults)


def test_hs_bounds_ordering():
    hs = HashinShtrikmanBounds(rho_matrix=1000.0, rho_fluid=25.0)
    lower, upper = hs.bounds(phi=0.25)
    assert lower <= upper


def test_hs_bounds_array_input():
    hs = HashinShtrikmanBounds(rho_matrix=1000.0, rho_fluid=25.0)
    lower, upper = hs.bounds(phi=np.array([0.1, 0.3, 0.5]))
    assert lower.shape == (3,)
    assert np.all(lower <= upper)


def test_hs_in_bounds_true_and_false():
    hs = HashinShtrikmanBounds(rho_matrix=1000.0, rho_fluid=25.0)
    lower, upper = hs.bounds(phi=0.25)
    mid = np.sqrt(lower * upper)
    assert bool(hs.in_bounds(mid, phi=0.25)) is True
    assert bool(hs.in_bounds(upper * 1000.0, phi=0.25)) is False


def test_hs_in_bounds_margin_widens_acceptance():
    hs = HashinShtrikmanBounds(rho_matrix=1000.0, rho_fluid=25.0)
    _, upper = hs.bounds(phi=0.25)
    just_outside = upper * 1.5
    assert bool(hs.in_bounds(just_outside, phi=0.25, margin=0.0)) is False
    assert bool(hs.in_bounds(just_outside, phi=0.25, margin=1.0)) is True


def test_hs_repr():
    hs = HashinShtrikmanBounds(rho_matrix=1000.0, rho_fluid=25.0)
    assert "HashinShtrikmanBounds" in repr(hs)


# ─────────────────────────────────────────────────────────────────────────────
# Hydraulic property functions
# ─────────────────────────────────────────────────────────────────────────────


def test_kozeny_carman_K_scalar_and_array():
    K = pp.kozeny_carman_K(0.30, d50_m=1e-3)
    assert isinstance(K, float)
    assert K > 0
    K_arr = pp.kozeny_carman_K(np.array([0.1, 0.3, 0.5]))
    assert K_arr.shape == (3,)
    assert np.all(np.diff(K_arr) > 0)  # K increases with porosity


def test_rho_to_hydraulic_conductivity():
    a = ArchieModel()
    K = pp.rho_to_hydraulic_conductivity(rho=100.0, archie=a, rho_w=20.0)
    assert K > 0


def test_transmissivity_scalar_and_broadcast():
    T = pp.transmissivity(K=1e-4, thickness=50.0)
    assert T == pytest.approx(5e-3)
    T_arr = pp.transmissivity(K=np.array([1e-4, 2e-4]), thickness=10.0)
    assert T_arr.shape == (2,)


def test_transmissivity_clips_negative_inputs():
    T = pp.transmissivity(K=-1.0, thickness=10.0)
    assert T == 0.0


def test_storativity_confined_and_unconfined():
    Sc, Su = pp.storativity(phi=0.3, thickness=50.0, specific_storage=1e-4)
    assert Sc == pytest.approx(50.0 * 1e-4)
    assert Su == pytest.approx(0.3)


# ─────────────────────────────────────────────────────────────────────────────
# Water chemistry
# ─────────────────────────────────────────────────────────────────────────────


def test_rho_w_to_tds_and_back_roundtrip():
    rho_w = 20.0
    tds = pp.rho_w_to_tds(rho_w)
    rho_w2 = pp.tds_to_rho_w(tds)
    assert rho_w2 == pytest.approx(rho_w, rel=1e-6)


def test_rho_w_to_tds_temperature_correction():
    tds_25 = pp.rho_w_to_tds(20.0, temp_c=25.0)
    tds_10 = pp.rho_w_to_tds(20.0, temp_c=10.0)
    assert tds_25 != tds_10


def test_ec_mscm_to_rho_and_back_roundtrip():
    ec = 0.5
    rho = pp.ec_mscm_to_rho(ec)
    ec2 = pp.rho_to_ec_mscm(rho)
    assert ec2 == pytest.approx(ec, rel=1e-6)


def test_ec_rho_array_inputs():
    ec = np.array([0.1, 0.5, 2.0])
    rho = pp.ec_mscm_to_rho(ec)
    assert rho.shape == (3,)


# ─────────────────────────────────────────────────────────────────────────────
# EM-specific helpers
# ─────────────────────────────────────────────────────────────────────────────


def test_skin_depth_and_bostick_depth_identical():
    d1 = pp.skin_depth(rho=100.0, freq=1.0)
    d2 = pp.bostick_depth(rho_a=100.0, freq=1.0)
    assert d1 == pytest.approx(d2)


def test_skin_depth_array():
    d = pp.skin_depth(rho=np.array([10.0, 100.0]), freq=1.0)
    assert d.shape == (2,)
    assert d[1] > d[0]  # higher rho -> deeper penetration


def test_aquifer_top_from_profile_low_direction_found():
    rho_log10 = np.log10(np.array([2000.0, 800.0, 50.0, 20.0]))
    z = np.array([5.0, 15.0, 30.0, 60.0])
    depth = pp.aquifer_top_from_profile(
        rho_log10, z, rho_threshold_ohm_m=300.0, direction="low"
    )
    assert depth == 30.0


def test_aquifer_top_from_profile_high_direction_found():
    rho_log10 = np.log10(np.array([20.0, 50.0, 800.0, 2000.0]))
    z = np.array([5.0, 15.0, 30.0, 60.0])
    depth = pp.aquifer_top_from_profile(
        rho_log10, z, rho_threshold_ohm_m=300.0, direction="high"
    )
    assert depth == 30.0


def test_aquifer_top_from_profile_not_found_returns_none():
    rho_log10 = np.log10(np.array([2000.0, 1800.0, 1500.0]))
    z = np.array([5.0, 15.0, 30.0])
    depth = pp.aquifer_top_from_profile(
        rho_log10, z, rho_threshold_ohm_m=300.0, direction="low"
    )
    assert depth is None


def test_aquifer_top_from_profile_min_depth_skips_shallow_crossing():
    rho_log10 = np.log10(np.array([50.0, 2000.0, 20.0]))
    z = np.array([5.0, 15.0, 30.0])
    depth = pp.aquifer_top_from_profile(
        rho_log10,
        z,
        rho_threshold_ohm_m=300.0,
        direction="low",
        min_depth=10.0,
    )
    assert depth == 30.0  # z=5 crossing skipped by min_depth


def test_water_table_from_profile_found():
    a = ArchieModel(m=1.8, n=2.0)
    rho_log10 = np.log10(np.array([2000.0, 800.0, 15.0, 10.0]))
    z = np.array([5.0, 15.0, 30.0, 60.0])
    depth = pp.water_table_from_profile(rho_log10, z, a, rho_w=20.0)
    assert depth is not None
    assert depth in (30.0, 60.0)


def test_water_table_from_profile_not_found_returns_none():
    a = ArchieModel(m=1.8, n=2.0)
    rho_log10 = np.log10(np.array([5000.0, 4000.0, 3000.0]))
    z = np.array([5.0, 15.0, 30.0])
    depth = pp.water_table_from_profile(rho_log10, z, a, rho_w=20.0)
    assert depth is None


def test_water_table_from_profile_skips_nan_and_nonpositive_and_shallow():
    a = ArchieModel(m=1.8, n=2.0)
    rho_log10 = np.array([np.nan, np.log10(10.0), np.log10(10.0)])
    z = np.array([0.1, 0.3, 30.0])  # first two below default min_depth=0.5
    depth = pp.water_table_from_profile(rho_log10, z, a, rho_w=20.0)
    assert depth == 30.0


def test_fractured_zone_K_scalar_and_array():
    K = pp.fractured_zone_K(rho=1000.0, rho_matrix=5000.0)
    assert 0.0 <= K <= 1e-1
    K_arr = pp.fractured_zone_K(rho=np.array([100.0, 4000.0]), rho_matrix=5000.0)
    assert K_arr.shape == (2,)
    assert K_arr[0] > K_arr[1]  # bigger contrast -> higher fracture K


def test_fractured_zone_K_at_matrix_resistivity_is_near_zero():
    K = pp.fractured_zone_K(rho=5000.0, rho_matrix=5000.0)
    assert K == pytest.approx(0.0, abs=1e-6)
