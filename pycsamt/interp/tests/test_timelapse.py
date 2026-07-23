# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp.timelapse — TimeLapseEM, assert_compatible_grids."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.interp._base import ResistivityModel
from pycsamt.interp.petrophysics import ArchieModel, WaxmanSmitsModel
from pycsamt.interp.timelapse import (
    TimeLapseEM,
    _to_archie,
    assert_compatible_grids,
)

_Z = np.array([10.0, 30.0, 60.0])
_X = np.array([0.0, 100.0])


def _survey(rho_linear, method="TDEM"):
    rho = np.log10(np.array(rho_linear))
    return ResistivityModel.from_array(
        rho, x_centers=_X, z_centers=_Z, method=method
    )


def _dry():
    return _survey([[1000.0, 1000.0], [800.0, 800.0], [600.0, 600.0]])


def _wet():
    return _survey([[500.0, 500.0], [100.0, 100.0], [50.0, 50.0]])


def _recharge():
    return _survey([[200.0, 200.0], [50.0, 50.0], [20.0, 20.0]])


# ─────────────────────────────────────────────────────────────────────────────
# assert_compatible_grids
# ─────────────────────────────────────────────────────────────────────────────


def test_assert_compatible_grids_ok():
    assert_compatible_grids([_dry(), _wet(), _recharge()])  # no raise


def test_assert_compatible_grids_shape_mismatch():
    other = ResistivityModel.from_array(
        np.zeros((2, 2)),
        x_centers=_X,
        z_centers=np.array([10.0, 30.0]),
        method="X",
    )
    with pytest.raises(ValueError, match="shape"):
        assert_compatible_grids([_dry(), other])


def test_assert_compatible_grids_x_mismatch():
    other = _survey([[1000.0, 1000.0], [800.0, 800.0], [600.0, 600.0]])
    other.x_centers = np.array([0.0, 999.0])
    with pytest.raises(ValueError, match="x_centers"):
        assert_compatible_grids([_dry(), other])


def test_assert_compatible_grids_z_mismatch():
    other = _survey([[1000.0, 1000.0], [800.0, 800.0], [600.0, 600.0]])
    other.z_centers = np.array([10.0, 30.0, 9999.0])
    with pytest.raises(ValueError, match="z_centers"):
        assert_compatible_grids([_dry(), other])


# ─────────────────────────────────────────────────────────────────────────────
# TimeLapseEM — construction
# ─────────────────────────────────────────────────────────────────────────────


def test_init_requires_at_least_two_surveys():
    with pytest.raises(ValueError, match="At least two"):
        TimeLapseEM([_dry()])


def test_init_default_times_and_labels():
    tl = TimeLapseEM([_dry(), _wet()])
    assert tl.times == [0, 1]
    assert tl.labels == ["T00", "T01"]


def test_init_custom_times_and_labels():
    tl = TimeLapseEM(
        [_dry(), _wet(), _recharge()],
        times=[0.0, 30.0, 90.0],
        labels=["dry", "wet", "recharge"],
    )
    assert tl.times == [0.0, 30.0, 90.0]
    assert tl.labels == ["dry", "wet", "recharge"]


def test_init_times_length_mismatch_raises():
    with pytest.raises(ValueError, match="len\\(times\\)"):
        TimeLapseEM([_dry(), _wet()], times=[0.0])


def test_init_labels_length_mismatch_raises():
    with pytest.raises(ValueError, match="len\\(labels\\)"):
        TimeLapseEM([_dry(), _wet()], labels=["only-one"])


def test_init_propagates_grid_mismatch():
    other = ResistivityModel.from_array(
        np.zeros((2, 2)), x_centers=_X, z_centers=np.array([1.0, 2.0])
    )
    with pytest.raises(ValueError, match="shape"):
        TimeLapseEM([_dry(), other])


def test_properties():
    tl = TimeLapseEM([_dry(), _wet(), _recharge()])
    assert tl.n_surveys == 3
    assert tl.n_x == 2
    assert tl.n_z == 3


def test_repr():
    tl = TimeLapseEM([_dry(), _wet()], labels=["a", "b"])
    r = repr(tl)
    assert "n_surveys=2" in r
    assert "'a'" in r


# ─────────────────────────────────────────────────────────────────────────────
# resistivity_change
# ─────────────────────────────────────────────────────────────────────────────


def test_resistivity_change_default_baseline():
    tl = TimeLapseEM([_dry(), _wet(), _recharge()])
    deltas = tl.resistivity_change()
    assert len(deltas) == 2
    # wetting -> resistivity decreases relative to dry baseline
    assert np.all(deltas[0] < 0)


def test_resistivity_change_custom_baseline():
    tl = TimeLapseEM([_dry(), _wet(), _recharge()])
    deltas = tl.resistivity_change(baseline_idx=1)
    assert len(deltas) == 2  # surveys 0 and 2, relative to survey 1


# ─────────────────────────────────────────────────────────────────────────────
# saturation_change / water_content_change
# ─────────────────────────────────────────────────────────────────────────────


def test_saturation_change_archie_scalar_phi():
    tl = TimeLapseEM([_dry(), _wet()])
    dsw = tl.saturation_change(ArchieModel(), phi=0.25, rho_w=0.025)
    assert len(dsw) == 1
    assert dsw[0].shape == (3, 2)
    # wetter survey -> saturation increase relative to dry baseline
    assert np.all(dsw[0] >= 0)


def test_saturation_change_waxman_smits_petro():
    tl = TimeLapseEM([_dry(), _wet()])
    dsw = tl.saturation_change(WaxmanSmitsModel(sigma_s=0.01), phi=0.25)
    assert dsw[0].shape == (3, 2)


def test_saturation_change_array_phi():
    tl = TimeLapseEM([_dry(), _wet()])
    phi_arr = np.full((3, 2), 0.3)
    dsw = tl.saturation_change(ArchieModel(), phi=phi_arr)
    assert dsw[0].shape == (3, 2)


def test_saturation_change_invalid_petro_type_raises():
    tl = TimeLapseEM([_dry(), _wet()])
    with pytest.raises(TypeError):
        tl.saturation_change("not-a-petro-model")


def test_water_content_change_equals_dsw_times_phi():
    tl = TimeLapseEM([_dry(), _wet()])
    dsw = tl.saturation_change(ArchieModel(), phi=0.25)
    dtheta = tl.water_content_change(ArchieModel(), phi=0.25)
    np.testing.assert_allclose(dtheta[0], dsw[0] * 0.25)


# ─────────────────────────────────────────────────────────────────────────────
# water_table_displacement / water_table_map
# ─────────────────────────────────────────────────────────────────────────────


def test_water_table_displacement_two_surveys_returns_single_row():
    tl = TimeLapseEM([_dry(), _wet()])
    disp = tl.water_table_displacement(ArchieModel(), rho_w=20.0)
    assert disp.shape == (2,)  # single row (n_x,), not stacked


def test_water_table_displacement_multi_survey_stacks_rows():
    tl = TimeLapseEM([_dry(), _wet(), _recharge()])
    disp = tl.water_table_displacement(ArchieModel(), rho_w=20.0)
    assert disp.shape == (2, 2)  # (n_surveys-1, n_x)


def test_water_table_map_shape_and_nan_where_undetected():
    resistive_only = _survey(
        [[5000.0, 5000.0], [4000.0, 4000.0], [3000.0, 3000.0]]
    )
    tl = TimeLapseEM([resistive_only, _wet()])
    wt = tl.water_table_map(ArchieModel(), rho_w=20.0)
    assert wt.shape == (2, 2)
    assert np.all(np.isnan(wt[0]))  # never crosses threshold


# ─────────────────────────────────────────────────────────────────────────────
# resistivity_stats
# ─────────────────────────────────────────────────────────────────────────────


def test_resistivity_stats_keys_and_shapes():
    tl = TimeLapseEM([_dry(), _wet(), _recharge()])
    stats = tl.resistivity_stats()
    assert set(stats) == {
        "mean_delta",
        "std_delta",
        "max_increase",
        "max_decrease",
    }
    for v in stats.values():
        assert v.shape == (3, 2)


# ─────────────────────────────────────────────────────────────────────────────
# _to_archie
# ─────────────────────────────────────────────────────────────────────────────


def test_to_archie_passthrough():
    a = ArchieModel(m=1.9)
    assert _to_archie(a) is a


def test_to_archie_converts_waxman_smits():
    ws = WaxmanSmitsModel(m=1.9, n=2.1, a=1.05, sigma_s=0.02)
    a = _to_archie(ws)
    assert isinstance(a, ArchieModel)
    assert a.m == pytest.approx(1.9)
    assert a.n == pytest.approx(2.1)
    assert a.a == pytest.approx(1.05)


def test_to_archie_invalid_type_raises():
    with pytest.raises(TypeError, match="petro must be"):
        _to_archie(object())
