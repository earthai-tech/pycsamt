# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp.calibrate — ModelCalibrator."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.interp._base import ResistivityModel
from pycsamt.interp.borehole import Borehole, Interval
from pycsamt.interp.calibrate import ModelCalibrator
from pycsamt.interp.lithology import RockDatabase


def _model():
    rho = np.log10(
        np.array(
            [
                [100.0, 12.0, 1500.0],
                [80.0, 25.0, 1200.0],
                [250.0, 70.0, 900.0],
                [1500.0, 1800.0, 2500.0],
            ]
        )
    )
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 500.0, 1000.0]),
        z_centers=np.array([25.0, 75.0, 150.0, 300.0]),
        station_x=np.array([0.0, 500.0, 1000.0]),
        station_names=["S1", "S2", "S3"],
        method="test",
    )


def _model_with_nan():
    rho = np.log10(
        np.array(
            [
                [100.0, 12.0, 1500.0],
                [80.0, 25.0, 1200.0],
                [250.0, 70.0, 900.0],
                [1500.0, 1800.0, 2500.0],
            ]
        )
    )
    rho[1, 0] = np.nan
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 500.0, 1000.0]),
        z_centers=np.array([25.0, 75.0, 150.0, 300.0]),
        station_x=np.array([0.0, 500.0, 1000.0]),
        station_names=["S1", "S2", "S3"],
        method="test",
    )


def _bh_partial(x=0.0):
    """Covers z=25 (soft-replace), leaves z=75 TRES-less (interval w/ None),
    leaves z=150 in a gap (no interval), covers z=300 far from CRM (step-2)."""
    return Borehole(
        "BoA",
        x=x,
        intervals=[
            Interval(top=0.0, bottom=50.0, lithology="topsoil", resistivity=105.0),
            Interval(top=50.0, bottom=100.0, lithology="clay", resistivity=None),
            Interval(
                top=200.0, bottom=400.0, lithology="basement", resistivity=3000.0
            ),
        ],
    )


# ─────────────────────────────────────────────────────────────────────────────
# Construction
# ─────────────────────────────────────────────────────────────────────────────


def test_default_init():
    cal = ModelCalibrator()
    assert cal.ptol == 0.10
    assert cal.max_borehole_distance == 500.0
    assert isinstance(cal.db, RockDatabase)
    assert cal.verbose is True
    assert repr(cal) == "ModelCalibrator(ptol=0.1, unfitted)"


def test_custom_db_passed_through():
    db = RockDatabase.default()
    cal = ModelCalibrator(db=db)
    assert cal.db is db


# ─────────────────────────────────────────────────────────────────────────────
# Guard rails before fit()
# ─────────────────────────────────────────────────────────────────────────────


def test_access_before_fit_raises():
    cal = ModelCalibrator()
    with pytest.raises(RuntimeError, match="Call fit"):
        cal.calibrated_model()
    with pytest.raises(RuntimeError, match="Call fit"):
        cal.misfit_map()
    with pytest.raises(RuntimeError, match="Call fit"):
        cal.stratigraphic_logs()


# ─────────────────────────────────────────────────────────────────────────────
# fit — no boreholes / far boreholes -> autolayer-only path
# ─────────────────────────────────────────────────────────────────────────────


def test_fit_without_boreholes_uses_autolayer_everywhere(capsys):
    cal = ModelCalibrator(verbose=True).fit(_model(), boreholes=None)
    assert cal._is_fitted is True
    nm = cal.calibrated_model()
    assert nm.method == "test+calibrated"
    assert nm.rho_2d.shape == (4, 3)

    db = RockDatabase.default()
    expected_00 = db.classify(100.0).log_rho_mid
    assert nm.rho_2d[0, 0] == pytest.approx(expected_00)

    out = capsys.readouterr().out
    assert "ModelCalibrator: fitted 3 columns" in out
    assert "0 borehole(s)" in out


def test_fit_with_empty_borehole_list_same_as_none():
    cal1 = ModelCalibrator(verbose=False).fit(_model(), boreholes=[])
    cal2 = ModelCalibrator(verbose=False).fit(_model(), boreholes=None)
    np.testing.assert_array_equal(
        cal1.calibrated_model().rho_2d, cal2.calibrated_model().rho_2d
    )


def test_fit_borehole_outside_max_distance_falls_back_to_autolayer():
    far_bh = _bh_partial(x=100_000.0)
    cal = ModelCalibrator(verbose=False).fit(_model(), boreholes=[far_bh])
    nm = cal.calibrated_model()
    db = RockDatabase.default()
    expected_00 = db.classify(100.0).log_rho_mid
    assert nm.rho_2d[0, 0] == pytest.approx(expected_00)


def test_fit_without_boreholes_skips_nan_crm_cells():
    """_autolayer_column must leave nan cells untouched (no borehole path)."""
    cal = ModelCalibrator(verbose=False).fit(_model_with_nan(), boreholes=None)
    nm = cal.calibrated_model()
    assert np.isnan(nm.rho_2d[1, 0])


def test_fit_verbose_false_suppresses_output(capsys):
    ModelCalibrator(verbose=False).fit(_model())
    out = capsys.readouterr().out
    assert out == ""


# ─────────────────────────────────────────────────────────────────────────────
# fit — borehole within range: soft-replace / autolayer branches
# ─────────────────────────────────────────────────────────────────────────────


def test_fit_borehole_soft_replace_branch():
    bh = _bh_partial(x=0.0)
    cal = ModelCalibrator(ptol=0.10, verbose=False).fit(_model(), [bh])
    nm = cal.calibrated_model()
    # z=25: TRES=105 vs CRM=100 -> xi_se=0.048 <= 0.10 -> soft replace
    assert nm.rho_2d[0, 0] == pytest.approx(np.log10(105.0))


def test_fit_borehole_tres_missing_uses_autolayer():
    bh = _bh_partial(x=0.0)
    cal = ModelCalibrator(verbose=False).fit(_model(), [bh])
    nm = cal.calibrated_model()
    db = RockDatabase.default()
    # z=75: interval present but resistivity=None -> TRES nan -> autolayer(CRM=80)
    assert nm.rho_2d[1, 0] == pytest.approx(db.classify(80.0).log_rho_mid)
    # z=150: falls in the interval gap -> TRES nan -> autolayer(CRM=250)
    assert nm.rho_2d[2, 0] == pytest.approx(db.classify(250.0).log_rho_mid)


def test_fit_borehole_step2_autolayer_when_outside_ptol():
    bh = _bh_partial(x=0.0)
    cal = ModelCalibrator(ptol=0.10, verbose=False).fit(_model(), [bh])
    nm = cal.calibrated_model()
    db = RockDatabase.default()
    # z=300: TRES=3000 vs CRM=1500 -> xi_se=0.5 > 0.10 -> autolayer(TRES=3000)
    assert nm.rho_2d[3, 0] == pytest.approx(db.classify(3000.0).log_rho_mid)


def test_fit_nan_crm_cell_is_skipped():
    bh = _bh_partial(x=0.0)
    cal = ModelCalibrator(verbose=False).fit(_model_with_nan(), [bh])
    nm = cal.calibrated_model()
    assert np.isnan(nm.rho_2d[1, 0])


def test_fit_selects_nearest_of_multiple_boreholes():
    near = _bh_partial(x=0.0)
    far = Borehole(
        "BoB",
        x=1000.0,
        intervals=[Interval(top=0.0, bottom=400.0, lithology="x", resistivity=999.0)],
    )
    cal = ModelCalibrator(verbose=False).fit(_model(), [far, near])
    nm = cal.calibrated_model()
    # station S1 (x=0) is nearest to `near` -> soft-replace with 105.0 applies
    assert nm.rho_2d[0, 0] == pytest.approx(np.log10(105.0))


# ─────────────────────────────────────────────────────────────────────────────
# misfit_map
# ─────────────────────────────────────────────────────────────────────────────


def test_misfit_map_shape_and_range():
    bh = _bh_partial(x=0.0)
    cal = ModelCalibrator(verbose=False).fit(_model(), [bh])
    g = cal.misfit_map()
    assert g.shape == (4, 3)
    assert np.all(g >= 0)
    # returned array must be a copy, not the internal buffer
    g[0, 0] = -1.0
    assert cal.misfit_map()[0, 0] != -1.0


def test_compute_misfit_zero_crm_column_no_crash():
    crm = np.zeros((3, 2))
    nm = np.full((3, 2), 0.01)
    g = ModelCalibrator._compute_misfit(crm, nm)
    assert np.all(np.isfinite(g))
    assert g.shape == (3, 2)


# ─────────────────────────────────────────────────────────────────────────────
# stratigraphic_logs
# ─────────────────────────────────────────────────────────────────────────────


def test_stratigraphic_logs_one_per_station():
    cal = ModelCalibrator(verbose=False).fit(_model())
    logs = cal.stratigraphic_logs()
    assert len(logs) == 3
    assert {log.station_name for log in logs} == {"S1", "S2", "S3"}


def test_stratigraphic_logs_crm_vs_nm_differ_after_calibration():
    bh = _bh_partial(x=0.0)
    cal = ModelCalibrator(verbose=False).fit(_model(), [bh])
    logs_nm = cal.stratigraphic_logs(model="nm")
    logs_crm = cal.stratigraphic_logs(model="crm")
    log_nm_s1 = next(l for l in logs_nm if l.station_name == "S1")
    log_crm_s1 = next(l for l in logs_crm if l.station_name == "S1")
    assert not np.allclose(log_nm_s1.rho_log10, log_crm_s1.rho_log10)


def test_stratigraphic_logs_custom_db_and_merge_tolerance():
    cal = ModelCalibrator(verbose=False).fit(_model())
    custom_db = RockDatabase.default()
    logs = cal.stratigraphic_logs(db=custom_db, merge_tolerance=0.5)
    assert len(logs) == 3


def test_stratigraphic_logs_fallback_to_all_columns_when_no_stations():
    rho = np.log10(np.array([[10.0, 100.0], [50.0, 500.0]]))
    model = ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 250.0]),
        z_centers=np.array([10.0, 30.0]),
        station_x=np.array([]),
        station_names=[],
        method="nostations",
    )
    cal = ModelCalibrator(verbose=False).fit(model)
    logs = cal.stratigraphic_logs()
    assert len(logs) == 2
    assert logs[0].station_name == "S000"
    assert logs[1].station_name == "S001"


# ─────────────────────────────────────────────────────────────────────────────
# repr
# ─────────────────────────────────────────────────────────────────────────────


def test_repr_fitted():
    cal = ModelCalibrator(verbose=False).fit(_model())
    assert repr(cal) == "ModelCalibrator(ptol=0.1, fitted)"
