# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

import numpy as np
import pytest

from pycsamt.utils.features import (
    find_closest_positions,
    find_nearest_indices,
    find_position_bounds,
    select_anomaly_peak,
)

# --------------------------- select_anomaly_peak --------------------------


def test_select_anomaly_peak_raw_with_positions():
    res = [100, 80, 90, 110]
    pos = [0, 1, 2, 3]
    out = select_anomaly_peak(res, pos)
    assert out["peak_pos"] == 1.0
    assert out["rho_peak"] == 80.0
    assert out["pos_bounds"] == (0.0, 3.0)
    assert out["rho_bounds"] == (80.0, 110.0)
    np.testing.assert_allclose(out["profile"], res)


def test_select_anomaly_peak_raw_positions_none_with_dipole():
    res = [100, 80, 90, 110]
    out = select_anomaly_peak(res, positions=None, dipole=5.0)
    assert out["peak_pos"] == 5.0
    assert out["pos_bounds"] == (0.0, 15.0)


def test_select_anomaly_peak_raw_neither_positions_nor_dipole():
    with pytest.raises(ValueError):
        select_anomaly_peak([100, 80, 90], positions=None, dipole=None)


def test_select_anomaly_peak_rank_greater_than_one():
    res = [100, 80, 90, 110]
    pos = [0, 1, 2, 3]
    out = select_anomaly_peak(res, pos, rank=2)
    assert out["rho_peak"] == 90.0
    assert out["peak_pos"] == 2.0


def test_select_anomaly_peak_anomaly_infos_matching_key():
    infos = {"1_pk_a": np.array([50.0, 20.0, 60.0])}
    out = select_anomaly_peak(
        [1, 2, 3], [0, 1, 2], rank=1, anomaly_infos=infos
    )
    assert out["rho_peak"] == 20.0
    assert out["peak_pos"] == 1.0


def test_select_anomaly_peak_anomaly_infos_no_matching_key():
    infos = {"2_pk_a": np.array([50.0, 20.0, 60.0])}
    with pytest.raises(KeyError):
        select_anomaly_peak([1, 2, 3], [0, 1, 2], rank=1, anomaly_infos=infos)


def test_select_anomaly_peak_anomaly_infos_positions_derived_from_dipole():
    infos = {"1_pk_a": np.array([50.0, 20.0, 60.0])}
    out = select_anomaly_peak(
        [1, 2, 3], positions=None, dipole=2.0, rank=1, anomaly_infos=infos
    )
    assert out["peak_pos"] == 2.0


def test_select_anomaly_peak_anomaly_infos_positions_none_dipole_none():
    infos = {"1_pk_a": np.array([50.0, 20.0, 60.0])}
    with pytest.raises(ValueError):
        select_anomaly_peak(
            [1, 2, 3],
            positions=None,
            dipole=None,
            rank=1,
            anomaly_infos=infos,
        )


def test_select_anomaly_peak_user_peak_inside_bounds():
    res = [100, 80, 90, 110]
    pos = [0, 1, 2, 3]
    out = select_anomaly_peak(res, pos, user_peak=2.5)
    assert out["peak_pos"] == 2.5


def test_select_anomaly_peak_user_peak_outside_bounds_warns():
    res = [100, 80, 90, 110]
    pos = [0, 1, 2, 3]
    with pytest.warns(UserWarning):
        out = select_anomaly_peak(res, pos, user_peak=50.0)
    assert out["peak_pos"] == 1.0


def test_select_anomaly_peak_return_bounds_false():
    res = [100, 80, 90, 110]
    pos = [0, 1, 2, 3]
    out = select_anomaly_peak(res, pos, return_bounds=False)
    assert "pos_bounds" not in out
    assert "rho_bounds" not in out


def test_select_anomaly_peak_empty_resistivity():
    with pytest.raises(ValueError):
        select_anomaly_peak([], [])


def test_select_anomaly_peak_positions_length_mismatch():
    with pytest.raises(ValueError):
        select_anomaly_peak([100, 80, 90], [0, 1])


# --------------------------- find_position_bounds -------------------------


def test_find_position_bounds_peak_as_float():
    out = find_position_bounds(1.0, 999, [100, 80, 90])
    assert out["pos_min"] == 1.0
    assert out["pos_max"] == 3.0
    assert out["rho_min"] == 80.0
    assert out["rho_max"] == 100.0


def test_find_position_bounds_peak_as_string():
    out = find_position_bounds("pk12.5", 80, [100, 80, 90])
    assert out["pos_min"] == 11.5
    assert out["pos_max"] == 13.5


def test_find_position_bounds_peak_string_parse_failure():
    with pytest.raises(ValueError):
        find_position_bounds("abc", 80, [100, 80, 90])


def test_find_position_bounds_positions_given():
    out = find_position_bounds(
        1.0, 80, [100, 80, 90], positions=[10, 20, 30]
    )
    assert out["pos_min"] == 10.0
    assert out["pos_max"] == 30.0


def test_find_position_bounds_positions_auto_from_dipole():
    out = find_position_bounds(5.0, 80, [100, 80, 90], dipole=2.0)
    assert out["pos_min"] == 3.0
    assert out["pos_max"] == 7.0


def test_find_position_bounds_empty_rho_range():
    with pytest.raises(ValueError):
        find_position_bounds(1.0, 80, [])


def test_find_position_bounds_positions_length_mismatch():
    with pytest.raises(ValueError):
        find_position_bounds(1.0, 80, [100, 80, 90], positions=[0, 1])


def test_find_position_bounds_rho_peak_not_present_falls_back_to_zero():
    out = find_position_bounds(5.0, 12345.0, [100, 80, 90], dipole=1.0)
    assert out["pos_min"] == 5.0
    assert out["pos_max"] == 7.0


# --------------------------- find_closest_positions ------------------------


def test_find_closest_positions_return_indices_only():
    ref = [0, 10, 20, 30]
    idxs = find_closest_positions(ref, [12, -5, 25])
    assert idxs == [1, 0, 2]


def test_find_closest_positions_return_values():
    ref = [0, 10, 20, 30]
    idxs, vals = find_closest_positions(ref, [12, -5], return_values=True)
    assert idxs == [1, 0]
    assert vals == [10.0, 0.0]


def test_find_closest_positions_non_1d_raises():
    with pytest.raises(ValueError):
        find_closest_positions([[1, 2], [3, 4]], [1, 2])


# --------------------------- find_nearest_indices ---------------------------


def test_find_nearest_indices_side_both():
    ref = [0, 10, 20, 30]
    assert find_nearest_indices(ref, [12, -5, 25]) == [1, 0, 2]


def test_find_nearest_indices_side_left():
    ref = [0, 10, 20, 30]
    assert find_nearest_indices(ref, [12, -5], side="left") == [1, 0]


def test_find_nearest_indices_side_right_return_values():
    ref = [0, 10, 20, 30]
    idxs, vals = find_nearest_indices(
        ref, [12, 35], side="right", return_values=True
    )
    assert idxs == [2, 3]
    assert vals == [20.0, 30.0]


def test_find_nearest_indices_invalid_side():
    with pytest.raises(ValueError):
        find_nearest_indices([0, 10, 20], [5], side="bogus")


def test_find_nearest_indices_left_no_candidate_falls_back_to_first():
    ref = [10, 20, 30]
    assert find_nearest_indices(ref, [5], side="left") == [0]


def test_find_nearest_indices_right_no_candidate_falls_back_to_last():
    ref = [10, 20, 30]
    assert find_nearest_indices(ref, [100], side="right") == [2]
