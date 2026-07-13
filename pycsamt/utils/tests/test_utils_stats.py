# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.utils.stats`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import StatsError
from pycsamt.utils.stats import (
    drawn_boundaries,
    get_confidence_ratio,
    remove_outliers,
    scale_position,
)

# ------------------------- get_confidence_ratio -----------------------


def test_confidence_ratio_counts_nans_per_column():
    arr = np.array(
        [
            [1.0, np.nan],
            [2.0, 5.0],
            [3.0, np.nan],
            [4.0, 6.0],
        ]
    )
    ratio = get_confidence_ratio(arr, axis=0)
    assert np.allclose(ratio, [1.0, 0.5])


def test_confidence_ratio_axis_names_and_percent():
    arr = np.array(
        [
            [1.0, np.nan],
            [2.0, 5.0],
        ]
    )
    # 'cols'/'x' -> axis 1: ratio per row
    per_row = get_confidence_ratio(arr, axis="cols", as_percent=True)
    assert np.allclose(per_row, [50.0, 100.0])
    per_col = get_confidence_ratio(arr, axis="rows")
    assert np.allclose(per_col, [1.0, 0.5])


def test_confidence_ratio_invalid_values_and_average():
    arr = np.array(
        [
            [0.0, 1.0],
            [0.0, 2.0],
        ]
    )
    ratio = get_confidence_ratio(arr, axis=0, invalid=[0.0])
    assert np.allclose(ratio, [0.0, 1.0])
    avg = get_confidence_ratio(arr, axis=0, invalid=[0.0], average=True)
    assert avg == pytest.approx(0.5)


def test_confidence_ratio_scalar_and_errors():
    assert get_confidence_ratio(np.float64(3.0)) == 1.0
    assert get_confidence_ratio(np.float64("nan")) == 0.0
    with pytest.raises(StatsError):
        get_confidence_ratio(np.ones((2, 2)), axis="diagonal")
    with pytest.raises(StatsError):
        get_confidence_ratio(np.ones((2, 2)), axis=5)


# ---------------------------- remove_outliers -------------------------


def _data_with_outlier():
    base = np.array([10.0, 11.0, 9.0, 10.5, 9.5, 10.2, 9.8, 1000.0])
    return base


def test_remove_outliers_iqr_1d_drops_value():
    out = remove_outliers(_data_with_outlier(), method="IQR")
    assert 1000.0 not in out
    assert out.size == 7


def test_remove_outliers_fill_value_nan_and_interpolate():
    arr = _data_with_outlier()
    filled = remove_outliers(arr, fill_value=np.nan)
    assert np.isnan(filled[-1])
    assert filled.size == arr.size

    interp = remove_outliers(
        arr.reshape(1, -1), fill_value=np.nan, interpolate=True
    )
    assert not np.isnan(interp).any()


def test_remove_outliers_zscore_threshold():
    arr = _data_with_outlier()
    out = remove_outliers(arr, method="z-score", threshold=2.0)
    assert 1000.0 not in out


def test_remove_outliers_2d_drops_rows():
    arr = np.vstack(
        [
            np.full(3, 10.0),
            np.full(3, 11.0),
            np.array([10.0, 10.0, 1e6]),
        ]
    )
    out = remove_outliers(arr, method="IQR")
    assert out.shape[0] == 2


def test_remove_outliers_dataframe_paths():
    df = pd.DataFrame(
        {
            "a": [10.0, 11.0, 9.0, 10.5, 9.5, 10.2, 9.8, 1000.0],
            "b": np.arange(8.0),
        }
    )
    dropped = remove_outliers(df, method="IQR")
    assert isinstance(dropped, pd.DataFrame)
    assert 1000.0 not in dropped["a"].values

    filled = remove_outliers(df, fill_value=0.0)
    assert filled.shape == df.shape
    assert filled.loc[7, "a"] == 0.0


def test_remove_outliers_bad_inputs_raise():
    with pytest.raises(StatsError):
        remove_outliers(_data_with_outlier(), method="mad")
    with pytest.raises(StatsError):
        remove_outliers(pd.DataFrame({"a": [1.0]}), method="mad")
    with pytest.raises(StatsError):
        remove_outliers(np.array(["x", "y"], dtype=object))


# ----------------------------- scale_position -------------------------


def test_scale_position_recovers_linear_model():
    x = np.arange(10.0)
    y = 2.0 * x + 1.0
    y_fit, popt, pcov = scale_position(y, x)
    assert np.allclose(y_fit, y, atol=1e-8)
    assert popt == pytest.approx([2.0, 1.0])
    assert pcov.shape == (2, 2)


def test_scale_position_defaults_to_index_and_series_output():
    y = pd.Series([1.0, 3.0, 5.0, 7.0])
    y_fit, popt, _ = scale_position(y, asarray=False)
    assert isinstance(y_fit, pd.Series)
    assert popt[0] == pytest.approx(2.0)


def test_scale_position_dataframe_requires_column():
    df = pd.DataFrame({"a": [1.0, 2.0, 3.0]})
    with pytest.raises(ValueError):
        scale_position(df)
    y_fit, popt, _ = scale_position(df, column="a")
    assert popt[0] == pytest.approx(1.0)


def test_scale_position_custom_func():
    x = np.linspace(0, 2, 20)
    y = 3.0 * x**2

    def quad(x, a):
        return a * x**2

    _, popt, _ = scale_position(y, x, func=quad, initial_params=[1.0])
    assert popt[0] == pytest.approx(3.0)


# ---------------------------- drawn_boundaries ------------------------


def test_drawn_boundaries_assembles_walls():
    profile = [5.0, 3.0, 1.0, 4.0, 6.0]
    peak, idx, bounds = drawn_boundaries(profile, 1.0, 2)
    assert peak == 1.0
    assert idx == 2
    assert np.allclose(bounds, [5.0, 3.0, 1.0, 4.0, 6.0])


def test_drawn_boundaries_peak_at_edge():
    profile = [1.0, 2.0, 3.0]
    peak, idx, bounds = drawn_boundaries(profile, 1.0, 0)
    assert np.allclose(bounds, [1.0, 2.0, 3.0])


def test_drawn_boundaries_out_of_bounds_raises():
    with pytest.raises(StatsError):
        drawn_boundaries([1.0, 2.0], 1.0, 5)
