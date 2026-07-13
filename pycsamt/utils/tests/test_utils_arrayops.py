# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.utils.arrayops`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.utils.arrayops import (
    assert_xy_in,
    concat_array_from_list,
    drop_nan_in,
    fill_nan,
    frameify,
    interpolate_grid,
    is_iterable,
    reshape,
)

# ------------------------ concat_array_from_list ----------------------


def test_concat_columns_with_padding_and_none():
    a = np.arange(3)
    result = concat_array_from_list([a, None, [5]])
    assert result.shape == (3, 3)
    assert np.allclose(result[:, 0], [0, 1, 2])
    assert np.isnan(result[:, 1]).all()
    assert result[0, 2] == 5.0
    assert np.isnan(result[1:, 2]).all()


def test_concat_rows_axis0():
    result = concat_array_from_list([[1, 2, 3], [4, 5]], concat_axis=0)
    assert result.shape == (2, 3)
    assert np.isnan(result[1, 2])


def test_concat_empty_and_invalid_axis():
    assert concat_array_from_list([]).shape == (0, 0)
    with pytest.raises(ValueError):
        concat_array_from_list([[1]], concat_axis=2)


def test_concat_flattens_higher_dims_and_scalars():
    result = concat_array_from_list([np.ones((2, 2)), 7])
    assert result.shape == (4, 2)
    assert result[0, 1] == 7.0


# ------------------------------ is_iterable ---------------------------


def test_is_iterable_bool_semantics():
    assert is_iterable([1, 2])
    assert is_iterable("abc")
    assert not is_iterable("abc", exclude_string=True)
    assert not is_iterable(42)


def test_is_iterable_transform_and_parse():
    assert is_iterable(5, transform=True) == [5]
    assert is_iterable("ab", exclude_string=True, transform=True) == ["ab"]
    parsed = is_iterable("a b", transform=True, parse_string=True)
    assert parsed == ["a", "b"]
    with pytest.raises(ValueError):
        is_iterable("a b", parse_string=True)


# -------------------------------- reshape -----------------------------


def test_reshape_1d_orientations():
    a = np.arange(4)
    assert reshape(a) is not None
    assert reshape(a, axis=0).shape == (4, 1)
    assert reshape(a, axis=1).shape == (1, 4)
    assert reshape(a, axis=None).shape == (4,)


def test_reshape_2d_squeeze_and_passthrough():
    col = np.arange(4).reshape(4, 1)
    row = np.arange(4).reshape(1, 4)
    full = np.arange(6).reshape(2, 3)
    assert reshape(col).shape == (4,)
    assert reshape(row).shape == (4,)
    assert reshape(full).shape == (2, 3)
    assert reshape(full, axis=0).shape == (2, 3)


def test_reshape_invalid_inputs():
    with pytest.raises(ValueError):
        reshape(np.zeros((2, 2, 2)))
    with pytest.raises(ValueError):
        reshape(np.arange(3), axis=2)


# ------------------------------- frameify -----------------------------


def test_frameify_array_with_columns_and_types():
    arr = np.array([[1.0, "x"], [2.0, "y"]], dtype=object)
    df, nf, cf = frameify(
        arr, columns=["num", "cat"], return_feature_types=True
    )
    assert list(df.columns) == ["num", "cat"]
    assert nf == ["num"]
    assert cf == ["cat"]
    assert df["num"].dtype == np.float64


def test_frameify_drops_all_nan_columns_and_rows():
    df_in = pd.DataFrame(
        {
            "a": [1.0, 2.0, np.nan],
            "b": [np.nan, np.nan, np.nan],
        }
    )
    df = frameify(df_in)
    assert "b" not in df.columns
    assert len(df) == 2  # last row was all-NaN after dropping b


def test_frameify_replaces_empty_strings():
    df = frameify(
        pd.DataFrame(
            {
                "a": ["1", " ", "3"],
                "b": ["4", "5", "6"],
            }
        )
    )
    assert df["a"].isna().sum() == 1
    assert len(df) == 3  # row survives thanks to column b


def test_frameify_pop_cat_features():
    df_in = pd.DataFrame({"n": [1, 2], "c": ["u", "v"]})
    df = frameify(df_in, pop_cat_features=True)
    assert list(df.columns) == ["n"]


def test_frameify_rejects_non_array():
    with pytest.raises(TypeError):
        frameify("not-an-array")


# ------------------------------ assert_xy_in --------------------------


def test_assert_xy_in_from_dataframe_columns():
    data = pd.DataFrame({"e": [1.0, 2.0], "n": [3.0, 4.0]})
    x, y = assert_xy_in("e", "n", data=data)
    assert isinstance(x, np.ndarray)
    assert np.allclose(x, [1.0, 2.0])
    assert np.allclose(y, [3.0, 4.0])


def test_assert_xy_in_series_output_and_frame():
    x, y = assert_xy_in([1, 2], [3, 4], asarray=False)
    assert isinstance(x, pd.Series)
    df = assert_xy_in(
        [1, 2], [3, 4], to_frame=True, columns=["east", "north"]
    )
    assert list(df.columns) == ["east", "north"]


def test_assert_xy_in_dropna_and_numeric():
    x, y = assert_xy_in([1.0, np.nan, 3.0], [4.0, 5.0, 6.0], dropna=True)
    assert x.size == y.size == 2
    x2, _ = assert_xy_in(["1", "2"], ["3", "4"], xy_numeric=True)
    assert x2.dtype.kind in "if"


def test_assert_xy_in_errors():
    with pytest.raises(TypeError):
        assert_xy_in("col", [1, 2])
    with pytest.raises(KeyError):
        assert_xy_in("missing", "n", data=pd.DataFrame({"n": [1]}))
    with pytest.raises(ValueError):
        assert_xy_in([1, 2, 3], [1, 2])
    with pytest.raises(ValueError):
        assert_xy_in([1], [2], to_frame=True, columns=["only"])


# ---------------------------- interpolate_grid ------------------------


def test_interpolate_grid_fills_all_nans():
    x = [28, np.nan, 50, 60]
    y = [np.nan, 1000, 2000, 3000]
    xy = np.vstack((x, y)).T
    out = interpolate_grid(xy)
    assert out.shape == xy.shape
    assert not np.isnan(out).any()
    # known values preserved
    assert out[0, 0] == pytest.approx(28.0)
    assert out[3, 1] == pytest.approx(3000.0)


def test_interpolate_grid_1d_input():
    arr = np.array([1.0, np.nan, 3.0, 4.0])
    out = interpolate_grid(arr, method="linear")
    assert out.ndim == 1
    assert not np.isnan(out).any()


def test_interpolate_grid_constant_fill_value():
    arr = np.array([[1.0, np.nan], [np.nan, 4.0]])
    out = interpolate_grid(arr, fill_value=0.0)
    assert not np.isnan(out).any()


# -------------------------------- fill_nan ----------------------------


def test_fill_nan_1d_directions():
    a = np.array([np.nan, 1.0, np.nan, 2.0])
    ff = fill_nan(a, method="ff")
    assert np.isnan(ff[0]) and ff[2] == 1.0
    bf = fill_nan(a, method="bf")
    assert bf[0] == 1.0 and bf[2] == 2.0
    both = fill_nan(a, method="both")
    assert not np.isnan(both).any()


def test_fill_nan_2d_axis0_and_axis1():
    m = np.array([[1.0, np.nan, 3.0], [np.nan, 5.0, np.nan]])
    by_col = fill_nan(m, method="both", axis=0)
    assert np.allclose(by_col, [[1.0, 5.0, 3.0], [1.0, 5.0, 3.0]])
    by_row = fill_nan(m, method="both", axis=1)
    assert np.allclose(by_row, [[1.0, 1.0, 3.0], [5.0, 5.0, 5.0]])


def test_fill_nan_validation():
    with pytest.raises(ValueError):
        fill_nan(np.zeros((2, 2, 2)))
    with pytest.raises(ValueError):
        fill_nan(np.array([1.0]), method="sideways")
    with pytest.raises(ValueError):
        fill_nan(np.array([[1.0]]), axis=3)


# ------------------------------ drop_nan_in ---------------------------


def test_drop_nan_in_aligns_predictions():
    yt = np.array([1.0, 2.0, np.nan, 4.0])
    yp = np.array([0.9, 1.8, 3.1, 4.2])
    with pytest.warns(UserWarning):
        yt_f, yp_f = drop_nan_in(yt, yp, error="warn")
    assert np.allclose(yt_f, [1.0, 2.0, 4.0])
    assert np.allclose(yp_f, [0.9, 1.8, 4.2])


def test_drop_nan_in_policies():
    yt = np.array([1.0, np.nan])
    yp = np.array([1.0, 2.0])
    with pytest.raises(ValueError):
        drop_nan_in(yt, yp)  # default error='raise'
    with pytest.raises(ValueError):
        drop_nan_in(yt, yp, nan_policy="raise")
    same_t, same_p = drop_nan_in(yt, yp, nan_policy="propagate")
    assert same_t.size == 2 and same_p.size == 2
    om_t, om_p = drop_nan_in(yt, yp, nan_policy="omit")
    assert om_t.size == om_p.size == 1
    with pytest.raises(ValueError):
        drop_nan_in(yt, yp, nan_policy="bogus")


def test_drop_nan_in_shape_mismatch():
    with pytest.raises(ValueError):
        drop_nan_in(np.ones(3), np.ones(4))
