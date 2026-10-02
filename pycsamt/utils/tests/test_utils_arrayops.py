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


def test_concat_rejects_non_iterable_list_of_arrays():
    with pytest.raises(TypeError, match="must be iterable"):
        concat_array_from_list(42)


def test_concat_rejects_non_array_like_item():
    class Unconvertible:
        def __array__(self, dtype=None):
            raise RuntimeError("cannot convert")

    with pytest.raises(TypeError, match="is not array-like"):
        concat_array_from_list([Unconvertible()])


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


def test_is_iterable_wraps_string_parsing_failure(monkeypatch):
    from pycsamt.utils import arrayops as arrayops_mod

    def _boom(_text):
        raise RuntimeError("parse failure")

    monkeypatch.setattr(arrayops_mod, "str2columns", _boom)
    with pytest.raises(TypeError, match="Error parsing string"):
        is_iterable("a b", transform=True, parse_string=True)


def test_is_iterable_transform_falls_back_when_list_conversion_fails():
    class BrokenIterable:
        def __iter__(self):
            def _gen():
                yield 1
                raise RuntimeError("boom mid-iteration")

            return _gen()

    obj = BrokenIterable()
    assert is_iterable(obj, transform=True) == [obj]


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


def test_reshape_2d_axis0_and_axis1_targeted():
    col = np.arange(4).reshape(4, 1)  # n=4, m=1
    row = np.arange(4).reshape(1, 4)  # n=1, m=4

    # axis=0 with m already 1 -> reshape(n, 1) branch
    assert reshape(col, axis=0).shape == (4, 1)
    # axis=1 with n already 1 -> reshape(1, m) branch
    assert reshape(row, axis=1).shape == (1, 4)
    # axis=1 with n != 1 -> passthrough branch
    assert reshape(col, axis=1).shape == (4, 1)
    # axis=0 with m != 1 -> passthrough branch (via a row vector)
    assert reshape(row, axis=0).shape == (1, 4)


def test_reshape_invalid_inputs():
    with pytest.raises(ValueError):
        reshape(np.zeros((2, 2, 2)))
    with pytest.raises(ValueError):
        reshape(np.arange(3), axis=2)


# ------------------------------- frameify -----------------------------


def test_frameify_array_with_columns_and_types():
    arr = np.array([[1.0, "x"], [2.0, "y"]], dtype=object)
    df, nf, cf = frameify(arr, columns=["num", "cat"], return_feature_types=True)
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


def test_frameify_pop_cat_features_verbose_with_and_without_cat(capsys):
    numeric_only = pd.DataFrame({"n": [1, 2]})
    frameify(numeric_only, pop_cat_features=True, verbose=True)
    assert "does not contain any categorial" in capsys.readouterr().out

    df_in = pd.DataFrame({"n": [1, 2], "c": ["u", "v"]})
    frameify(df_in, pop_cat_features=True, verbose=True)
    assert capsys.readouterr().out  # some listing was printed


def test_frameify_dataframe_input_replaces_columns_verbose(capsys):
    df_in = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]})
    df = frameify(df_in, columns=["x", "y"], verbose=True)
    assert list(df.columns) == ["x", "y"]
    assert "Columns should be replaced" in capsys.readouterr().out


def test_frameify_dataframe_input_replaces_columns_quiet(capsys):
    df_in = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]})
    df = frameify(df_in, columns=["x", "y"])  # verbose=False (default)
    assert list(df.columns) == ["x", "y"]
    assert capsys.readouterr().out == ""


def test_frameify_drop_nan_columns_disabled_keeps_all_nan_column():
    df_in = pd.DataFrame({"a": [1.0, 2.0], "b": [np.nan, np.nan]})
    df = frameify(df_in, drop_nan_columns=False)
    assert "b" in df.columns


def test_frameify_sanitize_columns_skips_already_numeric_headers():
    # Purely integer column labels have nothing to sanitize, so
    # sanitize_columns=True leaves them untouched (no str conversion).
    df_in = pd.DataFrame(np.ones((2, 2)))  # integer column labels 0, 1
    df = frameify(df_in, sanitize_columns=True)
    assert list(df.columns) == [0, 1]


def test_frameify_sanitize_columns_cleans_string_headers():
    df_in = pd.DataFrame({"a b": [1.0], "c-d": [2.0]})
    df = frameify(df_in, sanitize_columns=True)
    assert list(df.columns) == ["a_b", "c_d"]


def test_frameify_drop_nan_columns_verbose_reports_found_and_none(capsys):
    with_nan = pd.DataFrame({"a": [1.0, 2.0], "b": [np.nan, np.nan]})
    frameify(with_nan, verbose=True)
    assert "NaN columns found" in capsys.readouterr().out

    without_nan = pd.DataFrame({"a": [1.0, 2.0]})
    frameify(without_nan, verbose=True)
    assert "No NaN column found" in capsys.readouterr().out


def test_frameify_how_not_all_keeps_partially_nan_rows():
    df_in = pd.DataFrame({"a": [1.0, np.nan], "b": [np.nan, np.nan]})
    df = frameify(df_in, how="any")
    # column "b" is dropped (all-NaN); row 1 (now all-NaN in "a" only)
    # is kept since how != "all" skips the row-drop step entirely.
    assert "b" not in df.columns
    assert len(df) == 2


def test_frameify_reset_index():
    df_in = pd.DataFrame({"a": [1.0, 2.0]}, index=[5, 6])
    df = frameify(df_in, reset_index=True)
    assert list(df.index) == [0, 1]


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
    df = assert_xy_in([1, 2], [3, 4], to_frame=True, columns=["east", "north"])
    assert list(df.columns) == ["east", "north"]


def test_assert_xy_in_dropna_and_numeric():
    x, y = assert_xy_in([1.0, np.nan, 3.0], [4.0, 5.0, 6.0], dropna=True)
    assert x.size == y.size == 2
    x2, _ = assert_xy_in(["1", "2"], ["3", "4"], xy_numeric=True)
    assert x2.dtype.kind in "if"


def test_assert_xy_in_accepts_series_directly():
    x_in = pd.Series([1.0, 2.0], name="x")
    y_in = pd.Series([3.0, 4.0], name="y")
    x, y = assert_xy_in(x_in, y_in, asarray=False)
    assert x is x_in
    assert y is y_in


def test_assert_xy_in_wraps_unseriesable_scalar():
    # A 2-D array cannot become a 1-D pandas Series directly; the
    # fallback wraps it as a single-element Series instead of raising.
    x, y = assert_xy_in(np.ones((2, 2)), [1], asarray=False)
    assert len(x) == 1
    assert len(y) == 1


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


def test_interpolate_grid_accepts_plain_list_without_dunder_array():
    arr = [[1.0, np.nan], [np.nan, 4.0]]  # plain list, no __array__
    out = interpolate_grid(arr)
    assert not np.isnan(out).any()


def test_interpolate_grid_sparse_input_auto_fill():
    # Only 2 valid points (< 4): triggers the fill-only fallback with
    # the default fill_value="auto" (forward/backward fill).
    arr = np.array([[1.0, np.nan], [np.nan, 4.0]])
    out = interpolate_grid(arr)  # fill_value defaults to "auto"
    assert not np.isnan(out).any()


def test_interpolate_grid_sparse_1d_input_auto_fill():
    arr = np.array([1.0, np.nan, np.nan, np.nan])  # < 4 valid points, 1D
    out = interpolate_grid(arr)
    assert out.ndim == 1
    assert not np.isnan(out).any()


def test_interpolate_grid_dense_constant_fill_value():
    # >= 4 valid points: exercises the real griddata path with a
    # non-"auto" fill_value for the remaining NaNs at the edges.
    x = [28, np.nan, 50, 60]
    y = [np.nan, 1000, 2000, 3000]
    xy = np.vstack((x, y)).T
    out = interpolate_grid(xy, fill_value=0.0)
    assert not np.isnan(out).any()


def test_interpolate_grid_view_plots_without_error():
    import matplotlib
    matplotlib.use("Agg")
    x = [28, np.nan, 50, 60]
    y = [np.nan, 1000, 2000, 3000]
    xy = np.vstack((x, y)).T
    out = interpolate_grid(xy, view=True)
    assert not np.isnan(out).any()


# ------------------------------- _fill_nan ------------------------------


def test_private_fill_nan_1d_input_and_every_method():
    from pycsamt.utils.arrayops import _fill_nan

    a = np.array([np.nan, 1.0, np.nan, 2.0])
    ff = _fill_nan(a, method="ff")
    assert ff[2] == 1.0
    bf = _fill_nan(a, method="bf")
    assert np.isnan(bf).sum() == 0 or bf[0] == 1.0
    both = _fill_nan(a, method="both")
    assert not np.isnan(both).any()


def test_private_fill_nan_rejects_bad_ndim_and_method():
    from pycsamt.utils.arrayops import _fill_nan

    with pytest.raises(ValueError, match="only 1D or 2D"):
        _fill_nan(np.zeros((2, 2, 2)))
    with pytest.raises(ValueError, match="Unknown method"):
        _fill_nan(np.array([1.0]), method="sideways")


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


def test_drop_nan_in_no_nan_present_skips_error_handling():
    yt = np.array([1.0, 2.0, 3.0])
    yp = np.array([1.1, 2.1, 3.1])
    yt_f, yp_f = drop_nan_in(yt, yp)  # default error="raise", but no NaNs
    assert np.allclose(yt_f, yt)
    assert np.allclose(yp_f, yp)


def test_drop_nan_in_rejects_invalid_error_value():
    yt = np.array([1.0, np.nan])
    yp = np.array([1.0, 2.0])
    with pytest.raises(ValueError, match="error must be one of"):
        drop_nan_in(yt, yp, error="bogus")
