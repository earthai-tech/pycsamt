# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Additional unit tests for :mod:`pycsamt.utils.cleaner`."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest

from pycsamt.utils.cleaner import (
    fill_nan,
    impute_missing,
    ismissing,
    sanitize_frame_cols,
)


# -------------------------- sanitize_frame_cols --------------------------


def test_sanitize_frame_cols_precompiled_regex():
    splitter = re.compile(r"[\s_]+")
    out = sanitize_frame_cols(["a b", "c_d"], regex=splitter)
    assert out == ["ab", "cd"]


def test_sanitize_frame_cols_strip_false_keeps_whitespace_padding():
    out = sanitize_frame_cols([" a "], strip=False, regex=re.compile(r"x"))
    assert out == [" a "]


def test_sanitize_frame_cols_applies_custom_func():
    out = sanitize_frame_cols(["ab"], func=str.upper)
    assert out == ["AB"]


def test_sanitize_frame_cols_case_upper():
    out = sanitize_frame_cols(["ab"], case="upper")
    assert out == ["AB"]


def test_sanitize_frame_cols_dataframe_inplace():
    df = pd.DataFrame({"A B": [1], "C-D": [2]})
    result = sanitize_frame_cols(df, case="lower", inplace=True)
    assert result is df
    assert list(df.columns) == ["ab", "cd"]


def test_sanitize_frame_cols_series_inplace():
    ser = pd.Series([1, 2], name="A B")
    result = sanitize_frame_cols(ser, case="lower", inplace=True)
    assert result is ser
    assert ser.name == "ab"


def test_sanitize_frame_cols_rejects_unconvertible_input():
    class NotIterable:
        def __iter__(self):
            raise RuntimeError("not iterable after all")

    with pytest.raises(TypeError, match="must be DataFrame, Series"):
        sanitize_frame_cols(NotIterable())


# ----------------------------- impute_missing -----------------------------


def test_impute_missing_skips_non_numeric_for_mean_and_median():
    df = pd.DataFrame({"num": [1.0, np.nan, 3.0], "txt": ["a", None, "c"]})
    out = impute_missing(df, strategy="mean", columns=["num", "txt"])
    assert out["txt"].isna().sum() == 1  # untouched, still has the NaN
    assert out["num"].isna().sum() == 0


def test_impute_missing_median_strategy():
    df = pd.DataFrame({"A": [1.0, np.nan, 3.0]})
    out = impute_missing(df, strategy="median")
    assert out["A"].tolist() == [1.0, 2.0, 3.0]


def test_impute_missing_mode_strategy_with_and_without_mode():
    df = pd.DataFrame({"A": [1.0, 1.0, np.nan]})
    out = impute_missing(df, strategy="mode", columns=["A"])
    assert out["A"].iloc[-1] == 1.0

    empty_df = pd.DataFrame({"A": [np.nan, np.nan]})
    out2 = impute_missing(empty_df, strategy="mode", columns=["A"])
    assert out2["A"].isna().all()  # no mode exists -> stays NaN


# -------------------------------- fill_nan --------------------------------


def test_fill_nan_method_aliases():
    a = np.array([np.nan, 1.0, np.nan])
    assert np.allclose(fill_nan(a, method="forward")[1:], [1.0, 1.0])
    assert np.allclose(fill_nan(a, method="fwd")[1:], [1.0, 1.0])
    b = np.array([np.nan, np.nan, 1.0])
    assert np.allclose(fill_nan(b, method="backward"), [1.0, 1.0, 1.0])
    assert np.allclose(fill_nan(b, method="bwd"), [1.0, 1.0, 1.0])


def test_fill_nan_ff_mode_1d_and_2d():
    assert np.allclose(
        fill_nan(np.array([np.nan, 1.0, np.nan]), method="ff")[1:], [1.0, 1.0],
    )
    out = fill_nan(np.array([[np.nan, 1.0], [2.0, np.nan]]), method="ff")
    assert out.shape == (2, 2)


# -------------------------------- ismissing --------------------------------


def test_ismissing_return_index_true_and_index_aliases():
    for token in ("true", "index", "ix"):
        _, missing = ismissing([0, 1, 2, 3], [0, 2], return_index=token)
        assert missing == [1, 3]


def test_ismissing_rejects_bad_return_index():
    with pytest.raises(ValueError, match="Invalid return_index"):
        ismissing([0, 1], [0], return_index="bogus")
