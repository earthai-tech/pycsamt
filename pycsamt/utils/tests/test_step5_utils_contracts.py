from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.utils.cleaner import (
    drop_constant_columns,
    fill_nan,
    impute_missing,
    ismissing,
    sanitize_frame_cols,
)
from pycsamt.utils.conversion import (
    convert,
    convert_temperature,
    convert_time,
    convert_value,
)


def test_metric_conversion_dispatch_and_unit_validation():
    assert convert_value("3.5km", target_unit="m") == pytest.approx(3500.0)
    assert convert_value("1200g", target_unit="kg") == pytest.approx(1.2)
    assert convert("20mm", "mm", "m") == pytest.approx(0.02)
    assert convert(1.2345, "km", "m", round_result=1) == pytest.approx(1234.5)

    with pytest.raises(ValueError, match="Incompatible units"):
        convert_value("10m", target_unit="g")


def test_temperature_and_time_conversions():
    assert convert_temperature("32F", unit_to="C") == pytest.approx(0.0)
    assert convert_temperature(273.15, unit_from="K", unit_to="C") == pytest.approx(0.0)
    assert convert(100, "C", "F", category="temperature") == pytest.approx(212.0)

    assert convert_time("2h", unit_to="min") == pytest.approx(120.0)
    assert convert("3600s", "s", "h") == pytest.approx(1.0)

    with pytest.raises(ValueError, match="Unsupported unit_to"):
        convert_time(1, unit_from="s", unit_to="fortnight")


def test_sanitize_frame_cols_preserves_data_and_handles_series_and_lists():
    df = pd.DataFrame({" A-1 ": [1, 2], "B#2": [3, 4]})

    cleaned = sanitize_frame_cols(df, case="lower")

    assert list(cleaned.columns) == ["a1", "b2"]
    pd.testing.assert_frame_equal(cleaned, pd.DataFrame({"a1": [1, 2], "b2": [3, 4]}))
    assert list(df.columns) == [" A-1 ", "B#2"]

    series = pd.Series([1, 2], name=" Site Name ")
    assert (
        sanitize_frame_cols(series, fill_pattern="_", case="upper").name == "SITE_NAME"
    )
    assert sanitize_frame_cols([" A ", "B-C"], case="lower") == ["a", "bc"]


def test_drop_constant_columns_and_rows():
    df = pd.DataFrame(
        {
            "constant": [1, 1, 1],
            "mostly": [2, 2, 3],
            "varying": [1, 2, 3],
        }
    )

    assert list(drop_constant_columns(df).columns) == ["mostly", "varying"]
    assert list(drop_constant_columns(df, threshold=0.8).columns) == [
        "mostly",
        "varying",
    ]
    assert list(drop_constant_columns(df, threshold=0.6).columns) == ["varying"]

    rows = pd.DataFrame([[1, 1, 1], [1, 2, 3], [np.nan, np.nan, np.nan]])
    out = drop_constant_columns(rows, axis="rows")
    assert out.index.tolist() == [1]


def test_impute_missing_strategies_and_errors():
    df = pd.DataFrame({"a": [1.0, np.nan, 3.0], "b": ["x", None, "x"]})

    meaned = impute_missing(df, strategy="mean")
    assert meaned.loc[1, "a"] == pytest.approx(2.0)
    assert pd.isna(meaned.loc[1, "b"])

    constant = impute_missing(
        df, strategy="constant", fill_value="missing", columns=["b"]
    )
    assert constant.loc[1, "b"] == "missing"

    with pytest.raises(KeyError, match="not found"):
        impute_missing(df, columns=["missing"])
    with pytest.raises(ValueError, match="Unknown strategy"):
        impute_missing(df, strategy="nearest")


def test_fill_nan_and_ismissing_contracts():
    np.testing.assert_allclose(
        fill_nan([np.nan, 1.0, np.nan, 2.0], method="both"),
        np.array([1.0, 1.0, 1.0, 2.0]),
    )
    np.testing.assert_allclose(
        fill_nan([[np.nan, 1.0], [2.0, np.nan]], method="backward"),
        np.array([[1.0, 1.0], [2.0, np.nan]]),
        equal_nan=True,
    )

    filled, missing = ismissing([0, 1, 2, 3], [0, 2], fill_value=-1)
    np.testing.assert_array_equal(filled, np.array([0, -1, 2, -1]))
    np.testing.assert_array_equal(missing, np.array([1, 3]))

    _, mask = ismissing([0, 1, 2], [1], return_index="mask")
    np.testing.assert_array_equal(mask, np.array([True, False, True]))
