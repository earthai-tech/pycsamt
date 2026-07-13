# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
pycsamt.compat.pandas
---------------------

Compatibility helpers for pandas operations to ensure robustness
across different versions or usage patterns.
"""

from __future__ import annotations

from typing import Any, Callable

import pandas as pd


def safe_series_map(
    df: pd.DataFrame,
    col: str,
    mapper: Callable[[Any], Any] | dict[Any, Any],
) -> pd.Series:
    """
    Safely apply a mapping function or dict to a DataFrame column.

    This function ensures the target column is a pandas Series
    before applying the `.map()` method. It acts as a robust
    compatibility layer, preventing errors like "'DataFrame'
    object has no attribute 'map'".

    Parameters
    ----------
    df : pd.DataFrame
        The source DataFrame.
    col : str
        The name of the column (Series) to which the mapping
        should be applied.
    mapper : callable or dict
        The mapping function or dictionary to apply to the Series.

    Returns
    -------
    pd.Series
        A new Series containing the mapped values.

    Raises
    ------
    KeyError
        If the specified column does not exist in the DataFrame.
    TypeError
        If the selected column object is not a pandas Series.
    """
    if col not in df.columns:
        raise KeyError(f"Column '{col}' not found in DataFrame.")

    series_obj = df[col]
    series: pd.Series
    # Handle the case where a selection might return a
    # single-column DataFrame instead of a Series.
    if isinstance(series_obj, pd.DataFrame):
        if series_obj.shape[1] == 1:
            series = series_obj.iloc[:, 0]
        else:
            # Mapping is ambiguous for a multi-column DataFrame
            raise TypeError(
                f"Selection for '{col}' returned a multi-column "
                "DataFrame. Map is only for a single Series."
            )
    else:
        series = series_obj

    if not isinstance(series, pd.Series):
        raise TypeError(
            f"Object for column '{col}' is not a Series, "
            f"but a {type(series).__name__}."
        )
    return series.map(mapper)
