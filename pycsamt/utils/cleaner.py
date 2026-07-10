"""
Cleaning utilities for DataFrame and Series columns.
"""

import re
from collections.abc import Sequence
from typing import (
    Any,
    Callable,
    Optional,
    Union,
)

import numpy as np
import pandas as pd

from ..decorators import check_empty, isdf

__all__ = [
    'sanitize_frame_cols',
    'drop_constant_columns',
    'impute_missing',
    'fill_nan',
    'ismissing'
]

@check_empty
def sanitize_frame_cols(
    df_or_cols: Union[pd.DataFrame, pd.Series, Sequence[str]],
    *,
    regex: Optional[Union[str, re.Pattern]] = None,
    pattern: Optional[str] = None,
    fill_pattern: str = '',
    func: Optional[Callable[[str], str]] = None,
    strip: bool = True,
    case: Optional[str] = None,
    inplace: bool = False
) -> Union[pd.DataFrame, pd.Series, list]:
    """
    Sanitize column names or a list of strings by removing
    unwanted characters and applying optional transforms.

    Can operate on DataFrame, Series.name, or a list of names.

    Parameters
    ----------
    df_or_cols : DataFrame, Series, or list of str
        If DataFrame, its columns are sanitized. If Series,
        its name is sanitized. If list, list items are sanitized.
    regex : str or Pattern, optional
        Precompiled regex or pattern string for unwanted chars.
    pattern : str, optional
        Fallback pattern if ``regex`` not provided.
    fill_pattern : str, default ''
        Replacement for matches of the pattern.
    func : callable, optional
        Custom function applied to each name after regex step.
    strip : bool, default True
        Strip whitespace from ends of names before processing.
    case : {'lower','upper',None}, default None
        Convert names to lowercase or uppercase.
    inplace : bool, default False
        Modify DataFrame or Series in place.

    Returns
    -------
    DataFrame, Series, or list
        Sanitized object, matching input type.
    """
    # Determine regex to use
    if isinstance(regex, re.Pattern):
        splitter = regex
    else:
        pat = regex if isinstance(regex, str) else (
            pattern or r'[_#&.\)\(*@!_,;\s-]+'
        )
        splitter = re.compile(pat)

    def _clean(name: Any) -> str:
        s = str(name)
        if strip:
            s = s.strip()
        # remove unwanted
        s = splitter.sub(fill_pattern, s)
        # apply user func
        if func is not None and callable(func):
            s = func(s)
        # case transform
        if case == 'lower':
            s = s.lower()
        elif case == 'upper':
            s = s.upper()
        return s

    # Process DataFrame columns
    if isinstance(df_or_cols, pd.DataFrame):
        df = df_or_cols
        new_cols = [_clean(c) for c in df.columns]
        if inplace:
            df.columns = new_cols
            return df
        else:
            return pd.DataFrame(
                df.values, columns=new_cols,
                index=df.index
            )
    # Process Series name
    if isinstance(df_or_cols, pd.Series):
        ser = df_or_cols
        new_name = _clean(ser.name)
        if inplace:
            ser.name = new_name
            return ser
        else:
            return pd.Series(
                data=ser.values,
                index=ser.index,
                name=new_name,
                dtype=ser.dtype
            )
    # Process list of strings
    try:
        seq = list(df_or_cols)
    except Exception:
        raise TypeError(
            'Input must be DataFrame, Series, or list of names'
        )
    return [_clean(item) for item in seq]


@isdf
def drop_constant_columns(
    df: pd.DataFrame,
    *,
    axis: str = 'columns',
    threshold: float = 1.0,
    inplace: bool = False
) -> pd.DataFrame:
    """
    Drop constant columns or rows from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame to prune.
    axis : {'columns', 'rows'}, default='columns'
        Axis along which to drop constant values.
    threshold : float, default=1.0
        Proportion of identical values to consider constant.
        Must be in (0, 1].
    inplace : bool, default=False
        If True, modify `df` in place.

    Returns
    -------
    pd.DataFrame
        DataFrame with constant axis entries removed.

    Notes
    -----
    - A column (or row) is dropped if the most frequent value
      accounts for >= `threshold` of all entries.
    - NaNs count as a distinct value by default.

    Examples
    --------
    >>> import pandas as pd
    >>> from pycsamt.utils.cleaner import drop_constant_columns
    >>> df = pd.DataFrame({'A':[1,1,1], 'B':[1,2,3]})
    >>> drop_constant_columns(df)
       B
    0  1
    1  2
    2  3
    """
    # Validate parameters
    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f"`df` must be a DataFrame, got {type(df)!r}"
        )
    axis = axis.lower()
    if axis not in ('columns', 'rows'):
        raise ValueError(
            "axis must be 'columns' or 'rows'"
        )
    if not (0 < threshold <= 1):
        raise ValueError(
            "threshold must be within (0, 1], got {threshold}"
        )

    # Prepare DataFrame copy
    df2 = df if inplace else df.copy()

    if axis == 'columns':
        # Identify columns to drop
        drop_cols = []
        for col in df2.columns:
            vc = df2[col].value_counts(
                normalize=True, dropna=False
            )
            freq = vc.iloc[0] if not vc.empty else 0
            if freq >= threshold:
                drop_cols.append(col)
        # Drop columns
        df2.drop(columns=drop_cols, inplace=True)
    else:
        # Identify rows to drop
        mask = df2.apply(
            lambda row: (
                row.value_counts(
                    normalize=True, dropna=False
                ).iloc[0] >= threshold
            ),
            axis=1
        )
        # Keep non-constant rows
        df2 = df2.loc[~mask]

    return df2

@isdf
def impute_missing(
    df: pd.DataFrame,
    *,
    strategy: str = 'mean',
    fill_value: Any = None,
    columns: Optional[Sequence[str]] = None,
    inplace: bool = False
) -> pd.DataFrame:
    """
    Impute missing values in specified DataFrame columns.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing NaNs.
    strategy : {'mean', 'median', 'mode', 'constant'},
        Imputation strategy to apply.
    fill_value : Any, default=None
        Constant used when strategy='constant'.
    columns : list of str, optional
        Subset of columns to impute; defaults to all numeric columns.
    inplace : bool, default=False
        If True, modify `df` in place.

    Returns
    -------
    pd.DataFrame
        DataFrame with missing values filled.

    Notes
    -----
    - 'mode' uses the first mode; if no mode, leaves NaN.
    - Non-numeric columns are skipped for numeric strategies.

    Examples
    --------
    >>> import pandas as pd
    >>> from pycsamt.utils.cleaner import impute_missing
    >>> df = pd.DataFrame({
    ...     'A': [1, None, 3],
    ...     'B': [None, 2, 2]
    ... })
    >>> impute_missing(df, strategy='median')
         A    B
    0  1.0  2.0
    1  2.0  2.0
    2  3.0  2.0
    """
    # Validate inputs
    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f"`df` must be DataFrame, got {type(df)!r}"
        )
    strat = strategy.lower()
    valid_strat = {'mean', 'median', 'mode', 'constant'}
    if strat not in valid_strat:
        raise ValueError(
            f"Unknown strategy: {strategy!r}; must be one of"
            f" {valid_strat}"
        )
    # Prepare DataFrame
    df2 = df if inplace else df.copy()
    # Determine columns to process
    cols = columns or df2.select_dtypes(
        include=[np.number]
    ).columns.tolist()
    for col in cols:
        if col not in df2.columns:
            raise KeyError(
                f"Column {col!r} not found in DataFrame"
            )
        series = df2[col]
        # Compute fill value
        if strat == 'mean':
            val = series.mean()
        elif strat == 'median':
            val = series.median()
        elif strat == 'mode':
            modes = series.mode(dropna=True)
            val = modes.iloc[0] if not modes.empty else np.nan
        else:  # constant
            val = fill_value
        # Fill missing
        series.fillna(val, inplace=True)
    return df2

@check_empty(['arr'])
def fill_nan(
    arr: Union[Sequence[Any], np.ndarray],
    *,
    method: str = 'ff'
) -> np.ndarray:
    """
    Efficiently forward/backward fill NaNs in array.

    Parameters
    ----------
    arr : sequence or ndarray
        Input 1D or 2D array with NaNs.
    method : {'ff','bf','both'}, default 'ff'
        'ff': forward-fill; 'bf': backward-fill;
        'both': apply ff then bf.

    Returns
    -------
    filled : ndarray
        Array with NaNs filled. Preserves shape.

    Notes
    -----
    - Works row-wise for 2D arrays.
    - 'ff' fills NaNs after valid values;
      'bf' fills before valid values.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.utils.arrayops import fill_nan
    >>> a = np.array([np.nan,1,np.nan,2])
    >>> fill_nan(a, method='ff')
    array([nan, 1., 1., 2.])
    >>> fill_nan(a, method='both')
    array([1., 1., 1., 2.])
    """

    # convert to float array
    a = np.array(arr, dtype=float)
    # handle dimensions
    if a.ndim == 1:
        a = a.reshape(1, -1)
        orig1d = True
    elif a.ndim == 2:
        orig1d = False
    else:
        raise ValueError(
            f"Only 1D/2D arrays supported; got ndim={a.ndim}"
        )
    # normalize method
    m = method.lower().strip()
    if m in ('forward','ff','fwd'):
        m = 'ff'
    elif m in ('backward','bf','bwd'):
        m = 'bf'
    elif m in ('both','fb','bff','ffbf'):
        m = 'both'
    else:
        raise ValueError(
            f"Unknown method {method!r}; use 'ff','bf','both'"
        )
    # define forward fill
    def _ffill(row):
        mask = np.isnan(row)
        idx = np.where(~mask, np.arange(row.size), 0)
        np.maximum.accumulate(idx, out=idx)
        return row[idx]
    # define backward fill
    def _bfill(row):
        mask = np.isnan(row)
        idx = np.where(~mask, np.arange(row.size), row.size-1)
        idx = np.minimum.accumulate(idx[::-1])[::-1]
        return row[idx]
    # apply along axis=1
    if m == 'ff':
        out = np.vstack(_ffill(r) for r in a)
    elif m == 'bf':
        out = np.vstack(_bfill(r) for r in a)
    else:
        tmp = np.vstack(_ffill(r) for r in a)
        out = np.vstack(_bfill(r) for r in tmp)
    return out[0] if orig1d else out


@check_empty(['refarr', 'arr'])
def ismissing(
    refarr: Union[Sequence[Any], np.ndarray],
    arr: Union[Sequence[Any], np.ndarray],
    *,
    fill_value: Any = np.nan,
    return_index: Union[bool, str] = False
) -> Union[tuple[np.ndarray, np.ndarray],
           tuple[np.ndarray, Sequence[int]]]:
    """
    Identify missing entries of `arr` in `refarr` and fill.

    Parameters
    ----------
    refarr : sequence or ndarray
        Reference array with desired full length.
    arr : sequence or ndarray
        Array of present values (subset of refarr).
    fill_value : Any, default np.nan
        Value to insert where arr misses refarr.
    return_index : bool or {'values','index','mask'},
        If False or 'values', return missing values.
        If 'index', return indices of missing.
        If 'mask', return boolean mask.

    Returns
    -------
    filled : ndarray
        Array matching `refarr` where missing are filled.
    missing_info : ndarray or list
        Missing values or indices or mask.

    Notes
    -----
    - Assumes comparability via equality.
    - Duplicates in refarr propagate multiplicity.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.utils.arrayops import ismissing
    >>> ref = np.arange(5)
    >>> arr = np.array([0,2,4])
    >>> f, miss = ismissing(ref, arr)
    >>> f
    array([ 0., nan, 2., nan, 4.])
    >>> miss
    array([1., 3.])
    >>> _, idx = ismissing(ref, arr, return_index='index')
    >>> idx
    [1, 3]
    >>> _, mask = ismissing(ref, arr, return_index='mask')
    >>> mask
    array([False, True, False, True, False])
    """

    # convert inputs
    ref = np.array(refarr, copy=True)
    arr2 = np.array(arr)
    if ref.ndim != 1 or arr2.ndim != 1:
        raise ValueError(
            "Both refarr and arr must be 1D sequences"
        )
    # normalize return_index
    r = str(return_index).lower()
    if r in ('false','values','value'):
        mode = 'values'
    elif r in ('true','index','ix'):
        mode = 'index'
    elif r == 'mask':
        mode = 'mask'
    else:
        raise ValueError(
            f"Invalid return_index {return_index!r}"
        )
    # presence mask
    mask = np.isin(ref, arr2)
    # fill missing
    filled = np.where(mask, ref, fill_value)
    # missing data info
    if mode == 'values':
        missing = ref[~mask]
    elif mode == 'index':
        missing = list(np.nonzero(~mask)[0])
    else:  # mask
        missing = ~mask
    return filled, missing

