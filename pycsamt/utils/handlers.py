# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

import inspect
import re
import warnings
from collections.abc import Iterable
from typing import Optional, Union

from .text import str2columns


def columns_manager(
    columns: Optional[Union[str, list, tuple]],
    default: Optional[list] = None,
    regex: Optional[re.Pattern] = None,
    pattern: Optional[str] = r"[@&,;#]",
    separator: Optional[str] = None,
    to_upper: bool = False,
    empty_as_none: bool = ...,
    to_string: bool = False,
    wrap_dict: bool = False,
    error: str = "raise",
) -> list:
    r"""
    A function to handle various types of column inputs, convert them
    into a list, and optionally process them based on additional parameters
    like converting to uppercase, handling empty values, or ensuring all items
    are strings.

    Parameters
    ----------
    columns : str, list, tuple, or None
        The input column names, which can be:
        - A string: treated as a list of column names split by a separator
          or regex.
        - A list or tuple: directly converted to a list if not already.
        - None: returns the default list or an empty list
        (if `empty_as_none` is False).

    default : list, optional
        Default list of columns to return if `columns` is None.

    regex : re.Pattern, optional
        A custom compiled regular expression to use for splitting string input.
        If not provided, the `pattern` parameter will be used.

    pattern : str, optional, default=r'[@&,;#]'
        The default regex pattern used to split the `columns` string if no `regex`
        is provided.

    separator : str, optional
        If `columns` is a string, this defines the separator used to split the string
        into a list of column names.

    to_upper : bool, default=False
        If True, converts all column names to uppercase.

    empty_as_none : bool, default=True
        If True, returns `None` when `columns` is empty or None. If False, an
        empty list is returned.

    to_string : bool, default=False
        If True, converts all items in `columns` to strings.

    error : str, default='warn'
        Specifies how to handle errors:
        - 'warn': issues a warning if any error occurs.
        - 'raise': raises an exception.
        - 'ignore': silently ignores any errors.

    Returns
    -------
    list
        A list of column names after processing.

    Example
    -------
    >>> from pycsamt.utils.handlers import columns_manager
    >>> columns_manager("col1, col2, col3", separator=",")
    ['col1', 'col2', 'col3']

    >>> columns_manager(['col1', 'col2', 'col3'], to_upper=True)
    ['COL1', 'COL2', 'COL3']
    """
    # Handle None input
    if columns is None:
        return (
            default
            if default is not None
            else (None if empty_as_none else [])
        )

    # Handle case where a single numeric value is passed, convert it to list
    if isinstance(columns, (int, float)):
        columns = [columns]

    elif callable(columns) or (isinstance(columns, dict) and wrap_dict):
        columns = [columns]

    ## Use inspect to determine if it is a class.
    # Alternatively, if the object is not iterable (has no __iter__ attribute),
    # we assume it's a single model instance.
    if inspect.isclass(columns) or not hasattr(columns, "__iter__"):
        columns = [columns]

    # If columns is a string, split by separator or use regex
    elif isinstance(columns, str):
        if separator is not None:
            columns = columns.split(separator)
        else:
            columns = str2columns(columns, regex=regex, pattern=pattern)

    # If columns is any iterable object, convert it to a list
    elif isinstance(columns, Iterable):
        try:
            columns = list(columns)
        except Exception as e:
            if error == "raise":
                raise ValueError("Error converting columns to list") from e
            elif error == "warn":
                warnings.warn(
                    f"Could not convert columns to list: {e}", stacklevel=2
                )
            else:
                pass  # Ignore errors silently

    # Ensure columns is a list at this point
    if isinstance(columns, list):
        if to_upper:
            # Check whether all items are strings before calling 'upper'
            if all(isinstance(col, str) for col in columns):
                columns = [col.upper() for col in columns]
            elif error == "raise":
                raise TypeError(
                    "All column names must be strings to convert to uppercase."
                )
            elif error == "warn":
                warnings.warn(
                    "Warning: Not all column names are strings,"
                    " skipping 'upper' conversion.",
                    stacklevel=2,
                )

        # Convert all items to string if requested
        if to_string:
            columns = [str(col) for col in columns]
    else:
        # If 'columns' is not a string, list, or tuple,
        # then it might be a single object
        # (for example, an instance of RandomForestRegressor).
        # In such a case, we attempt to check if it is iterable.
        # Since an instance of RandomForestRegressor
        # is neither callable nor a class, nor is it iterable
        # (i.e., it has no __iter__ # attribute), we wrap it into a list.
        if not isinstance(columns, (str, list, tuple)):
            try:
                iter(columns)
            except Exception:
                columns = [columns]

    return columns
