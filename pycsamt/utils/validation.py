# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Validation utilities for pycsamt.
"""

from __future__ import annotations

import inspect
import warnings
from collections.abc import Iterable, Sequence
from functools import wraps
from typing import (
    Any,
    Callable,
)

import numpy as np
import pandas as pd

from ..exceptions import NotFittedError, NotReadError

__all__ = [
    "check_is_fitted",
    "_assert_all_types",
    "_isin",
    "assert_ratio",
    "_validate_name_in",
    "_is_numeric_dtype",
    "check_consistency_size",
    "_is_arraylike_1d",
    "isinstance_relaxed",
    "has_read",
    "check_has_read",
    "isin",
    "isin_if",
    "ensure_n_items",
]


def ensure_n_items(
    *values: Any,
    n: int = 2,
    items: Any = None,
    expect: str = "mixed",
    coerce: bool = False,
    to_string: bool = False,
    allow_none: bool = False,
    allow_nan: bool = False,
    bounds: tuple[float, float] | None = None,
    unique: bool = False,
    return_as: str = "tuple",
    dtype: Any = None,
    name: str | None = None,
    error: str = "raise",
):
    r"""
    Parse and validate ``n`` elements from flexible input.

    Accepts either:
    - a single iterable of length ``n`` (not a string), or
    - ``n`` separate positional values (``*values``).

    Enforces type policy (numeric / string / mixed), optional
    coercion, and optional bounds. Returns the values in the
    requested container type.

    Parameters
    ----------
    *values : Any
        Positional values. If provided, ``items`` is ignored.
    n : int, default 2
        Number of elements expected.
    items : Any, optional
        A single iterable-of-``n`` alternative to ``*values``.
    expect : {"numeric","string","mixed"}, default "mixed"
        Type policy for the elements.
    coerce : bool, default False
        If ``True`` and ``expect="numeric"``, try to coerce
        to ``float``.
    to_string : bool, default False
        Convert elements to ``str`` (applies to any policy).
    allow_none : bool, default False
        If ``True``, ``None`` is permitted.
    allow_nan : bool, default False
        If ``True``, ``nan`` is permitted (numeric policy).
    bounds : (float, float), optional
        Numeric lower/upper bounds (inclusive).
    unique : bool, default False
        Require that all elements are distinct.
    return_as : {"tuple","list","array"}, default "tuple"
        Container type for the result.
    dtype : Any, optional
        Numpy dtype used when ``return_as="array"``.
    name : str, optional
        Label used in error messages.
    error : {"raise","warn","ignore"}, default "raise"
        Behavior on validation failures.

    Returns
    -------
    tuple | list | ndarray
        Validated values in the requested container.

    Notes
    -----
    - Strings are *not* treated as iterables for splitting.
    - ``bool`` is not accepted as numeric.
    """
    label = name or "values"

    # Collect inputs into a flat list of length n
    if values:
        seq = list(values)
    else:
        if _iter_not_str(items):
            seq = list(items)  # type: ignore[arg-type]
        else:
            seq = [items]

    if len(seq) != n:
        msg = f"Expected {n} {label}; got {len(seq)}."
        if error == "raise":
            raise ValueError(msg)
        if error == "warn":
            warnings.warn(msg, stacklevel=2)

    # If extra values are supplied, trim but warn.
    if len(seq) > n:
        msg = f"Got {len(seq)} {label}; truncating to {n}."
        if error == "raise":
            raise ValueError(msg)
        if error == "warn":
            warnings.warn(msg, stacklevel=2)
        seq = seq[:n]

    # If fewer, pad with None to reach n (warn/ignore only).
    if len(seq) < n:
        pad = [None] * (n - len(seq))
        seq = list(seq) + pad

    # Element-wise normalization
    out: list[Any] = []
    for i, x in enumerate(seq):
        # None handling
        if x is None:
            if not allow_none:
                msg = f"{label}[{i}] is None; not allowed."
                if error == "raise":
                    raise TypeError(msg)
                if error == "warn":
                    warnings.warn(msg, stacklevel=2)
            out.append(None)
            continue

        # string conversion (global)
        if to_string:
            out.append(str(x))
            continue

        # policy: string
        if expect == "string":
            if isinstance(x, (str, bytes, bytearray)):
                out.append(str(x))
            else:
                msg = (
                    f"{label}[{i}] must be string; got {type(x).__name__!r}."
                )
                if error == "raise":
                    raise TypeError(msg)
                if error == "warn":
                    warnings.warn(msg, stacklevel=2)
                out.append(str(x))
            continue

        # policy: numeric
        if expect == "numeric":
            try:
                v = (
                    _to_float(x) if coerce else float(x)  # may raise
                )
            except Exception as exc:
                msg = f"{label}[{i}] not numeric: {x!r}."
                if error == "raise":
                    raise TypeError(msg) from exc
                if error == "warn":
                    warnings.warn(msg, stacklevel=2)
                v = np.nan

            if not allow_nan and np.isnan(v):
                msg = f"{label}[{i}] is NaN; not allowed."
                if error == "raise":
                    raise ValueError(msg)
                if error == "warn":
                    warnings.warn(msg, stacklevel=2)

            if bounds is not None and np.isfinite(v):
                lo, hi = bounds
                if v < lo or v > hi:
                    msg = f"{label}[{i}]={v} outside [{lo}, {hi}]."
                    if error == "raise":
                        raise ValueError(msg)
                    if error == "warn":
                        warnings.warn(msg, stacklevel=2)

            out.append(v)
            continue

        # policy: mixed (no checks)
        out.append(x)

    # Uniqueness check
    if unique:
        # stringify for hashability with stable compare
        sig = tuple(map(repr, out))
        if len(set(sig)) != len(out):
            msg = f"{label} must be unique."
            if error == "raise":
                raise ValueError(msg)
            if error == "warn":
                warnings.warn(msg, stacklevel=2)

    # Materialize in requested container
    if return_as == "list":
        return list(out)
    if return_as == "array":
        return np.array(out, dtype=dtype)
    # default: tuple
    return tuple(out)


def _iter_not_str(x: Any) -> bool:
    """True if iterable but not a (byte)string."""
    if isinstance(x, (str, bytes, bytearray)):
        return False
    return hasattr(x, "__iter__")


def _to_float(x: Any) -> float:
    """Best-effort float coercion (keeps nan)."""
    if x is None:
        return np.nan
    if isinstance(x, (float, np.floating)):
        return float(x)
    if isinstance(x, bool):
        # avoid treating bool as numeric
        raise TypeError("bool is not accepted as numeric.")
    if isinstance(x, (int, np.integer)):
        return float(x)
    try:
        v = float(x)
        return v
    except Exception as exc:
        raise TypeError(f"Cannot coerce {x!r} to float.") from exc


def has_read(
    obj: Any = None,
    *,
    attributes: str | list[str] | None = None,
    msg: str | None = None,
) -> bool:
    r"""Check if an object has been populated with data.

    This function provides a standardized way to verify that an
    object's `read` method (or equivalent) has been called
    before other methods are accessed. It uses a multi-stage
    checking process for maximum flexibility.

    Parameters
    ----------
    obj : object, optional
        The instance of the class to check. If ``None``, the
        function will attempt to find the instance (`self`) from
        the calling method's frame.
    attributes : str, list of str, optional
        Attribute names to verify for existence and content. This
        serves as a fallback check if no explicit read-state flag
        is found on the object.
    msg : str, optional
        A custom error message to raise if validation fails. If
        not provided, a default message is generated.

    Returns
    -------
    bool
        ``True`` if the object is considered "read". The function
        does not return ``False`` but raises an error instead.

    Raises
    ------
    NotReadError
        If the object has not been read, based on one of the
        checking mechanisms.
    ValueError
        If `obj` is ``None`` and `self` cannot be found in the
        caller's scope.

    Notes
    -----
    The check is performed with the following priority:

    1. **Custom Method**: Looks for a `__has_read__` method on
       the object. If it exists and returns ``False``, an error
       is raised.
    2. **Boolean Flag**: Looks for a `_has_read` attribute. If
       it is explicitly ``False``, an error is raised.
    3. **Attribute Inspection**: If neither of the above is found,
       it checks the `attributes` argument to ensure the specified
       attributes exist and are not empty.

    Examples
    --------
    >>> from pycsamt.exceptions import NotReadError
    >>> class MyReader:
    ...     def __init__(self):
    ...         self.data = None
    ...     def read_data(self):
    ...         self.data = [1, 2, 3]
    ...     def process(self):
    ...         has_read(self, attributes='data')
    ...         return sum(self.data)
    ...
    >>> reader = MyReader()
    >>> try:
    ...     reader.process()
    ... except NotReadError as e:
    ...     print(e)
    'MyReader' has not been read yet. Attribute 'data' is ...
    >>> reader.read_data()
    >>> reader.process()
    6

    See Also
    --------
    check_has_read : A decorator that uses this function.
    """
    if obj is None:
        # Introspect the caller's frame to find 'self'
        frame = inspect.currentframe()
        if frame and frame.f_back:
            obj = frame.f_back.f_locals.get("self")
        if obj is None:
            raise ValueError(
                "No object provided or found as 'self' to check."
            )

    # 1) Custom __has_read__ method
    if hasattr(obj, "__has_read__") and callable(obj.__has_read__):
        if not obj.__has_read__():
            custom_msg = msg or (
                f"'{obj.__class__.__name__}' reports not read via "
                "__has_read__()."
            )
            raise NotReadError(custom_msg)
        return True

    # 2) Boolean flag
    if hasattr(obj, "_has_read"):
        flag = obj._has_read
        if flag is False:
            custom_msg = msg or (
                f"'{obj.__class__.__name__}' has `_has_read` "
                "flag set to False."
            )
            raise NotReadError(custom_msg)
        if flag is True:
            return True

    # 3) Fallback check: verify specified attributes
    if attributes is None:
        attributes = []
    if isinstance(attributes, str):
        attributes = [attributes]

    if not attributes:
        custom_msg = msg or (
            f"'{obj.__class__.__name__}' has no explicit read flag "
            "and no attributes were specified for checking."
        )
        raise NotReadError(custom_msg)

    for attr in attributes:
        value = getattr(obj, attr, None)
        is_empty = False
        if value is None:
            is_empty = True
        elif isinstance(value, pd.DataFrame) and value.empty:
            is_empty = True
        elif hasattr(value, "__len__") and len(value) == 0:
            is_empty = True

        if is_empty:
            custom_msg = msg or (
                f"'{obj.__class__.__name__}' has not been read. "
                f"Attribute '{attr}' is missing or empty. Call a "
                "read or from_* method first."
            )
            raise NotReadError(custom_msg)

    return True


def check_has_read(
    attributes: str | list[str] | None = None,
    *,
    msg: str | None = None,
) -> Callable:
    r"""Decorator to verify an object has been read.

    This decorator uses the :func:`has_read` function to protect
    a method, ensuring the object instance (`self`) has been
    properly populated with data before the method is executed.

    Parameters
    ----------
    attributes : str or list of str, optional
        The attribute name(s) to pass to the `has_read` function
        for validation. This is typically the name of a DataFrame
        or a core data attribute.
    msg : str, optional
        A custom error message to raise if validation fails.

    Returns
    -------
    callable
        The decorated function.

    Examples
    --------
    >>> from pycsamt.exceptions import NotReadError
    >>> class MyDataLoader:
    ...     def __init__(self):
    ...         self.df = None
    ...     def load(self):
    ...         # In a real scenario, this would load data.
    ...         self.df = "some data"
    ...     @check_has_read(attributes="df")
    ...     def get_data(self):
    ...         return self.df
    ...
    >>> loader = MyDataLoader()
    >>> try:
    ...     loader.get_data()
    ... except NotReadError as e:
    ...     print(e)
    'MyDataLoader' has not been read yet. Attribute 'df' is ...
    >>> loader.load()
    >>> loader.get_data()
    'some data'

    See Also
    --------
    has_read : The underlying validation function.
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            # 'self' is the instance the decorated method is called on
            has_read(self, attributes=attributes, msg=msg)
            return func(self, *args, **kwargs)

        return wrapper

    return decorator


def check_is_fitted(
    obj: Any = None, attributes: list[str] | None = None
) -> bool:
    r"""
    Validate that an object is "fitted" before use.

    - If `obj` has a __is_fitted__ method, it is called.
    - Else if `attributes` is provided, each named attribute on
    `obj` must be non-None.
    - Else, if `obj` has an `_fitted` or `_is_fitted` attribute,
    it must be True.

    If no `obj` is passed, attempts to use `self` from the caller's frame.

    Raises
    ------
    NotFittedError
        If the fitted checks fail.
    ValueError
        If no object can be determined.

    Returns
    -------
    bool
        True if fitted, otherwise raises.
    """
    # Retrieve obj from caller if not explicitly provided
    if obj is None:
        frame = inspect.currentframe().f_back
        obj = frame.f_locals.get("self")  # type: ignore
    if obj is None:
        raise ValueError(
            "No object provided or found as 'self' to check fitting."
        )

    # 1) Custom __is_fitted__ method
    if hasattr(obj, "__is_fitted__") and callable(obj.__is_fitted__):
        if not obj.__is_fitted__():
            raise NotFittedError(
                f"{obj.__class__.__name__} reports"
                f" not fitted via __is_fitted__()."
            )
        return True

    # 2) Explicit attribute presence
    if attributes:
        missing = [
            attr for attr in attributes if getattr(obj, attr, None) is None
        ]
        if missing:
            raise NotFittedError(
                f"{obj.__class__.__name__} missing required"
                f" attribute(s) {missing} (not fitted)."
            )
        return True

    # 3) Boolean flag
    flag = None
    if hasattr(obj, "_is_fitted"):
        flag = obj._is_fitted  # type: ignore
    elif hasattr(obj, "_fitted"):
        flag = obj._fitted  # type: ignore

    if isinstance(flag, bool):
        if not flag:
            raise NotFittedError(
                f"{obj.__class__.__name__}._fitted flag is False."
            )
        return True

    # No fitting indicator found
    raise NotFittedError(
        f"Could not determine fitting state for {obj.__class__.__name__}."
    )


def _assert_all_types(
    obj: Any, *expected_types: type, objname: str = None
) -> Any:
    """
    Assert that an object is an instance of the given types.

    Parameters
    ----------
    obj : any
        The object to test.
    *expected_types : type
        One or more expected types for `obj`.
    objname : str, optional
        Name of the object for error messages.

    Returns
    -------
    obj : any
        The original object if type check passes.

    Raises
    ------
    TypeError
        If `obj` is not an instance of any of `expected_types`.
    """
    if not expected_types:
        raise TypeError("No expected types provided for type assertion.")
    # flatten nested tuples in expected_types
    types_tuple: tuple = ()
    for t in expected_types:
        if isinstance(t, tuple):
            types_tuple += t
        else:
            types_tuple += (t,)
    if not isinstance(obj, types_tuple):
        type_names = ", ".join(t.__name__ for t in types_tuple)
        prefix = f"'{objname}' " if objname else ""
        plural = "s" if len(types_tuple) > 1 else ""
        raise TypeError(
            f"{prefix}expected type{plural} {type_names}, "
            f"got {type(obj).__name__}"
        )
    return obj


def _isin(
    arr: np.ndarray | Any,
    subarr: np.ndarray | Any,
    *,
    return_mask: bool = False,
) -> bool | np.ndarray:
    r"""
    Test membership of ``subarr`` in ``arr``, or return mask.

    Parameters
    ----------
    arr : array-like
        Universe of values.
    subarr : array-like or scalar
        Values to test for presence in ``arr``.
    return_mask : bool, default False
        If ``True``, return a boolean mask over ``arr``.
        Else return ``True`` if *all* of ``subarr`` are in
        ``arr``.

    Returns
    -------
    mask : ndarray of bool
        When ``return_mask=True``.
    result : bool
        When ``return_mask=False``; ``True`` if all elements
        of ``subarr`` are present in ``arr``.

    Raises
    ------
    ValueError
        If inputs cannot be coerced to arrays.
    """
    try:
        a = np.asarray(arr, dtype=object)
        s = np.asarray(subarr, dtype=object)
    except Exception as exc:
        raise ValueError(f"Invalid inputs for membership: {exc}") from exc

    mask = np.isin(a, s)
    if return_mask:
        return mask

    s1 = np.atleast_1d(s).ravel()
    # unique values to avoid redundant checks
    try:
        s1 = np.unique(s1)
    except Exception:
        pass
    return bool(np.all(np.isin(s1, a)))


def isin(
    arr: Any,
    subarr: Any,
    *,
    return_mask: bool = False,
    match: str = "all",
    assume_unique: bool = False,
    equal_nan: bool = False,
    invert: bool = False,
    return_missing: bool = False,
    return_present: bool = False,
):
    r"""
    Improved membership test with flexible outputs.

    Parameters
    ----------
    arr : array-like
        Universe of values to test against.
    subarr : array-like or scalar
        Values to check for presence in ``arr``.
    return_mask : bool, default False
        If ``True``, return a boolean mask with the same
        shape as ``arr`` (membership of ``arr`` in ``subarr``).
    match : {"all","any","count"}, default "all"
        Reduction mode when ``return_mask`` is ``False``:

        - ``"all"``  → True if all of ``subarr`` in ``arr``.
        - ``"any"``  → True if any of ``subarr`` in ``arr``.
        - ``"count"``→ Number of elements of ``subarr`` in
          ``arr`` (unique-sensitive per NumPy rules).

    assume_unique : bool, default False
        Forwarded to :func:`numpy.isin` for speed when both
        inputs have unique values.
    equal_nan : bool, default False
        Treat NaNs in ``arr`` and ``subarr`` as equal.
    invert : bool, default False
        When ``return_mask=True``, invert membership on the
        elementwise mask (same as NumPy).
    return_missing : bool, default False
        When ``return_mask=False``, also return a list of
        items from ``subarr`` that are not in ``arr``.
    return_present : bool, default False
        When ``return_mask=False``, also return a list of
        items from ``subarr`` that are present in ``arr``.

    Returns
    -------
    mask : ndarray of bool
        If ``return_mask=True``.
    result : bool or int
        If ``return_mask=False``. ``bool`` for ``"all"`` and
        ``"any"``; ``int`` for ``"count"``.
    (result, missing) : tuple
        When ``return_missing=True``.
    (result, present) : tuple
        When ``return_present=True``.
    (result, missing, present) : tuple
        When both flags are ``True``.

    Notes
    -----
    - Uses ``dtype=object`` to robustly compare mixed types.
    - For reduction, membership is computed on flattened
      ``subarr`` values.

    Examples
    --------
    >>> isin([1, 2, 3], [2, 5], match="any")
    True
    >>> isin([1, 2, 3], [2, 5], match="all")
    False
    >>> isin([1, 2, 3], [2, 2, 3], match="count")
    2
    >>> isin([1, 2, 3], [2, 5], return_missing=True)
    (False, [5])
    """
    try:
        a = np.asarray(arr, dtype=object)
        s = np.asarray(subarr, dtype=object)
    except Exception as exc:
        raise ValueError(f"Invalid inputs for membership: {exc}") from exc

    # np.isin has no `equal_nan` support: NaN membership is
    # patched in manually below (pd.isna also handles object dtype)
    def _isin_with_nan(elements, test_elements, inv):
        m = np.isin(
            elements,
            test_elements,
            assume_unique=assume_unique,
            invert=inv,
        )
        if equal_nan and bool(np.any(pd.isna(test_elements))):
            nan_m = pd.isna(elements)
            m = (m & ~nan_m) if inv else (m | nan_m)
        return m

    # elementwise mask over arr against subarr
    mask = _isin_with_nan(a, s, invert)
    if return_mask:
        return mask

    # reduction over subarr membership in arr
    s1 = np.atleast_1d(s).ravel()
    in_s = _isin_with_nan(s1, a, False)

    if match not in {"all", "any", "count"}:
        raise ValueError("match must be 'all', 'any', or 'count'.")

    if match == "all":
        result = bool(np.all(in_s))
    elif match == "any":
        result = bool(np.any(in_s))
    else:
        # count *distinct* values of subarr present in arr, per the
        # documented example: isin([1, 2, 3], [2, 2, 3], "count") -> 2
        try:
            s_uniq = np.unique(s1)
        except Exception:
            s_uniq = np.asarray(
                list(dict.fromkeys(s1.tolist())), dtype=object
            )
        result = int(np.count_nonzero(_isin_with_nan(s_uniq, a, False)))

    # optional extras
    out = (result,)
    if return_missing:
        missing = s1[~in_s].tolist()
        out = out + (missing,)
    if return_present:
        present = s1[in_s].tolist()
        out = out + (present,)

    return out[0] if len(out) == 1 else out


def isin_if(
    o: Iterable,
    items: str | Iterable,
    *,
    error: str = "raise",
    return_diff: bool = False,
    return_intersect: bool = False,
) -> list | None:
    r"""
    Check presence of ``items`` inside iterable ``o``.

    When requested, return the missing items (difference) or
    the items found (intersection). Optionally raise or warn
    if some items are missing.

    Parameters
    ----------
    o : iterable
        Container to search in. If a string, it is treated
        as a single token (not a char sequence).
    items : str or iterable
        Item or collection to check. Strings are treated as a
        single token.
    error : {"raise","warn","ignore"}, default "raise"
        Behavior when missing items are detected (ignored if
        a return is requested).
    return_diff : bool, default False
        If ``True``, return the list of missing items.
    return_intersect : bool, default False
        If ``True``, return the list of intersecting items.

    Returns
    -------
    list or None
        Missing items (``return_diff``), or intersection
        (``return_intersect``). Otherwise ``None``.

    Raises
    ------
    TypeError
        If ``o`` is not iterable.
    ValueError
        If ``error`` is invalid or missing items with
        ``error='raise'``.
    """
    if not isinstance(o, Iterable):
        raise TypeError("`o` must be an iterable.")

    def _to_set(v, is_items: bool) -> set:
        # strings → single token
        if isinstance(v, (str, bytes)):
            return {v}
        try:
            return set(v)  # type: ignore[arg-type]
        except TypeError:
            # scalars (e.g., numbers)
            return {v}

    set_o = _to_set(o, is_items=False)
    set_items = _to_set(items, is_items=True)

    inter = list(set_o & set_items)
    missing = list(set_items - set_o)

    # request to return overrides error behavior
    if return_diff or return_intersect:
        error = "ignore"

    if missing:
        if error not in {"raise", "warn", "ignore"}:
            raise ValueError("error must be {'raise','warn','ignore'}.")
        msg = "Missing item(s): " + ", ".join(repr(x) for x in missing)
        if error == "raise":
            raise ValueError(msg)
        if error == "warn":
            warnings.warn(msg, stacklevel=2)

    if return_diff:
        return missing
    if return_intersect:
        return inter
    return None


def assert_ratio(
    v: Any,
    bounds: Sequence[float] | None = None,
    exclude_value: float | None = None,
    in_percent: bool = False,
    name: str = "rate",
) -> float:
    """
    Assert that a ratio value falls within given bounds and
    optionally exclude a specific value.

    Parameters
    ----------
    v : any
        Value to assert. If string with '%', treated as percent.
    bounds : sequence of two floats, optional
        (lower, upper) inclusive bounds for `v`.
    exclude_value : float, optional
        A value that `v` must not equal. If conversion fails,
        lower bound is used for exclusion if provided.
    in_percent : bool, default False
        Interpret and return `v` as fraction (divide by 100).
    name : str, default 'rate'
        Name for error messages.

    Returns
    -------
    float
        The validated (and possibly converted) value.

    Raises
    ------
    TypeError
        If `v` cannot be converted to float.
    ValueError
        If `v` is outside `bounds` or equals `exclude_value`.
    """
    # parse string with percent
    if isinstance(v, str):
        if "%" in v:
            in_percent = True
        v = v.replace("%", "")
    # convert to float
    try:
        val = float(v)
    except Exception:
        raise TypeError(
            f"Unable to convert {type(v).__name__!r} to float: {v!r}"
        )
    # percent conversion
    if in_percent:
        if 1 < val <= 100:
            val /= 100.0
    # unpack bounds
    low = up = None
    if bounds is not None:
        if len(bounds) != 2:
            raise ValueError(
                f"`bounds` must have two elements, got {bounds!r}"
            )
        low, up = bounds
        if low is not None and up is not None:
            if not (low <= val <= up):
                raise ValueError(
                    f"{name} must be between {low} and {up}, got {val}"
                )
    # handle exclusion
    if exclude_value is not None:
        try:
            excl = float(exclude_value)
        except Exception:
            excl = low
            if excl is None:
                warnings.warn(
                    "Cannot exclude value without valid bounds", stacklevel=2
                )
        if val == excl:
            raise ValueError(f"{name} excluding {excl}, got {val}")
    # post-check for percent
    if in_percent and val > 1.0:
        msg = f"{name} as percent must be <= 1.0, got {val}"
        raise ValueError(msg)
    return val


def _validate_name_in(
    name: str,
    defaults: Sequence[str] | str = "",
    expect_name: str | None = None,
    exception: Exception | None = None,
    deep: bool = False,
) -> bool | str:
    """
    Assert that `name` exists within `defaults`.

    Parameters
    ----------
    name : str
        Name to validate.
    defaults : sequence of str or str, default ''
        Allowed names for validation.
    expect_name : str, optional
        If provided, returned when validation passes.
    exception : Exception, optional
        Exception to raise on failure.
    deep : bool, default False
        If True, concatenates defaults and checks substring.

    Returns
    -------
    bool or str
        True if valid and no expect_name; otherwise expect_name or False.

    Raises
    ------
    Exception
        If not valid and `exception` is provided.
    """
    sname = name.lower().strip()
    # normalize defaults
    defs_seq = []
    if isinstance(defaults, str):
        defs_seq = [defaults.lower().strip()]
    else:
        try:
            defs_seq = [str(d).lower().strip() for d in defaults]
        except Exception:
            raise TypeError("`defaults` must be str or sequence of str")
    # check membership
    if deep:
        merged = "".join(defs_seq)
        valid = sname in merged
    else:
        valid = sname in defs_seq
    # determine return
    result = expect_name if (valid and expect_name) else valid
    if not valid and exception is not None:
        raise exception
    return result


def isinstance_relaxed(
    instance: Any,
    cls: type | tuple[type, ...],
) -> bool:
    r"""
    Robust ``isinstance`` variant tolerant to reloads/tuples.

    First performs a normal :func:`isinstance` check. If it
    fails, it compares class *names* and *module tails*
    (final path segment) to mitigate false negatives that can
    occur after module reloads or aliasing.

    Parameters
    ----------
    instance : object
        Object to test.
    cls : type or tuple of types
        Target class or a tuple of classes.

    Returns
    -------
    bool
        ``True`` if the object is an instance of any target
        class, either by direct check or by relaxed matching.

    Notes
    -----
    The relaxed path matches when both class names are equal
    and the last component of the module path matches, e.g.,
    ``"pkg.mod"`` and ``"other.mod"`` → ``"mod"``.

    Examples
    --------
    >>> class A: ...
    >>> a = A()
    >>> isinstance_relaxed(a, A)
    True
    >>> isinstance_relaxed(a, (int, A))
    True
    """
    # direct isinstance for fast path
    if isinstance(instance, cls):
        return True

    # normalize target to a tuple
    targets: tuple[type, ...]
    if isinstance(cls, tuple):
        targets = cls
    else:
        targets = (cls,)

    # safe accessors
    inst_type = getattr(instance, "__class__", type(instance))
    iname = getattr(inst_type, "__name__", None)
    imod = getattr(inst_type, "__module__", "") or ""

    imod_tail = imod.rsplit(".", 1)[-1]

    for t in targets:
        tname = getattr(t, "__name__", None)
        tmod = getattr(t, "__module__", "") or ""
        tmod_tail = tmod.rsplit(".", 1)[-1]

        if iname is None or tname is None:
            continue

        # relaxed: name match + module-tail match
        if iname == tname and imod_tail == tmod_tail:
            return True

    return False


def _is_numeric_dtype(o, /, to_array: bool = False) -> bool:
    r"""
    Check whether ``o`` has a numeric dtype as a NumPy array.

    Numeric kinds are booleans, unsigned/signed ints, floats,
    and complex numbers.

    Parameters
    ----------
    o : array-like
        Iterable to check.
    to_array : bool, default False
        If ``True``, coerce ``o`` to ``np.array`` first.

    Returns
    -------
    bool
        ``True`` if dtype is numeric, else ``False``.
    """
    _NUMERIC_KINDS = set("buifc")

    if not hasattr(o, "__iter__"):
        raise TypeError(f"'o' must be iterable. Got: {type(o).__name__!r}")

    if to_array:
        o = np.array(o)

    if not hasattr(o, "__array__"):
        raise ValueError(f"Expect array-like. Got: {type(o).__name__!r}")

    # prefer dtype.kind on ndarray/Series/DataFrame
    kind = (
        o.values.dtype.kind
        if (hasattr(o, "columns") or hasattr(o, "name"))
        else o.dtype.kind
    )
    return kind in _NUMERIC_KINDS


def _check_consistency_size(ar1, ar2, /, error: str = "raise") -> bool:
    r"""
    Check same length for two arrays. Optionally raise.

    Parameters
    ----------
    ar1, ar2 : array-like
        Arrays to compare.
    error : {"raise","ignore"}, default "raise"
        Raise on mismatch or return ``False``.

    Returns
    -------
    bool
        ``True`` if sizes match, else ``False``.
    """
    same = len(ar1) == len(ar2)
    if not same and error == "raise":
        msg = f"Array sizes must match: '{len(ar1)}' vs '{len(ar2)}'."
        raise AssertionError(msg)
    return same


def check_consistency_size(*arrays) -> None:
    r"""
    Ensure all arrays have the same number of samples.

    Parameters
    ----------
    *arrays : array-like
        Inputs to validate. ``None`` entries are ignored.

    Raises
    ------
    ValueError
        If lengths differ among provided arrays.
    """
    lengths = [len(X) for X in arrays if X is not None]
    uniques = np.unique(lengths)
    if len(uniques) > 1:
        raise ValueError(
            f"Inconsistent sample sizes: {[int(l) for l in lengths]}"
        )


def _is_arraylike_1d(x) -> bool:
    r"""
    Return ``True`` if input is 1-D array-like (not scalar).

    Raises
    ------
    TypeError
        If input is not array-like.
    """
    if not hasattr(x, "__array__"):
        raise TypeError(f"Expect a 1-D array. Got: {type(x).__name__!r}")
    if not _is_arraylike_not_scalar(x):
        return False
    nd = getattr(x, "ndim", None)
    shape = getattr(x, "shape", None)
    if nd is None or shape is None:
        # fallback: treat as 1-D if len works
        try:
            _ = len(x)
            return True
        except Exception:
            return False
    return (nd < 2) or (nd == 2 and shape[1] == 1)


def _is_arraylike(x) -> bool:
    r"""
    Return whether the input is array-like.
    """
    return (
        hasattr(x, "__len__")
        or hasattr(x, "shape")
        or hasattr(x, "__array__")
    )


def _is_arraylike_not_scalar(array) -> bool:
    r"""
    Return ``True`` if array-like and not a scalar.
    """
    return _is_arraylike(array) and not np.isscalar(array)
