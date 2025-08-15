# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Validation utilities for pycsamt.
"""

import inspect
import warnings
from functools import wraps 
from typing import ( 
    List, 
    Optional, 
    Any, 
    Union, 
    Sequence, 
    Callable
)

import numpy as np 
import pandas as pd 

from ..exceptions import ( 
    NotFittedError, 
    NotReadError 
)

__all__ = [
    'check_is_fitted', 
    '_assert_all_types',
    '_isin',
    'assert_ratio',
    '_validate_name_in', 
    '_is_numeric_dtype', 
    'check_consistency_size', 
    '_is_arraylike_1d', 
    'is_instance_extended', 
    "has_read", 
    "check_has_read"
]

def has_read(
    obj: Any = None,
    *,
    attributes: Union[str, List[str], None] = None,
    msg: Optional[str] = None,
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
            obj = frame.f_back.f_locals.get('self')
        if obj is None:
            raise ValueError(
                "No object provided or found as 'self' to check."
            )

    # 1) Custom __has_read__ method
    if hasattr(obj, '__has_read__') and callable(
        getattr(obj, '__has_read__')
    ):
        if not obj.__has_read__():
            custom_msg = msg or (
                f"'{obj.__class__.__name__}' reports not read via "
                "__has_read__()."
            )
            raise NotReadError(custom_msg)
        return True

    # 2) Boolean flag
    if hasattr(obj, '_has_read'):
        flag = getattr(obj, '_has_read')
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
    attributes: Union[str, List[str], None] = None,
    *,
    msg: Optional[str] = None,
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
        obj: Any = None, attributes: Optional[List[str]] = None
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
        obj = frame.f_locals.get('self')  # type: ignore
    if obj is None:
        raise ValueError(
            "No object provided or found as 'self' to check fitting.")

    # 1) Custom __is_fitted__ method
    if hasattr(obj, '__is_fitted__') and callable(obj.__is_fitted__):
        if not obj.__is_fitted__():
            raise NotFittedError(
                f"{obj.__class__.__name__} reports"
                f" not fitted via __is_fitted__().")
        return True

    # 2) Explicit attribute presence
    if attributes:
        missing = [attr for attr in attributes 
                   if getattr(obj, attr, None) is None]
        if missing:
            raise NotFittedError(
                f"{obj.__class__.__name__} missing required"
                f" attribute(s) {missing} (not fitted)."
            )
        return True

    # 3) Boolean flag
    flag = None
    if hasattr(obj, '_is_fitted'):
        flag = getattr(obj, '_is_fitted')  # type: ignore
    elif hasattr(obj, '_fitted'):
        flag = getattr(obj, '_fitted')  # type: ignore

    if isinstance(flag, bool):
        if not flag:
            raise NotFittedError(
                f"{obj.__class__.__name__}._fitted flag is False.")
        return True

    # No fitting indicator found
    raise NotFittedError(
        f"Could not determine fitting state for {obj.__class__.__name__}."
    )

def _assert_all_types(
    obj: Any,
    *expected_types: type,
    objname: str = None
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
        raise TypeError(
            "No expected types provided for type assertion."
        )
    # flatten nested tuples in expected_types
    types_tuple: tuple = ()
    for t in expected_types:
        if isinstance(t, tuple):
            types_tuple += t
        else:
            types_tuple += (t,)
    if not isinstance(obj, types_tuple):
        type_names = ", ".join(t.__name__ for t in types_tuple)
        prefix = f"'{objname}' " if objname else ''
        plural = 's' if len(types_tuple) > 1 else ''
        raise TypeError(
            f"{prefix}expected type{plural} {type_names}, "
            f"got {type(obj).__name__}"
        )
    return obj


def _isin(
    arr: Union[np.ndarray, Any],
    subarr: Union[np.ndarray, Any],
    return_mask: bool = False
) -> Union[bool, np.ndarray]:
    """
    Test whether all elements of `subarr` are in `arr`, or return mask.

    Parameters
    ----------
    arr : array-like
        Array of elements to test against.
    subarr : array-like or scalar
        Element(s) to test for membership in `arr`.
    return_mask : bool, default False
        If True, return boolean mask of same shape as `arr`.
        If False, return whether all elements of `subarr` are in `arr`.

    Returns
    -------
    mask : ndarray of bool
        Boolean mask if `return_mask` is True.
    result : bool
        True if all elements of `subarr` are in `arr`, False otherwise.

    Raises
    ------
    ValueError
        If inputs cannot be converted to NumPy arrays.
    """
    try:
        a = np.asarray(arr)
        s = np.asarray(subarr)
    except Exception as e:
        raise ValueError(
            f"Invalid input for membership test: {e}"
        )
    mask = np.isin(a, s)
    if return_mask:
        return mask
    # flatten unique subarr values
    try:
        s_vals = np.unique(s.ravel())
    except Exception:
        s_vals = np.array([s])
    return bool(np.all(np.isin(s_vals, a)))


def assert_ratio(
    v: Any,
    bounds: Optional[Sequence[float]] = None,
    exclude_value: Optional[float] = None,
    in_percent: bool = False,
    name: str = 'rate'
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
        if '%' in v:
            in_percent = True
        v = v.replace('%', '')
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
                    "Cannot exclude value without valid bounds"
                )
        if val == excl:
            raise ValueError(
                f"{name} excluding {excl}, got {val}"
            )
    # post-check for percent
    if in_percent and val > 1.0:
        msg = f"{name} as percent must be <= 1.0, got {val}"
        raise ValueError(msg)
    return val


def _validate_name_in(
    name: str,
    defaults: Union[Sequence[str], str] = '',
    expect_name: Optional[str] = None,
    exception: Optional[Exception] = None,
    deep: bool = False
) -> Union[bool, str]:
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
            raise TypeError(
                "`defaults` must be str or sequence of str"
            )
    # check membership
    if deep:
        merged = ''.join(defs_seq)
        valid = sname in merged
    else:
        valid = sname in defs_seq
    # determine return
    result = expect_name if (valid and expect_name) else valid
    if not valid and exception is not None:
        raise exception
    return result

def _is_numeric_dtype (o, / , to_array =False ): 
    """ Determine whether the argument has a numeric datatype, when
    converted to a NumPy array.

    Booleans, unsigned integers, signed integers, floats and complex
    numbers are the kinds of numeric datatype. 
    
    :param o: object, arraylike 
        Object presumed to be an array 
    :param to_array: bool, default=False 
        If `o` is passed as non-array like list or tuple or other iterable 
        object. Setting `to_array` to ``True`` will convert `o` to array. 
    :return: bool, 
        ``True`` if `o` has a numeric dtype and ``False`` otherwise. 
    """ 
    _NUMERIC_KINDS = set('buifc')
    if not hasattr (o, '__iter__'): 
        raise TypeError ("'o' is expected to be an iterable object."
                         f" got: {type(o).__name__!r}")
    if to_array : 
        o = np.array (o )
    if not hasattr(o, '__array__'): 
        raise ValueError (f"Expect type array, got: {type (o).__name__!r}")
    # use NUMERICKIND rather than # pd.api.types.is_numeric_dtype(arr) 
    # for series and dataframes
    return ( o.values.dtype.kind   
            if ( hasattr(o, 'columns') or hasattr (o, 'name'))
            else o.dtype.kind ) in _NUMERIC_KINDS 
        
def _check_consistency_size (ar1, ar2 , /  , error ='raise') :
    """ Check consistency of two arrays and raises error if both sizes 
    are differents. 
    Returns 'False' if sizes are not consistent and error is set to 'ignore'.
    """
    if error =='raise': 
        msg =("Array sizes must be consistent: '{}' and '{}' were given.")
        assert len(ar1)==len(ar2), msg.format(len(ar1), len(ar2))
        
    return len(ar1)==len(ar2) 

def check_consistency_size ( *arrays ): 
    """ Check consistency of array and raises error otherwise."""
    lengths = [len(X) for X in arrays if X is not None]
    uniques = np.unique(lengths)
    if len(uniques) > 1:
        raise ValueError(
            "Found input variables with inconsistent numbers of samples: %r"
            % [int(l) for l in lengths]
        )
        
def _is_arraylike_1d (x) :
    """ Returns whether the input is arraylike one dimensional and not a scalar"""
    if not hasattr (x, '__array__'): 
        raise TypeError ("Expects a one-dimensional array, "
                         f"got: {type(x).__name__!r}")
    _is_arraylike_not_scalar(x)
    return _is_arraylike_not_scalar(x) and  (  len(x.shape )< 2 or ( 
        len(x.shape ) ==2 and x.shape [1]==1 )) 

def _is_arraylike(x):
    """Returns whether the input is array-like."""
    return hasattr(x, "__len__") or hasattr(x, "shape") or hasattr(x, "__array__")


def _is_arraylike_not_scalar(array):
    """Return True if array is array-like and not a scalar"""
    return _is_arraylike(array) and not np.isscalar(array)

def is_instance_extended(instance, cls):
    """
    Performs an enhanced isinstance check that can gracefully handle a tuple 
    of classes and module reloading issues, facilitating a more robust type 
    checking, especially in environments where classes might be reloaded or 
    imported differently, potentially leading to false negatives with the 
    standard isinstance function.

    Parameters
    ----------
    instance : object
        The object to check.
    cls : type or tuple of types
        The target class, classes, or a tuple of classes to check against. 
        If `cls` is not a tuple, it will be converted to one for uniform handling.

    Returns
    -------
    bool
        True if `instance` is an instance of any class in `cls`, considering class 
        name and module path matches. False otherwise.

    Examples
    --------
    >>> class MyClass:
    ...     pass
    ...
    >>> obj = MyClass()
    >>> is_instance_extended(obj, MyClass)
    True

    # Demonstrating with module reloading issue
    >>> import importlib
    >>> importlib.reload(MyClass)
    <module 'MyClass' from '...'>
    >>> is_instance_extended(obj, MyClass)
    False  # This might vary based on how MyClass is defined and reloaded

    # Using a tuple of classes
    >>> class AnotherClass:
    ...     pass
    ...
    >>> is_instance_extended(obj, (MyClass, AnotherClass))
    True

    Note
    ----
    This function is particularly useful in dynamic environments where classes may 
    be reloaded or when dealing with complex import hierarchies that could lead to 
    situations where the standard `isinstance` check might erroneously return False 
    due to objects being instances of classes that have been reloaded or imported 
    under different namespaces.
    """
    if not isinstance(cls, tuple):
        cls = (cls,)  # Make cls a tuple if it isn't already, for uniform handling

    direct_check = any(isinstance(instance, single_cls) for single_cls in cls)
    if direct_check:
        return True

    for single_cls in cls:
        if instance.__class__.__name__ == single_cls.__name__:
            instance_module = instance.__class__.__module__.split('.')[-1]
            cls_module = single_cls.__module__.split('.')[-1]
            if instance_module == cls_module:
                return True
    return False

