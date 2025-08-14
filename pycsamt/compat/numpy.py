# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

"""
A compatibility layer so k-diagram works under both NumPy 1.x and 2.x.
Avoids removed aliases (np.bool, np.int, np.float, …).
"""

import numpy as np
from numpy.lib import NumpyVersion as _NV

_NVERSION = _NV(np.__version__)
IS_NP2 = _NVERSION >= "2.0.0"
IS_NP125_PLUS = _NVERSION >= "1.25.0"

# ---- Stable type aliases (prefer builtins) ----
# These work as dtype arguments and for isinstance checks in both 1.x and 2.x
int_ = int
float_ = float
bool_ = bool
complex_ = complex
str_ = str
object_ = object

# If you *need* NumPy scalar classes as well, expose them separately:
NP_INT = np.int32  # exists in 1.x and 2.x
NP_FLOAT = np.float64  # exists
NP_BOOL = np.bool_  # exists
# try:
#     NP_COMPLEX = np.complexfloating  # or np.complex128 if you need a concrete dtype
# except:
# NP_COMPLEX = np.complexfloating   # old
NP_COMPLEX = np.complex128  # new: concrete dtype, no deprecation

NP_STR = np.str_  # exists

# ---- Moved / renamed helpers ----
in1d = np.isin
row_stack = np.vstack
try:
    trapz = np.trapezoid  # new name
except AttributeError:
    trapz = np.trapz  # old name

# ---- AxisError location changed in 2.x ----
try:
    from numpy.exceptions import AxisError
except Exception:
    from numpy import AxisError

# Expose a stable “default int” dtype
default_int = np.intp


# === AxisError import compatibility ===
# In NumPy 2.x it lives under numpy.exceptions.
try:
    from numpy.exceptions import AxisError  # Noqa
except ImportError:
    from numpy import AxisError


# ---- Promotion warnings helper (NumPy 2.x only) ----
def set_promotion_warn(state: str = "weak_and_warn") -> None:
    """
    During testing you can enable warnings on changed type-promotion
    behavior (NumPy 2.0+).  E.g.:

        import warnings
        warnings.simplefilter('error')
        compat_numpy.set_promotion_warn()

    to turn those into errors and catch any unintended changes.
    """
    if hasattr(np, "_set_promotion_state"):
        np._set_promotion_state(state)


# ---- asarray signature wrapper (2.x adds copy=) ----
def asarray(x, dtype=None, copy=None):
    """
    During testing you can enable warnings on changed type-promotion
    behavior (NumPy 2.0+).  E.g.:

        import warnings
        warnings.simplefilter('error')
        compat_numpy.set_promotion_warn()

    to turn those into errors and catch any unintended changes.
    """
    # if "copy" in getattr(np.asarray, "__code__", ()).co_varnames:
    #     return np.asarray(x, dtype=dtype, copy=copy)
    # return np.asarray(x, dtype=dtype)

    try:
        # NumPy >= 2.0 supports copy=
        return np.asarray(x, dtype=dtype, copy=copy)
    except TypeError:
        # NumPy 1.x: no copy= kwarg
        return np.asarray(x, dtype=dtype)


# === Default integer dtype ===
# Expose what the “default” int type is in this NumPy build.
default_int = np.intp if IS_NP2 else np.int_


# --- find_common_type replacement ------------------------------------------


def _to_dtype_list(seq):
    return [np.dtype(x) for x in (seq or [])]


if hasattr(np, "find_common_type") and not IS_NP125_PLUS:
    # NumPy 1.x: keep behavior, no warnings from our shim
    # NumPy < 1.25: safe to call
    def find_common_type(array_dtypes=None, scalar_dtypes=None):
        return np.find_common_type(
            _to_dtype_list(array_dtypes), _to_dtype_list(scalar_dtypes)
        )

else:  # NumPy 2.x: use result_type
    # NumPy >= 1.25 (incl. 2.x): use result_type to avoid DeprecationWarning
    def find_common_type(array_dtypes=None, scalar_dtypes=None):
        arr = _to_dtype_list(array_dtypes)
        scl = _to_dtype_list(scalar_dtypes)
        if not arr and not scl:
            return np.dtype("float64")
        return np.result_type(*arr, *scl)


# === Public API ===
__all__ = [
    "IS_NP2",
    "IS_NP125_PLUS",
    "int_",
    "float_",
    "bool_",
    "object_",
    "complex_",
    "str_",
    "NP_INT",
    "NP_FLOAT",
    "NP_BOOL",
    "NP_COMPLEX",
    "NP_STR",
    "in1d",
    "row_stack",
    "trapz",
    "AxisError",
    "set_promotion_warn",
    "asarray",
    "default_int",
    "find_common_type",
]
