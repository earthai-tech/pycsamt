"""Shared typing helpers for PyCSAMT.

This module centralizes type aliases used across the v2 codebase.
It intentionally keeps the legacy names that appeared in v1
annotations, while mapping them to safer runtime objects.  The goal is
to let old annotations such as ``ArrayLike[DType[float]]`` continue to
evaluate, and to provide clearer aliases for new code.

Most aliases in this module are documentation and static-analysis
helpers.  They should not be used for runtime validation.  Use
validators from :mod:`pycsamt.utils.validation` or local parsing code
when input values must be checked.

Examples
--------
>>> from pycsamt.api.typing import ArrayLike, NDArray, PathLike
>>> def normalize(values: ArrayLike[float]) -> NDArray[float]:
...     ...

>>> def read_any(path: PathLike) -> str:
...     ...
"""

from __future__ import annotations

import sys
from collections.abc import (
    Callable,
    Iterable,
    Iterator,
    Mapping,
    MutableMapping,
    Sequence,
)
from os import PathLike as OSPathLike
from pathlib import Path
from typing import (
    Any,
    ClassVar,
    Generic,
    Literal,
    Optional,
    Protocol,
    SupportsFloat,
    SupportsInt,
    Text,
    TypeVar,
    Union,
    overload,
)

# ``TypeAlias`` and ``TypeGuard`` were only added to ``typing`` in Python 3.10;
# on 3.9 they live in ``typing_extensions``. Import compatibly so the package
# (requires-python >= 3.9) imports on every supported version.
if sys.version_info >= (3, 10):
    from typing import TypeAlias, TypeGuard
else:  # pragma: no cover - exercised on Python 3.9 CI
    from typing_extensions import TypeAlias, TypeGuard

import numpy as np
from numpy.typing import ArrayLike as NumpyArrayLike
from numpy.typing import DTypeLike
from numpy.typing import NDArray as NumpyNDArray

try:  # pragma: no cover - pandas is expected but kept optional.
    import pandas as _pd
except Exception:  # pragma: no cover
    _pd = None

__all__ = [
    "Any",
    "Array",
    "Array1D",
    "Array2D",
    "ArrayLike",
    "Callable",
    "ClassVar",
    "DataFrame",
    "DataFrameLike",
    "Dict",
    "DType",
    "DTypeLike",
    "EDIO",
    "F",
    "FloatArray",
    "Generic",
    "IndexLike",
    "IntArray",
    "Iterable",
    "Iterator",
    "K",
    "List",
    "Literal",
    "M",
    "Mapping",
    "MutableMapping",
    "N",
    "NDArray",
    "NumpyArrayLike",
    "NumpyNDArray",
    "Numeric",
    "Optional",
    "Path",
    "PathLike",
    "Protocol",
    "SP",
    "Scalar",
    "Sequence",
    "Series",
    "SeriesLike",
    "Shape",
    "Sub",
    "SupportsArray",
    "SupportsFloat",
    "SupportsInt",
    "T",
    "Text",
    "Type",
    "TypeAlias",
    "TypeGuard",
    "Tuple",
    "U",
    "Union",
    "V",
    "ZO",
    "overload",
]


# ---------------------------------------------------------------------------
# Standard-library compatibility exports
# ---------------------------------------------------------------------------

# The v1 module exported names from :mod:`typing`.  Keep those names so
# older modules can continue to import from ``pycsamt.api.typing``.
List = list
Tuple = tuple
Dict = dict
Type = type

T = TypeVar("T")
V = TypeVar("V")
K = TypeVar("K")
M = TypeVar("M", bound=int)
N = TypeVar("N", bound=int)
U = TypeVar("U")
F = TypeVar("F", bound=Callable[..., Any])


# ---------------------------------------------------------------------------
# Modern v2 aliases
# ---------------------------------------------------------------------------

PathLike: TypeAlias = Union[str, bytes, Path, OSPathLike[str]]
"""Path accepted by PyCSAMT readers and writers."""

Scalar: TypeAlias = Union[str, bytes, int, float, complex, bool]
"""Common scalar value accepted by lightweight utilities."""

Numeric: TypeAlias = Union[int, float, complex, np.number]
"""Python or NumPy numeric scalar."""

Array: TypeAlias = NumpyNDArray[Any]
"""Concrete NumPy array with arbitrary dtype and shape."""

Array1D: TypeAlias = NumpyNDArray[Any]
"""One-dimensional NumPy array by convention."""

Array2D: TypeAlias = NumpyNDArray[Any]
"""Two-dimensional NumPy array by convention."""

FloatArray: TypeAlias = NumpyNDArray[np.floating[Any]]
"""NumPy array with floating dtype."""

IntArray: TypeAlias = NumpyNDArray[np.integer[Any]]
"""NumPy array with integer dtype."""

IndexLike: TypeAlias = Union[int, slice, Sequence[int], NumpyNDArray[Any]]
"""Index selector accepted by array utilities."""

SeriesLike: TypeAlias = Any
"""Object behaving like a pandas Series."""

DataFrameLike: TypeAlias = Any
"""Object behaving like a pandas DataFrame."""


if _pd is not None:
    Series = _pd.Series
    DataFrame = _pd.DataFrame
else:  # pragma: no cover
    class Series:  # noqa: D101
        def __class_getitem__(cls, item: Any) -> type[Series]:
            return cls

    class DataFrame:  # noqa: D101
        def __class_getitem__(cls, item: Any) -> type[DataFrame]:
            return cls


# ---------------------------------------------------------------------------
# Legacy-compatible generic names
# ---------------------------------------------------------------------------

class _CompatAlias:
    """Runtime-safe base for legacy subscripted aliases.

    Older modules use aliases with several incompatible shapes, for
    example ``ArrayLike[DType[float]]`` and
    ``ArrayLike[float, DType[float]]``.  Returning the class from
    ``__class_getitem__`` preserves runtime evaluation without
    pretending that these aliases perform validation.
    """

    __args__: ClassVar[tuple[Any, ...]] = ()

    def __class_getitem__(cls, item: Any) -> type[_CompatAlias]:
        return cls


class Shape(_CompatAlias):
    """Shape marker kept for legacy annotations."""


class DType(_CompatAlias):
    """Dtype marker kept for legacy annotations.

    New code should prefer :class:`numpy.typing.DTypeLike` or a concrete
    NumPy dtype annotation.
    """


class ArrayLike(_CompatAlias):
    """Array-like input accepted by NumPy conversion routines.

    New code may use :class:`numpy.typing.ArrayLike` directly.  This
    compatibility alias remains subscriptable with one or two arguments
    for old PyCSAMT annotations.
    """


class NDArray(_CompatAlias):
    """NumPy array marker kept for legacy annotations.

    New code may use :class:`numpy.typing.NDArray` directly for stricter
    dtype annotations.
    """


class Sub(_CompatAlias):
    """Subset marker for legacy array annotations."""


class SP(_CompatAlias):
    """Station-position marker for legacy annotations."""


class EDIO(_CompatAlias):
    """Electrical Data Interchange object marker."""


class ZO(_CompatAlias):
    """Impedance tensor object marker."""


class SupportsArray(Protocol):
    """Protocol for objects convertible to NumPy arrays."""

    def __array__(self, dtype: DTypeLike | None = None) -> np.ndarray:
        """Return a NumPy representation of the object."""


def is_path_like(value: object) -> TypeGuard[PathLike]:
    """Return whether ``value`` can be treated as a path."""

    return isinstance(value, (str, bytes, Path, OSPathLike))
