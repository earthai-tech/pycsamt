# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0 

"""
`pycsamt` Type variables
========================

.. |CSAMT| replace:: Controlled-Source Audio-Magnetotellurics

.. _pycsamt: https://github.com/your-repo/pycsamt/
.. _pandas DataFrame: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
.. _Series: https://pandas.pydata.org/docs/reference/api/pandas.Series.html

Some customized type variables need to be explained for easy
understanding in the whole package. Indeed, customized type
hints are used to define the type of arguments.

**M**: Represents an integer variable `IntVar` to denote the
    number of rows in an ``Array``.

**N**: Similar to ``M``, *N* represents the number of columns
    in an ``Array``. It is bound to an integer variable.

**T**: A generic type, standing for `Any` type of variable.

**U**: Stands for a single dimension in an array. For instance:
    >>> import numpy as np
    >>> array = np.arange(4).shape
    ... (4, )

**S**: Indicates the `Shape` status. It is bound by `M`, `U`,
    `N`. 'U' stands for a single dimension in a 1D array.
    The `AddShape` class allows for extending arrays beyond two
    dimensions.

**D**: Stands for dtype object. It is bound with :class:`DType`.

**Array**: Defines a one-dimensional array where `DType` can be
    specified. For instance:
    >>> import numpy as np
    >>> from pycsamt.api.typing import TypeVar, Array, DType
    >>> T = TypeVar ('T', float)
    >>> A = TypeVar ('A', str, bytes )
    >>> arr1:Array[T, DType[T]] = np.arange(21) # dtype ='float'
    >>> arr2: Array[A, DType[A]] = arr1.astype ('str') # dtype ='str'

**NDArray**: Stands for multi-dimensional arrays (more than one
    dimension).

**Sub**: Represents a subset of an ``Array``. For example,
    extracting a specific zone from a |CSAMT| data line.

**SP**: Stands for Station Positions. The unit of position may
    vary, but the default is in meters, starting from position 0.

**Series**: A generic type hint for a `pandas Series`_ object.

**DataFrame**: A generic type hint for a `pandas DataFrame`_
    object.

**EDIO**: Stands for Electrical Data Interchange (EDI) Object.
    It is an object built from `pycsamt`_ or `MTpy`_ packages.

---
Additional definition for common arguments
===========================================

**data_array**: A data array, typically holding apparent
    resistivity values. The type hint is usually
    ``Array[float, DType[float]]`` or ``List[float]``.

**p**: Represents station location positions. The type hint is
    ``Array[int, DType[int]]`` or ``List[int]``. Values are
    expected to be integers.

**cz**: Stands for Conductive Zone. It is a subset of a data
    array and shares the same type hint. The ``Sub`` type is
    used for better demarcation.
"""
from __future__ import annotations

from typing import (
    List,
    Tuple,
    Sequence,
    Dict,
    Iterable,
    Callable,
    Union,
    Any,
    Generic,
    Optional,
    Type,
    Mapping,
    Text,
    TypeVar,
    Iterator,
    SupportsInt,
)

__all__ = [
    "List", "Tuple", "Sequence", "Dict", "Iterable",
    "Callable", "Any", "Generic", "Optional", "Union",
    "Type", "Mapping", "Text", "Shape", "DType", "NDArray",
    "ArrayLike", "EDIO", "Sub", "SP", "F", "T", "V",
    "Series", "Iterator", "SupportsInt", "ZO"
]

T = TypeVar('T')
V = TypeVar('V')
K = TypeVar('K')
M = TypeVar('M', bound=int)
N = TypeVar('N', bound=int)
U = TypeVar('U')
D = TypeVar('D', bound='DType')
S = TypeVar('S', bound='Shape')


class AddShape(Generic[S]):
    """ An extra bound to top the `Shape` for dimensions > 2."""


class Shape(Generic[M, S], AddShape[S]):
    """ Generic to construct a tuple shape for NDarray."""
    def __getitem__(self, M, N) -> S:
        ...


class DType(Generic[T]):
    """ DType can be Any Type so it holds 'T' type variable. """
    def __getitem__(self, T) -> T:
        ...


class ArrayLike(Generic[T, D]):
    """ Array Type here means the 1D array. """
    def __getitem__(self, T) -> Union['ArrayLike', T]:
        ...


class NDArray(ArrayLike[T, DType[T]], Generic[T, D]):
    """ NDarray has M-rows, N-columns, Shape and DType object."""
    def __getitem__(self, T) -> T:
        ...


class F(Generic[T]):
    """ Generic class for functions, methods and classes. """
    def __getitem__(
        self, item: Callable[..., T]
    ) -> Union['F', Callable[..., T], T, Any]:
        return self


class Sub(Generic[T]):
    """ Return subset of an Array. """
    ...


class SP(Generic[T, D]):
    """ Station position arrays hold integer values. """
    ...


class Series(DType[T], Generic[T]):
    """ To reference the pandas `Series`_ object. """
    def __getitem__(self, item: T) -> 'Series':
        return self


class EDIO(Generic[T]):
    """
    EDIO stands for Electrical Data Interchange (EDI) Object.
    It is an EDI object built from `pycsamt` or `MTpy`.
    """
    def __getitem__(self, T: str | T) -> object:
        ...


class ZO(Generic[T]):
    """
    ZO stands for Impedance tensor Object. It is an Impedance
    tensor object built from `pycsamt.core.z.Z` or
    `mtpy.core.z.Z`. It is a 3D data with dimensions
    (n_freq, 2, 2).
    """
    def __getitem__(self, T: str | T) -> object:
        ...


class DataFrame(Series[T], Generic[T]):
    """ Type hint for a `pandas DataFrame`_ object. """
    def __getitem__(self, item: T) -> 'DataFrame':
        return self
