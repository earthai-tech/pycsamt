# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Common base utilities for Z-related containers.

This module defines :class:`BaseEM`, a light-weight mixin-style
base class that provides:

- friendly ``__str__`` and ``__repr__`` summaries
- a robust way to infer the number of frequencies
- shallow / deep copy helpers
- safe subsetting along the frequency axis
- convenience predicates (e.g., ``has_freq``, ``has_errors``)
- a compact ``summary()`` string

It is intentionally minimal and non-intrusive so that
:class:`~pycsamt.z.resphase.ResPhase`,
:class:`~pycsamt.z.z.Z`, and
:class:`~pycsamt.z.tipper.Tipper`
can inherit from it without changes.

Notes
-----
``BaseEM`` does **not** enforce any attribute names. It relies on
convention if present:

- ``_freq``: 1-D array of frequencies.
- Arrays with first dimension equal to the number of frequencies
  (e.g., ``_z``, ``_z_err``, ``_tipper``, ``_resistivity``,
  ``_phase``) are detected for slicing.

If none of these exist, ``n_freq`` is reported as ``0``.

Examples
--------
>>> from pycsamt.z.base import BaseEM
>>> class Dummy(BaseEM):
...     def __init__(self, f, a):
...         super().__init__(name="dummy")
...         self._freq = f
...         self._z = a
...
>>> import numpy as np
>>> d = Dummy(np.array([10., 1.]), np.zeros((2, 2, 2), complex))
>>> print(d.n_freq)
2
>>> d2 = d.subset([1])   # take the second frequency only
>>> print(d2.n_freq)
1
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple, Union

import copy as _copy
import numpy as np

from ..log.logger import get_logger

_StrSeq = Sequence[str]
_Idx = Union[int, slice, np.ndarray, Sequence[int]]


@dataclass
class BaseEM:
    """
    Minimal, non-intrusive base for EM containers.

    Parameters
    ----------
    name : str, optional
        Human-friendly name used in ``__str__``.
    meta : dict, optional
        Arbitrary metadata attached to the object.

    Attributes
    ----------
    name : str or None
        Display name.
    meta : dict
        Arbitrary metadata.
    log
        Module logger (via :func:`get_logger`).

    Notes
    -----
    Subclasses can override ``_sliceable_predicate`` to customize
    which attributes should be sliced by :meth:`subset`.
    """
    name: Optional[str] = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.log = get_logger(type(self).__name__)


    # Introspection
    @property
    def n_freq(self) -> int:
        """
        Infer the number of frequencies.

        Returns
        -------
        int
            Length of ``_freq`` if present; otherwise the first
            dimension of the first array-like attribute found that
            looks frequency-stacked; ``0`` if unknown.
        """
        # Prefer explicit frequency vector
        f = getattr(self, "_freq", None)
        if isinstance(f, np.ndarray) and f.ndim == 1:
            return int(f.size)

        # Fallback: inspect known array fields
        for attr in ("_z", "_tipper", "_resistivity", "_phase"):
            a = getattr(self, attr, None)
            if isinstance(a, np.ndarray) and a.ndim >= 1:
                return int(a.shape[0])

        return 0

    @property
    def has_freq(self) -> bool:
        """Whether a valid ``_freq`` vector is available."""
        f = getattr(self, "_freq", None)
        return isinstance(f, np.ndarray) and f.ndim == 1 and f.size > 0

    @property
    def has_errors(self) -> bool:
        """Whether any error array appears to be present."""
        for attr in ("_z_err", "_tipper_err", "_resistivity_err", "_phase_err"):
            a = getattr(self, attr, None)
            if isinstance(a, np.ndarray):
                return True
        return False

    # Copy helpers

    def copy(self) -> "BaseEM":
        """Shallow copy of the instance."""
        return _copy.copy(self)

    def deepcopy(self) -> "BaseEM":
        """Deep copy of the instance."""
        return _copy.deepcopy(self)

    # Slicing (along frequency axis)
    def _sliceable_predicate(self, name: str, value: object, n: int) -> bool:
        """
        Decide if an attribute should be sliced by :meth:`subset`.

        Parameters
        ----------
        name : str
            Attribute name.
        value : object
            Attribute value.
        n : int
            Expected number of frequencies.

        Returns
        -------
        bool
            ``True`` if ``value`` should be sliced along axis-0.
        """
        if not isinstance(value, np.ndarray):
            return False
        if value.ndim == 0:
            return False
        # Slice arrays that match the frequency stack length
        return value.shape[0] == n

    def subset(self, indices: _Idx) -> "BaseEM":
        """
        Return a deep-copied subset along the frequency axis.

        Parameters
        ----------
        indices : int, slice, array-like of int
            Indices selecting along frequency axis (axis-0).

        Returns
        -------
        BaseEM
            A deep copy with all sliceable arrays subsetted.

        Notes
        -----
        - ``_freq`` is subset if present.
        - Arrays for which :meth:`_sliceable_predicate` returns
          ``True`` are subset along axis-0.
        """
        n = self.n_freq
        if n == 0:
            # nothing to slice
            return self.deepcopy()

        # Normalize indices to a numpy index
        idx = indices
        # Build new object by deep-copying, then slice selected attrs
        new = self.deepcopy()

        # Slice _freq, if any
        f = getattr(new, "_freq", None)
        if isinstance(f, np.ndarray) and f.ndim == 1 and f.size == n:
            new._freq = f[idx]

        # Generic slicing for frequency-aligned arrays
        for name, value in list(vars(new).items()):
            if self._sliceable_predicate(name, value, n):
                try:
                    setattr(new, name, value[idx])
                except Exception as exc:  # pragma: no cover
                    self.log.debug(
                        "Skipping slice for attr %s: %s", name, exc
                    )
        return new


    # String representation
    def _array_sig(self) -> Tuple[_StrSeq, _StrSeq]:
        """
        Collect array-shaped attribute signatures.

        Returns
        -------
        (names, sigs) : tuple of lists
            Names and compact ``shape@dtype`` strings.
        """
        names: list[str] = []
        sigs: list[str] = []
        for name, value in vars(self).items():
            if isinstance(value, np.ndarray):
                names.append(name)
                sigs.append(f"{tuple(value.shape)}@{value.dtype}")
        return names, sigs

    def summary(self) -> str:
        """
        Short, human-friendly summary.

        Returns
        -------
        str
            A compact multi-line description.
        """
        cls = type(self).__name__
        n = self.n_freq
        parts = [f"{cls}: {self.name or '-'}", f"n_freq: {n}"]
        if self.has_freq:
            f = getattr(self, "_freq")
            fmin = float(np.min(f)) if f.size else np.nan
            fmax = float(np.max(f)) if f.size else np.nan
            parts.append(f"f[Hz]: min={fmin:.6g}, max={fmax:.6g}")
        parts.append(f"errors: {'yes' if self.has_errors else 'no'}")
        an, asig = self._array_sig()
        if an:
            # Keep one-per-line but compact
            parts.append("arrays:")
            for name, sig in zip(an, asig):
                parts.append(f"  - {name}: {sig}")
        return "\n".join(parts)

    def __str__(self) -> str:  # pragma: no cover
        return self.summary()

    def __repr__(self) -> str:  # pragma: no cover
        cls = type(self).__name__
        return (
            f"{cls}(name={self.name!r}, n_freq={self.n_freq}, "
            f"errors={self.has_errors})"
        )
