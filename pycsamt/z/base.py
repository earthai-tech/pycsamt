# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple, Union, Any

import copy as _copy
import numpy as np

from ..log.logger import get_logger
from ..exceptions import ZError

_StrSeq = Sequence[str]
_Idx = Union[int, slice, np.ndarray, Sequence[int]]


@dataclass
class EMBase:
    r"""
    Minimal, mixin-style base class for EM containers.

    ``EMBase`` provides a small set of uniform behaviors used by
    higher-level containers such as :class:`~pycsamt.z.z.Z`,
    :class:`~pycsamt.z.resphase.ResPhase`, and
    :class:`~pycsamt.z.tipper.Tipper`.  It centralizes frequency
    handling, safe slicing along the frequency axis, shape checks,
    light verbosity control, and compact summaries.

    The class is deliberately conservative: it makes as few
    assumptions as possible about attribute names and shapes, while
    still enabling useful defaults for common containers.

    Parameters
    ----------
    name : str, optional
        Human-friendly display name.  Used by :meth:`summary` and
        representations.
    meta : dict, optional
        Arbitrary metadata attached to the instance.  Not used by
        core logic.
    verbose : int, default 0
        Controls the object logger level.  ``0`` → WARNING,
        ``1`` → INFO, ``>=2`` → DEBUG.  See :meth:`set_verbose`.

    Attributes
    ----------
    name : str or None
        Display name.
    meta : dict
        User metadata stored verbatim.
    verbose : int
        Current verbosity level.  See :meth:`set_verbose`.
    log : logging.Logger
        Per-class logger created via :func:`~pycsamt.log.logger.
        get_logger`.
    freq : ndarray or None
        One-dimensional frequency vector in Hz.  The property enforces
        finite, strictly positive values.  Stored internally as
        ``_freq`` when set.
    n_freq : int
        Inferred number of frequencies.  If ``_freq`` exists, its
        length is returned.  Otherwise, the first dimension of any
        known array attribute (e.g., ``_z``, ``_tipper``,
        ``_resistivity``, ``_phase``) is used.  Returns ``0`` if
        unknown.

    Notes
    -----
    **Design.**
    The class aims to be a small, composable layer.  It avoids
    imposing a strict attribute schema so that existing containers
    can adopt it without refactors.  Uniform frequency handling is
    offered through the :pyattr:`freq` property.

    **Slicing.**
    :meth:`subset` returns a deep-copied view where every attribute
    whose first dimension equals :pyattr:`n_freq` is sliced along
    axis-0.  This covers typical arrays like ``_z``, ``_tipper``,
    ``_phase``, and their error fields.  Subclasses can override
    :meth:`_sliceable_predicate` for fine-grained control.

    **Validation.**
    :meth:`validate_shapes` checks that all frequency-aligned arrays
    share the same leading dimension (:pyattr:`n_freq`).  It raises
    :class:`~pycsamt.exceptions.ZError` on mismatches.

    **Verbosity.**
    :meth:`set_verbose` updates the instance logger level only.  It
    does not alter global logging configuration.

    **Back-compat.**
    ``BaseEM`` remains available as an alias to ``EMBase`` for
    compatibility with older imports.

    Examples
    --------
    Minimal subclass that stores a frequency vector and a tensor:

    >>> import numpy as np
    >>> from pycsamt.z.base import BaseEM
    >>> class Dummy(BaseEM):
    ...     def __init__(self, f, z, *, name=None, verbose=0):
    ...         super().__init__(name=name, verbose=verbose)
    ...         self.freq = f
    ...         self._z = np.asarray(z, complex)
    ...
    >>> f = np.array([10.0, 1.0])
    >>> z = np.zeros((2, 2, 2), complex)
    >>> d = Dummy(f, z, name="site-A", verbose=1)
    >>> d.n_freq
    2
    >>> print(d.summary())  # doctest: +ELLIPSIS
    EMBase: site-A
    n_freq: 2
    f[Hz]: min=1, max=10
    errors: no
    arrays:
      - _z: (2, 2, 2)@complex128

    Frequency-axis subsetting with a boolean mask:

    >>> m = np.array([True, False])
    >>> d2 = d.subset(m)
    >>> d2.n_freq
    1

    Validating shapes (raises on mismatch):

    >>> d._phase = np.zeros((1, 2, 2))
    >>> from pycsamt.exceptions import ZError
    >>> _ = d.validate_shapes()  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
        ...
    ZError: freq-aligned arrays mismatch n_freq. ...

    See Also
    --------
    pycsamt.z.z.Z :
        Impedance tensor container building on resistivity/phase.
    pycsamt.z.resphase.ResPhase :
        Resistivity/phase container with Z linkage and errors.
    pycsamt.z.tipper.Tipper :
        Complex tipper container with rotation and arrow tools.

    References
    ----------
    .. [1] Chave, A. D., & Jones, A. G. (2012).
           *The Magnetotelluric Method: Theory and Practice*.
           Cambridge University Press.
    .. [2] Simpson, F., & Bahr, K. (2005).
           *Practical Magnetotellurics*. Cambridge University Press.
    """

    name: Optional[str] = None
    meta: dict[str, Any] = field(default_factory=dict)
    verbose: int = 0

    def __post_init__(self) -> None:
        # subclasses may not call super(); be defensive elsewhere
        self.log = get_logger(type(self).__name__)
        self.set_verbose(self.verbose)

    def set_verbose(self, level: int = 0) -> "EMBase":
        self.verbose = int(level)
        lg = getattr(self, "log", None)
        if lg is None:
            lg = get_logger(type(self).__name__)
            self.log = lg
        if self.verbose <= 0:
            lg.setLevel("WARNING")
        elif self.verbose == 1:
            lg.setLevel("INFO")
        else:
            lg.setLevel("DEBUG")
        return self

    @property
    def freq(self) -> Optional[np.ndarray]:
        f = getattr(self, "_freq", None)
        return None if f is None else f
    
    @freq.setter
    def freq(self, f: Optional[Sequence[float]]) -> None:
        if f is None:
            self._freq = None
            return
        ff = np.asarray(f, dtype=float)
        if ff.ndim != 1:
            raise ZError("freq must be 1-D")
        if not np.all(np.isfinite(ff)):
            raise ZError("freq must be finite")
        if np.any(ff <= 0.0):
            raise ZError("freq must be > 0")
        self._freq = ff
    

    @property
    def n_freq(self) -> int:
        f = getattr(self, "_freq", None)
        if isinstance(f, np.ndarray) and f.ndim == 1:
            return int(f.size)
        for attr in ("_z", "_tipper", "_resistivity", "_phase"):
            a = getattr(self, attr, None)
            if isinstance(a, np.ndarray) and a.ndim >= 1:
                return int(a.shape[0])
        return 0

    def __len__(self) -> int:
        return self.n_freq

    @property
    def has_freq(self) -> bool:
        f = getattr(self, "_freq", None)
        return isinstance(f, np.ndarray) and f.ndim == 1 and f.size > 0

    @property
    def has_errors(self) -> bool:
        for attr in (
            "_z_err",
            "_tipper_err",
            "_resistivity_err",
            "_phase_err",
        ):
            a = getattr(self, attr, None)
            if isinstance(a, np.ndarray):
                return True
        return False


    def copy(self) -> "EMBase":
        return _copy.copy(self)

    def deepcopy(self) -> "EMBase":
        return _copy.deepcopy(self)

    def _sliceable_predicate(
        self,
        name: str,
        value: object,
        n: int,
    ) -> bool:
        if not isinstance(value, np.ndarray):
            return False
        if value.ndim == 0:
            return False
        return value.shape[0] == n

    def subset(self, indices: _Idx) -> "EMBase":
        n = self.n_freq
        if n == 0:
            return self.deepcopy()
        idx = indices
        new = self.deepcopy()

        f = getattr(new, "_freq", None)
        if isinstance(f, np.ndarray) and f.ndim == 1 and f.size == n:
            new._freq = f[idx]

        for name, value in list(vars(new).items()):
            if self._sliceable_predicate(name, value, n):
                try:
                    setattr(new, name, value[idx])
                except Exception as exc:  # pragma: no cover
                    lg = getattr(
                        self,
                        "log",
                        get_logger(type(self).__name__),
                    )
                    lg.debug("skip slice %s: %s", name, exc)
        return new

    # small alias that reads nicely in notebooks
    def select(self, indices: _Idx) -> "EMBase":
        return self.subset(indices)

    def validate_shapes(self) -> None:
        n = self.n_freq
        if n == 0:
            return
        bad: list[str] = []
        for name, value in vars(self).items():
            if not isinstance(value, np.ndarray):
                continue
            if value.ndim == 0:
                continue
            if value.shape[0] != n:
                bad.append(f"{name}:{tuple(value.shape)}")
        if bad:
            raise ZError(
                (
                    "freq-aligned arrays mismatch n_freq. "
                    f"n={n}; bad={', '.join(bad)}"
                )
            )

    def _array_sig(self) -> Tuple[_StrSeq, _StrSeq]:
        names: list[str] = []
        sigs: list[str] = []
        for name, value in vars(self).items():
            if isinstance(value, np.ndarray):
                names.append(name)
                sigs.append(f"{tuple(value.shape)}@{value.dtype}")
        return names, sigs

    def summary(self) -> str:
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


# keep backward compat name
BaseEM = EMBase
