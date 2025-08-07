# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""pycsamt.zonge.tensor
Helpers that *promote* the four scalar component columns
(Zxx, Zxy, Zyx, Zyy – and their derived ρ / φ / σ) to true
``(n_freq, 2, 2)`` tensors.

Keeps all heavy maths out of the high-level objects.
"""
from __future__ import annotations

from dataclasses import dataclass #, field
from typing      import Dict, Sequence, Any, Final
import numpy as np

from ..constants import PI                                           # noqa: F401

class TensorBuildError (Exception):
    """ Tensor base exceptions"""
    pass 

@dataclass(slots=True)
class ImpedanceTensor:
    """Container holding *all* 2×2 stacks for a station / site.

    Attributes
    ----------
    z : ndarray, complex
        Shape ``(n_freq, 2, 2)``.
    rho, phase : ndarray, float
        Same shape as *z*.  Phase is **degrees**.
    z_err, rho_err, phase_err : ndarray | None
        Optional 1-σ uncertainties.  When *None* no error-prop
        was provided for the corresponding field.
    """
    z:          np.ndarray
    rho:        np.ndarray
    phase:      np.ndarray
    z_err:      np.ndarray | None = None
    rho_err:    np.ndarray | None = None
    phase_err:  np.ndarray | None = None


    def _slice(self, i: int, j: int, field: str = "z") -> np.ndarray:
        return getattr(self, field)[:, i, j]

    # impedance
    @property 
    def z_xx(self): return self._slice(0, 0, "z")
    @property 
    def z_xy(self): return self._slice(0, 1, "z")
    @property 
    def z_yx(self): return self._slice(1, 0, "z")
    @property 
    def z_yy(self): return self._slice(1, 1, "z")
    # resistivity
    @property 
    def rho_xx(self): return self._slice(0, 0, "rho")
    @property 
    def rho_xy(self): return self._slice(0, 1, "rho")
    @property 
    def rho_yx(self): return self._slice(1, 0, "rho")
    @property 
    def rho_yy(self): return self._slice(1, 1, "rho")
    # phase (deg)
    @property 
    def ph_xx(self): return self._slice(0, 0, "phase")
    @property 
    def ph_xy(self): return self._slice(0, 1, "phase")
    @property 
    def ph_yx(self): return self._slice(1, 0, "phase")
    @property 
    def ph_yy(self): return self._slice(1, 1, "phase")

    def __str__(self) -> str:                       # user-friendly
        blocks = ", ".join(
            k for k, v in (("Z", self.z),
                            ("ρ", self.rho),
                            ("ϕ", self.phase)) if v is not None
        )
        nf, ns = self.z.shape[:2] if self.z is not None else \
                  self.rho.shape[:2] if self.rho is not None else \
                  self.phase.shape[:2]
        return f"ImpedanceTensor({blocks}; n_freq={nf}, n_station={ns})"

    def __repr__(self) -> str:                      # dev-oriented
        cls = self.__class__.__name__
        nf, ns = (self.z or self.rho or self.phase).shape[:2]
        flags = "/".join(flag for flag in
                         ("Z"  if self.z     is not None else "-",
                          "ρ"  if self.rho   is not None else "-",
                          "ϕ"  if self.phase is not None else "-"))
        return (f"<{cls} {flags} "
                f"(n_freq={nf}, n_station={ns}) at 0x{id(self):x}>")



class TensorFactory:
    """
    Tiny helper that stacks four *xx/xy/yx/yy* component vectors into the
    canonical ``(n_freq, 1, 2, 2)`` shape expected by
    :class:`~pycsamt.zonge.Z`.
    """

    @classmethod
    def build(
        cls,
        *,
        # any **one** (or more) of the following can be supplied
        z:        Dict[str, Sequence[Any]] | None = None,
        rho:      Dict[str, Sequence[Any]] | None = None,
        phase:    Dict[str, Sequence[Any]] | None = None,
        # optional one-sigma errors
        z_err:    Dict[str, Sequence[Any]] | None = None,
        rho_err:  Dict[str, Sequence[Any]] | None = None,
        ph_err:   Dict[str, Sequence[Any]] | None = None,
        dtype_z:  Any = complex,
        dtype_f:  Any = float,
    ) -> "ImpedanceTensor":
        """
        Assemble an :class:`ImpedanceTensor` from the components that are
        **actually** available.  Any non-supplied block (and its error)
        will be set to ``None``.

        Example
        -------
        >>> t = TensorFactory.build(rho=rho_dict)   # only ρₐ known
        >>> assert t.rho is not None and t.z is None
        """
        if z is rho is phase is None:
            raise TensorBuildError("At least one of z / rho / phase is required")

        # figure out grid size from the first non-None block
        proto = next(b for b in (z, rho, phase) if b is not None)
        cls._check_dict(proto, "prototype")
        n_freq = len(next(iter(proto.values())))

        # helper to stack or return None
        def _maybe_stack(block, kind, dtype):
            if block is None:
                return None
            cls._check_dict(block, kind)
            cls._assert_equal_len(block, n_freq, kind)
            return cls._stack(block, dtype)

        z_t     = _maybe_stack(z,        "z",     dtype_z)
        rho_t   = _maybe_stack(rho,      "rho",   dtype_f)
        phase_t = _maybe_stack(phase,    "phase", dtype_f)
        ze_t    = _maybe_stack(z_err,    "z_err", dtype_f)
        re_t    = _maybe_stack(rho_err,  "rho_err", dtype_f)
        pe_t    = _maybe_stack(ph_err,   "ph_err", dtype_f)

        return ImpedanceTensor(z_t, rho_t, phase_t, ze_t, re_t, pe_t)


    _REQ_KEYS = ("xx", "xy", "yx", "yy")

    @classmethod
    def _check_dict(cls, d: Dict[str, Sequence[Any]] | None, label: str) -> None:
        if d is None:
            return
        missing = [k for k in cls._REQ_KEYS if k not in d]
        if missing:
            raise TensorBuildError(f"{label}: missing keys {missing}")

    @staticmethod
    def _assert_equal_len(d: Dict[str, Sequence[Any]], n: int, label: str) -> None:
        for k, v in d.items():
            if len(v) != n:
                raise TensorBuildError(
                    f"{label}.{k}: length {len(v)} ≠ reference {n}"
                )

    @staticmethod
    def _stack(block: Dict[str, Sequence[Any]] | None, dtype) -> np.ndarray:
        if block is None:
            return None     # type: ignore[return-value]
        xx, xy, yx, yy = (np.asarray(block[k], dtype=dtype) for k in (
            "xx", "xy", "yx", "yy"))
        return np.stack(  # (n_freq, 2, 2)
            (np.stack((xx, xy), axis=-1),
             np.stack((yx, yy), axis=-1)),
            axis=-2
        )[..., None, :, :]        # → (n_freq, 1, 2, 2)

class _TensorFactory:
    """Build :class:`ImpedanceTensor` **from flat column arrays**.

    All arrays must be *1-D* and of **equal length** *(n_freq)*.

    Examples
    --------
    ```python
    it = TensorFactory.build(
        z       = dict(xx=z_xx, xy=z_xy, yx=z_yx, yy=z_yy),
        rho     = dict(xx=rho_xx, xy=rho_xy, yx=rho_yx, yy=rho_yy),
        phase   = dict(xx=ph_xx, xy=ph_xy, yx=ph_yx, yy=ph_yy),
        z_err   = dict(xx=err_xx, xy=err_xy, yx=err_yx, yy=err_yy),
    )
    # → ImpedanceTensor ready for further processing
    ```
    """

    _COMP_ORDER: Final[tuple[str, str, str, str]] = ("xx", "xy", "yx", "yy")

    @classmethod
    def build(
        cls,
        *,
        z:        Dict[str, Sequence[Any]] | None = None,
        rho:      Dict[str, Sequence[Any]] | None = None,
        phase:    Dict[str, Sequence[Any]] | None = None,
        z_err:    Dict[str, Sequence[Any]] | None = None,
        rho_err:  Dict[str, Sequence[Any]] | None = None,
        ph_err:   Dict[str, Sequence[Any]] | None = None,
        freq:     Sequence[float] | np.ndarray | None = None,
        dtype_z:  Any = complex,
        dtype_f:  Any = float,
    ) -> "ImpedanceTensor":
        """
        Flexible constructor for :class:`ImpedanceTensor`.

        Every block (*z*, ρ, φ, their σ) is **optional**.  When a block is
        provided it must contain the four components ``{'xx','xy','yx','yy'}``
        and all components must share the same *n_freq*.

        The factory refuses to build an *empty* tensor (at least one of
        *z*, ρ or φ must be given).
        """
        supplied = {tag: blk for tag, blk in
                    dict(z=z, rho=rho, phase=phase,
                         z_err=z_err, rho_err=rho_err, ph_err=ph_err).items()
                    if blk is not None}

        if not supplied:
            raise TensorBuildError(
                "Nothing to build – provide at least one data block")

        # determine n_freq from the first block 
        tag0, first = next(iter(supplied.items()))
        cls._check_dict(first, tag0)
        n_freq = len(next(iter(first.values())))
        cls._assert_equal_len(first, n_freq, tag0)

        # validate all other supplied blocks 
        for tag, mapping in supplied.items():
            if tag == tag0:
                continue
            cls._check_dict(mapping, tag)
            cls._assert_equal_len(mapping, n_freq, tag)

        # stack everything (missing → None) 
        stack = lambda blk, dt: None if blk is None else cls._stack(blk, dt)

        return ImpedanceTensor(
            z        = stack(z,       dtype_z),
            rho      = stack(rho,     dtype_f),
            phase    = stack(phase,   dtype_f),
            z_err    = stack(z_err,   dtype_f),
            rho_err  = stack(rho_err, dtype_f),
            ph_err   = stack(ph_err,  dtype_f),
            freq     = None if freq is None else np.asarray(freq, dtype=float),
        )

    @staticmethod
    def _check_dict(d: Dict[str, Sequence[Any]], name: str) -> None:
        missing = [c for c in TensorFactory._COMP_ORDER if c not in d]
        if missing:
            raise TensorBuildError(f"{name!s} dict missing keys: {missing}")

    @staticmethod
    def _assert_equal_len(
        d: Dict[str, Sequence[Any]],
        n: int,
        name: str
    ) -> None:
        for k, v in d.items():
            if len(v) != n:
                raise TensorBuildError(
                    f"Inconsistent length in {name!s}.{k}: "
                    f"expected {n}, got {len(v)}"
                )

    @staticmethod
    def _stack(
        d: Dict[str, Sequence[Any]] | None,
        dtype: Any = float
    ) -> np.ndarray:
        if d is None:
            return None  # type: ignore[return-value]
        a = np.column_stack([np.asarray(d[k], dtype=dtype)
                             for k in TensorFactory._COMP_ORDER])
        n = a.shape[0]
        return a.reshape(n, 2, 2)
