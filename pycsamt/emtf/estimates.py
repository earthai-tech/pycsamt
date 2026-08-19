# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Statistical estimate objects used by the format-neutral EMTF core."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import PyCSAMTObject

__all__ = ["StatisticalEstimate"]


@dataclass(repr=False)
class StatisticalEstimate(PyCSAMTObject):
    """Represent a frequency-indexed statistical estimate.

    The object deliberately stores the estimate *as supplied* instead of
    coercing all uncertainty information into the legacy ``z_err`` model.
    This is essential for later support of full covariance matrices.

    Parameters
    ----------
    name : str
        Short code, e.g. ``"VAR"`` or ``"RESIDCOV"``.
    data : array-like
        Numerical estimate. Real and complex arrays are supported.
    kind : str
        Semantic name such as ``"variance"`` or
        ``"inverse_signal_covariance"``.
    units : str or None
        Units when defined.
    attrs : dict
        Additional non-serialization-specific metadata.
    """

    name: str
    data: Any
    kind: str
    units: str | None = None
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Normalize and validate the estimate in place."""
        self.name = str(self.name).strip().upper()
        self.kind = str(self.kind).strip().lower()
        if not self.name:
            raise ValueError("statistical estimate name must be non-empty")
        if not self.kind:
            raise ValueError("statistical estimate kind must be non-empty")
        arr = np.asarray(self.data)
        if arr.dtype.kind not in "biufc":
            raise TypeError("statistical estimate data must be numeric")
        self.data = arr
        self.attrs = dict(self.attrs or {})

    @property
    def shape(self) -> tuple[int, ...]:
        """Return the stored array shape."""
        return tuple(self.data.shape)

    @property
    def is_complex(self) -> bool:
        """Return whether the stored estimate has a complex dtype."""
        return bool(np.iscomplexobj(self.data))

    def copy(self) -> "StatisticalEstimate":
        """Return a detached copy of the estimate."""
        return StatisticalEstimate(
            name=self.name,
            data=np.array(self.data, copy=True),
            kind=self.kind,
            units=self.units,
            attrs=dict(self.attrs),
        )
