# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Mesh descriptions for EM inversion workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject

__all__ = ["InversionMesh"]


@dataclass
class InversionMesh(PyCSAMTObject, MetadataMixin):
    """Lightweight mesh/grid descriptor.

    Phase 1 uses this as metadata for external backends and result
    conversion. Numerical backends may attach native meshes in
    ``native`` while exposing common cell centres through ``x_centers``
    and ``z_centers``.
    """

    dimension: str = "1d"
    x_centers: Any = None
    z_centers: Any = None
    native: Any = field(default=None, repr=False)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.dimension = str(self.dimension).lower()
        self.x_centers = _array_or_none(self.x_centers)
        self.z_centers = _array_or_none(self.z_centers)
        self.validate()

    @classmethod
    def for_1d(cls, depths: Any) -> "InversionMesh":
        """Build a 1-D mesh descriptor from depth centres."""
        return cls(dimension="1d", x_centers=np.array([0.0]), z_centers=depths)

    def validate(self) -> None:
        if self.dimension not in {"1d", "2d", "3d"}:
            raise ValueError("dimension must be '1d', '2d', or '3d'.")
        if self.z_centers is not None and np.any(self.z_centers < 0):
            raise ValueError("z_centers must be positive downward.")


def _array_or_none(value: Any) -> np.ndarray | None:
    if value is None:
        return None
    return np.asarray(value, dtype=float)
