# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Model containers used by :mod:`pycsamt.inversion`."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject

__all__ = ["StartingModel", "ReferenceModel"]


@dataclass
class StartingModel(PyCSAMTObject, MetadataMixin):
    """Starting or recovered layered-earth model.

    Resistivities are linear ohm-m values. The final layer is a halfspace,
    so ``len(thicknesses)`` must equal ``len(resistivities) - 1``.
    """

    resistivities: Any
    thicknesses: Any
    name: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.resistivities = np.asarray(self.resistivities, dtype=float)
        self.thicknesses = np.asarray(self.thicknesses, dtype=float)
        self.validate()

    @classmethod
    def default(cls, n_layers: int = 4) -> "StartingModel":
        """Return a conservative layered starting model."""
        n_layers = int(n_layers)
        if n_layers < 2:
            raise ValueError("n_layers must be >= 2.")
        resistivities = np.full(n_layers, 100.0, dtype=float)
        thicknesses = np.geomspace(100.0, 2000.0, n_layers - 1)
        return cls(resistivities, thicknesses, name="default")

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "StartingModel":
        """Build from singular or plural key names."""
        return cls(
            data.get("resistivity", data.get("resistivities")),
            data.get("thickness", data.get("thicknesses")),
            name=str(data.get("name", "")),
            metadata=dict(data.get("metadata", {})),
        )

    @classmethod
    def coerce(cls, value: Any, *, n_layers: int = 4) -> "StartingModel":
        """Return *value* as a :class:`StartingModel`."""
        if value is None:
            return cls.default(n_layers=n_layers)
        if isinstance(value, cls):
            return value
        if isinstance(value, dict):
            return cls.from_dict(value)
        resistivity = getattr(value, "resistivity", None)
        if resistivity is None:
            resistivity = getattr(value, "resistivities", None)
        thickness = getattr(value, "thickness", None)
        if thickness is None:
            thickness = getattr(value, "thicknesses", None)
        return cls(resistivity, thickness, name=getattr(value, "name", ""))

    @property
    def n_layers(self) -> int:
        return int(self.resistivities.size)

    @property
    def depths(self) -> np.ndarray:
        """Top-of-layer depths in metres."""
        return np.r_[0.0, np.cumsum(self.thicknesses)]

    def to_layered_model(self):
        """Return the existing :class:`pycsamt.forward.LayeredModel`."""
        from ..forward import LayeredModel

        return LayeredModel(
            resistivity=self.resistivities.copy(),
            thickness=self.thicknesses.copy(),
            name=self.name,
        )

    def validate(self) -> None:
        if self.resistivities.ndim != 1 or self.thicknesses.ndim != 1:
            raise ValueError("resistivities and thicknesses must be 1-D.")
        if self.resistivities.size < 2:
            raise ValueError("at least two layers are required.")
        if self.thicknesses.size != self.resistivities.size - 1:
            raise ValueError("len(thicknesses) must be len(resistivities)-1.")
        if np.any(self.resistivities <= 0):
            raise ValueError("resistivities must be strictly positive.")
        if np.any(self.thicknesses <= 0):
            raise ValueError("thicknesses must be strictly positive.")


ReferenceModel = StartingModel
