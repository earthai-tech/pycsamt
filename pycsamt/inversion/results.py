# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unified inversion result container."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject

from .mesh import InversionMesh

__all__ = ["InversionResult"]


@dataclass
class InversionResult(PyCSAMTObject, MetadataMixin):
    """Backend-neutral post-inversion result.

    The result stores common outputs while retaining backend-native objects
    for advanced users.
    """

    method: str
    dimension: str
    backend: str
    status: str = "success"
    model: Any = None
    mesh: InversionMesh | None = None
    data: Any = None
    predicted: Any = None
    rms: float = float("nan")
    objective: float = float("nan")
    n_iter: int = 0
    workdir: str | None = None
    files: dict[str, str] = field(default_factory=dict)
    native: Any = field(default=None, repr=False)
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.method = str(self.method).lower()
        self.dimension = str(self.dimension).lower()
        self.backend = str(self.backend).lower()
        self.status = str(self.status)
        self.rms = float(self.rms)
        self.objective = float(self.objective)
        self.n_iter = int(self.n_iter)
        self.files = dict(self.files or {})
        self.warnings = [str(w) for w in self.warnings]
        self.metadata = dict(self.metadata or {})

    @property
    def converged(self) -> bool:
        """Whether the backend reported a successful inversion."""
        return self.status in {"success", "converged", "prepared", "loaded"}

    def to_resistivity_model(self):
        """Convert the result to :class:`pycsamt.interp.ResistivityModel`.

        2-D result arrays are passed through directly. A recovered 1-D
        layered model is expanded to one column so the interpretation API can
        consume it uniformly.
        """
        from ..interp import ResistivityModel

        if hasattr(self.native, "rho_2d"):
            try:
                return ResistivityModel.from_occam2d(self.native)
            except Exception:
                pass
        if isinstance(self.model, dict) and "rho_2d" in self.model:
            rho = np.asarray(self.model["rho_2d"], dtype=float)
            x = np.asarray(self.model.get("x_centers", np.arange(rho.shape[1])), dtype=float)
            z = np.asarray(self.model.get("z_centers", np.arange(rho.shape[0])), dtype=float)
            return ResistivityModel.from_array(
                rho,
                x,
                z,
                station_x=np.asarray(self.model.get("station_x", x), dtype=float),
                station_names=list(self.model.get("station_names", [])) or None,
                method=f"{self.backend}:{self.method}",
                rms=self.rms,
            )
        if self.model is None:
            raise ValueError("result has no model to convert.")

        resistivity = np.asarray(
            getattr(self.model, "resistivities", getattr(self.model, "resistivity", [])),
            dtype=float,
        )
        thickness = np.asarray(
            getattr(self.model, "thicknesses", getattr(self.model, "thickness", [])),
            dtype=float,
        )
        if resistivity.size == 0 or thickness.size != resistivity.size - 1:
            raise ValueError("cannot convert model to ResistivityModel.")

        tops = np.r_[0.0, np.cumsum(thickness)]
        bottoms = np.r_[tops[1:], tops[-1] + thickness[-1]]
        z = 0.5 * (tops + bottoms)
        rho_2d = np.log10(resistivity).reshape(-1, 1)
        return ResistivityModel.from_array(
            rho_2d,
            np.array([0.0]),
            z,
            station_x=np.array([0.0]),
            station_names=["S000"],
            method=f"{self.backend}:{self.method}:1d",
            rms=self.rms,
        )

    def summary(self, *, max_fields: int | None = None) -> str:
        status = self.status
        rms = "nan" if not np.isfinite(self.rms) else f"{self.rms:.3g}"
        return (
            f"InversionResult(method={self.method!r}, dimension={self.dimension!r}, "
            f"backend={self.backend!r}, status={status!r}, rms={rms})"
        )
