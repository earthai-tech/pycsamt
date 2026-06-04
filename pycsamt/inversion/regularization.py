# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regularization controls for inversion workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ..api.property import MetadataMixin, PyCSAMTObject

__all__ = ["Regularization"]


@dataclass
class Regularization(PyCSAMTObject, MetadataMixin):
    """Backend-neutral regularization settings."""

    kind: str = "smooth"
    alpha_s: float = 1.0
    alpha_x: float = 1.0
    alpha_z: float = 1.0
    reference_weight: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.kind = str(self.kind).lower()
        self.alpha_s = float(self.alpha_s)
        self.alpha_x = float(self.alpha_x)
        self.alpha_z = float(self.alpha_z)
        self.reference_weight = float(self.reference_weight)
        self.validate()

    def validate(self) -> None:
        if self.kind not in {"none", "smooth", "damped", "blocky"}:
            raise ValueError("regularization kind must be none/smooth/damped/blocky.")
        for name in ("alpha_s", "alpha_x", "alpha_z", "reference_weight"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be non-negative.")
