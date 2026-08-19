# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Matrix-oriented electromagnetic transfer-function scientific objects."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..core.base import MTBase
from .datatypes import DataTypeDefinition, get_emtf_datatype
from .estimates import StatisticalEstimate
from .validation import (
    normalize_channels,
    normalize_periods,
    normalize_tf_data,
)

__all__ = ["TransferFunction"]


@dataclass(repr=False)
class TransferFunction(MTBase):
    """Represent one matrix-valued electromagnetic transfer function.

    The canonical array layout is
    ``(n_period, n_output_channels, n_input_channels)``.  The matrix axes
    retain physical channel meaning, unlike EDI's disconnected component
    blocks.  No EDI or XML syntax is stored in this object.
    """

    name: str
    data: Any
    input_channels: tuple[str, ...] = field(default_factory=tuple)
    output_channels: tuple[str, ...] = field(default_factory=tuple)
    units: str | None = None
    periods: Any | None = None
    estimates: dict[str, StatisticalEstimate] = field(default_factory=dict)
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Normalize and validate matrix, channels, periods, and estimates."""
        raw_name = str(self.name).strip()
        if not raw_name:
            raise ValueError("transfer-function name must be non-empty")

        definition = get_emtf_datatype(raw_name)
        self.name = (
            definition.tag if definition is not None else raw_name.lower()
        )
        self.input_channels = normalize_channels(self.input_channels)
        self.output_channels = normalize_channels(self.output_channels)
        self.data = normalize_tf_data(
            self.data,
            n_output=len(self.output_channels),
            n_input=len(self.input_channels),
        )
        self.periods = normalize_periods(
            self.periods,
            n_periods=self.data.shape[0] if self.periods is not None else None,
        )
        self.attrs = dict(self.attrs or {})

        incoming = dict(self.estimates or {})
        self.estimates = {}
        for key, estimate in incoming.items():
            self.add_estimate(estimate, key=key, replace=True)

        if definition is not None:
            if self.units is None:
                self.units = definition.units
            self._validate_data_kind(definition)

    def _validate_data_kind(self, definition: DataTypeDefinition) -> None:
        if definition.data_kind == "real" and np.iscomplexobj(self.data):
            if np.any(np.imag(self.data) != 0.0):
                raise TypeError(
                    f"{definition.tag} is defined as real but complex data "
                    "with non-zero imaginary values were supplied"
                )
            self.data = np.real(self.data)
        elif definition.data_kind == "complex" and not np.iscomplexobj(
            self.data
        ):
            self.data = self.data.astype(complex)

    @property
    def definition(self) -> DataTypeDefinition | None:
        """Return the registry definition associated with this TF."""
        return get_emtf_datatype(self.name)

    @property
    def n_periods(self) -> int:
        """Number of frequency/period samples."""
        return int(self.data.shape[0])

    @property
    def n_output(self) -> int:
        """Number of output matrix rows."""
        return int(self.data.shape[1])

    @property
    def n_input(self) -> int:
        """Number of input matrix columns."""
        return int(self.data.shape[2])

    @property
    def shape(self) -> tuple[int, ...]:
        """Return the normalized TF data shape."""
        return tuple(self.data.shape)

    @property
    def frequency(self) -> np.ndarray | None:
        """Frequency vector in Hz derived from ``periods``."""
        if self.periods is None:
            return None
        return 1.0 / np.asarray(self.periods, dtype=float)

    def add_estimate(
        self,
        estimate: StatisticalEstimate,
        *,
        key: str | None = None,
        replace: bool = False,
    ) -> "TransferFunction":
        """Attach a statistical estimate to this transfer function."""
        if not isinstance(estimate, StatisticalEstimate):
            raise TypeError("estimate must be a StatisticalEstimate")
        est_key = str(key or estimate.kind).strip().lower()
        if not est_key:
            raise ValueError("estimate key must be non-empty")
        if est_key in self.estimates and not replace:
            raise ValueError(f"estimate already exists: {est_key}")

        if estimate.data.ndim > 0 and estimate.data.shape[0] != self.n_periods:
            raise ValueError(
                "estimate period axis does not match TF data: "
                f"{estimate.data.shape[0]} != {self.n_periods}"
            )
        self.estimates[est_key] = estimate
        return self

    def get_estimate(self, key: str) -> StatisticalEstimate | None:
        """Return an attached estimate by key, code, or semantic kind."""
        raw = str(key).strip()
        if not raw:
            return None
        direct = self.estimates.get(raw.lower())
        if direct is not None:
            return direct
        upper = raw.upper()
        lower = raw.lower()
        for estimate in self.estimates.values():
            if estimate.name == upper or estimate.kind == lower:
                return estimate
        return None

    def copy(self) -> "TransferFunction":
        """Return a detached copy of data, periods, estimates, and attrs."""
        return TransferFunction(
            name=self.name,
            data=np.array(self.data, copy=True),
            input_channels=tuple(self.input_channels),
            output_channels=tuple(self.output_channels),
            units=self.units,
            periods=(
                None
                if self.periods is None
                else np.array(self.periods, copy=True)
            ),
            estimates={
                key: estimate.copy()
                for key, estimate in self.estimates.items()
            },
            attrs=dict(self.attrs),
        )
