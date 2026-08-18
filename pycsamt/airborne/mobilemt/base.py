# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""MobileMT-specific metadata that does not duplicate EMTF mathematics."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ...core.base import CoreObject
from ...metadata import InstrumentMeta, SensorSpec, SiteMeta
from .constants import (
    MOBILEMT_INPUT_CHANNELS,
    MOBILEMT_NOMINAL_FREQUENCY_RANGE_HZ,
    MOBILEMT_NOMINAL_MAX_WINDOWS,
    MOBILEMT_NOMINAL_SAMPLING_RATE_HZ,
    MOBILEMT_OUTPUT_CHANNELS,
)

__all__ = ["MobileMTSystemSpec", "MobileMTReferenceStation"]


@dataclass(repr=False)
class MobileMTSystemSpec(CoreObject):
    """Published/declared MobileMT system characteristics.

    These values are descriptive metadata, not hard validation limits. A
    processed product may legitimately contain fewer windows or a subset of
    the nominal frequency range.
    """

    nominal_frequency_range_hz: tuple[float, float] = (
        MOBILEMT_NOMINAL_FREQUENCY_RANGE_HZ
    )
    nominal_max_frequency_windows: int = MOBILEMT_NOMINAL_MAX_WINDOWS
    nominal_sampling_rate_hz: float = MOBILEMT_NOMINAL_SAMPLING_RATE_HZ
    input_channels: tuple[str, ...] = MOBILEMT_INPUT_CHANNELS
    output_channels: tuple[str, ...] = MOBILEMT_OUTPUT_CHANNELS
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        low, high = (float(v) for v in self.nominal_frequency_range_hz)
        if not np.isfinite(low) or not np.isfinite(high):
            raise ValueError("nominal frequency range must be finite")
        if low <= 0.0 or high <= low:
            raise ValueError(
                "nominal frequency range must satisfy 0 < low < high"
            )
        self.nominal_frequency_range_hz = (low, high)

        nwin = int(self.nominal_max_frequency_windows)
        if nwin <= 0:
            raise ValueError("nominal_max_frequency_windows must be positive")
        self.nominal_max_frequency_windows = nwin

        rate = float(self.nominal_sampling_rate_hz)
        if not np.isfinite(rate) or rate <= 0.0:
            raise ValueError("nominal_sampling_rate_hz must be positive")
        self.nominal_sampling_rate_hz = rate

        self.input_channels = tuple(
            str(v).strip() for v in self.input_channels
        )
        self.output_channels = tuple(
            str(v).strip() for v in self.output_channels
        )
        if self.input_channels != MOBILEMT_INPUT_CHANNELS:
            raise ValueError(
                "MobileMT processed admittance inputs must be Ex/Ey"
            )
        if self.output_channels != MOBILEMT_OUTPUT_CHANNELS:
            raise ValueError(
                "MobileMT processed admittance outputs must be Hx/Hy/Hz"
            )
        self.attrs = dict(self.attrs or {})

    def nominal_frequency_mask(self, frequency: Any) -> np.ndarray:
        """Return a mask for values inside the published nominal range."""
        freq = np.asarray(frequency, dtype=float)
        low, high = self.nominal_frequency_range_hz
        return np.isfinite(freq) & (freq >= low) & (freq <= high)

    def to_instrument_meta(
        self,
        *,
        serial: str | None = None,
        software_version: str = "",
    ) -> InstrumentMeta:
        """Return a reusable :class:`InstrumentMeta` description.

        The sensor descriptions reflect only the published system-level
        geometry. Manufacturer model numbers and calibration details are not
        invented.
        """
        band = tuple(self.nominal_frequency_range_hz)
        return InstrumentMeta(
            system="MobileMT",
            serial=serial,
            magnetic_sensor=SensorSpec(
                sensor_type="induction_coil",
                frequency_range=band,
                notes="Three orthogonal airborne magnetic coils",
            ),
            electric_sensor=SensorSpec(
                sensor_type="electrode",
                frequency_range=band,
                notes="Ground horizontal electric reference station",
            ),
            software_version=str(software_version),
            notes=(
                "Nominal processed MobileMT specification; native delivery "
                "metadata should override these descriptive defaults."
            ),
        )


@dataclass(repr=False)
class MobileMTReferenceStation(CoreObject):
    """Ground electric reference-station metadata for MobileMT processing.

    This object intentionally describes only the processed reference station.
    It does not invent a raw electrode-file schema or model the additional
    internal reference electrode pairs used by proprietary processing.
    """

    station_id: str | None = None
    site: SiteMeta | None = None
    electric_channels: tuple[str, str] = MOBILEMT_INPUT_CHANNELS
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if self.station_id is not None:
            text = str(self.station_id).strip()
            self.station_id = text or None
        if self.site is not None and not isinstance(self.site, SiteMeta):
            raise TypeError("site must be a SiteMeta or None")
        channels = tuple(str(v).strip() for v in self.electric_channels)
        if channels != MOBILEMT_INPUT_CHANNELS:
            raise ValueError(
                "MobileMT reference electric channels must be Ex/Ey"
            )
        self.electric_channels = channels
        self.attrs = dict(self.attrs or {})

    @property
    def preferred_id(self) -> str | None:
        """Return explicit station ID, then SiteMeta preferred name."""
        if self.station_id:
            return self.station_id
        if self.site is not None:
            return self.site.preferred_name
        return None
