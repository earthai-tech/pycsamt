# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""ZTEM-specific metadata that reuses the common airborne/EMTF model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ...core.base import CoreObject
from ...metadata import InstrumentMeta, SensorSpec, SiteMeta
from .constants import (
    ZTEM_COMMON_PROCESSED_BANDS_HZ,
    ZTEM_INPUT_CHANNELS,
    ZTEM_NOMINAL_OUTPUT_RATE_HZ,
    ZTEM_NOMINAL_TIME_SERIES_SAMPLING_RATE_HZ,
    ZTEM_OUTPUT_CHANNELS,
    ZTEM_PRACTICAL_FREQUENCY_RANGE_HZ,
    ZTEM_TYPICAL_FREQUENCY_COUNT,
)

__all__ = ["ZTEMSystemSpec", "ZTEMReferenceStation"]


@dataclass(repr=False)
class ZTEMSystemSpec(CoreObject):
    """Published/declared ZTEM system characteristics.

    The values are informative defaults. ZTEM frequency selection depends on
    platform speed, sampling, signal strength, and noise, so processed survey
    products may legitimately use a different subset or frequency grid.
    """

    practical_frequency_range_hz: tuple[float, float] = (
        ZTEM_PRACTICAL_FREQUENCY_RANGE_HZ
    )
    common_processed_bands_hz: tuple[tuple[float, float], ...] = (
        ZTEM_COMMON_PROCESSED_BANDS_HZ
    )
    typical_frequency_count: tuple[int, int] = ZTEM_TYPICAL_FREQUENCY_COUNT
    time_series_sampling_rate_hz: float = (
        ZTEM_NOMINAL_TIME_SERIES_SAMPLING_RATE_HZ
    )
    nominal_output_rate_hz: float = ZTEM_NOMINAL_OUTPUT_RATE_HZ
    input_channels: tuple[str, ...] = ZTEM_INPUT_CHANNELS
    output_channels: tuple[str, ...] = ZTEM_OUTPUT_CHANNELS
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        low, high = (float(v) for v in self.practical_frequency_range_hz)
        if not np.isfinite(low) or not np.isfinite(high):
            raise ValueError("practical frequency range must be finite")
        if low <= 0.0 or high <= low:
            raise ValueError(
                "practical frequency range must satisfy 0 < low < high"
            )
        self.practical_frequency_range_hz = (low, high)

        bands: list[tuple[float, float]] = []
        for band in self.common_processed_bands_hz:
            if len(band) != 2:
                raise ValueError("each common ZTEM band must contain 2 values")
            start, stop = (float(v) for v in band)
            if not np.isfinite(start) or not np.isfinite(stop):
                raise ValueError("common ZTEM bands must be finite")
            if start <= 0.0 or stop <= start:
                raise ValueError(
                    "common ZTEM bands must satisfy 0 < low < high"
                )
            bands.append((start, stop))
        self.common_processed_bands_hz = tuple(bands)

        if len(self.typical_frequency_count) != 2:
            raise ValueError("typical_frequency_count must contain min/max")
        count_low, count_high = (int(v) for v in self.typical_frequency_count)
        if count_low <= 0 or count_high < count_low:
            raise ValueError("typical frequency count must be positive")
        self.typical_frequency_count = (count_low, count_high)

        rate = float(self.time_series_sampling_rate_hz)
        if not np.isfinite(rate) or rate <= 0.0:
            raise ValueError("time_series_sampling_rate_hz must be positive")
        self.time_series_sampling_rate_hz = rate

        output_rate = float(self.nominal_output_rate_hz)
        if not np.isfinite(output_rate) or output_rate <= 0.0:
            raise ValueError("nominal_output_rate_hz must be positive")
        self.nominal_output_rate_hz = output_rate

        self.input_channels = tuple(
            str(value).strip() for value in self.input_channels
        )
        self.output_channels = tuple(
            str(value).strip() for value in self.output_channels
        )
        if self.input_channels != ZTEM_INPUT_CHANNELS:
            raise ValueError("ZTEM reference inputs must be Hx/Hy")
        if self.output_channels != ZTEM_OUTPUT_CHANNELS:
            raise ValueError("ZTEM airborne output must be Hz")
        self.attrs = dict(self.attrs or {})

    def practical_frequency_mask(self, frequency: Any) -> np.ndarray:
        """Return a mask for values inside the published practical band."""
        freq = np.asarray(frequency, dtype=float)
        low, high = self.practical_frequency_range_hz
        return np.isfinite(freq) & (freq >= low) & (freq <= high)

    def to_instrument_meta(
        self,
        *,
        serial: str | None = None,
        software_version: str = "",
    ) -> InstrumentMeta:
        """Return reusable instrument metadata without vendor guesses."""
        return InstrumentMeta(
            system="ZTEM",
            serial=serial,
            magnetic_sensor=SensorSpec(
                sensor_type="induction_coil",
                frequency_range=tuple(self.practical_frequency_range_hz),
                notes=(
                    "Single-axis airborne vertical magnetic receiver; "
                    "horizontal Hx/Hy reference fields are measured at a "
                    "fixed ground station."
                ),
            ),
            software_version=str(software_version),
            notes=(
                "Published ZTEM system characteristics; actual survey "
                "delivery metadata should override descriptive defaults."
            ),
        )


@dataclass(repr=False)
class ZTEMReferenceStation(CoreObject):
    """Fixed ground magnetic reference station used by ZTEM processing."""

    station_id: str | None = None
    site: SiteMeta | None = None
    magnetic_channels: tuple[str, str] = ZTEM_INPUT_CHANNELS
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if self.station_id is not None:
            text = str(self.station_id).strip()
            self.station_id = text or None
        if self.site is not None and not isinstance(self.site, SiteMeta):
            raise TypeError("site must be a SiteMeta or None")
        channels = tuple(
            str(value).strip() for value in self.magnetic_channels
        )
        if channels != ZTEM_INPUT_CHANNELS:
            raise ValueError("ZTEM reference magnetic channels must be Hx/Hy")
        self.magnetic_channels = channels
        self.attrs = dict(self.attrs or {})

    @property
    def preferred_id(self) -> str | None:
        """Return explicit station ID, then the SiteMeta preferred name."""
        if self.station_id:
            return self.station_id
        if self.site is not None:
            return self.site.preferred_name
        return None
