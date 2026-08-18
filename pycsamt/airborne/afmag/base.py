# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""AFMAG-family metadata built on the common airborne model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ...core.base import CoreObject
from ...metadata import InstrumentMeta, SensorSpec, SiteMeta
from .constants import (
    AFMAG_ORIGINAL_COIL_COUNT,
    AFMAG_ORIGINAL_COIL_SEPARATION_DEG,
    AFMAG_ORIGINAL_COIL_TILT_DEG,
    AFMAG_ORIGINAL_HISTORICAL_BAND_HZ,
    AFMAG_ORIGINAL_TYPICAL_FREQUENCIES_HZ,
    AFMAG_TENSOR_INPUT_CHANNELS,
    AFMAG_TENSOR_NOMINAL_SAMPLING_RATE_HZ,
    AFMAG_TENSOR_OUTPUT_CHANNELS,
    AFMAG_TENSOR_PRACTICAL_FREQUENCY_RANGE_HZ,
    AFMAG_TENSOR_REFERENCE_CHANNELS,
    AFMAG_TENSOR_TYPICAL_FREQUENCY_COUNT,
)

__all__ = [
    "OriginalAFMAGSystemSpec",
    "AirMtSystemSpec",
    "AFMAGReferenceStation",
]


@dataclass(repr=False)
class OriginalAFMAGSystemSpec(CoreObject):
    """Descriptive characteristics of the historical comparator AFMAG.

    These values describe the original airborne implementation and are not
    parser constraints.  Historical systems reported a comparator deflection
    proportional to the polarization-plane tilt rather than a digital tensor.
    """

    historical_frequency_band_hz: tuple[float, float] = (
        AFMAG_ORIGINAL_HISTORICAL_BAND_HZ
    )
    typical_frequencies_hz: tuple[float, ...] = (
        AFMAG_ORIGINAL_TYPICAL_FREQUENCIES_HZ
    )
    coil_count: int = AFMAG_ORIGINAL_COIL_COUNT
    coil_tilt_deg: float = AFMAG_ORIGINAL_COIL_TILT_DEG
    coil_separation_deg: float = AFMAG_ORIGINAL_COIL_SEPARATION_DEG
    digital_recording: bool = False
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        low, high = (float(v) for v in self.historical_frequency_band_hz)
        if not np.isfinite(low) or not np.isfinite(high):
            raise ValueError("historical AFMAG frequency band must be finite")
        if low <= 0.0 or high <= low:
            raise ValueError("historical AFMAG frequency band is invalid")
        self.historical_frequency_band_hz = (low, high)

        freq = tuple(float(v) for v in self.typical_frequencies_hz)
        if not freq or any((not np.isfinite(v) or v <= 0.0) for v in freq):
            raise ValueError("typical AFMAG frequencies must be positive")
        self.typical_frequencies_hz = freq

        self.coil_count = int(self.coil_count)
        if self.coil_count != 2:
            raise ValueError("original comparator AFMAG uses two coils")
        for name in ("coil_tilt_deg", "coil_separation_deg"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive")
            setattr(self, name, value)
        self.digital_recording = bool(self.digital_recording)
        self.attrs = dict(self.attrs or {})

    def to_instrument_meta(
        self,
        *,
        serial: str | None = None,
        software_version: str = "",
    ) -> InstrumentMeta:
        """Return reusable instrument metadata for historical AFMAG."""
        return InstrumentMeta(
            system="AFMAG (original comparator)",
            serial=serial,
            magnetic_sensor=SensorSpec(
                sensor_type="induction_coil",
                frequency_range=self.historical_frequency_band_hz,
                notes=(
                    "Two airborne comparison coils; historical output was "
                    "proportional to magnetic polarization-plane tilt."
                ),
            ),
            software_version=str(software_version),
            notes=(
                "Historical AFMAG descriptive metadata; actual archival "
                "records may use instrument-specific deflection scales."
            ),
        )


@dataclass(repr=False)
class AirMtSystemSpec(CoreObject):
    """Published characteristics of the tensor AFMAG/AirMt generation."""

    practical_frequency_range_hz: tuple[float, float] = (
        AFMAG_TENSOR_PRACTICAL_FREQUENCY_RANGE_HZ
    )
    typical_frequency_count: tuple[int, int] = (
        AFMAG_TENSOR_TYPICAL_FREQUENCY_COUNT
    )
    time_series_sampling_rate_hz: float = (
        AFMAG_TENSOR_NOMINAL_SAMPLING_RATE_HZ
    )
    input_channels: tuple[str, ...] = AFMAG_TENSOR_INPUT_CHANNELS
    output_channels: tuple[str, ...] = AFMAG_TENSOR_OUTPUT_CHANNELS
    reference_channels: tuple[str, ...] = AFMAG_TENSOR_REFERENCE_CHANNELS
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        low, high = (float(v) for v in self.practical_frequency_range_hz)
        if not np.isfinite(low) or not np.isfinite(high):
            raise ValueError("AirMt frequency range must be finite")
        if low <= 0.0 or high <= low:
            raise ValueError("AirMt frequency range must satisfy low < high")
        self.practical_frequency_range_hz = (low, high)

        if len(self.typical_frequency_count) != 2:
            raise ValueError("typical_frequency_count must contain min/max")
        nmin, nmax = (int(v) for v in self.typical_frequency_count)
        if nmin <= 0 or nmax < nmin:
            raise ValueError("typical AirMt frequency count is invalid")
        self.typical_frequency_count = (nmin, nmax)

        rate = float(self.time_series_sampling_rate_hz)
        if not np.isfinite(rate) or rate <= 0.0:
            raise ValueError("AirMt sampling rate must be positive")
        self.time_series_sampling_rate_hz = rate

        self.input_channels = tuple(
            str(v).strip() for v in self.input_channels
        )
        self.output_channels = tuple(
            str(v).strip() for v in self.output_channels
        )
        self.reference_channels = tuple(
            str(v).strip() for v in self.reference_channels
        )
        if self.input_channels != AFMAG_TENSOR_INPUT_CHANNELS:
            raise ValueError("AirMt transfer inputs must be Hx/Hy")
        if self.output_channels != AFMAG_TENSOR_OUTPUT_CHANNELS:
            raise ValueError("AirMt airborne outputs must be Hx/Hy/Hz")
        if self.reference_channels != AFMAG_TENSOR_REFERENCE_CHANNELS:
            raise ValueError(
                "AirMt reference sensor channels must be Hx/Hy/Hz"
            )
        self.attrs = dict(self.attrs or {})

    def practical_frequency_mask(self, frequency: Any) -> np.ndarray:
        """Return a diagnostic mask for the descriptive frequency band."""
        freq = np.asarray(frequency, dtype=float)
        low, high = self.practical_frequency_range_hz
        return np.isfinite(freq) & (freq >= low) & (freq <= high)

    def to_instrument_meta(
        self,
        *,
        serial: str | None = None,
        software_version: str = "",
    ) -> InstrumentMeta:
        """Return reusable instrument metadata for tensor AFMAG/AirMt."""
        return InstrumentMeta(
            system="AirMt / tensor AFMAG",
            serial=serial,
            magnetic_sensor=SensorSpec(
                sensor_type="induction_coil",
                frequency_range=self.practical_frequency_range_hz,
                notes=(
                    "Three-component airborne magnetic receiver; transfer "
                    "functions reference horizontal ground Hx/Hy fields."
                ),
            ),
            software_version=str(software_version),
            notes=(
                "Published tensor-AFMAG/AirMt characteristics; actual "
                "delivery metadata should override descriptive defaults."
            ),
        )


@dataclass(repr=False)
class AFMAGReferenceStation(CoreObject):
    """Fixed magnetic reference station for tensor AFMAG/AirMt processing."""

    station_id: str | None = None
    site: SiteMeta | None = None
    measured_channels: tuple[str, ...] = AFMAG_TENSOR_REFERENCE_CHANNELS
    transfer_input_channels: tuple[str, ...] = AFMAG_TENSOR_INPUT_CHANNELS
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if self.station_id is not None:
            text = str(self.station_id).strip()
            self.station_id = text or None
        if self.site is not None and not isinstance(self.site, SiteMeta):
            raise TypeError("site must be a SiteMeta or None")
        measured = tuple(str(v).strip() for v in self.measured_channels)
        inputs = tuple(str(v).strip() for v in self.transfer_input_channels)
        if measured != AFMAG_TENSOR_REFERENCE_CHANNELS:
            raise ValueError("AirMt reference station must measure Hx/Hy/Hz")
        if inputs != AFMAG_TENSOR_INPUT_CHANNELS:
            raise ValueError("AirMt transfer inputs must be Hx/Hy")
        self.measured_channels = measured
        self.transfer_input_channels = inputs
        self.attrs = dict(self.attrs or {})

    @property
    def preferred_id(self) -> str | None:
        """Return explicit station ID, then the SiteMeta preferred name."""
        if self.station_id:
            return self.station_id
        if self.site is not None:
            return self.site.preferred_name
        return None
