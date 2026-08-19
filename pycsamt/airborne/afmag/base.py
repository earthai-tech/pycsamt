# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""AFMAG-family metadata built on the common airborne model.

Both generations' ``SystemSpec``/``ReferenceStation`` classes inherit
:class:`~pycsamt.core.base.CoreObject`, matching every sibling
technology's metadata
(:class:`~pycsamt.airborne.ztem.ZTEMSystemSpec`,
:class:`~pycsamt.airborne.mobilemt.MobileMTSystemSpec`, ...): they are
mutable-until-validated descriptive containers, not frozen registry
value objects, so :class:`~pycsamt.api.property.PyCSAMTObject` alone
would be the wrong base (see :mod:`pycsamt.airborne.registry` for
where that choice *is* the right one). Range/positivity/fixed-channel
normalization is delegated to :mod:`pycsamt.airborne.validation`
rather than reimplemented here field by field.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ...core.base import CoreObject
from ...metadata import InstrumentMeta, SensorSpec, SiteMeta
from ..validation import (
    normalize_count_range,
    normalize_fixed_channels,
    normalize_frequency,
    normalize_frequency_range,
    normalize_optional_identifier,
    normalize_positive_float,
)
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

    These values describe the original airborne implementation and are
    not parser constraints. Historical systems reported a comparator
    deflection proportional to the polarization-plane tilt rather than
    a digital tensor; see
    :class:`~pycsamt.airborne.afmag.adapter.AFMAGValidationError` and
    :func:`~pycsamt.airborne.afmag.adapter.build_original_afmag_emtf`
    for how that scalar response is built.

    Parameters
    ----------
    historical_frequency_band_hz : (float, float), optional
        Published historical operating band in Hz, ``(low, high)``
        with ``0 < low < high``.
    typical_frequencies_hz : tuple of float, optional
        Typical discrete operating frequencies in Hz historically used
        by comparator instruments (non-empty, finite, positive).
    coil_count : int, default 2
        Fixed at 2: the original comparator design used two crossed
        coils.
    coil_tilt_deg : float, optional
        Published nominal coil tilt angle in degrees.
    coil_separation_deg : float, optional
        Published nominal angular separation between the two coils in
        degrees.
    digital_recording : bool, default False
        Whether the archival record was digitized (``True``) or is a
        purely analogue comparator deflection (``False``).
    attrs : dict, optional
        Free-form extension metadata.

    Raises
    ------
    ValueError
        If *historical_frequency_band_hz* is not finite and correctly
        ordered, if *typical_frequencies_hz* is empty or not all
        finite/positive, if *coil_count* is not ``2``, or if
        *coil_tilt_deg*/*coil_separation_deg* is not finite/positive.
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
        """Normalize and range-check every descriptive field in place."""
        self.historical_frequency_band_hz = normalize_frequency_range(
            self.historical_frequency_band_hz,
            name="historical_frequency_band_hz",
        )
        self.typical_frequencies_hz = tuple(
            normalize_frequency(
                self.typical_frequencies_hz,
                name="typical_frequencies_hz",
            ).tolist()
        )

        self.coil_count = int(self.coil_count)
        if self.coil_count != 2:
            raise ValueError("original comparator AFMAG uses two coils")
        self.coil_tilt_deg = normalize_positive_float(
            self.coil_tilt_deg, name="coil_tilt_deg"
        )
        self.coil_separation_deg = normalize_positive_float(
            self.coil_separation_deg, name="coil_separation_deg"
        )
        self.digital_recording = bool(self.digital_recording)
        self.attrs = dict(self.attrs or {})

    def to_instrument_meta(
        self,
        *,
        serial: str | None = None,
        software_version: str = "",
    ) -> InstrumentMeta:
        """Return reusable instrument metadata for historical AFMAG.

        Only a ``magnetic_sensor`` is populated -- ``electric_sensor``
        stays ``None``, since original comparator AFMAG has no
        electric channel.
        """
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
    """Published characteristics of the tensor AFMAG/AirMt generation.

    Parameters
    ----------
    practical_frequency_range_hz : (float, float), optional
        Published practical frequency band in Hz, ``(low, high)`` with
        ``0 < low < high``. See :meth:`practical_frequency_mask`.
    typical_frequency_count : (int, int), optional
        Typical minimum/maximum count of processed frequency windows,
        ``(low, high)`` with ``0 < low <= high``.
    time_series_sampling_rate_hz : float, optional
        Published raw time-series sampling rate in Hz.
    input_channels : (str, str), default ("Hx", "Hy")
        Fixed transfer-function input channels; must equal
        ``("Hx", "Hy")``.
    output_channels : (str, str, str), default ("Hx", "Hy", "Hz")
        Fixed airborne transfer-function output channels; must equal
        ``("Hx", "Hy", "Hz")``.
    reference_channels : (str, str, str), default ("Hx", "Hy", "Hz")
        Fixed channels measured at the fixed ground reference station;
        must equal ``("Hx", "Hy", "Hz")``.
    attrs : dict, optional
        Free-form extension metadata.

    Raises
    ------
    ValueError
        If *practical_frequency_range_hz* is not finite and correctly
        ordered, if *typical_frequency_count* is not a valid
        ``0 < low <= high`` pair, if *time_series_sampling_rate_hz* is
        not finite/positive, or if *input_channels*/*output_channels*/
        *reference_channels* differs from its fixed value.
    """

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
        """Normalize and range-check every descriptive field in place."""
        self.practical_frequency_range_hz = normalize_frequency_range(
            self.practical_frequency_range_hz,
            name="practical_frequency_range_hz",
        )
        self.typical_frequency_count = normalize_count_range(
            self.typical_frequency_count,
            name="typical_frequency_count",
        )
        self.time_series_sampling_rate_hz = normalize_positive_float(
            self.time_series_sampling_rate_hz,
            name="time_series_sampling_rate_hz",
        )
        self.input_channels = normalize_fixed_channels(
            self.input_channels,
            expected=AFMAG_TENSOR_INPUT_CHANNELS,
            name="input_channels",
        )
        self.output_channels = normalize_fixed_channels(
            self.output_channels,
            expected=AFMAG_TENSOR_OUTPUT_CHANNELS,
            name="output_channels",
        )
        self.reference_channels = normalize_fixed_channels(
            self.reference_channels,
            expected=AFMAG_TENSOR_REFERENCE_CHANNELS,
            name="reference_channels",
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
        """Return reusable instrument metadata for tensor AFMAG/AirMt.

        Only a ``magnetic_sensor`` is populated -- ``electric_sensor``
        stays ``None``, since AirMt's interstation transfer function
        has no electric channel.
        """
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
    """Fixed magnetic reference station for tensor AFMAG/AirMt processing.

    Passed to
    :func:`~pycsamt.airborne.afmag.adapter.build_airmt_emtf` to
    populate :attr:`~pycsamt.emtf.EMTF.processing`'s remote-reference
    metadata; see
    :func:`~pycsamt.airborne.afmag.adapter._processing_for_reference`.

    Parameters
    ----------
    station_id : str, optional
        Explicit reference-station identifier. Falls back to
        ``site.preferred_name`` through :attr:`preferred_id` when
        omitted; see
        :func:`~pycsamt.airborne.validation.normalize_optional_identifier`.
    site : SiteMeta, optional
        Reference-station location/identity metadata.
    measured_channels : (str, str, str), default ("Hx", "Hy", "Hz")
        Fixed channels physically measured at the reference station;
        must equal ``("Hx", "Hy", "Hz")``.
    transfer_input_channels : (str, str), default ("Hx", "Hy")
        Fixed subset of *measured_channels* used as the AirMt transfer
        function's input; must equal ``("Hx", "Hy")``.
    attrs : dict, optional
        Free-form extension metadata.

    Raises
    ------
    TypeError
        If *site* is supplied and is not a
        :class:`~pycsamt.metadata.SiteMeta`.
    ValueError
        If *measured_channels* differs from ``("Hx", "Hy", "Hz")`` or
        *transfer_input_channels* differs from ``("Hx", "Hy")``.
    """

    station_id: str | None = None
    site: SiteMeta | None = None
    measured_channels: tuple[str, ...] = AFMAG_TENSOR_REFERENCE_CHANNELS
    transfer_input_channels: tuple[str, ...] = AFMAG_TENSOR_INPUT_CHANNELS
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Normalize the identifier/channels and type-check ``site``."""
        self.station_id = normalize_optional_identifier(self.station_id)
        if self.site is not None and not isinstance(self.site, SiteMeta):
            raise TypeError("site must be a SiteMeta or None")
        self.measured_channels = normalize_fixed_channels(
            self.measured_channels,
            expected=AFMAG_TENSOR_REFERENCE_CHANNELS,
            name="measured_channels",
        )
        self.transfer_input_channels = normalize_fixed_channels(
            self.transfer_input_channels,
            expected=AFMAG_TENSOR_INPUT_CHANNELS,
            name="transfer_input_channels",
        )
        self.attrs = dict(self.attrs or {})

    @property
    def preferred_id(self) -> str | None:
        """Return explicit station ID, then the SiteMeta preferred name."""
        if self.station_id:
            return self.station_id
        if self.site is not None:
            return self.site.preferred_name
        return None
