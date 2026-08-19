# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""ZTEM-specific metadata that reuses the common airborne/EMTF model.

Both classes inherit :class:`~pycsamt.core.base.CoreObject`, matching
:mod:`pycsamt.airborne.base` and every sibling technology's metadata
(:class:`~pycsamt.airborne.mobilemt.MobileMTSystemSpec`,
:class:`~pycsamt.airborne.afmag.AirMtSystemSpec`, ...): they are
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
    normalize_frequency_range,
    normalize_optional_identifier,
    normalize_positive_float,
)
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

    The values are informative defaults. ZTEM frequency selection
    depends on platform speed, sampling, signal strength, and noise,
    so processed survey products may legitimately use a different
    subset or frequency grid; nothing here constrains
    :func:`~pycsamt.airborne.ztem.build_ztem_emtf`'s accepted input.

    Parameters
    ----------
    practical_frequency_range_hz : (float, float), optional
        Published practical frequency band in Hz, ``(low, high)`` with
        ``0 < low < high``. See :meth:`practical_frequency_mask`.
    common_processed_bands_hz : tuple of (float, float), optional
        Commonly published processed sub-bands in Hz; each validated
        the same way as *practical_frequency_range_hz*.
    typical_frequency_count : (int, int), optional
        Typical minimum/maximum count of processed frequency windows,
        ``(low, high)`` with ``0 < low <= high``.
    time_series_sampling_rate_hz : float, optional
        Published raw time-series sampling rate in Hz.
    nominal_output_rate_hz : float, optional
        Published processed-output rate in Hz.
    input_channels : (str, str), default ("Hx", "Hy")
        Fixed ground-reference horizontal magnetic input channels;
        must equal ``("Hx", "Hy")``.
    output_channels : (str,), default ("Hz",)
        Fixed airborne vertical magnetic output channel; must equal
        ``("Hz",)``.
    attrs : dict, optional
        Free-form extension metadata.

    Raises
    ------
    ValueError
        If any frequency range/band/count is not finite and correctly
        ordered, if any rate is not finite and positive, or if
        *input_channels*/*output_channels* differ from their fixed
        values.
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
        """Normalize and range-check every descriptive field in place."""
        self.practical_frequency_range_hz = normalize_frequency_range(
            self.practical_frequency_range_hz,
            name="practical_frequency_range_hz",
        )
        self.common_processed_bands_hz = tuple(
            normalize_frequency_range(
                band,
                name="common_processed_bands_hz",
            )
            for band in self.common_processed_bands_hz
        )
        self.typical_frequency_count = normalize_count_range(
            self.typical_frequency_count,
            name="typical_frequency_count",
        )
        self.time_series_sampling_rate_hz = normalize_positive_float(
            self.time_series_sampling_rate_hz,
            name="time_series_sampling_rate_hz",
        )
        self.nominal_output_rate_hz = normalize_positive_float(
            self.nominal_output_rate_hz,
            name="nominal_output_rate_hz",
        )
        self.input_channels = normalize_fixed_channels(
            self.input_channels,
            expected=ZTEM_INPUT_CHANNELS,
            name="input_channels",
        )
        self.output_channels = normalize_fixed_channels(
            self.output_channels,
            expected=ZTEM_OUTPUT_CHANNELS,
            name="output_channels",
        )
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
        """Return reusable instrument metadata without vendor guesses.

        Only a ``magnetic_sensor`` is populated -- ``electric_sensor``
        stays ``None`` because ZTEM's tipper geometry has no electric
        channel, unlike MobileMT's admittance (see
        :meth:`~pycsamt.airborne.mobilemt.MobileMTSystemSpec.to_instrument_meta`).
        """
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
    """Fixed ground magnetic reference station used by ZTEM processing.

    Passed to :func:`~pycsamt.airborne.ztem.build_ztem_emtf` to
    populate :attr:`~pycsamt.emtf.EMTF.processing`'s remote-reference
    metadata; see
    :func:`~pycsamt.airborne.ztem.adapter._processing_for_reference`.

    Parameters
    ----------
    station_id : str, optional
        Explicit reference-station identifier. Falls back to
        ``site.preferred_name`` through :attr:`preferred_id` when
        omitted; see
        :func:`~pycsamt.airborne.validation.normalize_optional_identifier`.
    site : SiteMeta, optional
        Reference-station location/identity metadata.
    magnetic_channels : (str, str), default ("Hx", "Hy")
        Fixed horizontal magnetic channels measured at the reference
        station; must equal ``("Hx", "Hy")``.
    attrs : dict, optional
        Free-form extension metadata.

    Raises
    ------
    TypeError
        If *site* is supplied and is not a
        :class:`~pycsamt.metadata.SiteMeta`.
    ValueError
        If *magnetic_channels* differs from ``("Hx", "Hy")``.
    """

    station_id: str | None = None
    site: SiteMeta | None = None
    magnetic_channels: tuple[str, str] = ZTEM_INPUT_CHANNELS
    attrs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Normalize the identifier/channels and type-check ``site``."""
        self.station_id = normalize_optional_identifier(self.station_id)
        if self.site is not None and not isinstance(self.site, SiteMeta):
            raise TypeError("site must be a SiteMeta or None")
        self.magnetic_channels = normalize_fixed_channels(
            self.magnetic_channels,
            expected=ZTEM_INPUT_CHANNELS,
            name="magnetic_channels",
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
