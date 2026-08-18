"""ZTEM scientific adapter contract.

The module maps decoded ZTEM tipper values into the existing EMTF ``T``
transfer function and common airborne containers. It intentionally does not
invent a proprietary native file reader.
"""

from .adapter import (
    ZTEMValidationError,
    build_ztem_dataset,
    build_ztem_emtf,
    build_ztem_line,
    build_ztem_record,
    validate_ztem_transfer_function,
)
from .base import ZTEMReferenceStation, ZTEMSystemSpec
from .constants import (
    ZTEM_COMMON_PROCESSED_BANDS_HZ,
    ZTEM_INPUT_CHANNELS,
    ZTEM_NOMINAL_OUTPUT_RATE_HZ,
    ZTEM_NOMINAL_TIME_SERIES_SAMPLING_RATE_HZ,
    ZTEM_OUTPUT_CHANNELS,
    ZTEM_PRACTICAL_FREQUENCY_RANGE_HZ,
    ZTEM_TIPPER_COMPONENTS,
    ZTEM_TIPPER_CODE,
    ZTEM_TIPPER_TAG,
    ZTEM_TIPPER_UNITS,
    ZTEM_TYPICAL_FREQUENCY_COUNT,
)

__all__ = [
    "ZTEM_TIPPER_CODE",
    "ZTEM_TIPPER_TAG",
    "ZTEM_INPUT_CHANNELS",
    "ZTEM_OUTPUT_CHANNELS",
    "ZTEM_TIPPER_COMPONENTS",
    "ZTEM_TIPPER_UNITS",
    "ZTEM_PRACTICAL_FREQUENCY_RANGE_HZ",
    "ZTEM_COMMON_PROCESSED_BANDS_HZ",
    "ZTEM_TYPICAL_FREQUENCY_COUNT",
    "ZTEM_NOMINAL_TIME_SERIES_SAMPLING_RATE_HZ",
    "ZTEM_NOMINAL_OUTPUT_RATE_HZ",
    "ZTEMSystemSpec",
    "ZTEMReferenceStation",
    "ZTEMValidationError",
    "validate_ztem_transfer_function",
    "build_ztem_emtf",
    "build_ztem_record",
    "build_ztem_line",
    "build_ztem_dataset",
]
