"""MobileMT scientific adapter contract.

The module represents the published MobileMT 3x2 admittance relationship and
maps decoded arrays into the format-neutral :mod:`pycsamt.emtf` and airborne
containers. It intentionally does not guess a proprietary native file schema.
"""

from .adapter import (
    MobileMTValidationError,
    build_mobilemt_emtf,
    build_mobilemt_dataset,
    build_mobilemt_line,
    build_mobilemt_record,
    validate_mobilemt_transfer_function,
)
from .base import MobileMTReferenceStation, MobileMTSystemSpec
from .constants import (
    MOBILEMT_ADMITTANCE_CODE,
    MOBILEMT_ADMITTANCE_TAG,
    MOBILEMT_APPARENT_CONDUCTIVITY_FIELD,
    MOBILEMT_APPARENT_CONDUCTIVITY_UNITS,
    MOBILEMT_INPUT_CHANNELS,
    MOBILEMT_NOMINAL_FREQUENCY_RANGE_HZ,
    MOBILEMT_NOMINAL_MAX_WINDOWS,
    MOBILEMT_NOMINAL_SAMPLING_RATE_HZ,
    MOBILEMT_OUTPUT_CHANNELS,
)
from .datatypes import (
    MOBILEMT_ADMITTANCE_DEFINITION,
    register_mobilemt_datatypes,
)

# Importing the explicit technology subpackage may register its TF vocabulary;
# importing ``pycsamt.airborne`` alone remains technology-neutral.
register_mobilemt_datatypes()

__all__ = [
    "MOBILEMT_ADMITTANCE_CODE",
    "MOBILEMT_ADMITTANCE_TAG",
    "MOBILEMT_ADMITTANCE_DEFINITION",
    "MOBILEMT_INPUT_CHANNELS",
    "MOBILEMT_OUTPUT_CHANNELS",
    "MOBILEMT_APPARENT_CONDUCTIVITY_FIELD",
    "MOBILEMT_APPARENT_CONDUCTIVITY_UNITS",
    "MOBILEMT_NOMINAL_FREQUENCY_RANGE_HZ",
    "MOBILEMT_NOMINAL_MAX_WINDOWS",
    "MOBILEMT_NOMINAL_SAMPLING_RATE_HZ",
    "MobileMTSystemSpec",
    "MobileMTReferenceStation",
    "MobileMTValidationError",
    "register_mobilemt_datatypes",
    "validate_mobilemt_transfer_function",
    "build_mobilemt_emtf",
    "build_mobilemt_record",
    "build_mobilemt_line",
    "build_mobilemt_dataset",
]
