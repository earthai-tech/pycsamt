# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Stable scientific constants for the ZTEM adapter layer.

The values here describe the published ZTEM measurement relationship and
representative system characteristics. They are descriptive metadata only,
not native-file assumptions or hard acquisition limits.
"""

from __future__ import annotations

ZTEM_TIPPER_CODE = "T"
ZTEM_TIPPER_TAG = "tipper"
ZTEM_INPUT_CHANNELS = ("Hx", "Hy")
ZTEM_OUTPUT_CHANNELS = ("Hz",)
ZTEM_TIPPER_COMPONENTS = ("Tzx", "Tzy")
ZTEM_TIPPER_UNITS = "[]"

# Published practical/common processing characteristics. The usable band is
# survey dependent; these values must never be treated as parser limits.
ZTEM_PRACTICAL_FREQUENCY_RANGE_HZ = (22.0, 720.0)
ZTEM_COMMON_PROCESSED_BANDS_HZ = (
    (30.0, 720.0),
    (25.0, 600.0),
)
ZTEM_TYPICAL_FREQUENCY_COUNT = (5, 6)
ZTEM_NOMINAL_TIME_SERIES_SAMPLING_RATE_HZ = 2_000.0
ZTEM_NOMINAL_OUTPUT_RATE_HZ = 2.5

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
]
