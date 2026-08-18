# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Stable scientific constants for the MobileMT adapter layer.

The values in this module describe the published MobileMT measurement
relationship and nominal system specifications. They are deliberately not used
as hard file-format assumptions: real vendor delivery data may expose a
subset, a different processed frequency grid, or additional quality fields.
"""

from __future__ import annotations

MOBILEMT_ADMITTANCE_CODE = "Y"
MOBILEMT_ADMITTANCE_TAG = "mobilemt_admittance"
MOBILEMT_INPUT_CHANNELS = ("Ex", "Ey")
MOBILEMT_OUTPUT_CHANNELS = ("Hx", "Hy", "Hz")
MOBILEMT_APPARENT_CONDUCTIVITY_FIELD = "apparent_conductivity"
MOBILEMT_APPARENT_CONDUCTIVITY_UNITS = "mS/m"

# Published nominal specifications. They are informational defaults only and
# are not parser limits.
MOBILEMT_NOMINAL_FREQUENCY_RANGE_HZ = (19.0, 26_000.0)
MOBILEMT_NOMINAL_MAX_WINDOWS = 30
MOBILEMT_NOMINAL_SAMPLING_RATE_HZ = 73_728.0

__all__ = [
    "MOBILEMT_ADMITTANCE_CODE",
    "MOBILEMT_ADMITTANCE_TAG",
    "MOBILEMT_INPUT_CHANNELS",
    "MOBILEMT_OUTPUT_CHANNELS",
    "MOBILEMT_APPARENT_CONDUCTIVITY_FIELD",
    "MOBILEMT_APPARENT_CONDUCTIVITY_UNITS",
    "MOBILEMT_NOMINAL_FREQUENCY_RANGE_HZ",
    "MOBILEMT_NOMINAL_MAX_WINDOWS",
    "MOBILEMT_NOMINAL_SAMPLING_RATE_HZ",
]
