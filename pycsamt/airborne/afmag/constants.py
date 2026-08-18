# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Stable scientific constants for AFMAG-family adapter contracts.

The AFMAG family spans substantially different generations.  Original AFMAG
used an airborne magnetic-field comparator and reported a line-direction tilt
response.  Later tensor AFMAG/AirMt systems estimated interstation magnetic
transfer functions and rotation-invariant derived parameters.  These are kept
separate instead of being forced into one response schema.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Original comparator AFMAG
# ---------------------------------------------------------------------------
AFMAG_TILT_CODE = "AFTILT"
AFMAG_TILT_TAG = "afmag_tilt"
AFMAG_TILT_DEFAULT_UNITS = "deg"
AFMAG_ORIGINAL_TYPICAL_FREQUENCIES_HZ = (150.0, 510.0)
AFMAG_ORIGINAL_HISTORICAL_BAND_HZ = (1.0, 20_000.0)
AFMAG_ORIGINAL_COIL_COUNT = 2
AFMAG_ORIGINAL_COIL_TILT_DEG = 45.0
AFMAG_ORIGINAL_COIL_SEPARATION_DEG = 45.0

# ---------------------------------------------------------------------------
# Modern tensor AFMAG / AirMt
# ---------------------------------------------------------------------------
AFMAG_TENSOR_CODE = "TI"
AFMAG_TENSOR_TAG = "interstation_transfer_functions"
AFMAG_TENSOR_INPUT_CHANNELS = ("Hx", "Hy")
AFMAG_TENSOR_OUTPUT_CHANNELS = ("Hx", "Hy", "Hz")
AFMAG_TENSOR_UNITS = "[]"
AFMAG_TENSOR_COMPONENTS = (
    "Hxx",
    "Hxy",
    "Hyx",
    "Hyy",
    "Hzx",
    "Hzy",
)

AFMAG_AP_CODE = "AP"
AFMAG_AP_TAG = "airmt_amplification_parameter"
AFMAG_AP_UNITS = "[]"

# Published AirMt-style processing characteristics.  These are descriptive
# metadata only and must never become native-file parser limits.
AFMAG_TENSOR_PRACTICAL_FREQUENCY_RANGE_HZ = (20.0, 800.0)
AFMAG_TENSOR_TYPICAL_FREQUENCY_COUNT = (5, 6)
AFMAG_TENSOR_NOMINAL_SAMPLING_RATE_HZ = 2_000.0
AFMAG_TENSOR_REFERENCE_CHANNELS = ("Hx", "Hy", "Hz")

__all__ = [
    "AFMAG_TILT_CODE",
    "AFMAG_TILT_TAG",
    "AFMAG_TILT_DEFAULT_UNITS",
    "AFMAG_ORIGINAL_TYPICAL_FREQUENCIES_HZ",
    "AFMAG_ORIGINAL_HISTORICAL_BAND_HZ",
    "AFMAG_ORIGINAL_COIL_COUNT",
    "AFMAG_ORIGINAL_COIL_TILT_DEG",
    "AFMAG_ORIGINAL_COIL_SEPARATION_DEG",
    "AFMAG_TENSOR_CODE",
    "AFMAG_TENSOR_TAG",
    "AFMAG_TENSOR_INPUT_CHANNELS",
    "AFMAG_TENSOR_OUTPUT_CHANNELS",
    "AFMAG_TENSOR_UNITS",
    "AFMAG_TENSOR_COMPONENTS",
    "AFMAG_AP_CODE",
    "AFMAG_AP_TAG",
    "AFMAG_AP_UNITS",
    "AFMAG_TENSOR_PRACTICAL_FREQUENCY_RANGE_HZ",
    "AFMAG_TENSOR_TYPICAL_FREQUENCY_COUNT",
    "AFMAG_TENSOR_NOMINAL_SAMPLING_RATE_HZ",
    "AFMAG_TENSOR_REFERENCE_CHANNELS",
]
