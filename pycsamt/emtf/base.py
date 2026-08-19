# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Shared constants for the format-neutral EMTF scientific core."""

from __future__ import annotations

__all__ = [
    "IMPEDANCE_INPUT_CHANNELS",
    "IMPEDANCE_OUTPUT_CHANNELS",
    "TIPPER_INPUT_CHANNELS",
    "TIPPER_OUTPUT_CHANNELS",
    "LEGACY_STANDARD_ERROR_KIND",
    "LEGACY_STANDARD_ERROR_NAME",
]

IMPEDANCE_INPUT_CHANNELS = ("Hx", "Hy")
IMPEDANCE_OUTPUT_CHANNELS = ("Ex", "Ey")
TIPPER_INPUT_CHANNELS = ("Hx", "Hy")
TIPPER_OUTPUT_CHANNELS = ("Hz",)

# TFBundle's historical z_err/tipper_err fields are standard-error arrays,
# not EMTF VAR matrices. Keep that semantic distinction explicit.
LEGACY_STANDARD_ERROR_KIND = "standard_error"
LEGACY_STANDARD_ERROR_NAME = "STD"
