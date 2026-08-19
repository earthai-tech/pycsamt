from __future__ import annotations

from pycsamt.emtf.base import (
    IMPEDANCE_INPUT_CHANNELS,
    IMPEDANCE_OUTPUT_CHANNELS,
    LEGACY_STANDARD_ERROR_KIND,
    LEGACY_STANDARD_ERROR_NAME,
    TIPPER_INPUT_CHANNELS,
    TIPPER_OUTPUT_CHANNELS,
)


def test_impedance_channel_conventions():
    assert IMPEDANCE_INPUT_CHANNELS == ("Hx", "Hy")
    assert IMPEDANCE_OUTPUT_CHANNELS == ("Ex", "Ey")


def test_tipper_channel_conventions():
    assert TIPPER_INPUT_CHANNELS == ("Hx", "Hy")
    assert TIPPER_OUTPUT_CHANNELS == ("Hz",)


def test_legacy_standard_error_markers_are_distinct_from_emtf_var():
    # TFBundle's historical z_err/tipper_err are standard-error arrays,
    # not EMTF VAR (variance) matrices -- the naming must not collide
    # with the "VAR" estimate code used elsewhere in pycsamt.emtf.
    assert LEGACY_STANDARD_ERROR_KIND == "standard_error"
    assert LEGACY_STANDARD_ERROR_NAME == "STD"
    assert LEGACY_STANDARD_ERROR_NAME != "VAR"
