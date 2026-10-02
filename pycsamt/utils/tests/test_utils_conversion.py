# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unit tests for :mod:`pycsamt.utils.conversion`."""

from __future__ import annotations

import pytest

from pycsamt.utils.conversion import (
    convert,
    convert_temperature,
    convert_time,
    convert_value,
)


# ----------------------------- convert_value ----------------------------


def test_convert_value_with_prefix():
    assert convert_value("20mm") == pytest.approx(0.02)
    assert convert_value("3.5km", target_unit="m") == pytest.approx(3500.0)
    assert convert_value("1.2kg", target_unit="g") == pytest.approx(1200.0)


def test_convert_value_no_unit_suffix_returns_raw_number():
    assert convert_value("100") == 100.0
    assert convert_value(5, target_unit="m") == 5.0


def test_convert_value_rejects_unparseable_and_incompatible():
    with pytest.raises(ValueError, match="Cannot parse"):
        convert_value("not-a-value")
    with pytest.raises(ValueError, match="Incompatible units"):
        convert_value("10m", target_unit="g")


# --------------------------- convert_temperature -------------------------


def test_convert_temperature_all_targets():
    assert convert_temperature("0C", "C", "K") == pytest.approx(273.15)
    assert convert_temperature(100, "C", "F") == pytest.approx(212.0)
    assert convert_temperature(0, "C", "C") == pytest.approx(0.0)
    assert convert_temperature(273.15, "K", "C") == pytest.approx(0.0)
    assert convert_temperature(32, "F", "C") == pytest.approx(0.0)


def test_convert_temperature_rejects_bad_input():
    with pytest.raises(ValueError, match="Cannot parse temperature"):
        convert_temperature("bogus")
    with pytest.raises(ValueError, match="Unsupported unit_to"):
        convert_temperature(0, "C", "X")


# ------------------------------ convert_time -----------------------------


def test_convert_time_every_source_and_target_unit():
    assert convert_time("3600s", "s", "h") == pytest.approx(1.0)
    assert convert_time(2, "h", "min") == pytest.approx(120.0)
    assert convert_time(1, "min", "s") == pytest.approx(60.0)
    assert convert_time(1, "d", "h") == pytest.approx(24.0)
    assert convert_time(3600, "s", "min") == pytest.approx(60.0)
    assert convert_time(3600, "s", "d") == pytest.approx(3600 / 86400)


def test_convert_time_rejects_bad_input():
    with pytest.raises(ValueError, match="Cannot parse time"):
        convert_time("bogus")
    with pytest.raises(ValueError, match="Unsupported unit_to"):
        convert_time(1, "s", "x")


# -------------------------------- convert --------------------------------


def test_convert_dispatches_by_category_and_unit_hints():
    assert convert(1, "km", "m") == pytest.approx(1000.0)
    assert convert("20mm", "mm", "m") == pytest.approx(0.02)
    assert convert(100, "C", "F", category="temperature") == pytest.approx(212.0)
    assert convert("3600s", "s", "h") == pytest.approx(1.0)


def test_convert_rounds_result_when_requested():
    assert convert(1, "km", "m", round_result=1) == 1000.0
    assert convert(1, "km", "m", round_result=0) == 1000.0
