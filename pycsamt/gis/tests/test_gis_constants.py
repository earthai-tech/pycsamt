# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.gis.constants`.

Covers the UTM latitude-band letter designator (all branches) and
sanity checks on the static lookup tables.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.gis.constants import (
    _ECC_SQUARED_IDX,
    _EQUATORIAL_RADIUS_IDX,
    DEG2RAD,
    ELLIPSOIDS,
    EPSG_PROJ4,
    RAD2DEG,
    UTM_ZONE_DESIGNATOR,
    utm_letter_designator,
)


@pytest.mark.parametrize(
    "lat, letter",
    [
        (84, "X"),
        (72, "X"),
        (80, "X"),
        (70, "W"),
        (60, "V"),
        (50, "U"),
        (44, "T"),
        (36, "S"),
        (28, "R"),
        (20, "Q"),
        (12, "P"),
        (4, "N"),
        (-4, "M"),
        (-12, "L"),
        (-20, "K"),
        (-28, "J"),
        (-36, "H"),
        (-44, "G"),
        (-52, "F"),
        (-60, "E"),
        (-68, "D"),
        (-76, "C"),
        (-80, "C"),
        (90, "Z"),
        (-81, "Z"),
        (-90, "Z"),
    ],
)
def test_utm_letter_designator_all_bands(lat, letter):
    assert utm_letter_designator(lat) == letter


def test_utm_letter_designator_boundaries_are_half_open():
    # Each band is [lo, hi); the lower edge belongs to the band itself,
    # the value just below it belongs to the band below.
    assert utm_letter_designator(64) == "W"
    assert utm_letter_designator(63.999) == "V"
    assert utm_letter_designator(0) == "N"
    assert utm_letter_designator(-0.001) == "M"
    assert utm_letter_designator(-8) == "M"
    assert utm_letter_designator(-8.001) == "L"


def test_utm_zone_designator_table_shape():
    assert set(UTM_ZONE_DESIGNATOR["Z"]) == {-80, 84}
    for letter, bounds in UTM_ZONE_DESIGNATOR.items():
        assert len(bounds) == 2
        assert bounds[0] < bounds[1] or letter == "Z"


def test_deg_rad_conversion_constants():
    assert DEG2RAD == pytest.approx(np.pi / 180.0)
    assert RAD2DEG == pytest.approx(180.0 / np.pi)
    assert DEG2RAD * RAD2DEG == pytest.approx(1.0)


def test_ellipsoids_table_indices_and_wgs84_entry():
    wgs84 = [row for row in ELLIPSOIDS if row[1] == "WGS-84"]
    assert len(wgs84) == 1
    row = wgs84[0]
    assert row[0] == 23
    assert row[_EQUATORIAL_RADIUS_IDX] == pytest.approx(6378137.000)
    assert row[_ECC_SQUARED_IDX] == pytest.approx(0.006694380)


def test_epsg_proj4_table_has_zone_and_proj_string():
    assert 32611 in EPSG_PROJ4
    proj_str, zone = EPSG_PROJ4[32611]
    assert "+zone=11" in proj_str
    assert zone == 11
    # southern-hemisphere entries carry a "+south" flag
    proj_str_south, zone_south = EPSG_PROJ4[28356]
    assert "+south" in proj_str_south
    assert zone_south == 56
