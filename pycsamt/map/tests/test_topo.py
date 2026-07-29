# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
# ruff: noqa: E501
"""Tests for elevation sources (pycsamt.map.topo)."""

from __future__ import annotations

import base64

import numpy as np

from pycsamt.map import (
    apply_elevations,
    parse_elevation_file,
)
from pycsamt.map._core import MapData, StationRecord


class _Z:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _Edi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _Z()


def _data() -> MapData:
    return MapData(
        sites=[_Edi("A0"), _Edi("A1")],
        stations=(
            StationRecord("A0", 1.0, 2.0, None, "L1", 0),
            StationRecord("A1", 1.1, 2.1, None, "L1", 1),
        ),
    )


def test_apply_elevations_overrides_by_id() -> None:
    out = apply_elevations(_data(), {"A0": 500.0, "A1": 650.0})
    elevs = [s.elevation for s in out.stations]
    assert elevs == [500.0, 650.0]


def test_apply_elevations_ignores_unknown_and_nonfinite() -> None:
    out = apply_elevations(_data(), {"A0": 500.0, "ZZ": 9.0, "A1": float("nan")})
    elevs = [s.elevation for s in out.stations]
    assert elevs == [500.0, None]


def test_apply_elevations_empty_is_noop() -> None:
    data = _data()
    assert apply_elevations(data, {}) is data


def test_parse_csv_elevation_file() -> None:
    csv = "station,elevation\nA0,500\nA1,650\n"
    b64 = "data:text/csv;base64," + base64.b64encode(csv.encode()).decode()
    elev = parse_elevation_file(b64, "topo.csv")
    assert elev == {"A0": 500.0, "A1": 650.0}


def test_parse_csv_flexible_columns() -> None:
    csv = "ID,Z\nA0,12.5\nA1,13.0\n"
    b64 = base64.b64encode(csv.encode()).decode()
    elev = parse_elevation_file(b64, "x.csv")
    assert elev == {"A0": 12.5, "A1": 13.0}


def test_parse_unknown_columns_returns_empty() -> None:
    csv = "foo,bar\n1,2\n"
    b64 = base64.b64encode(csv.encode()).decode()
    assert parse_elevation_file(b64, "x.csv") == {}


def test_parse_bad_payload_returns_empty() -> None:
    assert parse_elevation_file("not-base64!!", "x.csv") == {}
