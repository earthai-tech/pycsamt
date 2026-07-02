# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for normalized map data and component handling."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.map import (
    MapData,
    StationRecord,
    component_spec,
    ensure_map_data,
    normalize_component,
    select_frequency,
    value_at_frequency_details,
)
from pycsamt.map._core import (
    component_values,
    frequency_axis,
    pseudosection_table,
)


class _Z:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.array(
            [
                [[10.0, 100.0], [400.0, 20.0]],
                [[30.0, 900.0], [1600.0, 40.0]],
            ]
        )
        self.phase = np.array(
            [
                [[1.0, 40.0], [50.0, 2.0]],
                [[3.0, 45.0], [55.0, 4.0]],
            ]
        )


class _Edi:
    station = "S00"
    Z = _Z()


class _Header:
    station = "H00"
    lat = "5.25"
    lon = "-4.15"
    elev = "48.0"
    line = "L18"


class _NestedEdi:
    Head = _Header()
    Z = _Z()


class _DictHeaderEdi:
    header = {
        "Station": "D00",
        "LATITUDE": 6.1,
        "LONGITUDE": -3.9,
        "ELEVATION": 52.0,
        "Profile": "L22",
    }
    Z = _Z()


class _CoordEdi:
    station = "C00"
    coords = (7.0, -2.0)
    Z = _Z()


class _Location:
    lat = 32.12
    lon = 119.13
    elev = 99.0


class _HeadSection:
    Location = _Location()


class _SectionEdi:
    """EDIFile-like: coords only reachable via get_section('head')."""

    station = "18-001A"
    Z = _Z()

    def get_section(self, name):
        return _HeadSection() if str(name).lower() == "head" else None


def test_map_data_contract_without_sites() -> None:
    data = MapData(
        sites=None,
        stations=(
            StationRecord(
                id=12,
                latitude="nan",
                longitude=2.0,
                line=18,
            ),
        ),
    )
    assert data.station_ids == ("12",)
    assert data.lines == ("18",)
    assert data.has_geo is False
    assert data.iter_edis() == ()


def test_ensure_map_data_reads_nested_header_coordinates(
    monkeypatch,
) -> None:
    from pycsamt.map import _core

    monkeypatch.setattr(
        _core,
        "ensure_sites",
        lambda data, recursive=True, verbose=0: data,
    )
    data = ensure_map_data([_NestedEdi(), _DictHeaderEdi()])
    first, second = data.stations

    assert data.has_geo is True
    assert first.id == "H00"
    assert first.latitude == 5.25
    assert first.longitude == -4.15
    assert first.elevation == 48.0
    assert first.line == "L18"
    assert second.id == "D00"
    assert second.latitude == 6.1
    assert second.longitude == -3.9
    assert second.elevation == 52.0
    assert second.line == "L22"


def test_ensure_map_data_reads_get_section_head_location(
    monkeypatch,
) -> None:
    """EDIFile exposes coords via get_section('head').Location, not attrs."""
    from pycsamt.map import _core

    monkeypatch.setattr(
        _core,
        "ensure_sites",
        lambda data, recursive=True, verbose=0: data,
    )
    data = ensure_map_data([_SectionEdi()])
    station = data.stations[0]
    assert data.has_geo is True
    assert station.latitude == 32.12
    assert station.longitude == 119.13
    assert station.elevation == 99.0


def test_ensure_map_data_accepts_coordinate_tuple(
    monkeypatch,
) -> None:
    from pycsamt.map import _core

    monkeypatch.setattr(
        _core,
        "ensure_sites",
        lambda data, recursive=True, verbose=0: data,
    )
    station = ensure_map_data([_CoordEdi()]).stations[0]

    assert station.latitude == 7.0
    assert station.longitude == -2.0


def test_frequency_axis_accepts_iterable_sites() -> None:
    data = MapData(
        sites=[_Edi()],
        stations=(StationRecord("S00"),),
    )
    assert frequency_axis(data).tolist() == [10.0, 1.0]


def test_select_frequency_reports_metadata() -> None:
    selected = select_frequency(
        [100.0, 10.0, 1.0],
        9.0,
        tolerance=2.0,
    )
    assert selected is not None
    assert selected.actual == 10.0
    assert selected.index == 1
    assert selected.within_tolerance is True


def test_frequency_details_honor_tolerance() -> None:
    data = MapData(
        sites=[_Edi()],
        stations=(StationRecord("S00"),),
    )
    details = value_at_frequency_details(
        data,
        frequency=9.0,
        tolerance=2.0,
    )
    assert details["S00"].selection.actual == 10.0
    skipped = value_at_frequency_details(
        data,
        frequency=9.0,
        tolerance=0.1,
    )
    assert skipped == {}


def test_component_aliases_and_derived_values() -> None:
    assert normalize_component("ZXY") == "xy"
    assert component_spec("det").mode == "det"
    vals = component_values(
        _Edi.Z.resistivity,
        "det",
        quantity="rho",
    )
    assert np.allclose(vals, [200.0, 1200.0])


def test_invalid_component_raises() -> None:
    with pytest.raises(ValueError, match="Unknown"):
        normalize_component("bad")


def test_pseudosection_table_uses_average_component() -> None:
    data = MapData(
        sites=[_Edi()],
        stations=(StationRecord("S00"),),
    )
    table = pseudosection_table(data, component="avg")
    assert table["value"].tolist() == [250.0, 1250.0]
