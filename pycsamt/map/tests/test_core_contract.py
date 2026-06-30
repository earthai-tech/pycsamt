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
