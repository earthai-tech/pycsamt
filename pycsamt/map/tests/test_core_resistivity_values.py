# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :func:`pycsamt.map.all_resistivity_values`."""

from __future__ import annotations

import numpy as np

from pycsamt.map import all_resistivity_values
from pycsamt.map._core import MapData, StationRecord


class _Z:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.array(
            [
                [[100.0, 100.0], [100.0, 100.0]],
                [[100.0, 10000.0], [100.0, -5.0]],  # a bad/negative cell
            ]
        )


class _Edi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _Z()


class _Sites:
    def as_list(self):
        return [_Edi("S00"), _Edi("S01")]


def test_edi_backed_data_flattens_finite_positive_resistivity():
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )
    vals = all_resistivity_values(data)
    assert vals.size > 0
    assert np.all(vals > 0)
    assert np.all(np.isfinite(vals))
    assert vals.max() == 10000.0


def test_section_backed_data_prefers_precomputed_sections():
    data = MapData(
        sites=None,
        stations=(),
        metadata={
            "sections": {
                "L1": {
                    "stations": np.array(["a", "b"]),
                    "rho": np.array([[10.0, 2000.0], [50.0, 3000.0]]),
                    "z": np.array([0.0, 100.0]),
                },
                "L2": {
                    "stations": np.array(["c"]),
                    "rho": np.array([[np.nan], [-1.0]]),
                    "z": np.array([0.0, 50.0]),
                },
            }
        },
    )
    vals = all_resistivity_values(data)
    assert sorted(vals.tolist()) == [10.0, 50.0, 2000.0, 3000.0]


def test_empty_data_returns_empty_array_not_an_error():
    data = MapData(sites=None, stations=(), metadata={})
    vals = all_resistivity_values(data)
    assert vals.size == 0
