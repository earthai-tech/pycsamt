# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Robustness tests for map edge cases."""

from __future__ import annotations

import builtins

import numpy as np
import pytest

from pycsamt.map import (
    ProfileMapOptions,
    StationMapOptions,
    VolumeMapOptions,
    _backends,
    build_3d_map,
    build_pseudosection,
    build_station_map,
)
from pycsamt.map._core import (
    MapData,
    StationRecord,
    pseudosection_table,
    skin_depth_at_frequency,
    value_at_frequency_details,
)


class _Z:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _ShortZ:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.ones((1, 2, 2)) * 50.0
        self.phase = np.ones((1, 2, 2)) * 20.0


class _Edi:
    def __init__(self, station: str, z_obj=None) -> None:
        self.station = station
        if z_obj is not None:
            self.Z = z_obj


class _Sites:
    def __init__(self, *edis) -> None:
        self._edis = list(edis)

    def as_list(self):
        return list(self._edis)


def test_empty_map_data_returns_valid_empty_figures() -> None:
    data = MapData(sites=None)
    station = build_station_map(data, StationMapOptions())
    assert station.data == ()
    assert build_pseudosection(
        data,
        ProfileMapOptions(components=("xy",)),
    ).data == ()
    assert build_3d_map(data, VolumeMapOptions()).data == ()


def test_missing_z_is_skipped_without_breaking_maps() -> None:
    data = MapData(
        sites=_Sites(_Edi("S00"), _Edi("S01", _Z())),
        stations=(
            StationRecord("S00", line="L1", index=0),
            StationRecord("S01", line="L1", index=1),
        ),
    )
    details = value_at_frequency_details(data, frequency=10.0)
    assert tuple(details) == ("S01",)
    table = pseudosection_table(data)
    assert table["station"].unique().tolist() == ["S01"]
    assert build_pseudosection(
        data,
        ProfileMapOptions(components=("xy",)),
    ).data


def test_no_geo_station_map_uses_profile_fallback() -> None:
    data = MapData(
        sites=_Sites(_Edi("S00", _Z()), _Edi("S01", _Z())),
        stations=(
            StationRecord("S00", index=0),
            StationRecord("S01", index=1),
        ),
    )
    fig = build_station_map(
        data,
        StationMapOptions(overlay="rho", frequency=10.0),
    )
    assert fig.data[0].type == "scatter"


def test_invalid_component_reaches_public_builders() -> None:
    data = MapData(
        sites=_Sites(_Edi("S00", _Z())),
        stations=(StationRecord("S00"),),
    )
    with pytest.raises(ValueError, match="Unknown"):
        build_station_map(
            data,
            StationMapOptions(
                overlay="rho",
                component="bad",
            ),
        )


def test_skin_depth_uses_actual_selected_frequency() -> None:
    data = MapData(
        sites=_Sites(_Edi("S00", _Z())),
        stations=(StationRecord("S00"),),
    )
    depth = skin_depth_at_frequency(
        data,
        frequency=9.0,
    )["S00"]
    assert np.isclose(depth, 503.0 * np.sqrt(100.0 / 10.0))


def test_short_response_array_is_skipped() -> None:
    data = MapData(
        sites=_Sites(_Edi("S00", _ShortZ())),
        stations=(StationRecord("S00"),),
    )
    details = value_at_frequency_details(data, frequency=1.0)
    assert details == {}


def test_missing_plotly_backend_message(monkeypatch) -> None:
    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.startswith("plotly"):
            raise ImportError("no plotly")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(ImportError, match="Plotly"):
        _backends.require_plotly()


def test_missing_matplotlib_backend(monkeypatch) -> None:
    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.startswith("matplotlib"):
            raise ImportError("no matplotlib")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(ImportError, match="Matplotlib"):
        _backends.require_matplotlib_pyplot()
