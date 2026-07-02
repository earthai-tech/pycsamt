# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the public map package API."""

from __future__ import annotations

import pycsamt.map as pcmap
from pycsamt.map._core import MapData, StationRecord


class _Z:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        import numpy as np

        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _Edi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _Z()


class _Sites:
    def as_list(self):
        return [_Edi("S00"), _Edi("S01")]


def _map_data() -> MapData:
    return MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
        ),
    )


def test_public_builders_are_exported() -> None:
    assert pcmap.build_station_map(_map_data()).data
    assert pcmap.build_profile_map(
        _map_data(),
        pcmap.ProfileMapOptions(components=("xy",)),
    ).data
    assert pcmap.build_pseudosection(
        _map_data(),
        pcmap.ProfileMapOptions(components=("xy",)),
    ).data
    assert pcmap.build_3d_map(
        _map_data(),
        pcmap.VolumeMapOptions(),
    ).data


def test_builder_chaining_public_api() -> None:
    data = _map_data()
    station = (
        pcmap.StationMap(data)
        .with_overlay("rho", frequency=10.0)
        .with_options(show_labels=False)
    )
    profile = (
        pcmap.ProfileMap(data)
        .with_quantity("phase")
        .with_component("xy")
    )
    volume = (
        pcmap.VolumeMap(data)
        .with_mode("surface")
        .with_quantity("phase")
        .with_component("xy")
    )
    assert station.figure().data
    assert profile.figure().data
    assert volume.figure().data


def test_plot_volume_map_alias() -> None:
    fig = pcmap.plot_volume_map(
        _map_data(),
        options=pcmap.VolumeMapOptions(mode="block"),
    )
    assert fig.data


def test_mapview_app_launchers_are_public(monkeypatch) -> None:
    import pycsamt.map._app as app_api

    calls = []

    def fake_launch(**kwargs):
        calls.append(("empty", kwargs))

    monkeypatch.setattr(app_api, "_launch_empty", fake_launch)
    pcmap.launch_mapview(open_browser=False, port=9001)

    assert pcmap.launch_app is app_api.launch_app
    assert pcmap.open_app is app_api.open_app
    assert calls == [
        (
            "empty",
            {
                "host": "127.0.0.1",
                "port": 9001,
                "debug": False,
                "open_browser": False,
            },
        )
    ]


def test_mapview_app_launcher_accepts_mapview(monkeypatch) -> None:
    calls = []
    view = pcmap.MapView(_map_data())

    def fake_launch(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr(view, "launch", fake_launch)
    pcmap.launch_app(view, open_browser=False, host="0.0.0.0")

    assert calls == [
        {
            "host": "0.0.0.0",
            "port": 8770,
            "debug": False,
            "open_browser": False,
        }
    ]
