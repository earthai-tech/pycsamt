# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the public map package API."""

from __future__ import annotations

from pathlib import Path

import pycsamt.map as pcmap
from pycsamt.map._core import MapData, StationRecord
from pycsamt.map.view import MapView


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
    profile = pcmap.ProfileMap(data).with_quantity("phase").with_component("xy")
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


def test_open_app_delegates_to_launch_app(monkeypatch) -> None:
    import pycsamt.map._app as app_api

    calls = []
    monkeypatch.setattr(
        app_api, "launch_app", lambda source, **kw: calls.append((source, kw))
    )
    pcmap.open_app("some-source", port=1234)

    assert calls == [("some-source", {"port": 1234})]


def test_as_map_view_passes_mapdata_through_directly() -> None:
    import pycsamt.map._app as app_api

    data = _map_data()
    view = app_api._as_map_view(
        data,
        theme="dark",
        backend="plotly",
        detect="folder",
        recursive=True,
    )

    assert isinstance(view, MapView)
    assert view.data is data
    assert view.theme == "dark"


def test_as_map_view_str_and_path_use_from_folder(monkeypatch, tmp_path) -> None:
    import pycsamt.map._app as app_api

    calls = []
    sentinel = object()

    def fake_from_folder(cls_or_path, path=None, **kw):
        # support classmethod call signature (cls, path, **kw)
        calls.append((path, kw))
        return sentinel

    monkeypatch.setattr(MapView, "from_folder", classmethod(fake_from_folder))

    result = app_api._as_map_view(
        str(tmp_path),
        theme="light",
        backend="plotly",
        detect="auto",
        recursive=False,
        line_map={"L1": ["S00"]},
    )
    assert result is sentinel

    result2 = app_api._as_map_view(
        tmp_path,
        theme="light",
        backend="plotly",
        detect="flat",
        recursive=True,
    )
    assert result2 is sentinel

    assert calls[0] == (
        str(tmp_path),
        {
            "detect": "auto",
            "recursive": False,
            "theme": "light",
            "backend": "plotly",
            "line_map": {"L1": ["S00"]},
        },
    )
    assert calls[1] == (
        tmp_path,
        {
            "detect": "flat",
            "recursive": True,
            "theme": "light",
            "backend": "plotly",
        },
    )


def test_as_map_view_fallback_uses_load_lines(monkeypatch) -> None:
    import pycsamt.map._app as app_api

    data = _map_data()
    calls = []

    def fake_load_lines(source, **kw):
        calls.append((source, kw))
        return data

    monkeypatch.setattr(app_api, "load_lines", fake_load_lines)

    source = {"L1": ["S00", "S01"]}
    view = app_api._as_map_view(
        source,
        theme="dark",
        backend="mpl",
        detect="folder",
        recursive=True,
        verbose=2,
    )

    assert isinstance(view, MapView)
    assert view.data is data
    assert view.backend == "mpl"
    assert calls == [(source, {"verbose": 2})]


def test_launch_empty_forwards_to_mapview_launch(monkeypatch) -> None:
    import pycsamt.map._app as app_api

    calls = []
    monkeypatch.setattr(
        "pycsamt.app.mapview.launch", lambda **kw: calls.append(kw)
    )

    app_api._launch_empty(
        host="127.0.0.1", port=1234, debug=True, open_browser=False
    )

    assert calls == [
        {
            "host": "127.0.0.1",
            "port": 1234,
            "debug": True,
            "open_browser": False,
        }
    ]
