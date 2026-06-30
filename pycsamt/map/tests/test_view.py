# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
# ruff: noqa: E501
"""Tests for the multi-line loader and the MapView façade."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import pycsamt.map.view as view_mod
from pycsamt.map import (
    MapView,
    load_lines,
    resolve_line_groups,
)
from pycsamt.map._core import MapData, StationRecord

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_THREE_EDIS = _PROJECT_ROOT / "data" / "3edis"
_WILLY = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA"


# ── fake survey (no file IO) ───────────────────────────


class _Z:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _Edi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _Z()


def _map_data() -> MapData:
    return MapData(
        sites=[_Edi("S00"), _Edi("S01"), _Edi("S02")],
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
            StationRecord("S02", 1.2, 2.2, 30.0, "L2", 2),
        ),
    )


# ── grouping (pure path logic) ─────────────────────────


def test_resolve_groups_by_subfolder() -> None:
    groups = resolve_line_groups(_WILLY, detect="folder")
    assert "L18PLT" in groups
    assert "L22PLT" in groups
    assert groups["L18PLT"]


def test_resolve_groups_flat_folder_is_one_line() -> None:
    groups = resolve_line_groups(_THREE_EDIS, detect="folder")
    assert list(groups) == [_THREE_EDIS.name]
    assert len(groups[_THREE_EDIS.name]) == 3


def test_resolve_groups_auto_prefix() -> None:
    groups = resolve_line_groups(
        _WILLY / "L18PLT",
        detect="auto",
    )
    assert "L18" in groups
    assert groups["L18"]


def test_resolve_groups_mapping_passthrough() -> None:
    groups = resolve_line_groups({"A": ["x.edi"], "B": ["y.edi"]})
    assert groups == {"A": ["x.edi"], "B": ["y.edi"]}


def test_resolve_groups_list_of_dirs() -> None:
    groups = resolve_line_groups([
        _WILLY / "L18PLT",
        _WILLY / "L22PLT",
    ])
    assert set(groups) == {"L18PLT", "L22PLT"}


def test_resolve_groups_missing_path_raises() -> None:
    with pytest.raises(FileNotFoundError):
        resolve_line_groups("does/not/exist")


def test_resolve_groups_unknown_detect_raises() -> None:
    with pytest.raises(ValueError, match="detect"):
        resolve_line_groups(_THREE_EDIS, detect="bogus")


# ── load_lines ─────────────────────────────────────────


def test_load_lines_passthrough_map_data() -> None:
    data = _map_data()
    assert load_lines(data) is data


# ── MapView façade (no file IO) ────────────────────────


def test_mapview_from_map_data_introspection() -> None:
    mv = MapView(_map_data())
    assert mv.lines == ("L1", "L2")
    assert mv.n_stations == 3
    assert mv.has_geo is True
    assert set(mv.stations) == {"S00", "S01", "S02"}
    assert "MapView" in repr(mv)


def test_mapview_renders_all_views() -> None:
    mv = MapView(_map_data())
    assert mv.station(overlay="rho", frequency=10.0).data
    assert mv.pseudosection(component="xy").data
    assert mv.map3d(mode="fence").data
    assert mv.figure("station").data


def test_mapview_table_columns() -> None:
    df = MapView(_map_data()).table()
    assert {"ID", "Latitude", "Longitude", "Line"} <= set(
        df.columns
    )
    assert len(df) == 3


def test_mapview_figure_rejects_unknown_view() -> None:
    with pytest.raises(ValueError, match="Unknown view"):
        MapView(_map_data()).figure("nope")


def test_mapview_theme_flows_into_options() -> None:
    mv = MapView(_map_data(), theme="dark")
    # publication theme override wins over the session default
    fig = mv.station(theme="publication")
    assert fig.data


def test_mapview_export_html(monkeypatch) -> None:
    written = {}

    def fake_export(fig, options):
        written["fig"] = fig
        written["path"] = Path(options.path)
        return written["path"]

    monkeypatch.setattr(view_mod, "export_figure", fake_export)
    mv = MapView(_map_data())
    out = mv.export(
        "station.html",
        view="station",
        overlay="index",
    )
    assert written["fig"].data
    assert out.suffix == ".html"


# ── integration with bundled EDIs (skip if absent) ─────


@pytest.mark.skipif(
    not _THREE_EDIS.exists(),
    reason="bundled data/3edis not available",
)
def test_load_lines_real_folder_flat() -> None:
    mv = MapView.from_folder(_THREE_EDIS, detect="flat")
    assert mv.n_stations == 3
    assert mv.lines == (_THREE_EDIS.name,)
    assert mv.pseudosection(component="xy") is not None
