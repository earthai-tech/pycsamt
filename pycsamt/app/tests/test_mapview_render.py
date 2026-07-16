# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unit tests for pycsamt.app.mapview._render — the MapView <-> Dash bridge."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.map._core import MapData, StationRecord
from pycsamt.map.view import MapView


class _Z:
    def __init__(self, freq=(10.0, 1.0)) -> None:
        self.freq = list(freq)
        self.resistivity = np.ones((len(self.freq), 2, 2)) * 100.0
        self.phase = np.ones((len(self.freq), 2, 2)) * 45.0


class _Edi:
    def __init__(self, station: str, freq=(10.0, 1.0)) -> None:
        self.station = station
        self.Z = _Z(freq)


def _map_data(**metadata) -> MapData:
    return MapData(
        sites=[_Edi("S00"), _Edi("S01"), _Edi("S02")],
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
            StationRecord("S02", 1.2, 2.2, 30.0, "L2", 2),
        ),
        metadata=metadata,
    )


def _view(**metadata) -> MapView:
    return MapView(_map_data(**metadata))


# ── store_from_view ────────────────────────────────────


class TestStoreFromView:
    def test_default_data_dir(self):
        from pycsamt.app.mapview._render import store_from_view

        store = store_from_view(_view())
        assert store["data_dir"] == "[browsed]"

    def test_custom_data_dir(self):
        from pycsamt.app.mapview._render import store_from_view

        store = store_from_view(_view(), data_dir="/some/folder")
        assert store["data_dir"] == "/some/folder"

    def test_station_records_and_counts(self):
        from pycsamt.app.mapview._render import store_from_view

        store = store_from_view(_view())
        assert len(store["station_records"]) == 3
        assert store["n_stations"] == 3
        assert store["n_lines"] == 2
        assert store["line_counts"] == {"L1": 2, "L2": 1}
        assert set(store["lines"]) == {"L1", "L2"}
        assert store["has_geo"] is True

    def test_frequencies_sorted_descending_and_positive(self):
        from pycsamt.app.mapview._render import store_from_view

        store = store_from_view(_view())
        freqs = store["frequencies"]
        assert freqs == sorted(freqs, reverse=True)
        assert all(f > 0 for f in freqs)


# ── _carried_metadata ──────────────────────────────────


class TestCarriedMetadata:
    def test_drops_counts(self):
        from pycsamt.app.mapview._render import _carried_metadata

        merged = _carried_metadata({"n_stations": 3, "n_profiles": 1, "a": 1})
        assert "n_stations" not in merged
        assert "n_profiles" not in merged
        assert merged["a"] == 1

    def test_later_source_wins(self):
        from pycsamt.app.mapview._render import _carried_metadata

        merged = _carried_metadata({"a": 1}, {"a": 2})
        assert merged["a"] == 2

    def test_sections_merge_across_sources(self):
        from pycsamt.app.mapview._render import _carried_metadata

        merged = _carried_metadata(
            {"sections": {"L1": "curtain1"}},
            {"sections": {"L2": "curtain2"}},
        )
        assert merged["sections"] == {"L1": "curtain1", "L2": "curtain2"}

    def test_no_sections_key_when_absent(self):
        from pycsamt.app.mapview._render import _carried_metadata

        merged = _carried_metadata({"a": 1}, None)
        assert "sections" not in merged

    def test_none_sources_ignored(self):
        from pycsamt.app.mapview._render import _carried_metadata

        assert _carried_metadata(None, None) == {}


# ── merge_views ─────────────────────────────────────────


class TestMergeViews:
    def test_new_station_wins_on_collision(self):
        from pycsamt.app.mapview._render import merge_views

        old = _view()
        new_data = MapData(
            sites=[_Edi("S00")],
            stations=(StationRecord("S00", 9.0, 9.0, 0.0, "L9", 0),),
        )
        new = MapView(new_data)
        merged = merge_views(old, new)
        by_id = {s.id: s for s in merged.data.stations}
        assert by_id["S00"].line == "L9"

    def test_order_preserved_new_stations_appended(self):
        from pycsamt.app.mapview._render import merge_views

        old = _view()
        new_data = MapData(
            sites=[_Edi("S99")],
            stations=(StationRecord("S99", 5.0, 5.0, 0.0, "L3", 0),),
        )
        new = MapView(new_data)
        merged = merge_views(old, new)
        ids = [s.id for s in merged.data.stations]
        assert ids == ["S00", "S01", "S02", "S99"]

    def test_theme_and_backend_from_old(self):
        from pycsamt.app.mapview._render import merge_views

        old = MapView(_map_data(), theme="dark", backend="plotly")
        new = MapView(_map_data())
        merged = merge_views(old, new)
        assert merged.theme == "dark"
        assert merged.backend == "plotly"

    def test_metadata_sections_carried(self):
        from pycsamt.app.mapview._render import merge_views

        old = MapView(_map_data(sections={"L1": "c1"}))
        new = MapView(_map_data(sections={"L2": "c2"}))
        merged = merge_views(old, new)
        assert merged.data.metadata["sections"] == {"L1": "c1", "L2": "c2"}


# ── exclude_stations / restrict_to_lines / apply_settings ─


class TestExcludeStations:
    def test_no_masked_returns_same_view(self):
        from pycsamt.app.mapview._render import exclude_stations

        view = _view()
        assert exclude_stations(view, None) is view
        assert exclude_stations(view, []) is view

    def test_masked_id_not_present_returns_same_view(self):
        from pycsamt.app.mapview._render import exclude_stations

        view = _view()
        assert exclude_stations(view, ["nope"]) is view

    def test_excludes_matching_station(self):
        from pycsamt.app.mapview._render import exclude_stations

        view = _view()
        result = exclude_stations(view, ["S00"])
        ids = {s.id for s in result.data.stations}
        assert ids == {"S01", "S02"}

    def test_masked_ids_stringified(self):
        from pycsamt.app.mapview._render import exclude_stations

        view = _view()
        result = exclude_stations(view, [0])  # not matching any real id
        assert {s.id for s in result.data.stations} == {"S00", "S01", "S02"}


class TestRestrictToLines:
    def test_empty_active_returns_same_view(self):
        from pycsamt.app.mapview._render import restrict_to_lines

        view = _view()
        assert restrict_to_lines(view, None) is view
        assert restrict_to_lines(view, []) is view

    def test_all_lines_active_returns_same_view(self):
        from pycsamt.app.mapview._render import restrict_to_lines

        view = _view()
        assert restrict_to_lines(view, ["L1", "L2"]) is view

    def test_restricts_to_single_line(self):
        from pycsamt.app.mapview._render import restrict_to_lines

        view = _view()
        result = restrict_to_lines(view, ["L1"])
        assert {s.id for s in result.data.stations} == {"S00", "S01"}


class TestApplySettings:
    def test_combines_line_and_mask_filters(self):
        from pycsamt.app.mapview._render import apply_settings

        view = _view()
        result = apply_settings(view, ["L1"], ["S00"])
        assert {s.id for s in result.data.stations} == {"S01"}

    def test_no_filters_returns_same_view(self):
        from pycsamt.app.mapview._render import apply_settings

        view = _view()
        assert apply_settings(view, None, None) is view


# ── reproject_view / project_to_crs / _source_epsg ─────


class TestSourceEpsg:
    def test_utm_north(self):
        from pycsamt.app.mapview._render import _source_epsg

        assert _source_epsg("utm", 32, None, "N") == 32632

    def test_utm_south(self):
        from pycsamt.app.mapview._render import _source_epsg

        assert _source_epsg("utm", 32, None, "S") == 32732

    def test_utm_invalid_zone_falls_back_to_50(self):
        from pycsamt.app.mapview._render import _source_epsg

        assert _source_epsg("utm", "abc", None, "N") == 32650

    def test_custom_epsg(self):
        from pycsamt.app.mapview._render import _source_epsg

        assert _source_epsg("custom", None, 3857, "N") == 3857

    def test_custom_invalid_epsg_falls_back_to_4326(self):
        from pycsamt.app.mapview._render import _source_epsg

        assert _source_epsg("custom", None, "not-a-number", "N") == 4326

    def test_custom_missing_epsg_falls_back_to_4326(self):
        from pycsamt.app.mapview._render import _source_epsg

        assert _source_epsg("custom", None, None, "N") == 4326


class TestReprojectView:
    def test_none_or_geo_mode_returns_same_view(self):
        from pycsamt.app.mapview._render import reproject_view

        view = _view()
        assert reproject_view(view, None, None, "N", None) is view
        assert reproject_view(view, "geo", None, "N", None) is view

    def test_code_4326_returns_same_view(self):
        from pycsamt.app.mapview._render import reproject_view

        view = _view()
        assert reproject_view(view, "custom", None, "N", 4326) is view

    def test_utm_reprojects_using_transform_xy(self, monkeypatch):
        import pycsamt.map.overlays as overlays_mod

        from pycsamt.app.mapview._render import reproject_view

        def fake_transform_xy(lon, lat, crs):
            return np.asarray(lon) + 100.0, np.asarray(lat) + 200.0

        monkeypatch.setattr(overlays_mod, "transform_xy", fake_transform_xy)

        view = _view()
        result = reproject_view(view, "utm", 32, "N", None)
        new_lon = [s.longitude for s in result.data.stations]
        old_lon = [s.longitude for s in view.data.stations]
        assert new_lon == [lo + 100.0 for lo in old_lon]

    def test_transform_failure_returns_same_view(self, monkeypatch):
        import pycsamt.map.overlays as overlays_mod

        from pycsamt.app.mapview._render import reproject_view

        def boom(*a, **k):
            raise RuntimeError("proj failure")

        monkeypatch.setattr(overlays_mod, "transform_xy", boom)

        view = _view()
        assert reproject_view(view, "utm", 32, "N", None) is view

    def test_no_finite_coords_returns_same_view(self):
        from pycsamt.app.mapview._render import reproject_view

        data = MapData(
            sites=[_Edi("S00")],
            stations=(StationRecord("S00", None, None, 0.0, "L1", 0),),
        )
        view = MapView(data)
        assert reproject_view(view, "utm", 32, "N", None) is view


class TestProjectToCrs:
    def test_geo_mode_returns_none_pair(self):
        from pycsamt.app.mapview._render import project_to_crs

        east, north, code = project_to_crs(
            [1.0], [2.0], "geo", None, "N", None
        )
        assert east is None and north is None
        assert code == 4326

    def test_utm_mode_returns_transformed_arrays(self, monkeypatch):
        import pycsamt.map.overlays as overlays_mod

        from pycsamt.app.mapview._render import project_to_crs

        def fake_transform_xy(lon, lat, crs):
            return np.asarray(lon) * 2, np.asarray(lat) * 2

        monkeypatch.setattr(overlays_mod, "transform_xy", fake_transform_xy)
        east, north, code = project_to_crs(
            [1.0, 2.0], [3.0, 4.0], "utm", 32, "N", None
        )
        assert list(east) == [2.0, 4.0]
        assert list(north) == [6.0, 8.0]
        assert code == 32632

    def test_transform_failure_returns_none_pair(self, monkeypatch):
        import pycsamt.map.overlays as overlays_mod

        from pycsamt.app.mapview._render import project_to_crs

        def boom(*a, **k):
            raise RuntimeError("boom")

        monkeypatch.setattr(overlays_mod, "transform_xy", boom)
        east, north, code = project_to_crs(
            [1.0], [2.0], "utm", 32, "N", None
        )
        assert east is None and north is None
        assert code == 32632


# ── figure_for ──────────────────────────────────────────


class TestFigureFor:
    def test_view_none_returns_empty_figure(self):
        from pycsamt.app.mapview._render import figure_for

        fig = figure_for("map", None, {})
        assert fig.data == () or fig.layout.annotations

    def test_zero_stations_returns_empty_figure(self):
        from pycsamt.app.mapview._render import figure_for

        empty_view = MapView(MapData(sites=[], stations=()))
        fig = figure_for("map", empty_view, {})
        assert fig.layout.annotations[0].text == "Load EDI lines to begin"

    def test_all_masked_returns_hidden_message(self):
        from pycsamt.app.mapview._render import figure_for

        view = _view()
        fig = figure_for(
            "map", view, {}, masked=["S00", "S01", "S02"]
        )
        assert "hidden or masked" in fig.layout.annotations[0].text

    def test_map_view_sets_uirevision(self):
        from pycsamt.app.mapview._render import figure_for

        view = _view()
        fig = figure_for("map", view, {}, fit=3)
        assert fig.layout.uirevision == "fit-3"

    def test_pseudosection_view_renders(self):
        from pycsamt.app.mapview._render import figure_for

        view = _view()
        fig = figure_for("pseudosection", view, {"component": "xy"})
        assert fig.data

    def test_map3d_view_sets_uirevision(self):
        from pycsamt.app.mapview._render import figure_for

        view = _view()
        fig = figure_for("map3d", view, {"mode3d": "fence"}, fit=2)
        assert fig.layout.uirevision == "fit-2"

    def test_unknown_view_name_returns_message(self):
        from pycsamt.app.mapview._render import figure_for

        view = _view()
        fig = figure_for("bogus", view, {})
        assert "Unknown view: bogus" in fig.layout.annotations[0].text

    def test_view_name_defaults_to_map_when_falsy(self):
        from pycsamt.app.mapview._render import figure_for

        view = _view()
        fig = figure_for("", view, {})
        assert fig.data

    def test_figure_background_is_transparent(self):
        from pycsamt.app.mapview._render import figure_for

        view = _view()
        fig = figure_for("map", view, {})
        assert fig.layout.paper_bgcolor == "rgba(0,0,0,0)"
        assert fig.layout.plot_bgcolor == "rgba(0,0,0,0)"


# ── _pair ───────────────────────────────────────────────


class TestPair:
    def test_valid_ordered_pair(self):
        from pycsamt.app.mapview._render import _pair

        assert _pair(1, 5) == (1.0, 5.0)

    def test_hi_not_greater_than_lo_returns_none(self):
        from pycsamt.app.mapview._render import _pair

        assert _pair(5, 5) is None
        assert _pair(5, 1) is None

    def test_negative_lo_returns_none(self):
        from pycsamt.app.mapview._render import _pair

        assert _pair(-1, 5) is None

    def test_non_numeric_returns_none(self):
        from pycsamt.app.mapview._render import _pair

        assert _pair("a", "b") is None
        assert _pair(None, None) is None


# ── _transparent / empty_figure ────────────────────────


class TestTransparentAndEmptyFigure:
    def test_transparent_sets_backgrounds(self):
        import plotly.graph_objects as go

        from pycsamt.app.mapview._render import _transparent

        fig = go.Figure()
        result = _transparent(fig)
        assert result.layout.paper_bgcolor == "rgba(0,0,0,0)"
        assert result.layout.plot_bgcolor == "rgba(0,0,0,0)"

    def test_transparent_clears_3d_scene(self):
        import plotly.graph_objects as go

        from pycsamt.app.mapview._render import _transparent

        fig = go.Figure(data=[go.Scatter3d(x=[0], y=[0], z=[0])])
        result = _transparent(fig)
        assert result.layout.scene.bgcolor == "rgba(0,0,0,0)"

    def test_empty_figure_default_message(self):
        from pycsamt.app.mapview._render import empty_figure

        fig = empty_figure()
        assert fig.layout.annotations[0].text == "Load EDI lines to begin"

    def test_empty_figure_custom_message_and_theme(self):
        from pycsamt.app.mapview._render import empty_figure

        fig = empty_figure(theme="dark", msg="Custom message")
        assert fig.layout.annotations[0].text == "Custom message"
        assert fig.layout.paper_bgcolor == "rgba(0,0,0,0)"
        assert fig.layout.xaxis.visible is False
        assert fig.layout.yaxis.visible is False
