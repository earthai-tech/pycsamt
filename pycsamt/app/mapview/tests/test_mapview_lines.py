# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.lines."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestColorsFor:
    def test_assigns_one_color_per_line(self):
        from pycsamt.app.mapview.callbacks.lines import _colors_for

        colors = _colors_for(["L1", "L2", "L3"])
        assert set(colors) == {"L1", "L2", "L3"}
        assert len(set(colors.values())) == 3

    def test_wraps_around_palette(self):
        from pycsamt.app.mapview.callbacks.lines import (
            _PALETTE,
            _colors_for,
        )

        lines = [f"L{i}" for i in range(len(_PALETTE) + 2)]
        colors = _colors_for(lines)
        assert colors["L0"] == colors[f"L{len(_PALETTE)}"]

    def test_empty_lines_returns_empty_dict(self):
        from pycsamt.app.mapview.callbacks.lines import _colors_for

        assert _colors_for([]) == {}


class TestRowMeta:
    def test_all_fields_present(self):
        from pycsamt.app.mapview.callbacks.lines import _row_meta

        rec = {"Latitude": 48.5, "Longitude": 7.75, "Elevation": 200}
        meta = _row_meta(rec)
        assert "48.5000" in meta or "+48.5000" in meta
        assert "200" in meta

    def test_missing_fields_returns_dash(self):
        from pycsamt.app.mapview.callbacks.lines import _row_meta

        assert _row_meta({}) == "—"

    def test_bad_values_are_swallowed(self):
        from pycsamt.app.mapview.callbacks.lines import _row_meta

        rec = {"Latitude": "not-a-number"}
        assert _row_meta(rec) == "—"


class TestFmtCoord:
    def test_degree_suffix_formats_five_decimals(self):
        from pycsamt.app.mapview.callbacks.lines import _fmt_coord

        assert _fmt_coord(7.75, "°") == "7.75000°"

    def test_meter_suffix_formats_with_commas(self):
        from pycsamt.app.mapview.callbacks.lines import _fmt_coord

        assert _fmt_coord(123456.789, " m") == "123,456.8 m"

    def test_non_numeric_returns_dash(self):
        from pycsamt.app.mapview.callbacks.lines import _fmt_coord

        assert _fmt_coord(None, "°") == "—"
        assert _fmt_coord("nope", "°") == "—"


class TestCoordinateFields:
    def test_geo_mode_returns_lat_lon_only(self):
        from pycsamt.app.mapview.callbacks.lines import _coordinate_fields

        rec = {"Latitude": 48.5, "Longitude": 7.75}
        rows = _coordinate_fields(
            rec, {"crs_mode": "geo"}, lambda label, value: (label, value)
        )
        assert rows == [("Latitude", "48.50000°"), ("Longitude", "7.75000°")]

    def test_missing_coords_returns_base_only(self):
        from pycsamt.app.mapview.callbacks.lines import _coordinate_fields

        rec = {"Latitude": None, "Longitude": None}
        rows = _coordinate_fields(
            rec, {"crs_mode": "utm"}, lambda label, value: (label, value)
        )
        assert len(rows) == 2

    def test_utm_mode_appends_easting_northing(self, monkeypatch):
        import pycsamt.app.mapview._render as render_mod
        import pycsamt.app.mapview.callbacks.lines as lines_mod

        def fake_project(lons, lats, mode, zone, hem, epsg):
            return [500000.0], [5370000.0], 32632

        monkeypatch.setattr(render_mod, "project_to_crs", fake_project)

        rec = {"Latitude": 48.5, "Longitude": 7.75}
        rows = lines_mod._coordinate_fields(
            rec,
            {"crs_mode": "utm", "utm_zone": 32, "utm_hem": "N"},
            lambda label, value: (label, value),
        )
        assert len(rows) == 4
        assert rows[2][0] == "Easting (EPSG:32632)"
        assert rows[3][0] == "Northing"


class TestAddProjectedColumns:
    def test_appends_easting_northing_columns(self, monkeypatch):
        import pycsamt.app.mapview._render as render_mod
        import pycsamt.app.mapview.callbacks.lines as lines_mod

        def fake_project(lons, lats, mode, zone, hem, epsg):
            return [1.0, 2.0], [3.0, 4.0], 32632

        monkeypatch.setattr(render_mod, "project_to_crs", fake_project)

        records = [
            {"Longitude": 7.75, "Latitude": 48.5},
            {"Longitude": 7.76, "Latitude": 48.51},
        ]
        cols = [{"name": "ID", "id": "ID"}]
        controls = {
            "crs_mode": "utm",
            "utm_zone": 32,
            "utm_hem": "N",
            "epsg": None,
        }
        new_records, new_cols = lines_mod._add_projected_columns(
            records, controls, cols
        )
        assert new_records[0]["E (32632)"] == 1.0
        assert new_records[0]["N"] == 3.0
        assert {"name": "E (32632)", "id": "E (32632)"} in new_cols
        assert {"name": "N", "id": "N"} in new_cols

    def test_projection_failure_returns_records_unchanged(self, monkeypatch):
        import pycsamt.app.mapview._render as render_mod
        import pycsamt.app.mapview.callbacks.lines as lines_mod

        monkeypatch.setattr(
            render_mod, "project_to_crs", lambda *a, **k: (None, None, 4326)
        )
        records = [{"Longitude": 7.75, "Latitude": 48.5}]
        cols = [{"name": "ID", "id": "ID"}]
        new_records, new_cols = lines_mod._add_projected_columns(
            records, {"crs_mode": "utm"}, cols
        )
        assert new_records == records
        assert new_cols == cols


class TestRegisterLines:
    def test_register_lines_is_callable(self):
        from pycsamt.app.mapview.callbacks.lines import register_lines

        assert callable(register_lines)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_LINES in cb_outputs
        assert IDs.LINE_PILLS in cb_outputs
        assert IDs.STORE_SELECTION in cb_outputs
        assert IDs.STATION_INSPECT in cb_outputs
