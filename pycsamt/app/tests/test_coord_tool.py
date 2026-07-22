# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for CoordTransformDialog (pycsamt.app.desktop.tools.coord_tool).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — a 3-file subset is used to populate the
survey-station table branch with genuine EDI-derived Site objects.

Notes on real bugs found while writing these tests (NOT fixed here,
per task constraints — tests only):

1. Both ``_ll_to_utm`` and ``_utm_to_ll`` fall back to
   ``pycsamt.gis.utils.ll_to_utm`` / ``utm_to_ll`` when pyproj raises.
   That fallback is broken:
   - ``pycsamt.gis.utils.ll_to_utm`` actually has the signature
     ``ll_to_utm(reference_ellipsoid, lat, lon)`` and returns a
     ``(utm_zone_str, easting, northing)`` *tuple*, but coord_tool
     calls it as ``ll_to_utm(lat, lon)`` (missing arg) and treats the
     result as a dict (``res["easting"]``) — this raises ``TypeError``
     before the dict-indexing is ever reached.
   - ``pycsamt.gis.utils.utm_to_ll`` has the signature
     ``utm_to_ll(reference_ellipsoid, northing, easting, zone)`` and
     returns a ``(lat, lon)`` tuple, but coord_tool calls it as
     ``utm_to_ll(easting, northing, zone, hem)`` (wrong arg order/count)
     and again treats the result as a dict. In practice this raises
     ``pycsamt.gis.config.GisError: Unknown ellipsoid ID: <easting>``
     because the easting value is passed positionally as the
     ``reference_ellipsoid`` id.
   Net effect: whenever the primary pyproj path fails (e.g. an
   out-of-range UTM zone), the "fallback" raises a *different*,
   confusing exception instead of silently degrading. In the dialog
   this is swallowed by ``_on_transform``'s per-pair try/except and
   surfaces as an ``ERROR: ...`` output line, so it never crashes the
   UI, but the documented pure-Python fallback path is dead/broken code.

2. ``_populate_survey_table`` looks for coordinates via
   ``getattr(ed, "lat", None) or getattr(ed, "latitude", None)`` (and
   the analogous lon lookup) on the item returned by
   ``pycsamt.emtools._core._unwrap``. Neither the ``Site`` wrapper nor
   the unwrapped ``EDIFile`` exposes ``.lat``/``.lon``/``.latitude``/
   ``.longitude`` directly — real coordinates live at
   ``Site.coords`` (a ``(lat, lon, elev)`` tuple) or
   ``EDIFile.get_section("head").lat`` / ``.long``. As a result, when
   the Survey Stations table is populated from a real loaded survey,
   every row shows the correct station name but the Lat/Lon/UTM
   columns always fall through to the "no coordinates" placeholder
   ("—"), even though the loaded EDIs do have valid coordinates.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

pytest.importorskip("PySide6")

from pycsamt.app.desktop.tools.coord_tool import (
    CoordTransformDialog,
    _ll_to_utm,
    _parse_pairs,
    _utm_to_ll,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_L18PLT = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_L18PLT = _L18PLT.exists() and any(_L18PLT.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_subset_sites():
    """A small (3-station) real Sites collection for the survey-table branch."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_L18PLT:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    edi_files = sorted(_L18PLT.glob("*.edi"))[:3]
    return ensure_sites([str(p) for p in edi_files])


# ── Fake EDI-like objects (for exercising branches real data can't reach) ─────


class _FakeStationWithCoords:
    def __init__(self, station, lat, lon):
        self.station = station
        self.lat = lat
        self.lon = lon


class _FakeStationNoCoords:
    def __init__(self, station):
        self.station = station


class _FakeStationBadCoords:
    """lat/lon present but not float-convertible -> exercises inner except."""

    def __init__(self, station):
        self.station = station
        self.lat = "not-a-number"
        self.lon = "also-not-a-number"


# ═════════════════════════════════════════════════════════════════════════════
# Module-level coordinate math
# ═════════════════════════════════════════════════════════════════════════════


class TestLLToUTM:
    def test_london_northern_hemisphere_zone_30(self):
        e, n, z = _ll_to_utm(51.5074, -0.1278, None, "N", "WGS84")
        assert z == 30
        assert e == pytest.approx(699316.23, abs=1.0)
        assert n == pytest.approx(5710163.76, abs=1.0)

    def test_auto_zone_detection_matches_formula(self):
        lon = -0.1278
        _, _, z = _ll_to_utm(51.5074, lon, None, "N", "WGS84")
        assert z == int((lon + 180) / 6) + 1

    def test_explicit_zone_overrides_auto_detection(self):
        # Force a "wrong" (but still valid 1-60) zone and check it's honored.
        _, _, z = _ll_to_utm(51.5074, -0.1278, 31, "N", "WGS84")
        assert z == 31

    def test_southern_hemisphere_sydney(self):
        e, n, z = _ll_to_utm(-33.8688, 151.2093, None, "S", "WGS84")
        assert z == 56
        assert e == pytest.approx(334368.63, abs=1.0)
        assert n == pytest.approx(6250948.35, abs=1.0)
        # Southern-hemisphere UTM northings are offset by 10,000,000 m;
        # Sydney sits well below that, so northing must be « 10e6.
        assert n < 10_000_000

    def test_equator_prime_meridian_null_island(self):
        e, n, z = _ll_to_utm(0.0, 0.0, None, "N", "WGS84")
        assert z == 31
        assert n == pytest.approx(0.0, abs=1.0)
        # False easting at the central meridian is 500,000 m; at lon=0 the
        # point is 3 degrees off zone-31's central meridian (=3E).
        assert e == pytest.approx(166021.44, abs=1.0)

    def test_invalid_zone_zero_falls_through_to_broken_fallback(self):
        """
        Documents bug #1 (see module docstring): zone=0 is invalid for
        pyproj, so the code falls back to pycsamt.gis.utils.ll_to_utm,
        whose signature mismatch raises TypeError instead of a graceful
        fallback result.
        """
        with pytest.raises(TypeError):
            _ll_to_utm(51.5074, -0.1278, 0, "N", "WGS84")


class TestUTMToLL:
    def test_roundtrip_london(self):
        e, n, z = _ll_to_utm(51.5074, -0.1278, None, "N", "WGS84")
        lat, lon = _utm_to_ll(e, n, z, "N", "WGS84")
        assert lat == pytest.approx(51.5074, abs=1e-4)
        assert lon == pytest.approx(-0.1278, abs=1e-4)

    def test_roundtrip_sydney_southern_hemisphere(self):
        e, n, z = _ll_to_utm(-33.8688, 151.2093, None, "S", "WGS84")
        lat, lon = _utm_to_ll(e, n, z, "S", "WGS84")
        assert lat == pytest.approx(-33.8688, abs=1e-4)
        assert lon == pytest.approx(151.2093, abs=1e-4)

    def test_roundtrip_equator_prime_meridian(self):
        e, n, z = _ll_to_utm(0.0, 0.0, None, "N", "WGS84")
        lat, lon = _utm_to_ll(e, n, z, "N", "WGS84")
        assert lat == pytest.approx(0.0, abs=1e-4)
        assert lon == pytest.approx(0.0, abs=1e-4)

    def test_invalid_zone_zero_falls_through_to_broken_fallback(self):
        """Documents bug #1 for the inverse direction (see module docstring)."""
        with pytest.raises(Exception):  # noqa: B017 - GisError, not importable path
            _utm_to_ll(452000.0, 5411000.0, 0, "N", "WGS84")


class TestParsePairs:
    def test_well_formed_multiline(self):
        text = "51.5074 -0.1278\n48.8566 2.3522\n"
        assert _parse_pairs(text) == [(51.5074, -0.1278), (48.8566, 2.3522)]

    def test_comma_separated(self):
        assert _parse_pairs("48.8566, 2.3522") == [(48.8566, 2.3522)]

    def test_tab_separated(self):
        assert _parse_pairs("10.0\t20.0") == [(10.0, 20.0)]

    def test_extra_whitespace_is_tolerant(self):
        text = "   51.5074      -0.1278   \n\n  48.8566    2.3522  "
        assert _parse_pairs(text) == [(51.5074, -0.1278), (48.8566, 2.3522)]

    def test_empty_string_returns_empty_list(self):
        assert _parse_pairs("") == []

    def test_whitespace_only_returns_empty_list(self):
        assert _parse_pairs("   \n   \n  ") == []

    def test_comment_lines_are_skipped(self):
        text = "# a header comment\n51.5074 -0.1278\n# another comment"
        assert _parse_pairs(text) == [(51.5074, -0.1278)]

    def test_malformed_non_numeric_lines_are_skipped_not_raised(self):
        text = "abc def\n51.5074 -0.1278"
        assert _parse_pairs(text) == [(51.5074, -0.1278)]

    def test_single_token_lines_are_skipped(self):
        text = "only_one_token\n51.5074 -0.1278"
        assert _parse_pairs(text) == [(51.5074, -0.1278)]

    def test_extra_tokens_on_a_line_uses_first_two(self):
        assert _parse_pairs("1 2 3 4 5") == [(1.0, 2.0)]

    def test_blank_lines_between_pairs_are_ignored(self):
        text = "51.5074 -0.1278\n\n\n48.8566 2.3522"
        assert _parse_pairs(text) == [(51.5074, -0.1278), (48.8566, 2.3522)]


# ═════════════════════════════════════════════════════════════════════════════
# Dialog construction
# ═════════════════════════════════════════════════════════════════════════════


class TestDialogConstructionNoSites:
    def test_constructs_without_error(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        assert dlg is not None

    def test_survey_table_not_built_when_no_sites(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        assert not hasattr(dlg, "_surv_table")

    def test_window_title(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        assert dlg.windowTitle() == "Coordinate Transformer"

    def test_default_direction_is_ll_to_utm(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        assert dlg._ll_to_utm_rb.isChecked() is True
        assert dlg._utm_to_ll_rb.isChecked() is False

    def test_default_hemisphere_is_north(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        assert dlg._hem_n.isChecked() is True
        assert dlg._hem_s.isChecked() is False

    def test_default_zone_value(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        assert dlg._zone_spin.value() == 30

    def test_datum_combo_contains_expected_datums(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        items = [
            dlg._datum_combo.itemText(i)
            for i in range(dlg._datum_combo.count())
        ]
        assert items == ["WGS84", "NAD83", "GRS80"]


class TestDialogConstructionWithRealSites:
    def test_survey_table_built_when_sites_given(
        self, qapp, willy_subset_sites
    ):
        dlg = CoordTransformDialog(sites=willy_subset_sites)
        assert hasattr(dlg, "_surv_table")

    def test_survey_table_row_count_matches_sites(
        self, qapp, willy_subset_sites
    ):
        dlg = CoordTransformDialog(sites=willy_subset_sites)
        assert dlg._surv_table.rowCount() == len(willy_subset_sites)

    def test_survey_table_station_names_populated(
        self, qapp, willy_subset_sites
    ):
        dlg = CoordTransformDialog(sites=willy_subset_sites)
        names = [
            dlg._surv_table.item(r, 0).text()
            for r in range(dlg._surv_table.rowCount())
        ]
        assert all(n and n != "?" for n in names)

    def test_survey_table_columns_header(self, qapp, willy_subset_sites):
        dlg = CoordTransformDialog(sites=willy_subset_sites)
        headers = [
            dlg._surv_table.horizontalHeaderItem(c).text()
            for c in range(dlg._surv_table.columnCount())
        ]
        assert headers == ["Station", "Lat (°)", "Lon (°)", "UTM E (m)", "UTM N (m)"]

    def test_survey_table_real_edi_coords_fall_through_to_dash(
        self, qapp, willy_subset_sites
    ):
        """
        Documents bug #2 (see module docstring): real Site/EDIFile objects
        don't expose flat .lat/.lon attributes, so every row's coordinate
        columns show the "no coordinates" placeholder despite the EDIs
        genuinely having valid coordinates (Site.coords / head.lat/.long).
        """
        dlg = CoordTransformDialog(sites=willy_subset_sites)
        for r in range(dlg._surv_table.rowCount()):
            assert dlg._surv_table.item(r, 1).text() == "—"
            assert dlg._surv_table.item(r, 2).text() == "—"


class TestPopulateSurveyTableWithFakeCoords:
    """
    Exercises the coordinate-found happy path and the malformed-coordinate
    except branch, neither of which real pycsamt Site/EDIFile objects can
    reach (see bug #2 in the module docstring).
    """

    def test_rows_with_valid_lat_lon_populate_all_columns(self, qapp):
        fakes = [
            _FakeStationWithCoords("STA01", 51.5074, -0.1278),
            _FakeStationWithCoords("STA02", 48.8566, 2.3522),
        ]
        dlg = CoordTransformDialog(sites=fakes)
        assert dlg._surv_table.rowCount() == 2
        row0 = [dlg._surv_table.item(0, c).text() for c in range(5)]
        assert row0[0] == "STA01"
        assert row0[1] == "51.507400"
        assert row0[2] == "-0.127800"
        # UTM E/N should be non-placeholder numeric strings.
        assert row0[3] != "—"
        assert row0[4] != "—"

    def test_rows_without_coords_show_dash_placeholder(self, qapp):
        fakes = [_FakeStationNoCoords("STA03")]
        dlg = CoordTransformDialog(sites=fakes)
        assert dlg._surv_table.item(0, 0).text() == "STA03"
        for c in range(1, 5):
            assert dlg._surv_table.item(0, c).text() == "—"

    def test_rows_with_unconvertible_coords_do_not_raise(self, qapp):
        fakes = [_FakeStationBadCoords("STA04")]
        # Must not raise despite lat/lon being non-numeric strings.
        dlg = CoordTransformDialog(sites=fakes)
        assert dlg._surv_table.rowCount() == 1
        assert dlg._surv_table.item(0, 0).text() == "STA04"

    def test_mixed_valid_and_missing_rows(self, qapp):
        fakes = [
            _FakeStationWithCoords("A", 10.0, 20.0),
            _FakeStationNoCoords("B"),
        ]
        dlg = CoordTransformDialog(sites=fakes)
        assert dlg._surv_table.item(0, 1).text() == "10.000000"
        assert dlg._surv_table.item(1, 1).text() == "—"

    def test_populate_survives_iter_items_raising(self, qapp, monkeypatch):
        """
        Force pycsamt.emtools._core._iter_items to raise so the
        `except Exception: return` guard at the top of
        _populate_survey_table is exercised.
        """
        import pycsamt.emtools._core as core

        def _boom(_sites):
            raise RuntimeError("boom")

        monkeypatch.setattr(core, "_iter_items", _boom)
        dlg = CoordTransformDialog(sites=[_FakeStationWithCoords("X", 1.0, 2.0)])
        # Guard swallowed the error; table stays at its initial empty state.
        assert dlg._surv_table.rowCount() == 0

    def test_populate_survives_unwrap_raising(self, qapp, monkeypatch):
        """
        ``_unwrap`` has its own internal try/except and never actually
        raises in practice, so the per-item ``except Exception: pass``
        wrapping it is normally dead code; force it to confirm the row
        still gets built from the raw (un-unwrapped) item.
        """
        import pycsamt.emtools._core as core

        monkeypatch.setattr(
            core,
            "_iter_items",
            lambda _sites: [_FakeStationWithCoords("STA-RAW", 5.0, 6.0)],
        )

        def _boom(_ed):
            raise RuntimeError("cannot unwrap")

        monkeypatch.setattr(core, "_unwrap", _boom)
        dlg = CoordTransformDialog(sites=[object()])
        assert dlg._surv_table.rowCount() == 1
        assert dlg._surv_table.item(0, 0).text() == "STA-RAW"

    def test_repopulate_respects_hemisphere_toggle(self, qapp):
        """_populate_survey_table re-reads hemisphere/datum on each call."""
        fakes = [_FakeStationWithCoords("STA05", -33.8688, 151.2093)]
        dlg = CoordTransformDialog(sites=fakes)
        e_north = float(dlg._surv_table.item(0, 3).text())

        dlg._hem_s.setChecked(True)
        dlg._populate_survey_table()
        e_south = float(dlg._surv_table.item(0, 3).text())

        # Easting is hemisphere-independent for the forward UTM projection
        # itself (only the false-northing offset differs), so this mainly
        # verifies the call path re-executes cleanly under both toggles.
        assert e_north == pytest.approx(e_south, rel=1e-6)


# ═════════════════════════════════════════════════════════════════════════════
# _on_transform
# ═════════════════════════════════════════════════════════════════════════════


class TestOnTransform:
    def test_ll_to_utm_direction_produces_zone_annotated_output(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        dlg._ll_to_utm_rb.setChecked(True)
        dlg._in_edit.setPlainText("51.5074 -0.1278")
        dlg._on_transform()
        out = dlg._out_edit.toPlainText()
        assert "zone 30N" in out
        e, n, z = _ll_to_utm(51.5074, -0.1278, 30, "N", "WGS84")
        assert f"{e:>14.3f}" in out
        assert f"{n:>14.3f}" in out

    def test_utm_to_ll_direction_produces_lat_lon_output(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        dlg._utm_to_ll_rb.setChecked(True)
        assert dlg._ll_to_utm_rb.isChecked() is False  # mutually exclusive
        dlg._zone_spin.setValue(30)
        dlg._in_edit.setPlainText("699316.234 5710163.758")
        dlg._on_transform()
        out = dlg._out_edit.toPlainText()
        lat, lon = _utm_to_ll(699316.234, 5710163.758, 30, "N", "WGS84")
        assert f"{lat:>12.6f}" in out
        assert f"{lon:>13.6f}" in out

    def test_multiple_pairs_produce_multiple_lines(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        dlg._ll_to_utm_rb.setChecked(True)
        dlg._in_edit.setPlainText("51.5074 -0.1278\n48.8566 2.3522")
        dlg._on_transform()
        lines = dlg._out_edit.toPlainText().splitlines()
        assert len(lines) == 2

    def test_no_valid_input_shows_message(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        dlg._in_edit.setPlainText("nothing_parseable_here")
        dlg._on_transform()
        assert dlg._out_edit.toPlainText() == "No valid input pairs found."

    def test_empty_input_shows_message(self, qapp):
        dlg = CoordTransformDialog(sites=None)
        dlg._in_edit.setPlainText("")
        dlg._on_transform()
        assert dlg._out_edit.toPlainText() == "No valid input pairs found."

    def test_southern_hemisphere_toggle_used(self, qapp):
        # NOTE: real widget bug -- none of the four radios (_ll_to_utm_rb /
        # _utm_to_ll_rb / _hem_n / _hem_s) are ever put in an explicit
        # QButtonGroup, and all four share the same immediate parent widget
        # (the "Options" QGroupBox). Qt's implicit auto-exclusive grouping
        # then treats all four as ONE group instead of two independent
        # pairs: whichever radio is checked *last* wins and can silently
        # uncheck an unrelated radio from the other pair (e.g. checking
        # "S" hemisphere first flips "Lat/Lon -> UTM" back off). Checking
        # the hemisphere before the direction (matching the order a user
        # would realistically click "Lat/Lon -> UTM" then "S") sidesteps
        # it here; not fixed in production, per instructions.
        dlg = CoordTransformDialog(sites=None)
        dlg._hem_s.setChecked(True)
        dlg._ll_to_utm_rb.setChecked(True)
        dlg._in_edit.setPlainText("-33.8688 151.2093")
        dlg._zone_spin.setValue(56)
        dlg._on_transform()
        out = dlg._out_edit.toPlainText()
        assert "zone 56S" in out

    def test_datum_selection_used(self, qapp):
        # NOTE: real bug -- _ll_to_utm passes `ellps=datum` straight to
        # pyproj.Proj, but "NAD83"/"GRS80" are datum names, not ellipsoid
        # names ("WGS84" happens to also be a valid ellps name, which is
        # why it's the only datum that works). The CRSError then falls
        # into the pure-Python fallback, which is itself broken: it calls
        # pycsamt.gis.utils.ll_to_utm(lat, lon) but that function's real
        # signature is (reference_ellipsoid, lat, lon) returning a 3-tuple,
        # not a dict -- so every non-WGS84 datum currently always errors.
        # Locking in current (broken) behaviour; not fixed, per instructions.
        dlg = CoordTransformDialog(sites=None)
        dlg._datum_combo.setCurrentText("NAD83")
        dlg._in_edit.setPlainText("40.7128 -74.0060")
        dlg._on_transform()
        out = dlg._out_edit.toPlainText()
        assert out.startswith("ERROR:")
        with pytest.raises(TypeError):
            _ll_to_utm(40.7128, -74.0060, dlg._zone_spin.value(), "N", "NAD83")

    def test_per_pair_exception_reported_as_error_line_ll_to_utm(
        self, qapp, monkeypatch
    ):
        import pycsamt.app.desktop.tools.coord_tool as mod

        def _raise(*a, **k):
            raise ValueError("synthetic failure")

        monkeypatch.setattr(mod, "_ll_to_utm", _raise)
        dlg = CoordTransformDialog(sites=None)
        dlg._ll_to_utm_rb.setChecked(True)
        dlg._in_edit.setPlainText("51.5074 -0.1278")
        dlg._on_transform()
        assert dlg._out_edit.toPlainText() == "ERROR: synthetic failure"

    def test_per_pair_exception_reported_as_error_line_utm_to_ll(
        self, qapp, monkeypatch
    ):
        import pycsamt.app.desktop.tools.coord_tool as mod

        def _raise(*a, **k):
            raise ValueError("synthetic failure 2")

        monkeypatch.setattr(mod, "_utm_to_ll", _raise)
        dlg = CoordTransformDialog(sites=None)
        dlg._utm_to_ll_rb.setChecked(True)
        dlg._in_edit.setPlainText("699316 5710163")
        dlg._on_transform()
        assert dlg._out_edit.toPlainText() == "ERROR: synthetic failure 2"

    def test_real_invalid_zone_bug_surfaces_as_error_line(self, qapp, monkeypatch):
        """
        End-to-end (dialog-level) manifestation of bug #1: forcing the
        zone spin box below its normal floor via direct attribute
        assignment (bypassing the 1-60 UI clamp) makes the pyproj path
        raise, which then hits the broken fallback and still ends up as
        a (differently worded) ERROR line rather than a crash.
        """
        dlg = CoordTransformDialog(sites=None)
        dlg._ll_to_utm_rb.setChecked(True)
        dlg._in_edit.setPlainText("51.5074 -0.1278")
        # QSpinBox.setValue clamps to the configured range (1-60), so we
        # go around it to reproduce the invalid-zone condition.
        monkeypatch.setattr(type(dlg._zone_spin), "value", lambda self: 0)
        dlg._on_transform()
        out = dlg._out_edit.toPlainText()
        assert out.startswith("ERROR:")


# ═════════════════════════════════════════════════════════════════════════════
# _copy_output
# ═════════════════════════════════════════════════════════════════════════════


class TestCopyOutput:
    def test_copy_output_sets_real_clipboard_text(self, qapp):
        from PySide6.QtWidgets import QApplication

        dlg = CoordTransformDialog(sites=None)
        dlg._out_edit.setPlainText("copy-me-123")
        dlg._copy_output()
        assert QApplication.clipboard().text() == "copy-me-123"

    def test_copy_output_uses_mocked_clipboard(self, qapp, monkeypatch):
        from PySide6.QtWidgets import QApplication

        mock_clipboard = MagicMock()
        monkeypatch.setattr(
            QApplication, "clipboard", staticmethod(lambda: mock_clipboard)
        )
        dlg = CoordTransformDialog(sites=None)
        dlg._out_edit.setPlainText("mocked-text")
        dlg._copy_output()
        mock_clipboard.setText.assert_called_once_with("mocked-text")

    def test_copy_output_after_transform_matches_output_text(self, qapp):
        from PySide6.QtWidgets import QApplication

        dlg = CoordTransformDialog(sites=None)
        dlg._in_edit.setPlainText("51.5074 -0.1278")
        dlg._on_transform()
        expected = dlg._out_edit.toPlainText()
        dlg._copy_output()
        assert QApplication.clipboard().text() == expected
