# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt.app.desktop.tools.elevation_tool``
(``_ElevWorker`` + ``ElevationEnrichDialog``).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — WILLY AMT EDIs; a 4-station subset is copied
into a session tmp dir and loaded via ``pycsamt.emtools.ensure_sites`` to
exercise ``_populate_stations`` against real EDI-backed ``Site`` objects.

Network safety
--------------
``_ElevWorker.run()`` calls ``pycsamt.gis.utils.get_elevation_from_api``
(imported *inside* the method body). Every test that exercises ``.run()``
monkeypatches that exact attribute on the ``pycsamt.gis.utils`` module
before calling it, so no test ever performs a real HTTP request.
``ElevationEnrichDialog`` orchestration tests never touch the real worker
at all: ``elevation_tool._ElevWorker`` is monkeypatched with a fake whose
``.start()`` synchronously emits through plain-callable signals — the same
idiom used in test_recompute_dlg.py / test_main_window_actions.py.

Known bug (found, NOT fixed — see final report)
-------------------------------------------------
``ElevationEnrichDialog._populate_stations`` unwraps each ``Site`` down to
its raw ``EDIFile`` (looking for an ``.edi``/``.site``/``.data`` attribute
that has ``.Z``) and then reads ``lat``/``latitude`` and ``lon``/
``longitude`` off *that* object. ``EDIFile`` exposes neither attribute
(coordinates live on the parsed ``>HEAD`` section, e.g.
``edi.get_section("head").lat``/``.long``, or via ``Site.coords``). As a
result, every real survey loaded through the normal app flow reports
"0 with coordinates" and the Fetch button is permanently disabled — the
dialog is non-functional against real data. See
``TestPopulateStationsRealData::test_real_sites_have_no_coords_bug``.
"""

from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.tools import elevation_tool
from pycsamt.app.desktop.tools.elevation_tool import (
    ElevationEnrichDialog,
    _ElevWorker,
)

# ── Paths / real-data fixtures ─────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_subset_sites(tmp_path_factory):
    """A 4-station Sites collection loaded from a small real-EDI subset."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    dest = tmp_path_factory.mktemp("elev_willy_subset")
    for f in sorted(_WILLY_L18.glob("*.edi"))[:4]:
        shutil.copy(f, dest / f.name)
    return ensure_sites(str(dest))


# ── Fake station-like items (bypass the real coordinate bug above) ─────────


class _FlatStationItem:
    """A pre-unwrapped item exposing station/lat/lon directly."""

    def __init__(self, station, lat, lon):
        self.station = station
        self.lat = lat
        self.lon = lon


def _flat_station_list():
    return [
        _FlatStationItem("S1", 10.0, 20.0),
        _FlatStationItem("S2", None, None),
        _FlatStationItem("S3", 30.0, 40.0),
    ]


class _RaisingZ:
    @property
    def Z(self):
        raise RuntimeError("boom-z")


class _RaisingEdiHolder:
    edi = _RaisingZ()


# ── Fake signal / fake worker for dialog-orchestration tests ───────────────


class _FakeSignal:
    def __init__(self):
        self._fn = None

    def connect(self, fn):
        self._fn = fn

    def emit(self, *a):
        if self._fn is not None:
            self._fn(*a)


def _fake_elev_worker_cls(
    *, progress_events=None, done_result=None, error_msg=None, auto_start=True
):
    class _FakeElevWorker:
        captured_stations = []
        captured_apis = []

        def __init__(self, stations, api):
            self.stations = stations
            self.api = api
            _FakeElevWorker.captured_stations.append(stations)
            _FakeElevWorker.captured_apis.append(api)
            self.progress = _FakeSignal()
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            self.started = False

        def start(self):
            self.started = True
            if not auto_start:
                return
            for ev in progress_events or []:
                self.progress.emit(*ev)
            if error_msg is not None:
                self.error.emit(error_msg)
            elif done_result is not None:
                self.done.emit(done_result)

    return _FakeElevWorker


# ── _ElevWorker — direct .run() calls, HTTP mocked ──────────────────────────


class TestElevWorkerRun:
    def test_full_success(self, qapp, monkeypatch):
        monkeypatch.setattr(
            "pycsamt.gis.utils.get_elevation_from_api",
            lambda lat, lon, api_name=None: lat + lon,
        )
        stations = [("S1", 10.0, 20.0), ("S2", 1.0, 2.0)]
        w = _ElevWorker(stations, api="open_meteo")

        progress_events = []
        done_results = []
        errors = []
        w.progress.connect(lambda *a: progress_events.append(a))
        w.done.connect(done_results.append)
        w.error.connect(errors.append)

        w.run()

        assert errors == []
        assert len(progress_events) == 2
        assert progress_events[0] == (1, 2, "S1", 30.0)
        assert progress_events[1] == (2, 2, "S2", 3.0)
        assert len(done_results) == 1
        assert done_results[0] == [
            ("S1", 10.0, 20.0, 30.0),
            ("S2", 1.0, 2.0, 3.0),
        ]

    def test_partial_failure_mixed_with_missing_coords(self, qapp, monkeypatch):
        def fake_api(lat, lon, api_name=None):
            if lat == 99.0:
                raise RuntimeError("station lookup failed")
            return 111.0

        monkeypatch.setattr(
            "pycsamt.gis.utils.get_elevation_from_api", fake_api
        )
        stations = [
            ("OK", 10.0, 20.0),
            ("NO_COORDS", None, None),
            ("API_FAILS", 99.0, 99.0),
        ]
        w = _ElevWorker(stations, api="open_topo_data")

        done_results = []
        w.done.connect(done_results.append)
        w.run()

        assert len(done_results) == 1
        results = done_results[0]
        assert results[0] == ("OK", 10.0, 20.0, 111.0)
        assert results[1][0] == "NO_COORDS"
        assert math.isnan(results[1][3])
        assert results[2][0] == "API_FAILS"
        assert math.isnan(results[2][3])

    def test_all_requests_fail(self, qapp, monkeypatch):
        def always_raises(lat, lon, api_name=None):
            raise ConnectionError("network down")

        monkeypatch.setattr(
            "pycsamt.gis.utils.get_elevation_from_api", always_raises
        )
        stations = [("S1", 1.0, 1.0), ("S2", 2.0, 2.0)]
        w = _ElevWorker(stations, api="open_meteo")

        progress_events = []
        done_results = []
        w.progress.connect(lambda *a: progress_events.append(a))
        w.done.connect(done_results.append)
        w.run()

        assert len(progress_events) == 2
        assert all(math.isnan(ev[3]) for ev in progress_events)
        assert len(done_results) == 1
        assert all(math.isnan(r[3]) for r in done_results[0])

    def test_import_error_emits_error_signal(self, qapp, monkeypatch):
        monkeypatch.delattr(
            "pycsamt.gis.utils.get_elevation_from_api", raising=True
        )
        w = _ElevWorker([("S1", 1.0, 1.0)], api="open_meteo")

        errors = []
        done_results = []
        w.error.connect(errors.append)
        w.done.connect(done_results.append)
        w.run()

        assert len(errors) == 1
        assert "Cannot import gis.utils" in errors[0]
        assert done_results == []

    def test_empty_station_list_emits_empty_done(self, qapp, monkeypatch):
        monkeypatch.setattr(
            "pycsamt.gis.utils.get_elevation_from_api",
            lambda lat, lon, api_name=None: 0.0,
        )
        w = _ElevWorker([], api="open_meteo")
        done_results = []
        w.done.connect(done_results.append)
        w.run()
        assert done_results == [[]]


# ── ElevationEnrichDialog — construction / _populate_stations ──────────────


class TestConstructionNoSites:
    def test_no_sites_status_text(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        assert dlg._status_lbl.text() == "No survey loaded."
        dlg.close()

    def test_no_sites_run_button_disabled(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        assert not dlg._run_btn.isEnabled()
        dlg.close()

    def test_no_sites_station_list_empty(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        assert dlg._station_list == []
        dlg.close()

    def test_export_button_disabled_initially(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        assert not dlg._export_btn.isEnabled()
        dlg.close()

    def test_api_combo_has_both_apis(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        items = [
            dlg._api_combo.itemText(i)
            for i in range(dlg._api_combo.count())
        ]
        assert items == ["open_meteo", "open_topo_data"]
        dlg.close()

    def test_table_has_four_columns(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        assert dlg._table.columnCount() == 4
        dlg.close()


class TestPopulateStationsFlatItems:
    """Items that already expose station/lat/lon directly (no unwrap needed)."""

    def test_station_count_and_coords(self, qapp):
        dlg = ElevationEnrichDialog(sites=_flat_station_list())
        assert dlg._station_list == [
            ("S1", 10.0, 20.0),
            ("S2", None, None),
            ("S3", 30.0, 40.0),
        ]
        dlg.close()

    def test_status_text_reports_counts(self, qapp):
        dlg = ElevationEnrichDialog(sites=_flat_station_list())
        assert dlg._status_lbl.text() == (
            "3 stations found, 2 with coordinates."
        )
        dlg.close()

    def test_run_button_enabled_when_some_have_coords(self, qapp):
        dlg = ElevationEnrichDialog(sites=_flat_station_list())
        assert dlg._run_btn.isEnabled()
        dlg.close()

    def test_run_button_disabled_when_none_have_coords(self, qapp):
        items = [_FlatStationItem("S1", None, None)]
        dlg = ElevationEnrichDialog(sites=items)
        assert not dlg._run_btn.isEnabled()
        assert "0 with coordinates" in dlg._status_lbl.text()
        dlg.close()

    def test_non_numeric_coords_become_none(self, qapp):
        items = [_FlatStationItem("S1", "not-a-number", "also-bad")]
        dlg = ElevationEnrichDialog(sites=items)
        assert dlg._station_list == [("S1", None, None)]
        assert "0 with coordinates" in dlg._status_lbl.text()
        dlg.close()

    def test_populate_exception_sets_status_and_disables_run(self, qapp):
        dlg = ElevationEnrichDialog(sites=[_RaisingEdiHolder()])
        assert dlg._status_lbl.text().startswith("Cannot read stations:")
        assert "boom-z" in dlg._status_lbl.text()
        assert not dlg._run_btn.isEnabled()
        dlg.close()


@pytest.mark.skipif(not _HAS_WILLY, reason="WILLY L18PLT data not available")
class TestPopulateStationsRealData:
    def test_real_sites_station_names_resolved(self, qapp, willy_subset_sites):
        dlg = ElevationEnrichDialog(sites=willy_subset_sites)
        assert len(dlg._station_list) == 4
        names = [n for n, _, _ in dlg._station_list]
        assert all(n for n in names)
        dlg.close()

    def test_real_sites_have_no_coords_bug(self, qapp, willy_subset_sites):
        """
        Documents a real bug (not fixed): _populate_stations unwraps each
        Site down to its raw EDIFile and reads .lat/.latitude/.lon/
        .longitude off it, but EDIFile exposes none of those — real
        coordinates live on the parsed >HEAD section
        (edi.get_section("head").lat/.long) or via Site.coords. So every
        station loaded through the normal app flow ends up with
        lat=lon=None here, and the Fetch button is disabled.
        """
        dlg = ElevationEnrichDialog(sites=willy_subset_sites)
        n_with_coords = sum(
            1 for _, la, lo in dlg._station_list if la is not None and lo is not None
        )
        assert n_with_coords == 0
        assert "0 with coordinates" in dlg._status_lbl.text()
        assert not dlg._run_btn.isEnabled()
        dlg.close()


# ── _on_run orchestration (fake worker) ─────────────────────────────────────


class TestOnRun:
    def test_constructs_worker_with_station_list_and_selected_api(
        self, qapp, monkeypatch
    ):
        fake_cls = _fake_elev_worker_cls(auto_start=False)
        monkeypatch.setattr(elevation_tool, "_ElevWorker", fake_cls)

        dlg = ElevationEnrichDialog(sites=_flat_station_list())
        dlg._api_combo.setCurrentText("open_topo_data")
        dlg._on_run()

        assert fake_cls.captured_stations[-1] == dlg._station_list
        assert fake_cls.captured_apis[-1] == "open_topo_data"
        assert dlg._worker.started
        dlg.close()

    def test_disables_buttons_and_resets_table_before_worker_runs(
        self, qapp, monkeypatch
    ):
        fake_cls = _fake_elev_worker_cls(auto_start=False)
        monkeypatch.setattr(elevation_tool, "_ElevWorker", fake_cls)

        dlg = ElevationEnrichDialog(sites=_flat_station_list())
        dlg._export_btn.setEnabled(True)  # simulate a prior successful run
        dlg._on_run()

        assert not dlg._run_btn.isEnabled()
        assert not dlg._export_btn.isEnabled()
        assert dlg._table.rowCount() == 0
        assert dlg._progress.value() == 0
        assert "Fetching elevations" in dlg._status_lbl.text()
        dlg.close()

    def test_full_success_flow_populates_table_and_reenables_buttons(
        self, qapp, monkeypatch
    ):
        progress_events = [
            (1, 2, "S1", 111.1),
            (2, 2, "S2", float("nan")),
        ]
        done_result = [
            ("S1", 10.0, 20.0, 111.1),
            ("S2", 30.0, 40.0, float("nan")),
        ]
        fake_cls = _fake_elev_worker_cls(
            progress_events=progress_events, done_result=done_result
        )
        monkeypatch.setattr(elevation_tool, "_ElevWorker", fake_cls)

        dlg = ElevationEnrichDialog(sites=_flat_station_list())
        dlg._on_run()

        assert dlg._table.rowCount() == 2
        assert dlg._table.item(0, 0).text() == "S1"
        assert dlg._table.item(0, 3).text() == "111.1"
        assert dlg._table.item(1, 3).text() == "ERROR"
        assert dlg._table.item(0, 1).text() == "10.000000"
        assert dlg._table.item(0, 2).text() == "20.000000"
        assert dlg._table.item(1, 1).text() == "30.000000"

        assert dlg._run_btn.isEnabled()
        assert dlg._export_btn.isEnabled()
        assert dlg._status_lbl.text() == "Done — 1 elevations fetched, 1 errors."
        assert dlg._results == done_result
        dlg.close()

    def test_error_flow_reenables_run_button_and_shows_message(
        self, qapp, monkeypatch
    ):
        fake_cls = _fake_elev_worker_cls(error_msg="Cannot import gis.utils: x")
        monkeypatch.setattr(elevation_tool, "_ElevWorker", fake_cls)

        dlg = ElevationEnrichDialog(sites=_flat_station_list())
        dlg._on_run()

        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: Cannot import gis.utils: x"
        # export stays disabled — no results ever arrived
        assert not dlg._export_btn.isEnabled()
        dlg.close()


# ── _on_progress / _on_done / _on_error direct unit tests ──────────────────


class TestOnProgress:
    def test_ok_row_uses_ok_background_and_formatted_value(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        dlg._on_progress(1, 3, "S1", 123.456)
        assert dlg._table.rowCount() == 1
        assert dlg._table.item(0, 0).text() == "S1"
        assert dlg._table.item(0, 3).text() == "123.5"
        assert dlg._table.item(0, 3).background().color() == elevation_tool._C_OK
        assert dlg._progress.maximum() == 3
        assert dlg._progress.value() == 1
        dlg.close()

    def test_error_row_uses_err_background_and_error_text(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        dlg._on_progress(2, 2, "S2", float("nan"))
        assert dlg._table.item(0, 3).text() == "ERROR"
        assert dlg._table.item(0, 3).background().color() == elevation_tool._C_ERR
        dlg.close()

    def test_multiple_progress_calls_append_rows(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        dlg._on_progress(1, 2, "S1", 1.0)
        dlg._on_progress(2, 2, "S2", 2.0)
        assert dlg._table.rowCount() == 2
        dlg.close()


class TestOnDone:
    def test_status_counts_ok_and_errors(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        results = [
            ("S1", 1.0, 2.0, 10.0),
            ("S2", 3.0, 4.0, float("nan")),
            ("S3", None, None, float("nan")),
        ]
        for i, (name, _, _, elev) in enumerate(results):
            dlg._on_progress(i + 1, len(results), name, elev)
        dlg._on_done(results)

        assert dlg._status_lbl.text() == "Done — 1 elevations fetched, 2 errors."
        assert dlg._results == results
        assert dlg._run_btn.isEnabled()
        assert dlg._export_btn.isEnabled()
        dlg.close()

    def test_fills_lat_lon_columns_including_missing(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        results = [("S1", 1.5, 2.5, 10.0), ("S2", None, None, float("nan"))]
        for i, (name, _, _, elev) in enumerate(results):
            dlg._on_progress(i + 1, len(results), name, elev)
        dlg._on_done(results)

        assert dlg._table.item(0, 1).text() == "1.500000"
        assert dlg._table.item(0, 2).text() == "2.500000"
        assert dlg._table.item(1, 1).text() == "—"
        assert dlg._table.item(1, 2).text() == "—"
        dlg.close()

    def test_empty_results_disables_export(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        dlg._on_done([])
        assert dlg._status_lbl.text() == "Done — 0 elevations fetched, 0 errors."
        assert not dlg._export_btn.isEnabled()
        dlg.close()


class TestOnError:
    def test_sets_status_and_reenables_run(self, qapp):
        dlg = ElevationEnrichDialog(sites=None)
        dlg._run_btn.setEnabled(False)
        dlg._on_error("timeout")
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: timeout"
        dlg.close()


# ── _on_export ───────────────────────────────────────────────────────────────


class TestOnExport:
    def test_no_results_does_not_open_save_dialog(self, qapp, monkeypatch):
        called = []
        monkeypatch.setattr(
            elevation_tool.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: called.append(1) or ("", "")),
        )
        dlg = ElevationEnrichDialog(sites=None)
        dlg._on_export()
        assert called == []
        dlg.close()

    def test_user_cancels_save_dialog_writes_nothing(self, qapp, monkeypatch, tmp_path):
        monkeypatch.setattr(
            elevation_tool.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        dlg = ElevationEnrichDialog(sites=None)
        dlg._results = [("S1", 1.0, 2.0, 10.0)]
        dlg._on_export()
        assert list(tmp_path.iterdir()) == []
        dlg.close()

    def test_writes_csv_with_expected_rows(self, qapp, monkeypatch, tmp_path):
        out_path = tmp_path / "out.csv"
        monkeypatch.setattr(
            elevation_tool.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(out_path), "CSV (*.csv)")),
        )
        dlg = ElevationEnrichDialog(sites=None)
        dlg._results = [
            ("S1", 1.5, 2.5, 10.25),
            ("S2", None, None, float("nan")),
        ]
        dlg._on_export()

        assert out_path.exists()
        with open(out_path, newline="", encoding="utf-8") as fh:
            rows = list(csv.reader(fh))
        assert rows[0] == ["station", "lat", "lon", "elevation_m"]
        assert rows[1] == ["S1", "1.500000", "2.500000", "10.25"]
        assert rows[2] == ["S2", "", "", ""]
        assert dlg._status_lbl.text() == f"Saved → {out_path}"
        dlg.close()

    def test_export_error_sets_status(self, qapp, monkeypatch, tmp_path):
        # A directory path makes open(path, "w") raise (IsADirectoryError on
        # POSIX, PermissionError on Windows) — either way an Exception the
        # handler must catch.
        monkeypatch.setattr(
            elevation_tool.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(tmp_path), "CSV (*.csv)")),
        )
        dlg = ElevationEnrichDialog(sites=None)
        dlg._results = [("S1", 1.0, 2.0, 10.0)]
        dlg._on_export()
        assert dlg._status_lbl.text().startswith("Export error:")
        dlg.close()
