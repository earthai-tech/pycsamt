# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.data (load modal, EDI upload, station
table/list).

Real data
---------
The ``simulated_edi`` session fixture (conftest.py) provides a real,
minimal, parseable EDI file; it's base64-encoded into the same
``data:<mime>;base64,<data>`` shape ``dcc.Upload`` produces, so
``load_data`` is exercised end to end through the real
``DataController``/``decode_upload_to_tempdir`` pipeline rather than
mocked.
"""

from __future__ import annotations

import base64

import pytest
from dash import no_update

from pycsamt.app.web.callbacks.data import (
    _detected_summary,
    _entry_line,
    _entry_names,
    _ext_counts,
    _filtered_upload_entries,
    _hex_to_rgb,
    _infer_line_counts,
    _line_from_path,
    _merge_store,
    _modal_title,
    _mode_hint,
    _normalise_names,
    _preflight_children,
    _sanitize_upload_name,
    _upload_entries,
)
from pycsamt.app.web.layout import IDs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(
        k for k in web_app.callback_map if all(s in k for s in substrings)
    )
    return _unwrap(web_app.callback_map[key])


def _cb_by_input(web_app, output_substr, input_id_substr):
    """Disambiguate callbacks sharing an (allow_duplicate) output by a
    substring of the real Input component id (including Dash's
    pattern-matching {"index":..., "type":...} ids, which stringify)."""
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(
            input_id_substr in str(i.get("id"))
            for i in v.get("inputs", [])
        ):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} "
        f"input~={input_id_substr!r}"
    )


def _set_triggered(prop_id):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(
        AttributeDict(triggered_inputs=[{"prop_id": prop_id}])
    )


def _clear_triggered():
    import dash._callback_context as cc

    cc.context_value.set({})


@pytest.fixture
def edi_data_uri(simulated_edi):
    raw = simulated_edi.read_bytes()
    b64 = base64.b64encode(raw).decode("ascii")
    return f"data:application/octet-stream;base64,{b64}"


# ── Pure helpers ─────────────────────────────────────────────────────────────


class TestNormaliseNames:
    def test_none_returns_empty(self):
        assert _normalise_names(None) == []

    def test_string_wrapped_in_list(self):
        assert _normalise_names("a.edi") == ["a.edi"]

    def test_list_stringified_and_falsy_dropped(self):
        assert _normalise_names(["a.edi", None, "", "b.edi"]) == [
            "a.edi",
            "b.edi",
        ]


class TestSanitizeUploadName:
    def test_empty_uses_fallback(self):
        assert _sanitize_upload_name("", "fallback.edi") == "fallback.edi"

    def test_strips_unsafe_chars(self):
        assert _sanitize_upload_name("a b?.edi", "x") == "a b_.edi"

    def test_preserves_relative_folders(self):
        assert (
            _sanitize_upload_name("L1/STA01.edi", "x") == "L1/STA01.edi"
        )

    def test_backslashes_normalised_to_slashes(self):
        assert (
            _sanitize_upload_name("L1\\STA01.edi", "x") == "L1/STA01.edi"
        )


class TestUploadEntries:
    def test_no_contents_returns_empty(self):
        assert _upload_entries(None, ["a.edi"]) == []

    def test_single_string_content_wrapped(self):
        entries = _upload_entries("data:...", "a.edi")
        assert len(entries) == 1
        assert entries[0]["filename"] == "a.edi"
        assert entries[0]["source"] == "upload"

    def test_multiple_contents_with_fallback_name(self):
        entries = _upload_entries(["c1", "c2"], ["a.edi"])
        assert entries[1]["original"] == "file_2.edi"


class TestExtCounts:
    def test_counts_recognised_and_other(self):
        counts = _ext_counts(["a.edi", "b.AVG", "c.j", "d.txt", "noext"])
        assert counts == {"edi": 1, "avg": 1, "j": 1, "other": 2}


class TestLineFromPath:
    def test_flat_file(self):
        assert _line_from_path("STA01.edi") == "flat files"

    def test_two_parts(self):
        assert _line_from_path("L1/STA01.edi") == "L1"

    def test_three_parts_uses_second(self):
        assert _line_from_path("survey/L1/STA01.edi") == "L1"

    def test_backslashes(self):
        assert _line_from_path("survey\\L1\\STA01.edi") == "L1"


class TestInferLineCounts:
    def test_groups_by_line(self):
        counts = _infer_line_counts(["L1/a.edi", "L1/b.edi", "L2/c.edi"])
        assert counts == {"L1": 2, "L2": 1}


class TestEntryLineAndFiltered:
    def test_entry_line_uses_filename_then_original(self):
        assert _entry_line({"filename": "L1/a.edi"}) == "L1"
        assert _entry_line({"original": "L2/b.edi"}) == "L2"

    def test_filtered_no_selection_returns_all(self):
        entries = [{"filename": "L1/a.edi"}, {"filename": "L2/b.edi"}]
        assert _filtered_upload_entries(entries, None) == entries

    def test_filtered_restricts_to_selected_lines(self):
        entries = [{"filename": "L1/a.edi"}, {"filename": "L2/b.edi"}]
        result = _filtered_upload_entries(entries, ["L1"])
        assert result == [{"filename": "L1/a.edi"}]


class TestEntryNames:
    def test_prefers_filename_over_original(self):
        assert _entry_names([{"filename": "a", "original": "b"}]) == ["a"]

    def test_falls_back_to_original(self):
        assert _entry_names([{"original": "b"}]) == ["b"]


class TestModalTitleAndHint:
    def test_append_title_and_hint(self):
        title = _modal_title("append")
        assert "Add Lines" in str(title)
        assert "merged" in _mode_hint("append")

    def test_replace_title_and_hint(self):
        title = _modal_title("replace")
        assert "Load Survey" in str(title)
        assert "replaced" in _mode_hint("replace")


class TestPreflightAndSummary:
    def test_preflight_empty(self):
        result = _preflight_children([], mode="replace", source="file")
        assert "Choose a folder" in str(result)

    def test_preflight_with_names(self):
        names = ["L1/a.edi", "L1/b.edi", "L2/c.avg"]
        result = _preflight_children(names, mode="replace", source="file")
        text = str(result)
        assert "EDI 2" in text
        assert "AVG 1" in text

    def test_preflight_hidden_lines_count(self):
        names = [f"L{i}/s.edi" for i in range(8)]
        result = _preflight_children(names, mode="append", source="folder")
        assert "+ 2 more line group" in str(result)

    def test_detected_summary_empty(self):
        result = _detected_summary([], mode="replace", source="file")
        assert "No data selected" in str(result)

    def test_detected_summary_with_names(self):
        result = _detected_summary(
            ["L1/a.edi", "L2/b.edi"], mode="append", source="folder"
        )
        text = str(result)
        assert "2 files" in text
        assert "2 line groups" in text
        assert "add lines" in text
        assert "folder" in text


class TestHexToRgb:
    def test_conversion(self):
        assert _hex_to_rgb("#ff0080") == (255, 0, 128)


class TestMergeStore:
    def test_merge_deduplicates_by_id_new_wins(self):
        existing = {
            "station_records": [
                {"ID": "S1", "Line": "L1"},
                {"ID": "S2", "Line": "L1"},
            ],
            "data_dir": "/old",
        }
        new_records = [{"ID": "S2", "Line": "L2"}, {"ID": "S3", "Line": "L2"}]

        class _FakeSites:
            pass

        import pycsamt.app.web.callbacks.data as data_mod

        merged_sites_obj = object()
        import unittest.mock as mock

        with mock.patch.object(
            data_mod, "cache_merge_sites", return_value=merged_sites_obj
        ) as m:
            store, sites = _merge_store(
                existing,
                new_records,
                _FakeSites(),
                "sess",
                new_dir="/new",
                n_new_lines=1,
                new_line_counts={"L2": 2},
            )
        assert store["n_stations"] == 3
        assert store["data_dir"] == "/old + /new"
        by_id = {r["ID"]: r for r in store["station_records"]}
        assert by_id["S2"]["Line"] == "L2"  # new wins
        assert sites is merged_sites_obj

    def test_merge_same_dir_not_duplicated(self):
        existing = {"station_records": [], "data_dir": "/same"}

        import unittest.mock as mock

        import pycsamt.app.web.callbacks.data as data_mod

        with mock.patch.object(data_mod, "cache_merge_sites"):
            store, _sites = _merge_store(
                existing,
                [],
                object(),
                "sess",
                new_dir="/same",
                n_new_lines=0,
                new_line_counts={},
            )
        assert store["data_dir"] == "/same"


# ── toggle_modal ─────────────────────────────────────────────────────────────


class TestToggleModal:
    def _fn(self, web_app):
        return _cb_multi(
            web_app, f"{IDs.MODAL_LOAD}.is_open", "load-mode-store.data"
        )

    def test_close_button_closes_modal(self, web_app):
        _set_triggered("modal-close-btn.n_clicks")
        try:
            result = self._fn(web_app)(
                0, 0, 1, 0, 0, 0, 0, 0, True
            )
        finally:
            _clear_triggered()
        is_open, *_rest = result
        assert is_open is False

    def test_add_line_button_opens_in_append_mode(self, web_app):
        _set_triggered(f"{IDs.BTN_ADD_LINE}.n_clicks")
        try:
            result = self._fn(web_app)(0, 1, 0, 0, 0, 0, 0, 0, False)
        finally:
            _clear_triggered()
        is_open, mode, replace_cls, append_cls, title, hint = result
        assert is_open is True
        assert mode == "append"
        assert "active" in append_cls
        assert "active" not in replace_cls

    def test_load_button_opens_in_replace_mode(self, web_app):
        _set_triggered(f"{IDs.BTN_LOAD}.n_clicks")
        try:
            result = self._fn(web_app)(1, 0, 0, 0, 0, 0, 0, 0, False)
        finally:
            _clear_triggered()
        is_open, mode, replace_cls, append_cls, title, hint = result
        assert mode == "replace"
        assert "active" in replace_cls


# ── capture_upload_selection ─────────────────────────────────────────────────


class TestCaptureUploadSelection:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, "load-upload-selection.data", IDs.UPLOAD
        )

    def test_from_folder_store(self, web_app, edi_data_uri):
        _set_triggered(f"{IDs.FOLDER_UPLOAD_STORE}.data")
        try:
            entries = self._fn(web_app)(
                None,
                {"contents": [edi_data_uri], "filenames": ["L1/a.edi"]},
                None,
            )
        finally:
            _clear_triggered()
        assert entries[0]["source"] == "folder"

    def test_from_direct_upload(self, web_app, edi_data_uri):
        _set_triggered(f"{IDs.UPLOAD}.contents")
        try:
            entries = self._fn(web_app)(edi_data_uri, None, "a.edi")
        finally:
            _clear_triggered()
        assert entries[0]["source"] == "upload"

    def test_unknown_trigger_returns_empty(self, web_app):
        _set_triggered("something-else.value")
        try:
            entries = self._fn(web_app)(None, None, None)
        finally:
            _clear_triggered()
        assert entries == []


# ── edit_upload_selection ────────────────────────────────────────────────────


class TestEditUploadSelection:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, "load-upload-selection.data", "load-file-remove"
        )

    def test_clear_initial_render_ignored(self, web_app):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        fn = self._fn(web_app)
        cc.context_value.set(
            AttributeDict(
                triggered_inputs=[
                    {
                        "prop_id": '{"index":"all","type":"load-file-clear"}.n_clicks',
                        "value": 0,
                    }
                ]
            )
        )
        try:
            result = fn([0], [0], [], [{"id": "e1"}], [])
        finally:
            cc.context_value.set({})
        assert result is no_update

    def test_clear_real_click_empties_selection(self, web_app):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        fn = self._fn(web_app)
        cc.context_value.set(
            AttributeDict(
                triggered_inputs=[
                    {
                        "prop_id": '{"index":"all","type":"load-file-clear"}.n_clicks',
                        "value": 1,
                    }
                ]
            )
        )
        try:
            result = fn([1], [1], [], [{"id": "e1"}], [])
        finally:
            cc.context_value.set({})
        assert result == []

    def test_remove_real_click_filters_entry(self, web_app):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        fn = self._fn(web_app)
        cc.context_value.set(
            AttributeDict(
                triggered_inputs=[
                    {
                        "prop_id": '{"index":"e1","type":"load-file-remove"}.n_clicks',
                        "value": 1,
                    }
                ]
            )
        )
        try:
            result = fn(
                [1],
                [0],
                [],
                [{"id": "e1"}, {"id": "e2"}],
                [],
            )
        finally:
            cc.context_value.set({})
        assert result == [{"id": "e2"}]

    def test_rename_updates_filename(self, web_app):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        fn = self._fn(web_app)
        cc.context_value.set(
            AttributeDict(
                triggered_inputs=[
                    {
                        "prop_id": '{"index":"e1","type":"load-file-name"}.value',
                        "value": "new_name.edi",
                    }
                ]
            )
        )
        try:
            result = fn(
                [],
                [],
                ["new_name.edi"],
                [{"id": "e1", "filename": "old.edi", "original": "old.edi"}],
                [{"index": "e1"}],
            )
        finally:
            cc.context_value.set({})
        assert result[0]["filename"] == "new_name.edi"

    def test_rename_unchanged_returns_no_update(self, web_app):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        fn = self._fn(web_app)
        cc.context_value.set(
            AttributeDict(
                triggered_inputs=[
                    {
                        "prop_id": '{"index":"e1","type":"load-file-name"}.value',
                        "value": "old.edi",
                    }
                ]
            )
        )
        try:
            result = fn(
                [],
                [],
                ["old.edi"],
                [{"id": "e1", "filename": "old.edi", "original": "old.edi"}],
                [{"index": "e1"}],
            )
        finally:
            cc.context_value.set({})
        assert result is no_update

    def test_default_trigger_returns_entries_unchanged(self, web_app):
        fn = self._fn(web_app)
        entries = [{"id": "e1"}]
        _set_triggered("")
        try:
            result = fn([], [], [], entries, [])
        finally:
            _clear_triggered()
        assert result == entries


# ── render_upload_manager ────────────────────────────────────────────────────


class TestRenderUploadManager:
    def _fn(self, web_app):
        return _cb_multi(
            web_app,
            "load-upload-file-manager.children",
            f"{IDs.UPLOAD_FILELIST}.children",
        )

    def test_empty_entries(self, web_app):
        manager, filelist = self._fn(web_app)([])
        assert manager == ""
        assert filelist == ""

    def test_renders_rows(self, web_app):
        entries = [
            {"id": "e1", "filename": "a.edi", "original": "a.edi"},
            {"id": "e2", "filename": "b", "original": "b"},
        ]
        manager, filelist = self._fn(web_app)(entries)
        assert manager != ""
        assert "2 file(s) ready" in filelist


# ── folder filter (populate) ─────────────────────────────────────────────────


class TestFolderFilterPopulate:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.LOAD_FOLDER_FILTER}.options")

    def test_empty_entries_hides_filter(self, web_app):
        opts, value, style = self._fn(web_app)([])
        assert opts == []
        assert style["display"] == "none"

    def test_single_line_hides_filter(self, web_app):
        entries = [{"filename": "L1/a.edi"}, {"filename": "L1/b.edi"}]
        opts, value, style = self._fn(web_app)(entries)
        assert style["display"] == "none"

    def test_multi_line_shows_filter(self, web_app):
        entries = [{"filename": "L1/a.edi"}, {"filename": "L2/b.edi"}]
        opts, value, style = self._fn(web_app)(entries)
        assert style["display"] == "block"
        assert {o["value"] for o in opts} == {"L1", "L2"}
        assert set(value) == {"L1", "L2"}


# ── update_preflight ─────────────────────────────────────────────────────────


class TestUpdatePreflight:
    def _fn(self, web_app):
        return _cb_multi(web_app, "load-preflight-preview.children")

    def test_no_entries(self, web_app):
        preview, summary, source = self._fn(web_app)([], None, "replace", None)
        assert source == "none"

    def test_with_entries(self, web_app):
        entries = [{"filename": "L1/a.edi", "source": "folder"}]
        preview, summary, source = self._fn(web_app)(
            entries, None, "replace", None
        )
        assert source == "folder"


# ── update_table ─────────────────────────────────────────────────────────────


class TestUpdateTable:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.STATION_TABLE}.data")

    def test_no_store_data(self, web_app):
        result = self._fn(web_app)(None)
        assert result[0] == []
        assert result[2] == "No data loaded"

    def test_with_records_full_kpis(self, web_app):
        store_data = {
            "station_records": [
                {
                    "ID": "S1",
                    "N_freq": 20,
                    "Latitude": 1.234,
                    "Line": "L1",
                },
                {
                    "ID": "S2",
                    "N_freq": 30,
                    "Latitude": 1.235,
                    "Line": "L1",
                },
            ],
            "n_stations": 2,
            "n_lines": 1,
            "data_dir": "/some/survey_dir",
        }
        (
            records,
            label,
            status,
            kpi_stations,
            kpi_freq,
            kpi_profiles,
            kpi_survey,
        ) = self._fn(web_app)(store_data)
        assert label == "2 stations loaded"
        assert kpi_stations == "2"
        assert kpi_freq == "25"
        assert kpi_profiles == "1"
        assert kpi_survey == "survey_dir"

    def test_uploaded_survey_label(self, web_app):
        store_data = {
            "station_records": [{"ID": "S1"}],
            "n_stations": 1,
            "data_dir": "[uploaded]",
        }
        result = self._fn(web_app)(store_data)
        assert result[-1] == "upload"

    def test_profiles_fallback_to_lat_clusters(self, web_app):
        store_data = {
            "station_records": [
                {"ID": "S1", "Latitude": 1.0},
                {"ID": "S2", "Latitude": 2.0},
            ],
            "n_stations": 2,
        }
        result = self._fn(web_app)(store_data)
        assert result[5] == "2"  # kpi_profiles


# ── render_station_list ──────────────────────────────────────────────────────


class TestRenderStationList:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.STATION_LINE_PILLS}.children")

    def test_no_store_data(self, web_app):
        pills, rows = self._fn(web_app)(None, None, None, None, "dark")
        assert pills == []

    def test_no_records(self, web_app):
        pills, rows = self._fn(web_app)(
            {"station_records": []}, None, None, None, "dark"
        )
        assert pills == []

    def test_renders_pills_and_rows(self, web_app):
        store_data = {
            "station_records": [
                {
                    "ID": "S1",
                    "Line": "L1",
                    "Latitude": 1.0,
                    "Longitude": 2.0,
                    "Elevation": 100,
                    "N_freq": 20,
                    "Tipper": True,
                },
            ]
        }
        active_lines = {"active": ["L1"], "all": ["L1"]}
        pills, rows = self._fn(web_app)(
            store_data, active_lines, {"station_id": "S1"}, None, "light"
        )
        assert len(pills) == 2  # "All" + "L1"
        assert len(rows) == 1

    def test_filter_excludes_other_lines(self, web_app):
        store_data = {
            "station_records": [
                {"ID": "S1", "Line": "L1"},
                {"ID": "S2", "Line": "L2"},
            ]
        }
        active_lines = {"active": ["L1", "L2"], "all": ["L1", "L2"]}
        pills, rows = self._fn(web_app)(
            store_data, active_lines, None, "L1", "dark"
        )
        assert len(rows) == 1

    def test_filter_matches_nothing_shows_message(self, web_app):
        store_data = {"station_records": [{"ID": "S1", "Line": "L1"}]}
        active_lines = {"active": ["L1"], "all": ["L1"]}
        pills, rows = self._fn(web_app)(
            store_data, active_lines, None, "L2", "dark"
        )
        assert "No stations match" in str(rows)

    def test_bad_numeric_metadata_swallowed(self, web_app):
        store_data = {
            "station_records": [
                {
                    "ID": "S1",
                    "Line": "L1",
                    "Latitude": "not-a-number",
                }
            ]
        }
        pills, rows = self._fn(web_app)(
            store_data, None, None, None, "dark"
        )
        assert len(rows) == 1


# ── toggle_line_filter ───────────────────────────────────────────────────────


class TestToggleLineFilter:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.STATION_LINE_FILTER}.data")

    def test_no_trigger_returns_no_update(self, web_app):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        cc.context_value.set(AttributeDict(triggered_inputs=[]))
        try:
            assert self._fn(web_app)([0], None) is no_update
        finally:
            _clear_triggered()

    def test_all_pill_clears_filter(self, web_app):
        _set_triggered('{"index":"__all__","type":"sta-line-pill"}.n_clicks')
        try:
            result = self._fn(web_app)([1], "L1")
        finally:
            _clear_triggered()
        assert result is None

    def test_line_pill_sets_filter(self, web_app):
        _set_triggered('{"index":"L1","type":"sta-line-pill"}.n_clicks')
        try:
            result = self._fn(web_app)([1], None)
        finally:
            _clear_triggered()
        assert result == "L1"

    def test_line_pill_toggles_off_when_already_active(self, web_app):
        _set_triggered('{"index":"L1","type":"sta-line-pill"}.n_clicks')
        try:
            result = self._fn(web_app)([1], "L1")
        finally:
            _clear_triggered()
        assert result is None


# ── load_data (real EDI end to end) ──────────────────────────────────────────


class TestLoadData:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.STORE_DATA}.data", IDs.BTN_LOAD_CONFIRM
        )

    def test_no_session_id(self, web_app):
        result = self._fn(web_app)(1, [], None, None, "replace", None)
        assert "Session not initialised" in result[1]

    def test_no_usable_entries(self, web_app):
        result = self._fn(web_app)(1, [], None, "sess", "replace", None)
        assert "Choose a survey folder" in result[1]

    def test_upload_success_replace_mode(self, web_app, edi_data_uri):
        entries = [
            {
                "id": "e1",
                "source": "upload",
                "filename": "SIM001.edi",
                "original": "SIM001.edi",
                "content": edi_data_uri,
            }
        ]
        store, feedback, is_open, filelist, selection = self._fn(web_app)(
            1, entries, None, "sess-upload", "replace", None
        )
        assert is_open is False
        assert store["n_stations"] == 1
        assert "Loaded 1 station" in feedback
        assert selection == []

    def test_upload_append_mode_merges(self, web_app, edi_data_uri, monkeypatch):
        import pycsamt.app.web.callbacks.data as data_mod

        entries = [
            {
                "id": "e1",
                "source": "upload",
                "filename": "SIM001.edi",
                "original": "SIM001.edi",
                "content": edi_data_uri,
            }
        ]
        existing_store = {
            "station_records": [{"ID": "OLD1", "Line": "flat files"}],
            "data_dir": "[uploaded]",
        }
        fake_merged_sites = object()
        monkeypatch.setattr(
            data_mod,
            "cache_merge_sites",
            lambda *a, **k: fake_merged_sites,
        )
        store, feedback, is_open, filelist, selection = self._fn(web_app)(
            1, entries, None, "sess-append", "append", existing_store
        )
        assert "Added" in feedback
        assert store["n_stations"] == 2

    def test_upload_no_recognised_files(self, web_app):
        entries = [
            {
                "id": "e1",
                "source": "upload",
                "filename": "readme.txt",
                "original": "readme.txt",
                "content": "data:text/plain;base64,aGVsbG8=",
            }
        ]
        result = self._fn(web_app)(1, entries, None, "sess", "replace", None)
        assert "No recognised files" in result[1]

    def test_upload_loader_exception_reported(
        self, web_app, edi_data_uri, monkeypatch
    ):
        import pycsamt.app.web.callbacks.data as data_mod

        def _boom(self, paths):
            raise RuntimeError("parse failed")

        monkeypatch.setattr(
            data_mod.DataController, "load", _boom
        )
        entries = [
            {
                "id": "e1",
                "source": "upload",
                "filename": "SIM001.edi",
                "original": "SIM001.edi",
                "content": edi_data_uri,
            }
        ]
        result = self._fn(web_app)(1, entries, None, "sess", "replace", None)
        assert "Error loading files" in result[1]

    def test_folder_success_replace_mode(self, web_app, edi_data_uri):
        entries = [
            {
                "id": "e1",
                "source": "folder",
                "filename": "L1/SIM001.edi",
                "original": "L1/SIM001.edi",
                "content": edi_data_uri,
            }
        ]
        store, feedback, is_open, filelist, selection = self._fn(web_app)(
            1, entries, None, "sess-folder", "replace", None
        )
        assert is_open is False
        assert store["n_stations"] == 1
        assert store["n_lines"] == 1
        assert "Loaded 1 station" in feedback

    def test_folder_append_mode(
        self, web_app, edi_data_uri, monkeypatch
    ):
        import pycsamt.app.web.callbacks.data as data_mod

        entries = [
            {
                "id": "e1",
                "source": "folder",
                "filename": "L2/SIM001.edi",
                "original": "L2/SIM001.edi",
                "content": edi_data_uri,
            }
        ]
        existing_store = {
            "station_records": [{"ID": "OLD1", "Line": "L1"}],
            "data_dir": "[browsed]",
        }
        monkeypatch.setattr(
            data_mod, "cache_merge_sites", lambda *a, **k: object()
        )
        store, feedback, is_open, filelist, selection = self._fn(web_app)(
            1, entries, None, "sess-folder-append", "append", existing_store
        )
        assert "Added" in feedback
        assert store["n_stations"] == 2

    def test_folder_no_filenames(self, web_app):
        entries = [
            {
                "id": "e1",
                "source": "folder",
                "filename": "",
                "original": "",
                "content": "",
            }
        ]
        # content is falsy so `usable` filters it out entirely -> "none" path
        result = self._fn(web_app)(1, entries, None, "sess", "replace", None)
        assert "Choose a survey folder" in result[1]

    def test_folder_loader_exception_reported(
        self, web_app, edi_data_uri, monkeypatch
    ):
        import pycsamt.app.web.callbacks.data as data_mod

        def _boom(self, paths, path_to_line=None):
            raise RuntimeError("folder parse failed")

        monkeypatch.setattr(data_mod.DataController, "load", _boom)
        entries = [
            {
                "id": "e1",
                "source": "folder",
                "filename": "L1/SIM001.edi",
                "original": "L1/SIM001.edi",
                "content": edi_data_uri,
            }
        ]
        result = self._fn(web_app)(1, entries, None, "sess", "replace", None)
        assert "Error" in result[1]
