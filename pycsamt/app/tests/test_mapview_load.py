# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.load."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestModalTitleAndHint:
    def test_append_mode_title_and_hint(self):
        from pycsamt.app.mapview.callbacks.load import (
            _mode_hint,
            _modal_title,
        )

        title = _modal_title("append")
        assert "Add Lines to Survey" in title
        assert "merged" in _mode_hint("append")

    def test_replace_mode_title_and_hint(self):
        from pycsamt.app.mapview.callbacks.load import (
            _mode_hint,
            _modal_title,
        )

        title = _modal_title("replace")
        assert "Load Survey Lines" in title
        assert "replaced" in _mode_hint("replace")


class TestSanitize:
    def test_empty_name_uses_fallback(self):
        from pycsamt.app.mapview.callbacks.load import _sanitize

        assert _sanitize("", "fallback.edi") == "fallback.edi"
        assert _sanitize(None, "fallback.edi") == "fallback.edi"

    def test_replaces_illegal_characters(self):
        from pycsamt.app.mapview.callbacks.load import _sanitize

        assert _sanitize("bad:name?.edi", "fb") == "bad_name_.edi"

    def test_preserves_nested_path_structure(self):
        from pycsamt.app.mapview.callbacks.load import _sanitize

        assert _sanitize("L1/site001.edi", "fb") == "L1/site001.edi"

    def test_backslashes_normalized_to_forward(self):
        from pycsamt.app.mapview.callbacks.load import _sanitize

        assert _sanitize("L1\\site001.edi", "fb") == "L1/site001.edi"


class TestEntries:
    def test_none_contents_returns_empty(self):
        from pycsamt.app.mapview.callbacks.load import _entries

        assert _entries(None, None, source="upload") == []

    def test_single_string_contents_wrapped_in_list(self):
        from pycsamt.app.mapview.callbacks.load import _entries

        out = _entries("data:...", "a.edi", source="upload")
        assert len(out) == 1
        assert out[0]["filename"] == "a.edi"
        assert out[0]["source"] == "upload"

    def test_non_edi_extension_is_skipped(self):
        from pycsamt.app.mapview.callbacks.load import _entries

        out = _entries(["data:..."], ["notes.txt"], source="upload")
        assert out == []

    def test_accepts_avg_and_j_extensions(self):
        from pycsamt.app.mapview.callbacks.load import _entries

        out = _entries(
            ["data:...", "data:..."], ["a.avg", "b.j"], source="folder"
        )
        assert len(out) == 2

    def test_missing_filename_gets_default(self):
        from pycsamt.app.mapview.callbacks.load import _entries

        out = _entries(["data:..."], [], source="upload")
        assert out[0]["original"] == "file_1.edi"


class TestInferLinesAndPath:
    def test_line_from_flat_path(self):
        from pycsamt.app.mapview.callbacks.load import _line_from_path

        assert _line_from_path("site001.edi") == "flat files"

    def test_line_from_single_folder(self):
        from pycsamt.app.mapview.callbacks.load import _line_from_path

        assert _line_from_path("L1/site001.edi") == "L1"

    def test_line_from_nested_folder_uses_second_part(self):
        from pycsamt.app.mapview.callbacks.load import _line_from_path

        assert _line_from_path("survey/L1/site001.edi") == "L1"

    def test_infer_lines_counts_by_group(self):
        from pycsamt.app.mapview.callbacks.load import _infer_lines

        counts = _infer_lines(
            ["L1/a.edi", "L1/b.edi", "L2/c.edi", "d.edi"]
        )
        assert counts == {"L1": 2, "L2": 1, "flat files": 1}

    def test_entry_line_prefers_original_over_filename(self):
        from pycsamt.app.mapview.callbacks.load import _entry_line

        entry = {"original": "L1/a.edi", "filename": "a.edi"}
        assert _entry_line(entry) == "L1"


class TestFilteredEntries:
    def test_no_selection_returns_all(self):
        from pycsamt.app.mapview.callbacks.load import _filtered_entries

        entries = [{"original": "L1/a.edi"}, {"original": "L2/b.edi"}]
        assert _filtered_entries(entries, None) == entries
        assert _filtered_entries(entries, []) == entries

    def test_filters_to_selected_lines(self):
        from pycsamt.app.mapview.callbacks.load import _filtered_entries

        entries = [{"original": "L1/a.edi"}, {"original": "L2/b.edi"}]
        result = _filtered_entries(entries, ["L1"])
        assert result == [{"original": "L1/a.edi"}]


class TestBuildView:
    def test_returns_view_on_first_success(self, monkeypatch):
        import pycsamt.app.mapview.callbacks.load as load_mod

        calls = {"n": 0}

        class _FakeMapView:
            @classmethod
            def from_lines(cls, line_map, theme):
                calls["n"] += 1
                return "built-view"

        monkeypatch.setattr(
            "pycsamt.map.MapView",
            _FakeMapView,
        )
        result = load_mod._build_view({"L1": ["a.edi"]}, "light")
        assert result == "built-view"
        assert calls["n"] == 1

    def test_retries_on_permission_error_then_succeeds(self, monkeypatch):
        import pycsamt.app.mapview.callbacks.load as load_mod

        calls = {"n": 0}

        class _FlakyMapView:
            @classmethod
            def from_lines(cls, line_map, theme):
                calls["n"] += 1
                if calls["n"] < 3:
                    raise PermissionError("locked")
                return "built-view"

        monkeypatch.setattr("pycsamt.map.MapView", _FlakyMapView)
        result = load_mod._build_view(
            {"L1": ["a.edi"]}, "light", attempts=4, delay=0.0
        )
        assert result == "built-view"
        assert calls["n"] == 3

    def test_raises_last_exception_after_all_attempts(self, monkeypatch):
        import pycsamt.app.mapview.callbacks.load as load_mod

        class _AlwaysLocked:
            @classmethod
            def from_lines(cls, line_map, theme):
                raise OSError("still locked")

        monkeypatch.setattr("pycsamt.map.MapView", _AlwaysLocked)
        with pytest.raises(OSError, match="still locked"):
            load_mod._build_view(
                {"L1": ["a.edi"]}, "light", attempts=2, delay=0.0
            )


class TestDecode:
    def test_upload_source_uses_decode_upload_to_tempdir(self, monkeypatch):
        import pycsamt.app.mapview.callbacks.load as load_mod
        import pycsamt.app.web.utils as web_utils

        monkeypatch.setattr(
            web_utils,
            "decode_upload_to_tempdir",
            lambda contents, filenames: (["a.edi", "b.edi"], "/tmp/x"),
        )
        line_map, tmpdir = load_mod._decode(
            [{"content": "c1", "filename": "a.edi"}], "upload"
        )
        assert line_map == {"uploaded": ["a.edi", "b.edi"]}
        assert tmpdir == "/tmp/x"

    def test_upload_source_empty_paths_returns_empty_map(self, monkeypatch):
        import pycsamt.app.mapview.callbacks.load as load_mod
        import pycsamt.app.web.utils as web_utils

        monkeypatch.setattr(
            web_utils,
            "decode_upload_to_tempdir",
            lambda contents, filenames: ([], "/tmp/y"),
        )
        line_map, tmpdir = load_mod._decode(
            [{"content": "c1", "filename": "a.edi"}], "upload"
        )
        assert line_map == {}
        assert tmpdir == "/tmp/y"

    def test_folder_source_uses_decode_folder_upload_to_tempdir(
        self, monkeypatch
    ):
        import pycsamt.app.mapview.callbacks.load as load_mod
        import pycsamt.app.web.utils as web_utils

        monkeypatch.setattr(
            web_utils,
            "decode_folder_upload_to_tempdir",
            lambda contents, filenames: (
                ["a.edi"],
                "/tmp/z",
                {"L1": ["a.edi"]},
            ),
        )
        line_map, tmpdir = load_mod._decode(
            [{"content": "c1", "filename": "L1/a.edi"}], "folder"
        )
        assert line_map == {"L1": ["a.edi"]}
        assert tmpdir == "/tmp/z"


class TestRegisterLoad:
    def test_register_load_is_callable(self):
        from pycsamt.app.mapview.callbacks.load import register_load

        assert callable(register_load)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.MODAL_LOAD in cb_outputs
        assert IDs.UPLOAD_SELECTION in cb_outputs
        assert IDs.BTN_LOAD_CONFIRM in cb_outputs
