# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.session."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestNowIso:
    def test_returns_iso_string_with_timezone(self):
        from pycsamt.app.mapview.callbacks.session import _now_iso

        stamp = _now_iso()
        assert "T" in stamp
        assert stamp.endswith("+00:00")


class TestIcons:
    def test_ok_icon_contains_message(self):
        from pycsamt.app.mapview.callbacks.session import _ok_icon

        span = _ok_icon("Saved!")
        assert "Saved!" in str(span)
        assert "bi-check-circle-fill" in str(span)

    def test_err_icon_contains_message(self):
        from pycsamt.app.mapview.callbacks.session import _err_icon

        span = _err_icon("Failed!")
        assert "Failed!" in str(span)
        assert "bi-x-circle-fill" in str(span)


class TestBuildSnapshot:
    def test_contains_all_expected_keys(self):
        from pycsamt.app.mapview.callbacks.session import _build_snapshot

        snap = _build_snapshot(
            "map",
            {"overlay": "index"},
            {"all": ["L1"]},
            None,
            ["S00"],
            "dark",
            {"n_stations": 3},
            "a note",
        )
        assert snap["version"] == "1.0"
        assert snap["app"] == "mapview"
        assert snap["view"] == "map"
        assert snap["controls"] == {"overlay": "index"}
        assert snap["masked"] == ["S00"]
        assert snap["theme"] == "dark"
        assert snap["note"] == "a note"
        assert "saved_at" in snap

    def test_missing_note_defaults_to_empty_string(self):
        from pycsamt.app.mapview.callbacks.session import _build_snapshot

        snap = _build_snapshot(
            "map", {}, {}, None, [], "light", {}, None
        )
        assert snap["note"] == ""


class TestValidateSnapshot:
    def test_valid_snapshot_passes(self):
        from pycsamt.app.mapview.callbacks.session import _validate_snapshot

        _validate_snapshot({"version": "1.0", "view": "map"})

    def test_non_dict_raises(self):
        from pycsamt.app.mapview.callbacks.session import _validate_snapshot

        with pytest.raises(ValueError, match="Not a valid session file"):
            _validate_snapshot("not-a-dict")

    def test_unsupported_version_raises(self):
        from pycsamt.app.mapview.callbacks.session import _validate_snapshot

        with pytest.raises(ValueError, match="Unsupported session version"):
            _validate_snapshot({"version": "2.0", "view": "map"})

    def test_missing_view_key_raises(self):
        from pycsamt.app.mapview.callbacks.session import _validate_snapshot

        with pytest.raises(ValueError, match="Missing 'view' key"):
            _validate_snapshot({"version": "1.0"})


class TestRegisterSession:
    def test_register_session_is_callable(self):
        from pycsamt.app.mapview.callbacks.session import register_session

        assert callable(register_session)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.CANVAS_SESSION in cb_outputs
        assert IDs.SESSION_SNAPSHOT in cb_outputs
        assert IDs.SESSION_DL in cb_outputs
