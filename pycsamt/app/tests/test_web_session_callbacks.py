# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.web.callbacks._session (session persistence)."""

from __future__ import annotations

import base64
import json

from pycsamt.app.web.callbacks._session import _validate_snapshot
from pycsamt.app.web.layout import IDs


def _unwrap(entry):
    """Return the plain undecorated callback function from a callback_map entry."""
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop: str):
    """Return the raw callback function registered for *output_id_prop*."""
    return _unwrap(web_app.callback_map[output_id_prop])


# ── pure helpers ─────────────────────────────────────────────────────────────


def test_validate_snapshot_accepts_well_formed_snapshot():
    _validate_snapshot({"version": "2.0", "store_data": {"n_stations": 3}})


def test_validate_snapshot_rejects_non_dict():
    import pytest

    with pytest.raises(ValueError):
        _validate_snapshot("not-a-dict")


def test_validate_snapshot_rejects_unsupported_version():
    import pytest

    with pytest.raises(ValueError, match="Unsupported session version"):
        _validate_snapshot({"version": "1.0", "store_data": {}})


def test_validate_snapshot_rejects_missing_store_data():
    import pytest

    with pytest.raises(ValueError, match="Missing 'store_data'"):
        _validate_snapshot({"version": "2.0"})


# ── _init_session_id ─────────────────────────────────────────────────────────


def test_init_session_id_generates_uuid_when_absent(web_app):
    fn = _cb(web_app, f"{IDs.SESSION_ID}.data")
    result = fn(None)
    assert isinstance(result, str) and len(result) == 36


def test_init_session_id_no_update_when_present(web_app):
    from dash import no_update

    fn = _cb(web_app, f"{IDs.SESSION_ID}.data")
    result = fn("already-set")
    assert result is no_update


# ── _toggle_session_canvas ───────────────────────────────────────────────────


def test_toggle_session_canvas_flips_state(web_app):
    fn = _cb(web_app, f"{IDs.SESSION_CANVAS}.is_open")
    assert fn(1, False) is True
    assert fn(1, True) is False


def test_toggle_session_canvas_no_update_without_click(web_app):
    from dash import no_update

    fn = _cb(web_app, f"{IDs.SESSION_CANVAS}.is_open")
    assert fn(0, False) is no_update


# ── _auto_snapshot ───────────────────────────────────────────────────────────


def _auto_snapshot_fn(web_app):
    matches = [k for k in web_app.callback_map if "session-autosave" in k]
    assert matches, "auto-snapshot callback not found"
    return _unwrap(web_app.callback_map[matches[0]])


def test_auto_snapshot_no_update_without_store_data(web_app):
    from dash import no_update

    fn = _auto_snapshot_fn(web_app)
    snap, chip = fn(None, None, None, None, None, None)
    assert snap is no_update
    assert chip is no_update


def test_auto_snapshot_builds_snapshot_and_chip(web_app):
    fn = _auto_snapshot_fn(web_app)
    store_data = {"n_stations": 4}
    snap, chip = fn(store_data, None, None, ["L1"], {"steps": []}, "note")
    assert snap["store_data"] == store_data
    assert snap["note"] == "note"
    assert snap["active_lines"] == ["L1"]
    assert snap["version"] == "2.0"
    assert chip is not None


# ── _session_download ────────────────────────────────────────────────────────


def _download_fn(web_app):
    matches = [k for k in web_app.callback_map if "session-download" in k]
    assert matches, "download callback not found"
    return _unwrap(web_app.callback_map[matches[0]])


def test_session_download_no_update_without_click(web_app):
    from dash import no_update

    fn = _download_fn(web_app)
    dl, feedback = fn(0, None, None, None, None, None, None)
    assert dl is no_update
    assert feedback is no_update


def test_session_download_errors_without_store_data(web_app):
    fn = _download_fn(web_app)
    dl, feedback = fn(1, None, None, None, None, None, None)
    from dash import no_update

    assert dl is no_update
    assert feedback is not None


def test_session_download_returns_json_payload(web_app):
    fn = _download_fn(web_app)
    store_data = {"n_stations": 2}
    dl, feedback = fn(1, store_data, ["L1"], {"steps": []}, None, None, "n")
    assert dl["filename"].startswith("pycsamt_session_")
    payload = json.loads(dl["content"])
    assert payload["store_data"] == store_data
    assert feedback is not None


# ── _session_upload_restore ──────────────────────────────────────────────────


def _upload_restore_fn(web_app):
    matches = [
        k
        for k in web_app.callback_map
        if IDs.SESSION_UL not in k and "store-data" in k and "corr-state" in k
    ]
    # narrow further: the upload-restore callback is keyed by the 6 outputs
    matches = [
        k
        for k in web_app.callback_map
        if k.count("store-data") == 1
        and "tool-elev-raw-store" in k
        and "session-feedback" in k
        and web_app.callback_map[k]["inputs"]
        and web_app.callback_map[k]["inputs"][0]["id"] == IDs.SESSION_UL
    ]
    assert len(matches) == 1, matches
    return _unwrap(web_app.callback_map[matches[0]])


def test_session_upload_restore_no_update_without_contents(web_app):
    fn = _upload_restore_fn(web_app)
    result = fn(None, "session.json")
    assert all(r is result[0] for r in result)


def test_session_upload_restore_success(web_app):
    fn = _upload_restore_fn(web_app)
    snap = {
        "version": "2.0",
        "store_data": {"n_stations": 7},
        "active_lines": ["L1"],
        "corr_store": {"steps": []},
        "elevation_raw": None,
        "elevation_corrected": None,
    }
    encoded = base64.b64encode(json.dumps(snap).encode()).decode()
    contents = f"data:application/json;base64,{encoded}"

    result = fn(contents, "session.json")
    store_data, active_lines, corr_store, elev_raw, elev_corr, feedback = (
        result
    )
    assert store_data == snap["store_data"]
    assert active_lines == ["L1"]
    assert feedback is not None


def test_session_upload_restore_handles_bad_payload(web_app):
    fn = _upload_restore_fn(web_app)
    contents = "data:application/json;base64," + base64.b64encode(
        b"not valid json"
    ).decode()
    result = fn(contents, "broken.json")
    from dash import no_update

    assert result[0] is no_update
    assert result[-1] is not None


# ── _session_browser_restore ─────────────────────────────────────────────────


def _browser_restore_fn(web_app):
    matches = [
        k
        for k in web_app.callback_map
        if k.count("store-data") == 1
        and "tool-elev-raw-store" in k
        and "session-feedback" in k
        and web_app.callback_map[k]["inputs"]
        and web_app.callback_map[k]["inputs"][0]["id"]
        == IDs.BTN_SESSION_RESTORE
    ]
    assert len(matches) == 1, matches
    return _unwrap(web_app.callback_map[matches[0]])


def test_session_browser_restore_no_update_without_click(web_app):
    fn = _browser_restore_fn(web_app)
    result = fn(0, {"version": "2.0", "store_data": {}})
    assert all(r is result[0] for r in result)


def test_session_browser_restore_errors_without_snapshot(web_app):
    fn = _browser_restore_fn(web_app)
    result = fn(1, None)
    from dash import no_update

    assert result[0] is no_update
    assert result[-1] is not None


def test_session_browser_restore_success(web_app):
    fn = _browser_restore_fn(web_app)
    snap = {
        "version": "2.0",
        "saved_at": "2026-01-01T00:00:00+00:00",
        "note": "hello",
        "store_data": {"n_stations": 5},
        "active_lines": None,
        "corr_store": None,
        "elevation_raw": None,
        "elevation_corrected": None,
    }
    result = fn(1, snap)
    store_data = result[0]
    corr_store = result[2]
    feedback = result[-1]
    assert store_data == snap["store_data"]
    assert corr_store == {"steps": []}
    assert feedback is not None


# ── _session_clear ───────────────────────────────────────────────────────────


def _clear_fn(web_app):
    matches = [k for k in web_app.callback_map if "btn-session-clear" in str(
        web_app.callback_map[k].get("inputs")
    )]
    assert len(matches) == 1, matches
    return _unwrap(web_app.callback_map[matches[0]])


def test_session_clear_no_update_without_click(web_app):
    fn = _clear_fn(web_app)
    result = fn(0)
    assert all(r is result[0] for r in result)


def test_session_clear_resets_snapshot(web_app):
    fn = _clear_fn(web_app)
    snapshot, autosave, feedback = fn(1)
    assert snapshot is None
    assert autosave == ""
    assert feedback is not None
