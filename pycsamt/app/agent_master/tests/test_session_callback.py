# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.session (save_session)."""

from __future__ import annotations

import json

import pytest

pytest.importorskip("dash", reason="dash required")

from dash.exceptions import PreventUpdate

from pycsamt.app.agent_master._ids import IDs
from pycsamt.app.agent_master.callbacks import session as session_mod


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _save_session_fn(agent_app):
    return _unwrap(agent_app.callback_map[f"{IDs.SAVE_STATUS}.children"])


def test_save_session_prevents_update_without_click(agent_app):
    fn = _save_session_fn(agent_app)
    with pytest.raises(PreventUpdate):
        fn(None, [], {}, {})


def test_save_session_writes_json_and_returns_status(agent_app, tmp_path, monkeypatch):
    monkeypatch.setattr(session_mod, "_SESSION_DIR", tmp_path)

    messages = [{"role": "user", "content": "hi"}]
    edi_store = {"path": "/data"}
    settings = {"provider": "offline", "key_openai": "secret"}

    fn = _save_session_fn(agent_app)
    result = fn(1, messages, edi_store, settings)

    saved = list(tmp_path.glob("session_*.json"))
    assert len(saved) == 1
    payload = json.loads(saved[0].read_text())
    assert payload["edi"] == edi_store
    assert payload["message_count"] == 1
    # settings keys prefixed "key_" (secrets) are stripped
    assert "key_openai" not in payload["settings"]
    assert payload["settings"]["provider"] == "offline"

    assert "Saved to" in str(result)


def test_session_path_creates_directory(tmp_path, monkeypatch):
    target = tmp_path / "nested" / "sessions"
    monkeypatch.setattr(session_mod, "_SESSION_DIR", target)

    path = session_mod._session_path()
    assert target.exists()
    assert path.parent == target
    assert path.name.startswith("session_")


def test_session_path_with_explicit_name(tmp_path, monkeypatch):
    monkeypatch.setattr(session_mod, "_SESSION_DIR", tmp_path)
    path = session_mod._session_path("mysession.json")
    assert path == tmp_path / "mysession.json"
