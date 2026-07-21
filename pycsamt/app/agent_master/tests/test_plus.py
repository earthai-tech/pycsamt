# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.plus."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")

from dash.exceptions import PreventUpdate

from pycsamt.app.agent_master._ids import IDs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _find_by_input(agent_app, input_id):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"] and entry["inputs"][0]["id"] == input_id
    ]
    assert len(matches) == 1, (input_id, matches)
    return _unwrap(agent_app.callback_map[matches[0]])


class TestOpenEdi:
    def test_prevents_update_without_click(self, agent_app):
        fn = _find_by_input(agent_app, IDs.PLUS_LOAD)
        with pytest.raises(PreventUpdate):
            fn(None)

    def test_opens_on_click(self, agent_app):
        fn = _find_by_input(agent_app, IDs.PLUS_LOAD)
        assert fn(1) is True


class TestOpenSettings:
    def test_prevents_update_without_click(self, agent_app):
        fn = _find_by_input(agent_app, IDs.PLUS_SETTINGS)
        with pytest.raises(PreventUpdate):
            fn(None)

    def test_opens_on_click(self, agent_app):
        fn = _find_by_input(agent_app, IDs.PLUS_SETTINGS)
        assert fn(1) is True


class TestPastePath:
    def test_prevents_update_without_click(self, agent_app):
        fn = _find_by_input(agent_app, IDs.PLUS_PATH)
        with pytest.raises(PreventUpdate):
            fn(None)

    def test_returns_placeholder_text(self, agent_app):
        fn = _find_by_input(agent_app, IDs.PLUS_PATH)
        assert fn(1) == "Load /path/to/your/EDIs, then "


class TestInjectWebLaunch:
    def test_prevents_update_without_click(self, agent_app):
        fn = _find_by_input(agent_app, IDs.PLUS_WEB)
        with pytest.raises(PreventUpdate):
            fn(None)

    def test_returns_open_web_app_text(self, agent_app):
        fn = _find_by_input(agent_app, IDs.PLUS_WEB)
        assert fn(1) == "open web app"
