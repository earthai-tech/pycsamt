# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.splash."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")

from dash.exceptions import PreventUpdate

from pycsamt.app.agent_master._ids import IDs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _card_load_opens_edi_fn(agent_app):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"] and entry["inputs"][0]["id"] == IDs.SPLASH_CARD_LOAD
    ]
    assert len(matches) == 1, matches
    return _unwrap(agent_app.callback_map[matches[0]])


def test_card_load_prevents_update_without_click(agent_app):
    fn = _card_load_opens_edi_fn(agent_app)
    with pytest.raises(PreventUpdate):
        fn(None)


def test_card_load_opens_edi_canvas_on_click(agent_app):
    fn = _card_load_opens_edi_fn(agent_app)
    assert fn(1) is True
