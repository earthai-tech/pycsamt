# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.help."""

from __future__ import annotations

import json

import pytest

pytest.importorskip("dash", reason="dash required")

from dash._utils import AttributeDict
from dash.exceptions import PreventUpdate

from pycsamt.app.agent_master._ids import IDs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _toggle_help_fn(agent_app):
    return _unwrap(agent_app.callback_map[f"{IDs.MODAL_HELP}.is_open"])


def _use_example_fn(agent_app):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"] and IDs.HELP_CHIP in entry["inputs"][0]["id"]
    ]
    assert len(matches) == 1, matches
    return _unwrap(agent_app.callback_map[matches[0]])


class TestToggleHelp:
    def test_raises_prevent_update_without_any_click(self, agent_app):
        fn = _toggle_help_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, None, False)

    def test_flips_open_state_on_open_click(self, agent_app):
        fn = _toggle_help_fn(agent_app)
        assert fn(1, None, False) is True
        assert fn(1, None, True) is False

    def test_flips_open_state_on_close_click(self, agent_app):
        fn = _toggle_help_fn(agent_app)
        assert fn(None, 1, True) is False


class TestUseExample:
    def _set_trigger(self, index):
        import dash._callback_context as cc

        prop_id = (
            json.dumps({"index": index, "type": IDs.HELP_CHIP}, sort_keys=True)
            + ".n_clicks"
        )
        cc.context_value.set(
            AttributeDict(
                triggered_inputs=[{"prop_id": prop_id, "value": 1}]
            )
        )

    def teardown_method(self, _method):
        import dash._callback_context as cc

        cc.context_value.set({})

    def test_raises_prevent_update_without_trigger(self, agent_app):
        import dash._callback_context as cc

        cc.context_value.set(AttributeDict(triggered_inputs=[]))
        fn = _use_example_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([None])

    def test_raises_prevent_update_for_out_of_range_index(self, agent_app):
        self._set_trigger(999)
        fn = _use_example_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([1])

    def test_fills_input_and_closes_modal(self, agent_app):
        from pycsamt.app.agent_master.layout import _HELP_EXAMPLES

        self._set_trigger(0)
        fn = _use_example_fn(agent_app)
        text, is_open = fn([1])
        assert text == _HELP_EXAMPLES[0]
        assert is_open is False
