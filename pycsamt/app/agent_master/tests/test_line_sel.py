# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Registration smoke tests for pycsamt.app.agent_master.callbacks.line_sel.

The callback bodies delegate to chat.py primitives (_new_job, _run_agent,
_quick_workflow, ...) which already have dedicated coverage in
test_chat_ui.py / test_dispatch_precedence.py; here we just confirm the
module wires up correctly.
"""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestRegisterLineSel:
    def test_register_line_sel_is_callable(self):
        from pycsamt.app.agent_master.callbacks.line_sel import (
            register_line_sel,
        )

        assert callable(register_line_sel)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.MODAL_LINE_SEL in cb_outputs
        assert IDs.LINE_SEL_LIST in cb_outputs
        assert IDs.BTN_LINE_RUN_SEL in cb_outputs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _open_line_modal_fn(agent_app):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"]
        and entry["inputs"][0]["id"] == "am-store-pending"
        and "am-modal-line-sel.is_open" in k
    ]
    assert len(matches) == 1, matches
    return _unwrap(agent_app.callback_map[matches[0]])


def _run_with_lines_fn(agent_app):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"]
        and entry["inputs"][0]["id"] == "am-btn-line-run-sel"
    ]
    assert len(matches) == 1, matches
    return _unwrap(agent_app.callback_map[matches[0]])


class TestOpenLineModal:
    def test_prevent_update_without_pending(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = _open_line_modal_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None)

    def test_prevent_update_for_other_disambiguation(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = _open_line_modal_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn({"disambiguation": "workflow"})

    def test_opens_with_line_options(self, agent_app):
        fn = _open_line_modal_fn(agent_app)
        is_open, opts, value = fn(
            {
                "disambiguation": "lines",
                "groups": {"L1": [1, 2], "L2": [3]},
            }
        )
        assert is_open is True
        assert opts == [
            {"label": "L1", "value": "L1"},
            {"label": "L2", "value": "L2"},
        ]
        assert value == []


class TestRunWithLines:
    def _set_no_trigger(self):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        cc.context_value.set(AttributeDict(triggered_inputs=[]))

    def _set_trigger(self, component_id):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        cc.context_value.set(
            AttributeDict(
                triggered_inputs=[
                    {"prop_id": f"{component_id}.n_clicks", "value": 1}
                ]
            )
        )

    def teardown_method(self, _method):
        import dash._callback_context as cc

        cc.context_value.set({})

    def test_prevent_update_without_trigger(self, agent_app):
        from dash.exceptions import PreventUpdate

        self._set_no_trigger()
        fn = _run_with_lines_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, None, [], {}, {}, {}, {}, [], [])

    def test_run_sel_without_selection_asks_to_select(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        self._set_trigger(IDs.BTN_LINE_RUN_SEL)
        fn = _run_with_lines_fn(agent_app)
        result = fn(
            1, None, [], {"text": "run it", "groups": {"L1": []}},
            {}, {}, {}, [], [],
        )
        from dash import no_update

        assert result[0] is no_update
        assert result[5] is True
        assert "Select at least one line" in str(result[6])

    def test_run_all_spawns_job_when_no_quick_workflow(
        self, agent_app, monkeypatch
    ):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.callbacks import line_sel as line_sel_mod

        monkeypatch.setattr(
            line_sel_mod, "_quick_workflow", lambda text: None
        )
        monkeypatch.setattr(
            line_sel_mod, "_new_job", lambda: "job-123"
        )
        monkeypatch.setattr(
            line_sel_mod, "_drop_workflow", lambda cfg: cfg
        )
        started = {}

        class _FakeThread:
            def __init__(self, target, args, daemon):
                started["target"] = target
                started["args"] = args

            def start(self):
                started["started"] = True

        monkeypatch.setattr(
            line_sel_mod.threading, "Thread", _FakeThread
        )

        self._set_trigger(IDs.BTN_LINE_RUN_ALL)
        fn = _run_with_lines_fn(agent_app)
        result = fn(
            None, 1, [], {"text": "run it", "groups": {"L1": [], "L2": []}},
            {}, {}, {}, [], [],
        )

        assert started.get("started") is True
        msgs, job, disabled, new_stored, pending, is_open, status = result
        assert job == {"jid": "job-123"}
        assert disabled is False
        assert pending == {}
        assert is_open is False

    def test_run_all_returns_waiting_bubble_when_quick_workflow_matches(
        self, agent_app, monkeypatch
    ):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.callbacks import line_sel as line_sel_mod

        monkeypatch.setattr(
            line_sel_mod, "_quick_workflow", lambda text: "ai_inversion"
        )

        self._set_trigger(IDs.BTN_LINE_RUN_ALL)
        fn = _run_with_lines_fn(agent_app)
        result = fn(
            None, 1, [], {"text": "run it", "groups": {"L1": []}},
            {}, {}, {}, [], [],
        )
        from dash import no_update

        msgs, job, disabled, new_stored, pending, is_open, status = result
        assert job is no_update
        assert pending["workflow"] == "ai_inversion"
        assert is_open is False

    def test_run_replaces_line_waiting_bubble_with_quick_workflow(
        self, agent_app, monkeypatch
    ):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.callbacks import line_sel as line_sel_mod

        monkeypatch.setattr(
            line_sel_mod, "_quick_workflow", lambda text: "ai_inversion"
        )

        waiting_bubble = {"props": {"id": "am-line-waiting-bubble"}}
        self._set_trigger(IDs.BTN_LINE_RUN_ALL)
        fn = _run_with_lines_fn(agent_app)
        result = fn(
            None,
            1,
            [],
            {"text": "run it", "groups": {"L1": []}},
            {},
            {},
            {},
            [],
            [waiting_bubble],
        )
        from dash import no_update

        msgs, job, disabled, new_stored, pending, is_open, status = result
        assert job is no_update
        assert pending["workflow"] == "ai_inversion"
        assert is_open is False
        assert msgs[0] != waiting_bubble  # replaced in place

    def test_run_replaces_line_waiting_bubble_without_quick_workflow(
        self, agent_app, monkeypatch
    ):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.callbacks import line_sel as line_sel_mod

        monkeypatch.setattr(
            line_sel_mod, "_quick_workflow", lambda text: None
        )
        monkeypatch.setattr(line_sel_mod, "_new_job", lambda: "job-456")
        monkeypatch.setattr(
            line_sel_mod, "_drop_workflow", lambda cfg: cfg
        )

        class _FakeThread:
            def __init__(self, target, args, daemon):
                pass

            def start(self):
                pass

        monkeypatch.setattr(line_sel_mod.threading, "Thread", _FakeThread)

        waiting_bubble = {"props": {"id": "am-line-waiting-bubble"}}
        self._set_trigger(IDs.BTN_LINE_RUN_ALL)
        fn = _run_with_lines_fn(agent_app)
        result = fn(
            None,
            1,
            [],
            {"text": "run it", "groups": {"L1": []}},
            {},
            {},
            {},
            [],
            [waiting_bubble],
        )

        msgs, job, disabled, new_stored, pending, is_open, status = result
        assert job == {"jid": "job-456"}
        assert msgs[0] != waiting_bubble  # replaced with thinking bubble
