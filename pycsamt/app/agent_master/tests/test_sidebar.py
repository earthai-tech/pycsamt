# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.sidebar."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


class TestRestoreBubble:
    def test_user_message_renders_user_row(self):
        from pycsamt.app.agent_master.callbacks.sidebar import (
            _restore_bubble,
        )

        div = _restore_bubble(
            {"role": "user", "content": "hi there", "ts": "10:00"}
        )
        text = str(div)
        assert "am-msg-row user" in div.className
        assert "hi there" in text
        assert "10:00" in text

    def test_assistant_message_renders_markdown_body(self):
        from pycsamt.app.agent_master.callbacks.sidebar import (
            _restore_bubble,
        )

        div = _restore_bubble(
            {"role": "assistant", "content": "**done**", "ts": "10:01"}
        )
        assert div.className == "am-msg-row"
        assert "am-bubble agent" in str(div)

    def test_assistant_message_empty_content_shows_placeholder(self):
        from pycsamt.app.agent_master.callbacks.sidebar import (
            _restore_bubble,
        )

        div = _restore_bubble({"role": "assistant", "content": "  "})
        assert "(no content)" in str(div)


class TestRegisterSidebar:
    def test_register_sidebar_is_callable(self):
        from pycsamt.app.agent_master.callbacks.sidebar import (
            register_sidebar,
        )

        assert callable(register_sidebar)

    def test_expected_outputs_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.STORE_HISTORY in cb_outputs
        assert IDs.SIDEBAR_HISTORY in cb_outputs
        assert IDs.SIDEBAR_FIGS in cb_outputs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _find(agent_app, input_id, output_hint):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"]
        and entry["inputs"][0]["id"] == input_id
        and output_hint in k
    ]
    assert len(matches) == 1, (input_id, output_hint, matches)
    return _unwrap(agent_app.callback_map[matches[0]])


def _set_pattern_trigger(component_id_dict, prop="n_clicks"):
    import json

    import dash._callback_context as cc
    from dash._utils import AttributeDict

    prop_id = json.dumps(component_id_dict, sort_keys=True) + f".{prop}"
    cc.context_value.set(
        AttributeDict(triggered_inputs=[{"prop_id": prop_id, "value": 1}])
    )


def _clear_trigger():
    import dash._callback_context as cc

    cc.context_value.set({})


class TestNewChat:
    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = _find(agent_app, "am-btn-new-chat", "am-chat-window.children")
        with pytest.raises(PreventUpdate):
            fn(None, None)

    def test_resets_stores_and_returns_welcome(self, agent_app, monkeypatch):
        from pycsamt.app.agent_master.callbacks import sidebar as sidebar_mod

        called = {}
        monkeypatch.setattr(
            sidebar_mod,
            "_chat_welcome",
            lambda: "welcome",
        )
        fn = _find(agent_app, "am-btn-new-chat", "am-chat-window.children")
        messages, edi, figs, job, children = fn(1, None)
        assert messages == []
        assert edi == {}
        assert figs == {}
        assert job == {}
        assert children == ["welcome"]


class TestAutoSave:
    def _fn(self, agent_app):
        return _find(agent_app, "am-store-messages", "am-store-history.data")

    def test_prevent_update_too_few_messages(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([{"role": "user", "content": "hi"}], [])

    def test_prevent_update_when_last_not_assistant(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        msgs = [
            {"role": "user", "content": "hi"},
            {"role": "user", "content": "again"},
        ]
        with pytest.raises(PreventUpdate):
            fn(msgs, [])

    def test_prevent_update_on_exact_restore(self, agent_app):
        from dash.exceptions import PreventUpdate

        msgs = [
            {"role": "user", "content": "hi"},
            {"role": "assistant", "content": "hello"},
        ]
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(msgs, [{"messages": msgs, "preview": "hi", "ts": "t"}])

    def test_inserts_new_entry_and_trims_to_20(self, agent_app):
        msgs = [
            {"role": "user", "content": "hi"},
            {"role": "assistant", "content": "hello"},
        ]
        history = [
            {"messages": [{"role": "x"}], "preview": f"old-{i}", "ts": "t"}
            for i in range(25)
        ]
        fn = self._fn(agent_app)
        result = fn(msgs, history)
        assert len(result) == 20
        assert result[0]["preview"] == "hi"
        assert result[0]["messages"] == msgs

    def test_replaces_existing_entry_with_same_preview(self, agent_app):
        msgs = [
            {"role": "user", "content": "hi"},
            {"role": "assistant", "content": "hello"},
        ]
        history = [{"messages": ["old"], "preview": "hi", "ts": "old-ts"}]
        fn = self._fn(agent_app)
        result = fn(msgs, history)
        assert len(result) == 1
        assert result[0]["messages"] == msgs


class TestSaveHistory:
    def _fn(self, agent_app):
        return _find(
            agent_app, "am-btn-save-session", "am-store-history.data"
        )

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, [{"role": "user", "content": "hi"}], [])

    def test_prevent_update_without_messages(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, [], [])

    def test_saves_entry(self, agent_app):
        fn = self._fn(agent_app)
        msgs = [{"role": "user", "content": "hello there"}]
        result = fn(1, msgs, [])
        assert result[0]["preview"] == "hello there"
        assert result[0]["messages"] == msgs


class TestRestoreSession:
    def teardown_method(self, _method):
        _clear_trigger()

    def _fn(self, agent_app):
        return _find(
            agent_app, '{"index":["ALL"],"type":"am-hist-item"}', ""
        )

    def test_prevent_update_without_clicks(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([None], [{"messages": [{"role": "user"}]}])

    def test_prevent_update_for_out_of_range_index(self, agent_app):
        from dash.exceptions import PreventUpdate

        _set_pattern_trigger({"index": 5, "type": "am-hist-item"})
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([1], [{"messages": [{"role": "user"}]}])

    def test_restores_messages_and_banner(self, agent_app):
        _set_pattern_trigger({"index": 0, "type": "am-hist-item"})
        history = [
            {
                "ts": "2026-01-01 10:00",
                "messages": [{"role": "user", "content": "hi"}],
            }
        ]
        fn = self._fn(agent_app)
        messages, children = fn([1], history)
        assert messages == history[0]["messages"]
        assert len(children) == 2  # banner + 1 restored bubble


class TestDeleteEntry:
    def teardown_method(self, _method):
        _clear_trigger()

    def _fn(self, agent_app):
        return _find(
            agent_app, '{"index":["ALL"],"type":"am-hist-del"}', ""
        )

    def test_prevent_update_without_clicks(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([None], [{"preview": "a"}])

    def test_deletes_entry_by_index(self, agent_app):
        _set_pattern_trigger({"index": 1, "type": "am-hist-del"})
        history = [{"preview": "a"}, {"preview": "b"}, {"preview": "c"}]
        fn = self._fn(agent_app)
        result = fn([1], history)
        assert result == [{"preview": "a"}, {"preview": "c"}]


class TestRenderHistory:
    def _fn(self, agent_app):
        return _find(agent_app, "am-store-history", "am-sidebar-history")

    def test_empty_state(self, agent_app):
        fn = self._fn(agent_app)
        result = fn([])
        assert "No saved sessions" in str(result)

    def test_renders_items(self, agent_app):
        fn = self._fn(agent_app)
        result = fn(
            [
                {"ts": "t1", "preview": "p1"},
                {"ts": "t2", "preview": "p2"},
            ]
        )
        assert len(result) == 2


class TestRenderFigs:
    def _fn(self, agent_app):
        return _find(agent_app, "am-store-figs", "am-sidebar-figs")

    def test_empty_state(self, agent_app):
        fn = self._fn(agent_app)
        result = fn({})
        assert "No figures yet" in str(result)

    def test_renders_thumbnails(self, agent_app):
        fn = self._fn(agent_app)
        result = fn({"fig1": {"b64": "abc", "title": "Fig 1"}})
        assert len(result) == 1
