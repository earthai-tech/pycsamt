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
