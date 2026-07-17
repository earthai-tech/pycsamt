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
