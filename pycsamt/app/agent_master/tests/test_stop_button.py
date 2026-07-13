# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the Send ⇄ Stop button behaviour.

While the assistant executes a task the polling interval is enabled, the
button becomes a Stop button, and clicking it cancels the active job
(drops the thinking bubble, halts polling, leaves the input untouched).

Requires Dash (the chat callbacks module imports it); skipped otherwise.
"""

from __future__ import annotations

import importlib.util
import unittest

_HAS_DASH = importlib.util.find_spec("dash") is not None


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestStripThinking(unittest.TestCase):
    def test_removes_thinking_bubble(self):
        import pycsamt.app.agent_master.callbacks.chat as C

        msgs = [
            {"props": {"id": "msg-1"}},
            {"props": {"id": "am-thinking-bubble"}},
        ]
        out = C._strip_thinking(msgs)
        ids = [m["props"]["id"] for m in out]
        self.assertEqual(ids, ["msg-1"])

    def test_handles_none(self):
        import pycsamt.app.agent_master.callbacks.chat as C

        self.assertEqual(C._strip_thinking(None), [])


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestStopJobResponse(unittest.TestCase):
    def test_stop_cancels_job_and_halts_polling(self):
        import pycsamt.app.agent_master.callbacks.chat as C

        jid = C._new_job()  # status == "running"
        msgs = [
            {"props": {"id": "user-1"}},
            {"props": {"id": "am-thinking-bubble"}},
        ]
        (
            out_msgs,
            store_job,
            interval_disabled,
            input_val,
            new_stored,
            pending,
        ) = C._stop_job_response(
            msgs,
            [{"role": "user", "content": "go"}],
            {"jid": jid},
        )

        # job marked cancelled
        self.assertEqual(C._get_job(jid)["status"], "cancelled")
        # polling halted + job store cleared
        self.assertTrue(interval_disabled)
        self.assertEqual(store_job, {})
        # thinking bubble gone (replaced by the stop notice)
        thinking = [
            m
            for m in out_msgs
            if isinstance(m, dict)
            and m.get("props", {}).get("id") == "am-thinking-bubble"
        ]
        self.assertEqual(thinking, [])
        # transcript records the stop
        self.assertEqual(new_stored[-1]["role"], "assistant")
        self.assertIn("stopped", new_stored[-1]["content"].lower())
        # input preserved (no_update sentinel, not a real value)
        from dash import no_update

        self.assertIs(input_val, no_update)

    def test_stop_without_job_is_safe(self):
        import pycsamt.app.agent_master.callbacks.chat as C

        out = C._stop_job_response([], [], {})
        # no job id → no crash, still returns a 6-tuple
        self.assertEqual(len(out), 6)
        self.assertTrue(out[2])  # interval disabled


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestToggleCallbackRegistered(unittest.TestCase):
    def test_button_toggle_callback_present(self):
        from pycsamt.app.agent_master import create_app

        app = create_app()
        keys = [k for k in app.callback_map if "am-btn-send" in k]
        # children/className/title are driven by one callback
        self.assertTrue(
            any("am-btn-send.children" in k for k in keys),
            msg=f"no toggle callback on send button: {keys}",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
