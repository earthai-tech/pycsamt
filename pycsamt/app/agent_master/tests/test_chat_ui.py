# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the executing-message redesign and the pin-message feature.

Requires Dash (the chat callbacks module imports it); skipped otherwise.
"""
from __future__ import annotations

import importlib.util
import unittest

_HAS_DASH = importlib.util.find_spec("dash") is not None


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestExecutingMessage(unittest.TestCase):

    def test_elapsed_format(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        self.assertEqual(C._fmt_elapsed(0), "0:00")
        self.assertEqual(C._fmt_elapsed(3), "0:03")
        self.assertEqual(C._fmt_elapsed(75), "1:15")
        self.assertEqual(C._fmt_elapsed(-5), "0:00")

    def test_thinking_bubble_has_timeline_and_id(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        steps = [
            {"label": "Parsing request...", "status": "done"},
            {"label": "Executing phase_analysis...",
             "status": "running"},
        ]
        b = C._thinking_bubble(
            steps, workflow="phase_analysis", elapsed=3.0
        )
        self.assertEqual(b.id, "am-thinking-bubble")

    def test_thinking_bubble_no_steps(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        b = C._thinking_bubble([])
        self.assertEqual(b.id, "am-thinking-bubble")

    def test_workflow_label_lookup(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        self.assertEqual(
            C._WF_RUNNING_LABEL["phase_analysis"],
            "phase tensor analysis",
        )
        self.assertIn("code_gen", C._WF_RUNNING_LABEL)


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestPinLogic(unittest.TestCase):

    def _messages(self):
        return [
            {"role": "user", "content": "run qc",
             "ts": "10:00", "mid": "am-msg-1"},
            {"role": "assistant",
             "content": "QC complete: 25 stations.",
             "ts": "10:01", "mid": "am-msg-2"},
        ]

    def test_pin_adds_entry(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        pins = C._apply_pin_toggle([], "am-msg-2", self._messages())
        self.assertEqual(len(pins), 1)
        self.assertEqual(pins[0]["mid"], "am-msg-2")
        self.assertEqual(pins[0]["role"], "assistant")
        self.assertIn("QC complete", pins[0]["snippet"])

    def test_pin_toggles_off(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        msgs = self._messages()
        pins = C._apply_pin_toggle([], "am-msg-2", msgs)
        pins = C._apply_pin_toggle(pins, "am-msg-2", msgs)
        self.assertEqual(pins, [])

    def test_pin_unknown_message_raises(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        with self.assertRaises(KeyError):
            C._apply_pin_toggle([], "am-msg-nope", self._messages())

    def test_remove_pin(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        pins = [{"mid": "a"}, {"mid": "b"}]
        self.assertEqual(
            C._remove_pin(pins, "a"), [{"mid": "b"}]
        )

    def test_snippet_truncates(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        long = "word " * 40
        snip = C._pin_snippet(long, limit=30)
        self.assertLessEqual(len(snip), 30)
        self.assertTrue(snip.endswith("…"))

    def test_bubbles_carry_mid_and_pin_button(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        ub = C._user_bubble("hi", mid="am-msg-1")
        ab = C._agent_bubble("done", mid="am-msg-2")
        self.assertEqual(ub.id, "am-msg-1")
        self.assertEqual(ab.id, "am-msg-2")
        # bubbles without a mid get no id (not pinnable)
        ab2 = C._agent_bubble("notice")
        self.assertIsNone(getattr(ab2, "id", None))


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestWorkflowNotForcedByStaleConfig(unittest.TestCase):
    """A persisted param-modal workflow must not hijack later requests."""

    def test_drop_workflow_strips_only_workflow(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        ic = {"workflow": "ai_inversion", "n_layers": 10,
              "epochs": 500}
        out = C._drop_workflow(ic)
        self.assertNotIn("workflow", out)
        self.assertEqual(out["n_layers"], 10)
        self.assertEqual(out["epochs"], 500)

    def test_drop_workflow_handles_none(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        self.assertEqual(C._drop_workflow(None), {})

    def test_request_routes_by_text_not_stale_workflow(self):
        import io, contextlib
        import pycsamt.app.agent_master.callbacks.chat as C
        import pycsamt.agents.orchestrator as O
        from pycsamt.agents._base import AgentResult

        captured = {}

        def fake_execute(self, input_data):
            captured["wf"] = (
                input_data.get("config", {}).get("workflow")
            )
            return AgentResult(
                "success", "stub",
                {"result": AgentResult("success", "i", {})},
            )

        orig = O.WorkflowOrchestratorAgent.execute
        O.WorkflowOrchestratorAgent.execute = fake_execute
        try:
            # Simulate the contaminated localStorage config, cleaned by
            # the caller (send_message / line_sel) before _run_agent.
            stale = {"workflow": "ai_inversion", "n_layers": 10}
            clean = C._drop_workflow(stale)
            jid = C._new_job()
            with contextlib.redirect_stdout(io.StringIO()):
                C._run_agent(
                    jid, "run static shift",
                    {"path": "/x", "groups": {}},
                    {"provider": "offline"}, clean, [],
                )
        finally:
            O.WorkflowOrchestratorAgent.execute = orig

        self.assertEqual(captured.get("wf"), "static_shift")


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestPinCallbacksRegistered(unittest.TestCase):

    def test_pin_callbacks_present(self):
        from pycsamt.app.agent_master import create_app
        app = create_app()
        keys = list(app.callback_map)
        self.assertTrue(
            any("am-store-pins" in k for k in keys),
            msg="no STORE_PINS callback registered",
        )
        self.assertTrue(
            any("am-sidebar-pins" in k for k in keys),
            msg="no sidebar pins render callback",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
