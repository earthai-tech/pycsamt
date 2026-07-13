# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Regression guard: the IntentRouter decides *intent only*; it must never
override the workflow *type* chosen by ContextInputAgent.

A terse-prompt LLM router can mislabel e.g. "phase tensor and dimensionality
analysis" as an inversion.  The ordered keyword registry inside
ContextInputAgent classifies it correctly, and that classification must win.

Requires Dash (the chat callbacks module imports it); skipped otherwise.
"""

from __future__ import annotations

import importlib.util
import unittest

_HAS_DASH = importlib.util.find_spec("dash") is not None


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestWorkflowTypeAuthority(unittest.TestCase):
    def test_router_workflow_slot_does_not_override(self):
        import pycsamt.agents.orchestrator as O
        import pycsamt.agents.router as R
        import pycsamt.app.agent_master.callbacks.chat as C
        from pycsamt.agents._base import AgentResult
        from pycsamt.agents.router import (
            WORKFLOW,
            RouterDecision,
        )

        # 1. Router returns the correct intent but a WRONG workflow slot.
        orig_route = R.IntentRouter.route
        R.IntentRouter.route = lambda self, text, history=None: (
            RouterDecision(
                intent=WORKFLOW,
                workflow="ai_inversion",
                confidence=0.9,
                source="llm",
            )
        )

        # 2. Capture the config the orchestrator actually receives,
        #    without running any real workflow.
        captured = {}

        def fake_execute(self, input_data):
            captured["workflow"] = input_data.get("config", {}).get(
                "workflow"
            )
            return AgentResult(
                status="success",
                summary="stub",
                data={"result": AgentResult("success", "inner", {})},
            )

        orig_exec = O.WorkflowOrchestratorAgent.execute
        O.WorkflowOrchestratorAgent.execute = fake_execute

        try:
            jid = C._new_job()
            C._run_agent(
                jid,
                "Run phase tensor and dimensionality analysis",
                {"path": "/some/edi/path", "groups": {}},
                {"provider": "offline"},
                {},
                [],
            )
        finally:
            R.IntentRouter.route = orig_route
            O.WorkflowOrchestratorAgent.execute = orig_exec

        # 3. ContextInputAgent's classification (phase_analysis) wins;
        #    the router's wrong ai_inversion slot is ignored.
        self.assertEqual(captured.get("workflow"), "phase_analysis")


if __name__ == "__main__":
    unittest.main(verbosity=2)
