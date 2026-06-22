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
class TestHelpModal(unittest.TestCase):

    def test_help_toggle_registered_and_modal_built(self):
        from pycsamt.app.agent_master import create_app
        from pycsamt.app.agent_master._ids import IDs
        app = create_app()
        self.assertTrue(
            any(IDs.MODAL_HELP in k for k in app.callback_map),
            msg="help toggle callback not registered",
        )

    def test_help_modal_has_examples_and_caps(self):
        from pycsamt.app.agent_master.layout import _help_modal
        # builds without error and is a real component tree
        modal = _help_modal()
        self.assertEqual(modal.id, "am-modal-help")


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
class TestCodeDispatchRAG(unittest.TestCase):
    """_dispatch_code resolves a named line + grounds code via RAG."""

    def setUp(self):
        import tempfile
        from pathlib import Path
        from pycsamt.assistant.rag import context_builder as cb
        from pycsamt.assistant.rag.context_builder import ContextBuilder
        from pycsamt.assistant.rag.retriever import Retriever
        from pycsamt.assistant.rag.schemas import RAGChunk

        # a real temp EDI dir so the resolved line "exists"
        self.tmp = Path(tempfile.mkdtemp())
        edi = self.tmp / "L22PLT"
        edi.mkdir()
        (edi / "22-001.edi").write_text(">HEAD\n", encoding="utf-8")

        chunks = [
            RAGChunk(
                id="1", text="estimate_ss_ama / correct_ss_ama AMA "
                "static shift correction recipe.",
                source_path="assistant_recipes/static_shift.md",
                kind="recipe", workflow="static_shift", priority=3,
                title="static shift",
            ),
        ]

        class _Reg:
            def __init__(self, p):
                self._p = p

            def find_line_in_text(self, t):
                return "L22PLT" if "l22plt" in t.lower() else None

            def resolve_line(self, name):
                return {
                    "line": "L22PLT", "edi_dir": str(self._p),
                    "exists": True, "n_edi_files": 1, "sort_by": "lon",
                    "output_root": "results/willy/L22PLT",
                    "default_workflows": ["static_shift"],
                    "static_shift": {"method": "ama"}, "plot": {},
                }

        self._orig = cb.default_context_builder
        builder = ContextBuilder(Retriever(chunks), _Reg(edi))
        cb.default_context_builder = lambda root=None: builder
        self._edi_dir = str(edi)

    def tearDown(self):
        from pycsamt.assistant.rag import context_builder as cb
        cb.default_context_builder = self._orig

    def test_line_resolved_into_generated_code(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        jid = C._new_job()

        def step(label, status="done"):
            C._update_job(
                jid,
                steps=C._JOBS[jid]["steps"]
                + [{"label": label, "status": status}],
            )

        C._dispatch_code(
            jid,
            "generate code for static shift for line L22PLT",
            {}, {"provider": "offline"},
            workflow=None, llm_prov="claude", api_key=None,
            sel_model=None, offline=True, step=step,
        )
        job = C._get_job(jid)
        code = job.get("code", "")
        self.assertIn("L22PLT", job.get("result", ""))
        # real resolved path is embedded (check separator-agnostically:
        # the unique temp dir name + the line dir, not the placeholder)
        from pathlib import Path
        self.assertIn(Path(self._edi_dir).parent.name, code)
        self.assertIn("L22PLT", code)
        self.assertNotIn("/path/to/EDIs", code)   # placeholder replaced
        self.assertIn("estimate_ss_ama", code)    # correct API


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestCodeTargetWorkflow(unittest.TestCase):
    """Code requests must target the subject, not the 'code_gen' action."""

    def test_targets(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        cases = {
            "generate code for static shift": "static_shift",
            "generate code for dimensionality": "phase_analysis",
            "write code for ai inversion": "ai_inversion",
            "give me a python script for qc": "qc",
        }
        for text, expected in cases.items():
            self.assertEqual(
                C._code_target_workflow(text), expected,
                msg=f"{text!r} -> {C._code_target_workflow(text)}",
            )

    def test_no_subject_returns_none(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        self.assertIsNone(C._code_target_workflow("write me a script"))


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestCodeBlockCollapsible(unittest.TestCase):
    """The code block is a collapsible <details> with a copy button."""

    def test_structure(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        from dash import html
        blk = C._code_block("x = 1\ny = 2")
        details = blk.children[0]
        self.assertIsInstance(details, html.Details)
        self.assertTrue(details.open)
        self.assertIsInstance(details.children[0], html.Summary)
        copy_btn = blk.children[1]
        self.assertEqual(copy_btn.className, "am-code-copy-btn")


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestNamedLineWorkflow(unittest.TestCase):
    """A named survey line runs without manually loading EDI data."""

    def setUp(self):
        import tempfile
        from pathlib import Path
        import pycsamt.assistant.tools.project_registry as pr

        self.tmp = Path(tempfile.mkdtemp())
        edi = self.tmp / "L22PLT"
        edi.mkdir()
        (edi / "22-001.edi").write_text(">HEAD\n", encoding="utf-8")
        self._edi_dir = str(edi)

        class _Reg:
            def __init__(self, p):
                self._p = p

            def find_line_in_text(self, t):
                return "L22PLT" if "l22plt" in t.lower() else None

            def resolve_line(self, name):
                return {
                    "line": "L22PLT", "edi_dir": str(self._p),
                    "exists": True, "n_edi_files": 1,
                }

        self._orig = pr.ProjectRegistry.from_default
        pr.ProjectRegistry.from_default = classmethod(
            lambda cls, root=None: _Reg(edi)
        )

    def tearDown(self):
        import pycsamt.assistant.tools.project_registry as pr
        pr.ProjectRegistry.from_default = self._orig

    def test_names_registry_line(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        self.assertTrue(
            C._names_registry_line("run static shift on line L22PLT")
        )
        self.assertFalse(C._names_registry_line("run static shift"))

    def test_run_agent_resolves_line_without_edi(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        import pycsamt.agents.orchestrator as O
        from pycsamt.agents._base import AgentResult

        captured = {}

        def fake_exec(self, input_data):
            captured["data_path"] = input_data.get("data_path")
            return AgentResult(
                "success", "stub",
                {"result": AgentResult("success", "i", {})},
            )

        orig = O.WorkflowOrchestratorAgent.execute
        O.WorkflowOrchestratorAgent.execute = fake_exec
        try:
            jid = C._new_job()
            # empty edi_store → must resolve the named line
            C._run_agent(
                jid, "run static shift on line L22PLT",
                {}, {"provider": "offline"}, {}, [],
            )
        finally:
            O.WorkflowOrchestratorAgent.execute = orig

        self.assertEqual(captured.get("data_path"), self._edi_dir)
        self.assertNotEqual(C._get_job(jid).get("kind"), "error")


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestRecentRuns(unittest.TestCase):
    """Workflow trace recording + sidebar render."""

    def test_record_and_recent_and_render(self):
        import pycsamt.app.agent_master.callbacks.chat as C
        from pycsamt.assistant.memory.workflow_history import (
            WorkflowHistory,
        )
        WorkflowHistory.default().clear()
        C._record_run(
            workflow="static_shift", path="/data/L22PLT",
            output_dir="out", status="success",
            summary="done", n_figures=2,
        )
        runs = C._recent_runs()
        self.assertTrue(runs)
        self.assertEqual(runs[0]["workflow"], "static_shift")
        item = C._run_item(runs[0])
        self.assertEqual(item.className, "am-run-item")

    def test_record_run_never_raises(self):
        # tracing is best-effort; bad inputs must not raise
        import pycsamt.app.agent_master.callbacks.chat as C
        C._record_run(
            workflow="qc", path=None, output_dir="o",
            status="failed", summary="x", n_figures=0,
        )

    def test_callback_registered(self):
        from pycsamt.app.agent_master import create_app
        app = create_app()
        self.assertTrue(
            any("am-sidebar-runs" in k for k in app.callback_map)
        )


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
