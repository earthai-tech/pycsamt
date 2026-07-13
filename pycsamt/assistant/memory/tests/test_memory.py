# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the assistant memory layer (RAG Step 9):
SessionState, ProjectState, WorkflowHistory (+ run_workflow tracing).
"""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import mkdtemp

from pycsamt.assistant.memory import (
    ProjectState,
    SessionState,
    WorkflowHistory,
    WorkflowRun,
)


class TestSessionState(unittest.TestCase):
    def test_record_and_roundtrip(self):
        s = SessionState()
        s.record_turn("user", "hi")
        s.record_turn("assistant", "hello")
        s.set_data(line="L22PLT", edi_path="/data/L22PLT")
        s.record_workflow("static_shift", "done")
        s.note("skew", "high")

        s2 = SessionState.from_dict(s.to_dict())
        self.assertEqual(s2.line, "L22PLT")
        self.assertEqual(s2.last_workflow, "static_shift")
        self.assertEqual(s2.facts["skew"], "high")
        self.assertEqual(len(s2.turns), 2)

    def test_recent_turns(self):
        s = SessionState()
        for i in range(10):
            s.record_turn("user", str(i))
        self.assertEqual(len(s.recent_turns(3)), 3)
        self.assertEqual(s.recent_turns(3)[-1]["content"], "9")

    def test_unique_ids(self):
        self.assertNotEqual(
            SessionState().session_id, SessionState().session_id
        )


class TestProjectState(unittest.TestCase):
    def setUp(self):
        self.path = Path(mkdtemp()) / "project_state.json"

    def test_remember_and_persist(self):
        ps = ProjectState(self.path)
        ps.remember_line("L22PLT", {"edi_dir": "/d/L22", "n_edi_files": 25})
        ps.set_preference("method", "ama")
        ps.save()

        ps2 = ProjectState(self.path)
        self.assertIn("L22PLT", ps2.lines_used())
        self.assertEqual(ps2.get_preference("method"), "ama")
        self.assertEqual(ps2.line_info("L22PLT")["n_edi_files"], 25)

    def test_missing_file_defaults(self):
        ps = ProjectState(self.path)  # not yet saved
        self.assertEqual(ps.lines_used(), [])
        self.assertIsNone(ps.get_preference("nope"))

    def test_lines_used_recency_order(self):
        ps = ProjectState(self.path)
        ps.remember_line("A", {})
        ps.remember_line("B", {})
        # set explicit timestamps (remember_line stamps "now", which is
        # the same second for both in a test)
        ps.data["lines"]["A"]["last_used"] = "2026-01-01T00:00:00"
        ps.data["lines"]["B"]["last_used"] = "2026-06-01T00:00:00"
        self.assertEqual(ps.lines_used()[0], "B")  # most recent first


class TestWorkflowHistory(unittest.TestCase):
    def setUp(self):
        self.path = Path(mkdtemp()) / "wh.jsonl"

    def test_record_and_read(self):
        h = WorkflowHistory(self.path)
        h.record(
            WorkflowRun(
                workflow="static_shift",
                status="success",
                line="L22PLT",
                summary="ok",
            )
        )
        h.record({"workflow": "qc", "status": "success"})
        runs = h.all()
        self.assertEqual(
            [r["workflow"] for r in runs], ["static_shift", "qc"]
        )
        self.assertEqual(runs[0]["line"], "L22PLT")
        self.assertIn("timestamp", runs[0])

    def test_recent(self):
        h = WorkflowHistory(self.path)
        for i in range(5):
            h.record({"workflow": f"w{i}"})
        self.assertEqual([r["workflow"] for r in h.recent(2)], ["w3", "w4"])

    def test_clear_and_empty(self):
        h = WorkflowHistory(self.path)
        self.assertEqual(h.all(), [])
        h.record({"workflow": "x"})
        h.clear()
        self.assertEqual(h.all(), [])


class TestRunWorkflowTracing(unittest.TestCase):
    def test_history_records_run(self):
        # run_workflow records to history when one is supplied
        from pycsamt.assistant.tools import (
            workflow_tools as wt,
        )

        # stub the agent so no real data is needed
        class _FakeResult:
            status = "success"
            summary = "stub"
            data = {"figures": {"f": object()}}

        class _FakeAgent:
            def __init__(self, **kw):
                pass

            def execute(self, _inp):
                return _FakeResult()

        import pycsamt.agents as A

        orig = A.DataQCAgent
        A.DataQCAgent = _FakeAgent
        try:
            reg = type(
                "R",
                (),
                {
                    "__contains__": lambda self, n: False,
                    "resolve_line": lambda self, n: {},
                },
            )()
            tmp = Path(mkdtemp())
            (tmp / "edi").mkdir()
            hist = WorkflowHistory(tmp / "wh.jsonl")
            res = wt.run_workflow(
                "qc",
                str(tmp / "edi"),
                output_dir=str(tmp / "out"),
                registry=reg,
                history=hist,
            )
            self.assertEqual(res["status"], "success")
            recorded = hist.all()
            self.assertEqual(len(recorded), 1)
            self.assertEqual(recorded[0]["workflow"], "qc")
            self.assertEqual(recorded[0]["n_figures"], 1)
        finally:
            A.DataQCAgent = orig


if __name__ == "__main__":
    unittest.main(verbosity=2)
