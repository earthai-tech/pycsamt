# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for conversational query rewriting (RAG Tier 2)."""
from __future__ import annotations

import unittest

from pycsamt.assistant.rag.rewrite import (
    is_follow_up,
    rewrite_query,
)


class TestIsFollowUp(unittest.TestCase):

    def test_anaphoric_followups(self):
        for q in ("do that again", "run it once more", "same for this one",
                  "now do that"):
            self.assertTrue(is_follow_up(q), q)

    def test_self_contained_is_not_followup(self):
        # Names its own subject → not context-dependent, even with "again".
        self.assertFalse(is_follow_up("run static shift again"))
        self.assertFalse(is_follow_up("how do I correct static shift"))

    def test_fresh_question_not_followup(self):
        self.assertFalse(is_follow_up("what is ensure_sites"))


class TestRewriteQuery(unittest.TestCase):

    def test_followup_inherits_last_workflow(self):
        out = rewrite_query("do that again", last_workflow="static_shift")
        self.assertNotEqual(out, "do that again")
        self.assertIn("static", out.lower())

    def test_selfcontained_query_unchanged(self):
        q = "how do I run quality control"
        self.assertEqual(rewrite_query(q, last_workflow="static_shift"), q)

    def test_line_reference_inherits_line(self):
        out = rewrite_query(
            "now do the same on that line",
            last_workflow="qc", last_line="L22PLT",
        )
        self.assertIn("L22PLT", out)

    def test_recent_turns_supply_workflow(self):
        turns = [
            {"role": "user", "content": "run phase tensor analysis"},
            {"role": "assistant", "content": "done"},
        ]
        out = rewrite_query("do it again", recent_turns=turns)
        self.assertIn("phase", out.lower())

    def test_no_context_returns_query(self):
        # A follow-up with nothing to inherit is returned unchanged.
        self.assertEqual(rewrite_query("do that again"), "do that again")


if __name__ == "__main__":
    unittest.main(verbosity=2)
