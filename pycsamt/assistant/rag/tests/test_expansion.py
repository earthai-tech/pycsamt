# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for deterministic query expansion (RAG Tier 1)."""
from __future__ import annotations

import unittest

from pycsamt.assistant.rag.expansion import (
    QUERY_EXPANSIONS,
    expand_query,
)


class TestExpandQuery(unittest.TestCase):

    def test_maps_user_language_to_domain_terms(self):
        # "spikes/outliers" → the denoise vocabulary the corpus uses
        terms = set(expand_query("clean up spikes and outliers"))
        self.assertTrue({"denoise", "rpca", "hampel"} & terms)

    def test_depth_of_investigation_phrasing(self):
        terms = set(expand_query("how far down can I trust my model"))
        self.assertIn("sensitivity", terms)
        self.assertIn("bostick", terms)

    def test_irrelevant_query_expands_to_nothing(self):
        # Critical: nonsense must not gain spurious terms (keeps retrieval
        # empty rather than drifting).
        self.assertEqual(expand_query("completely unrelated xyzzy query"), [])

    def test_no_duplicate_or_echoed_terms(self):
        # A term already present in the query is not re-added.
        out = expand_query("denoise the spikes")
        self.assertNotIn("denoise", out)          # already typed
        self.assertEqual(len(out), len(set(out)))  # no dupes

    def test_triggers_are_lowercased_substrings(self):
        # Case-insensitive matching.
        self.assertTrue(expand_query("Which DIRECTION does the geology run"))

    def test_expansion_table_shape(self):
        # Every entry is (trigger, tuple-of-terms) and non-empty.
        for trigger, terms in QUERY_EXPANSIONS:
            self.assertIsInstance(trigger, str)
            self.assertTrue(trigger)
            self.assertIsInstance(terms, tuple)
            self.assertTrue(terms)


if __name__ == "__main__":
    unittest.main(verbosity=2)
