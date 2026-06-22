# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Single-source-of-truth guard for the workflow registry.

The keyword table, alias map, and descriptions used to be duplicated in
:mod:`pycsamt.agents.context` and :mod:`pycsamt.agents.orchestrator`.  They are
now unified in :mod:`pycsamt.agents._workflows`.  These tests fail loudly if
the two consumers ever drift apart again.

Run::

    pytest pycsamt/agents/tests/test_workflow_registry.py -v
"""
from __future__ import annotations

import unittest

from pycsamt.agents._workflows import (
    WORKFLOW_ALIASES,
    WORKFLOW_DESCRIPTIONS,
    WORKFLOW_KEYWORDS,
    classify_workflow,
    normalise_workflow,
)

# Phrases that exercise the full table.
_PROBES = [
    "run a GCN 3D inversion",
    "U-Net 2D profile inversion",
    "ensemble uncertainty quantification",
    "joint inversion MT+TEM multi-modal",
    "prepare ModEM 3D inversion",
    "set up occam2d mesh and startup",
    "phase tensor and strike analysis",
    "denoise the survey",
    "write code for qc",
    "rotate tensors to strike",
    "bostick depth of investigation",
    "decimate frequencies",
    "tipper induction arrows",
    "run the full pipeline",
    "random gibberish with no keyword xyzzy",
]


class TestSharedRegistry(unittest.TestCase):
    """The registry module is well-formed."""

    def test_keyword_values_are_lists(self):
        for wf, kws in WORKFLOW_KEYWORDS.items():
            self.assertIsInstance(kws, list, msg=wf)
            self.assertGreater(len(kws), 0, msg=wf)

    def test_aliases_resolve_to_known_or_described(self):
        # Every alias target should be a real workflow id we describe
        # or that the keyword table can produce.
        producible = set(WORKFLOW_KEYWORDS) | set(WORKFLOW_DESCRIPTIONS)
        for alias, target in WORKFLOW_ALIASES.items():
            self.assertIn(
                target, producible,
                msg=f"alias {alias!r} -> unknown {target!r}",
            )

    def test_normalise_empty_defaults_qc(self):
        self.assertEqual(normalise_workflow(""), "qc")

    def test_normalise_resolves_aliases(self):
        self.assertEqual(normalise_workflow("PINN"), "pinn_inversion")
        self.assertEqual(normalise_workflow("1d-inversion"), "ai_inversion")
        self.assertEqual(normalise_workflow("codegen"), "code_gen")

    def test_classify_default_none(self):
        self.assertIsNone(classify_workflow("xyzzy no keywords here"))

    def test_classify_default_qc(self):
        self.assertEqual(
            classify_workflow("xyzzy no keywords here", default="qc"),
            "qc",
        )


class TestConsumersAgree(unittest.TestCase):
    """ContextInputAgent and the orchestrator must classify identically."""

    def test_context_and_orchestrator_agree(self):
        from pycsamt.agents.context import _regex_extract
        from pycsamt.agents.orchestrator import _keyword_classify

        disagreements = []
        for phrase in _PROBES:
            ctx = _regex_extract(phrase).get("workflow")
            ctx = normalise_workflow(ctx) if ctx else "qc"
            orch = _keyword_classify(phrase)[0]
            if ctx != orch:
                disagreements.append((phrase, ctx, orch))

        self.assertEqual(
            disagreements, [],
            msg=(
                "context vs orchestrator disagree (registry drift):\n"
                + "\n".join(
                    f"  {p!r}: context={c} orchestrator={o}"
                    for p, c, o in disagreements
                )
            ),
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
