# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Checkpoints must survive results that carry matplotlib figures.

Regression: AgentResult.data["figures"] holds Figure objects with
unpicklable closures (SecondaryAxis lambdas), so pickling the whole
result failed with "Could not save checkpoint … Can't pickle local
object". Figures are display-only; the checkpoint drops them.
"""
from __future__ import annotations

import pickle
import unittest

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

from pycsamt.agents.coordinator import AgentCoordinator
from pycsamt.agents._base import AgentResult


class TestCheckpointSafe(unittest.TestCase):

    def _fig_with_lambda(self):
        fig, ax = plt.subplots()
        # secondary axis uses lambdas → unpicklable
        ax.secondary_xaxis(
            "top", functions=(lambda x: x * 2, lambda x: x / 2)
        )
        return fig

    def test_figures_stripped_result_is_picklable(self):
        res = AgentResult(
            "success", "with fig",
            {
                "sites": [1, 2, 3],
                "figures": {"f": self._fig_with_lambda()},
                "single_fig": self._fig_with_lambda(),
                "fig_list": [self._fig_with_lambda()],
            },
        )
        safe = AgentCoordinator._checkpoint_safe(res)
        # all figure-bearing keys dropped, real data kept
        self.assertEqual(set(safe.data), {"sites"})
        self.assertEqual(safe.data["sites"], [1, 2, 3])
        # and it actually pickles now
        blob = pickle.dumps(safe)
        back = pickle.loads(blob)
        self.assertEqual(back.data["sites"], [1, 2, 3])
        self.assertEqual(back.status, "success")

    def test_no_figures_passthrough(self):
        res = AgentResult("success", "plain", {"a": 1, "b": "x"})
        safe = AgentCoordinator._checkpoint_safe(res)
        self.assertEqual(safe.data, {"a": 1, "b": "x"})

    def tearDown(self):
        plt.close("all")


if __name__ == "__main__":
    unittest.main(verbosity=2)
