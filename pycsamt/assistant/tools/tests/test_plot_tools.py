# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the line-aware plot tools (RAG Step 14).
"""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

from pycsamt.assistant.tools import plot_tools


class TestPlotToolsUnit(unittest.TestCase):
    def test_list_plots(self):
        kinds = plot_tools.list_plots()
        self.assertIn("snr", kinds)
        self.assertIn("confidence", kinds)

    def test_unknown_kind_raises(self):
        with self.assertRaises(ValueError):
            plot_tools.make_plot("not_a_kind", "/x")

    def test_missing_dir_raises(self):
        with self.assertRaises(FileNotFoundError):
            plot_tools.make_plot("snr", "/no/such/dir")

    def test_to_figure_coercion(self):
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        self.assertIs(plot_tools._to_figure(fig), fig)
        self.assertIs(plot_tools._to_figure(ax), fig)  # Axes
        self.assertIs(plot_tools._to_figure([ax]), fig)  # in a list
        self.assertIsNone(plot_tools._to_figure(None))
        self.assertIsNone(plot_tools._to_figure("nope"))
        plt.close(fig)


@unittest.skipUnless(Path("data/3edis").is_dir(), "bundled 3edis data not present")
class TestPlotToolsBundled(unittest.TestCase):
    def test_make_plot_saves_figure(self):
        out = tempfile.mkdtemp()
        r = plot_tools.make_plot("confidence", "data/3edis", output_dir=out)
        self.assertEqual(r["status"], "success")
        self.assertTrue(r["figure_path"])
        self.assertTrue(os.path.exists(r["figure_path"]))
        self.assertGreater(os.path.getsize(r["figure_path"]), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
