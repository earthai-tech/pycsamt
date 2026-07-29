# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the assistant project registry (RAG Step 2).

Hermetic: builds a temporary registry + EDI tree so it doesn't depend on
the bundled data. One test exercises the real bundled registry when
present.
"""

from __future__ import annotations

import textwrap
import unittest
from pathlib import Path
from tempfile import mkdtemp

from pycsamt.assistant.tools.project_registry import (
    ProjectRegistry,
    resolve_line,
)


def _make_registry(tmp: Path) -> Path:
    # a line that exists (with .edi files) and one that does not
    edi = tmp / "data" / "SURVEY" / "L22PLT"
    edi.mkdir(parents=True)
    for i in range(3):
        (edi / f"22-00{i}.edi").write_text(">HEAD\n", encoding="utf-8")

    reg = tmp / "projects" / "willy_project_registry.yml"
    reg.parent.mkdir(parents=True)
    reg.write_text(
        textwrap.dedent(
            """
            project:
              name: SURVEY
              type: AMT
              default_output_root: results/survey
            lines:
              L22PLT:
                edi_dir: data/SURVEY/L22PLT
                station_pattern: "22-*"
                sort_by: lon
                default_workflows: [load, qc, static_shift]
              L99ABS:
                edi_dir: /does/not/exist/L99
                station_pattern: "99-*"
            static_shift_defaults:
              method: ama
              half_window: 3
            plot_defaults:
              dpi: 150
            aliases:
              myline: L22PLT
            """
        ).strip(),
        encoding="utf-8",
    )
    return reg


class TestProjectRegistry(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(mkdtemp())
        self.reg_path = _make_registry(self.tmp)
        self.reg = ProjectRegistry(self.reg_path, root=self.tmp)

    def test_lines(self):
        self.assertEqual(self.reg.lines(), ["L22PLT", "L99ABS"])

    def test_resolve_existing_line(self):
        r = self.reg.resolve_line("L22PLT")
        self.assertEqual(r["line"], "L22PLT")
        self.assertTrue(r["exists"])
        self.assertEqual(r["n_edi_files"], 3)
        self.assertEqual(r["sort_by"], "lon")
        self.assertEqual(r["output_root"], "results/survey/L22PLT")
        self.assertEqual(r["static_shift"]["method"], "ama")
        self.assertTrue(Path(r["edi_dir"]).is_absolute())

    def test_messy_phrasing_and_case(self):
        for variant in ("line L22PLT", "l22plt", "  L22PLT  ", "line l22plt"):
            self.assertEqual(self.reg.resolve_line(variant)["line"], "L22PLT")

    def test_alias(self):
        self.assertEqual(self.reg.canonical("myline"), "L22PLT")
        self.assertEqual(self.reg.resolve_line("myline")["line"], "L22PLT")

    def test_missing_dir_marked_not_exists(self):
        r = self.reg.resolve_line("L99ABS")
        self.assertFalse(r["exists"])
        self.assertEqual(r["n_edi_files"], 0)

    def test_unknown_line_raises_with_known(self):
        with self.assertRaises(KeyError) as ctx:
            self.reg.resolve_line("L00XYZ")
        self.assertIn("L22PLT", str(ctx.exception))

    def test_find_line_in_text(self):
        self.assertEqual(
            self.reg.find_line_in_text(
                "write static shift code for line L22PLT please"
            ),
            "L22PLT",
        )
        self.assertIsNone(self.reg.find_line_in_text("just a general question"))

    def test_contains(self):
        self.assertIn("L22PLT", self.reg)
        self.assertIn("line l22plt", self.reg)
        self.assertNotIn("L00XYZ", self.reg)

    def test_relative_edi_dir_resolved_against_root(self):
        r = self.reg.resolve_line("L22PLT")
        self.assertTrue(r["edi_dir"].endswith("L22PLT") and "SURVEY" in r["edi_dir"])

    def test_resolve_line_convenience(self):
        r = resolve_line("L22PLT", registry_path=self.reg_path, root=self.tmp)
        self.assertEqual(r["line"], "L22PLT")


class TestBundledRegistry(unittest.TestCase):
    """Exercise the real registry shipped in projects/ when available."""

    def test_default_registry_resolves_willy(self):
        reg = ProjectRegistry.from_default()
        if reg is None:
            self.skipTest("no default registry present")
        self.assertIn("L22PLT", reg.lines())
        r = reg.resolve_line("L22PLT")
        # data is bundled in the repo, so it should resolve to an
        # existing dir with EDI files
        self.assertTrue(r["exists"])
        self.assertGreater(r["n_edi_files"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
