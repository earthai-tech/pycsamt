# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the executable assistant tools (RAG Step 7):
validate_generated_code + workflow runners.
"""
from __future__ import annotations

import os
import unittest
from pathlib import Path

from pycsamt.assistant.tools.validation_tools import validate_generated_code
from pycsamt.assistant.tools.workflow_tools import (
    resolve_target,
    run_workflow,
)


class TestValidateGeneratedCode(unittest.TestCase):

    def test_good_code(self):
        code = (
            "from pycsamt.emtools._core import ensure_sites\n"
            "from pycsamt.emtools.ss import estimate_ss_ama, correct_ss_ama\n"
            "sites = ensure_sites('/x')\n"
        )
        r = validate_generated_code(code)
        self.assertTrue(r["ok"])
        self.assertTrue(r["syntax_ok"])
        self.assertEqual(r["errors"], [])
        self.assertGreaterEqual(len(r["checked"]), 3)

    def test_hallucinated_symbol_flagged(self):
        code = (
            "from pycsamt import static_shift\n"
            "from pycsamt.emtools.ss import nonexistent_fn\n"
        )
        r = validate_generated_code(code)
        self.assertFalse(r["ok"])
        self.assertTrue(r["syntax_ok"])
        joined = " ".join(r["errors"])
        self.assertIn("static_shift", joined)
        self.assertIn("nonexistent_fn", joined)

    def test_syntax_error(self):
        r = validate_generated_code("def oops(:\n  pass\n")
        self.assertFalse(r["ok"])
        self.assertFalse(r["syntax_ok"])

    def test_empty(self):
        r = validate_generated_code("   ")
        self.assertFalse(r["ok"])

    def test_non_pycsamt_imports_ignored(self):
        code = "import numpy as np\nfrom math import sqrt\nx = sqrt(4)\n"
        r = validate_generated_code(code)
        self.assertTrue(r["ok"])
        self.assertEqual(r["checked"], [])  # nothing pycsamt to check

    def test_plain_module_import_checked(self):
        r = validate_generated_code("import pycsamt.emtools.ss\n")
        self.assertTrue(r["ok"])
        self.assertIn("pycsamt.emtools.ss", r["checked"])


class _FakeReg:
    def __contains__(self, name):
        return name == "L22PLT"

    def resolve_line(self, name):
        return {
            "line": "L22PLT", "edi_dir": "/data/L22PLT",
            "exists": True, "n_edi_files": 25,
            "output_root": "results/willy/L22PLT",
        }


class TestWorkflowTools(unittest.TestCase):

    def test_resolve_line(self):
        t = resolve_target("L22PLT", registry=_FakeReg())
        self.assertEqual(t["kind"], "line")
        self.assertEqual(t["line"], "L22PLT")

    def test_resolve_path(self):
        t = resolve_target("/some/path", registry=_FakeReg())
        self.assertEqual(t["kind"], "path")
        self.assertFalse(t["exists"])

    def test_unknown_workflow_raises(self):
        with self.assertRaises(ValueError):
            run_workflow("not_a_workflow", "/x", registry=_FakeReg())

    def test_missing_dir_raises(self):
        with self.assertRaises(FileNotFoundError):
            run_workflow(
                "static_shift", "/no/such/dir", registry=_FakeReg()
            )


@unittest.skipUnless(
    Path("data/3edis").is_dir(), "bundled 3edis data not present"
)
class TestRunStaticShiftBundled(unittest.TestCase):

    def test_run_on_bundled_data(self):
        import tempfile
        res = run_workflow(
            "static_shift", "data/3edis",
            output_dir=tempfile.mkdtemp(),
        )
        self.assertEqual(res["workflow"], "static_shift")
        self.assertEqual(res["status"], "success")
        self.assertEqual(res["path"], str(Path("data/3edis")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
