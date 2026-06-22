# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for the dash-free markdown parser used by the agent-master chat
to render typed responses (answers, generated code, capability cards).

The module is loaded directly by file path so these tests run even when the
GUI extras (Dash) are not installed — the parser only depends on ``re``.

Run::

    pytest pycsamt/app/agent_master/tests/test_markdown.py -v
"""
from __future__ import annotations

import importlib.util
import pathlib
import unittest

# Load _markdown.py in isolation (no Dash, no package __init__ chain).
_MOD_PATH = (
    pathlib.Path(__file__).resolve().parents[1] / "_markdown.py"
)
_spec = importlib.util.spec_from_file_location(
    "am_markdown_under_test", _MOD_PATH
)
md = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(md)  # type: ignore[union-attr]


class TestSplitInlineBold(unittest.TestCase):

    def test_basic(self):
        self.assertEqual(
            md.split_inline_bold("a **b** c"),
            [(False, "a "), (True, "b"), (False, " c")],
        )

    def test_plain(self):
        self.assertEqual(md.split_inline_bold("plain"), [(False, "plain")])

    def test_empty(self):
        self.assertEqual(md.split_inline_bold(""), [(False, "")])

    def test_leading_bold(self):
        self.assertEqual(
            md.split_inline_bold("**Name**: doc"),
            [(True, "Name"), (False, ": doc")],
        )


class TestParseMarkdown(unittest.TestCase):

    def test_heading_bullet_code(self):
        toks = md.parse_markdown("# Title\n- a\n```py\nx=1\n```")
        self.assertEqual(
            toks,
            [
                ("heading", "Title"),
                ("bullet", "a"),
                ("code", "py", "x=1"),
            ],
        )

    def test_fenced_code_multiline_default_lang(self):
        toks = md.parse_markdown("```\nline1\nline2\n```")
        self.assertEqual(toks, [("code", "python", "line1\nline2")])

    def test_paragraph_keeps_inline_markers(self):
        # inline bold is resolved at render time, not here
        self.assertEqual(
            md.parse_markdown("hello **world**"),
            [("para", "hello **world**")],
        )

    def test_blank_lines_tokenised(self):
        toks = md.parse_markdown("a\n\nb")
        self.assertEqual(
            toks,
            [("para", "a"), ("blank",), ("para", "b")],
        )

    def test_star_bullets(self):
        self.assertEqual(
            md.parse_markdown("* item"),
            [("bullet", "item")],
        )

    def test_unclosed_fence_consumes_rest(self):
        toks = md.parse_markdown("```py\nx=1\nstill code")
        self.assertEqual(toks, [("code", "py", "x=1\nstill code")])

    def test_empty_text(self):
        self.assertEqual(md.parse_markdown(""), [("blank",)])


if __name__ == "__main__":
    unittest.main(verbosity=2)
