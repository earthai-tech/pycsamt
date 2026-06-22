# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the RAG ingestion layer (Step 1).

Pure stdlib — no DL backend, no Dash, no network. Uses temporary files
so it doesn't depend on the repo layout.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from pycsamt.assistant.rag import config
from pycsamt.assistant.rag.ast_indexer import index_python_file, module_name
from pycsamt.assistant.rag.doc_indexer import index_doc_file, split_sections
from pycsamt.assistant.rag.ingest import (
    build_chunks,
    corpus_stats,
    load_chunks,
    save_chunks,
)
from pycsamt.assistant.rag.schemas import RAGChunk


_PY_SRC = '''\
"""Static-shift helpers (galvanic distortion)."""

def estimate_ss_ama(sites, *, half_window=3, sort_by="lon"):
    """Estimate AMA static-shift factors."""
    return sites


def _private_helper(x):
    return x


class StaticShiftAgent:
    """Detect and correct galvanic static shift."""

    def __init__(self, method="ama"):
        self.method = method

    def execute(self, input_data):
        """Run the correction."""
        return input_data

    def _hidden(self):
        return 1
'''


class TestConfig(unittest.TestCase):

    def test_is_excluded(self):
        self.assertTrue(config.is_excluded("pycsamt/__pycache__/x.pyc"))
        self.assertTrue(config.is_excluded("pycsamt/app/web/app.py"))
        self.assertFalse(config.is_excluded("pycsamt/emtools/ss.py"))

    def test_should_index(self):
        self.assertTrue(config.should_index("pycsamt/emtools/ss.py"))
        self.assertTrue(config.should_index("docs/source/x.rst"))
        self.assertTrue(config.should_index("README.md"))
        self.assertFalse(config.should_index("pycsamt/app/web/app.py"))
        self.assertFalse(config.should_index("setup.cfg"))
        self.assertFalse(config.should_index("random/elsewhere.py"))

    def test_priority(self):
        self.assertEqual(config.priority_for("pycsamt/agents/x.py"), 3)
        self.assertEqual(config.priority_for("docs/source/x.rst"), 2)

    def test_infer_workflow(self):
        self.assertEqual(
            config.infer_workflow("pycsamt/emtools/ss.py", "correct_ss_ama"),
            "static_shift",
        )
        self.assertEqual(
            config.infer_workflow("phase tensor strike"), "phase_analysis"
        )
        self.assertIsNone(config.infer_workflow("nothing relevant here"))


class TestAstIndexer(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.pkg = self.tmp / "pycsamt" / "emtools"
        self.pkg.mkdir(parents=True)
        self.f = self.pkg / "ss.py"
        self.f.write_text(_PY_SRC, encoding="utf-8")

    def test_module_name(self):
        self.assertEqual(
            module_name("pycsamt/emtools/ss.py"), "pycsamt.emtools.ss"
        )
        self.assertEqual(
            module_name("pycsamt/emtools/__init__.py"), "pycsamt.emtools"
        )

    def test_symbols_and_methods(self):
        chunks = index_python_file(self.f, self.tmp)
        symbols = {c.symbol for c in chunks}
        # module doc + public function + class + public methods + __init__
        self.assertIn("pycsamt.emtools.ss", symbols)            # module doc
        self.assertIn("pycsamt.emtools.ss.estimate_ss_ama", symbols)
        self.assertIn("pycsamt.emtools.ss.StaticShiftAgent", symbols)
        self.assertIn(
            "pycsamt.emtools.ss.StaticShiftAgent.execute", symbols
        )
        self.assertIn(
            "pycsamt.emtools.ss.StaticShiftAgent.__init__", symbols
        )
        # private symbols excluded
        self.assertNotIn("pycsamt.emtools.ss._private_helper", symbols)
        self.assertNotIn(
            "pycsamt.emtools.ss.StaticShiftAgent._hidden", symbols
        )

    def test_signature_and_workflow_captured(self):
        chunks = index_python_file(self.f, self.tmp)
        fn = next(
            c for c in chunks
            if c.symbol == "pycsamt.emtools.ss.estimate_ss_ama"
        )
        self.assertIn("half_window", fn.metadata["signature"])
        self.assertEqual(fn.workflow, "static_shift")
        self.assertEqual(fn.kind, "python_symbol")
        self.assertIn("Estimate AMA", fn.text)

    def test_syntax_error_returns_empty(self):
        bad = self.pkg / "bad.py"
        bad.write_text("def oops(:\n", encoding="utf-8")
        self.assertEqual(index_python_file(bad, self.tmp), [])


class TestDocIndexer(unittest.TestCase):

    def test_split_md(self):
        secs = split_sections(
            "# Title\nintro\n## Sub\nbody text here", ".md"
        )
        titles = [t for t, _ in secs]
        self.assertIn("Title", titles)
        self.assertIn("Sub", titles)

    def test_split_rst(self):
        rst = "Heading\n=======\nsome body content here\n"
        secs = split_sections(rst, ".rst")
        self.assertTrue(any(t == "Heading" for t, _ in secs))

    def test_index_doc_file(self):
        tmp = Path(tempfile.mkdtemp())
        d = tmp / "docs" / "source"
        d.mkdir(parents=True)
        f = d / "static_shift.rst"
        f.write_text(
            "Static shift\n============\n"
            "How to correct galvanic static shift in your data.\n",
            encoding="utf-8",
        )
        chunks = index_doc_file(f, tmp)
        self.assertTrue(chunks)
        self.assertEqual(chunks[0].kind, "doc_section")
        self.assertTrue(any(c.workflow == "static_shift" for c in chunks))


class TestIngest(unittest.TestCase):

    def _make_tree(self) -> Path:
        tmp = Path(tempfile.mkdtemp())
        (tmp / "pycsamt" / "emtools").mkdir(parents=True)
        (tmp / "pycsamt" / "emtools" / "ss.py").write_text(
            _PY_SRC, encoding="utf-8"
        )
        # excluded: app + pyc
        (tmp / "pycsamt" / "app").mkdir(parents=True)
        (tmp / "pycsamt" / "app" / "ui.py").write_text(
            "x = 1\n", encoding="utf-8"
        )
        (tmp / "README.md").write_text(
            "# pyCSAMT\nProcessing suite.\n", encoding="utf-8"
        )
        return tmp

    def test_build_chunks_filters_and_indexes(self):
        tmp = self._make_tree()
        files = [
            tmp / "pycsamt" / "emtools" / "ss.py",
            tmp / "pycsamt" / "app" / "ui.py",   # excluded by config
            tmp / "README.md",
        ]
        # build_chunks indexes any file passed; the exclusion lives in
        # iter_index_files, so test the walker path too:
        from pycsamt.assistant.rag.ingest import iter_index_files
        walked = {p.name for p in iter_index_files(tmp)}
        self.assertIn("ss.py", walked)
        self.assertIn("README.md", walked)
        self.assertNotIn("ui.py", walked)  # pycsamt/app excluded

        chunks = build_chunks(tmp)
        self.assertGreater(len(chunks), 0)
        self.assertTrue(
            any(c.source_path == "README.md" for c in chunks)
        )
        self.assertFalse(
            any("app/ui.py" in c.source_path for c in chunks)
        )

    def test_save_load_roundtrip(self):
        tmp = self._make_tree()
        chunks = build_chunks(tmp)
        out = tmp / "chunks.jsonl"
        save_chunks(chunks, out)
        back = load_chunks(out)
        self.assertEqual(len(back), len(chunks))
        self.assertEqual(back[0].symbol, chunks[0].symbol)
        self.assertIsInstance(back[0], RAGChunk)

    def test_corpus_stats(self):
        tmp = self._make_tree()
        st = corpus_stats(build_chunks(tmp))
        self.assertEqual(
            st["total"],
            sum(v for k, v in st.items() if k.startswith("kind:")),
        )


class TestSchemas(unittest.TestCase):

    def test_chunk_roundtrip(self):
        c = RAGChunk(
            id="x:1:2", text="hello", source_path="a/b.py",
            kind="python_symbol", symbol="a.b.f", workflow="qc",
            priority=3, metadata={"signature": "f()"},
        )
        d = c.to_dict()
        c2 = RAGChunk.from_dict(d)
        self.assertEqual(c2.symbol, "a.b.f")
        self.assertEqual(c2.metadata["signature"], "f()")


if __name__ == "__main__":
    unittest.main(verbosity=2)
