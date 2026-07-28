# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Symbol cross-reference graph + signature-aware metadata (RAG Tier 3)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from pycsamt.assistant.rag.ast_indexer import (
    index_python_file,
)
from pycsamt.assistant.rag.graph import (
    SymbolGraph,
    build_symbol_graph,
)
from pycsamt.assistant.rag.schemas import RAGChunk

_SRC = '''"""Static-shift helpers."""
import numpy as np
import pycsamt.emtools.qc as qc
from pycsamt.exceptions import ProcessingError
from pycsamt.site.base import Sites


def ensure_sites(path, recursive=False):
    """Load EDIs."""
    return path


def correct_ss_ama(sites, half_window: int = 3, *, method: str = "ama") -> dict:
    """Correct static shift."""
    s = ensure_sites(sites)                  # local def  -> edge
    frame = Sites(s)                         # from-import -> edge
    qc.build_qc_table(frame)                 # aliased module -> edge
    if not all(frame):                       # builtin -> NO edge
        raise ProcessingError("bad")         # exception -> NO edge
    np.log(frame)                            # foreign module -> NO edge
    frame.unique()                           # foreign method -> NO edge
    return estimate_ss_ama(s, half_window)


def estimate_ss_ama(sites, half_window=3):
    """Estimate factors."""
    return {}


class Plain:
    """No call refs of its own."""
'''


def _indexed() -> list[RAGChunk]:
    tmp = Path(tempfile.mkdtemp())
    (tmp / "pycsamt" / "emtools").mkdir(parents=True)
    f = tmp / "pycsamt" / "emtools" / "ss.py"
    f.write_text(_SRC, encoding="utf-8")
    return index_python_file(f, tmp)


def _by_symbol() -> dict[str, RAGChunk]:
    return {c.symbol: c for c in _indexed() if c.symbol}


class TestSignatureMetadata(unittest.TestCase):
    """Tier 3: structured call surface on code chunks."""

    def setUp(self):
        self.by_symbol = _by_symbol()

    def test_params_name_annotation_default(self):
        m = self.by_symbol["pycsamt.emtools.ss.correct_ss_ama"].metadata
        params = {p["name"]: p for p in m["params"]}
        self.assertIsNone(params["sites"]["default"])
        self.assertEqual(params["half_window"]["annotation"], "int")
        self.assertEqual(params["half_window"]["default"], "3")
        # keyword-only args are captured too
        self.assertEqual(params["method"]["annotation"], "str")

    def test_return_annotation(self):
        m = self.by_symbol["pycsamt.emtools.ss.correct_ss_ama"].metadata
        self.assertEqual(m["returns"], "dict")

    def test_class_has_no_params_or_refs(self):
        m = self.by_symbol["pycsamt.emtools.ss.Plain"].metadata
        self.assertEqual(m["params"], [])
        self.assertEqual(m["refs"], [])


class TestRefsExtraction(unittest.TestCase):
    """Import-aware resolution: only real pyCSAMT call targets become edges."""

    def setUp(self):
        self.refs = _by_symbol()["pycsamt.emtools.ss.correct_ss_ama"].metadata["refs"]

    def test_local_def_resolves_to_module_qualified(self):
        self.assertIn("pycsamt.emtools.ss.ensure_sites", self.refs)
        self.assertIn("pycsamt.emtools.ss.estimate_ss_ama", self.refs)

    def test_from_import_resolves(self):
        self.assertIn("pycsamt.site.base.Sites", self.refs)

    def test_aliased_module_call_resolves(self):
        self.assertIn("pycsamt.emtools.qc.build_qc_table", self.refs)

    def test_builtins_and_foreign_calls_excluded(self):
        joined = " ".join(self.refs)
        for noise in ("all", "np.log", "unique"):
            self.assertNotIn(f".{noise}", joined, msg=f"{noise} leaked in")

    def test_raised_exceptions_excluded(self):
        self.assertFalse(any("ProcessingError" in r for r in self.refs))

    def test_leaf_with_no_calls_has_no_refs(self):
        self.assertEqual(
            _by_symbol()["pycsamt.emtools.ss.estimate_ss_ama"].metadata["refs"],
            [],
        )


class TestSymbolGraph(unittest.TestCase):
    def test_related_resolves_qualified_refs(self):
        g = SymbolGraph(_indexed())
        rel = g.related("pycsamt.emtools.ss.correct_ss_ama")
        self.assertIn("pycsamt.emtools.ss.ensure_sites", rel)
        self.assertIn("pycsamt.emtools.ss.estimate_ss_ama", rel)

    def test_unknown_targets_dropped(self):
        # Sites / build_qc_table are not in this single-file corpus.
        g = SymbolGraph(_indexed())
        rel = g.related("pycsamt.emtools.ss.correct_ss_ama")
        self.assertNotIn("pycsamt.site.base.Sites", rel)

    def test_never_returns_self(self):
        g = SymbolGraph(_indexed())
        self.assertNotIn(
            "pycsamt.emtools.ss.correct_ss_ama",
            g.related("pycsamt.emtools.ss.correct_ss_ama"),
        )

    def test_unknown_symbol_returns_empty(self):
        self.assertEqual(SymbolGraph(_indexed()).related("pycsamt.no.Thing"), [])

    def test_unique_leaf_recovers_reexport(self):
        # ref points at a package re-export; the real symbol lives deeper.
        chunks = [
            RAGChunk(
                id="a",
                text="",
                source_path="pycsamt/a.py",
                kind="python_symbol",
                symbol="pycsamt.a.caller",
                module="pycsamt.a",
                metadata={"refs": ["pycsamt.emtools.ensure_sites"]},
            ),
            RAGChunk(
                id="b",
                text="",
                source_path="pycsamt/emtools/_core.py",
                kind="python_symbol",
                symbol="pycsamt.emtools._core.ensure_sites",
                module="pycsamt.emtools._core",
                metadata={},
            ),
        ]
        self.assertEqual(
            SymbolGraph(chunks).related("pycsamt.a.caller"),
            ["pycsamt.emtools._core.ensure_sites"],
        )

    def test_ambiguous_leaf_not_recovered(self):
        # Two symbols share the leaf name → refuse to guess.
        chunks = [
            RAGChunk(
                id="a",
                text="",
                source_path="pycsamt/a.py",
                kind="python_symbol",
                symbol="pycsamt.a.caller",
                module="pycsamt.a",
                metadata={"refs": ["pycsamt.x.run"]},
            ),
            RAGChunk(
                id="b",
                text="",
                source_path="pycsamt/b.py",
                kind="python_symbol",
                symbol="pycsamt.b.run",
                module="pycsamt.b",
                metadata={},
            ),
            RAGChunk(
                id="c",
                text="",
                source_path="pycsamt/c.py",
                kind="python_symbol",
                symbol="pycsamt.c.run",
                module="pycsamt.c",
                metadata={},
            ),
        ]
        self.assertEqual(SymbolGraph(chunks).related("pycsamt.a.caller"), [])

    def test_build_symbol_graph_memoises(self):
        chunks = _indexed()
        self.assertIs(build_symbol_graph(chunks), build_symbol_graph(chunks))


if __name__ == "__main__":
    unittest.main(verbosity=2)
