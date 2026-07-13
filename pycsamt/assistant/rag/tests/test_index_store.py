# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the persisted RAG index + CLI (RAG Step 5). Hermetic temp tree.
"""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import mkdtemp

from pycsamt.assistant.rag import cli
from pycsamt.assistant.rag.index_store import (
    build_index,
    index_exists,
    index_is_stale,
    load_index,
    read_manifest,
)
from pycsamt.assistant.rag.retriever import build_retriever

_PY = '''\
"""Static-shift helpers."""

def estimate_ss_ama(sites, half_window=3):
    """Estimate AMA static-shift factors."""
    return sites
'''


def _tree() -> Path:
    tmp = Path(mkdtemp())
    (tmp / "pycsamt" / "emtools").mkdir(parents=True)
    (tmp / "pycsamt" / "emtools" / "ss.py").write_text(_PY, encoding="utf-8")
    (tmp / "README.md").write_text(
        "# pyCSAMT\nProcessing suite.\n", encoding="utf-8"
    )
    return tmp


class TestIndexStore(unittest.TestCase):
    def setUp(self):
        self.root = _tree()

    def test_build_and_load(self):
        mf = build_index(root=self.root)
        self.assertGreater(mf["n_chunks"], 0)
        self.assertTrue(index_exists(root=self.root))
        chunks = load_index(root=self.root)
        self.assertEqual(len(chunks), mf["n_chunks"])
        self.assertTrue(
            any("estimate_ss_ama" in (c.symbol or "") for c in chunks)
        )

    def test_manifest(self):
        build_index(root=self.root)
        mf = read_manifest(root=self.root)
        self.assertIsNotNone(mf)
        self.assertEqual(mf["version"], 2)
        self.assertIn("stats", mf)
        self.assertIn("source_fingerprint", mf)

    def test_load_missing_returns_none(self):
        self.assertFalse(index_exists(root=self.root))
        self.assertIsNone(load_index(root=self.root))

    def test_build_no_save(self):
        mf = build_index(root=self.root, save=False)
        self.assertGreater(mf["n_chunks"], 0)
        self.assertFalse(index_exists(root=self.root))

    def test_retriever_prefers_persisted(self):
        build_index(root=self.root)
        # delete the source after building; persisted retriever still works
        (self.root / "pycsamt" / "emtools" / "ss.py").unlink()
        r = build_retriever(self.root, use_cache=False, prefer_persisted=True)
        ctx = r.search("estimate_ss_ama AMA static shift", k=3)
        self.assertTrue(ctx.chunks)  # came from the persisted index

    def test_retriever_fresh_ingest_when_disabled(self):
        # no index built; prefer_persisted ignored → fresh ingest
        r = build_retriever(
            self.root, use_cache=False, prefer_persisted=False
        )
        self.assertGreater(len(r.chunks), 0)


class TestFreshness(unittest.TestCase):
    """The stale-index guard keys on *content*, not mtime (Tier 3 fix)."""

    def setUp(self):
        self.root = _tree()
        build_index(root=self.root)

    def test_fresh_right_after_build(self):
        self.assertFalse(index_is_stale(root=self.root))

    def test_touching_mtime_does_not_invalidate(self):
        # Running the test suite bumps mtimes without changing content —
        # that must not report a false "stale index".
        import os
        import time

        f = self.root / "pycsamt" / "emtools" / "ss.py"
        future = time.time() + 120
        os.utime(f, (future, future))
        self.assertFalse(index_is_stale(root=self.root))

    def test_content_change_invalidates(self):
        f = self.root / "pycsamt" / "emtools" / "ss.py"
        f.write_text(
            f.read_text(encoding="utf-8") + "\n# edited\n", encoding="utf-8"
        )
        self.assertTrue(index_is_stale(root=self.root))

    def test_new_indexed_file_invalidates(self):
        (self.root / "pycsamt" / "emtools" / "qc.py").write_text(
            '"""QC."""\n', encoding="utf-8"
        )
        self.assertTrue(index_is_stale(root=self.root))

    def test_old_manifest_version_is_stale(self):
        import json

        mf_path = self.root / ".pycsamt_rag" / "manifest.json"
        mf = json.loads(mf_path.read_text(encoding="utf-8"))
        mf["version"] = 1  # pre-content-hash scheme
        mf_path.write_text(json.dumps(mf), encoding="utf-8")
        self.assertTrue(index_is_stale(root=self.root))


class TestCli(unittest.TestCase):
    def test_build_then_stats(self):
        root = _tree()
        self.assertEqual(cli.main(["--root", str(root), "build"]), 0)
        self.assertEqual(cli.main(["--root", str(root), "stats"]), 0)

    def test_stats_without_index(self):
        root = _tree()
        self.assertEqual(cli.main(["--root", str(root), "stats"]), 1)

    def test_query(self):
        root = _tree()
        rc = cli.main(
            ["--root", str(root), "query", "static shift", "-k", "3"]
        )
        self.assertEqual(rc, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
