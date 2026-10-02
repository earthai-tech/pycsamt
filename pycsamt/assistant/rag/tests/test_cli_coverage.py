# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Targeted coverage for :mod:`pycsamt.assistant.rag.cli`.

``test_index_store.py::TestCli`` covers the plain build/stats/query happy
paths. This file rounds out: the ``--embed`` branches of ``build``
(missing key, and a successful embed via a monkeypatched backend), the
embedded-manifest branch of ``stats``, the dense-query branches of
``query`` (no key -> inactive note, and a zero-result query), and the
``eval`` subcommand (both the clean and the violations outcome).
"""

from __future__ import annotations

from pathlib import Path
from tempfile import mkdtemp

from pycsamt.assistant.rag import cli

_PY = '''\
"""Static-shift helpers."""

def estimate_ss_ama(sites, half_window=3):
    """Estimate AMA static-shift factors."""
    return sites
'''


def _tree() -> Path:
    tmp = Path(mkdtemp())
    (tmp / "pycsamt" / "emtools").mkdir(parents=True)
    (tmp / "pycsamt" / "emtools" / "ss.py").write_text(
        _PY, encoding="utf-8"
    )
    (tmp / "README.md").write_text(
        "# pyCSAMT\nProcessing suite.\n", encoding="utf-8"
    )
    return tmp


def test_build_embed_without_key_fails(monkeypatch, capsys):
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    root = _tree()
    rc = cli.main(["--root", str(root), "build", "--embed"])
    assert rc == 1
    out = capsys.readouterr().out
    assert "--embed needs an API key" in out


def test_build_embed_with_key_succeeds(monkeypatch):
    import pycsamt.assistant.rag.cli as cli_mod

    def _fake_build_index(**kwargs):
        assert kwargs["embed"] is True
        assert kwargs["embed_api_key"] == "sk-test"
        return {
            "n_chunks": 3,
            "out_dir": "somewhere",
            "stats": {"total": 3, "kind:code": 3},
            "embedded": True,
            "embed_model": "fake:test",
            "embed_dim": 4,
        }

    monkeypatch.setattr(
        "pycsamt.assistant.rag.index_store.build_index", _fake_build_index
    )
    root = _tree()
    rc = cli_mod.main(
        ["--root", str(root), "build", "--embed", "--embed-key", "sk-test"]
    )
    assert rc == 0


def test_stats_reports_embedded_manifest(tmp_path, capsys):
    import json

    out_dir = tmp_path / ".idx"
    out_dir.mkdir()
    (out_dir / "chunks.jsonl").write_text("", encoding="utf-8")
    (out_dir / "manifest.json").write_text(
        json.dumps(
            {
                "version": 2,
                "created": "now",
                "n_chunks": 0,
                "stats": {"total": 0},
                "embedded": True,
                "embed_model": "fake:test",
                "embed_dim": 4,
            }
        ),
        encoding="utf-8",
    )
    rc = cli.main(["--out", str(out_dir), "stats"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "embeddings: fake:test" in out


def test_query_dense_without_key_prints_inactive_note(monkeypatch, capsys):
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    root = _tree()
    rc = cli.main(
        ["--root", str(root), "query", "static shift", "-k", "3", "--dense"]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "dense retrieval requested but inactive" in out


def test_query_no_results_prints_placeholder(tmp_path, capsys):
    empty_root = tmp_path / "empty"
    empty_root.mkdir()
    rc = cli.main(
        ["--root", str(empty_root), "query", "nonexistent topic", "-k", "3"]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "(no results)" in out


def test_eval_clean_run_returns_zero(monkeypatch, capsys):
    import pycsamt.assistant.evals.harness as harness

    class _Report:
        violations = []

        def summary(self):
            return "Eval over 1 records:"

    monkeypatch.setattr(harness, "load_suite", lambda p: [{"query": "q1"}])
    monkeypatch.setattr(
        harness, "evaluate", lambda records, k=10: _Report()
    )
    rc = cli.main(["eval", "--suite", "unused.jsonl"])
    assert rc == 0
    assert "Eval over 1 records" in capsys.readouterr().out


def test_eval_with_violations_returns_one(monkeypatch, capsys):
    import pycsamt.assistant.evals.harness as harness

    class _Report:
        violations = [{"query": "bad q", "found": ["x"]}]

        def summary(self):
            return "Eval over 1 records:"

    monkeypatch.setattr(harness, "load_suite", lambda p: [{"query": "q1"}])
    monkeypatch.setattr(
        harness, "evaluate", lambda records, k=10: _Report()
    )
    rc = cli.main(["eval", "--suite", "unused.jsonl"])
    assert rc == 1
    out = capsys.readouterr().out
    assert "bad q" in out
