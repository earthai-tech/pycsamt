# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.cli
=========================

Command-line helpers for the RAG index::

    python -m pycsamt.assistant.rag build          # build + save the index
    python -m pycsamt.assistant.rag stats          # show corpus stats
    python -m pycsamt.assistant.rag query "..."    # quick retrieval test
"""

from __future__ import annotations

import argparse
import os
import sys


def _cmd_build(args: argparse.Namespace) -> int:
    from .index_store import build_index

    key = None
    if args.embed:
        key = args.embed_key or os.environ.get("OPENAI_API_KEY")
        if not key:
            print(
                "--embed needs an API key: pass --embed-key or set "
                "OPENAI_API_KEY."
            )
            return 1
    mf = build_index(
        root=args.root,
        out_dir=args.out,
        embed=args.embed,
        embed_api_key=key,
        embed_provider=args.embed_provider,
        embed_model=args.embed_model,
    )
    print(f"Built index: {mf['n_chunks']} chunks -> {mf['out_dir']}")
    for k in sorted(mf["stats"]):
        if k.startswith("kind:") or k == "total":
            print(f"  {k}: {mf['stats'][k]}")
    if mf.get("embedded"):
        print(f"  embeddings: {mf['embed_model']} (dim={mf['embed_dim']})")
    return 0


def _cmd_stats(args: argparse.Namespace) -> int:
    from .index_store import index_exists, read_manifest

    if not index_exists(out_dir=args.out, root=args.root):
        print(
            "No persisted index. Run: python -m pycsamt.assistant.rag build"
        )
        return 1
    mf = read_manifest(out_dir=args.out, root=args.root) or {}
    print(
        f"Index created: {mf.get('created')}  ({mf.get('n_chunks')} chunks)"
    )
    if mf.get("embedded"):
        print(
            f"  embeddings: {mf.get('embed_model')} (dim={mf.get('embed_dim')})"
        )
    else:
        print("  embeddings: none (BM25 + expansion only)")
    for k in sorted(mf.get("stats", {})):
        print(f"  {k}: {mf['stats'][k]}")
    return 0


def _cmd_query(args: argparse.Namespace) -> int:
    from .retriever import build_retriever

    key = None
    if args.dense:
        key = args.embed_key or os.environ.get("OPENAI_API_KEY")
    retr = build_retriever(
        args.root, embed_api_key=key, embed_provider=args.embed_provider
    )
    if args.dense and not retr.dense_enabled:
        print(
            "(note: dense retrieval requested but inactive — build with "
            "--embed and provide a key)"
        )
    ctx = retr.search(args.text, k=args.k)
    if not ctx.chunks:
        print("(no results)")
        return 0
    for i, c in enumerate(ctx.chunks, 1):
        label = c.symbol or c.title or c.source_path
        print(f"{i:2d}. [{c.kind:13s} wf={c.workflow}] {label}")
        print(f"      {c.source_path}")
    return 0


def _cmd_eval(args: argparse.Namespace) -> int:
    from ..evals.harness import (
        evaluate,
        load_suite,
        suites_dir,
    )

    suite = args.suite or (suites_dir() / "rag_questions.jsonl")
    records = load_suite(suite)
    report = evaluate(records, k=args.k)
    print(report.summary())
    if report.violations:
        for v in report.violations:
            print(f"  ! {v['query']}: {v['found']}")
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pycsamt.assistant.rag",
        description="pyCSAMT RAG index tools",
    )
    parser.add_argument("--root", default=None, help="repo root")
    parser.add_argument("--out", default=None, help="index output dir")
    sub = parser.add_subparsers(dest="command", required=True)

    b = sub.add_parser("build", help="build and save the index")
    b.add_argument(
        "--embed",
        action="store_true",
        help="also compute + persist dense embeddings (needs an API key)",
    )
    b.add_argument("--embed-key", default=None, help="embedding API key")
    b.add_argument("--embed-provider", default="openai", help="embed backend")
    b.add_argument("--embed-model", default=None, help="embedding model id")
    sub.add_parser("stats", help="show stats of the persisted index")
    q = sub.add_parser("query", help="run a retrieval query")
    q.add_argument("text", help="query text")
    q.add_argument("-k", type=int, default=8, help="number of results")
    q.add_argument(
        "--dense",
        action="store_true",
        help="enable dense fusion (needs an embedded index + key)",
    )
    q.add_argument("--embed-key", default=None, help="embedding API key")
    q.add_argument("--embed-provider", default="openai", help="embed backend")

    e = sub.add_parser("eval", help="run an eval suite")
    e.add_argument("--suite", default=None, help="path to a .jsonl suite")
    e.add_argument("-k", type=int, default=10, help="retrieval depth")

    args = parser.parse_args(argv)
    return {
        "build": _cmd_build,
        "stats": _cmd_stats,
        "query": _cmd_query,
        "eval": _cmd_eval,
    }[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
