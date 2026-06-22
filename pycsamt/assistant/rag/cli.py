# -*- coding: utf-8 -*-
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
import sys


def _cmd_build(args: argparse.Namespace) -> int:
    from .index_store import build_index
    mf = build_index(root=args.root, out_dir=args.out)
    print(f"Built index: {mf['n_chunks']} chunks -> {mf['out_dir']}")
    for k in sorted(mf["stats"]):
        if k.startswith("kind:") or k == "total":
            print(f"  {k}: {mf['stats'][k]}")
    return 0


def _cmd_stats(args: argparse.Namespace) -> int:
    from .index_store import read_manifest, index_exists
    if not index_exists(out_dir=args.out, root=args.root):
        print("No persisted index. Run: python -m pycsamt.assistant.rag build")
        return 1
    mf = read_manifest(out_dir=args.out, root=args.root) or {}
    print(f"Index created: {mf.get('created')}  ({mf.get('n_chunks')} chunks)")
    for k in sorted(mf.get("stats", {})):
        print(f"  {k}: {mf['stats'][k]}")
    return 0


def _cmd_query(args: argparse.Namespace) -> int:
    from .retriever import build_retriever
    retr = build_retriever(args.root)
    ctx = retr.search(args.text, k=args.k)
    if not ctx.chunks:
        print("(no results)")
        return 0
    for i, c in enumerate(ctx.chunks, 1):
        label = c.symbol or c.title or c.source_path
        print(f"{i:2d}. [{c.kind:13s} wf={c.workflow}] {label}")
        print(f"      {c.source_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pycsamt.assistant.rag",
        description="pyCSAMT RAG index tools",
    )
    parser.add_argument("--root", default=None, help="repo root")
    parser.add_argument("--out", default=None, help="index output dir")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("build", help="build and save the index")
    sub.add_parser("stats", help="show stats of the persisted index")
    q = sub.add_parser("query", help="run a retrieval query")
    q.add_argument("text", help="query text")
    q.add_argument("-k", type=int, default=8, help="number of results")

    args = parser.parse_args(argv)
    return {
        "build": _cmd_build,
        "stats": _cmd_stats,
        "query": _cmd_query,
    }[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
