# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant
=================

Retrieval-augmented assistant layer for pyCSAMT v2.

This package adds the assistant's *technical memory* on top of the
existing agent layer (:mod:`pycsamt.agents`): a RAG index over the
package source, docstrings, docs, examples and recipes, plus executable
tools and project/session memory.

It is built **gradually** (see ``PYCSAMT-V2-RAG-IMPLEMENTATION.md``):

* ``rag/``    — ingestion, indexing, retrieval, context building
* ``tools/``  — real executable pyCSAMT tools (added later)
* ``memory/`` — session / project state (added later)

Nothing here is imported eagerly by the base package, so it adds no cost
to ``import pycsamt``.
"""

from __future__ import annotations

__all__: list[str] = []
