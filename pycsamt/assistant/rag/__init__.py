# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag
=====================

Retrieval-Augmented Generation for the pyCSAMT assistant.

Step 1 (this module set) covers **ingestion** — turning the package
source, docstrings and documentation into retrievable
:class:`~pycsamt.assistant.rag.schemas.RAGChunk` objects:

* :mod:`~pycsamt.assistant.rag.schemas`     — data types
* :mod:`~pycsamt.assistant.rag.config`      — what to index / exclude
* :mod:`~pycsamt.assistant.rag.ast_indexer` — Python → symbol chunks
* :mod:`~pycsamt.assistant.rag.doc_indexer` — .rst/.md → section chunks
* :mod:`~pycsamt.assistant.rag.ingest`      — walk the repo, build chunks

Later steps add the retriever, context builder, project registry and the
wiring into :class:`~pycsamt.agents.package_qa.PackageQAAgent`.
"""

from __future__ import annotations

from .schemas import RAGChunk, RetrievedContext

__all__ = ["RAGChunk", "RetrievedContext"]
