# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.graph
===========================

A lightweight **symbol cross-reference graph** over the corpus.

Retrieval finds the chunk that best matches a query, but a professional
answer also points to *what that symbol works with* — the functions it
calls, the helpers it depends on. The AST indexer records, per code
chunk, the leaf-names it calls (``metadata["refs"]``); this module
resolves those leaf-names back to the pyCSAMT symbols that define them and
exposes a ``related(symbol)`` lookup.

It is pure-Python and built from the same :class:`RAGChunk` list the
retriever uses — no extra parsing, no I/O. Because the indexer resolves
call targets through each module's *imports* (not by bare name), the edges
are precise: a builtin ``all(...)``, a ``numpy`` call or a DataFrame
method never becomes an edge. Only two lookups are needed here — the
qualified target, or a *uniquely* named symbol when the qualified path was
a package re-export.
"""

from __future__ import annotations

from .schemas import RAGChunk

__all__ = ["SymbolGraph", "build_symbol_graph"]

_CODE_KINDS = ("python_symbol", "python_method", "module_doc")


class SymbolGraph:
    """Resolve a symbol's call-refs to other corpus symbols ("See also")."""

    def __init__(self, chunks: list[RAGChunk]) -> None:
        # every symbol the corpus defines
        self._known: set[str] = set()
        # leaf-name -> defining symbols (used only for unique-name recovery)
        self._by_leaf: dict[str, list[str]] = {}
        # fully-qualified symbol -> the qualified symbols it calls
        self._refs: dict[str, list[str]] = {}
        for c in chunks:
            if not c.symbol or c.kind not in _CODE_KINDS:
                continue
            self._known.add(c.symbol)
            self._by_leaf.setdefault(c.symbol.rsplit(".", 1)[-1], []).append(c.symbol)
            refs = c.metadata.get("refs") or []
            if refs:
                self._refs[c.symbol] = list(refs)

    def _resolve(self, ref: str) -> str | None:
        """Qualified hit, else a *uniquely* named symbol (re-export case)."""
        if ref in self._known:
            return ref
        candidates = self._by_leaf.get(ref.rsplit(".", 1)[-1], [])
        return candidates[0] if len(candidates) == 1 else None

    def related(self, symbol: str, *, max_related: int = 6) -> list[str]:
        """Symbols directly called by *symbol*, in source order."""
        out: list[str] = []
        for ref in self._refs.get(symbol, []):
            target = self._resolve(ref)
            if target and target != symbol and target not in out:
                out.append(target)
                if len(out) >= max_related:
                    break
        return out

    def __len__(self) -> int:
        return len(self._refs)


# ── cached builder ──────────────────────────────────────────────────────────────
_CACHE: dict[int, SymbolGraph] = {}


def build_symbol_graph(chunks: list[RAGChunk]) -> SymbolGraph:
    """Build (and memoise, keyed by chunk-list identity) a SymbolGraph."""
    key = id(chunks)
    g = _CACHE.get(key)
    if g is None:
        g = SymbolGraph(chunks)
        _CACHE[key] = g
    return g
