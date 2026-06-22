# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.rag.ast_indexer
=================================

Turn a Python source file into :class:`~pycsamt.assistant.rag.schemas.RAGChunk`
objects using the AST — one chunk per public symbol (module docstring,
top-level function/class, and each public method of a class).

Plain text chunking loses the structure the assistant needs (which class,
which function, which parameters); indexing by symbol preserves it and
lets the retriever return precise, citable API entries.
"""

from __future__ import annotations

import ast
from pathlib import Path

from .config import infer_workflow, priority_for
from .schemas import RAGChunk

__all__ = ["index_python_file", "module_name"]

_MAX_CODE_CHARS = 4000  # keep chunks embeddable; trim very long bodies


def module_name(rel_path: str) -> str:
    """``pycsamt/emtools/ss.py`` → ``pycsamt.emtools.ss``."""
    p = rel_path.replace("\\", "/")
    if p.endswith("/__init__.py"):
        p = p[: -len("/__init__.py")]
    elif p.endswith(".py"):
        p = p[: -len(".py")]
    return p.strip("/").replace("/", ".")


def _signature(node: ast.AST) -> str:
    """Best-effort ``name(args) -> returns`` for a def/class node."""
    name = getattr(node, "name", "")
    if isinstance(node, ast.ClassDef):
        bases = []
        for b in node.bases:
            try:
                bases.append(ast.unparse(b))
            except Exception:
                pass
        return f"class {name}({', '.join(bases)})" if bases else f"class {name}"
    try:
        args = ast.unparse(node.args)  # type: ignore[attr-defined]
    except Exception:
        args = "..."
    ret = ""
    returns = getattr(node, "returns", None)
    if returns is not None:
        try:
            ret = f" -> {ast.unparse(returns)}"
        except Exception:
            ret = ""
    return f"{name}({args}){ret}"


def _source_segment(source_lines: list[str], node: ast.AST) -> str:
    start = getattr(node, "lineno", 1)
    end = getattr(node, "end_lineno", start)
    code = "\n".join(source_lines[start - 1:end])
    if len(code) > _MAX_CODE_CHARS:
        code = code[:_MAX_CODE_CHARS] + "\n# … (truncated)"
    return code


def _make_chunk(
    *,
    rel_path: str,
    module: str,
    name: str,
    qualifier: str,
    node: ast.AST,
    source_lines: list[str],
    kind: str,
) -> RAGChunk:
    doc = ast.get_docstring(node) or ""
    sig = _signature(node)
    code = _source_segment(source_lines, node)
    symbol = f"{module}.{qualifier}" if qualifier else module
    start = getattr(node, "lineno", 1)
    end = getattr(node, "end_lineno", start)
    text = (
        f"Symbol: {symbol}\n"
        f"Signature: {sig}\n"
        f"File: {rel_path}\n\n"
        f"Docstring:\n{doc}\n\n"
        f"Code:\n{code}"
    ).strip()
    return RAGChunk(
        id=f"{rel_path}:{start}:{end}",
        text=text,
        source_path=rel_path,
        kind=kind,
        symbol=symbol,
        module=module,
        title=name,
        workflow=infer_workflow(rel_path, symbol, doc),
        priority=priority_for(rel_path),
        metadata={
            "name": name,
            "signature": sig,
            "start_line": start,
            "end_line": end,
            "has_doc": bool(doc),
        },
    )


def index_python_file(path: Path, root: Path) -> list[RAGChunk]:
    """Index one ``.py`` file into a list of :class:`RAGChunk`.

    Captures the module docstring, each top-level function/class, and each
    public method of every top-level class. Returns ``[]`` on a syntax
    error (best-effort ingestion).
    """
    rel_path = path.relative_to(root).as_posix()
    try:
        source = path.read_text(encoding="utf-8", errors="ignore")
        tree = ast.parse(source)
    except (SyntaxError, ValueError):
        return []

    source_lines = source.splitlines()
    module = module_name(rel_path)
    chunks: list[RAGChunk] = []

    # module-level docstring
    mod_doc = ast.get_docstring(tree)
    if mod_doc:
        chunks.append(
            RAGChunk(
                id=f"{rel_path}:module",
                text=f"Module: {module}\nFile: {rel_path}\n\n{mod_doc}",
                source_path=rel_path,
                kind="module_doc",
                symbol=module,
                module=module,
                title=module,
                workflow=infer_workflow(rel_path, module, mod_doc),
                priority=priority_for(rel_path),
                metadata={"has_doc": True},
            )
        )

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith("_"):
                continue
            chunks.append(
                _make_chunk(
                    rel_path=rel_path, module=module, name=node.name,
                    qualifier=node.name, node=node,
                    source_lines=source_lines, kind="python_symbol",
                )
            )
        elif isinstance(node, ast.ClassDef):
            chunks.append(
                _make_chunk(
                    rel_path=rel_path, module=module, name=node.name,
                    qualifier=node.name, node=node,
                    source_lines=source_lines, kind="python_symbol",
                )
            )
            # public methods of the class
            for sub in node.body:
                if isinstance(
                    sub, (ast.FunctionDef, ast.AsyncFunctionDef)
                ):
                    if sub.name.startswith("_") and sub.name != "__init__":
                        continue
                    chunks.append(
                        _make_chunk(
                            rel_path=rel_path, module=module,
                            name=sub.name,
                            qualifier=f"{node.name}.{sub.name}",
                            node=sub, source_lines=source_lines,
                            kind="python_method",
                        )
                    )

    return chunks
