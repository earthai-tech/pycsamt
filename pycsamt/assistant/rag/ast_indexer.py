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

from .config import (
    infer_workflow,
    priority_for,
    workflow_for_path,
)
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
        return (
            f"class {name}({', '.join(bases)})" if bases else f"class {name}"
        )
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
    code = "\n".join(source_lines[start - 1 : end])
    if len(code) > _MAX_CODE_CHARS:
        code = code[:_MAX_CODE_CHARS] + "\n# … (truncated)"
    return code


def _unparse(node: ast.AST | None) -> str | None:
    if node is None:
        return None
    try:
        return ast.unparse(node)
    except Exception:  # noqa: BLE001
        return None


def _extract_params(node: ast.AST) -> list[dict[str, str | None]]:
    """Structured parameters (name / annotation / default) for a def node.

    Powers signature-aware answers and code generation — the caller gets
    the exact call surface, not just a flat signature string. ``self`` /
    ``cls`` are dropped; ``*args`` / ``**kwargs`` are kept with markers.
    Returns ``[]`` for classes and anything without an arg list.
    """
    a = getattr(node, "args", None)
    if not isinstance(a, ast.arguments):
        return []
    out: list[dict[str, str | None]] = []
    positional = list(getattr(a, "posonlyargs", [])) + list(a.args)
    defaults = list(a.defaults)
    first_default = len(positional) - len(defaults)
    for i, arg in enumerate(positional):
        if arg.arg in ("self", "cls"):
            continue
        default = (
            _unparse(defaults[i - first_default])
            if i - first_default >= 0
            else None
        )
        out.append(
            {
                "name": arg.arg,
                "annotation": _unparse(arg.annotation),
                "default": default,
            }
        )
    if a.vararg:
        out.append(
            {"name": f"*{a.vararg.arg}", "annotation": None, "default": None}
        )
    for arg, dflt in zip(a.kwonlyargs, a.kw_defaults):
        out.append(
            {
                "name": arg.arg,
                "annotation": _unparse(arg.annotation),
                "default": _unparse(dflt),
            }
        )
    if a.kwarg:
        out.append(
            {"name": f"**{a.kwarg.arg}", "annotation": None, "default": None}
        )
    return out


def _resolve_from_import(
    module: str, node: ast.ImportFrom, *, is_package: bool
) -> str | None:
    """Absolute dotted module for a ``from ... import`` (handles relatives)."""
    if not node.level:
        return node.module
    parts = module.split(".")
    # In a package's __init__, ``from .x`` is relative to the package itself;
    # in a submodule it is relative to the parent package.
    base = parts if is_package else parts[:-1]
    base = base[: len(base) - (node.level - 1)]
    if not base:
        return None
    return (
        f"{'.'.join(base)}.{node.module}" if node.module else ".".join(base)
    )


def _collect_imports(
    tree: ast.AST, module: str, *, is_package: bool
) -> tuple[dict[str, str], dict[str, str]]:
    """Map local names to the pyCSAMT symbols/modules they refer to.

    Returns ``(name_map, alias_map)``:

    * ``name_map`` — ``from pycsamt.emtools.ss import correct_ss_ama`` maps
      ``correct_ss_ama`` → the fully-qualified symbol.
    * ``alias_map`` — ``import pycsamt.emtools.ss as ss`` maps ``ss`` → the
      module, so ``ss.correct_ss_ama(...)`` resolves too.

    Only pyCSAMT targets are kept; stdlib/third-party imports are dropped so
    they can never become cross-reference edges.
    """
    name_map: dict[str, str] = {}
    alias_map: dict[str, str] = {}
    for sub in ast.walk(tree):
        if isinstance(sub, ast.ImportFrom):
            target = _resolve_from_import(module, sub, is_package=is_package)
            if not target or not target.startswith("pycsamt"):
                continue
            for a in sub.names:
                if a.name == "*":
                    continue
                name_map[a.asname or a.name] = f"{target}.{a.name}"
        elif isinstance(sub, ast.Import):
            for a in sub.names:
                if a.name.startswith("pycsamt") and a.asname:
                    alias_map[a.asname] = a.name
    return name_map, alias_map


def _is_exception_name(name: str) -> bool:
    return name.endswith(("Error", "Exception", "Warning"))


def _extract_refs(
    node: ast.AST,
    *,
    module: str,
    name_map: dict[str, str],
    alias_map: dict[str, str],
    local_defs: set[str],
    limit: int = 16,
) -> list[str]:
    """Fully-qualified pyCSAMT symbols *called* inside a function body.

    Resolution is import-aware rather than name-based, which is what keeps
    the graph honest: a bare ``all(...)`` / ``np.log(...)`` / ``df.unique()``
    is a builtin or a foreign method and yields nothing, while
    ``estimate_ss_ama(...)`` resolves through the module's imports (or its
    own top-level defs) to a real symbol. Raised exceptions are excluded —
    they are not part of a call surface worth suggesting.
    """
    refs: list[str] = []
    seen: set[str] = set()

    def _add(qualified: str, leaf: str) -> bool:
        if _is_exception_name(leaf) or len(leaf) < 3 or leaf.startswith("_"):
            return False
        if qualified in seen:
            return False
        seen.add(qualified)
        refs.append(qualified)
        return len(refs) >= limit

    for sub in ast.walk(node):
        if not isinstance(sub, ast.Call):
            continue
        fn = sub.func
        if isinstance(fn, ast.Name):
            if fn.id in name_map:
                if _add(name_map[fn.id], fn.id):
                    break
            elif fn.id in local_defs:
                if _add(f"{module}.{fn.id}", fn.id):
                    break
        elif isinstance(fn, ast.Attribute) and isinstance(fn.value, ast.Name):
            base = alias_map.get(fn.value.id) or name_map.get(fn.value.id)
            if base and _add(f"{base}.{fn.attr}", fn.attr):
                break
    return refs


def _make_chunk(
    *,
    rel_path: str,
    module: str,
    name: str,
    qualifier: str,
    node: ast.AST,
    source_lines: list[str],
    kind: str,
    refs_ctx: tuple[dict[str, str], dict[str, str], set[str]] | None = None,
) -> RAGChunk:
    doc = ast.get_docstring(node) or ""
    sig = _signature(node)
    code = _source_segment(source_lines, node)
    symbol = f"{module}.{qualifier}" if qualifier else module
    start = getattr(node, "lineno", 1)
    end = getattr(node, "end_lineno", start)
    is_func = isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
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
        workflow=infer_workflow(rel_path, symbol, doc)
        or workflow_for_path(rel_path),
        priority=priority_for(rel_path),
        metadata={
            "name": name,
            "signature": sig,
            "start_line": start,
            "end_line": end,
            "has_doc": bool(doc),
            # Tier 3: signature-aware + cross-reference metadata (funcs only;
            # classes carry a signature but no params/call-refs of their own).
            "params": _extract_params(node) if is_func else [],
            "returns": _unparse(getattr(node, "returns", None)),
            "refs": (
                _extract_refs(
                    node,
                    module=module,
                    name_map=refs_ctx[0],
                    alias_map=refs_ctx[1],
                    local_defs=refs_ctx[2],
                )
                if is_func and refs_ctx is not None
                else []
            ),
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

    # Import-aware cross-reference context, resolved once per file.
    is_package = rel_path.endswith("/__init__.py")
    name_map, alias_map = _collect_imports(
        tree, module, is_package=is_package
    )
    local_defs = {
        n.name
        for n in tree.body
        if isinstance(
            n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        )
    }
    refs_ctx = (name_map, alias_map, local_defs)

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
                workflow=infer_workflow(rel_path, module, mod_doc)
                or workflow_for_path(rel_path),
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
                    rel_path=rel_path,
                    module=module,
                    name=node.name,
                    qualifier=node.name,
                    node=node,
                    source_lines=source_lines,
                    kind="python_symbol",
                    refs_ctx=refs_ctx,
                )
            )
        elif isinstance(node, ast.ClassDef):
            chunks.append(
                _make_chunk(
                    rel_path=rel_path,
                    module=module,
                    name=node.name,
                    qualifier=node.name,
                    node=node,
                    source_lines=source_lines,
                    kind="python_symbol",
                    refs_ctx=refs_ctx,
                )
            )
            # public methods of the class
            for sub in node.body:
                if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    if sub.name.startswith("_") and sub.name != "__init__":
                        continue
                    chunks.append(
                        _make_chunk(
                            rel_path=rel_path,
                            module=module,
                            name=sub.name,
                            qualifier=f"{node.name}.{sub.name}",
                            node=sub,
                            source_lines=source_lines,
                            kind="python_method",
                            refs_ctx=refs_ctx,
                        )
                    )

    return chunks
