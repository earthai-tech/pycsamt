# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Bounded, static science API inspection for script generation."""

from __future__ import annotations

import ast
import re
from pathlib import Path

_SCIENCE_ROOTS = {"emtools", "site", "seg", "z", "forward", "ai", "models"}
_WORKFLOW_APIS = {
    "qc": ("pycsamt.emtools.qc.build_qc_table",),
    "static_shift": (
        "pycsamt.emtools.ss.estimate_ss_ama",
        "pycsamt.emtools.ss.correct_ss_ama",
    ),
    "forward": (
        "pycsamt.forward.synthetic.LayeredModel",
        "pycsamt.forward.em1d.MT1DForward",
    ),
}


def inspect_science_api(
    symbol: str, *, root: Path | None = None
) -> dict | None:
    """Read an exact module symbol without importing or executing its module."""
    parts = symbol.split(".")
    if (
        len(parts) < 4
        or parts[0] != "pycsamt"
        or parts[1] not in _SCIENCE_ROOTS
    ):
        return None
    if not all(p.isidentifier() for p in parts) or "tests" in parts:
        return None
    root = (root or Path(__file__).resolve().parents[3]).resolve()
    # Try module.function, then module.Class.method; never search arbitrary files.
    for cut in range(len(parts) - 1, 2, -1):
        path = root.joinpath(*parts[:cut]).with_suffix(".py")
        try:
            resolved = path.resolve()
            if not resolved.is_relative_to(root / "pycsamt" / parts[1]):
                return None
            if resolved.stat().st_size > 2_000_000:
                return None
            source = resolved.read_text(encoding="utf-8")
        except (OSError, UnicodeError):
            continue
        try:
            node = ast.parse(source)
            for name in parts[cut:]:
                node = next(
                    n
                    for n in node.body
                    if isinstance(
                        n,
                        (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
                    )
                    and n.name == name
                )
        except (SyntaxError, StopIteration, AttributeError):
            continue
        signature = symbol
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            signature += "(" + ast.unparse(node.args) + ")"
        else:
            init = next(
                (
                    n
                    for n in node.body
                    if isinstance(n, ast.FunctionDef) and n.name == "__init__"
                ),
                None,
            )
            fields = [
                ast.unparse(n).split(" = ", 1)[0]
                for n in node.body
                if isinstance(n, ast.AnnAssign)
            ]
            dataclass = any(
                isinstance(d, ast.Name) and d.id == "dataclass"
                for d in node.decorator_list
            )
            if init:
                signature += (
                    "(" + ast.unparse(init.args).removeprefix("self, ") + ")"
                )
            elif dataclass:
                signature += (
                    " [dataclass fields: "
                    + "; ".join(fields)
                    + "; consult documented defaults]"
                )
            else:
                signature += " [class; inspect constructor before calling]"
        methods = []
        if isinstance(node, ast.ClassDef):
            methods = [
                n.name + "(" + ast.unparse(n.args) + ")"
                for n in node.body
                if isinstance(n, ast.FunctionDef)
                and not n.name.startswith("_")
            ][:4]
        return {
            "symbol": symbol,
            "source_path": resolved.relative_to(root).as_posix(),
            "line": node.lineno,
            "signature": signature,
            "methods": methods,
            "doc": (ast.get_docstring(node) or "")[:500],
        }
    return None


def generation_api_evidence(
    request: str, workflow: str | None, candidates=(), *, root=None
) -> list[dict]:
    """Inspect explicit APIs first, then retrieved identities and workflow APIs."""
    explicit = re.findall(r"\bpycsamt(?:\.[A-Za-z_]\w*){2,}", request)
    names = [
        *explicit,
        *_WORKFLOW_APIS.get(workflow, ()),
        *candidates,
        "pycsamt.emtools._core.ensure_sites",
    ]
    cards = []
    for symbol in dict.fromkeys(names):
        if not isinstance(symbol, str):
            continue
        card = inspect_science_api(symbol, root=root)
        if card:
            cards.append(card)
        if len(cards) == 4:
            break
    return cards
