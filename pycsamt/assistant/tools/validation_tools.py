# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.tools.validation_tools
========================================

Deterministic checks on assistant-generated code — the "Validator" stage
of the assistant flow (RAG → tool → **validation** → answer).

:func:`validate_generated_code` parses the script and verifies that every
``pycsamt`` symbol it imports actually exists, so hallucinated APIs
(e.g. ``from pycsamt import static_shift``) are caught before the code
ever runs.
"""

from __future__ import annotations

import ast
import importlib
from pathlib import Path
from typing import Any

__all__ = ["validate_generated_code"]


def _pycsamt_imports(tree: ast.AST) -> list[tuple[str, list[str]]]:
    """Collect ``(module, names)`` for pycsamt imports in *tree*.

    ``names`` is empty for plain ``import pycsamt.x`` statements.
    """
    out: list[tuple[str, list[str]]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if mod == "pycsamt" or mod.startswith("pycsamt."):
                out.append((mod, [a.name for a in node.names]))
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name == "pycsamt" or a.name.startswith("pycsamt."):
                    out.append((a.name, []))
    return out


def _check_symbol(module: str, name: str) -> str | None:
    """Return an error string if ``module.name`` doesn't exist, else None.

    A failure to import the module is reported as *unverifiable* (None)
    rather than a hard error — we only flag a symbol we can prove absent.
    """
    try:
        mod = importlib.import_module(module)
    except Exception:  # noqa: BLE001 — can't verify (heavy/optional dep)
        return None
    if hasattr(mod, name):
        return None
    # maybe it's a submodule (e.g. from pycsamt import emtools)
    try:
        importlib.import_module(f"{module}.{name}")
        return None
    except Exception:  # noqa: BLE001
        return f"{module!r} has no attribute or submodule {name!r}"


def validate_generated_code(code: str) -> dict[str, Any]:
    """Validate a generated pyCSAMT script.

    Returns a report dict::

        {
          "ok": bool,            # no syntax errors and no missing symbols
          "syntax_ok": bool,
          "errors": [...],       # hard problems (syntax, missing symbols)
          "warnings": [...],     # unverifiable imports
          "checked": [...],      # pycsamt symbols verified present
        }
    """
    report: dict[str, Any] = {
        "ok": True,
        "syntax_ok": True,
        "errors": [],
        "warnings": [],
        "checked": [],
    }

    if not (code or "").strip():
        report.update(ok=False, syntax_ok=False)
        report["errors"].append("Empty code.")
        return report

    try:
        tree = ast.parse(code)
    except SyntaxError as exc:
        report.update(ok=False, syntax_ok=False)
        report["errors"].append(f"SyntaxError: {exc}")
        return report

    for module, names in _pycsamt_imports(tree):
        if not names:  # plain `import pycsamt.x`
            try:
                importlib.import_module(module)
                report["checked"].append(module)
            except Exception as exc:  # noqa: BLE001
                report["warnings"].append(
                    f"could not import {module!r}: {exc}"
                )
            continue
        for name in names:
            err = _check_symbol(module, name)
            if err:
                report["errors"].append(err)
            else:
                report["checked"].append(f"{module}.{name}")

    report["ok"] = report["syntax_ok"] and not report["errors"]
    return report


def validate_script_file(path: str | Path) -> dict[str, Any]:
    """Validate a script on disk (convenience wrapper)."""
    p = Path(path)
    try:
        code = p.read_text(encoding="utf-8")
    except OSError as exc:
        return {
            "ok": False,
            "syntax_ok": False,
            "errors": [f"cannot read {path!r}: {exc}"],
            "warnings": [],
            "checked": [],
        }
    return validate_generated_code(code)
