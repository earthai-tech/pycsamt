# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Resolve checkout APIs without importing their modules."""

from __future__ import annotations

import ast
import hashlib
import inspect
from pathlib import Path


class StaticAPI:
    def __init__(self, root=None):
        self.root = Path(root or Path(__file__).resolve().parents[3]).resolve()
        self.cache = {}

    def module(self, name):
        if not (name == "pycsamt" or name.startswith("pycsamt.")) or not all(
            p.isidentifier() for p in name.split(".")
        ):
            return "unverifiable", None, None
        for path in [
            self.root.joinpath(*name.split(".")).with_suffix(".py"),
            self.root.joinpath(*name.split("."), "__init__.py"),
        ]:
            if not path.resolve().is_relative_to(self.root / "pycsamt"):
                continue
            if path.is_file():
                try:
                    if path.stat().st_size > 2_000_000:
                        return "unverifiable", None, None
                    if name not in self.cache:
                        with path.open("rb") as stream:
                            raw = stream.read(2_000_001)
                        if len(raw) > 2_000_000:
                            return "unverifiable", None, None
                        self.cache[name] = (
                            ast.parse(raw.decode("utf-8")),
                            {
                                "path": path.relative_to(self.root).as_posix(),
                                "sha256": hashlib.sha256(raw).hexdigest(),
                            },
                        )
                    tree, evidence = self.cache[name]
                    return "passed", tree, evidence
                except (OSError, UnicodeError, SyntaxError):
                    return "unverifiable", None, None
        return (
            ("failed" if (self.root / "pycsamt").is_dir() else "unverifiable"),
            None,
            None,
        )

    def resolve(self, qualified, depth=0):
        if depth > 8:
            return "unverifiable", None, None
        parts = qualified.split(".")
        for cut in range(len(parts), 0, -1):
            module = ".".join(parts[:cut])
            state, tree, evidence = self.module(module)
            if state == "unverifiable":
                return state, None, None
            if tree is None:
                continue
            if cut == len(parts):
                return "passed", tree, evidence
            node = tree
            for offset, name in enumerate(parts[cut:]):
                found = next(
                    (
                        n
                        for n in getattr(node, "body", [])
                        if isinstance(
                            n,
                            (
                                ast.FunctionDef,
                                ast.AsyncFunctionDef,
                                ast.ClassDef,
                            ),
                        )
                        and n.name == name
                    ),
                    None,
                )
                if found is not None:
                    node = found
                    continue
                if not isinstance(node, ast.Module):
                    return "unverifiable", None, evidence
                for statement in tree.body:
                    if isinstance(statement, ast.ImportFrom):
                        for alias in statement.names:
                            if (alias.asname or alias.name) == name:
                                package = (
                                    module.split(".")
                                    if evidence["path"].endswith(
                                        "/__init__.py"
                                    )
                                    else module.split(".")[:-1]
                                )
                                prefix = (
                                    package[
                                        : len(package) - statement.level + 1
                                    ]
                                    if statement.level
                                    else []
                                )
                                target = ".".join(
                                    prefix
                                    + (
                                        [statement.module]
                                        if statement.module
                                        else []
                                    )
                                    + [alias.name]
                                    + parts[cut + offset + 1 :]
                                )
                                return self.resolve(target, depth + 1)
                    if isinstance(statement, (ast.Assign, ast.AnnAssign)):
                        targets = (
                            statement.targets
                            if isinstance(statement, ast.Assign)
                            else [statement.target]
                        )
                        if any(
                            isinstance(t, ast.Name) and t.id == name
                            for t in targets
                        ):
                            return (
                                (
                                    "passed"
                                    if cut + offset + 1 == len(parts)
                                    else "unverifiable"
                                ),
                                statement,
                                evidence,
                            )
                        if isinstance(statement.value, ast.Dict):
                            for key, value in zip(
                                statement.value.keys, statement.value.values
                            ):
                                if (
                                    isinstance(key, ast.Constant)
                                    and key.value == name
                                    and isinstance(value, ast.Constant)
                                    and isinstance(value.value, str)
                                    and value.value.startswith(".")
                                ):
                                    return self.resolve(
                                        module + value.value + "." + name,
                                        depth + 1,
                                    )
                conditional = any(
                    isinstance(n, (ast.ImportFrom, ast.Import))
                    and any(
                        (a.asname or a.name) == name or a.name == "*"
                        for a in n.names
                    )
                    for statement in tree.body
                    if isinstance(statement, (ast.If, ast.Try))
                    for n in ast.walk(statement)
                )
                if conditional:
                    return "unverifiable", None, evidence
                dynamic = any(
                    isinstance(n, ast.FunctionDef) and n.name == "__getattr__"
                    for n in tree.body
                )
                return (
                    (
                        "unverifiable"
                        if dynamic and module != "pycsamt"
                        else "failed"
                    ),
                    None,
                    evidence,
                )
            return "passed", node, {**evidence, "line": node.lineno}
        return "failed", None, None


def signature(node, *, bound=False):
    """Conservative signature: dynamic/decorated shapes remain unverifiable."""
    if isinstance(node, ast.ClassDef):
        init = next(
            (
                n
                for n in node.body
                if isinstance(n, ast.FunctionDef) and n.name == "__init__"
            ),
            None,
        )
        if init:
            return signature(init, bound=True)
        if not node.bases and any(
            isinstance(d, ast.Name) and d.id == "dataclass"
            for d in node.decorator_list
        ):
            params = []
            for field in node.body:
                if isinstance(field, ast.AnnAssign) and isinstance(
                    field.target, ast.Name
                ):
                    if (
                        isinstance(field.annotation, ast.Subscript)
                        and getattr(field.annotation.value, "id", "")
                        == "ClassVar"
                    ):
                        continue
                    if isinstance(field.value, ast.Call) and any(
                        k.arg in {"init", "kw_only"}
                        for k in field.value.keywords
                    ):
                        return None
                    params.append(
                        inspect.Parameter(
                            field.target.id,
                            inspect.Parameter.POSITIONAL_OR_KEYWORD,
                            default=inspect.Parameter.empty
                            if field.value is None
                            else object(),
                        )
                    )
            try:
                return inspect.Signature(params)
            except ValueError:
                return None
        return None
    if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        return None
    if any(
        not isinstance(d, ast.Name)
        or d.id not in {"staticmethod", "classmethod"}
        for d in node.decorator_list
    ):
        return None
    args = node.args
    if any(
        isinstance(d, ast.Name) and d.id == "staticmethod"
        for d in node.decorator_list
    ):
        bound = False
    positional = args.posonlyargs + args.args
    required = len(positional) - len(args.defaults)
    params = []
    for i, arg in enumerate(positional):
        if bound and i == 0 and arg.arg in {"self", "cls"}:
            continue
        params.append(
            inspect.Parameter(
                arg.arg,
                inspect.Parameter.POSITIONAL_ONLY
                if i < len(args.posonlyargs)
                else inspect.Parameter.POSITIONAL_OR_KEYWORD,
                default=inspect.Parameter.empty if i < required else object(),
            )
        )
    if args.vararg:
        params.append(
            inspect.Parameter(
                args.vararg.arg, inspect.Parameter.VAR_POSITIONAL
            )
        )
    for arg, default in zip(args.kwonlyargs, args.kw_defaults):
        params.append(
            inspect.Parameter(
                arg.arg,
                inspect.Parameter.KEYWORD_ONLY,
                default=inspect.Parameter.empty
                if default is None
                else object(),
            )
        )
    if args.kwarg:
        params.append(
            inspect.Parameter(args.kwarg.arg, inspect.Parameter.VAR_KEYWORD)
        )
    return inspect.Signature(params)
