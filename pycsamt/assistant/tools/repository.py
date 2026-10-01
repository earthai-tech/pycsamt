# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Bounded developer evidence, separate from the science retrieval index.

Every search reads current files. No source imports, execution, embeddings,
network access, or model-selected filesystem roots are involved.
"""

from __future__ import annotations

import ast
import hashlib
import math
import os
import re
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from functools import lru_cache
from importlib import metadata
from pathlib import Path

_ROOTS = (
    "pycsamt",
    "docs/source",
    "docs/examples",
    "examples",
    "assistant_recipes",
)
_SKIP = {
    "__pycache__",
    "build",
    "dist",
    "generated",
    "assets",
    "resources",
    "_static",
    "images",
    "node_modules",
}
_SECRET = re.compile(
    r"(?:secret|credential|password|token|api[_-]?key|private[_-]?key)", re.I
)
_STOP = set(
    "the and for how what where why does this that from with which explain implementation implemented source code function method pycsamt please find show tests test example examples"
    " not about into only can will would should could are was were has have its any all use used using".split()
)


def is_developer_question(text: str) -> bool:
    """Require implementation intent; ordinary science usage stays in RAG."""
    intent = re.search(
        r"\b(?:where|how|why|explain|find|show|inspect|trace|signature|docstring|implementation|implemented|source|tests?)\b",
        text,
        re.I,
    )
    scope = re.search(
        r"(?:pycsamt[/\\.]app|agent[_ ]master|callbacks?|repository|checkout|source code|implementation|unit tests?|test_[a-z_]+|_dispatch_[a-z_]+|pycsamt[/\\.](?:agents|assistant))",
        text,
        re.I,
    )
    action = re.match(
        r"\s*(?:write|generate|create|run|execute|apply|change|modify|delete)\b",
        text,
        re.I,
    )
    if action:
        return False
    if intent and scope:
        return True
    # A question naming an assistant-infrastructure function
    # ("does validate_generated_code prove ...?") is about the implementation.
    names = set(re.findall(r"\b[a-z]+(?:_[a-z0-9]+)+\b", text))
    asks = re.search(r"\b(?:does|do|is|are|can|what|how|why|where|which)\b|\?", text, re.I)
    return bool(asks and names & _infrastructure_names())


@lru_cache(maxsize=1)
def _infrastructure_names() -> frozenset[str]:
    """Function/class names defined in Agent Master's own infrastructure."""
    package = Path(__file__).resolve().parents[2]
    files = [
        *(package / "assistant").rglob("*.py"),
        *(package / "app" / "agent_master").rglob("*.py"),
        *(package / "agents").glob("_*.py"),
        *(package / "agents" / n for n in ("code_gen.py", "router.py", "package_qa.py")),
    ]
    names = set()
    for path in files:
        if "tests" in path.parts or not path.is_file():
            continue
        try:
            source = path.read_text(encoding="utf-8")
        except (OSError, UnicodeError):
            continue
        names.update(re.findall(r"^\s*(?:async\s+)?def\s+([a-z]\w*_\w+)", source, re.M))
    return frozenset(names)


class RepositoryTools:
    """Read only approved source roots; limits are hard upper bounds."""

    def __init__(self, root: Path | None = None):
        self.root = (root or Path(__file__).resolve().parents[3]).resolve()

    def _path(self, relative: str, include_tests: bool = False) -> Path:
        path = (self.root / relative).resolve()
        rel = path.relative_to(self.root)
        if not any(path.is_relative_to(self.root / p) for p in _ROOTS):
            raise ValueError("Outside approved source roots")
        if any(
            p.startswith(".") or p in _SKIP or _SECRET.search(p)
            for p in rel.parts
        ):
            raise ValueError("Excluded directory")
        if path.suffix not in {".py", ".rst", ".md"} or _SECRET.search(
            path.name
        ):
            raise ValueError("Excluded file type or sensitive filename")
        if not include_tests and (
            "tests" in rel.parts or path.name.startswith("test_")
        ):
            raise ValueError("Tests require explicit inclusion")
        return path

    def _read(
        self, relative: str, include_tests: bool = False
    ) -> tuple[Path, str, str]:
        path = self._path(relative, include_tests)
        if path.stat().st_size > 1_000_000:
            raise ValueError("Source exceeds one-megabyte limit")
        with path.open("rb") as stream:
            raw = stream.read(1_000_001)
        if len(raw) > 1_000_000:
            raise ValueError("Source exceeds one-megabyte limit")
        text = raw.decode("utf-8")
        # Mask quoted credential assignments without deleting definitions or
        # changing line numbers. In particular api_key=None is a signature,
        # not a secret, and must survive static inspection.
        assignment = re.compile(
            r"(?im)(\b\w*(?:secret|password|token|api_key|private_key)\w*\s*(?::[^\n=]+)?=\s*)"
            r"(?:[rubf]{0,2})(\"\"\"[\s\S]*?\"\"\"|'''[\s\S]*?'''|\"(?:\\.|[^\"\\\n])*\"|'(?:\\.|[^'\\\n])*')"
        )
        text = assignment.sub(
            lambda m: (
                m[1]
                + ('"""[credential omitted]' + "\n" * m[2].count("\n") + '"""')
            ),
            text,
        )
        text = re.sub(
            r"-----BEGIN [^-]*PRIVATE KEY-----.*?-----END [^-]*PRIVATE KEY-----",
            lambda m: "[private key omitted]" + "\n" * m[0].count("\n"),
            text,
            flags=re.S,
        )
        return path, text, hashlib.sha256(raw).hexdigest()

    def read_source(
        self,
        relative: str,
        *,
        start: int = 1,
        lines: int = 40,
        include_tests: bool = False,
    ) -> dict:
        path, text, digest = self._read(relative, include_tests)
        start = max(1, int(start))
        excerpt = "\n".join(
            text.splitlines()[
                start - 1 : start - 1 + min(80, max(1, int(lines)))
            ]
        )[:4000]
        rel = path.relative_to(self.root).as_posix()
        return {
            "path": rel,
            "line": start,
            "excerpt": excerpt,
            "sha256": digest,
            "kind": "test: expected behavior example, not public API documentation"
            if "tests" in path.parts or path.name.startswith("test_")
            else "source evidence",
        }

    def inspect_symbol(
        self, relative: str, symbol: str, *, include_tests: bool = False
    ) -> dict | None:
        _, source, _ = self._read(relative, include_tests)
        if not relative.endswith(".py"):
            return None
        try:
            tree = ast.parse(source)
        except SyntaxError:
            return None
        for node in ast.walk(tree):
            if (
                isinstance(
                    node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
                )
                and node.name == symbol
            ):
                card = self.read_source(
                    relative,
                    start=node.lineno,
                    lines=30,
                    include_tests=include_tests,
                )
                card.update(
                    symbol=node.name,
                    signature=node.name
                    + (
                        "(" + ast.unparse(node.args) + ")"
                        if hasattr(node, "args")
                        else " [class]"
                    ),
                    docstring=(ast.get_docstring(node) or "")[:1200],
                )
                return card
        return None

    def provenance(self) -> dict:
        checkout_version = "unknown"
        try:
            project = (self.root / "pyproject.toml").read_text(
                encoding="utf-8"
            )[:20000]
            match = re.search(r'^version\s*=\s*"([^"]+)"', project, re.M)
            if match:
                checkout_version = match[1]
        except OSError:
            pass
        try:
            version = metadata.version("pycsamt")
        except metadata.PackageNotFoundError:
            version = "not installed"
        loaded = getattr(sys.modules.get("pycsamt"), "__file__", None)
        matches = bool(
            loaded and Path(loaded).resolve().parent == self.root / "pycsamt"
        )
        try:
            revision = (
                subprocess.run(
                    ["git", "-C", str(self.root), "rev-parse", "HEAD"],
                    capture_output=True,
                    text=True,
                    timeout=2,
                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
                ).stdout.strip()
                or "unknown"
            )
        except (OSError, subprocess.TimeoutExpired):
            revision = "unknown"
        return {
            "root": str(self.root),
            "revision": revision,
            "installed_version": version,
            "checkout_version": checkout_version,
            "loaded_package_path": loaded,
            "checkout_matches_loaded_package": matches,
            "freshness": "live filesystem read; SHA-256 identifies each inspected file, including uncommitted changes",
            "read_at": datetime.now(timezone.utc).isoformat(),
            "limitation": "Checkout source may differ from the loaded/installed package; answers describe this checkout."
            if not matches or checkout_version != version
            else "Installed distribution version does not certify uncommitted checkout contents.",
        }

    def search(
        self, query: str, *, include_tests: bool = False, limit: int = 4
    ) -> dict:
        started = time.monotonic()
        terms = set(re.findall(r"[a-z_][a-z_0-9]{2,}", query.lower())) - _STOP
        if "agent master" in query.lower():
            terms.add("agent_master")
        # CamelCase names ("CodeGenerationAgent") are identifiers as much as
        # snake_case ones; common words must not outrank them.
        camel = {
            m.lower()
            for m in re.findall(r"\b[A-Z][a-z0-9]+(?:[A-Z][a-z0-9]+)+\b", query)
        } & terms
        paths, truncated = [], False
        for approved in _ROOTS:
            for directory, dirs, files in os.walk(
                self.root / approved, followlinks=False
            ):
                dirs[:] = sorted(
                    d
                    for d in dirs
                    if not d.startswith(".")
                    and d not in _SKIP
                    and (include_tests or d != "tests")
                    and (Path(directory) / d)
                    .resolve()
                    .is_relative_to(self.root / approved)
                )
                for name in sorted(files):
                    relative = (
                        (Path(directory) / name)
                        .relative_to(self.root)
                        .as_posix()
                    )
                    try:
                        self._path(relative, include_tests)
                    except ValueError:
                        continue
                    paths.append(relative)
                    if len(paths) >= 6000 or time.monotonic() - started > 5:
                        truncated = True
                        break
                if truncated:
                    break
            if truncated:
                break
        paths.sort(
            key=lambda p: (
                -sum(t in p.lower() for t in terms),
                -int(
                    include_tests
                    and ("/tests/" in p or Path(p).name.startswith("test_"))
                ),
                -int(p.endswith(".py")),
                p,
            )
        )
        candidates, size, scanned = [], 0, 0
        for relative in paths:
            if size >= 32_000_000 or time.monotonic() - started > 10:
                truncated = True
                break
            try:
                _, source, digest = self._read(relative, include_tests)
            except (OSError, ValueError, UnicodeError):
                continue
            scanned += 1
            size += len(source.encode("utf-8"))
            lower = source.lower()
            matched = {
                t
                for t in terms
                if re.search(r"(?<!\w)" + re.escape(t) + r"(?!\w)", lower)
                or t in relative.lower()
            }
            if not matched or not terms:
                continue
            # Require specific identifiers when supplied; do not substitute a
            # generic app hit for a missing named symbol.
            identifiers = {
                t for t in terms if "_" in t and t != "agent_master"
            } | camel
            if identifiers and not identifiers.intersection(matched):
                continue
            candidates.append((relative, matched, source, digest))
        # Rare terms ("rag", "cancel") identify the answer; words present
        # in most files ("app", "answer", "callbacks") must not outrank them.
        frequency = Counter(t for _, matched, _, _ in candidates for t in matched)
        weight = {t: math.log(1 + len(candidates) / n) for t, n in frequency.items()}
        # "Where is X defined / which code handles Y" asks for source, not prose.
        wants_source = re.search(
            r"\b(?:defined|implemented|handles?|handled|code|function|class)\b",
            query, re.I,
        )
        ranked = []
        for relative, matched, source, digest in candidates:
            score = sum(weight[t] for t in matched) + sum(
                2 * weight[t] for t in matched if t in relative.lower()
            )
            if wants_source and not relative.endswith(".py"):
                score *= 0.5
            ranked.append((score, relative, matched, source, digest))
        ranked.sort(key=lambda h: (-h[0], h[1]))
        hits = []
        for score, relative, matched, source, digest in ranked[:24]:
            best_line, best_score, symbol = 1, -1.0, ""
            for i, line in enumerate(source.splitlines(), 1):
                low = line.lower()
                value = sum(weight[t] for t in matched if t in low)
                definition = re.match(
                    r"\s*(?:async )?(?:def|class)\s+(\w+)", line
                )
                if definition and definition[1].lower() in matched:
                    value += 20
                    score += 5
                if value > best_score:
                    best_line, best_score = i, value
                    symbol = definition[1] if definition else ""
            hits.append((score, relative, best_line, symbol, digest))
        cards = []
        for _, relative, line, symbol, digest in sorted(
            hits, key=lambda h: (-h[0], h[1])
        )[: min(6, max(1, limit))]:
            try:
                card = (
                    self.inspect_symbol(
                        relative, symbol, include_tests=include_tests
                    )
                    if symbol
                    else None
                )
                card = card or self.read_source(
                    relative,
                    start=max(1, line - 2),
                    lines=25,
                    include_tests=include_tests,
                )
                card["changed_during_search"] = digest != card["sha256"]
                cards.append(card)
            except (ValueError, OSError, UnicodeError):
                truncated = True
        return {
            "scope": "developer",
            "sources": cards,
            "provenance": self.provenance(),
            "scanned_files": scanned,
            "truncated": truncated,
            "limitation": "Bounded search reached its limit; narrow the question to a path or symbol."
            if truncated
            else "Search is lexical; absent evidence is not proof that an implementation does not exist.",
        }


def evidence_text(result: dict) -> str:
    """Render citable source evidence; no inference of successful execution."""
    if not result["sources"]:
        return (
            "No matching developer evidence found in the approved source scope. "
            + result["limitation"]
        )
    parts = [
        "Developer source evidence (read only; not executed). Tests describe expected behavior, not public API contracts."
    ]
    for i, card in enumerate(result["sources"], 1):
        parts.append(
            f"[{i}] {card['path']}:{card['line']} — {card['kind']}\n{card.get('signature', '')}\n```python\n{card['excerpt'][:1200]}\n```"
        )
    info = result["provenance"]
    parts.append(
        f"Checkout revision: {info['revision']}; checkout version: {info['checkout_version']}; installed version: {info['installed_version']}. Sources read live at {info['read_at']}. {info['limitation']} {result['limitation']}"
    )
    return "\n\n".join(parts)
