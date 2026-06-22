# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.memory.project_state
======================================

Disk-persisted, per-project memory: which survey lines the user has
worked on, their last-used parameters, and preferences. Survives across
sessions so the assistant can offer sensible defaults and recall context.

Stored as JSON (default: ``<root>/.pycsamt_rag/project_state.json``).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

__all__ = ["ProjectState"]


def _default_path(root: Path | str | None = None) -> Path:
    from ..rag.ingest import repo_root
    base = Path(root) if root is not None else repo_root()
    return base / ".pycsamt_rag" / "project_state.json"


class ProjectState:
    """Persistent per-project memory (lines, params, preferences)."""

    def __init__(self, path: str | Path | None = None):
        self.path = Path(path) if path is not None else _default_path()
        self.data: dict[str, Any] = {
            "lines": {},          # line -> last-used info
            "preferences": {},    # arbitrary user prefs
            "updated": None,
        }
        if self.path.is_file():
            try:
                loaded = json.loads(
                    self.path.read_text(encoding="utf-8")
                )
                if isinstance(loaded, dict):
                    self.data.update(loaded)
            except (OSError, json.JSONDecodeError):
                pass

    # ── lines ────────────────────────────────────────────────────────────────
    def remember_line(self, line: str, info: dict[str, Any]) -> None:
        """Record that *line* was used, with its resolved info/params."""
        entry = dict(info)
        entry["last_used"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        self.data.setdefault("lines", {})[line] = entry

    def lines_used(self) -> list[str]:
        """Lines the user has worked on, most-recent first."""
        lines = self.data.get("lines", {})
        return sorted(
            lines, key=lambda k: lines[k].get("last_used", ""),
            reverse=True,
        )

    def line_info(self, line: str) -> dict[str, Any] | None:
        return self.data.get("lines", {}).get(line)

    # ── preferences ────────────────────────────────────────────────────────
    def set_preference(self, key: str, value: Any) -> None:
        self.data.setdefault("preferences", {})[key] = value

    def get_preference(self, key: str, default: Any = None) -> Any:
        return self.data.get("preferences", {}).get(key, default)

    # ── persistence ────────────────────────────────────────────────────────
    def save(self) -> Path:
        self.data["updated"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(
            json.dumps(self.data, indent=2), encoding="utf-8"
        )
        return self.path

    @classmethod
    def default(cls, root: Path | str | None = None) -> "ProjectState":
        return cls(_default_path(root))
