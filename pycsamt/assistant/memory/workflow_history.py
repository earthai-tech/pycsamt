# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.memory.workflow_history
=========================================

Append-only JSONL trace of executed workflows — observability for "what
did I run", debugging, and eval/telemetry. Each line is one
:class:`WorkflowRun`.

Default location: ``<root>/.pycsamt_rag/workflow_history.jsonl``.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

__all__ = ["WorkflowRun", "WorkflowHistory"]


@dataclass
class WorkflowRun:
    """One executed (or attempted) workflow."""

    workflow: str
    status: str = ""
    line: str | None = None
    path: str | None = None
    output_dir: str | None = None
    summary: str = ""
    n_figures: int = 0
    timestamp: str = field(default_factory=lambda: time.strftime("%Y-%m-%dT%H:%M:%S"))

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _default_path(root: Path | str | None = None) -> Path:
    from ..rag.ingest import repo_root

    base = Path(root) if root is not None else repo_root()
    return base / ".pycsamt_rag" / "workflow_history.jsonl"


class WorkflowHistory:
    """Append-only history of workflow runs (JSONL)."""

    def __init__(self, path: str | Path | None = None):
        self.path = Path(path) if path is not None else _default_path()

    def record(self, run: WorkflowRun | dict[str, Any]) -> None:
        """Append a run to the history file (creates it if needed)."""
        entry = run.to_dict() if isinstance(run, WorkflowRun) else dict(run)
        entry.setdefault("timestamp", time.strftime("%Y-%m-%dT%H:%M:%S"))
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    def all(self) -> list[dict[str, Any]]:
        """All recorded runs (oldest first)."""
        if not self.path.is_file():
            return []
        out: list[dict[str, Any]] = []
        with self.path.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        out.append(json.loads(line))
                    except json.JSONDecodeError:
                        continue
        return out

    def recent(self, n: int = 10) -> list[dict[str, Any]]:
        """The last *n* runs (most recent last)."""
        return self.all()[-n:]

    def clear(self) -> None:
        """Delete the history file."""
        if self.path.is_file():
            self.path.unlink()

    @classmethod
    def default(cls, root: Path | str | None = None) -> WorkflowHistory:
        return cls(_default_path(root))
