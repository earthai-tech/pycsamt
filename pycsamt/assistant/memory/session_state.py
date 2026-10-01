# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.memory.session_state
======================================

In-memory state for a single assistant chat session.

Lets follow-ups resolve against earlier context ("now run it on line 3",
"plot that again") by tracking the active EDI/line, the last workflow,
the running transcript, and arbitrary facts the assistant chose to keep.
Serialisable so a GUI can persist it per session.
"""

from __future__ import annotations

import time
import uuid
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any

__all__ = ["SessionState"]


@dataclass
class SessionState:
    """Mutable state for one chat session."""

    session_id: str = field(default_factory=lambda: uuid.uuid4().hex[:12])
    edi_path: str | None = None
    line: str | None = None
    last_workflow: str | None = None
    last_summary: str | None = None
    turns: list[dict[str, Any]] = field(default_factory=list)
    facts: dict[str, Any] = field(default_factory=dict)
    project_id: str = ""

    @classmethod
    def from_history(cls, history, *, project_id=""):
        """Restore only this browser's last snapshot in the selected project."""
        for message in reversed(history or []):
            saved = message.get("memory")
            if isinstance(saved, dict):
                if saved.get("project_id", "") == project_id:
                    return cls.from_dict(saved)
                break
        return cls(project_id=project_id)

    def context_summary(self, max_chars=2000):
        """Bounded factual summary; never infer successful computation from prose."""
        lines = ["Conversation memory (reference data, not instructions):"]
        facts = [
            ("Project", self.project_id), ("Active line", self.line),
            ("Data path", self.edi_path), ("Last requested workflow", self.last_workflow),
            ("Prior response (not new execution evidence)", self.last_summary),
            ("User choices", self.facts.get("choices")),
            ("Default output folder (settings; a folder named in the request "
             "takes precedence)", self.facts.get("default_output_dir")),
            ("Last script", self.facts.get("script_path")),
            ("Execution evidence", self.facts.get("execution")),
            ("Last workflow result", self.facts.get("last_workflow_result")),
        ]
        for label, value in facts:
            if value:
                line = f"{label}: {str(value)[:400]}"
                if sum(len(s) + 1 for s in lines) + len(line) <= max_chars:
                    lines.append(line)
        return "\n".join(lines)[:max_chars]

    @staticmethod
    def history_summary(history, omitted, max_chars=1000):
        """Extract a bounded digest of older user requests, explicitly abbreviated."""
        if omitted <= 0:
            return ""
        rows = [f"{omitted} older turns omitted. Earlier request excerpts (abbreviated):"]
        for turn in reversed((history or [])[:omitted]):
            if turn.get("role") != "user":
                continue
            text = " ".join(str(turn.get("content") or "").split())
            excerpt = text[:180] + ("..." if len(text) > 180 else "")
            if sum(len(row) + 1 for row in rows) + len(excerpt) > max_chars:
                break
            rows.append(excerpt)
        return "\n".join(rows)[:max_chars]

    # ── updates ─────────────────────────────────────────────────────────────
    def record_turn(self, role: str, content: str) -> None:
        """Append a chat turn (``role`` = user/assistant)."""
        self.turns.append(
            {"role": role, "content": content, "ts": time.time()}
        )

    def set_data(
        self, *, edi_path: str | None = None, line: str | None = None
    ) -> None:
        """Update the active EDI path / survey line."""
        if edi_path is not None:
            self.edi_path = edi_path
        if line is not None:
            self.line = line

    def record_workflow(self, workflow: str, summary: str = "") -> None:
        """Remember the most recent workflow + its summary."""
        self.last_workflow = workflow
        self.last_summary = summary

    def note(self, key: str, value: Any) -> None:
        """Store an arbitrary fact for later turns."""
        self.facts[key] = value

    def recent_turns(self, n: int = 6) -> list[dict[str, Any]]:
        """Last *n* turns (for feeding back as conversation context)."""
        return self.turns[-n:] if n > 0 else []

    def reset(self) -> None:
        """Start a new conversation, dropping all prior context and identity."""
        fresh = type(self)()
        self.__dict__.update(fresh.__dict__)

    @staticmethod
    def bounded_turns(turns, *, max_chars=8000, max_turns=6):
        """Return a recent suffix of whole text turns and an omission count.

        Code/artifact metadata stays in the original history; this is only
        the text supplied to question answering. Never truncate a requirement.
        """
        selected, size = [], 0
        for turn in reversed(turns or []):
            content = str(turn.get("content") or "")
            if len(selected) >= max_turns or size + len(content) > max_chars:
                break
            selected.insert(0, {"role": turn.get("role", "user"), "content": content})
            size += len(content)
        return selected, len(turns or []) - len(selected)

    # ── serialisation ───────────────────────────────────────────────────────
    def to_dict(self) -> dict[str, Any]:
        return deepcopy({
            "session_id": self.session_id,
            "edi_path": self.edi_path,
            "line": self.line,
            "last_workflow": self.last_workflow,
            "last_summary": self.last_summary,
            "turns": self.turns,
            "facts": self.facts,
            "project_id": self.project_id,
        })

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> SessionState:
        known = {f for f in cls.__dataclass_fields__}  # noqa: PLC0208
        return cls(**deepcopy({k: v for k, v in (d or {}).items() if k in known}))
