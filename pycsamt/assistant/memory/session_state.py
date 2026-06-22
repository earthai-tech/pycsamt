# -*- coding: utf-8 -*-
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
from dataclasses import dataclass, field
from typing import Any

__all__ = ["SessionState"]


@dataclass
class SessionState:
    """Mutable state for one chat session."""

    session_id: str = field(
        default_factory=lambda: uuid.uuid4().hex[:12]
    )
    edi_path: str | None = None
    line: str | None = None
    last_workflow: str | None = None
    last_summary: str | None = None
    turns: list[dict[str, Any]] = field(default_factory=list)
    facts: dict[str, Any] = field(default_factory=dict)

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
        return self.turns[-n:]

    # ── serialisation ───────────────────────────────────────────────────────
    def to_dict(self) -> dict[str, Any]:
        return {
            "session_id": self.session_id,
            "edi_path": self.edi_path,
            "line": self.line,
            "last_workflow": self.last_workflow,
            "last_summary": self.last_summary,
            "turns": self.turns,
            "facts": self.facts,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "SessionState":
        known = {f for f in cls.__dataclass_fields__}  # noqa: PLC0208
        return cls(**{k: v for k, v in (d or {}).items() if k in known})
