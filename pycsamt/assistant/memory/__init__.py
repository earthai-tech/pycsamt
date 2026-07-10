# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.memory
========================

State the assistant carries beyond a single retrieval call:

* :class:`~pycsamt.assistant.memory.session_state.SessionState` —
  in-memory state for one chat session (active line/EDI, last workflow,
  turns, facts);
* :class:`~pycsamt.assistant.memory.project_state.ProjectState` —
  disk-persisted per-project memory (lines worked on, preferences);
* :class:`~pycsamt.assistant.memory.workflow_history.WorkflowHistory` —
  append-only JSONL trace of executed workflows (observability).
"""

from __future__ import annotations

from .project_state import ProjectState
from .session_state import SessionState
from .workflow_history import WorkflowHistory, WorkflowRun

__all__ = [
    "SessionState",
    "ProjectState",
    "WorkflowHistory",
    "WorkflowRun",
]
