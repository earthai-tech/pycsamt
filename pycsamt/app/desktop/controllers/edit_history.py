# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Survey edit history for the desktop Edit menu (Qt-free).

Every change to the loaded survey -- a tool's result, a correction, the
pipeline, a manual edit -- becomes one :class:`EditStep` (label, the survey
before, the survey after).  :class:`EditHistory` gives undo / redo, a jump
to any step of the History list, "revert to as loaded" and a dirty flag
(edits not yet saved to disk).

Surveys are snapshotted with :func:`snapshot` (a deep copy), because some
tools edit sites in place; a snapshot that cannot be copied is kept by
reference.  The history is bounded (``max_steps``) so a long session does
not grow without limit.
"""

from __future__ import annotations

import copy
import datetime
from dataclasses import dataclass, field
from typing import Any

__all__ = ["EditHistory", "EditStep", "snapshot"]


def snapshot(sites: Any) -> Any:
    """An independent copy of *sites* (or *sites* itself if uncopyable)."""
    if sites is None:
        return None
    try:
        return copy.deepcopy(sites)
    except Exception:
        return sites


@dataclass
class EditStep:
    label: str
    before: Any
    after: Any
    time: str = field(
        default_factory=lambda: datetime.datetime.now().strftime("%H:%M:%S"))


class EditHistory:
    """Undo / redo stack over whole-survey snapshots."""

    def __init__(self, max_steps: int = 30) -> None:
        self.max_steps = int(max_steps)
        self._steps: list[EditStep] = []
        self._pos = 0  # number of applied steps
        self._original: Any = None
        self._saved_pos = 0  # position of the last save / load

    # ── lifecycle ─────────────────────────────────────────────────────
    def reset(self, loaded: Any) -> None:
        """A fresh survey was loaded: forget every step."""
        self._steps.clear()
        self._pos = 0
        self._original = snapshot(loaded)
        self._saved_pos = 0

    def record(self, label: str, before: Any, after: Any) -> None:
        """A change happened: *before* -> *after* (drops the redo branch)."""
        if self._original is None and before is not None:
            self._original = snapshot(before)
        del self._steps[self._pos:]
        if self._saved_pos > self._pos:
            self._saved_pos = -1  # the saved state is no longer reachable
        self._steps.append(EditStep(label or "Edit", snapshot(before),
                                    snapshot(after)))
        self._pos += 1
        overflow = len(self._steps) - self.max_steps
        if overflow > 0:
            del self._steps[:overflow]
            self._pos -= overflow
            self._saved_pos -= overflow

    def mark_saved(self) -> None:
        self._saved_pos = self._pos

    # ── navigation ────────────────────────────────────────────────────
    @property
    def can_undo(self) -> bool:
        return self._pos > 0

    @property
    def can_redo(self) -> bool:
        return self._pos < len(self._steps)

    def undo(self) -> Any:
        """The survey before the last applied step (moves back one)."""
        if not self.can_undo:
            raise IndexError("nothing to undo")
        self._pos -= 1
        return snapshot(self._steps[self._pos].before)

    def redo(self) -> Any:
        """The survey after the next step (moves forward one)."""
        if not self.can_redo:
            raise IndexError("nothing to redo")
        step = self._steps[self._pos]
        self._pos += 1
        return snapshot(step.after)

    def jump(self, n_applied: int) -> Any:
        """The survey with the first *n_applied* steps applied."""
        n = max(0, min(int(n_applied), len(self._steps)))
        self._pos = n
        if n == 0:
            first = self._steps[0].before if self._steps else self._original
            return snapshot(first)
        return snapshot(self._steps[n - 1].after)

    def original(self) -> Any:
        """The survey as it was loaded (not a history step itself)."""
        return snapshot(self._original)

    # ── state ─────────────────────────────────────────────────────────
    @property
    def position(self) -> int:
        return self._pos

    @property
    def dirty(self) -> bool:
        """Edits since the last load / save."""
        return self._pos != self._saved_pos

    @property
    def undo_label(self) -> str:
        return self._steps[self._pos - 1].label if self.can_undo else ""

    @property
    def redo_label(self) -> str:
        return self._steps[self._pos].label if self.can_redo else ""

    def entries(self) -> list[tuple[str, str, bool]]:
        """``(time, label, applied)`` for every step, oldest first."""
        return [(s.time, s.label, i < self._pos)
                for i, s in enumerate(self._steps)]

    def __len__(self) -> int:
        return len(self._steps)
