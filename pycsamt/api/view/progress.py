"""Progress helpers used by public pyCSAMT APIs."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import Iterable, Iterator, TypeVar

T = TypeVar("T")

__all__ = ["ProgressConfig", "iter_progress", "progress_enabled"]


@dataclass
class ProgressConfig:
    """Configuration for terminal progress display."""

    enabled: bool | str = "auto"
    desc: str | None = None
    leave: bool = False
    unit: str = "it"
    total: int | None = None


def progress_enabled(value: bool | str = "auto") -> bool:
    """Return whether progress should be displayed."""
    if isinstance(value, str):
        mode = value.lower()
        if mode == "auto":
            return bool(getattr(sys.stderr, "isatty", lambda: False)())
        if mode in {"1", "true", "yes", "on"}:
            return True
        if mode in {"0", "false", "no", "off", "none"}:
            return False
        raise ValueError("progress must be bool or 'auto'")
    return bool(value)


def iter_progress(
    iterable: Iterable[T],
    *,
    enabled: bool | str = "auto",
    desc: str | None = None,
    leave: bool = False,
    unit: str = "it",
    total: int | None = None,
) -> Iterator[T]:
    """Yield items, optionally wrapped by ``tqdm``."""
    if not progress_enabled(enabled):
        yield from iterable
        return
    try:
        from tqdm import tqdm
    except Exception:
        yield from iterable
        return
    yield from tqdm(
        iterable,
        desc=desc,
        leave=leave,
        unit=unit,
        total=total,
    )
