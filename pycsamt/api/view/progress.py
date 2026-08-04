"""Progress helpers used by public pyCSAMT APIs.

:mod:`pycsamt.utils.progress` implements the dependency-light progress
engine (:class:`~pycsamt.utils.progress.ProgressBar`, verbosity
normalization). This module layers package-wide *policy* on top of it
— a configurable singleton (:data:`PYCSAMT_PROGRESS`), following the
same ``configure()``/``context()``/``reset()`` shape as
:data:`pycsamt.api.pipe.PYCSAMT_PIPE` and
:data:`pycsamt.api.view.PYCSAMT_API_VIEW` — plus thin
backward-compatible wrappers (:func:`iter_progress`,
:func:`progress_enabled`) used elsewhere in the API layer.

Examples
--------
Silence every progress bar package-wide (e.g. for CI or batch jobs)::

    from pycsamt.api.view import configure_progress
    configure_progress(force_verbose=0)

Prefer plain log lines (no carriage-return redraw) in a non-tty run::

    from pycsamt.api.view import configure_progress
    configure_progress(default_verbose="log")
"""

from __future__ import annotations

import copy
import sys
from collections.abc import Generator, Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, TypeVar

from pycsamt.utils.progress import (
    DEFAULT_VERBOSE,
    ProgressBar,
    normalize_verbose,
)

T = TypeVar("T")

__all__ = [
    "PYCSAMT_PROGRESS",
    "ProgressAPIConfig",
    "ProgressConfig",
    "configure_progress",
    "get_progress_bar",
    "iter_progress",
    "progress_enabled",
    "reset_progress",
]

_STYLES = {"ascii", "unicode"}


@dataclass
class ProgressConfig:
    """Configuration for terminal progress display."""

    enabled: bool | str = "auto"
    desc: str | None = None
    leave: bool = False
    unit: str = "it"
    total: int | None = None


class ProgressAPIConfig:
    """Package-wide policy for :class:`~pycsamt.utils.progress.ProgressBar`.

    Attributes
    ----------
    default_verbose:
        Verbosity used when a call site's own ``verbose`` argument is
        ``None``. Does *not* override an explicit ``verbose=True/False``
        passed by the caller — use :attr:`force_verbose` for that.
    force_verbose:
        When set (``0``, ``1``, or ``2``), overrides every call site's
        ``verbose`` argument, explicit or not. ``None`` (default)
        leaves callers in control.
    style:
        Bar glyph set — ``"ascii"`` (default, safe on legacy Windows
        console code pages) or ``"unicode"``.
    min_interval:
        Minimum seconds between bar redraws.
    log_every:
        Default line-throttle for log-mode (verbose level ``2``)
        progress. ``None`` auto-computes ``~10`` lines per run.
    """

    def __init__(self) -> None:
        self._init_defaults()

    def _init_defaults(self) -> None:
        self.default_verbose: bool | int | str = DEFAULT_VERBOSE
        self.force_verbose: int | None = None
        self.style: str = "ascii"
        self.min_interval: float = 0.1
        self.log_every: int | None = None

    def configure(self, **kw: Any) -> None:
        """Set one or more configuration attributes by keyword."""
        if "style" in kw and kw["style"] not in _STYLES:
            raise ValueError(f"style must be one of {sorted(_STYLES)}.")
        for key, value in kw.items():
            if not hasattr(self, key):
                valid = sorted(vars(self))
                raise AttributeError(
                    f"Unknown progress config key {key!r}. Valid keys: {valid}"
                )
            setattr(self, key, value)

    @contextmanager
    def context(self, **kw: Any) -> Generator[ProgressAPIConfig, None, None]:
        """Temporarily override progress config, then restore it."""
        snapshot = copy.copy(vars(self))
        try:
            if kw:
                self.configure(**kw)
            yield self
        finally:
            self.__dict__.update(snapshot)

    def reset(self) -> None:
        """Reset to package defaults."""
        self._init_defaults()

    def resolve(self, verbose: bool | int | str | None = None) -> int:
        """Resolve a call site's ``verbose`` argument to ``0``/``1``/``2``."""
        if self.force_verbose is not None:
            return normalize_verbose(self.force_verbose)
        return normalize_verbose(verbose, default=normalize_verbose(
            self.default_verbose
        ))

    def summary(self) -> str:
        """Return a human-readable configuration summary."""
        lines = ["ProgressAPIConfig"]
        for key, value in vars(self).items():
            lines.append(f"  {key}: {value!r}")
        return "\n".join(lines)

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()


PYCSAMT_PROGRESS = ProgressAPIConfig()


def configure_progress(**kw: Any) -> None:
    """Configure :data:`PYCSAMT_PROGRESS` with keyword arguments."""
    PYCSAMT_PROGRESS.configure(**kw)


def reset_progress() -> None:
    """Reset :data:`PYCSAMT_PROGRESS` to package defaults."""
    PYCSAMT_PROGRESS.reset()


def get_progress_bar(
    total: int | None = None,
    *,
    desc: str | None = None,
    unit: str = "it",
    verbose: bool | int | str | None = None,
    leave: bool = False,
    log_every: int | None = None,
) -> ProgressBar:
    """Build a :class:`~pycsamt.utils.progress.ProgressBar` honouring
    :data:`PYCSAMT_PROGRESS` (style, throttle, forced verbosity).

    This is the recommended entry point for AI-inversion agents,
    dataset generation, and realization loops — it lets a user silence
    or reshape progress output package-wide via :func:`configure_progress`
    without touching every function's ``verbose`` default.
    """
    cfg = PYCSAMT_PROGRESS
    return ProgressBar(
        total=total,
        desc=desc,
        unit=unit,
        verbose=cfg.resolve(verbose),
        leave=leave,
        min_interval=cfg.min_interval,
        log_every=log_every if log_every is not None else cfg.log_every,
        style=cfg.style,
    )


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
    """Yield items, optionally wrapped by a :class:`ProgressBar`."""
    if not progress_enabled(enabled):
        yield from iterable
        return
    yield from ProgressBar.wrap(
        iterable,
        total=total,
        desc=desc,
        unit=unit,
        verbose=1,
        leave=leave,
        style=PYCSAMT_PROGRESS.style,
        min_interval=PYCSAMT_PROGRESS.min_interval,
    )
