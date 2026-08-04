# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt/utils/progress.py

Package-wide progress-reporting engine.

Mirrors the verbosity behaviour Keras' ``Progbar`` gives its training
loops, but as a single, reusable primitive for pyCSAMT: AI-inversion
agents, forward-dataset generation, and geological-realization
sampling all drive the same class instead of hand-rolling their own
``if verbose and i % step == 0: print(...)`` throttling.

Verbosity levels (see :func:`normalize_verbose`)
--------------------------------------------------
``0`` (``False``, ``"silent"``)
    No output at all.
``1`` (``True``, ``"bar"``)
    A single animated line: bar + ETA + rate + metrics, rendered by
    ``tqdm`` when available.
``2`` (``"log"``)
    One printed line per update (or per *log_every* updates) with no
    carriage-return redraw — friendly to log files, CI output, and
    non-interactive terminals.
``"auto"``
    ``1`` if the stream is a tty, else ``2``.

This module has no dependency on :mod:`pycsamt.api`; package-wide
defaults (style, throttle interval, forced verbosity) are layered on
top by :mod:`pycsamt.api.view.progress`.
"""

from __future__ import annotations

import sys
import time
from collections.abc import Iterable, Iterator
from typing import Any, TypeVar

__all__ = [
    "DEFAULT_VERBOSE",
    "ProgressBar",
    "normalize_verbose",
    "progress",
]

T = TypeVar("T")

DEFAULT_VERBOSE = 1

_SILENT = {"silent", "off", "none", "no", "false"}
_BAR = {"bar", "on", "yes", "true"}
_LOG = {"log", "line", "lines"}


def _stream_is_tty(stream: Any) -> bool:
    return bool(getattr(stream, "isatty", lambda: False)())


def normalize_verbose(
    verbose: bool | int | str | None,
    *,
    default: int = DEFAULT_VERBOSE,
    stream: Any = None,
) -> int:
    """Normalize a verbosity argument to ``0`` / ``1`` / ``2``.

    Accepts the historical ``bool`` convention used across pyCSAMT
    (``True`` -> ``1``, ``False`` -> ``0``) as well as Keras-style
    integer levels and a handful of string aliases, so existing
    ``verbose: bool = True`` call sites keep working unchanged.

    Parameters
    ----------
    verbose : bool, int, str, or None
        ``None`` resolves to *default*. ``"auto"`` resolves to ``1``
        if *stream* (default ``sys.stderr``) is a tty, else ``2``.
    default : int, default 1
        Value used when *verbose* is ``None``.
    stream : file-like, optional
        Stream used to resolve ``"auto"``. Defaults to ``sys.stderr``.

    Returns
    -------
    int
        ``0`` (silent), ``1`` (bar), or ``2`` (one line per update).
    """
    if verbose is None:
        return normalize_verbose(default, default=default, stream=stream)
    if isinstance(verbose, bool):
        return 1 if verbose else 0
    if isinstance(verbose, str):
        mode = verbose.strip().lower()
        if mode == "auto":
            return 1 if _stream_is_tty(stream or sys.stderr) else 2
        if mode in _SILENT:
            return 0
        if mode in _BAR:
            return 1
        if mode in _LOG:
            return 2
        raise ValueError(
            "verbose must be a bool, int, None, or one of "
            f"{sorted(_SILENT | _BAR | _LOG | {'auto'})!r}; got {verbose!r}."
        )
    level = int(verbose)
    if level <= 0:
        return 0
    if level == 1:
        return 1
    return 2


def _format_metrics(metrics: dict[str, Any] | None) -> str:
    if not metrics:
        return ""
    parts = []
    for key, value in metrics.items():
        if isinstance(value, float):
            parts.append(f"{key}: {value:.4f}")
        else:
            parts.append(f"{key}: {value}")
    return " - ".join(parts)


_BAR_FORMAT = (
    "{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} "
    "[{elapsed}<{remaining}, {rate_fmt}{postfix}]"
)


class ProgressBar:
    """A verbose-aware progress reporter shared across pyCSAMT.

    Supports two usage styles:

    Manual updates (epoch loops, multiprocessing-driven loops)::

        with ProgressBar(total=epochs, desc="Training", unit="epoch",
                          verbose=verbose) as bar:
            for epoch in range(epochs):
                ...
                bar.update(1, metrics={"loss": loss, "val_loss": val})

    Iterable wrapping — see :func:`progress` for the shorthand::

        for item in ProgressBar.wrap(items, desc="Generating dataset"):
            ...

    Parameters
    ----------
    total : int, optional
        Number of expected updates. ``None`` disables the ETA/percentage
        and, in log mode, disables throttling (every update is logged).
    desc : str, optional
        Short label shown before the bar / at the start of each line.
    unit : str, default "it"
        Unit name shown next to counts (``"epoch"``, ``"sample"``, ...).
    verbose : bool, int, str, or None, default None
        See :func:`normalize_verbose`.
    leave : bool, default False
        Keep the final bar line on screen once closed (bar mode only).
    min_interval : float, default 0.1
        Minimum seconds between re-renders in bar mode (throttles
        redraw rate; forwarded to ``tqdm`` as ``mininterval``).
    log_every : int, optional
        In log mode, print one line every *log_every* updates plus the
        final one. Defaults to ``max(1, total // 10)`` (~10 lines
        total), or ``1`` when *total* is unknown.
    style : {"unicode", "ascii"}, default "ascii"
        Bar glyph set (bar mode only). ``"ascii"`` is the safer default
        on legacy Windows console code pages.
    """

    def __init__(
        self,
        total: int | None = None,
        *,
        desc: str | None = None,
        unit: str = "it",
        verbose: bool | int | str | None = None,
        leave: bool = False,
        min_interval: float = 0.1,
        log_every: int | None = None,
        style: str = "ascii",
    ) -> None:
        self.total = total
        self.desc = desc
        self.unit = unit
        self.verbose = normalize_verbose(verbose)
        self.leave = leave
        self._n = 0
        self._closed = False
        self._t_start = time.perf_counter()
        self._tqdm = None

        if self.verbose == 1:
            try:
                from tqdm.auto import tqdm as _tqdm_cls
            except Exception:
                # tqdm should always be present (hard dependency), but
                # degrade gracefully rather than crash a training loop.
                self.verbose = 2
            else:
                self._tqdm = _tqdm_cls(
                    total=total,
                    desc=desc,
                    unit=unit,
                    leave=leave,
                    mininterval=min_interval,
                    ascii=(style == "ascii"),
                    bar_format=_BAR_FORMAT if total else None,
                )

        if log_every is not None:
            self._log_every = max(1, int(log_every))
        else:
            self._log_every = max(1, total // 10) if total else 1

    @classmethod
    def wrap(
        cls,
        iterable: Iterable[T],
        *,
        total: int | None = None,
        desc: str | None = None,
        unit: str = "it",
        verbose: bool | int | str | None = None,
        leave: bool = False,
        min_interval: float = 0.1,
        log_every: int | None = None,
        style: str = "ascii",
    ) -> Iterator[T]:
        """Yield items from *iterable*, reporting progress as they pass."""
        if total is None:
            try:
                total = len(iterable)  # type: ignore[arg-type]
            except TypeError:
                total = None
        with cls(
            total=total,
            desc=desc,
            unit=unit,
            verbose=verbose,
            leave=leave,
            min_interval=min_interval,
            log_every=log_every,
            style=style,
        ) as bar:
            for item in iterable:
                yield item
                bar.update(1)

    def update(self, n: int = 1, *, metrics: dict[str, Any] | None = None) -> None:
        """Advance the counter by *n* and refresh the display."""
        if self._closed:
            return
        self._n += n
        if self.verbose == 0:
            return
        if self._tqdm is not None:
            if metrics:
                self._tqdm.set_postfix(**metrics, refresh=False)
            self._tqdm.update(n)
            return
        # Log mode.
        at_end = self.total is not None and self._n >= self.total
        if self._n % self._log_every == 0 or at_end:
            elapsed = time.perf_counter() - self._t_start
            label = f"{self.desc}: " if self.desc else ""
            count = (
                f"{self._n}/{self.total}" if self.total else f"{self._n}"
            )
            msg = f"{label}{count} {self.unit} [{elapsed:.1f}s]"
            extra = _format_metrics(metrics)
            if extra:
                msg = f"{msg} - {extra}"
            print(msg)

    def set_description(self, desc: str) -> None:
        """Update the label shown before the bar / log lines."""
        self.desc = desc
        if self._tqdm is not None:
            self._tqdm.set_description(desc)

    def close(self) -> None:
        """Finalize the bar. Safe to call multiple times."""
        if self._closed:
            return
        self._closed = True
        if self._tqdm is not None:
            self._tqdm.close()

    def __enter__(self) -> ProgressBar:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()

    def __repr__(self) -> str:  # noqa: D105
        return (
            f"ProgressBar(n={self._n}, total={self.total}, "
            f"verbose={self.verbose})"
        )


def progress(
    iterable: Iterable[T],
    *,
    total: int | None = None,
    desc: str | None = None,
    unit: str = "it",
    verbose: bool | int | str | None = None,
    leave: bool = False,
    min_interval: float = 0.1,
    log_every: int | None = None,
    style: str = "ascii",
) -> Iterator[T]:
    """Shorthand for :meth:`ProgressBar.wrap`."""
    yield from ProgressBar.wrap(
        iterable,
        total=total,
        desc=desc,
        unit=unit,
        verbose=verbose,
        leave=leave,
        min_interval=min_interval,
        log_every=log_every,
        style=style,
    )
