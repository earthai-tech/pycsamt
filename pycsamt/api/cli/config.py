# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.api.cli.config
======================

Package-wide CLI configuration.  Same ``configure / reset / context``
pattern used by the rest of ``pycsamt.api``.

Quick start
-----------
**Read the current defaults**::

    from pycsamt.api.cli import PYCSAMT_CLI
    print(PYCSAMT_CLI)

**Change one setting for the whole session**::

    from pycsamt.api.cli import configure_cli
    configure_cli(log__level=1, output__format="json")

**Temporarily override, then restore**::

    from pycsamt.api.cli import PYCSAMT_CLI
    with PYCSAMT_CLI.context(log__level=2):
        ...   # debug verbosity active here only

**Reset to package defaults**::

    from pycsamt.api.cli import reset_cli
    reset_cli()
"""

from __future__ import annotations

import copy
import os
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Allowed values
# ---------------------------------------------------------------------------
_OUTPUT_FORMATS = {"text", "json", "csv"}
_VERBOSE_LEVELS = {0: "quiet", 1: "info", 2: "debug"}


# ---------------------------------------------------------------------------
# Sub-config dataclasses
# ---------------------------------------------------------------------------

@dataclass
class LogConfig:
    """Logging / verbosity settings."""

    level: int = 0
    color: bool = True
    file: Path | None = None

    def __post_init__(self) -> None:
        if self.level not in _VERBOSE_LEVELS:
            raise ValueError(
                f"log.level must be one of {sorted(_VERBOSE_LEVELS)}."
            )


@dataclass
class OutputConfig:
    """Output format and destination settings."""

    format: str = "text"
    dir: Path = field(default_factory=lambda: Path("."))
    overwrite: bool = False

    def __post_init__(self) -> None:
        if self.format not in _OUTPUT_FORMATS:
            raise ValueError(
                f"output.format must be one of {sorted(_OUTPUT_FORMATS)}."
            )
        self.dir = Path(self.dir)


@dataclass
class BuildConfig:
    """Build / processing execution settings."""

    n_jobs: int = 1
    cache: bool = True
    cache_dir: Path | None = None

    def __post_init__(self) -> None:
        if self.n_jobs < 1:
            raise ValueError("build.n_jobs must be >= 1.")
        if self.cache_dir is not None:
            self.cache_dir = Path(self.cache_dir)


# ---------------------------------------------------------------------------
# Container
# ---------------------------------------------------------------------------

class PyCSAMTCLI:
    """Package-wide CLI configuration container.

    Access the singleton via :data:`PYCSAMT_CLI`.  Do not instantiate
    directly unless you need an independent copy.
    """

    def __init__(self) -> None:
        self._init_defaults()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def configure(self, **kw: Any) -> None:
        """Set attributes using ``section__attribute`` dotted paths.

        Examples
        --------
        >>> PYCSAMT_CLI.configure(log__level=1, output__format="json")
        >>> PYCSAMT_CLI.configure(build__n_jobs=4)
        """
        for path, value in kw.items():
            parts = path.split("__")
            obj = self
            for part in parts[:-1]:
                obj = getattr(obj, part)
            setattr(obj, parts[-1], value)

    @contextmanager
    def context(self, **kw: Any) -> Generator[PyCSAMTCLI, None, None]:
        """Temporarily override settings, then restore on exit.

        Examples
        --------
        >>> with PYCSAMT_CLI.context(log__level=2):
        ...     run_command()
        """
        snapshot = self._snapshot()
        try:
            if kw:
                self.configure(**kw)
            yield self
        finally:
            self._restore(snapshot)

    def reset(self) -> None:
        """Reset all settings to package defaults."""
        self._init_defaults()

    def load_env(self) -> None:
        """Override settings from environment variables.

        Recognised variables
        --------------------
        PYCSAMT_VERBOSE   integer 0-2
        PYCSAMT_NO_COLOR  any non-empty value disables colour
        PYCSAMT_OUTPUT    text | json | csv
        PYCSAMT_OUTPUT_DIR  path
        PYCSAMT_JOBS      integer >= 1
        """
        if (v := os.environ.get("PYCSAMT_VERBOSE")) is not None:
            self.log.level = int(v)
        if os.environ.get("PYCSAMT_NO_COLOR"):
            self.log.color = False
        if (v := os.environ.get("PYCSAMT_OUTPUT")) is not None:
            self.output.format = v
        if (v := os.environ.get("PYCSAMT_OUTPUT_DIR")) is not None:
            self.output.dir = Path(v)
        if (v := os.environ.get("PYCSAMT_JOBS")) is not None:
            self.build.n_jobs = int(v)

    def summary(self) -> str:
        """Return a human-readable summary of current settings."""
        lines = [
            "PyCSAMTCLI",
            f"  log.level     = {self.log.level!r}  "
            f"({_VERBOSE_LEVELS.get(self.log.level, '?')})",
            f"  log.color     = {self.log.color}",
            f"  log.file      = {self.log.file!r}",
            f"  output.format = {self.output.format!r}",
            f"  output.dir    = {self.output.dir!r}",
            f"  output.overwrite = {self.output.overwrite}",
            f"  build.n_jobs  = {self.build.n_jobs}",
            f"  build.cache   = {self.build.cache}",
            f"  build.cache_dir = {self.build.cache_dir!r}",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _init_defaults(self) -> None:
        self.log = LogConfig()
        self.output = OutputConfig()
        self.build = BuildConfig()

    def _snapshot(self) -> dict[str, Any]:
        return {
            "log": copy.deepcopy(self.log),
            "output": copy.deepcopy(self.output),
            "build": copy.deepcopy(self.build),
        }

    def _restore(self, snapshot: dict[str, Any]) -> None:
        self.log = snapshot["log"]
        self.output = snapshot["output"]
        self.build = snapshot["build"]


# ---------------------------------------------------------------------------
# Singleton + module-level helpers
# ---------------------------------------------------------------------------

PYCSAMT_CLI = PyCSAMTCLI()
PYCSAMT_CLI.load_env()


def configure_cli(**kw: Any) -> None:
    """Configure :data:`PYCSAMT_CLI` with dotted-path arguments.

    Examples
    --------
    >>> from pycsamt.api.cli import configure_cli
    >>> configure_cli(log__level=1, output__format="json")
    """
    PYCSAMT_CLI.configure(**kw)


def reset_cli() -> None:
    """Reset :data:`PYCSAMT_CLI` to package defaults."""
    PYCSAMT_CLI.reset()


__all__ = [
    "BuildConfig",
    "LogConfig",
    "OutputConfig",
    "PYCSAMT_CLI",
    "PyCSAMTCLI",
    "configure_cli",
    "reset_cli",
]
