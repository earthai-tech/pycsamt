# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Foundation for all native Occam1D objects.

This module defines the strict contract shared by Occam1D file containers,
builders, runners, result objects, and plotting helpers. It centralizes:

* electromagnetic utilities inherited from :class:`MTBase`;
* validated verbosity and diagnostic logging;
* path normalization and file-boundary checks;
* explicit object lifecycle state;
* diagnostics, metadata, and warning provenance.

Logging and verbosity intentionally remain independent. Logging follows the
configured logger. Verbosity controls only concise user-facing messages.
"""

from __future__ import annotations

import logging
import sys
from collections.abc import Mapping
from enum import Enum
from os import PathLike as OSPathLike
from pathlib import Path
from types import MappingProxyType
from typing import Any, TextIO, Union

from ...core.base import MTBase
from ...log.logger import get_logger

PathLike = Union[str, OSPathLike]

__all__ = ["Occam1DBase", "Occam1DObjectState", "PathLike"]


class Occam1DObjectState(str, Enum):
    """Lifecycle states exposed by :class:`Occam1DBase`.

    ``NEW`` identifies an in-memory object. ``BOUND`` identifies an object
    associated with a filesystem path. ``INVALID`` records an explicit
    scientific or operational failure.
    """

    NEW = "new"
    BOUND = "bound"
    INVALID = "invalid"


class Occam1DBase(MTBase):
    r"""Provide the common contract for Occam1D API objects.

    This base is format-neutral. A data container may bind its ``path`` to an
    ``EMData_1.x`` file, while a runner may never bind a single file. Readers
    use :meth:`_require_input_file`; writers use
    :meth:`_prepare_output_file` and bind only after successful output.

    Parameters
    ----------
    verbose : int or bool, default=0
        User-facing verbosity. ``0`` is quiet, ``1`` reports major stages,
        and larger integers may report progressively finer detail. It does not
        configure diagnostic logging.
    logger : logging.Logger or logger-like object, optional
        Diagnostic logger providing callable ``debug``, ``info``, ``warning``,
        ``error``, and ``exception`` methods. A class-specific pyCSAMT logger
        is created when omitted.
    path : path-like, optional
        File associated with the object. It is normalized to an absolute
        :class:`pathlib.Path` without requiring existence.
    metadata : mapping, optional
        Free-form provenance copied into the object.
    stream : writable text stream, optional
        Destination used by :meth:`_emit`. Defaults to :data:`sys.stdout`.

    Attributes
    ----------
    verbose : int
        Validated user-facing verbosity level.
    logger : logging.Logger
        Diagnostic logger independent of ``verbose``.
    path : pathlib.Path or None
        Normalized associated path.
    metadata : dict
        Mutable, free-form provenance metadata.

    Raises
    ------
    TypeError
        If an argument has an unsupported type or unexpected constructor
        keywords are supplied.
    ValueError
        If verbosity is negative or a path is empty.

    Notes
    -----
    Inheritance from :class:`~pycsamt.core.base.MTBase` is deliberate.
    Occam1D objects operate on impedance, apparent resistivity, phase, and
    determinant soundings and should reuse canonical pyCSAMT conversions.

    Lifecycle state describes provenance, not inversion convergence. A parsed
    iteration file is ``BOUND`` even if its inversion did not converge.

    Examples
    --------
    Construct an in-memory object:

    >>> base = Occam1DBase(metadata={"station": "S01"})
    >>> base.state.value
    'new'
    >>> base.is_bound
    False

    Bind a path after successful I/O:

    >>> path = base._bind_path("run/S01/Startup")
    >>> path.name
    'Startup'
    >>> base.state.value
    'bound'

    User messages remain separate from logging:

    >>> import io
    >>> stream = io.StringIO()
    >>> base = Occam1DBase(verbose=1, stream=stream)
    >>> base._emit("input files prepared")
    True
    >>> stream.getvalue().strip()
    'input files prepared'
    """

    __repr_exclude__ = {"logger", "verbose", "metadata"}

    def __init__(
        self,
        *,
        verbose: int | bool = 0,
        logger: logging.Logger | Any | None = None,
        path: PathLike | None = None,
        metadata: Mapping[str, Any] | None = None,
        stream: TextIO | None = None,
        **kwargs: Any,
    ) -> None:
        """Initialize validated shared Occam1D state."""
        if kwargs:
            names = ", ".join(sorted(str(name) for name in kwargs))
            raise TypeError(f"Unexpected Occam1D base arguments: {names}")

        self.verbose = self._validate_verbose(verbose)
        self.logger = self._validate_logger(logger)
        self._stream = self._validate_stream(stream)
        self.metadata = self._validate_metadata(metadata)
        self._path: Path | None = None
        self._state = Occam1DObjectState.NEW
        self._state_reason: str | None = None
        self._warnings: list[str] = []
        if path is not None:
            self._bind_path(path)

    @staticmethod
    def _validate_verbose(value: int | bool) -> int:
        """Return a non-negative integer verbosity level."""
        if isinstance(value, bool):
            return int(value)
        if not isinstance(value, int):
            raise TypeError("verbose must be a non-negative integer or bool.")
        if value < 0:
            raise ValueError("verbose must be non-negative.")
        return value

    def _validate_logger(self, value):
        """Return a logger satisfying the diagnostic logging contract."""
        if value is None:
            name = f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            return get_logger(name)
        methods = ("debug", "info", "warning", "error", "exception")
        missing = [
            name
            for name in methods
            if not callable(getattr(value, name, None))
        ]
        if missing:
            names = ", ".join(missing)
            raise TypeError(f"logger is missing callable methods: {names}.")
        return value

    @staticmethod
    def _validate_stream(value: TextIO | None) -> TextIO:
        """Return a writable stream for user-facing messages."""
        stream = sys.stdout if value is None else value
        if not callable(getattr(stream, "write", None)):
            raise TypeError("stream must provide a callable write method.")
        return stream

    @staticmethod
    def _validate_metadata(value):
        """Return a defensive metadata copy."""
        if value is None:
            return {}
        if not isinstance(value, Mapping):
            raise TypeError("metadata must be a mapping or None.")
        return dict(value)

    @property
    def path(self) -> Path | None:
        """Normalized filesystem path associated with this object."""
        return self._path

    @path.setter
    def path(self, value: PathLike | None) -> None:
        if value is None:
            self._path = None
            if self._state is Occam1DObjectState.BOUND:
                self._state = Occam1DObjectState.NEW
                self._state_reason = None
            return
        self._bind_path(value)

    @property
    def state(self) -> Occam1DObjectState:
        """Current object lifecycle state."""
        return self._state

    @property
    def state_reason(self) -> str | None:
        """Reason recorded when the object was marked invalid."""
        return self._state_reason

    @property
    def is_bound(self) -> bool:
        """Whether the object is associated with a filesystem path."""
        return self._state is Occam1DObjectState.BOUND

    @property
    def is_valid(self) -> bool:
        """Whether the object has not been explicitly marked invalid."""
        return self._state is not Occam1DObjectState.INVALID

    @property
    def warnings(self) -> tuple[str, ...]:
        """Immutable snapshot of warnings recorded on the object."""
        return tuple(self._warnings)

    @property
    def metadata_view(self):
        """Read-only live view of object metadata."""
        return MappingProxyType(self.metadata)

    def _normalize_path(self, value: PathLike) -> Path:
        """Return an absolute path without checking existence."""
        if not isinstance(value, (str, OSPathLike)):
            raise TypeError("path must be a string or os.PathLike object.")
        if isinstance(value, str) and not value.strip():
            raise ValueError("path cannot be empty.")
        return Path(value).expanduser().resolve(strict=False)

    def _bind_path(self, value: PathLike) -> Path:
        """Associate the object with a path after successful I/O.

        This method does not require existence. Readers should first call
        :meth:`_require_input_file`; writers should bind after their stream
        closes successfully.
        """
        path = self._normalize_path(value)
        self._path = path
        if self._state is not Occam1DObjectState.INVALID:
            self._state = Occam1DObjectState.BOUND
            self._state_reason = None
        return path

    def _require_input_file(
        self,
        value: PathLike | None = None,
        *,
        purpose: str = "Occam1D input",
    ) -> Path:
        """Resolve and validate a readable regular input file.

        Parameters
        ----------
        value : path-like, optional
            Candidate file. The object's bound path is used when omitted.
        purpose : str, default="Occam1D input"
            Human-readable role included in exceptions.

        Returns
        -------
        pathlib.Path
            Existing readable regular file.

        Raises
        ------
        ValueError
            If no candidate path is available.
        FileNotFoundError
            If the candidate does not exist.
        IsADirectoryError
            If the candidate is a directory.
        PermissionError
            If the file cannot be opened for reading.
        """
        candidate = self.path if value is None else self._normalize_path(value)
        if candidate is None:
            raise ValueError(f"No path is available for {purpose}.")
        if not candidate.exists():
            raise FileNotFoundError(f"{purpose} file not found: {candidate}")
        if not candidate.is_file():
            raise IsADirectoryError(f"{purpose} is not a file: {candidate}")
        try:
            with candidate.open("r", encoding="utf8", errors="replace"):
                pass
        except PermissionError:
            raise PermissionError(f"{purpose} is not readable: {candidate}")
        return candidate

    def _prepare_output_file(
        self,
        value: PathLike,
        *,
        create_parent: bool = True,
    ) -> Path:
        """Validate an output target and optionally create its parent.

        Directories are rejected as file targets. Existing regular files are
        allowed; overwrite policy belongs to the public format writer.
        """
        target = self._normalize_path(value)
        if target.exists() and target.is_dir():
            raise IsADirectoryError(f"Output path is a directory: {target}")
        parent = target.parent
        if create_parent:
            parent.mkdir(parents=True, exist_ok=True)
        if not parent.is_dir():
            raise NotADirectoryError(
                f"Output parent is not a directory: {parent}"
            )
        return target

    def _emit(self, message: str, *, level: int = 1, end: str = "\n") -> bool:
        """Write a user-facing message when verbosity permits.

        Returns ``True`` when emitted and ``False`` when suppressed. This
        method never writes to the diagnostic logger.
        """
        if not isinstance(message, str):
            raise TypeError("message must be a string.")
        if not isinstance(level, int) or isinstance(level, bool) or level < 1:
            raise ValueError("level must be a positive integer.")
        if not isinstance(end, str):
            raise TypeError("end must be a string.")
        if self.verbose < level:
            return False
        self._stream.write(message + end)
        flush = getattr(self._stream, "flush", None)
        if callable(flush):
            flush()
        return True

    def add_warning(self, message: str, *, log: bool = True) -> None:
        """Record one deduplicated recoverable warning.

        ``log=True`` also calls ``logger.warning`` but never prints through the
        verbosity stream.
        """
        if not isinstance(message, str):
            raise TypeError("warning message must be a string.")
        text = message.strip()
        if not text:
            raise ValueError("warning message cannot be empty.")
        is_new = text not in self._warnings
        if is_new:
            self._warnings.append(text)
        if log and is_new:
            self.logger.warning("%s", text)

    def mark_invalid(self, reason: str) -> None:
        """Mark invalid state while preserving provenance diagnostics."""
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError("reason must be a non-empty string.")
        self._state = Occam1DObjectState.INVALID
        self._state_reason = reason.strip()
        self.logger.error("Occam1D object marked invalid: %s", reason.strip())

    def mark_valid(self) -> None:
        """Clear invalid state after successful explicit revalidation.

        The resulting state is ``BOUND`` when a path is associated with the
        object and ``NEW`` otherwise. Callers should invoke :meth:`validate`
        before this method when recovery depends on scientific invariants.
        """
        self._state = (
            Occam1DObjectState.BOUND
            if self.path is not None
            else Occam1DObjectState.NEW
        )
        self._state_reason = None

    def update_metadata(self, /, **values: Any):
        """Update metadata in place and return ``self``."""
        self.metadata.update(values)
        return self

    def diagnostics(self) -> dict[str, Any]:
        """Return JSON-friendly lifecycle and provenance diagnostics.

        Scientific subclasses may extend this mapping but should preserve its
        core keys.
        """
        return {
            "class": self.__class__.__qualname__,
            "module": self.__class__.__module__,
            "state": self.state.value,
            "state_reason": self.state_reason,
            "path": str(self.path) if self.path is not None else None,
            "is_bound": self.is_bound,
            "is_valid": self.is_valid,
            "verbose": self.verbose,
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }

    def validate(self) -> None:
        """Validate base invariants and invoke the subclass hook.

        Subclasses should override :meth:`_validate_state` instead of this
        method so shared invariants cannot be accidentally skipped.
        """
        self._validate_verbose(self.verbose)
        self._validate_logger(self.logger)
        self._validate_metadata(self.metadata)
        if self.path is not None:
            self._normalize_path(self.path)
        self._validate_state()

    def _validate_state(self) -> None:
        """Hook for scientific validation in concrete subclasses."""
        return None
