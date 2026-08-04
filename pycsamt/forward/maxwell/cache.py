# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Content-addressed cache for canonical Maxwell forward results.

Entries are keyed by ``MaxwellProblem.problem_hash``. Each result archive has
a SHA-256 sidecar, is written by atomic replacement, and is validated against
the requested problem when read. Lock files coordinate independent workers.
"""

from __future__ import annotations

import hashlib
import os
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .backends import MaxwellBackend
from .contracts import ForwardResult, MaxwellProblem
from .contracts_tri import TriProblem

__all__ = [
    "CacheCorruptionError",
    "CacheLockTimeoutError",
    "CacheEntry",
    "CacheStatistics",
    "MaxwellResultCache",
]


class CacheCorruptionError(RuntimeError):
    """Indicate that a cached archive failed integrity validation.

    Examples
    --------
    >>> isinstance(CacheCorruptionError("bad checksum"), RuntimeError)
    True
    """


class CacheLockTimeoutError(TimeoutError):
    """Indicate that a cache-key lock could not be acquired in time.

    Examples
    --------
    >>> isinstance(CacheLockTimeoutError("busy"), TimeoutError)
    True
    """


def _cache_key(value: str) -> str:
    key = str(value).strip().lower()
    hexadecimal = set("0123456789abcdef")
    if len(key) != 64 or any(char not in hexadecimal for char in key):
        raise ValueError(
            "cache key must be a 64-character hexadecimal SHA-256 digest."
        )
    return key


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class _KeyLock:
    def __init__(
        self,
        path: Path,
        timeout_s: float,
        poll_interval_s: float,
        stale_after_s: float,
    ) -> None:
        self.path = path
        self.timeout_s = timeout_s
        self.poll_interval_s = poll_interval_s
        self.stale_after_s = stale_after_s
        self._owned = False

    def __enter__(self) -> _KeyLock:
        deadline = time.monotonic() + self.timeout_s
        payload = f"pid={os.getpid()} created={time.time():.6f}\n".encode(
            "ascii"
        )
        while True:
            try:
                descriptor = os.open(
                    str(self.path),
                    os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                )
            except FileExistsError:
                self._remove_stale_lock()
                if time.monotonic() >= deadline:
                    raise CacheLockTimeoutError(
                        f"timed out waiting for cache lock {self.path.name!r}."
                    ) from None
                time.sleep(self.poll_interval_s)
                continue
            try:
                os.write(descriptor, payload)
            finally:
                os.close(descriptor)
            self._owned = True
            return self

    def __exit__(self, *args: Any) -> None:
        if self._owned:
            try:
                self.path.unlink()
            except FileNotFoundError:
                pass
            self._owned = False

    def _remove_stale_lock(self) -> None:
        try:
            age = time.time() - self.path.stat().st_mtime
        except FileNotFoundError:
            return
        if age <= self.stale_after_s:
            return
        try:
            self.path.unlink()
        except FileNotFoundError:
            pass


@dataclass(frozen=True)
class CacheEntry:
    """Describe one complete cache entry.

    Parameters
    ----------
    key : str
        Problem SHA-256 digest.
    archive_path : pathlib.Path
        Result archive location.
    size_bytes : int
        Combined archive and checksum size.
    modified_time_s : float
        Archive modification time as Unix seconds.

    Examples
    --------
    >>> entry = CacheEntry("0" * 64, Path("result.npz"), 10, 1.0)
    >>> entry.size_bytes
    10
    """

    key: str
    archive_path: Path
    size_bytes: int
    modified_time_s: float

    def __post_init__(self) -> None:
        key = _cache_key(self.key)
        path = Path(self.archive_path)
        if not isinstance(self.size_bytes, int) or self.size_bytes < 0:
            raise ValueError("size_bytes must be a non-negative integer.")
        modified = float(self.modified_time_s)
        if not modified >= 0 or not modified < float("inf"):
            raise ValueError(
                "modified_time_s must be finite and non-negative."
            )
        object.__setattr__(self, "key", key)
        object.__setattr__(self, "archive_path", path)
        object.__setattr__(self, "modified_time_s", modified)

    @property
    def checksum_path(self) -> Path:
        """Return the SHA-256 sidecar path.

        Returns
        -------
        pathlib.Path
            Archive path with ``.sha256`` appended.

        Examples
        --------
        >>> entry = CacheEntry("0" * 64, Path("r.npz"), 0, 0)
        >>> entry.checksum_path.name
        'r.npz.sha256'
        """
        suffix = self.archive_path.suffix + ".sha256"
        return self.archive_path.with_suffix(suffix)

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible entry metadata.

        Returns
        -------
        dict
            Key, path, size, and modification time.

        Examples
        --------
        >>> entry = CacheEntry("0" * 64, Path("r.npz"), 5, 2)
        >>> entry.to_dict()["size_bytes"]
        5
        """
        return {
            "key": self.key,
            "archive_path": str(self.archive_path),
            "size_bytes": self.size_bytes,
            "modified_time_s": self.modified_time_s,
        }


@dataclass(frozen=True)
class CacheStatistics:
    """Summarize the current on-disk cache state.

    Parameters
    ----------
    entry_count, total_bytes, orphan_count, corrupt_count : int
        Counts and storage for normal, incomplete, and quarantined files.

    Examples
    --------
    >>> CacheStatistics(2, 100, 0, 1).entry_count
    2
    """

    entry_count: int
    total_bytes: int
    orphan_count: int
    corrupt_count: int

    def __post_init__(self) -> None:
        names = (
            "entry_count",
            "total_bytes",
            "orphan_count",
            "corrupt_count",
        )
        for name in names:
            value = getattr(self, name)
            if not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a non-negative integer.")

    def to_dict(self) -> dict[str, int]:
        """Return JSON-compatible cache statistics.

        Returns
        -------
        dict
            Entry, byte, orphan, and corruption counts.

        Examples
        --------
        >>> CacheStatistics(1, 20, 0, 0).to_dict()["total_bytes"]
        20
        """
        return {
            "entry_count": self.entry_count,
            "total_bytes": self.total_bytes,
            "orphan_count": self.orphan_count,
            "corrupt_count": self.corrupt_count,
        }


class MaxwellResultCache:
    """Manage a validated, content-addressed result cache.

    Parameters
    ----------
    root : str or pathlib.Path
        Dedicated cache directory. It is created when absent.
    lock_timeout_s : float, default=300
        Maximum wait for another worker holding the same problem key.
    poll_interval_s : float, default=0.05
        Delay between lock acquisition attempts.
    stale_lock_s : float, default=3600
        Age after which an abandoned lock can be removed.
    quarantine_corrupt : bool, default=True
        Move corrupt files under ``root/corrupt``. When false, reads raise
        :class:`CacheCorruptionError` and leave the entry untouched.

    Examples
    --------
    >>> from tempfile import TemporaryDirectory
    >>> with TemporaryDirectory() as directory:
    ...     cache = MaxwellResultCache(directory)
    ...     cache.statistics().entry_count
    0
    """

    def __init__(
        self,
        root: str | Path,
        *,
        lock_timeout_s: float = 300.0,
        poll_interval_s: float = 0.05,
        stale_lock_s: float = 3600.0,
        quarantine_corrupt: bool = True,
    ) -> None:
        self._root = Path(root).expanduser().resolve()
        self._lock_timeout_s = self._positive_time(
            lock_timeout_s,
            "lock_timeout_s",
        )
        self._poll_interval_s = self._positive_time(
            poll_interval_s,
            "poll_interval_s",
        )
        self._stale_lock_s = self._positive_time(
            stale_lock_s,
            "stale_lock_s",
        )
        self._quarantine_corrupt = bool(quarantine_corrupt)
        self._entries = self._root / "entries"
        self._locks = self._root / "locks"
        self._corrupt = self._root / "corrupt"
        self._entries.mkdir(parents=True, exist_ok=True)
        self._locks.mkdir(parents=True, exist_ok=True)
        self._corrupt.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _positive_time(value: float, name: str) -> float:
        result = float(value)
        if not result > 0 or not result < float("inf"):
            raise ValueError(f"{name} must be positive and finite.")
        return result

    @property
    def root(self) -> Path:
        """Return the resolved cache root.

        Returns
        -------
        pathlib.Path
            Dedicated cache directory.

        Examples
        --------
        The returned path is always absolute.
        """
        return self._root

    def contains(self, problem: MaxwellProblem) -> bool:
        """Return whether a complete entry exists for a problem.

        Parameters
        ----------
        problem : MaxwellProblem
            Problem whose content hash identifies the entry.

        Returns
        -------
        bool
            True when both archive and checksum sidecar exist.

        Examples
        --------
        A newly created cache contains no problems.
        """
        self._require_problem(problem)
        archive, checksum = self._paths(problem.problem_hash)
        return archive.is_file() and checksum.is_file()

    def get(self, problem: MaxwellProblem) -> ForwardResult | None:
        """Load and validate a cached result.

        Parameters
        ----------
        problem : MaxwellProblem
            Exact problem expected by the caller.

        Returns
        -------
        ForwardResult or None
            Valid result, or None when no complete entry exists. Corruption is
            quarantined and treated as a miss when configured.

        Raises
        ------
        CacheCorruptionError
            If validation fails and quarantine is disabled.

        Examples
        --------
        Cache misses return None rather than raising ``KeyError``.
        """
        self._require_problem(problem)
        with self._lock(problem.problem_hash):
            return self._read_or_quarantine(problem)

    def _read_or_quarantine(
        self,
        problem: MaxwellProblem,
    ) -> ForwardResult | None:
        try:
            return self._read(problem)
        except CacheCorruptionError:
            if not self._quarantine_corrupt:
                raise
            self._quarantine(problem.problem_hash)
            return None

    def put(
        self,
        problem: MaxwellProblem,
        result: ForwardResult,
        *,
        overwrite: bool = False,
    ) -> CacheEntry:
        """Validate and atomically store one result.

        Parameters
        ----------
        problem : MaxwellProblem
            Exact simulation input.
        result : ForwardResult
            Canonical result matching ``problem``.
        overwrite : bool, default=False
            Replace an existing complete entry. Otherwise the validated
            existing entry is retained.

        Returns
        -------
        CacheEntry
            Metadata for the stored or retained entry.

        Examples
        --------
        Invalid problem/result pairs are rejected before writing an archive.
        """
        self._require_problem(problem)
        if not isinstance(result, ForwardResult):
            raise TypeError("result must be a ForwardResult.")
        result.validate_against(problem)
        with self._lock(problem.problem_hash):
            if not overwrite:
                try:
                    existing = self._read(problem)
                except CacheCorruptionError:
                    self._quarantine(problem.problem_hash)
                    existing = None
                if existing is not None:
                    return self.entry(problem.problem_hash)
            self._write(problem.problem_hash, result)
            return self.entry(problem.problem_hash)

    def get_or_solve(
        self,
        problem: MaxwellProblem,
        backend: MaxwellBackend,
    ) -> ForwardResult:
        """Return a hit or solve and cache one problem under a key lock.

        Parameters
        ----------
        problem : MaxwellProblem
            Simulation input and cache identity.
        backend : MaxwellBackend
            Conforming backend used only after a cache miss.

        Returns
        -------
        ForwardResult
            Valid cached or newly computed result.

        Notes
        -----
        The second read after locking prevents duplicate concurrent work.

        Examples
        --------
        Backend invocation is skipped whenever a validated hit exists.
        """
        self._require_problem(problem)
        if not isinstance(backend, MaxwellBackend):
            raise TypeError("backend must implement MaxwellBackend.")
        hit = self.get(problem)
        if hit is not None:
            return hit
        with self._lock(problem.problem_hash):
            hit = self._read_or_quarantine(problem)
            if hit is not None:
                return hit
            result = backend.solve(problem)
            if not isinstance(result, ForwardResult):
                raise TypeError("backend must return a ForwardResult.")
            result.validate_against(problem)
            self._write(problem.problem_hash, result)
            return result

    def entry(self, key: str) -> CacheEntry:
        """Return filesystem metadata for a complete entry.

        Parameters
        ----------
        key : str
            Problem SHA-256 digest.

        Returns
        -------
        CacheEntry
            Entry paths, size, and modification time.

        Raises
        ------
        KeyError
            If archive or checksum sidecar is missing.

        Examples
        --------
        This method does not deserialize the result archive.
        """
        normalized = _cache_key(key)
        archive, checksum = self._paths(normalized)
        if not archive.is_file() or not checksum.is_file():
            raise KeyError(
                f"cache entry {normalized!r} is incomplete or missing."
            )
        size = archive.stat().st_size + checksum.stat().st_size
        return CacheEntry(
            normalized,
            archive,
            size,
            archive.stat().st_mtime,
        )

    def entries(self) -> tuple[CacheEntry, ...]:
        """Return complete entries sorted by problem key.

        Returns
        -------
        tuple of CacheEntry
            Stable snapshot of complete cache entries.

        Examples
        --------
        An empty cache returns an empty tuple.
        """
        found = []
        for archive in self._entries.glob("*/*.npz"):
            try:
                found.append(self.entry(archive.stem))
            except (KeyError, ValueError):
                continue
        return tuple(sorted(found, key=lambda value: value.key))

    def remove(self, problem: MaxwellProblem) -> bool:
        """Remove one problem entry under its key lock.

        Parameters
        ----------
        problem : MaxwellProblem
            Exact problem identifying the entry.

        Returns
        -------
        bool
            True when at least one entry file was removed.

        Examples
        --------
        Removing an absent problem is a no-op returning False.
        """
        self._require_problem(problem)
        removed = False
        with self._lock(problem.problem_hash):
            for path in self._paths(problem.problem_hash):
                try:
                    path.unlink()
                    removed = True
                except FileNotFoundError:
                    pass
        return removed

    def prune(self, maximum_bytes: int) -> tuple[CacheEntry, ...]:
        """Remove oldest entries until storage is within a byte budget.

        Parameters
        ----------
        maximum_bytes : int
            Non-negative archive and checksum budget.

        Returns
        -------
        tuple of CacheEntry
            Entries removed, oldest first.

        Examples
        --------
        ``prune(0)`` removes every complete entry but leaves infrastructure.
        """
        if not isinstance(maximum_bytes, int) or maximum_bytes < 0:
            raise ValueError("maximum_bytes must be a non-negative integer.")
        entries = sorted(
            self.entries(),
            key=lambda value: (value.modified_time_s, value.key),
        )
        total = sum(value.size_bytes for value in entries)
        removed = []
        for value in entries:
            if total <= maximum_bytes:
                break
            archive, checksum = self._paths(value.key)
            with self._lock(value.key):
                for path in (archive, checksum):
                    try:
                        path.unlink()
                    except FileNotFoundError:
                        pass
            removed.append(value)
            total -= value.size_bytes
        return tuple(removed)

    def statistics(self) -> CacheStatistics:
        """Inspect complete, orphaned, and quarantined cache files.

        Returns
        -------
        CacheStatistics
            Current cache counts and complete-entry storage.

        Examples
        --------
        Statistics inspect metadata without deserializing archives.
        """
        entries = self.entries()
        archives = set(self._entries.glob("*/*.npz"))
        checksums = set(self._entries.glob("*/*.npz.sha256"))
        complete_archives = {value.archive_path for value in entries}
        complete_checksums = {value.checksum_path for value in entries}
        orphan_count = len(archives - complete_archives)
        orphan_count += len(checksums - complete_checksums)
        corrupt_count = sum(1 for path in self._corrupt.iterdir())
        return CacheStatistics(
            len(entries),
            sum(value.size_bytes for value in entries),
            orphan_count,
            corrupt_count,
        )

    @staticmethod
    def _require_problem(problem: MaxwellProblem) -> None:
        if not isinstance(problem, (MaxwellProblem, TriProblem)):
            raise TypeError("problem must be a MaxwellProblem or TriProblem.")

    def _paths(self, key: str) -> tuple[Path, Path]:
        normalized = _cache_key(key)
        directory = self._entries / normalized[:2]
        archive = directory / f"{normalized}.npz"
        checksum = directory / f"{normalized}.npz.sha256"
        return archive, checksum

    def _lock(self, key: str) -> _KeyLock:
        normalized = _cache_key(key)
        return _KeyLock(
            self._locks / f"{normalized}.lock",
            self._lock_timeout_s,
            self._poll_interval_s,
            self._stale_lock_s,
        )

    def _read(self, problem: MaxwellProblem) -> ForwardResult | None:
        archive, checksum = self._paths(problem.problem_hash)
        if not archive.exists() and not checksum.exists():
            return None
        if not archive.is_file() or not checksum.is_file():
            raise CacheCorruptionError("cache entry is incomplete.")
        try:
            expected = checksum.read_text(encoding="ascii").strip().lower()
        except OSError as exc:
            raise CacheCorruptionError(
                "cannot read checksum sidecar."
            ) from exc
        try:
            expected = _cache_key(expected)
        except ValueError as exc:
            raise CacheCorruptionError("checksum sidecar is invalid.") from exc
        if _file_digest(archive) != expected:
            raise CacheCorruptionError("result archive checksum mismatch.")
        try:
            result = ForwardResult.from_npz(archive)
            result.validate_against(problem)
        except Exception as exc:
            raise CacheCorruptionError(
                "result archive failed contract validation."
            ) from exc
        return result

    def _write(self, key: str, result: ForwardResult) -> None:
        archive, checksum = self._paths(key)
        archive.parent.mkdir(parents=True, exist_ok=True)
        token = uuid.uuid4().hex
        temporary_archive = archive.parent / (f".{archive.name}.{token}.npz")
        temporary_checksum = archive.parent / (f".{checksum.name}.{token}.tmp")
        try:
            result.to_npz(temporary_archive)
            digest = _file_digest(temporary_archive)
            temporary_checksum.write_text(
                digest + "\n",
                encoding="ascii",
            )
            os.replace(temporary_archive, archive)
            os.replace(temporary_checksum, checksum)
        finally:
            for path in (temporary_archive, temporary_checksum):
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass

    def _quarantine(self, key: str) -> None:
        archive, checksum = self._paths(key)
        token = uuid.uuid4().hex
        for path in (archive, checksum):
            if not path.exists():
                continue
            target = self._corrupt / f"{path.name}.{token}.corrupt"
            try:
                os.replace(path, target)
            except FileNotFoundError:
                pass
