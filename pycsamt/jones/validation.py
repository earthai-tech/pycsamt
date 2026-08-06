# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

import os
from abc import ABC, abstractmethod
from pathlib import Path

from ..exceptions import (
    FileHandlingError,
    JError,
)
from .config import (
    RE_BANNER,
    RE_BLANK,
    RE_COMMENT,
    RE_DATATYPE_UNITS,
    RE_INFO,
    RE_NPOINTS,
    RE_ROW_R,
    RE_ROW_TF,
    RE_STATION,
)

__all__ = ["IsJ", "is_j_file"]


def _read_lines(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    return text.splitlines()


def _is_comment(s: str) -> bool:
    return bool(RE_COMMENT.match(s))


def _is_blank(s: str) -> bool:
    return bool(RE_BLANK.match(s))


def _is_info(s: str) -> bool:
    return bool(RE_INFO.match(s))


def _first_nonblank(lines: list[str], i: int) -> int:
    n = len(lines)
    j = i
    while j < n and _is_blank(lines[j]):
        j += 1
    return j


def _find_dtype_after(lines: list[str], i: int) -> int | None:
    n = len(lines)
    j = _first_nonblank(lines, i + 1)
    if j >= n:
        return None
    # dtype first
    if RE_DATATYPE_UNITS.match(lines[j]):
        return j
    # or count first then dtype
    if RE_NPOINTS.match(lines[j]):
        k = j + 1
        while k < n:
            if RE_STATION.match(lines[k]):
                break
            if _is_blank(lines[k]):
                k += 1
                continue
            if RE_DATATYPE_UNITS.match(lines[k]):
                return k
            k += 1
    return None


def _find_count_after(lines: list[str], i: int) -> int | None:
    n = len(lines)
    j = _first_nonblank(lines, i + 1)
    if j >= n:
        return None
    # count first
    if RE_NPOINTS.match(lines[j]):
        return j
    # or dtype first then count
    if RE_DATATYPE_UNITS.match(lines[j]):
        k = j + 1
        while k < n and _is_blank(lines[k]):
            k += 1
        if k < n and RE_NPOINTS.match(lines[k]):
            return k
    return None


def _guess_kind(line: str) -> str | None:
    m = RE_DATATYPE_UNITS.match(line)
    if not m:
        return None
    return m.group("kind").upper()


def _has_at_least_one_row(lines: list[str], count_idx: int, kind: str) -> bool:
    n = len(lines)
    # next non-blank after count is the first row
    i = _first_nonblank(lines, count_idx + 1)
    if i >= n:
        return False
    # decide row regex by kind
    if kind in {"R", "S"}:
        return bool(RE_ROW_R.match(lines[i]))
    if kind in {"Z", "Q", "C", "T"}:
        return bool(RE_ROW_TF.match(lines[i]))
    # unknown kind: be conservative
    return False


def _scan_j_blocks(lines: list[str]) -> int:
    i, n = 0, len(lines)

    # Skip banner/comments/info/blank before first station
    while i < n and (
        _is_comment(lines[i])
        or _is_info(lines[i])
        or _is_blank(lines[i])
        or RE_BANNER.match(lines[i])
    ):
        i += 1

    n_ok = 0
    while i < n:
        m = RE_STATION.match(lines[i])
        if not m:
            i += 1
            continue

        # Find dtype and count in any accepted order
        j = _find_dtype_after(lines, i)
        k = _find_count_after(lines, i)
        if j is None or k is None:
            # advance to avoid infinite loop
            i += 1
            continue

        kind = _guess_kind(lines[j]) or ""
        if not kind:
            i += 1
            continue

        if not _has_at_least_one_row(lines, k, kind):
            i += 1
            continue

        n_ok += 1

        # advance roughly past this block using count
        try:
            cn = int(RE_NPOINTS.match(lines[k]).group("n"))  # type: ignore
        except Exception:
            cn = 0
        # station + (dtype|count) + (count|dtype) + rows
        i = max(k + 1 + max(cn, 0), i + 1)

    return n_ok


class IsJ(ABC):
    r"""
    Abstract base for A.G. Jones J-format validation helpers.

    Subclasses implement :attr:`IsJ.is_valid`.  The static method
    :py:meth:`_assert_j` provides a robust, file-level validator
    that accepts either a path or an existing :class:`IsJ`
    instance.  The check is heuristic, fast and tolerant of
    minor formatting issues (e.g., non-canonical head order).
    """

    @property
    @abstractmethod
    def is_valid(self) -> bool:
        r"""Return ``True`` when the instance is structurally valid."""
        raise NotImplementedError

    @staticmethod
    def _assert_j(
        file: str | os.PathLike | IsJ,
        *,
        deep: bool = True,
    ) -> bool:
        r"""
        Validate that ``file`` is a Jones J-format file (or obj).

        When ``deep`` is ``False``, only a light extension check
        is performed (``.j``, ``.jones``, ``.txt``,``.dat`` accepted).
        When ``deep`` is ``True`` (default), the function scans
        the text and requires at least one complete data block:
        a station line, a data-type line, a count line and a row
        that matches the expected numeric layout.
        """
        if file is None:
            raise FileHandlingError("NoneType not allowed")

        if isinstance(file, IsJ):
            if not file.is_valid:
                raise JError("J object reports invalid state")
            return True

        path = Path(file)  # type: ignore[arg-type]
        if not path.exists() or not path.is_file():
            raise FileNotFoundError(f"{str(path)!r} not a file")

        if not deep:
            ext = path.suffix.lower().lstrip(".")
            if ext in {"j", "jones", "txt", "dat"}:
                return True
            raise JError("Unexpected extension; expected .j/.jones/.txt/.dat")

        try:
            lines = _read_lines(path)
        except PermissionError:
            raise PermissionError("Permission denied while reading file")

        # Heuristic header checks (optional but helpful)
        # - zero or more info lines at the top are OK
        # - optional banner
        # We only require at least one valid block.
        n_blocks = _scan_j_blocks(lines)
        if n_blocks <= 0:
            raise JError(
                "Unrecognized Jones J file structure: found no "
                "valid data blocks (station + dtype + count + row)"
            )

        return True


def is_j_file(file: str | os.PathLike | IsJ, *, deep: bool = True) -> bool:
    """
    Convenience wrapper around :meth:`IsJ._assert_j`.
    """
    return IsJ._assert_j(file, deep=deep)
