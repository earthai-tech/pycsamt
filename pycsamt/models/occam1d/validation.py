# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Safe native-file recognition and validation for Occam1D.

Recognition is deliberately separate from scientific parsing. This module
performs bounded, side-effect-free inspection suitable for directory scans.
The owning format classes remain the authority for complete syntax and
scientific invariants.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from enum import Enum
from os import PathLike as OSPathLike
from pathlib import Path
from typing import Any, Union

PathLike = Union[str, OSPathLike]

__all__ = [
    "Occam1DFileReport",
    "Occam1DFileType",
    "Occam1DValidationStatus",
    "detect_file_type",
    "inspect_occam1d_file",
    "is_data_file",
    "is_iter_file",
    "is_log_file",
    "is_model_file",
    "is_response_file",
    "is_startup_file",
    "scan_occam1d_directory",
    "validate_occam1d_file",
]

_DEFAULT_MAX_LINES = 256
_DEFAULT_MAX_BYTES = 256 * 1024
_FORMAT_DATA = re.compile(r"\bemdata_1\.\d+\b", re.I)
_FORMAT_MODEL = re.compile(r"\bresistivity1dmod_1\.\d+\b", re.I)
_FORMAT_ITER = re.compile(r"\boccamiter_flex\b", re.I)
_HEADER = re.compile(r"^\s*([^:]+?)\s*:\s*(.*?)\s*$")
_LOG_ITERATION = re.compile(r"\*+\s*iteration\s+[0-9]+\s*\*+", re.I)
_LOG_RMS = re.compile(
    r"\b(?:(?:starting\s+)?r\.?m\.?s\.?|and\s+is)\b",
    re.I,
)
_KNOWN_RESPONSE_TYPES = frozenset({103, 104})


class Occam1DFileType(str, Enum):
    """Canonical categories returned by :func:`detect_file_type`.

    The enumeration derives from :class:`str`, preserving comparisons with
    historical values such as ``detect_file_type(path) == "data"``.
    """

    DATA = "data"
    MODEL = "model"
    STARTUP = "startup"
    ITER = "iter"
    RESPONSE = "response"
    LOG = "log"
    UNKNOWN = "unknown"


class Occam1DValidationStatus(str, Enum):
    """Outcome of bounded native-file inspection."""

    VALID = "valid"
    UNKNOWN = "unknown"
    AMBIGUOUS = "ambiguous"
    MISSING = "missing"
    NOT_FILE = "not_file"
    UNREADABLE = "unreadable"
    INVALID = "invalid"


@dataclass(frozen=True)
class Occam1DFileReport:
    """Immutable evidence from inspecting one filesystem candidate.

    Parameters
    ----------
    path : pathlib.Path
        Normalized absolute candidate path.
    file_type : Occam1DFileType
        Selected format category, or ``UNKNOWN``.
    status : Occam1DValidationStatus
        Inspection outcome.
    matches : tuple of Occam1DFileType
        Every category whose bounded signature matched.
    reasons : tuple of str
        Human-readable evidence or failure descriptions.
    iteration : int or None
        Native iteration header for startup/iteration files.
    lines_examined : int
        Number of decoded lines inspected.
    truncated : bool
        Whether configured line or byte limits bounded inspection.

    Notes
    -----
    ``valid`` means the signature is internally coherent. It does not replace
    full parsing by :class:`Occam1DData`, :class:`Occam1DModel`, or another
    native scientific object.
    """

    path: Path
    file_type: Occam1DFileType
    status: Occam1DValidationStatus
    matches: tuple[Occam1DFileType, ...] = ()
    reasons: tuple[str, ...] = ()
    iteration: int | None = None
    lines_examined: int = 0
    truncated: bool = False

    @property
    def valid(self) -> bool:
        """Whether one unambiguous native signature was recognized."""
        return self.status is Occam1DValidationStatus.VALID

    @property
    def ambiguous(self) -> bool:
        """Whether mutually incompatible signatures were recognized."""
        return self.status is Occam1DValidationStatus.AMBIGUOUS

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly independent representation."""
        return {
            "path": str(self.path),
            "file_type": self.file_type.value,
            "status": self.status.value,
            "valid": self.valid,
            "matches": [value.value for value in self.matches],
            "reasons": list(self.reasons),
            "iteration": self.iteration,
            "lines_examined": self.lines_examined,
            "truncated": self.truncated,
        }


def _normalize_path(path: PathLike) -> Path:
    """Return an absolute inspection path with explicit input validation."""
    if not isinstance(path, (str, OSPathLike)):
        raise TypeError("path must be a string or os.PathLike object.")
    if isinstance(path, str) and not path.strip():
        raise ValueError("path cannot be empty.")
    return Path(path).expanduser().resolve(strict=False)


def _positive_limit(name: str, value: int) -> int:
    """Validate a positive integer resource bound."""
    if not isinstance(value, int) or isinstance(value, bool):
        raise TypeError(f"{name} must be a positive integer.")
    if value < 1:
        raise ValueError(f"{name} must be strictly positive.")
    return value


def _read_bounded(
    path: Path,
    *,
    max_lines: int,
    max_bytes: int,
) -> tuple[list[str], bool]:
    """Read bounded UTF-8-compatible text without loading a large file."""
    lines = []
    consumed = 0
    truncated = False
    with path.open("rb") as stream:
        while len(lines) < max_lines:
            raw = stream.readline(max_bytes - consumed + 1)
            if not raw:
                break
            consumed += len(raw)
            if consumed > max_bytes:
                raw = raw[: max(0, len(raw) - (consumed - max_bytes))]
                truncated = True
            if b"\x00" in raw:
                raise ValueError("Candidate contains NUL bytes and is binary.")
            lines.append(raw.decode("utf8", errors="replace").rstrip("\r\n"))
            if truncated:
                break
        if len(lines) == max_lines and stream.read(1):
            truncated = True
    return lines, truncated


def _headers(lines: list[str]) -> dict[str, str]:
    """Return normalized native colon headers from bounded text."""
    values = {}
    for line in lines:
        match = _HEADER.match(line)
        if match:
            key = match.group(1).strip().lstrip("#! ").strip().lower()
            values[key] = match.group(2).strip()
    return values


def _integer_header(headers, name: str) -> int | None:
    """Return a strict integer header value, or ``None``."""
    value = headers.get(name)
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _numeric_fields(line: str) -> list[float] | None:
    """Parse a wholly numeric non-comment row."""
    text = line.strip()
    if not text or text.startswith(("!", "#")):
        return None
    try:
        values = [float(value.replace("D", "E").replace("d", "e"))
                  for value in text.split()]
    except ValueError:
        return None
    return values if values and all(math.isfinite(v) for v in values) else None


def _response_rows(lines: list[str]) -> int:
    """Count structurally plausible canonical or EMData response rows."""
    count = 0
    for line in lines:
        values = _numeric_fields(line)
        if values is None:
            continue
        if len(values) == 7:
            first = values[0]
            if first.is_integer() and int(first) in _KNOWN_RESPONSE_TYPES:
                frequency = values[1]
                code = first
            else:
                frequency = values[1]
                code = values[2]
            if (
                frequency.is_integer()
                and frequency >= 1
                and code.is_integer()
                and int(code) >= 100
            ):
                count += 1
        elif len(values) >= 8:
            code, frequency = values[:2]
            if (
                code.is_integer()
                and int(code) in _KNOWN_RESPONSE_TYPES
                and frequency.is_integer()
                and frequency >= 1
            ):
                count += 1
    return count


def _classify(lines: list[str]):
    """Return matched types, reasons, and an optional iteration number."""
    text = "\n".join(lines)
    headers = _headers(lines)
    matches = []
    reasons = []
    iteration = _integer_header(headers, "iteration")

    if _FORMAT_DATA.search(text):
        if "frequencies" in headers and "data" in headers:
            matches.append(Occam1DFileType.DATA)
            reasons.append("EMData format and required count headers found.")
        else:
            reasons.append("EMData token lacks Frequencies or Data header.")
    if _FORMAT_MODEL.search(text):
        if "layers" in headers:
            matches.append(Occam1DFileType.MODEL)
            reasons.append("Resistivity1DMod format and Layers header found.")
        else:
            reasons.append("Resistivity1DMod token lacks Layers header.")
    if _FORMAT_ITER.search(text):
        required = {"model file", "data file", "iteration"}
        if required.issubset(headers) and iteration is not None:
            kind = (
                Occam1DFileType.STARTUP
                if iteration == 0
                else Occam1DFileType.ITER
            )
            if iteration >= 0:
                matches.append(kind)
                reasons.append(
                    f"OCCAMITER_FLEX headers declare iteration {iteration}."
                )
        else:
            reasons.append("OCCAMITER_FLEX token lacks required headers.")

    response_count = _response_rows(lines)
    if response_count:
        matches.append(Occam1DFileType.RESPONSE)
        reasons.append(f"Found {response_count} plausible response row(s).")
    if _LOG_ITERATION.search(text) and _LOG_RMS.search(text):
        matches.append(Occam1DFileType.LOG)
        reasons.append("Iteration block and RMS log markers found.")
    return tuple(dict.fromkeys(matches)), tuple(reasons), iteration


def inspect_occam1d_file(
    path: PathLike,
    *,
    max_lines: int = _DEFAULT_MAX_LINES,
    max_bytes: int = _DEFAULT_MAX_BYTES,
) -> Occam1DFileReport:
    """Inspect one candidate using bounded native-format signatures.

    Missing paths, directories, binary content, and permission failures are
    represented in the returned report instead of raising. Invalid public
    arguments still raise immediately.
    """
    candidate = _normalize_path(path)
    max_lines = _positive_limit("max_lines", max_lines)
    max_bytes = _positive_limit("max_bytes", max_bytes)
    if not candidate.exists():
        return Occam1DFileReport(
            candidate,
            Occam1DFileType.UNKNOWN,
            Occam1DValidationStatus.MISSING,
            reasons=("Candidate does not exist.",),
        )
    if not candidate.is_file():
        return Occam1DFileReport(
            candidate,
            Occam1DFileType.UNKNOWN,
            Occam1DValidationStatus.NOT_FILE,
            reasons=("Candidate is not a regular file.",),
        )
    try:
        lines, truncated = _read_bounded(
            candidate,
            max_lines=max_lines,
            max_bytes=max_bytes,
        )
    except PermissionError:
        return Occam1DFileReport(
            candidate,
            Occam1DFileType.UNKNOWN,
            Occam1DValidationStatus.UNREADABLE,
            reasons=("Candidate is not readable.",),
        )
    except (OSError, ValueError) as error:
        return Occam1DFileReport(
            candidate,
            Occam1DFileType.UNKNOWN,
            Occam1DValidationStatus.INVALID,
            reasons=(str(error),),
        )
    matches, reasons, iteration = _classify(lines)
    if len(matches) == 1:
        kind = matches[0]
        status = Occam1DValidationStatus.VALID
    elif len(matches) > 1:
        kind = Occam1DFileType.UNKNOWN
        status = Occam1DValidationStatus.AMBIGUOUS
    else:
        kind = Occam1DFileType.UNKNOWN
        status = Occam1DValidationStatus.UNKNOWN
        if not reasons:
            reasons = ("No native Occam1D signature was recognized.",)
    return Occam1DFileReport(
        path=candidate,
        file_type=kind,
        status=status,
        matches=matches,
        reasons=reasons,
        iteration=iteration,
        lines_examined=len(lines),
        truncated=truncated,
    )


def validate_occam1d_file(
    path: PathLike,
    expected: Occam1DFileType | str | None = None,
    *,
    deep: bool = False,
) -> Occam1DFileReport:
    """Validate a signature and optionally invoke its full native parser.

    Parameters
    ----------
    path : path-like
        Candidate native text file.
    expected : Occam1DFileType or str, optional
        Required category. A mismatch returns an ``INVALID`` report.
    deep : bool, default=False
        Run the corresponding scientific parser after signature inspection.
        Parser failures are captured in the report.

    Returns
    -------
    Occam1DFileReport
        Immutable validation evidence. This function does not mutate files.
    """
    if not isinstance(deep, bool):
        raise TypeError("deep must be a bool.")
    report = inspect_occam1d_file(path)
    if expected is not None:
        try:
            expected_type = Occam1DFileType(expected)
        except (TypeError, ValueError):
            choices = ", ".join(value.value for value in Occam1DFileType)
            raise ValueError(f"expected must be one of: {choices}.") from None
        if report.file_type is not expected_type:
            return _invalid_report(
                report,
                f"Expected {expected_type.value}, detected "
                f"{report.file_type.value}.",
            )
    if not deep or not report.valid:
        return report
    try:
        _deep_parse(report.file_type, report.path)
    except (OSError, TypeError, ValueError) as error:
        return _invalid_report(report, f"Full parser rejected file: {error}")
    return report


def _invalid_report(
    report: Occam1DFileReport,
    reason: str,
) -> Occam1DFileReport:
    """Return ``report`` evidence with an invalid outcome."""
    return Occam1DFileReport(
        path=report.path,
        file_type=report.file_type,
        status=Occam1DValidationStatus.INVALID,
        matches=report.matches,
        reasons=report.reasons + (reason,),
        iteration=report.iteration,
        lines_examined=report.lines_examined,
        truncated=report.truncated,
    )


def _deep_parse(kind: Occam1DFileType, path: Path) -> None:
    """Dispatch locally imported authoritative parsers without cycles."""
    if kind is Occam1DFileType.DATA:
        from .data import Occam1DData

        Occam1DData.read(path)
    elif kind is Occam1DFileType.MODEL:
        from .model import Occam1DModel

        Occam1DModel.read(path)
    elif kind in {Occam1DFileType.STARTUP, Occam1DFileType.ITER}:
        from .startup import Occam1DStartup

        Occam1DStartup.read(path)
    elif kind is Occam1DFileType.RESPONSE:
        from .response import Occam1DResponse

        Occam1DResponse.read(path)
    elif kind is Occam1DFileType.LOG:
        from .log import Occam1DLog

        Occam1DLog.read(path)


def scan_occam1d_directory(
    directory: PathLike,
    *,
    include_unknown: bool = False,
    recursive: bool = False,
) -> tuple[Occam1DFileReport, ...]:
    """Return deterministic reports for regular files in a directory.

    Symbolic-link directories are never traversed. The default scan inspects
    direct children only, matching Occam1D run-directory semantics.
    """
    if not isinstance(include_unknown, bool):
        raise TypeError("include_unknown must be a bool.")
    if not isinstance(recursive, bool):
        raise TypeError("recursive must be a bool.")
    root = _normalize_path(directory)
    if not root.exists():
        raise FileNotFoundError(f"Occam1D directory not found: {root}")
    if not root.is_dir():
        raise NotADirectoryError(
            f"Occam1D scan target is not a directory: {root}"
        )
    candidates = root.rglob("*") if recursive else root.iterdir()
    reports = []
    for path in sorted(candidates, key=lambda item: str(item).lower()):
        if path.is_symlink() or not path.is_file():
            continue
        report = inspect_occam1d_file(path)
        if include_unknown or report.valid or report.ambiguous:
            reports.append(report)
    return tuple(reports)


def detect_file_type(path: PathLike) -> Occam1DFileType:
    """Return the unambiguous native category or ``UNKNOWN``.

    This compatibility function never raises for a missing, unreadable,
    binary, malformed, or ambiguous filesystem candidate.
    """
    return inspect_occam1d_file(path).file_type


def _is_type(path: PathLike, expected: Occam1DFileType) -> bool:
    """Return whether bounded inspection recognizes ``expected``."""
    try:
        return detect_file_type(path) is expected
    except (TypeError, ValueError):
        return False


def is_data_file(path: PathLike) -> bool:
    """Return whether ``path`` is an unambiguous native data file."""
    return _is_type(path, Occam1DFileType.DATA)


def is_model_file(path: PathLike) -> bool:
    """Return whether ``path`` is an unambiguous native layer-model file."""
    return _is_type(path, Occam1DFileType.MODEL)


def is_startup_file(path: PathLike) -> bool:
    """Return whether ``path`` declares OCCAMITER_FLEX iteration zero."""
    return _is_type(path, Occam1DFileType.STARTUP)


def is_iter_file(path: PathLike) -> bool:
    """Return whether ``path`` declares a positive solver iteration."""
    return _is_type(path, Occam1DFileType.ITER)


def is_response_file(path: PathLike) -> bool:
    """Return whether ``path`` contains plausible native response rows."""
    return _is_type(path, Occam1DFileType.RESPONSE)


def is_log_file(path: PathLike) -> bool:
    """Return whether ``path`` contains Occam iteration and RMS markers."""
    return _is_type(path, Occam1DFileType.LOG)
