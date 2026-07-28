"""Readers for Zonge TEMAVG processing ``.LOG`` files."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Union

from ..api.property import PyCSAMTObject

PathLike = Union[str, Path]

__all__ = [
    "TEMLog",
    "TEMLogRecord",
    "is_tem_log_file",
]


_VERSION_RE = re.compile(
    r"TEMAVG\s+(?P<version>[\d.]+),\s+Processed:\s+(?P<processed>.+)"
)
_READING_RE = re.compile(
    r'Reading\s+"(?P<source>[^"]+)"',
    re.IGNORECASE,
)
_CLOCK_RE = re.compile(
    r"Data uses\s+(?P<clock>[A-Z]+)\s+CLOCK\s+"
    r"\((?P<resolution>[-+\d.]+)\s*(?P<unit>[A-Za-z]+)\)",
    re.IGNORECASE,
)
_DATA_CLOCK_RE = re.compile(
    r"Data uses\s+(?P<frequency>[-+\d.]+)\s*MHz\s+"
    r"\((?P<clock>[A-Za-z]+)\)\s+clock\s+"
    r"\((?P<resolution>[-+\d.]+)\s*(?P<unit>[A-Za-z]+)\)",
    re.IGNORECASE,
)
_AVG_COUNT_RE = re.compile(
    r'File\s+"(?P<file>[^"]+)"\s+contains averaged data for\s+'
    r"(?P<count>\d+)\s+data sets",
    re.IGNORECASE,
)
_CLOSED_RE = re.compile(r'Log file\s+"(?P<file>[^"]+)"\s+closed', re.I)


@dataclass(frozen=True)
class TEMLogRecord:
    """One acquisition-summary row from a TEMAVG log.

    Parameters
    ----------
    station : float
        Station coordinate or station number.
    frequency : float
        Repetition frequency in hertz.
    loop : str
        Loop geometry or acquisition code, for example
        ``"INL"`` for in-loop data.
    component : str
        Measured component label.
    duty : str
        Duty-cycle string written by TEMAVG, for example
        ``"50%"``.
    first_window : float
        First gate centre time in microseconds, parsed from the
        ``Window1`` column.
    rx_moment : float
        Receiver moment or area value from the ``RxMoment``
        column.
    time_base : str
        Time-base or clock sample label from the ``Ts`` column.
    """

    station: float
    frequency: float
    loop: str
    component: str
    duty: str
    first_window: float
    rx_moment: float
    time_base: str


@dataclass(repr=False)
class TEMLog(PyCSAMTObject):
    """Parsed TEMAVG processing log.

    ``TEMLog`` stores processing provenance emitted by the
    Zonge TEMAVG program. It captures the stable ASCII metadata
    and acquisition-summary table while preserving the raw
    processing-mode text for later inspection.

    Parameters
    ----------
    path : pathlib.Path
        Source ``.LOG`` file.
    metadata : dict
        Parsed metadata such as source field file, clock type,
        clock resolution, output AVG file, data-set count, and
        close status.
    records : list of TEMLogRecord
        Parsed acquisition-summary rows.
    version : str, optional
        TEMAVG program version.
    processed : str, optional
        Processing date string.
    raw_modes : list of str
        Lines from the TEMAVG global and processing-mode
        sections. The original table uses legacy box-drawing
        characters, so it is preserved as text.
    """

    path: Path
    metadata: dict[str, Any] = field(default_factory=dict)
    records: list[TEMLogRecord] = field(default_factory=list)
    version: str | None = None
    processed: str | None = None
    raw_modes: list[str] = field(default_factory=list)
    verbose: int = 0
    logger: object | None = None

    @classmethod
    def read(
        cls,
        path: PathLike,
        *,
        verbose: int = 0,
        logger: object | None = None,
    ) -> TEMLog:
        """Read a TEMAVG ``.LOG`` processing file.

        Parameters
        ----------
        path : path-like
            TEMAVG log file produced while averaging and writing
            ``.AVG`` and ``.Z`` outputs.

        Returns
        -------
        TEMLog
            Parsed processing log.
        """
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(p)

        metadata: dict[str, Any] = {}
        records: list[TEMLogRecord] = []
        raw_modes: list[str] = []
        version = processed = None
        in_records = False
        in_modes = False

        with p.open("r", errors="replace") as fh:
            for line_no, raw in enumerate(fh, 1):
                line = raw.rstrip("\n")
                stripped = line.strip()
                if not stripped:
                    continue
                if line_no == 1:
                    match = _VERSION_RE.search(stripped)
                    if match:
                        version = match.group("version")
                        processed = match.group("processed").strip()
                    metadata["title"] = stripped
                    continue

                if stripped == "GLOBAL MODE LIST:":
                    in_modes = True
                if stripped.startswith("Reading "):
                    in_modes = False
                    _update_reading_metadata(stripped, metadata)
                    continue
                if in_modes:
                    raw_modes.append(stripped)
                    continue

                if stripped.startswith("[CLOCK]"):
                    _update_clock_metadata(stripped, metadata)
                    continue
                if stripped.startswith("Blk"):
                    in_records = True
                    continue
                if in_records and stripped.startswith("AVG"):
                    record = _parse_log_record(stripped, p, line_no)
                    records.append(record)
                    continue
                if stripped.startswith("Combine "):
                    steps = metadata.setdefault("combine_steps", [])
                    steps.append(stripped)
                    continue
                _update_output_metadata(stripped, metadata)

        return cls(
            path=p,
            metadata=metadata,
            records=records,
            version=version,
            processed=processed,
            raw_modes=raw_modes,
            verbose=verbose,
            logger=logger,
        )

    @property
    def n_records(self) -> int:
        """Number of acquisition-summary rows."""
        return len(self.records)

    @property
    def stations(self) -> list[float]:
        """Sorted station values represented in the log."""
        return sorted({rec.station for rec in self.records})

    def to_records(self) -> list[dict[str, Any]]:
        """Return acquisition rows as dictionaries."""
        rows: list[dict[str, Any]] = []
        for rec in self.records:
            rows.append(
                {
                    "source_file": self.path.name,
                    "station": rec.station,
                    "frequency": rec.frequency,
                    "loop": rec.loop,
                    "component": rec.component,
                    "duty": rec.duty,
                    "first_window_us": rec.first_window,
                    "rx_moment": rec.rx_moment,
                    "time_base": rec.time_base,
                    "field_file": self.metadata.get("field_file"),
                }
            )
        return rows

    def to_dataframe(self):
        """Return the acquisition-summary table as a DataFrame."""
        try:
            import pandas as pd
        except ImportError as exc:
            msg = "TEMLog.to_dataframe requires pandas."
            raise ImportError(msg) from exc
        return pd.DataFrame(self.to_records())


def is_tem_log_file(path: PathLike) -> bool:
    """Return ``True`` when ``path`` looks like a TEMAVG log."""
    p = Path(path)
    if not p.exists() or not p.is_file():
        return False
    try:
        with p.open("r", errors="replace") as fh:
            head = "".join(fh.readline() for _ in range(20))
    except OSError:
        return False
    return "TEMAVG" in head and "GLOBAL MODE LIST" in head


def _update_reading_metadata(
    line: str,
    metadata: dict[str, Any],
) -> None:
    """Parse the source field file line."""
    match = _READING_RE.search(line)
    if match:
        metadata["field_file"] = match.group("source")


def _update_clock_metadata(
    line: str,
    metadata: dict[str, Any],
) -> None:
    """Parse the first clock metadata line."""
    match = _CLOCK_RE.search(line)
    if match:
        metadata["clock_type"] = match.group("clock").upper()
        metadata["clock_resolution"] = float(match.group("resolution"))
        metadata["clock_unit"] = match.group("unit")


def _update_output_metadata(
    line: str,
    metadata: dict[str, Any],
) -> None:
    """Parse output and close-status metadata lines."""
    avg_count = _AVG_COUNT_RE.search(line)
    if avg_count:
        metadata["avg_file"] = avg_count.group("file")
        metadata["avg_dataset_count"] = int(avg_count.group("count"))
        return

    data_clock = _DATA_CLOCK_RE.search(line)
    if data_clock:
        frequency = data_clock.group("frequency")
        metadata["data_clock_mhz"] = float(frequency)
        metadata["data_clock_type"] = data_clock.group("clock").lower()
        metadata["data_clock_resolution"] = float(data_clock.group("resolution"))
        metadata["data_clock_unit"] = data_clock.group("unit")
        return

    closed = _CLOSED_RE.search(line)
    if closed:
        metadata["closed_file"] = closed.group("file")
        metadata["closed"] = True


def _parse_log_record(
    line: str,
    path: Path,
    line_no: int,
) -> TEMLogRecord:
    """Parse one TEMAVG log acquisition-summary row."""
    parts = line.split()
    if len(parts) < 9:
        msg = f"Malformed TEMAVG log row in {path!s} at line {line_no}."
        raise ValueError(msg)
    try:
        first_window = _parse_window_us(parts[6])
        return TEMLogRecord(
            station=float(parts[1]),
            frequency=float(parts[2]),
            loop=parts[3],
            component=parts[4],
            duty=parts[5],
            first_window=first_window,
            rx_moment=float(parts[7]),
            time_base=parts[8],
        )
    except ValueError as exc:
        msg = f"Cannot parse TEMAVG log row in {path!s} at line {line_no}."
        raise ValueError(msg) from exc


def _parse_window_us(value: str) -> float:
    """Parse a TEMAVG window value and return microseconds."""
    match = re.match(r"([-+\d.]+)\s*([A-Za-z]*)", value)
    if match is None:
        return float(value)
    number = float(match.group(1))
    unit = match.group(2).lower()
    if unit in {"u", "us", "usec"}:
        return number
    if unit in {"m", "ms", "msec"}:
        return number * 1000.0
    if unit in {"s", "sec"}:
        return number * 1_000_000.0
    return number
