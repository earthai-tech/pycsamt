"""Readers for Zonge TEMAVG contour ``.Z`` files."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Union

from ..api.property import PyCSAMTObject

PathLike = Union[str, Path]

__all__ = [
    "TEMZPlot",
    "TEMZRecord",
    "is_tem_z_file",
]


_FLOAT_RE = re.compile(
    r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"
)
_TAIL_RE = re.compile(
    r"^(?P<prefix>.*?)\s+"
    r"(?P<frequency>[-+]?\d+(?:\.\d+)?)\s*Hz\s+"
    r"W\s*(?P<window>\d+)\s*$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class TEMZRecord:
    """One row from a TEMAVG contour ``.Z`` file.

    Parameters
    ----------
    line : int
        Line or profile identifier stored in the first column.
    station : float
        Station coordinate or station number.
    time : float
        Gate centre time in milliseconds.
    magnitude : float
        Contoured transient magnitude, commonly in microvolts
        per ampere for TEMAVG output.
    frequency : float
        Repetition frequency in hertz.
    window : int
        Time-window number.
    """

    line: int
    station: float
    time: float
    magnitude: float
    frequency: float
    window: int

    @property
    def time_s(self) -> float:
        """Gate centre time in seconds."""
        return self.time * 1e-3


@dataclass(repr=False)
class TEMZPlot(PyCSAMTObject):
    """Parsed content of one TEMAVG contour ``.Z`` file.

    ``.Z`` files are plotting-oriented TEMAVG outputs. They
    store station, time gate, magnitude, frequency, and window
    values in a compact fixed-width text layout. They are useful
    for quick pseudo-section plotting and for checking that the
    processed ``.AVG`` magnitude column was exported correctly.

    Parameters
    ----------
    path : pathlib.Path
        Source ``.Z`` file.
    metadata : dict
        Parsed header metadata such as date, plot data type,
        window name, units, component, receiver area, contour
        mode, and profile mode when present.
    records : list of TEMZRecord
        Parsed plot rows.
    version : str, optional
        TEMAVG program version parsed from the first line.
    processed : str, optional
        Processing date parsed from the second line.

    Examples
    --------
    >>> from pycsamt.tdem.zplot import TEMZPlot
    >>> z = TEMZPlot.read("data/TEMAVG/JIANGSU/TEM100.Z")
    >>> z.n_records > 0
    True
    >>> z.records[0].window
    1
    """

    path: Path
    metadata: dict[str, Any] = field(default_factory=dict)
    records: list[TEMZRecord] = field(default_factory=list)
    version: str | None = None
    processed: str | None = None
    verbose: int = 0
    logger: object | None = None

    @classmethod
    def read(
        cls,
        path: PathLike,
        *,
        verbose: int = 0,
        logger: object | None = None,
    ) -> TEMZPlot:
        """Read a TEMAVG contour ``.Z`` file.

        Parameters
        ----------
        path : path-like
            TEMAVG contour file produced by the ``ZPLOT`` stage.

        Returns
        -------
        TEMZPlot
            Parsed contour records and header metadata.
        """
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(p)

        metadata: dict[str, Any] = {}
        records: list[TEMZRecord] = []
        version = processed = None
        in_table = False

        with p.open("r", errors="replace") as fh:
            for line_no, raw in enumerate(fh, 1):
                line = raw.strip()
                if not line:
                    continue
                if line_no == 1:
                    version = _parse_version(line)
                    metadata["title"] = line
                    continue
                if line.startswith("/*"):
                    processed = line.lstrip("/*").strip()
                    metadata["processed"] = processed
                    continue
                if line.startswith("$"):
                    _update_dollar_metadata(line, metadata)
                    continue
                if line.upper().startswith("WINDOW "):
                    metadata["window"] = line.split(None, 1)[1].strip()
                    continue
                if "VALUES IN" in line.upper():
                    metadata["units"] = line.strip()
                    continue
                if line.upper().startswith("COMPONENT:"):
                    _update_component_metadata(line, metadata)
                    continue
                if line.startswith("II"):
                    in_table = True
                    continue
                if in_table and _looks_like_z_record(line):
                    records.append(_parse_z_record(line, p, line_no))

        if not records:
            msg = f"No TEMAVG .Z data rows found in {p!s}."
            raise ValueError(msg)
        return cls(
            path=p,
            metadata=metadata,
            records=records,
            version=version,
            processed=processed,
            verbose=verbose,
            logger=logger,
        )

    @property
    def n_records(self) -> int:
        """Number of contour data rows."""
        return len(self.records)

    @property
    def stations(self) -> list[float]:
        """Sorted station values represented in the file."""
        return sorted({rec.station for rec in self.records})

    @property
    def windows(self) -> list[int]:
        """Sorted time-window numbers represented in the file."""
        return sorted({rec.window for rec in self.records})

    def rows_for_station(self, station: float) -> list[TEMZRecord]:
        """Return all contour rows for one station value."""
        return [rec for rec in self.records if rec.station == station]

    def to_records(self) -> list[dict[str, Any]]:
        """Return contour records as dictionaries."""
        rows: list[dict[str, Any]] = []
        for rec in self.records:
            rows.append(
                {
                    "source_file": self.path.name,
                    "line": rec.line,
                    "station": rec.station,
                    "time_ms": rec.time,
                    "time_s": rec.time_s,
                    "magnitude": rec.magnitude,
                    "frequency": rec.frequency,
                    "window": rec.window,
                    "component": self.metadata.get("component"),
                    "units": self.metadata.get("units"),
                    "rx_area": self.metadata.get("rx_area"),
                }
            )
        return rows

    def to_dataframe(self):
        """Return the parsed contour table as a DataFrame."""
        try:
            import pandas as pd
        except ImportError as exc:
            msg = "TEMZPlot.to_dataframe requires pandas."
            raise ImportError(msg) from exc
        return pd.DataFrame(self.to_records())


def is_tem_z_file(path: PathLike) -> bool:
    """Return ``True`` when ``path`` looks like a TEMAVG ``.Z`` file."""
    p = Path(path)
    if not p.exists() or not p.is_file():
        return False
    try:
        with p.open("r", errors="replace") as fh:
            head = "".join(fh.readline() for _ in range(14))
    except OSError:
        return False
    return "TEMAVG" in head and "ZPLOT" in head and "IIxxxxxxxx" in head


def _parse_version(line: str) -> str | None:
    """Parse a TEMAVG version from the first header line."""
    match = re.search(r"TEMAVG\s+([\d.]+)", line, re.IGNORECASE)
    return match.group(1) if match else None


def _update_dollar_metadata(
    line: str,
    metadata: dict[str, Any],
) -> None:
    """Update metadata from a dollar-prefixed header line."""
    content = line.lstrip("$").strip()
    if "=" not in content:
        return
    key, value = [part.strip() for part in content.split("=", 1)]
    if key.upper() == "DATE":
        metadata["date"] = value
    elif key.upper().startswith("ZPLOT:"):
        metadata["zplot_data"] = value
    else:
        metadata[key.lower()] = value


def _update_component_metadata(
    line: str,
    metadata: dict[str, Any],
) -> None:
    """Parse component and receiver area metadata."""
    component = re.search(r"Component:\s*([^,]+)", line, re.IGNORECASE)
    rx_area = re.search(r"Rxna:\s*([-\d.Ee+]+)", line, re.IGNORECASE)
    if component:
        metadata["component"] = component.group(1).strip()
    if rx_area:
        metadata["rx_area"] = float(rx_area.group(1))


def _looks_like_z_record(line: str) -> bool:
    """Return whether ``line`` begins like a contour data row."""
    return bool(re.match(r"^\d+\s+[-+]?\d", line))


def _parse_z_record(
    line: str,
    path: Path,
    line_no: int,
) -> TEMZRecord:
    """Parse one TEMAVG ``.Z`` contour row."""
    tail = _TAIL_RE.match(line)
    if tail is None:
        msg = f"Malformed TEMAVG .Z row in {path!s} at line {line_no}."
        raise ValueError(msg)
    values = _FLOAT_RE.findall(tail.group("prefix"))
    if len(values) < 4:
        msg = (
            f"Cannot parse TEMAVG .Z row in {path!s} "
            f"at line {line_no}."
        )
        raise ValueError(msg)
    return TEMZRecord(
        line=int(float(values[0])),
        station=float(values[1]),
        time=float(values[2]),
        magnitude=float(values[3]),
        frequency=float(tail.group("frequency")),
        window=int(tail.group("window")),
    )
