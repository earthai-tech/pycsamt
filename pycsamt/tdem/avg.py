"""Readers for Zonge TEMAVG processed ``.AVG`` files."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Union

import numpy as np

from ..api.property import PyCSAMTObject
from ._base import TEMSounding

PathLike = Union[str, Path]

__all__ = [
    "TEMAVG",
    "TEMAVGRecord",
    "is_temavg_file",
]


_META_RE = re.compile(r"^\$\s*TEM:\s*([^=]+)=\s*(.+?)\s*$")
_VERSION_RE = re.compile(
    r"TEMAVG\s+(?P<version>[\d.]+).*?"
    r"Dated\s+(?P<dated>[^,]+),\s+Processed\s+(?P<processed>.+)"
)


@dataclass(frozen=True)
class TEMAVGRecord:
    """One processed TEMAVG gate value.

    Parameters
    ----------
    skip : int
        TEMAVG skip or block flag stored in the first column.
    tx : float
        Transmitter identifier from the ``Tx`` column.
    station : float
        Station coordinate or station number along the survey
        profile.
    frequency : float
        Repetition frequency in hertz.
    component : str
        Measured component label, for example ``"Hz"``.
    current : float
        Transmitter current in amperes from the ``Amps``
        column.
    window : int
        Time-window number.
    time : float
        Gate centre time in milliseconds, as written by
        TEMAVG.
    magnitude : float
        Processed transient magnitude. For TEMAVG contour
        files this is commonly reported in microvolts per
        ampere.
    ramp_app_res : float
        Ramp-corrected apparent resistivity in ohm metres.
    depth : float
        TEMAVG depth estimate in metres.
    percent_magnitude : float
        Percent magnitude or percent error column.
    """

    skip: int
    tx: float
    station: float
    frequency: float
    component: str
    current: float
    window: int
    time: float
    magnitude: float
    ramp_app_res: float
    depth: float
    percent_magnitude: float

    @property
    def time_s(self) -> float:
        """Gate centre time in seconds."""
        return self.time * 1e-3


@dataclass(repr=False)
class TEMAVG(PyCSAMTObject):
    """Processed content of one Zonge TEMAVG ``.AVG`` file.

    ``TEMAVG`` reads the human-readable output produced by the
    Zonge TEMAVG processing program. The object stores global
    metadata such as transmitter ramp, loop dimensions, and
    receiver area together with one row per station/time gate.
    It is a processed time-domain EM container, not an EDI or
    frequency-domain impedance object.

    Parameters
    ----------
    path : pathlib.Path
        Source file path.
    metadata : dict
        Parsed header metadata. Keys include entries such as
        ``"TXramp"``, ``"TXdx"``, ``"TXdy"``, ``"TXarea"``,
        and ``"RXarea"`` when present.
    records : list of TEMAVGRecord
        Parsed processed data rows.
    version : str, optional
        TEMAVG program version parsed from the first line.
    dated : str, optional
        Field data date parsed from the first line.
    processed : str, optional
        Processing date parsed from the first line.

    Examples
    --------
    >>> from pycsamt.tdem.avg import TEMAVG
    >>> avg = TEMAVG.read("data/TEMAVG/JIANGSU/TEM100.AVG")
    >>> avg.n_records > 0
    True
    >>> avg.stations[:3]
    [100.0, 120.0, 140.0]
    """

    path: Path
    metadata: dict[str, Any] = field(default_factory=dict)
    records: list[TEMAVGRecord] = field(default_factory=list)
    version: str | None = None
    dated: str | None = None
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
    ) -> TEMAVG:
        """Read a Zonge TEMAVG ``.AVG`` file.

        Parameters
        ----------
        path : path-like
            Processed TEMAVG file. The reader expects the
            standard Zonge header followed by the table whose
            columns include station, frequency, component,
            current, window, time, magnitude, apparent
            resistivity, depth, and percent magnitude.

        Returns
        -------
        TEMAVG
            Parsed TEMAVG file.
        """
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(p)

        metadata: dict[str, Any] = {}
        records: list[TEMAVGRecord] = []
        version = dated = processed = None
        in_table = False

        with p.open("r", errors="replace") as fh:
            for line_no, raw in enumerate(fh, 1):
                line = raw.strip()
                if not line:
                    continue
                if line_no == 1:
                    match = _VERSION_RE.search(line)
                    if match:
                        version = match.group("version")
                        dated = match.group("dated").strip()
                        processed = match.group("processed").strip()
                    metadata["title"] = line.lstrip("\\").strip()
                    continue
                meta = _META_RE.match(line)
                if meta:
                    key = meta.group(1).strip()
                    metadata[key] = _parse_metadata_value(meta.group(2))
                    continue
                if line.startswith("skp"):
                    in_table = True
                    continue
                if not in_table or line.startswith("\\-"):
                    continue
                if _looks_like_record(line):
                    records.append(_parse_avg_record(line, p, line_no))

        if not records:
            msg = f"No TEMAVG data rows found in {p!s}."
            raise ValueError(msg)
        return cls(
            path=p,
            metadata=metadata,
            records=records,
            version=version,
            dated=dated,
            processed=processed,
            verbose=verbose,
            logger=logger,
        )

    @property
    def n_records(self) -> int:
        """Number of processed gate rows."""
        return len(self.records)

    @property
    def stations(self) -> list[float]:
        """Sorted station values represented in the file."""
        return sorted({rec.station for rec in self.records})

    @property
    def windows(self) -> list[int]:
        """Sorted time-window numbers represented in the file."""
        return sorted({rec.window for rec in self.records})

    @property
    def tx_area(self) -> float | None:
        """Transmitter loop area in square metres when present."""
        value = self.metadata.get("TXarea")
        return float(value) if value is not None else None

    @property
    def rx_area(self) -> float | None:
        """Receiver loop or coil area in square metres when present."""
        value = self.metadata.get("RXarea")
        return float(value) if value is not None else None

    @property
    def tx_dx(self) -> float | None:
        """Transmitter loop x dimension in metres when present."""
        value = self.metadata.get("TXdx")
        return float(value) if value is not None else None

    @property
    def tx_dy(self) -> float | None:
        """Transmitter loop y dimension in metres when present."""
        value = self.metadata.get("TXdy")
        return float(value) if value is not None else None

    def rows_for_station(self, station: float) -> list[TEMAVGRecord]:
        """Return all gate rows for one station value."""
        return [rec for rec in self.records if rec.station == station]

    def to_soundings(
        self,
        *,
        stations: list[float] | None = None,
        component: str = "Hz",
        frequency: float | None = None,
        data_column: str = "magnitude",
        magnitude_unit: str = "uV/A",
        data_type: str = "voltage",
        rx_turns: int = 1,
        tx_turns: int = 1,
        coordinate_table: Any | None = None,
        profile: float | None = None,
        station_name_template: str = "{stem}_{station:g}",
        min_gates: int = 1,
        verbose: int = 0,
        logger: object | None = None,
    ) -> list[TEMSounding]:
        """Build one :class:`TEMSounding` per station.

        Parameters
        ----------
        stations : list of float, optional
            Station values to export. When omitted, every station
            in the file is converted.
        component : str, default "Hz"
            Component label to select from the TEMAVG rows.
        frequency : float, optional
            Repetition frequency to select. When omitted, all
            rows matching ``component`` are used.
        data_column : {"magnitude"}, default "magnitude"
            TEMAVG data column used for the sounding decay.
            The processed magnitude column is currently the only
            supported transient column.
        magnitude_unit : str, default "uV/A"
            Unit of the TEMAVG magnitude column. ``"uV/A"``
            means microvolts per ampere. ``"V/A"``, ``"uV"``,
            ``"V"``, and ``"SI"`` are also accepted.
        data_type : str, default "voltage"
            Output ``TEMSounding`` data type. ``"voltage"``
            keeps the decay as receiver voltage. ``"dBdt"``
            divides by receiver turns and area. ``"dHdt"``
            additionally divides by :math:`\\mu_0`.
        rx_turns, tx_turns : int, default 1
            Receiver and transmitter turn counts passed to the
            resulting soundings.
        coordinate_table : object, optional
            Coordinate table exposing ``get(profile, point)``.
            Matching coordinates are copied to the sounding
            ``x``, ``y``, and ``elevation`` fields.
        profile : float, optional
            Profile id used with ``coordinate_table``. If omitted
            and no coordinate table is supplied, coordinates stay
            at their default values.
        station_name_template : str, default "{stem}_{station:g}"
            Format string used to create each station name.
            Available fields are ``stem``, ``station``,
            ``profile``, and ``component``.
        min_gates : int, default 1
            Minimum number of selected time gates required to
            create a sounding.

        Returns
        -------
        list of TEMSounding
            Station-wise soundings ready for late-time or Fourier
            conversion.
        """
        if data_column != "magnitude":
            msg = "Only data_column='magnitude' is currently supported."
            raise ValueError(msg)

        selected_stations = self.stations if stations is None else stations
        soundings: list[TEMSounding] = []
        for station in selected_stations:
            rows = [
                rec for rec in self.rows_for_station(station)
                if rec.component == component
            ]
            if frequency is not None:
                rows = [rec for rec in rows if rec.frequency == frequency]
            rows = sorted(rows, key=lambda rec: rec.time_s)
            if len(rows) < min_gates:
                continue

            current = float(rows[0].current)
            time_gates = np.array([rec.time_s for rec in rows], dtype=float)
            values = np.array([rec.magnitude for rec in rows], dtype=float)
            voltage = _scale_temavg_magnitude(
                values,
                current=current,
                magnitude_unit=magnitude_unit,
            )
            rx_area = self.rx_area or 1.0
            data = _temavg_data_for_type(
                voltage,
                current=current,
                tx_area=self._required_tx_area(),
                rx_area=rx_area,
                rx_turns=rx_turns,
                data_type=data_type,
            )
            errors = _temavg_percent_errors(data, rows)
            x = y = elevation = 0.0
            coord = None
            if coordinate_table is not None and profile is not None:
                coord = coordinate_table.get(profile, station)
            if coord is not None:
                x = coord.x
                y = coord.y
                elevation = coord.elevation

            name = station_name_template.format(
                stem=self.path.stem,
                station=float(station),
                profile=profile if profile is not None else "",
                component=component,
            )
            soundings.append(
                TEMSounding.from_arrays(
                    time_gates,
                    data,
                    current=current,
                    tx_area=self._required_tx_area(),
                    data_type=data_type,
                    tx_turns=tx_turns,
                    rx_area=rx_area,
                    rx_turns=rx_turns,
                    loop_shape=self._loop_shape(),
                    loop_dims=self._loop_dims(),
                    station_name=name,
                    x=x,
                    y=y,
                    elevation=elevation,
                    error=errors,
                    verbose=verbose,
                    logger=logger,
                )
            )
        return soundings

    def to_records(self) -> list[dict[str, Any]]:
        """Return records as dictionaries with file metadata."""
        rows: list[dict[str, Any]] = []
        for rec in self.records:
            row = {
                "source_file": self.path.name,
                "station": rec.station,
                "tx": rec.tx,
                "frequency": rec.frequency,
                "component": rec.component,
                "current": rec.current,
                "window": rec.window,
                "time_ms": rec.time,
                "time_s": rec.time_s,
                "magnitude": rec.magnitude,
                "ramp_app_res": rec.ramp_app_res,
                "depth": rec.depth,
                "percent_magnitude": rec.percent_magnitude,
                "tx_ramp": self.metadata.get("TXramp"),
                "tx_area": self.metadata.get("TXarea"),
                "rx_area": self.metadata.get("RXarea"),
            }
            rows.append(row)
        return rows

    def to_dataframe(self):
        """Return the parsed table as a :class:`pandas.DataFrame`."""
        try:
            import pandas as pd
        except ImportError as exc:
            msg = "TEMAVG.to_dataframe requires pandas."
            raise ImportError(msg) from exc
        return pd.DataFrame(self.to_records())

    def _required_tx_area(self) -> float:
        """Return transmitter area or raise a clear error."""
        if self.tx_area is None:
            msg = f"TEMAVG file {self.path!s} does not define TXarea."
            raise ValueError(msg)
        return self.tx_area

    def _loop_shape(self) -> str:
        """Infer transmitter loop shape from TEMAVG dimensions."""
        if self.tx_dx is None or self.tx_dy is None:
            return "square"
        if np.isclose(self.tx_dx, self.tx_dy):
            return "square"
        return "rectangular"

    def _loop_dims(self) -> tuple[float, ...]:
        """Infer transmitter loop dimensions from TEMAVG metadata."""
        if self.tx_dx is None or self.tx_dy is None:
            side = np.sqrt(self._required_tx_area())
            return (float(side),)
        if np.isclose(self.tx_dx, self.tx_dy):
            return (float(self.tx_dx),)
        return (float(self.tx_dx), float(self.tx_dy))


def is_temavg_file(path: PathLike) -> bool:
    """Return ``True`` when ``path`` looks like a TEMAVG file."""
    p = Path(path)
    if not p.exists() or not p.is_file():
        return False
    try:
        with p.open("r", errors="replace") as fh:
            head = "".join(fh.readline() for _ in range(8))
    except OSError:
        return False
    return "TEMAVG" in head and "TXarea" in head


def _parse_metadata_value(text: str) -> Any:
    """Parse a TEMAVG metadata value while preserving units."""
    value = text.strip()
    number = re.match(r"^([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)", value)
    if number:
        parsed = float(number.group(1))
        if parsed.is_integer():
            return int(parsed)
        return parsed
    return value


def _looks_like_record(line: str) -> bool:
    """Return whether ``line`` starts like a TEMAVG data row."""
    return bool(re.match(r"^\d+\s+[-+]?\d", line))


def _parse_avg_record(
    line: str,
    path: Path,
    line_no: int,
) -> TEMAVGRecord:
    """Parse one whitespace-delimited TEMAVG data row."""
    parts = line.split()
    if len(parts) < 12:
        msg = f"Malformed TEMAVG row in {path!s} at line {line_no}."
        raise ValueError(msg)
    try:
        return TEMAVGRecord(
            skip=int(parts[0]),
            tx=float(parts[1]),
            station=float(parts[2]),
            frequency=float(parts[3]),
            component=parts[4],
            current=float(parts[5]),
            window=int(parts[6]),
            time=float(parts[7]),
            magnitude=float(parts[8]),
            ramp_app_res=float(parts[9]),
            depth=float(parts[10]),
            percent_magnitude=_parse_float(parts[11]),
        )
    except ValueError as exc:
        msg = f"Cannot parse TEMAVG row in {path!s} at line {line_no}."
        raise ValueError(msg) from exc


def _parse_float(value: str) -> float:
    """Parse numeric TEMAVG tokens, mapping quality marks to nan."""
    if value.strip() in {"*", "-"}:
        return float("nan")
    return float(value)


def _scale_temavg_magnitude(
    values: np.ndarray,
    *,
    current: float,
    magnitude_unit: str,
) -> np.ndarray:
    """Scale TEMAVG magnitudes to the units expected downstream."""
    unit = magnitude_unit.strip().lower()
    if unit in {"uv/a", "microvolt/a", "microvolts/a"}:
        return values * current * 1e-6
    if unit in {"v/a", "volt/a", "volts/a"}:
        return values * current
    if unit in {"uv", "microvolt", "microvolts"}:
        return values * 1e-6
    if unit in {"v", "volt", "volts", "si"}:
        return values.copy()
    msg = (
        "magnitude_unit must be one of 'uV/A', 'V/A', "
        "'uV', 'V', or 'SI'."
    )
    raise ValueError(msg)


def _temavg_percent_errors(
    data: np.ndarray,
    rows: list[TEMAVGRecord],
) -> np.ndarray | None:
    """Build absolute errors from the TEMAVG percent column."""
    percent = np.array(
        [rec.percent_magnitude for rec in rows],
        dtype=float,
    )
    if not np.isfinite(percent).any() or np.nanmax(np.abs(percent)) == 0.0:
        return None
    return np.abs(data) * percent / 100.0


def _temavg_data_for_type(
    voltage: np.ndarray,
    *,
    current: float,
    tx_area: float,
    rx_area: float,
    rx_turns: int,
    data_type: str,
) -> np.ndarray:
    """Convert receiver voltage to a TEMSounding data convention."""
    if data_type == "voltage":
        return voltage
    if data_type == "dBdt":
        return voltage / (float(rx_turns) * float(rx_area))
    if data_type == "dHdt":
        mu0 = 4.0 * np.pi * 1e-7
        return voltage / (float(rx_turns) * float(rx_area) * mu0)
    if data_type == "normalized_voltage":
        return voltage / (float(current) * float(tx_area))
    msg = (
        "data_type must be one of 'voltage', 'dBdt', "
        "'dHdt', or 'normalized_voltage'."
    )
    raise ValueError(msg)
