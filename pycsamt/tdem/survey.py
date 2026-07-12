"""Folder-level readers for time-domain EM surveys."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Union

from ..api.property import PyCSAMTObject
from ._base import TEMSounding
from .avg import TEMAVG
from .coordinates import (
    TEMCoordinate,
    TEMCoordinateTable,
    read_tem_coordinates,
)
from .log import TEMLog
from .zplot import TEMZPlot

PathLike = Union[str, Path]

__all__ = [
    "TEMSurvey",
    "read_temavg_survey",
]


@dataclass(repr=False)
class TEMSurvey(PyCSAMTObject):
    """Collection of processed TEM files from one survey folder.

    Parameters
    ----------
    root : pathlib.Path
        Directory scanned for time-domain EM files.
    avg_files : dict of str to TEMAVG
        Parsed ``.AVG`` files keyed by file stem.
    z_files : dict of str to TEMZPlot
        Parsed ``.Z`` contour files keyed by file stem.
    log_files : dict of str to TEMLog
        Parsed ``.LOG`` processing files keyed by file stem.
    coordinates : TEMCoordinateTable, optional
        Profile/point coordinate table used to enrich exported
        survey records.
    companion_files : dict
        Available companion files grouped by stem. Each value
        may contain ``"log"``, ``"z"``, and other format keys
        that can be parsed in later processing stages.

    Examples
    --------
    >>> from pycsamt.tdem.survey import read_temavg_survey
    >>> survey = read_temavg_survey("data/TEMAVG/JIANGSU")
    >>> survey.n_avg_files > 0
    True
    """

    root: Path
    avg_files: dict[str, TEMAVG] = field(default_factory=dict)
    z_files: dict[str, TEMZPlot] = field(default_factory=dict)
    log_files: dict[str, TEMLog] = field(default_factory=dict)
    coordinates: TEMCoordinateTable | None = None
    companion_files: dict[str, dict[str, Path]] = field(
        default_factory=dict,
    )
    verbose: int = 0
    logger: object | None = None

    @property
    def n_avg_files(self) -> int:
        """Number of parsed ``.AVG`` files."""
        return len(self.avg_files)

    @property
    def n_z_files(self) -> int:
        """Number of parsed ``.Z`` contour files."""
        return len(self.z_files)

    @property
    def n_log_files(self) -> int:
        """Number of parsed ``.LOG`` processing files."""
        return len(self.log_files)

    @property
    def stems(self) -> list[str]:
        """Sorted file stems represented by parsed ``.AVG`` files."""
        return sorted(self.avg_files)

    @property
    def stations(self) -> list[float]:
        """Sorted union of all station values in parsed data."""
        values: set[float] = set()
        for avg in self.avg_files.values():
            values.update(avg.stations)
        return sorted(values)

    def to_records(self) -> list[dict[str, Any]]:
        """Return all parsed ``.AVG`` rows as dictionaries."""
        rows: list[dict[str, Any]] = []
        for stem in self.stems:
            profile = _profile_from_stem(stem)
            for row in self.avg_files[stem].to_records():
                row["profile"] = profile
                coord = self.coordinate_for(profile, row["station"])
                if coord is not None:
                    row.update(_coordinate_record(coord))
                rows.append(row)
        return rows

    def to_dataframe(self):
        """Return all parsed rows as a :class:`pandas.DataFrame`."""
        try:
            import pandas as pd
        except ImportError as exc:
            msg = "TEMSurvey.to_dataframe requires pandas."
            raise ImportError(msg) from exc
        return pd.DataFrame(self.to_records())

    def get(self, stem: str) -> TEMAVG:
        """Return one parsed ``.AVG`` file by stem."""
        return self.avg_files[stem]

    def get_z(self, stem: str) -> TEMZPlot:
        """Return one parsed ``.Z`` file by stem."""
        return self.z_files[stem]

    def get_log(self, stem: str) -> TEMLog:
        """Return one parsed ``.LOG`` file by stem."""
        return self.log_files[stem]

    def coordinate_for(
        self,
        profile: float,
        point: float,
    ) -> TEMCoordinate | None:
        """Return coordinate metadata for a profile/station pair."""
        if self.coordinates is None:
            return None
        return self.coordinates.get(profile, point)

    def to_soundings(
        self,
        *,
        stems: list[str] | None = None,
        component: str = "Hz",
        frequency: float | None = None,
        data_column: str = "magnitude",
        magnitude_unit: str = "uV/A",
        data_type: str = "voltage",
        rx_turns: int = 1,
        tx_turns: int = 1,
        min_gates: int = 1,
        verbose: int | None = None,
        logger: object | None = None,
    ) -> list[TEMSounding]:
        """Build station-wise :class:`TEMSounding` objects.

        Parameters
        ----------
        stems : list of str, optional
            File stems to export. When omitted, all parsed
            ``.AVG`` files are converted in sorted order.
        component : str, default "Hz"
            TEMAVG component to select.
        frequency : float, optional
            Repetition frequency to select. When omitted, all
            rows for ``component`` are used.
        data_column : {"magnitude"}, default "magnitude"
            TEMAVG transient column used for the sounding decay.
        magnitude_unit : str, default "uV/A"
            Unit of the TEMAVG magnitude column before conversion
            to the requested ``data_type`` representation.
        data_type : str, default "voltage"
            Data type stored in each sounding. The default is
            receiver voltage in volts.
        rx_turns, tx_turns : int, default 1
            Receiver and transmitter turn counts passed to each
            sounding.
        min_gates : int, default 1
            Minimum number of selected time gates required for a
            station to be exported.

        Returns
        -------
        list of TEMSounding
            Soundings suitable for :class:`LateTimeTransform` or
            :class:`TEMtoEDI`.
        """
        selected = self.stems if stems is None else stems
        soundings: list[TEMSounding] = []
        for stem in selected:
            avg = self.avg_files[stem]
            profile = _profile_from_stem(stem)
            soundings.extend(
                avg.to_soundings(
                    component=component,
                    frequency=frequency,
                    data_column=data_column,
                    magnitude_unit=magnitude_unit,
                    data_type=data_type,
                    rx_turns=rx_turns,
                    tx_turns=tx_turns,
                    coordinate_table=self.coordinates,
                    profile=profile,
                    min_gates=min_gates,
                    verbose=self.verbose if verbose is None else verbose,
                    logger=self.logger if logger is None else logger,
                )
            )
        return soundings


def read_temavg_survey(
    path: PathLike,
    *,
    pattern: str = "*.AVG",
    coordinate_file: PathLike | None = None,
    verbose: int = 0,
    logger: object | None = None,
) -> TEMSurvey:
    """Read a directory of Zonge TEMAVG processed files.

    Parameters
    ----------
    path : path-like
        Directory containing TEMAVG ``.AVG`` files and optional
        companion files such as ``.LOG``, ``.Z``, and ``.mde``.
    pattern : str, default "*.AVG"
        Glob pattern used to select processed average files.
    coordinate_file : path-like, optional
        Profile/point coordinate file. If omitted, the reader
        looks for common coordinate-table names in ``path``.

    Returns
    -------
    TEMSurvey
        Parsed survey collection.

    Examples
    --------
    >>> from pycsamt.tdem import read_temavg_survey
    >>> survey = read_temavg_survey("data/TEMAVG/JIANGSU")
    >>> len(survey.stations) > 0
    True
    """
    root = Path(path)
    if not root.exists():
        raise FileNotFoundError(root)
    if not root.is_dir():
        msg = f"TEMAVG survey path is not a directory: {root!s}"
        raise NotADirectoryError(msg)

    avg_files: dict[str, TEMAVG] = {}
    for avg_path in sorted(root.glob(pattern)):
        avg_files[avg_path.stem] = TEMAVG.read(
            avg_path,
            verbose=verbose,
            logger=logger,
        )

    if not avg_files:
        msg = f"No TEMAVG files matching {pattern!r} found in {root!s}."
        raise ValueError(msg)

    companions = _group_companion_files(root)
    z_files: dict[str, TEMZPlot] = {}
    log_files: dict[str, TEMLog] = {}
    for stem, files in companions.items():
        z_path = files.get("z")
        if z_path is not None:
            z_files[stem] = TEMZPlot.read(
                z_path,
                verbose=verbose,
                logger=logger,
            )
        log_path = files.get("log")
        if log_path is not None:
            log_files[stem] = TEMLog.read(
                log_path,
                verbose=verbose,
                logger=logger,
            )

    coordinates = _read_coordinates_if_available(
        root,
        coordinate_file,
        verbose=verbose,
        logger=logger,
    )

    return TEMSurvey(
        root=root,
        avg_files=avg_files,
        z_files=z_files,
        log_files=log_files,
        companion_files=companions,
        coordinates=coordinates,
        verbose=verbose,
        logger=logger,
    )


def _group_companion_files(root: Path) -> dict[str, dict[str, Path]]:
    """Group recognized companion files by stem."""
    grouped: dict[str, dict[str, Path]] = {}
    suffix_map = {
        ".log": "log",
        ".z": "z",
        ".mde": "mde",
    }
    for p in root.iterdir():
        key = suffix_map.get(p.suffix.lower())
        if key is None:
            continue
        grouped.setdefault(p.stem, {})[key] = p
    return grouped


def _read_coordinates_if_available(
    root: Path,
    coordinate_file: PathLike | None,
    *,
    verbose: int = 0,
    logger: object | None = None,
) -> TEMCoordinateTable | None:
    """Read a coordinate table when one can be found."""
    if coordinate_file is not None:
        coord_path = Path(coordinate_file)
        if not coord_path.is_absolute():
            coord_path = root / coord_path
        return read_tem_coordinates(
            coord_path,
            verbose=verbose,
            logger=logger,
        )

    candidates = [
        "Coordinate of measuring point.xls",
        "Coordinate of measuring point.xlsx",
        "coordinates.csv",
        "Coordinates.csv",
        "coordinate.csv",
    ]
    for name in candidates:
        candidate = root / name
        if candidate.exists():
            try:
                return read_tem_coordinates(
                    candidate,
                    verbose=verbose,
                    logger=logger,
                )
            except (ImportError, OSError, RuntimeError, ValueError):
                return None
    return None


def _profile_from_stem(stem: str) -> float:
    """Extract a numeric profile id from a TEM file stem."""
    digits = "".join(ch for ch in stem if ch.isdigit())
    return float(digits) if digits else float("nan")


def _coordinate_record(coord: TEMCoordinate) -> dict[str, Any]:
    """Return coordinate fields for a survey row."""
    return {
        "coord_profile": coord.profile,
        "coord_point": coord.point,
        "gauss_x": coord.gauss_x,
        "gauss_y": coord.gauss_y,
        "x": coord.x,
        "y": coord.y,
        "elevation": coord.elevation,
        "remark": coord.remark,
    }
