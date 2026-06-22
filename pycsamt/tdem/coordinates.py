"""Coordinate readers for TEM survey station tables."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from ..api.property import PyCSAMTObject

PathLike = str | Path

__all__ = [
    "TEMCoordinate",
    "TEMCoordinateTable",
    "read_tem_coordinates",
]


@dataclass(frozen=True)
class TEMCoordinate:
    """Coordinate metadata for one TEM station point.

    Parameters
    ----------
    profile : float
        Survey profile or line identifier.
    point : float
        Station point along the profile.
    gauss_x, gauss_y : float
        Projected Gauss coordinates in metres.
    x, y : float
        Relative local coordinates in metres.
    elevation : float
        Station elevation in metres.
    remark : str, optional
        Free-text field note from the coordinate table.
    """

    profile: float
    point: float
    gauss_x: float
    gauss_y: float
    x: float
    y: float
    elevation: float
    remark: str = ""


@dataclass(repr=False)
class TEMCoordinateTable(PyCSAMTObject):
    """Collection of TEM station coordinates.

    Parameters
    ----------
    path : pathlib.Path
        Source coordinate file.
    coordinates : dict
        Coordinate records keyed by ``(profile, point)``.
    """

    path: Path
    coordinates: dict[tuple[float, float], TEMCoordinate] = field(
        default_factory=dict,
    )
    verbose: int = 0
    logger: object | None = None

    @classmethod
    def read(
        cls,
        path: PathLike,
        *,
        verbose: int = 0,
        logger: object | None = None,
    ) -> TEMCoordinateTable:
        """Read a TEM coordinate table.

        Parameters
        ----------
        path : path-like
            Coordinate file. CSV and TSV files are read
            directly. Legacy ``.xls`` files are read with
            :func:`pandas.read_excel` when ``xlrd`` is
            installed, otherwise LibreOffice is used as a
            conversion fallback when available.

        Returns
        -------
        TEMCoordinateTable
            Parsed coordinate table.
        """
        return read_tem_coordinates(path, verbose=verbose, logger=logger)

    @property
    def n_points(self) -> int:
        """Number of coordinate records."""
        return len(self.coordinates)

    @property
    def profiles(self) -> list[float]:
        """Sorted profile identifiers."""
        return sorted({key[0] for key in self.coordinates})

    @property
    def points(self) -> list[float]:
        """Sorted station point identifiers."""
        return sorted({key[1] for key in self.coordinates})

    def get(
        self,
        profile: float,
        point: float,
    ) -> TEMCoordinate | None:
        """Return coordinate metadata for ``profile`` and ``point``."""
        return self.coordinates.get((float(profile), float(point)))

    def to_records(self) -> list[dict[str, Any]]:
        """Return coordinate records as dictionaries."""
        rows: list[dict[str, Any]] = []
        for coord in self.coordinates.values():
            rows.append(
                {
                    "profile": coord.profile,
                    "point": coord.point,
                    "gauss_x": coord.gauss_x,
                    "gauss_y": coord.gauss_y,
                    "x": coord.x,
                    "y": coord.y,
                    "elevation": coord.elevation,
                    "remark": coord.remark,
                }
            )
        return rows

    def to_dataframe(self):
        """Return coordinates as a :class:`pandas.DataFrame`."""
        try:
            import pandas as pd
        except ImportError as exc:
            msg = "TEMCoordinateTable.to_dataframe requires pandas."
            raise ImportError(msg) from exc
        return pd.DataFrame(self.to_records())


def read_tem_coordinates(
    path: PathLike,
    *,
    verbose: int = 0,
    logger: object | None = None,
) -> TEMCoordinateTable:
    """Read a TEM profile/point coordinate table.

    Parameters
    ----------
    path : path-like
        Coordinate file in CSV, TSV, XLS, or XLSX format.

    Returns
    -------
    TEMCoordinateTable
        Parsed coordinates keyed by ``(profile, point)``.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    suffix = p.suffix.lower()
    if suffix in {".csv", ".txt"}:
        rows = _read_text_rows(p, delimiter=",")
    elif suffix in {".tsv"}:
        rows = _read_text_rows(p, delimiter="\t")
    elif suffix in {".xls", ".xlsx"}:
        rows = _read_excel_rows(p)
    else:
        msg = f"Unsupported TEM coordinate file suffix: {suffix}"
        raise ValueError(msg)
    coords = _parse_coordinate_rows(rows)
    return TEMCoordinateTable(
        path=p,
        coordinates=coords,
        verbose=verbose,
        logger=logger,
    )


def _read_text_rows(path: Path, delimiter: str) -> list[list[Any]]:
    """Read rows from a delimited text coordinate file."""
    with path.open(newline="", errors="replace") as fh:
        return [row for row in csv.reader(fh, delimiter=delimiter)]


def _read_excel_rows(path: Path) -> list[list[Any]]:
    """Read Excel rows with pandas or LibreOffice fallback."""
    try:
        import pandas as pd

        frame = pd.read_excel(path, header=None)
        return frame.fillna("").values.tolist()
    except ImportError:
        pass
    except Exception:
        if path.suffix.lower() == ".xlsx":
            raise

    return _convert_xls_to_csv_rows(path)


def _convert_xls_to_csv_rows(path: Path) -> list[list[Any]]:
    """Convert a legacy Excel file to CSV using LibreOffice."""
    if shutil.which("libreoffice") is None:
        msg = (
            "Reading legacy .xls coordinate files requires xlrd "
            "or a LibreOffice executable on PATH."
        )
        raise ImportError(msg)
    with TemporaryDirectory() as tmp:
        outdir = Path(tmp)
        lo_profile = outdir / "lo-profile"
        cmd = [
            "libreoffice",
            "--headless",
            f"-env:UserInstallation=file://{lo_profile}",
            "--convert-to",
            "csv",
            "--outdir",
            str(outdir),
            str(path),
        ]
        env = os.environ.copy()
        env.update(
            {
                "HOME": tmp,
                "XDG_CACHE_HOME": str(outdir / "cache"),
                "XDG_CONFIG_HOME": str(outdir / "config"),
                "XDG_RUNTIME_DIR": str(outdir / "runtime"),
            }
        )
        runtime = Path(env["XDG_RUNTIME_DIR"])
        runtime.mkdir(mode=0o700, parents=True, exist_ok=True)
        try:
            subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                env=env,
            )
        except subprocess.CalledProcessError as exc:
            detail = (exc.stderr or exc.stdout or "").strip()
            msg = f"LibreOffice failed to convert {path!s}."
            if detail:
                msg = f"{msg} {detail}"
            raise RuntimeError(msg) from exc
        csv_path = outdir / f"{path.stem}.csv"
        if not csv_path.exists():
            matches = list(outdir.glob("*.csv"))
            if not matches:
                msg = f"LibreOffice did not create a CSV for {path!s}."
                raise FileNotFoundError(msg)
            csv_path = matches[0]
        rows = _read_text_rows(csv_path, delimiter=",")
    return rows


def _parse_coordinate_rows(
    rows: list[list[Any]],
) -> dict[tuple[float, float], TEMCoordinate]:
    """Parse rows from the known TEM coordinate table layout."""
    coords: dict[tuple[float, float], TEMCoordinate] = {}
    for row in rows:
        if len(row) < 7:
            continue
        try:
            profile = float(row[0])
            point = float(row[1])
            gauss_x = float(row[2])
            gauss_y = float(row[3])
            x = float(row[4])
            y = float(row[5])
            elevation = float(row[6])
        except (TypeError, ValueError):
            continue
        remark = str(row[7]).strip() if len(row) > 7 else ""
        key = (profile, point)
        coords[key] = TEMCoordinate(
            profile=profile,
            point=point,
            gauss_x=gauss_x,
            gauss_y=gauss_y,
            x=x,
            y=y,
            elevation=elevation,
            remark=remark,
        )
    if not coords:
        msg = "No coordinate rows found in TEM coordinate table."
        raise ValueError(msg)
    return coords
