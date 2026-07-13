# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Reader for the MARE2DEM group-level RMS log file.

Port of ``m2d_read_group_rms_log.m``.

MARE2DEM writes a CSV-like log with a header line naming each data
group, followed by rows of numeric RMS values one row per iteration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["GroupRMSLog", "read_group_rms_log"]


@dataclass
class GroupRMSLog:
    """Contents of one MARE2DEM group-RMS log file.

    Attributes
    ----------
    path : pathlib.Path or None
        Source file.
    headers : list of str
        Column names from the first line of the file.
    rms_log : numpy.ndarray, shape (n_iterations, n_groups)
        RMS values, one row per iteration.
    """

    path: Path | None = None
    headers: list[str] = field(default_factory=list)
    rms_log: np.ndarray = field(
        default_factory=lambda: np.empty((0, 0), dtype=float)
    )

    @property
    def n_iterations(self) -> int:
        return len(self.rms_log)

    @property
    def n_groups(self) -> int:
        return len(self.headers)

    def __repr__(self) -> str:
        return (
            f"GroupRMSLog(n_iterations={self.n_iterations}, "
            f"n_groups={self.n_groups})"
        )


def read_group_rms_log(path: str | Path) -> GroupRMSLog:
    """Read a MARE2DEM group-level RMS log file.

    Port of ``m2d_read_group_rms_log.m``.

    Parameters
    ----------
    path : path-like
        CSV-like log file written by MARE2DEM.

    Returns
    -------
    GroupRMSLog
        Parsed log with headers and numeric RMS table.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.group_rms import read_group_rms_log
    >>> log = read_group_rms_log("mare2dem_group_rms.log")
    >>> log.n_iterations
    42
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    result = GroupRMSLog(path=path.resolve())
    text = path.read_text(errors="replace")
    lines = text.splitlines()

    if not lines:
        return result

    # First line is the comma-separated header
    result.headers = [h.strip() for h in lines[0].split(",")]

    # Remaining lines are numeric data
    remaining = "\n".join(lines[1:])
    n_cols = len(result.headers)
    if n_cols > 0 and remaining.strip():
        vals = np.fromstring(
            remaining.replace("\n", " "), sep=" ", dtype=float
        )
        n_rows = len(vals) // n_cols
        if n_rows > 0:
            result.rms_log = vals[: n_rows * n_cols].reshape(n_rows, n_cols)

    return result
