# Author: L. Kouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt.utils.io

I/O utility functions for PyCSAMT v2.0.
"""

import os
import re
from typing import Optional, Union

import numpy as np
from scipy.interpolate import interp1d

from ..exceptions import PycsamtError, ValidationError

__all__ = [
    "stn_separation",
    "parse_stn_profile",
]


def stn_separation(
    easting: Union[np.ndarray, list, tuple],
    northing: Union[np.ndarray, list, tuple],
    interpolate: bool = False,
) -> tuple[np.ndarray, float]:
    """
    Compute station separations (distance between successive points).

    Parameters
    ----------
    easting : array_like
        Easting coordinates (m).
    northing : array_like
        Northing coordinates (m).
    interpolate : bool, optional
        If True, returns an array of length N matching the number of
        electrodes by extrapolating the first separation. If False,
        returns length N-1 (number of dipoles). Default is False.

    Returns
    -------
    separations : np.ndarray
        Distances between successive points.
    mean_sep : float
        Mean separation value.
    """
    # Convert inputs to numpy arrays of float
    try:
        e = np.asarray(easting, dtype=float)
        n = np.asarray(northing, dtype=float)
    except Exception as err:
        raise ValidationError(f"Invalid coordinate arrays: {err}")

    if e.shape != n.shape:
        raise ValidationError(
            f"Easting and northing shapes differ: {e.shape} vs {n.shape}"
        )

    count = e.size
    if count < 2:
        return np.array([], dtype=float), 0.0

    # Calculate pairwise distances
    deltas = np.sqrt((np.diff(e)) ** 2 + (np.diff(n)) ** 2)

    if interpolate:
        # Extrapolate separations to length = count
        indices = np.arange(1, count)
        f = interp1d(indices, deltas, fill_value="extrapolate")
        separations = f(np.arange(count))
    else:
        separations = deltas

    mean_sep = float(np.mean(separations))
    return separations, mean_sep


def parse_stn_profile(
    file_path: str, delimiter: Optional[str] = None
) -> dict[str, np.ndarray]:
    """
    Parse a station profile file (.stn) containing columns such as
    station position (dot), easting, northing, elevation.

    Parameters
    ----------
    file_path : str
        Path to the station profile file.
    delimiter : str, optional
        Field delimiter. If None, whitespace splitting is used.

    Returns
    -------
    result : Dict[str, np.ndarray]
        Dictionary with keys:
        - 'position': station position (float)
        - 'easting': easting coordinate (float)
        - 'northing': northing coordinate (float)
        - 'elevation': elevation (float)

    Raises
    ------
    PycsamtError
        If file cannot be read or parsed.
    """
    if not os.path.isfile(file_path):
        raise PycsamtError(f"File not found: {file_path}")

    with open(file_path, encoding="utf8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    # Skip comment lines starting with '>' or '!'
    data_lines = [ln for ln in lines if not re.match(r"^[>!+]", ln)]
    if not data_lines:
        raise PycsamtError(f"No data found in {file_path}")

    # First non-comment line is header
    header = data_lines[0]
    if delimiter:
        cols = header.split(delimiter)
    else:
        cols = header.split()
    # Normalize column names
    cols = [col.strip().strip('"') for col in cols]

    # Read data rows
    data = []
    for ln in data_lines[1:]:
        parts = ln.split(delimiter) if delimiter else ln.split()
        if len(parts) != len(cols):
            raise PycsamtError(
                f"Line has {len(parts)} fields but expected {len(cols)}: {ln}"
            )
        data.append([float(p) for p in parts])

    arr = np.array(data)

    result = {}
    for idx, name in enumerate(cols):
        key = name.lower()
        if key in ("dot", "station", "stn"):
            result["position"] = arr[:, idx]
        elif "east" in key or key == "e":
            result["easting"] = arr[:, idx]
        elif "north" in key or key == "n":
            result["northing"] = arr[:, idx]
        elif "elev" in key or "h" == key:
            result["elevation"] = arr[:, idx]
        else:
            # preserve any other columns under raw names
            result[key] = arr[:, idx]

    # Ensure required keys
    for req in ("position", "easting", "northing", "elevation"):
        if req not in result:
            raise PycsamtError(
                f"Missing required column '{req}' in {file_path}"
            )

    return result
