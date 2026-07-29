# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Survey profile line-orientation estimation.

Port of ``m2d_getLineOrientation.m``.

Fits a linear regression to a set of UTM northing/easting coordinates
and returns the orientation of the best-fit line in degrees clockwise
from geographic north (0–180°).
"""

from __future__ import annotations

import math

import numpy as np

__all__ = ["get_line_orientation", "project_onto_line"]


def get_line_orientation(
    northings: np.ndarray,
    eastings: np.ndarray,
) -> float:
    """Estimate the survey profile orientation from station UTM positions.

    Port of ``m2d_getLineOrientation.m``.

    The function fits a line to the input (northing, easting) pairs and
    returns its bearing.  The result is in the range 0–180°:

    * 0° / 180° → N–S profile
    * 90° → E–W profile

    Parameters
    ----------
    northings : array-like, shape (n,)
        UTM northing coordinates in metres.
    eastings : array-like, shape (n,)
        UTM easting coordinates in metres.

    Returns
    -------
    float
        Line orientation in degrees clockwise from geographic north
        (0 ≤ result ≤ 180).

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.models.mare2dem.geom.line_orientation import (
    ...     get_line_orientation,
    ... )
    >>> northings = np.array([0.0, 1000.0, 2000.0])
    >>> eastings = np.array([0.0, 0.0, 0.0])
    >>> get_line_orientation(northings, eastings)  # N-S profile
    0.0
    """
    n = np.asarray(northings, dtype=float).ravel()
    e = np.asarray(eastings, dtype=float).ravel()

    dN = n.max() - n.min()
    dE = e.max() - e.min()

    if dE > dN:
        # more E-W variation: fit N = a*E + b
        coeffs = np.polyfit(e, n, 1)
        orientation = math.atan(coeffs[0]) * 180.0 / math.pi
        orientation = 90.0 - orientation
    else:
        # more N-S variation: fit E = a*N + b
        coeffs = np.polyfit(n, e, 1)
        orientation = math.atan(coeffs[0]) * 180.0 / math.pi
        if orientation < 0:
            orientation += 180.0

    return float(orientation)


def project_onto_line(
    northings: np.ndarray,
    eastings: np.ndarray,
    utm0_north: float,
    utm0_east: float,
    line_orientation: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Project UTM (northing, easting) positions onto a survey profile.

    Parameters
    ----------
    northings : array-like
        Station UTM northings in metres.
    eastings : array-like
        Station UTM eastings in metres.
    utm0_north : float
        Profile origin northing in metres.
    utm0_east : float
        Profile origin easting in metres.
    line_orientation : float
        Profile orientation in degrees clockwise from north.

    Returns
    -------
    x : numpy.ndarray
        Cross-profile (x) offset in metres.
    y : numpy.ndarray
        Along-profile (y) offset in metres.
    """
    n = np.asarray(northings, dtype=float)
    e = np.asarray(eastings, dtype=float)

    theta = line_orientation - 90.0
    cc = math.cos(math.radians(theta))
    ss = math.sin(math.radians(theta))

    dN = n - utm0_north
    dE = e - utm0_east

    x = cc * dN + ss * dE
    y = -ss * dN + cc * dE
    return x, y
