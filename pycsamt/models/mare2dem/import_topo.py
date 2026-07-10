# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Import topography and project it onto a MARE2DEM survey profile.

Port of ``m2d_importTopo.m``.

:func:`import_topo` loads a topography / bathymetry file, converts
coordinates if needed, checks that the profile orientation is consistent
with the survey, and returns a (y, z) table aligned with the MARE2DEM
2-D model coordinate system.

Input file formats supported
-----------------------------
The function auto-detects the column layout from the
:class:`TopoConfig` specification:

* **GeoMapApp-style** (4 columns: Lon, Lat, [Elev or Dist], [Dist or Elev]).
* **Distance-based** (2 columns: distance-along-profile in km, depth in m).
* **Distance in metres** (2 columns: distance in m, depth in m).

Coordinate system
-----------------
MARE2DEM uses a right-handed system with ``y`` along the profile and
``z`` positive downward (sea-depth positive, elevation negative).
The returned ``z`` values are depths in metres (positive = below
surface).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .geom.line_orientation import get_line_orientation
from .geom.utm import lonlat_to_utm

__all__ = ["TopoConfig", "TopoProfile", "import_topo"]


@dataclass
class TopoConfig:
    """Configuration for loading one topography file.

    Attributes
    ----------
    topo_file : str or path-like
        Path to the topography file (whitespace-separated columns).
    col_longitude : int or None
        1-based column index for longitude. Required when the file has
        geographic coordinates.
    col_latitude : int or None
        1-based column index for latitude.
    col_elevation_m : int or None
        1-based column for elevation in metres (positive up). Mutually
        exclusive with *col_depth_m*.
    col_depth_m : int or None
        1-based column for depth in metres (positive down).
    col_distance_km : int or None
        1-based column for along-profile distance in km.
    col_distance_m : int or None
        1-based column for along-profile distance in metres.
    """

    topo_file: str | Path = ""
    col_longitude: int | None = None
    col_latitude: int | None = None
    col_elevation_m: int | None = None
    col_depth_m: int | None = None
    col_distance_km: int | None = None
    col_distance_m: int | None = None


@dataclass
class TopoProfile:
    """Loaded and projected topography profile.

    Attributes
    ----------
    y_topo : numpy.ndarray
        Along-profile position in metres (MARE2DEM y axis).
    z_topo : numpy.ndarray
        Depth in metres, positive down (MARE2DEM z axis).
    northings : numpy.ndarray or None
        UTM northing (m) — present when converted from Lon/Lat.
    eastings : numpy.ndarray or None
        UTM easting (m) — present when converted from Lon/Lat.
    """

    y_topo: np.ndarray
    z_topo: np.ndarray
    northings: np.ndarray | None = None
    eastings: np.ndarray | None = None


def import_topo(
    cfg: TopoConfig,
    *,
    utm_north0: float = 0.0,
    utm_east0: float = 0.0,
    utm_theta: float = 0.0,
    utm_zone: int | None = None,
    south_hemi: bool = False,
    orientation_tol: float = 5.0,
) -> TopoProfile:
    """Import a topography file and project it onto the survey profile.

    Port of ``m2d_importTopo.m``.

    Parameters
    ----------
    cfg : TopoConfig
        Column layout and file path.
    utm_north0 : float
        Profile UTM origin northing (metres).
    utm_east0 : float
        Profile UTM origin easting (metres).
    utm_theta : float
        Profile UTM theta angle (degrees) — the ``stUTM.theta`` field
        from the ``.emdata`` UTM block.  Note: the survey-line direction
        is ``theta + 90°``.
    utm_zone : int or None
        UTM zone number.  Required only when Lon/Lat input is used.
    south_hemi : bool, default False
        Southern-hemisphere flag for UTM conversion.
    orientation_tol : float, default 5.0
        Maximum allowed angle (degrees) between the topography profile
        orientation and the survey line orientation.

    Returns
    -------
    TopoProfile
        Projected topography (y, z) in MARE2DEM coordinates.

    Raises
    ------
    FileNotFoundError
        When the topography file does not exist.
    ValueError
        When the topography and survey orientations differ by more than
        *orientation_tol* degrees.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.import_topo import TopoConfig, import_topo
    >>> cfg = TopoConfig(topo_file="topo.txt", col_distance_km=1, col_depth_m=2)
    >>> prof = import_topo(cfg, utm_north0=0.0, utm_east0=0.0)
    >>> prof.y_topo
    array([0., ...])
    """
    path = Path(cfg.topo_file)
    if not path.exists():
        raise FileNotFoundError(f"Topography file not found: {path}")

    data = np.loadtxt(str(path))
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # ---- extract depth ----
    if cfg.col_depth_m is not None:
        z_raw = data[:, cfg.col_depth_m - 1]
    elif cfg.col_elevation_m is not None:
        z_raw = -data[:, cfg.col_elevation_m - 1]
    else:
        raise ValueError(
            "TopoConfig must specify either col_depth_m or col_elevation_m."
        )

    # ---- extract / compute y position ----
    northings: np.ndarray | None = None
    eastings: np.ndarray | None = None

    if cfg.col_longitude is not None and cfg.col_latitude is not None:
        lon = data[:, cfg.col_longitude - 1]
        lat = data[:, cfg.col_latitude - 1]

        east, north, zone_used, _ = lonlat_to_utm(
            lon, lat,
            zone=utm_zone,
            south_hemi=south_hemi,
        )
        northings = np.asarray(north, dtype=float)
        eastings = np.asarray(east, dtype=float)

        # check orientation consistency
        topo_orientation = get_line_orientation(northings, eastings)
        topo_orientation = topo_orientation % 180
        line_orientation = (utm_theta + 90.0) % 180
        if abs(topo_orientation - line_orientation) > orientation_tol:
            raise ValueError(
                f"Topography orientation ({topo_orientation:.1f}°) and survey "
                f"line orientation ({line_orientation:.1f}°) differ by more "
                f"than {orientation_tol}°.  Check the topography file."
            )

        # project onto the survey profile
        theta = utm_theta
        cc = math.cos(math.radians(theta))
        ss = math.sin(math.radians(theta))
        dN = northings - utm_north0
        dE = eastings - utm_east0
        # x = cross-profile, y = along-profile
        y_raw = -ss * dN + cc * dE

    elif cfg.col_distance_km is not None:
        y_raw = data[:, cfg.col_distance_km - 1] * 1000.0

    elif cfg.col_distance_m is not None:
        y_raw = data[:, cfg.col_distance_m - 1]

    else:
        raise ValueError(
            "TopoConfig must specify geographic columns (col_longitude + "
            "col_latitude) or a distance column (col_distance_km or "
            "col_distance_m)."
        )

    return TopoProfile(
        y_topo=np.asarray(y_raw, dtype=float),
        z_topo=np.asarray(z_raw, dtype=float),
        northings=northings,
        eastings=eastings,
    )
