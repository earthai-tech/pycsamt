# -*- coding: utf-8 -*-
# Author: Kouadio Laurent alias Daniel <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.gis.config

Configuration loader for the GIS subpackage.

- Verifies availability of GDAL or a PROJ-based fallback
  (``pyproj``).
- Loads EPSG definitions into :data:`EPSG_DICT` for use by
  lightweight transformers.

On import, the module attempts to locate GDAL resources.
If unavailable, it verifies that ``pyproj`` can be imported.
EPSG definitions are then loaded from ``pyproj`` data where
possible; otherwise a local ``epsg.npy`` mapping is used.

Exports
-------
HAS_GDAL : bool
    ``True`` when GDAL/OSR is available.
EPSG_DICT : dict[int, str]
    Mapping of EPSG codes to PROJ strings.
GDALMissingError : Exception
    Raised when a GDAL-dependent routine is used without GDAL.
GisError : Exception
    Base exception for GIS utilities.

Notes
-----
- Modern ``pyproj`` (>=3) ships a SQLite ``proj.db`` rather
  than a flat ``epsg`` file. This loader uses a best-effort
  approach to find legacy definitions; if none are found, the
  local ``epsg.npy`` is used as a fallback.
"""

from __future__ import annotations

import os
import re
from typing import Dict

import numpy as np

from pycsamt.decorators import GdalDataCheck
from pycsamt.log.logger import get_logger

logger = get_logger(__name__)


# Availability checks

_HAS_GDAL = GdalDataCheck._locate_gdal_data()
if not _HAS_GDAL:
    try:
        import pyproj  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "Either GDAL or pyproj must be installed for GIS "
            "features."
        ) from exc


# Build EPSG dictionary
EPSG_DICT: Dict[int, str] = {}

try:
    import pyproj  # noqa: F401

    # pyproj 3+ uses proj.db; older releases may ship 'epsg'
    data_dir = getattr(pyproj, "pyproj_datadir", None)
    epsg_path = (
        os.path.join(data_dir, "epsg") if data_dir else None
    )

    if epsg_path and os.path.isfile(epsg_path):
        with open(epsg_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip().startswith("#"):
                    continue
                codes = re.findall(r"<(\d+)>", line)
                if not codes:
                    continue
                code = int(codes[0])
                m = re.search(r">([^<]+)<", line)
                if m:
                    EPSG_DICT[code] = m.group(1).strip()
    else:
        raise FileNotFoundError

except Exception:
    # Fallback: load local numpy file of EPSG mappings
    here = os.path.dirname(os.path.abspath(__file__))
    epsg_file = os.path.join(here, "epsg.npy")
    try:
        EPSG_DICT = np.load(
            epsg_file,
            allow_pickle=True,
        ).item()
    except Exception as e:  # pragma: no cover
        logger.error(
            "Failed to load EPSG definitions: %s",
            str(e),
        )
        EPSG_DICT = {}

# Expose availability
HAS_GDAL = _HAS_GDAL


class GDALMissingError(RuntimeError):
    r"""
    Raised when a GDAL/OSR-dependent function is called while
    GDAL is not available.

    Notes
    -----
    Install GDAL (and the ``osgeo`` Python bindings) to use
    GDAL-backed transforms. Prefer high-level helpers that can
    fall back to PROJ, such as
    :func:`project_point_ll2utm` and
    :func:`project_point_utm2ll`, when GDAL is absent.
    """


class GisError(Exception):
    r"""
    Base exception for GIS utilities.
    """

# Create a configuration class for elevation APIs
class ElevationAPIConfig:
    """Configuration for elevation API services."""
    
    # Available elevation API endpoints
    APIS = {
        'open_meteo': {
            'url': 'https://api.open-meteo.com/v1/elevation',
            'params_format': 'comma_separated',
            'response_key': 'elevation'
        },
        'open_topo_data': {
            'url': 'https://api.opentopodata.org/v1/elevation',
            'params_format': 'pipe_separated',
            'response_key': 'results.elevation'
        },
        'usgs_ned': {
            'url': 'https://nationalmap.gov/epqs/pqs.php',
            'params_format': 'individual',
            'response_key': (
                'USGS_Elevation_Point_Query_Service.Elevation_Query'
            )
        }
    }
    
    # Default API to use
    DEFAULT_API = 'open_meteo'
    
    @classmethod
    def get_api_config(cls, api_name=None):
        """Get configuration for a specific API."""
        api_name = api_name or cls.DEFAULT_API
        return cls.APIS.get(api_name, cls.APIS[cls.DEFAULT_API])


__all__ = [
    "HAS_GDAL",
    "EPSG_DICT",
    "GDALMissingError",
    "GisError",
]
