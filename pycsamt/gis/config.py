# -*- coding: utf-8 -*-
# Author: Kouadio Laurent alias Daniel <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.gis.config

Configuration loader for GIS subpackage:
 - Verifies GDAL or PyProj availability
 - Loads EPSG definitions into EPSG_DICT
"""
import os

import re

import numpy as np
from pycsamt.decorators import GdalDataCheck
from pycsamt.log.logger import get_logger

logger = get_logger(__name__)

# Determine GDAL availability once
_HAS_GDAL = GdalDataCheck._locate_gdal_data()
if not _HAS_GDAL:
    try:
        import pyproj  # noqa: F401
    except ImportError:
        raise RuntimeError("Either GDAL or PyProj must be installed for GIS features.")

# Build EPSG dictionary
EPSG_DICT = {}

try:
    import pyproj #Noqa 

    # pyproj 3+ uses proj.db; older versions provide epsg file
    data_dir = getattr(pyproj, 'pyproj_datadir', None)
    epsg_path = os.path.join(data_dir, 'epsg') if data_dir else None

    if epsg_path and os.path.isfile(epsg_path):
        with open(epsg_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.strip().startswith('#'):
                    continue
                codes = re.findall(r'<(\d+)>', line)
                if not codes:
                    continue
                code = int(codes[0])
                proj4 = re.search(r'>([^<]+)<', line)
                if proj4:
                    EPSG_DICT[code] = proj4.group(1).strip()
    else:
        raise FileNotFoundError

except Exception:
    # Fallback: load local numpy file of EPSG mappings
    here = os.path.dirname(os.path.abspath(__file__))
    epsg_file = os.path.join(here, 'epsg.npy')
    try:
        EPSG_DICT = np.load(epsg_file, allow_pickle=True).item()
    except Exception as e:
        logger.error(f"Failed to load EPSG definitions: {e}")
        EPSG_DICT = {}

# Expose availability
HAS_GDAL = _HAS_GDAL

__all__ = ['HAS_GDAL', 'EPSG_DICT']

