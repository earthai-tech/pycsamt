# Author: Kouadio Laurent alias Daniel <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.gis.constants

Static GIS constants and lookup tables for pycsamt:
 - UTM_ZONE_DESIGNATOR: letter -> [lon_min, lon_max]
 - Ellipsoid parameters and conversion factors
 - utm_letter_designator: latitude -> UTM zone letter
"""

import numpy as np

#
# UTM zone designator ranges: letter -> [min_longitude, max_longitude]
UTM_ZONE_DESIGNATOR = {
    "X": [72, 84],
    "W": [64, 72],
    "V": [56, 64],
    "U": [48, 56],
    "T": [40, 48],
    "S": [32, 40],
    "R": [24, 32],
    "Q": [16, 24],
    "P": [8, 16],
    "N": [0, 8],
    "M": [-8, 0],
    "L": [-16, -8],
    "K": [-24, -16],
    "J": [-32, -24],
    "H": [-40, -32],
    "G": [-48, -40],
    "F": [-56, -48],
    "E": [-64, -56],
    "D": [-72, -64],
    "C": [-80, -72],
    "Z": [-80, 84],
}

# Conversion constants
# -
DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi

# Indices for ellipsoid parameter tuples
_EQUATORIAL_RADIUS_IDX = 2
_ECC_SQUARED_IDX = 3

# EPSG Proj4 string to zone mapping (spatialreference.org)
EPSG_PROJ4 = {
    28350: [
        "+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        50,
    ],
    28351: [
        "+proj=utm +zone=51 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        51,
    ],
    28352: [
        "+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        52,
    ],
    28353: [
        "+proj=utm +zone=53 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        53,
    ],
    28354: [
        "+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        54,
    ],
    28355: [
        "+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        55,
    ],
    28356: [
        "+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        56,
    ],
    3112: [
        "+proj=lcc +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        0,
    ],
    4326: ["+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 0],
    32609: [
        "+proj=utm +zone=9  +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        9,
    ],
    32610: [
        "+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        10,
    ],
    32611: [
        "+proj=utm +zone=11 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        11,
    ],
    32612: [
        "+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        12,
    ],
    32613: [
        "+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        13,
    ],
    32614: [
        "+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        14,
    ],
    32615: [
        "+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        15,
    ],
    32616: [
        "+proj=utm +zone=16 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        16,
    ],
    32617: [
        "+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        17,
    ],
    32618: [
        "+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        18,
    ],
    32619: [
        "+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        19,
    ],
}

# Reference ellipsoids: [id, name, equatorial_radius (m), eccentricity_squared]
ELLIPSOIDS = [
    [-1, "Placeholder", 0, 0],
    [1, "Airy", 6377563.396, 0.00667054],
    [2, "Australian National", 6378160.000, 0.006694542],
    [3, "Bessel 1841", 6377397.155, 0.006674372],
    [4, "Bessel 1841 (Namibia)", 6377484.000, 0.006674372],
    [5, "Clarke 1866", 6378206.400, 0.006768658],
    [6, "Clarke 1880", 6378249.145, 0.006803511],
    [7, "Everest", 6377276.345, 0.006637847],
    [8, "Fischer 1960 (Mercury)", 6378166.000, 0.006693422],
    [9, "Fischer 1968", 6378150.000, 0.006693422],
    [10, "GRS 1967", 6378160.000, 0.006694605],
    [11, "GRS 1980", 6378137.000, 0.006694380],
    [12, "Helmert 1906", 6378200.000, 0.006693422],
    [13, "Hough", 6378270.000, 0.006722670],
    [14, "International", 6378388.000, 0.006722670],
    [15, "Krassovsky", 6378245.000, 0.006693422],
    [16, "Modified Airy", 6377340.000, 0.006670540],
    [17, "Modified Everest", 6377304.000, 0.006637847],
    [18, "Modified Fischer 1960", 6378155.000, 0.006693422],
    [19, "South American 1969", 6378160.000, 0.006694542],
    [20, "WGS 60", 6378165.000, 0.006693422],
    [21, "WGS 66", 6378145.000, 0.006694542],
    [22, "WGS-72", 6378135.000, 0.006694318],
    [23, "WGS-84", 6378137.000, 0.006694380],
]


# UTM letter designator by latitude
def utm_letter_designator(lat: float) -> str:
    """
    Returns the UTM zone letter for a given latitude.
    'Z' if latitude outside 80°S to 84°N.
    """
    if 84 >= lat >= 72:
        return "X"
    elif 72 > lat >= 64:
        return "W"
    elif 64 > lat >= 56:
        return "V"
    elif 56 > lat >= 48:
        return "U"
    elif 48 > lat >= 40:
        return "T"
    elif 40 > lat >= 32:
        return "S"
    elif 32 > lat >= 24:
        return "R"
    elif 24 > lat >= 16:
        return "Q"
    elif 16 > lat >= 8:
        return "P"
    elif 8 > lat >= 0:
        return "N"
    elif 0 > lat >= -8:
        return "M"
    elif -8 > lat >= -16:
        return "L"
    elif -16 > lat >= -24:
        return "K"
    elif -24 > lat >= -32:
        return "J"
    elif -32 > lat >= -40:
        return "H"
    elif -40 > lat >= -48:
        return "G"
    elif -48 > lat >= -56:
        return "F"
    elif -56 > lat >= -64:
        return "E"
    elif -64 > lat >= -72:
        return "D"
    elif -72 > lat >= -80:
        return "C"
    else:
        return "Z"


__all__ = [
    "UTM_ZONE_DESIGNATOR",
    "DEG2RAD",
    "RAD2DEG",
    "ELLIPSOIDS",
    "EPSG_PROJ4",
    "_ECC_SQUARED_IDX",
    "_EQUATORIAL_RADIUS_IDX",
    "utm_letter_designator",
]
