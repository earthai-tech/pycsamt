# -*- coding: utf-8 -*-
# Author: L. Kouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt.gis.utils

GIS utilities for pycsamt (v2.0).
Provides coordinate transformation functions between
latitude/longitude and UTM, using GDAL (preferred) or PyProj.

Main functions:
  * project_point_ll2utm
  * project_point_utm2ll

Supports input formats for lat/lon:
  * 'DD:MM:SS.sss'
  * 'DD.decimal_degrees'
  * float(DD.decimal_degrees)

Examples:
  >>> from pycsamt.gis.utils import project_point_ll2utm
  >>> project_point_ll2utm(-118.34, 34.05)
"""

import numpy as np
from typing import Optional, Tuple

from ..decorators import Deprecated
from ..log.logger import get_logger
from .config import HAS_GDAL, EPSG_DICT
from .constants import (
    DEG2RAD,
    RAD2DEG,
    ELLIPSOIDS,
    _EQUATORIAL_RADIUS_IDX, 
    _ECC_SQUARED_IDX, 
    utm_letter_designator,
)

# GDAL or PyProj backend setup
if HAS_GDAL:
    from osgeo import osr
    from osgeo.ogr import OGRERR_NONE
else:
    import pyproj  # type: ignore

logger = get_logger(__name__)


__all__ = [
    'project_point_ll2utm',
    'project_point_utm2ll',
    'll_to_utm',
    'utm_to_ll',
    'get_utm_string_from_sr',
    'get_utm_zone',
    'utm_zone_to_epsg',
    'get_epsg',
    'assert_xy_coordinate_system', 
    'decimal_to_dms', 
    'dms_to_decimal'
    
]


class GisError(Exception):
    """Base exception for GIS utilities."""
    pass

def assert_xy_coordinate_system (x, y ): 
    """ Assert the coordinates system of x and y. 
    
    Parameters 
    ------------
    x, y : arraylike 1d 
       Array of position coordinates x and y  
    
    Returns 
    ---------
    cs: str
       Coordinates system  ['utm'| 'dms'|'ll'] as:
       
       -'ll' for longitude -latitude coordinate system 
       - 'dms' for degree-minute-second (DD:MM:SS) coordinate system.
       - 'utm' for 'UTM' coordinate system. 
      
    Note 
    -------
    Note that any other values that does not fit longitude-latitude ('ll') 
    or Degree-minute-seconds ('DD:MM:SS') should be considered as 
    'UTM' coordinate system. 
       
    Examples 
    --------
    >>> import numpy as np 
    >>> np.random.seed (42 ) 
    >>> from pycsamt.gis.utils import assert_xy_coordinate_system
    >>> x, y = np.random.rand(7 ), np.arange (7 )
    >>> assert_xy_coordinate_system (x, y)
    'll'
    >>> assert_xy_coordinate_system (x, y*180) 
    'utm'
    >>> x = ['28:24:43.08','28:24:42.69','28:24:42.31'] 
    >>> y = ['109:19:58.34','109:19:58.93','109:19:59.51'] 
    >>> assert_xy_coordinate_system (x, y) 
    'dms'
    """
    def isdms (a ): 
        return all ([':' in str (s) for s in a ] )
        
    cs ='utm'
    x = np.array (x ) ; y = np.array (y) # for consistency 
    
    if ( 
            isdms (x) or isdms (y) ) : 
        cs ='dms'
    
    elif (
             ( (x < 90).all () and  (y <180).all())
          or (( x <180).all() and (y <90).all())
          ): 
        cs ='ll'
    
        
    return cs 

# DMS (degrees, minutes, seconds) parsing helpers

def _assert_minutes(minutes: float) -> float:
    """
    Validate minutes component: 0 <= minutes < 60.
    """
    if not (0 <= minutes < 60):
        raise ValueError(f"Minutes must be in [0,60): {minutes}")
    return minutes


def _assert_seconds(seconds: float) -> float:
    """
    Validate seconds component: 0 <= seconds < 60.
    """
    if not (0 <= seconds < 60):
        raise ValueError(f"Seconds must be in [0,60): {seconds}")
    return seconds


def _rollover_dms(unit: float, value: float) -> Tuple[float, float]:
    """
    Handle rollover when seconds or minutes >= 60 by carrying
    over to higher unit.

    Returns (unit + carry, remainder).
    """
    carry = int(value // 60)
    remainder = value % 60
    if carry:
        return unit + carry, remainder
    return unit, remainder


def dms_to_decimal( position_str): 
    "Convert a DMS string ('DD:MM:SS.sss') or decimal degrees to a float."
    return convert_position_str2float (position_str) 

def convert_position_str2float(position_str: str) -> Optional[float]:
    """
    Convert a DMS string ('DD:MM:SS.sss') or decimal degrees to a float.

    Parameters
    ----------
    position_str : str
        DMS string or decimal degrees string/float.

    Returns
    -------
    float or None
        Decimal degrees value, or None for invalid input.
    """
    if position_str in (None, 'None'):
        return None
    # Try decimal
    try:
        return float(position_str)
    except (TypeError, ValueError):
        pass

    # Must be DMS format
    parts = position_str.split(':')
    if len(parts) != 3:
        raise ValueError(f"Invalid DMS format, expected DD:MM:SS: {position_str}")
    deg = float(parts[0])
    minutes = _assert_minutes(float(parts[1]))
    seconds = _assert_seconds(float(parts[2]))

    minutes, seconds = _rollover_dms(minutes, seconds)
    deg, minutes = _rollover_dms(deg, minutes)

    sign = -1 if deg < 0 else 1
    return sign * (abs(deg) + minutes/60.0 + seconds/3600.0)


def assert_lat_value(latitude):
    """
    make sure latitude is in decimal degrees
    """
    if latitude in [None, 'None']:
        return None
    try:
        lat_value = float(latitude)

    except TypeError:
        return None

    except ValueError:
        lat_value = convert_position_str2float(latitude)

    if abs(lat_value) >= 90:
        print("==> The lat_value =", lat_value)
        raise ValueError(f'|Latitude| > 90, unacceptable!: {str(lat_value)!r}')

    return lat_value


def assert_lon_value(longitude):
    """
    make sure longitude is in decimal degrees
    """
    
    if longitude in [None, 'None']:
        return None
    try:
        lon_value = float(longitude)

    except TypeError:
        return None

    except ValueError:
        lon_value = convert_position_str2float(longitude)

    if abs(lon_value) >= 180:
        raise ValueError(f'|Longitude| > 180, unacceptable!: {str(lon_value)!r}')

    return lon_value


def assert_elevation_value(elevation):
    """
    make sure elevation is a floating point number
    """

    try:
        elev_value = float(elevation)
    except (ValueError, TypeError):
        elev_value = 0.0
        logger.warn('{0} is not a number, setting elevation to 0'.format(elevation))

    return elev_value

def decimal_to_dms( position): 
    """convert position float to a string in the format of DD:MM:SS"""
    return convert_position_float2str (position )

def convert_position_float2str(position):
    """
    convert position float to a string in the format of DD:MM:SS
    
    Arguments
    -------------
        **position** : float
                       decimal degrees of latitude or longitude
                       
    Returns
    --------------
        **position_str** : string
                          latitude or longitude in format of DD:MM:SS.ms
                          
    Example
    -------------
        >>> import mtpy.utils.gis_tools as gis_tools
        >>> gis_tools.convert_position_float2str(-118.34563)
        
    """

    assert type(position) is float, 'Given value is not a float'

    deg = int(position)
    sign = 1
    if deg < 0:
        sign = -1

    deg = abs(deg)
    minutes = (abs(position) - deg) * 60.
    # need to round seconds to 4 decimal places otherwise machine precision
    # keeps the 60 second roll over and the string is incorrect.
    sec = np.round((minutes - int(minutes)) * 60., 4)
    if sec >= 60.:
        minutes += 1
        sec = 0

    if int(minutes) == 60:
        deg += 1
        minutes = 0
        
    position_str = '{0}:{1:02.0f}:{2:05.2f}'.format(sign * int(deg),
                                                    int(minutes),
                                                    sec)

    return position_str



@Deprecated(
    "Use of GDAL SpatialReference for UTM string is deprecated;"  
    "Refer to get_utm_zone for standard UTM formatting."
 )
def get_utm_string_from_sr(spatial_ref: 'osr.SpatialReference') -> str:
    """
    Return UTM zone string (e.g. '11N') from a GDAL SpatialReference.
    """
    zone = spatial_ref.GetUTMZone()
    if zone > 0:
        return f"{zone}N"
    elif zone < 0:
        return f"{abs(zone)}S"
    return str(zone)


def get_utm_zone(latitude: float, longitude: float) -> Tuple[int, bool, str]:
    """
    Determine the UTM zone number, hemisphere, and zone string for
    given latitude and longitude.

    Returns
    -------
    zone_number : int
    is_northern : bool
    zone_string : str (e.g. '11N')
    """
    # zone width = 6° starting at -180°
    zone_number = int((longitude + 180) / 6) + 1
    # handle wrap-around
    zone_number = ((zone_number - 1) % 60) + 1
    is_northern = latitude >= 0
    letter = utm_letter_designator(latitude)
    zone_string = f"{zone_number}{letter}"
    return zone_number, is_northern, zone_string


def utm_zone_to_epsg(
        zone_number: int, is_northern: bool
        ) -> Optional[int]:
    """
    Lookup EPSG code for WGS84 UTM zone.

    Parameters
    ----------
    zone_number : int
    is_northern : bool
    """
    datum_flag = '' if is_northern else '+south'
    for epsg_code, proj4 in EPSG_DICT.items():
        if ( 
                f"+zone={zone_number}" in proj4 
                and '+datum=WGS84' in proj4 
                and datum_flag in proj4
            ):
            return epsg_code
    return None


def get_epsg(latitude: float, longitude: float) -> Optional[int]:
    """
    Get EPSG code for the UTM projection (WGS84) of a latitude/longitude.
    """
    zone_number, is_northern, _ = get_utm_zone(latitude, longitude)
    return utm_zone_to_epsg(zone_number, is_northern)


def project_point_ll2utm(lat, lon, datum='WGS84', utm_zone=None, epsg=None):
    """
    Project a point that is in Lat, Lon (will be converted to decimal degrees)
    into UTM coordinates.

    Arguments:
    ---------------
        **lat** : float or string (DD:MM:SS.ms)
                  latitude of point

        **lon** : float or string (DD:MM:SS.ms)
                  longitude of point

        **datum** : string
                    well known datum ex. WGS84, NAD27, NAD83, etc.

        **utm_zone** : string
                       zone number and 'S' or 'N' e.g. '55S'

        **epsg** : int
                   epsg number defining projection (see
                   http://spatialreference.org/ref/ for moreinfo)
                   Overrides utm_zone if both are provided

    Returns:
    --------------
        **proj_point**: tuple(easting, northing, zone)
                        projected point in UTM in Datum

    """
    if lat is None or lon is None:
        return None, None, None

    # make sure the lat and lon are in decimal degrees
    if np.iterable(lon) and np.iterable(lat):
        lat = np.array([assert_lat_value(lat_value) for lat_value in lat])
        lon = np.array([assert_lon_value(lon_value) for lon_value in lon])
        assert lat.size == lon.size
    else:
        lat = np.array([assert_lat_value(lat)])
        lon = np.array([assert_lon_value(lon)])

    if HAS_GDAL:
        # set lat lon coordinate system
        ll_cs = osr.SpatialReference()
        if isinstance(datum, int):
            ogrerr = ll_cs.ImportFromEPSG(datum)
            if ogrerr != OGRERR_NONE:
                raise GisError("GDAL/osgeo ogr error code: {}".format(ogrerr))
        elif isinstance(datum, str):
            ogrerr = ll_cs.SetWellKnownGeogCS(datum)
            if ogrerr != OGRERR_NONE:
                raise GisError("GDAL/osgeo ogr error code: {}".format(ogrerr))
        else:
            raise GisError("""datum {0} not understood, needs to be EPSG as int
                               or a well known datum as a string""".format(datum))

        # set utm coordinate system
        utm_cs = osr.SpatialReference()
    # end if

    # project point on to EPSG coordinate system if given
    if isinstance(epsg, int):
        if HAS_GDAL:
            ogrerr = utm_cs.ImportFromEPSG(epsg)
            if ogrerr != OGRERR_NONE:
                raise GisError("GDAL/osgeo ogr error code: {}".format(ogrerr))
        else:
            pp = pyproj.Proj('+init=EPSG:%d'%(epsg))
        # end if
    # otherwise project onto given datum
    elif epsg is None:
        if HAS_GDAL:
            ogrerr = utm_cs.CopyGeogCSFrom(ll_cs)
            if ogrerr != OGRERR_NONE:
                raise GisError("GDAL/osgeo ogr error code: {}".format(ogrerr))
        # end if
        if utm_zone is None or not isinstance(None, str) or utm_zone.lower() == 'none':
            # get the UTM zone in the datum coordinate system, otherwise
            zone_number, is_northern, utm_zone = get_utm_zone(lat.mean(),
                                                              lon.mean())
        else:
            # get zone number and is_northern from utm_zone string
            zone_number = int(utm_zone[0:-1])
            is_northern = True if utm_zone[-1].lower() > 'n' else False

        if(HAS_GDAL):
            utm_cs.SetUTM(zone_number, is_northern)
        else:
            projstring = '+proj=utm +zone=%d +%s +datum=%s' % \
                         (zone_number, 'north' if is_northern else 'south', datum)
            pp = pyproj.Proj(projstring)
        # end if
    # end if

    # return different results depending on if lat/lon are iterable
    projected_point = np.zeros_like(lat, dtype=[('easting', np.float),
                                                ('northing', np.float),
                                                ('elev', np.float),
                                                ('utm_zone', 'U4')])

    if(HAS_GDAL):
        ll2utm = osr.CoordinateTransformation(ll_cs, utm_cs).TransformPoint
    else:
        ll2utm = pp
    # end if

    for ii in range(lat.size):
        point = ll2utm(lon[ii], lat[ii])
        projected_point['easting'][ii] = point[0]
        projected_point['northing'][ii] = point[1]
        if(HAS_GDAL): projected_point['elev'][ii] = point[2]
        projected_point['utm_zone'][ii] = utm_zone if utm_zone is not None else get_utm_zone(lat[ii], lon[ii])[2]
    # end for

    # if just projecting one point, then return as a tuple so as not to break
    # anything.  In the future we should adapt to just return a record array
    if len(projected_point) == 1:
        return (projected_point['easting'][0],
                projected_point['northing'][0],
                projected_point['utm_zone'][0])
    else:
        return np.rec.array(projected_point)

def project_point_utm2ll(easting, northing, utm_zone, datum='WGS84', epsg=3149):
    """
    Project a point that is in Lat, Lon (will be converted to decimal degrees)
    into UTM coordinates.
    
    Arguments:
    ---------------
        **easting** : float
                    easting coordinate in meters
                    
        **northing** : float
                    northing coordinate in meters
        
        **utm_zone** : string (##N or ##S)
                      utm zone in the form of number and North or South
                      hemisphere, 10S or 03N
        
        **datum** : string
                    well known datum ex. WGS84, NAD27, etc.
                    
    Returns:
    --------------
        **proj_point**: tuple(lat, lon)
                        projected point in lat and lon in Datum, as decimal
                        degrees.
                    
    """
    try:
        easting = float(easting)
    except ValueError:
        raise GisError("easting is not a float")
    try:
        northing = float(northing)
    except ValueError:
        raise GisError("northing is not a float")

    if HAS_GDAL:
        # set utm coordinate system
        utm_cs = osr.SpatialReference()
        utm_cs.SetWellKnownGeogCS(datum)
    # end if

    if epsg is not None:
        if HAS_GDAL:
            ogrerr = utm_cs.ImportFromEPSG(epsg)
            if ogrerr != OGRERR_NONE:
                raise Exception("GDAL/osgeo ogr error code: {}".format(ogrerr))
        else:
            pp = pyproj.Proj('+init=EPSG:%d'%(epsg))
        # end if
    elif isinstance(utm_zone, str) or isinstance(utm_zone, np.bytes_):
        # the isinstance(utm_zone, str) could be False in python3 due to numpy datatype change.
        # So FZ added  isinstance(utm_zone, np.bytes_) and convert the utm_zone into string
        if isinstance(utm_zone, np.bytes_):
            utm_zone = utm_zone.decode('UTF-8') # b'54J'
        try:
            zone_number = int(utm_zone[0:-1])  #b'54J'
            zone_letter = utm_zone[-1]
        except ValueError:
            raise ValueError('Zone number {0} is not a number'.format(utm_zone[0:-1]))
        is_northern = True if zone_letter.lower() >= 'n' else False
    elif isinstance(utm_zone, int):
        # std UTM code returned by gdal
        is_northern = False if utm_zone < 0 else True
        zone_number = abs(utm_zone)
    else:
        print("epsg and utm_zone", str(epsg), str(utm_zone))

        raise NotImplementedError("utm_zone type (%s, %s) not supported"%(type(utm_zone), str(utm_zone)))
    
    if epsg is None:
        if HAS_GDAL:
            utm_cs.SetUTM(zone_number, is_northern)
        else:
            projstring = '+proj=utm +zone=%d +%s +datum=%s' % \
                         (zone_number, 'north' if is_northern else 'south', datum)
            pp = pyproj.Proj(projstring)
        # end if
    # end if

    if HAS_GDAL:
        ll_cs = utm_cs.CloneGeogCS()
        utm2ll = osr.CoordinateTransformation(utm_cs, ll_cs).TransformPoint
        ll_point = list(utm2ll(easting, northing, 0.))
    else:
        ll_point = pp(easting, northing, inverse=True)
    # end if

    # be sure to round out the numbers to remove computing with floats
    return round(ll_point[1], 6), round(ll_point[0], 6)


def project_points_ll2utm(lat, lon, datum='WGS84', utm_zone=None, epsg=None):
    """
    Project a list of points that is in Lat, Lon (will be converted to decimal 
    degrees) into UTM coordinates.
    
    Arguments:
    ---------------
        **lat** : float or string (DD:MM:SS.ms)
                  latitude of point
                  
        **lon** : float or string (DD:MM:SS.ms)
                  longitude of point
        
        **datum** : string
                    well known datum ex. WGS84, NAD27, NAD83, etc.

        **utm_zone** : string
                       zone number and 'S' or 'N' e.g. '55S'. Defaults to the
                       centre point of the provided points
                       
        **epsg** : int
                   epsg number defining projection (see 
                   http://spatialreference.org/ref/ for moreinfo)
                   Overrides utm_zone if both are provided

    Returns:
    --------------
        **proj_point**: tuple(easting, northing, zone)
                        projected point in UTM in Datum
                    
    """

    lat = np.array(lat)
    lon = np.array(lon)

    # check length of arrays
    if np.shape(lat) != np.shape(lon):
        raise ValueError("latitude and longitude arrays are of different lengths")

    # flatten, if necessary
    flattened = False
    llshape = np.shape(lat)
    if llshape[0] > 1:
        flattened = True
        lat = lat.flatten()
        lon = lon.flatten()

    '''
    # check lat/lon values
    # this is incredibly slow; disabling for the time being
    for ii in range(len(lat)):
        lat[ii] = assert_lat_value(lat[ii])
        lon[ii] = assert_lon_value(lon[ii])
    '''
    
    if lat is None or lon is None:
        return None, None, None

    if HAS_GDAL:
        # set utm coordinate system
        utm_cs = osr.SpatialReference()
        utm_cs.SetWellKnownGeogCS(datum)

        # set lat, lon coordinate system
        ll_cs = utm_cs.CloneGeogCS()
        ll_cs.ExportToPrettyWkt()
    # end if

    # get zone number, north and zone name
    if epsg is not None:
        if HAS_GDAL:
            # set projection info
            ogrerr = utm_cs.ImportFromEPSG(epsg)
            if ogrerr != OGRERR_NONE:
                raise Exception("GDAL/osgeo ogr error code: {}".format(ogrerr))
            # get utm zone (for information) if applicable
            utm_zone = utm_cs.GetUTMZone()

            # Whilst some projections e.g. Geoscience Australia Lambert (epsg3112) do
            # not yield UTM zones, they provide eastings and northings for the whole
            # Australian region. We therefore set UTM zones, only when a valid UTM zone
            # is available
            if(utm_zone>0):
                # set projection info
                utm_cs.SetUTM(abs(utm_zone), utm_zone > 0)
        else:
            pp = pyproj.Proj('+init=EPSG:%d'%(epsg))
        # end if
    else:
        if utm_zone is not None:
            # get zone number and is_northern from utm_zone string
            zone_number = int(utm_zone[0:-1])
            is_northern = True if utm_zone[-1].lower() > 'n' else False
        else:
            # get centre point and get zone from that
            latc = (np.nanmax(lat) + np.nanmin(lat)) / 2.
            lonc = (np.nanmax(lon) + np.nanmin(lon)) / 2.
            zone_number, is_northern, utm_zone = get_utm_zone(latc, lonc)
        # set projection info

        if(HAS_GDAL):
            utm_cs.SetUTM(zone_number, is_northern)
        else:
            projstring = '+proj=utm +zone=%d +%s +datum=%s' % \
                         (zone_number, 'north' if is_northern else 'south', datum)
            pp = pyproj.Proj(projstring)
        # end if
    # end if
    
    if HAS_GDAL:
        ll2utm = osr.CoordinateTransformation(ll_cs, utm_cs).TransformPoints
        easting, northing, elev = np.array(ll2utm(np.array([lon, lat]).T)).T

    else:
        ll2utm = pp
        easting, northing = ll2utm(lon, lat)
    # end if

    projected_point = (easting, northing, utm_zone)

    # reshape back into original shape
    if flattened:
        lat = lat.reshape(llshape)
        lon = lon.reshape(llshape)

    return projected_point

def ll_to_utm(
    reference_ellipsoid: int,
    lat: float,
    lon: float
) -> Tuple[str, float, float]:
    """
    Convert latitude/longitude to UTM coordinates.
    Equations from USGS Bulletin 1532.

    Parameters
    ----------
    reference_ellipsoid : int
        ID of ellipsoid from ELLIPSOIDS (e.g. 23 for WGS-84).
    lat, lon : float
        Geographic coordinates in decimal degrees.

    Returns
    -------
    utm_zone : str
        UTM zone string (e.g. '11N').
    easting : float
    northing : float
    """
    # Get ellipsoid parameters
    try:
        ell = next(e for e in ELLIPSOIDS if e[0] == reference_ellipsoid)
    except StopIteration:
        raise GisError(f"Unknown ellipsoid ID: {reference_ellipsoid}")
        
    a = ell[_EQUATORIAL_RADIUS_IDX]
    ecc_sq = ell[_ECC_SQUARED_IDX]
    k0 = 0.9996

    # normalize longitude into [-180, 180)
    long_temp = ((lon + 180) % 360) - 180

    lat_rad = lat * DEG2RAD
    long_rad = long_temp * DEG2RAD

    # determine UTM zone number
    zone_number = int((long_temp + 180) / 6) + 1
    # Norway exception
    if 56.0 <= lat < 64.0 and 3.0 <= long_temp < 12.0:
        zone_number = 32
    # Svalbard exceptions
    if 72.0 <= lat < 84.0:
        if   0.0  <= long_temp <  9.0: zone_number = 31
        elif  9.0 <= long_temp < 21.0: zone_number = 33
        elif 21.0 <= long_temp < 33.0: zone_number = 35
        elif 33.0 <= long_temp < 42.0: zone_number = 37

    long_origin = (zone_number - 1) * 6 - 180 + 3
    long_origin_rad = long_origin * DEG2RAD

    utm_zone = f"{zone_number}{utm_letter_designator(lat)}"

    ecc_prime_sq = ecc_sq / (1 - ecc_sq)
    N = a / np.sqrt(1 - ecc_sq * np.sin(lat_rad) ** 2)
    T = np.tan(lat_rad) ** 2
    C = ecc_prime_sq * np.cos(lat_rad) ** 2
    A = np.cos(lat_rad) * (long_rad - long_origin_rad)

    # meridional arc
    M = a * (
        (1 - ecc_sq/4 - 3*ecc_sq**2/64 - 5*ecc_sq**3/256) * lat_rad
        - (3*ecc_sq/8 + 3*ecc_sq**2/32 + 45*ecc_sq**3/1024) * np.sin(2*lat_rad)
        + (15*ecc_sq**2/256 + 45*ecc_sq**3/1024) * np.sin(4*lat_rad)
        - (35*ecc_sq**3/3072) * np.sin(6*lat_rad)
    )

    # easting
    easting = (
        k0 * N * (
            A
            + (1 - T + C) * A**3 / 6
            + (5 - 18*T + T**2 + 72*C - 58*ecc_prime_sq) * A**5 / 120
        )
        + 500_000.0
    )

    # northing
    northing = k0 * (
        M + N * np.tan(lat_rad) * (
            A**2 / 2
            + (5 - T + 9*C + 4*C**2) * A**4 / 24
            + (61 - 58*T + T**2 + 600*C - 330*ecc_prime_sq) * A**6 / 720
        )
    )
    # southern hemisphere offset
    if lat < 0:
        northing += 10_000_000.0

    return utm_zone, easting, northing

def utm_to_ll(
    reference_ellipsoid: int,
    northing: float,
    easting: float,
    zone: str
) -> Tuple[float, float]:
    """
    Convert UTM coordinates to latitude/longitude.
    Inverse of ll_to_utm, equations from USGS Bulletin 1532.

    Parameters
    ----------
    reference_ellipsoid : int
        ID of ellipsoid from ELLIPSOIDS.
    northing, easting : float
        UTM coordinates in meters.
    zone : str
        UTM zone string (e.g. '11N').

    Returns
    -------
    lat, lon : float
        Geographic coordinates in decimal degrees.
    """
    # Get ellipsoid parameters
    try:
        ell = next(e for e in ELLIPSOIDS if e[0] == reference_ellipsoid)
    except StopIteration:
        raise GisError(f"Unknown ellipsoid ID: {reference_ellipsoid}")
    a = ell[_EQUATORIAL_RADIUS_IDX]
    ecc_sq = ell[_ECC_SQUARED_IDX]
    k0 = 0.9996

    # Compute e1 parameter
    e1 = (1 - np.sqrt(1 - ecc_sq)) / (1 + np.sqrt(1 - ecc_sq))

    # Remove offsets
    x = easting - 500_000.0
    y = northing

    # Determine hemisphere
    zone_letter = zone[-1]
    if zone_letter.upper() < 'N':
        # southern hemisphere
        y -= 10_000_000.0

    # Central meridian of zone
    zone_number = int(zone[:-1])
    long_origin = (zone_number - 1) * 6 - 180 + 3

    # Calculate footpoint latitude
    M = y / k0
    mu = M / (a * (1 - ecc_sq/4 - 3*ecc_sq**2/64 - 5*ecc_sq**3/256))
    phi1_rad = (
        mu
        + (3*e1/2 - 27*e1**3/32) * np.sin(2*mu)
        + (21*e1**2/16 - 55*e1**4/32) * np.sin(4*mu)
        + (151*e1**3/96) * np.sin(6*mu)
    )

    # Precompute
    ecc_prime_sq = ecc_sq / (1 - ecc_sq)
    n1 = a / np.sqrt(1 - ecc_sq * np.sin(phi1_rad)**2)
    t1 = np.tan(phi1_rad)**2
    c1 = ecc_prime_sq * np.cos(phi1_rad)**2
    r1 = a * (1 - ecc_sq) / (1 - ecc_sq * np.sin(phi1_rad)**2)**1.5
    d  = x / (n1 * k0)

    # Latitude
    lat_rad = (
        phi1_rad
        - (n1 * np.tan(phi1_rad) / r1) * (
            d**2/2
            - (5 + 3*t1 + 10*c1 - 4*c1**2 - 9*ecc_prime_sq) * d**4/24
            + (61 + 90*t1 + 298*c1 + 45*t1**2 - 252*ecc_prime_sq - 3*c1**2) * d**6/720
        )
    )
    lat = lat_rad * RAD2DEG

    # Longitude
    lon_rad = (
        d
        - (1 + 2*t1 + c1) * d**3/6
        + (5 - 2*c1 + 28*t1 - 3*c1**2 + 8*ecc_prime_sq + 24*t1**2) * d**5/120
    ) / np.cos(phi1_rad)
    lon = long_origin + lon_rad * RAD2DEG

    return lat, lon


def epsg_project(x, y, epsg_from, epsg_to):
    """
    project some xy points using the pyproj modules
    """

    try:
        import pyproj
    except ImportError:
        print("please install pyproj")
        return
    if epsg_from is not None:
        try:
            p1 = pyproj.Proj(EPSG_DICT[epsg_from])
            p2 = pyproj.Proj(EPSG_DICT[epsg_to])
        except KeyError:
            print("Surface or data epsg either not in dictionary or None")
            return

    return pyproj.transform(p1, p2, x, y)


def utm_wgs84_conv(lat, lon):
    """
    Bidirectional UTM-WGS84 converter https://github.com/Turbo87/utm/blob/master/utm/conversion.py
    :param lat:
    :param lon:
    :return: tuple(e, n, zone, lett)
    """

    import utm  # pip install utm
    tup = utm.from_latlon(lat, lon)

    (new_lat, new_lon) = utm.to_latlon(tup[0], tup[1], tup[2], tup[3])
    # print (new_lat,new_lon)  # should be same as the input param

    # checking correctess
    if abs(lat - new_lat) > 1.0 * np.e - 10:
        print("Warning: lat and new_lat should be equal!")

    if abs(lon - new_lon) > 1.0 * np.e - 10:
        print("Warning: lon and new_lon should be equal!")

    return tup


#@gdal_data_check
#@deprecated("This function may be removed in later release. mtpy.utils.gis_tools.project_point_utm2ll() should be "
#            "used instead.")
def transform_utm_to_ll(easting, northing, zone,
                        reference_ellipsoid='WGS84'):
    utm_coordinate_system = osr.SpatialReference()
    # Set geographic coordinate system to handle lat/lon
    utm_coordinate_system.SetWellKnownGeogCS(reference_ellipsoid)

    try:
        zone_number = int(zone[0:-1])
        zone_letter = zone[-1]
    except ValueError:
        raise ValueError('Zone number {0} is not a number'.format(zone[0:-1]))
    is_northern = True if zone_letter.lower() >= 'n' else False

    utm_coordinate_system.SetUTM(zone_number, is_northern)

    # Clone ONLY the geographic coordinate system
    ll_coordinate_system = utm_coordinate_system.CloneGeogCS()

    # create transform component
    utm_to_ll_geo_transform = osr.CoordinateTransformation(utm_coordinate_system,
                                                           ll_coordinate_system)
    # returns lon, lat, altitude
    return utm_to_ll_geo_transform.TransformPoint(easting, northing, 0)


#@gdal_data_check
@Deprecated("This function may be removed in later release. mtpy.utils.gis_tools.project_point_ll2utm() should be "
            "used instead.")
def transform_ll_to_utm(lon, lat, reference_ellipsoid='WGS84'):
    """
    transform a (lon,lat) to  a UTM coordinate.
    The UTM zone number will be determined by longitude. South-North will be determined by Lat.
    :param lon: degree
    :param lat: degree
    :param reference_ellipsoid:
    :return: utm_coordinate_system, utm_point
    """

    def get_utm_zone(longitude):
        return (int(1 + (longitude + 180.0) / 6.0))

    def is_northern(latitude):
        """
        Determines if given latitude is a northern for UTM
        """
        # if (latitude < 0.0):
        #     return 0
        # else:
        #     return 1
        return latitude >= 0

    utm_coordinate_system = osr.SpatialReference()
    # Set geographic coordinate system to handle lat/lon
    utm_coordinate_system.SetWellKnownGeogCS(reference_ellipsoid)
    utm_coordinate_system.SetUTM(get_utm_zone(lon), is_northern(lat))

    # Clone ONLY the geographic coordinate system
    ll_coordinate_system = utm_coordinate_system.CloneGeogCS()
    # create transform component
    ll_to_utm_geo_transform = osr.CoordinateTransformation(ll_coordinate_system,
                                                           utm_coordinate_system)

    utm_point = ll_to_utm_geo_transform.TransformPoint(lon, lat, 0)

    # returns easting, northing, altitude
    return utm_coordinate_system, utm_point





