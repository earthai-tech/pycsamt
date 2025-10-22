# -*- coding: utf-8 -*-
# Author: L. Kouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt.site.location

Station location handling: conversion between
latitude/longitude and UTM, elevation, and
computed azimuth.
"""

from typing import Optional, Union, Tuple

import numpy as np  

from ..api.typing import ArrayLike 
from ..log.logger import get_logger
from ..exceptions import LocationError
from ..gis.utils import ll_to_utm, utm_to_ll
from ..utils.arrayops import is_iterable, assert_xy_in 
from ..utils.validation import _check_consistency_size 

logger = get_logger(__name__)

__all__ = ['Location']

class Location:
    """
    Encapsulates a station's geographic and UTM coordinates.

    Attributes
    ----------
    latitude : float
    longitude : float
    elevation : float
    utm_zone : str
    easting : float
    northing : float
    """

    def __init__(
        self,
        latitude: Optional[float] = None,
        longitude: Optional[float] = None,
        elevation: float = 0.0,
        utm_zone: Optional[str] = None
    ):
        self._latitude: Optional[float] = None
        self._longitude: Optional[float] = None
        self._elevation: float = elevation
        self._utm_zone: Optional[str] = utm_zone
        self._easting: Optional[float] = None
        self._northing: Optional[float] = None

        if latitude is not None and longitude is not None:
            self.latitude = latitude
            self.longitude = longitude
            self.to_utm()
        elif utm_zone is not None:
            self.utm_zone = utm_zone

    @property
    def latitude(self) -> Optional[float]:
        return self._latitude

    @latitude.setter
    def latitude(self, val: Union[float, str]):
        try:
            self._latitude = float(val)
        except Exception:
            raise LocationError(f"Invalid latitude: {val}")

    @property
    def longitude(self) -> Optional[float]:
        return self._longitude

    @longitude.setter
    def longitude(self, val: Union[float, str]):
        try:
            self._longitude = float(val)
        except Exception:
            raise LocationError(f"Invalid longitude: {val}")

    @property
    def elevation(self) -> float:
        return self._elevation

    @elevation.setter
    def elevation(self, val: Union[float, str]):
        try:
            self._elevation = float(val)
        except Exception:
            raise LocationError(f"Invalid elevation: {val}")

    @property
    def utm_zone(self) -> Optional[str]:
        return self._utm_zone

    @utm_zone.setter
    def utm_zone(self, val: str):
        if not isinstance(val, str) or len(val) < 2:
            raise LocationError(f"Invalid UTM zone: {val}")
        self._utm_zone = val

    @property
    def easting(self) -> Optional[float]:
        return self._easting

    @property
    def northing(self) -> Optional[float]:
        return self._northing

    def to_utm(self) -> Tuple[str, float, float]:
        """
        Convert current lat/lon to UTM, storing easting, northing, utm_zone.
        """
        if self._latitude is None or self._longitude is None:
            raise LocationError(
                "Latitude and longitude must be set"
                " before UTM conversion.")

        zone, easting, northing = ll_to_utm(
            reference_ellipsoid=23,
            lat=self._latitude,
            lon=self._longitude
        )
        self._utm_zone = zone
        self._easting = easting
        self._northing = northing
        return zone, easting, northing

    def to_latlon(self) -> Tuple[float, float]:
        """
        Convert current easting/northing/utm_zone back to latitude/longitude.
        """
        if self._easting is None or self._northing is None or not self._utm_zone:
            raise LocationError(
                "UTM parameters must be set"
                " before lat/lon conversion.")

        lat, lon = utm_to_ll(
            reference_ellipsoid=23,
            northing=self._northing,
            easting=self._easting,
            zone=self._utm_zone
        )
        self._latitude = lat
        self._longitude = lon
        return lat, lon
    
    @staticmethod 
    def to_utm_in(
        lats:ArrayLike, 
        lons:ArrayLike, 
        *, 
        data=None, 
        utm_zone:str =None, 
        datum: str=None, 
        **kws 
        ):
        """ Convert array of longitude and latitude to array of X, y UTM 
        coordinates. 
        
        Parameters 
        ------------
        lats: arraylike 1d, 
           Array composed of latitude values 
        lons: ArrayLike 1d, 
           Array composed of longitude values. 
           
        data: pd.dataFrame
           Accepts dataframe containing the latitude and longitude coordinates
           by specifying the columns through 'lats' and 'lons' columns. 

        utm_zone: Optional, string
           zone number and 'S' or 'N' e.g. '55S'. Defaults to the centre point
           of the provided points, 
        datum: string
            well known datum ex. WGS84, NAD27, NAD83, etc.
            
        kws: dict, 
           Keywords argument passed to :meth:`~watex.site.Location.to_utm`. 
           
        Returns
        --------
        (easts, norths): Iterable object composed of easting and northing 
           coordinates. 
            
        """
        if ( 
                (isinstance (lats, str ) 
                or isinstance ( lons, str) )
                and data is None): 
            raise TypeError ("Data can't be None when latitude or longitude"
                             " has a string argumemt.")
        if data is not None: 
            lats, lons = assert_xy_in(lats, lons , data = data )      
        
        emsg = ("longitude is " if lons is None else 'latitude is') if (
            lats is None or lons is None) else "Both longitude and latitude are"
        
        if lats is None or lons is None: 
            raise TypeError (emsg + "missing.") 
        # For consistency, check iterable values 
        
        lats = is_iterable(lats, exclude_string= True, transform =True ) 
        lons = is_iterable(lons, exclude_string= True, transform= True )
        
        _check_consistency_size(lats, lons)
        lats = np.array(lats ) ; lons = np.array (lons )
        easts = np.zeros_like (lats , dtype = np.float64)
        norths = np.zeros_like (lons , dtype = np.float64)
        
        for ii, (lat, lon) in enumerate (zip (lats, lons )) : 
            Loc = Location()
            x, y  = Loc.to_utm (lat, lon , utm_zone= utm_zone , **kws )
            easts[ii] = x ; norths [ii] =y 
    
        return  easts, norths 
    
            
    @staticmethod 
    def to_latlon_in(
        easts:ArrayLike, 
        norths:ArrayLike, 
        *, 
        data=None,
        utm_zone:str=None, 
        datum: str=None, 
        **kws 
        ):
        """ 
        Convert array of longitude and latitude to array of X, y UTM 
        coordinates. 
        
        Parameters 
        ------------
        lats: arraylike 1d, 
           Array composed of latitude values 
        lons: ArrayLike 1d, 
           Array composed of longitude values. 
 
        utm_zone: Optional, string
           zone number and 'S' or 'N' e.g. '55S'. Defaults to the centre point
           of the provided points, 
           
        data: pd.dataFrame
           Accepts dataframe containing the easting and northing coordinates
           by specifying the columns through 'easts' and 'lats' columns. 
           

        datum: string
           well known datum ex. WGS84, NAD27, NAD83, etc.
            
        kws: dict, 
           Keywords argument passed to :meth:`~watex.site.Location.to_latlon`. 
           
        Returns 
        -------
        (lats, lons): Iterable object composed of latitude and longitude 
           coordinates. 
            
        """
        if (  ( isinstance (easts, str ) 
            or isinstance ( norths, str) )
            and data is None): 
            raise TypeError ("Data can't be None when easting or northing"
                             " has a string argumemt.")
        if data is not None: 
            easts , norths = assert_xy_in(easts, norths, data = data )
            
        emsg = ("easting is" if easts is None else 'northing is') if (
            easts is None or norths is None) else (
                "Both easting and northing are")
        
        if easts is None or norths is None: 
            raise TypeError (emsg + " missing.") 
            
        # For consistency, check iterable values 
        easts = is_iterable(easts, exclude_string= True,  transform =True ) 
        norths = is_iterable(norths, transform= True , exclude_string =True 
                             )
        
        _check_consistency_size(easts, norths)
        easts = np.array(easts ) ; norths = np.array (norths )
        lats_lons =[]
        for east, north in zip (easts, norths ) : 
           Loc = Location()
           lats_lons.append (
               Loc.to_latlon(east, north , utm_zone= utm_zone , **kws )
               )
           
        return tuple (zip ( *lats_lons))  
    
    def __repr__(self) -> str:
        return (
            f"Location(lat={self._latitude}, lon={self._longitude}, "
            f"elev={self._elevation}, utm_zone={self._utm_zone}, "
            f"easting={self._easting}, northing={self._northing})"
        )


