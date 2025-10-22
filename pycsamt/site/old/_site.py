# -*- coding: utf-8 -*-
# Author: L. Kouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt.site.site

Site: manages a collection of station locations
for easy lookup and transformation.
"""
import os
from typing import List, Dict, Optional, Union

import numpy as np
import pandas as pd

from pycsamt.site.profile import StationProfile
from pycsamt.site.location import Location
from pycsamt.utils.mathext import compute_azimuth 
from pycsamt.exceptions import SiteError


class Site:
    def __init__(
        self,
        data_fn: Optional[str] = None,
        /,
        easting: Optional[Union[List[float], np.ndarray]] = None,
        northing: Optional[Union[List[float], np.ndarray]] = None,
        lat: Optional[Union[List[float], np.ndarray]] = None,
        lon: Optional[Union[List[float], np.ndarray]] = None,
        stn_names: Optional[List[str]] = None,
        utm_zone: Optional[str] = None
    ):
        self.utm_zone = utm_zone
        self.stn_names: List[str] = []
        self.latitudes: Dict[str, float] = {}
        self.longitudes: Dict[str, float] = {}
        self.easting: Dict[str, float] = {}
        self.northing: Dict[str, float] = {}
        self.elevation: Dict[str, float] = {}
        self.azimuth: Dict[str, float] = {}

        if data_fn or (easting is not None and northing is not None):
            self.set_info(
                data_fn=data_fn,
                easting=easting,
                northing=northing,
                lat=lat,
                lon=lon,
                stn_names=stn_names,
                utm_zone=utm_zone
            )

    def set_info(
        self,
        data_fn: Optional[str] = None,
        easting: Optional[Union[List[float], np.ndarray]] = None,
        northing: Optional[Union[List[float], np.ndarray]] = None,
        lat: Optional[Union[List[float], np.ndarray]] = None,
        lon: Optional[Union[List[float], np.ndarray]] = None,
        stn_names: Optional[List[str]] = None,
        utm_zone: Optional[str] = None
    ) -> None:
        """
        Populate site from a .stn file or parallel coordinate arrays.

        Parameters
        ----------
        data_fn : str, optional
            Path to a station profile file (.stn).
        easting, northing : array_like, optional
            UTM coordinates. Must both be provided.
        lat, lon : array_like, optional
            Geographic coordinates. Must both be provided.
        stn_names : list of str, optional
            Station identifiers. If missing, auto-generated.
        utm_zone : str, optional
            UTM zone designator (e.g. '49N').
        """
        # override utm_zone if given
        if utm_zone:
            self.utm_zone = utm_zone

        # Case 1: load from file
        if data_fn:
            if not os.path.isfile(data_fn):
                raise SiteError(f"File not found: {data_fn}")
            profile = StationProfile(data_fn).fit()
            # extract arrays
            self.stn_names = [f'S{ii:02}' for ii in range(len(profile.easting))]
            self.easting = dict(zip(self.stn_names, profile.easting))
            self.northing = dict(zip(self.stn_names, profile.northing))
            self.elevation = dict(zip(self.stn_names, profile.elevation))
            # compute lat/lon if missing
            try:
                lat_arr, lon_arr = profile.to_latlon()
            except Exception:
                lat_arr = np.full_like(profile.easting, np.nan)
                lon_arr = np.full_like(profile.northing, np.nan)
            self.latitudes = dict(zip(self.stn_names, lat_arr))
            self.longitudes = dict(zip(self.stn_names, lon_arr))
            # compute azimuth
            try:
                azi = profile.compute_azimuth()
                # align length by dropping last station
                self.azimuth = dict(zip(self.stn_names[:-1], azi))
            except Exception:
                self.azimuth = {}
        # Case 2: load from coordinate arrays
        elif easting is not None and northing is not None:
            east = np.array(easting, dtype=float)
            north = np.array(northing, dtype=float)
            if east.shape != north.shape:
                raise SiteError("Easting and northing must match shapes.")
            n = east.size
            names = (stn_names if stn_names is not None
                     else [f'S{ii:02}' for ii in range(n)])
            if len(names) != n:
                raise SiteError("Station names length must equal number of points.")
            self.stn_names = names
            self.easting = dict(zip(names, east))
            self.northing = dict(zip(names, north))
            # optional lat/lon
            if lat is not None and lon is not None:
                lat_arr = np.array(lat, dtype=float)
                lon_arr = np.array(lon, dtype=float)
                self.latitudes = dict(zip(names, lat_arr))
                self.longitudes = dict(zip(names, lon_arr))
            else:
                # compute lat/lon via Location
                latitudes = []
                longitudes = []
                for e, n in zip(east, north):
                    loc = Location(
                        easting=e, northing=n,
                        utm_zone=self.utm_zone
                    )
                    la, lo = loc.to_latlon()
                    latitudes.append(la)
                    longitudes.append(lo)
                self.latitudes = dict(zip(names, latitudes))
                self.longitudes = dict(zip(names, longitudes))
            # elevation and azimuth default
            self.elevation = dict(zip(names, np.zeros(n)))
            # compute azimuth between stations
            az = compute_azimuth(easting=east, northing=north)
            self.azimuth = dict(zip(names[:-1], az))
        else:
            raise SiteError(
                "Insufficient data: provide data_fn or coordinate arrays.")

    def get_station(self, name: str) -> Location:
        """
        Return a Location object for the given station name.
        """
        if name not in self.stn_names:
            raise SiteError(f"Unknown station: {name}")
        loc = Location(
            latitude=self.latitudes.get(name),
            longitude=self.longitudes.get(name),
            elevation=self.elevation.get(name),
            utm_zone=self.utm_zone
        )
        # override easting/northing if available
        if name in self.easting and name in self.northing:
            loc._easting = self.easting[name]
            loc._northing = self.northing[name]
        return loc

    def as_dict(self) -> Dict[str, Dict[str, float]]:
        """
        Return all site data as nested dictionaries.
        """
        return {
            'latitude': self.latitudes,
            'longitude': self.longitudes,
            'easting': self.easting,
            'northing': self.northing,
            'elevation': self.elevation,
            'azimuth': self.azimuth
        }

    def __getitem__(self, key: str) -> Dict[str, float]:
        """
        Access site metrics by key: 'latitude','easting', etc.
        """
        data = self.as_dict().get(key)
        if data is None:
            raise SiteError(f"Invalid key: {key}")
        return data

    def to_dataframe(self) -> 'pd.DataFrame':
        """
        Return site data as a pandas DataFrame with columns:
        station, latitude, longitude, easting, northing,
        elevation, azimuth.
        """
        
        df = pd.DataFrame({
            'station': self.stn_names,
            'latitude': [self.latitudes.get(s) for s in self.stn_names],
            'longitude': [self.longitudes.get(s) for s in self.stn_names],
            'easting': [self.easting.get(s) for s in self.stn_names],
            'northing': [self.northing.get(s) for s in self.stn_names],
            'elevation': [self.elevation.get(s) for s in self.stn_names],
            'azimuth': [self.azimuth.get(s) for s in self.stn_names[:-1]] + [None]
        })
        return df

    def bounding_box(self) -> Dict[str, float]:
        """
        Compute bounding box of site in both geographic
        and UTM coordinates.

        Returns a dict with keys: lat_min, lat_max,
        lon_min, lon_max, east_min, east_max,
        north_min, north_max.
        """
        llats = np.array(list(self.latitudes.values()), float)
        llons = np.array(list(self.longitudes.values()), float)
        easts = np.array(list(self.easting.values()), float)
        norths = np.array(list(self.northing.values()), float)
        return {
            'lat_min': llats.min(), 'lat_max': llats.max(),
            'lon_min': llons.min(), 'lon_max': llons.max(),
            'east_min': easts.min(), 'east_max': easts.max(),
            'north_min': norths.min(), 'north_max': norths.max()
        }

    def filter_by_bbox(self, bbox: Dict[str, float]) -> 'Site':
        """
        Return a new Site containing only stations within
        the given bounding box. bbox requires keys: lat_min,
        lat_max, lon_min, lon_max.
        """
        mask = []
        for s in self.stn_names:
            lat = self.latitudes[s]
            lon = self.longitudes[s]
            if (bbox['lat_min'] <= lat <= bbox['lat_max'] and
                bbox['lon_min'] <= lon <= bbox['lon_max']):
                mask.append(s)
        new = Site()
        new.utm_zone = self.utm_zone
        new.stn_names = mask
        for attr in ('latitudes','longitudes','easting',
                     'northing','elevation','azimuth'):
            orig = getattr(self, attr)
            setattr(new, attr, {s: orig[s] for s in mask if s in orig})
        return new

    def nearest(self, lat: float = None, lon: float = None,
                easting: float = None, northing: float = None,
                n: int = 1) -> List[str]:
        """
        Find the n nearest stations to the given point.
        Specify either lat/lon (deg) or easting/northing (m).
        Returns list of station names.
        """
        if easting is not None and northing is not None:
            coords = np.array([self.easting[s], self.northing[s]]
                              for s in self.stn_names)
            target = np.array([easting, northing])
        elif lat is not None and lon is not None:
            # convert to UTM
            loc = Location(latitude=lat, longitude=lon)
            _, easting, northing = loc.to_utm()
            coords = np.array([self.easting[s], self.northing[s]]
                              for s in self.stn_names)
            target = np.array([easting, northing])
        else:
            raise SiteError('Provide lat/lon or easting/northing')
        dists = np.linalg.norm(coords - target, axis=1)
        idx = np.argsort(dists)[:n]
        return [self.stn_names[i] for i in idx]
    
    def __repr__(self) -> str:
        return (
            f"<Site n_stations={len(self.stn_names)} "
            f"utm_zone={self.utm_zone!r}>"
        )

Site.__doc__=r"""
Container for multiple station locations.

Parameters
----------
data_fn : str, optional
    Path to a station profile file (.stn).

    Must be provided if no coordinate arrays
    are given.
easting, northing : array-like, optional
    Parallel arrays of UTM coordinates (m).
lat, lon : array-like, optional
    Parallel arrays of geographic coordinates
    (degrees).
stn_names : list of str, optional
    Station identifiers; generated automatically
    if omitted.
utm_zone : str, optional
    UTM zone designator (e.g. '49N').

Attributes
----------
stn_names : List[str]
    Ordered list of station identifiers.
latitudes : Dict[str, float]
    Mapping station → latitude (deg).
longitudes : Dict[str, float]
    Mapping station → longitude (deg).
easting : Dict[str, float]
    Mapping station → UTM easting (m).
northing : Dict[str, float]
    Mapping station → UTM northing (m).
elevation : Dict[str, float]
    Mapping station → elevation (m).
azimuth : Dict[str, float]
    Mapping station → azimuth (deg).
utm_zone : Optional[str]
    Shared UTM zone for all stations.

Notes
-----
- Internally uses StationProfile to parse .stn files.
- If both file and arrays are supplied, file takes
  precedence.
- Missing values are filled with NaN or zeros.

Examples
--------
>>> site = Site(data_fn='station.stn')
>>> df = site.to_dataframe()
>>> bbox = site.bounding_box()
>>> nearby = site.nearest(lat=34.05,
...                      lon=-118.25, n=3)
"""