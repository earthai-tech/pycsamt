# -*- coding: utf-8 -*-
# Author: L. Kouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt.site.profile

Station: handles reading and transforming station profile files
(e.g. Zonge .stn or SEG-EDI station files), computes separations,
angles, dipole lengths, etc.

Topography : 
    
"""
import os
import time
import shutil
from typing import Optional, Tuple, List 

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

from ..api.typing import Optional, ArrayLike, F
from ..exceptions import NotFittedError, PycsamtError, TopographyError 
from ..utils.io import parse_stn_profile
from ..utils.func_utils import compute_azimuth, stn_separation
from ..utils.mathext import round_dipole_length
from .location import Location




# XXX TO DO : 
# -*- coding: utf-8 -*-
#   Licence:BSD-3-Clause
#   Author: LKouadio <etanoyau@gmail.com>
"""
Manage site data.
"""
from __future__ import annotations 
import os
import shutil
import copy 
import re

import warnings 
import numpy as np 

from ._watexlog import watexlog 

from .exceptions import ( 
    SiteError, 
    ProfileError,
    NotFittedError, 
    )
from .property import ( 
    UTM_DESIGNATOR
    )
from .utils.coreutils import makeCoords 
from .utils.exmath import ( 
    get_bearing, 
    scalePosition, 
    )

from .utils.funcutils import (
 
    interpolate_grid,
    reshape,
    
    ) 
from ..gis.utils import ( 
    assert_elevation_value, 
    assert_lat_value, 
    assert_lon_value, 
    ll_to_utm, 
    utm_to_ll, 
    project_point_ll2utm, 
    project_point_utm2ll, 
    assert_xy_coordinate_system, 
    convert_position_float2str, 
    convert_position_str2float, 
    )
# from .utils.gistools import (
#     assert_elevation_value, 
#     assert_lat_value, 
#     assert_lon_value , 
#     ll_to_utm, 
#     utm_to_ll, 
#     project_point_ll2utm, 
#     project_point_utm2ll,
#     assert_xy_coordinate_system,
#     convert_position_str2float, 
#     convert_position_float2str, 
#     )
from ..utils.arrayops import assert_xy_in, is_iterable , frameify
from ..utils.validation import ( 
    _check_consistency_size , 
    _is_numeric_dtype , 
    _assert_all_types,
    _validate_name_in
) 

__all__= ['Location', 'Profile']


__all__ = [
    'Station',
    'Topography', 
    'Profile'
]

class Station: # StationProfile renamed to Station
    """High-level interface for station profile data."""
    def __init__(self, profile_fn: str = None):
        self.profile_fn = profile_fn
        self._fitted = False

        # data containers
        self.positions = None
        self.easting = None
        self.northing = None
        self.latitude = None
        self.longitude = None
        self.elevation = None
        self.azimuth = None
        self.dipole_length = None
        self.profile_angle = None
        self.geoelectric_strike = None
        self.utm_zone = None

    def fit(self, profile_fn: str = None):
        """
        Read and parse the station profile file.
        Populates easting, northing, elevation, positions, utm_zone.
        """
        fn = profile_fn or self.profile_fn
        if fn is None:
            raise ValueError("No profile filename provided to fit().")
        if not os.path.isfile(fn):
            raise PycsamtError(f"File not found: {fn}")

        data = parse_stn_profile(fn)
        self.positions = np.asarray(data.get('position', []), dtype=float)
        self.easting = np.asarray(data.get('easting', []), dtype=float)
        self.northing = np.asarray(data.get('northing', []), dtype=float)
        self.elevation = np.asarray(data.get('elevation', []), dtype=float)
        self.utm_zone = data.get('utm_zone')

        # convert UTM to lat/lon
        loc = Location(utm_zone=self.utm_zone)
        # assign protected attributes for conversion
        loc._easting = self.easting
        loc._northing = self.northing
        lat, lon = loc.to_latlon()
        self.latitude = lat
        self.longitude = lon

        self._fitted = True
        return self

    def _check_fitted(self):
        if not self._fitted:
            raise NotFittedError(
                "StationProfile must be fitted before use.")

    def compute_azimuth(self, extrapolate: bool = True):
        self._check_fitted()
        self.azimuth = compute_azimuth(
            easting=self.easting,
            northing=self.northing,
            utm_zone=self.utm_zone,
            extrapolate=extrapolate
        )
        return self.azimuth

    def compute_separation(self, interpolate: bool = False):
        self._check_fitted()
        sep, _ = stn_separation(
            easting=self.easting,
            northing=self.northing,
            interpolate=interpolate
        )
        self.positions = sep
        return self.positions

    def compute_dipole_length(self):
        self._check_fitted()
        if self.positions is None:
            self.compute_separation()
        self.dipole_length = round_dipole_length(self.positions.mean())
        return self.dipole_length

    def compute_profile_angle(self, easting=None, northing=None):
        """
        Compute the geoelectric profile angle and strike.
        """
        self._check_fitted()
        if easting is not None:
            self.easting = np.asarray(easting, dtype=float)
        if northing is not None:
            self.northing = np.asarray(northing, dtype=float)
        if self.easting is None or self.northing is None:
            raise PycsamtError(
                "Easting and northing must be set before computing profile angle.")

        r1 = linregress(self.easting, self.northing)
        r2 = linregress(self.northing, self.easting)
        if r2.stderr < r1.stderr:
            m = 1.0 / r2.slope
        else:
            m = r1.slope

        angle = (90.0 - np.degrees(np.arctan(m))) % 180.0
        self.profile_angle = np.round(angle, 2)
        strike = (self.profile_angle + 90.0) % 180.0
        self.geoelectric_strike = np.floor(strike)
        return self.profile_angle, self.geoelectric_strike
    
    @staticmethod
    def adjust_coordinates(
        x: float = 0.0,
        y: float = 0.0,
        station_pk=None,
        easting=None,
        northing=None,
        elevation=None,
        rewrite: bool = False,
        savepath: str = None,
        profile_fn: str = None
    ):
        """
        Adjust station profile coordinates by offsets x (easting) and y (northing).
        Optionally rewrite a new .stn file if rewrite=True.

        Returns:
            station_pk, easting, northing, elevation arrays.
        """
        # read existing data if filename provided
        if profile_fn:
            data = parse_stn_profile(profile_fn)
            station_pk = np.asarray(data.get('position', []), dtype=float)
            easting = np.asarray(data.get('easting', []), dtype=float)
            northing = np.asarray(data.get('northing', []), dtype=float)
            elevation = np.asarray(data.get('elevation', []), dtype=float)
        else:
            station_pk = np.asarray(station_pk or [], dtype=float)
            easting = np.asarray(easting or [], dtype=float)
            northing = np.asarray(northing or [], dtype=float)
            elevation = np.asarray(elevation or [], dtype=float)

        # apply offsets
        easting = easting + float(x)
        northing = northing + float(y)

        if rewrite and profile_fn:
            # regenerate .stn file
            lines = ['"""dot""","""e""","""n""","""h"""\n']
            for pk, ex, no, el in zip(station_pk, easting, northing, elevation):
                lines.append(f"{int(pk):<7},{ex:12.3f},{no:12.3f},{el:8.2f}\n")
            base = os.path.splitext(os.path.basename(profile_fn))[0]
            newfn = f"{base}_adj.stn"
            with open(newfn, 'w') as f:
                f.writelines(lines)
            if savepath:
                shutil.move(newfn, os.path.join(savepath, newfn))
            print(f"Rewritten adjusted station file: {newfn}")

        return station_pk, easting, northing, elevation
    
    def straighten_profile(
        self,
        X=None,
        Y=None,
        mode: str = 'classic',
        adjust_offset: tuple = (0.0, 0.0),
        rewrite: bool = False,
        savepath: str = None
    ):
        """
        Straighten or redistribute station coordinates along a profile line.

        Parameters:
            X, Y           arrays of original easting/northing
            mode           'classic', 'equidistant', or 'natural'
            adjust_offset  (dx, dy) offsets applied before straightening
            rewrite        if True, write a new .stn file
            savepath       directory to save rewritten file
        Returns:
            new_easting, new_northing arrays
        """
        self._check_fitted()
        if X is not None: 
            self.easting = np.asarray(X, dtype=float)
        if Y is not None: 
            self.northing = np.asarray(Y, dtype=float)
        if self.easting is None or self.northing is None:
            raise PycsamtError("Coordinates must be set to straighten profile.")

        n = self.easting.size
        # apply initial offset
        dx, dy = map(float, adjust_offset)
        self.easting += dx
        self.northing += dy

        # compute reference dipole length if needed
        if self.dipole_length is None:
            self.compute_dipole_length()

        # store originals
        orig_e = self.easting.copy()
        orig_n = self.northing.copy()

        if mode.startswith('classic'):
            # linear interpolation between endpoints
            idx = np.arange(n) # noqa
            self.easting = np.linspace(orig_e[0], orig_e[-1], n)
            self.northing = np.linspace(orig_n[0], orig_n[-1], n)
            
        elif mode.startswith('equ'):
            # equal-distance spacing using dipole_length
            total = (n - 1) * self.dipole_length
            dist = np.linspace(0, total, n)
            # vector from start to end
            ve = orig_e[-1] - orig_e[0]
            vn = orig_n[-1] - orig_n[0]
            norm = np.hypot(ve, vn)
            dir_e = ve / norm
            dir_n = vn / norm
            self.easting = orig_e[0] + dir_e * dist
            self.northing = orig_n[0] + dir_n * dist
        else:
            # natural: maintain step distances but align direction
            steps_e = np.diff(orig_e)
            steps_n = np.diff(orig_n)
            avg_e = steps_e.mean()
            avg_n = steps_n.mean()
            new_e = [orig_e[0]]
            new_n = [orig_n[0]]
            for i in range(1, n):
                new_e.append(new_e[-1] + avg_e)
                new_n.append(new_n[-1] + avg_n)
            self.easting = np.array(new_e)
            self.northing = np.array(new_n)

        # optionally rewrite file
        if rewrite and self.profile_fn:
            self.adjust_coordinates(
                x=0, y=0,
                station_pk=self.positions,
                easting=self.easting,
                northing=self.northing,
                elevation=self.elevation,
                rewrite=True,
                savepath=savepath,
                profile_fn=self.profile_fn
            )

        return self.easting, self.northing
    
    def rewrite_profile(
        self,
        easting=None,
        northing=None,
        elevation=None,
        area_name=None,
        username=None,
        add_azimuth=False,
        **kwargs
    ):
        """
        Write a new station profile .stn file based on current coordinates.
        """
        # pull kwargs
        utm_zone = kwargs.pop('UTM_Zone', None)
        dipole_length = kwargs.pop('dipole_length', None)
        output_name = kwargs.pop('output_name', 'new_profile')
        savepath = kwargs.pop('savepath', None)

        if savepath:
            self.savepath = savepath
        if utm_zone:
            self.utm_zone = utm_zone
        if area_name is None:
            area_name = ''
        if username is None:
            username = ''
        if easting is not None:
            self.easting = np.asarray(easting, dtype=float)
        if northing is not None:
            self.northing = np.asarray(northing, dtype=float)
        if elevation is not None:
            self.elevation = np.asarray(elevation, dtype=float)

        if self.easting is None or self.northing is None:
            raise PycsamtError("Easting and northing required to rewrite profile.")

        # header
        lines = []
        lines.append(f">LOCATION: {area_name}\n")
        lines.append(f">USERNAME: {username}\n")
        lines.append(f">DATE: {time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())}\n")
        lines.append(f">UTM_ZONE: WGS84-{self.utm_zone}\n")
        lines.append(f">SOFTWARE: pyCSAMT\n") # noqa

        # compute azimuth if requested
        if add_azimuth:
            if self.azimuth is None:
                self.compute_azimuth()
            az = self.azimuth
        else:
            az = None

        # station positions
        if dipole_length is None:
            self.compute_dipole_length()
        stn_pos = np.arange(
            0, self.easting.size * self.dipole_length, self.dipole_length)

        # table header
        header = ['""sta""', '""len(m)""', '""east(m)""',
                  '""north(m)""', '""lat(°)""', '""long(°)""', '""elev(m)""']
        if add_azimuth:
            header.append('""azim(°)""')
        lines.append(','.join(header) + '\n')

        # rows
        loc = Location(utm_zone=self.utm_zone)
        for i, pk in enumerate(stn_pos):
            loc._easting = self.easting[i]
            loc._northing = self.northing[i]
            lat, lon = loc.to_latlon()
            row = [f"S{i:02}", f"{pk:.0f}",
                   f"{self.easting[i]:.3f}", f"{self.northing[i]:.3f}",
                   f"{lat:.6f}", f"{lon:.6f}", f"{self.elevation[i]:.2f}"]
            if add_azimuth:
                row.append(f"{az[i]:.2f}")
            lines.append(','.join(row) + '\n')

        # write file
        fn_out = f"{output_name}.stn"
        with open(fn_out, 'w', encoding='utf8') as f:
            f.writelines(lines)
        if self.savepath:
            shutil.move(fn_out, os.path.join(self.savepath, fn_out))
        return fn_out
    
    def to_latlon(self, utm_zone: str = None):
        self._check_fitted()
        loc = Location(utm_zone=utm_zone or self.utm_zone)
        loc._easting = self.easting
        loc._northing = self.northing
        lat, lon = loc.to_latlon()
        self.latitude = lat
        self.longitude = lon
        return lat, lon

    def to_utm(self, utm_zone: str = None):
        self._check_fitted()
        loc = Location(
            latitude=self.latitude,
            longitude=self.longitude
        )
        zone, east, north = loc.to_utm()
        self.utm_zone = zone
        self.easting = east
        self.northing = north
        return self.easting, self.northing, self.utm_zone
    
    def plot(
        self,
        kind: str = 'plan',
        ax=None,
        figsize: tuple = (8, 6),
        show: bool = True,
        **kwargs
    ):
        """
        Plot the station profile.

        Parameters
        ----------
        kind : {'plan','elevation','line'}, default 'plan'
            'plan'      : top-down view (easting vs northing)
            'elevation' : profile along distance vs elevation
            'line'      : easting & northing vs station index
        ax : matplotlib.axes.Axes, optional
            Axis to plot on; created if None.
        figsize : tuple, default (8,6)
            Figure size if ax is None.
        show : bool, default True
            Call plt.show() if True.
        **kwargs : passed to plotting function

        Returns
        -------
        ax : matplotlib.axes.Axes
        """
        self._check_fitted()
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        # plan view
        if kind == 'plan':
            ax.scatter(
                self.easting,
                self.northing,
                c=self.elevation,
                cmap=kwargs.pop('cmap', 'viridis'),
                **kwargs
            )
            ax.set_xlabel('Easting (m)')
            ax.set_ylabel('Northing (m)')
            ax.set_title('Plan View of Station Profile')
        # elevation vs distance
        elif kind == 'elevation':
            if self.positions is None:
                self.compute_separation()
            dist = np.insert(np.cumsum(self.positions), 0, 0)
            ax.plot(
                dist,
                self.elevation,
                **kwargs
            )
            ax.set_xlabel('Distance Along Profile (m)')
            ax.set_ylabel('Elevation (m)')
            ax.set_title('Elevation Profile')
        # line plot of coords
        elif kind == 'line':
            idx = np.arange(len(self.easting))
            ax.plot(idx, self.easting, label='Easting', **kwargs)
            ax.plot(idx, self.northing, label='Northing', **kwargs)
            ax.set_xlabel('Station Index')
            ax.set_ylabel('Coordinate (m)')
            ax.set_title('Easting & Northing vs Station Index')
            ax.legend()
        else:
            raise ValueError(f"Unknown plot kind: {kind}")
        if show:
            plt.show()
        return ax
    
    def __repr__(self):
        cls = self.__class__.__name__
        fn = os.path.basename(self.profile_fn) if self.profile_fn else 'None'
        return f"<{cls} fn={fn!r} fitted={self._fitted}>"


Station.__doc__=r"""
Handles reading and transforming station profile files
(e.g. Zonge .stn or SEG-EDI station files), computes separations,
angles, dipole lengths, etc.

Parameters
----------
profile_fn : str, optional
    Path to station profile file (.stn). If not provided,
    must be passed to :meth:`fit`.

Attributes
----------
positions : ndarray, shape (n,)
    Station separations between dipoles (m).
easting : ndarray, shape (n,)
    UTM easting coordinates (m).
northing : ndarray, shape (n,)
    UTM northing coordinates (m).
latitude : ndarray, shape (n,)
    Geographic latitude (degrees).
longitude : ndarray, shape (n,)
    Geographic longitude (degrees).
elevation : ndarray, shape (n,)
    Elevation at each station (m).
azimuth : ndarray, shape (n-1,)
    Azimuth angles between consecutive stations (degrees).
dipole_length : float
    Averaged dipole length used for spacing (m).
utm_zone : str
    UTM zone designator (e.g. '49N').

Notes
-----
- Call :meth:`fit` before computing attributes.
- Methods raise :exc:`NotFittedError` if not yet fitted.

Examples
--------
>>> from pycsamt.site.profile import StationProfile
>>> profile = StationProfile('data/K1.stn')
>>> profile.fit()
>>> profile.compute_separation()
array([ 36.02,  35....])
>>> profile.compute_azimuth()
array([123.45, ...])
>>> profile.compute_dipole_length()
50.0
>>> profile.to_latlon()
(array([...]), array([...])))
>>> profile.to_utm()
('49N', array([...]), array([...])))
"""


class Topography:
    """
    High-level handler for .bln or CSV elevation files.

    Parameters
    ----------
    path_bl : str, optional
        Path to a .bln or comma-delimited file with columns:
        offset, northing, elevation, [station_id].
    step : float, optional
        Station spacing (same units as offsets);
        if None, inferred from file.
    encoding : str, default 'utf-8'
        File encoding.
    """
    def __init__(
        self,
        path_bl: Optional[str] = None,
        step: Optional[float] = None,
        encoding: str = 'utf-8'
    ):
        self.path_bl = path_bl
        self.encoding = encoding
        self.step = step

        # parsed arrays
        self._offset: Optional[np.ndarray] = None  # shape (n,)
        self._north: Optional[np.ndarray] = None   # shape (n,)
        self._elev: Optional[np.ndarray] = None    # shape (n,)
        self._stations: Optional[List[str]] = None # optional station IDs
        # computed eastings
        self._east: Optional[np.ndarray] = None    # shape (n,)
        self._fitted: bool = False

        if path_bl:
            self.fit(path_bl, step)


    def fit(
        self,
        path_bl: Optional[str] = None,
        step: Optional[float] = None
    ) -> 'Topography':
        """
        Read and parse the topography file (.bln, .csv, or .txt).

        Parameters
        ----------
        path_bl : str, optional
            Override stored path to file (bln/csv/txt).
        step : float, optional
            Override station spacing.

        Returns
        -------
        self
        """
        fn = path_bl or self.path_bl
        if fn is None or not os.path.isfile(fn):
            raise TopographyError(f"File not found: {fn}")

        ext = os.path.splitext(fn)[1].lower()
        data = []
        with open(fn, 'r', encoding=self.encoding) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                # normalize delimiters for bln, csv, txt
                if ext in ['.bln', '.csv', '.txt']:
                    parts = [s for s in line.replace(',', ' ').split() if s]
                else:
                    parts = line.split()
                # expect at least offset, northing, elevation
                if len(parts) < 3:
                    continue
                try:
                    off = float(parts[0])
                    nor = float(parts[1])
                    elv = float(parts[2])
                except ValueError:
                    continue
                # optional station ID
                sid = parts[3] if len(parts) > 3 else None
                data.append((off, nor, elv, sid))

        if not data:
            raise TopographyError("No valid data found in file.")

        arr = np.array([d[:3] for d in data], dtype=float).T
        self._offset = arr[0]
        self._north = arr[1]
        self._elev = arr[2]
        self._stations = [d[3] for d in data if d[3] is not None]

        # infer step if not provided
        self.step = step if step is not None else _default_step_from_file(
            self._offset.size,
            float(self._offset[-1] - self._offset[0])
        )
        # treat east = relative offset
        self._east = self._offset.copy()

        self._fitted = True
        return self

    def _check_fitted(self):
        if not self._fitted:
            raise TopographyError("Topography.fit() must be called first.")
            
    def profile(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return raw offset and elevation arrays.

        Returns
        -------
        offset : ndarray, shape (n,)
        elevation : ndarray, shape (n,)
        """
        self._check_fitted()
        return self._offset.copy(), self._elev.copy()
    
    def interpolate(
        self,
        new_step: float,
        kind: str = 'linear'
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Resample profile at a new spacing.

        Parameters
        ----------
        new_step : float
            Desired spacing between points.
        kind : str, default 'linear'
            Interpolation kind for np.interp.

        Returns
        -------
        new_offset, new_elev : ndarrays
        """
        self._check_fitted()
        total = self._offset[-1] - self._offset[0]
        npt = int(np.floor(total / new_step)) + 1
        offs = np.linspace(self._offset[0], self._offset[-1], npt)
        elev_new = np.interp(offs, self._offset, self._elev)
        return offs, elev_new
    
    def model(self, utm_zone: str) -> np.ndarray:
        """
        Build full elevation model (Easting, Northing, Elevation).

        Parameters
        ----------
        utm_zone : str
            UTM zone designator for conversion.

        Returns
        -------
        model : ndarray, shape (3, n)
            Rows: [easting, northing, elevation]
        """
        self._check_fitted()
        # convert each relative point to absolute UTM
        eastings, northings = [], []
        for off, nor in zip(self._offset, self._north):
            loc = Location(utm_zone=utm_zone)
            loc._easting = off
            loc._northing = nor
            lat, lon = loc.to_latlon()
            _, ea, no = loc.to_utm()
            eastings.append(ea)
            northings.append(no)
        return np.vstack((np.array(eastings),
                          np.array(northings),
                          self._elev.copy()))
    def plot_profile(
        self,
        ax=None,
        show: bool = True,
        **plt_kwargs
    ) -> plt.Axes:
        """
        Plot offset vs elevation profile.

        Parameters
        ----------
        ax : matplotlib Axes, optional
            Axes to draw on. Creates one if None.
        show : bool, default True
            Call plt.show() if True.
        **plt_kwargs
            Passed to ax.plot().

        Returns
        -------
        ax : matplotlib Axes
        """
        self._check_fitted()
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(self._offset, self._elev, **plt_kwargs)
        ax.set_xlabel('Offset')
        ax.set_ylabel('Elevation')
        ax.set_title('Elevation Profile')
        if show:
            plt.show()
        return ax

    def plot_model(
        self,
        utm_zone: str,
        ax=None,
        show: bool = True,
        **plt_kwargs
    ) -> plt.Axes:
        """
        Plot UTM model easting vs northing, color-coded by elevation.

        Parameters
        ----------
        utm_zone : str
            UTM zone for conversion.
        ax : matplotlib Axes, optional
        show : bool, default True
        **plt_kwargs
            Passed to ax.scatter().

        Returns
        -------
        ax : matplotlib Axes
        """
        self._check_fitted()
        model = self.model(utm_zone)
        if ax is None:
            fig, ax = plt.subplots()
        sc = ax.scatter(model[0], model[1], c=model[2], **plt_kwargs)
        plt.colorbar(sc, ax=ax, label='Elevation')
        ax.set_xlabel('Easting')
        ax.set_ylabel('Northing')
        ax.set_title('Topography Model')
        if show:
            plt.show()
        return ax
    
    @property
    def elevation_profile(self) -> np.ndarray:
        """
        Elevation profile as 2×n array: [offsets; elevations].
        """
        self._check_fitted()
        return np.vstack((self._offset, self._elev))

    @property
    def elevation_model(self) -> np.ndarray:
        """
        Full model as 3×n array: [easting; northing; elevation].
        """
        self._check_fitted()
        return np.vstack((self._east, self._north, self._elev))

    @property
    def n_points(self) -> int:
        """Number of elevation points."""
        self._check_fitted()
        return self._offset.size

    def summary(self) -> dict:
        """
        Basic statistics.

        Returns
        -------
        info : dict
            Keys: n_points, min_elev, max_elev, total_length
        """
        self._check_fitted()
        return {
            'n_points': self.n_points,
            'min_elev': float(self._elev.min()),
            'max_elev': float(self._elev.max()),
            'total_length': float(self._offset[-1] - self._offset[0])
        }

    def __repr__(self) -> str:
        fn = os.path.basename(self.path_bl) if self.path_bl else 'None'
        pts = self.n_points if self._fitted else 0
        return f"<Topography fn={fn!r} points={pts}>"


def _default_step_from_file(n_points: int, total_length: float) -> float:
    """
    Infer step size from number of points and total length.
    """
    return total_length / (n_points - 1) if n_points > 1 else 0.0



class Profile: 
    
    def __init__(
        self , 
        *, 
        utm_zone =None, 
        coordinate_system=None, 
        datum = 'WGS84', 
        epsg = None, 
        reference_ellipsoid =23, 
        ): 
        self.epsg=epsg 
        self.datum =datum 
        self.reference_ellipsoid= reference_ellipsoid 
        self.utm_zone=utm_zone
        self.coordinate_system= coordinate_system

    def fit(
        self, 
        x:ArrayLike |str , 
        y:ArrayLike | str ,
        elev: Optional [ArrayLike] =None, 
        **fit_params
        ): 
        """ Populate profile attributes from x and y. 
        
        By default if the coordinate system is given as latitude/longitude 
        x, y are latitude and longitude respectively. 
 
        Parameters
        ------------
        x, y: ArrayLike 1d /str  
           One dimensional arrays. `x` can be consider as the abscissa of the  
           landmark and `y` as ordinates array.  If `x` or `y` is passed as 
           string argument, `data` must be passed as `fit_params` keyword 
           arguments and the name of `x` and `y` must be a column name of 
           the `data`. 
           By default `x` and `y` are considered as `longitude` and `latitude` 
           when ``dms`` or ``ll`` coordinate system is passed. 
  
        elev: ArrayLike 1d/str 
            Arraylike 1d of elevation at each positions. If not supplied 
            should be set to null at each points. If given, it must be  
            consistent with `x` and `y`. 
            
        data: pd.DataFrame 
           Data containing `x` and `y` values as series. Then if `x` and `y`
           are given as string argument, their names must be included in the 
           data columns. Otherwise an error will raise. 
        
        Return 
        ---------
        self :`watex.site.Profile` 
           Object of Profile. 
           
        Note 
        ------
        When `data` is supplied and `x` and `y` are given by their names 
        existing in the dataframe columns, by default, the non-numerical 
        data are removed. However, if `y` and `x` are given in DD:MM:SS in 
        the dataframe, the coordinate system must explicitly set to ``dms`
        to keep the non-numerical values  in the data. 
        
        """
        #15921 epsg io
        emsg =("'DD:MM:SS.ms' coordinate system is detected.")
        data = fit_params.pop('data', None) 
  
        self.coordinate_system = str(
            self.coordinate_system).lower().strip() 
        if  ( 
                self.coordinate_system.find ('dms')>=0 
            or 'dd:mm' in self.coordinate_system ) : 
            self.coordinate_system='dms'
    
        if data is not None: 
            # when coordinate_system is explicity passed 
            # and data is given, don't suppress the non-numerical 
            # data. Keep it untouch and try to convert x, y to 
            # float type later. 
            suppress_cf = True 
            if self.coordinate_system=='dms': 
               suppress_cf= False 
            # suppress non numerical values 
            data = frameify(data, pop_cat_features= suppress_cf)
            # get elevation if string is given 
            if isinstance (elev , str): 
                elev = elev if elev in data.columns else None 
                if elev is None: 
                    warnings.warn("Elevation not found in the dataframe")
                if elev is not None: 
                    elev = np.array (data[elev]) 
        
        self.x, self.y = assert_xy_in( x , y, data = data )
        
        # set x and y names useful for 
        # plotting labels. 
        self.xname_=None ; self.yname_=None 
        if ( 
                isinstance ( x, str ) 
                or hasattr (x, "name" )
                ) : 
            self.xname_= x if isinstance ( x, str ) else x.name 
        if ( 
                isinstance ( y, str ) 
                or hasattr (y, "name" )
                ) : 
           self.yname_= y if isinstance ( y, str ) else y.name 
        
        if self.coordinate_system in ('none', 'auto'): 
            self.coordinate_system = assert_xy_coordinate_system (
                self.x, self.y )
        
        if self.coordinate_system =='dms' : 
            try : 
                self.x, self.y = Profile.dms2ll(self.x, self.y) 
            except ValueError as e: 
                raise ValueError( emsg + str(e))
            else: 
                # initialize the system to ll if 
                # conversion succeeded. 
                self.coordinate_system =='ll' 
        # set 0. to elevation 
        self.elev = elev if elev is not None else np.zeros_like(
            self.x, dtype = np.float32 ) 
        self.elev = np.array (self.elev , dtype = np.float32 )
        
        if not _check_consistency_size (self.x , np.array( self.elev),
                                        error ='ignore'  ): 
            raise ProfileError ("Elevation and x or y must be consistent."
                              f" Got {len(x)} and {len(elev)}") 
        return self 
        
    def bearing (self, *, to_degree:bool = True ): 
        """
        Compute the bearing between two coordinates.
        
        A bearing is a direction of one point relative to another point, 
        usually given as an angle measured clockwise from north. 
        In navigation, bearings are often used to determine the direction 
        to a destination or to plot a course on a map. There are two main 
        types of bearings: absolute bearing and relative bearing.
        
        - Relative bearing 
          refers to the angle between the forward direction of a craft 
          (heading) and the location of another object. 
        - Absolute bearing 
          refers to the angle between the magnetic north (magnetic bearing) 
          or true north (true bearing) and an object
          
        The  bearing (:math:`\beta`) between two coordinates points 1(lon1, lat1) 
        and 2 (lon2, lat2) can becalculated as:
        
        .. math:: 
            
            \beta = atan2(sin(lon2-lon1)*cos(lat2), cos(lat1)*sin(lat2) – \
                          sin(lat1)*cos(lat2)*cos(lon2-lon1))
                
        By default, the first and last coordinates for points 1 and 2 are 
        used as `latlon1` and `latlon2` respectively. 
        
        Parameters 
        ------------
        to_degree: bool, default=True 
          Convert the bearing from radians to degree. 
          
        Returns 
        --------
        :math:`\beta` : float, 
           The value of bearing in degree ( default). 
           
        """
        msg =(
            "x, y are not in longitude/latitude format"
            " while 'utm_zone' is not set. The bearing"
            " should be less accurate. Provide the UTM"
            " zone to improve the accuracy.")
        
        self.inspect 
        
        if self.coordinate_system =='ll': 
            xs = np.array(copy.deepcopy(self.x)) 
            ys = np.array(copy.deepcopy(self.y))
        
        if ( 
                self.coordinate_system !='ll' 
                and self.utm_zone is None) : 
            warnings.warn(msg ) 
            
        self.utm_zone = self.utm_zone or '49R'    
        if self.coordinate_system !='ll': 
            # recompute value to lat/lon 
            # from easting/northing
            ys , xs = Location.to_latlon_in(
                self.x, self.y, utm_zone= self.utm_zone) 
        # compute bearing.     
        return get_bearing((ys[0], xs[0]) , ( ys[-1], xs[-1] ), 
                           to_deg= to_degree 
                           ) 
    
    def distance (self, *, average_distance: bool = True): 
        """Compute the distance between profile coordinates points.
        
        Preferably, it uses the UTM coordinates positions. By default 
        the coordinate system is automatically detected. 
        
        Parameters 
        -----------
        average_distance: bool, default =True, 
           Returns the average value of the distance between different points. 

        Returns 
        ---------
        d: Arraylike of shape (N-1) 
          Is the distance between points or the average of all distances.  
        
        See Also 
        ---------
        watex.utils.exmath.get_distance: 
            Get the distance from longitude/latitude or easting/northing 
            coordinates. 
            
        Examples 
        --------- 
        >>> import numpy as np 
        >>> from watex.site import Profile 
        >>> posx = np.random.rand (7) *10 
        >>> posx = np.abs ( np.random.randn (7) * 12 ) 
        >>> po= Profile().fit(posx, posy )
        >>> # convert data to UTM and compute distance becuase 
        >>> # our toy example has value less than 90 and 180.
        >>> po.distance () 
        251053.3287093233
        >>> po.coordinate_system = 'UTM' 
        >>> 6.451210308544236
        """
        self.inspect 
        
        x, y = self.x , self.y 
        
        if self.coordinate_system =='ll': 
            x , y = Location.to_utm_in(
                lats= self.y, lons =self.y, 
                epsg = self.epsg , 
                datum = self.datum , 
                reference_ellipsoid= self.reference_ellipsoid
            ) 
    
        d = np.sqrt( np.diff (x) **2 + np.diff (y)**2 ) 
        
        return d.mean()  if average_distance else d 
    
    def scale_positions (self, **sp_kws): 
        """
        Scale the position coordinates along x and y.
        
        Parameters
        -----------
        sp_kws: dict 
           Keyword arguments passed to :func:`~watex.utils.scalePosition`. 
           
        Returns 
        --------
        x, y : Arraylikes 
           Scaled positions 
           
        See also 
        ---------
        watex.utils.scalePositions: 
            Scale positions using the `scipy` curve fit. 
        watex.utils.exmath.scale_positions: 
            Scale and shift positions using bearing. 
            
        """
        self.inspect 
        x, *_ = scalePosition(self.x , **sp_kws) 
        y, *_ = scalePosition(self.y , **sp_kws)
        
        return x, y 
    
    def shift_positions (
        self, 
        *, 
        step:float= None, 
        use_average_dist:bool=False,
        angle:Optional[float]=None,
        ): 
        """
        Shift the x and y  position coordinates from the step and angle. 
         
        By default, it consider `x` and `y` as easting/latitude and 
        northing/longitude coordinates respectively. It latitude and longitude 
        are given, specify the parameter `is_latlon` to ``True``. 
        
        Parameters
        ------------
        step: float, Optional 
           The positions separation. If not given, the average distance between 
           all positions should be used instead. 
        use_average_dist: bool, default=False, 
           Use the distance computed between positions for the correction. 
           
        angle: float, Optional 
           Bearing angle in degree to shift the profile line . If ``None``, 
           the ``bearing`` of `x` and `y` is used instead.
           
        Returns 
        --------
        xx, yy: Arraylike 1d, 
           The arrays of position correction from `x` and `y` using the 
           bearing. 
           
        See Also 
        ---------
        watex.utils.exmath.get_bearing: 
            Compute the  direction of one point relative to another point. 
          
        Examples
        ---------
        >>> from watex.utils.exmath import scale_positions 
        >>> east = [336698.731, 336714.574, 336730.305] 
        >>> north = [3143970.128, 3143957.934, 3143945.76]
        >>> p= Profile().fit (east, north)
        >>> east_c , north_c= p.scale_positions (east, north, step =20) 
        >>> east_c , north_c
        (array([336686.69198337, 336702.53498337, 336718.26598337]),
         array([3143986.09866306, 3143973.90466306, 3143961.73066306]))
        """
        cs ={'ll': 'longitude/latitude', 
             'dms': 'dd:mm:ss'}
        self.inspect 
        
        if self.coordinate_system !='utm': 
            cse = cs.get('ll').title () if self.coordinate_system =='ll' else (
                cs.get('dms').upper()  if self.coordinate_system =='dms' else
              self.coordinate_system .upper() )
            warnings.warn((
                "{0!r} coordinates system is detected. Shifting positions" 
                " expects UTM coordinates. It is recommended to convert"
                " positions data to UTM (ref:`watex.site.to_utm_in`) before"
                " processing. The following with {0} coordinates"
                " might lead to invalid results. Use at your own risk." 
                ).format(cse)
                          )
            
        if ( not use_average_dist
            and step is None 
            ): 
            warnings.warn("Step is not given. The mean-distance of"
                          " positions should be used instead.")
            use_average_dist =True 
            
        if use_average_dist: 
            step = self.distance (
                average_distance= use_average_dist) 
 
        step  = float (_assert_all_types(step, float, int, objname ='Step'))

        if angle is not None: 
            angle  = float (_assert_all_types(
                step, float, int, objname ='Bearing angle (in degree)'))
            
            angle = np.deg2rad (angle %360)  
            
        angle = angle or self.bearing(to_degree =False ) # return bearing in rad.
     
        xx = self.x + ( step * np.cos (angle))
        yy = self.y +  (step  * np.sin(angle))

        return xx, yy 
    
    @staticmethod 
    def f_ (ar , /,  func: str | F  = 'dms->ll'): 
        """
        Converter position function from dms to longitude/latitude degree 
        decimal or vice versa.
        
        Convert position from str (DD:MM:SS) to float (latitude/longitude)
        and vice versa.
        
        Parameters 
        -----------
        ar: ArrayLike 1d, 
           Array containing the position coordinates for conversion. 
        func: Callable or str, default ='dms->ll'
            Converter functions. They can be: 
               
            - :func:`~watex.utils.gistools.convert_position_str2float` 
              for ``:dd:mm:ss``  to foat(long, lat) coordinates. If 
              string is passed it should be ['dms2ll'|'dmstoll'|'dms-<ll']. 
            - :func:`~watex.utils.gistools.convert_position_float2str` 
              from float (long, lat) in decimal degree coordinates 
              to ``dd:mm:ss``. When string is passed, it should be 
              ['ll2dms'|'lltodms'|'ll->dms']
               
        Returns
        --------
        generator obj. 
           Map object composed of value converted. 
        
        """
        if isinstance (func, str): 
            f = _validate_name_in(func, defaults = (
                'll2dms', 'lltodms', 'll->dms'), 
                expect_name= convert_position_float2str 
                ) 
        if not callable (f): 
            f = convert_position_str2float 
            
        # faster than 
        # x = np.array ( [ convert_position_str2float(i) for i in x ], 
        #               dtype = np.float64)
        # while apply_along not possible due to the string dtype during 
        # the loop
        # x = np.apply_along_axis(convert_position_str2float, 0, 
        #                         np.array (x, dtype =str ))
        return map ( lambda i: f(i) , ar)
    
    @staticmethod 
    def dms2ll (x:ArrayLike|str  , y:ArrayLike|str , *,  data=None): 
        """ Convert array x and y from DD:MM:SS to degree decimal -longitude 
        (x) and latitude (y). 
        
        Parameters
        -----------
        x, y: ArrayLike containing the degree-minutes-seconds (DMS) coordinates
           positions.
        Returns 
        --------
        x, y: Arraylike or str
           ArrayLike in degree decimal coordinates format. By default `x` and 
           `y` are longitude and latitude respectively. 
           If `x` and `y` has a string argument, data must be provided, 
           otherwise an error raises. 
           
        data: pd.DataFrame 
          DataFrame containg the position x and y. 
          
          .. versionadded:: 0.2.5 
          
        Example
        ---------
        >>> from watex.site import Profile 
        >>> x=['20:15:35'] ; y =['7:45:8.5'] 
        >>> Profile.dms2ll (x, y)
        Out[83]: (array([20.25972222]), array([7.75236111]))
        """
        if ( 
                (isinstance (x, str) 
                or isinstance (y, str)) 
                and data is None
                ): 
            raise TypeError(
                "Data can't be None when x or y has a string argument.")
        if data is not None: 
            x, y = assert_xy_in(x, y, data=data )
            
        if not _is_numeric_dtype(x , to_array =True ): 
           # reconvert object to str for consistency  
           x = np.array ( list ( Profile.f_(x)), dtype = np.float64 )
        if not _is_numeric_dtype(y , to_array =True): 
           y = np.array ( list (Profile.f_ (y)), dtype = np.float64 )
          
        return x, y
    
    @staticmethod 
    def ll2dms (x:ArrayLike|str  , y:ArrayLike|str , *, data=None ): 
        """
        Convert array x and y from degree decimal  to 
        degree-minutes-seconds (DMS)
        
        Parameters 
        -----------
        x, y: ArrayLike or str 
           Arrays containing the degree decimal position coordinates.
           If `x` and `y` has a string argument, data must be provided,
           otherwise an error raises
        
        Returns 
        --------
        x, y: Arraylike 
           ArrayLike in DD:MM:SS coordinates format.
          
        data: pd.DataFrame 
          DataFrame containg the position x and y.
          
          .. versionadded:: 0.2.5
          
        Example
        ---------
        >>> from watex.site import Profile 
        >>> x =[15.18 ] ; y =[19.60]
        >>> Profile.ll2dms (x, y)
        Out[84]: (array(['15:10:48.00'], dtype='<U11'), 
                  array(['19:36:00.00'], dtype='<U11'))
        """
        if ( (isinstance (x, str) 
              or isinstance (y, str)) 
              and data is None): 
            raise TypeError(
                "Data can't be None when x or y has a string argument.")
            
        if data is not None: 
            x, y = assert_xy_in(x, y, data=data )
        if _is_numeric_dtype(x , to_array =True ):  
           x = np.array ( list (Profile.f_ (x, 'll->dms')), dtype = str )
        if _is_numeric_dtype(y , to_array =True): 
           y = np.array ( list (Profile.f_ (y, 'll->dms')), dtype = str)
          
        return x, y
    
    def make_xy_coordinates (
        self, 
        *, 
        step =None,
        r = None, 
        todms:bool =False, 
        **kws 
        ): 
        """
        Generate synthetic coordinates from references latitude and 
        longitude from x and y. 
        
        Parameters 
        -----------
        step: float or str 
            Offset or the distance of seperation between different sites 
            in meters. If the value is given as string type, except 
            the ``km``, it should be considered as a ``m`` value. Only 
            meters and kilometers are accepables.
            
        r: float or int 
            The rotate angle in degrees. Rotate the angle features 
            toward the direction of the projection profile. 
            Default value use the :meth:`~.bearing` value in degrees. 
               
        todms: bool, Default=False
            Reconvert the longitude/latitude degree decimal values into 
            the DD:MM:SS. 
     
        kws: dict, 
           Additional keywords of :func:`~watex.utils.exmath.makeCoords`.   
           
        Returns 
        ----------
        lon, lat : ArrayLike 
           ArrayLike of synthetic coordinates latitude and longitude. 
          
        See Also 
        --------
        watex.utils.coreutils.makeCoords: 
            Create mutiple coordinates with different strategies using 
            references latitude and longitude. 
            
        Examples 
        ----------
        >>> import watex as wx 
        >>> # get two coordinates from random stations
        >>> data = wx.make_erp (n_stations =2 ).frame 
        >>> p = Profile ().fit(data.longitude, data.latitude)
        >>> longs, lats = p.make_xy_coordinates(step =100 , todms =True )
        >>> longs
        Out[40]: array(['-99:13:48.11', '-99:13:48.11'], dtype='<U12')
        >>> lats
        array(['-85:31:37.57', '-85:31:37.57'], dtype='<U12')
        >>> p = Profile ().fit(data.easting, data.northing)
        >>> longs, lats = p.make_xy_coordinates(step =100 , todms =True )
        >>> longs, lats
        Out[43]: 
        (array(['110:58:55.00', '110:58:55.00'], dtype='<U12'),
         array(['4:32:10.96', '4:32:10.96'], dtype='<U10'))
        """
        def _auto ( v): 
            if str (v).lower().strip() in ('none', 'auto'): 
                return None 
            
        self.inspect 
        step = _auto (step ); r = _auto (r )
        # use distance by default and bearing as r angle 
        step = step or self.distance (average_distance= True )
        r = r or self.bearing() 
  
        reflat = ( self.y[0], self.y[-1]) 
        reflong = (self.x[0], self.x[-1])
  
        isutm = False if self.coordinate_system =='ll' else True 
        utm_zone =  kws.pop ('utm_zone', None ) or self.utm_zone 
        
        return makeCoords(
            reflong, 
            reflat, 
            nsites = len(self.x ), 
            r= r ,  
            step =step , 
            todms=todms, 
            utm_zone= utm_zone, 
            is_utm= isutm, 
            datum=self.datum, 
            espg=self.epsg,
            **kws
            ) 
    
    def interpolate(
        self, 
        method ='linear', 
        inplace =True, 
        **kws
        ): 
        """
        Interpolate x, y and elev ( if applicable).
        
        
        Parameters 
        -----------
        inplace: bool, default=True 
           Erase existing value of x , y and elev with the interpolated one. 
           If ``False`` , it returns interpolated x, y and elev and 
           keep the previous attributes x, y and elev untouchable.

        method: bool, default='nearest', 
           Method of interpolation. One of ['nearest'|'linear'|'cubic'] 
         
        kws: dict,
           Additional keywords arguments passed to
           :func:`~watex.utils.funcutils.interpolate_grid`. 
           
        Returns
        --------
        self|x, y, elev: :class:`~watex.site.Profile` or Arraylikes 
           `:class:`.Profile` objects if `inplace` is ``True`` or 
           interpolated x, y and elev. 
           
        See Also 
        --------
        watex.utils.funcutils.interpolate_grid: 
            Interpolate two dimensional array. 
            
        Examples 
        --------
        >>> import numpy as np 
        >>> from watex.site import Profile 
        >>> x = [ 28, np.nan, 50, 60 ]
        >>> y =[ np.nan, 1000, 2000, 3000]
        >>> elev=[ 0, 1 , np.nan, np.nan]
        >>> po = Profile().fit(x, y, elev )
        >>> po.interpolate () 
        >>> po.x 
        array([28., 39., 50., 60.])
        >>> po.y 
        array([1000., 1000., 2000., 3000.])
        >>> po.elev 
        array([0., 1., 1., 1.])

        """ 
        emsg = ("Interpolation expects x, y and elev to have a consistent"
                " size. sizes x, y and elev are {}, {} and {}.")
        
        self.inspect 
        
        try :
            _check_consistency_size(self.x, self.y) 
            _check_consistency_size(self.y, self.elev)
        except: 
            raise ProfileError(emsg.format(
                len(self.x), len(self.y), len(self.elev )) )
            
        # make a grid data along axis =0 
        ar = np.vstack ((self.x, self.y, self.elev ))
        # interpolate is done along axis =0 so 
        # for x, y and elev we may transpose 
        # the data first and do the reverse 
        # processing back for new x, y and z 
        ar = interpolate_grid(ar, method = method , **kws) 
        
        x , y, elev = np.vsplit (ar  , 3 )
        x, y, elev = reshape (x) , reshape (y), reshape (elev)
        if inplace : 
            self.x, self.y, self.elev = x, y, elev 
            
            return self 
        
        return x, y, elev 

    @staticmethod 
    def out (
        x:ArrayLike , 
        y:ArrayLike, 
        /, 
        elev:ArrayLike =None, 
        filename:str=None, 
        basename:str=None , 
        counter:bool= False , 
        sep:str =None , 
        savepath:str =None, 
        extension:str ='.txt', 
        how:str='py', 
        verbose:int=0, 
        ):
        """Export coordinates data. 
        
        File to export should be ['name' | 'x' |'y' |'elev' ]. 
        
        Parameters 
        -----------
        x: Arraylike 
           X-coordinate values to export. 
           
        y: ArrayLike, 
          Y -coordinate values to export. 
          
        elev: ArrayLike, 
           Elevation values to export. If not given should be null. 
           
        filename: str, 
           Name of output file. 
           
        basename: str, default='Site'
           Station/site name of each coordinate (x, y)  position. If not given 
           should be ``site``. 
           
        counter: bool, default=False 
            If ``True``, count the number site and prefix it. For instance of 
            `42` sites with the ``basename=AMT``, the first counter site 
            should be ``00_AMT for Python indexing (``how='py'``)
            
        sep: str, default =' ' 
           the separator between site, x, y, elev. 
           
        savepath: str, 
           Output file destination. Default is current directory. 
           
        extension: str, default='.txt'
           The format of output coordinates files. If the `extension` is 
           suffixed with the `filename`, it should be used instead. 
           
        how: str, default='py'
          Way for counting the site/station. Any other value should starting 
          counting the site from 1. 
             
        verbose:int, default=False 
           show message to where the file is saved. 
           
        .. versionadded:: 0.2.1 
        
        
        Examples 
        ---------
        >>> import numpy as np 
        >>> import watex as wx 
        >>> data = wx.make_erp (n_stations=24).frame 
        >>> wx.site.Profile.out (data.easting, data.northing, verbose =True  )
        >>> elev = np.random.permutation  (len(data.easting)) *100 
        >>> wx.site.Profile.out (data.easting, data.northing, elev, 
                                 filename ='coordinates.txt', basename ='AMT-E1', 
                                 counter =True )
        """
       
        x = np.array( is_iterable(
            x, exclude_string =True, transform =True )) 
        y = np.array( is_iterable(y , exclude_string =True, transform =True )) 
        elev = np.array( is_iterable(elev , exclude_string =True, 
                                     transform =True )) 
        
        _check_consistency_size(x, y )
        
        elev = elev if elev is not None else np.repeat([0.], len(x))
        # for consitency 
        _check_consistency_size(y, elev )
        
        basename = basename or 'Site'
        basename = str (basename )
        filename= filename or 'site'
        filename = str(filename )
        
        sep = sep or " "
        sep = str(sep )
 
        # make site 
        sites = [(( f"{ii:02}_" if how=='py' else f"{ii+1:02}_") 
                  if counter else '') + basename 
                 for ii in range (len(x )) ]
        writelines =[]
        #x = x.astype (str) ; y = y.astype (str) ; elev= elev.astype (str)
        for site, east, north , elev in zip ( sites, x, y, elev ): 
            writelines.extend(
                [site + sep + 
                 f"{east:.3f}"  +sep +  
                 f"{north:.3f}"+ sep + 
                 f"{elev:.3f}"  + '\n' ]
                )

        filename, ex  = os.path.splitext(filename )
        extension = extension if ex =='' else ex
        
        extension ='.'+ extension.replace ('.', '') # for consistency 
        
        with open (filename + extension, 'w',encoding='utf8') as df:
            df.writelines(writelines)

        if savepath : 
            shutil.move (filename + extension, savepath )
            
        if verbose: 
            print(f"--> File {filename+extension} is successfull written" + 
                  (f"to {savepath}." if savepath else '.'))
        
    
    def plot (
        self, 
        sites: ArrayLike[str]=None , 
        basename:str ='S', 
        coerce:bool=True,  
        **kws): 
        """Base plot of the profile. 
        
        Parameters 
        -----------
        
        sites: str, Arraylike 
            The name of the sites that fits each position x, y
 
        basename: str, default='S' 
           the text to prefix the `site` position when the text is not given. 
        
        coerce:bool, default=True 
           Force the plot despite the given sites do not match the number of  
           positions `x` and `y`. If ``False``, number of positions must be 
           consistent with x and y, otherwise error raises. 
           
        kws: dict, 
           Additional keyword arguments passed to
           :func:`~watex.utils.plotutils.plot_text`
           
        .. versionadded:: 0.2.1 
        
        
        """
        from .utils.plotutils import plot_text 
        
        self.inspect 
        
        return plot_text (
            self.x, 
            self.y , 
            text = sites, 
            coerce =coerce, 
            basename=basename,
            xlabel= self.xname_ or '', 
            ylabel = self.yname_ or '', 
            **kws 
            )
 
    @property 
    def inspect (self): 
        """ Inspect object whether is fitted or not"""
        msg = ( "{obj.__class__.__name__} instance is not fitted yet."
               " Call 'fit' with appropriate arguments before using"
               " this method"
               )
        if not hasattr (self, 'x'): 
            raise NotFittedError(msg.format(
                obj=self)
            )
        return 1 
    
    def __repr__(self): 
        """ Represent the output class format """
        t_=("utm_zone" ,"coordinate_system","datum" ,
            "epsg" ,"reference_ellipsoid" ) 
        
        outm = ( '<{!r}:' + ', '.join(
            [f"{k}={getattr(self, k)!r}" for k in t_]) + '>' 
            ) 
        return  outm.format(self.__class__.__name__)
    
    
Profile.__doc__="""\
Profile class deals with the positions collected in the survey area. 

In principle, there is no exact solution  between longitude/latitude to 
x/y coordinates. Indeed, there is no isometric map from the sphere to the 
plane. when you convert lon/lat coordinates from the sphere to x/y 
coordinates in the plane, we cannot hope that all lengths will be preserved 
by this operation. Therefore, we have to accept some kind of deformation. 
For smallish parts of earth's surface, transverse Mercator is quite common.

By default, we use ``x`` for `longitude` and ``y`` for `latitude`. This is 
useful when data is passed as longitude-latitude (``ll``) or degree-minutes-
seconds (``dms``) in the `fit` method.

Parameters 
----------
utm_zone: Optional, string
   zone number and 'S' or 'N' e.g. '55S'. Default to the centre point
   of coordinates points in the survey area. It should be a string (##N or ##S)
   in the form of number and North or South hemisphere, 10S or 03N
coordinate_system: str, ['utm'|'dms'|'ll'] 
   The coordinate system in which the data points for the profile is collected. 
   If not given, the auto-detection will be triggered and find the  suitable 
   coordinate system. However, it is recommended to provide it for consistency. 
   Note that if `x` and `y` are composed of value less than 180 degrees 
   for longitude and 90 degrees for latitude, it should be considered as  
   longitude-latitude (``ll``) coordinates system. If `x` and `y` are 
   degree-minutes-second (``dms`` or ``dd:mm:ss``) data, they must be 
   specify as coordinate system in order to accept the non-numerical data 
   before transforming to ``ll``. If ``data`` is passed to the :meth:`.fit`
   method and ``dms`` is not specify, `x` and `y` values should be discarded.
   
datum: string, default = 'WGS84'
   well known datum ex. WGS84, NAD27, NAD83, etc.

epsg: Optional, int
   epsg number defining projection (
        see http://spatialreference.org/ref/ for moreinfo)
   Overrides utm_zone if both are provided. or 
   https://epsg.io/?q=China%20kind%3APROJCRS

reference_ellipsoid: int, default=23 
   reference ellipsoids is derived from Peter H. Dana's website-
   http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
   Department of Geography, University of Texas at Austin
   Internet: pdana@mail.utexas.edu . Default is ``23`` constrained to 
   WGS84. 
    
Examples 
--------- 
>>> from watex.datasets import load_edis
>>> from watex.site import Profile  
>>> xy = load_edis (samples =17 , as_frame =True , key ='latitude longitude') 
>>> xy.head (2)
    longitude   latitude
0  110.485833  26.051390
1  110.486153  26.051794 
>>> pro= Profile ().fit( xy.longitude, xy.latitude) 
>>> pro.distance ()
62.890276656978244
>>> pro.bearing () 
35.4252016495945
>>> pro.make_xy_coordinates( ) 
(array([110.48583316, 110.48615319, 110.48647322, 110.48679325,
       110.48711328, 110.48743331, 110.48775334]), 
 array([26.05138954, 26.05179371, 26.05219788, 26.05260205, 26.05300622,
       26.05341038, 26.05381455]))
"""
    
class Location (object): 

    def __init__(
        self, 
        lat= None, 
        lon=None,
        elev=None, 
        datum='WGS84', 
        epsg=None,
        reference_ellipsoid=23,
        utm_zone=None, 
        **kwds
        ) :
        self._logging = watexlog.get_watex_logger(
            self.__class__.__name__)
        self._lat = lat 
        self._lon = lon 
        self._utm_zone=utm_zone
        self._elev=elev
        
        self.datum=datum
        self.epsg=epsg 
        self.reference_ellipsoid=reference_ellipsoid
        
        self._east= None 
        self._north =None 
        
        for key in list(kwds.keys()) : 
            setattr (self, key, kwds[key])
            
    @property 
    def utm_zone (self): 
        return self._utm_zone
        
    @property 
    def lat(self): 
        return self._lat
    @property 
    def lon(self) : 
        return self._lon
    @property 
    def east(self ): 
        return self._east
    @property 
    def north(self):
        return self._north
    @property 
    def elev(self): 
        return self._elev

    @utm_zone.setter 
    def utm_zone (self, utm_zone): 
        utm_xone = copy.deepcopy(utm_zone) 
        utm_zone = str(utm_zone).upper().strip() 
        str_ = f"{'|'.join([ i for i in UTM_DESIGNATOR.keys()]).lower()}"
        regex= re.compile(rf'{str_}', flags=re.IGNORECASE)
        if regex.search(utm_zone) is None: 
           raise SiteError (f"Wrong UTM zone value!: {utm_xone!r} ")
        self._utm_zone =utm_zone.upper() 
    
    @lat.setter 
    def lat (self, lat): 
        self._lat = assert_lat_value(lat)
        
        
    @lon.setter 
    def lon(self, lon) : 
        self._lon = assert_lon_value(lon)
        
    @elev.setter 
    def elev(self, elev) : 
        self._elev = assert_elevation_value(elev ) 
        
    @east.setter 
    def east (self, east): 
        self._east = np.array(east, dtype = float )

    @north.setter 
    def north (self, north): 
        self._north = np.array(north, dtype =float)
        
    def to_utm (
        self, 
        lat:float=None, 
        lon:float=None, 
        datum:str= None ,
        utm_zone:str=None, 
        epsg:int=None, 
        reference_ellipsoid:int= None
        ): 
        """
        Project coordinates to utm if coordinates are in degrees at  
        given reference ellipsoid constrained to WGS84 by default. 
         
        Parameters 
        -----------
        lat: float or string (DD:MM:SS.ms)
            latitude of point
                  
        lon: float or string (DD:MM:SS.ms)
            longitude of point
        
        datum: string, default='WGS84'
            well known datum ex. WGS84, NAD27, NAD83, etc.

        utm_zone: Optional, string
            zone number and 'S' or 'N' e.g. '55S'. Defaults to the centre point
            of the provided points
                       
        epsg: Optional, int
            epsg number defining projection (see http://spatialreference.org/ref/ for moreinfo)
            Overrides utm_zone if both are provided

        reference_ellipsoid: Optional, int 
            reference ellipsoids is derived from Peter H. Dana's website-
            http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
            Department of Geography, University of Texas at Austin
            Internet: pdana@mail.utexas.edu . Default is ``23`` constrained to 
            WGS84. 
             
        Returns
        --------
        proj_point: tuple(easting, northing, zone)
            projected point in UTM in Datum
            
        """
        
        if lat is not None : 
            self.lat = lat
        if lon is not None : 
            self.lon =lon
            
        if (self.lat or self.lon) is None : 
            raise SiteError (" Latitude and longitude must not be None")
            
        self.epsg = epsg or self.epsg 
        self.datum = datum  or self.datum 
        
        self.reference_ellipsoid = reference_ellipsoid or \
            self.reference_ellipsoid
     
        try : 
            self.east, self.north, self.utm_zone= project_point_ll2utm(
                lat = self.lat ,
                lon= self.lon,
                datum = self.datum , 
                utm_zone = self.utm_zone ,
                epsg= self.epsg 
                ) 
        except :
            self.utm_zone, self.east, self.north= ll_to_utm(
                reference_ellipsoid=self.reference_ellipsoid,
                lat = self.lat, 
                lon= self.lon)

        return float(self.east), float(self.north)
    
        
    def to_latlon(
        self, 
        east:float=None, 
        north:float= None, 
        utm_zone:str=None, 
        reference_ellipsoid:int=None , 
        datum:str = None, 
        epsg: int= None, 
        ): 
        """
        Project coodinate on longitude latitude once  data are utm at  given
        reference ellispoid constrained to WGS-84 by default. 
        
        Parameters 
        -----------
        east: float
            easting coordinate in meters
                    
        north: float
            northing coordinate in meters
        
        utm_zone: Optional, string (##N or ##S)
            utm zone in the form of number and North or South hemisphere, 10S or 03N
        
        datum: string, default ='WGS84'
            well known datum ex. WGS84, NAD27, etc.
            
        epsg: Optional, int
            epsg number defining projection (see http://spatialreference.org/ref/ for moreinfo)
            Overrides utm_zone if both are provided    
            
        reference_ellipsoid: Optional, int 
            reference ellipsoids is derived from Peter H. Dana's website-
            http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
            Department of Geography, University of Texas at Austin
            Internet: pdana@mail.utexas.edu . Default is ``23`` constrained to 
            WGS84.
            
        Returns
        -------
        proj_point: tuple(lat, lon)
            projected point in lat and lon in Datum, as decimal degrees.
                        
        """
        if east is not None: 
            self.east = east 
        if north is not None: 
            self.north = north 
        
        if (self.east or self.north) is None : 
            raise SiteError(" Easting and northing must not be None")
            
        if utm_zone is not None : 
            self._utm_zone =utm_zone 
            
        self.epsg = epsg or self.epsg 
        self.datum = datum or self.datum 
        
        self.reference_ellipsoid = reference_ellipsoid or \
            self.reference_ellipsoid
                          
        try : 
            self.lat , self.lon = project_point_utm2ll(
                easting = self.east,
                northing= self.north, 
                utm_zone = self.utm_zone, 
                datum= self.datum, 
                epsg= self.epsg 
                ) 
        except : 
             self.lat, self.lon  = utm_to_ll(
                reference_ellipsoid=self.reference_ellipsoid, 
                northing = self.north , 
                easting= self.east, 
                zone= self.utm_zone 
                                  ) 
            
        return self.lat, self.lon 
        
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
           
           .. versionadded:: 0.2.5
           
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
           
        .. versionadded:: 0.1.8 
        
        See Also
        ----------
        watex.site.Location.to_utm: 
            Convert longitude and latitude value to easting and northing 
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
           
           .. versionadded:: 0.2.5
           
        datum: string
           well known datum ex. WGS84, NAD27, NAD83, etc.
            
        kws: dict, 
           Keywords argument passed to :meth:`~watex.site.Location.to_latlon`. 
           
        Returns 
        -------
        (lats, lons): Iterable object composed of latitude and longitude 
           coordinates. 
          
        .. versionadded:: 0.1.8 
        
        See Also
        ------------
        watex.site.Location.to_latlon: 
           Convert easting and northing value to latitude and  longitude 
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
                   
Location.__doc__="""\
Location class

Assert, convert coordinates lat/lon , east/north  
to appropriate formats. 

Location does not follow the `scikit-learn` API in order to encompass the 
the syntax of :term:`pycsamt` and :term:`mtpy`. The latter follows the 
Generci Mapping Tool (`GMT <ttps://www.generic-mapping-tools.org>`__)  API 
format.   
  
Parameters  
--------------
lat: float or string (DD:MM:SS.ms)
    latitude of point
          
lon: float or string (DD:MM:SS.ms)
    longitude of point
    
east: float
    easting coordinate in meters
            
north: float
    northing coordinate in meters
     
datum: string
    well known datum ex. WGS84, NAD27, NAD83, etc.

utm_zone: Optional, string
    zone number and 'S' or 'N' e.g. '55S'. Defaults to the centre point
    of the provided points
               
epsg: Optional, int
    epsg number defining projection (see http://spatialreference.org/ref/ for moreinfo)
    Overrides utm_zone if both are provided

reference_ellipsoid: Optional, int 
    reference ellipsoids is derived from Peter H. Dana's website-
    http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
    Department of Geography, University of Texas at Austin
    Internet: pdana@mail.utexas.edu . Default is ``23`` constrained to 
    WGS84. 
"""     
        
       

     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        