# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Survey layout.

- Station: robust line-geometry container that understands modern/legacy
  AVG frames, exposes unique station positions, IDs, and handy header
  ($Stn.*) derivations for round-trips.
- Topography: robust container for .stn station location files.

"""
from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from dataclasses import field
from pathlib import Path
from typing import (
    Any,
    Literal,
)

import numpy as np
import pandas as pd
from scipy.interpolate import griddata

from ..compat.python import dc
from ..exceptions import ProcessingError, StationError
from ..gis.utils import (
    assert_xy_coordinate_system,
    calculate_azimuth,
    get_elevation_from_api,
    get_elevation_from_utm,
    normalize_lat_lon,
    to_ll,
    to_utm,
)
from ..log.logger import get_logger
from ..utils.deps import import_optional_dependency
from ..utils.validation import ensure_n_items
from .base import AVGComponentBase
from .utils import (
    find_and_rename_column,
    number_stations,
    read_stn,
)

__all__ = ["Station", "Topography"]


logger = get_logger(__name__)


class Topography(AVGComponentBase):
    r"""A container for station topography and location data.

    This class is designed to read, manage, and process spatial
    information from Zonge ``.stn`` files, as described in the
    ASTATIC manual [1]_. It handles both legacy (space-delimited)
    and modern (comma-delimited) formats and provides a suite of
    tools for coordinate conversion, regularization, and gridding.

    Parameters
    ----------
    data : pandas.DataFrame, optional
        A pre-loaded DataFrame containing station location data.
        If provided, it will be standardized upon initialization.
    meta : mapping, optional
        An optional dictionary for metadata, consistent with the
        :class:`~.base.AVGComponentBase` API.
    verbose : bool, default False
        Controls the level of detail in logging output.

    Attributes
    ----------
    stations : numpy.ndarray
        An array of the station numbers or identifiers.
    eastings, northings : numpy.ndarray
        Arrays of the station UTM coordinates.
    elevations : numpy.ndarray
        An array of the station elevations.

    Methods
    -------
    read(source)
        Reads and parses a ``.stn`` file or a DataFrame.
    generate(...)
        A static method to create a synthetic survey line.
    correct_coords(...)
        Regularizes station locations to a best-fit straight line.
    convert_coords(...)
        Converts coordinates between UTM and geographic (lat/lon).
    to_grid(...)
        Interpolates the scattered station data onto a regular 2D
        grid for contouring.
    get_step()
        Calculates the distance between consecutive stations.

    Notes
    -----
    The `read` method is designed to be robust, automatically
    detecting the delimiter and header row of ``.stn`` files.
    Upon reading, it standardizes column names to a canonical
    schema (e.g., 'easting', 'northing', 'elevation'), making
    the data consistent for all subsequent processing steps.

    Examples
    --------
    >>> from pycsamt.zonge.survey import Topography
    >>> # Load topography from a .stn file
    >>> topo = Topography().read('data/avg/K1.stn')
    >>>
    >>> # Generate a synthetic survey line
    >>> synthetic_topo = Topography.generate(
    ...     start_coord=(500000, 4000000),
    ...     n_stations=20,
    ...     step=50,
    ...     azimuth=45
    ... )
    >>> # Get the average station spacing
    >>> avg_step = synthetic_topo.get_step().mean()

    References
    ----------
    .. [1] Zonge International, Inc. (2014). *ASTATIC v3.70
           User Manual*, "STN Files" section, p. 36.
    """
    def __init__(
        self,
        data: pd.DataFrame | None = None,
        meta: Mapping[str, Any] | None = None,
        *,
        utm_zone: str =None,
        epsg: int =None,
        name: str | None = None,
        verbose: bool = False
    ) -> None:
        super().__init__(
            data=data, meta=meta, name=name or "Topography",
            verbose=verbose
        )
        self.utm_zone =utm_zone
        self.epsg= epsg
        self._azimuths: np.ndarray | None = None

    def read(
        self,
        source: str | Path | pd.DataFrame,
        meta: Mapping[str, Any] | None = None,
        **kws: Any
    ) -> Topography:
        r"""Read topography data from a file path or DataFrame.

        This is the primary data ingestion method for the Topography
        class. It is designed to be robust, handling various ``.stn``
        file formats and standardizing the data into a consistent
        internal structure.

        Parameters
        ----------
        source : str, pathlib.Path, or pandas.DataFrame
            The data source to load. This can be:
            - A string or `pathlib.Path` pointing to a Zonge ``.stn``
              file.
            - A `pandas.DataFrame` containing station location data.
        meta : mapping, optional
            An optional dictionary for metadata, consistent with the
            :class:`~.base.AVGComponentBase` API. This is typically
            not used for ``.stn`` files.

        Returns
        -------
        self : Topography
            The method returns the instance of the class, allowing
            for convenient method chaining.

        Raises
        ------
        TypeError
            If the `source` is of an unsupported type.
        StationError
            If the ``.stn`` file is malformed or a valid header row
            cannot be found.

        Notes
        -----
        The file parser is designed to be flexible:
        - It automatically detects whether the file is comma-delimited
          or space-delimited.
        - It intelligently searches for a header row by looking for
          common keywords (e.g., 'station', 'easting'), rather than
          assuming a fixed position.
        - All column names are normalized to a canonical schema
          (e.g., 'easting', 'northing', 'elevation') after reading.

        Examples
        --------
        >>> from pycsamt.zonge.survey import Topography
        >>> # Load topography from a .stn file
        >>> topo = Topography().read('data/avg/K1.stn')
        >>> print(topo.stations[:5])
        [150. 200. 250. 300. 350.]
        """
        if isinstance(source, (str, Path)):
            df = read_stn(source)
        elif isinstance(source, pd.DataFrame):
            df = source.copy()
        else:
            raise TypeError("Source must be a file path or DataFrame.")

        self._frame = self._normalize_stn_columns(df)

        # initialize longitude and latitude
        self._longitude = pd.Series (
            np.zeros((self._frame.shape[0],), dtype=float ))
        self._latitude = pd.Series (
            np.zeros((self._frame.shape[0],), dtype=float )
        )

        # reset the cache when new data is loaded
        self._azimuths = None

        return self

    def _normalize_stn_columns(
            self, df: pd.DataFrame
        ) -> pd.DataFrame:
        """
        Find and rename columns to a canonical schema.
        """
        rename_map = {}
        for col in df.columns:
            col_lower = str(col).lower()
            if "station" in col_lower or "dot" in col_lower:
                rename_map[col] = "station"
            elif "east" in col_lower or col_lower == 'e':
                rename_map[col] = "easting"
            elif "north" in col_lower or col_lower == 'n':
                rename_map[col] = "northing"
            elif "elev" in col_lower or col_lower == 'h':
                rename_map[col] = "elevation"
            elif "head" in col_lower: rename_map[col] = "heading"
            elif "pitch" in col_lower: rename_map[col] = "pitch"
            elif "roll" in col_lower: rename_map[col] = "roll"

        df = df.rename(columns=rename_map)
        # Convert all columns to numeric
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        required = ["station", "easting", "northing", "elevation"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ProcessingError(
                f"STN data is missing required columns: {missing}"
            )
        return df.dropna(subset=required)

    def convert_coords(
        self,
        to: Literal["utm", "ll", "auto"] = "auto",
        *,
        inplace: bool = True
    ) -> pd.DataFrame | None:
        r"""Convert between UTM and geographic (lat/lon) coordinates.

        This method provides a high-level interface for coordinate
        system conversion, leveraging the underlying utilities in the
        :mod:`~pycsamt.gis.utils` module.

        Parameters
        ----------
        to : {'utm', 'll', 'auto'}, default 'auto'
            The target coordinate system.

            - 'utm': Convert to Universal Transverse Mercator.
            - 'll': Convert to latitude/longitude decimal degrees.
            - 'auto': Automatically detect the current system and
              convert to the other.

        update_inplace : bool, default True
            If ``True``, the internal DataFrame of the `Topography`
            object is updated with the new coordinate columns.
            If ``False``, a new DataFrame containing both the original
            and converted coordinates is returned.

        Returns
        -------
        pandas.DataFrame or None
            - If `update_inplace` is ``True``, returns ``None`` and
              modifies the object's internal frame.
            - If `update_inplace` is ``False``, returns a new DataFrame
              with the added coordinate columns.

        Notes
        -----
        The function first uses
        :func:`~.gis.utils.assert_xy_coordinate_system` to determine
        the current coordinate system of the data. It then calls either
        :func:`~.gis.utils.to_utm` or :func:`~.gis.utils.to_ll` to
        perform the conversion.

        See Also
        --------
        pycsamt.gis.utils.to_utm : The underlying UTM conversion utility.
        pycsamt.gis.utils.to_ll : The underlying lat/lon conversion utility.
        """
        if self._frame.empty:
            raise ProcessingError("Topography data has not been loaded.")

        x_coords = self.easting
        y_coords = self.northing

        current_system = assert_xy_coordinate_system(x_coords, y_coords)

        if to == "auto":
            to = "ll" if current_system == "utm" else "utm"

        if to == current_system:
            if self.verbose:
                self._logger.info(
                    f"Coordinates are already in '{to}' system. "
                    "No conversion performed."
                )
            return

        if to == "ll":
            if self.utm_zone is None and self.epsg is None:
                self._logger.error (
                    "UTM zone or 'epsg' needs to be set for"
                    " converting easting/northing"
                    " to longitude/latitude"
                )
                new_df = pd.DataFrame(
                    {"latitude": self._latitude,
                     "longitude": self._longitude
                     },
                    index=self._frame.index
                )
            else:
                lat, lon = to_ll(
                    x_coords, y_coords,
                    zone =self.utm_zone,
                    epsg = self.epsg,
                    as_frame=False
                )
                new_df = pd.DataFrame(
                    {"latitude": lat, "longitude": lon},
                    index=self._frame.index
                )
            self._longitude = new_df['longitude']
            self._latitude = new_df['latitude']

        elif to == "utm":
            east, north, zone = to_utm(
                x_coords, y_coords,
                epsg= self.epsg,
                utm_zone= self.utm_zone,
                as_frame=False
            )
            new_df = pd.DataFrame(
                {"easting": east, "northing": north, "utm_zone": zone},
                index=self._frame.index
            )

        else:
            raise ValueError(
                f"Invalid target system '{to}'. Must be "
                "'utm', 'll', or 'auto'."
            )

        if inplace:
            for col in new_df.columns:
                self._frame[col] = new_df[col]
            if self.verbose:
                self._logger.info(
                    f"Coordinates converted to '{to}' and updated "
                    "in place."
                )

            return None

        return pd.concat([self._frame, new_df], axis=1)

    def to_grid(
        self,
        resolution: int = 100,
        method: Literal['linear', 'cubic', 'nearest'] = 'cubic'
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r"""Interpolate scattered station data onto a regular 2D grid.

        This method takes the scattered station locations (easting,
        northing) and their corresponding elevations and
        interpolates them onto a regular grid, which is essential
        for creating contour maps or surface plots.

        Parameters
        ----------
        resolution : int, default 100
            The number of points to create in each dimension (x and
            y) of the output grid. A higher number results in a
            finer grid.
        method : {'linear', 'cubic', 'nearest'}, default 'cubic'
            The interpolation method to be used by the underlying
            :func:`scipy.interpolate.griddata` function.
            - 'linear': Performs linear interpolation.
            - 'cubic': Performs cubic interpolation for a smoother
              surface.
            - 'nearest': Uses the value of the nearest data point.

        Returns
        -------
        grid_x : numpy.ndarray
            A 2D array of the X (easting) coordinates of the grid.
        grid_y : numpy.ndarray
            A 2D array of the Y (northing) coordinates of the grid.
        grid_z : numpy.ndarray
            A 2D array of the interpolated Z (elevation) values on
            the grid.

        Raises
        ------
        ProcessingError
            If the topography data has not been loaded first.

        See Also
        --------
        scipy.interpolate.griddata : The core interpolation
            function used by this method.
        pycsamt.zonge.plot.Plot.plot_location_map : A method for
            visualizing the gridded data.
        """

        if self._frame.empty:
            raise ProcessingError("Topography data has not been loaded.")

        points = self._frame[['easting', 'northing']].values
        values = self.elevations

        # Create the target grid
        grid_x, grid_y = np.mgrid[
            self.easting.min():self.easting.max():complex(resolution),
            self.northing.min():self.northing.max():complex(resolution)
        ]

        # Interpolate the data
        grid_z = griddata(
            points, values, (grid_x, grid_y), method=method)

        return grid_x, grid_y, grid_z

    def get_step(self) -> pd.Series:
        r"""Calculate the distance between consecutive stations.

        This method computes the Euclidean distance between each
        station and the next one along the survey line, based on
        their easting and northing coordinates.

        Returns
        -------
        pandas.Series
            A Series containing the calculated step distances in
            meters. The first value is always 0. The length of the
            Series matches the number of stations.

        Notes
        -----
        The calculation assumes a Cartesian coordinate system (like
        UTM) where the Pythagorean theorem can be applied to find
        the distance between points. The result is useful for
        assessing the regularity of station spacing and for
        providing a default step size for other methods like
        :meth:`correct_coords`.

        Examples
        --------
        >>> from pycsamt.zonge.survey import Topography
        >>> topo = Topography().read('data/avg/K1.stn')
        >>> steps = topo.get_step()
        >>> print(f"Average station spacing: {steps.mean():.2f} m")
        """
        if self._frame.empty or len(self._frame) < 2:
            return pd.Series(dtype=float)

        dx = np.diff(self.easting)
        dy = np.diff(self.northing)
        steps = np.hypot(dx, dy)

        return pd.Series(
            np.concatenate(([0], steps)),
            index=self._frame.index
        )

    def correct_coords(
        self,
        step: float | None = None,
        *,
        inplace: bool = True
    ) -> pd.DataFrame | None:
        r"""Regularize station coordinates to a best-fit straight line.

        This processing tool corrects for minor deviations in survey
        line geometry by projecting all station locations onto a
        best-fit straight line and re-spacing them at a uniform
        interval. This is a common step for preparing data for 2D
        inversion or gridding.

        Parameters
        ----------
        step : float, optional
            The desired uniform distance between stations along the
            corrected line. If ``None``, the average step distance is
            calculated automatically using the :meth:`get_step`
            method.
        update_inplace : bool, default True
            If ``True``, the internal DataFrame of the `Topography`
            object is updated in place with the new, corrected
            coordinates. If ``False``, a new DataFrame with the
            corrected coordinates is returned.

        Returns
        -------
        pandas.DataFrame or None
            - If `update_inplace` is ``True``, returns ``None``.
            - If `update_inplace` is ``False``, returns a new
              DataFrame with the corrected 'easting' and 'northing'
              columns.

        Notes
        -----
        The correction process involves two main steps:
        1.  A first-degree polynomial (a straight line) is fitted
            to the original easting and northing coordinates using a
            least-squares regression.
        2.  New station locations are generated along this ideal
            line at a constant spacing defined by the `step`
            parameter.

        See Also
        --------
        get_step : The method used to automatically determine the
            average station spacing.
        """
        if self._frame.empty or len(self._frame) < 2:
            warnings.warn(
                "Not enough data to correct coordinates.", stacklevel=2)
            return

        x = self.easting
        y = self.northing

        # 1. Determine the best-fit line
        m, c = np.polyfit(x, y, 1) # y = mx + c

        # 2. Calculate the uniform step distance
        if step is None:
            step = self.get_step().mean()
            if self.verbose:
                self._logger.info(
                    f"Using auto-detected average step of {step:.2f} m."
                )

        # 3. Project the first station onto the line
        x0, y0 = x[0], y[0]
        x_proj_start = (x0 + m * (y0 - c)) / (1 + m**2)
        y_proj_start = m * x_proj_start + c

        # 4. Generate new points along the line
        line_direction = np.array([1, m]) / np.sqrt(1 + m**2)
        distances = np.arange(len(x)) * step

        new_eastings = x_proj_start + distances * line_direction[0]
        new_northings = y_proj_start + distances * line_direction[1]

        if inplace:
            self._frame['easting'] = new_eastings
            self._frame['northing'] = new_northings
            if self.verbose:
                self._logger.info(
                    "Coordinates have been corrected and updated in place."
                )

            return None

        df_corrected = self._frame.copy()
        df_corrected['easting'] = new_eastings
        df_corrected['northing'] = new_northings

        return df_corrected

    @staticmethod
    def generate(
        start_coord: tuple[float, float],
        n_stations: int,
        step: float,
        azimuth: float,
        *,
        initial_station_name: float = 0.,
        initial_elevation: float = 0.,
        elevation_gradient: float = 0.,
        coord_type: Literal['utm', 'll'] = 'utm'
    ) -> Topography:
        r"""Generate a synthetic station topography dataset.

        This static method acts as a factory for creating a new
        :class:`Topography` object based on survey design
        parameters. It is useful for creating test data, planning
        surveys, or generating a regularized coordinate set for
        modeling.

        Parameters
        ----------
        start_coord : tuple[float, float]
            The starting coordinate of the survey line. The format
            depends on `coord_type`:

            - For 'utm': (easting, northing) in meters.
            - For 'll': (latitude, longitude) in decimal degrees.

        n_stations : int
            The total number of stations to generate along the line.
        step : float
            The distance between consecutive stations, in meters.
        azimuth : float
            The azimuth of the survey line in degrees clockwise from
            North (e.g., 0 for North, 90 for East).
        initial_station_name : float, default 0.0
            The name or number of the first station. Subsequent
            station names are incremented by `step`.
        initial_elevation : float, default 0.0
            The elevation of the first station, in meters.
        elevation_gradient : float, default 0.0
            The change in elevation per meter along the survey line.
            A positive value creates an upward slope.
        coord_type : {'utm', 'll'}, default 'utm'
            The coordinate system of the `start_coord`.

        Returns
        -------
        Topography
            A new, fully populated instance of the `Topography` class.

        Notes
        -----
        When `coord_type` is 'll', the function uses geodetic
        calculations to accurately generate new points along the
        great-circle path and then converts the final lat/lon
        coordinates to UTM for storage.
        """
        if n_stations <= 0:
            raise ValueError("Number of stations must be positive.")

        station_names = np.arange(
            n_stations) * step + initial_station_name
        distances = np.arange(n_stations) * step
        elevations = initial_elevation + distances * elevation_gradient

        start_coord =ensure_n_items (
            items=start_coord,
            name ='(latitude, longitude)',
            expect="numeric",
            coerce=True,
            n=2
        )
        if coord_type == 'll':
            import_optional_dependency(
                'geopy', extra=(
                    "'geopy' is required for geodetectic"
                    " position calculations")
               )
            # For lat/lon, we must calculate geodetic positions
            from geopy.distance import geodesic
            start_lat, start_lon = normalize_lat_lon (
                *start_coord, assume='latlon')
            # start_lat, start_lon = start_coord

            points = [start_coord]
            for _i in range(1, n_stations):
                new_point = geodesic(
                    meters=step).destination(points[-1], bearing=azimuth)
                points.append((new_point.latitude, new_point.longitude))

            latitudes, longitudes = zip(*points)
            eastings, northings, _ = to_utm(latitudes, longitudes)

        else: # UTM coordinates
            azimuth_rad = np.deg2rad(90 - azimuth) # Convert from bearing
            dx = distances * np.cos(azimuth_rad)
            dy = distances * np.sin(azimuth_rad)
            start_x, start_y = start_coord
            eastings = start_x + dx
            northings = start_y + dy

        df = pd.DataFrame({
            'station': station_names,
            'easting': eastings,
            'northing': northings,
            'elevation': elevations
        })

        return Topography(data=df)

    def write(self) -> Sequence[str]:
        r"""Serialize the topography data to .stn format lines.

        This method converts the internal DataFrame of topography
        data into a list of strings that conform to the standard
        Zonge ``.stn`` file format. This is useful for exporting
        processed or generated topography data.

        Returns
        -------
        list[str]
            A list of strings, where the first string is the comma-
            separated header and subsequent strings are the data
            rows.

        Notes
        -----
        The output is designed to be a "modern" ``.stn`` file, using
        comma-separated values, which is compatible with the ASTATIC
        program and other Zonge software.

        Examples
        --------
        >>> from pycsamt.zonge.survey import Topography
        >>> topo = Topography.generate(
        ...     start_coord=(500000, 4000000), n_stations=3,
        ...     step=100, azimuth=90
        ... )
        >>> stn_lines = topo.write()
        >>> for line in stn_lines:
        ...     print(line)
        station,easting,northing,elevation
        0.0,500000.0,4000000.0,0.0
        100.0,500100.0,4000000.0,0.0
        200.0,500200.0,4000000.0,0.0
        """
        if self._frame.empty:
            return []

        header = ",".join(self._frame.columns)
        data_lines = self._frame.to_csv(
            index=False, header=False, lineterminator='\n'
        ).strip('\n').split('\n')

        return [header] + data_lines

    def get_elevation_from(
        self,
        from_: Literal["utm", "api"] = "utm",
        zone: str | None = None,
        datum: str = "WGS84",
    ) -> np.ndarray:
        r"""Fetches station elevations from an external API.

        This method retrieves elevation data for each station in
        the survey using an online service. It can operate using
        either the existing UTM coordinates or geographic (lat/lon)
        coordinates.

        Parameters
        ----------
        from_ : {'utm', 'api'}, default 'utm'
            Specifies the coordinate system to use for the API
            query.
            - 'utm': Uses the ``easting`` and ``northing`` columns.
              The ``zone`` parameter is required.
            - 'api': Uses latitude and longitude. If these are not
              present, they are automatically calculated from the
              UTM coordinates.
        zone : str, optional
            The UTM zone designator (e.g., "11S", "32N"). This is
            **required** when ``from_`` is 'utm'.
        datum : str, default "WGS84"
            The geodetic datum to assume for coordinate
            conversions.

        Returns
        -------
        numpy.ndarray
            An array of elevation values in meters corresponding
            to each station.

        Raises
        ------
        ProcessingError
            If the topography data has not been loaded first.
        ValueError
            If ``from_`` is 'utm' but no ``zone`` is provided, or
            if an invalid ``from_`` value is given.

        See Also
        --------
        get_elevation_from_utm : Fetches elevation from UTM data.
        get_elevation_from_api : Fetches elevation from lat/lon data.
        convert_coords : Converts between coordinate systems.
        """
        # Ensure the topography data frame has been loaded first
        if self._frame.empty:
            raise ProcessingError(
                "Topography data has not been loaded. "
                "Call the .read() method first."
            )

        elevations = np.array([], dtype=float)

        zone = zone or self.utm_zone
        if from_ == "utm":
            # The UTM zone is required for this operation
            if zone is None:
                raise ValueError(
                    "A UTM 'zone' must be provided when "
                    "from_='utm'."
                )
            # Fetch elevations using existing easting/northing
            elevations = get_elevation_from_utm(
                easting=self.easting,
                northing=self.northing,
                zone=zone,
                datum=datum,
            )

        elif from_ == "api":
            # Check if lat/lon data exists; if not, create it
            if 'latitude' not in self._frame.columns or (
                    'longitude' not in self._frame.columns):
                if self.verbose:
                    self._logger.info(
                        "Latitude/Longitude not found. "
                        "Converting from UTM."
                    )
                self.convert_coords(to='ll', inplace=True)

            # Fetch elevations using latitude and longitude
            elevations = get_elevation_from_api(
                latitude=self.latitude,
                longitude=self.longitude,
            )
        else:
            raise ValueError(
                f"Invalid 'from_' argument: '{from_}'. "
                "Must be 'utm' or 'api'."
            )

        return elevations

    def get_azimuth(
        self,
        *,
        mode: Literal["mean", "median"] | None = None
    ) -> float | np.ndarray:
        r"""Returns the survey line azimuth(s).

        This method can return the full array of segment azimuths
        or a single aggregated value (mean or median) representing
        the overall trend of the survey line.

        Parameters
        ----------
        mode : {'mean', 'median'}, optional
            Determines the return format.
            - ``None`` (default): Returns the full NumPy array of
              azimuths for each station segment.
            - 'mean': Returns the circular mean of all segment
              azimuths, providing the average direction.
            - 'median': Returns the median of all segment azimuths.

        Returns
        -------
        float or numpy.ndarray
            - If ``mode`` is ``None``, returns the full array of
              azimuths.
            - If ``mode`` is 'mean' or 'median', returns a single
              float representing the aggregated value.

        Notes
        -----
        The 'mean' is calculated using circular statistics to
        correctly average angles (e.g., the mean of 350° and 10°
        is 0°, not 180°). The 'median' is calculated linearly and
        is a good approximation for relatively straight lines.
        """
        # Retrieve the azimuths using the caching property
        azimuths = self.azimuth

        if mode is None:
            return azimuths

        # Filter out NaN values before aggregating
        valid_azimuths = azimuths[~np.isnan(azimuths)]

        if valid_azimuths.size == 0:
            return np.nan if mode else np.array([])

        if mode == "mean":
            # Correctly average angles using circular statistics
            rads = np.deg2rad(valid_azimuths)
            sin_mean = np.mean(np.sin(rads))
            cos_mean = np.mean(np.cos(rads))
            mean_rad = np.arctan2(sin_mean, cos_mean)
            mean_deg = np.rad2deg(mean_rad)
            # Normalize to 0-360 range
            return (mean_deg + 360) % 360

        if mode == "median":
            return np.median(valid_azimuths)

        raise ValueError(
            f"Invalid mode '{mode}'. Choose from 'mean', 'median',"
            " or None."
        )

    @property
    def azimuth(self) -> np.ndarray:
        r"""Calculates and returns the forward azimuths.

        This property computes the azimuth for each line segment
        (from one station to the next) along the survey line. The
        result is cached for efficiency.

        Returns
        -------
        numpy.ndarray
            An array of azimuths in degrees (0-360), where the
            azimuth at index ``i`` corresponds to the direction
            from station ``i`` to ``i+1``.
        """
        # Check if the azimuths have already been calculated
        if self._azimuths is None:
            # Ensure there's enough data to calculate a direction
            if self.easting.size < 2:
                self._azimuths = np.array([], dtype=float)
            else:
                # If not cached, calculate and store the result
                self._azimuths = calculate_azimuth(
                    self.easting, self.northing
                )
        return self._azimuths

    @property
    def bearing(self) -> np.ndarray:
        r"""Alias for the azimuth property.

        In this context, bearing is treated as synonymous with
        azimuth (0-360 degrees clockwise from North).
        """
        return self.azimuth

    @property
    def stations(self) -> np.ndarray:
        # This will be deprecated and
        return self._frame.get(
            "station", pd.Series(dtype=float)).values

    @property
    def easting(self) -> np.ndarray:
        return self._frame.get(
            "easting", pd.Series(dtype=float)).values

    @property
    def northing(self) -> np.ndarray:
        # north for single value of northing
        return self._frame.get(
            "northing", pd.Series(dtype=float)).values

    @property
    def elevations(self) -> np.ndarray:
        return self._frame.get(
            "elevation", pd.Series(dtype=float)).values
    @property
    def elevation (self):
        # use elev for single station and elevation
        # for pd.Series or array
        # of multiple elevations
        return self.elevations

    @property
    def latitude(self) -> np.ndarray:
        """Exposes the latitude data as a NumPy array.

        Returns
        -------
        numpy.ndarray
            An array of station latitudes, or an empty array if
            latitude data is not available in the frame.
        """
        return self._frame.get(
            "latitude", self._latitude).values

    @property
    def longitude(self) -> np.ndarray:
        """Exposes the longitude data as a NumPy array.

        Returns
        -------
        numpy.ndarray
            An array of station longitudes, or an empty array if
            longitude data is not available in the frame.
        """
        return self._frame.get(
            "longitude", self._longitude).values


def _to_float_series(s: pd.Series) -> pd.Series:
    """Coerce a Series to float, preserving NaN for non-numeric."""
    return pd.to_numeric(s, errors="coerce").astype(float)


def _median_inc(vals: np.ndarray) -> float | None:
    """Median increment if grid-like, else None."""
    if vals.size < 2:
        return None
    diffs = np.diff(vals)
    inc = np.median(diffs)
    if np.allclose(diffs, inc, rtol=0.0, atol=1e-6):
        return float(inc)
    return None


@dc(slots=True)
class Station(AVGComponentBase):
    """
    One-dimensional survey-line geometry container.

    This component reads the *station* coordinate from a tidy AVG frame
    (usually column ``'station'``), validates and stores unique positions,
    and (optionally) normalizes the origin to zero. It also derives a
    stable set of header keys ($Stn.*) for round-tripping modern headers.

    Key features
    -----------
    • Accepts either a DataFrame (from AVG parsers) or raw arrays.
    • Handles ragged frequency coverage (no equal-count requirement).
    • Tracks units; internal helpers can work in metres when desired.
    • Derives $Stn.Beg / $Stn.Inc / $Stn.Left / $Stn.Right if grid-like.
    • Exposes dictionary mappings useful for downstream selection.

    Notes
    -----
    We **do not** modify the caller's AVG frame; we keep our own view
    in ``self._frame`` that at minimum contains a numeric ``'station'``
    column and (if relevant) a ``'station_m'`` column when the input
    length unit is feet ('ft') and conversion to metres is helpful.
    """

    # configuration
    unit: str = "m"
    normalize: bool = False
    allow_ragged: bool = True

    # derived state
    values: np.ndarray = field(default_factory=lambda: np.empty(0))
    names: list[str] = field(default_factory=list)
    index_by_value: dict[float, np.ndarray] = field(default_factory=dict)
    index_by_name: dict[str, np.ndarray] = field(default_factory=dict)


    def read(
        self,
        source: pd.DataFrame | Sequence[float] | np.ndarray,
        meta: Mapping[str, Any] | None = None,
        *,
        unit: str | None = None,
        names: Sequence[str] | None = None,
        normalize: bool | None = None,
        allow_ragged: bool | None = None,
    ) -> None:
        """
        Populate the component from *source*.

        Parameters
        ----------
        source
            DataFrame with a ``'station'`` column **or** a 1-D array of
            station values (possibly repeated across frequencies).
        meta
            Optional file metadata. If provided and ``unit`` is None, we
            try ``meta['Unit.Length']`` for ('m'|'ft').
        unit
            Length unit of *source* values. Defaults to 'm' or to
            ``meta['Unit.Length']`` when available.
        names
            Optional labels for *unique* station positions. If omitted,
            we auto-generate: 'S00', 'S01', ...
        normalize
            When True, shift the origin so the first station is 0.0.
        allow_ragged
            When False, all stations must have identical row counts
            (one per frequency). When True (default), ragged coverage
            is allowed.
        """
        # resolve options
        munit = (meta or {}).get("Unit.Length") if meta else None
        unit = (unit or munit or self.unit or "m").lower()
        normalize = self.normalize if normalize is None else normalize
        allow_ragged = (
            self.allow_ragged if allow_ragged is None else allow_ragged
        )

        # build a minimal frame with a numeric 'station' column
        if isinstance(source, pd.DataFrame):
            source = find_and_rename_column(source.copy(), 'station')

            if "station" not in source.columns:
                raise StationError("column 'station' missing in frame")

            frame = source[["station"]].copy()
            frame["station"] = _to_float_series(frame["station"])
        else:
            arr = np.asarray(source, dtype=float).ravel()
            if arr.size == 0:
                raise StationError("empty station array")
            frame = pd.DataFrame({"station": arr})

        # drop rows where station cannot be parsed (all we can do)
        frame = frame.loc[frame["station"].notna()].reset_index(drop=True)
        if frame.empty:
            raise StationError("no valid station values")

        # normalize if requested
        if normalize:
            frame["station"] = frame["station"] - frame["station"].min()

        # optionally provide a metres view if unit is feet
        if unit == "ft":
            frame["station_m"] = frame["station"] / 3.280839895
        elif unit == "km":
            frame["station_m"] = frame["station"] * 1e3
        else:
            # mirror for symmetry; callers can ignore it
            frame["station_m"] = frame["station"]

        # sort by station coord for stable introspection
        frame.sort_values("station", kind="mergesort", inplace=True)
        frame.reset_index(drop=True, inplace=True)

        # unique positions and index maps
        uniq = frame["station"].dropna().unique()
        uniq.sort()
        self.values = uniq.astype(float)

        # build group index maps
        groups = frame.groupby("station", sort=True, dropna=False).indices
        idx_by_val: dict[float, np.ndarray] = {}
        for k, v in groups.items():
            try:
                key = float(k)
            except Exception:
                # groupby won't give NaN key here (we dropped NaNs)
                continue
            idx_by_val[key] = np.asarray(sorted(v), dtype=int)
        self.index_by_value = idx_by_val

        # ragged check if requested
        if not allow_ragged and self.values.size > 0:
            counts = np.array([len(idx_by_val[val]) for val in self.values])
            if not np.all(counts == counts[0]):
                raise StationError("inconsistent rows per station")

        # station names (one per unique position)
        if names is not None:
            if len(names) != self.values.size:
                raise StationError(
                    "`names` length must match unique station count"
                )
            self.names = list(names)
        else:
            ids, _ = number_stations(self.values.size, 1, prefix="S")
            self.names = ids

        # name→index map (same row indices as value map)
        self.index_by_name = {
            nm: self.index_by_value[val]
            for nm, val in zip(self.names, self.values, strict=True)
        }

        # stash the frame + a small meta view in the component payload
        self._frame = frame
        self._meta = dict(meta or {})
        # remember configuration for write/pretty
        self.unit = unit
        self.normalize = normalize
        self.allow_ragged = allow_ragged

    def write(self) -> Sequence[str]:
        """
        Emit a small, human-friendly header preamble plus a CSV of the
        station coordinate. This is mainly useful for diagnostics and
        tests; full file writing should be orchestrated at a higher
        level.
        """
        if self._frame.empty:
            return []

        # Try to derive canonical $Stn.* keys for the header
        keys = self.to_keywords()

        lines: list[str] = []
        lines.append(r"\ Station Geometry")
        for k, v in keys.items():
            lines.append(f"${k}={v}")
        lines.append("")  # spacer before table

        # CSV of station and metres view (if present)
        cols = ["station"] + (
            ["station_m"] if "station_m" in self._frame.columns else [])
        csv = self._frame.loc[:, cols].to_csv(index=False, float_format="%.6g")
        lines.extend(csv.splitlines())
        return lines

    @property
    def n_unique(self) -> int:
        """Number of unique station positions."""
        return int(self.values.size)

    @property
    def span(self) -> tuple[float, float] | None:
        """(min, max) of the station coordinate; None if empty."""
        if self.values.size == 0:
            return None
        return (float(self.values.min()), float(self.values.max()))

    @property
    def increment(self) -> float | None:
        """Grid increment if evenly spaced, else None."""
        return _median_inc(self.values)

    def label_map(self) -> dict[float, str]:
        """Map station numeric value → generated label (e.g., 'S03')."""
        return {float(v): n for v, n in zip(
            self.values, self.names, strict=True)}

    def to_keywords(self) -> dict[str, Any]:
        """
        Derive a stable set of $Stn.* keys from current geometry:

        - Stn.Beg   : first station
        - Stn.Inc   : increment (when grid-like)
        - Stn.Left  : left edge (min)
        - Stn.Right : right edge (max)

        For convenience, we also mirror GDP-station versions when
        increment is grid-like: Stn.GdpBeg / Stn.GdpInc.
        """
        out: dict[str, Any] = {}
        if self.values.size == 0:
            return out

        left, right = float(self.values.min()), float(self.values.max())
        out["Stn.Left"] = left
        out["Stn.Right"] = right
        out["Stn.Beg"] = left

        inc = self.increment
        if inc is not None:
            out["Stn.Inc"] = inc
            out["Stn.GdpBeg"] = left
            out["Stn.GdpInc"] = inc
        return out


    def __str__(self) -> str:
        if self.values.size == 0:
            return "Station(empty)"
        span = self.span
        inc = self.increment
        span_txt = f"{span[0]:g}–{span[1]:g}" if span else "?"
        inc_txt = f"{inc:g}" if inc is not None else "ragged"
        return (
            f"Station(n={self.n_unique}, span={span_txt} {self.unit}, "
            f"inc={inc_txt})"
        )

    __repr__ = __str__
