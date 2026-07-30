# Author: L. Kouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import math
from collections.abc import Iterable
from typing import Any, Union

import numpy as np

from .api.property import PyCSAMTObject
from .api.typing import ArrayLike
from .exceptions import LocationError
from .gis.utils import (
    assert_lat_value,
    assert_lon_value,
    ll_to_utm,
    utm_to_ll,
)
from .log.logger import get_logger
from .utils.arrayops import assert_xy_in, is_iterable
from .utils.validation import _check_consistency_size

logger = get_logger(__name__)

__all__ = ["Location"]


_Number = Union[float, int, str]


__all__ = ["BaseLoc", "Location", "UTMZone", "Bounds", "GeoPath"]


class BaseLoc(PyCSAMTObject):
    """
    Base mixin for point-like location objects with lat/lon and
    UTM representations.

    This class does not store data by itself. Subclasses are
    expected to define the following private attributes:

    - _latitude : float or None
    - _longitude : float or None
    - _elevation : float
    - _utm_zone : str or None
    - _easting : float or None
    - _northing : float or None

    It provides robust validation, conversion, vectorized
    helpers, and convenience operations such as distance and
    approximate equality.

    Notes
    -----
    - Latitude and longitude are validated and coerced using
      :func:`assert_lat_value` and :func:`assert_lon_value`.
    - UTM transformations are delegated to
      :func:`ll_to_utm` and :func:`utm_to_ll`.
    - Accessors ``has_latlon`` and ``has_utm`` report whether
      the minimal fields for each representation are present.
    - Methods prefixed with ``_require_*`` assert presence and
      return normalized values or raise :class:`LocationError`.

    Examples
    --------
    >>> class Station(BaseLoc):
    ...     def __init__(self, lat, lon):
    ...         self._latitude = float(lat)
    ...         self._longitude = float(lon)
    ...         self._elevation = 0.0
    ...         self._utm_zone = None
    ...         self._easting = None
    ...         self._northing = None
    >>> s = Station(10.0, 5.0)
    >>> s.has_latlon
    True
    >>> s.ensure_utm()  # doctest: +ELLIPSIS
    ('31N', ..., ...)
    >>> s2 = Station(10.0001, 5.0001)
    >>> round(s.distance_to(s2), 1) >= 0.0
    True

    See Also
    --------
    Location
        Concrete implementation that exposes public properties
        for latitude, longitude, elevation, and UTM fields.
    ll_to_utm, utm_to_ll
        Coordinate converters.
    assert_lat_value, assert_lon_value
        Validators for latitude and longitude.

    References
    ----------
    .. [1] Snyder, J. P. 1987. Map Projections, A Working
           Manual. USGS Professional Paper 1395.
    .. [2] NGA. 2014. The Universal Grids and the Transverse
           Mercator and Polar Stereographic Map Projections.
    """

    __slots__ = ()
    # expects subclass attrs:
    # _latitude,_longitude,_elevation,_utm_zone,_easting,_northing

    @property
    def has_latlon(self) -> bool:
        r"""
        bool: True if both latitude and longitude are present.

        Notes
        -----
        This does not validate numerical ranges. Use
        ``_require_latlon`` to retrieve normalized values.
        """
        la = getattr(self, "_latitude", None)
        lo = getattr(self, "_longitude", None)
        return la is not None and lo is not None

    @property
    def has_utm(self) -> bool:
        r"""
        bool: True if UTM zone, easting, and northing are present.

        Notes
        -----
        The zone string is not normalized here. Use
        ``_require_utm`` to normalize and validate.
        """

        z = getattr(self, "_utm_zone", None)
        e = getattr(self, "_easting", None)
        n = getattr(self, "_northing", None)
        return bool(z) and e is not None and n is not None

    def _require_latlon(self) -> tuple[float, float]:
        r"""
        Return normalized latitude and longitude.

        Returns
        -------
        tuple of float
            ``(latitude, longitude)`` in decimal degrees.

        Raises
        ------
        LocationError
            If either component is missing or invalid.

        Notes
        -----
        Validation and coercion are delegated to
        :func:`assert_lat_value` and :func:`assert_lon_value`.
        """

        la = getattr(self, "_latitude", None)
        lo = getattr(self, "_longitude", None)
        if la is None or lo is None:
            raise LocationError("lat/lon not set")
        return _coerce_latlon(la, lo)

    def _require_utm(self) -> tuple[str, float, float]:
        r"""
        Return normalized UTM triple.

        Returns
        -------
        tuple
            ``(zone, easting, northing)`` where ``zone`` is a
            normalized ``"NNH"`` string and coordinates are floats.

        Raises
        ------
        LocationError
            If any component is missing or cannot be coerced.

        Notes
        -----
        The zone string is normalized by ``_norm_zone``.
        Numeric fields are coerced by ``_as_float``.
        """

        z = getattr(self, "_utm_zone", None)
        e = getattr(self, "_easting", None)
        n = getattr(self, "_northing", None)
        if not z or e is None or n is None:
            raise LocationError("UTM not set")
        z = _norm_zone(z)
        e = _as_float(e, "easting")
        n = _as_float(n, "northing")
        return z, float(e), float(n)

    def _set_latlon(self, la: Any, lo: Any) -> None:
        r"""
        Set latitude and longitude after coercion and validation.

        Parameters
        ----------
        la : Any
            Latitude value or string parsable by
            :func:`assert_lat_value`.
        lo : Any
            Longitude value or string parsable by
            :func:`assert_lon_value`.

        Raises
        ------
        LocationError
            If inputs are missing or invalid.
        """

        la_f, lo_f = _coerce_latlon(la, lo)
        self._latitude = la_f
        self._longitude = lo_f

    def _set_utm(
        self,
        zone: Any,
        easting: Any,
        northing: Any,
    ) -> None:
        r"""
        Set UTM components after normalization and coercion.

        Parameters
        ----------
        zone : Any
            Zone designator. Integers map to hemisphere by sign.
            Strings are normalized to ``"NNH"`` via ``_norm_zone``.
        easting : Any
            Easting in meters.
        northing : Any
            Northing in meters.

        Raises
        ------
        LocationError
            If zone is invalid or coordinates cannot be coerced.
        """

        z = _norm_zone(zone)
        e = _as_float(easting, "easting")
        n = _as_float(northing, "northing")
        self._utm_zone = z
        self._easting = float(e)
        self._northing = float(n)

    def _clear_ll(self) -> None:
        r"""
        Clear latitude and longitude fields by setting them to
        ``None``.
        """

        self._latitude = None
        self._longitude = None

    def _clear_utm(self) -> None:
        r"""
        Clear UTM fields by setting zone, easting, and northing to
        ``None``.
        """

        self._utm_zone = None
        self._easting = None
        self._northing = None

    def ensure_utm(self) -> tuple[str, float, float]:
        r"""
        Ensure UTM is available, converting from lat/lon if needed.

        Returns
        -------
        tuple
            ``(zone, easting, northing)`` as normalized values.

        Raises
        ------
        LocationError
            If neither lat/lon nor UTM is sufficiently defined.

        Notes
        -----
        - If UTM is already present it is returned unchanged.
        - Otherwise lat/lon are required and converted via
          :func:`ll_to_utm` with ``reference_ellipsoid=23``.

        Examples
        --------
        >>> class P(BaseLoc):
        ...     def __init__(self, la, lo):
        ...         self._latitude = la
        ...         self._longitude = lo
        ...         self._elevation = 0.0
        ...         self._utm_zone = None
        ...         self._easting = None
        ...         self._northing = None
        >>> P(10.0, 5.0).ensure_utm()[0]
        '31N'
        """

        if self.has_utm:
            return self._require_utm()
        la, lo = self._require_latlon()
        z, e, n = ll_to_utm(
            reference_ellipsoid=23,
            lat=la,
            lon=lo,
        )
        self._set_utm(z, e, n)
        return self._require_utm()

    def ensure_latlon(self) -> tuple[float, float]:
        r"""
        Ensure lat/lon are available, converting from UTM if needed.

        Returns
        -------
        tuple of float
            ``(latitude, longitude)`` in decimal degrees.

        Raises
        ------
        LocationError
            If neither UTM nor lat/lon is sufficiently defined.

        Notes
        -----
        - If lat/lon are already present they are returned.
        - Otherwise UTM must be present and is converted via
          :func:`utm_to_ll` with ``reference_ellipsoid=23``.
        """

        if self.has_latlon:
            return self._require_latlon()
        z, e, n = self._require_utm()
        la, lo = utm_to_ll(
            reference_ellipsoid=23,
            northing=n,
            easting=e,
            zone=z,
        )
        self._set_latlon(la, lo)
        return self._require_latlon()

    @staticmethod
    def to_utm_batch(
        lats: Iterable[Any],
        lons: Iterable[Any],
    ) -> tuple[np.ndarray, np.ndarray, list[str]]:
        r"""
        Vectorized conversion from lat/lon arrays to UTM arrays.

        Parameters
        ----------
        lats : Iterable[Any]
            Iterable of latitude values or parsable strings.
        lons : Iterable[Any]
            Iterable of longitude values or parsable strings.

        Returns
        -------
        tuple
            ``(eastings, northings, zones)`` where the first two are
            ``numpy.ndarray`` and zones is a ``list[str]``.

        Raises
        ------
        LocationError
            If sizes mismatch or any item cannot be coerced.

        Notes
        -----
        - Each input pair is validated with
          :func:`assert_lat_value` and :func:`assert_lon_value`.
        - Per-point conversion uses :func:`ll_to_utm`.

        See Also
        --------
        to_latlon_batch
            Inverse vectorized conversion.
        """

        lats = np.asarray(list(lats), dtype=object)
        lons = np.asarray(list(lons), dtype=object)
        if lats.size != lons.size:
            raise LocationError("size mismatch")
        es = np.zeros(lats.size, dtype=float)
        ns = np.zeros(lats.size, dtype=float)
        zs: list[str] = []
        for i, (la, lo) in enumerate(zip(lats, lons)):
            la_f, lo_f = _coerce_latlon(la, lo)
            z, e, n = ll_to_utm(
                reference_ellipsoid=23,
                lat=la_f,
                lon=lo_f,
            )
            zs.append(str(z))
            es[i] = float(e)
            ns[i] = float(n)
        return es, ns, zs

    @staticmethod
    def to_latlon_batch(
        easts: Iterable[Any],
        norths: Iterable[Any],
        zones: Iterable[Any],
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Vectorized conversion from UTM arrays to lat/lon arrays.

        Parameters
        ----------
        easts : Iterable[Any]
            Easting values in meters.
        norths : Iterable[Any]
            Northing values in meters.
        zones : Iterable[Any]
            Zone designators, normalized by ``_norm_zone``.

        Returns
        -------
        tuple
            ``(lats, lons)`` as ``numpy.ndarray`` of floats.

        Raises
        ------
        LocationError
            If sizes mismatch or zone normalization fails.

        Notes
        -----
        Per-point conversion uses :func:`utm_to_ll`.
        """

        easts = np.asarray(list(easts), dtype=float)
        norths = np.asarray(list(norths), dtype=float)
        zones = list(zones)
        if easts.size != norths.size or easts.size != len(zones):
            raise LocationError("size mismatch")
        las = np.zeros_like(easts, dtype=float)
        los = np.zeros_like(easts, dtype=float)
        for i, (e, n, z) in enumerate(zip(easts, norths, zones)):
            z_s = _norm_zone(z)
            la, lo = utm_to_ll(
                reference_ellipsoid=23,
                northing=float(n),
                easting=float(e),
                zone=z_s,
            )
            las[i] = float(la)
            los[i] = float(lo)
        return las, los

    def distance_to(self, other: Any) -> float:
        r"""
        Great-circle distance between this object and another.

        Parameters
        ----------
        other : Any
            Another point-like object. If it exposes ``latitude``
            and ``longitude`` attributes, they are used. If it has
            ``ensure_latlon()``, that method is called.

        Returns
        -------
        float
            Distance in meters computed via the haversine formula.

        Raises
        ------
        LocationError
            If the other object cannot provide lat/lon.

        Notes
        -----
        The haversine distance is given by

        .. math::

           d = 2 R \\arcsin\\Big(\\sqrt{a}\\Big)

        with

        .. math::

           a = \\sin^2(\\Delta\\varphi/2) +
               \\cos\\varphi_1\\cos\\varphi_2\\sin^2(\\Delta\\lambda/2)

        where angles are in radians and :math:`R` is the mean Earth
        radius.

        References
        ----------
        .. [1] Sinnott, R. W. 1984. Virtues of the Haversine.
               Sky and Telescope 68(2), 159.
        """

        la1, lo1 = self.ensure_latlon()
        la2 = getattr(other, "latitude", None)
        lo2 = getattr(other, "longitude", None)
        if la2 is None or lo2 is None:
            if hasattr(other, "ensure_latlon"):
                la2, lo2 = other.ensure_latlon()  # type: ignore
            else:
                raise LocationError("other has no lat/lon")
        la2_f, lo2_f = _coerce_latlon(la2, lo2)
        return _haversine_m(
            float(la1),
            float(lo1),
            float(la2_f),
            float(lo2_f),
        )

    def nearly_equals(
        self,
        other: Any,
        *,
        tol_m: float = 0.25,
    ) -> bool:
        r"""
        Approximate equality based on geodesic tolerance.

        Parameters
        ----------
        other : Any
            Another point-like object.
        tol_m : float, optional
            Distance tolerance in meters. Default is 0.25 m.

        Returns
        -------
        bool
            True if ``distance_to(other) <= tol_m``.

        Notes
        -----
        This is useful to compare coordinates that differ by small
        rounding errors or mixed conversions.
        """

        try:
            d = self.distance_to(other)
        except LocationError:
            return False
        return d <= tol_m

    def to_dict(self) -> dict[str, Any]:
        r"""
        Serialize the object to a dictionary of primitive values.

        Returns
        -------
        dict
            Dictionary with keys
            ``latitude``, ``longitude``, ``elevation``,
            ``utm_zone``, ``easting``, ``northing``.

        Notes
        -----
        Fields may be ``None`` if not set.
        The result is suitable for JSON or YAML dumps.
        """
        out: dict[str, Any] = {}
        la = getattr(self, "_latitude", None)
        lo = getattr(self, "_longitude", None)
        el = getattr(self, "_elevation", None)
        z = getattr(self, "_utm_zone", None)
        e = getattr(self, "_easting", None)
        n = getattr(self, "_northing", None)
        out.update(
            latitude=la,
            longitude=lo,
            elevation=el,
            utm_zone=z,
            easting=e,
            northing=n,
        )
        return out

    @classmethod
    def from_dict(cls, d: dict[str, Any]):
        r"""
        Construct an instance from a dictionary.

        Parameters
        ----------
        d : dict
            Mapping with optional keys
            ``latitude``, ``longitude``, ``elevation``,
            ``utm_zone``, ``easting``, ``northing``.

        Returns
        -------
        BaseLoc
            Instance of ``cls`` with fields set.

        Notes
        -----
        - Missing fields default to ``None`` or ``0.0`` for
          elevation.
        - When both lat/lon and UTM are given, both are stored
          without cross-validation. Use ``ensure_*`` to sync.
        """

        obj = cls.__new__(cls)  # type: ignore[misc]
        obj._latitude = None
        obj._longitude = None
        obj._elevation = float(d.get("elevation", 0.0))
        obj._utm_zone = None
        obj._easting = None
        obj._northing = None
        la = d.get("latitude", None)
        lo = d.get("longitude", None)
        z = d.get("utm_zone", None)
        e = d.get("easting", None)
        n = d.get("northing", None)
        if la is not None and lo is not None:
            obj._set_latlon(la, lo)
        if z is not None and e is not None and n is not None:
            obj._set_utm(z, e, n)
        return obj

    def __eq__(self, other: object) -> bool:
        r"""
        Semantic equality for location objects.

        Parameters
        ----------
        other : object
            Object to compare against.

        Returns
        -------
        bool
            - If both have lat/lon, compare rounded to 7 decimals.
            - Else if both have UTM, compare zone and rounded meter
              integers.
            - Else fall back to ``nearly_equals`` with a very small
              tolerance.

        Notes
        -----
        Returns ``NotImplemented`` for non-BaseLoc objects, which
        lets Python try reversed comparisons.
        """

        if not isinstance(other, BaseLoc):
            return NotImplemented
        if self.has_latlon and other.has_latlon:
            la1, lo1 = self._require_latlon()
            la2, lo2 = other._require_latlon()
            return round(la1, 7) == round(la2, 7) and round(lo1, 7) == round(
                lo2, 7
            )
        if self.has_utm and other.has_utm:
            z1, e1, n1 = self._require_utm()
            z2, e2, n2 = other._require_utm()
            return (
                z1 == z2
                and int(round(e1)) == int(round(e2))
                and int(round(n1)) == int(round(n2))
            )
        try:
            return self.nearly_equals(  # type: ignore[arg-type]
                other,
                tol_m=0.01,
            )
        except Exception:
            return False

    def __hash__(self) -> int:
        r"""
        Hash compatible with equality semantics.

        Returns
        -------
        int
            - If lat/lon are present, hash rounded lat/lon.
            - Else if UTM is present, hash zone and rounded meters.
            - Else hash a unique placeholder with object id.

        Notes
        -----
        This enables use in sets and as dict keys when locations
        are effectively immutable.
        """

        if self.has_latlon:
            la, lo = self._require_latlon()
            return hash((round(la, 7), round(lo, 7)))
        if self.has_utm:
            z, e, n = self._require_utm()
            return hash((z, int(round(e)), int(round(n))))
        return hash(("empty", id(self)))


class Location(BaseLoc):
    r"""
    Concrete point-like location with lat/lon and UTM fields.

    This class stores one station or point in two coordinate
    systems: geographic (latitude, longitude) and UTM
    (zone, easting, northing). It inherits robust validation,
    conversion, and small utilities from :class:`BaseLoc`.

    Initialization favors lat/lon. If both lat/lon are given,
    UTM fields passed at the same time are ignored and UTM is
    computed from lat/lon. If only UTM is provided, a valid
    UTM zone string is required.

    Notes
    -----
    - Validation for latitude and longitude is delegated to
      :func:`assert_lat_value` and :func:`assert_lon_value`.
      Any error is raised as :class:`LocationError`.
    - Conversions use :func:`ll_to_utm` and :func:`utm_to_ll`
      with ``reference_ellipsoid=23``.
    - The zone string is normalized to ``"NNH"`` by
      the internal helper ``_norm_zone``.
    - Use :meth:`to_utm` and :meth:`to_latlon` to perform
      on-demand conversions. These methods call the base
      ``ensure_*`` utilities and update cached fields.

    Examples
    --------
    Create from lat/lon and compute UTM:

    >>> loc = Location(latitude=10.0, longitude=5.0)
    >>> z, e, n = loc.to_utm()  # doctest: +ELLIPSIS
    >>> isinstance(z, str) and isinstance(e, float)
    True

    Create from UTM:

    >>> loc = Location(utm_zone="31N", easting=500000.0, northing=4649776.0)
    >>> la, lo = loc.to_latlon()  # doctest: +ELLIPSIS
    >>> isinstance(la, float) and isinstance(lo, float)
    True

    See Also
    --------
    BaseLoc
        Mixin offering validation, conversion, batch helpers.
    UTMZone
        Zone parser and formatter.
    Bounds, GeoPath
        Light geometry helpers useful around locations.

    References
    ----------
    .. [1] Snyder, J. P. 1987. Map Projections, A Working
           Manual. USGS Professional Paper 1395.
    """

    def __init__(
        self,
        latitude: _Number | None = None,
        longitude: _Number | None = None,
        *,
        elevation: _Number = 0.0,
        utm_zone: str | None = None,
        easting: _Number | None = None,
        northing: _Number | None = None,
    ):
        self._latitude: float | None = None
        self._longitude: float | None = None
        self._elevation: float = 0.0
        self._utm_zone: str | None = None
        self._easting: float | None = None
        self._northing: float | None = None

        self.elevation = elevation

        if utm_zone is not None:
            # normalize early to keep invariants
            self.utm_zone = utm_zone

        if latitude is not None and longitude is not None:
            self.latitude = latitude
            self.longitude = longitude
            if easting is not None or northing is not None:
                logger.warning(
                    "Ignoring provided UTM easting/northing; "
                    "lat/lon takes precedence."
                )
            # keep old behavior: auto-compute UTM
            self.to_utm()
        elif easting is not None and northing is not None:
            if self._utm_zone is None:
                raise LocationError(
                    "utm_zone is required when initializing "
                    "from UTM easting/northing."
                )
            self._set_utm(self._utm_zone, easting, northing)
        else:
            logger.debug(
                "Initialized empty Location; set coordinates "
                "before converting."
            )

    @property
    def latitude(self) -> float | None:
        return self._latitude

    @latitude.setter
    def latitude(self, val: _Number) -> None:
        r"""
        float or None : Geographic latitude in decimal degrees.

        Setting this property validates and coerces the input using
        :func:`assert_lat_value`. Any error is raised as
        :class:`LocationError`.

        Notes
        -----
        - Accepts float-like numbers or sexagesimal strings as
          supported by the validator.
        - Setting latitude does not automatically update UTM.
          Call :meth:`to_utm` or :meth:`ensure_utm` when needed.

        Examples
        --------
        >>> loc = Location()
        >>> loc.latitude = 12.5
        >>> loc.latitude
        12.5
        """

        try:
            v = assert_lat_value(val)
        except (TypeError, ValueError) as exc:
            raise LocationError(f"Invalid latitude: {val!r}") from exc
        if v is None:
            raise LocationError("Invalid latitude: None")
        self._latitude = float(v)

    @property
    def longitude(self) -> float | None:
        return self._longitude

    @longitude.setter
    def longitude(self, val: _Number) -> None:
        r"""
        float or None : Geographic longitude in decimal degrees.

        Setting this property validates and coerces the input using
        :func:`assert_lon_value`. Any error is raised as
        :class:`LocationError`.

        Notes
        -----
        - Accepts float-like numbers or sexagesimal strings as
          supported by the validator.
        - Setting longitude does not automatically update UTM.
          Call :meth:`to_utm` or :meth:`ensure_utm` when needed.

        Examples
        --------
        >>> loc = Location()
        >>> loc.longitude = -3.25
        >>> loc.longitude
        -3.25
        """

        try:
            v = assert_lon_value(val)
        except (TypeError, ValueError) as exc:
            raise LocationError(f"Invalid longitude: {val!r}") from exc
        if v is None:
            raise LocationError("Invalid longitude: None")
        self._longitude = float(v)

    @property
    def elevation(self) -> float:
        return self._elevation

    @elevation.setter
    def elevation(self, val: _Number) -> None:
        r"""
        float : Elevation in meters above reference datum.

        The value is coerced to float. Invalid inputs raise
        :class:`LocationError`.

        Notes
        -----
        Elevation is not involved in UTM or lat/lon conversions,
        but is stored alongside the position for convenience.
        """

        self._elevation = _as_float(val, "elevation")

    @property
    def utm_zone(self) -> str | None:
        return self._utm_zone

    @utm_zone.setter
    def utm_zone(self, val: str) -> None:
        r"""
        str or None : UTM zone designator.

        Setting this property normalizes the value to ``"NNH"``,
        where ``NN`` is the zone number in ``[1, 60]`` and
        ``H`` is the hemisphere letter ``"N"`` or ``"S"``.
        Invalid inputs raise :class:`LocationError`.

        Notes
        -----
        Normalization is performed by the internal ``_norm_zone``.
        """

        self._utm_zone = _norm_zone(val)

    @property
    def easting(self) -> float | None:
        r"""
        float or None : UTM easting in meters.

        Notes
        -----
        This is a read-only view of the internal field. Use
        :meth:`to_utm` or :meth:`ensure_utm` to populate it.
        """

        return self._easting

    @property
    def northing(self) -> float | None:
        r"""
        float or None : UTM northing in meters.

        Notes
        -----
        This is a read-only view of the internal field. Use
        :meth:`to_utm` or :meth:`ensure_utm` to populate it.
        """
        return self._northing

    def to_utm(self) -> tuple[str, float, float]:
        r"""
        Convert to UTM and return the normalized triple.

        Returns
        -------
        tuple
            ``(zone, easting, northing)`` where zone is ``"NNH"`` and
            coordinates are floats.

        Raises
        ------
        LocationError
            If geographic coordinates are not set or invalid.

        Notes
        -----
        - Calls :meth:`BaseLoc.ensure_utm`. If UTM is already
          present, it is returned. Otherwise, lat/lon are required.
        - Conversion uses :func:`ll_to_utm` with
          ``reference_ellipsoid=23``.

        Examples
        --------
        >>> loc = Location(latitude=10.0, longitude=5.0)
        >>> z, e, n = loc.to_utm()  # doctest: +ELLIPSIS
        >>> isinstance(z, str) and isinstance(e, float)
        True
        """

        z, e, n = self.ensure_utm()
        logger.debug(
            "Converted to UTM zone=%s, E=%.6f, N=%.6f",
            z,
            e,
            n,
        )
        return z, e, n

    def to_latlon(self) -> tuple[float, float]:
        r"""
        Convert to geographic lat/lon and return the pair.

        Returns
        -------
        tuple of float
            ``(latitude, longitude)`` in decimal degrees.

        Raises
        ------
        LocationError
            If UTM fields are not set or invalid.

        Notes
        -----
        - Calls :meth:`BaseLoc.ensure_latlon`. If lat/lon are
          already present, they are returned. Otherwise, UTM must
          be present.
        - Conversion uses :func:`utm_to_ll` with
          ``reference_ellipsoid=23``.

        Examples
        --------
        >>> loc = Location(
        ...     utm_zone="31N", easting=500000.0, northing=4649776.0
        ... )
        >>> la, lo = loc.to_latlon()  # doctest: +ELLIPSIS
        >>> isinstance(la, float) and isinstance(lo, float)
        True
        """

        la, lo = self.ensure_latlon()
        logger.debug(
            "Converted to lat=%.8f, lon=%.8f",
            la,
            lo,
        )
        return la, lo

    @staticmethod
    def to_utm_in(
        lats: ArrayLike,
        lons: ArrayLike,
        *,
        data=None,
        utm_zone: str | None = None,
        datum: str | None = None,
        **kws,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Vectorized conversion of lat/lon arrays to UTM arrays.

        Parameters
        ----------
        lats : ArrayLike
            Iterable or column name for latitudes. When a string is
            given, ``data`` must be provided and will be used with
            :func:`assert_xy_in`.
        lons : ArrayLike
            Iterable or column name for longitudes. When a string is
            given, ``data`` must be provided and will be used with
            :func:`assert_xy_in`.
        data : Any, optional
            Tabular container used when ``lats`` or ``lons`` are
            column names.
        utm_zone : str, optional
            Requested zone for consistency checking. A warning is
            logged if a point falls in a different zone.
        datum : str, optional
            Kept for backward compatibility. Not used.
        **kws :
            Accepted for compatibility. Ignored.

        Returns
        -------
        tuple
            ``(eastings, northings)`` as ``numpy.ndarray`` of float.

        Raises
        ------
        TypeError
            If ``data`` is missing while inputs are column names,
            or if arrays are not the same size.
        LocationError
            If an element cannot be validated or coerced.

        Notes
        -----
        - Each point is validated via the scalar validators and
          converted with :func:`ll_to_utm`.
        - For batch conversion that also returns zones, see
          :meth:`BaseLoc.to_utm_batch`.

        See Also
        --------
        BaseLoc.to_utm_batch
            Returns eastings, northings, and zones.

        References
        ----------
        .. [1] Snyder, J. P. 1987. Map Projections, A Working
               Manual. USGS Professional Paper 1395.
        """

        # Keep signature for compatibility. 'datum' unused.
        if (isinstance(lats, str) or isinstance(lons, str)) and data is None:
            raise TypeError(
                "Data can't be None when latitude or longitude "
                "is a string (column name)."
            )
        if data is not None:
            lats, lons = assert_xy_in(lats, lons, data=data)
        if lats is None or lons is None:
            raise TypeError("Both longitude and latitude are required.")

        lats = is_iterable(lats, exclude_string=True, transform=True)
        lons = is_iterable(lons, exclude_string=True, transform=True)
        _check_consistency_size(lats, lons)

        lats = np.asarray(lats, dtype=object)
        lons = np.asarray(lons, dtype=object)

        easts = np.zeros(lats.size, dtype=float)
        norths = np.zeros(lats.size, dtype=float)

        for i, (la, lo) in enumerate(zip(lats, lons)):
            vla = assert_lat_value(la)
            vlo = assert_lon_value(lo)
            if vla is None or vlo is None:
                raise LocationError(f"Invalid lat/lon at index {i}")
            z, e, n = ll_to_utm(
                reference_ellipsoid=23,
                lat=float(vla),
                lon=float(vlo),
            )
            if utm_zone is not None and str(z) != str(utm_zone):
                logger.warning(
                    "Point %d in zone %s differs from requested %s.",
                    i,
                    z,
                    utm_zone,
                )
            easts[i] = float(e)
            norths[i] = float(n)

        return easts, norths

    @staticmethod
    def to_latlon_in(
        easts: ArrayLike,
        norths: ArrayLike,
        *,
        data=None,
        utm_zone: str | None = None,
        datum: str | None = None,
        **kws,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Vectorized conversion of UTM arrays to lat/lon arrays.

        Parameters
        ----------
        easts : ArrayLike
            Easting values in meters, or a column name if ``data``
            is provided.
        norths : ArrayLike
            Northing values in meters, or a column name if
            ``data`` is provided.
        data : Any, optional
            Tabular container used when inputs are column names.
        utm_zone : str
            UTM zone designator applied to all input rows.
        datum : str, optional
            Kept for backward compatibility. Not used.
        **kws :
            Accepted for compatibility. Ignored.

        Returns
        -------
        tuple
            ``(lats, lons)`` as ``numpy.ndarray`` of float.

        Raises
        ------
        TypeError
            If ``data`` is missing while inputs are column names,
            or if arrays are not the same size.
        LocationError
            If ``utm_zone`` is missing or cannot be normalized.

        Notes
        -----
        - The zone string is normalized once and applied to all
          points.
        - Per-point conversion uses :func:`utm_to_ll`.
        - For batch conversion with per-row zones, use
          :meth:`BaseLoc.to_latlon_batch`.

        See Also
        --------
        BaseLoc.to_latlon_batch
            Inverse vectorized conversion that accepts per-row
            zones.

        References
        ----------
        .. [1] Snyder, J. P. 1987. Map Projections, A Working
               Manual. USGS Professional Paper 1395.
        """

        # Keep signature for compatibility. 'datum' unused.
        if (
            isinstance(easts, str) or isinstance(norths, str)
        ) and data is None:
            raise TypeError(
                "Data can't be None when easting or northing "
                "is a string (column name)."
            )
        if data is not None:
            easts, norths = assert_xy_in(easts, norths, data=data)
        if easts is None or norths is None:
            raise TypeError("Both easting and northing are required.")

        easts = is_iterable(easts, exclude_string=True, transform=True)
        norths = is_iterable(norths, exclude_string=True, transform=True)
        _check_consistency_size(easts, norths)

        easts = np.asarray(easts, dtype=float)
        norths = np.asarray(norths, dtype=float)

        if utm_zone is None:
            raise LocationError(
                "utm_zone is required to convert UTM array to lat/lon."
            )
        z = _norm_zone(utm_zone)

        lats: list[float] = []
        lons: list[float] = []

        for e, n in zip(easts, norths):
            la, lo = utm_to_ll(
                reference_ellipsoid=23,
                northing=float(n),
                easting=float(e),
                zone=z,
            )
            lats.append(float(la))
            lons.append(float(lo))

        return np.asarray(lats), np.asarray(lons)

    def __repr__(self) -> str:
        """
        Return a concise representation for debugging.

        Returns
        -------
        str
            A one-line string with current field values.

        Notes
        -----
        The representation is not a strict round-trip and may be
        changed in future for clarity.
        """

        return (
            "Location("
            f"lat={self._latitude}, "
            f"lon={self._longitude}, "
            f"elev={self._elevation}, "
            f"utm_zone={self._utm_zone}, "
            f"easting={self._easting}, "
            f"northing={self._northing}"
            ")"
        )


class UTMZone(PyCSAMTObject):
    r"""
    Lightweight UTM zone parser and formatter.

    This class normalizes a zone designator to the canonical
    ``"NNH"`` form and exposes the numeric zone and hemisphere.

    Notes
    -----
    - Normalization is delegated to the internal helper
      ``_norm_zone`` which enforces ``1 <= NN <= 60`` and
      hemisphere in ``{"N","S"}``.
    - Instances are immutable in practice; use :meth:`as_str`
      to render the canonical string.

    Examples
    --------
    >>> UTMZone("31n").as_str()
    '31N'
    >>> UTMZone("-31").as_str()  # negative implies south
    '31S'

    See Also
    --------
    Location
        Concrete point holding both lat/lon and UTM values.

    References
    ----------
    .. [1] NGA. 2014. The Universal Grids and the Transverse
           Mercator and Polar Stereographic Map Projections.
    """

    def __init__(self, zone: str):
        s = _norm_zone(zone)
        self._num = int(s[:-1])
        self._hem = s[-1]

    @property
    def number(self) -> int:
        """int : Zone number in ``[1, 60]``."""
        return int(self._num)

    @property
    def hemisphere(self) -> str:
        """str : Hemisphere letter, either ``"N"`` or ``"S"``."""
        return str(self._hem)

    def as_str(self) -> str:
        r"""
        Return the canonical ``"NNH"`` string.

        Returns
        -------
        str
            Canonical zone string composed of the zone number and
            hemisphere letter.

        Examples
        --------
        >>> UTMZone("31n").as_str()
        '31N'
        """

        return f"{self.number}{self.hemisphere}"

    @classmethod
    def from_latlon(cls, lat: float, lon: float) -> UTMZone:
        r"""
        Infer a UTM zone from geographic coordinates.

        Parameters
        ----------
        lat : float
            Latitude in decimal degrees.
        lon : float
            Longitude in decimal degrees.

        Returns
        -------
        UTMZone
            New instance for the inferred zone.

        Notes
        -----
        - Latitude and longitude are validated via ``_lat_ok`` and
          ``_lon_ok``. Errors are raised as :class:`LocationError`.
        - The zone is obtained via :func:`ll_to_utm` using
          ``reference_ellipsoid=23``.

        Examples
        --------
        >>> UTMZone.from_latlon(10.0, 5.0).as_str()
        '31N'
        """

        lat = _lat_ok(lat, "latitude")
        lon = _lon_ok(lon, "longitude")
        z, _, _ = ll_to_utm(reference_ellipsoid=23, lat=lat, lon=lon)
        return cls(z)

    def __repr__(self) -> str:
        return f"UTMZone({self.as_str()!r})"


class Bounds(PyCSAMTObject):
    r"""
    Axis-aligned geographic bounding box.

    The box is defined in degrees by ``min_lat``, ``min_lon``,
    ``max_lat``, and ``max_lon``. Inputs are validated and
    swapped as needed so that minimums are less than maximums.

    Notes
    -----
    - Coordinates are validated and coerced with ``_lat_ok``
      and ``_lon_ok`` which raise :class:`LocationError` on
      failure.
    - Methods such as :meth:`buffer_m`, :meth:`union`,
      :meth:`intersection`, and :meth:`to_utm_rect` provide
      convenient spatial operations.

    Examples
    --------
    >>> b = Bounds(-2.0, 9.5, 1.0, 12.0)
    >>> b.center()
    (-0.5, 10.75)
    >>> b.contains(0.0, 10.0)
    True

    See Also
    --------
    Location
        Point-like object that can be tested against a box.
    GeoPath
        Polyline helper; its :meth:`bbox` returns a ``Bounds``.
    """

    def __init__(
        self,
        min_lat: float,
        min_lon: float,
        max_lat: float,
        max_lon: float,
    ):
        a = _lat_ok(min_lat, "min_lat")
        c = _lat_ok(max_lat, "max_lat")
        b = _lon_ok(min_lon, "min_lon")
        d = _lon_ok(max_lon, "max_lon")
        if a > c:
            a, c = c, a
        if b > d:
            b, d = d, b
        self.min_lat = a
        self.min_lon = b
        self.max_lat = c
        self.max_lon = d

    @classmethod
    def from_points(cls, lats, lons, *, data=None) -> Bounds:
        r"""
        Create a bounding box from arrays of points.

        Parameters
        ----------
        lats : array-like or str
            Latitudes or a column name when ``data`` is provided.
        lons : array-like or str
            Longitudes or a column name when ``data`` is provided.
        data : Any, optional
            Tabular container used with column names. Resolved via
            :func:`assert_xy_in`.

        Returns
        -------
        Bounds
            The minimal box spanning the given points.

        Raises
        ------
        TypeError
            If ``data`` is required but not provided, or if sizes
            mismatch.
        LocationError
            If any element cannot be validated.

        Notes
        -----
        - Elements are validated with ``_lat_ok`` and ``_lon_ok``.
        - NaN-aware min and max are computed with NumPy.

        Examples
        --------
        >>> Bounds.from_points([0, 1], [10, 12]).to_tuple()
        (0.0, 10.0, 1.0, 12.0)
        """

        if (isinstance(lats, str) or isinstance(lons, str)) and data is None:
            raise TypeError("data required when lats/lons are column names")
        if data is not None:
            lats, lons = assert_xy_in(lats, lons, data=data)
        lats = is_iterable(lats, exclude_string=True, transform=True)
        lons = is_iterable(lons, exclude_string=True, transform=True)
        _check_consistency_size(lats, lons)
        la = np.asarray(lats, dtype=object)
        lo = np.asarray(lons, dtype=object)
        la = np.array([_lat_ok(x) for x in la], dtype=float)
        lo = np.array([_lon_ok(x) for x in lo], dtype=float)
        return cls(
            float(np.nanmin(la)),
            float(np.nanmin(lo)),
            float(np.nanmax(la)),
            float(np.nanmax(lo)),
        )

    def to_tuple(self) -> tuple[float, float, float, float]:
        r"""
        Return the box as a 4-tuple.

        Returns
        -------
        tuple of float
            ``(min_lat, min_lon, max_lat, max_lon)``.
        """

        return (
            self.min_lat,
            self.min_lon,
            self.max_lat,
            self.max_lon,
        )

    def center(self) -> tuple[float, float]:
        r"""
        Return the center point of the box.

        Returns
        -------
        tuple of float
            ``(latitude, longitude)`` at the box center.
        """

        return (
            (self.min_lat + self.max_lat) / 2.0,
            (self.min_lon + self.max_lon) / 2.0,
        )

    def contains(self, lat: float, lon: float) -> bool:
        r"""
        Test whether a point lies inside the box (inclusive).

        Parameters
        ----------
        lat : float
            Latitude in degrees.
        lon : float
            Longitude in degrees.

        Returns
        -------
        bool
            True if the point lies inside or on the edges.

        Notes
        -----
        Inputs are validated with ``_lat_ok`` and ``_lon_ok``.
        """

        la = _lat_ok(lat)
        lo = _lon_ok(lon)
        return (
            self.min_lat <= la <= self.max_lat
            and self.min_lon <= lo <= self.max_lon
        )

    def buffer_m(self, m: float) -> Bounds:
        r"""
        Return a new box expanded by a metric buffer.

        Parameters
        ----------
        m : float
            Buffer distance in meters added to all sides.

        Returns
        -------
        Bounds
            Expanded bounding box.

        Notes
        -----
        - Conversion meters-to-degrees uses a simple local scale:
          ``1 deg lat = 111320 m`` and
          ``1 deg lon = 111320 m * cos(lat)`` at the box center.
        - This is an approximation intended for small buffers.

        Examples
        --------
        >>> b = Bounds(0.0, 0.0, 1.0, 1.0)
        >>> b2 = b.buffer_m(1000.0)
        >>> b2.min_lat < b.min_lat and b2.max_lat > b.max_lat
        True
        """

        la_c, lo_c = self.center()
        d_la = m * _deg_per_m_lat()
        d_lo = m * _deg_per_m_lon(la_c)
        return Bounds(
            self.min_lat - d_la,
            self.min_lon - d_lo,
            self.max_lat + d_la,
            self.max_lon + d_lo,
        )

    def union(self, other: Bounds) -> Bounds:
        r"""
        Return a new box expanded by a metric buffer.

        Parameters
        ----------
        m : float
            Buffer distance in meters added to all sides.

        Returns
        -------
        Bounds
            Expanded bounding box.

        Notes
        -----
        - Conversion meters-to-degrees uses a simple local scale:
          ``1 deg lat = 111320 m`` and
          ``1 deg lon = 111320 m * cos(lat)`` at the box center.
        - This is an approximation intended for small buffers.

        Examples
        --------
        >>> b = Bounds(0.0, 0.0, 1.0, 1.0)
        >>> b2 = b.buffer_m(1000.0)
        >>> b2.min_lat < b.min_lat and b2.max_lat > b.max_lat
        True
        """

        return Bounds(
            min(self.min_lat, other.min_lat),
            min(self.min_lon, other.min_lon),
            max(self.max_lat, other.max_lat),
            max(self.max_lon, other.max_lon),
        )

    def intersection(self, other: Bounds) -> Bounds | None:
        r"""
        Return the minimal box containing both boxes.

        Parameters
        ----------
        other : Bounds
            Another bounding box.

        Returns
        -------
        Bounds
            The union of the two boxes.
        """

        mi_la = max(self.min_lat, other.min_lat)
        mi_lo = max(self.min_lon, other.min_lon)
        ma_la = min(self.max_lat, other.max_lat)
        ma_lo = min(self.max_lon, other.max_lon)
        if mi_la > ma_la or mi_lo > ma_lo:
            return None
        return Bounds(mi_la, mi_lo, ma_la, ma_lo)

    def to_utm_rect(
        self,
        zone: str | None = None,
    ) -> tuple[str, float, float, float, float]:
        r"""
        Map the geographic box to a UTM-aligned rectangle.

        Parameters
        ----------
        zone : str, optional
            Target UTM zone. If not given, the zone at the box
            center is used. The value is normalized to ``"NNH"``.

        Returns
        -------
        tuple
            ``(zone, min_easting, min_northing, max_easting,
            max_northing)``.

        Notes
        -----
        - Each geographic corner is converted via :func:`ll_to_utm`.
        - If a corner falls in a different zone than the requested
          one, a warning is logged and the requested zone is kept.
        - The returned rectangle is axis-aligned in the chosen UTM
          zone. It may over-approximate the true projected shape.

        References
        ----------
        .. [1] Snyder, J. P. 1987. Map Projections, A Working
               Manual. USGS Professional Paper 1395.
        """

        la_c, lo_c = self.center()
        if zone is None:
            zone, _, _ = ll_to_utm(reference_ellipsoid=23, lat=la_c, lon=lo_c)
        zone = _norm_zone(zone)
        pts = [
            (self.min_lat, self.min_lon),
            (self.min_lat, self.max_lon),
            (self.max_lat, self.min_lon),
            (self.max_lat, self.max_lon),
        ]
        es: list[float] = []
        ns: list[float] = []
        for i, (la, lo) in enumerate(pts):
            z, e, n = ll_to_utm(reference_ellipsoid=23, lat=la, lon=lo)
            if _norm_zone(z) != zone:
                logger.warning(
                    "Corner %d zone %s != %s; using %s.",
                    i,
                    z,
                    zone,
                    zone,
                )
            es.append(float(e))
            ns.append(float(n))
        return (
            zone,
            float(np.min(es)),
            float(np.min(ns)),
            float(np.max(es)),
            float(np.max(ns)),
        )

    def __repr__(self) -> str:
        a, b, c, d = self.to_tuple()
        return f"Bounds({a:.6f},{b:.6f},{c:.6f},{d:.6f})"


class GeoPath(PyCSAMTObject):
    r"""
    Lightweight polyline over geographic points.

    A ``GeoPath`` stores an ordered list of points as pairs of
    ``(latitude, longitude)`` in decimal degrees. It provides
    basic geometry such as bounding box, cumulative length on
    the sphere (haversine), UTM export, and path simplification
    via the Ramer Douglas Peucker algorithm.

    Notes
    -----
    - Points are stored as floats in degrees in the order they
      were appended.
    - :meth:`length_m` uses the haversine distance between
      consecutive vertices with a fixed Earth radius. This is
      a good approximation for short to medium paths.
    - :meth:`simplify` applies RDP in a local planar frame
      derived from the path bounding box.

    Examples
    --------
    >>> g = GeoPath([(0.0, 0.0), (0.0, 0.01), (0.01, 0.01)])
    >>> len(g)
    3
    >>> bb = g.bbox()
    ... bb.to_tuple()  # doctest: +ELLIPSIS
    (0.0, 0.0, 0.01, 0.01)
    >>> L = g.length_m()
    ... L >= 0.0
    True

    See Also
    --------
    Bounds
        Axis aligned geographic bounding box.
    Location
        Single point object with conversion helpers.
    BaseLoc.distance_to
        Great circle distance used internally.

    References
    ----------
    .. [1] Ramer, U. 1972. An iterative procedure for the
           polygonal approximation of plane curves.
    .. [2] Douglas, D. H., Peucker, T. K. 1973. Algorithms for
           the reduction of the number of points required to
           represent a digitized line or its caricature.
    """

    def __init__(self, pts=None):
        self._pts: list[tuple[float, float]] = []
        if pts is not None:
            self.extend(pts)

    def append(self, lat: _Number, lon: _Number) -> None:
        r"""
        Append one point to the end of the path.

        Parameters
        ----------
        lat : float or int or str
            Latitude in decimal degrees or parsable string.
        lon : float or int or str
            Longitude in decimal degrees or parsable string.

        Raises
        ------
        LocationError
            If the input cannot be validated.

        Notes
        -----
        Validation uses ``_lat_ok`` and ``_lon_ok`` which wrap the
        project latitude and longitude validators.
        """

        la = _lat_ok(lat)
        lo = _lon_ok(lon)
        self._pts.append((la, lo))

    def extend(self, pts) -> None:
        r"""
        Append multiple points to the path.

        Parameters
        ----------
        pts : iterable
            Iterable of ``(lat, lon)`` pairs or :class:`Location`
            objects. Locations are converted with
            ``ensure_latlon()``.

        Raises
        ------
        LocationError
            If any element fails validation.

        Examples
        --------
        >>> g = GeoPath()
        >>> g.extend([(0.0, 0.0), (0.0, 0.01)])
        >>> len(g)
        2
        """

        for p in pts:
            if isinstance(p, Location):
                la, lo = p.ensure_latlon()
                self.append(la, lo)
            else:
                la, lo = p
                self.append(la, lo)

    def __len__(self) -> int:
        r"""
        Return the number of vertices in the path.

        Returns
        -------
        int
            Count of stored points.
        """
        return len(self._pts)

    def __iter__(self):
        r"""
        Iterate over points as ``(lat, lon)`` tuples.

        Yields
        ------
        tuple of float
            Latitude and longitude in degrees, in append order.
        """

        return iter(self._pts)

    def bbox(self) -> Bounds:
        r"""
        Compute the minimal axis aligned bounding box.

        Returns
        -------
        Bounds
            Bounding box that contains all points.

        Raises
        ------
        LocationError
            If the path has no points.

        Notes
        -----
        The result is computed with NaN aware min and max over the
        latitude and longitude arrays.
        """
        if not self._pts:
            raise LocationError("Empty GeoPath")
        a = np.asarray([p[0] for p in self._pts], dtype=float)
        b = np.asarray([p[1] for p in self._pts], dtype=float)
        return Bounds(
            float(np.nanmin(a)),
            float(np.nanmin(b)),
            float(np.nanmax(a)),
            float(np.nanmax(b)),
        )

    def length_m(self) -> float:
        r"""
        Total great circle length of the path in meters.

        Returns
        -------
        float
            Sum of haversine distances between consecutive points.

        Notes
        -----
        - If the path has fewer than two points, the length is 0.
        - The haversine distance is computed with a fixed Earth
          radius. See :meth:`BaseLoc.distance_to` for the formula.

        Examples
        --------
        >>> GeoPath([(0, 0), (0, 0.01)]).length_m() >= 0.0
        True
        """

        if len(self._pts) < 2:
            return 0.0
        tot = 0.0
        pa = self._pts[0]
        for pb in self._pts[1:]:
            tot += _haversine_m(pa[0], pa[1], pb[0], pb[1])
            pa = pb
        return tot

    def to_utm(
        self,
        *,
        ensure_single_zone: bool = False,
    ):
        r"""
        Project all vertices to UTM.

        Parameters
        ----------
        ensure_single_zone : bool, optional
            If True, raise an error when points span multiple UTM
            zones. Default is False.

        Returns
        -------
        tuple
            ``(zones, eastings, northings)`` where ``zones`` is a
            ``list[str]`` and coordinates are ``numpy.ndarray``.

        Raises
        ------
        LocationError
            If ``ensure_single_zone`` is True and more than one UTM
            zone is found.

        Notes
        -----
        - Each vertex is converted via :func:`ll_to_utm`.
        - Zone strings are normalized to the ``"NNH"`` form.
        """
        zones: list[str] = []
        easts: list[float] = []
        norths: list[float] = []
        for la, lo in self._pts:
            z, e, n = ll_to_utm(reference_ellipsoid=23, lat=la, lon=lo)
            zones.append(_norm_zone(z))
            easts.append(float(e))
            norths.append(float(n))
        if ensure_single_zone and len(set(zones)) > 1:
            raise LocationError("Multiple UTM zones in path.")
        return zones, np.asarray(easts), np.asarray(norths)

    def simplify(self, eps_m: float) -> GeoPath:
        r"""
        Return a simplified path using RDP in a local planar frame.

        Parameters
        ----------
        eps_m : float
            Tolerance in meters. Larger values produce fewer points.

        Returns
        -------
        GeoPath
            New path containing a subset of the original vertices.

        Notes
        -----
        - The path is projected to a local planar metric frame by
          scaling degrees using the center latitude:
          ``1 deg lat = 111320 m`` and
          ``1 deg lon = 111320 m * cos(lat)``.
        - The RDP algorithm keeps endpoints and removes interior
          points whose perpendicular distance to the current segment
          is at most ``eps_m``.

        Examples
        --------
        >>> g = GeoPath([(0.0, 0.0), (0.0, 0.001), (0.0, 0.01)])
        >>> g2 = g.simplify(eps_m=100.0)
        >>> len(g2) <= len(g)
        True

        References
        ----------
        .. [1] Ramer, U. 1972. An iterative procedure for the
               polygonal approximation of plane curves.
        .. [2] Douglas, D. H., Peucker, T. K. 1973. Algorithms for
               the reduction of the number of points required to
               represent a digitized line or its caricature.
        """

        if len(self._pts) <= 2:
            return GeoPath(self._pts)

        la_c, _ = self.bbox().center()
        sx = _deg_per_m_lon(la_c)
        sy = _deg_per_m_lat()

        def proj(p):
            return (p[1] / sx, p[0] / sy)

        def perp(a, b, p):
            (x1, y1), (x2, y2) = proj(a), proj(b)
            x0, y0 = proj(p)
            dx = x2 - x1
            dy = y2 - y1
            if dx == 0 and dy == 0:
                return math.hypot(x0 - x1, y0 - y1)
            t = (x0 - x1) * dx + (y0 - y1) * dy
            t /= dx * dx + dy * dy
            t = max(0.0, min(1.0, t))
            xt = x1 + t * dx
            yt = y1 + t * dy
            return math.hypot(x0 - xt, y0 - yt)

        def rdp(pts, eps):
            if len(pts) <= 2:
                return pts
            a, b = pts[0], pts[-1]
            idx, dmax = 0, -1.0
            for i in range(1, len(pts) - 1):
                d = perp(a, b, pts[i])
                if d > dmax:
                    idx, dmax = i, d
            if dmax * 1.0 <= eps:
                return [a, b]
            left = rdp(pts[: idx + 1], eps)
            right = rdp(pts[idx:], eps)
            return left[:-1] + right

        out = rdp(self._pts, eps_m)
        return GeoPath(out)

    def __repr__(self) -> str:
        """Return a concise representation for debugging."""
        return f"GeoPath(n={len(self._pts)})"


def _lat_ok(x, name="lat") -> float:
    v = assert_lat_value(x)
    if v is None:
        raise LocationError(f"Invalid {name}: {x!r}")
    return float(v)


def _lon_ok(x, name="lon") -> float:
    v = assert_lon_value(x)
    if v is None:
        raise LocationError(f"Invalid {name}: {x!r}")
    return float(v)


def _as_float(x: Any, name: str) -> float:
    try:
        return float(x)  # type: ignore[arg-type]
    except Exception as exc:
        raise LocationError(f"Invalid {name}: {x!r}") from exc


def _haversine_m(
    la1: float,
    lo1: float,
    la2: float,
    lo2: float,
) -> float:
    r = 6_371_008.8
    dlat = math.radians(la2 - la1)
    dlon = math.radians(lo2 - lo1)
    a = (
        math.sin(dlat / 2.0) ** 2
        + math.cos(math.radians(la1))
        * math.cos(math.radians(la2))
        * math.sin(dlon / 2.0) ** 2
    )
    c = 2.0 * math.asin(min(1.0, math.sqrt(a)))
    return r * c


def _coerce_latlon(la: Any, lo: Any) -> tuple[float, float]:
    try:
        la_v = assert_lat_value(la)
        lo_v = assert_lon_value(lo)
    except (TypeError, ValueError) as exc:
        raise LocationError("lat/lon missing or invalid") from exc
    if la_v is None or lo_v is None:
        raise LocationError("lat/lon missing or invalid")
    return float(la_v), float(lo_v)


def _norm_zone(z: Any) -> str:
    if isinstance(z, bytes):
        z = z.decode("utf-8")
    if isinstance(z, int):
        n = abs(z)
        hem = "N" if z >= 0 else "S"
        s = f"{n}{hem}"
    else:
        s = str(z).strip().upper()
    if len(s) < 2:
        raise LocationError(f"Invalid UTM zone: {z!r}")
    num = s[:-1]
    hem = s[-1]
    try:
        n = int(num)
    except Exception as exc:
        raise LocationError(f"Invalid UTM zone: {z!r}") from exc
    if not (1 <= n <= 60) or hem not in "NS":
        raise LocationError(f"Invalid UTM zone: {z!r}")
    return f"{n}{hem}"


def _deg_per_m_lat() -> float:
    # ~111.32 km per degree latitude
    return 1.0 / 111_320.0


def _deg_per_m_lon(lat: float) -> float:
    # scale by cos(lat)
    c = math.cos(math.radians(lat))
    return 1.0 / (111_320.0 * max(c, 1e-12))


def _to_float(name: str, value: _Number) -> float:
    try:
        v = float(value)  # type: ignore[arg-type]
    except Exception as exc:
        raise LocationError(f"Invalid {name}: {value!r}") from exc
    return v


def _in_range(name: str, v: float, lo: float, hi: float) -> float:
    if np.isnan(v) or v < lo or v > hi:
        raise LocationError(f"{name} out of range [{lo}, {hi}]: {v}")
    return v
