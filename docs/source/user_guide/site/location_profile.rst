.. _site_location_profile:

Location And Profiles
=====================

.. currentmodule:: pycsamt.site.location

The site location tools handle station coordinates, topography assignment,
coordinate projection, distance and bearing calculations, and survey-line
chainage. The profile tools in :mod:`pycsamt.site.profile` build a 1-D line
model from station locations so that a survey can be sorted, sliced, checked
for gaps, and passed cleanly into processing or inversion preparation.

Use this page when you need to:

* parse latitude, longitude, or elevation values from field tables;
* normalize EDI ``HEAD`` coordinate fields;
* apply corrected GPS or topography tables to EDI sites;
* project lon/lat coordinates into another CRS;
* compute station spacing, bearing, or chainage along a survey line;
* infer the dominant line orientation from station coordinates;
* build a :class:`pycsamt.site.profile.Profile` for line-ordered workflows.

Location Tool Map
-----------------

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Object or function
     - Scope
     - Main purpose
   * - :class:`Coord`
     - One coordinate.
     - Store latitude, longitude, and elevation as a small dataclass.
   * - :func:`parse_lat`, :func:`parse_lon`, :func:`parse_elev`
     - Field values.
     - Parse decimal or DMS-like coordinate text into numeric values.
   * - :func:`ensure_head_coords`
     - One EDI-like object.
     - Ensure the EDI ``HEAD`` section has numeric ``lat``, ``lon``/``long``,
       and ``elev`` fields.
   * - :func:`apply_topography`
     - One or many sites.
     - Update site coordinates from a table matched by station identifier.
   * - :func:`project`
     - Coordinate arrays.
     - Transform points between coordinate reference systems.
   * - :func:`distance`
     - Two points.
     - Compute geodetic, flat, or projected distance in meters.
   * - :func:`bearing`
     - Two points.
     - Compute azimuth from one point to another.
   * - :func:`chainage_along`
     - One profile axis.
     - Project station positions onto a survey-line axis.

Profile Tool Map
----------------

.. currentmodule:: pycsamt.site.profile

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Object or function
     - Main purpose
   * - :class:`Profile`
     - Store profile origin, azimuth, per-station chainage, spacing
       statistics, and detected large gaps.
   * - :meth:`Profile.from_sites`
     - Build a profile from EDI-like sites or ``Site`` wrappers.
   * - :meth:`Profile.sort_sites`
     - Return sites ordered by increasing chainage.
   * - :meth:`Profile.slice`
     - Return stations whose chainage lies inside a distance window.
   * - :meth:`Profile.resample`
     - Build a regular chainage grid between minimum and maximum station
       chainage.
   * - :meth:`Profile.summary`
     - Return spacing and gap summary statistics.
   * - :func:`infer_line_orientation`
     - Estimate survey-line azimuth using PCA on local coordinate offsets.

Coordinate Containers
---------------------

Use :class:`pycsamt.site.location.Coord` when a small explicit coordinate
object is clearer than passing loose tuples.

.. code-block:: python
   :linenos:

   from pycsamt.site.location import Coord

   station = Coord(lat=10.25, lon=20.75, elev=640.0)

   print(station.lat)
   print(station.lon)
   print(station.elev)

``Coord`` does not validate values at construction time. Validation happens in
helpers such as :func:`ensure_head_coords`, :func:`distance`,
:func:`bearing`, and :func:`chainage_along`.

Parsing Coordinates
-------------------

Field coordinate tables are not always consistent. The parser accepts numeric
values and common degree/minute/second style strings.

.. code-block:: python
   :linenos:

   from pycsamt.site.location import parse_elev, parse_lat, parse_lon

   lat1 = parse_lat("45N")
   lat2 = parse_lat("45 30 0 S")
   lon1 = parse_lon("123W")
   lon2 = parse_lon("123 15 30 E")
   elev = parse_elev("1200")

   print(lat1, lat2)
   print(lon1, lon2)
   print(elev)

Parsing conventions:

* north and east are positive;
* south and west are negative;
* signed decimal values are accepted;
* unparseable values return ``NaN`` instead of raising;
* elevation is interpreted in meters.

Accepted examples include ``"10"``, ``"-3.2"``, ``"12 30 00 N"``,
``"12:30:00W"``, ``"3.5W"``, and numeric floats.

Normalizing EDI Head Coordinates
--------------------------------

:func:`ensure_head_coords` creates or updates coordinate fields on an EDI-like
object. It guarantees that the ``HEAD`` section carries numeric ``lat``,
``lon``, ``long``, and ``elev`` values.

.. code-block:: python
   :linenos:

   from pycsamt.seg.edi import EDIFile
   from pycsamt.site.location import ensure_head_coords

   edi = EDIFile("data/edi/S01.edi")

   head = ensure_head_coords(
       edi,
       lat="35 07 30 N",
       lon="12 45 00 E",
       elev="1234",
   )

   print(head.lat, head.lon, head.long, head.elev)

If a coordinate is missing or invalid, the function uses ``empty`` as the
fallback sentinel. By default ``empty`` is ``0.0``.

.. code-block:: python
   :linenos:

   head = ensure_head_coords(edi, empty=-9999.0)

Use this helper early when downstream code expects numeric coordinates.

Applying Topography Tables
--------------------------

:func:`apply_topography` updates one site, a list of sites, or a container with
``._items`` from a table matched by station identifier.

.. code-block:: python
   :linenos:

   import pandas as pd

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.location import apply_topography

   collection = EDICollection.from_sources("data/edi")

   topo = pd.DataFrame(
       {
           "station": ["S01", "S02", "S03"],
           "latitude": [35.125, 35.200, 35.275],
           "longitude": [12.750, 12.900, 13.050],
           "elevation": [1234.0, 1180.0, 1215.0],
       }
   )

   updated = apply_topography(collection, topo, inplace=False)

Station identifiers are matched case-insensitively against common columns:
``station``, ``site``, ``dataid``, ``id``, and ``name``. Coordinate columns are
searched among:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Coordinate
     - Accepted columns
   * - Latitude
     - ``latitude``, ``lat``.
   * - Longitude
     - ``longitude``, ``lon``, ``long``.
   * - Elevation
     - ``elevation``, ``elev``, ``alt``.

Use ``inplace=False`` when you want a copied structure and ``inplace=True``
when the input objects should be updated directly.

Coordinate Projection
---------------------

:func:`project` transforms coordinate pairs between coordinate reference
systems. When :mod:`pyproj` is available, pyCSAMT uses it with
``always_xy=True``. GDAL is used as a fallback when available.

.. code-block:: python
   :linenos:

   from pycsamt.site.location import project

   x, y = project(
       [(12.5, 35.25), (12.6, 35.30)],
       crs_from="EPSG:4326",
       crs_to="EPSG:32633",
   )

   print(x)
   print(y)

Important convention: projected coordinate pairs are interpreted as ``(x, y)``.
For geographic CRS values, that means ``(lon, lat)`` when using
``EPSG:4326`` with this function.

If neither pyproj nor GDAL is available, :func:`project` raises
``RuntimeError``.

Distance And Bearing
--------------------

:func:`distance` computes the distance between two coordinates in meters.
Three modes are available:

``mode="geodetic"``
   Haversine approximation on a spherical Earth. This is the default.

``mode="flat"``
   Local small-extent planar approximation using latitude-scaled longitude
   distances.

``mode="utm"``
   Project both points to UTM, or to ``crs_to`` when supplied, and compute
   Euclidean distance.

.. code-block:: python
   :linenos:

   from pycsamt.site.location import Coord, bearing, distance

   a = Coord(0.0, 0.0, 0.0)
   b = Coord(0.0, 1.0, 0.0)

   d = distance(a, b, mode="geodetic")
   az = bearing(a, b)

   print(d)   # about 111 km at the equator
   print(az)  # about 90 degrees, east

:func:`bearing` uses the same mode names. It returns azimuth in degrees with
0 degrees pointing north and 90 degrees pointing east.

The geodetic distance is based on:

.. math::

   d =
   2R \arcsin
   \left(
      \sqrt{
         \sin^2(\Delta \phi / 2)
         + \cos \phi_1 \cos \phi_2
           \sin^2(\Delta \lambda / 2)
      }
   \right).

The geodetic bearing uses the spherical forward azimuth:

.. math::

   \theta =
   \operatorname{atan2}
   \left(
      \sin\Delta\lambda \cos\phi_2,\,
      \cos\phi_1 \sin\phi_2
      - \sin\phi_1 \cos\phi_2 \cos\Delta\lambda
   \right).

Chainage Along A Line
---------------------

:func:`chainage_along` projects station coordinates onto a survey-line axis.
The axis is defined by an origin and azimuth. Azimuth follows the geophysical
line convention used in pyCSAMT: 0 degrees is north and 90 degrees is east.

.. code-block:: python
   :linenos:

   from pycsamt.site.location import chainage_along

   origin = (0.0, 0.0)
   station = (0.0, 1.0 / 111.0)  # about 1 km east near the equator

   s = chainage_along(origin, azimuth=90.0, pts=station)
   print(s)

For local offsets :math:`dx` east and :math:`dy` north, chainage is:

.. math::

   s = dx \sin A + dy \cos A,

where :math:`A` is the profile azimuth in radians.

The function returns a scalar for one point and a NumPy array for a sequence of
points.

.. code-block:: python
   :linenos:

   points = [
       (0.0, 0.0),
       (0.0, 1.0 / 111.0),
       (0.0, 2.0 / 111.0),
   ]

   chainages = chainage_along(origin, 90.0, points)
   print(chainages)

Inferring Line Orientation
--------------------------

.. currentmodule:: pycsamt.site.profile

:func:`infer_line_orientation` estimates the dominant survey-line axis from a
set of site coordinates. It builds local Cartesian offsets and applies PCA.
The principal component with the largest variance defines the line direction.

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.profile import infer_line_orientation

   collection = EDICollection.from_sources("data/edi")

   azimuth = infer_line_orientation(collection)
   print(azimuth)

The returned azimuth is in degrees, with 0 degrees north and 90 degrees east.
Because a line axis has no intrinsic direction, the result should be
interpreted modulo 180 degrees. For example, 45 degrees and 225 degrees
describe the same line axis.

Building Profiles
-----------------

:class:`Profile` stores a survey-line representation:

* ``origin`` as a :class:`pycsamt.site.location.Coord`;
* ``azimuth`` in degrees;
* ``chainages`` as ``{station_name: chainage_m}``;
* ``spacing_stats`` for line spacing;
* ``gaps`` for large station-spacing intervals.

Build one from sites:

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.profile import Profile

   collection = EDICollection.from_sources("data/edi")

   profile = Profile.from_sites(collection)

   print(profile.azimuth)
   print(profile.chainages)
   print(profile.spacing_stats)
   print(profile.gaps)

If no ``origin`` is supplied, the first finite site coordinate is used. If no
``azimuth`` is supplied, :func:`infer_line_orientation` is used.

You can also force both:

.. code-block:: python
   :linenos:

   from pycsamt.site.location import Coord

   profile = Profile.from_sites(
       collection,
       origin=Coord(35.0, 12.0, 0.0),
       azimuth=90.0,
   )

Sorting And Slicing Profiles
----------------------------

Use :meth:`Profile.sort_sites` to return sites ordered by increasing chainage.
Sites without finite coordinates are dropped.

.. code-block:: python
   :linenos:

   ordered = profile.sort_sites(collection)
   print([site.name for site in ordered if hasattr(site, "name")])

Use :meth:`Profile.slice` to select a chainage window:

.. code-block:: python
   :linenos:

   window = profile.slice(500.0, 2500.0)
   print(window)

The returned mapping is ordered by chainage.

Resampling And Summary
----------------------

:meth:`Profile.resample` builds a regular chainage grid between the minimum and
maximum station chainage.

.. code-block:: python
   :linenos:

   grid = profile.resample(step=250.0)
   print(grid)

:meth:`Profile.summary` returns spacing and gap information:

.. code-block:: python
   :linenos:

   summary = profile.summary()

   print(summary["n_sites"])
   print(summary.get("spacing_mean"))
   print(summary.get("spacing_med"))
   print(summary["n_gaps"])

Spacing statistics include:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``spacing_mean``
     - Mean station spacing in meters.
   * - ``spacing_med``
     - Median station spacing in meters.
   * - ``spacing_min``
     - Minimum station spacing in meters.
   * - ``spacing_max``
     - Maximum station spacing in meters.
   * - ``s_min``, ``s_max``
     - Minimum and maximum finite chainage.
   * - ``n_sites``
     - Number of stations with finite chainage.
   * - ``n_gaps``
     - Number of large spacing gaps.

Gap Detection
-------------

After chainages are sorted, :class:`Profile` computes station spacings with
``numpy.diff``. A large gap is recorded when spacing exceeds
``1.5 * median_spacing``. Gaps are stored as ``(s_left, s_right)`` chainage
intervals.

.. code-block:: python
   :linenos:

   if profile.gaps:
       for left, right in profile.gaps:
           print(f"Gap from {left:.1f} m to {right:.1f} m")

This is a simple QC rule, not a geological interpretation. Use it to find
survey-line acquisition holes, missing stations, or irregular spacing before
building inversions or pseudo-sections.

End-To-End Example
------------------

The following example applies a topography table, builds a profile, sorts the
survey, and checks spacing.

.. code-block:: python
   :linenos:

   import pandas as pd

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.location import apply_topography
   from pycsamt.site.profile import Profile

   collection = EDICollection.from_sources("data/edi")

   topo = pd.read_csv("data/station_coordinates.csv")
   collection = apply_topography(collection, topo, inplace=False)

   profile = Profile.from_sites(collection)
   ordered = profile.sort_sites(collection)
   summary = profile.summary()

   print(profile.azimuth)
   print(summary)
   print([site.name for site in ordered if hasattr(site, "name")])

Common Mistakes
---------------

Swapping latitude and longitude in projection
   :func:`project` uses ``always_xy=True`` when pyproj is available. For
   ``EPSG:4326``, pass ``(lon, lat)`` pairs, not ``(lat, lon)`` pairs.

Treating local flat chainage as high-precision geodesy
   :func:`chainage_along` uses a local flat approximation. It is suitable for
   survey-line ordering and spacing checks, not cadastral-grade positioning.

Forgetting the 180-degree ambiguity of line orientation
   A line azimuth of 45 degrees and 225 degrees describes the same physical
   axis. Choose the direction that matches acquisition convention when
   reporting chainage.

Using uncorrected EDI header coordinates
   Raw EDI coordinates can be missing, rounded, or copied from acquisition
   templates. Apply the authoritative station table before building profiles.

Ignoring missing coordinates
   Profile methods drop stations without finite coordinates. Check
   ``profile.summary()`` and compare the number of profiled stations with the
   number of loaded sites.

Next Pages
----------

Continue with:

* :doc:`containers` for the site and collection objects that carry locations;
* :doc:`editing` for assigning coordinates from tables and projected values;
* :doc:`selection` for selecting stations by chainage, bounding box, or
  frequency coverage;
* :doc:`export_reporting` for writing cleaned sites and reporting profile-ready
  station sets.

