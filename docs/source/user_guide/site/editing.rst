.. _site_editing:

Site Editing
============

.. currentmodule:: pycsamt.site.edit

The :mod:`pycsamt.site.edit` module contains practical editing helpers for
EDI-like station objects and site collections. These helpers are designed for
survey preparation: rotate tensors, reduce frequency ranges, normalize station
names, assign coordinates, fill missing arrays, and recompute derived
resistivity/phase quantities after edits.

The editing functions are deliberately tolerant. When an optional section is
missing or has an incompatible shape, most functions skip that part and keep
processing the rest of the object. This is helpful for real field data, where
some stations may lack tipper, errors, or derived arrays.

For a complete survey-level cleanup that reads EDI files, applies edits,
recomputes derived quantities, preserves line folders, and writes new
pyCSAMT-authored EDI files, use :doc:`recompute`. The functions on this page
are the lower-level building blocks used by that workflow.

Editing Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Function
     - Scope
     - Main purpose
   * - :func:`rotate`
     - One site
     - Rotate impedance tensors and tipper values by an azimuthal angle.
   * - :func:`select_freq`
     - One site
     - Subset frequency-indexed arrays by range, indices, or boolean mask.
   * - :func:`rename`
     - One site
     - Set an explicit station name or apply a naming policy.
   * - :func:`set_coords`
     - One site
     - Update latitude, longitude, and elevation in the EDI header.
   * - :func:`fill_missing`
     - One site
     - Fill or allocate missing Z and tipper arrays.
   * - :func:`recompute_res_phase`
     - One site
     - Recompute apparent resistivity and phase from the impedance tensor.
   * - :func:`rotate_all`
     - Collection
     - Rotate every site in a collection.
   * - :func:`select_freq_all`
     - Collection
     - Apply frequency subsetting to every site.
   * - :func:`rename_all`
     - Collection
     - Rename stations across a collection.
   * - :func:`set_coords_all`
     - Collection
     - Assign coordinates from a callable, mapping, or table holder.
   * - :func:`set_coords_from_table`
     - Collection
     - Assign coordinates from CSV, DataFrame, structured array, or row list.
   * - :func:`set_coords_from_en`
     - One site
     - Project easting/northing coordinates and store lon/lat on a site.

Copy Versus In-Place Editing
----------------------------

Most editing helpers accept ``inplace``. The default is ``False``:

.. code-block:: python
   :linenos:

   from pycsamt.seg.edi import EDIFile
   from pycsamt.site.edit import rename

   edi = EDIFile("data/edi/S01.edi")

   edited = rename(edi, name="LINE01_S01", inplace=False)

   print(edited is edi)  # False

Use ``inplace=True`` when you intentionally want to mutate the provided object:

.. code-block:: python
   :linenos:

   rename(edi, name="LINE01_S01", inplace=True)

For collection helpers, ``inplace=False`` returns a new
:class:`pycsamt.site.base.Sites` wrapper. The original EDI objects are left
untouched as far as the helper can manage. Use ``inplace=True`` only when the
calling workflow owns the input objects.

Tensor Rotation
---------------

:func:`rotate` rotates the impedance tensor and, when present, the tipper.
For the impedance tensor, pyCSAMT applies:

.. math::

   Z' = R Z R^{-1},

where :math:`R` is the horizontal rotation matrix built from ``angle_deg``.

.. code-block:: python
   :linenos:

   from pycsamt.seg.edi import EDIFile
   from pycsamt.site.edit import rotate

   edi = EDIFile("data/edi/S01.edi")
   rotated = rotate(edi, angle_deg=30.0)

The helper checks common tensor attribute names such as ``z``, ``impedance``,
and ``_z``. Tipper arrays may live under ``T``, ``TIP``, ``Tip``, or sometimes
as a tipper-like attribute under ``Z``.

Error arrays are rotated with a magnitude-only approximation using absolute
values of the rotation matrices. Treat this as a practical uncertainty
propagation, not a full covariance rotation.

Rotate a whole collection with :func:`rotate_all`:

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.edit import rotate_all

   collection = EDICollection.from_sources("data/edi")
   rotated_sites = rotate_all(collection, angle_deg=30.0)

Frequency Subsetting
--------------------

:func:`select_freq` subsets a station along its frequency axis and slices
aligned arrays together. It handles Z, Z errors, apparent resistivity, phase,
tipper, tipper errors, and several common aliases.

Keep a frequency range:

.. code-block:: python
   :linenos:

   from pycsamt.site.edit import select_freq

   band = select_freq(
       edi,
       fmin=1.0,
       fmax=1000.0,
       inplace=False,
   )

Keep explicit row indices:

.. code-block:: python
   :linenos:

   edges = select_freq(edi, keep=[0, -1])

Use a boolean mask:

.. code-block:: python
   :linenos:

   import numpy as np

   freq = np.asarray(edi.Z.freq)
   mask = freq >= 10.0
   high = select_freq(edi, keep=mask)

When ``keep`` is provided, ``fmin`` and ``fmax`` are ignored. Use
:func:`select_freq_all` for collections:

.. code-block:: python
   :linenos:

   from pycsamt.site.edit import select_freq_all

   trimmed = select_freq_all(
       collection,
       fmin=1.0,
       fmax=1000.0,
   )

Station Renaming
----------------

:func:`rename` updates common station identifiers in the EDI header and mirrors
the result to ``edi.name`` when possible. This makes later lookup through
:class:`pycsamt.site.base.Site` and :class:`pycsamt.site.base.Sites` more
consistent.

Set an explicit name:

.. code-block:: python
   :linenos:

   from pycsamt.site.edit import rename

   renamed = rename(edi, name="LINE01_S01")

Apply a policy to the current station name:

.. code-block:: python
   :linenos:

   renamed = rename(
       edi,
       policy=lambda old: f"LINE01_{old}",
   )

If both ``name`` and ``policy`` are provided, the explicit ``name`` wins.
Renaming does not change the on-disk filename. Use export templates later if
you want output filenames to follow the new station names.

For collections, use :func:`rename_all`:

.. code-block:: python
   :linenos:

   from pathlib import Path
   from pycsamt.site.edit import rename_all

   renamed_sites = rename_all(
       collection,
       name_fn=lambda edi: f"L01_{Path(getattr(edi, 'path', '')).stem}",
   )

Prefer ``name_fn`` when the original EDI files have duplicate station names.
It receives the EDI object and can use file paths, metadata, or external state
to generate unique identifiers.

Coordinate Editing
------------------

:func:`set_coords` updates latitude, longitude, and elevation on one site.
Only the values provided are changed.

.. code-block:: python
   :linenos:

   from pycsamt.site.edit import set_coords

   moved = set_coords(
       edi,
       lat=35.125,
       lon=12.750,
       elev=1234.0,
   )

To apply coordinates to many sites, use :func:`set_coords_all`.

From a mapping keyed by station name:

.. code-block:: python
   :linenos:

   from pycsamt.site.edit import set_coords_all

   coords = {
       "S01": (35.125, 12.750, 1234.0),
       "S02": (35.200, 12.900, 1180.0),
   }

   updated = set_coords_all(collection, coords)

From a callable:

.. code-block:: python
   :linenos:

   def lookup(edi):
       name = getattr(edi, "name", "")
       if name == "S01":
           return 35.125, 12.750, 1234.0
       return None

   updated = set_coords_all(collection, lookup)

From an object exposing a ``.frame`` attribute:

.. code-block:: python
   :linenos:

   class CoordinateTable:
       def __init__(self, frame):
           self.frame = frame

   updated = set_coords_all(collection, CoordinateTable(frame))

Coordinate Tables
-----------------

:func:`set_coords_from_table` is the high-level table loader. It accepts:

* path to a CSV file;
* path to a whitespace-delimited text file;
* :class:`pandas.DataFrame`;
* NumPy structured array;
* list of dictionaries;
* list of row-like tuples that pandas can convert to a frame.

With standard geographic columns:

.. code-block:: python
   :linenos:

   import pandas as pd
   from pycsamt.site.edit import set_coords_from_table

   table = pd.DataFrame(
       {
           "station": ["S01", "S02"],
           "lat": [35.125, 35.200],
           "lon": [12.750, 12.900],
           "elev": [1234.0, 1180.0],
       }
   )

   updated = set_coords_from_table(collection, table)

The resolver understands common aliases:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Canonical field
     - Accepted names
   * - ``station``
     - ``station``, ``name``, ``site``, ``id``.
   * - ``lat``
     - ``lat``, ``latitude``.
   * - ``lon``
     - ``lon``, ``long``, ``longitude``.
   * - ``elev``
     - ``elev``, ``elevation``, ``z``.
   * - ``easting``
     - ``easting``, ``x``.
   * - ``northing``
     - ``northing``, ``y``.

Use an explicit column map when a field table uses non-standard names:

.. code-block:: python
   :linenos:

   updated = set_coords_from_table(
       collection,
       table,
       columns={
           "station": "name",
           "lat": "latitude",
           "lon": "long",
           "elev": "elevation",
       },
   )

When both geographic coordinates and easting/northing are present, latitude and
longitude are preferred.

Easting/Northing Conversion
---------------------------

If a table provides projected coordinates, pass ``crs_from``:

.. code-block:: python
   :linenos:

   table = pd.DataFrame(
       {
           "station": ["S10", "S11"],
           "easting": [400000.0, 401250.0],
           "northing": [5750000.0, 5750400.0],
           "elev": [250.0, 252.0],
       }
   )

   updated = set_coords_from_table(
       collection,
       table,
       crs_from="EPSG:32631",
   )

Projection uses :mod:`pyproj`. If ``pyproj`` is not installed and projected
coordinates must be converted, the helper raises ``ImportError``.

For one site, use :func:`set_coords_from_en`:

.. code-block:: python
   :linenos:

   from pycsamt.site.edit import set_coords_from_en

   projected = set_coords_from_en(
       edi,
       easting=400000.0,
       northing=5750000.0,
       crs_from="EPSG:32631",
       elev=250.0,
   )

Missing Data Filling
--------------------

:func:`fill_missing` fills or allocates missing Z and tipper arrays. It is
often used before diagnostics or before recomputing derived quantities.

.. code-block:: python
   :linenos:

   from pycsamt.site.edit import fill_missing

   filled = fill_missing(
       edi,
       how="zero",
       components=("Z",),
   )

Available fill policies are:

``how="zero"``
   Replace non-finite values with numeric zeros.

``how="nan"``
   Replace non-finite values with ``NaN``.

The ``components`` argument is case-insensitive and accepts ``"Z"`` and
``"Tip"``. Z arrays are expected to have shape ``(n_freq, 2, 2)``. Tipper
arrays are expected to have shape ``(n_freq, 2)``. The number of rows is
inferred from the frequency vector.

Use this function carefully. Zero-filled tensors are convenient for tests and
some robust workflows, but zeros are not measured values. Keep a note in the
processing log when missing data have been filled.

Recomputing Resistivity And Phase
---------------------------------

:func:`recompute_res_phase` calls the available Z-section
``compute_resistivity_phase()`` method. It is best used after changing
frequency rows, tensor values, or masks.

.. code-block:: python
   :linenos:

   from pycsamt.site.edit import recompute_res_phase

   edited = select_freq(edi, fmin=1.0, fmax=1000.0)
   edited = recompute_res_phase(edited)

The function is best-effort. If the Z section is missing, incompatible, or does
not expose a recomputation method, the object is returned unchanged.

Practical Preparation Workflow
------------------------------

The following pattern is common before diagnostics or inversion preparation:

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.edit import (
       fill_missing,
       recompute_res_phase,
       rename_all,
       rotate_all,
       select_freq_all,
       set_coords_from_table,
   )

   collection = EDICollection.from_sources("data/raw_edi")

   sites = rename_all(
       collection,
       name_fn=lambda edi: f"L01_{getattr(edi, 'name', 'site')}",
   )
   sites = set_coords_from_table(sites, "data/station_coordinates.csv")
   sites = rotate_all(sites, angle_deg=30.0)
   sites = select_freq_all(sites, fmin=1.0, fmax=1000.0)

   prepared = []
   for edi in sites.as_list():
       edi = fill_missing(edi, how="nan", components=("Z",))
       edi = recompute_res_phase(edi)
       prepared.append(edi)

Each stage should be recorded in a processing log, especially rotation angles,
frequency bands, coordinate sources, and missing-data fill policies.

Common Mistakes
---------------

Using zero fill as if it were measured data
   ``fill_missing(..., how="zero")`` is useful for deterministic tests and
   some defensive algorithms. It should not silently replace failed field
   measurements in scientific interpretation.

Renaming without checking duplicates
   A policy such as ``lambda name: "X_" + name`` preserves duplicates if the
   original station names were duplicated. Use ``name_fn`` with file stems or a
   station table when uniqueness matters.

Mixing frequency masks from different stations
   ``keep=[...]`` and boolean masks are applied row-wise. A mask built from one
   station may not represent the same frequencies on another station if their
   frequency axes differ.

Forgetting to recompute derived arrays
   After changing Z or frequency rows, apparent resistivity and phase may be
   stale. Run :func:`recompute_res_phase` when downstream steps use derived
   arrays.

Ignoring CRS metadata
   Easting/northing values are meaningless without their source CRS. Always
   provide the correct ``crs_from`` when using projected coordinates.

Next Pages
----------

Continue with:

* :doc:`containers` for the object model behind edited sites;
* :doc:`selection` for filtering stations before editing;
* :doc:`recompute` for end-to-end EDI recomputation and pyCSAMT rewriting;
* :doc:`location_profile` for coordinate, distance, and chainage tools;
* :doc:`computed_diagnostics` for checking edited sites before export or
  inversion.
