3-D Quick-Look Maps
===================

The volume tools build 3-D quick-look visualizations from profile
pseudosection data.  They estimate pseudo-depth from apparent
resistivity and period, then render the values as fence, block,
depth-slice, or isosurface views.

.. warning::

   These maps are not inversion models.  Use them to inspect trends,
   compare lines, and communicate quick-look targets.  Geological
   interpretation should be checked against inversion and QC products.

What The 3-D Map Represents
---------------------------

The 3-D map module does not read an inversion mesh.  It starts from the
same impedance-derived pseudosection table used by the profile tools,
then places each period sample at an approximate skin-depth scale:

.. math::

   z \approx 503 \sqrt{\rho_a T}

where ``rho_a`` is apparent resistivity and ``T`` is period.  This is a
quick-look depth proxy.  It is useful for comparing survey lines and
screening targets, but it should not be treated as a recovered earth
model.

The displayed color can be apparent resistivity or phase:

``quantity="resistivity"`` or ``quantity="rho"``
   Color by apparent resistivity.

``quantity="phase"``
   Color by phase, while apparent resistivity is still used to derive
   pseudo-depth and to apply ``rho_range`` filters.

Data Preparation
----------------

For a single line, pass the EDI folder directly.  For multi-line 3-D
views, load all lines first so line names and station ordering are
stable across every mode.

.. code-block:: python

   from pycsamt.map import load_lines

   data = load_lines(
       "data/AMT/WILLY_DATA",
       detect="folder",
       recursive=True,
   )

   print(data.lines)
   print(data.station_ids[:5])

The volume builder groups stations by ``StationRecord.line``.  If no
line metadata is available, every station is placed into a single line
named ``"line"``.

Function API
------------

Use :func:`pycsamt.map.plot_volume_map` or the equivalent
:func:`pycsamt.map.plot_3d_map` for one-shot figures.

.. code-block:: python

   from pycsamt.map import VolumeMapOptions, plot_volume_map

   fig = plot_volume_map(
       data,
       options=VolumeMapOptions(
           mode="fence",
           quantity="resistivity",
           component="xy",
           depth_range=(0.0, 2000.0),
           period_range=(0.001, 10.0),
           show_stations=True,
       ),
   )

The returned object is a Plotly figure.  Use ``fig.show()`` in a
notebook or export it with :doc:`export`.

Builder API
-----------

Use :class:`pycsamt.map.VolumeMap` when you want to reuse normalized
data and switch modes or quantities without reloading files.
``VolumeMap`` is an alias of :class:`pycsamt.map.Map3D`.

.. code-block:: python

   from pycsamt.map import VolumeMap

   fig = (
       VolumeMap("data/AMT/WILLY_DATA/L18PLT")
       .with_mode("surface")
       .with_quantity("phase")
       .with_component("xy")
       .figure()
   )

The builder methods are immutable: each call returns a new builder that
shares the same ``MapData`` but carries different options.

.. code-block:: python

   base = VolumeMap(data).with_options(
       depth_range=(0.0, 2500.0),
       component="xy",
   )

   fence = base.with_mode("fence").figure()
   slices = base.with_mode("depth").with_options(n_slices=6).figure()

Modes
-----

``fence``
   One pseudo-depth surface per survey line.  This is the best first
   view for multi-line surveys because it keeps each profile readable.

``block``
   Sparse volume rendering from all finite pseudo-depth samples.  It is
   useful for a compact 3-D impression, but can hide line structure on
   sparse surveys.

``depth``
   Horizontal pseudo-depth slices.  Values are interpolated at the
   slice depths generated from ``depth_range`` or from the available
   pseudo-depth span.

``surface``
   Isosurfaces across the pseudo-depth point cloud.  Use
   ``iso_range`` and ``surface_count`` to control which value shells
   are visible.

Mode Examples
-------------

Fence view with one draped surface per line:

.. code-block:: python

   fig = plot_volume_map(
       data,
       options=VolumeMapOptions(
           mode="fence",
           quantity="resistivity",
           component="xy",
           show_labels=True,
           show_contours=True,
       ),
   )

Block volume view:

.. code-block:: python

   fig = plot_volume_map(
       data,
       options=VolumeMapOptions(
           mode="block",
           opacity=0.35,
           surface_count=18,
       ),
   )

Depth slices:

.. code-block:: python

   fig = plot_volume_map(
       data,
       options=VolumeMapOptions(
           mode="depth",
           depth_range=(0.0, 3000.0),
           n_slices=7,
           show_contours=True,
       ),
   )

Isosurface view:

.. code-block:: python

   fig = plot_volume_map(
       data,
       options=VolumeMapOptions(
           mode="surface",
           iso_range=(1.0, 3.0),
           surface_count=5,
           opacity=0.55,
       ),
   )

Filtering
---------

``depth_range`` clips the pseudo-depth axis.  ``period_range`` filters
the periods before grid construction.  ``rho_range`` masks samples by
apparent resistivity, even when the displayed quantity is phase.

``iso_range`` controls the value range used by isosurface rendering.

Use ``value_range`` to keep colorbars comparable across multiple
figures:

.. code-block:: python

   options = VolumeMapOptions(
       mode="fence",
       quantity="resistivity",
       log_color=True,
       value_range=(10.0, 10000.0),
       rho_range=(10.0, 10000.0),
       period_range=(0.001, 5.0),
       depth_range=(0.0, 2500.0),
   )

``rho_range`` filters in physical apparent-resistivity units.  When
``log_color=True``, ``value_range`` is converted to log10 color space
for resistivity colorbars.

Components
----------

The ``component`` option selects the impedance component used to build
the pseudosection table:

``"xy"``, ``"yx"``, ``"xx"``, ``"yy"``
   Individual tensor components.

``"avg"``
   Average of ``xy`` and ``yx``.

``"det"``
   Determinant-style derived value for resistivity, or average phase.

Use the same component for volume maps that you use in profile
pseudosections when you want the 2-D and 3-D views to compare directly.

Line Spacing And Azimuth
------------------------

``line_spacing`` controls the offset between profile lines in the 3-D
scene.  ``azimuth`` rotates those line offsets.

.. code-block:: python

   fig = plot_volume_map(
       data,
       options=VolumeMapOptions(
           mode="fence",
           line_spacing=1.5,
           azimuth=30.0,
       ),
   )

With ``azimuth=0``, line offsets appear along the scene ``y`` axis.
With ``azimuth=90``, offsets are shifted into the ``x`` direction.

Topography And Terrain
----------------------

By default, depth is plotted downward from a flat surface.  Enable
``topography`` to use station elevations from the loaded EDI metadata:

.. code-block:: python

   fig = plot_volume_map(
       data,
       options=VolumeMapOptions(
           mode="fence",
           topography=True,
           show_terrain=True,
       ),
   )

When ``topography=True``, the vertical axis is labelled
``Elevation - depth (m)`` and each pseudo-depth surface is shifted by
the station elevations.  ``show_terrain=True`` adds a terrain trace at
the top of each line.  If elevations are missing, zeros are used for
those stations.

Station Markers
---------------

Set ``show_stations=True`` to add station markers at the survey
surface.

.. code-block:: python

   fig = plot_volume_map(
       data,
       options=VolumeMapOptions(
           mode="fence",
           show_stations=True,
           station_symbol="diamond",
           station_size=5,
           station_color="#111827",
       ),
   )

Markers use the same line offsets and optional topography shift as the
volume surfaces.

Color And Theme Controls
------------------------

Volume maps support the shared map themes and Plotly color scales:

.. code-block:: python

   fig = plot_volume_map(
       data,
       options=VolumeMapOptions(
           mode="depth",
           theme="dark",
           cmap="Turbo",
           opacity=0.75,
           title="Depth slices: XY resistivity",
       ),
   )

For resistivity, ``log_color=True`` is the default.  For phase, values
are shown linearly and the colorbar title becomes ``Phase (deg)``.

Exporting 3-D Views
-------------------

HTML is the safest export for 3-D Plotly figures because it preserves
rotation, zoom, hover labels, and all surfaces:

.. code-block:: python

   from pycsamt.map import write_html

   write_html(fig, "outputs/volume_fence.html")

Static image export is possible when a Plotly image backend such as
Kaleido is installed:

.. code-block:: python

   from pycsamt.map import save_png

   save_png(fig, "outputs/volume_fence.png", width=1600, height=1000)

Troubleshooting
---------------

Empty 3-D figure
   No pseudosection rows could be built.  Check that stations have a
   valid ``Z`` object with frequency, resistivity, and phase arrays.

Only one line appears
   Line metadata may be missing.  Use :func:`pycsamt.map.load_lines`
   with an explicit mapping or ``detect="folder"`` before plotting.

Depth range removes everything
   The pseudo-depth estimate may be outside your requested
   ``depth_range``.  Temporarily remove the range and inspect the full
   extent.

Phase map still responds to ``rho_range``
   This is expected.  Apparent resistivity is still used to estimate
   pseudo-depth and to apply resistivity masks.

Isosurfaces look empty
   ``iso_range`` may not overlap the color-space values.  For
   resistivity with ``log_color=True``, use log10 values in
   ``iso_range``.

Terrain is flat
   Elevations may be missing or non-finite.  The terrain shift uses
   station elevations from the normalized ``StationRecord`` objects.
