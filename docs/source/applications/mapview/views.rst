Views And Controls
==================

MapView organises everything around four views, switched from the left
rail: **Station**, **Profile**, **Pseudosection**, and **3-D**. All four
read the same loaded survey and respect the line picker, so switching
views never changes *what* is selected — only how it is drawn.

Station Map
-----------

The 2-D basemap view of the survey:

* station markers with optional labels, coloured per survey line;
* profile-line overlays connecting the stations of each line;
* an optional Surfer-style contour overlay draped over the basemap;
* the map toolbar offers fit-to-survey, label and profile toggles,
  contour toggle, and basemap selection.

Profile
-------

A single survey line seen along-profile — station positions and
spacing, with the profile unit (metres or kilometres) selectable in the
sidebar. Use it to sanity-check station ordering and gaps before
reading the pseudosection.

Pseudosection
-------------

Apparent-resistivity / phase pseudosections along the selected line:

* *Smooth sections* interpolates between stations for a continuous
  image; leave it off to see the raw station columns;
* *Section resolution* controls the interpolation grid density;
* a log-scale toggle for the colour axis.

3-D View
--------

The full survey in one interactive scene:

* vertical **fence sections** along each survey line;
* **draped topography** (on by default) built from the station
  elevations, with an optional terrain line per profile;
* station markers and labels in the scene;
* contour and smoothing controls shared with the pseudosection;
* grouped accordion controls keep the 3-D-only options (resolution,
  vertical exaggeration, scene styling) in one place.

The 3-D scene is WebGL — see :doc:`troubleshooting` if it renders
blank.

.. seealso::

   :ref:`3-D maps gallery <sphx_glr_examples_map_3d>`
       The same station-map and 3-D figures produced from code, executed
       at build time.
