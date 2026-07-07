.. _site:

Site Tools
==========

The :mod:`pycsamt.site` package provides station-centric tools for working
with EDI sites and survey-line collections. It sits between low-level EDI
parsing and higher-level processing, inversion, plotting, and reporting.

Use this section when you need to:

* wrap one EDI file as a single site object;
* manage many sites as a survey collection;
* inspect station names, coordinates, frequency axes, impedance arrays, phase,
  apparent resistivity, and tipper values;
* select subsets of stations by name, index, frequency coverage, coordinates,
  chainage, or data quality;
* edit station metadata, rotate tensors, subset frequencies, fill missing
  arrays, or recompute derived resistivity and phase;
* normalize foreign or raw EDI files and rewrite them with pyCSAMT while
  preserving survey-line folder structure;
* build survey-line profiles and chainages from coordinates;
* export cleaned sites, write manifests, create zip packages, and generate
  human-readable site reports.

Package Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Page
     - Modules
     - Main purpose
   * - :doc:`containers`
     - ``pycsamt.site.base``
     - ``SiteMixin``, ``Site``, ``Sites``, ``to_sites``, and ``to_edis``
       wrappers around EDI files and collections.
   * - :doc:`location_profile`
     - ``pycsamt.site.location``, ``pycsamt.site.profile``
     - Coordinate parsing, topography, projection, distances, bearing,
       chainage, and survey profile inference.
   * - :doc:`selection`
     - ``pycsamt.site.selection``
     - Station filtering by name, position, chainage, frequency window,
       bounding box, custom predicates, and quality checks.
   * - :doc:`editing`
     - ``pycsamt.site.edit``
     - Rotation, renaming, coordinate assignment, frequency subsetting,
       missing-value filling, and resistivity/phase recomputation.
   * - :doc:`recompute`
     - ``pycsamt.site.recompute``
     - Survey-level EDI recomputation, pyCSAMT EDI rewriting, line-folder
       output, flattening, progress display, and manifests.
   * - :doc:`computed_diagnostics`
     - ``pycsamt.site.compute``
     - Strike estimate, apparent resistivity at a frequency, phase slope, and
       tipper magnitude diagnostics.
   * - :doc:`export_reporting`
     - ``pycsamt.site.export``, ``pycsamt.site.report``
     - Site writing, batch export, zip packaging, manifest files, and terminal
       or tabular reports.
   * - :doc:`utilities`
     - ``pycsamt.site.utils``
     - Shared helpers for EDI detection, collection coercion, names,
       coordinates, frequency matching, and angle conversion.

Typical Workflow
----------------

The site layer is commonly used as an early survey preparation stage.

.. code-block:: python
   :linenos:

   from pycsamt.site import Sites
   from pycsamt.site.selection import by_freq, keep_finite_z
   from pycsamt.site.edit import rotate_all
   from pycsamt.site.profile import Profile
   from pycsamt.site.recompute import recompute_edis

   sites = Sites.from_path("data/edi")
   sites = keep_finite_z(sites)
   sites = by_freq(sites, fmin=1.0, fmax=1000.0)
   sites = rotate_all(sites, angle_deg=30.0)

   profile = Profile.from_sites(sites)
   print(profile.spacing_stats)

For a complete EDI rewrite workflow, use :func:`pycsamt.site.recompute.recompute_edis`:

.. code-block:: python
   :linenos:

   result = recompute_edis(
       "data/willy",
       rotate_angle=30.0,
       template="{source_stem}.edi",
       overwrite=True,
       progress=True,
   )

   print(result.output_root)

Documentation Contents
----------------------

.. toctree::
   :maxdepth: 2

   containers
   location_profile
   selection
   editing
   recompute
   computed_diagnostics
   ../map/index
   export_reporting
   utilities
