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

Guide Sections
--------------

.. grid:: 1 1 2 2
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Containers
      :link: containers
      :link-type: doc
      :img-top: ../../_static/icons/container.svg
      :class-card: pycsamt-card sd-text-center

      Work with ``SiteMixin``, ``Site``, ``Sites``, ``to_sites``, and
      ``to_edis`` wrappers around EDI files and survey collections.

   .. grid-item-card:: Location & Profile
      :link: location_profile
      :link-type: doc
      :img-top: ../../_static/icons/location.svg
      :class-card: pycsamt-card sd-text-center

      Parse coordinates, infer survey profiles, compute distances, bearings,
      chainage, elevation, and line geometry.

   .. grid-item-card:: Selection
      :link: selection
      :link-type: doc
      :img-top: ../../_static/icons/selection.svg
      :class-card: pycsamt-card sd-text-center

      Filter stations by name, index, chainage, frequency coverage,
      coordinates, bounding box, and custom quality predicates.

   .. grid-item-card:: Editing
      :link: editing
      :link-type: doc
      :img-top: ../../_static/icons/edit.svg
      :class-card: pycsamt-card sd-text-center

      Rotate tensors, rename stations, assign coordinates, subset
      frequencies, fill missing arrays, and recompute response products.

   .. grid-item-card:: Recompute
      :link: recompute
      :link-type: doc
      :img-top: ../../_static/icons/recompute.svg
      :class-card: pycsamt-card sd-text-center

      Rewrite EDI files with pyCSAMT, preserve line-folder structure, flatten
      outputs, show progress, and generate manifests.

   .. grid-item-card:: Diagnostics
      :link: computed_diagnostics
      :link-type: doc
      :img-top: ../../_static/icons/diagnostic.svg
      :class-card: pycsamt-card sd-text-center

      Compute strike estimates, apparent resistivity at a frequency, phase
      slope, and tipper magnitude diagnostics.

   .. grid-item-card:: Export & Reporting
      :link: export_reporting
      :link-type: doc
      :img-top: ../../_static/icons/export.svg
      :class-card: pycsamt-card sd-text-center

      Write cleaned sites, batch exports, zip packages, manifests, and
      human-readable site reports.

   .. grid-item-card:: Utilities
      :link: utilities
      :link-type: doc
      :img-top: ../../_static/icons/controls.svg
      :class-card: pycsamt-card sd-text-center

      Use shared helpers for EDI detection, collection coercion, station
      names, coordinates, frequency matching, and angle conversion.

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

.. toctree::
   :maxdepth: 2
   :hidden:

   containers
   location_profile
   selection
   editing
   recompute
   computed_diagnostics
   export_reporting
   utilities
