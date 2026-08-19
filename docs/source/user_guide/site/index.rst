.. _site:

Site Tools
==========

The :mod:`pycsamt.site` package provides station-centric tools for working
with EDI or EMTF-XML sites and survey-line collections. It sits between
low-level EDI/EMTF-XML parsing and higher-level processing, inversion,
plotting, and reporting.

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

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   containers
   metadata
   location_profile
   selection
   editing
   recompute
   computed_diagnostics
   export_reporting
   utilities
