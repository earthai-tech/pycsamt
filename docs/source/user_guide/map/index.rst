Map Tools
=========

The :mod:`pycsamt.map` package is the code-first mapping layer for
pyCSAMT surveys.  It gives script and notebook users direct control
over station maps, profile pseudosections, 3-D quick-look maps,
overlays, multi-line survey loading, and figure export.

Use this section when you want to work with maps from Python code.  The
web and desktop applications can build on the same rendering layer, but
UI state and callbacks do not live here.

.. warning::

   3-D maps in :mod:`pycsamt.map` are pseudo-depth visualizations
   derived from apparent resistivity and period.  They are useful for
   QC, screening, and communication, but they are not constrained
   inversion models.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   loading
   station
   profile
   volume
   overlays
   export
   mapview
