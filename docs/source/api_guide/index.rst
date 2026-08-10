.. _api-guide:

API Guide
=========

:mod:`pycsamt.api` is the package's common entry point: the readers that
turn field data into survey objects, the opt-in tabular view layer behind
``api=True``, the shared mesh-drawing layer behind every resistivity
section and inversion figure, and the ``configure_*``/``reset_*``
pattern that controls every one of those behaviours at runtime — output
directories, plot styles, axis conventions, CLI defaults, and AI-agent
budgets alike. Learning that one pattern once means every family below
already behaves the way you expect.

Use this guide when you need to decide what a dataframe-returning
function gives back, draw or restyle a computational/inversion mesh, or
change any runtime setting — from a single figure's DPI to the station
ordering used package-wide — without hunting through individual function
signatures. For the complete generated reference of every
``pycsamt.api`` function and class, see :doc:`../api/api`.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   views
   mesh
   style
   contour
   interpretation
   agent_config
   section
   station_rendering
   view_controls
   configuration
